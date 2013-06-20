/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file BamLoader.hh
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_BAM_LOADER_HH
#define iSAAC_IO_BAM_LOADER_HH

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "bgzf/BgzfReader.hh"
#include "flowcell/ReadMetadata.hh"
#include "bam/BamParser.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace io
{

//#pragma GCC push_options
//#pragma GCC optimize ("0")

struct BamLoaderException : common::IoException
{
    BamLoaderException(const std::string &message) : common::IoException(EINVAL, message){}
};

class BamLoader
{
    std::string bamPath_;
    bgzf::ParallelBgzfReader bgzfReader_;

    // previous pass content of currentPassBam_. Assumed that processor is interested to
    // be able to access both buffers at the same time in order to deal with pairs (or reads) that
    // overlap the buffer boundary
    std::vector<char> lastPassBam_;
    unsigned lastUnparsedBytes_;
    boost::array<std::vector<char>, 2> decompressionBuffers_;
    boost::array<unsigned, 2> unparsedBytes_;

    bam::BamParser bamParser_;

    common::ThreadVector decompressParseParallelizationThreads_;
    boost::mutex mutex_;
    boost::condition_variable stateChangedCondition_;
    unsigned nextDecompressorThread_;
    unsigned nextParserThread_;

public:
    BamLoader(
        std::size_t maxPathLength,
        common::ThreadVector &threads,
        const unsigned coresMax);

    BamLoader(
        std::size_t maxPathLength,
        common::ThreadVector &threads);

    void open(const boost::filesystem::path &bamPath)
    {
        if (bamPath_ != bamPath.string())
        {
            bamPath_ = bamPath.c_str();// ensure actual copying, prevent path buffer sharing
            lastUnparsedBytes_ = 0;
            unparsedBytes_[0] = 0;
            unparsedBytes_[1] = 0;
            nextDecompressorThread_ = 0;
            nextParserThread_ = 0;
            bgzfReader_.open(bamPath);
            bamParser_.reset();
        }
        else
        {
            ISAAC_THREAD_CERR << "Keeping bam stream open for " << bamPath_ << std::endl;
        }
    }

    template <typename ProcessorT>
    void load(ProcessorT processor);

private:
    static const unsigned UNPARSED_BYTES_MAX = 1024*100;

    template <typename RecordProcessorT>
    void parallelLoad(const unsigned threadNumber,
        bool &wantMoreData, bool &exception, RecordProcessorT processor);
    template <typename RecordProcessorRemoveOldT>
    void swapBuffers(RecordProcessorRemoveOldT removeOld, std::vector<char> &current, unsigned &unparsedBytes);
    bool waitForLoadSlot(boost::unique_lock<boost::mutex> &lock, bool &exception, const unsigned threadNumber);
    void returnLoadSlot();
    bool waitForParseSlot(boost::unique_lock<boost::mutex> &lock, bool &wantMoreData, bool &exception, const unsigned threadNumber);
    void returnParseSlot(bool wantMoreData);
    void reserveBuffers(std::size_t maxPathLength);
};


/**
 * \brief Exchanges last and current buffers. Resets current buffer, removes all index entries that
 *        have been processed by the previous pass
 */
template <typename RecordProcessorRemoveOldT>
void BamLoader::swapBuffers(RecordProcessorRemoveOldT removeOld, std::vector<char> &current, unsigned &unparsedBytes)
{
    removeOld(lastPassBam_.begin().base(), lastPassBam_.end().base());
    lastPassBam_.swap(current);
    lastUnparsedBytes_ = unparsedBytes;
    unparsedBytes = 0;
    current.clear();
}

inline bool BamLoader::waitForLoadSlot(boost::unique_lock<boost::mutex> &lock, bool &exception, const unsigned threadNumber)
{
    while(!exception && nextDecompressorThread_ != threadNumber)
    {
        stateChangedCondition_.wait(lock);
    }
    return !exception;
}

inline void BamLoader::returnLoadSlot()
{
    nextDecompressorThread_ = (nextDecompressorThread_ + 1) % decompressParseParallelizationThreads_.size();
    stateChangedCondition_.notify_all();
}

inline bool BamLoader::waitForParseSlot(boost::unique_lock<boost::mutex> &lock, bool &wantMoreData, bool &exception, const unsigned threadNumber)
{
    while(!exception && wantMoreData && nextParserThread_ != threadNumber)
    {
        stateChangedCondition_.wait(lock);
    }

    return !exception && wantMoreData;
}

inline void BamLoader::returnParseSlot(bool suspending)
{
    if (!suspending)
    {
        nextParserThread_ = (nextParserThread_ + 1) % decompressParseParallelizationThreads_.size();
    }
    else
    {
        ISAAC_THREAD_CERR << "BamLoader::returnParseSlot nextParserThread_:" << nextParserThread_ << std::endl;
    }
    stateChangedCondition_.notify_all();
    // otherwise, the thread will need to repeat parsing as some of the data did not fit into processor
}

template <typename RecordProcessorT>
void BamLoader::parallelLoad(
    const unsigned threadNumber,
    bool &wantMoreData, bool &exception,
    RecordProcessorT processor)
{
    boost::unique_lock<boost::mutex> lock(mutex_);

    try
    {
        while (!exception && wantMoreData)
        {
            if(decompressionBuffers_[threadNumber].empty())
            {
                if (!waitForLoadSlot(lock, exception, threadNumber))
                {
                    break;
                }
                ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&BamLoader::returnLoadSlot, this))
                {
                    decompressionBuffers_[threadNumber].resize(UNPARSED_BYTES_MAX);
                    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                    if (!bgzfReader_.readMoreData(decompressionBuffers_[threadNumber]))
                    {
                        decompressionBuffers_[threadNumber].clear();
                        ISAAC_THREAD_CERR << "no more data on thread " << threadNumber << std::endl;
                        return;
                    }
//                    ISAAC_THREAD_CERR << "read " << decompressionBuffers_[threadNumber].size() << " bytes, done:" << std::endl;
                }
            }
            else
            {
                // This means we've filled decompressionBuffer_ during last parallelLoad call, but the parser was not able to finish
                // processing currentPassBam_. Just wait until it is done
                ISAAC_THREAD_CERR << "Thread " << threadNumber << " already has " << decompressionBuffers_[threadNumber].size() << " bytes" << std::endl;
            }

            if (!waitForParseSlot(lock, wantMoreData, exception, threadNumber))
            {
                break;
            }

            bool suspending = false;
            ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&BamLoader::returnParseSlot, this, boost::ref(suspending)))
            {
                if (!decompressionBuffers_[threadNumber].empty())
                {
                    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                    std::vector<char>::const_iterator unparsedBegin;
                    if (unparsedBytes_[threadNumber])
                    {
                        ISAAC_ASSERT_MSG(!lastUnparsedBytes_, "lastUnparsedBytes_ must be 0 on resume. Got: " << lastUnparsedBytes_);
                        if (-1U == unparsedBytes_[threadNumber])
                        {
                            ISAAC_THREAD_CERR << "force-resuming on thread " << threadNumber << " unparsed: " << unparsedBytes_[threadNumber] << std::endl;
                            unparsedBegin = decompressionBuffers_[threadNumber].end();
                        }
                        else
                        {
                            unparsedBegin = decompressionBuffers_[threadNumber].end() - unparsedBytes_[threadNumber];
                            ISAAC_THREAD_CERR << "resuming on thread " << threadNumber << " unparsed: " << unparsedBytes_[threadNumber] << std::endl;
                        }
                    }
                    else
                    {
                        ISAAC_ASSERT_MSG(lastUnparsedBytes_ <= UNPARSED_BYTES_MAX, "Too much unparsed from previous step: " << lastUnparsedBytes_);
                        std::copy(lastPassBam_.end() - lastUnparsedBytes_, lastPassBam_.end(), decompressionBuffers_[threadNumber].begin() + UNPARSED_BYTES_MAX - lastUnparsedBytes_);
                        unparsedBegin = decompressionBuffers_[threadNumber].begin() + UNPARSED_BYTES_MAX - lastUnparsedBytes_;
                        lastUnparsedBytes_ = 0;
                    }
                    wantMoreData = bamParser_.parse(unparsedBegin, decompressionBuffers_[threadNumber].end(), boost::get<0>(processor));
                    // TODO: here the assumption is that unparsedBytes_ will point at the last bam block that did not end in the buffer.
                    // Currently processor does not stop in the middle of indexing, but if it did that, we'd be forced to
                    // potentially copy more than UNPARSED_BYTES_MAX of data into the next buffer.
                    unparsedBytes_[threadNumber] = std::distance<std::vector<char>::const_iterator>(unparsedBegin, decompressionBuffers_[threadNumber].end());

                    if (wantMoreData)
                    {
                        swapBuffers(boost::get<1>(processor), decompressionBuffers_[threadNumber], unparsedBytes_[threadNumber]);
                        ISAAC_ASSERT_MSG(lastUnparsedBytes_ <= UNPARSED_BYTES_MAX, "Too much unparsed from this step: " << lastUnparsedBytes_);
                        suspending = false;
                    }
                    else
                    {
                        if(!unparsedBytes_[threadNumber])
                        {
                            ISAAC_THREAD_CERR << "force-suspending on thread " << threadNumber << " wantMoreData: " << wantMoreData << " lastUnparsedBytes_:" << lastUnparsedBytes_ << std::endl;
                            unparsedBytes_[threadNumber] = -1U;
                        }
                        else
                        {
                            ISAAC_THREAD_CERR << "suspending on thread " << threadNumber << " unparsed: " << unparsedBytes_[threadNumber] << " lastUnparsedBytes_:" << lastUnparsedBytes_ << std::endl;
                        }
                        suspending = true;
                    }
                }
                else
                {
                    ISAAC_ASSERT_MSG(bgzfReader_.isEof(), "Expected end of compressed bam data stream");
                    if (lastUnparsedBytes_)
                    {
                        BOOST_THROW_EXCEPTION(BamLoaderException(
                            (boost::format("Reached the end of the bam file with %d bytes unparsed. Truncated Bam?") % lastUnparsedBytes_).str()));
                    }
                    // no more data to process
                    // ensure processor has a chance to deal with the last batch of blocks
                    swapBuffers(boost::get<1>(processor), decompressionBuffers_[threadNumber], unparsedBytes_[threadNumber]);
                }
            }
        }
    }
    catch(...)
    {
        exception = true;
        stateChangedCondition_.notify_all();
        throw;
    }
}

template <typename ProcessorT>
void BamLoader::load(ProcessorT processor)
{
    bool wantMoreData = true;
    bool exception = false;

    decompressParseParallelizationThreads_.execute(
        boost::bind(&BamLoader::parallelLoad<ProcessorT>, this,
                    _1, boost::ref(wantMoreData), boost::ref(exception),
                    processor));

}
//#pragma GCC pop_options


} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_BAM_LOADER_HH
