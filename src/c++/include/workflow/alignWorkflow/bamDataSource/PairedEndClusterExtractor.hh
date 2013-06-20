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
 ** \file PairedEndClusterExtractor.hh
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_PAIRED_END_CLUSTER_EXTRACTOR_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_PAIRED_END_CLUSTER_EXTRACTOR_HH

#include <cmath>

#include <boost/crc.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/back.hpp>

#include "bam/Bam.hh"
#include "bam/BamParser.hh"
#include "common/FastIo.hh"
#include "common/FileSystem.hh"
#include "flowcell/ReadMetadata.hh"
#include "io/FileBufCache.hh"
#include "reference/ReferencePosition.hh"

//#pragma GCC push_options
//#pragma GCC optimize ("0")

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{
namespace bamDataSource
{
struct ClusterExtractorException : common::IoException
{
    ClusterExtractorException(const std::string &message) : common::IoException(EINVAL, message){}
};

struct IndexRecord
{
    typedef unsigned long NameHashType;
    NameHashType nameHash_;
    const bam::BamBlockHeader * bamRecordPointer_;
    static const NameHashType EXTRACTED = -1UL;

    IndexRecord() : nameHash_(0), bamRecordPointer_(0)
    {
    }

    IndexRecord(
        const bam::BamBlockHeader &block) : bamRecordPointer_(&block)
    {
        ISAAC_ASSERT_MSG(sizeof(nameHash_) <= block.getReadNameLength(), "Name too short. Make the code that copes with it.");
        // last few bytes of names are usually quite differernt.
        const char *readNameEnd = block.read_name + block.getReadNameLength() - sizeof(nameHash_);
        nameHash_ = *reinterpret_cast<const NameHashType*>(readNameEnd);
        if (EXTRACTED == nameHash_)
        {
            nameHash_ = 0;
        }
    }

    const bam::BamBlockHeader &getBlock() const {return *bamRecordPointer_;}
    void markExtracted() {nameHash_ = EXTRACTED;}
    bool isExtracted() const {return EXTRACTED == nameHash_;}

    friend std::ostream& operator << (std::ostream &os, const IndexRecord &idx)
    {
        return os << "BamBufferIndexRecord(" << idx.nameHash_ << ", " << idx.getBlock() << ")";
    }

    friend bool operator <(const IndexRecord &left, const IndexRecord &right)
    {
        // comparing by name is enough, comparing by r1 position is much faster
        if (left.nameHash_ < right.nameHash_)
        {
            return true;
        }
        else if (left.nameHash_ == right.nameHash_)
        {
            const bam::BamBlockHeader &leftBlock = left.getBlock();
            const bam::BamBlockHeader &rightBlock = right.getBlock();
            const int namecmp = strcmp(leftBlock.read_name, rightBlock.read_name);
            if (0 > namecmp)
            {
                return true;
            }
            else if (0 == namecmp)
            {
                // put read one on top so that it's simpler to extract data
                return leftBlock.isReadOne();
            }
        }
        return false;
    }
};

class TempFileClusterExtractor : std::vector<char>
{
    static const bool allowUnpairedReads_ = true;
    std::string tempFilePath_;
    io::FileBufWithReopen unpairedReadsFile_;

    std::vector<std::vector<char>::const_iterator> recordIndex_;
    std::vector<std::vector<char>::const_iterator>::const_iterator firstUnextracted_;

public:
    typedef unsigned char FlagsType;
    static const FlagsType READ_ONE_FLAG = 1;
    static const FlagsType PASS_FILTER_FLAG = 2;

    TempFileClusterExtractor(
        const std::size_t maxTempFilePathLength,
        const std::size_t bufferSize,
        const std::size_t minClusterLength) :
            unpairedReadsFile_(std::ios_base::in | std::ios_base::binary),
            firstUnextracted_(recordIndex_.end())
    {
        reserve(bufferSize);
        // ignore the fact that each cluster will have more than just the bcl data
        recordIndex_.reserve(bufferSize / minClusterLength);
        tempFilePath_.reserve(maxTempFilePathLength);
    }

    bool isEmpty() const {return recordIndex_.end() == firstUnextracted_;}

    void open(const boost::filesystem::path &tempFilePath, std::streamsize expectedFileSize);

    template <typename ClusterInsertIt, typename PfInsertIt>
    unsigned extractClusters(
        const unsigned r1Length,
        const unsigned r2Length,
        unsigned clusterCount,
        ClusterInsertIt &clusterIt,
        PfInsertIt &pfIt)
    {
        ISAAC_THREAD_CERR << "TempFileClusterExtractor::extractClusters: " << clusterCount << std::endl;
        while (recordIndex_.end() != firstUnextracted_ && clusterCount)
        {
            std::vector<std::vector<char>::const_iterator>::const_iterator prev = firstUnextracted_++;
            if (recordIndex_.end() == firstUnextracted_ || strcmp(getReadName(*prev), getReadName(*firstUnextracted_)))
            {
                if (!allowUnpairedReads_)
                {
                    BOOST_THROW_EXCEPTION(ClusterExtractorException(
                        (boost::format("No pair for read name %s in %s") % getReadName(*prev) % tempFilePath_).str()));
                }
                else
                {
                    if (isReadOne(*prev))
                    {
                        clusterIt = std::copy(getBclBegin(*prev), getBclEnd(*prev), clusterIt);
                    }
                    const unsigned toFill = isReadOne(*prev) ? r2Length : r1Length;
                    clusterIt = std::fill_n(clusterIt, toFill, 0);
                    if (!isReadOne(*prev))
                    {
                        clusterIt = std::copy(getBclBegin(*prev), getBclEnd(*prev), clusterIt);
                    }
                    *pfIt++ = isPf(*prev);
                }
            }
            else
            {
                ISAAC_ASSERT_MSG(!isReadOne(*firstUnextracted_), "Out of two reads, second one was expected to be read 2 " <<
                                 getReadName(*prev) << ":" << isReadOne(*prev) <<
                                 " " << getReadName(*firstUnextracted_) << ":" << isReadOne(*firstUnextracted_));
                ISAAC_ASSERT_MSG(isReadOne(*prev), "Out of two reads, first one was expected to be read 1 " <<
                                 getReadName(*prev) << ":" << isReadOne(*prev) <<
                                 " " << getReadName(*firstUnextracted_) << ":" << isReadOne(*firstUnextracted_));
                clusterIt = std::copy(getBclBegin(*prev), getBclEnd(*prev), clusterIt);
                ISAAC_ASSERT_MSG(isPf(*prev) == isPf(*firstUnextracted_), "Pf flag must be the same for both reads of the cluster " << getReadName(*prev));
                *pfIt++ = isPf(*prev);
                clusterIt = std::copy(getBclBegin(*firstUnextracted_), getBclEnd(*firstUnextracted_), clusterIt);
                if (recordIndex_.end() != firstUnextracted_)
                {
                    ++firstUnextracted_;
                }
            }

            --clusterCount;
        }
        ISAAC_THREAD_CERR << "TempFileClusterExtractor::extractClusters done: " << clusterCount << std::endl;
        return clusterCount;
    }

private:
    static const std::vector<char>::const_iterator getNextRecord(const std::vector<char>::const_iterator it)
        {return it + reinterpret_cast<const unsigned &>(*it);}

    static const unsigned char *getBclBegin(const std::vector<char>::const_iterator it)
        {return reinterpret_cast<const unsigned char *>(std::find(getReadName(it), getNextRecord(it).base(), 0) + 1);}

    static const unsigned char *getBclEnd(const std::vector<char>::const_iterator it)
        {return reinterpret_cast<const unsigned char *>(getNextRecord(it).base());}

    static std::size_t getReadLength(const std::vector<char>::const_iterator it)
        {return std::distance(getBclBegin(it), getBclEnd(it));}

    static const char *getReadName(const std::vector<char>::const_iterator it)
    {
        // name is a null-terminated string following the unsigned record length and a bool flag
        return it.base() + sizeof(unsigned) + sizeof(bool);
    }

    static bool isReadOne(const std::vector<char>::const_iterator it)
        {return getFlags(it) & READ_ONE_FLAG;}

    static bool isPf(const std::vector<char>::const_iterator it)
        {return getFlags(it) & PASS_FILTER_FLAG;}

    static unsigned char getFlags(const std::vector<char>::const_iterator it)
        {return *reinterpret_cast<const unsigned char*>(it.base() + sizeof(unsigned));}

    static bool compareNameAndRead(const std::vector<char>::const_iterator left, const std::vector<char>::const_iterator right)
    {
        const int namecmp = strcmp(getReadName(left), getReadName(right));
        return 0 > namecmp || (0 == namecmp && isReadOne(left) > isReadOne(right));
    }
};

class UnpairedReadsCache
{
    template <unsigned requiredWidth>
    struct CrcSelector
    {
        template <unsigned minVal, unsigned maxVal, unsigned val>
        struct MinMax
        {
            BOOST_STATIC_CONSTANT(unsigned, value
             = (val < minVal) ? minVal : (val > maxVal) ? maxVal : val);
        };

        typedef boost::mpl::vector<
            boost::crc_optimal<5, 0x09>,
            boost::crc_optimal<6, 0x03>,
            boost::crc_optimal<7, 0x09> > AllowedCrcs;

        typedef typename boost::mpl::front<AllowedCrcs>::type FrontCrc;
        typedef typename boost::mpl::back<AllowedCrcs>::type BackCrc;

        typedef typename boost::mpl::at_c<AllowedCrcs, (MinMax<FrontCrc::bit_count, BackCrc::bit_count, requiredWidth>::value - FrontCrc::bit_count)  >::type Selection;

        static unsigned getStringCrc(const char *string, const std::size_t stringLength)
        {
            Selection crc;
            crc.process_bytes(string, stringLength);
            return crc.checksum();
        }
    };

    const unsigned crcWidth_;
    const boost::filesystem::path &tempDirectoryPath_;
    const std::string directorySeparator_;
    std::vector<std::string> tempFilePaths_;
    std::vector<std::size_t> tempFileSizes_;
    std::vector<std::string>::const_iterator extractorFileIterator_;
    bool extracting_;
    std::vector<io::FileBufHolder<io::FileBufWithReopen> > tempFiles_;
    std::string tempFilePathBuffer_;
    TempFileClusterExtractor extractor_;
    // aim to have ~3 gigabyte temp files assuming none of the input reads pair
    static const std::size_t UNPAIRED_BUFFER_SIZE = 1024UL * 1024UL * 1024UL * 3UL;
public:
    UnpairedReadsCache(
        const boost::filesystem::path &tempDirectoryPath,
        const std::size_t maxBamFileSize,
        const std::size_t maxFlowcellIdLength,
        const std::size_t minClusterLength) :
            crcWidth_(log2(std::max<std::size_t>(1, maxBamFileSize / UNPAIRED_BUFFER_SIZE))),
            tempDirectoryPath_(tempDirectoryPath),
            directorySeparator_(common::getDirectorySeparator()),
            tempFilePaths_(1 << getEffectiveCrcWidth<7>(crcWidth_)),
            tempFileSizes_(tempFilePaths_.size(), 0),
            extractorFileIterator_(tempFilePaths_.end()),
            extracting_(false),
            tempFiles_(
                1 << getEffectiveCrcWidth<7>(crcWidth_),
                io::FileBufHolder<io::FileBufWithReopen>(std::ios_base::out | std::ios_base::app | std::ios_base::binary,
                                                         getMaxTempFilePathLength(maxFlowcellIdLength))),
            extractor_(getMaxTempFilePathLength(maxFlowcellIdLength), UNPAIRED_BUFFER_SIZE, minClusterLength)
    {
        tempFilePathBuffer_.reserve(getMaxTempFilePathLength(maxFlowcellIdLength));
        BOOST_FOREACH(std::string &tempPath, tempFilePaths_)
        {
            tempPath.reserve(tempFilePathBuffer_.capacity());
        }
    }

    bool extractingUnpaired() const {return extracting_;}

    void open(const std::string &flowcellId)
    {
        unsigned i = 0;
        BOOST_FOREACH(std::string &tempPath, tempFilePaths_)
        {
            tempPath = makeTempFilePath(flowcellId, i);
            unlink(tempPath.c_str());
            tempFiles_[i].reopen(tempPath.c_str(), io::FileBufWithReopen::SequentialOnce);
            tempFileSizes_[i] = 0;
            ++i;
        }
        extracting_ = false;
    }

    void startExtractingUnpaired(const std::string &flowcellId)
    {
        ISAAC_THREAD_CERR << "startExtractingUnpaired " << std::endl;
        std::for_each(tempFiles_.begin(), tempFiles_.end(), boost::bind(&io::FileBufHolder<io::FileBufWithReopen>::flush, _1));

        extractorFileIterator_ = tempFilePaths_.begin();
        extractor_.open(
            *extractorFileIterator_, tempFileSizes_[extractorFileIterator_ - tempFilePaths_.begin()]);
        extracting_ = true;
    }

    template <typename IteratorT>
    void storeUnpaired(
        IteratorT unpairedBegin,
        IteratorT unpairedEnd,
        const flowcell::ReadMetadataList &readMetadataList);

    template <typename ClusterInsertIt, typename PfInsertIt>
    unsigned extractClusters(
        const unsigned r1Length,
        const unsigned r2Length,
        unsigned clusterCount,
        ClusterInsertIt &clusterIt,
        PfInsertIt &pfIt)
    {
        while (tempFilePaths_.end() != extractorFileIterator_ && clusterCount)
        {
            clusterCount = extractor_.extractClusters(r1Length, r2Length, clusterCount, clusterIt, pfIt);
            if (clusterCount)
            {
                ++extractorFileIterator_;
                if (tempFilePaths_.end() != extractorFileIterator_)
                {
                    extractor_.open(
                        *extractorFileIterator_,
                        tempFileSizes_[extractorFileIterator_ - tempFilePaths_.begin()]);
                }
            }
        }
        return clusterCount;
    }


private:
    const char* makeTempFilePath(const std::string &flowcellId, unsigned crc)
    {
        return makeTempFilePath(flowcellId, crc, tempFilePathBuffer_).c_str();
    }

    const std::string& makeTempFilePath(const std::string &flowcellId, unsigned crc, std::string &buffer)
    {
        buffer = tempDirectoryPath_.c_str();
        buffer += directorySeparator_;
        buffer += flowcellId;
        buffer += "-unpaired-";
        common::appendUnsignedInteger(buffer, crc);
        buffer += ".tmp";
        return buffer;
    }

    std::size_t getMaxTempFilePathLength(const std::size_t maxFlowcellIdLength)
    {
        std::string buffer;
        return makeTempFilePath(std::string(maxFlowcellIdLength, 'a'), 9999, buffer).size();
    }

    template <unsigned N>
    static unsigned getNameCrc(unsigned crcWidth, const char *name, const std::size_t nameLength)
    {
        ISAAC_ASSERT_MSG(crcWidth <= N, "Requested crcWidth is too big: " << crcWidth << " max: " << N);
        if (N == crcWidth)
        {
            return CrcSelector<N>::getStringCrc(name, nameLength);
        }

        return getNameCrc<N-1>(crcWidth, name, nameLength);
    }

    template <unsigned N>
    static unsigned getEffectiveCrcWidth(unsigned requiredWidth)
    {
        ISAAC_ASSERT_MSG(requiredWidth <= N, "Requested crcWidth is too big: " << requiredWidth << " max: " << N);
        if (N == requiredWidth)
        {
            return CrcSelector<N>::Selection::bit_count;
        }

        return getEffectiveCrcWidth<N-1>(requiredWidth);
    }

};

template <>
inline unsigned UnpairedReadsCache::getNameCrc<0>(unsigned crcWidth, const char *name, const std::size_t nameLength)
{
    return CrcSelector<0>::getStringCrc(name, nameLength);
}

template <>
inline unsigned UnpairedReadsCache::getEffectiveCrcWidth<0>(unsigned requiredWidth)
{
    return CrcSelector<0>::Selection::bit_count;
}


class PairedEndClusterExtractor :
    // each bgzf buffer is unlikely to expand to over 65535 bytes. At least ones produced by iSAAC will not.
    // make a crazy assumption that each bam block is 1 byte long, two buffers...
    common::FiniteCapacityVector<IndexRecord, 65535*2>
{
    typedef common::FiniteCapacityVector<IndexRecord, 65535*2> BaseT;
    iterator firstUnextracted_;

    std::string flowcellId_;

    UnpairedReadsCache unpairedReadCache_;
public:
    using BaseT::size;
    PairedEndClusterExtractor(
        const boost::filesystem::path &tempDirectoryPath,
        const std::size_t maxBamFileLength,
        const std::size_t maxFlowcellIdLength,
        const std::size_t minClusterLength) :
            firstUnextracted_(end()),
//            unpairedBegin_(begin()),
//            unpairedEnd_(begin()),
            unpairedReadCache_(
                tempDirectoryPath,
                maxBamFileLength,
                maxFlowcellIdLength,
                minClusterLength)
    {
    }

    void open(const std::string &flowcellId)
    {
        if (flowcellId_ != flowcellId)
        {
            reset();
            unpairedReadCache_.open(flowcellId);
            flowcellId_ = flowcellId;
        }
    }

    bool extractingUnpaired() const {return unpairedReadCache_.extractingUnpaired();}
    bool isEmpty() const {return end() == firstUnextracted_;}


    template <typename ClusterInsertIt, typename PfInserIt>
    bool append(
        const bam::BamBlockHeader &block, const bool lastBlock,
        unsigned &clusterCount,
        const flowcell::ReadMetadataList &readMetadataList,
        ClusterInsertIt &clustersIt,
        PfInserIt &pfIt)
    {
        push_back(IndexRecord(block));

        if (lastBlock)
        {
    //  ISAAC_THREAD_CERR << "parsed " << index_.size() << " entries, unparsedBytes_:" << unparsedBytes_ << std::endl;
            sort();
    //  ISAAC_THREAD_CERR << "sorted " << index_.size() << " entries" << std::endl;
            //                     unsigned before = clusterCount;
            clusterCount = extractPairedReads(clusterCount, clustersIt, pfIt, readMetadataList);
    //                    ISAAC_THREAD_CERR << "extracted " << before - clusterCount << " clusters" << std::endl;
        }

        // if clusterCount is non-zero, we can extract more clusters into the clustersIt. Say "we want more data"
        return !!clusterCount;
    }

    void removeOld(
        const char* rangeStart,
        const char* rangeEnd,
        const flowcell::ReadMetadataList &readMetadataList);

    template <typename ClusterInsertIt, typename PfInsertIt>
    unsigned extractUnpaired(
        const unsigned r1Length,
        const unsigned r2Length,
        unsigned clusterCount,
        ClusterInsertIt &clusterIt,
        PfInsertIt &pfIt)
    {
        return unpairedReadCache_.extractClusters(r1Length, r2Length, clusterCount, clusterIt, pfIt);
    }


    /**
     * \brief For index entries that have all the reads needed, copies bcl data into the output and removes the index entry
     *
     * \return number of clusters extracted
     */
    template <typename ClusterAccessIt, typename PfInsertIt>
    unsigned extractPairedReads(
        unsigned clusterCount,
        ClusterAccessIt &clusterIt,
        PfInsertIt &pfIt,
        const flowcell::ReadMetadataList &readMetadataList)
    {
        if (2 > std::distance(firstUnextracted_, end()))
        {
            return clusterCount;
        }
        iterator it = firstUnextracted_ + 1;
        for (; clusterCount && end() != it;)
        {
    //            ISAAC_THREAD_CERR << *(it-1) << " " << clusterCount << std::endl;
    //            ISAAC_THREAD_CERR << *it << " " << clusterCount << std::endl;
            if (readNamesMatch(*(it - 1), *it))
            {
                const bam::BamBlockHeader &r1Block = (it - 1)->getBlock();
                const bam::BamBlockHeader &r2Block = it->getBlock();
                ISAAC_ASSERT_MSG(r1Block.isReadOne(), "Sort order must put r1 block before r2 block " << r1Block  << " " << r2Block);
                ISAAC_ASSERT_MSG(!r2Block.isReadOne(), "Sort order must put r1 block before r2 block " << r1Block  << " " << r2Block);
                ISAAC_ASSERT_MSG(r2Block.isPf() == r1Block.isPf(), "Pf flag must be the same for both reads of the cluster " << r2Block << " " << r1Block);

                const flowcell::ReadMetadata &r1Metadata = readMetadataList.at(!r1Block.isReadOne());
                clusterIt = extractBcl(r1Block, clusterIt, r1Metadata);
                (it - 1)->markExtracted();

                const flowcell::ReadMetadata &r2Metadata = readMetadataList.at(!r2Block.isReadOne());
                clusterIt = extractBcl(r2Block, clusterIt, r2Metadata);
                it->markExtracted();

                *pfIt++ = r2Block.isPf();

                --clusterCount;
                ++it;
                if (end() != it)
                {
                    ++it;
                    if (end() == it)
                    {
                        // this is the case when the last line in the index has its mate in the next buffer and this is
                        // the first unpaired entry in the current buffer.
                        break;
                    }
                }
            }
            else
            {
    //                ISAAC_THREAD_CERR << "unpaired " << " " << *it << std::endl;
                ++it;
            }
        }
        if (clusterCount)
        {
            erase(std::remove_if(begin(), end(), boost::bind(&IndexRecord::isExtracted, _1)), end());
            reset();
        }
        else
        {
            firstUnextracted_ = it - 1;
            ISAAC_THREAD_CERR << "Out of clusters: " <<
                std::distance(begin(), firstUnextracted_) << " " <<
                std::distance(begin(), it) << " " <<
                std::distance(it, end()) << " " <<
                std::endl;
        }
        return clusterCount;
    }

    void startExtractingUnpaired()
    {
        unpairedReadCache_.startExtractingUnpaired(flowcellId_);
    }

private:
    void sort();

    static bool readNamesMatch(const IndexRecord &left, const IndexRecord &right)
    {
        const bam::BamBlockHeader &leftBlock = left.getBlock();
        const bam::BamBlockHeader &rightBlock = right.getBlock();

        return (leftBlock.getReadNameLength() == rightBlock.getReadNameLength() &&
            leftBlock.read_name + leftBlock.getReadNameLength() ==
                std::mismatch(leftBlock.read_name, leftBlock.read_name + leftBlock.getReadNameLength(), rightBlock.read_name).first);
    }

    void storeUnpaired(
        BaseT::const_iterator unpairedBegin,
        BaseT::const_iterator unpairedEnd,
        const flowcell::ReadMetadataList &readMetadataList)
    {
        unpairedReadCache_.storeUnpaired(unpairedBegin, unpairedEnd, readMetadataList);
    }
    void reset();
};

} // namespace bamDataSource
} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

//#pragma GCC pop_options


#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_PAIRED_END_CLUSTER_EXTRACTOR_HH