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
 ** \file PairedEndClusterExtractor.cpp
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/integer/static_min_max.hpp>

#include "common/Debug.hh"
#include "workflow/alignWorkflow/bamDataSource/PairedEndClusterExtractor.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{
namespace bamDataSource
{

void TempFileClusterExtractor::open(const boost::filesystem::path &tempFilePath, std::streamsize expectedFileSize)
{
    if (tempFilePath_ != tempFilePath.string())
    {
        recordIndex_.clear();
        unpairedReadsFile_.reopen(tempFilePath.c_str(), io::FileBufWithReopen::SequentialOnce);
        std::istream is(&unpairedReadsFile_);
        if (!is.eof())
        {
            if (!is)
            {
                BOOST_THROW_EXCEPTION(isaac::common::IoException(
                    errno, (boost::format("Unable to open file for reading %s") % tempFilePath.string()).str()));
            }

            resize(expectedFileSize + 1);
            if (!is.read(&front(), size()) && !is.eof())
            {
                BOOST_THROW_EXCEPTION(isaac::common::IoException(
                    errno, (boost::format("Unable to read %d bytes from %s") % size() % tempFilePath.string()).str()));
            }

            ISAAC_ASSERT_MSG(is.gcount() == expectedFileSize, "Read mismatching number of bytes from the file: " << expectedFileSize <<
                " from " << tempFilePath);

            resize(expectedFileSize);

            for (std::vector<char>::const_iterator it = begin(); end() != it; it = getNextRecord(it))
            {
                recordIndex_.push_back(it);
            }

            std::sort(recordIndex_.begin(), recordIndex_.end(), compareNameAndRead);
            ISAAC_THREAD_CERR << "TempFileClusterExtractor::open: " << recordIndex_.size() << " " << tempFilePath << std::endl;
        }
        else
        {
            ISAAC_THREAD_CERR << "TempFileClusterExtractor::open: empty " << tempFilePath << std::endl;
        }
    }
    firstUnextracted_ = recordIndex_.begin();
    tempFilePath_ = tempFilePath.c_str();
}

template <typename IteratorT>
void UnpairedReadsCache::storeUnpaired(
    IteratorT unpairedBegin,
    IteratorT unpairedEnd,
    const flowcell::ReadMetadataList &readMetadataList)
{
    common::FiniteCapacityVector<char, 10240> bcl;
    BOOST_FOREACH(const IndexRecord &idx, std::make_pair(unpairedBegin, unpairedEnd))
    {
        const flowcell::ReadMetadata &readMetadata = readMetadataList.at(!idx.getBlock().isReadOne());
        const unsigned effectiveReadLength = bam::getEffectiveReadLength(idx.getBlock(), readMetadata);

        const bam::BamBlockHeader &block = idx.getBlock();

        const unsigned nameCrc = getNameCrc<7>(crcWidth_, block.read_name, block.getReadNameLength());
        std::ostream os(tempFiles_[nameCrc].get());

        const unsigned nameLength = block.getReadNameLength();
        const unsigned recordLength = sizeof(recordLength) + sizeof(bool) + nameLength + effectiveReadLength;
        if (!os.write(reinterpret_cast<const char*>(&recordLength), sizeof(unsigned)))
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(
                errno, (boost::format("Failed to write: %d bytes into %s") % sizeof(unsigned) % tempFilePaths_[nameCrc]).str()));
        }
        tempFileSizes_[nameCrc] += sizeof(unsigned);

        const TempFileClusterExtractor::FlagsType flags =
            (block.isReadOne() ? TempFileClusterExtractor::READ_ONE_FLAG : 0) |
                (block.isPf() ? TempFileClusterExtractor::PASS_FILTER_FLAG : 0) ;
        if (!os.write(reinterpret_cast<const char*>(&flags), sizeof(flags)))
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(
                errno, (boost::format("Failed to write: %d bytes into %s") % nameLength % tempFilePaths_[nameCrc]).str()));
        }
        tempFileSizes_[nameCrc] += sizeof(flags);

        if (!os.write(block.read_name, nameLength))
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(
                errno, (boost::format("Failed to write: %d bytes into %s") % nameLength % tempFilePaths_[nameCrc]).str()));
        }
        tempFileSizes_[nameCrc] += nameLength;

        bcl.resize(effectiveReadLength);
        bam::extractBcl(idx.getBlock(), bcl.begin(), readMetadata);

        if (!os.write(&bcl.front(), bcl.size()))
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(
                errno, (boost::format("Failed to write: %d bytes into %s") % bcl.size() % tempFilePaths_[nameCrc]).str()));
        }
        tempFileSizes_[nameCrc] += bcl.size();
    }
}

void PairedEndClusterExtractor::reset()
{
    firstUnextracted_ = end();
}

inline bool inRange(const IndexRecord &idx, const char* rangeStart, const char* rangeEnd)
{
    return reinterpret_cast<const char*>(idx.bamRecordPointer_) >= rangeStart &&
        reinterpret_cast<const char*>(idx.bamRecordPointer_) < rangeEnd;
}

/**
 * \brief stores bam records in [rangeStart,rangeEnd) into temporary files and frees memory associated with them
 * \precondition All records in the extractor are assumed to be unpaired
 */
void PairedEndClusterExtractor::removeOld(
    const char* rangeStart,
    const char* rangeEnd,
    const flowcell::ReadMetadataList &readMetadataList)
{
    BaseT::iterator oldBegin =
        std::partition(begin(), end(), !boost::bind(&inRange, _1, rangeStart, rangeEnd));

//    ISAAC_THREAD_CERR << "remove old " << std::size_t(rangeStart) << "-" << std::size_t(rangeEnd) << " oldBegin!=begin():" << (oldBegin != begin()) << std::endl;
    storeUnpaired(oldBegin, end(), readMetadataList);

    erase(oldBegin, end());
    reset();
}

void PairedEndClusterExtractor::sort()
{
    std::sort(begin(), end());
    firstUnextracted_ = begin();
}

template
unsigned PairedEndClusterExtractor::extractUnpaired<std::vector<char>::iterator, std::vector<bool>::iterator >(
    const unsigned r1Length,
    const unsigned r2Length,
    unsigned clusterCount, std::vector<char>::iterator &clusterIt, std::vector<bool>::iterator &pfIt );

template
unsigned PairedEndClusterExtractor::extractPairedReads<std::vector<char>::iterator, std::vector<bool>::iterator >(
    unsigned clusterCount, std::vector<char>::iterator &clusterIt, std::vector<bool>::iterator &pfIt,
    const flowcell::ReadMetadataList &readMetadataList);

} // namespace bamDataSource
} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
