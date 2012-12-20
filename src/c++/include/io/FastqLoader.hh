/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/downloads/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file FastqLoader.hh
 **
 ** Component to read FASTQ files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FASTQ_LOADER_HH
#define iSAAC_IO_FASTQ_LOADER_HH

#include "io/FastqReader.hh"

namespace isaac
{
namespace io
{

//#pragma GCC push_options
//#pragma GCC optimize ("0")

class FastqLoader
{
    FastqReader read1Reader_;
    FastqReader read2Reader_;
    bool paired_;
public:
    /**
     * \brief initializes paired fastq loader.
     */

    FastqLoader(
        const bool allowVariableLength,
        const boost::filesystem::path &read1Path,
        const boost::filesystem::path &read2Path) :
        read1Reader_(allowVariableLength, read1Path), read2Reader_(allowVariableLength, read2Path), paired_(true)
    {
    }


    /**
     * \brief Creates uninitialized fastq loader
     */
    FastqLoader(const bool allowVariableLength) :
        read1Reader_(allowVariableLength), read2Reader_(allowVariableLength), paired_(false)
    {
    }

    void reservePathBuffers(std::size_t maxPathLength)
    {
        read1Reader_.reservePathBuffers(maxPathLength);
        read2Reader_.reservePathBuffers(maxPathLength);
    }

    void open(
        const boost::filesystem::path &read1Path,
        const boost::filesystem::path &read2Path)
    {
        read1Reader_.open(read1Path);
        read2Reader_.open(read2Path);
        paired_ = true;
    }

    /**
     * \param clusterCount  Maximum number of clusters to load
     * \param clusterLength Expected cluster length. Will fail if the read cluster length does not match
     * \param it            Insert iterator for the buffer that is sufficient to load the clusterCount
     *                      clusters of clusterLength
     *
     * \return Actual number of loaded clusters
     */
    template <typename InsertIt>
    unsigned loadClusters(unsigned clusterCount, const flowcell::ReadMetadataList &readMetadataList, InsertIt &it)
    {
        ISAAC_ASSERT_MSG(paired_, "Single-ended fastq data is not supported yet");
        unsigned readClusters = 0;
        while (readClusters < clusterCount && read1Reader_.hasData())
        {
            read1Reader_.getBcl(readMetadataList.at(0), it);
            if (!read2Reader_.hasData())
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Read2 fastq file ended before read1 is over. Read clusters %d, %s") %
                    readClusters % read2Reader_.getPath()).str()));
            }
            read2Reader_.getBcl(readMetadataList.at(1), it);
            read1Reader_.next();
            read2Reader_.next();
            ++readClusters;
        }
        if (!read1Reader_.hasData() && read2Reader_.hasData())
        {
//            ISAAC_ASSERT_MSG(false, "Read1 fastq file ended before read2 is over. Read clusters" << readClusters);
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Read1 fastq file ended before read2 is over. Read clusters %d, file: %s") %
                readClusters % read1Reader_.getPath()).str()));
        }
        return readClusters;
    }
};

//#pragma GCC pop_options

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FASTQ_LOADER_HH
