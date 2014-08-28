/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file FastqLoader.hh
 **
 ** Component to read FASTQ files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FASTQ_LOADER_HH
#define iSAAC_IO_FASTQ_LOADER_HH

#include "common/Threads.hpp"
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
    common::ThreadVector &threads_;
    const unsigned inputLoadersMax_;

public:
    /**
     * \brief initializes paired fastq loader.
     */
    /**
     * \brief Creates uninitialized fastq loader
     */
    FastqLoader(
        const bool allowVariableLength,
        std::size_t maxPathLength,
        common::ThreadVector &threads,
        const unsigned inputLoadersMax) :
        read1Reader_(allowVariableLength),
        read2Reader_(allowVariableLength),
        paired_(false),
        threads_(threads),
        inputLoadersMax_(inputLoadersMax)
    {
        read1Reader_.reservePathBuffers(maxPathLength);
        read2Reader_.reservePathBuffers(maxPathLength);
    }

    void open(
        const boost::filesystem::path &read1Path)
    {
        read1Reader_.open(read1Path);
        paired_ = false;
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
        if (1 == readMetadataList.size())
        {
            return loadSingleRead(read1Reader_, clusterCount, readMetadataList.at(0), 0, it);
        }
        else
        {
            ISAAC_ASSERT_MSG(2 == readMetadataList.size(), "Only paired and single-ended data is supported");
            unsigned readClusters[2] = {0,0};
            InsertIt it1 = it;
            if (2 <= inputLoadersMax_)
            {
                it += readMetadataList.at(0).getLength();
                boost::reference_wrapper<InsertIt> insertIterators[] = {boost::ref(it1), boost::ref(it)};
                threads_.execute(boost::bind(&FastqLoader::threadLoadPairedReads<InsertIt>, this,
                                            clusterCount, boost::ref(readMetadataList),
                                            boost::ref(readClusters), insertIterators, _1),
                                 2);
            }
            else
            {
                ISAAC_ASSERT_MSG(1 == inputLoadersMax_, "At least one thread is expected for IO")
                readClusters[0] = loadSingleRead(read1Reader_, clusterCount, readMetadataList.at(0), readMetadataList.at(1).getLength(), it1);
                it += readMetadataList.at(0).getLength();
                readClusters[1] = loadSingleRead(read2Reader_, clusterCount, readMetadataList.at(1), readMetadataList.at(0).getLength(), it);
            }

            if (readClusters[0] != readClusters[1])
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Mismatching number of cluster read for r1/r2 = %d/%d, files: %s/%s") %
                    readClusters[0] % readClusters[1] % read1Reader_.getPath() % read2Reader_.getPath()).str()));
            }

            return readClusters[0];
        }

    }
private:
    template <typename InsertIt>
    static unsigned loadSingleRead(FastqReader &reader, unsigned clusterCount,
                            const flowcell::ReadMetadata &readMetadata,
                            const unsigned step, InsertIt &it)
    {
        unsigned clustersToRead = clusterCount;
        for (;clustersToRead && reader.hasData();)
        {
            it = reader.extractBcl(readMetadata, it);
            reader.next();
            // avoid debug glibc complaining about advancing iterator past the end of the container
            if (--clustersToRead)
            {
                std::advance(it, step);
            }
        }
        return clusterCount - clustersToRead;
    }

    template <typename InsertIt>
    void threadLoadPairedReads(unsigned clusterCount,
                              const flowcell::ReadMetadataList &readMetadataList,
                              unsigned readClusters[2], boost::reference_wrapper<InsertIt> insertIterators[2], int threadNumber)
    {
        readClusters[threadNumber] = loadSingleRead(
            threadNumber ? read2Reader_ : read1Reader_,
            clusterCount, readMetadataList.at(threadNumber),
            readMetadataList.at((threadNumber + 1) % 2).getLength(), insertIterators[threadNumber].get());
    }
};

//#pragma GCC pop_options

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FASTQ_LOADER_HH
