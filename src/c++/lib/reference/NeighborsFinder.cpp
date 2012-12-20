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
 ** \file NeighborsFinder.cpp
 **
 ** \brief See NeighborsFinder.hh.
 **
 ** \author Come Raczy
 **/

#include <fstream>
#include <cerrno>
#include <cstring>
#include <ctime>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "common/Debug.hh"
#include "common/ParallelSort.hpp"
#include "reference/NeighborsFinder.hh"
#include "reference/SortedReferenceXml.hh"
#include "reference/ReferenceKmer.hh"
#include "common/Exceptions.hh"
#include "oligo/Kmer.hh"
#include "oligo/Permutate.hh"

namespace isaac
{
namespace reference
{

namespace bfs = boost::filesystem;

NeighborsFinder::NeighborsFinder(
    const bfs::path &inputFile,
    const bfs::path &outputDirectory,
    const bfs::path &outputFile,
    const bfs::path &tempFile,
    const unsigned jobs)
    : inputFile_(inputFile)
    , outputDirectory_(outputDirectory)
    , outputFile_(outputFile)
    , tempFile_(tempFile)
    , jobs_(jobs)
{
}

/// Helper function to load the list of kmers
//NeighborsFinder::KmerList getKmerList(const boost::filesystem::path &sortedReferenceXml);
NeighborsFinder::KmerList getKmerList(const SortedReferenceXml &sortedReferenceXml);


void NeighborsFinder::run() const
{
    using boost::format;
    using common::IoException;
    std::ifstream is(inputFile_.string().c_str());
    SortedReferenceXml sortedReferenceXml;
    if (!is || !(is >> sortedReferenceXml))
    {
        const format message = format("Failed to open XML sorted reference: %s") % strerror(errno);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }

    std::ofstream os(outputFile_.string().c_str());
    if (!os)
    {
        const format message = format("Failed to open output file file %s for writing: %s") % outputFile_ % strerror(errno);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }

    generateNeighbors(sortedReferenceXml);
    std::vector<MaskFile> maskFileList = sortedReferenceXml.getMaskFileList("ABCD");
    updateSortedReference(maskFileList);

    BOOST_FOREACH(MaskFile &maskFile, maskFileList)
    {
        sortedReferenceXml.addMaskFile("ABCD", maskFile.maskWidth, maskFile.mask, maskFile.path,
                                       maskFile.kmers, maskFile.maxPrefixRangeCount);
    }

    if (!(os << sortedReferenceXml))
    {
        const format message = format("Failed to write output file file %s: %s") % outputFile_ % strerror(errno);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
}

oligo::Kmer reverseComplement(oligo::Kmer kmer)
{
    oligo::Kmer reversed = 0;
    kmer = ~kmer; // complement all the bases
    for (unsigned i = 0; oligo::kmerLength > i; ++i)
    {
        reversed <<= 2;
        reversed |= (kmer & 3);
        kmer >>= 2;
    }
    return reversed;
}

void NeighborsFinder::updateSortedReference(std::vector<MaskFile> &maskFileList) const
{
    using boost::format;
    using common::IoException;
    std::ifstream neighbors(tempFile_.string().c_str());
    if (!neighbors)
    {
        const format message = format("Failed to open neighbors file %s for reading: %s") % tempFile_ % strerror(errno);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
    clock_t start;
    BOOST_FOREACH(MaskFile &maskFile, maskFileList)
    {
        start = clock();
        const bfs::path oldMaskFile = maskFile.path; //bfs::path(maskFile.path).replace_extension(".orig");
        ISAAC_THREAD_CERR << "Annotating " << oldMaskFile << std::endl;
        if (!exists(oldMaskFile))
        {
            const format message = format("Mask file %s does not exist: %s") % oldMaskFile;
            BOOST_THROW_EXCEPTION(IoException(ENOENT, message.str()));
        }
/*
        boost::system::error_code errorCode;
        rename(maskFile.path, oldMaskFile, errorCode);
        if (errorCode)
        {
            const format message = format("Failed to rename mask file %s to %s: %s") % maskFile.path % oldMaskFile % errorCode.message();
            BOOST_THROW_EXCEPTION(IoException(errorCode.value(), message.str()));
        }
*/
        maskFile.path = outputDirectory_ / maskFile.path.filename();
        std::ifstream maskInput(oldMaskFile.string().c_str());
        if (!maskInput)
        {
            const format message = format("Failed to open mask file %s for reading: %s") % oldMaskFile % strerror(errno);
            BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
        }
        std::ofstream maskOutput(maskFile.path.string().c_str());
        if (!maskOutput)
        {
            const format message = format("Failed to open mask file %s for writing: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
        }
        oligo::Kmer currentNeighbor = 0;
        neighbors.read(reinterpret_cast<char *>(&currentNeighbor), sizeof(currentNeighbor));
        while(maskInput && maskOutput)
        {
            ReferenceKmer referenceKmer;
            if (maskInput.read(reinterpret_cast<char *>(&referenceKmer), sizeof(referenceKmer)))
            {
                while (neighbors && currentNeighbor < referenceKmer.getKmer())
                {
                    neighbors.read(reinterpret_cast<char *>(&currentNeighbor), sizeof(currentNeighbor));
                }
                if (!neighbors && !neighbors.eof())
                {
                    const format message = format("Failed to read neighbor from %s: %s") % tempFile_ % strerror(errno);
                    BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
                }

                referenceKmer.getReferencePosition().setNeighbors(neighbors && currentNeighbor == referenceKmer.getKmer());

                if (!maskOutput.write(reinterpret_cast<const char *>(&referenceKmer), sizeof(referenceKmer)))
                {
                    const format message = format("Failed to write reference k-mer into %s: %s") % maskFile.path % strerror(errno);
                    BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
                }
            }
        }
        if (!maskInput.eof() && !neighbors.eof())
        {
            const format message = format("Failed to update %s with neighbors information: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
        }
        ISAAC_THREAD_CERR << "Adding neighbors information done in " << (clock() - start) / 1000 << " ms for " << maskFile.path << std::endl;
    }
}

inline bool compareAnnotatedKmer(const NeighborsFinder::AnnotatedKmer &lhs, const NeighborsFinder::AnnotatedKmer &rhs)
{
    return lhs.value < rhs.value;
}

inline bool compareAnnotatedKmerMask(const NeighborsFinder::AnnotatedKmer &lhs, const NeighborsFinder::AnnotatedKmer &rhs)
{
    return (lhs.value >> oligo::kmerLength) < (rhs.value >> oligo::kmerLength);
}

inline bool isAnnotatedKmerEqual(const NeighborsFinder::AnnotatedKmer &lhs, const NeighborsFinder::AnnotatedKmer &rhs)
{
    return lhs.value == rhs.value;
}

void NeighborsFinder::generateNeighbors(const SortedReferenceXml &sortedReferenceXml) const
{
    KmerList kmerList = getKmerList(sortedReferenceXml);
    std::vector<oligo::Permutate> permutateList = oligo::getPermutateList(4);
    // iterate over all possible permutations
    clock_t start = clock();
    BOOST_FOREACH(const oligo::Permutate &permutate, permutateList)
    {
        start = clock();
        ISAAC_THREAD_CERR << "Permuting all k-mers (" << kmerList.size() << " k-mers) " << permutate.toString() << std::endl;
        BOOST_FOREACH(AnnotatedKmer &kmer, kmerList)
        {
            //std::cerr << (boost::format("    %016x:%d") % kmer.value % (int)kmer.count).str();
            kmer.value = permutate(kmer.value);
            //std::cerr << (boost::format(" %016x:%d") % kmer.value % (int)kmer.count).str() << std::endl;
        }
        ISAAC_THREAD_CERR << "Permuting all k-mers done (" << kmerList.size() << " k-mers) " << permutate.toString() << " in " << (clock() - start) / 1000 << " ms" << std::endl;
        start = clock();
        ISAAC_THREAD_CERR << "Sorting all k-mers (" << kmerList.size() << " k-mers)" << std::endl;
        //std::sort(kmerList.begin(), kmerList.end());
        common::parallelSort(kmerList, compareAnnotatedKmerMask);
        ISAAC_THREAD_CERR << "Sorting all k-mers done in " << (clock() - start) / 1000 << " ms" << std::endl;
        start = clock();
        ISAAC_THREAD_CERR << "Finding neighbors" << std::endl;
        // find neighbors
        findNeighbors(kmerList, jobs_);
        unsigned count = 0;
        ISAAC_THREAD_CERR << "Counting neighbors" << std::endl;
        BOOST_FOREACH(AnnotatedKmer &kmer, kmerList)
        {
            if (kmer.hasNeighbors)
            {
                ++count;
            }
        }
        ISAAC_THREAD_CERR << "Found " << count
                          << " neighbors in " << kmerList.size() << " kmers" << std::endl;
        ISAAC_THREAD_CERR << "Finding neighbors done in " << (clock() - start) / 1000 << " ms" << std::endl;
    }
    start = clock();
    ISAAC_THREAD_CERR << "Reordering all k-mers (" << kmerList.size() << " k-mers)" << std::endl;
    BOOST_FOREACH(AnnotatedKmer &kmer, kmerList)
    {
        kmer.value = permutateList.back().reorder(kmer.value);
    }
    ISAAC_THREAD_CERR << "Reordering all k-mers done in " << (clock() - start) / 1000 << " ms" << std::endl;
    start = clock();
    ISAAC_THREAD_CERR << "Sorting all k-mers (" << kmerList.size() << " k-mers)" << std::endl;
    std::sort(kmerList.begin(), kmerList.end());
    ISAAC_THREAD_CERR << "Sorting all k-mers done in " << (clock() - start) / 1000 << " ms" << std::endl;

    storeNeighborKmers(kmerList);
}

void NeighborsFinder::storeNeighborKmers(const KmerList &kmerList) const
{
    std::ofstream neighbors(tempFile_.string().c_str());
    if (!neighbors)
    {
        using boost::format;
        using common::IoException;
        const format message = format("Failed to open neighbors file %s for writing: %s") % tempFile_ % strerror(errno);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
    const clock_t start = clock();
    ISAAC_THREAD_CERR << "Saving all k-mers with non-equal neighbors (" << kmerList.size() << " k-mers)" << std::endl;
    unsigned neighborsCount = 0;
    BOOST_FOREACH(const AnnotatedKmer &annotatedKmer, kmerList)
    {
        if (annotatedKmer.hasNeighbors)
        {
            ++neighborsCount;
            //std::cerr << (boost::format(" Saving %016x:%d") % annotatedKmer.value % (int)annotatedKmer.count).str() << std::endl;
            if (!neighbors.write(reinterpret_cast<const char *>(&annotatedKmer.value), sizeof(annotatedKmer.value)))
            {
                using boost::format;
                using common::IoException;
                const format message = format("Failed to write neighborhood into file %s: %s") % tempFile_ % strerror(errno);
                BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
            }
        }
    }
    ISAAC_THREAD_CERR << "Saving all k-mers with non-equal neighbors done in " << (clock() - start) / 1000 << " ms (" 
                      << neighborsCount << "/" << kmerList.size() << " have non-equal neighbors)" << std::endl;
}

void NeighborsFinder::findNeighbors(KmerList &kmerList, const unsigned jobs)
{
    KmerList::iterator kmerListBegin = kmerList.begin();
    boost::thread_group threads;
    while (kmerList.end() != kmerListBegin)
    {
        assert(threads.size() < jobs);
        const unsigned remainingThreads = jobs - threads.size();
        const size_t tailSize = kmerList.end() - kmerListBegin;
        KmerList::iterator kmerListEnd = kmerListBegin + tailSize / remainingThreads;
        const oligo::Kmer currentPrefix = (kmerListEnd->value) >> oligo::kmerLength;
        while (kmerList.end() != kmerListEnd && currentPrefix == (kmerListEnd->value) >> oligo::kmerLength)
        {
            ++kmerListEnd;
        }
        threads.create_thread(boost::bind(&NeighborsFinder::findNeighborsParallel, kmerListBegin, kmerListEnd));
        kmerListBegin = kmerListEnd;
    }
    threads.join_all();
}

void NeighborsFinder::findNeighborsParallel(const KmerList::iterator kmerListBegin, const KmerList::iterator kmerListEnd)
{
    const clock_t start = clock();
    ISAAC_THREAD_CERR << "findNeighborsParallel: " << (unsigned)(kmerListEnd - kmerListBegin) << " kmers" << std::endl;
    KmerList::iterator blockBegin = kmerListBegin;
    while (kmerListEnd != blockBegin)
    {
        const oligo::Kmer currentPrefix = (blockBegin->value) >> oligo::kmerLength;
        KmerList::iterator blockEnd = blockBegin;
        while (kmerListEnd != blockEnd && currentPrefix == (blockEnd->value >> oligo::kmerLength))
        {
            ++blockEnd;
        }

        for (KmerList::iterator i = blockBegin; blockEnd > i; ++i)
        {
            markNeighbors(i, blockEnd);
        }
        blockBegin = blockEnd;
    }
    ISAAC_THREAD_CERR << "findNeighborsParallel done in " << (clock() - start) / 1000 << " ms: " 
                      << (unsigned)(kmerListEnd - kmerListBegin) << " kmers" << std::endl;
}

/**
 * \brief Marks all kmers within 4 mismatches of *blockBegin.
 *        Also marks *blockBegin if it has any neighbors.
 *
 */
void NeighborsFinder::markNeighbors(
    const KmerList::iterator blockBegin,
    const KmerList::const_iterator blockEnd)
{
    ISAAC_ASSERT_MSG(blockBegin <= blockEnd, "Improper range");

    KmerList::iterator current = blockBegin;
    AnnotatedKmer &kmer = *blockBegin;
    while (blockEnd != current)
    {
        if (!kmer.hasNeighbors || !current->hasNeighbors)
        {
            oligo::Kmer values[] = {kmer.value, current->value};
            unsigned width = oligo::kmerLength / 2;
            unsigned mismatchCount = 0;
            while (4 >= mismatchCount && width--)
            {
                const oligo::Kmer x = values[0] ^ values[1];
                if (!x)
                {
                    // If there is no difference, no point to continue counting.
                    // Also ensures that matching kmers don't get marked as neighbors
                    break;
                }
                if (3 & x)
                {
                    ++mismatchCount;
                }
                values[0] >>= 2;
                values[1] >>= 2;
            }
            if (mismatchCount && 4 >= mismatchCount)
            {
                kmer.hasNeighbors = true;
                current->hasNeighbors = true;
            }
        }
        ++current;
    }
}

//NeighborsFinder::KmerList getKmerList(const boost::filesystem::path &sortedReferenceXml)

NeighborsFinder::AnnotatedKmer reverseComplementAnnotatedKmer(NeighborsFinder::AnnotatedKmer ak)
{
    ak.value = reverseComplement(ak.value);
    return ak;
}

NeighborsFinder::KmerList getKmerList(const SortedReferenceXml &sortedReferenceXml)
{
    const std::vector<MaskFile> maskFileList = sortedReferenceXml.getMaskFileList("ABCD");
    typedef NeighborsFinder::KmerList KmerList;
    typedef NeighborsFinder::AnnotatedKmer AnnotatedKmer;
    NeighborsFinder::KmerList kmerList;
    ISAAC_THREAD_CERR << "reserving memory for " << sortedReferenceXml.getTotalKmers() * 2 << " kmers" << std::endl;
    kmerList.reserve(sortedReferenceXml.getTotalKmers() * 2);
    ISAAC_THREAD_CERR << "reserving memory done for " << kmerList.capacity() << " kmers" << std::endl;
    // load all the kmers found in all the ABCD files
    //const std::vector<MaskFile> maskFileList = getMaskFileList(sortedReferenceXml, "ABCD");
    BOOST_FOREACH(const MaskFile &maskFile, maskFileList)
    {
        std::ifstream is(maskFile.path.string().c_str());
        if (!is)
        {
            using boost::format;
            using common::IoException;
            const format message = format("Failed to open sorted reference file %s: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
        }
        ReferenceKmer referenceKmer;
        while (is.read(reinterpret_cast<char*>(&referenceKmer), sizeof(referenceKmer)))
        {
            if (kmerList.empty() || referenceKmer.getKmer() != kmerList.back().value)
            {
                kmerList.push_back(AnnotatedKmer(referenceKmer.getKmer(), false));
            }
        }
        if (!is.eof())
        {
            using boost::format;
            using common::IoException;
            const format message = format("Failed to read sorted reference file %s: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
        }
    }
    ISAAC_THREAD_CERR << "loading done for " << kmerList.size() << " unique forward kmers" << std::endl;
    std::transform(kmerList.begin(), kmerList.end(), std::back_inserter(kmerList), &reverseComplementAnnotatedKmer);
    ISAAC_THREAD_CERR << "generating reverse complements done for " << kmerList.size() / 2 << " unique forward kmers" << std::endl;

    common::parallelSort(kmerList, compareAnnotatedKmer);
    kmerList.erase(std::unique(kmerList.begin(), kmerList.end(), isAnnotatedKmerEqual), kmerList.end());
    ISAAC_THREAD_CERR << "removing complement duplicates done for " << kmerList.size() << " unique kmers (forward + reverse) " << std::endl;

    // We have the memory to produce the full copy of data because the current implementation of the parallelSort
    // requires it. Reduce the footprint by getting rid of the chunk of kmerList, containing the duplicate kmers.
    NeighborsFinder::KmerList ret(kmerList.begin(), kmerList.end());
    return ret;
}

} // namespace reference
} //namespace isaac
