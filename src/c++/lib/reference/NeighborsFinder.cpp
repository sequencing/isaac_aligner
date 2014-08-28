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

template <typename KmerT>
NeighborsFinder<KmerT>::NeighborsFinder(
    const bool parallelSort,
    const bfs::path &inputFile,
    const bfs::path &outputDirectory,
    const bfs::path &outputFile,
    const bfs::path &tempFile,
    const unsigned jobs)
    : parallelSort_(parallelSort)
    , inputFile_(inputFile)
    , outputDirectory_(outputDirectory)
    , outputFile_(outputFile)
    , tempFile_(tempFile)
    , jobs_(jobs)
{
}


template <typename KmerT>
void NeighborsFinder<KmerT>::run() const
{
    SortedReferenceMetadata sortedReferenceMetadata = loadSortedReferenceXml(inputFile_);

    generateNeighbors(sortedReferenceMetadata);
    updateSortedReference(sortedReferenceMetadata.getMaskFileList(oligo::KmerTraits<KmerT>::KMER_BASES));
    saveSortedReferenceXml(outputFile_, sortedReferenceMetadata);
}

template <typename KmerT>
KmerT reverseComplement(KmerT kmer)
{
    KmerT reversed = 0;
    kmer = ~kmer; // complement all the bases
    for (unsigned i = 0; oligo::KmerTraits<KmerT>::KMER_BASES > i; ++i)
    {
        reversed <<= 2;
        reversed |= (kmer & 3);
        kmer >>= 2;
    }
    return reversed;
}

template <typename KmerT>
void NeighborsFinder<KmerT>::updateSortedReference(SortedReferenceMetadata::MaskFiles &maskFileList) const
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
    KmerT currentNeighbor = 0;
    neighbors.read(reinterpret_cast<char *>(&currentNeighbor), sizeof(currentNeighbor));
    BOOST_FOREACH(SortedReferenceMetadata::MaskFile &maskFile, maskFileList)
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
        while(maskInput && maskOutput)
        {
            ReferenceKmer<KmerT> referenceKmer;
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

                referenceKmer.setNeighbors(neighbors && currentNeighbor == referenceKmer.getKmer());

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

template <typename KmerT>
inline bool compareAnnotatedKmer(
    const typename NeighborsFinder<KmerT>::AnnotatedKmer &lhs,
    const typename NeighborsFinder<KmerT>::AnnotatedKmer &rhs)
{
    return lhs.value < rhs.value;
}

template <typename KmerT>
inline bool compareAnnotatedKmerMask(
    const typename NeighborsFinder<KmerT>::AnnotatedKmer &lhs,
    const typename NeighborsFinder<KmerT>::AnnotatedKmer &rhs)
{
    return (lhs.value >> oligo::KmerTraits<KmerT>::KMER_BASES) < (rhs.value >> oligo::KmerTraits<KmerT>::KMER_BASES);
}

template <typename KmerT>
inline bool isAnnotatedKmerEqual(
    const typename NeighborsFinder<KmerT>::AnnotatedKmer &lhs,
    const typename NeighborsFinder<KmerT>::AnnotatedKmer &rhs)
{
    return lhs.value == rhs.value;
}

template <typename KmerT>
void NeighborsFinder<KmerT>::generateNeighbors(const SortedReferenceMetadata &sortedReferenceMetadata) const
{
    KmerList kmerList = getKmerList(sortedReferenceMetadata);
    std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<KmerT>(4);
    // iterate over all possible permutations
    clock_t start = clock();
    BOOST_FOREACH(const oligo::Permutate &permutate, permutateList)
    {
        start = clock();
        ISAAC_THREAD_CERR << "Permuting all k-mers (" << kmerList.size() << " k-mers) " << permutate.toString() << std::endl;
        BOOST_FOREACH(AnnotatedKmer &kmer, kmerList)
        {
//            std::cerr << (boost::format("    %016x:%d") % oligo::traceKmer(kmer.value) % (int)kmer.hasNeighbors).str();
            kmer.value = permutate(kmer.value);
//            std::cerr << (boost::format(" %016x:%d") % oligo::traceKmer(kmer.value) % (int)kmer.hasNeighbors).str() << std::endl;
        }
        ISAAC_THREAD_CERR << "Permuting all k-mers done (" << kmerList.size() << " k-mers) " << permutate.toString() << " in " << (clock() - start) / 1000 << " ms" << std::endl;
        start = clock();
        ISAAC_THREAD_CERR << "Sorting all k-mers (" << kmerList.size() << " k-mers)" << std::endl;
        if (parallelSort_)
        {
            common::parallelSort(kmerList, &compareAnnotatedKmerMask<KmerT>);
        }
        else
        {
            std::sort(kmerList.begin(), kmerList.end());
        }
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

template <typename KmerT>
void NeighborsFinder<KmerT>::storeNeighborKmers(const KmerList &kmerList) const
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

template <typename KmerT>
void NeighborsFinder<KmerT>::findNeighbors(KmerList &kmerList, const unsigned jobs)
{
    typename KmerList::iterator kmerListBegin = kmerList.begin();
    boost::thread_group threads;
    while (kmerList.end() != kmerListBegin)
    {
        assert(threads.size() < jobs);
        const unsigned remainingThreads = jobs - threads.size();
        const size_t tailSize = std::distance(kmerListBegin, kmerList.end());
        typename KmerList::iterator kmerListEnd = kmerListBegin + tailSize / remainingThreads;
        if (kmerList.end() != kmerListEnd)
        {
            const KmerT currentPrefix = (kmerListEnd->value) >> oligo::KmerTraits<KmerT>::KMER_BASES;
            while (kmerList.end() != kmerListEnd && currentPrefix == (kmerListEnd->value) >> oligo::KmerTraits<KmerT>::KMER_BASES)
            {
                ++kmerListEnd;
            }
        }
        threads.create_thread(boost::bind(&NeighborsFinder::findNeighborsParallel, kmerListBegin, kmerListEnd));
        kmerListBegin = kmerListEnd;
    }
    threads.join_all();
}

template <typename KmerT>
void NeighborsFinder<KmerT>::findNeighborsParallel(
    const typename KmerList::iterator kmerListBegin,
    const typename KmerList::iterator kmerListEnd)
{
    const clock_t start = clock();
    ISAAC_THREAD_CERR << "findNeighborsParallel: " << (unsigned)(kmerListEnd - kmerListBegin) << " kmers" << std::endl;
    typename KmerList::iterator blockBegin = kmerListBegin;
    while (kmerListEnd != blockBegin)
    {
        const KmerT currentPrefix = (blockBegin->value) >> oligo::KmerTraits<KmerT>::KMER_BASES;
        typename KmerList::iterator blockEnd = blockBegin;
        while (kmerListEnd != blockEnd && currentPrefix == (blockEnd->value >> oligo::KmerTraits<KmerT>::KMER_BASES))
        {
            ++blockEnd;
        }

        for (typename KmerList::iterator i = blockBegin; blockEnd > i; ++i)
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
template <typename KmerT>
void NeighborsFinder<KmerT>::markNeighbors(
    const typename KmerList::iterator blockBegin,
    const typename KmerList::const_iterator blockEnd)
{
    ISAAC_ASSERT_MSG(blockBegin <= blockEnd, "Improper range");

    typename KmerList::iterator current = blockBegin;
    AnnotatedKmer &kmer = *blockBegin;
    while (blockEnd != current)
    {
        if (!kmer.hasNeighbors || !current->hasNeighbors)
        {
            KmerT values[] = {kmer.value, current->value};
            unsigned width = oligo::KmerTraits<KmerT>::KMER_BASES / 2;
            unsigned mismatchCount = 0;
            while (4 >= mismatchCount && width--)
            {
                const KmerT x = values[0] ^ values[1];
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

//NeighborsFinder::KmerList getKmerList(const boost::filesystem::path &sortedReferenceMetadata)

template <typename KmerT>
typename NeighborsFinder<KmerT>::AnnotatedKmer reverseComplementAnnotatedKmer(typename NeighborsFinder<KmerT>::AnnotatedKmer ak)
{
    ak.value = reverseComplement(ak.value);
    return ak;
}

template <typename KmerT>
typename NeighborsFinder<KmerT>::KmerList NeighborsFinder<KmerT>::getKmerList(const SortedReferenceMetadata &sortedReferenceMetadata) const
{
    const std::vector<SortedReferenceMetadata::MaskFile> &maskFileList =
        sortedReferenceMetadata.getMaskFileList(oligo::KmerTraits<KmerT>::KMER_BASES);
    KmerList kmerList;
    ISAAC_THREAD_CERR << "reserving memory for " <<
        sortedReferenceMetadata.getTotalKmers(oligo::KmerTraits<KmerT>::KMER_BASES) * 2 << " kmers" << std::endl;
    kmerList.reserve(sortedReferenceMetadata.getTotalKmers(oligo::KmerTraits<KmerT>::KMER_BASES) * 2);
    ISAAC_THREAD_CERR << "reserving memory done for " << kmerList.capacity() << " kmers" << std::endl;
    // load all the kmers found in all the ABCD files
    BOOST_FOREACH(const SortedReferenceMetadata::MaskFile &maskFile, maskFileList)
    {
        std::ifstream is(maskFile.path.string().c_str());
        if (!is)
        {
            using boost::format;
            using common::IoException;
            const format message = format("Failed to open sorted reference file %s: %s") % maskFile.path % strerror(errno);
            BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
        }
        ReferenceKmer<KmerT> referenceKmer;
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
    std::transform(kmerList.begin(), kmerList.end(), std::back_inserter(kmerList), &reverseComplementAnnotatedKmer<KmerT>);
    ISAAC_THREAD_CERR << "generating reverse complements done for " << kmerList.size() / 2 << " unique forward kmers" << std::endl;

    if (parallelSort_)
    {
        common::parallelSort(kmerList, &compareAnnotatedKmer<KmerT>);
    }
    else
    {
        std::sort(kmerList.begin(), kmerList.end(), &compareAnnotatedKmer<KmerT>);
    }
    kmerList.erase(std::unique(kmerList.begin(), kmerList.end(), &isAnnotatedKmerEqual<KmerT>), kmerList.end());
    ISAAC_THREAD_CERR << "removing complement duplicates done for " << kmerList.size() << " unique kmers (forward + reverse) " << std::endl;

    return kmerList;
}

template class NeighborsFinder<oligo::ShortKmerType>;
template class NeighborsFinder<oligo::KmerType>;
template class NeighborsFinder<oligo::LongKmerType>;

} // namespace reference
} //namespace isaac
