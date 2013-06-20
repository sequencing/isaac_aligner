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
 ** \file ReferenceSorter.cpp
 **
 ** Top level component to produce a sorted reference.
 **
 ** \author Come Raczy
 **/

#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/io/ios_state.hpp>

#include "common/Exceptions.hh"
#include "common/SystemCompatibility.hh"
#include "io/BitsetLoader.hh"
#include "io/BitsetSaver.hh"
#include "io/FastaReader.hh"
#include "oligo/Nucleotides.hh"
#include "oligo/Mask.hh"
#include "reference/ReferencePosition.hh"
#include "reference/ReferenceSorter.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

template <typename KmerT>
ReferenceSorter<KmerT>::ReferenceSorter (
    const unsigned int maskWidth,
    const unsigned mask,
    const boost::filesystem::path &genomeFile,
    const boost::filesystem::path &genomeNeighborsFile,
    const boost::filesystem::path &outputFile,
    const unsigned repeatThreshold
    )
    : repeatThreshold_(repeatThreshold)
    , maskWidth_(maskWidth)
    , mask_(mask)
    , msbMask_(~((~(KmerT(0)) >> maskWidth)))
    , maskBits_(KmerT(mask_) << (oligo::KmerTraits<KmerT>::KMER_BASES * oligo::BITS_PER_BASE - maskWidth_))
    , unpermutatedMsbMask_(msbMask_)
    , unpermutatedMaskBits_(maskBits_)
    , genomeFile_(genomeFile)
    , genomeNeighborsFile_(genomeNeighborsFile)
    , outputFile_(boost::filesystem::absolute(outputFile))
{
    boost::io::ios_flags_saver svr(std::cerr);
    std::cerr <<
            "Constructing ReferenceSorter: for " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers " <<
            " mask width: " << maskWidth_ <<
            " msbMask_: " << oligo::TraceKmer<KmerT>(msbMask_) <<
            " maskBits_: " << oligo::TraceKmer<KmerT>(maskBits_) <<
            " unpermutatedMsbMask_:" << std::hex << oligo::TraceKmer<KmerT>(unpermutatedMsbMask_) <<
            " unpermutatedMaskBits_:" << std::hex << oligo::TraceKmer<KmerT>(unpermutatedMaskBits_) <<
            " genomeFile_: " << genomeFile_ <<
            " outputFile_: " << outputFile_ <<
            std::endl;

    BOOST_ASSERT(mask_ < isaac::oligo::getMaskCount(maskWidth_) && "Mask value cannot exceed the allowed bit width");
}


template <typename KmerT>
void ReferenceSorter<KmerT>::run()
{
    std::vector<unsigned long> contigOffsets;

    unsigned long genomeLength = 0;
    {
        genomeLength = loadReference(contigOffsets);
        ISAAC_THREAD_CERR << "Loaded genome from " << genomeFile_ << " found " << genomeLength << " bases" << std::endl;
    }

    sortReference();

    std::vector<bool> neighbors;
    if (!genomeNeighborsFile_.empty())
    {
        io::BitsetLoader loader(genomeNeighborsFile_);
        const unsigned long neighborsCount = loader.load(genomeLength, neighbors);
        ISAAC_THREAD_CERR << "Scanning " << genomeNeighborsFile_ << " found " << neighborsCount << " neighbors among " << genomeLength << " bases" << std::endl;
    }
    saveReference(contigOffsets, neighbors, genomeLength);
}

/**
 * \brief Load kmers matching the unpermutatedMaskBits_.
 *
 * \return vector of contig base offsets in the order found in fasta file
 */
template <typename KmerT>
unsigned long ReferenceSorter<KmerT>::loadReference(
    std::vector<unsigned long> &retContigOffsets)
{
    std::cerr << "Loading " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << std::endl;
    const clock_t start = clock();

    retContigOffsets.clear();
    isaac::io::MultiFastaReader multiFastaReader(std::vector<boost::filesystem::path>(1, genomeFile_));
    char base;
    KmerT forward = 0;
    KmerT reverse = 0;
    unsigned int badKmer = oligo::KmerTraits<KmerT>::KMER_BASES;
    unsigned long position = 0;
    bool newContig = false;
    unsigned long lastContigOffset = 0UL;
    while (multiFastaReader.get(base, newContig))
    {
        if (newContig)
        {
            std::cerr << "New contig: " << multiFastaReader.getContigId() << " found at position " << position << std::endl;
            position = 0;
            badKmer = oligo::KmerTraits<KmerT>::KMER_BASES;
            retContigOffsets.push_back(lastContigOffset);
        }
        if (badKmer)
        {
            --badKmer;
        }
        const KmerT baseValue = oligo::getValue(base);
        if (baseValue >> oligo::BITS_PER_BASE)
        {
            badKmer = oligo::KmerTraits<KmerT>::KMER_BASES;
        }
        forward <<= oligo::BITS_PER_BASE;
        forward |= baseValue;
        reverse >>= oligo::BITS_PER_BASE;
        reverse |= (((~baseValue) & oligo::BITS_PER_BASE_MASK) << (oligo::BITS_PER_BASE * oligo::KmerTraits<KmerT>::KMER_BASES - oligo::BITS_PER_BASE));
        //std:: cerr << forward << '\t' << position << '\t' << reverse << '\t' << position << '\t' << bad16mer << '\t' << base << '\t' << baseValue << '\n';
        const unsigned long kmerPosition = (position + 1) - oligo::KmerTraits<KmerT>::KMER_BASES;
        if (0 == badKmer)
        {
            ISAAC_ASSERT_MSG(position + 1 >= oligo::KmerTraits<KmerT>::KMER_BASES, "Kmers at the start of contig are always bad");
            addToReference(forward, ReferencePosition(multiFastaReader.getContigId(), kmerPosition, false));
            // Reverse kmers are added only to be able to properly count repeats.
            // Mark them as having neighbors to be able to filter them out before storing the results.
            addToReference(reverse, ReferencePosition(multiFastaReader.getContigId(), kmerPosition, true));
        }
        ++position;
        ++lastContigOffset;
    }
    std::cerr << "Loading " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << " done in " << (clock() - start) / 1000 << "ms" << std::endl;
    return lastContigOffset;
}

template <typename KmerT>
void ReferenceSorter<KmerT>::addToReference(const KmerT kmer, const ReferencePosition &referencePosition)
{
    if (unpermutatedMaskBits_ == (kmer & unpermutatedMsbMask_))
    {
//        std::cerr << std::hex << kmer << '\t' << referencePosition << '\n';
        reference_.push_back(ReferenceKmer<KmerT>(kmer, referencePosition));
    }
}

template <typename KmerT>
void ReferenceSorter<KmerT>::sortReference()
{
    std::cerr << "Sorting " << reference_.size() << " " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << std::endl;
    const clock_t start = clock();
    std::sort(reference_.begin(), reference_.end(), &compareKmer<KmerT>);
    std::cerr << "Sorting " << reference_.size() << " " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << " done in " << (clock() - start) / 1000 << "ms" << std::endl;
}

template <typename KmerT>
void ReferenceSorter<KmerT>::saveReference(
    const std::vector<unsigned long> &contigOffsets,
    const std::vector<bool> &neighbors,
    const unsigned long genomeLength)
{
    std::cerr << "Saving " << reference_.size() << " " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers" << std::endl;
    const clock_t start = clock();

    std::ofstream os(outputFile_.c_str());
    if (!os)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno,"Failed to create file " + outputFile_.string()));
    }

    typename std::vector<ReferenceKmer<KmerT> >::iterator current(reference_.begin());
    std::size_t neighborKmers = 0;
    std::size_t storedKmers = 0;
    while(reference_.end() != current)
    {
        std::pair<typename std::vector<ReferenceKmer<KmerT> >::iterator,
                  typename std::vector<ReferenceKmer<KmerT> >::iterator> sameKmerRange =
                std::equal_range(current, reference_.end(), *current, &compareKmer<KmerT>);
        const std::size_t kmerMatches = std::distance(sameKmerRange.first, sameKmerRange.second);

        // the kmers we want to store are those that don't have the neighbors flag set by loadReference.
        const typename std::vector<ReferenceKmer<KmerT> >::const_iterator firstToStore =
            std::find_if(sameKmerRange.first, sameKmerRange.second,
                         boost::bind(&ReferenceKmer<KmerT>::hasNoNeighbors, _1));

        if (firstToStore != sameKmerRange.second)
        {
            if (repeatThreshold_ < kmerMatches)
            {
                std::cerr << "Skipping kmer " << oligo::bases(current->getKmer()) << " as it generates " << kmerMatches << "matches\n";

                static const ReferencePosition tooManyMatchPosition(ReferencePosition::TooManyMatch);
                const ReferenceKmer<KmerT> tooManyMatchKmer(sameKmerRange.first->getKmer(), tooManyMatchPosition);
                //std::cerr << std::hex << referenceKmer.first << '\t' << referenceKmer.second << '\n';
                if (!os.write(reinterpret_cast<const char*>(&tooManyMatchKmer), sizeof(tooManyMatchKmer)))
                {
                    BOOST_THROW_EXCEPTION(common::IoException(errno,"Failed to write toomanymatch reference kmer into " + outputFile_.string()));
                }
                ++storedKmers;
            }
            else
            {
                const bool kmerHasNeighbors =
                    // neighborhood annotation is available
                    !neighbors.empty() &&
                    neighbors.at(contigOffsets.at(firstToStore->getReferencePosition().getContigId()) +
                                 firstToStore->getReferencePosition().getPosition());

                BOOST_FOREACH(ReferenceKmer<KmerT> referenceKmer, sameKmerRange)
                {
                    // the kmers we want to store are those that don't have the neighbors flag set by loadReference.
                    if (referenceKmer.hasNoNeighbors())
                    {
                        if (kmerHasNeighbors)
                        {
                            referenceKmer.setNeighbors();
                        }
                        if (!os.write(reinterpret_cast<const char*>(&referenceKmer), sizeof(referenceKmer)))
                        {
                            BOOST_THROW_EXCEPTION(common::IoException(errno,"Failed to write reference kmer into " + outputFile_.string()));
                        }
                        ++storedKmers;
                    }
                }
            }

        }

        current = sameKmerRange.second;
    }
    os.flush();
    os.close();
    std::cerr << "Saving " << storedKmers << " " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers with " <<
        neighborKmers << " neighbors done in " << (clock() - start) / 1000 << "ms" << std::endl;

    SortedReferenceMetadata sortedReference;
    sortedReference.addMaskFile(oligo::KmerTraits<KmerT>::KMER_BASES, maskWidth_, mask_, outputFile_, storedKmers);
    saveSortedReferenceXml(std::cout, sortedReference);
}

template class ReferenceSorter<oligo::ShortKmerType>;
template class ReferenceSorter<oligo::KmerType>;
template class ReferenceSorter<oligo::LongKmerType>;

} // namespace reference
} // namespace isaac
