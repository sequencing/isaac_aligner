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
 ** \file ReferenceSorter.cpp
 **
 ** Top level component to produce a sorted reference.
 **
 ** \author Come Raczy
 **/

#include <boost/assert.hpp>
#include <boost/io/ios_state.hpp>

#include "io/FastaReader.hh"
#include "oligo/Nucleotides.hh"
#include "oligo/Mask.hh"
#include "oligo/Permutations.hh"
#include "reference/ReferenceSorter.hh"
#include "reference/SortedReferenceXml.hh"
#include "common/SystemCompatibility.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace reference
{

ReferenceSorter::ReferenceSorter (
    const unsigned int maskWidth,
    const isaac::oligo::Kmer mask,
    const boost::filesystem::path &genomeFile,
    const boost::filesystem::path &genomeNeighborsFile,
    const std::string &permutationName,
    const boost::filesystem::path &outputFile,
    const unsigned repeatThreshold
    )
    : maskWidth_(maskWidth)
    , mask_(mask)
    , msbMask_(~((~(oligo::Kmer(0)) >> maskWidth)))
    , maskBits_(mask_ << (oligo::kmerLength * 2 - maskWidth))
    , unpermutatedMsbMask_(oligo::getReversePermutation(permutationName)(msbMask_))
    , unpermutatedMaskBits_(oligo::getReversePermutation(permutationName)(maskBits_))
    , permutation_(oligo::getPermutation(permutationName))
    , genomeFile_(genomeFile)
    , genomeNeighborsFile_(genomeNeighborsFile)
    , permutationName_(permutationName)
    , outputFile_(outputFile)
    , repeatThreshold_(repeatThreshold)
{
    boost::io::ios_flags_saver svr(std::cerr);
    std::cerr <<
            "Constructing ReferenceSorter: for " << isaac::oligo::kmerLength << "-mers " <<
            " mask width: " << maskWidth_ <<
            " msbMask_: 0x" << std::hex << msbMask_ <<
            " maskBits_: 0x" << std::hex << maskBits_ <<
            " unpermutatedMsbMask_: 0x" << std::hex << unpermutatedMsbMask_ <<
            " unpermutatedMaskBits_: 0x" << std::hex << unpermutatedMaskBits_ <<
            " genomeFile_: " << genomeFile_ <<
            " permutationName_: " << permutationName_ <<
            " outputFile_: " << outputFile_ <<
            std::endl;

    BOOST_ASSERT(mask_ < isaac::oligo::getMaskCount(maskWidth_) && "Mask value cannot exceed the allowed bit width");
}



void ReferenceSorter::run()
{
    const std::vector<unsigned long> contigOffsets = loadReference();
    sortReference();

    std::vector<bool> neighbors;
    if (!genomeNeighborsFile_.empty())
    {
        const unsigned genomeSize = boost::filesystem::file_size(genomeNeighborsFile_) * 8;
        neighbors.reserve(genomeSize);

        std::ifstream genomeNeighborsInput(genomeNeighborsFile_.c_str());
        if (!genomeNeighborsInput)
        {
            const boost::format message = boost::format("Failed to open genome neighbors file %s for reading: %s") % genomeNeighborsFile_ % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }

        unsigned genomeNeighbors = 0;
        while(genomeNeighborsInput)
        {
            char flags = 0;
            if (genomeNeighborsInput.read(reinterpret_cast<char *>(&flags), sizeof(flags)))
            {
                for (int i = 0; i < 8; ++i)
                {
                    const bool positionHasNeighbors = (flags >> i) & 0x01;
                    neighbors.push_back(positionHasNeighbors);
                    genomeNeighbors += positionHasNeighbors;
                }
            }
        }

        ISAAC_ASSERT_MSG(neighbors.size() == genomeSize, "Incorrect number of flags stored");

        if (!genomeNeighborsInput.eof())
        {
            const boost::format message = boost::format("Failed to scan %s to the end") % genomeNeighborsFile_ % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        else
        {
            ISAAC_THREAD_CERR << "Scanning " << genomeNeighborsFile_ << " found " << genomeNeighbors << " neighbors " << std::endl;
        }
    }
    saveReference(contigOffsets, neighbors);
}

/**
 * \brief Load kmers matching the unpermutatedMaskBits_.
 *
 * \return vector of contig base offsets in the order found in fasta file
 */
std::vector<unsigned long> ReferenceSorter::loadReference()
{
    std::vector<unsigned long> retContigOffsets;
    std::cerr << "Loading " << isaac::oligo::kmerLength << "-mers" << std::endl;
    const clock_t start = clock();

    isaac::io::MultiFastaReader multiFastaReader(std::vector<boost::filesystem::path>(1, genomeFile_));
    char base;
    oligo::Kmer forward = 0;
    //Pair::first_type reverse = 0;
    unsigned int badKmer = isaac::oligo::kmerLength;
    unsigned long position = 0;
    bool newContig = false;
    unsigned long lastContigOffset = 0UL;
    while (multiFastaReader.get(base, newContig))
    {
        if (newContig)
        {
            std::cerr << "New contig: " << multiFastaReader.getContigId() << " found at position " << position << std::endl;
            position = 0;
            badKmer = oligo::kmerLength;
            retContigOffsets.push_back(lastContigOffset);
        }
        if (badKmer)
        {
            --badKmer;
        }
        const oligo::Kmer baseValue = oligo::getValue(base);
        if (baseValue >> 2)
        {
            badKmer = oligo::kmerLength;
        }
        forward <<= 2;
        forward |= baseValue;
        //reverse >>= 2;
        //reverse |= (((~baseValue) & 3) << (2 * kmerLength - 2));
        //std:: cerr << forward << '\t' << position << '\t' << reverse << '\t' << position << '\t' << bad16mer << '\t' << base << '\t' << baseValue << '\n';
        if (0 == badKmer)
        {
            assert(position + 1 >= oligo::kmerLength);
            const unsigned long kmerPosition = (position + 1) - oligo::kmerLength;
            //std::cerr << (boost::format("%2d: %5d: %s") % multiFastaReader.getContigId() % kmerPosition % oligo::bases(forward)).str() << std::endl;
            addToReference(forward, ReferencePosition(multiFastaReader.getContigId(), kmerPosition, false));
            //addToReference(reverse, -position);
            //if (!os)
            //{
            //    std::cerr << "Failed to write 16-mers for position " << position << std::endl;
            //    exit(1);
            //}
        }
        ++position;
        ++lastContigOffset;
    }
    std::cerr << "Loading " << isaac::oligo::kmerLength << "-mers" << " done in " << (clock() - start) / 1000 << "ms" << std::endl;
    return retContigOffsets;
}

void ReferenceSorter::addToReference(const oligo::Kmer kmer, const ReferencePosition &referencePosition)
{
    if (unpermutatedMaskBits_ == (kmer & unpermutatedMsbMask_))
    {
//        std::cerr << std::hex << kmer << '\t' << referencePosition << '\n';
        reference_.push_back(ReferenceKmer(permutation_(kmer), referencePosition));
    }
}

void ReferenceSorter::sortReference()
{
    std::cerr << "Sorting " << reference_.size() << " " << isaac::oligo::kmerLength << "-mers" << std::endl;
    const clock_t start = clock();
    std::sort(reference_.begin(), reference_.end(), compareKmer);
    std::cerr << "Sorting " << reference_.size() << " " << isaac::oligo::kmerLength << "-mers" << " done in " << (clock() - start) / 1000 << "ms" << std::endl;
}

void ReferenceSorter::saveReference(
    const std::vector<unsigned long> &contigOffsets,
    const std::vector<bool> &neighbors)
{
    std::cerr << "Saving " << reference_.size() << " " << isaac::oligo::kmerLength << "-mers" << std::endl;
    const clock_t start = clock();

    std::ofstream os(outputFile_.c_str());
    if (!os)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno,"Failed to create file " + outputFile_.string()));
    }

    std::vector<ReferenceKmer>::iterator current(reference_.begin());
    size_t skippedKmers(0);
    unsigned maxPrefixRangeCount = 0;
    unsigned currentPrefixRangecount = 0;
    oligo::Kmer currentPrefix = 0;
    const unsigned suffixBits = oligo::kmerBitLength / 2;
    unsigned neighborKmers = 0;
    while(reference_.end() != current)
    {
        if (currentPrefix != (current->getKmer() >> suffixBits))
        {
            maxPrefixRangeCount = std::max(maxPrefixRangeCount, currentPrefixRangecount);
            currentPrefixRangecount = 0;
            currentPrefix = current->getKmer() >> suffixBits;
        }
        std::pair<std::vector<ReferenceKmer>::iterator, std::vector<ReferenceKmer>::iterator> sameKmerRange =
                std::equal_range(current, reference_.end(), *current, compareKmer);
        size_t kmerMatches(sameKmerRange.second - sameKmerRange.first);

        if (!neighbors.empty())
        {
            const reference::ReferencePosition pos = sameKmerRange.first->getReferencePosition();
            const bool kmerHasNeighbors = neighbors.at(contigOffsets.at(pos.getContigId()) + pos.getPosition());
            if (kmerHasNeighbors)
            {
                std::for_each(sameKmerRange.first, sameKmerRange.second,
                              boost::bind(&ReferenceKmer::setNeighbors, _1));
                neighborKmers += std::distance(sameKmerRange.first, sameKmerRange.second);
            }
        }
        if (repeatThreshold_ < kmerMatches)
        {
            skippedKmers += kmerMatches;
            std::cerr << "Skipping kmer " << oligo::bases(current->getKmer()) << " as it generates " << kmerMatches << "matches\n";

            ++currentPrefixRangecount;
            static const ReferencePosition tooManyMatchPosition(ReferencePosition::TooManyMatch);
            const ReferenceKmer tooManyMatchKmer(sameKmerRange.first->getKmer(), tooManyMatchPosition);
            //std::cerr << std::hex << referenceKmer.first << '\t' << referenceKmer.second << '\n';
            if (!os.write(reinterpret_cast<const char*>(&tooManyMatchKmer), sizeof(tooManyMatchKmer)))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno,"Failed to write reference kmer into " + outputFile_.string()));
            }
        }
        else
        {
            BOOST_FOREACH(const ReferenceKmer &referenceKmer, sameKmerRange)
            {
                ++currentPrefixRangecount;
                //std::cerr << std::hex << referenceKmer.first << '\t' << referenceKmer.second << '\n';
                if (!os.write(reinterpret_cast<const char*>(&referenceKmer), sizeof(referenceKmer)))
                {
                    BOOST_THROW_EXCEPTION(common::IoException(errno,"Failed to write reference kmer into " + outputFile_.string()));
                }
            }
        }
        current = sameKmerRange.second;
    }
    os.flush();
    os.close();
    std::cerr << "Saving " << reference_.size() - skippedKmers << " " << isaac::oligo::kmerLength << "-mers with " <<
        neighborKmers << " neighbors done in " << (clock() - start) / 1000 << "ms" << std::endl;

    SortedReferenceXml sortedReference;
    maxPrefixRangeCount = std::max(maxPrefixRangeCount, currentPrefixRangecount);
    sortedReference.addMaskFile(permutationName_, maskWidth_, mask_, outputFile_, reference_.size() - skippedKmers, maxPrefixRangeCount);
    using namespace isaac::reference;
    std::cout << sortedReference;
}

} // namespace reference
} // namespace isaac
