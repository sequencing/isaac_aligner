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
 ** \file ReferenceSorter.hh
 **
 ** Top level component to produce a sorted reference.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_SORTER_HH
#define iSAAC_REFERENCE_REFERENCE_SORTER_HH

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>

#include "oligo/Kmer.hh"
#include "oligo/Permutations.hh"
#include "reference/ReferenceKmer.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace reference
{

class ReferenceSorter: boost::noncopyable
{
public:
    ReferenceSorter(
        const unsigned int maskWidth,
        const oligo::Kmer mask,
        const boost::filesystem::path &genomeFile,
        const boost::filesystem::path &genomeNeighborsFile,
        const std::string &permutationName,
        const boost::filesystem::path &outputFile,
        const unsigned repeatThreshold);
    void run();
private:
    const unsigned int maskWidth_;
    const oligo::Kmer mask_;

    // the mask highlight bits in original kmer (ABCD)
    const oligo::Kmer msbMask_;
    // the mask value in the original kmer (ABCD) shifted to the topmost position
    const oligo::Kmer maskBits_;

    // the mask highlight bits as if the permutated kmer having them all set was unpermutated back into ABCD
    const oligo::Kmer unpermutatedMsbMask_;
    // the mask as if the permutated kemer containing it was unpermutated back into ABCD
    const oligo::Kmer unpermutatedMaskBits_;

    const oligo::Permutation permutation_;

    const boost::filesystem::path genomeFile_;
    const boost::filesystem::path genomeNeighborsFile_;
    const std::string permutationName_;

    const boost::filesystem::path outputFile_;
    std::vector<ReferenceKmer> reference_;

    const unsigned repeatThreshold_;

    std::vector<unsigned long> loadReference();
    void sortReference();
    void saveReference(
        const std::vector<unsigned long> &contigOffsets,
        const std::vector<bool> &neighbors);
    void addToReference(const oligo::Kmer kmer, const ReferencePosition &referencePosition);
};

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERNECE_SORTER_HH
