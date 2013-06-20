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
 ** \file OverlappingEndsClipper.hh
 **
 ** \brief Utility classes for detecting and removing overlapping parts of the reads when template is so short
 **        that the bit in the middle gets sequenced twice.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_OVERLAPPING_ENDS_CLIPPER_HH
#define iSAAC_ALIGNMENT_OVERLAPPING_ENDS_CLIPPER_HH

#include "alignment/BamTemplate.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class OverlappingEndsClipper
{
public:
    OverlappingEndsClipper()
    {
        reserve();
    }

    // this is needed to keep clippers in a vector. There is no case of copying them with state preservation.
    OverlappingEndsClipper(const OverlappingEndsClipper &that)
    {
        reserve();
    }

    void clip(
        const std::vector<reference::Contig> &contigList,
        BamTemplate &bamTemplate);

    void reset()
    {
        cigarBuffer_.clear();
    }

    /// storing OverlappingEndsClipper objects in the vector requires this operator although it is not expected to be executed at runtime
    OverlappingEndsClipper &operator=(const OverlappingEndsClipper &read)
    {
        ISAAC_ASSERT_MSG(false, "OverlappingEndsClipper objects are not supposed to be reassigned");
        return *this;
    }

private:
    Cigar cigarBuffer_;

    int diffBaseQualities(
        std::vector<char>::const_iterator left,
        std::vector<char>::const_iterator right,
        unsigned length);

    void reserve()
    {
        // should be enough for two reads
        cigarBuffer_.reserve(10000);
    }

};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_OVERLAPPING_ENDS_CLIPPER_HH
