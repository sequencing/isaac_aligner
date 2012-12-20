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
 ** \file SemialignedEndsClipper.hh
 **
 ** \brief Utility classes for detecting and removing fragment parts that contain
 **        sequences of the adapters
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_SEMIALIGNED_ENDS_CLIPPER_HH
#define iSAAC_ALIGNMENT_SEMIALIGNED_ENDS_CLIPPER_HH

#include "alignment/BamTemplate.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class SemialignedEndsClipper
{
    static const unsigned CONSECUTIVE_MATCHES_MIN = 5;
public:
    SemialignedEndsClipper()
    {
        reserve();
    }

    // this is needed to keep clippers in a vector. There is no case of copying them with state preservation.
    SemialignedEndsClipper(const SemialignedEndsClipper &that)
    {
        reserve();
    }

    void clip(
        const std::vector<reference::Contig> &barcodeContigList,
        BamTemplate &bamTemplate);

    bool clip(
        const std::vector<reference::Contig> &barcodeContigList,
        FragmentMetadata &fragmentMetadata);

    void reset()
    {
        cigarBuffer_.clear();
    }

    void reserve()
    {
        // should be enough for two reads
        cigarBuffer_.reserve(10000);
    }

    /// storing SemialignedEndsClipper objects in the vector requires this operator although it is not expected to be executed at runtime
    SemialignedEndsClipper &operator=(const SemialignedEndsClipper &read)
    {
        ISAAC_ASSERT_MSG(false, "SemialignedEndsClipper objects are not supposed to be reassigned");
        return *this;
    }

private:
    Cigar cigarBuffer_;

    bool clipLeftSide(
        const std::vector<reference::Contig> &contigList,
        FragmentMetadata &fragmentMetadata);

    void clipRightSide(
        const std::vector<reference::Contig> &contigList,
        FragmentMetadata &fragmentMetadata);

};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEMIALIGNED_ENDS_CLIPPER_HH
