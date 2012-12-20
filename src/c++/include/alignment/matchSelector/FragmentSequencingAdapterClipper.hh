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
 ** \file FragmentSequencingAdapterClipper.hh
 **
 ** \brief Utility classes for detecting and removing fragment parts that contain
 **        sequences of the adapters
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_SEQUENCING_ADAPTER_CLIPPER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_SEQUENCING_ADAPTER_CLIPPER_HH

#include "alignment/matchSelector/SequencingAdapter.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class FragmentSequencingAdapterClipper: public boost::noncopyable
{
    // percent of mismatching bases below which the flank is assumed
    // to be too good for the real adaptor-containing read.
    static const unsigned TOO_GOOD_READ_MISMATCH_PERCENT = 40;

    const matchSelector::SequencingAdapterList &sequencingAdapters_;
public:
    explicit FragmentSequencingAdapterClipper(
        const matchSelector::SequencingAdapterList &sequencingAdapters):
            sequencingAdapters_(sequencingAdapters)
    {
    }

    void checkInitStrand(
        const FragmentMetadata &fragmentMetadata,
        const reference::Contig &contig);

    void clip(
        const reference::Contig &contig,
        FragmentMetadata &fragment,
        std::vector<char>::const_iterator &sequenceBegin,
        std::vector<char>::const_iterator &sequenceEnd) const;
private:
    struct SequencingAdapterRange
    {
        SequencingAdapterRange() : initialized_(false), empty_(true), unbounded_(false){}
        bool initialized_;
        bool empty_;
        bool unbounded_;
        std::vector<char>::const_iterator adapterRangeBegin_;
        std::vector<char>::const_iterator adapterRangeEnd_;
    };

    struct StrandSequencingAdapterRange
    {
        // 0 - forward range, 1 - reverse range
        matchSelector::FragmentSequencingAdapterClipper::SequencingAdapterRange strandRange_[2];
    };

    StrandSequencingAdapterRange strandAdapters_;

    static bool decideWhichSideToClip(
        const reference::Contig &contig,
        const long contigPosition,
        const std::vector<char>::const_iterator sequenceBegin,
        const std::vector<char>::const_iterator sequenceEnd,
        const SequencingAdapterRange &adapterRange,
        bool &clipBackwards);

};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_SEQUENCING_ADAPTER_CLIPPER_HH
