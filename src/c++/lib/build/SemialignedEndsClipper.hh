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
#include "build/PackedFragmentBuffer.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace build
{

class SemialignedEndsClipper : boost::noncopyable
{
    static const unsigned CONSECUTIVE_MATCHES_MIN = 5;
public:
    SemialignedEndsClipper(alignment::Cigar &cigarBuffer) : cigarBuffer_(cigarBuffer)
    {
    }

    void clip(
        const std::vector<reference::Contig> &contigs,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment);

private:
    alignment::Cigar &cigarBuffer_;

    void clipLeftSide(
        const std::vector<reference::Contig> &contigList,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment);

    void clipRightSide(
        const std::vector<reference::Contig> &contigList,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment);

};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEMIALIGNED_ENDS_CLIPPER_HH
