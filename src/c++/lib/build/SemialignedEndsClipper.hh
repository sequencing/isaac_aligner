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
#include "io/Fragment.hh"
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

    bool clipLeftSide(
        const std::vector<reference::Contig> &contigList,
        const reference::ReferencePosition binEndPos,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment);

    bool clipRightSide(
        const std::vector<reference::Contig> &contigList,
        PackedFragmentBuffer::Index &index,
        io::FragmentAccessor &fragment);

};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEMIALIGNED_ENDS_CLIPPER_HH
