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
 ** \file BamLayout.hh
 **
 ** Specialization of Layout for bam flowcell.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_BAM_LAYOUT_HH
#define iSAAC_FLOWCELL_BAM_LAYOUT_HH

#include "flowcell/Layout.hh"

namespace isaac
{
namespace flowcell
{
struct BamFilePathAttributeTag
{
    typedef boost::filesystem::path value_type;
    friend std::ostream &operator << (std::ostream &os, const BamFilePathAttributeTag &tag){return os << "BamFilePathAttributeTag";}
};

template<>
const boost::filesystem::path & Layout::getAttribute<Layout::Bam, BamFilePathAttributeTag>(
    boost::filesystem::path &result) const;

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_BAM_LAYOUT_HH
