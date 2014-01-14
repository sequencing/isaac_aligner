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
