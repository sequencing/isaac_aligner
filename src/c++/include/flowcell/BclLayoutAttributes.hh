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
 ** \file BclLayoutAttributes.hh
 **
 ** Attribute tag definitions common to all bcl flowcells
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_BCL_LAYOUT_ATTRIBUTES_HH
#define iSAAC_FLOWCELL_BCL_LAYOUT_ATTRIBUTES_HH

#include <boost/filesystem.hpp>

namespace isaac
{
namespace flowcell
{

struct BclFilePathAttributeTag
{
    typedef boost::filesystem::path value_type;
    friend std::ostream &operator << (std::ostream &os, const BclFilePathAttributeTag &tag){return os << "BclFilePathAttributeTag";}
};
struct FiltersFilePathAttributeTag
{
    typedef boost::filesystem::path value_type;
    friend std::ostream &operator << (std::ostream &os, const FiltersFilePathAttributeTag &tag){return os << "FiltersFilePathAttributeTag";}
};
struct PositionsFilePathAttributeTag
{
    typedef boost::filesystem::path value_type;
    friend std::ostream &operator << (std::ostream &os, const PositionsFilePathAttributeTag &tag){return os << "PositionsFilePathAttributeTag";}
};

inline bool isClocsPath(const boost::filesystem::path& path)
{
    static const char dotClocs[] = {".clocs"};
    static const size_t dotClocsLength = sizeof(dotClocs) - 1;
    return path.string().length() > dotClocsLength &&
        0 == path.string().compare(path.string().size() - dotClocsLength, dotClocsLength, dotClocs);
}

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_BCL_LAYOUT_ATTRIBUTES_HH
