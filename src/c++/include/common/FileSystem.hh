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
 ** \file FileSystem.hh
 **
 ** file system helper utilities.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_FILE_SYSTEM_HH
#define iSAAC_COMMON_FILE_SYSTEM_HH

#include <vector>
#include <boost/filesystem/path.hpp>

#include "common/Debug.hh"

namespace isaac
{
namespace common
{

void createDirectories(const std::vector<boost::filesystem::path> &createList);


inline bool isDotGzPath(const boost::filesystem::path& path)
{
    static const char dotGz[] = {".gz"};
    static const size_t dotGzLength = sizeof(dotGz) - 1;
    return path.string().length() > dotGzLength &&
        0 == path.string().compare(path.string().size() - dotGzLength, dotGzLength, dotGz);
}

} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_FILE_SYSTEM_HH
