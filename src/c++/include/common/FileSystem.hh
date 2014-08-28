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

void createDirectories(std::vector<boost::filesystem::path> createList);


inline bool isDotGzPath(const boost::filesystem::path& path)
{
    static const char dotGz[] = {'.', 'g', 'z'};
    return path.string().length() > sizeof(dotGz) &&
        0 == path.string().compare(path.string().size() - sizeof(dotGz), sizeof(dotGz), dotGz);
}

char getDirectorySeparatorChar();

} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_FILE_SYSTEM_HH
