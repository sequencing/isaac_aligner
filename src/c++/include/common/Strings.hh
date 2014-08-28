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

#ifndef iSAAC_COMMON_STRINGS_HH
#define iSAAC_COMMON_STRINGS_HH

#include <string>

namespace isaac
{
namespace common
{

inline void replaceSubstring(std::string& str, const std::string& substr, const std::string& newstr)
{
  size_t pos = 0;
  while((pos = str.find(substr, pos)) != std::string::npos)
  {
     str.replace(pos, substr.length(), newstr);
     pos += newstr.length();
  }
}


} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_STRINGS_HH
