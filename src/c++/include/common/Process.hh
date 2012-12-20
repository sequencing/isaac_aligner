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
 ** \file Process.hh
 **
 ** process management helper utilities.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_PROCESS_HPP
#define iSAAC_COMMON_PROCESS_HPP

#include "common/Debug.hh"

namespace isaac
{
namespace common
{

void executeCommand(const std::string &cmd);

} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_PROCESS_HPP
