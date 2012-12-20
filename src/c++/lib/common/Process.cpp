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
 ** \file Process.cpp
 **
 ** Process management helper utilities.
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "common/Process.hh"

namespace isaac
{
namespace common
{

void executeCommand(const std::string &cmd)
{
    ISAAC_THREAD_CERR << "Executing : " << cmd << std::endl;
    if (system(cmd.c_str()))
    {
        BOOST_THROW_EXCEPTION(
            common::IsaacException(
                errno, (boost::format("Couldn't execute the following command: %s: %s") % cmd % strerror(errno)).str()));
    }
    ISAAC_THREAD_CERR << "Executing done: " << cmd << std::endl;
}

} // namespace common
} // namespace isaac
