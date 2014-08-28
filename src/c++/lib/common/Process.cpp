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
