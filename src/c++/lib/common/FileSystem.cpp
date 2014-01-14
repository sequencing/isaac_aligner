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
 ** \file FileSystem.cpp
 **
 ** file system helper utilities.
 **
 ** \author Roman Petrovski
 **/
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>

#include "common/Exceptions.hh"
#include "common/FileSystem.hh"

namespace isaac
{
namespace common
{

void createDirectories(std::vector<boost::filesystem::path> createList)
{
    std::sort(createList.begin(), createList.end());
    createList.erase(std::unique(createList.begin(), createList.end()), createList.end());
    BOOST_FOREACH(const boost::filesystem::path &directory, createList)
    {
        if (!boost::filesystem::exists(directory))
        {
            ISAAC_THREAD_CERR << "creating directory " << directory << std::endl;
            boost::system::error_code errorCode;
            if(!boost::filesystem::create_directory(directory, errorCode) &&!exists(directory))
            {
                using common::IoException;
                BOOST_THROW_EXCEPTION(IoException(errorCode.value(), "Failed to create directory " + directory.string()));
            }
        }
    }
}


char makeDirectorySeparatorChar()
{
    boost::filesystem::path slash("/");
    slash.make_preferred();
    ISAAC_ASSERT_MSG(1 == slash.string().size(), "Unexpected directory separator char length greater than 1: " << slash)
    return *slash.native().begin();
}

static const char DIRECTORY_SEPARATOR_CHAR = makeDirectorySeparatorChar();

char getDirectorySeparatorChar()
{
    return  DIRECTORY_SEPARATOR_CHAR;
}


} // namespace common
} // namespace isaac
