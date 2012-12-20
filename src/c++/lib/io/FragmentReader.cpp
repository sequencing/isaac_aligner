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
 ** \file FragmentReader.cpp
 **
 ** \brief See FragmentReader.hh
 ** 
 ** \author Come Raczy
 **/

#include <cerrno>
#include <boost/format.hpp>

#include "io/FragmentReader.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace io
{

FragmentReader::FragmentReader(const bfs::path filePath)
    : std::ifstream(filePath.string().c_str())
    , filePath_(filePath)
{
    if (!*this)
    {
        using common::IoException;
        using boost::format;
        const format message = format("Failed to open %s") % filePath;
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
}

FragmentReader::~FragmentReader()
{
}

std::istream &FragmentReader::read(char *buffer)
{
    unsigned *dataLength = reinterpret_cast<unsigned*>(buffer);
    if (std::istream::read(buffer, sizeof(*dataLength)))
    {
        std::istream::read(buffer + sizeof(*dataLength), *dataLength - sizeof(*dataLength));
    }
    return *this;
}

unsigned FragmentReader::getDataLength(char *buffer)
{
    return *reinterpret_cast<unsigned*>(buffer);
}
reference::ReferencePosition FragmentReader::getTemplatePosition(char *buffer)
{
    return reference::ReferencePosition(*reinterpret_cast<unsigned long*>(buffer + 4));
}

//unsigned getFragmentLength(char *buffer);
//unsigned short getReadLength();
//unsigned short getCigarLength();

} // namespace io
} // namespace isaac
