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
 ** \file BitsetSaver.cpp
 **
 ** Helper for saving neighbor flags and such from a binary file
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/BitsetSaver.hh"

namespace isaac
{
namespace io
{

BitsetSaver::BitsetSaver(const boost::filesystem::path filePath) :
    filePath_(filePath),
    os_(filePath.c_str())
{
    if (!os_)
    {
        const boost::format message = boost::format("Failed to open file %s for writing: %s") % filePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}

void BitsetSaver::save(
    const std::vector<bool> &bits)
{
    char byte = 0;
    for (size_t i = 0; i < bits.size(); ++ i)
    {
        const unsigned shift = i % 8;
        const char positionHasNeighbors = bits[i];
        byte |= positionHasNeighbors << shift;
        if (7 == shift)
        {
            if (!os_.write(&byte, sizeof(byte)))
            {
                const boost::format message = boost::format("Failed to write bits into %s: %s") % filePath_.string() % strerror(errno);
                BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
            }
            byte = 0;
        }
    }

    if (!os_.write(&byte, sizeof(byte)))
    {
        const boost::format message = boost::format("Failed to write final bits byte into %s: %s") % filePath_.string() % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}

} // namespace reference
} // namespace isaac
