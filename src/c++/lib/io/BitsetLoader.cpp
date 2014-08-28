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
 ** \file BitsetLoader.cpp
 **
 ** Helper for loading neighbor flags and such from a binary file
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/BitsetLoader.hh"

namespace isaac
{
namespace io
{

BitsetLoader::BitsetLoader(
    const boost::filesystem::path &filePath) :
    filePath_(filePath),
    is_(filePath_.c_str())
{
    if (!is_)
    {
        const boost::format message = boost::format("Failed to open bitset file %s for reading: %s") % filePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}

std::size_t BitsetLoader::load(
    const unsigned long genomeSize,
    std::vector<bool> &bits)
{
    unsigned long genomeSizeLeft = genomeSize;
    bits.reserve(genomeSize);
    bits.clear();
    unsigned long bitsSet = 0;
    while(is_)
    {
        char flags = 0;
        if (is_.read(reinterpret_cast<char *>(&flags), sizeof(flags)))
        {
            for (int i = 0; i < 8 && genomeSizeLeft--; ++i)
            {
                const bool positionHasNeighbors = (flags >> i) & 0x01;
                bits.push_back(positionHasNeighbors);
                bitsSet += positionHasNeighbors;
            }
        }
    }

    ISAAC_ASSERT_MSG(bits.size() == genomeSize, "Incorrect number of flags loaded from " << filePath_ << " expected:" << genomeSize << " got: " << bits.size());

    return bitsSet;
}

std::size_t BitsetLoader::load(
    std::vector<bool> &bits)
{
    const unsigned long genomeSize = boost::filesystem::file_size(filePath_) * 8;

    return load(genomeSize, bits);
}


} // namespace reference
} // namespace isaac
