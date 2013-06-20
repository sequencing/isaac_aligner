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
 ** \file BitsetLoader.hh
 **
 ** \brief Helper for loading neighbor flags and such from a binary file
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_IO_BITSET_LOADER_HH
#define ISAAC_IO_BITSET_LOADER_HH

#include <fstream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>

namespace isaac
{
namespace io
{

class BitsetLoader: boost::noncopyable
{
    boost::filesystem::path filePath_;
    std::ifstream is_;
public:
    BitsetLoader(const boost::filesystem::path &filePath);

    std::size_t load(
        const unsigned long genomeSize,
        std::vector<bool> &bits);

    std::size_t load(
        std::vector<bool> &bits);
};

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_IO_BITSET_LOADER_HH
