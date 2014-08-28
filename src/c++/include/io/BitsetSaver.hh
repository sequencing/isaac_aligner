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
 ** \file BitsetSaver.hh
 **
 ** \brief Helper for saving neighbor flags and such from a binary file
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_REFERENCE_BITSET_SAVER_HH
#define ISAAC_REFERENCE_BITSET_SAVER_HH

#include <iostream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>

namespace isaac
{
namespace io
{

class BitsetSaver: boost::noncopyable
{
    boost::filesystem::path filePath_;
    std::ofstream os_;
public:
    BitsetSaver(const boost::filesystem::path filePath);
    void save(
        const std::vector<bool> &bits);
};

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_BITSET_SAVER_HH
