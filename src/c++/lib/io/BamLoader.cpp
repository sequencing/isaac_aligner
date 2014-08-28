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
 ** \file BamLoader.cpp
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#include "common/Debug.hh"
#include "io/BamLoader.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace io
{

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BamLoader::BamLoader(
    std::size_t maxPathLength,
    common::ThreadVector &threads,
    const unsigned coresMax) :
    bgzfReader_(threads, coresMax),
    lastUnparsedBytes_(0),
    decompressParseParallelizationThreads_(2),
    nextDecompressorThread_(0),
    nextParserThread_(0)
{
    reserveBuffers(maxPathLength);
}

void BamLoader::reserveBuffers(std::size_t maxPathLength)
{
    // give more than one thread a chance to unpack a bgzf block. Don't give too much as then the L3 cache
    // gets trashed and things take longer.
    // IMPORTANT!!!, don't use coresMax in this calculation. MatchFinder and MatchSelector are likely to
    // provide different numbers there. If the buffer sizes differ between processing stages, the
    // read pairing will be different between two attempts of reading the same tile. Use the constant
    // that will result in the same buffer size through the run.
    const std::size_t bufferSize = UNPARSED_BYTES_MAX + UNCOMPRESSED_BGZF_BLOCK_SIZE * BGZF_BLOCKS_PER_CLUSTER_BLOCK;
    lastPassBam_.reserve(bufferSize);
    decompressionBuffers_[0].reserve(bufferSize);
    decompressionBuffers_[1].reserve(bufferSize);
}



} // namespace io
} // namespace isaac
