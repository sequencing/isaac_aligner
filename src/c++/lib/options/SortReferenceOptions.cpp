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
 ** \file SortReferenceOptions.cpp
 **
 ** Command line options for 'sortReference'
 **
 ** \author Come Raczy
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "options/SortReferenceOptions.hh"

namespace isaac
{
namespace options
{

namespace bpo = boost::program_options;

SortReferenceOptions::SortReferenceOptions()
    : seedLength(32)
    , maskWidth(6)
    , mask(0)
    , repeatThreshold(1000)
{
     namedOptions_.add_options()
        ("genome-file,g",       bpo::value<std::string>(&genomeFile),
                                "Path to the reference genome")
        ("genome-neighbors,n",  bpo::value<boost::filesystem::path>(&genomeNeighborsFile),
                                "Path to the file containing neighbor flags (one bit per genome file position)")
        ("mask,m",              bpo::value<unsigned long>(&mask),
                                "mask used to filter the k-mers counted by this process (must be strictly less than 2^mask-width")
        ("mask-width,w",        bpo::value<unsigned int>(&maskWidth)->default_value(maskWidth),
                                "Width in bits of the mask used to split the sorted files")
        ("repeat-threshold",    bpo::value<unsigned int>(&repeatThreshold)->default_value(repeatThreshold),
                                "Maximum number of k-mer occurrences in genome for it to be counted as repeat")
        ("output-file,o",       bpo::value<boost::filesystem::path>(&outFile), "Output file path.")
        ("seed-length,s",       bpo::value<unsigned int>(&seedLength)->default_value(seedLength),
                                "Length of reference k-mer in bases. 64 or 32 is supported.")
        ;
}

void SortReferenceOptions::postProcess(bpo::variables_map &vm)
{
    if(vm.count("help"))
    {
        return;
    }
    using isaac::common::InvalidOptionException;
    using boost::format;
    const std::vector<std::string> requiredOptions = boost::assign::list_of("mask")("genome-file")("output-file");
    BOOST_FOREACH(const std::string &required, requiredOptions)
    {
        if(!vm.count(required))
        {
            const format message = format("\n   *** The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
    const unsigned int maskCount = (1 << maskWidth);
    if(maskCount <= mask)
    {
        const format message = format("\n   *** The mask must be strictly less than %d: mask = %d ***\n") % maskCount % mask;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }

    if (16 != seedLength && 32 != seedLength && 64 != seedLength)
    {
        const format message = format("\n   *** The seed-length must be either 16, 32 or 64. Got: %d ***\n") % seedLength;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}

} //namespace option
} // namespace isaac
