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

static const std::vector<std::string> permutationNameList = boost::assign::list_of
    ("ABCD")("ACBD")("ADBC")("BCDA")("BDAC")("CDAB");

SortReferenceOptions::SortReferenceOptions()
    : maskWidth(6)
    , mask(0)
    , repeatThreshold(1000)
{
    std::vector<std::string> permutationList = permutationNameList;
    BOOST_FOREACH(std::string &name, permutationList)
    {
        name = std::string("""") + name +  std::string("""");
    }
    const std::string permutationNamesString = boost::algorithm::join(permutationList, ", ");
    namedOptions_.add_options()
        ("mask-width,w",        bpo::value<unsigned int>(&maskWidth)->default_value(maskWidth),
                                "Width in bits of the mask used to split the sorted files")
        ("mask,m",              bpo::value<unsigned long>(&mask),
                                "mask used to filter the k-mers counted by this process (must be strictly less than 2^mask-width")
        ("repeat-threshold",    bpo::value<unsigned int>(&repeatThreshold)->default_value(repeatThreshold),
                                "Maximum number of kmer occurrences in genome for it to be counted as repeat")
        ("genome-file,g",       bpo::value<std::string>(&genomeFile),
                                "Path to the reference genome")
        ("genome-neighbors,n",  bpo::value<boost::filesystem::path>(&genomeNeighborsFile),
                                "Path to the file containing neighbor flags (one bit per genome file position)")
        ("permutation-name,p",  bpo::value<std::string>(&permutationName),
                                (boost::format("Name of the permutation to apply: %s") %
                                        permutationNamesString).str().c_str())
        ("output-file,o",       bpo::value<boost::filesystem::path>(&outFile), "Output file path.")
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
    const std::vector<std::string> requiredOptions = boost::assign::list_of("mask")("genome-file")("permutation-name")("output-file");
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
    if (permutationNameList.end() == std::find(permutationNameList.begin(), permutationNameList.end(), permutationName))
    {
        const format message = boost::format("unknown permutation name %d") % permutationName;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
}

} //namespace option
} // namespace isaac
