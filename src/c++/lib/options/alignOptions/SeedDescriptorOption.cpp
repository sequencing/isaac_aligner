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
 ** \file SeedDescriptorOption.cpp
 **
 ** seeds option parsing
 **
 ** \author Roman Petrovski
 **/
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/SeedMetadata.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"

#include "SeedDescriptorOption.hh"


namespace isaac
{
namespace options
{
namespace alignOptions
{

alignment::SeedMetadataList parseSeedDescriptor(
    const std::vector<flowcell::ReadMetadata> &readMetadataList,
    const std::string &seedDescriptor)
{
    using boost::algorithm::split;
    using boost::algorithm::is_any_of;
    if (seedDescriptor.empty())
    {
        const boost::format message = boost::format("\n   *** The seed descriptor is empty. At least one seed is needed ***\n");
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    const unsigned SEED_LENGTH = 32;
    alignment::SeedMetadataList seedMetadataList;
    std::vector<std::string> seedDescriptorList; // split by read
    split(seedDescriptorList, seedDescriptor,  is_any_of(","));
    if (readMetadataList.size() < seedDescriptorList.size())
    {
        const boost::format message = boost::format("\n   *** Too many lists-of-seeds in seed-descriptor '%s': found %d: %d reads only ***\n") %
            seedDescriptor % seedDescriptorList.size() % readMetadataList.size();
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    // extend the last list-of-seeds to all subsequent reads if needed
    seedDescriptorList.resize(readMetadataList.size(), seedDescriptorList.back());
    // create all the seeds
    std::vector<flowcell::ReadMetadata>::const_iterator readMetadata = readMetadataList.begin();
    BOOST_FOREACH(const std::string &descriptor, seedDescriptorList)
    {
        const unsigned readIndex = readMetadata->getIndex();
        if (descriptor.empty())
        {
            const boost::format message = boost::format("\n   *** The seed descriptor for read index %d is empty. At least one seed is needed ***\n") %
                readIndex;
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
        std::vector<std::string> offsetList;
        split(offsetList, descriptor,  is_any_of(":"));
        if (offsetList.empty())
        {
            const boost::format message = boost::format("\n   *** The list of seed offsets read index %d is empty. At least one seed is needed for each read ***\n") %
                readIndex;
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
        BOOST_FOREACH(const std::string &offsetString, offsetList)
        {
            try
            {
                const unsigned offset = boost::lexical_cast<unsigned>(offsetString);
                alignment::SeedMetadata seedMetadata(offset, SEED_LENGTH, readIndex, seedMetadataList.size());
                if (offset + SEED_LENGTH > readMetadata->getLength())
                {
                    ISAAC_THREAD_CERR << "WARNING: ignored " << seedMetadata <<
                        " as it stretches beyond the read " << (readIndex + 1) <<
                        " which is " << readMetadata->getLength() << " bases long" << std::endl;
                }
                else
                {
                    seedMetadataList.push_back(seedMetadata);
                    ISAAC_THREAD_CERR << "constructed " << seedMetadataList.back() << std::endl;
                }
            }
            catch(boost::bad_lexical_cast &)
            {
                const boost::format message = boost::format("\n   *** Invalid seed offset '%s' found in '%s' ***\n") %
                    offsetString % seedDescriptor;
                BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
            }
        }
        ++readMetadata;
    }

    return seedMetadataList;
}

} // namespace alignOptions
} // namespace option
} // namespace isaac
