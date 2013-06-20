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
 ** \file BamLoader.cpp
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>

#include "common/Debug.hh"
#include "bam/BamParser.hh"

namespace isaac
{
namespace bam
{

bool BamParser::skipHeader(
    std::vector<char>::const_iterator &it,
    std::vector<char>::const_iterator end)
{
    if (-1U == headerBytesToSkip_)
    {
        char magic[4] = {'B', 'A', 'M', '\1'};
        if (std::size_t(std::distance(it, end)) < sizeof(magic) ||
            it + sizeof(magic) != std::mismatch(it, it + sizeof(magic), magic).first)
        {
            BOOST_THROW_EXCEPTION(BamParserException("Bam magic first 4 bytes of the data are not present "));
        }

        it += sizeof(magic);
        unsigned l_text = 0;
        if (std::size_t(std::distance(it, end)) < sizeof(l_text))
        {
            BOOST_THROW_EXCEPTION(BamParserException("Bam header l_text is not present"));
        }
        it = common::extractLittleEndian(it, l_text);

        headerBytesToSkip_ = l_text;
    }
    const unsigned increment = std::min<unsigned>(headerBytesToSkip_, std::distance(it, end));
    it += increment;
    headerBytesToSkip_ -= increment;
    if (!headerBytesToSkip_)
    {
        ISAAC_THREAD_CERR << "BAM Header skipped" << std::endl;
    }
    return true;
}

bool BamParser::skipReferences(
    std::vector<char>::const_iterator &it,
    std::vector<char>::const_iterator end)
{
    if (end == it)
    {
        return true;
    }

    if (-1 == referenceSequencesToSkip_)
    {
        int n_ref = 0;

        if (std::size_t(std::distance(it, end)) < sizeof(n_ref))
        {
            return true;
        }
        it = common::extractLittleEndian(it, n_ref);
        referenceSequencesToSkip_ = n_ref;
    }

    ISAAC_THREAD_CERR << "Will need to skip " << referenceSequencesToSkip_ << " reference sequences" << std::endl;

    while (referenceSequencesToSkip_)
    {
        int l_name = 0;
        if (std::size_t(std::distance(it, end)) < sizeof(l_name))
        {
            return true;
        }

        const std::vector<char>::const_iterator itRef = it;
        it = common::extractLittleEndian(it, l_name);
        if (std::size_t(std::distance(it, end)) < l_name + sizeof(int))
        {
            it = itRef;
            return true;
        }
//        ISAAC_THREAD_CERR << "SQ:" << std::string(it, it + l_name)  << std::endl;

        it += l_name + sizeof(int);
        --referenceSequencesToSkip_;
    }

    if (!referenceSequencesToSkip_)
    {
        ISAAC_THREAD_CERR << "BAM References skipped" << std::endl;
    }

    return true;
}


} // namespace io
} // namespace isaac
