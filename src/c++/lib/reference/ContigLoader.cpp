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
 ** \file ContigLoader.cpp
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "reference/ContigLoader.hh"

namespace isaac
{
namespace reference
{

void loadContig(
    const reference::SortedReferenceMetadata::Contig &xmlContig,
    std::vector<char> &forward)
{
    forward.clear();
    forward.reserve(xmlContig.totalBases_);
    std::ifstream is(xmlContig.filePath_.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open reference file " + xmlContig.filePath_.string()));
    }
    if (!is.seekg(xmlContig.offset_))
    {
        using common::IoException;
        using boost::format;
        const format message = (boost::format("Failed to reach offset %d in reference file % s") % xmlContig.offset_ % xmlContig.filePath_);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
//        ISAAC_THREAD_CERR << (boost::format("Contig seek %s (%3d:%8d): %s") % xmlContig.name_ % xmlContig.index_ % xmlContig.totalBases_ % xmlContig.filePath_).str() << std::endl;
    static const oligo::Translator translator = oligo::getTranslator(true, oligo::invalidOligo);
    char base = 0;
    while(is && (forward.size() < xmlContig.totalBases_) && is.get(base))
    {
        ISAAC_ASSERT_MSG(base >= 0, "Valid fasta characters can't be negative");
        ISAAC_ASSERT_MSG(base < static_cast<long>(translator.size()), "Translator is expected to handle 256 values");
        if (std::isalpha(base))
        {
            forward.push_back(oligo::getBase(translator[base], true));
        }
    }
    if (xmlContig.totalBases_ != forward.size())
    {
        using common::IoException;
        using boost::format;
        const format message = (format("Failed to read %d bases from reference file % s: %d") % xmlContig.totalBases_ % xmlContig.filePath_.string() % forward.size());
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
}

struct DummyFilter {bool operator() (const unsigned contigIdx) const {return true;}} dummyFilter;
/**
 * \brief loads all the fasta file contigs into memory on multiple threads
 */
std::vector<reference::Contig> loadContigs(
    const reference::SortedReferenceMetadata::Contigs &xmlContigs,
    common::ThreadVector &loadThreads)
{
    return loadContigs(xmlContigs, dummyFilter, loadThreads);
}


} // namespace reference
} // namespace isaac
