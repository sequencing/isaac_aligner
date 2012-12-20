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
    const reference::SortedReferenceXml::Contig &xmlContig,
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
    const std::vector<unsigned int> translator = oligo::getTranslator(true, oligo::invalidOligo);
    char base = 0;
    while(is && (forward.size() < xmlContig.totalBases_) && is.get(base))
    {
        ISAAC_ASSERT_MSG(base >= 0, "Valid fasta characters can't be negative");
        ISAAC_ASSERT_MSG(base < static_cast<long>(translator.size()), "Translator is expected to handle 256 values");
        if (std::isalpha(base))
        {
            forward.push_back(oligo::getBase(translator.at(base), true));
        }
    }
    if (xmlContig.totalBases_ != forward.size())
    {
        using common::IoException;
        using boost::format;
        const format message = (format("Failed to read %d bases from reference file % s: %d") % xmlContig.totalBases_ % xmlContig.filePath_.string() % forward.size());
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
    ISAAC_THREAD_CERR << (boost::format("Contig %s (%3d:%8d): %s\n") % xmlContig.name_ % xmlContig.index_ % xmlContig.totalBases_ % xmlContig.filePath_).str();
    static const std::vector<char>::size_type maxBasesToPrintFromEachEnd(35);
    ISAAC_THREAD_CERR
        << std::string(forward.begin(), forward.begin() + std::min(forward.size()/2, maxBasesToPrintFromEachEnd)) +
            (forward.size() <= maxBasesToPrintFromEachEnd * 2 ? "" : " ... ") +
            std::string(forward.end() - std::min(forward.size()/2, maxBasesToPrintFromEachEnd), forward.end())
        << std::endl;
}

void loadContigsParallel(
    std::vector<reference::SortedReferenceXml::Contig>::const_iterator &nextContigToLoad,
    const std::vector<reference::SortedReferenceXml::Contig>::const_iterator contigsEnd,
    std::vector<reference::Contig> &contigList,
    boost::mutex &mutex)
{
    boost::lock_guard<boost::mutex> lock(mutex);
    while (contigsEnd != nextContigToLoad)
    {
        const std::vector<reference::SortedReferenceXml::Contig>::const_iterator ourContig = nextContigToLoad++;
        ISAAC_ASSERT_MSG(contigList[ourContig->karyotypeIndex_].index_ == ourContig->karyotypeIndex_, "Unexpected order of preallocated contigs or index collision");
        {
            common::unlock_guard<boost::mutex> unlock(mutex);
            loadContig(*ourContig, contigList[ourContig->karyotypeIndex_].forward_);
        }
        contigList[ourContig->karyotypeIndex_].index_ = ourContig->index_;
    }
}

/**
 * \brief loads the fasta file contigs into memory on multiple threads
 */
std::vector<reference::Contig> loadContigs(
    const reference::SortedReferenceXml &sortedReferenceXml,
    common::ThreadVector &loadThreads)
{
    const std::vector<reference::SortedReferenceXml::Contig> xmlContigs = sortedReferenceXml.getContigs();
    std::vector<reference::Contig> ret;
    ret.reserve(xmlContigs.size());
    BOOST_FOREACH(const reference::SortedReferenceXml::Contig &xmlContig, xmlContigs)
    {
        ISAAC_ASSERT_MSG(ret.size() == xmlContig.index_ , "Expected sequentially ordered starting with 0");

        ret.push_back(reference::Contig(xmlContig.index_, xmlContig.name_));
    }

    std::vector<reference::SortedReferenceXml::Contig>::const_iterator nextContigToLoad = xmlContigs.begin();
    boost::mutex mutex;
    loadThreads.execute(boost::bind(&loadContigsParallel,
                                    boost::ref(nextContigToLoad),
                                    xmlContigs.end(),
                                    boost::ref(ret),
                                    boost::ref(mutex)));

    return ret;
}

/**
 * \brief loads the fasta file contigs into memory on multiple threads
 */
std::vector<std::vector<reference::Contig> > loadContigs(
    const reference::SortedReferenceXmlList &sortedReferenceXmlList,
    common::ThreadVector &loadThreads)
{
    std::vector<std::vector<reference::Contig> > ret(sortedReferenceXmlList.size());

    BOOST_FOREACH(const reference::SortedReferenceXml &sortedReferenceXml, sortedReferenceXmlList)
    {
        std::vector<reference::Contig> contigList = loadContigs(sortedReferenceXml, loadThreads);
        const unsigned referenceIndex = &sortedReferenceXml - &sortedReferenceXmlList.front();
        ret.at(referenceIndex).swap(contigList);
    }
    return ret;
}

} // namespace reference
} // namespace isaac
