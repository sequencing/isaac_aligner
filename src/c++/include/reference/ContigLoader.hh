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
 ** \file ContigsLoader.hh
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_CONTIGS_PRINTER_HH
#define iSAAC_REFERENCE_CONTIGS_PRINTER_HH

#include <boost/format.hpp>

#include "common/Threads.hpp"
#include "reference/Contig.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{

std::vector<Contig> loadContigs(
    const reference::SortedReferenceMetadata::Contigs &xmlContigs,
    common::ThreadVector &loadThreads);

void loadContig(
    const reference::SortedReferenceMetadata::Contig &xmlContig,
    std::vector<char> &forward);

template <typename ShouldLoadF> void loadContigsParallel(
    ShouldLoadF &shouldLoad,
    std::vector<reference::SortedReferenceMetadata::Contig>::const_iterator &nextContigToLoad,
    const std::vector<reference::SortedReferenceMetadata::Contig>::const_iterator contigsEnd,
    std::vector<reference::Contig> &contigList,
    boost::mutex &mutex)
{
    const unsigned traceStep = pow(10, int(log10((contigList.size() + 99) / 100)));
    boost::lock_guard<boost::mutex> lock(mutex);
    while (contigsEnd != nextContigToLoad)
    {
        const std::vector<reference::SortedReferenceMetadata::Contig>::const_iterator ourContig = nextContigToLoad++;
        ISAAC_ASSERT_MSG(contigList[ourContig->karyotypeIndex_].index_ == ourContig->karyotypeIndex_, "Unexpected order of preallocated contigs or index collision");
        if (shouldLoad(ourContig->karyotypeIndex_))
        {
            common::unlock_guard<boost::mutex> unlock(mutex);
            const reference::SortedReferenceMetadata::Contig &xmlContig = *ourContig;
            std::vector<char> &forward = contigList[ourContig->karyotypeIndex_].forward_;
            loadContig(xmlContig, forward);
            if (!(xmlContig.index_ % traceStep))
            {
                ISAAC_THREAD_CERR << (boost::format("Contig %s (%3d:%8d): %s\n") % xmlContig.name_ % xmlContig.index_ % xmlContig.totalBases_ % xmlContig.filePath_).str();
                static const std::vector<char>::size_type maxBasesToPrintFromEachEnd(35);
                ISAAC_THREAD_CERR
                    << std::string(forward.begin(), forward.begin() + std::min(forward.size()/2, maxBasesToPrintFromEachEnd)) +
                        (forward.size() <= maxBasesToPrintFromEachEnd * 2 ? "" : " ... ") +
                        std::string(forward.end() - std::min(forward.size()/2, maxBasesToPrintFromEachEnd), forward.end())
                    << std::endl;
            }
        }
        contigList[ourContig->karyotypeIndex_].index_ = ourContig->index_;
    }
}

/**
 * \brief loads the fasta file contigs into memory on multiple threads unless shouldLoad(contig->index_) returns false
 */
template <typename ShouldLoadF> std::vector<reference::Contig> loadContigs(
    const reference::SortedReferenceMetadata::Contigs &xmlContigs,
    ShouldLoadF shouldLoad,
    common::ThreadVector &loadThreads)
{
    std::vector<reference::Contig> ret;
    ret.reserve(xmlContigs.size());
    BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &xmlContig, xmlContigs)
    {
        ISAAC_ASSERT_MSG(ret.size() == xmlContig.index_ , "Expected sequentially ordered starting with 0");

        ret.push_back(reference::Contig(xmlContig.index_, xmlContig.name_));
    }

    std::vector<reference::SortedReferenceMetadata::Contig>::const_iterator nextContigToLoad = xmlContigs.begin();
    boost::mutex mutex;
    loadThreads.execute(boost::bind(&loadContigsParallel<ShouldLoadF>,
                                    boost::ref(shouldLoad),
                                    boost::ref(nextContigToLoad),
                                    xmlContigs.end(),
                                    boost::ref(ret),
                                    boost::ref(mutex)));

    return ret;
}

/**
 * \brief loads the fasta file contigs into memory on multiple threads
 */
template <typename FilterT> std::vector<std::vector<reference::Contig> > loadContigs(
    const reference::SortedReferenceMetadataList &SortedReferenceMetadataList,
    const FilterT &loadedContigFilter,
    common::ThreadVector &loadThreads)
{
    ISAAC_TRACE_STAT("loadContigs ");
    std::vector<std::vector<reference::Contig> > ret(SortedReferenceMetadataList.size());

    BOOST_FOREACH(const reference::SortedReferenceMetadata &SortedReferenceMetadata, SortedReferenceMetadataList)
    {
        const unsigned referenceIndex = &SortedReferenceMetadata - &SortedReferenceMetadataList.front();
        std::vector<reference::Contig> contigList =
            loadContigs(SortedReferenceMetadata.getContigs(),
                        boost::bind(&FilterT::isMapped, loadedContigFilter, referenceIndex, _1),
                        loadThreads);
        ret.at(referenceIndex).swap(contigList);
    }

    ISAAC_TRACE_STAT("loadContigs done ");

    return ret;
}


} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_CONTIGS_PRINTER_HH
