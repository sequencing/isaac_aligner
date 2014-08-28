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
 ** \file BuildStatsXml.cpp
 **
 ** \brief Xml Serialization of Build statistics.
 **
 ** \author Roman Petrovski
 **/

#include <functional>
#include <boost/foreach.hpp>

#include "BuildStatsXml.hh"
#include "xml/XmlWriter.hh"

namespace isaac
{
namespace build
{

inline bool orderByProjectSample(const flowcell::BarcodeMetadata &left, const flowcell::BarcodeMetadata &right)
{
    return
        left.getProject() < right.getProject() ||
        (left.getProject() == right.getProject() && (left.getSampleName() < right.getSampleName()));
}

struct GetBinContigId : public std::unary_function<const alignment::BinMetadataCRef, unsigned>
{
    result_type operator()(argument_type& binRef) const
    {
        const alignment::BinMetadata &bin = binRef;
        return bin.getBinStart().getContigId();
    }
};

BuildStatsXml::BuildStatsXml(
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
    const alignment::BinMetadataCRefList &bins,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const BuildStats &buildStats) :
    sortedReferenceMetadataList_(sortedReferenceMetadataList),
    bins_(bins),
    orderedBarcodeMetadataList_(barcodeMetadataList),
    buildStats_(buildStats)
{
    std::sort(orderedBarcodeMetadataList_.begin(), orderedBarcodeMetadataList_.end(), orderByProjectSample);
}

void BuildStatsXml::dumpContigs(
    xml::XmlWriter &xmlWriter,
    const std::vector<reference::SortedReferenceMetadata::Contig> &contigs,
    flowcell::BarcodeMetadataList::const_iterator sampleBarcodesBegin,
    flowcell::BarcodeMetadataList::const_iterator sampleBarcodesEnd)
{
    // skip unaligned as they cannot be asked for contig id and mess up the sort order
    alignment::BinMetadataCRefList::const_iterator binsBegin = std::find_if(
        bins_.begin(), bins_.end(),!boost::bind(&alignment::BinMetadata::isUnalignedBin, _1));
    alignment::BinMetadataCRefList::const_iterator binsEnd = std::find_if(
        std::reverse_iterator<alignment::BinMetadataCRefList::const_iterator>(bins_.end()),
        std::reverse_iterator<alignment::BinMetadataCRefList::const_iterator>(binsBegin),
        !boost::bind(&alignment::BinMetadata::isUnalignedBin, _1)).base();
    BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &contig, contigs)
    {
        typedef boost::transform_iterator<GetBinContigId, alignment::BinMetadataCRefList::const_iterator> ContigIdIterator;
        std::pair<ContigIdIterator, ContigIdIterator> binRange =
            std::equal_range(ContigIdIterator(binsBegin), ContigIdIterator(binsEnd), contig.index_);
        if (binRange.first != binRange.second)
        {
            ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Contig")
            {
                xmlWriter.writeAttribute("name", contig.name_);
                xmlWriter.writeElement("ReferenceTotalBases", contig.totalBases_);
                for(alignment::BinMetadataCRefList::const_iterator binIt = binRange.first.base();
                    binRange.second.base() != binIt; ++binIt)
                {
                    const alignment::BinMetadata &bin = *binIt;
                    if (bin.isUnalignedBin())
                    {
                        // skip unaligned reads stats. TODO: figure out if those can make useful part of build stats
                        continue;
                    }
                    if (bin.getBinStart().getContigId() != contig.index_)
                    {
                        continue;
                    }
                    const std::size_t binStatsIndex = std::distance(bins_.begin(), binIt);
                    const unsigned long totalFragments = std::accumulate(
                        sampleBarcodesBegin, sampleBarcodesEnd, 0UL,
                        bind(std::plus<unsigned long>(), _1,
                             boost::bind(&BuildStats::getTotalFragments, &buildStats_,
                                         binStatsIndex, boost::bind(&flowcell::BarcodeMetadata::getIndex, _2))));
                    const unsigned long uniqueFragments = std::accumulate(
                        sampleBarcodesBegin, sampleBarcodesEnd, 0UL,
                        bind(std::plus<unsigned long>(), _1,
                             boost::bind(&BuildStats::getUniqueFragments, &buildStats_,
                                         binStatsIndex, boost::bind(&flowcell::BarcodeMetadata::getIndex, _2))));

                    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Bin")
                    {
                        xmlWriter.writeAttribute("offset", bin.getBinStart().getPosition());
                        xmlWriter.writeElement("TotalFragments", totalFragments);
                        xmlWriter.writeElement("UniqueFragments", uniqueFragments);
                    }
                }
            }
        }
    }
}

void BuildStatsXml::serialize(std::ostream &os)
{
    ISAAC_THREAD_CERR << "Generating Build statistics" << std::endl;
    xml::XmlWriter xmlWriter(os);
    ISAAC_XML_WRITER_ELEMENT_BLOCK(xmlWriter, "Stats")
    {
        // multiple bins can belong to the same contig. We need to avoid creating duplicate Contig elements
        std::string lastProject, lastSample;
        flowcell::BarcodeMetadataList::const_iterator lastBarcodeIt = orderedBarcodeMetadataList_.begin();
        for(flowcell::BarcodeMetadataList::const_iterator barcodeIt = orderedBarcodeMetadataList_.begin();
            orderedBarcodeMetadataList_.end() != barcodeIt; ++barcodeIt)
        {
            const flowcell::BarcodeMetadata &barcode = *barcodeIt;
            if (!barcode.isUnmappedReference())
            {
                if (lastSample != barcode.getSampleName() || lastProject != barcode.getProject())
                {
                    if (!lastProject.empty())
                    {
                        dumpContigs(
                            xmlWriter, sortedReferenceMetadataList_.at(lastBarcodeIt->getReferenceIndex()).getKaryotypeOrderedContigs(),
                            lastBarcodeIt, barcodeIt);

                        if (lastProject != barcode.getProject())
                        {
                            xmlWriter.endElement(); //close Sample
                            xmlWriter.endElement(); //close Project

                            xmlWriter.startElement("Project");
                            xmlWriter.writeAttribute("name", barcode.getProject());
                        }
                        else
                        {
                            ISAAC_ASSERT_MSG(lastSample != barcode.getSampleName(), "Sample name was expected to be different here");
                            xmlWriter.endElement(); //close Sample
                        }

                        xmlWriter.startElement("Sample");
                        xmlWriter.writeAttribute("name", barcode.getSampleName());
                    }
                    else
                    {
                        xmlWriter.startElement("Project");
                        xmlWriter.writeAttribute("name", barcode.getProject());
                        xmlWriter.startElement("Sample");
                        xmlWriter.writeAttribute("name", barcode.getSampleName());
                    }
                    lastProject = barcode.getProject();
                    lastSample = barcode.getSampleName();
                    lastBarcodeIt = barcodeIt;
                }
            }
        }

        if (!lastProject.empty())
        {
            dumpContigs(
                xmlWriter, sortedReferenceMetadataList_.at(lastBarcodeIt->getReferenceIndex()).getKaryotypeOrderedContigs(),
                lastBarcodeIt, orderedBarcodeMetadataList_.end());
            xmlWriter.endElement(); //close Sample
            xmlWriter.endElement(); //close Project
        }
    }
    ISAAC_THREAD_CERR << "Generating Build statistics done" << std::endl;
}

} //namespace build
} //namespace isaac

