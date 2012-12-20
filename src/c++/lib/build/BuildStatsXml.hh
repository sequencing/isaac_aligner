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
 ** \file BuildStatsXml.hh
 **
 ** \brief Xml Serialization of Build statistics.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_BUILD_BUILD_STATS_XML_H
#define ISAAC_BUILD_BUILD_STATS_XML_H

#include "io/PtreeXml.hh"

namespace isaac
{
namespace build
{

class BuildStatsXml : public boost::property_tree::ptree
{
    static bool orderByProjectSample(const flowcell::BarcodeMetadata &left, const flowcell::BarcodeMetadata &right)
    {
        return
            left.getProject() < right.getProject() ||
                (left.getProject() == right.getProject() && (left.getSampleName() < right.getSampleName()));
    }

    void addProjectSampleReferenceContigBin(
        const std::string &projectName,
        const std::string &sampleName,
        const std::vector<reference::SortedReferenceXml::Contig> &contigs,
        const alignment::BinMetadata &bin,
        unsigned long totalFragments,
        unsigned long uniqueFragments)
    {
        if (contigs.size() > bin.getBinStart().getContigId())
        {
            // Some references have this bin and some are too short for that.
            const boost::property_tree::path projectSampleValuePrefix(
                "Stats/<indexed>Project/<name>" + projectName +
                "/<indexed>Sample/<name>" + sampleName, '/');

            const reference::SortedReferenceXml::Contig &contig = contigs.at(bin.getBinStart().getContigId());
            const boost::property_tree::path projectSampleContigValuePrefix(
                projectSampleValuePrefix / "<indexed>Contig" / boost::property_tree::path("<name>" + contig.name_, '/'));

            put(projectSampleContigValuePrefix / "ReferenceTotalBases", contig.totalBases_);


            const boost::property_tree::path projectSampleContigBinValuePrefix(
                projectSampleContigValuePrefix /
                "<indexed>Bin" / ("<offset>" + boost::lexical_cast<std::string>(bin.getBinStart().getPosition())).c_str());

            add(projectSampleContigBinValuePrefix / "TotalFragments", totalFragments);
            add(projectSampleContigBinValuePrefix / "UniqueFragments", uniqueFragments);
        }
    }

public:
    BuildStatsXml(
        const reference::SortedReferenceXmlList &sortedReferenceXmlList,
        const alignment::BinMetadataList &bins,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const BuildStats &buildStats)
    {
        flowcell::BarcodeMetadataList orderedBarcodeMetadataList(barcodeMetadataList);
        std::sort(orderedBarcodeMetadataList.begin(), orderedBarcodeMetadataList.end(), orderByProjectSample);

        BOOST_FOREACH(const alignment::BinMetadata &bin, bins)
        {
            if (bin.isUnalignedBin())
            {
                // skip unaligned reads stats. TODO: figure out if those can make useful part of build stats
                continue;
            }
            std::vector<reference::SortedReferenceXml::Contig> contigs;
            std::string lastProject, lastSample;
            unsigned long totalFragments = 0UL;
            unsigned long uniqueFragments = 0UL;
            BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, orderedBarcodeMetadataList)
            {
                if (!barcode.isUnmappedReference())
                {
                    if (lastProject != barcode.getProject() || lastSample != barcode.getSampleName())
                    {
                        if (!lastProject.empty())
                        {
                            addProjectSampleReferenceContigBin(lastProject, lastSample,
                                                               contigs, bin,
                                                               totalFragments, uniqueFragments);
                        }
                        lastProject = barcode.getProject();
                        lastSample = barcode.getSampleName();
                        contigs = sortedReferenceXmlList.at(barcode.getReferenceIndex()).getKaryotypeOrderedContigs();

                        totalFragments = 0UL;
                        uniqueFragments = 0UL;
                    }

                    totalFragments += buildStats.getTotalFragments(bin.getIndex(), barcode.getIndex());
                    uniqueFragments += buildStats.getUniqueFragments(bin.getIndex(), barcode.getIndex());
                }
            }
            if (!lastProject.empty())
            {
                addProjectSampleReferenceContigBin(lastProject, lastSample,
                                                   contigs, bin,
                                                   totalFragments, uniqueFragments);
            }
        }
    }
};

inline std::ostream &operator << (std::ostream &os, const BuildStatsXml &tree)
{
    return isaac::io::serializeAsXml(os, tree);
}

} //namespace build
} //namespace isaac

#endif //ISAAC_BUILD_BUILD_STATS_XML_H
