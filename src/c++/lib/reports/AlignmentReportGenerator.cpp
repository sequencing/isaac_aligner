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
 ** \file AlignmentReportGenerator.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 **
 ** \author Roman Petrovski
 **/

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

//#include <libxml/xmlmemory.h>
//#include <libxml/debugXML.h>
//#include <libxml/HTMLtree.h>
//#include <libxml/xmlIO.h>
//#include <libxml/DOCBparser.h>
//#include <libxml/xinclude.h>
//#include <libxml/catalog.h>
#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>
#include <libexslt/exslt.h>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FileSystem.hh"
#include "common/Process.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "package/InstallationPaths.hh"
#include "reports/AlignmentReportGenerator.hh"

extern int xmlLoadExtDtdDefaultValue;

namespace isaac
{
namespace reports
{

AlignmentReportGenerator::AlignmentReportGenerator(
    const std::vector<flowcell::Layout> &flowcellLayoutList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const boost::filesystem::path &alignmentStatsXmlPath,
    const boost::filesystem::path &demultiplexingStatsXmlPath,
    const boost::filesystem::path &tempDirectory,
    const boost::filesystem::path &outputDirectory,
    const ImageFileFormat imageFileFormat)
    :flowcellLayoutList_(flowcellLayoutList),
     alignmentStatsXmlPath_(alignmentStatsXmlPath),
     demultiplexingStatsXmlPath_(demultiplexingStatsXmlPath),
     tempDirectory_(tempDirectory),
     gnuplotScriptPath_(tempDirectory_/"gnuplot-generate-images.dat"),
     outputDirectoryHtml_(outputDirectory/"html"),
     outputDirectoryImages_(outputDirectory/"gif"),
     imageFileFormat_(imageFileFormat)

{
    std::vector<boost::filesystem::path> createList =
        boost::assign::list_of(outputDirectoryHtml_)(outputDirectoryImages_);
    createList.push_back(outputDirectoryHtml_/"all");
    createList.push_back(outputDirectoryHtml_/"all"/"all");
    createList.push_back(outputDirectoryHtml_/"all"/"all"/"all");
    createList.push_back(outputDirectoryHtml_/"all"/"all"/"all"/"all");
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcodeMetadata, barcodeMetadataList)
    {
        createList.push_back(outputDirectoryHtml_/"all"/barcodeMetadata.getProject());
        createList.push_back(outputDirectoryHtml_/"all"/barcodeMetadata.getProject()/"all");
        createList.push_back(outputDirectoryHtml_/"all"/barcodeMetadata.getProject()/"all"/"all");
        createList.push_back(outputDirectoryHtml_/"all"/barcodeMetadata.getProject()/barcodeMetadata.getSampleName());
        createList.push_back(outputDirectoryHtml_/"all"/barcodeMetadata.getProject()/barcodeMetadata.getSampleName()/"all");
    }
    BOOST_FOREACH(const flowcell::Layout& flowcell, flowcellLayoutList_)
    {
        const boost::filesystem::path flowcellHtmlPath = outputDirectoryHtml_ / flowcell.getFlowcellId();
        createList.push_back(flowcellHtmlPath);
        const boost::filesystem::path flowcellImagesPath = outputDirectoryImages_ / flowcell.getFlowcellId();
        createList.push_back(flowcellImagesPath);
        createList.push_back(flowcellHtmlPath/"all");
        createList.push_back(flowcellHtmlPath/"all"/"all");
        createList.push_back(flowcellHtmlPath/"all"/"all"/"all");
        createList.push_back(flowcellImagesPath/"all");
        createList.push_back(flowcellImagesPath/"all"/"all");
        createList.push_back(flowcellImagesPath/"all"/"all"/"all");
        BOOST_FOREACH(const flowcell::BarcodeMetadata &barcodeMetadata, barcodeMetadataList)
        {
            createList.push_back(flowcellHtmlPath/barcodeMetadata.getProject());
            createList.push_back(flowcellHtmlPath/barcodeMetadata.getProject()/"all");
            createList.push_back(flowcellHtmlPath/barcodeMetadata.getProject()/"all"/"all");
            createList.push_back(flowcellHtmlPath/barcodeMetadata.getProject()/barcodeMetadata.getSampleName());
            createList.push_back(flowcellHtmlPath/barcodeMetadata.getProject()/barcodeMetadata.getSampleName()/"all");
            createList.push_back(flowcellHtmlPath/barcodeMetadata.getProject()/barcodeMetadata.getSampleName()/barcodeMetadata.getName());

//            createList.push_back(flowcellImagesPath/barcodeMetadata.getProject());
//            createList.push_back(flowcellImagesPath/barcodeMetadata.getProject()/barcodeMetadata.getSampleName());
//            createList.push_back(flowcellImagesPath/barcodeMetadata.getProject()/barcodeMetadata.getSampleName()/barcodeMetadata.getName());
        }
    }

    common::createDirectories(createList);

    boost::filesystem::remove(gnuplotScriptPath_);
    if (boost::filesystem::exists(gnuplotScriptPath_))
    {
        BOOST_THROW_EXCEPTION(
            common::IoException(errno, "Failed to remove the following temporary file:" + gnuplotScriptPath_.string()));
    }

}

void AlignmentReportGenerator::run()
{
    const char *params[] = {"OUTPUT_DIRECTORY_HTML_PARAM", "''",
                            "OUTPUT_DIRECTORY_IMAGES_PARAM", "''",
                            "TEMP_DIRECTORY_PARAM", "''",
                            "DEMULTIPLEXING_STATS_XML_PARAM", "''",
                            "GNUPLOT_SCRIPT_PATH_PARAM", "''",
                            0};
    std::string quotedOutputHtmlDirectory = "'" + outputDirectoryHtml_.string() + "'";
    params[1] = quotedOutputHtmlDirectory.c_str();
    std::string quotedOutputImagesDirectory = "'" + outputDirectoryImages_.string() + "'";
    params[3] = quotedOutputImagesDirectory.c_str();
    std::string quotedTempDirectory = "'" + tempDirectory_.string() + "'";
    params[5] = quotedTempDirectory.c_str();
    std::string quotedDemultiplexingStatsXml= "'" + demultiplexingStatsXmlPath_.string() + "'";
    params[7] = quotedDemultiplexingStatsXml.c_str();
    std::string quotedGnuplotScriptPath = "'" + gnuplotScriptPath_.string() + "'";
    params[9] = quotedGnuplotScriptPath.c_str();

    xmlSubstituteEntitiesDefault(1);
    xmlLoadExtDtdDefaultValue = 1;

    exsltRegisterAll();

    xsltStylesheetPtr cur = xsltParseStylesheetFile((const xmlChar *)iSAAC_FULL_DATADIR_XSL_ALIGNMENT_NON_MULTIPLEXED_REPORT_XSL.c_str());
    xmlDocPtr doc = xmlParseFile(alignmentStatsXmlPath_.c_str());
    xmlDocPtr res = xsltApplyStylesheet(cur, doc, params);
//    int xsltprocRes = xsltSaveResultToFile(0, res, cur);

    xsltFreeStylesheet(cur);
    xmlFreeDoc(res);
    xmlFreeDoc(doc);

    xsltCleanupGlobals();
    xmlCleanupParser();
    if (!res)
    {
        BOOST_THROW_EXCEPTION(common::LibXsltException());
    }
//    if (-1 == xsltprocRes)
//    {
//        BOOST_THROW_EXCEPTION(common::LibXsltException());
//    }

    // Only generate output plots if not turned off
    if(none != imageFileFormat_)
    {
        // And the gnuplot script path exists
        if (boost::filesystem::exists(gnuplotScriptPath_))
        {
            std::string gnuplotCmd = "gnuplot " + gnuplotScriptPath_.string();
            common::executeCommand(gnuplotCmd);
        }
    }
}


} // namespace reports
} // namespace isaac
