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
 ** \file AlignmentReportGenerator.hh
 **
 ** Reorders alingments and stores them in bam file.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REPORTS_ALIGNMENT_REPORT_GENERATOR_HH
#define iSAAC_REPORTS_ALIGNMENT_REPORT_GENERATOR_HH

#include <boost/filesystem.hpp>

#include "flowcell/Layout.hh"

namespace isaac
{
namespace reports
{

class AlignmentReportGenerator
{
public:
    /**
     * \brief Stats image file format choices
     */
    enum ImageFileFormat
    {
        gif,       // Produce .gif files for plots
        none       // Do not produce any plot files
    };


    const std::vector<flowcell::Layout> &flowcellLayoutList_;
    const boost::filesystem::path alignmentStatsXmlPath_;
    const boost::filesystem::path demultiplexingStatsXmlPath_;
    const boost::filesystem::path tempDirectory_;
    const boost::filesystem::path gnuplotScriptPath_;
    const boost::filesystem::path outputDirectoryHtml_;
    const boost::filesystem::path outputDirectoryImages_;
    const ImageFileFormat imageFileFormat_;

public:
    AlignmentReportGenerator(
        const std::vector<flowcell::Layout> &flowcellLayoutList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const boost::filesystem::path &alignmentStatsXmlPath,
        const boost::filesystem::path &demultiplexingStatsXmlPath,
        const boost::filesystem::path &tempDirectory,
        const boost::filesystem::path &outputDirectory,
        const ImageFileFormat imageFileFormat);

    void run();

};

} // namespace reports
} // namespace isaac

#endif // #ifndef iSAAC_REPORTS_ALIGNMENT_REPORT_GENERATOR_HH
