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
 ** \file CasavaIntegration.cpp
 **
 ** \brief see CasavaIntegration.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include "casava/CasavaIntegration.hh"
#include "common/config.h"
#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/Process.hh"

namespace isaac
{
namespace casava
{

CasavaIntegration::CasavaIntegration(
    const boost::filesystem::path &tempDirectory,
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const build::BarcodeBamMapping &barcodeBamMapping,
    const reference::SortedReferenceXmlList &sortedReferenceXmlList,
    const flowcell::FlowcellLayoutList &flowcellLayoutList) :
    tempDirectoryCasava_(tempDirectory / "Casava"),
    glueMakefile_(tempDirectoryCasava_ / "Makefile"),
    outputDirectory_(outputDirectory),
    barcodeMetadataList_(barcodeMetadataList),
    barcodeBamMapping_(barcodeBamMapping),
    sortedReferenceXmlList_(sortedReferenceXmlList),
    flowcellLayoutList_(flowcellLayoutList)
{

}

void CasavaIntegration::reset(const std::string &casavaArgv)
{
    std::vector<boost::filesystem::path> createDirList;
    createDirList.push_back(tempDirectoryCasava_);

    ISAAC_THREAD_CERR << "Erasing " << tempDirectoryCasava_ << std::endl;
    boost::filesystem::remove_all(tempDirectoryCasava_);
    ISAAC_THREAD_CERR << "Erasing done " << tempDirectoryCasava_ << std::endl;

    unsigned casavaBuildIndexCreate = 0;

    std::vector<std::pair<boost::filesystem::path, boost::filesystem::path> > buildPathList;
    std::vector<boost::filesystem::path> resultsPathList;
    std::vector<std::string> configCmds;
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList_)
    {
        if (casavaBuildIndexCreate == barcodeBamMapping_.getFileIndex(barcode))
        {
            const boost::filesystem::path &bamPath = barcodeBamMapping_.getFilePath(barcode);
            if (!barcode.isUnmappedReference())
            {
                ISAAC_THREAD_CERR << "Configuring CASAVA for " << bamPath << std::endl;

                const boost::filesystem::path projectDir = tempDirectoryCasava_ / barcode.getProject();
                const boost::filesystem::path tmpDir = projectDir / barcode.getSampleName();
                const boost::filesystem::path sampleOutputDir =
                    outputDirectory_ / barcode.getProject() / barcode.getSampleName() / "Casava";

                configCmds.push_back(makeConfigureCasavaCmd(flowcellLayoutList_.at(barcode.getFlowcellIndex()),
                                                            sortedReferenceXmlList_.at(barcode.getReferenceIndex()),
                                                            bamPath, tmpDir, casavaArgv));
                createDirList.push_back(projectDir);
                createDirList.push_back(sampleOutputDir);
                buildPathList.push_back(std::make_pair(tmpDir, sampleOutputDir));

                ISAAC_THREAD_CERR << "Configuring CASAVA done for " << bamPath << std::endl;
            }
            ++casavaBuildIndexCreate;
        }
    }
    ISAAC_ASSERT_MSG(barcodeBamMapping_.getTotalFiles() == casavaBuildIndexCreate, "must process all bam files");


    common::createDirectories(createDirList);
    std::for_each(configCmds.begin(), configCmds.end(), &common::executeCommand);

    spitOutTheGlueMakeFile(buildPathList);
}

boost::filesystem::path getFastaFilePath(
    const reference::SortedReferenceXml &sortedReferenceXml)
{
    const reference::SortedReferenceXml::Contigs contigs = sortedReferenceXml.getContigs();
    reference::SortedReferenceXml::Contigs::const_iterator mismatchingFasta =
        std::find_if(contigs.begin(), contigs.end(),
                     boost::bind(&reference::SortedReferenceXml::Contig::filePath_, _1) != contigs.at(0).filePath_);
    if(contigs.end() != mismatchingFasta)
    {
        BOOST_THROW_EXCEPTION(common::FeatureNotAvailable(
            (boost::format("Casava integration requires single reference genome file. Found: %s and %s") %
                contigs.at(0).filePath_ % mismatchingFasta->filePath_).str()));
    }

    return contigs.at(0).filePath_;
}

void CasavaIntegration::execute(const unsigned coresMax)
{
    common::executeCommand(std::string("make CASAVA_LOG_LEVEL:=3 ") +
        " -j " + boost::lexical_cast<std::string>(coresMax) +
        " -C " + tempDirectoryCasava_.string());
}

void CasavaIntegration::spitOutTheGlueMakeFile(
    const std::vector<PairTempFinal > &buildPathList)
{
    std::ofstream os(glueMakefile_.c_str());
    if (!os)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to create file:" + glueMakefile_.string()));
    }

    static const std::string finalFile = (boost::filesystem::path("tmp") / "casava.finished").string();
    if (!(os << "all: transfer.finished" << std::endl << std::endl << "%/" << finalFile <<
        ":\n\t$(MAKE) -C $(dir $@).. && touch $@"  <<
        std::endl << std::endl))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to append to file:" + glueMakefile_.string()));
    }

    if (!(os << "transfer.finished: "))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write rule to file:" + glueMakefile_.string()));
    }

    BOOST_FOREACH(const PairTempFinal &tempFinal, buildPathList)
    {
        const boost::filesystem::path casavaFinishedFilePath = tempFinal.first / finalFile;
        if (!(os << " " << casavaFinishedFilePath.string()))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write dependency to file:" + glueMakefile_.string()));
        }
    }

    BOOST_FOREACH(const PairTempFinal &tempFinal, buildPathList)
    {
        if (!(os <<
            "\n\t" << iSAAC_FULL_LIBEXECDIR_CASAVA_EXTRACT_RESULTS_SH <<
                " '" << tempFinal.first.string() <<
                "' '" << tempFinal.second.string() + "' && \\"))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write transfer cmd to file:" + glueMakefile_.string()));
        }
    }

    if (!(os << "\n\ttouch $@" << std::endl << std::endl))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to finish file:" + glueMakefile_.string()));
    }
}

std::string CasavaIntegration::makeConfigureCasavaCmd(
    const flowcell::Layout &flowcell,
    const reference::SortedReferenceXml &sortedReferenceXml,
    const boost::filesystem::path &bamPath,
    const boost::filesystem::path &tmpDir,
    const std::string &casavaArgv) const
{
    const std::string configureBuildCmd(
        std::string(iSAAC_FULL_BINDIR_CONFIGURE_BUILD_SH) +
        ((1 == flowcell.getReadMetadataList().size()) ? " --readMode single " : " --readMode paired ") +
        std::string(" --indelsReadStartDepthCutoff 10") +
// CNVseg has hard-coded bin size of 50000000 and segfaults otherwise
//        " --binSizeProject 5000000"
//        " --binSizeBuild 5000000" +
        casavaArgv +
        " --samtoolsRefFile '" + getFastaFilePath(sortedReferenceXml).string() + "'"
        " --inSortedBam '" + bamPath.string() + "'" +
        " --outDir '" + tmpDir.string() + "'" +
        " --make ");

    return configureBuildCmd;
}

} // namespace casava
} // namespace isaac
