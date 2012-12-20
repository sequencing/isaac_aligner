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
 ** \file CasavaIntegration.hh
 **
 ** \brief Helper for configuring and executing CASAVA variant calling.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_CASAVA_CASAVA_INTEGRATION_H
#define ISAAC_CASAVA_CASAVA_INTEGRATION_H

#include <boost/filesystem.hpp>

#include "build/BarcodeBamMapping.hh"
#include "flowcell/Layout.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac {

namespace casava {


class CasavaIntegration: boost::noncopyable
{
    const boost::filesystem::path tempDirectoryCasava_;
    const boost::filesystem::path glueMakefile_;

    const boost::filesystem::path &outputDirectory_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const build::BarcodeBamMapping &barcodeBamMapping_;
    const reference::SortedReferenceXmlList &sortedReferenceXmlList_;
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;

public:
    CasavaIntegration(
        const boost::filesystem::path &tempDirectory,
        const boost::filesystem::path &outputDirectory,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const build::BarcodeBamMapping &barcodeBamMapping,
        const reference::SortedReferenceXmlList &sortedReferenceXmlList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList);

    void reset(const std::string &casavaArgv);
    void execute(const unsigned maxCores);

private:
    std::string makeConfigureCasavaCmd(
        const flowcell::Layout &flowcell,
        const reference::SortedReferenceXml &sortedReferenceXml,
        const boost::filesystem::path &bamPath,
        const boost::filesystem::path &tmpDir,
        const std::string &casavaArgv) const;

    typedef std::pair<boost::filesystem::path, boost::filesystem::path> PairTempFinal;
    void spitOutTheGlueMakeFile(
        const std::vector<PairTempFinal > &buildPathList);

};

} //namespace casava
} //namespace isaac

#endif //ISAAC_CASAVA_CASAVA_INTEGRATION_H
