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
 ** \file Layout.cpp
 **
 ** Flowcell file locations and such
 **
 ** \author Roman Petrovski
 **/

#include "rta/ConfigXml.hh"
#include "flowcell/Layout.hh"

namespace isaac
{
namespace flowcell
{


Layout::Layout(const boost::filesystem::path &baseCallsDirectory,
        const Format format,
        const FormatSpecificData &formatSpecificData,
        const unsigned laneNumberMax,
        const std::vector<unsigned> &barcodeCycles,
        const flowcell::ReadMetadataList &readMetadataList,
        const alignment::SeedMetadataList &seedMetadataList,
        const std::string &flowcellId)
     : baseCallsPath_(baseCallsDirectory)
     , format_(format)
     , formatSpecificData_(formatSpecificData)
     , laneNumberMax_(laneNumberMax)
     , barcodeCycles_(barcodeCycles)
     , flowcellId_(flowcellId)
     , readMetadataList_(readMetadataList)
     , seedMetadataList_(seedMetadataList)
     , dataCycles_(flowcell::getAllCycleNumbers(readMetadataList_))
     , index_(0)
{
}


} // namespace flowcell
} // namespace isaac
