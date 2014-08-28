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
 ** \file BarcodeBamMapping.hh
 **
 ** Helper class for mapping barcodes to output files.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BARCODE_BAM_HH
#define iSAAC_BUILD_BARCODE_BAM_HH

#include <iterator>
#include <boost/filesystem.hpp>

#include "flowcell/BarcodeMetadata.hh"

namespace isaac
{
namespace build
{

class BarcodeBamMapping
{
public:
    typedef std::vector<unsigned> BarcodeSampleIndexMap;
    typedef std::vector<unsigned> BarcodeProjectIndexMap;

    BarcodeBamMapping() : projectIndexMax_(-1U){}
    /**
     * \param projectIds one entry per barcode index mapping it to the corresponding project id
     * \param sampleIds one entry per barcode index mapping it to the corresponding sample
     * \param samplePaths one entry per sample id
     */
    BarcodeBamMapping(
        const BarcodeProjectIndexMap &projectIds,
        const BarcodeSampleIndexMap &sampleIds,
        const std::vector<boost::filesystem::path> &samplePaths):
            barcodeProjectIndex_(projectIds),
            projectIndexMax_(std::distance(barcodeProjectIndex_.begin(), std::max_element(barcodeProjectIndex_.begin(), barcodeProjectIndex_.end()))),
            barcodeSampleIndex_(sampleIds), samplePaths(samplePaths){}
    /// Each position in the vector contains unique index of the project-sample
    const BarcodeSampleIndexMap &getSampleIndexMap() const {return barcodeSampleIndex_;}
    const std::vector<boost::filesystem::path> &getPaths() const {return samplePaths;}
    unsigned getTotalBarcodes() const {return barcodeSampleIndex_.size();}
    unsigned getTotalSamples() const {return samplePaths.size();}
    unsigned getProjectIndex(const unsigned barcodeIndex) const {return barcodeProjectIndex_.at(barcodeIndex);}
    unsigned getMaxProjectIndex() const {return projectIndexMax_;}
    unsigned getSampleIndex(const unsigned barcodeIndex) const {return barcodeSampleIndex_.at(barcodeIndex);}
    const boost::filesystem::path &getFilePath(const flowcell::BarcodeMetadata &barcode) const
    {
        return samplePaths.at(getSampleIndex(barcode.getIndex()));
    }

private:
    template<class Archive> friend void serialize(Archive & ar, BarcodeBamMapping &, const unsigned int file_version);
    BarcodeProjectIndexMap barcodeProjectIndex_;
    unsigned projectIndexMax_;
    BarcodeSampleIndexMap barcodeSampleIndex_;
    std::vector<boost::filesystem::path> samplePaths;
};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BARCODE_BAM_HH
