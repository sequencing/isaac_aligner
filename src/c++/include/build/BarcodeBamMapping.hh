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
 ** \file BarcodeBamMapping.hh
 **
 ** Helper class for mapping barcodes to output files.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BARCODE_BAM_HH
#define iSAAC_BUILD_BARCODE_BAM_HH

#include <boost/filesystem.hpp>

#include "flowcell/BarcodeMetadata.hh"

namespace isaac
{
namespace build
{

class BarcodeBamMapping
{
public:
    BarcodeBamMapping(){}
    BarcodeBamMapping(
        const std::vector<unsigned> &sampleIds,
        const std::vector<boost::filesystem::path> &samplePaths):
            first(sampleIds), second(samplePaths){}
    /// Each position in the vector contains unique index of the project-sample
    typedef std::vector<unsigned> BarcodeSampleIndexMap;
    const BarcodeSampleIndexMap &getIndexMap() const {return first;}
    const std::vector<boost::filesystem::path> &getPaths() const {return second;}
    unsigned getTotalBarcodes() const {return first.size();}
    unsigned getTotalFiles() const {return second.size();}
    unsigned getFileIndex(const unsigned barcodeIndex) const {return first.at(barcodeIndex);}
    unsigned getFileIndex(const flowcell::BarcodeMetadata &barcode) const {return getFileIndex(barcode.getIndex());}
    const boost::filesystem::path &getFilePath(const flowcell::BarcodeMetadata &barcode) const
    {
        return second.at(getFileIndex(barcode));
    }
    void mapToNew(const flowcell::BarcodeMetadata &barcode, const boost::filesystem::path bamPath)
    {
        ISAAC_ASSERT_MSG(barcode.getIndex() == first.size(), "this implementation expects barcodes to arrive sequentially without gaps");
        first.push_back(second.size());
        second.push_back(bamPath);
    }
    void mapToExisting(const flowcell::BarcodeMetadata &barcode, const unsigned pathIndex)
    {
        ISAAC_ASSERT_MSG(pathIndex < second.size(), "pathIndex outside of range");
        ISAAC_ASSERT_MSG(barcode.getIndex() == first.size(), "this implementation expects barcodes to arrive sequentially without gaps");
        first.push_back(pathIndex);
    }
private:
    template<class Archive> friend void serialize(Archive & ar, BarcodeBamMapping &, const unsigned int file_version);
    std::vector<unsigned> first;
    std::vector<boost::filesystem::path> second;
};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BARCODE_BAM_HH
