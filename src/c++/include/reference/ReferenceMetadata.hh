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
 ** \file ReferenceMetadata.hh
 **
 ** Packaging of the metadata associated to a sorted reference.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_REFERENCE_METADATA_HH
#define iSAAC_FLOWCELL_REFERENCE_METADATA_HH

#include <iostream>
#include <vector>

namespace isaac
{
namespace reference
{

/**
 ** \brief Read-only interface to the metadata associated to a reference.
 **
 ** The intended usage is for barcode management in ordered collections (the index
 ** in the collectionis is associated to each barcode metadata instance).
 ** Index 0 is reserved for mapping the barcode sequences that don't match the known ones
 **
 **/
class ReferenceMetadata
{
public:
    ReferenceMetadata(const std::string &name,
                      const boost::filesystem::path &xmlPath,
                      const unsigned index) :
        name_(name), xmlPath_(xmlPath), index_(index)
    {}

    const std::string &getName() const {return name_;}
    const boost::filesystem::path &getXmlPath() const {return xmlPath_;}
    unsigned getIndex() const {return index_;}

private:
    std::string name_;
    boost::filesystem::path xmlPath_;
    unsigned index_;
};

typedef std::vector<ReferenceMetadata> ReferenceMetadataList;

inline std::ostream &operator<<(std::ostream &os, const ReferenceMetadata &referenceMetadata)
{
    return os << "ReferenceMetadata(" << referenceMetadata.getName() << "," <<
        referenceMetadata.getXmlPath() << "," << referenceMetadata.getIndex() << ")";
}

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_REFERENCE_METADATA_HH

