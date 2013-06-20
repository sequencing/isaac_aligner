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
 ** \file SeedMetadata.hh
 **
 ** \brief Metadata associated to a seed (offset, length, read, index)
 ** 
 ** A seed is defined as a contiguous sequence of cycles, starting at a specific
 ** offset on a given read.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_SEED_METADATA_HH
#define iSAAC_ALIGNMENT_SEED_METADATA_HH

#include <utility>
#include <iostream>

#include "common/Debug.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Trivial representation of a seed as a contiguous set of cycles.
 **
 ** The intended usage is for seed management in ordered collections (the index
 ** in the collection is associated to each instance).
 **/
class SeedMetadata
{
public:
    /**
     ** \brief Constructor for the metadata relevant to a seed
     **
     ** \param offset offset of the seed in the read
     ** \param length length of the seed
     ** \param read index of the read in the associated list of ReadMetadata
     ** \param index index of this instance in the associated list of SeedMetadata
     **/
    SeedMetadata(unsigned int offset, unsigned length, unsigned int readIndex, unsigned int index)
        : offset_(offset)
        , length_(length)
        , readIndex_(readIndex)
        , index_(index)
    {
        ISAAC_ASSERT_MSG(offset_ < std::numeric_limits<unsigned short>::max(), "Unexpectedly large seed offset");
        ISAAC_ASSERT_MSG(16 == length || 32 == length_ || 64 == length_, "Unexpected seed length. Only seed length 32 and 64 are supported");
    }
    SeedMetadata(const SeedMetadata &seed)
        : offset_(seed.offset_)
        , length_(seed.length_)
        , readIndex_(seed.readIndex_)
        , index_(seed.index_)
    {}
    SeedMetadata &operator=(const SeedMetadata &seed)
    {
        if (this != &seed)
        {
            offset_ = seed.offset_;
            length_ = seed.length_;
            readIndex_ = seed.readIndex_;
            index_ = seed.index_;
        }
        return *this;
    }
    virtual ~SeedMetadata() {}
    unsigned short getOffset() const {return offset_;}
    unsigned short getLength() const {return length_;}
    unsigned int getReadIndex() const {return readIndex_;}
    unsigned int getIndex() const {return index_;}
protected:
    void setOffset(unsigned int offset) {offset_ = offset;}
    void setLength(unsigned int length) {length_ = length;}
    void setReadIndex(unsigned int readIndex) {readIndex_ = readIndex;}
    void setIndex(unsigned int index) {index_ = index;}
private:
    unsigned short offset_;
    unsigned short length_;
    unsigned int readIndex_;
    unsigned int index_;
};

typedef std::vector<SeedMetadata> SeedMetadataList;

inline std::ostream &operator<<(std::ostream &os, const SeedMetadata &seedMetadata)
{
    return os << "SeedMetadata("
              << seedMetadata.getOffset() << ", "
              << seedMetadata.getLength() << ", "
              << seedMetadata.getReadIndex() << ", "
              << seedMetadata.getIndex()
              << ")";
}

inline std::ostream &operator<<(std::ostream &os, const std::vector<SeedMetadata> &seedMetadataList)
{
    std::copy(seedMetadataList.begin(), seedMetadataList.end(), std::ostream_iterator<SeedMetadata>(os, " "));
    return os;
}


} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEED_METADATA_HH
