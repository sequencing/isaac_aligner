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
 ** \file readMetadata.hh
 **
 ** Packaging of the metadata associated to a read.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_FLOWCELL_READ_METADATA_HH
#define iSAAC_FLOWCELL_READ_METADATA_HH

#include <iostream>
#include <vector>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

namespace isaac
{
namespace flowcell
{

/**
 ** \brief Read-only interface to the metadata associated to a read.
 **
 ** The intended usage is for read management in ordered collections (the index
 ** in the collectionis is associated to each read metadata instance).
 **
 ** \todo Provide a robust and flexible indexing mechanism for the reads.
 **/
class ReadMetadata
{
public:
    /**
     ** \brief Detailed constructor
     **
     ** \param number user-facing read number. Will appear in report and file names
     **
     ** \param cycleList list of cycles for that read
     **
     ** \param index 0-based index of the read - used to reverse-map to the read list
     **
     ** \param firstReadCycle first cycle that belongs to the read without use base masking applied
     **
     ** \param offset -1U for index reads, partial sum of the lengths of the previous reads for data reads
     **/
    ReadMetadata(const unsigned number, const std::vector<unsigned> &cycleList, unsigned index, unsigned offset, unsigned firstReadCycle);
    // constructor for unit tests, don't use elsewhere!
    ReadMetadata(unsigned firstCycle, unsigned lastCycle, unsigned index, unsigned offset);
    /// \brief Allow inheritance
    virtual ~ReadMetadata() {}
    unsigned getLength() const {return cycleList_.size();}
    unsigned getFirstReadCycle() const {return firstReadCycle_;}
    unsigned getFirstCycle() const {return cycleList_.front();}
    unsigned getLastCycle() const {return cycleList_.back();}
    const std::vector<unsigned> &getCycles() const {return cycleList_;}
    unsigned getIndex() const {return index_;}
    unsigned getNumber() const {return number_;}
    unsigned getOffset() const {return offset_;}
    bool operator==(const ReadMetadata &rhs) const;
    bool operator!=(const ReadMetadata &rhs) const
    {
        return !(rhs == *this);
    }

    void setNumber(unsigned number){number_ = number;}
private:
    unsigned number_;
    std::vector<unsigned> cycleList_;
    unsigned index_;
    unsigned offset_;
    unsigned firstReadCycle_;
};

class ReadMetadataList : public std::vector<ReadMetadata>
{
public:
    ReadMetadataList(const std::vector<ReadMetadata> &that):
        std::vector<ReadMetadata>(that){}
    ReadMetadataList(){}
};


unsigned getTotalReadLength(const ReadMetadataList &readMetadataList);
unsigned getMaxCycleNumber(const ReadMetadataList &readMetadataList);

std::vector<unsigned> getAllCycleNumbers(const ReadMetadataList &readMetadataList);

inline std::ostream &operator<<(std::ostream &os, const ReadMetadata &readMetadata)
{
    return os << "ReadMetadata(" 
              << readMetadata.getNumber() << ", "
              << readMetadata.getLength() << " [" 
              << readMetadata.getFirstCycle() << ", " 
              << readMetadata.getLastCycle() << "], " 
              << readMetadata.getIndex() << "id, "
              << readMetadata.getOffset() << "off,"
              << readMetadata.getFirstReadCycle() << "frc)";
}

inline unsigned getMaxReadLength(const flowcell::ReadMetadataList &readMetadataList)
{
    return std::max_element(readMetadataList.begin(), readMetadataList.end(),
                            boost::bind(&flowcell::ReadMetadata::getLength, _1)<
                            boost::bind(&flowcell::ReadMetadata::getLength, _2))->getLength();
}


} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_READ_METADATA_HH

