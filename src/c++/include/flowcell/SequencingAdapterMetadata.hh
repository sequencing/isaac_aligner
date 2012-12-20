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
 ** \file SequencingAdapterMetadata.hh
 **
 ** Packaging of the metadata associated to a barcode.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_SEQUENCING_ADAPTER_METADATA_HH
#define iSAAC_FLOWCELL_SEQUENCING_ADAPTER_METADATA_HH

#include <iostream>
#include <iterator>
#include <vector>

namespace isaac
{
namespace flowcell
{

/**
 ** \brief Description of a sequencing adapter.
 **/
class SequencingAdapterMetadata
{
public:
    /**
     * \brief default constructor creates 'unknown index' barcode
     */
    SequencingAdapterMetadata() : reverse_(false), clipLength_(0){;}

    SequencingAdapterMetadata(const std::string &sequence,
                              const bool reverse) :
                                  sequence_(sequence), reverse_(reverse), clipLength_(sequence_.size()){}

    SequencingAdapterMetadata(const std::string &sequence,
                              const bool reverse,
                              const unsigned clipLength) :
                                  sequence_(sequence), reverse_(reverse), clipLength_(clipLength){}

    const std::string &getSequence() const {return sequence_;}
    bool isReverse() const { return reverse_; }
    unsigned getClipLength() const {return clipLength_;}
    bool isUnbounded() const {return !clipLength_;}

    bool operator == (const SequencingAdapterMetadata &right) const
    {
        return sequence_ == right.sequence_ && reverse_ == right.reverse_ && clipLength_ == right.clipLength_;
    }
    bool operator != (const SequencingAdapterMetadata &right) const
    {
        return !(*this == right);
    }
private:
    // adapter sequence in the direction of the reference
    std::string sequence_;
    // Indicates the direction in which the adapter is expected to be sequenced in respect to reference
    bool reverse_;
    // number of bases to be clipped starting from the beginning of the adapter. 0 means all bases.
    int clipLength_;
};

typedef std::vector<SequencingAdapterMetadata> SequencingAdapterMetadataList;
extern const SequencingAdapterMetadataList NEXTERA_STANDARD_ADAPTERS;
extern const SequencingAdapterMetadataList NEXTERA_MATEPAIR_ADAPTERS;

inline std::ostream &operator<<(std::ostream &os, const SequencingAdapterMetadata &adapter)
{
    return os << "SequencingAdapterMetadata(" <<
        adapter.getSequence() << "," << adapter.isReverse() << "," << adapter.getClipLength() << ")";
}

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_SEQUENCING_ADAPTER_METADATA_HH

