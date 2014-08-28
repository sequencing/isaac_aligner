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
 ** \file FragmentCollector.hh
 **
 ** \brief Buffering of fragments with indexing information for BAM generation.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_COLLECTOR_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_COLLECTOR_HH

#include <boost/filesystem.hpp>
#include <boost/integer/static_min_max.hpp>

#include "alignment/FragmentMetadata.hh"
#include "alignment/BamTemplate.hh"
#include "alignment/BinMetadata.hh"
#include "alignment/matchSelector/BinIndexMap.hh"
#include "io/Fragment.hh"


namespace isaac
{
namespace alignment
{
namespace matchSelector
{

/**
 * \brief buffer of fixed-size records capable of holding all cluster fragments
 */
class FragmentBuffer
{
public:
    FragmentBuffer(const unsigned long reserveClusters,
                   const flowcell::FlowcellLayoutList &flowcellLayoutList):
        recordLength_(getRecordLength(flowcellLayoutList)),
        readOffsets_(getReadOffsets(flowcellLayoutList)), // offset of the first read is 0
        clusters_(reserveClusters)
    {
        // data buffer is pre-allocated as the fragments will be put in there by mutiple threads using cluster id
        // as target location
        reserve(reserveClusters);
        ISAAC_THREAD_CERR << "Constructed FragmentBuffer for " << clusters_ << " clusters. Record length: "
            << recordLength_ << " read offsets : " << readOffsets_[0] << "," <<
            (READS_MAX == readOffsets_.size() ? readOffsets_[1] : 0)  <<  std::endl;
    }

    class IndexRecord
    {
        // this constructor is only used for searching
        IndexRecord(reference::ReferencePosition fStrandPos, char *dataBytes)
            : fStrandPos_(fStrandPos), dataBytes_(dataBytes)
        {
        }
        friend class FragmentBuffer;

    public:
        IndexRecord() : fStrandPos_(0), dataBytes_(0)
        {
        }

        reference::ReferencePosition fStrandPos_;
        char *dataBytes_;

        bool initialized() const
        {
            return !!dataBytes_;
        }

        io::FragmentHeader &fragmentHeader()
        {
            return *reinterpret_cast<io::FragmentHeader*>(dataBytes_);
        }

        const io::FragmentHeader &fragmentHeader() const
        {
            return *reinterpret_cast<const io::FragmentHeader*>(dataBytes_);
        }

        io::FragmentAccessor &fragment()
        {
            return *reinterpret_cast<io::FragmentAccessor*>(dataBytes_);
        }

        const io::FragmentAccessor &fragment() const
        {
            return *reinterpret_cast<const io::FragmentAccessor*>(dataBytes_);
        }

        char *fragmentData()
        {
            return dataBytes_ + sizeof(io::FragmentHeader);
        }
    };

    typedef std::vector<IndexRecord>::const_iterator IndexConstIterator;
    typedef std::vector<IndexRecord>::iterator IndexIterator;

    void resize(const unsigned long clusters)
    {
        clusters_ = clusters;
        index_.resize(clusters_ * readOffsets_.size());
        data_.resize(clusters_ * recordLength_);
    }

    void reserve(const unsigned long clusters)
    {
        index_.reserve(clusters * readOffsets_.size());
        data_.reserve(clusters * recordLength_);
    }

    void clear()
    {
        index_.clear();
        data_.clear();
    }

    void unreserve()
    {
        std::vector<char>().swap(data_);
        std::vector<IndexRecord>().swap(index_);
    }

    bool empty() const
    {
        ISAAC_ASSERT_MSG(index_.empty() == data_.empty(), "index_ and data_ must be in sync");
        return data_.empty();
    }

    void swap(FragmentBuffer &another)
    {
        ISAAC_ASSERT_MSG(recordLength_ == another.recordLength_, "Buffers must be formatted identically");
        ISAAC_ASSERT_MSG(readOffsets_.size() == another.readOffsets_.size(), "Buffers must be formatted identically");
        ISAAC_ASSERT_MSG(readOffsets_[0] == another.readOffsets_[0], "Read offsets must match");
        ISAAC_ASSERT_MSG(readOffsets_.size() == 1 || readOffsets_[1] == another.readOffsets_[1], "Read offsets must match");
        using std::swap;
        swap(clusters_, another.clusters_);
        index_.swap(another.index_);
        data_.swap(another.data_);
    }

    IndexRecord &initialize(const unsigned long clusterId, const unsigned readIndex)
    {
        IndexRecord& ret = index_[clusterId * readOffsets_.size() + readIndex];

        ISAAC_ASSERT_MSG(!ret.dataBytes_, "Did not expect the same cluster read to be initialized twice");
        ret.dataBytes_ = &data_[clusterId * recordLength_ + readOffsets_[readIndex]];
        return ret;
    }

    const IndexRecord &getRecordStart(const unsigned long clusterId, const unsigned readIndex) const
    {
        const IndexRecord& ret = index_.at(clusterId * readOffsets_.size() + readIndex);

        const void * expected = &data_.at(clusterId * recordLength_ + readOffsets_[readIndex]);
        if(ret.dataBytes_ && ret.dataBytes_ != expected)
        {
            ISAAC_THREAD_CERR << "Improperly initialized index record for cluster: " << clusterId << " read " << readIndex <<
               " expected: " << expected << " actual: " << (void *)ret.dataBytes_ << std::endl;
            ISAAC_ASSERT_MSG(false, "Improperly initialized index record");
        }
        return ret;
    }

    unsigned long getClusters() const
    {
        return clusters_;
    }

    void sortIndex(const BinIndexMap &binIndexMap)
    {
        std::sort(index_.begin(), index_.end(), boost::bind(orderIndexByBin, _1, _2, boost::ref(binIndexMap)));
    }

    IndexIterator binBegin(
        size_t bin, const BinIndexMap &binIndexMap)
    {
        if (!bin)
        {
            // Bin 0 contains unaligned fragments and does not work in terms of positions.
            return index_.begin();
        }
        else if (bin <= binIndexMap.getHighestBinIndex())
        {
            // set lowest non-0 pointer to make the searched object look initialized() and sort to the start of the bin
            char *blah = 0;
            const IndexRecord binBegin(binIndexMap.getBinFirstPos(bin), ++blah);

            return lower_bound(index_.begin(), index_.end(), binBegin,
                               boost::bind(orderIndexByBin, _1, _2, boost::ref(binIndexMap)));
        }
        else
        {
            //ISAAC_THREAD_CERR << "requesting bin begin for " << bin << std::endl;
            ISAAC_ASSERT_MSG(binIndexMap.getHighestBinIndex() + 1 == bin, "Only one bin past the end can be requested here");
            return indexEnd();
        }
    }

    IndexIterator indexEnd()
    {
        return index_.end();
    }

private:
    static const unsigned READS_MAX = 2;
    const unsigned recordLength_;
    const common::FiniteCapacityVector<unsigned, 2> readOffsets_;
    unsigned long clusters_;
    std::vector<IndexRecord> index_;
    std::vector<char> data_;

    static bool orderIndexByBin(const IndexRecord &left, const IndexRecord &right, const BinIndexMap &binIndexMap)
    {
        if (!left.initialized())
        {
            if (right.initialized())
            {
                // push the uninitialized to the bottom. They will not pass to build stage
                return false;
            }
            else
            {
                // no ordering between two uninitialized
                return false;
            }
        }
        else if (!right.initialized())
        {
            // push the uninitialized to the bottom. They will not pass to build stage
            return true;
        }

        if (!left.fStrandPos_.isNoMatch())
        {
            if (right.fStrandPos_.isNoMatch())
            {
                // push unaligned (but not the shadow) ones to the top, they go to bin 0
                return false;
            }
            else //if both aligned
            {
                const size_t leftBin = binIndexMap.getBinIndex(left.fStrandPos_);
                const size_t rightBin = binIndexMap.getBinIndex(right.fStrandPos_);
                if (leftBin == rightBin)
                {
                    // If mates stay in the same bin, their offsets need to point at
                    // each other for realignment to work.
                    // this ensures the mates stay together which is important to
                    // updating offsets when they get loaded for bam generation
                    return &left.fragmentHeader() < &right.fragmentHeader();
                }
                return leftBin < rightBin;
            }
        }
        else if (!right.fStrandPos_.isNoMatch())
        {
            // push unaligned (but not the shadow) ones to the top, they go to bin 0
            return true;
        }
        else
        {
            // no ordering between two isNoMatch
            return false;
        }
    }

    static unsigned getRecordLength(const flowcell::FlowcellLayoutList &flowcellLayoutList)
    {
        const unsigned read2MaxLength = flowcell::getMaxReadLength(flowcellLayoutList, 1);
        return
            io::FragmentHeader::getMaxTotalLength(flowcell::getMaxReadLength(flowcellLayoutList, 0)) +
            (read2MaxLength ? (io::FragmentHeader::getMaxTotalLength(read2MaxLength)) : 0);
    }

    /**
     * \brief First read is located at the beginning, Second is at io::FragmentHeader::getMaxTotalLength of the
     *        first read length. Only two reads are supported
     */
    static common::FiniteCapacityVector<unsigned, 2> getReadOffsets(const flowcell::FlowcellLayoutList &flowcellLayoutList)
    {
        common::FiniteCapacityVector<unsigned, 2> ret;
        ret.push_back(0);
        const unsigned readsMax = flowcell::getMaxReadCount(flowcellLayoutList);
        if (READS_MAX == readsMax)
        {
            ret.push_back(io::FragmentHeader::getMaxTotalLength(flowcell::getMaxReadLength(flowcellLayoutList, 0)));
        }
        else
        {
            ISAAC_ASSERT_MSG(1 == readsMax, "Unexpected reads count " << readsMax);
        }
        return ret;
    }

};

/**
 ** \brief Buffering of fragments
 **
 ** - dataLength: total length on the block in bytes [unsigned int]
 ** - templatePosition: leftmost position of template [unsigned long]
 ** - fragmentLength: observed on reference [unsigned int]
 ** - readLength: number of bases in the read [unsigned short]
 ** - cigarLength: number of operations in cigar [unsigned short]
 ** - fragmentCount: number of fragments in template [unsigned byte]
 ** - fragmentIndex: 0-based index of this fragment [unsigned byte]
 ** - bclData: the read, in BCL format [readLength * unsigned byte]
 ** - cigar: the list of cigar operations [cigarLength * unsigned int]
 ** - allFragments: metadata for all fragments [fragmentCount*fragmentMetadata]
 **
 ** The fragment metadata is as follows:
 ** - alignmentQuality: [unsigned short] (see note below)
 ** - fragmentPosition: leftmost position on reference [unsigned long]
 **
 ** Note: the alignmentQuality encodes the orientation on the most significant
 ** bit.
 **/
class FragmentCollector
{
public:
    FragmentCollector(
        const BinIndexMap &binIndexMap,
        const unsigned long reserveClusters,
        const flowcell::FlowcellLayoutList &flowcellLayoutList)
    : binIndexMap_(binIndexMap),
      buffer_(reserveClusters, flowcellLayoutList)
    {
    }
    ~FragmentCollector();

    void resize(const unsigned long clusters)
    {
        buffer_.resize(clusters);
    }

    void add(
        const alignment::BamTemplate &bamTemplate,
        unsigned fragmentIndex, const unsigned barcodeIdx);

    void swapBuffer(FragmentBuffer &newBuffer)
    {
        buffer_.swap(newBuffer);
    }

private:
    const BinIndexMap &binIndexMap_;
    FragmentBuffer buffer_;
    void storeBclAndCigar(const alignment::FragmentMetadata & fragment, FragmentBuffer::IndexRecord & recordStart);
};

inline std::ostream &operator <<(std::ostream &os, const FragmentBuffer::IndexRecord &idxr)
{
    os << "FragmentBuffer::IndexRecord(" <<
        idxr.fStrandPos_ << "," << static_cast<const void*>(idxr.dataBytes_) << ")";
    return os;
}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_COLLECTOR_HH
