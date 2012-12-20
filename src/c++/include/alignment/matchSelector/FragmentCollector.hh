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
#include "io/FragmentIndex.hh"
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
            << recordLength_ << " read offsets : " << readOffsets_.at(0) << "," << readOffsets_.at(1) <<  std::endl;
    }

    class IndexRecord
    {
        // this constructor is only used for searching
        IndexRecord(reference::ReferencePosition fStrandPos, char *dataBytes)
            : fStrandPos_(fStrandPos), dataBytes_(dataBytes), mateDataBytes_(dataBytes)
        {
        }
        friend class FragmentBuffer;


        io::FragmentIndex &indexBase()
        {
            return *reinterpret_cast<io::FragmentIndex*>(dataBytes_);
        }

        io::FragmentIndex &mateIndexBase()
        {
            return *reinterpret_cast<io::FragmentIndex*>(mateDataBytes_);
        }

    public:
        IndexRecord() : fStrandPos_(0), dataBytes_(0), mateDataBytes_(0)
        {
        }

        static const size_t indexSize_ =
            boost::static_unsigned_max<boost::static_unsigned_max<sizeof(io::FStrandFragmentIndex),
                                       sizeof(io::RStrandOrShadowFragmentIndex)>::value,
                                       sizeof(io::SeFragmentIndex)>::value;

        reference::ReferencePosition fStrandPos_;
        char *dataBytes_;
        char *mateDataBytes_;

        bool initialized() const
        {
            return !!dataBytes_;
        }

        io::NmFragmentIndex &nmIndex()
        {
            return *reinterpret_cast<io::NmFragmentIndex*>(dataBytes_);
        }

        const io::NmFragmentIndex &nmIndex() const
        {
            return *reinterpret_cast<const io::NmFragmentIndex*>(dataBytes_);
        }

        io::SeFragmentIndex &seIndex()
        {
            return *reinterpret_cast<io::SeFragmentIndex*>(dataBytes_);
        }

        const io::SeFragmentIndex &seIndex() const
        {
            return *reinterpret_cast<const io::SeFragmentIndex*>(dataBytes_);
        }

        const io::PairEndIndex &peIndex() const
        {
            return *reinterpret_cast<io::PairEndIndex*>(dataBytes_);
        }

        io::FStrandFragmentIndex &fIndex()
        {
            return *reinterpret_cast<io::FStrandFragmentIndex*>(dataBytes_);
        }

        const io::FStrandFragmentIndex &fIndex() const
        {
            return *reinterpret_cast<const io::FStrandFragmentIndex*>(dataBytes_);
        }

        io::RStrandOrShadowFragmentIndex &rsIndex()
        {
            return *reinterpret_cast<io::RStrandOrShadowFragmentIndex*>(dataBytes_);
        }

        const io::RStrandOrShadowFragmentIndex &rsIndex() const
        {
            return *reinterpret_cast<const io::RStrandOrShadowFragmentIndex*>(dataBytes_);
        }

        io::FragmentHeader &fragmentHeader()
        {
            return *reinterpret_cast<io::FragmentHeader*>(dataBytes_ + indexSize_);
        }

        const io::FragmentHeader &fragmentHeader() const
        {
            return *reinterpret_cast<const io::FragmentHeader*>(dataBytes_ + indexSize_);
        }

        io::FragmentAccessor &fragment()
        {
            return *reinterpret_cast<io::FragmentAccessor*>(dataBytes_ + indexSize_);
        }

        const io::FragmentAccessor &fragment() const
        {
            return *reinterpret_cast<const io::FragmentAccessor*>(dataBytes_ + indexSize_);
        }

        char *fragmentData()
        {
            return dataBytes_ + indexSize_ + sizeof(io::FragmentHeader);
        }

        /**
         * \param updateMate if true, mateDataOffset_ in the mate is set too. Currently it is expected that
         *        updateMate is false for single-ended fragments or if the mate belongs to a different bin because
         *          a) bins are processed asynchronously on separate threads.
         *          b) Gap realignment uses offsets to update mate assuming it is within the same bin.
         */
        void setDataOffset(const unsigned long offset, const bool updateMate)
        {
            indexBase().dataOffset_ = offset;
            if (updateMate)
            {
                mateIndexBase().mateDataOffset_ = offset;
            }
            else
            {
                // otherwise, pretend that mate information is not accessible. Same will happen when setDataOffset is called for the mate.
                indexBase().mateDataOffset_ = offset;
            }
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
        ISAAC_ASSERT_MSG(readOffsets_.size() == 2, "Buffers must be for paired-ended data");
        ISAAC_ASSERT_MSG(readOffsets_.at(0) == another.readOffsets_.at(0), "Read offsets must match");
        ISAAC_ASSERT_MSG(readOffsets_.at(1) == another.readOffsets_.at(1), "Read offsets must match");
        using std::swap;
        swap(clusters_, another.clusters_);
        index_.swap(another.index_);
        data_.swap(another.data_);
    }

    IndexRecord &initialize(const unsigned long clusterId, const unsigned readIndex)
    {
        IndexRecord& ret = index_[clusterId * readOffsets_.size() + readIndex];

        ISAAC_ASSERT_MSG(!ret.dataBytes_, "Did not expect the same cluster read to be initialized twice");
        ret.dataBytes_ = &data_[clusterId * recordLength_ + readOffsets_.at(readIndex)];
        ret.mateDataBytes_ = &data_[clusterId * recordLength_ + readOffsets_.at((readIndex + 1) % readsMax_)];
        return ret;
    }

    const IndexRecord &getRecordStart(const unsigned long clusterId, const unsigned readIndex) const
    {
        const IndexRecord& ret = index_.at(clusterId * readOffsets_.size() + readIndex);

        const void * expected = &data_.at(clusterId * recordLength_ + readOffsets_.at(readIndex));
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
    unsigned getReads() const
    {
        return readOffsets_.size();
    }
//    unsigned long getDataStartOffset(const RecordLayout &record) const
//    {
//        return record.dataBytes_ - &data_.front();
//    }

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
        else if (bin < binIndexMap.getBinCount())
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
            ISAAC_ASSERT_MSG(binIndexMap.getBinCount() == bin, "Only one bin past the end can be requested here");
            return indexEnd();
        }
    }

    IndexIterator indexEnd()
    {
        return index_.end();
    }

private:
    static const unsigned readsMax_ = 2;
    const unsigned recordLength_;
    const std::vector<unsigned> readOffsets_;
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
            IndexRecord::indexSize_ + io::FragmentHeader::getMaxTotalLength(flowcell::getMaxReadLength(flowcellLayoutList, 0)) +
            (read2MaxLength ? (IndexRecord::indexSize_ + io::FragmentHeader::getMaxTotalLength(read2MaxLength)) : 0);
    }

    /**
     * \brief First read is located at the beginning, Second is at io::FragmentHeader::getMaxTotalLength of the
     *        first read length. Only two reads are supported
     */
    static std::vector<unsigned> getReadOffsets(const flowcell::FlowcellLayoutList &flowcellLayoutList)
    {
        std::vector<unsigned> ret = boost::assign::list_of
            (0u)
            (IndexRecord::indexSize_ + io::FragmentHeader::getMaxTotalLength(flowcell::getMaxReadLength(flowcellLayoutList, 0)));
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
        idxr.fStrandPos_ << "," <<
        static_cast<const void*>(idxr.dataBytes_) << "," <<
        static_cast<const void*>(idxr.mateDataBytes_) << ")";
    return os;
}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_COLLECTOR_HH
