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
 ** \file Cigar.hh
 **
 ** \brief Tools for creation, handling, management of BAM CIGAR
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_CIGAR_HH
#define iSAAC_ALIGNMENT_CIGAR_HH

#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>
#include <stdint.h>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "common/FastIo.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace alignment
{

class Cigar: public std::vector<uint32_t>
{
public:
    Cigar(const std::vector<uint32_t> &that) : std::vector<uint32_t>(that)
    {
    }

    Cigar(unsigned reservedSize = 0)
    {
        reserve(reservedSize);
    }
    enum OpCode {
        ALIGN = 0, // 'M'
        INSERT = 1, // 'I'
        DELETE = 2, // 'D'
        SKIP = 3, // 'N'
        SOFT_CLIP = 4, // 'S'
        HARD_CLIP = 5, // 'H'
        PAD = 6, // 'P'
        MATCH = 7, // '='
        MISMATCH = 8, // 'X'
        UNKNOWN = 9 // '?'
    };

    // we're in read length of hundreds. assume read length of thousands plus the operation code char
    static const unsigned OPERATION_CHARS_MAX = 5;
    void addOperation(unsigned length, OpCode opCode)
    {
        push_back(encode(length, opCode));
    }
    std::string toString() const;
    std::string toString(unsigned offset, unsigned length) const;
    static std::string toString(const std::vector<uint32_t> &cigarBuffer, unsigned offset, unsigned length);
    static char opCodeToChar(const Cigar::OpCode opCode)
    {
        static const char opCodes[] = {'M','I','D','N','S','H','P','=','X','?'};
        ISAAC_ASSERT_MSG(sizeof(opCodes) > std::size_t(opCode), "Unexpected CIGAR op code: " << opCode);
        return opCodes[opCode];
    }
    /**
     * \brief Serializes cigar to a container. Does not push terminating zero
     *
     * \return const reference to result
     */
    template <typename IteratorT, typename ContainerT>
    static const ContainerT &toString(IteratorT begin, IteratorT end, ContainerT &result)
    {
        for (IteratorT v = begin; end != v; ++v)
        {
            const std::pair<unsigned, Cigar::OpCode> d = Cigar::decode(*v);
            common::appendUnsignedInteger(result, d.first);
            result.push_back(opCodeToChar(d.second));
        }
        return result;
    }

    template <typename IteratorT>
    static std::string toString(IteratorT begin, IteratorT end)
    {
        std::string result;
        return toString(begin, end, result);
    }

    template <typename IteratorT>
    static std::ostream& toStream(IteratorT begin, IteratorT end, std::ostream &os)
    {
        for (IteratorT v = begin; end != v; ++v)
        {
            const std::pair<unsigned, Cigar::OpCode> d = Cigar::decode(*v);
            os << d.first << opCodeToChar(d.second);
        }
        return os;
    }

    template <typename IteratorT>
    static unsigned getReadLength(IteratorT begin, IteratorT end)
    {
        unsigned ret = 0;
        BOOST_FOREACH(const value_type v, std::make_pair(begin, end))
        {
            const Component d = Cigar::decode(v);
            switch (d.second)
            {
            case ALIGN:
            case INSERT:
            case SOFT_CLIP:
                ret += d.first;
                break;
            default:
                break;
            }
        }
        return ret;
    }

    template <typename IteratorT>
    static unsigned getMappedLength(IteratorT begin, IteratorT end)
    {
        unsigned ret = 0;
        BOOST_FOREACH(const value_type v, std::make_pair(begin, end))
        {
            const Component d = Cigar::decode(v);
            switch (d.second)
            {
            case ALIGN:
                ret += d.first;
                break;
            default:
                break;
            }
        }
        return ret;
    }


    static uint32_t encode(unsigned length, OpCode opCode)
    {
        assert(ALIGN <= opCode);
        assert(UNKNOWN >= opCode);
        return (length<<4) | opCode;
    }
    typedef std::pair<unsigned, OpCode> Component;
    static Component decode(const uint32_t value)
    {
        const unsigned length = value>>4;
        const unsigned code = std::min(value & 0x0F, static_cast<unsigned>(UNKNOWN));
        return std::pair<unsigned, OpCode>(length, static_cast<OpCode>(code));
    }


    static unsigned getMaxLength(const unsigned readLength)
    {
        return getMaxOpeations(readLength) * sizeof(value_type);
    }

    static unsigned getMinLength()
    {
        return 1 * sizeof(value_type);
    }

    static unsigned getMaxOpeations(const unsigned readLength)
    {
        const unsigned minBasesToFindAnIndel = 10; //assuming at least 10 bases are needed to identify an indel
        const unsigned maxCigarIndels = (readLength / minBasesToFindAnIndel);
        const unsigned cigarOpsPerIndel = 2; //assuming one indel requries one match
        const unsigned oneMatchOp = 1;
        const unsigned maxHardClipOps = 2; //one at each end
        const unsigned maxSoftClipOps = 2; //one at each end
        const unsigned maxCigarOperations =
            maxSoftClipOps + maxHardClipOps + oneMatchOp + maxCigarIndels * cigarOpsPerIndel;

        return maxCigarOperations;
    }

    static unsigned getMaxOpeationsForRead(const flowcell::ReadMetadata &readMetadata)
    {
        return getMaxLength(readMetadata.getLength());
    }

    static unsigned getMaxOperationsForReads(const flowcell::ReadMetadataList &readMetadataList)
    {
        return std::accumulate(readMetadataList.begin(), readMetadataList.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&Cigar::getMaxOpeationsForRead, _2)));
    }

    static unsigned getMaxOperationsForReads(const flowcell::FlowcellLayoutList &flowcellLayoutList)
    {
        unsigned ret = 0;
        BOOST_FOREACH(const flowcell::Layout &flowcell, flowcellLayoutList)
        {
            ret = std::max(ret, getMaxOperationsForReads(flowcell.getReadMetadataList()));
        }
        return ret;
    }
};

inline std::ostream &operator <<(std::ostream &os, const Cigar::Component &cigarComponent)
{
    return os << "CigarComponent(" << cigarComponent.first << "," << cigarComponent.second << ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_CIGAR_HH
