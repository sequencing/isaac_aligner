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
 ** \file FiltersMapper.hh
 **
 ** Helper class for mapping Filter files into memory.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FILTERS_MAPPER_HH
#define iSAAC_IO_FILTERS_MAPPER_HH

#include "common/Debug.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace io
{

class FiltersMapper
{
public:
    FiltersMapper(const bool ignoreMissingFilterFiles) :
        ignoreMissingFilterFiles_(ignoreMissingFilterFiles),
        fileBufCache_(1, std::ios_base::in | std::ios_base::binary)
    {
    }
    void mapTile(const boost::filesystem::path &filtersFilePath, const unsigned clusterCount)
    {
        clusterCount_ = clusterCount;
        tileData_.clear();
        version_ = load(filtersFilePath, Autodetect);
    }

    unsigned getClusterCount() const {return clusterCount_;}

    template <typename InsertIteratorT>
    void getPf(InsertIteratorT it) const
    {
        switch (version_)
        {
        case V0:
            getPf(reinterpret_cast<const V0Header &>(tileData_.front()), it, clusterCount_);
            break;
        case V1:
            getPf(reinterpret_cast<const V1Header &>(tileData_.front()), it, clusterCount_);
            break;
        case V2:
            getPf(reinterpret_cast<const V2Header &>(tileData_.front()), it, clusterCount_);
            break;
        case V3:
            getPf(reinterpret_cast<const V3Header &>(tileData_.front()), it, clusterCount_);
            break;
        default:
            assert(false && "Autodetection is expected to give only the values listed above");
            break;
        }
    }

    void reservePathBuffers(const size_t reservePathLength)
    {
        fileBufCache_.reservePathBuffers(reservePathLength);
    }

    void reserveBuffer(const unsigned maxClusterCount)
    {
        tileData_.reserve(getMaxPossibleExpectedFileSize(maxClusterCount));
    }

    void unreserve()
    {
        std::vector<char>().swap(tileData_);
        fileBufCache_.unreserve();
    }
private:
    const bool ignoreMissingFilterFiles_;
    io::FileBufCache<io::FileBufWithReopen> fileBufCache_;
    unsigned clusterCount_;
    std::vector<char> tileData_;
    enum Version
    {
        Autodetect = -1,
        V0 = 0,
        FirstSupported = V0,
        V1 = 1,
        V2 = 2,
        V3 = 3,
        FirstUnsupported = 4
    };
    Version version_;


    struct V0Header
    {
        struct Header
        {
            unsigned clusters; //unsigned 32bits little endian integer: number N of clusters
        } header;

        typedef char value_type;
        value_type values[1]; //    * unsigned 8bits integer:
        //    * bit 0 is pass or failed
    };//__attribute__ ((packed));

    struct V1Header
    {
        typedef short value_type;
        struct Header
        {
            unsigned clusters; //unsigned 32bits little endian integer: number N of clusters
        } header;
        value_type values[1]; //bytes 4..(2*N+3): unsigned 16bits integer:
        //    * bit 0: pass or failed
        //    * bit 1: was the read identified as a control?
        //    * bit 2: was the match ambiguous?
        //    * bit 3: did the read match the phiX tag?
        //    * bit 4: did the read align to match the phiX tag?
        //    * bit 5: did the read match the control index sequence? (specified in controls.fata, TGTCACA)
        //    * bits 6,7: reserved for future use
        //    * bits 8..15: the report key for the matched record in the controls.fasta file (specified by the REPORT_KEY metadata)

    };//__attribute__ ((packed));

    struct V2Header
    {
        struct Header
        {
            unsigned header; //0x00000000 for disambiguation
            unsigned version; //version number = 0x00000002
            unsigned clusters; //unsigned 32bits little endian integer: number N of clusters
        } header;
        typedef short value_type;
        value_type values[1]; //unsigned 16bits integer:
        //    *  bit 0: pass or failed
        //    * bit 1: was the read identified as a control?
        //    * bit 2: was the match ambiguous?
        //    * bit 3: did the read match the phiX tag?
        //    * bit 4: did the read align to match the phiX tag?
        //    * bit 5: did the read match the control index sequence? (specified in controls.fata, TGTCACA)
        //    * bits 6,7: reserved for future use
        //    * bits 8..15: the report key for the matched record in the controls.fasta file (specified by the REPORT_KEY metadata)
    };//__attribute__ ((packed));

    struct V3Header
    {
        struct Header
        {
            unsigned header; //0x00000000 for disambiguation
            unsigned version; //version number = 0x00000003
            unsigned clusters; //unsigned 32bits little endian integer: number N of clusters
        } header;
        typedef char value_type;
        value_type values[1]; //unsigned 8bits integer:
        //    * bit 0: pass or failed
    };//__attribute__ ((packed));


    template <typename HeaderT>
    static unsigned getExpectedFileSize(const unsigned clusters) {
        return sizeof(typename HeaderT::Header) + clusters * sizeof(typename HeaderT::value_type);
    }

    template <typename HeaderT, typename InsertIteratorT>
    static void getPf(const HeaderT &header, InsertIteratorT it, const unsigned clusters) {
//        ISAAC_THREAD_CERR << "getPf " << clusters << " actual " << header.clusters << std::endl;

        assert(header.header.clusters == clusters);
        std::copy(header.values, header.values + clusters, it);
    }

    Version detectVersion(const boost::filesystem::path &filterFilePath, std::istream &is) const
    {
        Version assumedVersion(V0);
        unsigned int clusterCount(0);
        if (!is.read(reinterpret_cast<char *>(&clusterCount), 4))
        {
            BOOST_THROW_EXCEPTION(
                common::IoException(errno, (boost::format("Failed to read cluster count from filters file %s: %s") % filterFilePath % strerror(errno)).str()));
        }

        if (!clusterCount)
        {
            // V2 or V3
            unsigned int version(0);
            if (!is.read(reinterpret_cast<char *>(&version), 4))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read version from %s: %s") % filterFilePath % strerror(errno)).str()));
            }
            if (V2 != version && V3 != version)
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Unexpected filter file version (%d). from %s:") % version % filterFilePath ).str()));
            }
            assumedVersion = static_cast<Version>(version);
        }
        else
        {
            // V0 or V1
            const std::istream::pos_type valuesStartPos = is.tellg();
            if (!is.seekg(0, std::ios::end))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to seek to end of " + filterFilePath.string()));
            }
            const std::istream::pos_type fileEndPos = is.tellg();
            const unsigned long valueBytes = fileEndPos - valuesStartPos;
            if (clusterCount_ == valueBytes)
            {
                assumedVersion = V0;
            }
            else if (clusterCount_ == valueBytes / 2)
            {
                assumedVersion = V1;
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Unexpected file length (%d +%d) when detecting filter file format version for % clusters. from %s:") %
                    valuesStartPos % fileEndPos % clusterCount_ % filterFilePath ).str()));
            }
        }
        if (!is.seekg(0, std::ios::beg))
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to seek back to beginning of file " + filterFilePath.string()));
        }
        return assumedVersion;
    }

    Version load(const boost::filesystem::path &filterFilePath, Version assumedVersion)
    {
        if (!boost::filesystem::exists(filterFilePath))
        {
            if (!ignoreMissingFilterFiles_)
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, std::string("File does not exist: ") + filterFilePath.string()));
            }
            else
            {
                ISAAC_THREAD_CERR << "WARNING: ignoring missing filter file: " << filterFilePath << std::endl;

                if (Autodetect == assumedVersion)
                {
                    assumedVersion = V0;
                }
                unsigned expectedFileSize = getVersionExpectedFileSize(assumedVersion, clusterCount_);
                tileData_.resize(expectedFileSize);
                std::fill(tileData_.begin(), tileData_.end(), 1);
                reinterpret_cast<V0Header &>(tileData_.front()).header.clusters = clusterCount_;

                ISAAC_THREAD_CERR << "Assuming " << clusterCount_ << " clusters pass filter due to missing " << filterFilePath << std::endl;
            }
        }
        else
        {
            std::istream is(fileBufCache_.get(filterFilePath));
            if (!is)
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open %s: %s") % filterFilePath % strerror(errno)).str()));
            }

            if (Autodetect == assumedVersion)
            {
                assumedVersion = detectVersion(filterFilePath, is);
            }

            unsigned expectedFileSize = getVersionExpectedFileSize(assumedVersion, clusterCount_);

            tileData_.resize(expectedFileSize);
            if (!is.read(&tileData_.front(), expectedFileSize))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d bytes from %s") %
                    expectedFileSize % filterFilePath ).str()));
            }

            ISAAC_THREAD_CERR << "Read " << clusterCount_ << " filter values from filter file version " << assumedVersion << ": " << filterFilePath << std::endl;

        }
        return assumedVersion;
    }

    static unsigned getVersionExpectedFileSize(const Version assumedVersion, unsigned maxClusterCount)
    {
        unsigned expectedFileSize =
            V0 == assumedVersion
            ? getExpectedFileSize<V0Header>(maxClusterCount) :
            V1 == assumedVersion
            ? getExpectedFileSize<V1Header>(maxClusterCount) :
            V2 == assumedVersion
            ? getExpectedFileSize<V2Header>(maxClusterCount) : getExpectedFileSize<V3Header>(maxClusterCount);
        return expectedFileSize;
    }

    static unsigned getMaxPossibleExpectedFileSize(unsigned maxClusterCount)
    {
        unsigned expectedFileSize = 0;
        for (Version v = FirstSupported; v < FirstUnsupported; v = Version(v+1))
        {
            expectedFileSize = std::max(expectedFileSize, getVersionExpectedFileSize(v, maxClusterCount));
        }
        return expectedFileSize;
    }

};


} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FILTERS_MAPPER_HH
