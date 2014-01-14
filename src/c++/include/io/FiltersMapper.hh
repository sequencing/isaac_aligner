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
        fileBufCache_(1, std::ios_base::in | std::ios_base::binary),
        clusterCount_(0),
        version_(FirstUnsupported)
    {
    }

    /**
     * \param clusterOffset non-zero value allows extracting tile information from correct offset in lane filter file.
     *                      Tile filter file requests must have it set to 0
     */
    void mapTile(
        const boost::filesystem::path &filtersFilePath,
        const unsigned clusterCount,
        const unsigned long clusterOffset = ONE_TILE_PER_FILE)
    {
        clusterCount_ = clusterCount;
        tileData_.clear();
        version_ = load(filtersFilePath, clusterOffset, Autodetect);
    }

    template <typename InsertIteratorT>
    void getPf(InsertIteratorT it) const
    {
        versionSpecific<GetPfAction>(version_, boost::cref(tileData_), it, clusterCount_, UNUSED);
    }

    void reserveBuffers(const size_t reservePathLength, const unsigned maxClusterCount)
    {
        fileBufCache_.reservePathBuffers(reservePathLength);
        tileData_.reserve(getMaxPossibleExpectedFileSize(maxClusterCount));
    }

private:
    typedef boost::error_info<struct tag_errmsg, std::string> errmsg_info;
    static const unsigned long ONE_TILE_PER_FILE = -1UL;
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

    typedef const unsigned UnusedT;
    static UnusedT UNUSED = 0;
    template <template <typename HeaderT> class ActionT,
        typename Arg0T, typename Arg1T, typename Arg2T, typename Arg3T>
    static void versionSpecific(unsigned version, Arg0T arg0, Arg1T arg1, Arg2T arg2, Arg3T arg3)
    {
        if (V0 == version)
        {
            boost::bind(ActionT<V0Header>(), arg0 ,arg1, arg2, arg3)();
        }
        else if (V1 == version)
        {
            boost::bind(ActionT<V1Header>(), arg0 ,arg1, arg2, arg3)();
        }
        else if (V2 == version)
        {
            boost::bind(ActionT<V2Header>(), arg0 ,arg1, arg2, arg3)();
        }
        else if (V3 == version)
        {
            boost::bind(ActionT<V3Header>(), arg0 ,arg1, arg2, arg3)();
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected versionSpecific call for unknown version " << version);
        }
    }


    template <typename HeaderT>
    static std::size_t getExpectedFileSize(const unsigned clusters)
    {
        return sizeof(typename HeaderT::Header) + clusters * sizeof(typename HeaderT::value_type);
    }

    template <typename HeaderT>
    struct GetExpectedFileSizeAction
    {
        typedef void result_type;
        void operator()(const unsigned clusters, unsigned &size, UnusedT, UnusedT) const
        {
            size = getExpectedFileSize<HeaderT>(clusters);
        }
    };

    static unsigned getVersionExpectedFileSize(const Version assumedVersion, unsigned maxClusterCount)
    {
        unsigned expectedFileSize = 0;
        versionSpecific<GetExpectedFileSizeAction>(assumedVersion, maxClusterCount, boost::ref(expectedFileSize), UNUSED, UNUSED);
        return expectedFileSize;
    }

    template <typename HeaderT>
    struct ReadDataAction
    {
        typedef void result_type;
        void operator()(
            std::istream &is,
            std::vector<char> &tileData,
            const unsigned clusterCount,
            const unsigned long clusterOffset) const
        {
            unsigned expectedFileSize = getExpectedFileSize<HeaderT>(clusterCount);
            tileData.resize(expectedFileSize);

            if (!is.read(&tileData.front(), sizeof(HeaderT)))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d header bytes") %
                    sizeof(HeaderT)).str()));
            }

            if (ONE_TILE_PER_FILE != clusterOffset)
            {
                // patch the clusters number as multitile filter files contain the total number of clusters
                HeaderT &header = reinterpret_cast<HeaderT&>(tileData.front());
                header.header.clusters = clusterCount;
                const unsigned long clusterByteOffset = clusterOffset * sizeof(typename HeaderT::value_type);
                if (!is.seekg(clusterByteOffset, is.cur))
                {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to seek %d bytes") %
                        clusterByteOffset).str()));
                }
            }

            if (!is.read(&tileData.front() + sizeof(HeaderT), expectedFileSize - sizeof(HeaderT)))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %d bytes") %
                    expectedFileSize).str()));
            }
        }
    };

    template <typename HeaderT>
    struct GetPfAction
    {
        typedef void result_type;
        template<typename InsertIteratorT>
        void operator()(const std::vector<char> &tileData, InsertIteratorT it, const unsigned clusters, UnusedT) const
        {
            const HeaderT &header = reinterpret_cast<const HeaderT&>(tileData.front());
            ISAAC_ASSERT_MSG(header.header.clusters == clusters, "Requested number of pf values (" << clusters <<
                ") does not match the loaded:" << unsigned(header.header.clusters));
            std::copy(header.values, header.values + clusters, it);
        }
    };

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


    Version load(const boost::filesystem::path &filterFilePath, const unsigned long clusterOffset, Version assumedVersion)
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

            try
            {
                versionSpecific<ReadDataAction>(assumedVersion, boost::ref(is), boost::ref(tileData_), clusterCount_, clusterOffset);
            }
            catch (boost::exception &e)
            {
                e << errmsg_info(" While reading from " + filterFilePath.string());
                throw;
            }
            ISAAC_THREAD_CERR << "Read " << clusterCount_ << " filter values from filter file version " << assumedVersion << ": " << filterFilePath << std::endl;
        }
        return assumedVersion;
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
