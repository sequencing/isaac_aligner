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
 ** \file FileSinkWithMd5.hh
 **
 ** \brief produces .md5 checksum file next to the data file.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FILE_SINK_WITH_MD5_HH
#define iSAAC_IO_FILE_SINK_WITH_MD5_HH

#include <boost/filesystem/path.hpp>
#include <boost/iostreams/device/file.hpp>

#include "common/MD5Sum.hh"

namespace isaac
{
namespace io
{

template<typename Ch>
struct BasicFileSinkWithMd5 : private boost::iostreams::basic_file<Ch> {
    typedef Ch char_type;
    struct category
        : boost::iostreams::output_seekable,
          boost::iostreams::device_tag,
          boost::iostreams::closable_tag,
          boost::iostreams::flushable_tag
        { };
    using boost::iostreams::basic_file<Ch>::seek;
    using boost::iostreams::basic_file<Ch>::is_open;
    using boost::iostreams::basic_file<Ch>::close;
    using boost::iostreams::basic_file<Ch>::flush;
    BasicFileSinkWithMd5( const std::string& path,
                     BOOST_IOS::openmode mode = BOOST_IOS::out )
        : boost::iostreams::basic_file<Ch>(path, mode & ~BOOST_IOS::in, BOOST_IOS::out),
          filePath_(path)
        { }
    ~BasicFileSinkWithMd5()
    {
    }
    void open( const std::string& path,
               BOOST_IOS::openmode mode = BOOST_IOS::out )
    {
        filePath_ = path;
        md5Sum_.clear();
        boost::iostreams::basic_file<Ch>::open(path, mode & ~BOOST_IOS::in, BOOST_IOS::out);
    }

    std::streamsize write(const char_type* s, std::streamsize n)
    {
        const std::streamsize ret = boost::iostreams::basic_file<Ch>::write(s, n);
        md5Sum_.update(s, n);
        return ret;
    }

    void close()
    {
        boost::iostreams::basic_file<Ch>::close();
        const std::string md5String = md5Sum_.getHexStringDigest();
        std::ofstream md5File((filePath_.string() + ".md5").c_str(), BOOST_IOS::out);
        md5File << md5String << " *" << filePath_.filename().string() << std::endl;
        if (!md5File)
        {
            BOOST_THROW_EXCEPTION(
                common::IoException(errno, (boost::format("Failed to write bytes into md5 stream for %s") % md5String.size() % filePath_.string()).str()));
        }
        ISAAC_THREAD_CERR << "md5 checksum for "  << filePath_.string() << ":" << md5String << std::endl;
    }

private:
    boost::filesystem::path filePath_;
    common::MD5Sum md5Sum_;
};

typedef BasicFileSinkWithMd5<char> FileSinkWithMd5;


} // namespace io
} // namespace isaac


#endif // iSAAC_IO_FILE_SINK_WITH_MD5_HH
