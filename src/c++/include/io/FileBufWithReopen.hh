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
 ** \file FileBufWithReopen.hh
 **
 ** Same as std::filebuf but has reopen which is useful in parts of the code that are not
 ** allowed to do resource allocations.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FILE_BUF_WITH_REOPEN_HH
#define iSAAC_IO_FILE_BUF_WITH_REOPEN_HH

#include <fcntl.h>
#include <fstream>
#include <vector>

#include "common/Debug.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace io
{

extern const std::vector<const char*> iosBaseToStdioOpenModesTranslationTable;

template<typename _CharT, typename _Traits = std::char_traits<_CharT> >
class basic_FileBufWithReopen : public std::basic_filebuf<_CharT, _Traits>
{
    const std::ios_base::openmode mode_;

public:
    basic_FileBufWithReopen(std::ios_base::openmode mode) : mode_(mode)
    {
        if (!reserve())
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to allocate a file handle"));
        }
    }
    typedef std::basic_filebuf<_CharT, _Traits> base_type;
    typedef typename base_type::__filebuf_type __filebuf_type;

    enum FadviseFlags
    {
        normal = 0,
        sequential = 1,
        random = 2,
        noreuse = 4,
        willneed = 8,
        dontneed = 16,
        SequentialOnce = sequential|dontneed,
        SequentialOften = sequential|willneed,
    };
    /**
     * \brief opens the specified file with the access mode of the file that is currently open
     *        Supposedly does not fail due to memory allocation failure as it does not do any
     *        Use reserve on a new object before you can execute reopen
     */
    __filebuf_type* reopen(const char* s, FadviseFlags fadvise = normal)
    {
        const char * const openMode = iosFlagsToStdioMode(mode_);
        ISAAC_ASSERT_MSG(openMode, "Combination of open mode flags is invalid");
//        std::cerr << "reopening " << __s << " with mode " << openMode << "\n";
        ISAAC_ASSERT_MSG(this->is_open(), "The file must be already open for reopen to be possible");

        if ((mode_ & std::ios_base::out))
        {
            //flush any pending data before reopening a (usually) different file, reset eof flag
            ISAAC_ASSERT_MSG(0 == this->pubsync(), "pubsync failed with errno: " << errno << " " << strerror(errno));
        }


        if ((fadvise & noreuse) && posix_fadvise(fileno(this->_M_file.file()), 0, 0, POSIX_FADV_NOREUSE) &&
            errno)
        {
            // && errno check above is reqired since POSIX_FADV_NOREUSE fails with errno 0 on /dev/null
            ISAAC_THREAD_CERR << "WARNING: posix_fadvise failed for POSIX_FADV_NOREUSE with " << errno << "(" <<
                strerror(errno) << ")" << " file: " << s << std::endl;
        }
        if ((fadvise & willneed) && posix_fadvise(fileno(this->_M_file.file()), 0, 0, POSIX_FADV_WILLNEED) &&
            errno)
        {
            // && errno check above is reqired since POSIX_FADV_NOREUSE fails with errno 0 on /dev/null
            ISAAC_THREAD_CERR << "WARNING: posix_fadvise failed for POSIX_FADV_WILLNEED with " << errno << "(" <<
                strerror(errno) << ")" << " file: " << s << std::endl;
        }
        if ((fadvise & dontneed) && posix_fadvise(fileno(this->_M_file.file()), 0, 0, POSIX_FADV_DONTNEED) &&
            errno)
        {
            // && errno check above is reqired since POSIX_FADV_NOREUSE fails with errno 0 on /dev/null
            ISAAC_THREAD_CERR << "WARNING: posix_fadvise failed for POSIX_FADV_DONTNEED with " << errno << "(" <<
                strerror(errno) << ")" << " file: " << s << std::endl;
        }


        FILE* result = freopen(s, openMode, this->_M_file.file());


        if (result)
        {

            ISAAC_ASSERT_MSG(this->_M_file.file() == result, "According to specs, returned pointer must be the same as the one passed to freopen");
            if (!(mode_ & std::ios_base::app))
            {
                if(0 != this->seekpos(0, mode_))
                {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, s));
                }
            }


            if ((fadvise & sequential) && posix_fadvise(fileno(result), 0, 0, POSIX_FADV_SEQUENTIAL))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, s));
            }
            if ((fadvise & random) && posix_fadvise(fileno(result), 0, 0, POSIX_FADV_RANDOM))
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, s));
            }
            return this;
        }
        this->close();
        return 0;
    }

    std::ios_base::openmode mode() const {return mode_;}
    /**
     * \brief Reserves a file handle in a specified mode. This mode will be used during any subsequent reopen.
     *        Currently by opening /dev/null.
     */
    bool reserve()
    {
        // TODO: make it windows-compatible
        static const char * const fileThatAlwaysExists = "/dev/null";
        return !!this->open(fileThatAlwaysExists, mode_);
    }

    void flush()
    {
        if ((mode_ & std::ios_base::out))
        {
            //flush any pending data before reopening a (usually) different file
            //TODO: make some sensible error reporting...
            ISAAC_ASSERT_MSG(0 == this->pubsync(), "TODO: make some sensible error reporting...");
        }
    }

private:
    static const char * iosFlagsToStdioMode(std::ios_base::openmode mode)
    {
        const unsigned openModeIndex =
            !!(mode & std::ios_base::binary) << 4 |
            !!(mode & std::ios_base::in    ) << 3 |
            !!(mode & std::ios_base::out   ) << 2 |
            !!(mode & std::ios_base::trunc ) << 1 |
            !!(mode & std::ios_base::app   ) << 0
            ;

        return iosBaseToStdioOpenModesTranslationTable.at(openModeIndex);
    }
};

typedef basic_FileBufWithReopen<char> FileBufWithReopen;

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FILE_BUF_WITH_REOPEN_HH
