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
 ** \file MatchReader.cpp
 **
 ** Abstract component to read matches.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_IO_MATCH_READER_HH
#define iSAAC_IO_MATCH_READER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>

#include "alignment/Match.hh"
#include "common/Exceptions.hh"
#include "io/FileBufCache.hh"
#include "io/FileBufWithReopen.hh"

namespace isaac
{
namespace io
{

namespace bfs = boost::filesystem;
/**
 ** \brief a component that reads the matches from a single file.
 **
 **/
class MatchReader
{
public:
    MatchReader() : fileBuf_(1, std::ios_base::binary|std::ios_base::in)
    {}

    // throws on any failure (including eof())
    void read(const bfs::path &matchFilePath, alignment::Match *destination, unsigned long count)
    {
        if (!std::istream(fileBuf_.get(matchFilePath, FileBufWithReopen::SequentialOnce)).read(reinterpret_cast<char *>(destination), sizeof(alignment::Match) * count))
        {
            BOOST_THROW_EXCEPTION(common::IoException(
                errno, (boost::format("Failed to read %u (%u bytes) matches from file %s") %
                    count % (sizeof(alignment::Match) * count) % matchFilePath.string()).str()));
        }
    }

    void reservePathBuffers(const size_t reservePathLength)
    {
        fileBuf_.reservePathBuffers(reservePathLength);
    }

private:
    FileBufCache<FileBufWithReopen> fileBuf_;
};

} //namespace io
} //namespace isaac

#endif // #ifndef iSAAC_IO_MATCH_READER_HH
