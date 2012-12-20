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
 ** \file FragmentReader.hh
 **
 ** \brief Reading of fragments data containing appropriate information for BAM generation.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_IO_FRAGMENT_READER_HH
#define iSAAC_IO_FRAGMENT_READER_HH

#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace io
{

namespace bfs = boost::filesystem;

class FragmentReader: public std::ifstream
{
public:
    FragmentReader(const bfs::path filePath);
    virtual ~FragmentReader();
    /**
     ** \brief Read the raw data for a single fragment and store it into the indicated buffer
     **
     ** Note: it is the responsibility of the calling code to ensure that enough
     ** spece is available in the indicated buffer.
     **
     ** TODO: make it a template to allow both pointers and iterators
     **/
    std::istream &read(char *buffer);
    /// Get the total storage size used by the raw fragment at the indicated buffer
    static unsigned getDataLength(char *buffer);
    /// Get the leftmost reference position of the associated BAM template
    static reference::ReferencePosition getTemplatePosition(char *buffer);
    /// Get the length of the fragment observed on the reference (as opposed to the read lenght)
    static unsigned getFragmentLength(char *buffer);
    /// Get the length of the read
    static unsigned short getReadLength();
    /// Get the CIGAR length
    static unsigned short getCigarLength();
    
private:
    const bfs::path filePath_;
};

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FRAGMENT_READER_HH
