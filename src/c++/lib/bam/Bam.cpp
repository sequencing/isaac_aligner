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
 ** \file Bam.cpp
 **
 ** \brief collection of helpers for bam serialization.
 **
 ** \author Roman Petrovski
 **/

#include "bam/Bam.hh"


namespace isaac
{
namespace bam
{

void serialize(std::ostream &os, const char* bytes, size_t size)
{
//    std::cerr << "writing: " << size << " bytes\n";
    if (!os.write(bytes, size))
    {
        BOOST_THROW_EXCEPTION(
            common::IoException(errno, (boost::format("Failed to write %d bytes into bam stream") % size).str()));
    }
}

void serializeBgzfFooter(std::ostream &os)
{
    // For some strange reason, samtools wants an empty block at the very end of a compressed bam.
    // They call it 'magic'. Note, last \0 is removed (comparing to the oritinal bgzf.c) because
    // the C++ compiler (rightfully) complains.
    const static char magic[28] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0";
    serialize(os, magic, sizeof(magic));
}

} //namespace bam
} // namespace isaac
