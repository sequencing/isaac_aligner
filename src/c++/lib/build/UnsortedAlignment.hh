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
 ** \file UnsortedAlignment.hh
 **
 ** Helper class for accessing the unsorted alignment data
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_UNSORTED_ALIGNMENT_HH
#define iSAAC_BUILD_UNSORTED_ALIGNMENT_HH

namespace isaac
{
namespace build
{

namespace bip=boost::interprocess;
namespace bfs=boost::filesystem;

struct UnsortedAlignment
{
    unsigned long alignmentPos_;
    int templateLength_;
    unsigned short pras_;
    unsigned short sras_;
    unsigned short mateSras_;

    bool isSingleton() const {
        return sras_ && !mateSras_;
    }

    bool isShadow() const {
        return !sras_ && mateSras_;
    }

    char read_name[11];
    unsigned cigar[1];
    unsigned char seq[4];
    char qual[8];
    const UnsortedAlignment *next() const
    {
        return reinterpret_cast<const UnsortedAlignment *>(reinterpret_cast<const char*>(this) + sizeof(UnsortedAlignment));
    }
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_UNSORTED_ALIGNMENT_HH
