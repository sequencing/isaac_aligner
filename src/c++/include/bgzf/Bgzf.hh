/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file Bgzf.hh
 **
 ** \brief structures and constants required tow work with bgzf-compressed data.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BGZF_BGZF_HH
#define iSAAC_BGZF_BGZF_HH

#include "common/Endianness.hh"

namespace isaac
{
namespace bgzf
{

struct BAM_XFIELD
{
    unsigned char XLEN[2];

    unsigned char SI1;
    unsigned char SI2;
    unsigned char SLEN[2];
    unsigned char BSIZE[2];

    unsigned getXLEN() const {return unsigned(XLEN[0]) + XLEN[1] * 256;}
    unsigned getBSIZE() const {return unsigned(BSIZE[0]) + BSIZE[1] * 256;}

} __attribute__ ((packed));
BOOST_STATIC_ASSERT(8 == sizeof(BAM_XFIELD));

struct Header
{
    unsigned char ID1;
    unsigned char ID2;
    unsigned char CM;
    unsigned char FLG;
    unsigned char MTIME[4];
    unsigned char XFL;
    unsigned char OS;

    BAM_XFIELD xfield;

    unsigned getCDATASize() const
    {
        return xfield.getBSIZE() - xfield.getXLEN() - 19;
    }

} __attribute__ ((packed));
BOOST_STATIC_ASSERT(18 == sizeof(Header));

struct Footer
{
    unsigned char CRC32[4];
    unsigned char ISIZE[4];

    unsigned getISIZE() const
    {unsigned ret; common::extractLittleEndian(ISIZE, ret); return ret;
    /*unsigned(ISIZE[0]) + ISIZE[1] * 256 + ISIZE[2] * 256 * 256 + ISIZE[3] * 256 * 256 * 256;*/}

} __attribute__ ((packed));
BOOST_STATIC_ASSERT(8 == sizeof(Footer));


} // namespace bgzf
} // namespace isaac


#endif // iSAAC_BGZF_BGZF_HH
