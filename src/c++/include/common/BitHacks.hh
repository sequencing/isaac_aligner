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
 ** \file BitHacks.hh
 **
 ** \brief http://graphics.stanford.edu/~seander/bithacks.html
 **
 ** \author Roman Petrovski
 **/

/**
 * \return next power of two value that is >= v
 */
inline unsigned long upperPowerOfTwo(unsigned long v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
}

inline unsigned countBitsSet(unsigned v)
{
    v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
    const unsigned c = (((v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24; // count
    return c;
}

// find the number of trailing zeros in 32-bit v
inline int lsbSet(const unsigned v)
{
    static const int MultiplyDeBruijnBitPosition[32] =
    {
      0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
      31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    };
    const int r = MultiplyDeBruijnBitPosition[((uint32_t)((v & -v) * 0x077CB531U)) >> 27];
    return r;
}
