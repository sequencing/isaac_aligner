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
 ** \file BandedSmithWaterman.cpp
 **
 ** \brief See BandedSmithWaterman.hh
 **
 ** \author Come Raczy
 **/

#include <iostream>
#include <cassert>
#include <algorithm>
#include <xmmintrin.h> 
//#include <smmintrin.h> // sse4 instructions for _mm_insert_epi8

#include <boost/format.hpp>

#include "alignment/BandedSmithWaterman.hh"

namespace isaac
{
namespace alignment
{

BandedSmithWaterman::BandedSmithWaterman(const int matchScore, const int mismatchScore,
                                         const int gapOpenScore, const int gapExtendScore,
                                         const size_t maxReadLength)
    : matchScore_(matchScore)
    , mismatchScore_(mismatchScore)
    , gapOpenScore_(gapOpenScore)
    , gapExtendScore_(gapExtendScore)
    , maxReadLength_(maxReadLength)
    , T_((char *)_mm_malloc (maxReadLength * 3 *sizeof(__m128i), 16))
{
}

BandedSmithWaterman::~BandedSmithWaterman()
{
    _mm_free(T_);
}

// insert in register 0 only -- workaround for missing sse4 instruction set
inline __m128i _mm_insert_epi8(__m128i v, char c, int)
{
    const unsigned int tmp = _mm_extract_epi16(v, 0) & 0xffffff00u;
    const unsigned s = c & 0xff;
    return _mm_insert_epi16(v, tmp | s, 0);
}

// helper functions for debugging
std::string epi8(__m128i v);
std::string epi8s(__m128i v);
std::string epi8c(__m128i v);
std::ostream &operator<<(std::ostream &os, const __m128i &H);

unsigned BandedSmithWaterman::align(
    const std::vector<char> &query,
    const std::vector<char>::const_iterator databaseBegin,
    const std::vector<char>::const_iterator databaseEnd,
    Cigar &cigar) const
{
    return align(query.begin(), query.end(), databaseBegin, databaseEnd, cigar);
}

unsigned BandedSmithWaterman::align(
    const std::vector<char>::const_iterator queryBegin,
    const std::vector<char>::const_iterator queryEnd,
    const std::vector<char>::const_iterator databaseBegin,
    const std::vector<char>::const_iterator databaseEnd,
    Cigar &cigar) const
{
    assert(databaseEnd > databaseBegin);
    const size_t querySize = std::distance(queryBegin, queryEnd);
    assert(querySize + widestGapSize - 1 == (unsigned long)(databaseEnd - databaseBegin));
    assert(querySize <= maxReadLength_);
    const size_t originalCigarSize = cigar.size();
    //std::cerr << matchScore_ << " " << mismatchScore_ << " " << gapOpenScore_ << " " << gapExtendScore_  << std::endl;
    //std::cerr << "   " << query << std::endl;
    //std::cerr << database << std::endl;
    //__m128i *nextH = allH_;
    //allH_.clear();
    //allH_.reserve(query.length());
    //const __m128i GapOpenScore = _mm_insert_epi16(_mm_set1_epi16(gapOpenScore_), 0, 0);
    //const __m128i GapExtendScore = _mm_insert_epi16(_mm_set1_epi16(gapExtendScore_), 0, 0);
    __m128i *t = (__m128i *)T_;
    const __m128i GapOpenScore = _mm_set1_epi16(gapOpenScore_);
    const __m128i GapExtendScore = _mm_set1_epi16(gapExtendScore_);
    // Initialize E, F and G
    __m128i E[2], F[2], G[2];
    for (unsigned int i = 0; 2 > i; ++i)
    {
        E[i] = _mm_set1_epi16(-999);
        F[i] = _mm_set1_epi16(0); // Should this be -10000?
        G[i] = _mm_set1_epi16(-999);
    }
    G[0] = _mm_insert_epi16(G[0], 0, 0);
    // initialize D -- Note that we must leave the leftmost position empty (else it will be discarded before use)
    __m128i D = _mm_setzero_si128();
    for (unsigned int i = 0; widestGapSize > i + 1; ++i)
    {
        D = _mm_slli_si128(D, 1);
        D = _mm_insert_epi8(D, *(databaseBegin + i), 0);
    }
    //std::cerr << "    D:" << epi8(D) << "   " << database << std::endl;
    // iterate over all bases in the query
    //std::cerr << std::endl << database << std::endl << query << std::endl;
    std::vector<char>::const_iterator queryCurrent = queryBegin;
    for (unsigned queryOffset = 0; queryEnd != queryCurrent; ++queryOffset, ++queryCurrent)
    {
        __m128i tmp0[2], tmp1[2], tmp2[2];
        // F[i, j] = max(G[i-1, j] - open, E[i-1, j] - open, F[i-1, j] - extend)
        // G[i-1, j] and E[i-1, j]
        tmp0[1] = _mm_slli_si128(G[1], 2);
        tmp0[1] = _mm_insert_epi16(tmp0[1], _mm_extract_epi16(G[0], 7), 0);
        tmp0[0] = _mm_slli_si128(G[0], 2); // is the 0 initialisation alright?
        tmp1[1] = _mm_slli_si128(E[1], 2);
        tmp1[1] = _mm_insert_epi16(tmp1[1], _mm_extract_epi16(E[0], 7), 0);
        tmp1[0] = _mm_slli_si128(E[0], 2); // is the 0 initialisation alright?

        for (unsigned j = 0; 2 > j; ++j)
        {
            // identify which matrix provided the max (hack: true is 0xffff)
            tmp2[j] = _mm_cmplt_epi16(tmp0[j], tmp1[j]); // 0xffff if G[i-1, j] < E[i-1, j], 0 otherwise
            tmp2[j] = _mm_srli_epi16(tmp2[j], 15); // shift 15 bits to have only 1 for true
        }
        __m128i TF = _mm_packs_epi16(tmp2[0], tmp2[1]);
        // get max(G[i-1, j] - open, E[i-1, j] - open)
        __m128i newF[2];
        for (unsigned int j = 0; 2 > j; ++j)
        {
            newF[j] = _mm_max_epi16(tmp0[j], tmp1[j]);// newF = max(G[i-1, j], E[i-1, j])
            newF[j] = _mm_sub_epi16(newF[j], GapOpenScore); // newF = max(G[i-1, j] - open, E[i-1, j] - open)
        }
        // get F[i-1, j] - extend
        tmp0[1] = _mm_slli_si128(F[1], 2);
        tmp0[1] = _mm_insert_epi16(tmp0[1], _mm_extract_epi16(F[0], 7), 0);
        tmp0[1] = _mm_sub_epi16(tmp0[1], GapExtendScore);
        tmp0[0] = _mm_slli_si128(F[0], 2);
        tmp0[0] = _mm_sub_epi16(tmp0[0], GapExtendScore);
        // correct TF
        for (unsigned j = 0; 2 > j; ++j)
        {
            tmp2[j] = _mm_cmplt_epi16(newF[j], tmp0[j]); // 0xffff if max(G[i-1, j], E[i-1,j] - open < F[i-1, j] - extend
            //tmp2[j] = _mm_srli_epi16(tmp2[j], 14); // shift 14 bits to have 3 for true
            tmp2[j] = _mm_slli_epi16(_mm_srli_epi16(tmp2[j], 15), 1); // shift right 15 bits and left 1 bit to have 2 for true
        }
        TF = _mm_max_epu8(_mm_packs_epi16(tmp2[0], tmp2[1]), TF); // 0, 1, or 2 for G, E or F
        TF = _mm_insert_epi8(TF, 0, 0); // 0, 1, or 2 for G, E or F
        // correct F according to (F[i-1, j] - extend)
        for (unsigned int j = 0; 2 > j; ++j)
        {
            newF[j] = _mm_max_epi16(newF[j], tmp0[j]);
        }
        newF[0] = _mm_insert_epi16(newF[0], -999, 0);
        // G[i, j] = max(G[i-1, j-1], E[i-1, j-1], F[i-1, j-1]
        __m128i newG[2];
        for (unsigned int j = 0; 2 > j; ++j)
        {
            tmp2[j] = _mm_cmplt_epi16(G[j], E[j]); // 0xffff if G[i-1, j-1] < E[i-1, j-1], 0 otherwise
            tmp2[j] = _mm_srli_epi16(tmp2[j], 15); // shift 15 bits to have only 1 for true
            newG[j] = _mm_max_epi16(G[j], E[j]);// newG = max(G[i-1, j-1], E[i-1, j-1])
        }
        __m128i TG = _mm_packs_epi16(tmp2[0], tmp2[1]);
        // correct G and TG
        //std::cerr << "newG: " << newG[1] << newG[0] << std::endl;
        for (unsigned int j = 0; 2 > j; ++j)
        {
            tmp2[j] =  _mm_cmplt_epi16(newG[j], F[j]); // 0xffff if max(G[i-1, j-1], E[i-1,j-1] < F[i-1, j-1]
            tmp2[j] = _mm_slli_epi16(_mm_srli_epi16(tmp2[j], 15), 1); // shift right 15 bits and left 1 bit to have 2 for true
            newG[j] = _mm_max_epi16(newG[j], F[j]);
        }
        //std::cerr << "newG: " << newG[1] << newG[0] << std::endl;
        //std::cerr << "   G: " << G[1] << G[0] << std::endl;
        //std::cerr << "   E: " << E[1] <<  E[0] << std::endl;
        //std::cerr << "   F: " << F[1] <<  F[0] << std::endl;
        //std::cerr << "tmp0: " << tmp0[1] << tmp0[0] << std::endl;
        //std::cerr << "tmp1: " << tmp1[1] << tmp1[0] << std::endl;
        TG = _mm_max_epi16(_mm_packs_epi16(tmp2[0], tmp2[1]), TG); // 0, 1, or 2 for G, E or F
        // add the match/mismatch score
        // load the query base in all 8 values of the register
        __m128i Q = _mm_set1_epi8(*queryCurrent);
        // shift the database by 1 byte to the left and add the new base
        D = _mm_slli_si128(D, 1);
        D = _mm_insert_epi8(D, *(databaseBegin + queryOffset + widestGapSize - 1), 0);
        // compare query and database. 0xff if different (that also the sign bits)
        const __m128i B = ~_mm_cmpeq_epi8(Q, D);
#if 0
        //std::cerr << std::endl << database << std::endl << query << std::endl;;
        std::cerr << (boost::format("%2d  D:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8c(D) << std::endl;

        std::cerr << (boost::format("%2d  Q:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8c(Q) << std::endl;
        //std::cerr << (boost::format("%2d  B:") % (i+1)).str();
        //for (unsigned j = 0; j < i; ++j)
        //{
        //    std::cerr << "   0";
        //}
        //std::cerr << epi8(B) << std::endl;
        //exit(1);
#endif
        // set match/mismatch scores, according to comparison
        const __m128i Match = _mm_andnot_si128(B, _mm_set1_epi8(matchScore_));
        const __m128i Mismatch = _mm_and_si128(B, _mm_set1_epi8(mismatchScore_));
        // add the match/mismatch scored to HH
        const __m128i W = _mm_add_epi8(Match, Mismatch);

#if 0
        std::cerr << (boost::format("%2d  W:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8s(W) << std::endl;
#endif
        newG[0] = _mm_add_epi16(newG[0], _mm_unpacklo_epi8(W, B));
        newG[1] = _mm_add_epi16(newG[1], _mm_unpackhi_epi8(W, B));
        // E[i,j] = max(G[i, j-1] - open, E[i, j-1] - extend, F[i, j-1] - open)
        __m128i TE = _mm_setzero_si128();
        // E should never be the maximum in the leftmost side of the window
        short g = -999; // -gapOpenScore_;
        short e = -999; // -gapExtendScore_;
        short f = -999; // -gapOpenScore_;
        for (unsigned j = 0; 2 > j; ++j)
        {
            tmp0[j] = _mm_sub_epi16(newG[j], GapOpenScore);
            tmp1[j] = _mm_sub_epi16(newF[j], GapOpenScore);
        }
        //std::cerr << "newG: " << newG[1] << newG[0] << std::endl;
        //std::cerr << "tmp0: " << tmp0[1] << tmp0[0] << std::endl;
        //std::cerr << "newF: " << newF[1] << newF[0] << std::endl;
        //std::cerr << "tmp1: " << tmp1[1] << tmp1[0] << std::endl;
        //std::cerr <<  "WIDEST_GAP_SIZE = " << WIDEST_GAP_SIZE << std::endl;
        for (unsigned int j = 0; j < widestGapSize; ++j)
        {
            //std::cerr << "---- j = " << j << std::endl;
            if (8 == j)
            {
                //std::cerr << "   E: " << E[1] << E[0] << std::endl;
                //std::cerr << "tmp0: " << tmp0[1] << tmp0[0] << std::endl;
                //std::cerr << "tmp1: " << tmp1[1] << tmp1[0] << std::endl;
                tmp0[1] = tmp0[0];
                tmp1[1] = tmp1[0];
                E[1] = E[0];
            }
            short max = g;
            short tMax = 0;
            if (e > g && e > f)
            {
                max = e;
                tMax = 1;
            }
            else if (f > g)
            {
                max = f;
                tMax = 2;
            }
            //std::cerr << (boost::format("g = %d, e = %d, f = %d, max = %d, tMax = %d") % g % e % f % max % tMax ).str() << std::endl;
            TE = _mm_slli_si128(TE, 1);
            TE = _mm_insert_epi8(TE, tMax, 0);
            E[0] = _mm_slli_si128(E[0], 2);
            E[0] = _mm_insert_epi16(E[0], max, 0);
            g = _mm_extract_epi16(tmp0[1], 7);
            e = max - gapExtendScore_;
            f = _mm_extract_epi16(tmp1[1], 7);
            tmp0[1] = _mm_slli_si128(tmp0[1], 2);
            tmp1[1] = _mm_slli_si128(tmp1[1], 2);
            //std::cerr << "   E: " << E[1] << E[0] << std::endl;
            //std::cerr << "  TE: " << epi8(TE) << std::endl;
        }
        //E[0] = _mm_slli_si128(E[0], 2);
        for (unsigned j = 0; 2 > j; ++j)
        {
            G[j] = newG[j];
            F[j] = newF[j];
        }
        // TODO: add support for databases shorter than query + widestGapSize
        // store the matrix types
        _mm_store_si128(t++, TG);
        _mm_store_si128(t++, TE);
        _mm_store_si128(t++, TF);
#if 0
        std::cerr << (boost::format("%2d  G:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << G[1] << G[0] << std::endl;
        std::cerr << (boost::format("%2d TG:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8(TG) << std::endl;
        std::cerr << (boost::format("%2d  E:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << E[1] << E[0] << std::endl;
        std::cerr << (boost::format("%2d TE:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8(TE) << std::endl;
        std::cerr << (boost::format("%2d  F:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << F[1] << F[0] << std::endl;
        std::cerr << (boost::format("%2d TF:") % (queryOffset+1)).str();
        for (unsigned j = 0; j < queryOffset; ++j)
        {
            std::cerr << " 0";
        }
        std::cerr << epi8(TF) << std::endl;
#endif
    }
    // find the max of E, F and G at the end
    short max = _mm_extract_epi16(G[1], 7) - 1;
    int ii = querySize - 1;
    int jj = ii;
    unsigned maxType = 0;
    __m128i *TT[] = {G, E, F};
    for (unsigned int k = 0; 2 > k; ++k)
    {
        const unsigned int kk = (k + 1) % 2;
        //std:: cerr << "kk = " << kk << ", max = " << max << std::endl;
        //std::cerr << "type 0: " << TT[0][kk] << std::endl;
        //std::cerr << "type 1: " << TT[1][kk] << std::endl;
        //std::cerr << "type 2: " << TT[2][kk] << std::endl;
        for (unsigned j = 8; 0 < j; --j)
        {
            //std::cerr << "j = " << j;
            for (unsigned type = 0; 3 > type; ++type)
            {
                const short value = _mm_extract_epi16(TT[type][kk], 7);
                //std::cerr << ", type = " << type << ", value = " << value;
                TT[type][kk] = _mm_slli_si128(TT[type][kk], 2);
                if (value > max)
                {
                    //std::cerr << std::endl << "max = " << max << ", value = " << value<< ", j = "  << j<< ", kk = "  << kk << std::endl;
                    max = value;
                    jj = 8 * kk + j - 1;
                    maxType = type;
                }
            }
            //std::cerr << std::endl;
        }
    }
    //std::cerr << (boost::format("ii = %d, jj = %d, max = %d, maxType = %d") % (ii + 1) % jj % max % maxType).str() << std::endl;
    const int jjIncrement[] = {0, 1, -1};
    const int iiIncrement[] = {-1, 0, -1};
    const Cigar::OpCode opCodes[] = {Cigar::ALIGN, Cigar::DELETE, Cigar::INSERT};
    unsigned opLength = 0;
    //std::string newDatabase;
    //std::string newQuery;
    //std::cerr << std::endl << "building CIGAR" << std::endl;
    if (jj > 0)
    {
        cigar.addOperation(jj, Cigar::DELETE);
    }
    while(ii >= 0 && jj >= 0 && jj <= 15)
    {
#if 0
        if (0 == maxType)
        {
            newDatabase.push_back(database[ii + 15 - jj]);
            newQuery.push_back(query[ii]);
        }
        else if(1 == maxType)
        {
            newDatabase.push_back(database[ii + 15 - jj]);
            newQuery.push_back('-');
        }
        else if(2 == maxType)
        {
            newDatabase.push_back('-');
            newQuery.push_back(query[ii]);
        }
#endif
        ++opLength;
        const unsigned nextMaxType = T_[(ii * 3 + maxType) * sizeof(__m128i) + jj];
        //std::cerr << (boost::format("ii = %d, jj = %d, maxType = %d, nextMaxType = %d, opLength = %d") %
        //              ii % jj % maxType %nextMaxType % opLength ).str() << std::endl;
        if (nextMaxType != maxType)
        {
            cigar.addOperation(opLength, opCodes[maxType]);
            opLength = 0;
        }
        ii += iiIncrement[maxType];
        jj += jjIncrement[maxType];
        maxType = nextMaxType;
    }
    assert(-1 == ii);
    if (1 != maxType && opLength)
    {
        cigar.addOperation(opLength, opCodes[maxType]);
        opLength = 0;
    }
    if (15 > jj)
    {
        cigar.addOperation(opLength + 15 - jj, Cigar::DELETE);
        opLength = 0;
    }
    assert(0 == opLength);
    //descriptor.push_back(d[maxType]);
    const std::pair<unsigned, Cigar::OpCode> firstCigar = Cigar::decode(cigar.back());
    if(Cigar::DELETE == firstCigar.second)
    {
        //CASAVA does not like CIGAR beginning with a deletion in the data
        cigar.pop_back();
        ISAAC_ASSERT_MSG(Cigar::DELETE != Cigar::decode(cigar.back()).second, "two Cigar::DELETE cannot be next to each other");
    }
    std::reverse(cigar.begin() + originalCigarSize, cigar.end());
    if(Cigar::DELETE == Cigar::decode(cigar.back()).second)
    {
        //CASAVA does not like CIGAR ending with a deletion in the data
        cigar.pop_back();
        ISAAC_ASSERT_MSG(Cigar::DELETE != Cigar::decode(cigar.back()).second, "two Cigar::DELETE cannot be next to each other");
    }
    return firstCigar.first;
    //std::cerr << std::endl << "CIGAR: " << cigar.toString() << std::endl;
    //std::reverse(newDatabase.begin(), newDatabase.end());
    //std::reverse(newQuery.begin(), newQuery.end());
    //std::cerr << "       " << descriptor << std::endl;
    //std::cerr << database << std::endl;
    //std::cerr << "       " << query << std::endl;
    //std::cerr << "       " << newQuery << std::endl;
    //std::cerr << "       " << newDatabase << std::endl;
}

std::string epi8(__m128i v)
{
    std::string result;
    for (unsigned int i = 0; 16 > i; ++i)
    {
        result += (boost::format("%3d") % (_mm_extract_epi16(v, 7) >> 8)).str();
        v = _mm_slli_si128(v, 1);
    }
    return result;
}

std::string epi8s(__m128i v)
{
    std::string result;
    for (unsigned int i = 0; 16 > i; ++i)
    {
        result += (boost::format("%3d") % ((int)(_mm_extract_epi16(v, 7)) >> 8)).str();
        v = _mm_slli_si128(v, 1);
    }
    return result;
}

std::string epi8c(__m128i v)
{
    std::string result;
    for (unsigned int i = 0; 16 > i; ++i)
    {
        result += (boost::format("  %c") % char(_mm_extract_epi16(v, 7) >> 8)).str();
        v = _mm_slli_si128(v, 1);
    }
    return result;
}

std::ostream &operator<<(std::ostream &os, const __m128i &H)
{   
    __m128i tmp0 = H;
    for (unsigned int i = 0; i < 8; ++i)
    {
        short v = _mm_extract_epi16(tmp0, 7);
        std::cerr << (boost::format("%3d") % std::max(short(-99), v)).str();
        tmp0 = _mm_slli_si128(tmp0, 2);
    }
    return os;
}

} // namespace alignment
} // namespace isaac
