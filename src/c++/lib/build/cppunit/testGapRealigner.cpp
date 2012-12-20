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
 ** \file testGapRealigner.cpp
 **
 ** Test cases for gap realignment.
 **
 ** \author Roman Petrovski
 **/

#include <string>

#include "build/GapRealigner.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"

using namespace isaac;


#include "RegistryName.hh"
#include "testGapRealigner.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestGapRealigner, registryName("TestGapRealigner"));

TestGapRealigner::TestGapRealigner()
{


}

void TestGapRealigner::setUp()
{
}

void TestGapRealigner::tearDown()
{
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

template <typename GapRealignerT>
void setInsertions(const std::string &ref, GapRealignerT &realigner)
{
    typename GapRealignerT::GapType currentGap(reference::ReferencePosition(0, 0), 0);
    BOOST_FOREACH(const char refChar, ref)
    {
        if ('*' == refChar)
        {
            currentGap.length_ -= 1;
        }
        else
        {
            if (currentGap.length_)
            {
                realigner.addGap(0, currentGap);
                currentGap.length_ = 0;
            }
            currentGap.pos_ += 1;
//            ISAAC_THREAD_CERR << "currentGap=" << currentGap << std::endl;
        }
    }

    if (currentGap.length_)
    {
        realigner.addGap(0, currentGap);
    }
}

template <typename GapRealignerT>
void setDeletions(const std::string &read, GapRealignerT &realigner)
{
    typename GapRealignerT::GapType currentGap(reference::ReferencePosition(0, 0), 0);
    BOOST_FOREACH(const char readChar, read)
    {
        if ('-' != readChar)
        {
            if (currentGap.length_)
            {
                typename GapRealignerT::GapType deletion(currentGap);
                deletion.pos_ -= currentGap.length_;
                realigner.addGap(0, deletion);
                currentGap.length_ = 0;
            }
        }
        else
        {
            currentGap.length_ += 1;
        }
        currentGap.pos_ += 1;
    }

    if (currentGap.length_)
    {
        typename GapRealignerT::GapType deletion(currentGap);
        deletion.pos_ -= currentGap.length_;
        realigner.addGap(0, deletion);
    }
}

struct TestFragmentAccessor : public io::FragmentAccessor
{
    static const unsigned maxReadLength_ = 1000;
    static const unsigned maxCigarLength_ = 1000;
    unsigned char buffer_[maxReadLength_ + maxCigarLength_ * sizeof(unsigned)];
    TestFragmentAccessor(
        const alignment::FragmentMetadata &fragment,
        const reference::ReferencePosition fStrandPosition,
        const std::string &read,
        const alignment::Cigar &cigar,
        const unsigned short editDistance)
    {
        CPPUNIT_ASSERT(maxCigarLength_ > unsigned(cigar.size()));
        CPPUNIT_ASSERT(maxReadLength_ > unsigned(read.size()));

        fStrandPosition_ = fStrandPosition;

        unsigned char *b = const_cast<unsigned char*>(basesBegin());
        BOOST_FOREACH(unsigned char base,
                      std::make_pair(read.begin(), read.end()))
        {
            if (std::string::npos != std::string("ACGTN").find(base))
            {
                *b++ = oligo::getValue(base);
                ++readLength_;
            }
        }
        std::copy(cigar.begin(), cigar.end(), const_cast<unsigned*>(cigarBegin()));
        cigarLength_ = cigar.size();
        editDistance_ = editDistance;
        flags_.mateBinTheSame_ = true;
        alignmentScore_ = 1;
        lowClipped_ = fragment.lowClipped;
        highClipped_ = fragment.highClipped;
    }

    const char *begin() const {return reinterpret_cast<const char *>(this);}
    const char *end() const {return reinterpret_cast<const char *>(cigarEnd());}
};


struct RealignResult
{
    RealignResult(
        const reference::ReferencePosition fStrandPosition,
        const std::string &originalCigar, const unsigned short editDistance):
            originalPos_(fStrandPosition),
            originalCigar_(originalCigar),
            originalEditDistance_(editDistance){}
    reference::ReferencePosition originalPos_;
    reference::ReferencePosition realignedPos_;
    std::string originalCigar_;
    std::string realignedCigar_;
    unsigned short originalEditDistance_;
    unsigned short realignedEditDistance_;
};

TestFragmentAccessor initFragment(
    const alignment::FragmentMetadata &fragment,
    const std::string &read,
    const std::string &ref)
{
    size_t unclippedPos = std::distance(read.begin(),
                                        std::find_if(read.begin(), read.end(),
                                                     boost::bind(&boost::cref<char>, _1) != ' '));

    reference::ReferencePosition fStrandPos(0, fragment.leftClipped() + unclippedPos);

    ISAAC_ASSERT_MSG(fStrandPos.getPosition() < ref.size(), "Reference too short");

    using alignment::Cigar;
    Cigar originalCigar;
    std::string::const_iterator readIterator(read.begin() + unclippedPos);
    std::string::const_iterator refIterator(ref.begin() + unclippedPos);

    unsigned short editDistance = 0;
    if (fragment.leftClipped())
    {
        originalCigar.push_back(Cigar::encode(fragment.leftClipped(), Cigar::SOFT_CLIP));
        for (;
            (read.begin() + fStrandPos.getPosition()) != readIterator && ref.end() != refIterator;
            ++readIterator, ++refIterator)
        {
            editDistance += *refIterator != *readIterator;
        }
    }

    Cigar::Component cigarBit(0, Cigar::ALIGN);
    for (;
        (read.end() - fragment.rightClipped()) != readIterator && ref.end() != refIterator;
        ++readIterator, ++refIterator)
    {
        if ('-' == *readIterator)
        {
            if (Cigar::DELETE == cigarBit.second)
            {
                ++cigarBit.first;
            }
            else
            {
                originalCigar.push_back(Cigar::encode(cigarBit.first, cigarBit.second));
                cigarBit.first = 1;
                cigarBit.second = Cigar::DELETE;
            }
            ISAAC_ASSERT_MSG('*' != *refIterator, "Overlap between insertion and deletion is not allowed");
            ++editDistance;
        }
        else if ('*' == *refIterator)
        {
            if (Cigar::INSERT == cigarBit.second)
            {
                ++cigarBit.first;
            }
            else
            {
                if (cigarBit.first)
                {
                    originalCigar.push_back(Cigar::encode(cigarBit.first, cigarBit.second));
                }
                cigarBit.first = 1;
                cigarBit.second = Cigar::INSERT;
            }
            ++editDistance;
        }
        else
        {
            if (Cigar::ALIGN == cigarBit.second)
            {
                ++cigarBit.first;
            }
            else
            {
                originalCigar.push_back(Cigar::encode(cigarBit.first, cigarBit.second));
                cigarBit.first = 1;
                cigarBit.second = Cigar::ALIGN;
            }
            editDistance += *refIterator != *readIterator;
        }
    }

    if (cigarBit.first)
    {
        originalCigar.push_back(Cigar::encode(cigarBit.first, cigarBit.second));
    }

    if (fragment.rightClipped())
    {
        originalCigar.push_back(Cigar::encode(fragment.rightClipped(), Cigar::SOFT_CLIP));
        for (;
            read.end() != readIterator && ref.end() != refIterator;
            ++readIterator, ++refIterator)
        {
            editDistance += *refIterator != *readIterator;
        }
    }

    TestFragmentAccessor ret(fragment, fStrandPos, read, originalCigar, editDistance);
    return ret;
}

template<typename GapRealignerT>
void addGaps(const std::string &ref, const std::string &gaps, GapRealignerT &realigner)
{
    typename GapRealignerT::GapType currentGap(reference::ReferencePosition(0,0), 0);
    for (std::string::const_iterator gapsIterator(gaps.begin()), refIterator(ref.begin());
        gaps.end() != gapsIterator;
        ++gapsIterator, ++refIterator)
    {
        if (ref.end() == refIterator)
        {
            if (currentGap.length_)
            {
                //if gaps string is longer than reference, wrap around
                typename GapRealignerT::GapType gap(currentGap.pos_ - std::abs(currentGap.length_), currentGap.length_);
                realigner.addGap(0, gap);
            }
            currentGap.pos_ = reference::ReferencePosition(0,0);
            currentGap.length_ = 0;
            refIterator = ref.begin();
        }
        if ('*' == *refIterator)
        {
            ISAAC_ASSERT_MSG(' ' == *gapsIterator, "Spacers in gap string must be in place of fragment insertions");
            continue;
        }

        if ('*' != *gapsIterator)
        {
            if (currentGap.isInsertion())
            {
                typename GapRealignerT::GapType insertion(currentGap.pos_ + currentGap.length_, currentGap.length_);
                realigner.addGap(0, insertion);
                currentGap.length_ = 0;
            }
        }
        else
        {
            --currentGap.length_;
        }

        if('-' != *gapsIterator)
        {
            if (currentGap.isDeletion())
            {
                typename GapRealignerT::GapType deletion(currentGap.pos_ - currentGap.length_, currentGap.length_);
                realigner.addGap(0, deletion);
                currentGap.length_ = 0;
            }
        }
        else
        {
            ++currentGap.length_;
        }

        currentGap.pos_ += 1;
    }

    if (currentGap.length_)
    {
        typename GapRealignerT::GapType gap(currentGap.pos_ - std::abs(currentGap.length_), currentGap.length_);
        realigner.addGap(0, gap);
    }
}

RealignResult realign(
    const std::string &read,
    const std::string &ref,
    const std::string &gaps,
    const isaac::alignment::FragmentMetadata &init,
    const reference::ReferencePosition binStartPos = reference::ReferencePosition(0, 0),
    const reference::ReferencePosition binEndPos = reference::ReferencePosition(reference::ReferencePosition::NoMatch))
{

    TestFragmentAccessor fragment = initFragment(init, read, ref);

    std::vector<std::vector<reference::Contig> > contigList(
        1, std::vector<reference::Contig>(1, reference::Contig(0, "testContig")));
    std::remove_copy_if(ref.begin(), ref.end(), std::back_inserter(contigList.at(0).at(0).forward_),
                        (boost::bind(&boost::cref<char>, _1) == '*' ||
                            boost::bind(&boost::cref<char>, _1) == ' '));

    const build::BarcodeBamMapping barcodeBamMapping(std::vector<unsigned>(1,0),
                                                     std::vector<boost::filesystem::path>(1,""));

    isaac::flowcell::BarcodeMetadataList barcodeMetadataList(1);
    barcodeMetadataList.at(0).setUnknown();
    barcodeMetadataList.at(0).setIndex(0);
    barcodeMetadataList.at(0).setReferenceIndex(0);

    alignment::BinMetadata bin(barcodeMetadataList.size(), 0, reference::ReferencePosition(0,0), 1000000, "tada");

    build::GapRealigner realigner(false, barcodeMetadataList, contigList, barcodeBamMapping);

    addGaps(ref, gaps, realigner);
// this does not work due to * considered    addGaps(ref, ref, realigner, false);
//    addGaps(ref, read, realigner, false);
    realigner.finalizeGaps();

    RealignResult ret(fragment.fStrandPosition_,
                      alignment::Cigar::toString(fragment.cigarBegin(), fragment.cigarEnd()),
                      fragment.editDistance_);

    build::PackedFragmentBuffer dataBuffer;
    alignment::BinMetadata realBin(bin);
    realBin.incrementDataSize(isaac::reference::ReferencePosition(0,0), sizeof(fragment));
    realBin.incrementCigarLength(isaac::reference::ReferencePosition(0,0), 1024, 0);
    realigner.reserve(realBin);
    dataBuffer.resize(realBin);
    std::copy(fragment.begin(), fragment.end(), dataBuffer.begin());

    build::PackedFragmentBuffer::Index index(fragment.fStrandPosition_, 0, 0,
                                             fragment.cigarBegin(), fragment.cigarEnd());

    realigner.realign(binStartPos, binEndPos, index, dataBuffer);

    ret.realignedPos_ = index.pos_;
    ret.realignedCigar_ = alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_);
    const io::FragmentAccessor &updatedFragment = dataBuffer.getFragment(index);
    CPPUNIT_ASSERT_EQUAL(updatedFragment.fStrandPosition_, index.pos_);

    ret.realignedEditDistance_ = updatedFragment.editDistance_;

    return ret;
}

RealignResult realign(
    const std::string &read,
    const std::string &ref,
    const std::string &gaps)
{
    return realign(read, ref, gaps, isaac::alignment::FragmentMetadata());
}


void TestGapRealigner::testAllTogether()
{
    using alignment::Cigar;

    {
        const RealignResult result = realign(
            "      CCCGGAAATTCACACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCG",
            "AAGGAACTCGGCAAACCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGAAAAAA",
            "      *  *                                                                          ");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
            "GACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "  *                                              *              *      ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("2M1I47M1I15M1I7M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
             "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTCAAGAACCACTACTTCCTTAAAAAAAA",
             "                                                                          ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("74M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
             "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGCCATACCATTCTCAAGAACCACTACTTCCTTAAAAAAAA",
             "                                                                 *       ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("65M1I8M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(4, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
             "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "AAAGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTCAAGAACCACTACTTCCTTAAAAAAAA",
             "                                                                             ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("74M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,3), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.realignedEditDistance_));
    }


    {
        const RealignResult result = realign(
             "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "CTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAAAA",
             "***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("3S97M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }

    ISAAC_THREAD_CERR << "----------------------------------------------------------------------------" << std::endl << std::endl << std::endl;

    { // ensuring realignment does not move read into the next bin
        { // first ensure the realignment works for this example when bin boundary is not crossed
            const RealignResult result = realign(
                 "TATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
                 "AAAAAAAAAAAAAAAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAAAA",
                 "----------------",
                 isaac::alignment::FragmentMetadata(),
                 reference::ReferencePosition(0,0), reference::ReferencePosition(0,20));
            CPPUNIT_ASSERT_EQUAL(std::string("84M"), result.originalCigar_);
            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
            CPPUNIT_ASSERT_EQUAL(57, int(result.originalEditDistance_));

            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,16), result.realignedPos_);
            CPPUNIT_ASSERT_EQUAL(std::string("84M"), result.realignedCigar_);
            CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
        }
        { // now check that it does not
            const RealignResult result = realign(
                 "TATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
                 "AAAAAAAAAAAAAAAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAAAA",
                 "----------------",
                 isaac::alignment::FragmentMetadata(),
                 reference::ReferencePosition(0,0), reference::ReferencePosition(0,10));
            CPPUNIT_ASSERT_EQUAL(std::string("84M"), result.originalCigar_);
            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
            CPPUNIT_ASSERT_EQUAL(57, int(result.originalEditDistance_));

            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
            CPPUNIT_ASSERT_EQUAL(std::string("84M"), result.realignedCigar_);
            CPPUNIT_ASSERT_EQUAL(57, int(result.realignedEditDistance_));
        }
    }

    {
        const RealignResult result = realign(
             "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "CTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAAAA",
             "***                                            *              *      ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("3S47M1I15M1I7M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));
    }
/* Seems to be illegal test case. The 4 base insertion beginning 1 base before the read cannot move read by 3 bases!
    {
        const RealignResult result = realign(
             " GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "ACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAAAA",
             "****");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("3S97M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }


    {
        const RealignResult result = realign(
             " GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTACAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "ACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGACAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAAAA",
             "****                                            *              *      ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("3S47M1I15M1I7M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
    }
 */

    // first three bases collapse into insertion, rest moves left
    {
        const RealignResult result = realign(
             " GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "ACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAAAA",
             " ***                                            *              *      ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("3S47M1I15M1I7M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));
    }

    // fifth base collapses into insertion, first four bases move right
    {
        const RealignResult result = realign(
              "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
              "AGACCCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAA",
              "     *                                            *              *      ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("4M1I45M1I15M1I7M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.realignedEditDistance_));
    }

    // first 5 bases collapse into insertion
    {
        const RealignResult result = realign(
              "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
              "ATTGACAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAA",
              "     *****                                        *              *      ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(45, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("5S45M1I15M1I7M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,5), result.realignedPos_);
        // soft clip does not count as edit distance
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
            "   GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCGTTCTTCTCAAGAACCACTACTTC",
            "AGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCACA",
            "                                                                                                   --");
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(78, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
    }


    {
        const RealignResult result = realign(
            "   GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCGTTCTTCTCAAGAACCACTACTTC",
            "AAAGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCACA",
            "                                                                                                   --  ");
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,3), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTT",
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
            "                                                                                                ---");
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("96M3D1M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
            "   GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTC",
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
            "                                                                                                 ---");
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,3), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(67, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
            "GACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTT",
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
            "  -                                               -               -       ---");
        CPPUNIT_ASSERT_EQUAL(std::string("94M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(72, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("2M1D47M1D15M1D7M3D23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.realignedEditDistance_));
    }


    {
        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTC",
            "AAAGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
            "---");
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(67, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,3), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }


    {
        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTC",
            "AAAAGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
            "----                                                                                                   "
            " *");
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(74, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
            " GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCT-G",
            "AGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTTGG",
            "                                                                                                    --  "
            "                                                                                                     -  "
            "                                                                                                    -   "
            );
        CPPUNIT_ASSERT_EQUAL(std::string("99M1D1M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.originalEditDistance_));
        CPPUNIT_ASSERT(std::string("99M2D1M") == result.realignedCigar_||
                       std::string("99M1D1D1M") == result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
             "      GA--TCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "AAAGAAAACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
             "    --- --");
        CPPUNIT_ASSERT_EQUAL(std::string("2M2D96M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("1M3D1M2D96M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,3), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));

    }

    {
        const RealignResult result = realign(
             "      GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "AAAGAAGACCTCAATCAGGCAATATGAAGTT*CAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCAT*CTTCTCAAGAACCACTACTTCCTTAAAAAA",
             "                                *                                                 *");
        CPPUNIT_ASSERT_EQUAL(std::string("25M1I49M1I24M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("25M1I49M1I24M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.realignedEditDistance_));

    }

    {
        const RealignResult result = realign(
             "      GA--TCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "AAAGAAAACCTCAATCAGGCAATATGAAGTT*CAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCAT*CTTCTCAAGAACCACTACTTCCTTAAAAAA",
             "    ---                      ** **                                                                              "
             "        --                      *                                                 *                             "
             );
        CPPUNIT_ASSERT_EQUAL(std::string("2M2D21M1I49M1I24M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("1M3D1M2D21M1I49M1I24M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,3), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(7, int(result.realignedEditDistance_));

    }


    {
        // the shortest is expected to win
        const RealignResult result = realign(
             "  GTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "GGATCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
             " --                                                                                                      "
             "  -                                                                                                      "
             " --                                                                                                      "
         );
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,2), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("1M1D96M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));

    }

    {
        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTAAAAA  ",
            "                                                                                               *******");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("95M5S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTC-TTT",
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCAATT",
            "                                                                                                  ***"
            "                                                                                                 -   "
            );
        CPPUNIT_ASSERT_EQUAL(std::string("97M1D3M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(std::string("97M3S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }


/*  Currently disabled because of SAAC-253 (see test for it below)
    {
        const RealignResult result = realign(
            "GACCCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCATT",
            "AGAATCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCATT",
            "    ****                                                                                            "
            "    ---                                                                                             "
            "       ------                                                                                       "
            );

        CPPUNIT_ASSERT_EQUAL(std::string("91M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(64, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("4S87M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,13), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }
*/

    // overlapping gaps are not merged anymore
    {
        const RealignResult result = realign(
            " GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCG-G",
            "AGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCGTTGGGGGGG",
            "                                                                                                   ----      "
            "                                                                                                    -        "
            );
        CPPUNIT_ASSERT_EQUAL(std::string("99M1D1M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("98M4D2M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(4, int(result.realignedEditDistance_));
    }

    {// test for cigar compacting     * SAAC-247 GATK 'Adjacent I/D events in read' error
        const RealignResult result = realign(
              "   GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
              "AAAGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
              "                                                                  **                                      "
              "                                                                  *                                       "
              );
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(23, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("63M3I34M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,3), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.realignedEditDistance_));
    }


    { //SAAC-251    Gap realigner moves perfectly aligning read by one base
        const RealignResult result = realign(
              "    GAATCATCGAATGGACTCGAATGGAATAATCCTTGAACGGAATCGATTGGAATCATCATCGGATGGATACGANTGGAATCATCATTGANTGGAATCGAAT",
              "AATGGAATCATCGAATGGACTCGAATGGAATAATCCTTGAACGGAATCGATTGGAATCATCATCGGATGGATACGAATGGAATCATCATTGAATGGAATCGAATGGAA",
              "   -                                                                                                        "
              );
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
    }

    {// SAAC-253   Gap realigner produces a pair of insertion and deletion of an equal size in place of a mismatch
        const RealignResult result = realign(
              "    GAATCATCGAATGGACTCGAATGGAATAATCCTTGAACGGAATCGATTGGAACCATCATCGGATGGATACGAATGGAATCATCATTGAATGGAATCGAAT",
              "AATGGAATCATCGAATGGACTCGAATGGAATAATCCTTGAACGGAATCGATTGGAATCATCATCGGATGGATACGAATGGAATCATCATTGAATGGAATCGAATGGAA",
              "                                                        -                                                   "
              "                                                        *"
              );
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));
    }

    {// Test for preservation of alignment-independent clipping
        isaac::alignment::FragmentMetadata fragmentMetadata;
        fragmentMetadata.lowClipped = 6;
        fragmentMetadata.highClipped = 24;

        const RealignResult result = realign(
            " GACTCAAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTT",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "        *",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S65M24S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,7), result.originalPos_);

        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,7), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("6S1M1I63M24S"), result.realignedCigar_);
    }

    {// Test for preservation of alignment-independent clipping
        isaac::alignment::FragmentMetadata fragmentMetadata;
        fragmentMetadata.lowClipped = 6;
        fragmentMetadata.highClipped = 24;

        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGCCATACCATTCTCAAGAACCACTACTTCCTT",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "   *                                                                  ",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S65M24S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(4, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("6S65M24S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(4, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);
    }

    {// Test for preservation of alignment-independent clipping (right side clipping prevents the introduction of the insertion at the end)
        isaac::alignment::FragmentMetadata fragmentMetadata;
        fragmentMetadata.lowClipped = 6;
        fragmentMetadata.highClipped = 24;

        const RealignResult result = realign(
        //    SSSSSS                                                                     SSSSSSSSSSSSSSSSSSSSSSSS
            " GACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
//            ||||||||||||||||||||||||||||||||||||||||||||||||||xxx|xx|xxxx|x|xxx|xxxxxx||xxxxxxxxxxx|xxxxx|xxxxx
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "                                                  *              *      ***",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S69M24S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,7), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(40, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,7), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("6S43M1I15M1I9M24S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(15, int(result.realignedEditDistance_));
    }

    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
        isaac::alignment::FragmentMetadata fragmentMetadata;
        fragmentMetadata.lowClipped = 6;
        fragmentMetadata.highClipped = 20;

        const RealignResult result = realign(
        //   SSSSSS                                                                          SSSSSSSSSSSSSSSSSSSS
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
        //   xxx||||||||||||||||||||||||||||||||||||||||||||||||xxx|xx|xxxx|x|xxx|xxxxxx||xxxxxxxxxxx|xxxxx|xxxxx
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "   *                                              *              *      ***",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S74M20S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(43, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("6S44M1I15M1I7M3I3M20S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(8, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);
    }

}

void TestGapRealigner::testFull()
{
    testAllTogether();

    // this is too complex for fast realigner
    {
        const RealignResult result = realign(
            "      CCCGGAAATTCACACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCG",
            "AAGGAACTCGGCAAACCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGAAAAAA",
            "      *  *     *   *                                                                ");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(std::string("1M1I6M1I91M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,8), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(7, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
              "GACCACAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
              "ATTGTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAA",
              "    ****");

        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("4S96M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.realignedPos_);
        // soft clip does not count as edit distance but first unclipped base is mismatching
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));
    }

    {
        const RealignResult result = realign(
              "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
              "ATTGACAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAA",
              "    ****                                          *              *      ***");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(45, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("4S46M1I15M1I7M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.realignedPos_);
        // soft clip does not count as edit distance but first unclipped base is mismatching
        CPPUNIT_ASSERT_EQUAL(6, int(result.realignedEditDistance_));
    }

}
