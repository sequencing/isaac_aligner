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
 ** \file testGapRealigner.cpp
 **
 ** Test cases for gap realignment.
 **
 ** \author Roman Petrovski
 **/

#include <string>

#include "build/gapRealigner/OverlappingGapsFilter.hh"
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
        const io::FragmentHeader &fragment,
        const reference::ReferencePosition fStrandPosition,
        const std::string &read,
        const alignment::Cigar &cigar,
        const unsigned observedLength,
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
                *b++ = (0x03 & oligo::getValue(base)) | ('N' == base ? 0 : 0x20);
                ++readLength_;
            }
        }
        std::copy(cigar.begin(), cigar.end(), const_cast<unsigned*>(cigarBegin()));
        cigarLength_ = cigar.size();
        observedLength_ = observedLength;
        editDistance_ = editDistance;
        mateFStrandPosition_ = fStrandPosition;
        alignmentScore_ = 1;
        lowClipped_ = fragment.lowClipped_;
        highClipped_ = fragment.highClipped_;
        clusterId_ = fragment.clusterId_;
    }

    const char *begin() const {return reinterpret_cast<const char *>(this);}
    const char *end() const {return reinterpret_cast<const char *>(cigarEnd());}
};


struct RealignResult
{
    RealignResult(
        const reference::ReferencePosition fStrandPosition,
        const std::string &originalCigar, const unsigned short editDistance,
        const build::gapRealigner::GapsRange range):
            originalPos_(fStrandPosition),
            originalCigar_(originalCigar),
            originalEditDistance_(editDistance),
            overlappingGapsFilter_(range){}
    reference::ReferencePosition originalPos_;
    reference::ReferencePosition realignedPos_;
    std::string originalCigar_;
    std::string realignedCigar_;
    unsigned short originalEditDistance_;
    unsigned short realignedEditDistance_;
    build::gapRealigner::OverlappingGapsFilter overlappingGapsFilter_;
};

TestFragmentAccessor initFragment(
    const io::FragmentHeader &fragment,
    const std::string &read,
    const std::string &ref)
{
    size_t unclippedPos = std::distance(read.begin(),
                                        std::find_if(read.begin(), read.end(),
                                                     boost::bind(&boost::cref<char>, _1) != ' '));

    size_t referenceLeftOverhang = std::distance(ref.begin(),
                                        std::find_if(ref.begin(), ref.end(),
                                                     boost::bind(&boost::cref<char>, _1) != ' '));

    reference::ReferencePosition fStrandPos(0, fragment.leftClipped() + unclippedPos);

    ISAAC_ASSERT_MSG(fStrandPos.getPosition() < ref.size(), "Reference too short");

    using alignment::Cigar;
    Cigar originalCigar;
    std::string::const_iterator readIterator(read.begin() + unclippedPos);
    std::string::const_iterator refIterator(ref.begin() + unclippedPos);

    unsigned short editDistance = 0;
    if (fragment.leftClipped() + referenceLeftOverhang)
    {
        originalCigar.push_back(Cigar::encode(fragment.leftClipped() + referenceLeftOverhang, Cigar::SOFT_CLIP));
        while((read.begin() + fStrandPos.getPosition() + referenceLeftOverhang) != readIterator && ref.end() != refIterator)
        {
            ++readIterator, ++refIterator;
        }
    }

    unsigned observedLength = 0;
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
            ++observedLength;
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
            ++observedLength;
        }
    }

    if (cigarBit.first)
    {
        originalCigar.push_back(Cigar::encode(cigarBit.first, cigarBit.second));
    }

    if (ref.end() == refIterator)
    {
        if (read.end() != readIterator)
        {
            originalCigar.push_back(Cigar::encode(std::distance(readIterator, read.end()), Cigar::SOFT_CLIP));
        }
    }
    else if (fragment.rightClipped())
    {
        originalCigar.push_back(Cigar::encode(fragment.rightClipped(), Cigar::SOFT_CLIP));
        while(read.end() != readIterator && ref.end() != refIterator)
        {
            ++readIterator, ++refIterator;
        }
    }

    TestFragmentAccessor ret(fragment, fStrandPos, read, originalCigar, observedLength, editDistance);
    return ret;
}

template<typename RealignerGapsT>
void addGaps(const std::string &ref, const std::string &gaps, RealignerGapsT &realigner)
{
    typename RealignerGapsT::GapType currentGap(reference::ReferencePosition(0,0), 0);
    for (std::string::const_iterator gapsIterator(gaps.begin()), refIterator(ref.begin());
        gaps.end() != gapsIterator;
        ++gapsIterator, ++refIterator)
    {
        if (ref.end() == refIterator)
        {
            if (currentGap.length_)
            {
                //if gaps string is longer than reference, wrap around
                typename RealignerGapsT::GapType gap(currentGap.pos_ - std::abs(currentGap.length_), currentGap.length_);
                realigner.addGap(gap);
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
                typename RealignerGapsT::GapType insertion(currentGap.pos_ + currentGap.length_, currentGap.length_);
                realigner.addGap(insertion);
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
                typename RealignerGapsT::GapType deletion(currentGap.pos_ - currentGap.length_, currentGap.length_);
                realigner.addGap(deletion);
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
        typename RealignerGapsT::GapType gap(currentGap.pos_ - std::abs(currentGap.length_), currentGap.length_);
        realigner.addGap(gap);
    }
}

RealignResult realign(
    const unsigned mismatchCost,
    const unsigned gapOpenCost,
    const std::string &read,
    const std::string &ref,
    const std::string &gaps,
    const io::FragmentHeader &init,
    const reference::ReferencePosition binStartPos = reference::ReferencePosition(0, 0),
    reference::ReferencePosition binEndPos = reference::ReferencePosition(reference::ReferencePosition::NoMatch))
{

//    ISAAC_THREAD_CERR << "Initialized " <<
    TestFragmentAccessor fragment = initFragment(init, read, ref);
//        oligo::bclToString(fragment.basesBegin(), fragment.basesEnd() - fragment.basesBegin()) << " " << fragment << std::endl;

    std::vector<std::vector<reference::Contig> > contigList(
        1, std::vector<reference::Contig>(1, reference::Contig(0, "testContig")));
    std::remove_copy_if(ref.begin(), ref.end(), std::back_inserter(contigList.at(0).at(0).forward_),
                        (boost::bind(&boost::cref<char>, _1) == '*' ||
                            boost::bind(&boost::cref<char>, _1) == ' '));

    if (binEndPos.isNoMatch())
    {
        binEndPos = binStartPos + contigList.at(0).at(0).forward_.size();
        ISAAC_THREAD_CERR << "binEndPos:" << binEndPos << std::endl;
    }


    isaac::flowcell::BarcodeMetadataList barcodeMetadataList(1);
    barcodeMetadataList.at(0).setUnknown();
    barcodeMetadataList.at(0).setIndex(0);
    barcodeMetadataList.at(0).setReferenceIndex(0);

    alignment::BinMetadata bin(barcodeMetadataList.size(), 0, reference::ReferencePosition(0,0), 1000000, "tada", 0);

    std::vector<alignment::TemplateLengthStatistics> templateLengthStatistics(1);
    build::GapRealigner realigner(true, false, 8, mismatchCost, gapOpenCost, 0, false, barcodeMetadataList, templateLengthStatistics, contigList);
    build::RealignerGaps realignerGaps;
    addGaps(ref, gaps, realignerGaps);
// this does not work due to * considered    addGaps(ref, ref, realigner, false);
//    addGaps(ref, read, realigner, false);
    realignerGaps.finalizeGaps();

    static build::gapRealigner::Gaps foundGaps;
    foundGaps.reserve(100000);
    RealignResult ret(fragment.fStrandPosition_,
                      alignment::Cigar::toString(fragment.cigarBegin(), fragment.cigarEnd()),
                      fragment.editDistance_,
                      realignerGaps.findGaps(0, binStartPos, binStartPos, binEndPos, foundGaps));

    build::PackedFragmentBuffer dataBuffer;
    alignment::BinMetadata realBin(bin);
    realBin.incrementDataSize(isaac::reference::ReferencePosition(0,0), sizeof(fragment));
    realBin.incrementCigarLength(isaac::reference::ReferencePosition(0,0), 1024, 0);
    realigner.reserve(realBin);
    dataBuffer.resize(realBin);
    std::copy(fragment.begin(), fragment.end(), dataBuffer.begin());

    build::PackedFragmentBuffer::Index index(fragment.fStrandPosition_, 0, 0,
                                             fragment.cigarBegin(), fragment.cigarEnd());

    realigner.realign(realignerGaps, binStartPos, binEndPos, index, dataBuffer.getFragment(index), dataBuffer);

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
    const std::string &gaps,
    const io::FragmentHeader &init,
    const reference::ReferencePosition binStartPos = reference::ReferencePosition(0, 0),
    const reference::ReferencePosition binEndPos = reference::ReferencePosition(reference::ReferencePosition::NoMatch))
{
    return realign(1, 0, read, ref, gaps, init, binStartPos, binEndPos);
}

RealignResult realign(
    const unsigned mismatchCost,
    const unsigned gapOpenCost,
    const std::string &read,
    const std::string &ref,
    const std::string &gaps)
{
    return realign(mismatchCost, gapOpenCost, read, ref, gaps, io::FragmentHeader());
}

RealignResult realign(
    const std::string &read,
    const std::string &ref,
    const std::string &gaps)
{
    return realign(1, 0, read, ref, gaps, io::FragmentHeader());
}


void TestGapRealigner::testFull()
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    { // ensuring realignment does not move read into the next bin
        { // first ensure the realignment works for this example when bin boundary is not crossed
            const RealignResult result = realign(
                 "TATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
                 "AAAAAAAAAAAAAAAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAAAA",
                 "----------------",
                 io::FragmentHeader(),
                 reference::ReferencePosition(0,0), reference::ReferencePosition(0,20));
            CPPUNIT_ASSERT_EQUAL(std::string("84M"), result.originalCigar_);
            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
            CPPUNIT_ASSERT_EQUAL(57, int(result.originalEditDistance_));

            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,16), result.realignedPos_);
            CPPUNIT_ASSERT_EQUAL(std::string("84M"), result.realignedCigar_);
            CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));

            CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
        }
        { // now check that it does not
            const RealignResult result = realign(
                 "TATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
                 "AAAAAAAAAAAAAAAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAAAA",
                 "----------------",
                 io::FragmentHeader(),
                 reference::ReferencePosition(0,0), reference::ReferencePosition(0,10));
            CPPUNIT_ASSERT_EQUAL(std::string("84M"), result.originalCigar_);
            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
            CPPUNIT_ASSERT_EQUAL(57, int(result.originalEditDistance_));

            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
            CPPUNIT_ASSERT_EQUAL(std::string("84M"), result.realignedCigar_);
            CPPUNIT_ASSERT_EQUAL(57, int(result.realignedEditDistance_));

            CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
        }
    }

    { // ensuring realignment does not move read into the next bin
        { // first ensure the realignment works for this example when bin boundary is not crossed
            io::FragmentHeader fragment;
            fragment.lowClipped_ = 16;

            const RealignResult result = realign(
                 "TATGAAGTTGCAGGAACTGGAAGAGGAGAGAT",
                 "AAAAAAAAAAAAAAAATATGAAGTTGCAGGAAC",
                 "----------------",
                 fragment);
            CPPUNIT_ASSERT_EQUAL(std::string("16S16M"), result.originalCigar_);
            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,16), result.originalPos_);
            CPPUNIT_ASSERT_EQUAL(9, int(result.originalEditDistance_));

            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,32), result.realignedPos_);
            CPPUNIT_ASSERT_EQUAL(std::string("16S16M"), result.realignedCigar_);
            CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));

            CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
        }
        { // now check that it does not
            io::FragmentHeader fragment;
            fragment.lowClipped_ = 16;

            const RealignResult result = realign(
                 "TATGAAGTTGCAGGAACTGGAAGAGGAGAGAT",
                 "AAAAAAAAAAAAAAAATATGAAGTTGCAGGAA",
                 "----------------",
                 fragment);
            CPPUNIT_ASSERT_EQUAL(std::string("16S16M"), result.originalCigar_);
            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,16), result.originalPos_);
            CPPUNIT_ASSERT_EQUAL(9, int(result.originalEditDistance_));

            CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,16), result.realignedPos_);
            CPPUNIT_ASSERT_EQUAL(std::string("16S16M"), result.realignedCigar_);
            CPPUNIT_ASSERT_EQUAL(9, int(result.realignedEditDistance_));

            CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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
        CPPUNIT_ASSERT_EQUAL(75, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("3S47M1I15M1I7M3I23M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {
        const RealignResult result = realign(3,4,
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTTT",
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTT",
            "                                                                                                ---");
        CPPUNIT_ASSERT_EQUAL(std::string("98M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("96M3D2M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }


    {
        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTC",
            "AAAAGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
            "----                                                                                                          "
            " *");
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(74, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("97M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(1U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(3U, result.overlappingGapsFilter_.overlap(0));
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

        CPPUNIT_ASSERT_EQUAL(2U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(3U, result.overlappingGapsFilter_.overlap(0));
        CPPUNIT_ASSERT_EQUAL(6U, result.overlappingGapsFilter_.overlap(1));
    }

    {
        const RealignResult result = realign(3,4,
            " GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCT-GT",
            "AGACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTTGT",
            "                                                                                                    --  "
            "                                                                                                     -  "
            "                                                                                                    -   "
            );
        CPPUNIT_ASSERT_EQUAL(std::string("99M1D2M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("99M2D2M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(2U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(3U, result.overlappingGapsFilter_.overlap(0));
        CPPUNIT_ASSERT_EQUAL(6U, result.overlappingGapsFilter_.overlap(1));
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {
        const RealignResult result = realign(3,4,
             "    GGGA--TCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "AGGGAAAACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
             "    --- --");
        CPPUNIT_ASSERT_EQUAL(std::string("4M2D96M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("3M3D1M2D96M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(1U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(3U, result.overlappingGapsFilter_.overlap(0));

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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(1U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(3U, result.overlappingGapsFilter_.overlap(0));
    }

/* disabled due to SAAC-381 Gap realigner fails when read is realigned against multiple insertion gaps starting at the same position
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
*/
    testMore();
}


void TestGapRealigner::testMore()
{
    using alignment::Cigar;

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
        CPPUNIT_ASSERT_EQUAL(2, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    { // N-containing reads get realigned and their edit distance broken for no reason
      // Turned out to be a bug in TestFragmentAccessor fragment initialization
        const RealignResult result = realign(
              "    GAATCATCGAATGGACTCGAATGGAATAATCCTTGAACGGAATCGATTGGAATCATCATCGGATGGATACGANTGGAATCATCATTGANTGGAATCGAAT",
              "AATGGAATCATCGAATGGACTCGAATGGAATAATCCTTGAACGGAATCGATTGGAATCATCATCGGATGGATACGAATGGAATCATCATTGAATGGAATCGAATGGAA",
              "    -                                                                                                       "
              );
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.originalEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,4), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(1U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(3U, result.overlappingGapsFilter_.overlap(0));
    }

    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 1;
        fragmentMetadata.highClipped_ = 0;

        const RealignResult result = realign(
           //S
            "TTGCTCAGGAAGCTTCCTTCAAAATGTCTACTG",
            "TTGACGGCCCAGCTTCCTTCAAAATGTCTACTG",
            "            ************",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("1S32M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,1), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(7, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("12S21M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,12), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }


    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 128;
        fragmentMetadata.highClipped_ = 0;

        const RealignResult result = realign(
         //  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
            "GTCTTGTCTTGTTTTTGGTTTTTTTTGATTTTTCCTAGTTGTATTTTTTTTTTTTTTTTTTTTTTTTTTTAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGGTTGAACACATTGGCCTCAGGAAGCTTCCTTCAAAATGTCTACTGTTCACGAAATCCTGTGCAAGCTCAGCTTGGAGGGTGATCACTCTACACCCCCAAGTGCATATGGGTCTGTCAAAG",
            "GTCTTGTCTTGTTTTTGGTTTTTTTTGATTTTTCCTAGTTGTATTTTTTTTTTTTTTTTTTTTTTTTTTTAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGGTTGAACACATTGGCACGGCCCAGCTTCCTTCAAAATGTCTACTGTTCACGAAATCCTGTGCAAGCTCAGCTTGGAGGGTGATCACTCTACACCCCCAAGTGCATATGGGTCTGTCAAAG",
            "                                                                                                                                                          *********************************************************************************",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("128S122M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,128), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(7, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("154S96M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,154), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    // notice that this one is not supposed to work unless the insertions are allowed to eat bases beyond the sequence start
//    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
//
//        isaac::alignment::FragmentMetadata fragmentMetadata;
//        fragmentMetadata.lowClipped = 0;
//        fragmentMetadata.highClipped = 0;//128;
//
//        const RealignResult result = realign(
//         //                                                                                                                            SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
//            "TGGTTGAACACATTGGCCTCAGGAAGCTTCCTTCAAAATGTCTACTGTTCACGAAATCCTGTGCAAGCTCAGCTTGGAGGGTGATCACTCTACACCCCCAAGTGCATATGGGTCTGTCAAAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTTACCAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAACAAAAACAAAGATAAAACCTAAAAAAATAATACAAAAT",
//            "TGGTTGAACACATTGGCACGGCCCAGCTTCCTTCAAAATGTCTACTGTTCACGAAATCCTGTGCAAGCTCAGCTTGGAGGGTGATCACTCTACACCCCCAAGTGCATATGGGTCTGTCAAAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTTACCAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAACAAAAACAAAGATAAAACCTAAAAAAATAATACAAAAT",
//            "                          *********************************************************************************",
//            fragmentMetadata);
////        CPPUNIT_ASSERT_EQUAL(std::string("122M128S"), result.originalCigar_);
//        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
//        CPPUNIT_ASSERT_EQUAL(7, int(result.originalEditDistance_));
//
//        CPPUNIT_ASSERT_EQUAL(std::string("81S41M128S"), result.realignedCigar_);
//        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
//        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,26), result.realignedPos_);
//
//        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
//    }


    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 56;
        fragmentMetadata.highClipped_ = 0;

        const RealignResult result = realign(
        //   SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS                                         SSSSSSSS
            "AGACTCAATCAGGCAATATGAAGTTGCAGGGATTACACATAGACAGAT"
                                                                                                     "ACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGGATTACACATAGACAGATCACAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "                                                -----------------------------------------                                 ",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("56S25M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,56), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(18, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("56S25M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,97), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 56;
        fragmentMetadata.highClipped_ = 6;

        const RealignResult result = realign(
        //   SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS                                         SSSSSSSS                   SSSSSS
            "AGACTCAATCAGGCAATATGAAGTTGCAGGGATTACACATAGACAGAT"
                                                                                                     "ACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGGATTACACATAGACAGATCACAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "                                                -----------------------------------------                                 ",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("56S19M6S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,56), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(13, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("56S19M6S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,97), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 0;
        fragmentMetadata.highClipped_ = 56;

        const RealignResult result = realign(
        //                            SSSSSSSSSSSSSSSSSSSSSSS                                         SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
            "                                         "
            "AGACTCAATCAGGCAATATGAAGTTGCAGGGATTACACATAGACAGAT"
                                                                                                     "ACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGGATTACACATAGACAGATCACAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "                                                -----------------------------------------                                 ",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("25M56S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,41), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(22, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("25M56S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 6;
        fragmentMetadata.highClipped_ = 56;

        const RealignResult result = realign(
        //   SSSSSS                   SSSSSSSSSSSSSSSSSSSSSSS                                         SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
            "                                         "
            "AGACTCAATCAGGCAATATGAAGTTGCAGGGATTACACATAGACAGAT"
                                                                                                     "ACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGGATTACACATAGACAGATCACAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "                                                -----------------------------------------                                 ",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S19M56S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,47), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(16, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("6S19M56S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for preservation of alignment-independent clipping
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 6;
        fragmentMetadata.highClipped_ = 24;

        const RealignResult result = realign(
            " GACTCAAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTT",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "        *",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S65M24S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,7), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(50, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,7), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("6S1M1I63M24S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Testing that alignment-independent clipping does not prevent realignment (SAAC-446)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 6;
        fragmentMetadata.highClipped_ = 0;

        const RealignResult result = realign(1,2,
            "                    TCTGAG"
                                      "GGTCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "AGACTCGGTCAGGAAAAAAAAAAAAAAAAAAAACAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "             --------------------",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S95M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,26), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(20, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("6S7M20D88M"), result.realignedCigar_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {   // Test for correct compacting of XXXIYYS CIGAR
        // (make sure realignment does not occur as in this case the entire reads gets soft-clipped away)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 0;
        fragmentMetadata.highClipped_ = 8;

        const RealignResult result = realign(
            "GACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTT",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "**************************************************************************************",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("86M8S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(67, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(67, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(std::string("86M8S"), result.realignedCigar_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for preservation of alignment-independent clipping
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 6;
        fragmentMetadata.highClipped_ = 24;

        const RealignResult result = realign(
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGCCATACCATTCTCAAGAACCACTACTTCCTT",
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "   *                                                                  ",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S65M24S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("6S65M24S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for preservation of alignment-independent clipping (right side clipping prevents the introduction of the insertion at the end)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 6;
        fragmentMetadata.highClipped_ = 24;

        const RealignResult result = realign(
        //    SSSSSS                                                                     SSSSSSSSSSSSSSSSSSSSSSSS
            " GACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
//            ||||||||||||||||||||||||||||||||||||||||||||||||||xxx|xx|xxxx|x|xxx|xxxxxx||xxxxxxxxxxx|xxxxx|xxxxx
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "                                                  *              *      ***",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S69M24S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,7), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(19, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,7), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(std::string("6S43M1I15M1I9M24S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }


    {// Test for preservation of alignment-independent clipping (right-side clipping begins after the three-base insertion introduced)
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 6;
        fragmentMetadata.highClipped_ = 20;

        const RealignResult result = realign(
        //   SSSSSS                                                                          SSSSSSSSSSSSSSSSSSSS
            "GACCTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
        //   xxx||||||||||||||||||||||||||||||||||||||||||||||||xxx|xx|xxxx|x|xxx|xxxxxx||xxxxxxxxxxx|xxxxx|xxxxx
            "AGACTCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTCAGGCTTATCTTGGCATACCATTCTCAAGAACCACTACTTCCTTAAAAAA",
            "   *                                              *              *      ***",
            fragmentMetadata);
        CPPUNIT_ASSERT_EQUAL(std::string("6S74M20S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(22, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("6S44M1I15M1I7M3I3M20S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {
        const RealignResult result = realign(
            "      CCCGGAAATTCACACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCG",
            "AAGGAACTCCGGAAACCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGAAAAAA",
            "      *  *     *   *                                                                                            "
            );
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(7, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("1M1I6M1I91M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,8), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(5, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
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

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for preservation of original alignment if the realignment yield equivalent mismatches edit distance
        const RealignResult result = realign(
            "CAGAAACCAGATTTTTATTCG-TTGATGAAAGTCCTTGCAGTTTTTCCCATGGTCTATTTGGAGAACCACTACATACTAGAAAGCTAGTATGACAAAATTT",
        //   |||||||||||||||||||||x|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            "CAGAAACCAGATTTTTATTCGCTAGATGAAAGTCCTTGCAGTTTTTCCCATGGTCTATTTGGAGAACCACTACATACTAGAAAGCTAGTATGACAAAATTTT",
            "                       -");
        CPPUNIT_ASSERT_EQUAL(std::string("21M1D79M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("21M1D79M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(2, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for preservation of original alignment if the realignment requires to introduce two separate gaps to get rid of one mismatch
        const RealignResult result = realign(3, 4,
            "GGCTGTTGGTTTTTTTTTTCCTTCTGTAAAATAAAACTACCTCCTCCTACCTGTTTCTGCTTCCAAGAAGCTTAAGTGAACATGATTTCCAAACTGTCTT",
        //   ||||||||x|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            "GGCTGTTGTTTTTTTTTTTCCTTCTGTAAAATAAAACTACCTCCTCCTACCTGTTTCTGCTTCCAAGAAGCTTAAGTGAACATGATTTCCAAACTGTCTT",
            "       *                                                                                            "
            "        -");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(1, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
}

    {// Test for introduction of long deletion
        const RealignResult result = realign(3, 4,
            "                                            TAGTAACAGTTGGTGGGCCGGCACTGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAAGGCAGGTGG",
                                                    //   ||xxx|xx|xxxxx|x||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            "GGTGCAGACTAGTAACAGTTGGTGGGCCGGCACTGATCCATTAATATATATTGCCCAGGTGCCGGCACTGATCCATTAATATATATTGCCCAGGTGCCGTGGCTCACCTATAATCCCAGCACTTTGAGAGGCCAAGGCAGGTGGATTGCTTGAGCCCAGGAG",
            "                         -----------------------------------");
        CPPUNIT_ASSERT_EQUAL(std::string("100M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,44), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(11, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("16M35D84M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(35, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,9), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for what was causing assertion failure in GapRealigner::findStartPos
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.flags_.reverse_ = true;

        const RealignResult result = realign(3, 4,
            "A-----C",
            "AGATCAG",
            "   **");
        CPPUNIT_ASSERT_EQUAL(std::string("1M5D1M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("1M5D1M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);

        // unit tests don't throw sequence gaps into the mix
        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    {// Test for what was causing assertion failure in GapRealigner::findStartPos
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.flags_.reverse_ = true;

        const RealignResult result = realign(3, 4,
            "A-----C",
            "AGATCAG",
            " ----- "
            "   **");
        CPPUNIT_ASSERT_EQUAL(std::string("1M5D1M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("1M5D1M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(6, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(1U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(0x3U, result.overlappingGapsFilter_.overlap(0));
    }

    {   // SAAC-381 Gap realigner fails when read is realigned against multiple insertion gaps starting at the same position
        const RealignResult result = realign(3, 4,
            "GGCTGTTGGTTTTTTTTTTCCTTCTGTAAAATAAAACTACC",
        //   xxxxxxxxxxxx|||||||||||||||||||||||||||||
            "AAAAAAAAAAAATTTTTTTCCTTCTGTAAAATAAAACTACC",
            "            ***********************      "
            "            **********************       "
            );
        CPPUNIT_ASSERT_EQUAL(std::string("41M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(12, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("41M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(12, int(result.realignedEditDistance_));
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,0), result.realignedPos_);

        CPPUNIT_ASSERT_EQUAL(1U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(0x3U, result.overlappingGapsFilter_.overlap(0));
    }

    {
        // the multiple overlapping deletions
        const RealignResult result = realign(
             "  ATCAATCAAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTT",
             "GGATCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
             "          ---                                                                                            "
             "         ----                                                                                            "
             "         -----                                                                                           "
             "         ------                                                                                          "
             "         -------                                                                                         "
             "         --------                                                                                        "
             "         ---------                                                                                       "
             "         ----------                                                                                      "
             "         -----------                                                                                     "
             "         ------------                                                                                    "
             "         -------------                                                                                   "
             "         --------------                                                                                  "
             "                        --------------                                                                   "
         );
        CPPUNIT_ASSERT_EQUAL(std::string("94M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,2), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(58, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("8M3D86M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,2), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(1U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(0xfffU, result.overlappingGapsFilter_.overlap(0));
    }

    { // check to prevent cases where long insertion reduces mismatches but misplaces the entire read
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 16;
        fragmentMetadata.highClipped_ = 0;

        const RealignResult result = realign(
             "TTGCATTCGTCCTGGCATGAAGTACGTCTGGCTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
             "GGATCAATCAGGCAATATGAAGTTGCAGGAACTGGAAGAGGAGAGATAGTTCAGGCTTATCTTGGCCATACCATTCTTCTCAAGAACCACTACTTCCTTAAAAAA",
             "                ********************************************************************************         ",
             fragmentMetadata
         );
        CPPUNIT_ASSERT_EQUAL(std::string("16S89M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,16), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(8, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("16S89M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,16), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(8, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

    { // *SAAC-514 gap realigner does not soft-clip the end of the read hanging outside the chromosome
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 6;
        fragmentMetadata.highClipped_ = 0;
        const RealignResult result = realign(
             "NNNNNNACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGATCACAGGTCTA",
             "CCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATG",
             "                                                                        -                                       "
             "                                                                          *****",
             fragmentMetadata
         );
        CPPUNIT_ASSERT_EQUAL(std::string("6S106M13S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(39, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("6S66M1D1M5I38M9S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,6), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(41, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

/*
    {
        // the multiple overlapping deletions
        const RealignResult result = realign(
            // offset spaces
            "                                                                                                         "
            "                                                                                                         "
            "                                                                                                         "
            "                                                                      "
            // original alignment position
            "ACTCGTTCCTCTACCCTCCCCTACCTCCGTTTTCTTTTTTTTTTTAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCCT"
            "GTCTCTCCCTTCTGCCCACCCTCAGAACGGGCGCTCGATCAGCAAGCGGGGCGTGAGCTGTCAGTTTGGGCCTGACGTCACCAAGGCCTTCTTGGAAGAGAACAAC"
            "CTGGACTATATCATCCGCAGCCACGAAGTCAAGGCCGG",
            // reference
            "TCCAAGCTGAGCACGCTCGTGGAAACCACACTCAAAGAGACAGAGAAGATTACAGTATGTGGGGACACCCATGGCCAGTTCTATGACCTCCTCAACATATTCGAG"
            "CTCAACGGTTTACCCTCGGAGACCAACCCCTATATATTTAATGGTGACTTTGTGGACCGAGGCTCCTTCTCTGTAGAAGTGATCCTCACCCTTTTCGGCTTCAAG"
            "CTCCTGTACCCAGATCACTTTCACCTCCTTCGAGGCAACCACGAGACAGACAACATGAACCAGATCTACGGTTTCGAGGGTGAGGTGAAGGCCAAGTACACAGCC"
            "CAGATGTACGAGCTCTTTAGCGAGGTGTTCGAGTGGCTCCCGTTGGCCCAGTGCATCAACGGCAAAGTGC"
            // original alignment position
            "TGATCATGCACGGAGGCCTGTTCAGTGAAGACGGTGTCACCCTGGATGACATCCGGAAAATTGAGCGGAATCGACAACCCCCAGATTCAGGGCCCATGTGTGACC"
            "TGCTCTGGTCAGATCCACAGCCACAGAACGGGCGCTCGATCAGCAAGCGGGGCGTGAGCTGTCAGTTTGGGCCTGACGTCACCAAGGCCTTCTTGGAAGAGAACA"
            "ACCTGGACTATATCATCCGCAGCCACGAAGTCAAGGCCGAGGGCTACGAGGTGGCTCACGGAGGCCGCTGTGTCACCGTCTTCTCTGCCCCCAACTACTGCGACCAGATGGGGAACAAAGCCTCC",
            "          ---                                                                                            "
         );
        CPPUNIT_ASSERT_EQUAL(std::string("94M"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,2), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(58, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("8M3D86M"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,2), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(3, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(1U, result.overlappingGapsFilter_.overlapsCount());
        CPPUNIT_ASSERT_EQUAL(0xfffU, result.overlappingGapsFilter_.overlap(0));
    }*/



    {   // As the deletion pushes the part of the read outside the contig, the following ugly cigar used to get formed: 222S3M6D0M25S
        io::FragmentHeader fragmentMetadata;
        fragmentMetadata.lowClipped_ = 222;
        fragmentMetadata.highClipped_ = 2;
        const RealignResult result = realign(3,4,
             "GAATGGAATAAAATGGAATCATCCAATGGAAGAGAATTGAATCATCATCGAATGGAATGGAATAGAATCATCGAAAGGAATCGAATAGAATCATCAAATGAAA"
             "TCGAATGGAATCATCATTGTATAGAATCGAATAGAATCAACATCAAATGGAATCAAATGGAATCATCATCGAATGGAATCGAATGGAATCATCATCAAATGGA"
             "ATCGAATGGAATCATCATCAAATGGAATCAAAAATAACCATCAT",
             "GAATGGAATGAAATGGAATCATCCAATGGAAGAGAATTGAATCATCATCGAATGGAATCGAATAGAATCATCGAAAGGAATCTACTAGAATCATCAAATGAAA"
             "TCGAATGGAATCATCATTGAATAGAATCCAATGGAATCAACATCAAATGGAATCAAATGGAATCATCATCGAATGGAATCGAATGGAATCATCATCAAATGGA"
             "ATCGAATGGAATCATCCAATATAATAGAATTGAATCACCATCGAACGG",
             //             "AATCGAATAGAATCATCGAATGAACTCAAATGGAATCATCATTGAATGGAATCGA"
             //             "ATCATCAACGAATGGAATTGAATGGAATCATAGAATGGAATCCAATGTAACCATCATCGAATTGAACCCAATGGAATCATTAAATGGACTCGAATAGAATCAT"
             //             "CGAATGGACTCGAATGGAATCATCATCGAATGGAATAGAATGGAGTTGGAATCGAATGGAATCATCATCA",
             "                                                              -----------------------"
                                                                           "                                         "
             "                                                                                                       "
             "                   ------",
             fragmentMetadata
         );
        CPPUNIT_ASSERT_EQUAL(std::string("222S26M2S"), result.originalCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,222), result.originalPos_);
        CPPUNIT_ASSERT_EQUAL(15, int(result.originalEditDistance_));

        CPPUNIT_ASSERT_EQUAL(std::string("222S3M25S"), result.realignedCigar_);
        CPPUNIT_ASSERT_EQUAL(reference::ReferencePosition(0,245), result.realignedPos_);
        CPPUNIT_ASSERT_EQUAL(0, int(result.realignedEditDistance_));

        CPPUNIT_ASSERT_EQUAL(0U, result.overlappingGapsFilter_.overlapsCount());
    }

}

