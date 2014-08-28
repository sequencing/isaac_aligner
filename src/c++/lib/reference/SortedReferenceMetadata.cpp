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
 ** \file SortedReferenceMetadata.cpp
 **
 ** Information about the pre-processed reference data files.
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{

const unsigned SortedReferenceMetadata::OLDEST_SUPPORTED_REFERENCE_FORMAT_VERSION;
const unsigned SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION;


void SortedReferenceMetadata::putContig(
    const unsigned long genomicOffset,
    const std::string& name,
    const boost::filesystem::path &sequencePath,
    const unsigned long byteOffset,
    const unsigned long byteSize,
    const unsigned long totalBases,
    const unsigned long acgtBases,
    const unsigned index,
    const unsigned karyotypeIndex,
    const std::string &bamSqAs,
    const std::string &bamSqUr,
    const std::string &bamM5
    )
{
    contigs_.push_back(Contig(
        index, karyotypeIndex,
        name,
        sequencePath,
        byteOffset,
        byteSize,
        genomicOffset,
        totalBases, acgtBases,
        bamSqAs,
        bamSqUr,
        bamM5));
}
void SortedReferenceMetadata::addMaskFile(
    const unsigned seedLength,
    const unsigned int maskWidth,
    const unsigned mask, const boost::filesystem::path &filePath,
    const size_t kmers)
{
    if (maskFiles_.empty())
    {
        defaultMaskWidth_ = maskWidth;
    }
    else
    {
        ISAAC_ASSERT_MSG(maskWidth == defaultMaskWidth_, "Mask width must match");
    }
    maskFiles_[seedLength].push_back(MaskFile(filePath, maskWidth, mask, kmers));
}

SortedReferenceMetadata::Contigs SortedReferenceMetadata::getKaryotypeOrderedContigs() const
{
    Contigs ret = getContigs();
    std::sort(ret.begin(), ret.end(),
              boost::bind(&reference::SortedReferenceMetadata::Contig::karyotypeIndex_, _1) <
              boost::bind(&reference::SortedReferenceMetadata::Contig::karyotypeIndex_, _2));
    return ret;
}

unsigned long SortedReferenceMetadata::getTotalKmers(const unsigned seedLength) const
{
    const unsigned maskWidth = getDefaultMaskWidth();
    unsigned long ret = 0;
    BOOST_FOREACH(const MaskFile &mask, maskFiles_.at(seedLength))
    {
        if (maskWidth == mask.maskWidth)
        {
            ret += mask.kmers;
        }
    }
    return ret;
}

void SortedReferenceMetadata::merge(SortedReferenceMetadata &that)
{
    ISAAC_ASSERT_MSG(formatVersion_ == that.formatVersion_, "Incompatible formats: " << formatVersion_ << " and " << that.formatVersion_ << " cannot be merged");
    if(!contigs_.empty() && !that.contigs_.empty())
    {
        if (contigs_ != that.contigs_)
        {
            BOOST_THROW_EXCEPTION(common::FeatureNotAvailable("Cannot merge references with different contig lists"));
        }
    }
    else
    {
        contigs_.insert(contigs_.end(), that.contigs_.begin(), that.contigs_.end());
    }

    if (maskFiles_.empty())
    {
        defaultMaskWidth_ = that.defaultMaskWidth_;
    }
    ISAAC_ASSERT_MSG(defaultMaskWidth_ == that.defaultMaskWidth_, "Incompatible mask widths: " << defaultMaskWidth_ << " and " << that.defaultMaskWidth_ << " cannot be merged");

    BOOST_FOREACH(AllMaskFiles::value_type &seedMaskFiles, that.maskFiles_)
    {
        MaskFiles &maskFiles = getMaskFileList(seedMaskFiles.first);
        maskFiles.insert(maskFiles.end(), seedMaskFiles.second.begin(), seedMaskFiles.second.end());
    }
}

bool SortedReferenceMetadata::singleFileReference() const {
    ISAAC_ASSERT_MSG(0 != getContigsCount(), "Contigs list must not be empty.");
    Contigs::const_iterator different = std::find_if(contigs_.begin() + 1, contigs_.end(),
                                                     boost::bind(&Contig::filePath_, _1) != contigs_.front().filePath_);
    return contigs_.end() == different;
}

} // namespace reference
} // namespace isaac

