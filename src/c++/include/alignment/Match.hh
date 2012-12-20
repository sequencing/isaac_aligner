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
 ** \file Match.hh
 **
 ** \brief Abstract representation of the match of a seed to a reference position.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_HH
#define iSAAC_ALIGNMENT_MATCH_HH

#include <iostream>

#include "oligo/Kmer.hh"
#include "alignment/SeedId.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief a component that represents a match from a seed to a given reference location.
 **
 **/
struct Match
{
    Match(alignment::SeedId seedId = alignment::SeedId(0),
          reference::ReferencePosition location = reference::ReferencePosition(0)) : seedId(seedId), location(location) {}
    alignment::SeedId seedId;
    reference::ReferencePosition location;

    alignment::SeedId getSeedId() const
    {
        return seedId;
    }

    unsigned long getTile() const
    {
        return seedId.getTile();
    }

    unsigned long getTileBarcode() const
    {
        return seedId.getTileBarcode();
    }

    unsigned long getCluster() const
    {
        return seedId.getCluster();
    }

    unsigned long getBarcode() const
    {
        return seedId.getBarcode();
    }

    bool isNoMatch() const {return location.isNoMatch();}
    bool isTooManyMatch() const {return location.isTooManyMatch();}

};

inline std::ostream &operator<<(std::ostream &os, const Match &match)
{
    return os << "Match(" << match.seedId << ", " << match.location << ")";
}

} //namespace alignment
} //namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_HH
