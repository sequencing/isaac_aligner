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
 ** \file AlignWorkflowSerialization.hh
 **
 ** \brief Xml Serialization of Workflow state.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_WORKFLOW_SERIALIZATION_H
#define ISAAC_WORKFLOW_SERIALIZATION_H

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_free.hpp>

#include "alignment/MatchTally.hh"
#include "workflow/AlignWorkflow.hh"

BOOST_SERIALIZATION_SPLIT_FREE(boost::filesystem::path)

/**
 * \brief serialization implementation types that don't require private member access
 * can be placed in the boost::serialization namespace
 */
namespace boost {
namespace serialization {

template <class Archive>
void save(Archive &ar, const boost::filesystem::path &p, const unsigned int version)
{
    ar << boost::serialization::make_nvp("path", p.string());
}

template <class Archive>
void load(Archive &ar, boost::filesystem::path &p, const unsigned int version)
{
    std::string tmp;
    ar >> boost::serialization::make_nvp("path", tmp);
    p = tmp;
}

template <class Archive>
void serialize(Archive &ar, std::pair<std::string, std::string> &ps, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(ps.first);
    ar & BOOST_SERIALIZATION_NVP(ps.second);
}

template <class Archive>
void serialize(Archive &ar, isaac::alignment::MatchTally::FileTally &ft, const unsigned int version)
{
    ar & boost::serialization::make_nvp("path", ft.first);
    ar & boost::serialization::make_nvp("count", ft.second);
//    ar & boost::serialization::make_nvp("barcodeTally_",
//                                        boost::serialization::base_object<std::vector<unsigned long> >(
//                                            ft.barcodeTally_));
    ar & boost::serialization::make_nvp("barcodeTally_", ft.barcodeTally_);
}


} //namespace serialization
} //namespace boost


/**
 * \brief serialization for types that declare serialize as friend, has to be in the same namespace as
 * the type
 */
namespace isaac {

namespace reference {

template <class Archive>
void serialize(Archive &ar, ReferencePosition &pos, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(pos.value_);
}

}
namespace alignment {

template <class Archive>
void serialize(Archive &ar, MatchTally &mt, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(mt.allTallies_);
}

template <class Archive>
void serialize(Archive &ar, BarcodeCounts &bc, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(bc.elements_);
    ar & BOOST_SERIALIZATION_NVP(bc.gaps_);
    ar & BOOST_SERIALIZATION_NVP(bc.cigarLength_);
}

template <class Archive>
void serialize(Archive &ar, BinChunk &bch, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(bch.barcodeBreakdown_);
    ar & BOOST_SERIALIZATION_NVP(bch.dataSize_);
}

template <class Archive>
void serialize(Archive &ar, BinDataDistribution &bdd, const unsigned int version)
{
    ar & boost::serialization::make_nvp<std::vector<BinChunk> >("dataDistribution_", bdd);
    ar & BOOST_SERIALIZATION_NVP(bdd.chunkSize_);
}

template <class Archive>
void serialize(Archive &ar, BinMetadata &bm, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(bm.binIndex_);
    ar & BOOST_SERIALIZATION_NVP(bm.binStart_);
    ar & BOOST_SERIALIZATION_NVP(bm.length_);
    ar & BOOST_SERIALIZATION_NVP(bm.binFilePath_);
    ar & BOOST_SERIALIZATION_NVP(bm.fIdxFilePath_);
    ar & BOOST_SERIALIZATION_NVP(bm.rIdxFilePath_);
    ar & BOOST_SERIALIZATION_NVP(bm.seIdxFilePath_);
    ar & BOOST_SERIALIZATION_NVP(bm.dataSize_);
    ar & BOOST_SERIALIZATION_NVP(bm.dataOffset_);
    ar & BOOST_SERIALIZATION_NVP(bm.seIdxElements_);
    ar & BOOST_SERIALIZATION_NVP(bm.rIdxElements_);
    ar & BOOST_SERIALIZATION_NVP(bm.fIdxElements_);
    ar & BOOST_SERIALIZATION_NVP(bm.nmElements_);
    ar & BOOST_SERIALIZATION_NVP(bm.dataDistribution_);
}

template <class Archive>
void serialize(Archive &ar, alignment::BinMetadataList &bml, const unsigned int version)
{
    ar & boost::serialization::make_nvp(
        "bml",
        boost::serialization::base_object<std::vector<alignment::BinMetadata> >(bml));
}

} //namespace alignment

namespace flowcell {

template <class Archive>
void serialize(Archive &ar, TileMetadata &tm, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(tm.flowcellId_);
    ar & BOOST_SERIALIZATION_NVP(tm.flowcellIndex_);
    ar & BOOST_SERIALIZATION_NVP(tm.tile_);
    ar & BOOST_SERIALIZATION_NVP(tm.tileString_);
    ar & BOOST_SERIALIZATION_NVP(tm.lane_);
    ar & BOOST_SERIALIZATION_NVP(tm.laneString_);
    ar & BOOST_SERIALIZATION_NVP(tm.baseCallsPath_);
    ar & BOOST_SERIALIZATION_NVP(tm.clusterCount_);
    ar & BOOST_SERIALIZATION_NVP(tm.compression_);
    ar & BOOST_SERIALIZATION_NVP(tm.index_);
}

template <class Archive>
void serialize(Archive &ar, flowcell::TileMetadataList &tml, const unsigned int version)
{
    ar & boost::serialization::make_nvp(
        "tml",
        boost::serialization::base_object<std::vector<flowcell::TileMetadata> >(tml));
}

} //namespace flowcell

namespace build {

template <class Archive>
void serialize(Archive &ar, BarcodeBamMapping &bbm, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(bbm.first);
    ar & BOOST_SERIALIZATION_NVP(bbm.second);
}

}

namespace workflow {

namespace alignWorkflow {
template <class Archive>
void serialize(Archive &ar, FoundMatchesMetadata &fmm, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(fmm.tileMetadataList_);
    ar & BOOST_SERIALIZATION_NVP(fmm.matchTally_);
    ar & boost::serialization::make_nvp("fmm.matchDistribution_",
        boost::serialization::base_object<std::vector<std::vector<unsigned> > >(fmm.matchDistribution_));
}

} //namespace alignWorkflow

template <class Archive>
void serialize(Archive &ar, AlignWorkflow &a, const unsigned int version)
{
//    ar & BOOST_SERIALIZATION_NVP(a.repeatThreshold_);
    ar & BOOST_SERIALIZATION_NVP(a.state_);
    if (AlignWorkflow::MatchFinderDone <= a.state_)
    {
        ar & BOOST_SERIALIZATION_NVP(a.foundMatchesMetadata_);

        if (AlignWorkflow::MatchSelectorDone <= a.state_)
        {
            ar & BOOST_SERIALIZATION_NVP(a.selectedMatchesMetadata_);
            if (AlignWorkflow::BamDone <= a.state_)
            {
                ar & BOOST_SERIALIZATION_NVP(a.barcodeBamMapping_);
            }
        }
    }
}

inline void save(const boost::filesystem::path &stateFilePath, const AlignWorkflow &aligner)
{
    std::ofstream ofs(stateFilePath.string().c_str());

    boost::archive::xml_oarchive oa(ofs);
    oa << BOOST_SERIALIZATION_NVP(aligner);
}

inline void load(const boost::filesystem::path &stateFilePath, AlignWorkflow &aligner)
{
    std::ifstream ifs(stateFilePath.string().c_str());

    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(aligner);
}

} //namespace workflow
} //namespace isaac

#endif //ISAAC_WORKFLOW_SERIALIZATION_H
