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
 ** \file AlignWorkflowSerialization.hh
 **
 ** \brief Xml Serialization of Workflow state.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_WORKFLOW_SERIALIZATION_H
#define ISAAC_WORKFLOW_SERIALIZATION_H

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/serialization/vector.hpp>

#include "alignment/MatchTally.hh"
#include "common/BoostArchiveHelpers.hh"
#include "workflow/AlignWorkflow.hh"

/**
 * \brief serialization implementation types that don't require private member access
 * can be placed in the boost::serialization namespace
 */
namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, isaac::alignment::MatchTally::FileTally &ft, const unsigned int version)
{
    ar & boost::serialization::make_nvp("path", ft.path_);
    ar & boost::serialization::make_nvp("count", ft.matchCount_);
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


template <class Archive>
void serialize(Archive &ar, TemplateLengthStatistics &tls, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(tls.min_);
    ar & BOOST_SERIALIZATION_NVP(tls.max_);
    ar & BOOST_SERIALIZATION_NVP(tls.median_);
    ar & BOOST_SERIALIZATION_NVP(tls.lowStdDev_);
    ar & BOOST_SERIALIZATION_NVP(tls.highStdDev_);
    ar & BOOST_SERIALIZATION_NVP(tls.bestModels_[0]);
    ar & BOOST_SERIALIZATION_NVP(tls.bestModels_[1]);
    ar & BOOST_SERIALIZATION_NVP(tls.stable_);
    ar & BOOST_SERIALIZATION_NVP(tls.mateMin_);
    ar & BOOST_SERIALIZATION_NVP(tls.mateMin_);
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
    ar & BOOST_SERIALIZATION_NVP(tm.clusterCount_);
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
    ar & BOOST_SERIALIZATION_NVP(bbm.barcodeProjectIndex_);
    ar & BOOST_SERIALIZATION_NVP(bbm.projectIndexMax_);
    ar & BOOST_SERIALIZATION_NVP(bbm.barcodeSampleIndex_);
    ar & BOOST_SERIALIZATION_NVP(bbm.samplePaths);
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
            ar & BOOST_SERIALIZATION_NVP(a.barcodeTemplateLengthStatistics_);
            if (AlignWorkflow::BamDone <= a.state_)
            {
                ar & BOOST_SERIALIZATION_NVP(a.barcodeBamMapping_);
            }
        }
    }
}

inline void save(const boost::filesystem::path &stateFilePath, const AlignWorkflow &aligner)
{
    const boost::filesystem::path tmp = stateFilePath.string() + ".tmp";

    ISAAC_THREAD_CERR << "Saving workflow state to " << stateFilePath << std::endl;

    {
        std::ofstream ofs(tmp.string().c_str());
        boost::archive::text_oarchive oa(ofs);
        oa << BOOST_SERIALIZATION_NVP(aligner);
    }
    boost::filesystem::rename(tmp, stateFilePath);

    ISAAC_THREAD_CERR << "Saving workflow state done to " << stateFilePath << std::endl;
}

inline void load(const boost::filesystem::path &stateFilePath, AlignWorkflow &aligner)
{
    ISAAC_THREAD_CERR << "Loading workflow state from " << stateFilePath << std::endl;
    std::ifstream ifs(stateFilePath.string().c_str());

    boost::archive::text_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(aligner);

    ISAAC_THREAD_CERR << "Loading workflow state done from " << stateFilePath << std::endl;
}

} //namespace workflow
} //namespace isaac

#endif //ISAAC_WORKFLOW_SERIALIZATION_H
