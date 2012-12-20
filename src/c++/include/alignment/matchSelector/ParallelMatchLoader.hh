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
 ** \file ParallelMatchLoader.hh
 **
 ** Helper component to combat the directory access latency when loading matches.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALILGNMENT_MATCH_SELECTOR_PARALLEL_MATCH_LOADER_HH
#define iSAAC_ALILGNMENT_MATCH_SELECTOR_PARALLEL_MATCH_LOADER_HH

#include <boost/lambda/bind.hpp>
#include <boost/noncopyable.hpp>

#include "alignment/Match.hh"
#include "common/Debug.hh"
#include "common/Threads.hpp"
#include "io/MatchReader.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

/**
 ** \brief a component that reads the matches from a multiple files in parallel.
 **
 **/
class ParallelMatchLoader: boost::noncopyable
{
    mutable boost::mutex mutex_;
    common::ThreadVector &threads_;
    std::vector<io::MatchReader> threadMatchReaders_;
public:
    ParallelMatchLoader( common::ThreadVector &threads):
        threads_(threads), threadMatchReaders_(threads_.size())
    {
    }
    void load(const std::vector<MatchTally::FileTally> &fileTallyList,
              std::vector<Match> &matches)
    {
        // get the list of files and the match count from the match tally
        // Calculate the total number of matches for the tile
        using boost::lambda::_1;
        using boost::lambda::_2;
        using boost::lambda::bind;
        const unsigned long totalMatchCount = std::accumulate(
            fileTallyList.begin(), fileTallyList.end(), 0UL,
            bind(std::plus<unsigned long>(),
                 _1,
                 bind<unsigned long>(&MatchTally::FileTally::second, _2)));

        // allocate the storage for the matches
        matches.resize(totalMatchCount);
        // load all the matches


        std::vector<Match>::iterator destination = matches.begin();
        std::vector<MatchTally::FileTally>::const_iterator filesBegin = fileTallyList.begin();

        threads_.execute(boost::bind(
                &ParallelMatchLoader::threadLoadMatches, this, _1,
                boost::ref(destination),
                boost::ref(filesBegin),
                fileTallyList.end()));

        assert(matches.end() == destination);
    }

    void reservePathBuffers(const size_t reservePathLength)
    {
        std::for_each(threadMatchReaders_.begin(), threadMatchReaders_.end(),
                      boost::bind(&io::MatchReader::reservePathBuffers, _1, reservePathLength));
    }

    void unreserve()
    {
        std::vector<io::MatchReader>().swap(threadMatchReaders_);
    }
private:

    void threadLoadMatches(const unsigned threadNumber, std::vector<Match>::iterator &destination,
                           std::vector<MatchTally::FileTally>::const_iterator &file,
                           const std::vector<MatchTally::FileTally>::const_iterator filesEnd)
    {
        while(true)
        {
            std::vector<MatchTally::FileTally>::const_iterator ourFile = filesEnd;
            std::vector<Match>::iterator ourDestination;
            {
                boost::lock_guard<boost::mutex> lock(mutex_);
                while (filesEnd != file && !file->second) ++file;
                if (filesEnd != file)
                {
                    ourDestination = destination;
                    destination += file->second;
                    ourFile = file++;
                }
            }

            if (filesEnd != ourFile)
            {
                threadMatchReaders_[threadNumber].read(ourFile->first, &*ourDestination, ourFile->second);
//                ISAAC_THREAD_CERR << " loaded " << ourFile->first << " : " << ourFile->second << std::endl;
            }
            else
            {
                break;
            }
        }
    }

};

} //namspace matchSelector
} //namespace alignment
} //namespace isaac

#endif // #ifndef iSAAC_ALILGNMENT_MATCH_SELECTOR_PARALLEL_MATCH_LOADER_HH
