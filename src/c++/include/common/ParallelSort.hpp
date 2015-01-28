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
 ** \file ParallelSort.hpp
 **
 ** Formerly, a homemade implementation of parallel sort. Currently just a place
 ** where the redirection to gnu parallel sort is hidden.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_PARALLEL_SORT_HPP
#define iSAAC_COMMON_PARALLEL_SORT_HPP

#include <vector>
#include <parallel/algorithm>

#include "common/Debug.hh"
#include "common/Threads.hpp"

namespace isaac
{
namespace common
{

//template <typename T, typename P> const T& ma(const T& left, const T& right, P p)
//{
//    return p(left, right) ? right : left;
//}

//template <typename T, typename P> const T& mi(const T& left, const T& right, P p)
//{
//    return p(left, right) ? left : right;
//}

template <typename IteratorT, class Compare>
class ParallelSorter
{
    boost::mutex m_;
    boost::condition_variable c_;
    int partitioningJobs_;

    struct Subjob
    {
        Subjob(){;}
        Subjob(IteratorT begin, IteratorT end) :
            begin_(begin), end_(end){}
        IteratorT begin_;
        IteratorT end_;
        bool operator <(const Subjob &that) const {return size() > that.size();}
        std::size_t size() const {return std::distance(begin_, end_);}
    };
    typedef std::vector<Subjob> Subjobs;

    void thread(Subjobs &subjobs, const unsigned long minsize, const Compare &comp)
    {
        boost::unique_lock<boost::mutex> lock(m_);
        while (true)
        {
            if (subjobs.empty())
            {
                if (!partitioningJobs_)
                {
                    break;
                }
                // there is a job partitioning a chunk of data. Wait for it produce the results
                c_.wait(lock);
                continue;
            }
            std::pop_heap(subjobs.begin(), subjobs.end());

            Subjob ourJob = subjobs.back();
            subjobs.pop_back();

            while (true)
            {
                if (ourJob.size() <= minsize)
                {
                    isaac::common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                    std::sort(ourJob.begin_, ourJob.end_, comp);
                    break;
                }
                else
                {
// for some reason, median of 3 pivot makes sort run longer on the test data...
//                    typename IteratorT::value_type pivots[] = {*(sj.begin_ + sj.size() / 2), *(sj.begin_), *(sj.end_ - 1)};
//                    using std::swap;
//                    if (comp(pivots[1], pivots[0]))
//                        swap(pivots[0], pivots[1]);
//                    if (comp(pivots[2], pivots[0]))
//                        swap(pivots[0], pivots[2]);
//                    if (comp(pivots[2], pivots[1]))
//                        swap(pivots[1], pivots[2]);
//                    const typename IteratorT::value_type pivot = pivots[1];

//                    const typename IteratorT::value_type pivot = ma(mi(pivots[0], pivots[1], comp), mi(ma(pivots[0], pivots[1], comp), pivots[2], comp), comp);

                    const typename IteratorT::value_type pivot = *(ourJob.begin_ + ourJob.size() / 2);

                    ++partitioningJobs_;
                    IteratorT midpoint;
                    {
                        isaac::common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                        midpoint = std::partition(ourJob.begin_, ourJob.end_, boost::bind(comp, _1, pivot));
                    }
                    --partitioningJobs_;
                    if (midpoint != ourJob.begin_ && midpoint != ourJob.end_)
                    {
    //                  std::cerr << threadNumber << " partitioned " << std::distance(begin, end) << " to " <<  std::distance(begin, midpoint) << " and " << std::distance(midpoint, end) << " minsize " << minsize << std::endl;
                        // prioritize partitioning above sorting
                        if (std::distance(ourJob.begin_, midpoint) < std::distance(midpoint, ourJob.end_))
                        {
                            subjobs.push_back(Subjob(ourJob.begin_, midpoint));
                            ourJob.begin_ = midpoint;
                        }
                        else
                        {
                            subjobs.push_back(Subjob(midpoint, ourJob.end_));
                            ourJob.end_ = midpoint;
                        }
                        std::push_heap(subjobs.begin(), subjobs.end());
                    }
                    else
                    {
                        // failed to partition, just let everybody else know this job is not going to partition, sort and go back to checking the subjobs queue
                        c_.notify_all();
                        isaac::common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                        std::sort(ourJob.begin_, ourJob.end_, comp);
                        break;
                    }
                    c_.notify_all();
                }
            }
        }
    }
public:
    ParallelSorter() : partitioningJobs_(0){}

    /**
     * \brief performs in-place sort on multiple threads. Should handle well randomly distributed data and sorted data.
     *        Uses a bit of dynamic memory for jobs priority queue.
     */

    void sort(IteratorT begin, IteratorT end, const Compare &comp, isaac::common::ThreadVector &threads, const unsigned threadsMax)
    {
        Subjobs subjobs(1, Subjob(begin, end));
        threads.execute(boost::bind(
            &ParallelSorter::thread, this,
            boost::ref(subjobs),
            // no reason to make single stretch shorter than size/threads except for
            // when some sort quicker than others. / 10/ allows a bit of rebalancing
            // when amount of work turns out to be inequal.
            std::distance(begin, end) / threads.size() / 100,
            boost::ref(comp)),
                        threadsMax);
    }
};

template <class Iterator, class Compare>
void parallelSort (Iterator begin, Iterator end, const Compare &comp, isaac::common::ThreadVector &threads, const unsigned threadsMax)
{
    ParallelSorter<Iterator, Compare> sorter;
    sorter.sort(begin, end, comp, threads, threadsMax);
}

template <class Iterator, class Compare>
void parallelSort (Iterator begin, Iterator end, const Compare &comp)
{
    isaac::common::ThreadVector threads(boost::thread::hardware_concurrency());
    parallelSort(begin, end, comp, threads, threads.size());
}

/**
 ** \brief Calling to the gnu parallel sort implementation. NOTE!!!
 **        Seems to require as much extra dynamic memory as there is data to sort in gcc 4.6.
 **/
template <class T, class Compare>
void parallelSort (std::vector<T> &v, const Compare &comp)
{
    parallelSort(v.begin(), v.end(), comp);
}

} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_PARALLEL_SORT_HPP
