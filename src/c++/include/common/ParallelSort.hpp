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

namespace isaac
{
namespace common
{

template <class Iterator, class Compare>
void parallelSort (Iterator begin, Iterator end, const Compare &comp)
{
    // TODO: note that http://gcc.gnu.org/bugzilla/show_bug.cgi?id=40852 has been spotted with gcc 4.3
    // so far the only well-tested gcc versions are 4.6.1 and 4.6.2
    // TODO: gcc (GCC) 4.4.5 20110214 (Red Hat 4.4.5-6) parallel sort implemenatation seems to leak memory
//    ISAAC_THREAD_CERR << "__gnu_parallel:sort for " << end - begin << " elements " << std::endl;
    __gnu_parallel::sort(begin, end, comp);
//    ISAAC_THREAD_CERR << "__gnu_parallel:sort done for " << end - begin << " elements " << std::endl;
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
