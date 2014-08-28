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
//    const std::size_t valueSize = sizeof(typename std::iterator_traits<Iterator>::value_type);
//    const std::size_t elements = std::distance(begin, end);
//    const std::size_t bytes = elements * valueSize;
//    ISAAC_THREAD_CERR << "__gnu_parallel:sort for " << elements << " elements (" << bytes << " bytes)" << std::endl;
//    ISAAC_TRACE_STAT("parallelSort before ")
    __gnu_parallel::sort(begin, end, comp);
//    ISAAC_TRACE_STAT("parallelSort after ")
//    ISAAC_THREAD_CERR << "__gnu_parallel:sort done for " << elements << " elements (" << bytes << " bytes)" << std::endl;
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
