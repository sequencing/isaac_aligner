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
 ** generic parallel sort.
 **
 ** \author Come Raczy
 **/

#include "common/ParallelSort.hpp"

/*
 * std::inplace_merge attempts to create a _Temporary_buffer with enough space
 * for the whole range. This is greedy and is likely to starve other threads or
 * processes. To prevent this, we must make sure that the constructor of
 * _Temporary_buffer shows the restraint that we need.  Note that inplace_merge
 * is an adaptive algorithm and that by default it would only allocate as much
 * memory is available, which is usually fine, but miserably fails when there is
 * a large number of threads merging in parallel (the first thread hogs the
 * whole memory and the others starve).
 *
 * The workaround is to reimplement the template specialization for the
 * temporary buffer, limitting the size of the buffer to our needs.
 */
namespace std
{

template<>
_Temporary_buffer<std::vector<int>::iterator, int>::
_Temporary_buffer(std::vector<int>::iterator first, std::vector<int>::iterator last)
    : _M_original_len(std::distance(first, last) / 100)
    , _M_len(0), _M_buffer(0)
{
    // Workaround for a __type_traits bug in the pre-7.3 compiler.
    typedef std::__is_scalar<int>::__type _Trivial;

    try
    {
        pair<pointer, size_type> __p(get_temporary_buffer<
                                     value_type>(_M_original_len));
        _M_buffer = __p.first;
        _M_len = __p.second;
        if (_M_len > 0)
          std::uninitialized_fill_n(_M_buffer, _M_len, 0); //_M_initialize_buffer(*first, _Trivial());
    }
    catch(...)
    {
        std::return_temporary_buffer(_M_buffer);
        _M_buffer = 0;
        _M_len = 0;
        __throw_exception_again;
    }
}

} // namespace std
