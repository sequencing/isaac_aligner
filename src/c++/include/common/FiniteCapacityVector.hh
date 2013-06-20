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
 ** \file FiniteCapacityVector.hh
 **
 ** Something that behaves more or less like std::vector but does not use dynamic memory
 ** and thus has fixed capacity().
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_FINITE_CAPACITY_VECTOR_HH
#define iSAAC_COMMON_FINITE_CAPACITY_VECTOR_HH

#include <boost/numeric/ublas/storage.hpp>

namespace isaac
{
namespace common
{

// Bounded array - with allocator for size_type and difference_type
template<class T, std::size_t N, class ALLOC = std::allocator<T> >
class FiniteCapacityVector: boost::numeric::ublas::bounded_array<T, N, ALLOC>
{
    typedef boost::numeric::ublas::bounded_array<T, N, ALLOC> BaseType;
public:
    using BaseType::begin;
    using BaseType::empty;
    using BaseType::iterator;
    using BaseType::const_iterator;
    using BaseType::const_reference;
    using BaseType::end;
    using BaseType::operator [];
    using BaseType::resize;
    using BaseType::size;

    FiniteCapacityVector()
    {

    }

    FiniteCapacityVector(std::size_t s, const T& v)
    {
        resize(s, v);
    }

    void clear()
    {
        BaseType::resize(0);
    }

    static std::size_t capacity() {return N;}

    void push_back(const T& x)
    {
        BaseType::resize(BaseType::size() + 1, x);
    }

    T pop_back()
    {
        const T ret = back();
        BaseType::resize(BaseType::size() - 1);
        return ret;
    }

    T& front()
    {
        BOOST_UBLAS_CHECK (0 < size(), boost::numeric::ublas::bad_index ());
        return *begin();
    }

    const T& front() const
    {
        BOOST_UBLAS_CHECK (0 < BaseType::size(), boost::numeric::ublas::bad_index ());
        return *begin();
    }

    T& back()
    {
        BOOST_UBLAS_CHECK (0 < size(), boost::numeric::ublas::bad_index ());
        return *BaseType::rbegin();
    }

    const T& back() const
    {
        BOOST_UBLAS_CHECK (0 < BaseType::size(), boost::numeric::ublas::bad_index ());
        return *BaseType::rbegin();
    }

    void erase(typename BaseType::iterator b, typename BaseType::iterator e)
    {
        resize(std::distance(begin(), std::copy(e, end(), b)));
    }
};

} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_FINITE_CAPACITY_VECTOR_HH
