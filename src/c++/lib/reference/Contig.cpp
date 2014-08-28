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
 ** \file Contig.
 **
 ** \brief See Contig.hh
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "reference/Contig.hh"

namespace isaac
{
namespace reference
{

size_t genomeLength(const std::vector<Contig> &contigList)
{
    using boost::lambda::_1;
    using boost::lambda::_2;
    using boost::lambda::bind;
    return std::accumulate(
        contigList.begin(), contigList.end(),
        size_t(0), bind<size_t>(std::plus<size_t>(), _1, bind(&Contig::getLength, _2)));
}

} // namespace reference
} // namespace isaac
