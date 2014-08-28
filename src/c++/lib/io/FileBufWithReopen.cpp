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
 ** \file FileBufWithReopen.cpp
 **
 ** Same as std::filebuf but has reopen.
 **
 ** \author Roman Petrovski
 **/

#include <vector>

#include <boost/assign.hpp>

namespace isaac
{
namespace io
{

extern const std::vector<const char*> iosBaseToStdioOpenModesTranslationTable = boost::assign::list_of<const char*>
                        // inspired by /usr/include/c++/4.3.2/fstream
                        // | ios_base Flag combination      |
                        // |binary  in  out  trunc  app     |
            (0      )   // |                                |
            ("a"    )   // |                         +      |
            (0      )   // |                   +            |
            (0      )   // |                   +     +      |
            ("w"    )   // |             +                  |
            ("a"    )   // |             +           +      |
            ("w"    )   // |             +     +            |
            (0      )   // |             +     +     +      |
            ("r"    )   // |         +                      |
            ("a+"   )   // |         +               +      |
            (0      )   // |         +         +            |
            (0      )   // |         +         +     +      |
            ("r+"   )   // |         +   +                  |
            ("a+"   )   // |         +   +           +      |
            ("w+"   )   // |         +   +     +            |
            (0      )   // |         +   +     +     +      |
            (0      )   // |   +                            |
            ("ab"   )   // |   +                     +      |
            (0      )   // |   +               +            |
            (0      )   // |   +               +     +      |
            ("wb"   )   // |   +         +                  |
            ("ab"   )   // |   +         +           +      |
            ("wb"   )   // |   +         +     +            |
            (0      )   // |   +         +     +     +      |
            ("rb"   )   // |   +     +                      |
            ("a+b"  )   // |   +     +               +      |
            (0      )   // |   +     +         +            |
            (0      )   // |   +     +         +     +      |
            ("r+b"  )   // |   +     +   +                  |
            ("a+b"  )   // |   +     +   +           +      |
            ("w+b"  )   // |   +     +   +     +            |
            (0      )   // |   +     +   +     +     +      |
            ;

} // namespace io
} // namespace isaac
