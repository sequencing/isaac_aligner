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
 ** \file SampleSheetConstraints.hh
 **
 ** SampleSheet.csv grammar definition
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_MISEQ_SAMPLE_SHEET_CONSTRAINTS_HH
#define iSAAC_DEMULTIPLEXING_MISEQ_SAMPLE_SHEET_CONSTRAINTS_HH

namespace isaac
{
namespace demultiplexing
{

/**
 * \brief returns the input string or throws an exception if illegal characters are present
 */
inline const std::string &checkIllegalCharacters(const std::string &columnName, const std::string &str)
{
    static const char illegalCharacters[]           = "\t\n\r/,";
    static const std::string illegalCharactersEscaped("\\t\\n\\r/,");

    if (std::string::npos != str.find_first_of(illegalCharacters))
    {
        BOOST_THROW_EXCEPTION(
            common::IoException(errno, "Value '" + str + "' is not allowed." +
                                " The following characters are not allowed in sample sheet " + columnName +
                                " column: " + illegalCharactersEscaped));
    }
    return str;
}

} // namespace demultiplexing
} // namespace isaac

#endif //iSAAC_DEMULTIPLEXING_MISEQ_SAMPLE_SHEET_CONSTRAINTS_HH

