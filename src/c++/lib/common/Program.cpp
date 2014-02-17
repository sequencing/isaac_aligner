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
 ** \file Program.cpp
 **
 ** Implementation of the skeleton of all c++ programs.
 **
 ** \author Come Raczy
 **/

#include <sys/ioctl.h>

#include <boost/foreach.hpp>

#include "common/Program.hh"

namespace isaac
{
namespace common
{

#define HELP_STR "help"
#define HELP_MD_STR "help-md"
#define HELP_DEFAULTS_STR "help-defaults"

Options::Options()
{
    namedOptions_.add_options()(HELP_STR ",h", "produce help message and exit");
    namedOptions_.add_options()(HELP_MD_STR, "produce help message pre-formatted as a markdown file section and exit");
    namedOptions_.add_options()(HELP_DEFAULTS_STR, "produce tab-delimited list of command line options and their default values");
    namedOptions_.add_options()("version,v", "print program version information");
}

Options::Action Options::parse(int argc, const char * const argv[])
{
    try
    {
        bpo::options_description allOptions("Allowed options");
        allOptions.add(namedOptions_).add(unnamedOptions_);
        vm_.clear();
        bpo::store(
                bpo::command_line_parser(argc, argv). options(allOptions).positional(
                        positionalOptions_).run(), vm_);
        bpo::notify(vm_);
        if (vm_.count(HELP_STR) || vm_.count(HELP_MD_STR) || vm_.count(HELP_DEFAULTS_STR))
        {
            return HELP;
        }
        else if (vm_.count("version"))
        {
            return VERSION;
        }
        postProcess(vm_);
    }
    catch (const boost::program_options::multiple_values &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
        return ABORT;
    }
    catch (const boost::program_options::multiple_occurrences &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
        return ABORT;
    }
    catch (const boost::program_options::required_option &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
        return ABORT;
    }
    catch (const boost::exception &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << boost::diagnostic_information(e) << std::endl;
        return ABORT;
    }
    catch (const std::exception &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << e.what() << std::endl;
        return ABORT;
    }
    return RUN;
}

bool compareOptionName(
    const boost::shared_ptr<boost::program_options::option_description> &left,
    const boost::shared_ptr<boost::program_options::option_description> &right)
{
    return left->long_name() < right->long_name();
}

std::string Options::helpDefaults(const OptionDescriptionPtrs &sortedOptions) const
{
    std::string ret;
    BOOST_FOREACH(const OptionDescriptionPtr &odPtr, sortedOptions)
    {
        ret += odPtr->long_name() + "\t" + odPtr->format_parameter() + "\n";
    }
    return ret;
}

static winsize ioctlSTDERR_TIOCGWINSZ()
{
    winsize ret = {0,0,0,0};
    if (-1 == ioctl(STDERR_FILENO, TIOCGWINSZ, &ret))
    {
        ret.ws_col = bpo::options_description::m_default_line_length;
    }
    return ret;
}

std::string Options::help(const OptionDescriptionPtrs &sortedOptions, const bool markdown) const
{
    std::ostringstream os;
    if (!markdown)
    {
        os << this->usagePrefix() << std::endl << std::endl;
    }

    //markdown lines are prepended by two spaces
    const unsigned effectiveLineLength = markdown ? MARKDOWN_LINE_LENGTH - 2 : ioctlSTDERR_TIOCGWINSZ().ws_col;
    bpo::options_description printedDescriptions(!markdown ? "Command line options" : "", effectiveLineLength, effectiveLineLength - 50);
    BOOST_FOREACH(const OptionDescriptionPtr &odPtr, sortedOptions)
    {
        printedDescriptions.add(odPtr);
    }
    os << printedDescriptions << std::endl;

    if (markdown)
    {
        std::vector<std::string> lines;
        std::string str = os.str();
        boost::algorithm::split(lines, str, boost::algorithm::is_any_of("\n\r"));
        os.str("");
        BOOST_FOREACH(const std::string &line, lines)
        {
            // pre-pend two spaces to the two spaces that boost adds so that we get 4 spaces for markdown to
            // recognise pre-formatted text.
            os << "  " << line << "\n";
        }
    }
    return os.str();
}

std::string Options::usage() const
{
    OptionDescriptionPtrs sortedOptions = namedOptions_.options();
    std::sort(sortedOptions.begin(), sortedOptions.end(), compareOptionName);

    if (vm_.count(HELP_DEFAULTS_STR))
    {
        return helpDefaults(sortedOptions);
    }

    return help(sortedOptions, vm_.count(HELP_MD_STR));
}

} // namespace common
} // namespace isaac
