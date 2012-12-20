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
 ** \file FastaReader.cpp
 **
 ** Component to read FASTA files.
 **
 ** \author Come Raczy
 **/

#include "io/FastaReader.hh"

namespace isaac
{
namespace io
{

FastaReader &FastaReader::get(char &base, bool &newContig)
{
    char tmp;
    newContig = false;
    while (std::ifstream::get(tmp) && ('\r' == tmp || '\n' == tmp)) ;

    while (*this && '>' == tmp)
    {
        ignore(1000000000, '\n');
        get(base, newContig);
        newContig = true;
        return *this;
    }

    if (*this)
    {
        base = tmp;
    }
    assert('\n' != base && '\r' != base);
    return *this;
}

MultiFastaReader::MultiFastaReader(const std::vector<boost::filesystem::path> &fastaPathList)
    : fastaPathList_(fastaPathList)
    , nextFastaPath_(fastaPathList_.begin())
    , contigId_(-1)
{
    std::cerr << *nextFastaPath_ << std::endl;
    open(nextFastaPath_->string().c_str());
    ++nextFastaPath_;
}

MultiFastaReader &MultiFastaReader::get(char &base, bool &newContig)
{
    char tmp;
    FastaReader::get(tmp, newContig);
    if (!*this)
    {
        close();
        if (fastaPathList_.end() != nextFastaPath_)
        {
            clear();
            std::cerr << *nextFastaPath_ << std::endl;
            open(nextFastaPath_->string().c_str());
            ++nextFastaPath_;
            if (get(base, newContig))
            {
                newContig = true;
                ++contigId_;
            }
            return *this;
        }
    }
    else
    {
        base = tmp;
        if (newContig)
        {
            ++contigId_;
        }
    }
    return *this;
}

} // namespace io
} // namespace isaac
