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
 ** \file FastaReader.cpp
 **
 ** Component to read FASTA files.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_IO_FASTA_READER_HH
#define iSAAC_IO_FASTA_READER_HH

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/utility.hpp>

namespace isaac
{
namespace io
{

class FastaReader: public std::ifstream, boost::noncopyable
{
public:
    /// read a base and set the 'newContig' flag when starting a new contig
    FastaReader &get(char &base, bool &newContig);
private:
};

class MultiFastaReader: public FastaReader
{
public:
    MultiFastaReader(const std::vector<boost::filesystem::path> &fastaPathList);
    MultiFastaReader &get(char &base, bool &newContig);
    int getContigId() {return contigId_;}
private:
    MultiFastaReader();
    const std::vector<boost::filesystem::path> fastaPathList_;
    std::vector<boost::filesystem::path>::const_iterator nextFastaPath_;
    int contigId_;
};

} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FASTA_READER_HH
