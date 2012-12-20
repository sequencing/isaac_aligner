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
 ** \file sortReference.cpp
 **
 ** The main for the contigs metadata printer.
 **
 ** \author Come Raczy
 **/

#include "options/PrintContigsOptions.hh"
#include "reference/ContigsPrinter.hh"

void printContigs(const isaac::options::PrintContigsOptions &options);

int main(int argc, char *argv[])
{
    isaac::common::run(printContigs, argc, argv);
}

void printContigs(const isaac::options::PrintContigsOptions &options)
{
    isaac::reference::ContigsPrinter contigsPrinter(
        options.originalMetadataPath,
        options.genomeFile);
    contigsPrinter.run();
}
