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
