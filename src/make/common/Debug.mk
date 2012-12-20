################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2012 Illumina, Inc.
##
## This software is provided under the terms and conditions of the
## Illumina Open Source Software License 1.
##
## You should have received a copy of the Illumina Open Source
## Software License 1 along with this program. If not, see
## <https://github.com/downloads/sequencing/licenses/>.
##
## The distribution includes the code libraries listed below in the
## 'redist' sub-directory. These are distributed according to the
## licensing terms governing each library.
##
################################################################################
##
## file Debug.mk
##
## brief Partial makefile providing basic debugging support
##
## - Redefines the SHELL
## - provides print finctionalities
##
## author Come Raczy
##
################################################################################

PRINTCMDGOALS:=$(filter print-%, $(MAKECMDGOALS))

ifneq (,$(PRINTCMDGOALS))
# Target to print the value of a variable and exit
print-%: ; @$(error $* is: $($*))

# This will translate the command line "make all print-X" into the dependency
# all: print-X
# as a result, X will get the target-specific value
$(filter-out print-%, $(MAKECMDGOALS)): $(filter print-%, $(MAKECMDGOALS))
endif
