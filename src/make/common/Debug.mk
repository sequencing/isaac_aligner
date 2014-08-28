################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2014 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## BSD 2-Clause License
##
## You should have received a copy of the BSD 2-Clause License
## along with this program. If not, see
## <https://github.com/sequencing/licenses/>.
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
