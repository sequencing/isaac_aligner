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
## file Log.mk
##
## brief Partial makefile providing basic folder hierarchy creation support
##
## author Roman Petrovski
##
################################################################################

# Pattern rule to create folder hierarchies. If first mkdir fails due to a race
# condition, sleep 5 seconds and retry. Note that patterns dont' work with 
# .SECONDARY on qmake, therefore all sentinels have to be declared as .PRECIOUS.
# As result, they have to be kept when cleaning up Temp or else, everything will
# be rebuilt once the Temp is destroyed.

.PRECIOUS: %/.sentinel
%/.sentinel:
	$(MKDIR) -p $* $(OR) ( $(SLEEP) 5 $(AND) $(MKDIR) -p $* ); $(TOUCH) $@

