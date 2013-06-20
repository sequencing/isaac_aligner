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
## <https://github.com/sequencing/licenses/>.
##
## The distribution includes the code libraries listed below in the
## 'redist' sub-directory. These are distributed according to the
## licensing terms governing each library.
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

