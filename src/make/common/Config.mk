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
## file Config.mk
##
## brief Common configuration file for all makefiles.
##
## Defines paths, file names, extensions, etc.
##
## author Roman Petrovski
##
################################################################################

# The shell pipefail option is required for accurate workflow control
SHELL = /bin/bash -o pipefail

# iSAAC installation directories
TOOLS_DIR:=@iSAAC_FULL_BINDIR@
BIN_DIR:=@iSAAC_FULL_BINDIR@
DATA_DIR:=@iSAAC_FULL_DATADIR@
LIBEXEC_DIR:=@iSAAC_FULL_LIBEXECDIR@

# System tools
AWK = awk
CAT = cat
CD = cd
CHMOD = chmod
# install should be used instead of cp to have the proper permission bits on the target files.
#CP = cp
#CPDIR = cp -R
CUT = cut
DATE = date
DIFF = diff
ECHO = echo
EGREP = egrep
EXIT = exit
FGREP = fgrep
FIND = find
GREP = grep
GUNZIP = gunzip
MKDIR = mkdir
MV = mv
PASTE = paste
PRINTF = printf
PWD = pwd
RM = rm -f
SED = sed
SET=set
SEQ=seq
SORT = sort -T $(SORT_TMPDIR)
TAR = tar
TEE = tee
TOUCH = touch
UNIQ = uniq
XARGS = xargs
XSLTPROC = xsltproc
WC = wc
SLEEP=sleep

AND=&&
OR=||

# Logging

# level 0 - no logging except for error messages. All recipe stderr is suppressed
# level 1 - only text preffixed with known keywords such as (INFO, ERROR, WARNING and similar) is allowed
#           to escape recipes stderr
# level 2 - unfiltered stderr
iSAAC_LOG_LEVEL:=2

LOG_DATE_FORMAT:=%F %T
LOG_DEBUG=[[ $(iSAAC_LOG_LEVEL) != 2 ]] $(OR) 1>&2 $(ECHO)
LOG_WARNING=[[ $(iSAAC_LOG_LEVEL) == 0 ]] $(OR) 1>&2 $(ECHO) -e "WARNING:"
LOG_INFO=[[ $(iSAAC_LOG_LEVEL) == 0 ]] $(OR) 1>&2 $(ECHO) -e "INFO:"

#Condition checking
BOOLEAN_TRUE_WORDS:=Y y yes YES Yes on ON On true TRUE True 1 ok OK Ok
BOOLEAN_FALSE_WORDS:=N n no NO No off OFF Off false FALSE False 0
CHECK_VALID_BOOLEAN=$(if $(filter $(BOOLEAN_TRUE_WORDS) $(BOOLEAN_FALSE_WORDS), $(1)),$(1),\
                                $(error "Incorrect boolean value '$(1)'. Allowed values: $(BOOLEAN_TRUE_WORDS) $(BOOLEAN_FALSE_WORDS)"))
IS_BOOLEAN_TRUE=$(if $(strip $(1)),$(filter $(BOOLEAN_TRUE_WORDS), $(call CHECK_VALID_BOOLEAN,$(1))))

# Global macros
UNPROTECTED_TARGET=$@
SAFEPIPETARGET = $@.tmp && mv $@.tmp $@
SAFEPIPETARGET2=$(subst .,.tmp.,$@) && mv $(subst .,.tmp.,$@) $@
SAFEPIPETARGET3=$(CURDIR)/$@.tmp && mv $(CURDIR)/$@.tmp $(CURDIR)/$@
TOUCH_TARGET=$(TOUCH) $@
CHECK_TARGET=@if [ ! -e $@ ]; then echo "Error: $@ does not exist."; exit 1; fi

# Structure of the analysis folder
TEMP_DIR:=Temp

#########################################################
# tools

# ... findNeighbors
EXTRACT_NEIGHBORS:=$(LIBEXEC_DIR)/extractNeighbors
# ...

# ... findNeighbors
FIND_NEIGHBORS:=$(LIBEXEC_DIR)/findNeighbors
# ...

# ... printContigs
PRINT_CONTIGS:=$(LIBEXEC_DIR)/printContigs
# ...

# ... isaac-reorder-reference
REORDER_REFERENCE:=$(BIN_DIR)/isaac-reorder-reference
# ...

# ... sortReference
SORT_REFERENCE:=$(LIBEXEC_DIR)/sortReference
# ...

MERGE_XML_DOCUMENTS_XSL:=$(DATA_DIR)/xsl/common/MergeXmlDocuments.xsl 
