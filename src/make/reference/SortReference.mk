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
## file SortReference.mk
##
## brief Defines appropriate rules
##
## author Roman Petrovski
##
################################################################################

# first target needs to be defined in the beginning. Ohterwise includes such as
# Log.mk cause unexpected behavior
firsttarget: all

MAKEFILES_DIR:=@iSAAC_FULL_DATADIR@/makefiles

# Import the global configuration
include $(MAKEFILES_DIR)/common/Config.mk

include $(MAKEFILES_DIR)/common/Sentinel.mk

# Import the logging functionalities
include $(MAKEFILES_DIR)/common/Log.mk

# Import the debug functionalities
include $(MAKEFILES_DIR)/common/Debug.mk

include $(MAKEFILES_DIR)/reference/Config.mk

include config.mk

ifeq (,$(MASK_WIDTH))
$(error "MASK_WIDTH is not defined")
endif

ifeq (,$(GENOME_FILE))
$(error "GENOME_FILE is not defined")
endif


MASK_COUNT:=$(shell $(AWK) 'BEGIN{print 2^$(MASK_WIDTH)}')
MASK_LIST:=$(wordlist 1, $(MASK_COUNT), $(shell $(SEQ) --equal-width 0 $(MASK_COUNT)))

GENOME_NAME:=$(if $(GENOME_NAME),$(GENOME_NAME),$(notdir $(GENOME_FILE)))

MASK_FILE_PREFIX:=$(GENOME_NAME)-32mer-$(MASK_WIDTH)bit-
SORTED_REFERENCE_XML:=sorted-reference.xml
CONTIGS_XML:=contigs.xml

permutation=$(word 1,$(subst -, ,$(@:$(TEMP_DIR)/$(MASK_FILE_PREFIX)%$(MASK_FILE_XML_SUFFIX)=%)))
mask=$(word 2,$(subst -, ,$(@:$(TEMP_DIR)/$(MASK_FILE_PREFIX)%$(MASK_FILE_XML_SUFFIX)=%)))
mask_file=$(@:$(TEMP_DIR)/%$(MASK_FILE_XML_SUFFIX)=%$(MASK_FILE_SUFFIX))

ALL_MASK_XMLS:=$(foreach p, $(PERMUTATION_NAME_LIST), $(foreach m, $(MASK_LIST), $(TEMP_DIR)/$(MASK_FILE_PREFIX)$(p)-$(m)$(MASK_FILE_XML_SUFFIX)))
ALL_MASKS:=$(foreach p, $(PERMUTATION_NAME_LIST), $(foreach m, $(MASK_LIST), $(MASK_FILE_PREFIX)$(p)-$(m)$(MASK_FILE_SUFFIX)))
ALL_MASKS_TMP:=$(foreach p, $(PERMUTATION_NAME_LIST), $(foreach m, $(MASK_LIST), $(MASK_FILE_PREFIX)$(p)-$(m)$(MASK_TMP_FILE_SUFFIX)))

$(ALL_MASK_XMLS): $(GENOME_FILE) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(SORT_REFERENCE) -g $(GENOME_FILE) --mask-width $(MASK_WIDTH) --mask $(mask) \
		--permutation-name $(permutation) --output-file $(TEMP_DIR)/$(mask_file) --repeat-threshold $(REPEAT_THRESHOLD) >$(SAFEPIPETARGET)

$(TEMP_DIR)/$(CONTIGS_XML): $(GENOME_FILE) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(PRINT_CONTIGS) -g $(GENOME_FILE) >$(SAFEPIPETARGET)

# The order of prerequisites matters. We want the Contigs to come before Permutations in merged file
$(TEMP_DIR)/$(SORTED_REFERENCE_XML): $(TEMP_DIR)/$(CONTIGS_XML) $(ALL_MASK_XMLS)
	$(CMDPREFIX) $(ECHO) '<?xml version="1.0"?><SortedReference/>' $(foreach part, $^, \
		| $(XSLTPROC) --param with "'$(part)'" $(MERGE_XML_DOCUMENTS_XSL) -) \
	> $(SAFEPIPETARGET)

$(SORTED_REFERENCE_XML):$(TEMP_DIR)/$(SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(FIND_NEIGHBORS) -i $< -t $(TEMP_DIR)/neighbors.dat --output-directory $(CURDIR) $(FIND_NEIGHBORS_OPTIONS) -o $(SAFEPIPETARGET)

all: $(SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(LOG_INFO) "All done!"

.PHONY: clean_intermediary.done
clean_intermediary.done:
	$(CMDPREFIX) $(RM) $(ALL_MASK_XMLS) $(ALL_MASKS_TMP) $(TEMP_DIR)/$(CONTIGS_XML) $(TEMP_DIR)/$(SORTED_REFERENCE_XML)

.PHONY: clean
clean: clean_intermediary.done
	$(CMDPREFIX) $(RM) $(ALL_MASKS) $(SORTED_REFERENCE_XML)
