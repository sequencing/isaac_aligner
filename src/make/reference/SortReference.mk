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

ifeq (,$(ISAAC_HOME))
MAKEFILES_DIR:=@iSAAC_FULL_DATADIR@/makefiles
else
MAKEFILES_DIR:=$(ISAAC_HOME)/@iSAAC_PARTIAL_DATADIR@/makefiles
endif

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

MASK_FILE_PREFIX:=$(GENOME_NAME)-$(SEED_LENGTH)mer-$(MASK_WIDTH)bit-
SORTED_REFERENCE_XML:=sorted-reference.xml
CONTIGS_XML:=contigs.xml
# actual kmers that have neighbors in the genome 
NEIGHBORS_DAT:=neighbors.dat
GENOME_NEIGHBORS_DAT:=genome-neighbors.1bpb
GENOME_NEIGHBORS_DAT_PATTERN:=genome-neighbors%1bpb
HIGH_REPEATS_DAT:=repeats-$(REPEAT_THRESHOLD).1bpb
HIGH_REPEATS_DAT_PATTERN:=repeats-$(REPEAT_THRESHOLD)%1bpb

mask=$(word 1,$(subst -, ,$(@:$(TEMP_DIR)/$(MASK_FILE_PREFIX)%$(MASK_FILE_XML_SUFFIX)=%)))
mask_file=$(@:$(TEMP_DIR)/%$(MASK_FILE_XML_SUFFIX)=%$(MASK_FILE_SUFFIX))

ALL_MASK_XMLS:=$(foreach m, $(MASK_LIST), $(TEMP_DIR)/$(MASK_FILE_PREFIX)$(m)$(MASK_FILE_XML_SUFFIX))
ALL_MASKS:=$(foreach m, $(MASK_LIST), $(MASK_FILE_PREFIX)$(m)$(MASK_FILE_SUFFIX))
ALL_MASKS_TMP:=$(foreach m, $(MASK_LIST), $(MASK_FILE_PREFIX)$(m)$(MASK_TMP_FILE_SUFFIX))

$(ALL_MASK_XMLS): $(GENOME_FILE) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(SORT_REFERENCE) -g $(GENOME_FILE) --mask-width $(MASK_WIDTH) --mask $(mask) \
		--seed-length $(SEED_LENGTH) \
		--output-file $(TEMP_DIR)/$(mask_file) \
		--repeat-threshold $(REPEAT_THRESHOLD) >$(SAFEPIPETARGET)

$(TEMP_DIR)/$(CONTIGS_XML): $(GENOME_FILE) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(PRINT_CONTIGS) -g $(GENOME_FILE) >$(SAFEPIPETARGET)

$(TEMP_DIR)/$(SORTED_REFERENCE_XML): $(TEMP_DIR)/$(CONTIGS_XML) $(ALL_MASK_XMLS)
	$(CMDPREFIX) $(MERGE_REFERENCES) $(foreach part, $^, -i '$(part)') -o $(SAFEPIPETARGET)

ifeq (false,$(DONT_ANNOTATE))
$(SORTED_REFERENCE_XML):$(TEMP_DIR)/$(SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(FIND_NEIGHBORS) -i $< -t $(TEMP_DIR)/$(NEIGHBORS_DAT) \
		--seed-length $(SEED_LENGTH) \
		--output-directory $(CURDIR) $(FIND_NEIGHBORS_OPTIONS) -o $(SAFEPIPETARGET)

$(GENOME_NEIGHBORS_DAT_PATTERN) $(HIGH_REPEATS_DAT_PATTERN): $(SORTED_REFERENCE_XML) $(ALL_MASK_XMLS)
	$(CMDPREFIX) $(EXTRACT_NEIGHBORS) --reference-genome $< \
		--seed-length $(SEED_LENGTH) \
		--output-file $(GENOME_NEIGHBORS_DAT).tmp --high-repeats-file $(HIGH_REPEATS_DAT).tmp && \
	$(MV) $(GENOME_NEIGHBORS_DAT).tmp $(GENOME_NEIGHBORS_DAT) && \
	$(MV) $(HIGH_REPEATS_DAT).tmp $(HIGH_REPEATS_DAT)

all: $(SORTED_REFERENCE_XML) $(GENOME_NEIGHBORS_DAT) $(HIGH_REPEATS_DAT)
	$(CMDPREFIX) $(LOG_INFO) "All done!"
else
$(SORTED_REFERENCE_XML):$(TEMP_DIR)/$(SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(CP) $< $(SAFEPIPETARGET)

all: $(SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(LOG_INFO) "All done!"
endif

