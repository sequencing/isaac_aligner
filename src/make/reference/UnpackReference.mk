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
## file UnpackReference.mk
##
## brief Defines appropriate rules
##
## author Roman Petrovski
##
################################################################################

# first target needs to be defined in the beginning. Ohterwise includes such as
# Log.mk cause unexpected behavior
.PHONY: firsttarget
firsttarget: all

THIS_MAKEFILE:=$(word $(words $(MAKEFILE_LIST)), $(MAKEFILE_LIST))


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

ifeq (,$(MASK_WIDTH))
$(error "MASK_WIDTH is not defined")
endif

ifeq (,$(INPUT_ARCHIVE))
$(error "INPUT_ARCHIVE is not defined")
endif

UNPACKED_SORTED_REFERENCE_XML:=$(TEMP_DIR)/sorted-reference.xml

UNPACKED_GENOME_FILE:=$(TEMP_DIR)/genome.fa
GENOME_FILE:=$(CURDIR)/genome.fa


MASK_COUNT:=$(shell $(AWK) 'BEGIN{print 2^$(MASK_WIDTH)}')
MASK_LIST:=$(wordlist 1, $(MASK_COUNT), $(shell $(SEQ) --equal-width 0 $(MASK_COUNT)))

GENOME_NAME:=genome.fa

GENOME_NEIGHBORS_PREFIX:=genome-neighbors-
GENOME_NEIGHBORS_SUFFIX:=mer.dat
GENOME_NEIGHBORS_PATTERN:=$(GENOME_NEIGHBORS_PREFIX)%$(GENOME_NEIGHBORS_SUFFIX)

HIGH_REPEATS_PREFIX:=repeats-$(REPEAT_THRESHOLD)-
HIGH_REPEATS_SUFFIX:=mer.1bpb
HIGH_REPEATS_PATTERN:=$(HIGH_REPEATS_PREFIX)%$(HIGH_REPEATS_SUFFIX)

MASK_FILE_PREFIX:=$(GENOME_NAME)-
MASK_FILE_MIDDLE:=mer-$(MASK_WIDTH)bit-

SORTED_REFERENCE_XML:=sorted-reference.xml
CONTIGS_XML:=contigs.xml

mask=$(word 2,$(subst $(MASK_FILE_MIDDLE), ,$(@:$(TEMP_DIR)/$(MASK_FILE_PREFIX)%$(MASK_FILE_XML_SUFFIX)=%)))
seed_length=$(word 1,$(subst $(MASK_FILE_MIDDLE), ,$(@:$(TEMP_DIR)/$(MASK_FILE_PREFIX)%$(MASK_FILE_XML_SUFFIX)=%)))

get_format_version=$(shell $(XSLTPROC) $(GET_FORMAT_VERSION_XSL) $(UNPACKED_SORTED_REFERENCE_XML))
get_seed_length_list=$(shell $(XSLTPROC) $(GET_SUPPORTED_SEED_LENGTHS_XSL) $(UNPACKED_SORTED_REFERENCE_XML))

get_all_genome_neighbors=$(foreach s, $(get_seed_length_list), $(GENOME_NEIGHBORS_PREFIX)$(s)$(GENOME_NEIGHBORS_SUFFIX))
get_all_high_repeats=$(foreach s, $(get_seed_length_list), $(HIGH_REPEATS_PREFIX)$(s)$(HIGH_REPEATS_SUFFIX))

ifneq (,$(UNPACK_REFERENCE_SUBMAKE))

ifneq ($(CURRENT_REFERENCE_FORMAT_VERSION),$(get_format_version))
$(error "Unsupported packed reference format $(get_format_version). Version $(CURRENT_REFERENCE_FORMAT_VERSION) is required")
endif

ALL_MASK_XMLS:=$(foreach sl, $(get_seed_length_list), $(foreach m, $(MASK_LIST), $(TEMP_DIR)/$(MASK_FILE_PREFIX)$(sl)$(MASK_FILE_MIDDLE)$(m)$(MASK_FILE_XML_SUFFIX)))
ALL_MASKS:=$(foreach sl, $(get_seed_length_list), $(foreach m, $(MASK_LIST), $(TEMP_DIR)/$(MASK_FILE_PREFIX)$(sl)$(MASK_FILE_MIDDLE)$(m)$(MASK_FILE_SUFFIX)))
$(ALL_MASK_XMLS): $(GENOME_FILE)
	$(CMDPREFIX) $(SORT_REFERENCE) -g $(GENOME_FILE) --mask-width $(MASK_WIDTH) --mask $(mask) \
		--output-file $(notdir $(@:%$(MASK_FILE_XML_SUFFIX)=%$(MASK_FILE_SUFFIX))) \
		--repeat-threshold $(REPEAT_THRESHOLD) \
		--seed-length $(seed_length) \
		--genome-neighbors $(TEMP_DIR)/$(GENOME_NEIGHBORS_PREFIX)$(seed_length)$(GENOME_NEIGHBORS_SUFFIX) >$(SAFEPIPETARGET)

$(TEMP_DIR)/$(CONTIGS_XML): $(GENOME_FILE) $(TEMP_DIR)/.sentinel $(UNPACKED_SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(PRINT_CONTIGS) -g $(GENOME_FILE) --original-metadata $(UNPACKED_SORTED_REFERENCE_XML) >$(SAFEPIPETARGET)

$(SORTED_REFERENCE_XML): $(TEMP_DIR)/$(CONTIGS_XML) $(ALL_MASK_XMLS)
	$(CMDPREFIX) $(MERGE_REFERENCES) $(foreach part, $^, -i '$(part)') -o $(SAFEPIPETARGET)

$(GENOME_NEIGHBORS_PATTERN) $(HIGH_REPEATS_PATTERN): $(SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(EXTRACT_NEIGHBORS) --reference-genome $< \
	--seed-length $* \
	--output-file $(GENOME_NEIGHBORS_PREFIX)$*$(GENOME_NEIGHBORS_SUFFIX).tmp --high-repeats-file $(HIGH_REPEATS_PREFIX)$*$(HIGH_REPEATS_SUFFIX).tmp && \
	$(MV) $(GENOME_NEIGHBORS_PREFIX)$*$(GENOME_NEIGHBORS_SUFFIX).tmp $(GENOME_NEIGHBORS_PREFIX)$*$(GENOME_NEIGHBORS_SUFFIX) && \
	$(MV) $(HIGH_REPEATS_PREFIX)$*$(HIGH_REPEATS_SUFFIX).tmp $(HIGH_REPEATS_PREFIX)$*$(HIGH_REPEATS_SUFFIX)
endif


$(UNPACKED_SORTED_REFERENCE_XML): $(INPUT_ARCHIVE) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(TAR) -C $(TEMP_DIR) -xvf $(INPUT_ARCHIVE) |xargs -l -I blah $(TOUCH) $(TEMP_DIR)/blah && \
	$(MV) $(UNPACKED_GENOME_FILE) $(GENOME_FILE)

.PHONY: all
all: SHELL:=$(SHELL_LOG_OLD)
all: $(UNPACKED_SORTED_REFERENCE_XML)
	$(MAKE) -f $(THIS_MAKEFILE) $(SORTED_REFERENCE_XML) $(get_all_genome_neighbors) $(get_all_high_repeats) UNPACK_REFERENCE_SUBMAKE:=yes && \
	$(CMDPREFIX) $(LOG_INFO) "All done!"


