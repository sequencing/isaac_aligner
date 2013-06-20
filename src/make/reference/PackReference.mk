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
## file PackReference.mk
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

ifeq (,$(REFERENCE_GENOME))
$(error "REFERENCE_GENOME is not defined")
endif

ifeq (,$(OUTPUT_FILE))
$(error "OUTPUT_FILE is not defined")
endif

GENOME_NEIGHBORS_PREFIX:=$(TEMP_DIR)/genome-neighbors-
GENOME_NEIGHBORS_SUFFIX:=mer.dat
GENOME_FASTA:=$(TEMP_DIR)/genome.fa
SORTED_REFERENCE_XML:=$(TEMP_DIR)/sorted-reference.xml

SUPPORTED_SEED_LENGTHS:=$(shell xsltproc $(GET_SUPPORTED_SEED_LENGTHS_XSL) $(REFERENCE_GENOME))

ALL_GENOME_NEIGHBORS:=$(foreach s, $(SUPPORTED_SEED_LENGTHS), $(GENOME_NEIGHBORS_PREFIX)$(s)$(GENOME_NEIGHBORS_SUFFIX))

genome_neighbors_seed_length=$(@:$(GENOME_NEIGHBORS_PREFIX)%$(GENOME_NEIGHBORS_SUFFIX)=%)

$(GENOME_FASTA): $(REFERENCE_GENOME) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(REORDER_REFERENCE) --reference-genome $< --output-xml $(SORTED_REFERENCE_XML).tmp --output-fasta $(SAFEPIPETARGET) &&\
	$(MV) $(SORTED_REFERENCE_XML).tmp $(SORTED_REFERENCE_XML)

$(ALL_GENOME_NEIGHBORS): $(REFERENCE_GENOME) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(EXTRACT_NEIGHBORS) --reference-genome $< --seed-length $(genome_neighbors_seed_length) --output-file $(SAFEPIPETARGET)

$(OUTPUT_FILE): $(ALL_GENOME_NEIGHBORS) $(GENOME_FASTA)
	$(CMDPREFIX) $(TAR) -czvO -C $(TEMP_DIR) \
		$(foreach gn, $(ALL_GENOME_NEIGHBORS), $(notdir $(gn)))  \
		$(notdir $(GENOME_FASTA)) \
		$(notdir $(SORTED_REFERENCE_XML)) >$(SAFEPIPETARGET) 

all: $(OUTPUT_FILE)

