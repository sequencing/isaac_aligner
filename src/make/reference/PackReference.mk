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

MAKEFILES_DIR:=@iSAAC_FULL_DATADIR@/makefiles

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

GENOME_NEIGHBORS_DAT:=$(TEMP_DIR)/genome-neighbors.dat
GENOME_FASTA:=$(TEMP_DIR)/genome.fa
SORTED_REFERENCE_XML:=$(TEMP_DIR)/sorted-reference.xml

$(GENOME_FASTA): $(REFERENCE_GENOME) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(REORDER_REFERENCE) --reference-genome $< --output-xml $(SORTED_REFERENCE_XML).tmp --output-fasta $(SAFEPIPETARGET) &&\
	$(MV) $(SORTED_REFERENCE_XML).tmp $(SORTED_REFERENCE_XML)

$(GENOME_NEIGHBORS_DAT): $(REFERENCE_GENOME) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(EXTRACT_NEIGHBORS) --reference-genome $< --output-file $(SAFEPIPETARGET)

$(OUTPUT_FILE): $(GENOME_NEIGHBORS_DAT) $(GENOME_FASTA)
	$(CMDPREFIX) $(TAR) -czvO \
		-C $(dir $(GENOME_NEIGHBORS_DAT)) $(notdir $(GENOME_NEIGHBORS_DAT))  \
		$(notdir $(GENOME_FASTA)) \
		$(notdir $(SORTED_REFERENCE_XML)) >$(SAFEPIPETARGET) 

all: $(OUTPUT_FILE)

