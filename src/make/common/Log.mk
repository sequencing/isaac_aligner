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
## brief Partial makefile providing basic logging support
##
## - Redefines the SHELL
##
## author Roman Petrovski
##
################################################################################

SHELL_LOG_OLD := $(SHELL)

# Important: any command that has word qmake in it is executed locally on machine that is
# running qmake (that's how qmake has coded it in). Make sure below that if qmake is the
# current make, the word is renamed before it is passed to the shell.
# Important: some targets have a very big number of prerequisites. Using $(wordlist... to
# avoid 'make: execvp: .../loggingshell.sh: Argument list too long'

# in make -n mode use debug sheell to be able to see the information on the prerequisites 
ifeq (n,$(filter n,$(MAKEFLAGS)))
override SHELL = $(if $@,$(warning [$@ ($?)]))$(SHELL_LOG_OLD)
override SHELL_LOG = $(if $@,$(warning [$@ ($?)]))$(SHELL_LOG_OLD)
else

ifneq (,$(JOB_NAME))

NAME_SHELL_FILE:=$(JOB_NAME)
# make sure it does not contain word qmake or else everyting will be done locally to qmake node 
NAME_SHELL_FILE:=$(subst qmake,q-make,$(NAME_SHELL_FILE))
# Don't put shell in $(TEMP_DIR). It gets cleaned up!
# Make sure the shell path is absolute.
NAME_SHELL_FOLDER:=$(CURDIR)/NameShell
NAME_SHELL:=$(NAME_SHELL_FOLDER)/$(NAME_SHELL_FILE)

# SGE uses the name of shell script for recipe job names. Generate one and use it
# so that user can do something like:
#  qmake -now no -cwd -v PATH -N blah3  -- -j 48 all
# note: the #$ -N sge trick in the job script is flaky on the sge side. Some recipes see it
# and some not. this results in different recipe job names for the same workflow
$(NAME_SHELL): SHELL:=$(SHELL)
$(NAME_SHELL): $(dir $(NAME_SHELL)).sentinel
	echo -e 'script=$$1\nshift\nexec $$script "$$@"' >$@.tmp \
	$(AND) $(CHMOD) 755 $(SAFEPIPETARGET)

# depend this makefile of name shell so that the rule gets executed before anything
$(MAKEFILES_DIR)/common/Log.mk: $(NAME_SHELL)
endif


LOG_SHELL:= $(NAME_SHELL) $(LIBEXEC_DIR)/loggingShell.sh
override SHELL = $(if $@,$(LOG_SHELL) '$(iSAAC_LOG_LEVEL)' '$(LOG_DATE_FORMAT)' '$(subst qmake,q-make,$(MAKE))' '$@' '$(wordlist 1, 50, $^)' '$(wordlist 1, 50, $?)' '$(SHELL_LOG_OLD)',$(SHELL_LOG_OLD))
override SHELL_LOG = $(if $@,$(LOG_SHELL) '$(iSAAC_LOG_LEVEL)' '$(LOG_DATE_FORMAT)' '$(subst qmake,q-make,$(MAKE))' '$@' '$(wordlist 1, 50, $^)' '$(wordlist 1, 50, $?)' '$(SHELL_LOG_OLD)',$(SHELL_LOG_OLD))
endif
