CXX ?= g++

RELEASE_OPT = -O2
DEBUG_OPT = -g -O0
PROFILE_OPT = -O2 -pg -g

SRC_PATH = src
BUILD_PATH = pamir-obj
BIN_PATH = pamir
UTIL_PATH = $(BIN_PATH)/util

SCRIPT_SOURCE = scripts
SCRIPT_PATH = $(BIN_PATH)/scripts
SRC_EXT = cc

SOURCE_FILES =  pamir.cc aligner.cc bam_parser.cc common.cc genome.cc logger.cc partition.cc record.cc assembler.cc sam_parser.cc sort.cc extractor.cc 
UTIL_SRC_FILES = extract_support.cc smoother.cc clean_megablast.cc recalibrate.cc
TAMIR_SRC_FILES = process_reads.cc process_orphans.cc process_range.cc edlib.cc sam_processing.cc

SCRIPT_FILES = merge_refs.py contig_graph.py  filter_by_setcover.py  filtering.py  generate_setcover_input.py  genotyping.py  prep-ctgs.py  remove_contaminations.py  sort_vcf.py  version_check.py


SOURCES = $(patsubst %, $(SRC_PATH)/%, $(SOURCE_FILES))
OBJECTS = $(SOURCES:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)


UTIL_SRC = $(patsubst %, $(SRC_PATH)/%, $(UTIL_SRC_FILES))
UTIL_OBJ = $(UTIL_SRC:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)

.SECONDARY: $(UTIL_OBJ)

TAMIR_SRC = $(patsubst %, $(SRC_PATH)/%, $(TAMIR_SRC_FILES))
TAMIR_OBJ = $(TAMIR_SRC:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)


SCRIPTS = $(patsubst %, $(SCRIPT_SOURCE)/%, $(SCRIPT_FILES))
COPIED_SCRIPTS = $(patsubst %, $(SCRIPT_PATH)/%, $(SCRIPT_FILES))


EXE=pamir
TAMIR_EXE=process
UTIL_EXE = $(UTIL_SRC_FILES:%.cc=$(UTIL_PATH)/%)


DEPS = $(OBJECTS:.o=.d) $(UTIL_OBJ:.o=.d) $(TAMIR_OBJ:.o=.d)

COMPILE_FLAGS = -std=c++14 #-Wall -Wextra #Removed because pamir gives way too many warnings 
INCLUDES = -I $(SRC_PATH)/include/

LFLAGS = $(LDFLAGS) -lm -lz

.PHONY: default_make
default_make: release

.PHONY: release
.PHONY: debug
.PHONY: profile
.PHONY: install


release: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(RELEASE_OPT)
release: dirs
	@$(MAKE) all
	
debug: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(DEBUG_OPT)
debug: dirs
	@$(MAKE) all

profile: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(PROFILE_OPT)
profile: dirs
	@$(MAKE) all

install:  export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(RELEASE_OPT)
install: dirs
	@$(MAKE) install_helper
	@$(MAKE) build_clean

.PHONY: dirs
dirs:  
	@mkdir -p $(BUILD_PATH)
	@mkdir -p $(BIN_PATH)
	@mkdir -p $(UTIL_PATH)




clean: build_clean bin_clean

.PHONY: clean

build_clean:
	@$(RM)  $(OBJECTS)
	@$(RM)  $(UTIL_OBJ)
	@$(RM)  $(TAMIR_OBJ)
	@$(RM)  $(DEPS)
	@$(RM) -d $(BUILD_PATH)

.PHONY: build_clean

bin_clean:
	$(RM)  $(UTIL_EXE)
	$(RM)  $(BIN_PATH)/$(EXE)
	$(RM)	$(BIN_PATH)/.snakemake/ -r	
	$(RM)  $(BIN_PATH)/Snakefile
	$(RM)  $(UTIL_PATH)/$(TAMIR_EXE)
	$(RM) -d $(UTIL_PATH)
	$(RM) -dr $(SCRIPT_PATH)
	$(RM) -d $(BIN_PATH)
	
$(COPIED_SCRIPTS): dirs $(SCRIPTS)
	@mkdir -p $(BIN_PATH)/scripts
	@cp scripts/* $(BIN_PATH)/scripts
	@cp Snakefile $(BIN_PATH)

-include $(DEPS)


.PHONY: all
all: $(BIN_PATH)/$(EXE) $(UTIL_EXE) $(UTIL_PATH)/$(TAMIR_EXE) $(COPIED_SCRIPTS)

.PHONY: install_helper

install_helper: $(BIN_PATH)/$(EXE) $(UTIL_EXE) $(UTIL_PATH)/$(TAMIR_EXE) $(COPIED_SCRIPTS)
	@echo -e "\nPlease add `pwd` to your path;\n\nOR\n\nmv pamir.sh /usr/bin\nmv pamir /usr/bin\n"

$(BIN_PATH)/$(EXE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $@

$(UTIL_PATH)/%: $(BUILD_PATH)/%.o $(BUILD_PATH)/common.o
	$(CXX) $(CXXFLAGS) $(LFLAGS) -o $@ $< $(BUILD_PATH)/common.o

$(UTIL_PATH)/$(TAMIR_EXE): $(TAMIR_OBJ)
	$(CXX)  $(LFAGS) -o $@ $(TAMIR_OBJ)

$(BUILD_PATH)/%.o: $(SRC_PATH)/%.$(SRC_EXT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@
