CXX ?= g++

RELEASE_OPT = -O2
DEBUG_OPT = -g -O0
PROFILE_OPT = -O2 -pg -g

EXT_PATH = ext
SRC_PATH = src
BUILD_PATH = pamir-obj
BIN_PATH = pamir


LOGGER_EXT_PATH = $(EXT_PATH)/util-logger/include/logger.h
LOGGER_HED_PATH  = $(SRC_PATH)/include/logger.h


EDLIB_EXT_HED_PATH = $(EXT_PATH)/edlib/edlib/include/edlib.h
EDLIB_EXT_SRC_PATH = $(EXT_PATH)/edlib/edlib/src/edlib.cpp

EDLIB_SRC_PATH = $(SRC_PATH)/edlib.cc
EDLIB_HED_PATH = $(SRC_PATH)/include/edlib.h

SCRIPT_SOURCE = scripts
SCRIPT_PATH = $(BIN_PATH)/scripts
SRC_EXT = cc

SOURCE_FILES =  pamir.cc aligner.cc common.cc genome.cc partition.cc assembler.cc sam_parser.cc extractor.cc  smoother.cc  process_orphans.cc process_range.cc edlib.cc sam_processing.cc recalibrate.cc


SCRIPT_FILES = merge_refs.py contig_graph.py  filter_by_setcover.py  filtering.py  generate_setcover_input.py  prep-ctgs.py  remove_contaminations.py  sort_vcf.py  version_check.py process_repeats.py process_unique.py new_generate_setcover.py

SOURCES = $(patsubst %, $(SRC_PATH)/%, $(SOURCE_FILES))
OBJECTS = $(SOURCES:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)

.SECONDARY: $(UTIL_OBJ)

SCRIPTS = $(patsubst %, $(SCRIPT_SOURCE)/%, $(SCRIPT_FILES))
COPIED_SCRIPTS = $(patsubst %, $(SCRIPT_PATH)/%, $(SCRIPT_FILES))

EXE=pamir

DEPS = $(OBJECTS:.o=.d) $(UTIL_OBJ:.o=.d) $(PROCESSING_OBJ:.o=.d)

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
#	@$(MAKE) build_clean

.PHONY: dirs
dirs:  
	@mkdir -p $(BUILD_PATH)
	@mkdir -p $(BIN_PATH)

clean: build_clean bin_clean

.PHONY: clean

build_clean:
	@$(RM)  $(OBJECTS)
	@$(RM)  $(DEPS)
	@$(RM) -d $(BUILD_PATH)
	@$(RM) $(SRC_PATH)/edlib.cc
	@$(RM) $(SRC_PATH)/include/edlib.h
	@$(RM) $(SRC_PATH)/include/logger.h
.PHONY: build_clean

bin_clean:
	$(RM)  $(BIN_PATH)/$(EXE)
	$(RM)	$(BIN_PATH)/.snakemake/ -r	
	$(RM)  $(BIN_PATH)/Snakefile
	$(RM)  $(BIN_PATH)/version
	$(RM)  $(BIN_PATH)/environment.yaml
	$(RM) -dr $(SCRIPT_PATH)
	$(RM) -d $(BIN_PATH)
	
$(COPIED_SCRIPTS): dirs $(SCRIPTS)
	@mkdir -p $(BIN_PATH)/scripts
	@cp scripts/* $(BIN_PATH)/scripts
	@cp Snakefile $(BIN_PATH)
	@cp environment.yaml $(BIN_PATH)
-include $(DEPS)

$(LOGGER_EXT_PATH): 
	@echo Please clone the repository with --recursive option!; exit 1;

$(EDLIB_EXT_SRC_PATH):
	@echo Please clone the repository with --recursive option!; exit 1;

$(EDLIB_EXT_HED_PATH):
	@echo Please clone the repository with --recursive option!; exit 1;

$(LOGGER_HED_PATH): $(LOGGER_EXT_PATH)
	@cp $(LOGGER_EXT_PATH) $(LOGGER_HED_PATH)

$(EDLIB_SRC_PATH): $(EDLIB_EXT_SRC_PATH)
	@cp $(EDLIB_EXT_SRC_PATH) $(EDLIB_SRC_PATH)

$(EDLIB_HED_PATH): $(EDLIB_EXT_HED_PATH)
	@cp $(EDLIB_EXT_HED_PATH) $(EDLIB_HED_PATH)


.PHONY: all
all: $(EDLIB_HED_PATH) $(EDLIB_SRC_PATH) $(LOGGER_HED_PATH) $(BIN_PATH)/$(EXE)  $(COPIED_SCRIPTS)
	@git rev-parse  --short HEAD > $(BIN_PATH)/version
.PHONY: install_helper

install_helper: $(BIN_PATH)/$(EXE) $(COPIED_SCRIPTS)
	@printf "\nPlease add pamir to your path:\n\nexport PATH=\044PATH:`pwd`\n\nOR\n\nsudo mv pamir.sh /usr/bin\nsudo mv pamir /usr/bin\n\n"

$(BIN_PATH)/$(EXE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $@

$(BUILD_PATH)/%.o: $(SRC_PATH)/%.$(SRC_EXT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@
