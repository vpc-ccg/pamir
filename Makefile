CXX=g++
FLAGS= -O3 -std=c++11 
CFLAGS= -c $(FLAGS) -Wfatal-errors
SOURCE_FILES=partition.cc pamir.cc assembler.cc genome.cc aligner.cc extractor.cc common.cc bam_parser.cc sam_parser.cc record.cc sort.cc logger.cc



EXT = cc

LDFLAGS=-lm -lz
SRC_PATH = src
INCLUDES = -I $(SRC_PATH)/include

OBJ_PATH = obj
SOURCES = $(patsubst %, $(SRC_PATH)/%, $(SOURCE_FILES))
OBJECTS=$(SOURCES:$(SRC_PATH)/%.cc=$(OBJ_PATH)/%.o) 
EXECUTABLE=pamir
all: dirs snp rc es sm mf cm
basic: dirs snp rc es mf cm


rc: 
	g++ -O3 $(INCLUDES)  -o recalibrate $(SRC_PATH)/recalibrate.cc
cm: 
	g++ -O3 $(INCLUDES)  -o cleanmega $(SRC_PATH)/clean_megablast.cc
es: 
	g++ -O3 $(INCLUDES)  -o extract_support $(SRC_PATH)/extract_support.cc $(SRC_PATH)/common.cc
sm:
	g++ -O3 $(INCLUDES)  -o smoother -g -std=c++1y -Wfatal-errors $(SRC_PATH)/smoother.cc
mf: 
	make -C mrsfast

dirs:
	mkdir -p $(OBJ_PATH)

snp: $(SOURCES) $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(OBJECTS) $(INCLUDES) $(LDFLAGS) -o $@


$(OBJ_PATH)/%.o: $(SRC_PATH)/%.$(EXT)
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@ 

clean:
	rm -f obj/*.o pamir recalibrate extract_support smoother cleanmega
	rm -d $(OBJ_PATH)

.PHONY: rc

.PHONY: cm
.PHONY: es
.PHONY: sm
.PHONY: mf
.PHONY: dirs
.PHONY: snp
.PHONY: all
.PHONY: basic
.PHONY: clean
