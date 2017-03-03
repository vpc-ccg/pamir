CC=g++
FLAGS= -O3 -std=c++11 
CFLAGS= -c $(FLAGS) -Wfatal-errors
SOURCES=partition.cc pamir.cc assembler.cc genome.cc aligner.cc  assembler_ext.cc extractor.cc common.cc bam_parser.cc sam_parser.cc record.cc sort.cc logger.cc
LDFLAGS=-lm -lz
OBJECTS=$(SOURCES:.cc=.o) 
EXECUTABLE=pamir
all: snp pp rc es sm mf cm
basic: snp pp rc es mf cm

pp: 
	g++ -std=gnu++0x -O3 -o partition_processor partition_processor.cc common.cc
rc: 
	g++ -O3 -o recalibrate recalibrate.cc
cm: 
	g++ -O3 -o clean clean_megablast.cc
es: 
	g++ -O3 -o extract_support extract_support.cc common.cc
sm:
	g++ -O3 -o smoother -g -std=c++1y -Wfatal-errors -I$(BOOST_INCLUDE) smoother.cc
mf: 
	make -C mrsfast

snp: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm -f *.o
	rm -f pamir
	rm -f partition_processor
	rm -f recalibrate
	rm -f extract_support
	rm -f smoother
	rm -f clean
