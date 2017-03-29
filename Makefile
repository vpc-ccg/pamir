CC=g++
FLAGS= -O3 -std=c++11 
CFLAGS= -c $(FLAGS) -Wfatal-errors
SOURCES=partition.cpp pamir.cpp assembler.cpp genome.cpp aligner.cpp extractor.cpp common.cpp bam_parser.cpp sam_parser.cpp record.cpp sort.cpp logger.cpp
LDFLAGS=-lm -lz
OBJECTS=$(SOURCES:.cc=.o) 
EXECUTABLE=pamir
all: snp rc es sm mf cm
basic: snp rc es mf cm

rc: 
	g++ -O3 -o recalibrate recalibrate.cpp
cm: 
	g++ -O3 -o clean clean_megablast.cpp
es: 
	g++ -O3 -o extract_support extract_support.cpp common.cpp
sm:
	g++ -O3 -o smoother -g -std=c++1y -Wfatal-errors smoother.cpp
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
	rm -f recalibrate
	rm -f extract_support
	rm -f smoother
	rm -f clean
