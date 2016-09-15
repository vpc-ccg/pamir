CC=g++
CFLAGS= -c -g -std=gnu++0x -O3  
SOURCES=partition.cc sniper.cc assembler.cc genome.cc aligner.cc  assembler_ext.cc extractor.cc common.cc bam_parser.cc sam_parser.cc record.cc sort.cc
LDFLAGS=-lm -lz
OBJECTS=$(SOURCES:.cc=.o) 
EXECUTABLE=sniper
all: snp pp rc rd es

pp: 
	g++ -std=gnu++0x -O3 -o partition_processor partition_processor.cc common.cc
rc: 
	g++ -O3 -o recalibrate recalibrate.cc
rd: 
	g++ -O3 -o remove_duplicate_insertions remove_duplicate_insertions.cc
es: 
	g++ -O3 -o extract_support extract_support.cc common.cc
snp: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm -f *.o
	rm -f sniper
	rm -f partition_processor
	rm -f recalibrate
	rm -f remove_duplicate_insertions
	rm -f extract_support
