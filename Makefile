CC=g++
CFLAGS= -c -g -std=gnu++0x -I ../dz 
SOURCES=partition.cc sniper.cc assembler.cc genome.cc aligner.cc logger.cc assembler_ext.cc extractor.cc 
LDFLAGS=-lm -L ../zlib-1.2.8/ -lz -lpthread -pthread ../dz/libdeez.a
#LDFLAGS=-lm -lz -lpthread -pthread
OBJECTS=$(SOURCES:.cc=.o) 
EXECUTABLE=sniper
all: snp pp rc rd
pp: 
	g++ -O3 -o partition_processor partition_processor.cc
rc: 
	g++ -O3 -o recalibrate recalibrate.cc
rd: 
	g++ -O3 -o remove_duplicate_insertions remove_duplicate_insertions.cc

snp: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	#$(CC) $(OBJECTS) -lz -pthread -o $@
	$(CC) $(OBJECTS) $(LDFLAGS) ../dz/libdeez.a -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@ -Idz

dzlib:
	make -B -j lib -C ../dz

clean:
	rm -f *.o
	rm -f sniper
	rm -f partition_processor
	rm -f recalibrate
	rm -f remove_duplicate_insertions
