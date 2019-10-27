CC=gcc
CFLAGS=-O3 -Wall

all: 	ssdraw

ssdraw2d.o:  ssdraw2d.c ReadCommandLine.c InputFileReader.c RenderCirclePolyGraph.c utils.c SecStructDrawHeader.h
	$(CC) $(CFLAGS) -c ssdraw2d.c

ReadCommandLine.o: ReadCommandLine.c SecStructDrawHeader.h
	$(CC) $(CFLAGS) -c ReadCommandLine.c

InputFileReader.o: InputFileReader.c utils.c SecStructDrawHeader.h
	$(CC) $(CFLAGS) -c InputFileReader.c

utils.o: utils.c SecStructDrawHeader.h
	$(CC) $(CFLAGS) -c utils.c

RenderCirclePolyGraph.o: RenderCirclePolyGraph.c utils.c SecStructDrawHeader.h
	$(CC) $(CFLAGS) -c RenderCirclePolyGraph.c

DrawCanonicalSS.o: DrawCanonicalSS.c utils.c SecStructDrawHeader.h
	$(CC) $(CFLAGS) -c DrawCanonicalSS.c

ssdraw: ssdraw2d.o RenderCirclePolyGraph.o DrawCanonicalSS.o ReadCommandLine.o InputFileReader.o utils.o
	$(CC) $(CFLAGS) ssdraw2d.o RenderCirclePolyGraph.o DrawCanonicalSS.o ReadCommandLine.o InputFileReader.o utils.o -o ssdraw2d -lm

clean: 
	rm -f *.o
	rm -f ssdraw2d
	rm -f *~
