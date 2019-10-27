CC=gcc
CFLAGS=-g -Wall

all: 	CirclePolyGraph

CirclePolyGraph.o:  CirclePolyGraph.c ReadCommandLine.c InputFileReader.c CirclePolyGraphHeaderFile.h
	$(CC) $(CFLAGS) -c CirclePolyGraph.c

ReadCommandLine.o: ReadCommandLine.c CirclePolyGraphHeaderFile.h
	$(CC) $(CFLAGS) -c ReadCommandLine.c

InputFileReader.o: InputFileReader.c utils.c CirclePolyGraphHeaderFile.h
	$(CC) $(CFLAGS) -c InputFileReader.c

utils.o: utils.c CirclePolyGraphHeaderFile.h
	$(CC) $(CFLAGS) -c utils.c

CirclePolyGraph: CirclePolyGraph.o ReadCommandLine.o InputFileReader.o utils.o
	$(CC) $(CFLAGS) CirclePolyGraph.o ReadCommandLine.o InputFileReader.o utils.o -o circlepolygraph -lm


clean: 
	rm -f *.o
	rm -f circlepolygraph
	rm -f *~
