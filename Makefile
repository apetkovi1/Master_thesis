CC=g++
CFLAGS=-I. $(shell root-config --cflags)
LDFLAGS = $(shell root-config --libs)
DEPS = Analyzer.h cConstants.h 

%.o: %.cpp $(DEPS)
	$(CC) $(LDFLAGS) -c -o $@ $< $(CFLAGS)

analyzer: analyze.o Analyzer.o
	$(CC) $(LDFLAGS) -o analyze analyze.o Analyzer.o

clean:
	rm -rf *.o analyze
