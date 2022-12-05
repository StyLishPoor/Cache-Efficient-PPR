CC = g++
CPPFLAGS = -Wno-deprecated -O3 -c -m64 -march=native -std=c++17 -DSFMT_MEXP=19937 -I SFMT-src-1.5.1/
LDFLAGS = -static -O3 -m64 -DSFMT_MEXP=19937 -I SFMT-src-1.5.1/
SOURCES = main.cpp graph.cpp random.cpp helper.cpp util.cpp SFMT-src-1.5.1/SFMT.c
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = SSPPR

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o :
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o
