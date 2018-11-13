CC=g++
COPTS=`root-config --cflags` -I/usr/local/root/include -g
LDOPTS=`root-config --libs` -g 
CFLAGS=-Wall -pedantic -std=c++11 -O3
AUTHOR= Adrian Setniewski
#C++ Files
SOURCES =  main.cpp  Dict.cpp Vec3.cpp
OBJECTS = $(SOURCES:.cpp=.o)
#Dictionary classes
HEADERS = Vec3.hpp
EXECUTABLE=program
all: $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $^ $(LDOPTS) $(CFLAGS)
	rm -f *.o
#C++ files
.cpp.o:
	$(CC) -o $@ $^ -c $(COPTS) $(CFLAGS)
#Dictionary for ROOT classes
Dict.cpp: $(HEADERS)
	@echo "Generating dictionary ..."
	@rootcint -f  Dict.cpp -c -P -I$ROOTSYS  $(HEADERS)
clean:
	@rm -f $(OBJECTS)  $(EXECUTABLE) *.o *.d Dict.cpp Dict.h
info:
	echo Compiler - $(CC) Flags - $(CFLAGS) Copyright $(AUTHOR)


