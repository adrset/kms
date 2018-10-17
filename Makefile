CC=g++
CFLAGS=-Wall -pedantic -std=c++11
AUTHOR= Adrian Setniewski
all: pr clean
pr: main.o Vec3.o
	$(CC) -o program main.o Vec3.o $(CFLAGS)
main.o: main.cpp
	$(CC) -o main.o -c main.cpp $(CFLAGS)
Vec3.o: Vec3.cpp
	$(CC) -o Vec3.o -c Vec3.cpp $(CFLAGS)
clean:
	rm -f *.o
info:
	echo Compiler - $(CC) Flags - $(CFLAGS) Copyright $(AUTHOR)

