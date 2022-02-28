CC=g++
CFLAGS= -std=c++11 -c

all: randmst

randmst:
	$(CC) $(CFLAGS) justprims.cpp

clean: 
	rm *.o justprims.cpp
	