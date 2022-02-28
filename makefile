CC=g++
CFLAGS= -std=c++11 -c

all: randmst

randmst:
	$(CC) $(CFLAGS) indexpriorityqueue.cpp

new:
	$(CC) $(CFLAGS) justprims.cpp

clean: 
	rm *.o indexpriorityqueue.cpp justprims.cpp
	