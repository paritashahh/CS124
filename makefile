CC=g++
CFLAGS= -std=c++11

all: randmst

randmst:
	$(CC) $(CFLAGS) indexpriorityqueue.cpp

clean: 
	rm *.o randmst
	g++ -c indexpriorityqueue.cpp
	