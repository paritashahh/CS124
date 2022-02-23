CC=clang
CFLAGS= -Wall -g

all: randmst

randmst: pa1.o indexpriorityqueue.o
	$(CC) pa1.o indexpriorityqueue.o -o randmst

test: indexpriorityqueue.o
	$(CC) indexpriorityqueue.o -o randmst

pa1.o: pa1.cc
	$(CC) $(CFLAGS) pa1.cc

indexpriorityqueue.o: indexpriorityqueue.cc
	$(CC) $(CFLAGS) indexpriorityqueue.cc

clean: 
	rm *.o randmst
	g++ -c pa1.cc indexpriorityqueue.cc
	