CC=clang
CFLAGS= -Wall -g

all: randmst

randmst: mst.o creategraph.o node.o
	$(CC) mst.o creategraph.o node.o -o randmst

pa1.o: pa1.cc
	$(CC) $(CFLAGS) pa1.cc


clean: 
	rm *.o randmst

	g++ -c pa1.cc
	