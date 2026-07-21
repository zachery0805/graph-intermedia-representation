CC = gcc -std=c99 -fopenmp

bic: bic.o graph.o
	$(CC) -o bic bic.c graph.o

bic-serial: bic-serial.o graph.o
	$(CC) -o bic-serial bic-serial.c graph.o

graph.o: graph.c graph.h
	$(CC) -c graph.c