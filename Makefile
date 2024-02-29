CC=gcc
CFLAGS=-g

all: main

tsp.o: tsp.h tsp.c
	$(CC) $(CFLAGS) -c tsp.c -o tsp.o

main: main.c tsp.h tsp.o
	$(CC) $(CFLAGS) tsp.o main.c -o main

