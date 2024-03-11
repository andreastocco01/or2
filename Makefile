CC=gcc
CFLAGS=-g -DDEBUG
LINK=-lm
# CFLAGS= -lm

all: main

format:
	find . -iname "*.h" -o -iname "*.c" | xargs clang-format -i

tsp.o: tsp.h tsp.c
	$(CC) $(CFLAGS) -c tsp.c -o tsp.o

util.o: util.h util.c
	$(CC) $(CFLAGS) -c util.c -o util.o

tsp_greedy.o: tsp_greedy.h tsp_greedy.c
	$(CC) $(CFLAGS) -c tsp_greedy.c -o tsp_greedy.o

main: main.c tsp.h tsp.o util.h util.o tsp_greedy.h tsp_greedy.o
	$(CC) $(CFLAGS) tsp.o util.o tsp_greedy.o main.c $(LINK) -o main

clean:
	rm -f tsp.o util.o tsp_greedy.o main
