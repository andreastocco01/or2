CC=gcc
CFLAGS=-g -DDEBUG
LINK=-lm
# CFLAGS= -g -O3

all: main

format:
	find . -iname "*.h" -o -iname "*.c" | xargs clang-format -i

tsp.o: tsp.h tsp.c
	$(CC) $(CFLAGS) -c tsp.c -o tsp.o

util.o: util.h util.c
	$(CC) $(CFLAGS) -c util.c -o util.o

tsp_greedy.o: tsp_greedy.h tsp_greedy.c
	$(CC) $(CFLAGS) -c tsp_greedy.c -o tsp_greedy.o

tsp_tabu.o: tsp_tabu.h tsp_tabu.c
	$(CC) $(CFLAGS) -c tsp_tabu.c -o tsp_tabu.o

tsp_vns.o: tsp_vns.h tsp_vns.c
	$(CC) $(CFLAGS) -c tsp_vns.c -o tsp_vns.o

eventlog.o: eventlog.h eventlog.c
	$(CC) $(CFLAGS) -c eventlog.c -o eventlog.o

main: main.c tsp.h tsp.o util.h util.o tsp_greedy.h tsp_greedy.o tsp_tabu.h tsp_tabu.o tsp_vns.h tsp_vns.o eventlog.o
	$(CC) $(CFLAGS) tsp.o util.o tsp_greedy.o tsp_tabu.o tsp_vns.o eventlog.o main.c $(LINK) -o main

clean:
	rm -f tsp.o util.o tsp_greedy.o tsp_tabu.o tsp_vns.o main
