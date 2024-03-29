CC=gcc
CFLAGS=-g -DDEBUG
LINK=-lm
# CFLAGS= -g -O3

# from command line
CPLEX_PATH ?= /opt/ibm/ILOG/CPLEX_Studio2211

CPLEX_INC_DIR = $(CPLEX_PATH)/cplex/include
CPLEX_LIB_DIR = $(CPLEX_PATH)/cplex/lib/x86-64_linux/static_pic
CPLEX_CONCERT_INC_DIR = $(CPLEX_PATH)/concert/include
CPLEX_CONCERT_LIB_DIR = $(CPLEX_PATH)/concert/lib/x86-64_linux/static_pic
CPLEX_LIBS = -lcplex -lilocplex -lconcert -lpthread

LINKFLAGS = -L$(CPLEX_CONCERT_LIB_DIR) -L$(CPLEX_LIB_DIR)
INCFLAGS = -I$(CPLEX_CONCERT_INC_DIR) -I$(CPLEX_INC_DIR)

all: main

format:
	find . -iname "*.h" -o -iname "*.c" | xargs clang-format -i

tsp.o: tsp.h tsp.c
	$(CC) $(CFLAGS) -c tsp.c -o tsp.o

tsp_cplex.o: tsp_cplex.h tsp_cplex.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c tsp_cplex.c -o tsp_cplex.o

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

main: main.c tsp.h tsp.o util.h util.o tsp_greedy.h tsp_greedy.o tsp_tabu.h tsp_tabu.o tsp_vns.h tsp_vns.o eventlog.o tsp_cplex.o
	$(CC) $(CFLAGS) $(LINKFLAGS) $(INCFLAGS) tsp.o util.o tsp_greedy.o tsp_tabu.o tsp_vns.o eventlog.o tsp_cplex.o main.c $(LINK) $(CPLEX_LIBS) -o main

clean:
	rm -f tsp.o util.o tsp_greedy.o tsp_tabu.o tsp_vns.o eventlog.o tsp_cplex.o main
