CC=gcc
CFLAGS=-g -DDEBUG -Wall
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

BUILD_OBJS = tsp.o \
	tsp_cplex.o \
	util.o \
	tsp_greedy.o \
	tsp_tabu.o \
	tsp_vns.o \
	eventlog.o \
	tsp_instance.o

all: main

%.o : %.c %.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

main: main.c $(BUILD_OBJS)
	$(CC) $(CFLAGS) $(LINKFLAGS) $(INCFLAGS) $(BUILD_OBJS) main.c $(LINK) $(CPLEX_LIBS) -o main

clean:
	rm -f $(BUILD_OBJS) main

format:
	find . -iname "*.h" -o -iname "*.c" | xargs clang-format -i
