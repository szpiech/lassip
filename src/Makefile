# Makefile for lassip

#debugging
#G++FLAG = -g -Wall -Werror
#For dynamic linking
#LINK_OPTS = -L$(L_PATH) -lgsl -lgslcblas

#for osx systems
CC = g++
G++FLAG = -O3 -m64 -mmmx -msse -msse2
#LINK_OPTS2 = -pthread -lz -lgsl -lgslcblas
LINK_OPTS2 = -lz -pthread
#-lgsl -lgslcblas
I_PATH = -I../include
L_PATH = ../lib/osx

#for linux systems
#CC = g++
#G++FLAG = -O3 -m64 -mmmx -msse -msse2
#LINK_OPTS2 = -pthread -lz
#I_PATH = -I../include
#L_PATH = ../lib/linux

#for windows, using MinGW build environment
#CC = g++.exe
#G++FLAG = -DPTW32_STATIC_LIB -O3 -static-libgcc -static-libstdc++
#LINK_OPTS2 = ../lib/win32/libpthreadGC2.a ../lib/win32/libz.a 
#I_PATH = -I../include -I../include/win32 
#L_PATH = ../lib/win32

#For static linking of norm program to gsl libs
LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a

all : lassip

lassip : lassip-main.o lassip-data.o binom.o param_t.o gzstream.o lassip-winstats.o lassip-wintools.o
	$(CC) $(G++FLAG) -o lassip lassip-main.o lassip-data.o binom.o param_t.o gzstream.o lassip-winstats.o lassip-wintools.o $(LINK_OPTS2) $(LINK_OPTS)

lassip-main.o : lassip-main.cpp
	$(CC) $(G++FLAG) -c lassip-main.cpp $(I_PATH)

lassip-data.o : lassip-data.cpp
	$(CC) $(G++FLAG) -c lassip-data.cpp $(I_PATH)

lassip-winstats.o : lassip-winstats.cpp
	$(CC) $(G++FLAG) -c lassip-winstats.cpp $(I_PATH)

lassip-wintools.o : lassip-wintools.cpp
	$(CC) $(G++FLAG) -c lassip-wintools.cpp $(I_PATH)

binom.o : binom.cpp
	$(CC) $(G++FLAG) -c binom.cpp

param_t.o : param_t.cpp
	$(CC) $(G++FLAG) -c param_t.cpp

gzstream.o : gzstream.cpp
	$(CC) $(G++FLAG) -c gzstream.cpp $(I_PATH)

clean :
	rm *.o
