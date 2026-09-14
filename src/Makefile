# Makefile for lassip

#debugging
#G++FLAG = -g -Wall -Werror
#For dynamic linking

CC = g++

# Optimisation, architecture and warning flags are separate variables so that any
# of them can be overridden from the command line. ARCHFLAGS is empty by default:
# the -m64 -mmmx -msse -msse2 this Makefile used to pass unconditionally are
# x86-only, and clang on arm64 accepts but ignores them.
#
#   make                                        # portable -O3 build
#   make ARCHFLAGS="-m64 -mmmx -msse -msse2"    # x86-64, as lassip <= 1.2.2 built
#   make ARCHFLAGS="-march=native"              # tuned for this machine
#   make OPTFLAGS="-O0 -g" WARNFLAGS="-Wall -Wextra"   # debugging
#   make check                                  # build, then run tests/run_tests.sh
STDFLAG   ?= -std=c++11
OPTFLAGS  ?= -O3
ARCHFLAGS ?=
WARNFLAGS ?=
G++FLAG   ?= $(STDFLAG) $(OPTFLAGS) $(ARCHFLAGS) $(WARNFLAGS)

# -MMD -MP makes the compiler emit a .d file per object listing the headers it
# included, so that editing a header rebuilds what depends on it. Without this,
# a header-only change left stale objects to be linked against new declarations.
DEPFLAGS = -MMD -MP

LINK_OPTS2 = -lz -pthread


#for windows, using MinGW build environment
#CC = g++.exe
#G++FLAG = -DPTW32_STATIC_LIB -O3 -static-libgcc -static-libstdc++
#LINK_OPTS2 = ../lib/win32/libpthreadGC2.a ../lib/win32/libz.a 
#I_PATH = -I../include -I../include/win32 
#L_PATH = ../lib/win32


SRCS = lassip-main.cpp lassip-data.cpp param_t.cpp gzstream.cpp \
       lassip-winstats.cpp lassip-wintools.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.d)

.PHONY : all check clean

all : lassip

lassip : $(OBJS)
	$(CC) $(G++FLAG) -o $@ $(OBJS) $(LINK_OPTS2) $(LINK_OPTS)

%.o : %.cpp
	$(CC) $(G++FLAG) $(DEPFLAGS) -c $< $(I_PATH)

check : lassip
	../tests/run_tests.sh --bin ./lassip

clean :
	rm -f $(OBJS) $(DEPS) lassip

-include $(DEPS)
