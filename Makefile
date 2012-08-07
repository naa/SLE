# Customize below to fit your system
# includes and libs
INCS = -I. -I/usr/include 
LIBS = -L/usr/lib -lc -lm

# flags
CFLAGS = -g -std=c99 -pedantic -Wall -O0 ${INCS} ${CPPFLAGS}
LDFLAGS = -g ${LIBS}

# compiler and linker
CC = cc

SRC = ising.c
OBJ = ${SRC:.c=.o}

all: ising

.c.o:
	@echo CC $<
	@${CC} -c ${CFLAGS} $<

ising: ${OBJ}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ} ${LDFLAGS}
