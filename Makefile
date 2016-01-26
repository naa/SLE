# Customize below to fit your system
# includes and libs
INCS = -I. -I/usr/include 
LIBS = -L/usr/lib -lc -lm -lgmp

# flags
CFLAGS = -g -std=c99 -pedantic -Wall -O0 ${INCS} ${CPPFLAGS}
LDFLAGS = -g ${LIBS}

# compiler and linker
CC = cc

SRC = ising.c
SRC2 = wang.c
OBJ = ${SRC:.c=.o}
OBJ2 = ${SRC2:.c=.o}

all: ising wang

.c.o:
	@echo CC $<
	@${CC} -c ${CFLAGS} $<

ising: ${OBJ}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ} ${LDFLAGS}

wang: ${OBJ2}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ2} ${LDFLAGS}
