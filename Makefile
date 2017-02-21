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
SRC3 = tricritical.c
SRC4 = tricritical-wang.c
SRC5 = interface.c
OBJ = ${SRC:.c=.o}
OBJ2 = ${SRC2:.c=.o}
OBJ3 = ${SRC3:.c=.o}
OBJ4 = ${SRC4:.c=.o}
OBJ5 = ${SRC5:.c=.o}

all: ising wang tricritical tricritical-wang interface stats

.c.o:
	@echo CC $<
	@${CC} -c ${CFLAGS} $<

ising: ${OBJ}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ} ${LDFLAGS}

wang: ${OBJ2}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ2} ${LDFLAGS}

tricritical: ${OBJ3}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ3} ${LDFLAGS}


tricritical-wang: ${OBJ4}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ4} ${LDFLAGS}


interface: ${OBJ5}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ5} ${LDFLAGS}

stats:
	gcc -lm -o stats stats.c
