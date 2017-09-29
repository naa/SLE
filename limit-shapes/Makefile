# Customize below to fit your system
# includes and libs
INCS = -I. -I/usr/include 
LIBS = -L/usr/lib -lc -lm -lgmp

# flags
CFLAGS = -g -std=gnu99 -pedantic -Wall -O0 ${INCS} ${CPPFLAGS}
LDFLAGS = -g ${LIBS}

# compiler and linker
CC = cc

SRC = wang.c
OBJ = ${SRC:.c=.o}

all: wang

.c.o:
	@echo CC $<
	@${CC} -c ${CFLAGS} $<

wang: ${OBJ}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ} ${LDFLAGS}

