UNAME := $(shell uname)
CFLAGS = -Wall

SOURCES=$(shell ls *.c)
OBJECTS=$(SOURCES:.c=.o)

ifeq ($(UNAME), Linux)
LDFLAGS = -lgsl -lgslcblas -lnlopt
INCDIRS =
endif

ifeq ($(UNAME), Darwin)
LDFLAGS = -L/opt/local/lib -lgsl 
INCDIRS = -I/opt/local/include
endif

CFLAGS += ${LDFLAGS} ${INCDIRS} -g 

all: solver

solver: ${OBJECTS}
	${CC} -o $@ ${OBJECTS} ${CFLAGS}

%.o: %.c
	${CC} -c $< ${INCDIRS}

clean:
	rm *.o solver
