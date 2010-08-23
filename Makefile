UNAME := $(shell uname)
CFLAGS = -Wall

SOURCES=$(shell ls *.c)
OBJECTS=$(SOURCES:.c=.o)

ifeq ($(UNAME), Linux)
LDFLAGS = -lgsl -lgslcblas -lpthread
INCDIRS =
endif

ifeq ($(UNAME), Darwin)
LDFLAGS = -L/opt/local/lib -lgsl -lpthread
INCDIRS = -I/opt/local/include
endif

ifeq ($(NETSERVER), true)
	LDFLAGS = -L/shared/home-05/guna/installs/lib -lgsl -lgslcblas -lpthread
	INCDIRS = -I/shared/home-05/guna/installs/include
endif

CFLAGS += ${INCDIRS} -g

all: solver

solver: ${OBJECTS}
	${CC} -o $@ ${OBJECTS} ${LDFLAGS} ${CFLAGS}

%.o: %.c
	${CC} -c $< ${INCDIRS} ${CFLAGS}

clean:
	rm *.o solver
