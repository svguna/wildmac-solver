UNAME := $(shell uname)
CFLAGS = -Wall

PROB_SOURCES=chain.c hashtable.c probability_chain.c energy.c pdf_chain.c pthread_sem.c hashkeys.c probability.c prob-solver.c
PROB_OBJECTS=$(PROB_SOURCES:.c=.o)

DET_SOURCES=det-solver.c
DET_OBJECTS=$(DET_SOURCES:.c=.o)

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

CFLAGS += ${INCDIRS} -O3

all: det-solver prob-solver

prob-solver: ${PROB_OBJECTS}
	${CC} -o $@ ${PROB_OBJECTS} ${LDFLAGS} ${CFLAGS}

det-solver: ${DET_OBJECTS}
	${CC} -o $@ ${DET_OBJECTS}

%.o: %.c
	${CC} -c $< ${INCDIRS} ${CFLAGS}

clean:
	rm *.o prob-solver det-solver
