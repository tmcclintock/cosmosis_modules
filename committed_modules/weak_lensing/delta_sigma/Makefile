include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS = -I ${GSL_LIB} 

USER_LDFLAGS= -L ${GSL_LIB} -lgsl -lgslcblas -lcosmosis

all: delta_sigma.so

delta_sigma.so: calc_delta_sigma.c delta_sigma.c 
	${CC} ${CFLAGS} -shared -o delta_sigma.so calc_delta_sigma.c delta_sigma.c $(LDFLAGS)

.PHONY: clean

clean:
	rm -f delta_sigma.so delta_sigma.so.dSYM delta_sigma.o