include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS = -I ${GSL_LIB}

USER_LDFLAGS= -L ${GSL_LIB} -lgsl -lgslcblas -lcosmosis

all: xi_mm.so

xi_mm.so: calc_xi_mm.c xi_mm.c 
	${CC} ${CFLAGS} -shared -o xi_mm.so calc_xi_mm.c xi_mm.c $(LDFLAGS)

.PHONY: clean

clean:
	rm -f xi_mm.so xi_mm.so.dSYM xi_mm.o