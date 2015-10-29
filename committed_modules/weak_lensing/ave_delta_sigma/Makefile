include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS = -I ${GSL_LIB} 

USER_LDFLAGS= -L ${GSL_LIB} -lgsl -lgslcblas -lcosmosis

all: ave_delta_sigma.so

ave_delta_sigma.so: calc_ave_delta_sigma.c ave_delta_sigma.c 
	${CC} ${CFLAGS} -shared -o ave_delta_sigma.so calc_ave_delta_sigma.c ave_delta_sigma.c $(LDFLAGS)

.PHONY: clean

clean:
	rm -f ave_delta_sigma.so ave_delta_sigma.so.dSYM ave_delta_sigma.o