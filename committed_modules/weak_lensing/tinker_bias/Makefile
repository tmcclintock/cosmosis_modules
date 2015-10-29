include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS = -I ${GSL_LIB}

USER_LDFLAGS= -L ${GSL_LIB} -lgsl -lgslcblas -lcosmosis

all: tinker_bias.so

tinker_bias.so: calc_tinker_bias.c tinker_bias.c 
	${CC} ${CFLAGS} -shared -o tinker_bias.so calc_tinker_bias.c tinker_bias.c $(LDFLAGS)

.PHONY: clean

clean:
	rm -f tinker_bias.so tinker_bias.so.dSYM tinker_bias.o