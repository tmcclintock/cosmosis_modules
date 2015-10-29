include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS = -I ${GSL_LIB} 

USER_LDFLAGS= -L ${GSL_LIB} -lgsl -lgslcblas -lcosmosis

all: sigma_r.so

sigma_r.so: calc_sigma_r.c sigma_r.c 
	${CC} ${CFLAGS} -shared -o sigma_r.so calc_sigma_r.c sigma_r.c $(LDFLAGS)

.PHONY: clean

clean:
	rm -f sigma_r.so sigma_r.so.dSYM sigma_r.o