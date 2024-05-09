FC:=gfortran
FCFLAGS=-fdefault-real-8 -fimplicit-none

pso: pso.f90
	${FC} ${FCFLAGS} -fopenmp $^ -o $@

noparallel: pso.f90
	${FC} ${FCFLAGS} $^ -o $@

clean:
	-rm *.com *.com_energy OUT
	-rm pso noparallel

pso-methane: pso-methane.f90
	${FC} ${FCFLAGS} -fopenmp $^ -o $@

parse.o: parse.f90
	${FC} ${FCFLAGS} -c $^

pso-ethane: parse.o pso-ethane.f90
	${FC} ${FCFLAGS} $^ -o $@

.PHONY: clean
