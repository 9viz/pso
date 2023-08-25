FC:=gfortran

pso: pso.f90
	${FC} -fopenmp $< -o $@

noparallel: pso.f90
	${FC} $< -o $@

clean:
	-rm *.com *.com_energy OUT
	-rm pso noparallel

.PHONY: clean
