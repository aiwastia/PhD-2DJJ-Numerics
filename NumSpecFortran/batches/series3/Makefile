FC = gfortran
FFLAGS = -llapack

OBJS = selectedsubr.o io.o zufall.o scans.o JJ_gap.o

VPATH= src

dir: datatest
	python dirmaker.py

datatest: JJ_gap infile_JJ_bands_1 infile_JJ_2
	$(shell ./$^)

%.mod %.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@


JJ_gap: $(OBJS)
	$(FC) $(OBJS) $(FFLAGS) -o $@

clean:
	rm *.mod *.o

pplot: dir
	$(shell gnuplot pltpng)
