
  PROG = ~/bin/pdb2angles


  OBJS = pdb2angles.o
#--------------------------------------------------------------------------
# Testing for 1. gfortran 2. ifort
FC := $(shell which gfortran 2>> error.dat)
FC := $(shell which ifort 2> error.dat) 
ifndef FC
#ifeq ($(FC),)
$(warning ifort not found.)
$(warning gfortran not found.)
$(warning Please adjust FC= yourself!)
$(error aborting...)
endif

#--------------------------------------------------------------------------
  FC = ifort 
#  FC = gfortran 
  LINKER = $(FC) -static
FFLAGS =    -O 
         
.PHONY: all
.PHONY: clean
.PHONY: sync


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@


# link
$(PROG):$(OBJS) 
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)


clean:
	rm -f *.o *.mod $(PROG) 
	rm -f $(patsubst %.F, %.f, $(wildcard *.F))


sync:
	rsync -vP $(PROG) kruse@alpha:~/bin/
	rsync -vP $(PROG) kruse@tripura:~/bin/
	rsync -vP $(PROG) kruse@195.178.69.202:~/bin/
#	cp $(PROG) ~/bin/
