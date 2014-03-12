OBJECT  = 3d_gpe
OBJS    = 3d_gpe.o
#FC     = sunf95
FC      = gfortran
#FFLAGS = -g -fast
FFLAGS =  -mcmodel=medium -fopenmp#-fno-range-check
#LIB = -llapack
#FFLAGS = -pg 
#LDFLAGS =
#INCLUDE        =
#-----------------------------------------------------------------------
%.o : %.f90
        $(FC) $(FFLAGS) $(INCLUDE) -c $*.f90

all:    $(OBJECT)

clean :
        rm -f $(OBJECT) *.o *.mod
#-----------------------------------------------------------------------
$(OBJECT): $(OBJS)
        $(FC) $(FFLAGS) $(LIB) $(INCLUDE) $(OBJS) $(LDFLAGS)  -o $@

#diffusion.o: ic_bc.o io.o parameters.o solve.o transform.o variables.o
~                                                                             
