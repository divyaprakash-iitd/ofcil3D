CC	= gcc
FORT	= gfortran
OPTS    = -O3 -fbounds-check

OBJECTS= mod_io.o fem2d.o mod_solid.o soft_particles.o

INC_DIR=

LIB_L= 

#ibmc: $(OBJECTS) 
#	$(FORT) $(OPTS) -o $@ $(OBJECTS) $(LIB_DIR) $(LIB_L)

fem2d.o: fem2d.f90 mod_solid.o Makefile
	$(FORT) -c $(OPTS) $<

mod_io.o: mod_io.f90 Makefile
	$(FORT) -c $(OPTS) $<

mod_solid.o: mod_solid.f90 Makefile
	$(FORT) -c $(OPTS) $<

soft_particles.o: soft_particles.f90 mod_io.o fem2d.o Makefile
	$(FORT) -c -fPIC $(OPTS) $<

clean:
	rm -f *.mod *.bak *~ $(OBJECTS) 
