gfortran -c -fPIC fem3d.f90 -o fem3d.o

# Compile soft_particles.f90 and link with object files to create shared library
gfortran -c -fPIC soft_particles.f90 -o soft_particles.o
gfortran -shared -o libsoftparticles3D.so soft_particles.o fem3d.o

