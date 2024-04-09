module soft_particles
    use mod_io
    use fem2d
    implicit none

    ! Parameters
    real(8), parameter :: PI = 3.141592653589793
    
    ! Computational Domain
    real(8)    :: Lx = 25 !50 !0.3
    real(8)    :: Ly = 25 !50 !0.1
    
    ! Particle information
    integer         :: nvp ! Number of vertices in the particle
    real(real64)    :: aa, bb ! Major and minor axis
    real(8)         :: Kp, Bp
    integer(int32)  :: itnum ! iteration number

    ! FEM data
    type(festruct), allocatable :: particles(:)
    integer(int32)              :: ntri, npp, nparticle
    real(real64), allocatable   :: pb(:,:,:)
    integer,      allocatable   :: mp(:,:)
    real(real64), allocatable   :: paelem(:)
    real(real64), allocatable   :: pp(:,:)
    real(real64), allocatable   :: FN(:,:)
    real(real64), allocatable   :: UN(:,:)
    logical, allocatable        :: pboundary(:,:)
    integer(int32)              :: femdata(2)
    integer(int32)              :: i, j, err

    ! Namelists for input
    namelist /particleprops/ Kp, nvp, Bp, aa, bb

    !---------------------- Begin Calculations ------------------------------------!

    contains 

    subroutine sayhello() bind(C)
        use iso_c_binding, only: C_INT, C_CHAR
        implicit none
        print *, "Hello!"
    end subroutine sayhello

    subroutine generateellipse(noelpts) bind(C)
    use iso_c_binding, only: C_INT, C_CHAR
    implicit none

    integer(C_INT), intent(inout) :: noelpts

    ! Read input data from file
    open(1004,file="input_params.dat",form='formatted')
    READ(unit=1004,nml=particleprops,iostat=err)
    close(1004)


    ! Allocate fem variables
    !------- femdata----------!
    OPEN(33, FILE="femdata.bin",&
     FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
    READ(33) femdata
    close(33)
    !---------------------!
    npp         = femdata(1)
    ntri        = femdata(2)
    nparticle   = 1
    allocate(mp(ntri,3), pb(ntri,2,3), pp(npp,2), paelem(ntri), & 
                FN(npp,2), UN(npp,2), pboundary(npp,2), particles(nparticle))

    print *, npp, ntri
     
    ! Create fem particle
    call read_fem_data(mp,paelem,pb,pp)
    pboundary = .FALSE.
    PP(:,1) = PP(:,1) + Lx/2.0d0 
    PP(:,2) = PP(:,2) + Ly/2.0d0  
    particles(1) = festruct(MP,PP,PB,pboundary,paelem,FN,UN,bp,kp,1.0d0) ! kp = kval, co = bp , dl = 1.0d0

    ! Open the file for writing
    OPEN(UNIT=10, FILE="MP.txt", STATUS='replace', ACTION='write')
    ! Loop through the matrix and write its elements to the file
    DO i = 1, ntri
        WRITE(10, '(3I5)') (particles(1)%M(i, j), j = 1, 3)
    ENDDO
    ! Close the file
    CLOSE(10)

    itnum = 1
    call write_field(particles(1)%XE,'P',itnum)

    noelpts = npp

    end subroutine generateellipse
   
    subroutine arraycheck(pxyz,n) bind(C)
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none
        
        integer(c_int), intent(in) :: n
        real(c_double), intent(inout)   :: pxyz(n)
        ! integer(c_int), intent(inout)   :: pxyz(:)
        ! integer(c_int), intent(inout)   :: pxyz(n)
        
        integer(int32) :: i, npoints
        npoints = size(particles(1)%XE,1)
       
        ! print *, "SIZE: ", size(pxyz)
        ! print *, "pxyz1: ", pxyz(5)
        do i = 1,5!npoints
            ! pxyz(i)   = particles(1)%XE(i,1)
            pxyz(i)   = 23651491.23
        end do


    end subroutine arraycheck

    subroutine getpositions(XC,YC,ZC,nn) bind(C)
        ! It takes in the position arrays defined in openfoam and fills
        ! it with the particle's position values
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none

        integer(c_int), intent(in)      :: nn
        real(c_double), intent(inout)   :: XC(nn),YC(nn),ZC(nn)

        integer(int32) :: i, nparticles, npoints

        ! print *, "Size of XC: ", size(XC)
        nparticles = 1

        npoints = size(particles(1)%XE,1)

        do i = 1,npoints
            XC(i)   = particles(1)%XE(i,1)
            YC(i)   = particles(1)%XE(i,2)
            ZC(i)   = 0.005d0 ! Fill it with the value of the z-component of the mesh cell center
        end do
    end subroutine getpositions
   
    subroutine calculateforces(FXC,FYC,FZC,nn) bind(C)
        ! Calculates the forces in the particle
        ! Transfers those forces to the arrays passed in by openfoam
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none

        integer(c_int), intent(in)      :: nn
        real(c_double), intent(inout)   :: FXC(nn),FYC(nn),FZC(nn)

        integer(int32) :: i, nparticles, npoints

        nparticles = 1

        npoints = size(particles(1)%XE,1)

        do i = 1,nparticles
            call particles(i)%calculate_forces()
        end do

        do i = 1,npoints
            FXC(i)  = particles(1)%fden(i,1)
            FYC(i)  = particles(1)%fden(i,2)
            FZC(i)  = 0.0
        end do
    end subroutine calculateforces


    subroutine updatepositions(U,V,W,dt,nn) bind(C)
        ! Take in the empty position arrays from openfoam
        ! Fill it up with values
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none

        integer(c_int), intent(in)      :: nn
        real(c_double), intent(inout)   :: U(nn), V(nn) ,W(nn)
        real(c_double), intent(in)      :: dt

        integer(int32) :: i, nparticles, npoints

        nparticles = 1

        npoints = size(particles(1)%XE,1)
        
        do i = 1,npoints
            particles(1)%U(i,1) = U(i)
            particles(1)%U(i,2) = V(i)
        end do

        itnum = itnum + 1
        
        if (mod(itnum,200).eq.0) then
            call write_field(particles(1)%XE,'P',itnum)
        end if

        do i = 1,nparticles
            call particles(i)%update_position(dt)
        end do

    end subroutine updatepositions

    ! Create 3 subroutines
    ! 1. Creates the ellipse and it's coordinates and connectivity. Basically reads it form the python generated file.
    ! 2. Calls the force calculation and fills up the force vector.
    ! 3. Takes the velocity vector from the C program and uses it to update the particle's nodes positions.

end module soft_particles
