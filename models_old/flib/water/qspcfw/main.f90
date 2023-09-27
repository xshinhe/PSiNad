program main
	implicit none

	integer::natoms,i
    double precision::pbox(3)
    double precision::e,cutoff,ewald_parm
    double precision,allocatable::x(:), f(:), q(:)

    !ordered by O-H1-H2, x-y-z
    natoms = 801
    pbox = (/20.008000d0,  20.079000d0,  20.057000d0/)
!     pbox = (/22.000000d0,  22.000000d0,  22.000000d0/)
    cutoff = 10.0d0
    ewald_parm = 0.34864
    allocate(x(3*natoms))
    allocate(f(3*natoms))
    allocate(q(natoms))
!     x=(/0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, -0.6d0, 0.8d0, 0.0d0/)

    open(1,file='H2O267.xyz')
    read(1,*) x
+
!        do i = 1, natoms
!          read(1,*) x(3*i-2), x(3*i-1), x(3*i), q(i)
!        enddo

    
    do i=1,natoms
    	if(mod(i,3).eq.1) then
    		q(i)=-0.82d0
    	else
    		q(i)=0.41d0
    	end if
    end do

	call get_force(natoms, pbox, x, q, f, e, cutoff, ewald_parm)

! 	write(*,*) 'check, finish get_force'
    write(*,*) e
!     write(*,*) f

! 	write(*,*) f
! 	write(*,*) e

end program