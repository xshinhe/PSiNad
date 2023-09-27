module diffsizes
   integer,parameter :: nbdirsmax = 12   ! 3*4 for H2O2
end module diffsizes

!--- PUSH REAL8 values (x, xd(:))to memory
subroutine pushreal8_dv(x,xd,nbdirs)
   use diffsizes
   integer :: nbdirs
   real(8) :: x, xd(nbdirsmax)
   integer :: i
   call pushreal8(x)
   call pushreal8array(xd,nbdirs)
!   do i = 1, nbdirsmax
!      call pushreal8(xd(i))
!   end do
end subroutine pushreal8_dv

!--- POP REAL8 values (x, xd(:)) from memory
subroutine popreal8_dv(x,xd,nbdirs)
   use diffsizes
   integer :: nbdirs
   real(8) :: x, xd(nbdirsmax)
   integer :: i
!   do i = nbdirsmax,1,-1
!      call popreal8(xd(i))
!   end do
call popreal8array(xd,nbdirs)
   call popreal8(x)
end subroutine popreal8_dv

!--- PUSH REAL8 values (x(:), xd(:,:)) to memory
subroutine pushreal8array_dv(x,xd,n,nbdirs)
   use diffsizes
   integer :: n,nbdirs
   real(8) :: x(n), xd(nbdirsmax,n)
   integer :: i
   do i = 1, n
      call pushreal8_dv(x(i),xd(:,i),nbdirs)
   end do
end subroutine pushreal8array_dv

!--- POP REAL8 values (x(:),xd(:,:)) from memory
subroutine popreal8array_dv(x,xd,n,nbdirs)
   use diffsizes
   integer :: n,nbdirs
   real(8) :: x(n), xd(nbdirsmax,n)
   integer :: i
   do i = n, 1, -1
      call popreal8_dv(x(i),xd(:,i),nbdirs)
   end do
end subroutine popreal8array_dv
