module constant
    implicit none
    integer :: P = 200
    integer, parameter :: constP = 200
    integer, parameter :: ndim=3! number of beads
    integer, parameter :: w=2 
    !real*8,parameter ::convert_fs=1.d-15  !s to fs
    !real*8,parameter ::dt=10.d0*convert_fs
    !    real*8, parameter :: t(0:2)=1.d0/3.d0, v(0:2)=1.d0/3.d0
    !t1=2t_0,t2=t3=t_1
    !real*8, parameter :: v(0:w-1)=0.5d0, t(0:w-1)=0.5d0
    integer,parameter :: period=6,Noa=216, ntot = ndim * noa
    real*8,parameter ::m=3.349339d-26
    real*8, parameter :: a0=5.29177d-11
    real*8, parameter :: Na=6.0221415d23
    real*8, parameter :: Eh=2.6255d6/Na
    real*8, parameter :: kb=1.3806488d-23
    real*8, parameter :: Jtoau=0.229371d18
    !...............
    !Lenard Jones
    !        real*8, parameter :: elj=35.6d0*kb
    real*8, parameter :: dlj=2.83834d-10!749d-10,
    real*8, parameter :: dlj6=dlj**6
    real*8, parameter :: a1 = 1.48735474d2, a2 = 2.3520917d0
    real*8, parameter :: a3 = 3.749754d4, a4 = 5.160348d3, a5 = 8.3724779d0, a6 = 6.8028878d0
    real*8,parameter :: dens=1.224d3
    real*8, parameter :: length=(m/dens*Noa)**(1.d0/dble(ndim))
    real*8, parameter :: le=length/period
    real*8, parameter :: pi=3.14159265359d0
    !real*8 :: a(P),b(P),c(P),q(P)!,coln(P)  ! coefficients for staging transform
    !real*8 :: mj(P) 
    real(8), parameter :: rcut=2.5d0*dlj,rlist=2.9d0*dlj
    real(8), parameter :: rcut2=rcut**2,rlist2=rlist**2
    integer, parameter :: nblmax=2000*(Noa+1),rfnbl=20,rgrid=400
end module

module params
    use constant
    implicit none
    real*8 :: dist2(nblmax, constP), dist(nblmax, constP)
    integer :: list(nblmax, constP), point(noa, constP), tcount
    integer :: rcount(rgrid)
    real*8 :: yj0(ndim, noa, constP) 
    !real*8, allocatable :: dist2(:, :), dist(:, :)
    !integer, allocatable :: list(:, :), point(:, :)
    !integer, allocatable :: rcount(:)
    !real*8, allocatable :: yj0(:, :, :) 
    integer :: cnt = 0
end module

subroutine init(cor, pntot) bind(c, name="liquidne_init_ccall")
    use constant
    use params
    use random
    implicit none

    integer pntot
    real*8 :: cor(pntot)
    real*8 rand
    real*8 fy(ndim, noa, P), y(ndim, noa, P)
    real*8 qn(ndim, noa)
    !real*8, allocatable :: fy(:, :, :), y(:, :, :)
    !real*8, allocatable :: qn(:, :)

    integer i, j, k, l

    P = pntot / ntot

    !allocate(dist2(nblmax, P))
    !allocate(dist(nblmax, P))
    !allocate(list(nblmax, P))
    !allocate(point(noa, P))
    !allocate(rcount(rgrid))
    !allocate(yj0(ndim, noa, P))

    !allocate(fy(ndim, noa, P))
    !allocate(y(ndim, noa, P))
    !allocate(qn(ndim, noa))

    !print *, "2 ",  P
    do j = 1, noa
    k = j
    do l = ndim, 1, -1
    i = k / period**(l-1)
    qn(l, j) = le*i
    k = k - i * period**(l-1)
    end do
    end do

    !print *, "3 ",  P
    do i = 1, P
    do j = 1, noa
    do k = 1, ndim
    call random_number(rand)
    y(k, j, i) = qn(k, j) + le*0.1d0*(rand-0.5)
    cor((i-1) * ntot + (j-1)*ndim + k) = y(k, j, i)
    end do
    end do
    !print *, "4 ",  P
    call nblist(y(:,:,i),yj0(:,:,i),list(:,i),point(:,i),rcount)
    call calfy(fy(:,:,i),y(:,:,i),dist2(:,i),dist(:,i),list(:,i),point(:,i))

    !if (allocated(fy)) deallocate(fy)
    !if (allocated(y)) deallocate(y)
    !if (allocated(qn)) deallocate(qn)
    end do
end subroutine

subroutine finalize() bind(c, name="liquidne_finalize_ccall")
    use constant
    use params
    implicit none

    !if (allocated(dist2)) deallocate(dist2)
    !if (allocated(dist)) deallocate(dist)
    !if (allocated(list)) deallocate(list)
    !if (allocated(point)) deallocate(point)
    !if (allocated(rcount)) deallocate(rcount)
    !if (allocated(yj0)) deallocate(yj0)

end subroutine

!function getenergy(cor) bind(c, name="liquidne_energy_ccall")
!    use constant
!    use params
!    implicit none
!
!    real*8 :: cor(:)
!    real*8 getenergy
!
!    !real*8 y(ndim ,noa, P), eptmp(P)
!    real*8, allocatable :: y(:, :, :), eptmp(:)
!
!    integer i, j, k
!
!    allocate(y(ndim, noa, P))
!    allocate(eptmp(P))
!
!    do k = 1, P
!    do i = 1, noa
!    do j = 1, ndim
!    y(j, i, k) = cor((k-1) * ntot + (i-1) * ndim + j)
!    end do
!    end do
!    end do
!    call calPotential(eptmp,y,dist2,dist,list,point)
!
!    getenergy = 0.0
!    do i = 1, P
!    getenergy = getenergy + eptmp(i)
!    end do
!
!    if (mod(cnt, 20) == 0) then
!        call checklist(y, yj0, list, point, tcount, rcount)
!        !print *, "check"
!    end if
!
!    if (allocated(y)) deallocate(y)
!    if (allocated(eptmp)) deallocate(eptmp)
!
!    return
!end function

subroutine calforce(cor, potential, force) bind(c, name="liquidne_force_ccall")
    use constant
    use params
    implicit none

    real*8 :: cor(ntot*P)
    real*8 :: force(ntot*P)
    real*8 :: potential(P)

    real*8 fy(ndim, noa), y(ndim, noa, P), eptmp(P)
    !real*8 fy(ndim, noa), eptmp(P)
    !real*8, allocatable :: y(:, :, :)
    integer i, j, k

    !allocate(y(ndim, noa, P))

    cnt = cnt + 1

    do k = 1, P
    do i = 1, noa
    do j = 1, ndim
    y(j, i, k) = cor((k-1) * ntot + (i-1)*ndim + j)
    !y(j, i, k) = yj0(j, i, k)
    end do
    end do


    call calfy(fy, y(:, :, k), dist2(:, k), dist(:, k), list(:, k), point(:, k))

    do i = 1, noa
    do j = 1, ndim
    force((k-1) * ntot + (i-1)*ndim + j) = fy(j, i) * Jtoau
    end do
    end do

    end do



    call calPotential(eptmp,y,dist2,dist,list,point)
    !call calPotential1(getforce,y,dist2(:, n),dist(:, n),list(:, n),point(:, n))

    do i = 1, P
    potential(i) = eptmp(i) * Jtoau
    end do

    if (mod(cnt, 20) == 0) then
        call checklist(y, yj0, list, point, tcount, rcount)
    end if


    !if (allocated(y)) deallocate(y)

    return 
end subroutine

subroutine calfy(fy,y,dist2,dist,list,point)!,test)
    use constant
    implicit none
    real(8) :: fy(ndim,Noa),y(ndim,Noa),dist2(nblmax),dist(nblmax)
    integer :: list(nblmax),point(Noa)
    !        integer,optional :: test
    real(8) :: force(ndim),vector(ndim)
    real(8) :: pow,ra,dist6,ra2
    integer :: i,j,n
    fy=0.0d0
    do i=1,Noa-1
    do j=point(i),point(i+1)-1
    n=list(j)
    call calDist(y(:,i),y(:,n),vector(:),dist2(j))
    dist(j)=sqrt(dist2(j))
    ra=dist(j)/a0
    ra2=ra**2
    dist6=dist2(j)**3/dlj6
    if (dist2(j)<rcut2) then
        force(:)=-Eh*(-a0*a1*a2*exp(-a2*ra)/dist(j)+((((-12.d0*a3/ra2+10.d0*a4)/ra2+&
            &8.d0*a5)/ra2+6.d0*a6)/ra2**4))/a0**2*vector(:)
        fy(:,i)=fy(:,i)+force(:)
        fy(:,n)=fy(:,n)-force(:)
    end if
    end do
    end do
end subroutine

subroutine restrictCoord(y)
    use constant
    implicit none
    real*8 :: y(ndim,Noa,P)
    !real*8 :: y(:,:,:)
    integer :: i,j,k
    do k=1,P
    do j=1,Noa
    do i=1,ndim
    y(i,j,k)=y(i,j,k)-dnint(y(i,j,k)/length)*length
    end do
    end do
    end do
end subroutine

subroutine calDist(y1,y2,vec,dist2)
    use constant
    implicit none
    real*8 :: y1(ndim),y2(ndim),vec(ndim),dist2
    integer :: i,j
    integer :: l
    real(8) :: mp,r
    do l=1,ndim
    mp=y1(l)-y2(l)
    vec(l)=mp-dnint(mp/length)*length
    if (abs(vec(l))>abs(vec(l)-length)) then
        !write(*,*) "small",vec(l),vec(l)-length,dnint(mp/length)
        vec(l)=vec(l)-length
    else !end if
        if (abs(vec(l))>abs(vec(l)+length)) then
            !write(*,*) "big",vec(l),vec(l)+length,dnint(mp/length)
            vec(l)=vec(l)+length
        end if
    end if
    end do
    dist2=sum(vec(:)**2)
end subroutine

subroutine calPotential1(eptmp,y,dist2,dist,list,point)
    use constant
    implicit none
    real(8) :: eptmp(1),y(ndim,Noa,1),dist2(nblmax,1),dist(nblmax,1)
    integer :: list(nblmax,1),point(Noa,1)
    real(8) :: force(ndim),vector(ndim)
    real(8) :: pow,ra,ra2
    integer :: i,j,k,l
    eptmp=0.0d0
    do k=1,1
    do i=1,Noa-1
    do j=point(i,k),point(i+1,k)-1
    l=list(j,k)
    !       call calDist(y(:,i,k),y(:,l,k),vector,dist2(j,k))
    ra=dist(j,k)/a0
    ra2=ra**2
    if (dist2(j,k)<rcut2) then
        eptmp(k)=eptmp(k)+Eh*(a1*exp(-a2*ra)+(((a3/ra2-a4)/ra2-a5)/ra2-a6)/ra2**3)
    end if
    end do
    end do
    end do
end subroutine


subroutine calPotential(eptmp,y,dist2,dist,list,point)
    use constant
    implicit none
    real(8) :: eptmp(P),y(ndim,Noa,P),dist2(nblmax,P),dist(nblmax,P)
    integer :: list(nblmax,P),point(Noa,P)
    !real(8) :: eptmp(:),y(:,:,:),dist2(:,:),dist(:,:)
    !integer :: list(:,:),point(:,:)
    real(8) :: force(ndim),vector(ndim)
    real(8) :: pow,ra,ra2
    integer :: i,j,k,l
    eptmp=0.0d0
    do k=1,P
    do i=1,Noa-1
    do j=point(i,k),point(i+1,k)-1
    l=list(j,k)
    !       call calDist(y(:,i,k),y(:,l,k),vector,dist2(j,k))
    ra=dist(j,k)/a0
    ra2=ra**2
    if (dist2(j,k)<rcut2) then
        eptmp(k)=eptmp(k)+Eh*(a1*exp(-a2*ra)+(((a3/ra2-a4)/ra2-a5)/ra2-a6)/ra2**3)
    end if
    end do
    end do
    end do
end subroutine

subroutine nblist(yj,yj0,list,point,rcount)
    use constant
    implicit none
    real*8 :: yj(ndim,Noa),yj0(ndim,Noa)!,pj(ndim,Noa)
    integer :: list(nblmax), point(Noa),rcount(rgrid)
    integer :: i,j,nlist,n,k
    real(8) :: dist,dist2,vector(ndim)
    list=0
    nlist=0
    point(1)=1
    do i=1,Noa-1
    do j=i+1,Noa
    call calDist(yj(:,i),yj(:,j),vector,dist2)
    if(dist2<rlist2) then
        nlist=nlist+1
        list(nlist)=j
        dist=sqrt(dist2)
        n=int(100.0d0*dist/dlj)+1
        if(n<=rgrid) then
            rcount(n)=rcount(n)+1
        end if
    end if
    end do
    if(nlist>nblmax-100) then
        write(*,*) "list is too small"
        !                        call savepq(yj,pj)
        stop
    end if
    point(i+1)=nlist+1
    end do
    yj0=yj
end subroutine
subroutine checklist(yj,yj0,list,point,tcount,rcount)
    use constant
    implicit none
    real*8 :: yj(ndim,Noa,P),yj0(ndim,Noa,P)
    integer :: list(nblmax,P),point(Noa,P)
    !real*8 :: yj(:,:,:),yj0(:,:,:)
    !integer :: list(:,:),point(:,:)
    integer :: tcount, rcount(rgrid)
    integer :: l,k
    real*8 :: maxa,maxb,d2
    real*8 :: vec(ndim)
    do k=1,P
    maxb=0.d0
    maxa=0.d0
    do l=1,Noa
    call calDist(yj(:,l,k),yj0(:,l,k),vec,d2)
    if (d2>maxb) then
        maxb=d2
        if(d2>maxa) then
            maxb=maxa
            maxa=d2
        end if
    end if
    end do
    if (sqrt(maxa)+sqrt(maxb)>rlist-rcut) then
        call nblist(yj(:,:,k),yj0(:,:,k),list(:,k),point(:,k),rcount)
        tcount=tcount+1
    end if
    end do
end subroutine
