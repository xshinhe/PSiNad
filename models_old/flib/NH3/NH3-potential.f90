!
! ================================================================================
!
        ! Etot   (energy)
        ! force  (force vector)
        ! d2v    (Hessian Matrix)
        ! xtemp  (Coordinate)
        ! 1.  Only force
        ! 2.  Force and Hessian
!
      !   call nh3_pot_cart(Etot, force, d2v, xtemp, 1)    !  Force
      !   call nh3_pot_cart(Etot, force, d2v, xtemp, 2)    !  Hessian 
!
! ================================================================================
!
! =====================================================
!
      Module NH3_INFO
        implicit none
        save
!
! -----------------------
!  Unit conversion
! -----------------------
!
        real*8,parameter:: ctime = 2.41887d-17    ! Time (s/au)
        real*8,parameter:: cal_to_Jole = 4.184d0   ! cal to Jole
        real*8,parameter:: n_aov = 6.022d23        ! Avogadro constant
        real*8,parameter:: m_atom = 1.661d-27      ! atomic mass constant to kilogram
        real*8,parameter:: ceng  = 4.35975d-18    ! Energy (J/au)
        real*8,parameter:: au_to_a0 = 0.5291772d0 ! A/au
        real*8,parameter:: autocmwave = 2.19474d5 ! Energy (cm^-1/au)
!
        real*8,parameter:: kb = 1.38066d-23    ! Boltzman constant (J/K)
        real*8,Parameter:: pi = 3.141592653589793d0
        real*8,parameter:: hbar = 1.0d0
        complex*16,parameter:: c1 = (0.0d0,1.0d0)
!        
        integer,parameter:: nps = 4
        real*8:: fnps
!
        real*8:: fmass(3*nps)
!
    end Module NH3_INFO
!
! ============================================================================
!
! ================================================================================
!
!************************************************************************
!> module that contains subroutines for PES and its derivatives for NH3.
!! self-contained, no need external module.
!! modified from pes_nh3.f90 and related input file nh3_cbs_ss_5.inp
!! which are in the supplymentary materials of [Yurchenko2005].
!! @see [Yurchenko2005] J. Chem. Phys. 123, 134308 (2005).
!! @author zhangzj
!! @date 2015/3
!************************************************************************
module nh3_mod

   implicit none
   ! constants
   real(8), parameter :: PI = 3.141592653589793d0
   real(8), parameter :: SQRT2=sqrt(2.d0)
   real(8), parameter :: SQRT3=sqrt(3.d0)
   real(8), parameter :: SQRT6=sqrt(6.d0)
 
   ! status of whether potential parameters are saved
   logical, save :: param_saved = .false.
   ! potential parameters to be read from file and saved
   double precision, save :: &      
      Rhoedg,re14,aa1,ve  ,  &
      f1a,f2a,f3a,f4a,f5a,f6a,f7a,f8a, &
      f1a1,f2a1,f3a1,f4a1,f5a1,f6a1,  &
      f0a11,f1a11,f2a11,f3a11,f4a11, &
      f0a12,f1a12,f2a12,f3a12,f4a12, &
      f0a14,f1a14,f2a14,f3a14,f4a14, &
      f0a44,f1a44,f2a44,f3a44,f4a44, &
      f0a111,f1a111,f2a111,f3a111  , &
      f0a112,f1a112,f2a112,f3a112  , &
      f0a114,f1a114,f2a114,f3a114  , &
      f0a123,f1a123,f2a123,f3a123  , &
      f0a124,f1a124,f2a124,f3a124  , &
      f0a144,f1a144,f2a144,f3a144  , &
      f0a155,f1a155,f2a155,f3a155  , &
      f0a455,f1a455,f2a455,f3a455  , &
      f0a1111,f1a1111,f2a1111      , &
      f0a1112,f1a1112,f2a1112      , &
      f0a1114,f1a1114,f2a1114      , &
      f0a1122,f1a1122,f2a1122      , &
      f0a1123,f1a1123,f2a1123      , &
      f0a1124,f1a1124,f2a1124      , &
      f0a1125,f1a1125,f2a1125      , &
      f0a1144,f1a1144,f2a1144      , &
      f0a1155,f1a1155,f2a1155      , &
      f0a1244,f1a1244,f2a1244      , &
      f0a1255,f1a1255,f2a1255      , &
      f0a1444,f1a1444,f2a1444      , &
      f0a1455,f1a1455,f2a1455      , &
      f0a4444,f1a4444,f2a4444      , &
      f0a44444 ,f1a44444 ,           &
      f2a44444 ,f0a33455 ,f1a33455 ,f2a33455 ,f0a33445 ,f1a33445 ,&
      f2a33445 ,f0a33345 ,f1a33345 ,f2a33345 ,f0a33344 ,f1a33344 ,&
      f2a33344 ,f0a33334 ,f1a33334 ,f2a33334 ,f0a33333 ,f1a33333 ,&
      f2a33333 ,f0a25555 ,f1a25555 ,f2a25555 ,f0a24455 ,f1a24455 ,&
      f2a24455 ,f0a24445 ,f1a24445 ,f2a24445 ,f0a23333 ,f1a23333 ,&
      f2a23333 ,f0a13455 ,f1a13455 ,f2a13455 ,f0a13445 ,f1a13445 ,&
      f2a13445 ,f0a13345 ,f1a13345 ,f2a13345 ,f0a12355 ,f1a12355 ,&
      f2a12355 ,f0a11334 ,f1a11334 ,f2a11334 ,f0a11333 ,f1a11333 ,&
      f2a11333 ,f0a11255 ,f1a11255 ,f2a11255 ,f0a11245 ,f1a11245 ,&
      f2a11245 ,f0a11234 ,f1a11234 ,f2a11234 ,f0a11233 ,f1a11233 ,&
      f2a11233 ,f0a11135 ,f1a11135 ,f2a11135 ,f0a11134 ,f1a11134 ,&
      f2a11134 ,f0a11123 ,f1a11123 ,f2a11123 ,f0a555555,f1a555555,&
      f2a555555,f0a444444,f1a444444,f2a444444,f0a335555,f1a335555,&
      f2a335555,f0a334455,f1a334455,f2a334455,f0a334445,f1a334445,&
      f2a334445,f0a333555,f1a333555,f2a333555,f0a333333,f1a333333,&
      f2a333333,f0a244555,f1a244555,f2a244555,f0a244455,f1a244455,&
      f2a244455,f0a233445,f1a233445,f2a233445,f0a233444,f1a233444,&
      f2a233444,f0a233345,f1a233345,f2a233345,f0a233344,f1a233344,&
      f2a233344,f0a233335,f1a233335,f2a233335,f0a223355,f1a223355,&
      f2a223355,f0a222335,f1a222335,f2a222335,f0a222334,f1a222334,&
      f2a222334,f0a222333,f1a222333,f2a222333,f0a222255,f1a222255,&
      f2a222255,f0a222245,f1a222245,f2a222245,f0a222233,f1a222233,&
      f2a222233,f0a222224,f1a222224,f2a222224,f0a145555,f1a145555,&
      f2a145555,f0a134444,f1a134444,f2a134444,f0a133444,f1a133444,&
      f2a133444,f0a133345,f1a133345,f2a133345,f0a133334,f1a133334,&
      f2a133334,f0a133333,f1a133333,f2a133333,f0a124555,f1a124555,&
      f2a124555,f0a124455,f1a124455,f2a124455,f0a123455,f1a123455,&
      f2a123455,f0a123345,f1a123345,f2a123345,f0a113555,f1a113555,&
      f2a113555,f0a113345,f1a113345,f2a113345,f0a112355,f1a112355,&
      f2a112355,f0a112335,f1a112335,f2a112335,f0a112233,f1a112233,&
      f2a112233,f0a111444,f1a111444,f2a111444,f0a111234,f1a111234,&
      f2a111234,f0a111233,f1a111233,f2a111233,f0a111123,f1a111123,&
      f2a111123
contains
   !---------------------------------------------------------------------
   !> calculate potential and derivatives with cartesian coordinates
   !! @see [Yurchenko2005] J. Chem. Phys. 123, 134308 (2005).
   !! @param [out] v, energy, in cm-1
   !!              output, convert it back to a.u.
   !! @param [out] dv, gradients, units are cm-1, angstrom and radian
   !!              output, convert it back to a.u.
   !! @param [out] ddv, hessian, units are cm-1, angstrom and radian
   !!              output, convert it back to a.u.
   !! @param [in] cart, cartesian coordinates
   !! @param [in] flag, control out values, 0: v, 1: v,dv, 2: v,dv,ddv
   !---------------------------------------------------------------------
   subroutine nh3_pot_cart(v, dvdx, ddvdxdx, cart, flag) bind(c, name="nh3_pot_cart_ccall")
      use NH3_INFO, only : au_to_a0,  autocmwave
      implicit none
      real(8), intent(out) :: v, dvdx(12), ddvdxdx(12,12)
      real(8), intent(inout)  :: cart(12)
      integer, intent(in)  :: flag
 
      real(8) :: r(6), drdx(6,12), ddrdxdx(6,12,12)
      real(8) :: dvdr(6), ddvdrdr(6,6)
 
      if (flag < 0) then
         write(*,*) 'flag should be 0 for v, or 1 for v and dv, or 2 for v,dv and ddv.'
         stop
      end if
 
      ! convert a.u. to Angstrom, length
      cart = cart * au_to_a0     
 
      ! convert cartesian to internal coordinates
      call nh3_cart_to_inter(r, drdx, ddrdxdx, cart, flag) 
 
      ! calculate potential and derivatives with internal coordinates
      call nh3_pot_inter(v, dvdr, ddvdrdr, r, flag)
 
      ! derivatives
      if (flag > 0) then
         call calc_dvdx_from_dvdy(dvdx, ddvdxdx, dvdr, ddvdrdr, drdx, ddrdxdx, 12, 6, flag)
      end if
!
      ! convert Angstrom to a.u., length
      cart = cart / au_to_a0
      if (flag == 0)then
         v = v / autocmwave   ! convert cm^-1 to a.u.
      elseif(flag == 1)then
         v = v / autocmwave   ! convert cm^-1 to a.u.
         dvdx = dvdx / autocmwave * au_to_a0
      elseif(flag == 2)then 
         v = v / autocmwave   ! convert cm^-1 to a.u.
         dvdx = dvdx / autocmwave * au_to_a0
         ddvdxdx = ddvdxdx / autocmwave * au_to_a0**2
      endif
!
   end subroutine nh3_pot_cart
!
   !---------------------------------------------------------------------
   !> calculate potential and derivatives with internal coordinates
   !! @see [Yurchenko2005] J. Chem. Phys. 123, 134308 (2005).
   !! @param [out] v, energy, in cm-1
   !! @param [out] dvdr, gradients, units are cm-1, angstrom and radian
   !! @param [out] ddvdrdr, hessian, units are cm-1, angstrom and radian
   !! @param [in] r, internal coordinates
   !! @param [in] flag, control out values, 0: v, 1: v,dv, 2: v,dv,ddv
   !---------------------------------------------------------------------
   subroutine nh3_pot_inter(v, dvdr, ddvdrdr, r, flag)
      implicit none
      real(8), intent(in)  :: r(6)
      integer, intent(in)  :: flag
      real(8), intent(out) :: v, dvdr(6), ddvdrdr(6,6)
 
      real(8) :: work(6), dvdw(6), ddvdwdw(6,6)
      real(8) :: dwdr(6,6), ddwdrdr(6,6,6)
 
      if (flag < 0) then
         write(*,*) 'flag should be 0 for v, or 1 for v and dv, or 2 for v,dv and ddv.'
         stop
      end if
 
      ! convert internal coordinates to working coordinates
      call nh3_inter_to_work(work, dwdr, ddwdrdr, r, flag)
 
      ! calculate potential energy and derivatives with working coordinates
      call nh3_pot_work(v,dvdw,ddvdwdw,work,flag)
 
      ! derivatives
      if (flag > 0) then
         call calc_dvdx_from_dvdy(dvdr, ddvdrdr, dvdw, ddvdwdw, dwdr, ddwdrdr, 6, 6, flag)
      end if
 
      return
   end subroutine nh3_pot_inter 
 
   !---------------------------------------------------------------------
   !> calculate potential and derivatives with working coordinates
   !! @param [out] v, potential energy, in cm-1
   !! @param [out] dv, 1st derivatives of v to working coordinates
   !! @param [out] ddv, 2nd derivatives of v to working coordinates
   !! @param [in] work, working coordinates (y1,y2,y3,y4,y5,y6)
   !! @param [in] flag, control out values, 0: v, 1: v,dv, 2: v,dv,ddv
   !---------------------------------------------------------------------
   subroutine nh3_pot_work(v,dv,ddv,work,flag)
      implicit none
      real(8), intent(out) :: v, dv(6), ddv(6,6)
      real(8), intent(in)  :: work(6)
      integer, intent(in)  :: flag
 
      real(8) :: y1,y2,y3,y4,y5,coro
      real(8) :: v0,dv0,ddv0
      real(8) :: v1,v2,v3,v4,v5,v6
      real(8) :: s1,s2,s3,s4,s5
      real(8),dimension(6) :: dv1,dv2,dv3,dv4,dv5,dv6
      real(8),dimension(6) :: ds1,ds2,ds3,ds4,ds5
      real(8),dimension(6,6) :: ddv1,ddv2,ddv3,ddv4,ddv5,ddv6
      real(8),dimension(6,6) :: dds1,dds2,dds3,dds4,dds5
 
      ! sinrho functions multiplied with parameter
      double precision   :: &
         fea1  ,                                    &
         fea11  ,fea12  ,fea14  ,fea44  ,                   &
         fea111 ,fea112 ,fea114 ,fea123 ,                   &
         fea124 ,fea144 ,fea155 ,fea455 ,                   &
         fea1111,fea1112,fea1114,fea1122,                   &
         fea1123,fea1124,fea1125,fea1144,                   &
         fea1155,fea1244,fea1255,fea1444,                   &
         fea1455,fea4444
      double precision ::   &
         fea44444 ,fea33455 ,fea33445 ,fea33345 ,fea33344 ,&
         fea33334 ,fea33333 ,fea25555 ,fea24455 ,fea24445 ,fea23333 ,&
         fea13455 ,fea13445 ,fea13345 ,fea12355 ,fea11334 ,fea11333 ,&
         fea11255 ,fea11245 ,fea11234 ,fea11233 ,fea11135 ,fea11134 ,&
         fea11123 ,fea555555,fea444444,fea335555,fea334455,fea334445,&
         fea333555,fea333333,fea244555,fea244455,fea233445,fea233444,&
         fea233345,fea233344,fea233335,fea223355,fea222335,fea222334,&
         fea222333,fea222255,fea222245,fea222233,fea222224,fea145555,&
         fea134444,fea133444,fea133345,fea133334,fea133333,fea124555,&
         fea124455,fea123455,fea123345,fea113555,fea113345,fea112355,&
         fea112335,fea112233,fea111444,fea111234,fea111233,fea111123
 
      double precision   :: &
         dfea1  ,                                    &
         dfea11  ,dfea12  ,dfea14  ,dfea44  ,                   &
         dfea111 ,dfea112 ,dfea114 ,dfea123 ,                   &
         dfea124 ,dfea144 ,dfea155 ,dfea455 ,                   &
         dfea1111,dfea1112,dfea1114,dfea1122,                   &
         dfea1123,dfea1124,dfea1125,dfea1144,                   &
         dfea1155,dfea1244,dfea1255,dfea1444,                   &
         dfea1455,dfea4444
      double precision ::   &
         dfea44444 ,dfea33455 ,dfea33445 ,dfea33345 ,dfea33344 ,&
         dfea33334 ,dfea33333 ,dfea25555 ,dfea24455 ,dfea24445 ,dfea23333 ,&
         dfea13455 ,dfea13445 ,dfea13345 ,dfea12355 ,dfea11334 ,dfea11333 ,&
         dfea11255 ,dfea11245 ,dfea11234 ,dfea11233 ,dfea11135 ,dfea11134 ,&
         dfea11123 ,dfea555555,dfea444444,dfea335555,dfea334455,dfea334445,&
         dfea333555,dfea333333,dfea244555,dfea244455,dfea233445,dfea233444,&
         dfea233345,dfea233344,dfea233335,dfea223355,dfea222335,dfea222334,&
         dfea222333,dfea222255,dfea222245,dfea222233,dfea222224,dfea145555,&
         dfea134444,dfea133444,dfea133345,dfea133334,dfea133333,dfea124555,&
         dfea124455,dfea123455,dfea123345,dfea113555,dfea113345,dfea112355,&
         dfea112335,dfea112233,dfea111444,dfea111234,dfea111233,dfea111123
 
      double precision   :: &
         ddfea1  ,                                    &
         ddfea11  ,ddfea12  ,ddfea14  ,ddfea44  ,                   &
         ddfea111 ,ddfea112 ,ddfea114 ,ddfea123 ,                   &
         ddfea124 ,ddfea144 ,ddfea155 ,ddfea455 ,                   &
         ddfea1111,ddfea1112,ddfea1114,ddfea1122,                   &
         ddfea1123,ddfea1124,ddfea1125,ddfea1144,                   &
         ddfea1155,ddfea1244,ddfea1255,ddfea1444,                   &
         ddfea1455,ddfea4444
      double precision ::   &
         ddfea44444 ,ddfea33455 ,ddfea33445 ,ddfea33345 ,ddfea33344 ,&
         ddfea33334 ,ddfea33333 ,ddfea25555 ,ddfea24455 ,ddfea24445 ,ddfea23333 ,&
         ddfea13455 ,ddfea13445 ,ddfea13345 ,ddfea12355 ,ddfea11334 ,ddfea11333 ,&
         ddfea11255 ,ddfea11245 ,ddfea11234 ,ddfea11233 ,ddfea11135 ,ddfea11134 ,&
         ddfea11123 ,ddfea555555,ddfea444444,ddfea335555,ddfea334455,ddfea334445,&
         ddfea333555,ddfea333333,ddfea244555,ddfea244455,ddfea233445,ddfea233444,&
         ddfea233345,ddfea233344,ddfea233335,ddfea223355,ddfea222335,ddfea222334,&
         ddfea222333,ddfea222255,ddfea222245,ddfea222233,ddfea222224,ddfea145555,&
         ddfea134444,ddfea133444,ddfea133345,ddfea133334,ddfea133333,ddfea124555,&
         ddfea124455,ddfea123455,ddfea123345,ddfea113555,ddfea113345,ddfea112355,&
         ddfea112335,ddfea112233,ddfea111444,ddfea111234,ddfea111233,ddfea111123
 
      integer :: i,j
 
      if (flag < 0) then
         write(*,*) 'flag should be 0 for v, or 1 for v and dv, or 2 for v,dv and ddv.'
         stop
      end if
 
      y1=work(1)
      y2=work(2)
      y3=work(3)
      y4=work(4)
      y5=work(5)
      coro=work(6)
 
      ! set parameters (only for the first time) 
      call set_param()
 
      ! sinrho functions 
 
      fea1= f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4+f5a1*coro**5+f6a1*coro**6
 
      fea11=   f0a11+f1a11*coro+f2a11*coro**2+f3a11*coro**3+f4a11*coro**4
      fea12=   f0a12+f1a12*coro+f2a12*coro**2+f3a12*coro**3+f4a12*coro**4
      fea14=   f0a14+f1a14*coro+f2a14*coro**2+f3a14*coro**3+f4a14*coro**4
      fea44=   f0a44+f1a44*coro+f2a44*coro**2+f3a44*coro**3+f4a44*coro**4
 
      fea111= f0a111+f1a111*coro+f2a111*coro**2+f3a111*coro**3
      fea112= f0a112+f1a112*coro+f2a112*coro**2+f3a112*coro**3
      fea114= f0a114+f1a114*coro+f2a114*coro**2+f3a114*coro**3
      fea123= f0a123+f1a123*coro+f2a123*coro**2+f3a123*coro**3
      fea124= f0a124+f1a124*coro+f2a124*coro**2+f3a124*coro**3
      fea144= f0a144+f1a144*coro+f2a144*coro**2+f3a144*coro**3
      fea155= f0a155+f1a155*coro+f2a155*coro**2+f3a155*coro**3
      fea455= f0a455+f1a455*coro+f2a455*coro**2+f3a455*coro**3
 
      fea1111= f0a1111+f1a1111*coro+f2a1111*coro**2
      fea1112= f0a1112+f1a1112*coro+f2a1112*coro**2
      fea1114= f0a1114+f1a1114*coro+f2a1114*coro**2
      fea1122= f0a1122+f1a1122*coro+f2a1122*coro**2
      fea1123= f0a1123+f1a1123*coro+f2a1123*coro**2
      fea1124= f0a1124+f1a1124*coro+f2a1124*coro**2
      fea1125= f0a1125+f1a1125*coro+f2a1125*coro**2
      fea1144= f0a1144+f1a1144*coro+f2a1144*coro**2
      fea1155= f0a1155+f1a1155*coro+f2a1155*coro**2
      fea1244= f0a1244+f1a1244*coro+f2a1244*coro**2
      fea1255= f0a1255+f1a1255*coro+f2a1255*coro**2
      fea1444= f0a1444+f1a1444*coro+f2a1444*coro**2
      fea1455= f0a1455+f1a1455*coro+f2a1455*coro**2
      fea4444= f0a4444+f1a4444*coro+f2a4444*coro**2
 
      fea44444 = f0a44444  + f1a44444 *coro+ f2a44444 *coro**2
      fea33455 = f0a33455  + f1a33455 *coro+ f2a33455 *coro**2
      fea33445 = f0a33445  + f1a33445 *coro+ f2a33445 *coro**2
      fea33345 = f0a33345  + f1a33345 *coro+ f2a33345 *coro**2
      fea33344 = f0a33344  + f1a33344 *coro+ f2a33344 *coro**2
      fea33334 = f0a33334  + f1a33334 *coro+ f2a33334 *coro**2
      fea33333 = f0a33333  + f1a33333 *coro+ f2a33333 *coro**2
      fea25555 = f0a25555  + f1a25555 *coro+ f2a25555 *coro**2
      fea24455 = f0a24455  + f1a24455 *coro+ f2a24455 *coro**2
      fea24445 = f0a24445  + f1a24445 *coro+ f2a24445 *coro**2
      fea23333 = f0a23333  + f1a23333 *coro+ f2a23333 *coro**2
      fea13455 = f0a13455  + f1a13455 *coro+ f2a13455 *coro**2
      fea13445 = f0a13445  + f1a13445 *coro+ f2a13445 *coro**2
      fea13345 = f0a13345  + f1a13345 *coro+ f2a13345 *coro**2
      fea12355 = f0a12355  + f1a12355 *coro+ f2a12355 *coro**2
      fea11334 = f0a11334  + f1a11334 *coro+ f2a11334 *coro**2
      fea11333 = f0a11333  + f1a11333 *coro+ f2a11333 *coro**2
      fea11255 = f0a11255  + f1a11255 *coro+ f2a11255 *coro**2
      fea11245 = f0a11245  + f1a11245 *coro+ f2a11245 *coro**2
      fea11234 = f0a11234  + f1a11234 *coro+ f2a11234 *coro**2
      fea11233 = f0a11233  + f1a11233 *coro+ f2a11233 *coro**2
      fea11135 = f0a11135  + f1a11135 *coro+ f2a11135 *coro**2
      fea11134 = f0a11134  + f1a11134 *coro+ f2a11134 *coro**2
      fea11123 = f0a11123  + f1a11123 *coro+ f2a11123 *coro**2
      fea555555= f0a555555 + f1a555555*coro+ f2a555555*coro**2
      fea444444= f0a444444 + f1a444444*coro+ f2a444444*coro**2
      fea335555= f0a335555 + f1a335555*coro+ f2a335555*coro**2
      fea334455= f0a334455 + f1a334455*coro+ f2a334455*coro**2
      fea334445= f0a334445 + f1a334445*coro+ f2a334445*coro**2
      fea333555= f0a333555 + f1a333555*coro+ f2a333555*coro**2
      fea333333= f0a333333 + f1a333333*coro+ f2a333333*coro**2
      fea244555= f0a244555 + f1a244555*coro+ f2a244555*coro**2
      fea244455= f0a244455 + f1a244455*coro+ f2a244455*coro**2
      fea233445= f0a233445 + f1a233445*coro+ f2a233445*coro**2
      fea233444= f0a233444 + f1a233444*coro+ f2a233444*coro**2
      fea233345= f0a233345 + f1a233345*coro+ f2a233345*coro**2
      fea233344= f0a233344 + f1a233344*coro+ f2a233344*coro**2
      fea233335= f0a233335 + f1a233335*coro+ f2a233335*coro**2
      fea223355= f0a223355 + f1a223355*coro+ f2a223355*coro**2
      fea222335= f0a222335 + f1a222335*coro+ f2a222335*coro**2
      fea222334= f0a222334 + f1a222334*coro+ f2a222334*coro**2
      fea222333= f0a222333 + f1a222333*coro+ f2a222333*coro**2
      fea222255= f0a222255 + f1a222255*coro+ f2a222255*coro**2
      fea222245= f0a222245 + f1a222245*coro+ f2a222245*coro**2
      fea222233= f0a222233 + f1a222233*coro+ f2a222233*coro**2
      fea222224= f0a222224 + f1a222224*coro+ f2a222224*coro**2
      fea145555= f0a145555 + f1a145555*coro+ f2a145555*coro**2
      fea134444= f0a134444 + f1a134444*coro+ f2a134444*coro**2
      fea133444= f0a133444 + f1a133444*coro+ f2a133444*coro**2
      fea133345= f0a133345 + f1a133345*coro+ f2a133345*coro**2
      fea133334= f0a133334 + f1a133334*coro+ f2a133334*coro**2
      fea133333= f0a133333 + f1a133333*coro+ f2a133333*coro**2
      fea124555= f0a124555 + f1a124555*coro+ f2a124555*coro**2
      fea124455= f0a124455 + f1a124455*coro+ f2a124455*coro**2
      fea123455= f0a123455 + f1a123455*coro+ f2a123455*coro**2
      fea123345= f0a123345 + f1a123345*coro+ f2a123345*coro**2
      fea113555= f0a113555 + f1a113555*coro+ f2a113555*coro**2
      fea113345= f0a113345 + f1a113345*coro+ f2a113345*coro**2
      fea112355= f0a112355 + f1a112355*coro+ f2a112355*coro**2
      fea112335= f0a112335 + f1a112335*coro+ f2a112335*coro**2
      fea112233= f0a112233 + f1a112233*coro+ f2a112233*coro**2
      fea111444= f0a111444 + f1a111444*coro+ f2a111444*coro**2
      fea111234= f0a111234 + f1a111234*coro+ f2a111234*coro**2
      fea111233= f0a111233 + f1a111233*coro+ f2a111233*coro**2
      fea111123= f0a111123 + f1a111123*coro+ f2a111123*coro**2
 
      ! calculate potential
      v0=-f1a*coro+f2a*coro**2+f3a*coro**3+f4a*coro**4+f5a*coro**5 &
         +f6a*coro**6+f7a*coro**7+f8a*coro**8
 
      v1 = (y3+y2+y1)*FEA1
 
      v2 = (y2*y3+y1*y3+y1*y2)*FEA12                                                    &
         +(y2**2+y3**2+y1**2)*FEA11                                                        &
         +(-sqrt(3.D0)*y3*y5/2.D0-y3*y4/2.D0+y1*y4+sqrt(3.D0)*y2*y5/2.D0-y2*y4/2.D0)*FEA14 &
         +(y5**2+y4**2)*FEA44
 
      v3 = (y1*y3*y4+y1*y2*y4-2.D0*y2*y3*y4+sqrt(3.D0)*y1*y2*y5-sqrt(3.D0)*y1*y3*y5)*FEA124        &
         +(3.D0/4.D0*y3*y4**2-sqrt(3.D0)*y3*y4*y5/2.D0+y1*y5**2+y2*y5**2/4.D0+3.D0/4.D0*y2*y4**2 &
         +sqrt(3.D0)*y2*y4*y5/2.D0+y3*y5**2/4.D0)*FEA155  &
         +(y2*y3**2+y1*y3**2+y1**2*y3+y1*y2**2+y2**2*y3+y1**2*y2)*FEA112+                             &
         (-y4**3/3.D0+y4*y5**2)*FEA455+FEA123*y1*y2*y3                                                &
         +(y1*y4**2+3.D0/4.D0*y3*y5**2+3.D0/4.D0*y2*y5**2+y2*y4**2/4.D0-sqrt(3.D0)*y2*y4*y5/2.D0 &
         +sqrt(3.D0)*y3*y4*y5/2.D0+y3*y4**2/4.D0)*FEA144 &
         +(y3**3+y2**3+y1**3)*FEA111                                                                  &
         +(-y2**2*y4/2.D0-y3**2*y4/2.D0+sqrt(3.D0)*y2**2*y5/2.D0+y1**2*y4-sqrt(3.D0)*y3**2*y5/2.D0)*FEA114
 

      s2 = (y4**4+y5**4+2.D0*y4**2*y5**2)*FEA4444+(3.D0/8.D0*sqrt(3.D0)*&
         &y2*y5**3-3.D0/8.D0*sqrt(3.D0)*y3*y4**2*y5-3.D0/8.D0*sqrt(3.D0)*y3*&
         &y5**3-9.D0/8.D0*y2*y4*y5**2-y3*y4**3/8.D0-y2*y4**3/8.D0-9.D0/8.D0*&
         &y3*y4*y5**2+y1*y4**3+3.D0/8.D0*sqrt(3.D0)*y2*y4**2*y5)*FEA1444+(3.&
         &D0/4.D0*y2**2*y4**2+3.D0/4.D0*y3**2*y4**2+y1**2*y5**2+y3**2*y5**2/&
         &4.D0-sqrt(3.D0)*y3**2*y4*y5/2.D0+sqrt(3.D0)*y2**2*y4*y5/2.D0+y2**2&
         &*y5**2/4.D0)*FEA1155 
      s1 = s2+(y3**2*y4**2/4.D0+3.D0/4.D0*y3**2*y5**2+y1**2*y4**2+y2**2*&
         &y4**2/4.D0+sqrt(3.D0)*y3**2*y4*y5/2.D0-sqrt(3.D0)*y2**2*y4*y5/2.D0&
         &+3.D0/4.D0*y2**2*y5**2)*FEA1144+(y1**3*y4+sqrt(3.D0)*y2**3*y5/2.D0&
         &-sqrt(3.D0)*y3**3*y5/2.D0-y2**3*y4/2.D0-y3**3*y4/2.D0)*FEA1114+(y2&
         &**4+y1**4+y3**4)*FEA1111+(sqrt(3.D0)*y1*y3*y4*y5+3.D0/2.D0*y2*y3*y&
         &5**2-y2*y3*y4**2/2.D0+y1*y2*y4**2-sqrt(3.D0)*y1*y2*y4*y5+y1*y3*y4*&
         &*2)*FEA1244 
      s2 = s1+(y1*y3*y5**2+y1*y2*y5**2-sqrt(3.D0)*y1*y3*y4*y5-y2*y3*y5**&
         &2/2.D0+3.D0/2.D0*y2*y3*y4**2+sqrt(3.D0)*y1*y2*y4*y5)*FEA1255+(-y1*y&
         &3**2*y4/2.D0+y1**2*y3*y4-sqrt(3.D0)*y1*y3**2*y5/2.D0-sqrt(3.D0)*y2&
         &*y3**2*y5/2.D0+y1**2*y2*y4+sqrt(3.D0)*y2**2*y3*y5/2.D0-y2**2*y3*y4&
         &/2.D0+sqrt(3.D0)*y1*y2**2*y5/2.D0-y2*y3**2*y4/2.D0-y1*y2**2*y4/2.D&
         &0)*FEA1124+(y1**2*y2*y5+sqrt(3.D0)*y1*y3**2*y4/2.D0+sqrt(3.D0)*y1*&
         &y2**2*y4/2.D0-sqrt(3.D0)*y2*y3**2*y4/2.D0-sqrt(3.D0)*y2**2*y3*y4/2&
         &.D0-y2**2*y3*y5/2.D0+y2*y3**2*y5/2.D0-y1*y3**2*y5/2.D0+y1*y2**2*y5&
         &/2.D0-y1**2*y3*y5)*FEA1125 
      v4 = s2+(y2*y3**3+y1**3*y3+y1**3*y2+y1*y2**3+y1*y3**3+y2**3*y3)*FE&
         &A1112+(y2**2*y3**2+y1**2*y3**2+y1**2*y2**2)*FEA1122+(y1*y2**2*y3+y&
         &1**2*y2*y3+y1*y2*y3**2)*FEA1123+(5.D0/8.D0*y2*y4*y5**2+sqrt(3.D0)*&
         &y2*y5**3/8.D0-sqrt(3.D0)*y3*y4**2*y5/8.D0+sqrt(3.D0)*y2*y4**2*y5/8&
         &.D0-3.D0/8.D0*y2*y4**3+y1*y4*y5**2-sqrt(3.D0)*y3*y5**3/8.D0+5.D0/8&
         &.D0*y3*y4*y5**2-3.D0/8.D0*y3*y4**3)*FEA1455
 

      s3 = (y4**5-2.D0*y4**3*y5**2-3.D0*y4*y5**4)*FEA44444+(-4.D0*y3*y4*&
         &y5**3*sqrt(3.D0)+9.D0*y1*y4**2*y5**2-3.D0/2.D0*y1*y4**4+4.D0*y2*y4&
         &*y5**3*sqrt(3.D0)+3.D0*y2*y4**4+5.D0/2.D0*y1*y5**4+3.D0*y3*y4**4+y&
         &2*y5**4+y3*y5**4)*FEA25555+(-y2*y4**4+y3*y4**2*y5**2-2.D0*y2*y4*y5&
         &**3*sqrt(3.D0)-y3*y4**4-7.D0/2.D0*y1*y4**2*y5**2-3.D0/4.D0*y1*y5**&
         &4+2.D0*y3*y4*y5**3*sqrt(3.D0)+y2*y4**2*y5**2+5.D0/4.D0*y1*y4**4)*F&
         &EA24455 
      s2 = s3+(y2*y4**3*y5-3.D0*y3*y4*y5**3+2.D0/3.D0*y3*y4**4*sqrt(3.D0&
         &)+3.D0/4.D0*y1*y5**4*sqrt(3.D0)+3.D0*y2*y4*y5**3-7.D0/12.D0*y1*y4*&
         &*4*sqrt(3.D0)+3.D0/2.D0*y1*y4**2*y5**2*sqrt(3.D0)-y3*y4**3*y5+2.D0&
         &/3.D0*y2*y4**4*sqrt(3.D0))*FEA24445+(-y2**2*y5**3+y3**2*y4**2*y5+y&
         &3**2*y5**3+4.D0/9.D0*y2**2*y4**3*sqrt(3.D0)-5.D0/9.D0*y1**2*y4**3*&
         &sqrt(3.D0)+4.D0/9.D0*y3**2*y4**3*sqrt(3.D0)-y2**2*y4**2*y5-y1**2*y&
         &4*y5**2*sqrt(3.D0))*FEA33445+(y3**2*y4*y5**2-y1**2*y4**3/3.D0-y3**&
         &2*y4**3/3.D0+y1**2*y4*y5**2+y2**2*y4*y5**2-y2**2*y4**3/3.D0)*FEA33&
         &455 
      s1 = s2+(-y2**3*y4*y5+y3**3*y4*y5+y2**3*y5**2*sqrt(3.D0)/3.D0+y1**&
         &3*y4**2*sqrt(3.D0)/2.D0+y3**3*y5**2*sqrt(3.D0)/3.D0-y1**3*y5**2*sq&
         &rt(3.D0)/6.D0)*FEA33345+(y3**3*y4**2+y3**3*y5**2+y2**3*y4**2+y2**3&
         &*y5**2+y1**3*y5**2+y1**3*y4**2)*FEA33344+(y3**4*y4+sqrt(3.D0)*y3**&
         &4*y5+y2**4*y4-2.D0*y1**4*y4-sqrt(3.D0)*y2**4*y5)*FEA33334+(y2**5+y&
         &3**5+y1**5)*FEA33333+(-4.D0/9.D0*y1*y2*y4**3*sqrt(3.D0)-y1*y2*y5**&
         &3+y1*y3*y4**2*y5+y2*y3*y4*y5**2*sqrt(3.D0)-y1*y2*y4**2*y5+5.D0/9.D&
         &0*y2*y3*y4**3*sqrt(3.D0)-4.D0/9.D0*y1*y3*y4**3*sqrt(3.D0)+y1*y3*y5&
         &**3)*FEA13445+(y2*y3*y4*y5**2+y1*y2*y4*y5**2-y2*y3*y4**3/3.D0-y1*y&
         &2*y4**3/3.D0-y1*y3*y4**3/3.D0+y1*y3*y4*y5**2)*FEA13455 
      s3 = s1+(y1**2*y3*y5**2+y2**2*y3*y4**2+y2**2*y3*y5**2+y1*y2**2*y5*&
         &*2+y1**2*y2*y5**2+y1*y2**2*y4**2+y2*y3**2*y4**2+y1*y3**2*y4**2+y1*&
         &*2*y3*y4**2+y1**2*y2*y4**2+y1*y3**2*y5**2+y2*y3**2*y5**2)*FEA11255&
         &+(2.D0/3.D0*y1**2*y3*y4**2*sqrt(3.D0)+y1*y3**2*y5**2*sqrt(3.D0)/2.&
         &D0+y1*y2**2*y5**2*sqrt(3.D0)/2.D0+y2**2*y3*y5**2*sqrt(3.D0)/2.D0-y&
         &1*y2**2*y4*y5+y2*y3**2*y4*y5+y1*y3**2*y4*y5-y2**2*y3*y4*y5+y2*y3**&
         &2*y4**2*sqrt(3.D0)/6.D0+y1*y3**2*y4**2*sqrt(3.D0)/6.D0+y1*y2**2*y4&
         &**2*sqrt(3.D0)/6.D0+2.D0/3.D0*y1**2*y2*y4**2*sqrt(3.D0)+y2*y3**2*y&
         &5**2*sqrt(3.D0)/2.D0+y2**2*y3*y4**2*sqrt(3.D0)/6.D0)*FEA13345 
      s4 = s3+(y1**2*y2*y4*y5+y1**2*y3*y4**2*sqrt(3.D0)/3.D0+y1**2*y2*y4&
         &**2*sqrt(3.D0)/3.D0-y1*y2**2*y4**2*sqrt(3.D0)/6.D0+y2*y3**2*y4*y5-&
         &y2**2*y3*y4*y5-y1**2*y3*y4*y5+y2*y3**2*y4**2*sqrt(3.D0)/3.D0+y1*y2&
         &**2*y5**2*sqrt(3.D0)/2.D0-y1*y3**2*y4**2*sqrt(3.D0)/6.D0+y2**2*y3*&
         &y4**2*sqrt(3.D0)/3.D0+y1*y3**2*y5**2*sqrt(3.D0)/2.D0)*FEA11245 
      s2 = s4+(-y1**3*y2*y5+y1**3*y3*y5+y2**3*y3*y5/2.D0-y1*y2**3*y4*sqr&
         &t(3.D0)/2.D0-y1*y2**3*y5/2.D0-y2*y3**3*y5/2.D0+y1*y3**3*y5/2.D0+y2&
         &**3*y3*y4*sqrt(3.D0)/2.D0+y2*y3**3*y4*sqrt(3.D0)/2.D0-y1*y3**3*y4*&
         &sqrt(3.D0)/2.D0)*FEA11135+(y1**3*y3*y4-y2**3*y3*y4/2.D0+y1**3*y2*y&
         &4-y2*y3**3*y4/2.D0-y1*y3**3*y4/2.D0+y1*y2**3*y5*sqrt(3.D0)/2.D0+y2&
         &**3*y3*y5*sqrt(3.D0)/2.D0-y2*y3**3*y5*sqrt(3.D0)/2.D0-y1*y2**3*y4/&
         &2.D0-y1*y3**3*y5*sqrt(3.D0)/2.D0)*FEA11134 
      v5 = s2+(y1*y2**4+y1**4*y3+y1**4*y2+y2**4*y3+y2*y3**4+y1*y3**4)*FE&
         &A23333+(-2.D0*y2**2*y3**2*y4+y1**2*y2**2*y4-sqrt(3.D0)*y1**2*y3**2&
         &*y5+sqrt(3.D0)*y1**2*y2**2*y5+y1**2*y3**2*y4)*FEA11334+(y1**2*y3**&
         &3+y1**3*y3**2+y2**2*y3**3+y1**2*y2**3+y1**3*y2**2+y2**3*y3**2)*FEA&
         &11333+(y1*y2*y3*y4**2+y1*y2*y3*y5**2)*FEA12355+(-y1*y2*y3**2*y4/2.&
         &D0-y1*y2**2*y3*y4/2.D0-sqrt(3.D0)*y1*y2*y3**2*y5/2.D0+y1**2*y2*y3*&
         &y4+sqrt(3.D0)*y1*y2**2*y3*y5/2.D0)*FEA11234+(y1*y2**3*y3+y1*y2*y3*&
         &*3+y1**3*y2*y3)*FEA11123+(y1**2*y2**2*y3+y1*y2**2*y3**2+y1**2*y2*y&
         &3**2)*FEA11233
 
      s3 = (y2**3*y4**3*sqrt(3.D0)-y2**3*y4**2*y5+y3**3*y4**2*y5-5.D0/3.&
         &D0*y2**3*y4*y5**2*sqrt(3.D0)+y3**3*y4**3*sqrt(3.D0)-5.D0/3.D0*y3**&
         &3*y4*y5**2*sqrt(3.D0)-y2**3*y5**3+y3**3*y5**3-8.D0/3.D0*y1**3*y4*y&
         &5**2*sqrt(3.D0))*FEA333555+(y1**4*y5**2*sqrt(3.D0)/2.D0+y2**4*y4*y&
         &5+y2**4*y4**2*sqrt(3.D0)/3.D0+y3**4*y4**2*sqrt(3.D0)/3.D0-y3**4*y4&
         &*y5-y1**4*y4**2*sqrt(3.D0)/6.D0)*FEA222245+(y1*y3**5+y1*y2**5+y2**&
         &5*y3+y1**5*y3+y1**5*y2+y2*y3**5)*FEA133333+(y1**4*y3*y4-2.D0*y2**4&
         &*y3*y4+y1**4*y2*y4+y1*y2**4*y5*sqrt(3.D0)+y1*y3**4*y4-2.D0*y2*y3**&
         &4*y4+y1**4*y2*y5*sqrt(3.D0)-y1*y3**4*y5*sqrt(3.D0)-y1**4*y3*y5*sqr&
         &t(3.D0)+y1*y2**4*y4)*FEA133334+(-y1*y2*y3*y4**3/3.D0+y1*y2*y3*y4*y&
         &5**2)*FEA123455 
      s4 = s3+(2.D0/3.D0*sqrt(3.D0)*y1*y2**2*y3**2*y4-y1**2*y2**2*y3*y5-&
         &sqrt(3.D0)*y1**2*y2**2*y3*y4/3.D0+y1**2*y2*y3**2*y5-sqrt(3.D0)*y1*&
         &*2*y2*y3**2*y4/3.D0)*FEA112335+(y1*y2**2*y3*y5**2+y1*y2*y3**2*y5**&
         &2+y1*y2*y3**2*y4**2+y1*y2**2*y3*y4**2+y1**2*y2*y3*y4**2+y1**2*y2*y&
         &3*y5**2)*FEA112355 
      s2 = s4+(y2**3*y3**2*y5-y1**3*y2**2*y5/2.D0-y1**2*y3**3*y5/2.D0-y2&
         &**2*y3**3*y5+y1**3*y2**2*y4*sqrt(3.D0)/2.D0-y1**2*y2**3*y4*sqrt(3.&
         &D0)/2.D0+y1**3*y3**2*y5/2.D0+y1**2*y2**3*y5/2.D0+y1**3*y3**2*y4*sq&
         &rt(3.D0)/2.D0-y1**2*y3**3*y4*sqrt(3.D0)/2.D0)*FEA222335+(-y1**2*y2&
         &**2*y5**2*sqrt(3.D0)/2.D0-y1**2*y3**2*y5**2*sqrt(3.D0)/2.D0-y1**2*&
         &y2**2*y4**2*sqrt(3.D0)/6.D0-y1**2*y2**2*y4*y5-2.D0/3.D0*y2**2*y3**&
         &2*y4**2*sqrt(3.D0)+y1**2*y3**2*y4*y5-y1**2*y3**2*y4**2*sqrt(3.D0)/&
         &6.D0)*FEA113345+(y2**2*y3**2*y5**2+y2**2*y3**2*y4**2+y1**2*y2**2*y&
         &5**2+y1**2*y3**2*y4**2+y1**2*y3**2*y5**2+y1**2*y2**2*y4**2)*FEA223&
         &355 
      s3 = s2+(y1*y2*y3**2*y4**2*sqrt(3.D0)/6.D0+y1*y2*y3**2*y4*y5+y1*y2&
         &*y3**2*y5**2*sqrt(3.D0)/2.D0+2.D0/3.D0*y1**2*y2*y3*y4**2*sqrt(3.D0&
         &)-y1*y2**2*y3*y4*y5+y1*y2**2*y3*y4**2*sqrt(3.D0)/6.D0+y1*y2**2*y3*&
         &y5**2*sqrt(3.D0)/2.D0)*FEA123345+(-y1**3*y2**2*y5*sqrt(3.D0)/2.D0-&
         &y1**3*y2**2*y4/2.D0-y1**3*y3**2*y4/2.D0-y1**2*y2**3*y4/2.D0+y1**3*&
         &y3**2*y5*sqrt(3.D0)/2.D0-y1**2*y3**3*y4/2.D0+y2**3*y3**2*y4-y1**2*&
         &y2**3*y5*sqrt(3.D0)/2.D0+y2**2*y3**3*y4+y1**2*y3**3*y5*sqrt(3.D0)/&
         &2.D0)*FEA222334+(3.D0*y3**2*y4**4+5.D0/2.D0*y1**2*y5**4+y2**2*y5**&
         &4+3.D0*y2**2*y4**4-4.D0*y3**2*y4*y5**3*sqrt(3.D0)+y3**2*y5**4+9.D0&
         &*y1**2*y4**2*y5**2-3.D0/2.D0*y1**2*y4**4+4.D0*y2**2*y4*y5**3*sqrt(&
         &3.D0))*FEA335555+(y1**3*y2**3+y1**3*y3**3+y2**3*y3**3)*FEA222333 
      s4 = s3+(y3*y4**5/5.D0-y2*y4**4*y5*sqrt(3.D0)/2.D0-2.D0/5.D0*y1*y4&
         &**5-2.D0*y1*y4**3*y5**2-3.D0/10.D0*y2*y5**5*sqrt(3.D0)+y3*y4**3*y5&
         &**2+y3*y4**4*y5*sqrt(3.D0)/2.D0+y2*y4**3*y5**2+3.D0/10.D0*y3*y5**5&
         &*sqrt(3.D0)+y2*y4**5/5.D0)*FEA244455+(y2**5*y4-2.D0*y1**5*y4-sqrt(&
         &3.D0)*y2**5*y5+y3**5*y4+sqrt(3.D0)*y3**5*y5)*FEA222224 
      s5 = s4+(-y3*y5**5*sqrt(3.D0)/5.D0+y2*y5**5*sqrt(3.D0)/5.D0+y1*y4*&
         &y5**4-7.D0/15.D0*y2*y4**5+y2*y4**4*y5*sqrt(3.D0)/3.D0-y3*y4**4*y5*&
         &sqrt(3.D0)/3.D0+y3*y4*y5**4+y2*y4*y5**4+2.D0*y1*y4**3*y5**2-7.D0/1&
         &5.D0*y3*y4**5-y1*y4**5/15.D0)*FEA145555 
      s1 = s5+(-sqrt(3.D0)*y1*y2*y3**3*y5/2.D0+y1**3*y2*y3*y4+sqrt(3.D0)&
         &*y1*y2**3*y3*y5/2.D0-y1*y2**3*y3*y4/2.D0-y1*y2*y3**3*y4/2.D0)*FEA1&
         &11234+(y3*y4**4*y5/3.D0+y3*y4**5*sqrt(3.D0)/18.D0-y2*y4**4*y5/3.D0&
         &-y2*y4*y5**4*sqrt(3.D0)/2.D0-y3*y4**2*y5**3+2.D0/9.D0*y1*y4**5*sqr&
         &t(3.D0)+y2*y4**5*sqrt(3.D0)/18.D0+y2*y4**2*y5**3-2.D0/3.D0*y1*y4**&
         &3*y5**2*sqrt(3.D0)-y3*y4*y5**4*sqrt(3.D0)/2.D0)*FEA244555+(y1*y2*y&
         &4**2*y5**2-3.D0/4.D0*y2*y3*y4**4-y1*y2*y5**4-y1*y3*y5**4+5.D0/4.D0&
         &*y2*y3*y5**4+y1*y3*y4**2*y5**2-7.D0/2.D0*y2*y3*y4**2*y5**2-2.D0*y1&
         &*y2*y4**3*y5*sqrt(3.D0)+2.D0*y1*y3*y4**3*y5*sqrt(3.D0))*FEA124455 
      s3 = s1+(y2**6+y1**6+y3**6)*FEA333333+(y1*y2**4*y3+y1**4*y2*y3+y1*&
         &y2*y3**4)*FEA111123+FEA112233*y1**2*y2**2*y3**2+(y1**4*y4**2+y2**4&
         &*y4**2+y2**4*y5**2+y3**4*y4**2+y1**4*y5**2+y3**4*y5**2)*FEA222255 
      s4 = s3+(3.D0*y1*y3*y5**4+y1*y3*y4**4+9.D0*y2*y3*y4**2*y5**2-3.D0/&
         &2.D0*y2*y3*y5**4-4.D0*y1*y3*y4**3*y5*sqrt(3.D0)+y1*y2*y4**4+4.D0*y&
         &1*y2*y4**3*y5*sqrt(3.D0)+3.D0*y1*y2*y5**4+5.D0/2.D0*y2*y3*y4**4)*F&
         &EA134444+(-y1*y3**2*y5**3*sqrt(3.D0)/3.D0-7.D0/3.D0*y1**2*y3*y4*y5&
         &**2+5.D0/3.D0*y1*y2**2*y4**2*y5*sqrt(3.D0)-13.D0/3.D0*y2**2*y3*y4*&
         &y5**2-4.D0/3.D0*y2*y3**2*y5**3*sqrt(3.D0)-7.D0/3.D0*y1**2*y2*y4*y5&
         &**2-16.D0/3.D0*y1*y3**2*y4*y5**2+4.D0/3.D0*y1**2*y3*y4**2*y5*sqrt(&
         &3.D0)+4.D0/3.D0*y2**2*y3*y5**3*sqrt(3.D0)+3.D0*y1**2*y2*y4**3+y2*y&
         &3**2*y4**3+y1*y2**2*y5**3*sqrt(3.D0)/3.D0+y2**2*y3*y4**3-13.D0/3.D&
         &0*y2*y3**2*y4*y5**2-5.D0/3.D0*y1*y3**2*y4**2*y5*sqrt(3.D0)-4.D0/3.&
         &D0*y1**2*y2*y4**2*y5*sqrt(3.D0)+3.D0*y1**2*y3*y4**3-16.D0/3.D0*y1*&
         &y2**2*y4*y5**2)*FEA233444 
      s5 = s4+(2.D0*y1*y3**2*y5**3+4.D0*y2*y3**2*y5**3+4.D0*y2**2*y3*y4*&
         &y5**2*sqrt(3.D0)-2.D0*y1*y2**2*y5**3+y1**2*y3*y4*y5**2*sqrt(3.D0)+&
         &6.D0*y1*y3**2*y4**2*y5-6.D0*y1*y2**2*y4**2*y5-3.D0*y1**2*y3*y4**2*&
         &y5+y1**2*y2*y4*y5**2*sqrt(3.D0)+4.D0*y1*y3**2*y4*y5**2*sqrt(3.D0)-&
         &3.D0*y1**2*y2*y4**3*sqrt(3.D0)-4.D0*y2**2*y3*y5**3+3.D0*y1**2*y2*y&
         &4**2*y5-y1**2*y2*y5**3+y1**2*y3*y5**3-3.D0*y1**2*y3*y4**3*sqrt(3.D&
         &0)+4.D0*y2*y3**2*y4*y5**2*sqrt(3.D0)+4.D0*y1*y2**2*y4*y5**2*sqrt(3&
         &.D0))*FEA113555 
      s2 = s5+(-2.D0/3.D0*y3**2*y4**4*sqrt(3.D0)-3.D0/2.D0*y1**2*y4**2*y&
         &5**2*sqrt(3.D0)-3.D0/4.D0*y1**2*y5**4*sqrt(3.D0)-y2**2*y4**3*y5+7.&
         &D0/12.D0*y1**2*y4**4*sqrt(3.D0)+y3**2*y4**3*y5+3.D0*y3**2*y4*y5**3&
         &-2.D0/3.D0*y2**2*y4**4*sqrt(3.D0)-3.D0*y2**2*y4*y5**3)*FEA334445+(&
         &-3.D0*y1*y3*y4**3*y5+2.D0/3.D0*y1*y2*y5**4*sqrt(3.D0)-y1*y3*y4*y5*&
         &*3+2.D0/3.D0*y1*y3*y5**4*sqrt(3.D0)+3.D0*y1*y2*y4**3*y5-7.D0/12.D0&
         &*y2*y3*y5**4*sqrt(3.D0)+3.D0/2.D0*y2*y3*y4**2*y5**2*sqrt(3.D0)+y1*&
         &y2*y4*y5**3+3.D0/4.D0*y2*y3*y4**4*sqrt(3.D0))*FEA124555+(2.D0*y3**&
         &2*y4*y5**3*sqrt(3.D0)-7.D0/2.D0*y1**2*y4**2*y5**2+y2**2*y4**2*y5**&
         &2-y2**2*y4**4-y3**2*y4**4-2.D0*y2**2*y4*y5**3*sqrt(3.D0)-3.D0/4.D0&
         &*y1**2*y5**4+5.D0/4.D0*y1**2*y4**4+y3**2*y4**2*y5**2)*FEA334455 
      s3 = s2+(-6.D0*y4**2*y5**4+9.D0*y4**4*y5**2+y5**6)*FEA555555+(y2*y&
         &3**3*y4**2+y2*y3**3*y5**2+y1*y3**3*y4**2+y1*y2**3*y4**2+y1**3*y2*y&
         &4**2+y1*y2**3*y5**2+y1**3*y3*y5**2+y1**3*y3*y4**2+y1**3*y2*y5**2+y&
         &2**3*y3*y4**2+y1*y3**3*y5**2+y2**3*y3*y5**2)*FEA233344+(y1*y2**3*y&
         &5**2*sqrt(3.D0)/6.D0-y2**3*y3*y5**2*sqrt(3.D0)/3.D0-y2*y3**3*y5**2&
         &*sqrt(3.D0)/3.D0+y1**3*y2*y4*y5-y1**3*y2*y5**2*sqrt(3.D0)/3.D0-y1*&
         &*3*y3*y4*y5-y1**3*y3*y5**2*sqrt(3.D0)/3.D0-y1*y3**3*y4**2*sqrt(3.D&
         &0)/2.D0+y1*y3**3*y5**2*sqrt(3.D0)/6.D0-y2**3*y3*y4*y5+y2*y3**3*y4*&
         &y5-y1*y2**3*y4**2*sqrt(3.D0)/2.D0)*FEA233345+(-3.D0*y2**3*y4*y5**2&
         &+y3**3*y4**3-3.D0*y3**3*y4*y5**2-3.D0*y1**3*y4*y5**2+y2**3*y4**3+y&
         &1**3*y4**3)*FEA111444+(y1*y2**3*y3**2+y1**3*y2**2*y3+y1**2*y2**3*y&
         &3+y1*y2**2*y3**3+y1**2*y2*y3**3+y1**3*y2*y3**2)*FEA111233 
      s4 = s3+(9.D0*y4**2*y5**4-6.D0*y4**4*y5**2+y4**6)*FEA444444+(-5.D0&
         &/3.D0*y1*y2**2*y4**2*y5*sqrt(3.D0)+y1*y2**2*y4**3-4.D0/3.D0*y1**2*&
         &y3*y4**2*y5*sqrt(3.D0)-2.D0*y1**2*y2*y4**3-y1*y2**2*y5**3*sqrt(3.D&
         &0)/3.D0+4.D0/3.D0*y2**2*y3*y4*y5**2-4.D0/3.D0*y2**2*y3*y5**3*sqrt(&
         &3.D0)-2.D0*y1**2*y3*y4**3+7.D0/3.D0*y1*y2**2*y4*y5**2-2.D0/3.D0*y1&
         &**2*y3*y4*y5**2+y1*y3**2*y4**3+4.D0/3.D0*y2*y3**2*y5**3*sqrt(3.D0)&
         &+y1*y3**2*y5**3*sqrt(3.D0)/3.D0+4.D0/3.D0*y1**2*y2*y4**2*y5*sqrt(3&
         &.D0)+4.D0/3.D0*y2*y3**2*y4*y5**2+5.D0/3.D0*y1*y3**2*y4**2*y5*sqrt(&
         &3.D0)-2.D0/3.D0*y1**2*y2*y4*y5**2+7.D0/3.D0*y1*y3**2*y4*y5**2)*FEA&
         &133444 
      s5 = s4+(-y1**3*y2*y4*y5+2.D0/3.D0*y2**3*y3*y5**2*sqrt(3.D0)+y1*y3&
         &**3*y4**2*sqrt(3.D0)/2.D0+y1**3*y3*y4**2*sqrt(3.D0)/2.D0+y1**3*y3*&
         &y5**2*sqrt(3.D0)/6.D0+y1**3*y2*y5**2*sqrt(3.D0)/6.D0+y1**3*y3*y4*y&
         &5+y1*y2**3*y5**2*sqrt(3.D0)/6.D0+y1**3*y2*y4**2*sqrt(3.D0)/2.D0+2.&
         &D0/3.D0*y2*y3**3*y5**2*sqrt(3.D0)-y1*y2**3*y4*y5+y1*y2**3*y4**2*sq&
         &rt(3.D0)/2.D0+y1*y3**3*y5**2*sqrt(3.D0)/6.D0+y1*y3**3*y4*y5)*FEA13&
         &3345 
      v6 = s5+(-y2**2*y3*y4**2*y5+y1**2*y3*y4*y5**2*sqrt(3.D0)/3.D0+y2*y&
         &3**2*y4**2*y5+y2*y3**2*y5**3-y1*y2**2*y5**3+4.D0/3.D0*y2**2*y3*y4*&
         &y5**2*sqrt(3.D0)+4.D0/3.D0*y2*y3**2*y4*y5**2*sqrt(3.D0)-y1*y2**2*y&
         &4**2*y5+4.D0/3.D0*y1*y3**2*y4*y5**2*sqrt(3.D0)-y2**2*y3*y5**3+y1*y&
         &3**2*y5**3+y1**2*y2*y4*y5**2*sqrt(3.D0)/3.D0-y1**2*y2*y4**3*sqrt(3&
         &.D0)+y1*y3**2*y4**2*y5-y1**2*y3*y4**3*sqrt(3.D0)+4.D0/3.D0*y1*y2**&
         &2*y4*y5**2*sqrt(3.D0))*FEA233445+(y2*y3**4*y4*sqrt(3.D0)-y1**4*y2*&
         &y5+y2**4*y3*y4*sqrt(3.D0)-y1**4*y3*y4*sqrt(3.D0)+y2*y3**4*y5-2.D0*&
         &y1*y2**4*y5+2.D0*y1*y3**4*y5-y1**4*y2*y4*sqrt(3.D0)+y1**4*y3*y5-y2&
         &**4*y3*y5)*FEA233335+(y2**2*y3**4+y1**4*y3**2+y1**2*y2**4+y2**4*y3&
         &**2+y1**2*y3**4+y1**4*y2**2)*FEA222233
 
      v = v0+v1+v2+v3+v4+v5+v6 ! + ve -- add it in test if compared to ab initio energy 
 
      if (flag > 0) then
         ! dfea
         dfea1= f1a1+f2a1*2d0*coro+f3a1*3d0*coro**2+f4a1*4d0*coro**3+f5a1*5d0*coro**4+f6a1*6d0*coro**5
 
         dfea11=   f1a11+f2a11*2d0*coro+f3a11*3d0*coro**2+f4a11*4d0*coro**3
         dfea12=   f1a12+f2a12*2d0*coro+f3a12*3d0*coro**2+f4a12*4d0*coro**3
         dfea14=   f1a14+f2a14*2d0*coro+f3a14*3d0*coro**2+f4a14*4d0*coro**3
         dfea44=   f1a44+f2a44*2d0*coro+f3a44*3d0*coro**2+f4a44*4d0*coro**3
 
         dfea111= f1a111+f2a111*2d0*coro+f3a111*3d0*coro**2
         dfea112= f1a112+f2a112*2d0*coro+f3a112*3d0*coro**2
         dfea114= f1a114+f2a114*2d0*coro+f3a114*3d0*coro**2
         dfea123= f1a123+f2a123*2d0*coro+f3a123*3d0*coro**2
         dfea124= f1a124+f2a124*2d0*coro+f3a124*3d0*coro**2
         dfea144= f1a144+f2a144*2d0*coro+f3a144*3d0*coro**2
         dfea155= f1a155+f2a155*2d0*coro+f3a155*3d0*coro**2
         dfea455= f1a455+f2a455*2d0*coro+f3a455*3d0*coro**2
 
         dfea1111= f1a1111+f2a1111*2d0*coro
         dfea1112= f1a1112+f2a1112*2d0*coro
         dfea1114= f1a1114+f2a1114*2d0*coro
         dfea1122= f1a1122+f2a1122*2d0*coro
         dfea1123= f1a1123+f2a1123*2d0*coro
         dfea1124= f1a1124+f2a1124*2d0*coro
         dfea1125= f1a1125+f2a1125*2d0*coro
         dfea1144= f1a1144+f2a1144*2d0*coro
         dfea1155= f1a1155+f2a1155*2d0*coro
         dfea1244= f1a1244+f2a1244*2d0*coro
         dfea1255= f1a1255+f2a1255*2d0*coro
         dfea1444= f1a1444+f2a1444*2d0*coro
         dfea1455= f1a1455+f2a1455*2d0*coro
         dfea4444= f1a4444+f2a4444*2d0*coro
 
         dfea44444 = f1a44444 + f2a44444 *2d0*coro
         dfea33455 = f1a33455 + f2a33455 *2d0*coro
         dfea33445 = f1a33445 + f2a33445 *2d0*coro
         dfea33345 = f1a33345 + f2a33345 *2d0*coro
         dfea33344 = f1a33344 + f2a33344 *2d0*coro
         dfea33334 = f1a33334 + f2a33334 *2d0*coro
         dfea33333 = f1a33333 + f2a33333 *2d0*coro
         dfea25555 = f1a25555 + f2a25555 *2d0*coro
         dfea24455 = f1a24455 + f2a24455 *2d0*coro
         dfea24445 = f1a24445 + f2a24445 *2d0*coro
         dfea23333 = f1a23333 + f2a23333 *2d0*coro
         dfea13455 = f1a13455 + f2a13455 *2d0*coro
         dfea13445 = f1a13445 + f2a13445 *2d0*coro
         dfea13345 = f1a13345 + f2a13345 *2d0*coro
         dfea12355 = f1a12355 + f2a12355 *2d0*coro
         dfea11334 = f1a11334 + f2a11334 *2d0*coro
         dfea11333 = f1a11333 + f2a11333 *2d0*coro
         dfea11255 = f1a11255 + f2a11255 *2d0*coro
         dfea11245 = f1a11245 + f2a11245 *2d0*coro
         dfea11234 = f1a11234 + f2a11234 *2d0*coro
         dfea11233 = f1a11233 + f2a11233 *2d0*coro
         dfea11135 = f1a11135 + f2a11135 *2d0*coro
         dfea11134 = f1a11134 + f2a11134 *2d0*coro
         dfea11123 = f1a11123 + f2a11123 *2d0*coro
         dfea555555= f1a555555+ f2a555555*2d0*coro
         dfea444444= f1a444444+ f2a444444*2d0*coro
         dfea335555= f1a335555+ f2a335555*2d0*coro
         dfea334455= f1a334455+ f2a334455*2d0*coro
         dfea334445= f1a334445+ f2a334445*2d0*coro
         dfea333555= f1a333555+ f2a333555*2d0*coro
         dfea333333= f1a333333+ f2a333333*2d0*coro
         dfea244555= f1a244555+ f2a244555*2d0*coro
         dfea244455= f1a244455+ f2a244455*2d0*coro
         dfea233445= f1a233445+ f2a233445*2d0*coro
         dfea233444= f1a233444+ f2a233444*2d0*coro
         dfea233345= f1a233345+ f2a233345*2d0*coro
         dfea233344= f1a233344+ f2a233344*2d0*coro
         dfea233335= f1a233335+ f2a233335*2d0*coro
         dfea223355= f1a223355+ f2a223355*2d0*coro
         dfea222335= f1a222335+ f2a222335*2d0*coro
         dfea222334= f1a222334+ f2a222334*2d0*coro
         dfea222333= f1a222333+ f2a222333*2d0*coro
         dfea222255= f1a222255+ f2a222255*2d0*coro
         dfea222245= f1a222245+ f2a222245*2d0*coro
         dfea222233= f1a222233+ f2a222233*2d0*coro
         dfea222224= f1a222224+ f2a222224*2d0*coro
         dfea145555= f1a145555+ f2a145555*2d0*coro
         dfea134444= f1a134444+ f2a134444*2d0*coro
         dfea133444= f1a133444+ f2a133444*2d0*coro
         dfea133345= f1a133345+ f2a133345*2d0*coro
         dfea133334= f1a133334+ f2a133334*2d0*coro
         dfea133333= f1a133333+ f2a133333*2d0*coro
         dfea124555= f1a124555+ f2a124555*2d0*coro
         dfea124455= f1a124455+ f2a124455*2d0*coro
         dfea123455= f1a123455+ f2a123455*2d0*coro
         dfea123345= f1a123345+ f2a123345*2d0*coro
         dfea113555= f1a113555+ f2a113555*2d0*coro
         dfea113345= f1a113345+ f2a113345*2d0*coro
         dfea112355= f1a112355+ f2a112355*2d0*coro
         dfea112335= f1a112335+ f2a112335*2d0*coro
         dfea112233= f1a112233+ f2a112233*2d0*coro
         dfea111444= f1a111444+ f2a111444*2d0*coro
         dfea111234= f1a111234+ f2a111234*2d0*coro
         dfea111233= f1a111233+ f2a111233*2d0*coro
         dfea111123= f1a111123+ f2a111123*2d0*coro
 
         ! calculate 1st derivatives of potential energy (gradients)
 
         dv0=-f1a+f2a*2d0*coro+f3a*3d0*coro**2+f4a*4d0*coro**3+f5a*5d0*coro**4 &
            +f6a*6d0*coro**5+f7a*7d0*coro**6+f8a*8d0*coro**7
 
         dv1(1)=( +1d0 )*FEA1 
 
         dv2(1)=( +y3 +y2 )*FEA12 +( +2d0*y1 )*FEA11 +( &
            +y4 )*FEA14 +0d0 
 
         dv3(1)=( y3*y4 +y2*y4 &
            +SQRT3*y2*y5 -SQRT3*y3*y5 )*FEA124 +( +y5**2 )*FEA155 +( +y3**2 &
            +2d0*y1*y3 +y2**2 +2d0*y1*y2 )*FEA112 +0d0 +y2*y3*FEA123 &
            +( y4**2 )*FEA144 +( +3d0*y1**2 )*FEA111 +( &
            +2d0*y1*y4 )*FEA114 
 
         ds2(1)=0d0 +( +y4**3 )*FEA1444 +( &
            +2d0*y1*y5**2 )*FEA1155 
 
         ds1(1)=ds2(1) +( +2d0*y1*y4**2 )*FEA1144 &
            +( 3d0*y1**2*y4 )*FEA1114 +( +4d0*y1**3 )*FEA1111 +( SQRT3*y3*y4*y5 &
            +y2*y4**2 -SQRT3*y2*y4*y5 +y3*y4**2 )*FEA1244 
 
         ds2(1)=ds1(1) &
            +( y3*y5**2 +y2*y5**2 -SQRT3*y3*y4*y5 +SQRT3*y2*y4*y5 )*FEA1255 &
            +( -y3**2*y4/2.D0 +2d0*y1*y3*y4 -SQRT3*y3**2*y5/2.D0 +2d0*y1*y2*y4 &
            +SQRT3*y2**2*y5/2.D0 -y2**2*y4/2.D0 )*FEA1124 +( 2d0*y1*y2*y5 &
            +SQRT3*y3**2*y4/2.D0 +SQRT3*y2**2*y4/2.D0 -y3**2*y5/2.D0 &
            +y2**2*y5/2.D0 -2d0*y1*y3*y5 )*FEA1125 
 
         dv4(1)=ds2(1) +( +3d0*y1**2*y3 &
            +3d0*y1**2*y2 +y2**3 +y3**3 )*FEA1112 +( +2d0*y1*y3**2 &
            +2d0*y1*y2**2 )*FEA1122 +( y2**2*y3 +2d0*y1*y2*y3 +y2*y3**2 )*FEA1123 &
            +( +y4*y5**2 )*FEA1455 
 
         ds3(1)=0d0 +( &
            +9.D0*y4**2*y5**2 -3.D0/2.D0*y4**4 +5.D0/2.D0*y5**4 )*FEA25555 &
            +( -7.D0/2.D0*y4**2*y5**2 -3.D0/4.D0*y5**4 &
            +5.D0/4.D0*y4**4 )*FEA24455 
 
         ds2(1)=ds3(1) +( &
            +3.D0/4.D0*y5**4*SQRT3 -7.D0/12.D0*y4**4*SQRT3 &
            +3.D0/2.D0*y4**2*y5**2*SQRT3 )*FEA24445 &
            +( -5.D0/9.D0*2d0*y1*y4**3*SQRT3 -2d0*y1*y4*y5**2*SQRT3 )*FEA33445 &
            +( -2d0*y1*y4**3/3.D0 +2d0*y1*y4*y5**2 )*FEA33455 
 
         ds1(1)=ds2(1) +( &
            +3d0*y1**2*y4**2*SQRT3/2.D0 -3d0*y1**2*y5**2*SQRT3/6.D0 )*FEA33345 +( &
            +3d0*y1**2*y5**2 +3d0*y1**2*y4**2 )*FEA33344 &
            +( -2.D0*4d0*y1**3*y4 )*FEA33334 +( +5d0*y1**4 )*FEA33333 &
            +( -4.D0/9.D0*y2*y4**3*SQRT3 -y2*y5**3 &
            +y3*y4**2*y5 -y2*y4**2*y5 -4.D0/9.D0*y3*y4**3*SQRT3 &
            +y3*y5**3 )*FEA13445 +( +y2*y4*y5**2 -y2*y4**3/3.D0 -y3*y4**3/3.D0 &
            +y3*y4*y5**2 )*FEA13455 
 
         ds3(1)=ds1(1) +( 2d0*y1*y3*y5**2 +y2**2*y5**2 &
            +2d0*y1*y2*y5**2 +y2**2*y4**2 +y3**2*y4**2 +2d0*y1*y3*y4**2 &
            +2d0*y1*y2*y4**2 +y3**2*y5**2 )*FEA11255 &
            +( 2.D0/3.D0*2d0*y1*y3*y4**2*SQRT3 +y3**2*y5**2*SQRT3/2.D0 &
            +y2**2*y5**2*SQRT3/2.D0 -y2**2*y4*y5 +y3**2*y4*y5 &
            +y3**2*y4**2*SQRT3/6.D0 +y2**2*y4**2*SQRT3/6.D0 &
            +2.D0/3.D0*2d0*y1*y2*y4**2*SQRT3 )*FEA13345 
 
         ds4(1)=ds3(1) &
            +( 2d0*y1*y2*y4*y5 +2d0*y1*y3*y4**2*SQRT3/3.D0 &
            +2d0*y1*y2*y4**2*SQRT3/3.D0 -y2**2*y4**2*SQRT3/6.D0 -2d0*y1*y3*y4*y5 &
            +y2**2*y5**2*SQRT3/2.D0 -y3**2*y4**2*SQRT3/6.D0 &
            +y3**2*y5**2*SQRT3/2.D0 )*FEA11245 
 
         ds2(1)=ds4(1) +( -3d0*y1**2*y2*y5 &
            +3d0*y1**2*y3*y5 -y2**3*y4*SQRT3/2.D0 -y2**3*y5/2.D0 &
            +y3**3*y5/2.D0 -y3**3*y4*SQRT3/2.D0 )*FEA11135 +( 3d0*y1**2*y3*y4 &
            +3d0*y1**2*y2*y4 -y3**3*y4/2.D0 &
            +y2**3*y5*SQRT3/2.D0 -y2**3*y4/2.D0 &
            -y3**3*y5*SQRT3/2.D0 )*FEA11134 
 
         dv5(1)=ds2(1) &
            +( y2**4 +4d0*y1**3*y3 +4d0*y1**3*y2 +y3**4 )*FEA23333 +( &
            +2d0*y1*y2**2*y4 -SQRT3*2d0*y1*y3**2*y5 +SQRT3*2d0*y1*y2**2*y5 &
            +2d0*y1*y3**2*y4 )*FEA11334 +( 2d0*y1*y3**3 +3d0*y1**2*y3**2 &
            +2d0*y1*y2**3 +3d0*y1**2*y2**2 )*FEA11333 +( y2*y3*y4**2 &
            +y2*y3*y5**2 )*FEA12355 &
            +( -y2*y3**2*y4/2.D0 -y2**2*y3*y4/2.D0 -SQRT3*y2*y3**2*y5/2.D0 &
            +2d0*y1*y2*y3*y4 +SQRT3*y2**2*y3*y5/2.D0 )*FEA11234 +( y2**3*y3 &
            +y2*y3**3 +3d0*y1**2*y2*y3 )*FEA11123 +( 2d0*y1*y2**2*y3 +y2**2*y3**2 &
            +2d0*y1*y2*y3**2 )*FEA11233 
 
         ds3(1)=( &
            -8.D0/3.D0*3d0*y1**2*y4*y5**2*SQRT3 )*FEA333555 &
            +( 4d0*y1**3*y5**2*SQRT3/2.D0 -4d0*y1**3*y4**2*SQRT3/6.D0 )*FEA222245 &
            +( y3**5 +y2**5 +5d0*y1**4*y3 +5d0*y1**4*y2 )*FEA133333 &
            +( 4d0*y1**3*y3*y4 +4d0*y1**3*y2*y4 +y2**4*y5*SQRT3 +y3**4*y4 &
            +4d0*y1**3*y2*y5*SQRT3 -y3**4*y5*SQRT3 -4d0*y1**3*y3*y5*SQRT3 &
            +y2**4*y4 )*FEA133334 +( -y2*y3*y4**3/3.D0 &
            +y2*y3*y4*y5**2 )*FEA123455 
 
         ds4(1)=ds3(1) &
            +( 2.D0/3.D0*SQRT3*y2**2*y3**2*y4 -2d0*y1*y2**2*y3*y5 &
            -SQRT3*2d0*y1*y2**2*y3*y4/3.D0 &
            +2d0*y1*y2*y3**2*y5 -SQRT3*2d0*y1*y2*y3**2*y4/3.D0 )*FEA112335 &
            +( y2**2*y3*y5**2 +y2*y3**2*y5**2 +y2*y3**2*y4**2 +y2**2*y3*y4**2 &
            +2d0*y1*y2*y3*y4**2 +2d0*y1*y2*y3*y5**2 )*FEA112355 
 
         ds2(1)=ds4(1) &
            +( -3d0*y1**2*y2**2*y5/2.D0 -2d0*y1*y3**3*y5/2.D0 &
            +3d0*y1**2*y2**2*y4*SQRT3/2.D0 -2d0*y1*y2**3*y4*SQRT3/2.D0 &
            +3d0*y1**2*y3**2*y5/2.D0 +2d0*y1*y2**3*y5/2.D0 &
            +3d0*y1**2*y3**2*y4*SQRT3/2.D0 -2d0*y1*y3**3*y4*SQRT3/2.D0 )*FEA222335 &
            +( -2d0*y1*y2**2*y5**2*SQRT3/2.D0 -2d0*y1*y3**2*y5**2*SQRT3/2.D0 &
            -2d0*y1*y2**2*y4**2*SQRT3/6.D0 -2d0*y1*y2**2*y4*y5 &
            +2d0*y1*y3**2*y4*y5 -2d0*y1*y3**2*y4**2*SQRT3/6.D0 )*FEA113345 +( &
            +2d0*y1*y2**2*y5**2 +2d0*y1*y3**2*y4**2 +2d0*y1*y3**2*y5**2 &
            +2d0*y1*y2**2*y4**2 )*FEA223355 
 
         ds3(1)=ds2(1) &
            +( y2*y3**2*y4**2*SQRT3/6.D0 +y2*y3**2*y4*y5 &
            +y2*y3**2*y5**2*SQRT3/2.D0 &
            +2.D0/3.D0*2d0*y1*y2*y3*y4**2*SQRT3 -y2**2*y3*y4*y5 &
            +y2**2*y3*y4**2*SQRT3/6.D0 +y2**2*y3*y5**2*SQRT3/2.D0 )*FEA123345 &
            +( -3d0*y1**2*y2**2*y5*SQRT3/2.D0 -3d0*y1**2*y2**2*y4/2.D0 &
            -3d0*y1**2*y3**2*y4/2.D0 -2d0*y1*y2**3*y4/2.D0 &
            +3d0*y1**2*y3**2*y5*SQRT3/2.D0 -2d0*y1*y3**3*y4/2.D0 &
            -2d0*y1*y2**3*y5*SQRT3/2.D0 &
            +2d0*y1*y3**3*y5*SQRT3/2.D0 )*FEA222334 +( +5.D0/2.D0*2d0*y1*y5**4 &
            +9.D0*2d0*y1*y4**2*y5**2 -3.D0/2.D0*2d0*y1*y4**4 )*FEA335555 &
            +( 3d0*y1**2*y2**3 +3d0*y1**2*y3**3 )*FEA222333 
 
         ds4(1)=ds3(1) &
            +( -2.D0/5.D0*y4**5 -2.D0*y4**3*y5**2 )*FEA244455 &
            +( -2.D0*5d0*y1**4*y4 )*FEA222224 
 
         ds5(1)=ds4(1) +( +y4*y5**4 &
            +2.D0*y4**3*y5**2 -y4**5/15.D0 )*FEA145555 
 
         ds1(1)=ds5(1) &
            +( -SQRT3*y2*y3**3*y5/2.D0 +3d0*y1**2*y2*y3*y4 &
            +SQRT3*y2**3*y3*y5/2.D0 -y2**3*y3*y4/2.D0 -y2*y3**3*y4/2.D0 )*FEA111234 &
            +( +2.D0/9.D0*y4**5*SQRT3 -2.D0/3.D0*y4**3*y5**2*SQRT3 )*FEA244555 &
            +( y2*y4**2*y5**2 -y2*y5**4 -y3*y5**4 &
            +y3*y4**2*y5**2 -2.D0*y2*y4**3*y5*SQRT3 &
            +2.D0*y3*y4**3*y5*SQRT3 )*FEA124455 
 
         ds3(1)=ds1(1) +( &
            +6d0*y1**5 )*FEA333333 +( y2**4*y3 +4d0*y1**3*y2*y3 &
            +y2*y3**4 )*FEA111123 +2d0*y1*y2**2*y3**2*FEA112233 +( 4d0*y1**3*y4**2 &
            +4d0*y1**3*y5**2 )*FEA222255 
 
         ds4(1)=ds3(1) +( 3.D0*y3*y5**4 &
            +y3*y4**4 -4.D0*y3*y4**3*y5*SQRT3 +y2*y4**4 +4.D0*y2*y4**3*y5*SQRT3 &
            +3.D0*y2*y5**4 )*FEA134444 &
            +( -y3**2*y5**3*SQRT3/3.D0 -7.D0/3.D0*2d0*y1*y3*y4*y5**2 &
            +5.D0/3.D0*y2**2*y4**2*y5*SQRT3 -7.D0/3.D0*2d0*y1*y2*y4*y5**2 &
            -16.D0/3.D0*y3**2*y4*y5**2 &
            +4.D0/3.D0*2d0*y1*y3*y4**2*y5*SQRT3 +3.D0*2d0*y1*y2*y4**3 &
            +y2**2*y5**3*SQRT3/3.D0 -5.D0/3.D0*y3**2*y4**2*y5*SQRT3 &
            -4.D0/3.D0*2d0*y1*y2*y4**2*y5*SQRT3 &
            +3.D0*2d0*y1*y3*y4**3 &
            -16.D0/3.D0*y2**2*y4*y5**2 )*FEA233444 
 
         ds5(1)=ds4(1) &
            +( 2.D0*y3**2*y5**3 -2.D0*y2**2*y5**3 +2d0*y1*y3*y4*y5**2*SQRT3 &
            +6.D0*y3**2*y4**2*y5 -6.D0*y2**2*y4**2*y5 -3.D0*2d0*y1*y3*y4**2*y5 &
            +2d0*y1*y2*y4*y5**2*SQRT3 &
            +4.D0*y3**2*y4*y5**2*SQRT3 -3.D0*2d0*y1*y2*y4**3*SQRT3 &
            +3.D0*2d0*y1*y2*y4**2*y5 -2d0*y1*y2*y5**3 &
            +2d0*y1*y3*y5**3 -3.D0*2d0*y1*y3*y4**3*SQRT3 &
            +4.D0*y2**2*y4*y5**2*SQRT3 )*FEA113555 
 
         ds2(1)=ds5(1) &
            +( -3.D0/2.D0*2d0*y1*y4**2*y5**2*SQRT3 -3.D0/4.D0*2d0*y1*y5**4*SQRT3 &
            +7.D0/12.D0*2d0*y1*y4**4*SQRT3 )*FEA334445 +( -3.D0*y3*y4**3*y5 &
            +2.D0/3.D0*y2*y5**4*SQRT3 -y3*y4*y5**3 +2.D0/3.D0*y3*y5**4*SQRT3 &
            +3.D0*y2*y4**3*y5 +y2*y4*y5**3 )*FEA124555 &
            +( -7.D0/2.D0*2d0*y1*y4**2*y5**2 -3.D0/4.D0*2d0*y1*y5**4 &
            +5.D0/4.D0*2d0*y1*y4**4 )*FEA334455 
 
         ds3(1)=ds2(1) +0d0 +( +y3**3*y4**2 &
            +y2**3*y4**2 +3d0*y1**2*y2*y4**2 +y2**3*y5**2 +3d0*y1**2*y3*y5**2 &
            +3d0*y1**2*y3*y4**2 +3d0*y1**2*y2*y5**2 +y3**3*y5**2 )*FEA233344 &
            +( y2**3*y5**2*SQRT3/6.D0 &
            +3d0*y1**2*y2*y4*y5 -3d0*y1**2*y2*y5**2*SQRT3/3.D0 -3d0*y1**2*y3*y4*y5 &
            -3d0*y1**2*y3*y5**2*SQRT3/3.D0 -y3**3*y4**2*SQRT3/2.D0 &
            +y3**3*y5**2*SQRT3/6.D0 -y2**3*y4**2*SQRT3/2.D0 )*FEA233345 &
            +( -3.D0*3d0*y1**2*y4*y5**2 +3d0*y1**2*y4**3 )*FEA111444 &
            +( y2**3*y3**2 +3d0*y1**2*y2**2*y3 +2d0*y1*y2**3*y3 +y2**2*y3**3 &
            +2d0*y1*y2*y3**3 +3d0*y1**2*y2*y3**2 )*FEA111233 
 
         ds4(1)=ds3(1) +0d0 &
            +( -5.D0/3.D0*y2**2*y4**2*y5*SQRT3 &
            +y2**2*y4**3 -4.D0/3.D0*2d0*y1*y3*y4**2*y5*SQRT3 -2.D0*2d0*y1*y2*y4**3 &
            -y2**2*y5**3*SQRT3/3.D0 -2.D0*2d0*y1*y3*y4**3 &
            +7.D0/3.D0*y2**2*y4*y5**2 -2.D0/3.D0*2d0*y1*y3*y4*y5**2 +y3**2*y4**3 &
            +y3**2*y5**3*SQRT3/3.D0 +4.D0/3.D0*2d0*y1*y2*y4**2*y5*SQRT3 &
            +5.D0/3.D0*y3**2*y4**2*y5*SQRT3 -2.D0/3.D0*2d0*y1*y2*y4*y5**2 &
            +7.D0/3.D0*y3**2*y4*y5**2 )*FEA133444 
 
         ds5(1)=ds4(1) &
            +( -3d0*y1**2*y2*y4*y5 +y3**3*y4**2*SQRT3/2.D0 &
            +3d0*y1**2*y3*y4**2*SQRT3/2.D0 +3d0*y1**2*y3*y5**2*SQRT3/6.D0 &
            +3d0*y1**2*y2*y5**2*SQRT3/6.D0 +3d0*y1**2*y3*y4*y5 &
            +y2**3*y5**2*SQRT3/6.D0 +3d0*y1**2*y2*y4**2*SQRT3/2.D0 -y2**3*y4*y5 &
            +y2**3*y4**2*SQRT3/2.D0 +y3**3*y5**2*SQRT3/6.D0 &
            +y3**3*y4*y5 )*FEA133345 
 
         dv6(1)=ds5(1) +( &
            +2d0*y1*y3*y4*y5**2*SQRT3/3.D0 -y2**2*y5**3 -y2**2*y4**2*y5 &
            +4.D0/3.D0*y3**2*y4*y5**2*SQRT3 +y3**2*y5**3 &
            +2d0*y1*y2*y4*y5**2*SQRT3/3.D0 -2d0*y1*y2*y4**3*SQRT3 &
            +y3**2*y4**2*y5 -2d0*y1*y3*y4**3*SQRT3 &
            +4.D0/3.D0*y2**2*y4*y5**2*SQRT3 )*FEA233445 &
            +( -4d0*y1**3*y2*y5 -4d0*y1**3*y3*y4*SQRT3 -2.D0*y2**4*y5 &
            +2.D0*y3**4*y5 -4d0*y1**3*y2*y4*SQRT3 +4d0*y1**3*y3*y5 )*FEA233335 +( &
            +4d0*y1**3*y3**2 +2d0*y1*y2**4 +2d0*y1*y3**4 &
            +4d0*y1**3*y2**2 )*FEA222233 
 
         dv(1)= +dv1(1) +dv2(1) +dv3(1) &
            +dv4(1) +dv5(1) +dv6(1)
 

         dv1(2)=( +1d0 )*FEA1 
 
         dv2(2)=( y3 +y1 )*FEA12 +( 2d0*y2 )*FEA11 +( &
            +SQRT3*y5/2.D0 -y4/2.D0 )*FEA14 +0d0 
 
         dv3(2)=( +y1*y4 -2.D0*y3*y4 &
            +SQRT3*y1*y5 )*FEA124 +( +y5**2/4.D0 +3.D0/4.D0*y4**2 &
            +SQRT3*y4*y5/2.D0 )*FEA155 +( y3**2 +y1*2d0*y2 +2d0*y2*y3 &
            +y1**2 )*FEA112 +0d0 +y1*y3*FEA123 +( +3.D0/4.D0*y5**2 &
            +y4**2/4.D0 -SQRT3*y4*y5/2.D0 )*FEA144 +( +3d0*y2**2 )*FEA111 &
            +( -2d0*y2*y4/2.D0 +SQRT3*2d0*y2*y5/2.D0 )*FEA114 
 
         ds2(2)=0d0 &
            +( 3.D0/8.D0*SQRT3*y5**3 -9.D0/8.D0*y4*y5**2 -y4**3/8.D0 &
            +3.D0/8.D0*SQRT3*y4**2*y5 )*FEA1444 +( 3.D0/4.D0*2d0*y2*y4**2 &
            +SQRT3*2d0*y2*y4*y5/2.D0 +2d0*y2*y5**2/4.D0 )*FEA1155 
 
         ds1(2)=ds2(2) +( &
            +2d0*y2*y4**2/4.D0 -SQRT3*2d0*y2*y4*y5/2.D0 &
            +3.D0/4.D0*2d0*y2*y5**2 )*FEA1144 +( &
            +SQRT3*3d0*y2**2*y5/2.D0 -3d0*y2**2*y4/2.D0 )*FEA1114 &
            +( 4d0*y2**3 )*FEA1111 +( +3.D0/2.D0*y3*y5**2 -y3*y4**2/2.D0 &
            +y1*y4**2 -SQRT3*y1*y4*y5 )*FEA1244 
 
         ds2(2)=ds1(2) +( &
            +y1*y5**2 -y3*y5**2/2.D0 +3.D0/2.D0*y3*y4**2 +SQRT3*y1*y4*y5 )*FEA1255 &
            +( -SQRT3*y3**2*y5/2.D0 +y1**2*y4 &
            +SQRT3*2d0*y2*y3*y5/2.D0 -2d0*y2*y3*y4/2.D0 &
            +SQRT3*y1*2d0*y2*y5/2.D0 -y3**2*y4/2.D0 -y1*2d0*y2*y4/2.D0 )*FEA1124 &
            +( y1**2*y5 &
            +SQRT3*y1*2d0*y2*y4/2.D0 -SQRT3*y3**2*y4/2.D0 -SQRT3*2d0*y2*y3*y4/2.D0 &
            -2d0*y2*y3*y5/2.D0 &
            +y3**2*y5/2.D0 +y1*2d0*y2*y5/2.D0 )*FEA1125 
 
         dv4(2)=ds2(2) +( y3**3 &
            +y1**3 +y1*3d0*y2**2 +3d0*y2**2*y3 )*FEA1112 +( 2d0*y2*y3**2 &
            +y1**2*2d0*y2 )*FEA1122 +( y1*2d0*y2*y3 +y1**2*y3 +y1*y3**2 )*FEA1123 &
            +( 5.D0/8.D0*y4*y5**2 +SQRT3*y5**3/8.D0 &
            +SQRT3*y4**2*y5/8.D0 -3.D0/8.D0*y4**3 )*FEA1455 
 
         ds3(2)=0d0 +( &
            +4.D0*y4*y5**3*SQRT3 +3.D0*y4**4 +y5**4 )*FEA25555 &
            +( -y4**4 -2.D0*y4*y5**3*SQRT3 +y4**2*y5**2 )*FEA24455 
 
         ds2(2)=ds3(2) &
            +( y4**3*y5 +3.D0*y4*y5**3 +2.D0/3.D0*y4**4*SQRT3 )*FEA24445 &
            +( -2d0*y2*y5**3 &
            +4.D0/9.D0*2d0*y2*y4**3*SQRT3 -2d0*y2*y4**2*y5 )*FEA33445 +( &
            +2d0*y2*y4*y5**2 -2d0*y2*y4**3/3.D0 )*FEA33455 
 
         ds1(2)=ds2(2) &
            +( -3d0*y2**2*y4*y5 +3d0*y2**2*y5**2*SQRT3/3.D0 )*FEA33345 +( &
            +3d0*y2**2*y4**2 +3d0*y2**2*y5**2 )*FEA33344 +( &
            +4d0*y2**3*y4 -SQRT3*4d0*y2**3*y5 )*FEA33334 +( 5d0*y2**4 )*FEA33333 &
            +( -4.D0/9.D0*y1*y4**3*SQRT3 -y1*y5**3 +y3*y4*y5**2*SQRT3 -y1*y4**2*y5 &
            +5.D0/9.D0*y3*y4**3*SQRT3 )*FEA13445 +( y3*y4*y5**2 &
            +y1*y4*y5**2 -y3*y4**3/3.D0 -y1*y4**3/3.D0 )*FEA13455 
 
         ds3(2)=ds1(2) +( &
            +2d0*y2*y3*y4**2 +2d0*y2*y3*y5**2 +y1*2d0*y2*y5**2 +y1**2*y5**2 &
            +y1*2d0*y2*y4**2 +y3**2*y4**2 +y1**2*y4**2 +y3**2*y5**2 )*FEA11255 +( &
            +y1*2d0*y2*y5**2*SQRT3/2.D0 &
            +2d0*y2*y3*y5**2*SQRT3/2.D0 -y1*2d0*y2*y4*y5 &
            +y3**2*y4*y5 -2d0*y2*y3*y4*y5 +y3**2*y4**2*SQRT3/6.D0 &
            +y1*2d0*y2*y4**2*SQRT3/6.D0 +2.D0/3.D0*y1**2*y4**2*SQRT3 &
            +y3**2*y5**2*SQRT3/2.D0 &
            +2d0*y2*y3*y4**2*SQRT3/6.D0 )*FEA13345 
 
         ds4(2)=ds3(2) +( y1**2*y4*y5 &
            +y1**2*y4**2*SQRT3/3.D0 -y1*2d0*y2*y4**2*SQRT3/6.D0 &
            +y3**2*y4*y5 -2d0*y2*y3*y4*y5 +y3**2*y4**2*SQRT3/3.D0 &
            +y1*2d0*y2*y5**2*SQRT3/2.D0 &
            +2d0*y2*y3*y4**2*SQRT3/3.D0 )*FEA11245 
 
         ds2(2)=ds4(2) +( -y1**3*y5 &
            +3d0*y2**2*y3*y5/2.D0 -y1*3d0*y2**2*y4*SQRT3/2.D0 &
            -y1*3d0*y2**2*y5/2.D0 -y3**3*y5/2.D0 &
            +3d0*y2**2*y3*y4*SQRT3/2.D0 +y3**3*y4*SQRT3/2.D0 )*FEA11135 &
            +( -3d0*y2**2*y3*y4/2.D0 +y1**3*y4 -y3**3*y4/2.D0 &
            +y1*3d0*y2**2*y5*SQRT3/2.D0 &
            +3d0*y2**2*y3*y5*SQRT3/2.D0 -y3**3*y5*SQRT3/2.D0 &
            -y1*3d0*y2**2*y4/2.D0 )*FEA11134 
 
         dv5(2)=ds2(2) &
            +( y1*4d0*y2**3 +y1**4 +4d0*y2**3*y3 +y3**4 )*FEA23333 &
            +( -2.D0*2d0*y2*y3**2*y4 +y1**2*2d0*y2*y4 &
            +SQRT3*y1**2*2d0*y2*y5 )*FEA11334 +( +2d0*y2*y3**3 +y1**2*3d0*y2**2 &
            +y1**3*2d0*y2 +3d0*y2**2*y3**2 )*FEA11333 +( y1*y3*y4**2 &
            +y1*y3*y5**2 )*FEA12355 &
            +( -y1*y3**2*y4/2.D0 -y1*2d0*y2*y3*y4/2.D0 -SQRT3*y1*y3**2*y5/2.D0 &
            +y1**2*y3*y4 +SQRT3*y1*2d0*y2*y3*y5/2.D0 )*FEA11234 +( y1*3d0*y2**2*y3 &
            +y1*y3**3 +y1**3*y3 )*FEA11123 +( y1**2*2d0*y2*y3 +y1*2d0*y2*y3**2 &
            +y1**2*y3**2 )*FEA11233 
 
         ds3(2)=( 3d0*y2**2*y4**3*SQRT3 &
            -3d0*y2**2*y4**2*y5 -5.D0/3.D0*3d0*y2**2*y4*y5**2*SQRT3 &
            -3d0*y2**2*y5**3 )*FEA333555 &
            +( +4d0*y2**3*y4*y5 +4d0*y2**3*y4**2*SQRT3/3.D0 )*FEA222245 +( &
            +y1*5d0*y2**4 +5d0*y2**4*y3 +y1**5 +y3**5 )*FEA133333 &
            +( -2.D0*4d0*y2**3*y3*y4 +y1**4*y4 &
            +y1*4d0*y2**3*y5*SQRT3 -2.D0*y3**4*y4 +y1**4*y5*SQRT3 &
            +y1*4d0*y2**3*y4 )*FEA133334 +( -y1*y3*y4**3/3.D0 &
            +y1*y3*y4*y5**2 )*FEA123455 
 
         ds4(2)=ds3(2) &
            +( 2.D0/3.D0*SQRT3*y1*2d0*y2*y3**2*y4 -y1**2*2d0*y2*y3*y5 &
            -SQRT3*y1**2*2d0*y2*y3*y4/3.D0 &
            +y1**2*y3**2*y5 -SQRT3*y1**2*y3**2*y4/3.D0 )*FEA112335 &
            +( y1*2d0*y2*y3*y5**2 +y1*y3**2*y5**2 +y1*y3**2*y4**2 &
            +y1*2d0*y2*y3*y4**2 +y1**2*y3*y4**2 &
            +y1**2*y3*y5**2 )*FEA112355 
 
         ds2(2)=ds4(2) &
            +( 3d0*y2**2*y3**2*y5 -y1**3*2d0*y2*y5/2.D0 -2d0*y2*y3**3*y5 &
            +y1**3*2d0*y2*y4*SQRT3/2.D0 -y1**2*3d0*y2**2*y4*SQRT3/2.D0 &
            +y1**2*3d0*y2**2*y5/2.D0 )*FEA222335 &
            +( -y1**2*2d0*y2*y5**2*SQRT3/2.D0 -y1**2*2d0*y2*y4**2*SQRT3/6.D0 &
            -y1**2*2d0*y2*y4*y5 -2.D0/3.D0*2d0*y2*y3**2*y4**2*SQRT3 )*FEA113345 &
            +( 2d0*y2*y3**2*y5**2 +2d0*y2*y3**2*y4**2 +y1**2*2d0*y2*y5**2 &
            +y1**2*2d0*y2*y4**2 )*FEA223355 
 
         ds3(2)=ds2(2) &
            +( y1*y3**2*y4**2*SQRT3/6.D0 +y1*y3**2*y4*y5 &
            +y1*y3**2*y5**2*SQRT3/2.D0 &
            +2.D0/3.D0*y1**2*y3*y4**2*SQRT3 -y1*2d0*y2*y3*y4*y5 &
            +y1*2d0*y2*y3*y4**2*SQRT3/6.D0 &
            +y1*2d0*y2*y3*y5**2*SQRT3/2.D0 )*FEA123345 &
            +( -y1**3*2d0*y2*y5*SQRT3/2.D0 -y1**3*2d0*y2*y4/2.D0 &
            -y1**2*3d0*y2**2*y4/2.D0 &
            +3d0*y2**2*y3**2*y4 -y1**2*3d0*y2**2*y5*SQRT3/2.D0 &
            +2d0*y2*y3**3*y4 )*FEA222334 +( +2d0*y2*y5**4 +3.D0*2d0*y2*y4**4 &
            +4.D0*2d0*y2*y4*y5**3*SQRT3 )*FEA335555 +( y1**3*3d0*y2**2 &
            +3d0*y2**2*y3**3 )*FEA222333 
 
         ds4(2)=ds3(2) &
            +( -y4**4*y5*SQRT3/2.D0 -3.D0/10.D0*y5**5*SQRT3 +y4**3*y5**2 &
            +y4**5/5.D0 )*FEA244455 &
            +( 5d0*y2**4*y4 -SQRT3*5d0*y2**4*y5 )*FEA222224 
 
         ds5(2)=ds4(2) +( &
            +y5**5*SQRT3/5.D0 -7.D0/15.D0*y4**5 +y4**4*y5*SQRT3/3.D0 &
            +y4*y5**4 )*FEA145555 
 
         ds1(2)=ds5(2) +( -SQRT3*y1*y3**3*y5/2.D0 &
            +y1**3*y3*y4 &
            +SQRT3*y1*3d0*y2**2*y3*y5/2.D0 -y1*3d0*y2**2*y3*y4/2.D0 &
            -y1*y3**3*y4/2.D0 )*FEA111234 &
            +( -y4**4*y5/3.D0 -y4*y5**4*SQRT3/2.D0 +y4**5*SQRT3/18.D0 &
            +y4**2*y5**3 )*FEA244555 &
            +( y1*y4**2*y5**2 -3.D0/4.D0*y3*y4**4 -y1*y5**4 &
            +5.D0/4.D0*y3*y5**4 -7.D0/2.D0*y3*y4**2*y5**2 &
            -2.D0*y1*y4**3*y5*SQRT3 )*FEA124455 
 
         ds3(2)=ds1(2) &
            +( 6d0*y2**5 )*FEA333333 +( y1*4d0*y2**3*y3 +y1**4*y3 &
            +y1*y3**4 )*FEA111123 +y1**2*2d0*y2*y3**2*FEA112233 +( &
            +4d0*y2**3*y4**2 +4d0*y2**3*y5**2 )*FEA222255 
 
         ds4(2)=ds3(2) +( &
            +9.D0*y3*y4**2*y5**2 -3.D0/2.D0*y3*y5**4 +y1*y4**4 &
            +4.D0*y1*y4**3*y5*SQRT3 +3.D0*y1*y5**4 +5.D0/2.D0*y3*y4**4 )*FEA134444 &
            +( &
            +5.D0/3.D0*y1*2d0*y2*y4**2*y5*SQRT3 -13.D0/3.D0*2d0*y2*y3*y4*y5**2 &
            -4.D0/3.D0*y3**2*y5**3*SQRT3 -7.D0/3.D0*y1**2*y4*y5**2 &
            +4.D0/3.D0*2d0*y2*y3*y5**3*SQRT3 +3.D0*y1**2*y4**3 +y3**2*y4**3 &
            +y1*2d0*y2*y5**3*SQRT3/3.D0 &
            +2d0*y2*y3*y4**3 -13.D0/3.D0*y3**2*y4*y5**2 &
            -4.D0/3.D0*y1**2*y4**2*y5*SQRT3 &
            -16.D0/3.D0*y1*2d0*y2*y4*y5**2 )*FEA233444 
 
         ds5(2)=ds4(2) &
            +( +4.D0*y3**2*y5**3 &
            +4.D0*2d0*y2*y3*y4*y5**2*SQRT3 -2.D0*y1*2d0*y2*y5**3 &
            -6.D0*y1*2d0*y2*y4**2*y5 &
            +y1**2*y4*y5**2*SQRT3 -3.D0*y1**2*y4**3*SQRT3 -4.D0*2d0*y2*y3*y5**3 &
            +3.D0*y1**2*y4**2*y5 -y1**2*y5**3 +4.D0*y3**2*y4*y5**2*SQRT3 &
            +4.D0*y1*2d0*y2*y4*y5**2*SQRT3 )*FEA113555 
 
         ds2(2)=ds5(2) &
            +( -2d0*y2*y4**3*y5 -2.D0/3.D0*2d0*y2*y4**4*SQRT3 &
            -3.D0*2d0*y2*y4*y5**3 )*FEA334445 &
            +( +2.D0/3.D0*y1*y5**4*SQRT3 &
            +3.D0*y1*y4**3*y5 -7.D0/12.D0*y3*y5**4*SQRT3 &
            +3.D0/2.D0*y3*y4**2*y5**2*SQRT3 +y1*y4*y5**3 &
            +3.D0/4.D0*y3*y4**4*SQRT3 )*FEA124555 +( &
            +2d0*y2*y4**2*y5**2 -2d0*y2*y4**4 &
            -2.D0*2d0*y2*y4*y5**3*SQRT3 )*FEA334455 
 
         ds3(2)=ds2(2) &
            +0d0 +( y3**3*y4**2 +y3**3*y5**2 +y1*3d0*y2**2*y4**2 +y1**3*y4**2 &
            +y1*3d0*y2**2*y5**2 +y1**3*y5**2 +3d0*y2**2*y3*y4**2 &
            +3d0*y2**2*y3*y5**2 )*FEA233344 &
            +( y1*3d0*y2**2*y5**2*SQRT3/6.D0 -3d0*y2**2*y3*y5**2*SQRT3/3.D0 &
            -y3**3*y5**2*SQRT3/3.D0 &
            +y1**3*y4*y5 -y1**3*y5**2*SQRT3/3.D0 -3d0*y2**2*y3*y4*y5 &
            +y3**3*y4*y5 -y1*3d0*y2**2*y4**2*SQRT3/2.D0 )*FEA233345 &
            +( -3.D0*3d0*y2**2*y4*y5**2 +3d0*y2**2*y4**3 )*FEA111444 &
            +( y1*3d0*y2**2*y3**2 +y1**3*2d0*y2*y3 +y1**2*3d0*y2**2*y3 &
            +y1*2d0*y2*y3**3 +y1**2*y3**3 +y1**3*y3**2 )*FEA111233 
 
         ds4(2)=ds3(2) &
            +0d0 +( -5.D0/3.D0*y1*2d0*y2*y4**2*y5*SQRT3 &
            +y1*2d0*y2*y4**3 -2.D0*y1**2*y4**3 -y1*2d0*y2*y5**3*SQRT3/3.D0 &
            +4.D0/3.D0*2d0*y2*y3*y4*y5**2 -4.D0/3.D0*2d0*y2*y3*y5**3*SQRT3 &
            +7.D0/3.D0*y1*2d0*y2*y4*y5**2 +4.D0/3.D0*y3**2*y5**3*SQRT3 &
            +4.D0/3.D0*y1**2*y4**2*y5*SQRT3 &
            +4.D0/3.D0*y3**2*y4*y5**2 &
            -2.D0/3.D0*y1**2*y4*y5**2 )*FEA133444 
 
         ds5(2)=ds4(2) &
            +( -y1**3*y4*y5 +2.D0/3.D0*3d0*y2**2*y3*y5**2*SQRT3 &
            +y1**3*y5**2*SQRT3/6.D0 +y1*3d0*y2**2*y5**2*SQRT3/6.D0 &
            +y1**3*y4**2*SQRT3/2.D0 &
            +2.D0/3.D0*y3**3*y5**2*SQRT3 -y1*3d0*y2**2*y4*y5 &
            +y1*3d0*y2**2*y4**2*SQRT3/2.D0 )*FEA133345 
 
         dv6(2)=ds5(2) &
            +( -2d0*y2*y3*y4**2*y5 +y3**2*y4**2*y5 +y3**2*y5**3 -y1*2d0*y2*y5**3 &
            +4.D0/3.D0*2d0*y2*y3*y4*y5**2*SQRT3 &
            +4.D0/3.D0*y3**2*y4*y5**2*SQRT3 -y1*2d0*y2*y4**2*y5 -2d0*y2*y3*y5**3 &
            +y1**2*y4*y5**2*SQRT3/3.D0 -y1**2*y4**3*SQRT3 &
            +4.D0/3.D0*y1*2d0*y2*y4*y5**2*SQRT3 )*FEA233445 &
            +( y3**4*y4*SQRT3 -y1**4*y5 +4d0*y2**3*y3*y4*SQRT3 &
            +y3**4*y5 -2.D0*y1*4d0*y2**3*y5 -y1**4*y4*SQRT3 &
            -4d0*y2**3*y3*y5 )*FEA233335 &
            +( 2d0*y2*y3**4 +y1**2*4d0*y2**3 +4d0*y2**3*y3**2 &
            +y1**4*2d0*y2 )*FEA222233 
 
         dv(2)= +dv1(2) +dv2(2) +dv3(2) +dv4(2) &
            +dv5(2) +dv6(2)
 

         dv1(3)=( 1d0 )*FEA1 
 
         dv2(3)=( y2 +y1 )*FEA12 +( +2d0*y3 )*FEA11 &
            +( -SQRT3*y5/2.D0 -y4/2.D0 )*FEA14 &
            +0d0 
 
         dv3(3)=( y1*y4 -2.D0*y2*y4 -SQRT3*y1*y5 )*FEA124 &
            +( 3.D0/4.D0*y4**2 -SQRT3*y4*y5/2.D0 +y5**2/4.D0 )*FEA155 +( y2*2d0*y3 &
            +y1*2d0*y3 +y1**2 +y2**2 )*FEA112 +0d0 +y1*y2*FEA123 +( &
            +3.D0/4.D0*y5**2 +SQRT3*y4*y5/2.D0 +y4**2/4.D0 )*FEA144 &
            +( 3d0*y3**2 )*FEA111 &
            +( -2d0*y3*y4/2.D0 -SQRT3*2d0*y3*y5/2.D0 )*FEA114 
 
         ds2(3)=0d0 &
            +( -3.D0/8.D0*SQRT3*y4**2*y5 -3.D0/8.D0*SQRT3*y5**3 -y4**3/8.D0 &
            -9.D0/8.D0*y4*y5**2 )*FEA1444 &
            +( +3.D0/4.D0*2d0*y3*y4**2 &
            +2d0*y3*y5**2/4.D0 -SQRT3*2d0*y3*y4*y5/2.D0 )*FEA1155 
 
         ds1(3)=ds2(3) &
            +( 2d0*y3*y4**2/4.D0 +3.D0/4.D0*2d0*y3*y5**2 &
            +SQRT3*2d0*y3*y4*y5/2.D0 )*FEA1144 &
            +( -SQRT3*3d0*y3**2*y5/2.D0 -3d0*y3**2*y4/2.D0 )*FEA1114 +( &
            +4d0*y3**3 )*FEA1111 +( SQRT3*y1*y4*y5 &
            +3.D0/2.D0*y2*y5**2 -y2*y4**2/2.D0 +y1*y4**2 )*FEA1244 
 
         ds2(3)=ds1(3) &
            +( y1*y5**2 -SQRT3*y1*y4*y5 -y2*y5**2/2.D0 &
            +3.D0/2.D0*y2*y4**2 )*FEA1255 +( -y1*2d0*y3*y4/2.D0 &
            +y1**2*y4 -SQRT3*y1*2d0*y3*y5/2.D0 -SQRT3*y2*2d0*y3*y5/2.D0 &
            +SQRT3*y2**2*y5/2.D0 -y2**2*y4/2.D0 -y2*2d0*y3*y4/2.D0 )*FEA1124 +( &
            +SQRT3*y1*2d0*y3*y4/2.D0 -SQRT3*y2*2d0*y3*y4/2.D0 -SQRT3*y2**2*y4/2.D0 &
            -y2**2*y5/2.D0 &
            +y2*2d0*y3*y5/2.D0 -y1*2d0*y3*y5/2.D0 -y1**2*y5 )*FEA1125 
 
         dv4(3)=ds2(3) &
            +( y2*3d0*y3**2 +y1**3 +y1*3d0*y3**2 +y2**3 )*FEA1112 +( y2**2*2d0*y3 &
            +y1**2*2d0*y3 )*FEA1122 +( y1*y2**2 +y1**2*y2 +y1*y2*2d0*y3 )*FEA1123 &
            +( -SQRT3*y4**2*y5/8.D0 -SQRT3*y5**3/8.D0 &
            +5.D0/8.D0*y4*y5**2 -3.D0/8.D0*y4**3 )*FEA1455 
 
         ds3(3)=0d0 &
            +( -4.D0*y4*y5**3*SQRT3 +3.D0*y4**4 +y5**4 )*FEA25555 +( &
            +y4**2*y5**2 -y4**4 +2.D0*y4*y5**3*SQRT3 )*FEA24455 
 
         ds2(3)=ds3(3) &
            +( -3.D0*y4*y5**3 +2.D0/3.D0*y4**4*SQRT3 -y4**3*y5 )*FEA24445 +( &
            +2d0*y3*y4**2*y5 +2d0*y3*y5**3 &
            +4.D0/9.D0*2d0*y3*y4**3*SQRT3 )*FEA33445 &
            +( 2d0*y3*y4*y5**2 -2d0*y3*y4**3/3.D0 )*FEA33455 
 
         ds1(3)=ds2(3) +( &
            +3d0*y3**2*y4*y5 +3d0*y3**2*y5**2*SQRT3/3.D0 )*FEA33345 &
            +( 3d0*y3**2*y4**2 +3d0*y3**2*y5**2 )*FEA33344 +( 4d0*y3**3*y4 &
            +SQRT3*4d0*y3**3*y5 )*FEA33334 +( +5d0*y3**4 )*FEA33333 +( &
            +y1*y4**2*y5 +y2*y4*y5**2*SQRT3 &
            +5.D0/9.D0*y2*y4**3*SQRT3 -4.D0/9.D0*y1*y4**3*SQRT3 &
            +y1*y5**3 )*FEA13445 +( y2*y4*y5**2 -y2*y4**3/3.D0 -y1*y4**3/3.D0 &
            +y1*y4*y5**2 )*FEA13455 
 
         ds3(3)=ds1(3) +( y1**2*y5**2 +y2**2*y4**2 &
            +y2**2*y5**2 +y2*2d0*y3*y4**2 +y1*2d0*y3*y4**2 +y1**2*y4**2 &
            +y1*2d0*y3*y5**2 +y2*2d0*y3*y5**2 )*FEA11255 &
            +( 2.D0/3.D0*y1**2*y4**2*SQRT3 +y1*2d0*y3*y5**2*SQRT3/2.D0 &
            +y2**2*y5**2*SQRT3/2.D0 +y2*2d0*y3*y4*y5 +y1*2d0*y3*y4*y5 -y2**2*y4*y5 &
            +y2*2d0*y3*y4**2*SQRT3/6.D0 +y1*2d0*y3*y4**2*SQRT3/6.D0 &
            +y2*2d0*y3*y5**2*SQRT3/2.D0 &
            +y2**2*y4**2*SQRT3/6.D0 )*FEA13345 
 
         ds4(3)=ds3(3) +( &
            +y1**2*y4**2*SQRT3/3.D0 +y2*2d0*y3*y4*y5 -y2**2*y4*y5 -y1**2*y4*y5 &
            +y2*2d0*y3*y4**2*SQRT3/3.D0 -y1*2d0*y3*y4**2*SQRT3/6.D0 &
            +y2**2*y4**2*SQRT3/3.D0 &
            +y1*2d0*y3*y5**2*SQRT3/2.D0 )*FEA11245 
 
         ds2(3)=ds4(3) +( +y1**3*y5 &
            +y2**3*y5/2.D0 -y2*3d0*y3**2*y5/2.D0 +y1*3d0*y3**2*y5/2.D0 &
            +y2**3*y4*SQRT3/2.D0 &
            +y2*3d0*y3**2*y4*SQRT3/2.D0 -y1*3d0*y3**2*y4*SQRT3/2.D0 )*FEA11135 &
            +( y1**3*y4 -y2**3*y4/2.D0 -y2*3d0*y3**2*y4/2.D0 -y1*3d0*y3**2*y4/2.D0 &
            +y2**3*y5*SQRT3/2.D0 -y2*3d0*y3**2*y5*SQRT3/2.D0 &
            -y1*3d0*y3**2*y5*SQRT3/2.D0 )*FEA11134 
 
         dv5(3)=ds2(3) &
            +( +y1**4 +y2**4 +y2*4d0*y3**3 +y1*4d0*y3**3 )*FEA23333 &
            +( -2.D0*y2**2*2d0*y3*y4 -SQRT3*y1**2*2d0*y3*y5 &
            +y1**2*2d0*y3*y4 )*FEA11334 +( y1**2*3d0*y3**2 +y1**3*2d0*y3 &
            +y2**2*3d0*y3**2 +y2**3*2d0*y3 )*FEA11333 +( y1*y2*y4**2 &
            +y1*y2*y5**2 )*FEA12355 &
            +( -y1*y2*2d0*y3*y4/2.D0 -y1*y2**2*y4/2.D0 -SQRT3*y1*y2*2d0*y3*y5/2.D0 &
            +y1**2*y2*y4 +SQRT3*y1*y2**2*y5/2.D0 )*FEA11234 +( y1*y2**3 &
            +y1*y2*3d0*y3**2 +y1**3*y2 )*FEA11123 +( y1**2*y2**2 +y1*y2**2*2d0*y3 &
            +y1**2*y2*2d0*y3 )*FEA11233 
 
         ds3(3)=( +3d0*y3**2*y4**2*y5 &
            +3d0*y3**2*y4**3*SQRT3 -5.D0/3.D0*3d0*y3**2*y4*y5**2*SQRT3 &
            +3d0*y3**2*y5**3 )*FEA333555 +( &
            +4d0*y3**3*y4**2*SQRT3/3.D0 -4d0*y3**3*y4*y5 )*FEA222245 &
            +( y1*5d0*y3**4 +y2**5 +y1**5 +y2*5d0*y3**4 )*FEA133333 &
            +( y1**4*y4 -2.D0*y2**4*y4 &
            +y1*4d0*y3**3*y4 -2.D0*y2*4d0*y3**3*y4 -y1*4d0*y3**3*y5*SQRT3 &
            -y1**4*y5*SQRT3 )*FEA133334 &
            +( -y1*y2*y4**3/3.D0 +y1*y2*y4*y5**2 )*FEA123455 
 
         ds4(3)=ds3(3) &
            +( 2.D0/3.D0*SQRT3*y1*y2**2*2d0*y3*y4 -y1**2*y2**2*y5 &
            -SQRT3*y1**2*y2**2*y4/3.D0 &
            +y1**2*y2*2d0*y3*y5 -SQRT3*y1**2*y2*2d0*y3*y4/3.D0 )*FEA112335 &
            +( y1*y2**2*y5**2 +y1*y2*2d0*y3*y5**2 +y1*y2*2d0*y3*y4**2 &
            +y1*y2**2*y4**2 +y1**2*y2*y4**2 &
            +y1**2*y2*y5**2 )*FEA112355 
 
         ds2(3)=ds4(3) &
            +( y2**3*2d0*y3*y5 -y1**2*3d0*y3**2*y5/2.D0 -y2**2*3d0*y3**2*y5 &
            +y1**3*2d0*y3*y5/2.D0 &
            +y1**3*2d0*y3*y4*SQRT3/2.D0 -y1**2*3d0*y3**2*y4*SQRT3/2.D0 )*FEA222335 &
            +( -y1**2*2d0*y3*y5**2*SQRT3/2.D0 -2.D0/3.D0*y2**2*2d0*y3*y4**2*SQRT3 &
            +y1**2*2d0*y3*y4*y5 -y1**2*2d0*y3*y4**2*SQRT3/6.D0 )*FEA113345 &
            +( y2**2*2d0*y3*y5**2 +y2**2*2d0*y3*y4**2 +y1**2*2d0*y3*y4**2 &
            +y1**2*2d0*y3*y5**2 )*FEA223355 
 
         ds3(3)=ds2(3) &
            +( y1*y2*2d0*y3*y4**2*SQRT3/6.D0 +y1*y2*2d0*y3*y4*y5 &
            +y1*y2*2d0*y3*y5**2*SQRT3/2.D0 &
            +2.D0/3.D0*y1**2*y2*y4**2*SQRT3 -y1*y2**2*y4*y5 &
            +y1*y2**2*y4**2*SQRT3/6.D0 +y1*y2**2*y5**2*SQRT3/2.D0 )*FEA123345 &
            +( -y1**3*2d0*y3*y4/2.D0 &
            +y1**3*2d0*y3*y5*SQRT3/2.D0 -y1**2*3d0*y3**2*y4/2.D0 +y2**3*2d0*y3*y4 &
            +y2**2*3d0*y3**2*y4 +y1**2*3d0*y3**2*y5*SQRT3/2.D0 )*FEA222334 &
            +( 3.D0*2d0*y3*y4**4 -4.D0*2d0*y3*y4*y5**3*SQRT3 &
            +2d0*y3*y5**4 )*FEA335555 +( +y1**3*3d0*y3**2 &
            +y2**3*3d0*y3**2 )*FEA222333 
 
         ds4(3)=ds3(3) +( y4**5/5.D0 +y4**3*y5**2 &
            +y4**4*y5*SQRT3/2.D0 +3.D0/10.D0*y5**5*SQRT3 )*FEA244455 +( &
            +5d0*y3**4*y4 +SQRT3*5d0*y3**4*y5 )*FEA222224 
 
         ds5(3)=ds4(3) &
            +( -y5**5*SQRT3/5.D0 -y4**4*y5*SQRT3/3.D0 &
            +y4*y5**4 -7.D0/15.D0*y4**5 )*FEA145555 
 
         ds1(3)=ds5(3) &
            +( -SQRT3*y1*y2*3d0*y3**2*y5/2.D0 +y1**3*y2*y4 &
            +SQRT3*y1*y2**3*y5/2.D0 -y1*y2**3*y4/2.D0 &
            -y1*y2*3d0*y3**2*y4/2.D0 )*FEA111234 &
            +( y4**4*y5/3.D0 &
            +y4**5*SQRT3/18.D0 -y4**2*y5**3 -y4*y5**4*SQRT3/2.D0 )*FEA244555 &
            +( -3.D0/4.D0*y2*y4**4 -y1*y5**4 +5.D0/4.D0*y2*y5**4 &
            +y1*y4**2*y5**2 -7.D0/2.D0*y2*y4**2*y5**2 &
            +2.D0*y1*y4**3*y5*SQRT3 )*FEA124455 
 
         ds3(3)=ds1(3) +( &
            +6d0*y3**5 )*FEA333333 +( y1*y2**4 +y1**4*y2 &
            +y1*y2*4d0*y3**3 )*FEA111123 +y1**2*y2**2*2d0*y3*FEA112233 +( &
            +4d0*y3**3*y4**2 +4d0*y3**3*y5**2 )*FEA222255 
 
         ds4(3)=ds3(3) &
            +( 3.D0*y1*y5**4 +y1*y4**4 &
            +9.D0*y2*y4**2*y5**2 -3.D0/2.D0*y2*y5**4 -4.D0*y1*y4**3*y5*SQRT3 &
            +5.D0/2.D0*y2*y4**4 )*FEA134444 &
            +( -y1*2d0*y3*y5**3*SQRT3/3.D0 -7.D0/3.D0*y1**2*y4*y5**2 &
            -13.D0/3.D0*y2**2*y4*y5**2 -4.D0/3.D0*y2*2d0*y3*y5**3*SQRT3 &
            -16.D0/3.D0*y1*2d0*y3*y4*y5**2 &
            +4.D0/3.D0*y1**2*y4**2*y5*SQRT3 +4.D0/3.D0*y2**2*y5**3*SQRT3 &
            +y2*2d0*y3*y4**3 &
            +y2**2*y4**3 -13.D0/3.D0*y2*2d0*y3*y4*y5**2 &
            -5.D0/3.D0*y1*2d0*y3*y4**2*y5*SQRT3 &
            +3.D0*y1**2*y4**3 )*FEA233444 
 
         ds5(3)=ds4(3) +( 2.D0*y1*2d0*y3*y5**3 &
            +4.D0*y2*2d0*y3*y5**3 +4.D0*y2**2*y4*y5**2*SQRT3 +y1**2*y4*y5**2*SQRT3 &
            +6.D0*y1*2d0*y3*y4**2*y5 -3.D0*y1**2*y4**2*y5 &
            +4.D0*y1*2d0*y3*y4*y5**2*SQRT3 -4.D0*y2**2*y5**3 &
            +y1**2*y5**3 -3.D0*y1**2*y4**3*SQRT3 &
            +4.D0*y2*2d0*y3*y4*y5**2*SQRT3 )*FEA113555 
 
         ds2(3)=ds5(3) &
            +( -2.D0/3.D0*2d0*y3*y4**4*SQRT3 +2d0*y3*y4**3*y5 &
            +3.D0*2d0*y3*y4*y5**3 )*FEA334445 +( -3.D0*y1*y4**3*y5 -y1*y4*y5**3 &
            +2.D0/3.D0*y1*y5**4*SQRT3 -7.D0/12.D0*y2*y5**4*SQRT3 &
            +3.D0/2.D0*y2*y4**2*y5**2*SQRT3 +3.D0/4.D0*y2*y4**4*SQRT3 )*FEA124555 &
            +( 2.D0*2d0*y3*y4*y5**3*SQRT3 -2d0*y3*y4**4 &
            +2d0*y3*y4**2*y5**2 )*FEA334455 
 
         ds3(3)=ds2(3) +0d0 &
            +( y2*3d0*y3**2*y4**2 +y2*3d0*y3**2*y5**2 +y1*3d0*y3**2*y4**2 &
            +y1**3*y5**2 +y1**3*y4**2 +y2**3*y4**2 +y1*3d0*y3**2*y5**2 &
            +y2**3*y5**2 )*FEA233344 &
            +( -y2**3*y5**2*SQRT3/3.D0 -y2*3d0*y3**2*y5**2*SQRT3/3.D0 -y1**3*y4*y5 &
            -y1**3*y5**2*SQRT3/3.D0 -y1*3d0*y3**2*y4**2*SQRT3/2.D0 &
            +y1*3d0*y3**2*y5**2*SQRT3/6.D0 -y2**3*y4*y5 &
            +y2*3d0*y3**2*y4*y5 )*FEA233345 +( &
            +3d0*y3**2*y4**3 -3.D0*3d0*y3**2*y4*y5**2 )*FEA111444 &
            +( y1*y2**3*2d0*y3 +y1**3*y2**2 +y1**2*y2**3 +y1*y2**2*3d0*y3**2 &
            +y1**2*y2*3d0*y3**2 +y1**3*y2*2d0*y3 )*FEA111233 
 
         ds4(3)=ds3(3) +0d0 &
            +( -4.D0/3.D0*y1**2*y4**2*y5*SQRT3 &
            +4.D0/3.D0*y2**2*y4*y5**2 -4.D0/3.D0*y2**2*y5**3*SQRT3 &
            -2.D0*y1**2*y4**3 -2.D0/3.D0*y1**2*y4*y5**2 &
            +y1*2d0*y3*y4**3 +4.D0/3.D0*y2*2d0*y3*y5**3*SQRT3 &
            +y1*2d0*y3*y5**3*SQRT3/3.D0 +4.D0/3.D0*y2*2d0*y3*y4*y5**2 &
            +5.D0/3.D0*y1*2d0*y3*y4**2*y5*SQRT3 &
            +7.D0/3.D0*y1*2d0*y3*y4*y5**2 )*FEA133444 
 
         ds5(3)=ds4(3) +( &
            +2.D0/3.D0*y2**3*y5**2*SQRT3 +y1*3d0*y3**2*y4**2*SQRT3/2.D0 &
            +y1**3*y4**2*SQRT3/2.D0 +y1**3*y5**2*SQRT3/6.D0 +y1**3*y4*y5 &
            +2.D0/3.D0*y2*3d0*y3**2*y5**2*SQRT3 +y1*3d0*y3**2*y5**2*SQRT3/6.D0 &
            +y1*3d0*y3**2*y4*y5 )*FEA133345 
 
         dv6(3)=ds5(3) +( -y2**2*y4**2*y5 &
            +y1**2*y4*y5**2*SQRT3/3.D0 +y2*2d0*y3*y4**2*y5 +y2*2d0*y3*y5**3 &
            +4.D0/3.D0*y2**2*y4*y5**2*SQRT3 +4.D0/3.D0*y2*2d0*y3*y4*y5**2*SQRT3 &
            +4.D0/3.D0*y1*2d0*y3*y4*y5**2*SQRT3 -y2**2*y5**3 +y1*2d0*y3*y5**3 &
            +y1*2d0*y3*y4**2*y5 -y1**2*y4**3*SQRT3 )*FEA233445 &
            +( y2*4d0*y3**3*y4*SQRT3 +y2**4*y4*SQRT3 -y1**4*y4*SQRT3 &
            +y2*4d0*y3**3*y5 +2.D0*y1*4d0*y3**3*y5 +y1**4*y5 -y2**4*y5 )*FEA233335 &
            +( y2**2*4d0*y3**3 +y1**4*2d0*y3 +y2**4*2d0*y3 &
            +y1**2*4d0*y3**3 )*FEA222233 
 
         dv(3)= +dv1(3) +dv2(3) +dv3(3) &
            +dv4(3) +dv5(3) +dv6(3)
 

         dv1(4)=0d0 
 
         dv2(4)=0d0 +0d0 +( -y3*1d0/2.D0 +y1 -y2*1d0/2.D0 )*FEA14 +( &
            +2d0*y4 )*FEA44 
 
         dv3(4)=( y1*y3 +y1*y2 -2.D0*y2*y3 )*FEA124 &
            +( 3.D0/4.D0*y3*2d0*y4 -SQRT3*y3*y5/2.D0 +3.D0/4.D0*y2*2d0*y4 &
            +SQRT3*y2*y5/2.D0 )*FEA155 +0d0 +( -3d0*y4**2/3.D0 +y5**2 )*FEA455 &
            +( y1*2d0*y4 +y2*2d0*y4/4.D0 -SQRT3*y2*y5/2.D0 +SQRT3*y3*y5/2.D0 &
            +y3*2d0*y4/4.D0 )*FEA144 +0d0 +( -y2**2*1d0/2.D0 -y3**2*1d0/2.D0 &
            +y1**2 )*FEA114 
 
         ds2(4)=( 4d0*y4**3 +2.D0*2d0*y4*y5**2 )*FEA4444 &
            +( -3.D0/8.D0*SQRT3*y3*2d0*y4*y5 -9.D0/8.D0*y2*y5**2 &
            -y3*3d0*y4**2/8.D0 -y2*3d0*y4**2/8.D0 -9.D0/8.D0*y3*y5**2 &
            +y1*3d0*y4**2 +3.D0/8.D0*SQRT3*y2*2d0*y4*y5 )*FEA1444 &
            +( 3.D0/4.D0*y2**2*2d0*y4 +3.D0/4.D0*y3**2*2d0*y4 -SQRT3*y3**2*y5/2.D0 &
            +SQRT3*y2**2*y5/2.D0 )*FEA1155 
 
         ds1(4)=ds2(4) +( y3**2*2d0*y4/4.D0 &
            +y1**2*2d0*y4 +y2**2*2d0*y4/4.D0 &
            +SQRT3*y3**2*y5/2.D0 -SQRT3*y2**2*y5/2.D0 )*FEA1144 &
            +( y1**3 -y2**3*1d0/2.D0 -y3**3*1d0/2.D0 )*FEA1114 +0d0 &
            +( SQRT3*y1*y3*y5 -y2*y3*2d0*y4/2.D0 +y1*y2*2d0*y4 -SQRT3*y1*y2*y5 &
            +y1*y3*2d0*y4 )*FEA1244 
 
         ds2(4)=ds1(4) +( -SQRT3*y1*y3*y5 &
            +3.D0/2.D0*y2*y3*2d0*y4 +SQRT3*y1*y2*y5 )*FEA1255 &
            +( -y1*y3**2*1d0/2.D0 +y1**2*y3 &
            +y1**2*y2 -y2**2*y3*1d0/2.D0 -y2*y3**2*1d0/2.D0 &
            -y1*y2**2*1d0/2.D0 )*FEA1124 &
            +( +SQRT3*y1*y3**2*1d0/2.D0 &
            +SQRT3*y1*y2**2*1d0/2.D0 -SQRT3*y2*y3**2*1d0/2.D0 &
            -SQRT3*y2**2*y3*1d0/2.D0 )*FEA1125 
 
         dv4(4)=ds2(4) &
            +0d0 +0d0 +0d0 +( 5.D0/8.D0*y2*y5**2 -SQRT3*y3*2d0*y4*y5/8.D0 &
            +SQRT3*y2*2d0*y4*y5/8.D0 -3.D0/8.D0*y2*3d0*y4**2 +y1*y5**2 &
            +5.D0/8.D0*y3*y5**2 &
            -3.D0/8.D0*y3*3d0*y4**2 )*FEA1455 
 
         ds3(4)=( 5d0*y4**4 &
            -2.D0*3d0*y4**2*y5**2 -3.D0*y5**4 )*FEA44444 &
            +( -4.D0*y3*y5**3*SQRT3 +9.D0*y1*2d0*y4*y5**2 -3.D0/2.D0*y1*4d0*y4**3 &
            +4.D0*y2*y5**3*SQRT3 +3.D0*y2*4d0*y4**3 +3.D0*y3*4d0*y4**3 )*FEA25555 &
            +( -y2*4d0*y4**3 &
            +y3*2d0*y4*y5**2 -2.D0*y2*y5**3*SQRT3 -y3*4d0*y4**3 &
            -7.D0/2.D0*y1*2d0*y4*y5**2 &
            +2.D0*y3*y5**3*SQRT3 +y2*2d0*y4*y5**2 &
            +5.D0/4.D0*y1*4d0*y4**3 )*FEA24455 
 
         ds2(4)=ds3(4) &
            +( y2*3d0*y4**2*y5 -3.D0*y3*y5**3 +2.D0/3.D0*y3*4d0*y4**3*SQRT3 &
            +3.D0*y2*y5**3 -7.D0/12.D0*y1*4d0*y4**3*SQRT3 &
            +3.D0/2.D0*y1*2d0*y4*y5**2*SQRT3 -y3*3d0*y4**2*y5 &
            +2.D0/3.D0*y2*4d0*y4**3*SQRT3 )*FEA24445 +( +y3**2*2d0*y4*y5 &
            +4.D0/9.D0*y2**2*3d0*y4**2*SQRT3 -5.D0/9.D0*y1**2*3d0*y4**2*SQRT3 &
            +4.D0/9.D0*y3**2*3d0*y4**2*SQRT3 -y2**2*2d0*y4*y5 &
            -y1**2*y5**2*SQRT3 )*FEA33445 &
            +( y3**2*y5**2 -y1**2*3d0*y4**2/3.D0 -y3**2*3d0*y4**2/3.D0 &
            +y1**2*y5**2 &
            +y2**2*y5**2 -y2**2*3d0*y4**2/3.D0 )*FEA33455 
 
         ds1(4)=ds2(4) &
            +( -y2**3*y5 +y3**3*y5 +y1**3*2d0*y4*SQRT3/2.D0 )*FEA33345 &
            +( y3**3*2d0*y4 +y2**3*2d0*y4 +y1**3*2d0*y4 )*FEA33344 +( y3**4 &
            +y2**4 -2.D0*y1**4 )*FEA33334 +0d0 +( -4.D0/9.D0*y1*y2*3d0*y4**2*SQRT3 &
            +y1*y3*2d0*y4*y5 +y2*y3*y5**2*SQRT3 -y1*y2*2d0*y4*y5 &
            +5.D0/9.D0*y2*y3*3d0*y4**2*SQRT3 &
            -4.D0/9.D0*y1*y3*3d0*y4**2*SQRT3 )*FEA13445 &
            +( y2*y3*y5**2 &
            +y1*y2*y5**2 -y2*y3*3d0*y4**2/3.D0 -y1*y2*3d0*y4**2/3.D0 &
            -y1*y3*3d0*y4**2/3.D0 &
            +y1*y3*y5**2 )*FEA13455 
 
         ds3(4)=ds1(4) +( +y2**2*y3*2d0*y4 &
            +y1*y2**2*2d0*y4 +y2*y3**2*2d0*y4 +y1*y3**2*2d0*y4 +y1**2*y3*2d0*y4 &
            +y1**2*y2*2d0*y4 )*FEA11255 &
            +( 2.D0/3.D0*y1**2*y3*2d0*y4*SQRT3 -y1*y2**2*y5 +y2*y3**2*y5 &
            +y1*y3**2*y5 -y2**2*y3*y5 +y2*y3**2*2d0*y4*SQRT3/6.D0 &
            +y1*y3**2*2d0*y4*SQRT3/6.D0 +y1*y2**2*2d0*y4*SQRT3/6.D0 &
            +2.D0/3.D0*y1**2*y2*2d0*y4*SQRT3 &
            +y2**2*y3*2d0*y4*SQRT3/6.D0 )*FEA13345 
 
         ds4(4)=ds3(4) +( y1**2*y2*y5 &
            +y1**2*y3*2d0*y4*SQRT3/3.D0 &
            +y1**2*y2*2d0*y4*SQRT3/3.D0 -y1*y2**2*2d0*y4*SQRT3/6.D0 &
            +y2*y3**2*y5 -y2**2*y3*y5 -y1**2*y3*y5 &
            +y2*y3**2*2d0*y4*SQRT3/3.D0 -y1*y3**2*2d0*y4*SQRT3/6.D0 &
            +y2**2*y3*2d0*y4*SQRT3/3.D0 )*FEA11245 
 
         ds2(4)=ds4(4) &
            +( -y1*y2**3*SQRT3/2.D0 +y2**3*y3*SQRT3/2.D0 &
            +y2*y3**3*SQRT3/2.D0 -y1*y3**3*SQRT3/2.D0 )*FEA11135 &
            +( y1**3*y3 -y2**3*y3*1d0/2.D0 &
            +y1**3*y2 -y2*y3**3*1d0/2.D0 -y1*y3**3*1d0/2.D0 &
            -y1*y2**3*1d0/2.D0 )*FEA11134 
 
         dv5(4)=ds2(4) &
            +0d0 +( -2.D0*y2**2*y3**2 +y1**2*y2**2 +y1**2*y3**2 )*FEA11334 +0d0 &
            +( y1*y2*y3*2d0*y4 )*FEA12355 &
            +( -y1*y2*y3**2*1d0/2.D0 -y1*y2**2*y3*1d0/2.D0 +y1**2*y2*y3 )*FEA11234 &
            +0d0 +0d0 
 
         ds3(4)=( y2**3*3d0*y4**2*SQRT3 -y2**3*2d0*y4*y5 &
            +y3**3*2d0*y4*y5 -5.D0/3.D0*y2**3*y5**2*SQRT3 &
            +y3**3*3d0*y4**2*SQRT3 -5.D0/3.D0*y3**3*y5**2*SQRT3 &
            -8.D0/3.D0*y1**3*y5**2*SQRT3 )*FEA333555 &
            +( +y2**4*y5 +y2**4*2d0*y4*SQRT3/3.D0 &
            +y3**4*2d0*y4*SQRT3/3.D0 -y3**4*y5 -y1**4*2d0*y4*SQRT3/6.D0 )*FEA222245 &
            +0d0 +( y1**4*y3 -2.D0*y2**4*y3 +y1**4*y2 +y1*y3**4 -2.D0*y2*y3**4 &
            +y1*y2**4 )*FEA133334 +( -y1*y2*y3*3d0*y4**2/3.D0 &
            +y1*y2*y3*y5**2 )*FEA123455 
 
         ds4(4)=ds3(4) &
            +( 2.D0/3.D0*SQRT3*y1*y2**2*y3**2 -SQRT3*y1**2*y2**2*y3*1d0/3.D0 &
            -SQRT3*y1**2*y2*y3**2*1d0/3.D0 )*FEA112335 &
            +( +y1*y2*y3**2*2d0*y4 +y1*y2**2*y3*2d0*y4 &
            +y1**2*y2*y3*2d0*y4 )*FEA112355 
 
         ds2(4)=ds4(4) +( &
            +y1**3*y2**2*SQRT3/2.D0 -y1**2*y2**3*SQRT3/2.D0 &
            +y1**3*y3**2*SQRT3/2.D0 -y1**2*y3**3*SQRT3/2.D0 )*FEA222335 &
            +( -y1**2*y2**2*2d0*y4*SQRT3/6.D0 -y1**2*y2**2*y5 &
            -2.D0/3.D0*y2**2*y3**2*2d0*y4*SQRT3 &
            +y1**2*y3**2*y5 -y1**2*y3**2*2d0*y4*SQRT3/6.D0 )*FEA113345 +( &
            +y2**2*y3**2*2d0*y4 +y1**2*y3**2*2d0*y4 &
            +y1**2*y2**2*2d0*y4 )*FEA223355 
 
         ds3(4)=ds2(4) &
            +( y1*y2*y3**2*2d0*y4*SQRT3/6.D0 +y1*y2*y3**2*y5 &
            +2.D0/3.D0*y1**2*y2*y3*2d0*y4*SQRT3 -y1*y2**2*y3*y5 &
            +y1*y2**2*y3*2d0*y4*SQRT3/6.D0 )*FEA123345 &
            +( -y1**3*y2**2*1d0/2.D0 -y1**3*y3**2*1d0/2.D0 -y1**2*y2**3*1d0/2.D0 &
            -y1**2*y3**3*1d0/2.D0 &
            +y2**3*y3**2 +y2**2*y3**3 )*FEA222334 +( 3.D0*y3**2*4d0*y4**3 &
            +3.D0*y2**2*4d0*y4**3 -4.D0*y3**2*y5**3*SQRT3 &
            +9.D0*y1**2*2d0*y4*y5**2 -3.D0/2.D0*y1**2*4d0*y4**3 &
            +4.D0*y2**2*y5**3*SQRT3 )*FEA335555 +0d0 
 
         ds4(4)=ds3(4) &
            +( y3*5d0*y4**4/5.D0 -y2*4d0*y4**3*y5*SQRT3/2.D0 &
            -2.D0/5.D0*y1*5d0*y4**4 -2.D0*y1*3d0*y4**2*y5**2 &
            +y3*3d0*y4**2*y5**2 +y3*4d0*y4**3*y5*SQRT3/2.D0 +y2*3d0*y4**2*y5**2 &
            +y2*5d0*y4**4/5.D0 )*FEA244455 +( y2**5 -2.D0*y1**5 &
            +y3**5 )*FEA222224 
 
         ds5(4)=ds4(4) +( +y1*y5**4 -7.D0/15.D0*y2*5d0*y4**4 &
            +y2*4d0*y4**3*y5*SQRT3/3.D0 -y3*4d0*y4**3*y5*SQRT3/3.D0 +y3*y5**4 &
            +y2*y5**4 &
            +2.D0*y1*3d0*y4**2*y5**2 -7.D0/15.D0*y3*5d0*y4**4 &
            -y1*5d0*y4**4/15.D0 )*FEA145555 
 
         ds1(4)=ds5(4) &
            +( &
            +y1**3*y2*y3 -y1*y2**3*y3*1d0/2.D0 -y1*y2*y3**3*1d0/2.D0 )*FEA111234 &
            +( y3*4d0*y4**3*y5/3.D0 &
            +y3*5d0*y4**4*SQRT3/18.D0 -y2*4d0*y4**3*y5/3.D0 -y2*y5**4*SQRT3/2.D0 &
            -y3*2d0*y4*y5**3 &
            +2.D0/9.D0*y1*5d0*y4**4*SQRT3 +y2*5d0*y4**4*SQRT3/18.D0 &
            +y2*2d0*y4*y5**3 -2.D0/3.D0*y1*3d0*y4**2*y5**2*SQRT3 &
            -y3*y5**4*SQRT3/2.D0 )*FEA244555 &
            +( y1*y2*2d0*y4*y5**2 -3.D0/4.D0*y2*y3*4d0*y4**3 &
            +y1*y3*2d0*y4*y5**2 -7.D0/2.D0*y2*y3*2d0*y4*y5**2 &
            -2.D0*y1*y2*3d0*y4**2*y5*SQRT3 &
            +2.D0*y1*y3*3d0*y4**2*y5*SQRT3 )*FEA124455 
 
         ds3(4)=ds1(4) +0d0 +0d0 &
            +( y1**4*2d0*y4 +y2**4*2d0*y4 +y3**4*2d0*y4 )*FEA222255 
 
         ds4(4)=ds3(4) &
            +( +y1*y3*4d0*y4**3 &
            +9.D0*y2*y3*2d0*y4*y5**2 -4.D0*y1*y3*3d0*y4**2*y5*SQRT3 &
            +y1*y2*4d0*y4**3 +4.D0*y1*y2*3d0*y4**2*y5*SQRT3 &
            +5.D0/2.D0*y2*y3*4d0*y4**3 )*FEA134444 +( -7.D0/3.D0*y1**2*y3*y5**2 &
            +5.D0/3.D0*y1*y2**2*2d0*y4*y5*SQRT3 -13.D0/3.D0*y2**2*y3*y5**2 &
            -7.D0/3.D0*y1**2*y2*y5**2 -16.D0/3.D0*y1*y3**2*y5**2 &
            +4.D0/3.D0*y1**2*y3*2d0*y4*y5*SQRT3 +3.D0*y1**2*y2*3d0*y4**2 &
            +y2*y3**2*3d0*y4**2 &
            +y2**2*y3*3d0*y4**2 -13.D0/3.D0*y2*y3**2*y5**2 &
            -5.D0/3.D0*y1*y3**2*2d0*y4*y5*SQRT3 &
            -4.D0/3.D0*y1**2*y2*2d0*y4*y5*SQRT3 &
            +3.D0*y1**2*y3*3d0*y4**2 &
            -16.D0/3.D0*y1*y2**2*y5**2 )*FEA233444 
 
         ds5(4)=ds4(4) &
            +( +4.D0*y2**2*y3*y5**2*SQRT3 +y1**2*y3*y5**2*SQRT3 &
            +6.D0*y1*y3**2*2d0*y4*y5 -6.D0*y1*y2**2*2d0*y4*y5 &
            -3.D0*y1**2*y3*2d0*y4*y5 &
            +y1**2*y2*y5**2*SQRT3 &
            +4.D0*y1*y3**2*y5**2*SQRT3 -3.D0*y1**2*y2*3d0*y4**2*SQRT3 &
            +3.D0*y1**2*y2*2d0*y4*y5 -3.D0*y1**2*y3*3d0*y4**2*SQRT3 &
            +4.D0*y2*y3**2*y5**2*SQRT3 &
            +4.D0*y1*y2**2*y5**2*SQRT3 )*FEA113555 
 
         ds2(4)=ds5(4) &
            +( -2.D0/3.D0*y3**2*4d0*y4**3*SQRT3 &
            -3.D0/2.D0*y1**2*2d0*y4*y5**2*SQRT3 -y2**2*3d0*y4**2*y5 &
            +7.D0/12.D0*y1**2*4d0*y4**3*SQRT3 +y3**2*3d0*y4**2*y5 &
            +3.D0*y3**2*y5**3 -2.D0/3.D0*y2**2*4d0*y4**3*SQRT3 &
            -3.D0*y2**2*y5**3 )*FEA334445 &
            +( -3.D0*y1*y3*3d0*y4**2*y5 -y1*y3*y5**3 +3.D0*y1*y2*3d0*y4**2*y5 &
            +3.D0/2.D0*y2*y3*2d0*y4*y5**2*SQRT3 +y1*y2*y5**3 &
            +3.D0/4.D0*y2*y3*4d0*y4**3*SQRT3 )*FEA124555 &
            +( 2.D0*y3**2*y5**3*SQRT3 -7.D0/2.D0*y1**2*2d0*y4*y5**2 &
            +y2**2*2d0*y4*y5**2 -y2**2*4d0*y4**3 -y3**2*4d0*y4**3 &
            -2.D0*y2**2*y5**3*SQRT3 &
            +5.D0/4.D0*y1**2*4d0*y4**3 &
            +y3**2*2d0*y4*y5**2 )*FEA334455 
 
         ds3(4)=ds2(4) +( -6.D0*2d0*y4*y5**4 &
            +9.D0*4d0*y4**3*y5**2 )*FEA555555 +( y2*y3**3*2d0*y4 +y1*y3**3*2d0*y4 &
            +y1*y2**3*2d0*y4 +y1**3*y2*2d0*y4 +y1**3*y3*2d0*y4 &
            +y2**3*y3*2d0*y4 )*FEA233344 +( &
            +y1**3*y2*y5 -y1**3*y3*y5 -y1*y3**3*2d0*y4*SQRT3/2.D0 -y2**3*y3*y5 &
            +y2*y3**3*y5 -y1*y2**3*2d0*y4*SQRT3/2.D0 )*FEA233345 &
            +( -3.D0*y2**3*y5**2 &
            +y3**3*3d0*y4**2 -3.D0*y3**3*y5**2 -3.D0*y1**3*y5**2 +y2**3*3d0*y4**2 &
            +y1**3*3d0*y4**2 )*FEA111444 +0d0 
 
         ds4(4)=ds3(4) &
            +( 9.D0*2d0*y4*y5**4 -6.D0*4d0*y4**3*y5**2 +6d0*y4**5 )*FEA444444 &
            +( -5.D0/3.D0*y1*y2**2*2d0*y4*y5*SQRT3 &
            +y1*y2**2*3d0*y4**2 -4.D0/3.D0*y1**2*y3*2d0*y4*y5*SQRT3 &
            -2.D0*y1**2*y2*3d0*y4**2 &
            +4.D0/3.D0*y2**2*y3*y5**2 -2.D0*y1**2*y3*3d0*y4**2 &
            +7.D0/3.D0*y1*y2**2*y5**2 -2.D0/3.D0*y1**2*y3*y5**2 &
            +y1*y3**2*3d0*y4**2 +4.D0/3.D0*y1**2*y2*2d0*y4*y5*SQRT3 &
            +4.D0/3.D0*y2*y3**2*y5**2 &
            +5.D0/3.D0*y1*y3**2*2d0*y4*y5*SQRT3 -2.D0/3.D0*y1**2*y2*y5**2 &
            +7.D0/3.D0*y1*y3**2*y5**2 )*FEA133444 
 
         ds5(4)=ds4(4) +( -y1**3*y2*y5 &
            +y1*y3**3*2d0*y4*SQRT3/2.D0 +y1**3*y3*2d0*y4*SQRT3/2.D0 +y1**3*y3*y5 &
            +y1**3*y2*2d0*y4*SQRT3/2.D0 -y1*y2**3*y5 +y1*y2**3*2d0*y4*SQRT3/2.D0 &
            +y1*y3**3*y5 )*FEA133345 
 
         dv6(4)=ds5(4) +( -y2**2*y3*2d0*y4*y5 &
            +y1**2*y3*y5**2*SQRT3/3.D0 +y2*y3**2*2d0*y4*y5 &
            +4.D0/3.D0*y2**2*y3*y5**2*SQRT3 &
            +4.D0/3.D0*y2*y3**2*y5**2*SQRT3 -y1*y2**2*2d0*y4*y5 &
            +4.D0/3.D0*y1*y3**2*y5**2*SQRT3 &
            +y1**2*y2*y5**2*SQRT3/3.D0 -y1**2*y2*3d0*y4**2*SQRT3 &
            +y1*y3**2*2d0*y4*y5 -y1**2*y3*3d0*y4**2*SQRT3 &
            +4.D0/3.D0*y1*y2**2*y5**2*SQRT3 )*FEA233445 +( y2*y3**4*SQRT3 &
            +y2**4*y3*SQRT3 -y1**4*y3*SQRT3 -y1**4*y2*SQRT3 )*FEA233335 &
            +0d0 
 
         dv(4)= +dv1(4) +dv2(4) +dv3(4) +dv4(4) +dv5(4) +dv6(4)
 

         dv1(5)=0d0 
 
         dv2(5)=0d0 +0d0 +( -SQRT3*y3*1d0/2.D0 &
            +SQRT3*y2*1d0/2.D0 )*FEA14 +( 2d0*y5 )*FEA44 
 
         dv3(5)=( &
            +SQRT3*y1*y2 -SQRT3*y1*y3 )*FEA124 +( -SQRT3*y3*y4*1d0/2.D0 +y1*2d0*y5 &
            +y2*2d0*y5/4.D0 +SQRT3*y2*y4*1d0/2.D0 +y3*2d0*y5/4.D0 )*FEA155 +0d0 +( &
            +y4*2d0*y5 )*FEA455 +( +3.D0/4.D0*y3*2d0*y5 &
            +3.D0/4.D0*y2*2d0*y5 -SQRT3*y2*y4*1d0/2.D0 &
            +SQRT3*y3*y4*1d0/2.D0 )*FEA144 +0d0 +( &
            +SQRT3*y2**2*1d0/2.D0 -SQRT3*y3**2*1d0/2.D0 )*FEA114 
 
         ds2(5)=( &
            +4d0*y5**3 +2.D0*y4**2*2d0*y5 )*FEA4444 &
            +( 3.D0/8.D0*SQRT3*y2*3d0*y5**2 -3.D0/8.D0*SQRT3*y3*y4**2 &
            -3.D0/8.D0*SQRT3*y3*3d0*y5**2 -9.D0/8.D0*y2*y4*2d0*y5 &
            -9.D0/8.D0*y3*y4*2d0*y5 &
            +3.D0/8.D0*SQRT3*y2*y4**2 )*FEA1444 +( +y1**2*2d0*y5 &
            +y3**2*2d0*y5/4.D0 -SQRT3*y3**2*y4*1d0/2.D0 +SQRT3*y2**2*y4*1d0/2.D0 &
            +y2**2*2d0*y5/4.D0 )*FEA1155 
 
         ds1(5)=ds2(5) +( +3.D0/4.D0*y3**2*2d0*y5 &
            +SQRT3*y3**2*y4*1d0/2.D0 -SQRT3*y2**2*y4*1d0/2.D0 &
            +3.D0/4.D0*y2**2*2d0*y5 )*FEA1144 +( &
            +SQRT3*y2**3*1d0/2.D0 -SQRT3*y3**3*1d0/2.D0 )*FEA1114 +0d0 &
            +( SQRT3*y1*y3*y4 &
            +3.D0/2.D0*y2*y3*2d0*y5 -SQRT3*y1*y2*y4 )*FEA1244 
 
         ds2(5)=ds1(5) &
            +( y1*y3*2d0*y5 +y1*y2*2d0*y5 -SQRT3*y1*y3*y4 -y2*y3*2d0*y5/2.D0 &
            +SQRT3*y1*y2*y4 )*FEA1255 &
            +( -SQRT3*y1*y3**2*1d0/2.D0 -SQRT3*y2*y3**2*1d0/2.D0 &
            +SQRT3*y2**2*y3*1d0/2.D0 +SQRT3*y1*y2**2*1d0/2.D0 )*FEA1124 &
            +( y1**2*y2 -y2**2*y3*1d0/2.D0 +y2*y3**2*1d0/2.D0 -y1*y3**2*1d0/2.D0 &
            +y1*y2**2*1d0/2.D0 -y1**2*y3 )*FEA1125 
 
         dv4(5)=ds2(5) +0d0 +0d0 +0d0 &
            +( 5.D0/8.D0*y2*y4*2d0*y5 &
            +SQRT3*y2*3d0*y5**2/8.D0 -SQRT3*y3*y4**2*1d0/8.D0 &
            +SQRT3*y2*y4**2*1d0/8.D0 +y1*y4*2d0*y5 -SQRT3*y3*3d0*y5**2/8.D0 &
            +5.D0/8.D0*y3*y4*2d0*y5 )*FEA1455 
 
         ds3(5)=( -2.D0*y4**3*2d0*y5 &
            -3.D0*y4*4d0*y5**3 )*FEA44444 &
            +( -4.D0*y3*y4*3d0*y5**2*SQRT3 +9.D0*y1*y4**2*2d0*y5 &
            +4.D0*y2*y4*3d0*y5**2*SQRT3 +5.D0/2.D0*y1*4d0*y5**3 +y2*4d0*y5**3 &
            +y3*4d0*y5**3 )*FEA25555 +( &
            +y3*y4**2*2d0*y5 -2.D0*y2*y4*3d0*y5**2*SQRT3 &
            -7.D0/2.D0*y1*y4**2*2d0*y5 -3.D0/4.D0*y1*4d0*y5**3 &
            +2.D0*y3*y4*3d0*y5**2*SQRT3 +y2*y4**2*2d0*y5 )*FEA24455 
 
         ds2(5)=ds3(5) &
            +( y2*y4**3 -3.D0*y3*y4*3d0*y5**2 +3.D0/4.D0*y1*4d0*y5**3*SQRT3 &
            +3.D0*y2*y4*3d0*y5**2 &
            +3.D0/2.D0*y1*y4**2*2d0*y5*SQRT3 -y3*y4**3 )*FEA24445 &
            +( -y2**2*3d0*y5**2 +y3**2*y4**2 &
            +y3**2*3d0*y5**2 -y2**2*y4**2 -y1**2*y4*2d0*y5*SQRT3 )*FEA33445 &
            +( y3**2*y4*2d0*y5 +y1**2*y4*2d0*y5 &
            +y2**2*y4*2d0*y5 )*FEA33455 
 
         ds1(5)=ds2(5) +( -y2**3*y4 +y3**3*y4 &
            +y2**3*2d0*y5*SQRT3/3.D0 &
            +y3**3*2d0*y5*SQRT3/3.D0 -y1**3*2d0*y5*SQRT3/6.D0 )*FEA33345 +( &
            +y3**3*2d0*y5 +y2**3*2d0*y5 +y1**3*2d0*y5 )*FEA33344 +( &
            +SQRT3*y3**4 -SQRT3*y2**4 )*FEA33334 +0d0 +( -y1*y2*3d0*y5**2 &
            +y1*y3*y4**2 +y2*y3*y4*2d0*y5*SQRT3 -y1*y2*y4**2 &
            +y1*y3*3d0*y5**2 )*FEA13445 +( y2*y3*y4*2d0*y5 +y1*y2*y4*2d0*y5 &
            +y1*y3*y4*2d0*y5 )*FEA13455 
 
         ds3(5)=ds1(5) +( y1**2*y3*2d0*y5 &
            +y2**2*y3*2d0*y5 +y1*y2**2*2d0*y5 +y1**2*y2*2d0*y5 +y1*y3**2*2d0*y5 &
            +y2*y3**2*2d0*y5 )*FEA11255 +( +y1*y3**2*2d0*y5*SQRT3/2.D0 &
            +y1*y2**2*2d0*y5*SQRT3/2.D0 +y2**2*y3*2d0*y5*SQRT3/2.D0 -y1*y2**2*y4 &
            +y2*y3**2*y4 +y1*y3**2*y4 -y2**2*y3*y4 &
            +y2*y3**2*2d0*y5*SQRT3/2.D0 )*FEA13345 
 
         ds4(5)=ds3(5) +( y1**2*y2*y4 &
            +y2*y3**2*y4 -y2**2*y3*y4 -y1**2*y3*y4 +y1*y2**2*2d0*y5*SQRT3/2.D0 &
            +y1*y3**2*2d0*y5*SQRT3/2.D0 )*FEA11245 
 
         ds2(5)=ds4(5) +( -y1**3*y2 &
            +y1**3*y3 +y2**3*y3*1d0/2.D0 -y1*y2**3*1d0/2.D0 -y2*y3**3*1d0/2.D0 &
            +y1*y3**3*1d0/2.D0 )*FEA11135 +( +y1*y2**3*SQRT3/2.D0 &
            +y2**3*y3*SQRT3/2.D0 -y2*y3**3*SQRT3/2.D0 &
            -y1*y3**3*SQRT3/2.D0 )*FEA11134 
 
         dv5(5)=ds2(5) &
            +0d0 +( -SQRT3*y1**2*y3**2 +SQRT3*y1**2*y2**2 )*FEA11334 +0d0 +( &
            +y1*y2*y3*2d0*y5 )*FEA12355 +( -SQRT3*y1*y2*y3**2*1d0/2.D0 &
            +SQRT3*y1*y2**2*y3*1d0/2.D0 )*FEA11234 +0d0 +0d0 
 
         ds3(5)=( -y2**3*y4**2 &
            +y3**3*y4**2 -5.D0/3.D0*y2**3*y4*2d0*y5*SQRT3 &
            -5.D0/3.D0*y3**3*y4*2d0*y5*SQRT3 -y2**3*3d0*y5**2 &
            +y3**3*3d0*y5**2 -8.D0/3.D0*y1**3*y4*2d0*y5*SQRT3 )*FEA333555 &
            +( y1**4*2d0*y5*SQRT3/2.D0 +y2**4*y4 -y3**4*y4 )*FEA222245 +0d0 +( &
            +y1*y2**4*SQRT3 &
            +y1**4*y2*SQRT3 -y1*y3**4*SQRT3 -y1**4*y3*SQRT3 )*FEA133334 +( &
            +y1*y2*y3*y4*2d0*y5 )*FEA123455 
 
         ds4(5)=ds3(5) +( -y1**2*y2**2*y3 &
            +y1**2*y2*y3**2 )*FEA112335 +( y1*y2**2*y3*2d0*y5 +y1*y2*y3**2*2d0*y5 &
            +y1**2*y2*y3*2d0*y5 )*FEA112355 
 
         ds2(5)=ds4(5) &
            +( y2**3*y3**2 -y1**3*y2**2*1d0/2.D0 -y1**2*y3**3*1d0/2.D0 -y2**2*y3**3 &
            +y1**3*y3**2*1d0/2.D0 +y1**2*y2**3*1d0/2.D0 )*FEA222335 &
            +( -y1**2*y2**2*2d0*y5*SQRT3/2.D0 -y1**2*y3**2*2d0*y5*SQRT3/2.D0 &
            -y1**2*y2**2*y4 &
            +y1**2*y3**2*y4 )*FEA113345 +( y2**2*y3**2*2d0*y5 +y1**2*y2**2*2d0*y5 &
            +y1**2*y3**2*2d0*y5 )*FEA223355 
 
         ds3(5)=ds2(5) +( +y1*y2*y3**2*y4 &
            +y1*y2*y3**2*2d0*y5*SQRT3/2.D0 -y1*y2**2*y3*y4 &
            +y1*y2**2*y3*2d0*y5*SQRT3/2.D0 )*FEA123345 +( -y1**3*y2**2*SQRT3/2.D0 &
            +y1**3*y3**2*SQRT3/2.D0 -y1**2*y2**3*SQRT3/2.D0 &
            +y1**2*y3**3*SQRT3/2.D0 )*FEA222334 +( +5.D0/2.D0*y1**2*4d0*y5**3 &
            +y2**2*4d0*y5**3 -4.D0*y3**2*y4*3d0*y5**2*SQRT3 +y3**2*4d0*y5**3 &
            +9.D0*y1**2*y4**2*2d0*y5 +4.D0*y2**2*y4*3d0*y5**2*SQRT3 )*FEA335555 &
            +0d0 
 
         ds4(5)=ds3(5) &
            +( -y2*y4**4*SQRT3/2.D0 -2.D0*y1*y4**3*2d0*y5 &
            -3.D0/10.D0*y2*5d0*y5**4*SQRT3 &
            +y3*y4**3*2d0*y5 +y3*y4**4*SQRT3/2.D0 +y2*y4**3*2d0*y5 &
            +3.D0/10.D0*y3*5d0*y5**4*SQRT3 )*FEA244455 +( -SQRT3*y2**5 &
            +SQRT3*y3**5 )*FEA222224 
 
         ds5(5)=ds4(5) +( -y3*5d0*y5**4*SQRT3/5.D0 &
            +y2*5d0*y5**4*SQRT3/5.D0 +y1*y4*4d0*y5**3 &
            +y2*y4**4*SQRT3/3.D0 -y3*y4**4*SQRT3/3.D0 +y3*y4*4d0*y5**3 &
            +y2*y4*4d0*y5**3 +2.D0*y1*y4**3*2d0*y5 )*FEA145555 
 
         ds1(5)=ds5(5) &
            +( -SQRT3*y1*y2*y3**3*1d0/2.D0 +SQRT3*y1*y2**3*y3*1d0/2.D0 )*FEA111234 &
            +( y3*y4**4*1d0/3.D0 -y2*y4**4*1d0/3.D0 -y2*y4*4d0*y5**3*SQRT3/2.D0 &
            -y3*y4**2*3d0*y5**2 &
            +y2*y4**2*3d0*y5**2 -2.D0/3.D0*y1*y4**3*2d0*y5*SQRT3 &
            -y3*y4*4d0*y5**3*SQRT3/2.D0 )*FEA244555 &
            +( y1*y2*y4**2*2d0*y5 -y1*y2*4d0*y5**3 -y1*y3*4d0*y5**3 &
            +5.D0/4.D0*y2*y3*4d0*y5**3 &
            +y1*y3*y4**2*2d0*y5 -7.D0/2.D0*y2*y3*y4**2*2d0*y5 &
            -2.D0*y1*y2*y4**3*SQRT3 &
            +2.D0*y1*y3*y4**3*SQRT3 )*FEA124455 
 
         ds3(5)=ds1(5) +0d0 +0d0 +( &
            +y2**4*2d0*y5 +y1**4*2d0*y5 +y3**4*2d0*y5 )*FEA222255 
 
         ds4(5)=ds3(5) &
            +( 3.D0*y1*y3*4d0*y5**3 &
            +9.D0*y2*y3*y4**2*2d0*y5 -3.D0/2.D0*y2*y3*4d0*y5**3 &
            -4.D0*y1*y3*y4**3*SQRT3 &
            +4.D0*y1*y2*y4**3*SQRT3 +3.D0*y1*y2*4d0*y5**3 )*FEA134444 &
            +( -y1*y3**2*3d0*y5**2*SQRT3/3.D0 -7.D0/3.D0*y1**2*y3*y4*2d0*y5 &
            +5.D0/3.D0*y1*y2**2*y4**2*SQRT3 -13.D0/3.D0*y2**2*y3*y4*2d0*y5 &
            -4.D0/3.D0*y2*y3**2*3d0*y5**2*SQRT3 -7.D0/3.D0*y1**2*y2*y4*2d0*y5 &
            -16.D0/3.D0*y1*y3**2*y4*2d0*y5 &
            +4.D0/3.D0*y1**2*y3*y4**2*SQRT3 +4.D0/3.D0*y2**2*y3*3d0*y5**2*SQRT3 &
            +y1*y2**2*3d0*y5**2*SQRT3/3.D0 -13.D0/3.D0*y2*y3**2*y4*2d0*y5 &
            -5.D0/3.D0*y1*y3**2*y4**2*SQRT3 -4.D0/3.D0*y1**2*y2*y4**2*SQRT3 &
            -16.D0/3.D0*y1*y2**2*y4*2d0*y5 )*FEA233444 
 
         ds5(5)=ds4(5) &
            +( 2.D0*y1*y3**2*3d0*y5**2 +4.D0*y2*y3**2*3d0*y5**2 &
            +4.D0*y2**2*y3*y4*2d0*y5*SQRT3 -2.D0*y1*y2**2*3d0*y5**2 &
            +y1**2*y3*y4*2d0*y5*SQRT3 &
            +6.D0*y1*y3**2*y4**2 -6.D0*y1*y2**2*y4**2 -3.D0*y1**2*y3*y4**2 &
            +y1**2*y2*y4*2d0*y5*SQRT3 &
            +4.D0*y1*y3**2*y4*2d0*y5*SQRT3 -4.D0*y2**2*y3*3d0*y5**2 &
            +3.D0*y1**2*y2*y4**2 -y1**2*y2*3d0*y5**2 +y1**2*y3*3d0*y5**2 &
            +4.D0*y2*y3**2*y4*2d0*y5*SQRT3 &
            +4.D0*y1*y2**2*y4*2d0*y5*SQRT3 )*FEA113555 
 
         ds2(5)=ds5(5) &
            +( -3.D0/2.D0*y1**2*y4**2*2d0*y5*SQRT3 &
            -3.D0/4.D0*y1**2*4d0*y5**3*SQRT3 -y2**2*y4**3 &
            +y3**2*y4**3 &
            +3.D0*y3**2*y4*3d0*y5**2 -3.D0*y2**2*y4*3d0*y5**2 )*FEA334445 &
            +( -3.D0*y1*y3*y4**3 &
            +2.D0/3.D0*y1*y2*4d0*y5**3*SQRT3 -y1*y3*y4*3d0*y5**2 &
            +2.D0/3.D0*y1*y3*4d0*y5**3*SQRT3 &
            +3.D0*y1*y2*y4**3 -7.D0/12.D0*y2*y3*4d0*y5**3*SQRT3 &
            +3.D0/2.D0*y2*y3*y4**2*2d0*y5*SQRT3 +y1*y2*y4*3d0*y5**2 )*FEA124555 &
            +( 2.D0*y3**2*y4*3d0*y5**2*SQRT3 -7.D0/2.D0*y1**2*y4**2*2d0*y5 &
            +y2**2*y4**2*2d0*y5 -2.D0*y2**2*y4*3d0*y5**2*SQRT3 &
            -3.D0/4.D0*y1**2*4d0*y5**3 &
            +y3**2*y4**2*2d0*y5 )*FEA334455 
 
         ds3(5)=ds2(5) +( -6.D0*y4**2*4d0*y5**3 &
            +9.D0*y4**4*2d0*y5 +6d0*y5**5 )*FEA555555 +( +y2*y3**3*2d0*y5 &
            +y1*y2**3*2d0*y5 +y1**3*y3*2d0*y5 +y1**3*y2*2d0*y5 +y1*y3**3*2d0*y5 &
            +y2**3*y3*2d0*y5 )*FEA233344 &
            +( y1*y2**3*2d0*y5*SQRT3/6.D0 -y2**3*y3*2d0*y5*SQRT3/3.D0 &
            -y2*y3**3*2d0*y5*SQRT3/3.D0 &
            +y1**3*y2*y4 -y1**3*y2*2d0*y5*SQRT3/3.D0 -y1**3*y3*y4 &
            -y1**3*y3*2d0*y5*SQRT3/3.D0 &
            +y1*y3**3*2d0*y5*SQRT3/6.D0 -y2**3*y3*y4 +y2*y3**3*y4 )*FEA233345 &
            +( -3.D0*y2**3*y4*2d0*y5 -3.D0*y3**3*y4*2d0*y5 &
            -3.D0*y1**3*y4*2d0*y5 )*FEA111444 &
            +0d0 
 
         ds4(5)=ds3(5) &
            +( 9.D0*y4**2*4d0*y5**3 -6.D0*y4**4*2d0*y5 )*FEA444444 &
            +( -5.D0/3.D0*y1*y2**2*y4**2*SQRT3 -4.D0/3.D0*y1**2*y3*y4**2*SQRT3 &
            -y1*y2**2*3d0*y5**2*SQRT3/3.D0 &
            +4.D0/3.D0*y2**2*y3*y4*2d0*y5 -4.D0/3.D0*y2**2*y3*3d0*y5**2*SQRT3 &
            +7.D0/3.D0*y1*y2**2*y4*2d0*y5 -2.D0/3.D0*y1**2*y3*y4*2d0*y5 &
            +4.D0/3.D0*y2*y3**2*3d0*y5**2*SQRT3 +y1*y3**2*3d0*y5**2*SQRT3/3.D0 &
            +4.D0/3.D0*y1**2*y2*y4**2*SQRT3 +4.D0/3.D0*y2*y3**2*y4*2d0*y5 &
            +5.D0/3.D0*y1*y3**2*y4**2*SQRT3 -2.D0/3.D0*y1**2*y2*y4*2d0*y5 &
            +7.D0/3.D0*y1*y3**2*y4*2d0*y5 )*FEA133444 
 
         ds5(5)=ds4(5) &
            +( -y1**3*y2*y4 +2.D0/3.D0*y2**3*y3*2d0*y5*SQRT3 &
            +y1**3*y3*2d0*y5*SQRT3/6.D0 +y1**3*y2*2d0*y5*SQRT3/6.D0 +y1**3*y3*y4 &
            +y1*y2**3*2d0*y5*SQRT3/6.D0 &
            +2.D0/3.D0*y2*y3**3*2d0*y5*SQRT3 -y1*y2**3*y4 &
            +y1*y3**3*2d0*y5*SQRT3/6.D0 +y1*y3**3*y4 )*FEA133345 
 
         dv6(5)=ds5(5) &
            +( -y2**2*y3*y4**2 +y1**2*y3*y4*2d0*y5*SQRT3/3.D0 +y2*y3**2*y4**2 &
            +y2*y3**2*3d0*y5**2 -y1*y2**2*3d0*y5**2 &
            +4.D0/3.D0*y2**2*y3*y4*2d0*y5*SQRT3 &
            +4.D0/3.D0*y2*y3**2*y4*2d0*y5*SQRT3 -y1*y2**2*y4**2 &
            +4.D0/3.D0*y1*y3**2*y4*2d0*y5*SQRT3 -y2**2*y3*3d0*y5**2 &
            +y1*y3**2*3d0*y5**2 +y1**2*y2*y4*2d0*y5*SQRT3/3.D0 +y1*y3**2*y4**2 &
            +4.D0/3.D0*y1*y2**2*y4*2d0*y5*SQRT3 )*FEA233445 +( -y1**4*y2 &
            +y2*y3**4 -2.D0*y1*y2**4 +2.D0*y1*y3**4 &
            +y1**4*y3 -y2**4*y3 )*FEA233335 +0d0 
 
         dv(5)= +dv1(5) +dv2(5) &
            +dv3(5) +dv4(5) +dv5(5) +dv6(5)
 

         dv1(6)=( y3 +y2 +y1 )*DFEA1  
 
         dv2(6)=( y2*y3 +y1*y3 +y1*y2 )*DFEA12 &
            +( y2**2 +y3**2 +y1**2 )*DFEA11 +(  -SQRT3*y3*y5/2.D0 -y3*y4/2.D0 &
            +y1*y4 +SQRT3*y2*y5/2.D0 -y2*y4/2.D0 )*DFEA14 +( y5**2 &
            +y4**2 )*DFEA44  
 
         dv3(6)=( y1*y3*y4 +y1*y2*y4 -2.D0*y2*y3*y4 &
            +SQRT3*y1*y2*y5 -SQRT3*y1*y3*y5 )*DFEA124 &
            +( 3.D0/4.D0*y3*y4**2 -SQRT3*y3*y4*y5/2.D0 +y1*y5**2 +y2*y5**2/4.D0 &
            +3.D0/4.D0*y2*y4**2 +SQRT3*y2*y4*y5/2.D0 +y3*y5**2/4.D0 )*DFEA155 &
            +( y2*y3**2 +y1*y3**2 +y1**2*y3 +y1*y2**2 +y2**2*y3 &
            +y1**2*y2 )*DFEA112 +(  -y4**3/3.D0 +y4*y5**2 )*DFEA455 &
            +y1*y2*y3*DFEA123 +( y1*y4**2 +3.D0/4.D0*y3*y5**2 +3.D0/4.D0*y2*y5**2 &
            +y2*y4**2/4.D0 -SQRT3*y2*y4*y5/2.D0 +SQRT3*y3*y4*y5/2.D0 &
            +y3*y4**2/4.D0 )*DFEA144 +( y3**3 +y2**3 +y1**3 )*DFEA111 &
            +(  -y2**2*y4/2.D0 -y3**2*y4/2.D0 +SQRT3*y2**2*y5/2.D0 &
            +y1**2*y4 -SQRT3*y3**2*y5/2.D0 )*DFEA114   
 
         ds2(6)=( y4**4 +y5**4 &
            +2.D0*y4**2*y5**2 )*DFEA4444 &
            +( 3.D0/8.D0*SQRT3*y2*y5**3 -3.D0/8.D0*SQRT3*y3*y4**2*y5 &
            -3.D0/8.D0*SQRT3*y3*y5**3 -9.D0/8.D0*y2*y4*y5**2 -y3*y4**3/8.D0 &
            -y2*y4**3/8.D0 -9.D0/8.D0*y3*y4*y5**2 &
            +y1*y4**3 +3.D0/8.D0*SQRT3*y2*y4**2*y5 )*DFEA1444 &
            +( 3.D0/4.D0*y2**2*y4**2 +3.D0/4.D0*y3**2*y4**2 +y1**2*y5**2 &
            +y3**2*y5**2/4.D0 -SQRT3*y3**2*y4*y5/2.D0 +SQRT3*y2**2*y4*y5/2.D0 &
            +y2**2*y5**2/4.D0 )*DFEA1155 
 
         ds1(6)=ds2(6) +( y3**2*y4**2/4.D0 &
            +3.D0/4.D0*y3**2*y5**2 +y1**2*y4**2 +y2**2*y4**2/4.D0 &
            +SQRT3*y3**2*y4*y5/2.D0 -SQRT3*y2**2*y4*y5/2.D0 &
            +3.D0/4.D0*y2**2*y5**2 )*DFEA1144 +( y1**3*y4 &
            +SQRT3*y2**3*y5/2.D0 -SQRT3*y3**3*y5/2.D0 -y2**3*y4/2.D0 &
            -y3**3*y4/2.D0 )*DFEA1114 &
            +( y2**4 +y1**4 +y3**4 )*DFEA1111 +( SQRT3*y1*y3*y4*y5 &
            +3.D0/2.D0*y2*y3*y5**2 -y2*y3*y4**2/2.D0 &
            +y1*y2*y4**2 -SQRT3*y1*y2*y4*y5 +y1*y3*y4**2 )*DFEA1244 
 
         ds2(6)=ds1(6) &
            +( y1*y3*y5**2 +y1*y2*y5**2 -SQRT3*y1*y3*y4*y5 -y2*y3*y5**2/2.D0 &
            +3.D0/2.D0*y2*y3*y4**2 +SQRT3*y1*y2*y4*y5 )*DFEA1255 &
            +(  -y1*y3**2*y4/2.D0 &
            +y1**2*y3*y4 -SQRT3*y1*y3**2*y5/2.D0 -SQRT3*y2*y3**2*y5/2.D0 &
            +y1**2*y2*y4 +SQRT3*y2**2*y3*y5/2.D0 -y2**2*y3*y4/2.D0 &
            +SQRT3*y1*y2**2*y5/2.D0 -y2*y3**2*y4/2.D0 -y1*y2**2*y4/2.D0 )*DFEA1124 &
            +( y1**2*y2*y5 +SQRT3*y1*y3**2*y4/2.D0 &
            +SQRT3*y1*y2**2*y4/2.D0 -SQRT3*y2*y3**2*y4/2.D0 &
            -SQRT3*y2**2*y3*y4/2.D0 -y2**2*y3*y5/2.D0 &
            +y2*y3**2*y5/2.D0 -y1*y3**2*y5/2.D0 &
            +y1*y2**2*y5/2.D0 -y1**2*y3*y5 )*DFEA1125 
 
         dv4(6)=ds2(6) +( y2*y3**3 &
            +y1**3*y3 +y1**3*y2 +y1*y2**3 +y1*y3**3 +y2**3*y3 )*DFEA1112 &
            +( y2**2*y3**2 +y1**2*y3**2 +y1**2*y2**2 )*DFEA1122 +( y1*y2**2*y3 &
            +y1**2*y2*y3 +y1*y2*y3**2 )*DFEA1123 +( 5.D0/8.D0*y2*y4*y5**2 &
            +SQRT3*y2*y5**3/8.D0 -SQRT3*y3*y4**2*y5/8.D0 &
            +SQRT3*y2*y4**2*y5/8.D0 -3.D0/8.D0*y2*y4**3 &
            +y1*y4*y5**2 -SQRT3*y3*y5**3/8.D0 &
            +5.D0/8.D0*y3*y4*y5**2 -3.D0/8.D0*y3*y4**3 )*DFEA1455   
 
         ds3(6)=( y4**5 &
            -2.D0*y4**3*y5**2 -3.D0*y4*y5**4 )*DFEA44444 &
            +(  -4.D0*y3*y4*y5**3*SQRT3 +9.D0*y1*y4**2*y5**2 -3.D0/2.D0*y1*y4**4 &
            +4.D0*y2*y4*y5**3*SQRT3 +3.D0*y2*y4**4 +5.D0/2.D0*y1*y5**4 &
            +3.D0*y3*y4**4 +y2*y5**4 +y3*y5**4 )*DFEA25555 +(  -y2*y4**4 &
            +y3*y4**2*y5**2 -2.D0*y2*y4*y5**3*SQRT3 -y3*y4**4 &
            -7.D0/2.D0*y1*y4**2*y5**2 -3.D0/4.D0*y1*y5**4 &
            +2.D0*y3*y4*y5**3*SQRT3 +y2*y4**2*y5**2 &
            +5.D0/4.D0*y1*y4**4 )*DFEA24455 
 
         ds2(6)=ds3(6) &
            +( y2*y4**3*y5 -3.D0*y3*y4*y5**3 +2.D0/3.D0*y3*y4**4*SQRT3 &
            +3.D0/4.D0*y1*y5**4*SQRT3 +3.D0*y2*y4*y5**3 -7.D0/12.D0*y1*y4**4*SQRT3 &
            +3.D0/2.D0*y1*y4**2*y5**2*SQRT3 -y3*y4**3*y5 &
            +2.D0/3.D0*y2*y4**4*SQRT3 )*DFEA24445 +(  -y2**2*y5**3 +y3**2*y4**2*y5 &
            +y3**2*y5**3 +4.D0/9.D0*y2**2*y4**3*SQRT3 -5.D0/9.D0*y1**2*y4**3*SQRT3 &
            +4.D0/9.D0*y3**2*y4**3*SQRT3 -y2**2*y4**2*y5 &
            -y1**2*y4*y5**2*SQRT3 )*DFEA33445 &
            +( y3**2*y4*y5**2 -y1**2*y4**3/3.D0 -y3**2*y4**3/3.D0 +y1**2*y4*y5**2 &
            +y2**2*y4*y5**2 -y2**2*y4**3/3.D0 )*DFEA33455 
 
         ds1(6)=ds2(6) &
            +(  -y2**3*y4*y5 +y3**3*y4*y5 +y2**3*y5**2*SQRT3/3.D0 &
            +y1**3*y4**2*SQRT3/2.D0 &
            +y3**3*y5**2*SQRT3/3.D0 -y1**3*y5**2*SQRT3/6.D0 )*DFEA33345 &
            +( y3**3*y4**2 +y3**3*y5**2 +y2**3*y4**2 +y2**3*y5**2 +y1**3*y5**2 &
            +y1**3*y4**2 )*DFEA33344 +( y3**4*y4 +SQRT3*y3**4*y5 &
            +y2**4*y4 -2.D0*y1**4*y4 -SQRT3*y2**4*y5 )*DFEA33334 +( y2**5 +y3**5 &
            +y1**5 )*DFEA33333 +(  -4.D0/9.D0*y1*y2*y4**3*SQRT3 -y1*y2*y5**3 &
            +y1*y3*y4**2*y5 +y2*y3*y4*y5**2*SQRT3 -y1*y2*y4**2*y5 &
            +5.D0/9.D0*y2*y3*y4**3*SQRT3 -4.D0/9.D0*y1*y3*y4**3*SQRT3 &
            +y1*y3*y5**3 )*DFEA13445 +( y2*y3*y4*y5**2 &
            +y1*y2*y4*y5**2 -y2*y3*y4**3/3.D0 -y1*y2*y4**3/3.D0 -y1*y3*y4**3/3.D0 &
            +y1*y3*y4*y5**2 )*DFEA13455 
 
         ds3(6)=ds1(6) +( y1**2*y3*y5**2 &
            +y2**2*y3*y4**2 +y2**2*y3*y5**2 +y1*y2**2*y5**2 +y1**2*y2*y5**2 &
            +y1*y2**2*y4**2 +y2*y3**2*y4**2 +y1*y3**2*y4**2 +y1**2*y3*y4**2 &
            +y1**2*y2*y4**2 +y1*y3**2*y5**2 +y2*y3**2*y5**2 )*DFEA11255 &
            +( 2.D0/3.D0*y1**2*y3*y4**2*SQRT3 +y1*y3**2*y5**2*SQRT3/2.D0 &
            +y1*y2**2*y5**2*SQRT3/2.D0 +y2**2*y3*y5**2*SQRT3/2.D0 -y1*y2**2*y4*y5 &
            +y2*y3**2*y4*y5 +y1*y3**2*y4*y5 -y2**2*y3*y4*y5 &
            +y2*y3**2*y4**2*SQRT3/6.D0 +y1*y3**2*y4**2*SQRT3/6.D0 &
            +y1*y2**2*y4**2*SQRT3/6.D0 +2.D0/3.D0*y1**2*y2*y4**2*SQRT3 &
            +y2*y3**2*y5**2*SQRT3/2.D0 &
            +y2**2*y3*y4**2*SQRT3/6.D0 )*DFEA13345 
 
         ds4(6)=ds3(6) +( y1**2*y2*y4*y5 &
            +y1**2*y3*y4**2*SQRT3/3.D0 &
            +y1**2*y2*y4**2*SQRT3/3.D0 -y1*y2**2*y4**2*SQRT3/6.D0 &
            +y2*y3**2*y4*y5 -y2**2*y3*y4*y5 -y1**2*y3*y4*y5 &
            +y2*y3**2*y4**2*SQRT3/3.D0 &
            +y1*y2**2*y5**2*SQRT3/2.D0 -y1*y3**2*y4**2*SQRT3/6.D0 &
            +y2**2*y3*y4**2*SQRT3/3.D0 &
            +y1*y3**2*y5**2*SQRT3/2.D0 )*DFEA11245 
 
         ds2(6)=ds4(6) +(  -y1**3*y2*y5 &
            +y1**3*y3*y5 &
            +y2**3*y3*y5/2.D0 -y1*y2**3*y4*SQRT3/2.D0 -y1*y2**3*y5/2.D0 &
            -y2*y3**3*y5/2.D0 &
            +y1*y3**3*y5/2.D0 +y2**3*y3*y4*SQRT3/2.D0 &
            +y2*y3**3*y4*SQRT3/2.D0 -y1*y3**3*y4*SQRT3/2.D0 )*DFEA11135 &
            +( y1**3*y3*y4 -y2**3*y3*y4/2.D0 &
            +y1**3*y2*y4 -y2*y3**3*y4/2.D0 -y1*y3**3*y4/2.D0 &
            +y1*y2**3*y5*SQRT3/2.D0 &
            +y2**3*y3*y5*SQRT3/2.D0 -y2*y3**3*y5*SQRT3/2.D0 -y1*y2**3*y4/2.D0 &
            -y1*y3**3*y5*SQRT3/2.D0 )*DFEA11134 
 
         dv5(6)=ds2(6) &
            +( y1*y2**4 +y1**4*y3 +y1**4*y2 +y2**4*y3 +y2*y3**4 &
            +y1*y3**4 )*DFEA23333 +(  -2.D0*y2**2*y3**2*y4 &
            +y1**2*y2**2*y4 -SQRT3*y1**2*y3**2*y5 +SQRT3*y1**2*y2**2*y5 &
            +y1**2*y3**2*y4 )*DFEA11334 +( y1**2*y3**3 +y1**3*y3**2 +y2**2*y3**3 &
            +y1**2*y2**3 +y1**3*y2**2 +y2**3*y3**2 )*DFEA11333 +( y1*y2*y3*y4**2 &
            +y1*y2*y3*y5**2 )*DFEA12355 &
            +(  -y1*y2*y3**2*y4/2.D0 -y1*y2**2*y3*y4/2.D0 &
            -SQRT3*y1*y2*y3**2*y5/2.D0 &
            +y1**2*y2*y3*y4 +SQRT3*y1*y2**2*y3*y5/2.D0 )*DFEA11234 +( y1*y2**3*y3 &
            +y1*y2*y3**3 +y1**3*y2*y3 )*DFEA11123 +( y1**2*y2**2*y3 &
            +y1*y2**2*y3**2 &
            +y1**2*y2*y3**2 )*DFEA11233  
 
         ds3(6)=( y2**3*y4**3*SQRT3 -y2**3*y4**2*y5 &
            +y3**3*y4**2*y5 -5.D0/3.D0*y2**3*y4*y5**2*SQRT3 &
            +y3**3*y4**3*SQRT3 -5.D0/3.D0*y3**3*y4*y5**2*SQRT3 -y2**3*y5**3 &
            +y3**3*y5**3 -8.D0/3.D0*y1**3*y4*y5**2*SQRT3 )*DFEA333555 &
            +( y1**4*y5**2*SQRT3/2.D0 +y2**4*y4*y5 +y2**4*y4**2*SQRT3/3.D0 &
            +y3**4*y4**2*SQRT3/3.D0 -y3**4*y4*y5 &
            -y1**4*y4**2*SQRT3/6.D0 )*DFEA222245 &
            +( y1*y3**5 +y1*y2**5 +y2**5*y3 +y1**5*y3 +y1**5*y2 &
            +y2*y3**5 )*DFEA133333 +( y1**4*y3*y4 -2.D0*y2**4*y3*y4 +y1**4*y2*y4 &
            +y1*y2**4*y5*SQRT3 +y1*y3**4*y4 -2.D0*y2*y3**4*y4 &
            +y1**4*y2*y5*SQRT3 -y1*y3**4*y5*SQRT3 -y1**4*y3*y5*SQRT3 &
            +y1*y2**4*y4 )*DFEA133334 +(  -y1*y2*y3*y4**3/3.D0 &
            +y1*y2*y3*y4*y5**2 )*DFEA123455 
 
         ds4(6)=ds3(6) &
            +( 2.D0/3.D0*SQRT3*y1*y2**2*y3**2*y4 -y1**2*y2**2*y3*y5 &
            -SQRT3*y1**2*y2**2*y3*y4/3.D0 &
            +y1**2*y2*y3**2*y5 -SQRT3*y1**2*y2*y3**2*y4/3.D0 )*DFEA112335 &
            +( y1*y2**2*y3*y5**2 +y1*y2*y3**2*y5**2 +y1*y2*y3**2*y4**2 &
            +y1*y2**2*y3*y4**2 +y1**2*y2*y3*y4**2 &
            +y1**2*y2*y3*y5**2 )*DFEA112355 
 
         ds2(6)=ds4(6) &
            +( y2**3*y3**2*y5 -y1**3*y2**2*y5/2.D0 -y1**2*y3**3*y5/2.D0 &
            -y2**2*y3**3*y5 &
            +y1**3*y2**2*y4*SQRT3/2.D0 -y1**2*y2**3*y4*SQRT3/2.D0 &
            +y1**3*y3**2*y5/2.D0 +y1**2*y2**3*y5/2.D0 &
            +y1**3*y3**2*y4*SQRT3/2.D0 -y1**2*y3**3*y4*SQRT3/2.D0 )*DFEA222335 &
            +(  -y1**2*y2**2*y5**2*SQRT3/2.D0 -y1**2*y3**2*y5**2*SQRT3/2.D0 &
            -y1**2*y2**2*y4**2*SQRT3/6.D0 -y1**2*y2**2*y4*y5 &
            -2.D0/3.D0*y2**2*y3**2*y4**2*SQRT3 &
            +y1**2*y3**2*y4*y5 -y1**2*y3**2*y4**2*SQRT3/6.D0 )*DFEA113345 &
            +( y2**2*y3**2*y5**2 +y2**2*y3**2*y4**2 +y1**2*y2**2*y5**2 &
            +y1**2*y3**2*y4**2 +y1**2*y3**2*y5**2 &
            +y1**2*y2**2*y4**2 )*DFEA223355 
 
         ds3(6)=ds2(6) &
            +( y1*y2*y3**2*y4**2*SQRT3/6.D0 +y1*y2*y3**2*y4*y5 &
            +y1*y2*y3**2*y5**2*SQRT3/2.D0 &
            +2.D0/3.D0*y1**2*y2*y3*y4**2*SQRT3 -y1*y2**2*y3*y4*y5 &
            +y1*y2**2*y3*y4**2*SQRT3/6.D0 &
            +y1*y2**2*y3*y5**2*SQRT3/2.D0 )*DFEA123345 &
            +(  -y1**3*y2**2*y5*SQRT3/2.D0 -y1**3*y2**2*y4/2.D0 &
            -y1**3*y3**2*y4/2.D0 -y1**2*y2**3*y4/2.D0 &
            +y1**3*y3**2*y5*SQRT3/2.D0 -y1**2*y3**3*y4/2.D0 &
            +y2**3*y3**2*y4 -y1**2*y2**3*y5*SQRT3/2.D0 +y2**2*y3**3*y4 &
            +y1**2*y3**3*y5*SQRT3/2.D0 )*DFEA222334 +( 3.D0*y3**2*y4**4 &
            +5.D0/2.D0*y1**2*y5**4 +y2**2*y5**4 &
            +3.D0*y2**2*y4**4 -4.D0*y3**2*y4*y5**3*SQRT3 +y3**2*y5**4 &
            +9.D0*y1**2*y4**2*y5**2 -3.D0/2.D0*y1**2*y4**4 &
            +4.D0*y2**2*y4*y5**3*SQRT3 )*DFEA335555 +( y1**3*y2**3 +y1**3*y3**3 &
            +y2**3*y3**3 )*DFEA222333 
 
         ds4(6)=ds3(6) &
            +( y3*y4**5/5.D0 -y2*y4**4*y5*SQRT3/2.D0 -2.D0/5.D0*y1*y4**5 &
            -2.D0*y1*y4**3*y5**2 -3.D0/10.D0*y2*y5**5*SQRT3 &
            +y3*y4**3*y5**2 +y3*y4**4*y5*SQRT3/2.D0 +y2*y4**3*y5**2 &
            +3.D0/10.D0*y3*y5**5*SQRT3 +y2*y4**5/5.D0 )*DFEA244455 &
            +( y2**5*y4 -2.D0*y1**5*y4 -SQRT3*y2**5*y5 +y3**5*y4 &
            +SQRT3*y3**5*y5 )*DFEA222224 
 
         ds5(6)=ds4(6) +(  -y3*y5**5*SQRT3/5.D0 &
            +y2*y5**5*SQRT3/5.D0 +y1*y4*y5**4 -7.D0/15.D0*y2*y4**5 &
            +y2*y4**4*y5*SQRT3/3.D0 -y3*y4**4*y5*SQRT3/3.D0 +y3*y4*y5**4 &
            +y2*y4*y5**4 &
            +2.D0*y1*y4**3*y5**2 -7.D0/15.D0*y3*y4**5 &
            -y1*y4**5/15.D0 )*DFEA145555 
 
         ds1(6)=ds5(6) &
            +(  -SQRT3*y1*y2*y3**3*y5/2.D0 +y1**3*y2*y3*y4 &
            +SQRT3*y1*y2**3*y3*y5/2.D0 -y1*y2**3*y3*y4/2.D0 &
            -y1*y2*y3**3*y4/2.D0 )*DFEA111234 &
            +( y3*y4**4*y5/3.D0 &
            +y3*y4**5*SQRT3/18.D0 -y2*y4**4*y5/3.D0 -y2*y4*y5**4*SQRT3/2.D0 &
            -y3*y4**2*y5**3 &
            +2.D0/9.D0*y1*y4**5*SQRT3 +y2*y4**5*SQRT3/18.D0 &
            +y2*y4**2*y5**3 -2.D0/3.D0*y1*y4**3*y5**2*SQRT3 &
            -y3*y4*y5**4*SQRT3/2.D0 )*DFEA244555 &
            +( y1*y2*y4**2*y5**2 -3.D0/4.D0*y2*y3*y4**4 -y1*y2*y5**4 -y1*y3*y5**4 &
            +5.D0/4.D0*y2*y3*y5**4 &
            +y1*y3*y4**2*y5**2 -7.D0/2.D0*y2*y3*y4**2*y5**2 &
            -2.D0*y1*y2*y4**3*y5*SQRT3 &
            +2.D0*y1*y3*y4**3*y5*SQRT3 )*DFEA124455 
 
         ds3(6)=ds1(6) +( y2**6 +y1**6 &
            +y3**6 )*DFEA333333 +( y1*y2**4*y3 +y1**4*y2*y3 &
            +y1*y2*y3**4 )*DFEA111123 +y1**2*y2**2*y3**2*DFEA112233 +( y1**4*y4**2 &
            +y2**4*y4**2 +y2**4*y5**2 +y3**4*y4**2 +y1**4*y5**2 &
            +y3**4*y5**2 )*DFEA222255 
 
         ds4(6)=ds3(6) +( 3.D0*y1*y3*y5**4 &
            +y1*y3*y4**4 &
            +9.D0*y2*y3*y4**2*y5**2 -3.D0/2.D0*y2*y3*y5**4 &
            -4.D0*y1*y3*y4**3*y5*SQRT3 &
            +y1*y2*y4**4 +4.D0*y1*y2*y4**3*y5*SQRT3 +3.D0*y1*y2*y5**4 &
            +5.D0/2.D0*y2*y3*y4**4 )*DFEA134444 &
            +(  -y1*y3**2*y5**3*SQRT3/3.D0 -7.D0/3.D0*y1**2*y3*y4*y5**2 &
            +5.D0/3.D0*y1*y2**2*y4**2*y5*SQRT3 -13.D0/3.D0*y2**2*y3*y4*y5**2 &
            -4.D0/3.D0*y2*y3**2*y5**3*SQRT3 -7.D0/3.D0*y1**2*y2*y4*y5**2 &
            -16.D0/3.D0*y1*y3**2*y4*y5**2 &
            +4.D0/3.D0*y1**2*y3*y4**2*y5*SQRT3 +4.D0/3.D0*y2**2*y3*y5**3*SQRT3 &
            +3.D0*y1**2*y2*y4**3 +y2*y3**2*y4**3 +y1*y2**2*y5**3*SQRT3/3.D0 &
            +y2**2*y3*y4**3 -13.D0/3.D0*y2*y3**2*y4*y5**2 &
            -5.D0/3.D0*y1*y3**2*y4**2*y5*SQRT3 -4.D0/3.D0*y1**2*y2*y4**2*y5*SQRT3 &
            +3.D0*y1**2*y3*y4**3 &
            -16.D0/3.D0*y1*y2**2*y4*y5**2 )*DFEA233444 
 
         ds5(6)=ds4(6) &
            +( 2.D0*y1*y3**2*y5**3 +4.D0*y2*y3**2*y5**3 &
            +4.D0*y2**2*y3*y4*y5**2*SQRT3 -2.D0*y1*y2**2*y5**3 &
            +y1**2*y3*y4*y5**2*SQRT3 &
            +6.D0*y1*y3**2*y4**2*y5 -6.D0*y1*y2**2*y4**2*y5 -3.D0*y1**2*y3*y4**2*y5 &
            +y1**2*y2*y4*y5**2*SQRT3 &
            +4.D0*y1*y3**2*y4*y5**2*SQRT3 -3.D0*y1**2*y2*y4**3*SQRT3 &
            -4.D0*y2**2*y3*y5**3 &
            +3.D0*y1**2*y2*y4**2*y5 -y1**2*y2*y5**3 &
            +y1**2*y3*y5**3 -3.D0*y1**2*y3*y4**3*SQRT3 &
            +4.D0*y2*y3**2*y4*y5**2*SQRT3 &
            +4.D0*y1*y2**2*y4*y5**2*SQRT3 )*DFEA113555 
 
         ds2(6)=ds5(6) &
            +(  -2.D0/3.D0*y3**2*y4**4*SQRT3 -3.D0/2.D0*y1**2*y4**2*y5**2*SQRT3 &
            -3.D0/4.D0*y1**2*y5**4*SQRT3 -y2**2*y4**3*y5 &
            +7.D0/12.D0*y1**2*y4**4*SQRT3 +y3**2*y4**3*y5 &
            +3.D0*y3**2*y4*y5**3 -2.D0/3.D0*y2**2*y4**4*SQRT3 &
            -3.D0*y2**2*y4*y5**3 )*DFEA334445 &
            +(  -3.D0*y1*y3*y4**3*y5 +2.D0/3.D0*y1*y2*y5**4*SQRT3 -y1*y3*y4*y5**3 &
            +2.D0/3.D0*y1*y3*y5**4*SQRT3 &
            +3.D0*y1*y2*y4**3*y5 -7.D0/12.D0*y2*y3*y5**4*SQRT3 &
            +3.D0/2.D0*y2*y3*y4**2*y5**2*SQRT3 +y1*y2*y4*y5**3 &
            +3.D0/4.D0*y2*y3*y4**4*SQRT3 )*DFEA124555 &
            +( 2.D0*y3**2*y4*y5**3*SQRT3 -7.D0/2.D0*y1**2*y4**2*y5**2 &
            +y2**2*y4**2*y5**2 -y2**2*y4**4 -y3**2*y4**4 &
            -2.D0*y2**2*y4*y5**3*SQRT3 -3.D0/4.D0*y1**2*y5**4 &
            +5.D0/4.D0*y1**2*y4**4 +y3**2*y4**2*y5**2 )*DFEA334455 
 
         ds3(6)=ds2(6) &
            +(  -6.D0*y4**2*y5**4 +9.D0*y4**4*y5**2 +y5**6 )*DFEA555555 &
            +( y2*y3**3*y4**2 +y2*y3**3*y5**2 +y1*y3**3*y4**2 +y1*y2**3*y4**2 &
            +y1**3*y2*y4**2 +y1*y2**3*y5**2 +y1**3*y3*y5**2 +y1**3*y3*y4**2 &
            +y1**3*y2*y5**2 +y2**3*y3*y4**2 +y1*y3**3*y5**2 &
            +y2**3*y3*y5**2 )*DFEA233344 &
            +( y1*y2**3*y5**2*SQRT3/6.D0 -y2**3*y3*y5**2*SQRT3/3.D0 &
            -y2*y3**3*y5**2*SQRT3/3.D0 &
            +y1**3*y2*y4*y5 -y1**3*y2*y5**2*SQRT3/3.D0 -y1**3*y3*y4*y5 &
            -y1**3*y3*y5**2*SQRT3/3.D0 -y1*y3**3*y4**2*SQRT3/2.D0 &
            +y1*y3**3*y5**2*SQRT3/6.D0 -y2**3*y3*y4*y5 &
            +y2*y3**3*y4*y5 -y1*y2**3*y4**2*SQRT3/2.D0 )*DFEA233345 &
            +(  -3.D0*y2**3*y4*y5**2 &
            +y3**3*y4**3 -3.D0*y3**3*y4*y5**2 -3.D0*y1**3*y4*y5**2 +y2**3*y4**3 &
            +y1**3*y4**3 )*DFEA111444 +( y1*y2**3*y3**2 +y1**3*y2**2*y3 &
            +y1**2*y2**3*y3 +y1*y2**2*y3**3 +y1**2*y2*y3**3 &
            +y1**3*y2*y3**2 )*DFEA111233 
 
         ds4(6)=ds3(6) &
            +( 9.D0*y4**2*y5**4 -6.D0*y4**4*y5**2 +y4**6 )*DFEA444444 &
            +(  -5.D0/3.D0*y1*y2**2*y4**2*y5*SQRT3 &
            +y1*y2**2*y4**3 -4.D0/3.D0*y1**2*y3*y4**2*y5*SQRT3 &
            -2.D0*y1**2*y2*y4**3 -y1*y2**2*y5**3*SQRT3/3.D0 &
            +4.D0/3.D0*y2**2*y3*y4*y5**2 -4.D0/3.D0*y2**2*y3*y5**3*SQRT3 &
            -2.D0*y1**2*y3*y4**3 &
            +7.D0/3.D0*y1*y2**2*y4*y5**2 -2.D0/3.D0*y1**2*y3*y4*y5**2 &
            +y1*y3**2*y4**3 +4.D0/3.D0*y2*y3**2*y5**3*SQRT3 &
            +y1*y3**2*y5**3*SQRT3/3.D0 +4.D0/3.D0*y1**2*y2*y4**2*y5*SQRT3 &
            +4.D0/3.D0*y2*y3**2*y4*y5**2 &
            +5.D0/3.D0*y1*y3**2*y4**2*y5*SQRT3 -2.D0/3.D0*y1**2*y2*y4*y5**2 &
            +7.D0/3.D0*y1*y3**2*y4*y5**2 )*DFEA133444 
 
         ds5(6)=ds4(6) &
            +(  -y1**3*y2*y4*y5 +2.D0/3.D0*y2**3*y3*y5**2*SQRT3 &
            +y1*y3**3*y4**2*SQRT3/2.D0 +y1**3*y3*y4**2*SQRT3/2.D0 &
            +y1**3*y3*y5**2*SQRT3/6.D0 +y1**3*y2*y5**2*SQRT3/6.D0 +y1**3*y3*y4*y5 &
            +y1*y2**3*y5**2*SQRT3/6.D0 +y1**3*y2*y4**2*SQRT3/2.D0 &
            +2.D0/3.D0*y2*y3**3*y5**2*SQRT3 -y1*y2**3*y4*y5 &
            +y1*y2**3*y4**2*SQRT3/2.D0 +y1*y3**3*y5**2*SQRT3/6.D0 &
            +y1*y3**3*y4*y5 )*DFEA133345 
 
         dv6(6)=ds5(6) +(  -y2**2*y3*y4**2*y5 &
            +y1**2*y3*y4*y5**2*SQRT3/3.D0 +y2*y3**2*y4**2*y5 &
            +y2*y3**2*y5**3 -y1*y2**2*y5**3 +4.D0/3.D0*y2**2*y3*y4*y5**2*SQRT3 &
            +4.D0/3.D0*y2*y3**2*y4*y5**2*SQRT3 -y1*y2**2*y4**2*y5 &
            +4.D0/3.D0*y1*y3**2*y4*y5**2*SQRT3 -y2**2*y3*y5**3 +y1*y3**2*y5**3 &
            +y1**2*y2*y4*y5**2*SQRT3/3.D0 -y1**2*y2*y4**3*SQRT3 &
            +y1*y3**2*y4**2*y5 -y1**2*y3*y4**3*SQRT3 &
            +4.D0/3.D0*y1*y2**2*y4*y5**2*SQRT3 )*DFEA233445 &
            +( y2*y3**4*y4*SQRT3 -y1**4*y2*y5 &
            +y2**4*y3*y4*SQRT3 -y1**4*y3*y4*SQRT3 +y2*y3**4*y5 -2.D0*y1*y2**4*y5 &
            +2.D0*y1*y3**4*y5 -y1**4*y2*y4*SQRT3 &
            +y1**4*y3*y5 -y2**4*y3*y5 )*DFEA233335 +( y2**2*y3**4 +y1**4*y3**2 &
            +y1**2*y2**4 +y2**4*y3**2 +y1**2*y3**4 &
            +y1**4*y2**2 )*DFEA222233   
 
         dv(6)=dv0 +dv1(6) +dv2(6) +dv3(6) &
            +dv4(6) +dv5(6) +dv6(6)
 
         if (flag > 1) then
            ! 2nd derivatives of FEA
 
            ddfea1= f2a1*2d0+f3a1*6d0*coro+f4a1*12d0*coro**2+f5a1*20d0*coro**3+f6a1*30d0*coro**4
 
            ddfea11= f2a11*2d0+f3a11*6d0*coro+f4a11*12d0*coro**2
            ddfea12= f2a12*2d0+f3a12*6d0*coro+f4a12*12d0*coro**2
            ddfea14= f2a14*2d0+f3a14*6d0*coro+f4a14*12d0*coro**2
            ddfea44= f2a44*2d0+f3a44*6d0*coro+f4a44*12d0*coro**2
 
            ddfea111= f2a111*2d0+f3a111*6d0*coro
            ddfea112= f2a112*2d0+f3a112*6d0*coro
            ddfea114= f2a114*2d0+f3a114*6d0*coro
            ddfea123= f2a123*2d0+f3a123*6d0*coro
            ddfea124= f2a124*2d0+f3a124*6d0*coro
            ddfea144= f2a144*2d0+f3a144*6d0*coro
            ddfea155= f2a155*2d0+f3a155*6d0*coro
            ddfea455= f2a455*2d0+f3a455*6d0*coro
 
            ddfea1111= f2a1111*2d0
            ddfea1112= f2a1112*2d0
            ddfea1114= f2a1114*2d0
            ddfea1122= f2a1122*2d0
            ddfea1123= f2a1123*2d0
            ddfea1124= f2a1124*2d0
            ddfea1125= f2a1125*2d0
            ddfea1144= f2a1144*2d0
            ddfea1155= f2a1155*2d0
            ddfea1244= f2a1244*2d0
            ddfea1255= f2a1255*2d0
            ddfea1444= f2a1444*2d0
            ddfea1455= f2a1455*2d0
            ddfea4444= f2a4444*2d0
 
            ddfea44444 = f2a44444 *2d0
            ddfea33455 = f2a33455 *2d0
            ddfea33445 = f2a33445 *2d0
            ddfea33345 = f2a33345 *2d0
            ddfea33344 = f2a33344 *2d0
            ddfea33334 = f2a33334 *2d0
            ddfea33333 = f2a33333 *2d0
            ddfea25555 = f2a25555 *2d0
            ddfea24455 = f2a24455 *2d0
            ddfea24445 = f2a24445 *2d0
            ddfea23333 = f2a23333 *2d0
            ddfea13455 = f2a13455 *2d0
            ddfea13445 = f2a13445 *2d0
            ddfea13345 = f2a13345 *2d0
            ddfea12355 = f2a12355 *2d0
            ddfea11334 = f2a11334 *2d0
            ddfea11333 = f2a11333 *2d0
            ddfea11255 = f2a11255 *2d0
            ddfea11245 = f2a11245 *2d0
            ddfea11234 = f2a11234 *2d0
            ddfea11233 = f2a11233 *2d0
            ddfea11135 = f2a11135 *2d0
            ddfea11134 = f2a11134 *2d0
            ddfea11123 = f2a11123 *2d0
            ddfea555555= f2a555555*2d0
            ddfea444444= f2a444444*2d0
            ddfea335555= f2a335555*2d0
            ddfea334455= f2a334455*2d0
            ddfea334445= f2a334445*2d0
            ddfea333555= f2a333555*2d0
            ddfea333333= f2a333333*2d0
            ddfea244555= f2a244555*2d0
            ddfea244455= f2a244455*2d0
            ddfea233445= f2a233445*2d0
            ddfea233444= f2a233444*2d0
            ddfea233345= f2a233345*2d0
            ddfea233344= f2a233344*2d0
            ddfea233335= f2a233335*2d0
            ddfea223355= f2a223355*2d0
            ddfea222335= f2a222335*2d0
            ddfea222334= f2a222334*2d0
            ddfea222333= f2a222333*2d0
            ddfea222255= f2a222255*2d0
            ddfea222245= f2a222245*2d0
            ddfea222233= f2a222233*2d0
            ddfea222224= f2a222224*2d0
            ddfea145555= f2a145555*2d0
            ddfea134444= f2a134444*2d0
            ddfea133444= f2a133444*2d0
            ddfea133345= f2a133345*2d0
            ddfea133334= f2a133334*2d0
            ddfea133333= f2a133333*2d0
            ddfea124555= f2a124555*2d0
            ddfea124455= f2a124455*2d0
            ddfea123455= f2a123455*2d0
            ddfea123345= f2a123345*2d0
            ddfea113555= f2a113555*2d0
            ddfea113345= f2a113345*2d0
            ddfea112355= f2a112355*2d0
            ddfea112335= f2a112335*2d0
            ddfea112233= f2a112233*2d0
            ddfea111444= f2a111444*2d0
            ddfea111234= f2a111234*2d0
            ddfea111233= f2a111233*2d0
            ddfea111123= f2a111123*2d0
 
            ! calculate 2nd derivatives of potential energy (hessian)
            ddv1(1,1)=0d0 
            ddv2(1,1)=0d0 +( +2d0 )*FEA11 +0d0 +0d0 
            ddv3(1,1)=0d0 &
               +0d0 +( +2d0*y3 +2d0*y2 )*FEA112 +0d0 +0d0 +( +3d0*2d0*y1 )*FEA111 +( &
               +2d0*y4 )*FEA114 
            dds2(1,1)=0d0 +0d0 +( &
               +2d0*y5**2 )*FEA1155 
            dds1(1,1)=dds2(1,1) +( +2d0*y4**2 )*FEA1144 &
               +( 3d0*2d0*y1*y4 )*FEA1114 +( +4d0*3d0*y1**2 )*FEA1111 &
               +0d0 
            dds2(1,1)=dds1(1,1) +0d0 +( +2d0*y3*y4 +2d0*y2*y4 )*FEA1124 &
               +( 2d0*y2*y5 -2d0*y3*y5 )*FEA1125 
            ddv4(1,1)=dds2(1,1) +( &
               +3d0*2d0*y1*y3 +3d0*2d0*y1*y2 )*FEA1112 +( +2d0*y3**2 &
               +2d0*y2**2 )*FEA1122 +( +2d0*y2*y3 )*FEA1123 +0d0 
            dds3(1,1)=0d0 +0d0 &
               +0d0 
            dds2(1,1)=dds3(1,1) +0d0 &
               +( -5.D0/9.D0*2d0*y4**3*SQRT3 -2d0*y4*y5**2*SQRT3 )*FEA33445 &
               +( -2d0*y4**3/3.D0 +2d0*y4*y5**2 )*FEA33455 
            dds1(1,1)=dds2(1,1) +( &
               +3d0*2d0*y1*y4**2*SQRT3/2.D0 -3d0*2d0*y1*y5**2*SQRT3/6.D0 )*FEA33345 &
               +( +3d0*2d0*y1*y5**2 +3d0*2d0*y1*y4**2 )*FEA33344 &
               +( -2.D0*4d0*3d0*y1**2*y4 )*FEA33334 +( +5d0*4d0*y1**3 )*FEA33333 +0d0 &
               +0d0 
            dds3(1,1)=dds1(1,1) +( 2d0*y3*y5**2 +2d0*y2*y5**2 +2d0*y3*y4**2 &
               +2d0*y2*y4**2 )*FEA11255 +( 2.D0/3.D0*2d0*y3*y4**2*SQRT3 &
               +2.D0/3.D0*2d0*y2*y4**2*SQRT3 )*FEA13345 
            dds4(1,1)=dds3(1,1) &
               +( 2d0*y2*y4*y5 +2d0*y3*y4**2*SQRT3/3.D0 &
               +2d0*y2*y4**2*SQRT3/3.D0 -2d0*y3*y4*y5 )*FEA11245 
            dds2(1,1)=dds4(1,1) &
               +( -3d0*2d0*y1*y2*y5 +3d0*2d0*y1*y3*y5 )*FEA11135 +( 3d0*2d0*y1*y3*y4 &
               +3d0*2d0*y1*y2*y4 )*FEA11134 
            ddv5(1,1)=dds2(1,1) +( +4d0*3d0*y1**2*y3 &
               +4d0*3d0*y1**2*y2 )*FEA23333 +( +2d0*y2**2*y4 -SQRT3*2d0*y3**2*y5 &
               +SQRT3*2d0*y2**2*y5 +2d0*y3**2*y4 )*FEA11334 +( 2d0*y3**3 &
               +3d0*2d0*y1*y3**2 +2d0*y2**3 +3d0*2d0*y1*y2**2 )*FEA11333 +0d0 +( &
               +2d0*y2*y3*y4 )*FEA11234 +( +3d0*2d0*y1*y2*y3 )*FEA11123 &
               +( 2d0*y2**2*y3 &
               +2d0*y2*y3**2 )*FEA11233 
            dds3(1,1)=( &
               -8.D0/3.D0*3d0*2d0*y1*y4*y5**2*SQRT3 )*FEA333555 &
               +( 4d0*3d0*y1**2*y5**2*SQRT3/2.D0 &
               -4d0*3d0*y1**2*y4**2*SQRT3/6.D0 )*FEA222245 &
               +( +5d0*4d0*y1**3*y3 +5d0*4d0*y1**3*y2 )*FEA133333 &
               +( 4d0*3d0*y1**2*y3*y4 +4d0*3d0*y1**2*y2*y4 &
               +4d0*3d0*y1**2*y2*y5*SQRT3 -4d0*3d0*y1**2*y3*y5*SQRT3 )*FEA133334 &
               +0d0 
            dds4(1,1)=dds3(1,1) &
               +( -2d0*y2**2*y3*y5 -SQRT3*2d0*y2**2*y3*y4/3.D0 &
               +2d0*y2*y3**2*y5 -SQRT3*2d0*y2*y3**2*y4/3.D0 )*FEA112335 +( &
               +2d0*y2*y3*y4**2 +2d0*y2*y3*y5**2 )*FEA112355 
            dds2(1,1)=dds4(1,1) &
               +( -3d0*2d0*y1*y2**2*y5/2.D0 -2d0*y3**3*y5/2.D0 &
               +3d0*2d0*y1*y2**2*y4*SQRT3/2.D0 -2d0*y2**3*y4*SQRT3/2.D0 &
               +3d0*2d0*y1*y3**2*y5/2.D0 +2d0*y2**3*y5/2.D0 &
               +3d0*2d0*y1*y3**2*y4*SQRT3/2.D0 -2d0*y3**3*y4*SQRT3/2.D0 )*FEA222335 &
               +( -2d0*y2**2*y5**2*SQRT3/2.D0 -2d0*y3**2*y5**2*SQRT3/2.D0 &
               -2d0*y2**2*y4**2*SQRT3/6.D0 -2d0*y2**2*y4*y5 &
               +2d0*y3**2*y4*y5 -2d0*y3**2*y4**2*SQRT3/6.D0 )*FEA113345 +( &
               +2d0*y2**2*y5**2 +2d0*y3**2*y4**2 +2d0*y3**2*y5**2 &
               +2d0*y2**2*y4**2 )*FEA223355 
            dds3(1,1)=dds2(1,1) +( &
               +2.D0/3.D0*2d0*y2*y3*y4**2*SQRT3 )*FEA123345 &
               +( -3d0*2d0*y1*y2**2*y5*SQRT3/2.D0 -3d0*2d0*y1*y2**2*y4/2.D0 &
               -3d0*2d0*y1*y3**2*y4/2.D0 -2d0*y2**3*y4/2.D0 &
               +3d0*2d0*y1*y3**2*y5*SQRT3/2.D0 -2d0*y3**3*y4/2.D0 &
               -2d0*y2**3*y5*SQRT3/2.D0 &
               +2d0*y3**3*y5*SQRT3/2.D0 )*FEA222334 +( +5.D0/2.D0*2d0*y5**4 &
               +9.D0*2d0*y4**2*y5**2 -3.D0/2.D0*2d0*y4**4 )*FEA335555 &
               +( 3d0*2d0*y1*y2**3 +3d0*2d0*y1*y3**3 )*FEA222333 
            dds4(1,1)=dds3(1,1) &
               +0d0 +( -2.D0*5d0*4d0*y1**3*y4 )*FEA222224 
            dds5(1,1)=dds4(1,1) &
               +0d0 
            dds1(1,1)=dds5(1,1) +( +3d0*2d0*y1*y2*y3*y4 )*FEA111234 +0d0 &
               +0d0 
            dds3(1,1)=dds1(1,1) +( +6d0*5d0*y1**4 )*FEA333333 +( &
               +4d0*3d0*y1**2*y2*y3 )*FEA111123 +2d0*y2**2*y3**2*FEA112233 &
               +( 4d0*3d0*y1**2*y4**2 &
               +4d0*3d0*y1**2*y5**2 )*FEA222255 
            dds4(1,1)=dds3(1,1) +0d0 &
               +( -7.D0/3.D0*2d0*y3*y4*y5**2 -7.D0/3.D0*2d0*y2*y4*y5**2 &
               +4.D0/3.D0*2d0*y3*y4**2*y5*SQRT3 &
               +3.D0*2d0*y2*y4**3 -4.D0/3.D0*2d0*y2*y4**2*y5*SQRT3 &
               +3.D0*2d0*y3*y4**3 )*FEA233444 
            dds5(1,1)=dds4(1,1) +( &
               +2d0*y3*y4*y5**2*SQRT3 -3.D0*2d0*y3*y4**2*y5 &
               +2d0*y2*y4*y5**2*SQRT3 -3.D0*2d0*y2*y4**3*SQRT3 &
               +3.D0*2d0*y2*y4**2*y5 -2d0*y2*y5**3 &
               +2d0*y3*y5**3 -3.D0*2d0*y3*y4**3*SQRT3 )*FEA113555 
            dds2(1,1)=dds5(1,1) &
               +( -3.D0/2.D0*2d0*y4**2*y5**2*SQRT3 -3.D0/4.D0*2d0*y5**4*SQRT3 &
               +7.D0/12.D0*2d0*y4**4*SQRT3 )*FEA334445 +0d0 &
               +( -7.D0/2.D0*2d0*y4**2*y5**2 -3.D0/4.D0*2d0*y5**4 &
               +5.D0/4.D0*2d0*y4**4 )*FEA334455 
            dds3(1,1)=dds2(1,1) +0d0 +( &
               +3d0*2d0*y1*y2*y4**2 +3d0*2d0*y1*y3*y5**2 +3d0*2d0*y1*y3*y4**2 &
               +3d0*2d0*y1*y2*y5**2 )*FEA233344 +( &
               +3d0*2d0*y1*y2*y4*y5 -3d0*2d0*y1*y2*y5**2*SQRT3/3.D0 &
               -3d0*2d0*y1*y3*y4*y5 -3d0*2d0*y1*y3*y5**2*SQRT3/3.D0 )*FEA233345 &
               +( -3.D0*3d0*2d0*y1*y4*y5**2 +3d0*2d0*y1*y4**3 )*FEA111444 +( &
               +3d0*2d0*y1*y2**2*y3 +2d0*y2**3*y3 +2d0*y2*y3**3 &
               +3d0*2d0*y1*y2*y3**2 )*FEA111233 
            dds4(1,1)=dds3(1,1) +0d0 &
               +( -4.D0/3.D0*2d0*y3*y4**2*y5*SQRT3 -2.D0*2d0*y2*y4**3 &
               -2.D0*2d0*y3*y4**3 -2.D0/3.D0*2d0*y3*y4*y5**2 &
               +4.D0/3.D0*2d0*y2*y4**2*y5*SQRT3 &
               -2.D0/3.D0*2d0*y2*y4*y5**2 )*FEA133444 
            dds5(1,1)=dds4(1,1) &
               +( -3d0*2d0*y1*y2*y4*y5 +3d0*2d0*y1*y3*y4**2*SQRT3/2.D0 &
               +3d0*2d0*y1*y3*y5**2*SQRT3/6.D0 +3d0*2d0*y1*y2*y5**2*SQRT3/6.D0 &
               +3d0*2d0*y1*y3*y4*y5 &
               +3d0*2d0*y1*y2*y4**2*SQRT3/2.D0 )*FEA133345 
            ddv6(1,1)=dds5(1,1) +( &
               +2d0*y3*y4*y5**2*SQRT3/3.D0 &
               +2d0*y2*y4*y5**2*SQRT3/3.D0 -2d0*y2*y4**3*SQRT3 &
               -2d0*y3*y4**3*SQRT3 )*FEA233445 &
               +( -4d0*3d0*y1**2*y2*y5 -4d0*3d0*y1**2*y3*y4*SQRT3 &
               -4d0*3d0*y1**2*y2*y4*SQRT3 &
               +4d0*3d0*y1**2*y3*y5 )*FEA233335 +( +4d0*3d0*y1**2*y3**2 +2d0*y2**4 &
               +2d0*y3**4 +4d0*3d0*y1**2*y2**2 )*FEA222233 
 
            ddv(1,1)= +ddv1(1,1) &
               +ddv2(1,1) +ddv3(1,1) +ddv4(1,1) +ddv5(1,1) +ddv6(1,1)
 
            ddv1(1,2)=0d0 
            ddv2(1,2)=( +1d0 )*FEA12 +0d0 +0d0 +0d0 
            ddv3(1,2)=( +y4 &
               +SQRT3*y5 )*FEA124 +0d0 +( +2d0*y2 +2d0*y1 )*FEA112 +0d0 +y3*FEA123 &
               +0d0 +0d0 +0d0 
            dds2(1,2)=0d0 +0d0 +0d0 
            dds1(1,2)=dds2(1,2) +0d0 +0d0 &
               +0d0 +( +y4**2 -SQRT3*y4*y5 )*FEA1244 
            dds2(1,2)=dds1(1,2) +( +y5**2 &
               +SQRT3*y4*y5 )*FEA1255 +( +2d0*y1*y4 &
               +SQRT3*2d0*y2*y5/2.D0 -2d0*y2*y4/2.D0 )*FEA1124 +( 2d0*y1*y5 &
               +SQRT3*2d0*y2*y4/2.D0 +2d0*y2*y5/2.D0 )*FEA1125 
            ddv4(1,2)=dds2(1,2) +( &
               +3d0*y1**2 +3d0*y2**2 )*FEA1112 +( +2d0*y1*2d0*y2 )*FEA1122 &
               +( 2d0*y2*y3 +2d0*y1*y3 +y3**2 )*FEA1123 +0d0 
            dds3(1,2)=0d0 +0d0 &
               +0d0 
            dds2(1,2)=dds3(1,2) +0d0 +0d0 +0d0 
            dds1(1,2)=dds2(1,2) +0d0 +0d0 &
               +0d0 +0d0 +( -4.D0/9.D0*y4**3*SQRT3 -y5**3 -y4**2*y5 )*FEA13445 +( &
               +y4*y5**2 -y4**3/3.D0 )*FEA13455 
            dds3(1,2)=dds1(1,2) +( +2d0*y2*y5**2 &
               +2d0*y1*y5**2 +2d0*y2*y4**2 +2d0*y1*y4**2 )*FEA11255 +( &
               +2d0*y2*y5**2*SQRT3/2.D0 -2d0*y2*y4*y5 +2d0*y2*y4**2*SQRT3/6.D0 &
               +2.D0/3.D0*2d0*y1*y4**2*SQRT3 )*FEA13345 
            dds4(1,2)=dds3(1,2) &
               +( 2d0*y1*y4*y5 +2d0*y1*y4**2*SQRT3/3.D0 -2d0*y2*y4**2*SQRT3/6.D0 &
               +2d0*y2*y5**2*SQRT3/2.D0 )*FEA11245 
            dds2(1,2)=dds4(1,2) &
               +( -3d0*y1**2*y5 -3d0*y2**2*y4*SQRT3/2.D0 -3d0*y2**2*y5/2.D0 )*FEA11135 &
               +( +3d0*y1**2*y4 &
               +3d0*y2**2*y5*SQRT3/2.D0 &
               -3d0*y2**2*y4/2.D0 )*FEA11134 
            ddv5(1,2)=dds2(1,2) &
               +( 4d0*y2**3 +4d0*y1**3 )*FEA23333 +( +2d0*y1*2d0*y2*y4 &
               +SQRT3*2d0*y1*2d0*y2*y5 )*FEA11334 +( +2d0*y1*3d0*y2**2 &
               +3d0*y1**2*2d0*y2 )*FEA11333 +( y3*y4**2 +y3*y5**2 )*FEA12355 &
               +( -y3**2*y4/2.D0 -2d0*y2*y3*y4/2.D0 -SQRT3*y3**2*y5/2.D0 &
               +2d0*y1*y3*y4 +SQRT3*2d0*y2*y3*y5/2.D0 )*FEA11234 +( 3d0*y2**2*y3 &
               +y3**3 +3d0*y1**2*y3 )*FEA11123 +( 2d0*y1*2d0*y2*y3 +2d0*y2*y3**2 &
               +2d0*y1*y3**2 )*FEA11233 
            dds3(1,2)=0d0 +0d0 +( +5d0*y2**4 &
               +5d0*y1**4 )*FEA133333 +( +4d0*y1**3*y4 +4d0*y2**3*y5*SQRT3 &
               +4d0*y1**3*y5*SQRT3 +4d0*y2**3*y4 )*FEA133334 +( -y3*y4**3/3.D0 &
               +y3*y4*y5**2 )*FEA123455 
            dds4(1,2)=dds3(1,2) &
               +( 2.D0/3.D0*SQRT3*2d0*y2*y3**2*y4 -2d0*y1*2d0*y2*y3*y5 &
               -SQRT3*2d0*y1*2d0*y2*y3*y4/3.D0 &
               +2d0*y1*y3**2*y5 -SQRT3*2d0*y1*y3**2*y4/3.D0 )*FEA112335 &
               +( 2d0*y2*y3*y5**2 +y3**2*y5**2 +y3**2*y4**2 +2d0*y2*y3*y4**2 &
               +2d0*y1*y3*y4**2 +2d0*y1*y3*y5**2 )*FEA112355 
            dds2(1,2)=dds4(1,2) &
               +( -3d0*y1**2*2d0*y2*y5/2.D0 &
               +3d0*y1**2*2d0*y2*y4*SQRT3/2.D0 -2d0*y1*3d0*y2**2*y4*SQRT3/2.D0 &
               +2d0*y1*3d0*y2**2*y5/2.D0 )*FEA222335 &
               +( -2d0*y1*2d0*y2*y5**2*SQRT3/2.D0 -2d0*y1*2d0*y2*y4**2*SQRT3/6.D0 &
               -2d0*y1*2d0*y2*y4*y5 )*FEA113345 &
               +( +2d0*y1*2d0*y2*y5**2 &
               +2d0*y1*2d0*y2*y4**2 )*FEA223355 
            dds3(1,2)=dds2(1,2) &
               +( y3**2*y4**2*SQRT3/6.D0 +y3**2*y4*y5 +y3**2*y5**2*SQRT3/2.D0 &
               +2.D0/3.D0*2d0*y1*y3*y4**2*SQRT3 -2d0*y2*y3*y4*y5 &
               +2d0*y2*y3*y4**2*SQRT3/6.D0 +2d0*y2*y3*y5**2*SQRT3/2.D0 )*FEA123345 &
               +( -3d0*y1**2*2d0*y2*y5*SQRT3/2.D0 -3d0*y1**2*2d0*y2*y4/2.D0 &
               -2d0*y1*3d0*y2**2*y4/2.D0 -2d0*y1*3d0*y2**2*y5*SQRT3/2.D0 )*FEA222334 &
               +0d0 +( 3d0*y1**2*3d0*y2**2 )*FEA222333 
            dds4(1,2)=dds3(1,2) +0d0 &
               +0d0 
            dds5(1,2)=dds4(1,2) +0d0 
            dds1(1,2)=dds5(1,2) &
               +( -SQRT3*y3**3*y5/2.D0 +3d0*y1**2*y3*y4 &
               +SQRT3*3d0*y2**2*y3*y5/2.D0 -3d0*y2**2*y3*y4/2.D0 &
               -y3**3*y4/2.D0 )*FEA111234 &
               +0d0 &
               +( y4**2*y5**2 -y5**4 &
               -2.D0*y4**3*y5*SQRT3 )*FEA124455 
            dds3(1,2)=dds1(1,2) &
               +0d0 +( 4d0*y2**3*y3 +4d0*y1**3*y3 +y3**4 )*FEA111123 &
               +2d0*y1*2d0*y2*y3**2*FEA112233 +0d0 
            dds4(1,2)=dds3(1,2) +( +y4**4 &
               +4.D0*y4**3*y5*SQRT3 +3.D0*y5**4 )*FEA134444 +( &
               +5.D0/3.D0*2d0*y2*y4**2*y5*SQRT3 -7.D0/3.D0*2d0*y1*y4*y5**2 &
               +3.D0*2d0*y1*y4**3 &
               +2d0*y2*y5**3*SQRT3/3.D0 -4.D0/3.D0*2d0*y1*y4**2*y5*SQRT3 &
               -16.D0/3.D0*2d0*y2*y4*y5**2 )*FEA233444 
            dds5(1,2)=dds4(1,2) &
               +( -2.D0*2d0*y2*y5**3 -6.D0*2d0*y2*y4**2*y5 &
               +2d0*y1*y4*y5**2*SQRT3 -3.D0*2d0*y1*y4**3*SQRT3 &
               +3.D0*2d0*y1*y4**2*y5 -2d0*y1*y5**3 &
               +4.D0*2d0*y2*y4*y5**2*SQRT3 )*FEA113555 
            dds2(1,2)=dds5(1,2) +0d0 +( &
               +2.D0/3.D0*y5**4*SQRT3 +3.D0*y4**3*y5 +y4*y5**3 )*FEA124555 &
               +0d0 
            dds3(1,2)=dds2(1,2) +0d0 +( +3d0*y2**2*y4**2 +3d0*y1**2*y4**2 &
               +3d0*y2**2*y5**2 +3d0*y1**2*y5**2 )*FEA233344 &
               +( 3d0*y2**2*y5**2*SQRT3/6.D0 &
               +3d0*y1**2*y4*y5 -3d0*y1**2*y5**2*SQRT3/3.D0 &
               -3d0*y2**2*y4**2*SQRT3/2.D0 )*FEA233345 &
               +0d0 +( 3d0*y2**2*y3**2 +3d0*y1**2*2d0*y2*y3 +2d0*y1*3d0*y2**2*y3 &
               +2d0*y2*y3**3 +2d0*y1*y3**3 &
               +3d0*y1**2*y3**2 )*FEA111233 
            dds4(1,2)=dds3(1,2) +0d0 &
               +( -5.D0/3.D0*2d0*y2*y4**2*y5*SQRT3 &
               +2d0*y2*y4**3 -2.D0*2d0*y1*y4**3 -2d0*y2*y5**3*SQRT3/3.D0 &
               +7.D0/3.D0*2d0*y2*y4*y5**2 &
               +4.D0/3.D0*2d0*y1*y4**2*y5*SQRT3 &
               -2.D0/3.D0*2d0*y1*y4*y5**2 )*FEA133444 
            dds5(1,2)=dds4(1,2) &
               +( -3d0*y1**2*y4*y5 +3d0*y1**2*y5**2*SQRT3/6.D0 &
               +3d0*y2**2*y5**2*SQRT3/6.D0 &
               +3d0*y1**2*y4**2*SQRT3/2.D0 -3d0*y2**2*y4*y5 &
               +3d0*y2**2*y4**2*SQRT3/2.D0 )*FEA133345 
            ddv6(1,2)=dds5(1,2) &
               +( -2d0*y2*y5**3 -2d0*y2*y4**2*y5 &
               +2d0*y1*y4*y5**2*SQRT3/3.D0 -2d0*y1*y4**3*SQRT3 &
               +4.D0/3.D0*2d0*y2*y4*y5**2*SQRT3 )*FEA233445 &
               +( -4d0*y1**3*y5 -2.D0*4d0*y2**3*y5 -4d0*y1**3*y4*SQRT3 )*FEA233335 +( &
               +2d0*y1*4d0*y2**3 +4d0*y1**3*2d0*y2 )*FEA222233 
 
            ddv(1,2)= &
               +ddv1(1,2) +ddv2(1,2) +ddv3(1,2) +ddv4(1,2) +ddv5(1,2) +ddv6(1,2)
 
            ddv1(1,3)=0d0 
            ddv2(1,3)=( +1d0 )*FEA12 +0d0 +0d0 &
               +0d0 
            ddv3(1,3)=( y4 -SQRT3*y5 )*FEA124 +0d0 +( +2d0*y3 &
               +2d0*y1 )*FEA112 +0d0 +y2*FEA123 +0d0 +0d0 +0d0 
            dds2(1,3)=0d0 +0d0 &
               +0d0 
            dds1(1,3)=dds2(1,3) +0d0 +0d0 +0d0 +( SQRT3*y4*y5 &
               +y4**2 )*FEA1244 
            dds2(1,3)=dds1(1,3) +( y5**2 -SQRT3*y4*y5 )*FEA1255 &
               +( -2d0*y3*y4/2.D0 +2d0*y1*y4 -SQRT3*2d0*y3*y5/2.D0 )*FEA1124 +( &
               +SQRT3*2d0*y3*y4/2.D0 -2d0*y3*y5/2.D0 &
               -2d0*y1*y5 )*FEA1125 
            ddv4(1,3)=dds2(1,3) &
               +( +3d0*y1**2 +3d0*y3**2 )*FEA1112 +( +2d0*y1*2d0*y3 )*FEA1122 &
               +( y2**2 +2d0*y1*y2 +y2*2d0*y3 )*FEA1123 +0d0 
            dds3(1,3)=0d0 +0d0 &
               +0d0 
            dds2(1,3)=dds3(1,3) +0d0 +0d0 +0d0 
            dds1(1,3)=dds2(1,3) +0d0 +0d0 &
               +0d0 +0d0 +( +y4**2*y5 -4.D0/9.D0*y4**3*SQRT3 +y5**3 )*FEA13445 &
               +( -y4**3/3.D0 +y4*y5**2 )*FEA13455 
            dds3(1,3)=dds1(1,3) &
               +( 2d0*y1*y5**2 +2d0*y3*y4**2 +2d0*y1*y4**2 +2d0*y3*y5**2 )*FEA11255 &
               +( 2.D0/3.D0*2d0*y1*y4**2*SQRT3 +2d0*y3*y5**2*SQRT3/2.D0 +2d0*y3*y4*y5 &
               +2d0*y3*y4**2*SQRT3/6.D0 )*FEA13345 
            dds4(1,3)=dds3(1,3) +( &
               +2d0*y1*y4**2*SQRT3/3.D0 -2d0*y1*y4*y5 -2d0*y3*y4**2*SQRT3/6.D0 &
               +2d0*y3*y5**2*SQRT3/2.D0 )*FEA11245 
            dds2(1,3)=dds4(1,3) +( &
               +3d0*y1**2*y5 +3d0*y3**2*y5/2.D0 -3d0*y3**2*y4*SQRT3/2.D0 )*FEA11135 &
               +( 3d0*y1**2*y4 -3d0*y3**2*y4/2.D0 &
               -3d0*y3**2*y5*SQRT3/2.D0 )*FEA11134 
            ddv5(1,3)=dds2(1,3) &
               +( +4d0*y1**3 +4d0*y3**3 )*FEA23333 +( -SQRT3*2d0*y1*2d0*y3*y5 &
               +2d0*y1*2d0*y3*y4 )*FEA11334 +( 2d0*y1*3d0*y3**2 &
               +3d0*y1**2*2d0*y3 )*FEA11333 +( y2*y4**2 +y2*y5**2 )*FEA12355 &
               +( -y2*2d0*y3*y4/2.D0 -y2**2*y4/2.D0 -SQRT3*y2*2d0*y3*y5/2.D0 &
               +2d0*y1*y2*y4 +SQRT3*y2**2*y5/2.D0 )*FEA11234 +( y2**3 +y2*3d0*y3**2 &
               +3d0*y1**2*y2 )*FEA11123 +( 2d0*y1*y2**2 +y2**2*2d0*y3 &
               +2d0*y1*y2*2d0*y3 )*FEA11233 
            dds3(1,3)=0d0 +0d0 +( 5d0*y3**4 &
               +5d0*y1**4 )*FEA133333 +( 4d0*y1**3*y4 &
               +4d0*y3**3*y4 -4d0*y3**3*y5*SQRT3 -4d0*y1**3*y5*SQRT3 )*FEA133334 &
               +( -y2*y4**3/3.D0 +y2*y4*y5**2 )*FEA123455 
            dds4(1,3)=dds3(1,3) &
               +( 2.D0/3.D0*SQRT3*y2**2*2d0*y3*y4 -2d0*y1*y2**2*y5 &
               -SQRT3*2d0*y1*y2**2*y4/3.D0 &
               +2d0*y1*y2*2d0*y3*y5 -SQRT3*2d0*y1*y2*2d0*y3*y4/3.D0 )*FEA112335 &
               +( y2**2*y5**2 +y2*2d0*y3*y5**2 +y2*2d0*y3*y4**2 +y2**2*y4**2 &
               +2d0*y1*y2*y4**2 +2d0*y1*y2*y5**2 )*FEA112355 
            dds2(1,3)=dds4(1,3) &
               +( -2d0*y1*3d0*y3**2*y5/2.D0 +3d0*y1**2*2d0*y3*y5/2.D0 &
               +3d0*y1**2*2d0*y3*y4*SQRT3/2.D0 &
               -2d0*y1*3d0*y3**2*y4*SQRT3/2.D0 )*FEA222335 &
               +( -2d0*y1*2d0*y3*y5**2*SQRT3/2.D0 &
               +2d0*y1*2d0*y3*y4*y5 -2d0*y1*2d0*y3*y4**2*SQRT3/6.D0 )*FEA113345 +( &
               +2d0*y1*2d0*y3*y4**2 &
               +2d0*y1*2d0*y3*y5**2 )*FEA223355 
            dds3(1,3)=dds2(1,3) &
               +( y2*2d0*y3*y4**2*SQRT3/6.D0 +y2*2d0*y3*y4*y5 &
               +y2*2d0*y3*y5**2*SQRT3/2.D0 &
               +2.D0/3.D0*2d0*y1*y2*y4**2*SQRT3 -y2**2*y4*y5 +y2**2*y4**2*SQRT3/6.D0 &
               +y2**2*y5**2*SQRT3/2.D0 )*FEA123345 +( -3d0*y1**2*2d0*y3*y4/2.D0 &
               +3d0*y1**2*2d0*y3*y5*SQRT3/2.D0 -2d0*y1*3d0*y3**2*y4/2.D0 &
               +2d0*y1*3d0*y3**2*y5*SQRT3/2.D0 )*FEA222334 +0d0 +( &
               +3d0*y1**2*3d0*y3**2 )*FEA222333 
            dds4(1,3)=dds3(1,3) +0d0 &
               +0d0 
            dds5(1,3)=dds4(1,3) +0d0 
            dds1(1,3)=dds5(1,3) &
               +( -SQRT3*y2*3d0*y3**2*y5/2.D0 +3d0*y1**2*y2*y4 &
               +SQRT3*y2**3*y5/2.D0 -y2**3*y4/2.D0 -y2*3d0*y3**2*y4/2.D0 )*FEA111234 &
               +0d0 +( -y5**4 +y4**2*y5**2 &
               +2.D0*y4**3*y5*SQRT3 )*FEA124455 
            dds3(1,3)=dds1(1,3) +0d0 +( y2**4 &
               +4d0*y1**3*y2 +y2*4d0*y3**3 )*FEA111123 +2d0*y1*y2**2*2d0*y3*FEA112233 &
               +0d0 
            dds4(1,3)=dds3(1,3) +( 3.D0*y5**4 &
               +y4**4 -4.D0*y4**3*y5*SQRT3 )*FEA134444 &
               +( -2d0*y3*y5**3*SQRT3/3.D0 -7.D0/3.D0*2d0*y1*y4*y5**2 &
               -16.D0/3.D0*2d0*y3*y4*y5**2 &
               +4.D0/3.D0*2d0*y1*y4**2*y5*SQRT3 -5.D0/3.D0*2d0*y3*y4**2*y5*SQRT3 &
               +3.D0*2d0*y1*y4**3 )*FEA233444 
            dds5(1,3)=dds4(1,3) &
               +( 2.D0*2d0*y3*y5**3 +2d0*y1*y4*y5**2*SQRT3 &
               +6.D0*2d0*y3*y4**2*y5 -3.D0*2d0*y1*y4**2*y5 &
               +4.D0*2d0*y3*y4*y5**2*SQRT3 &
               +2d0*y1*y5**3 -3.D0*2d0*y1*y4**3*SQRT3 )*FEA113555 
            dds2(1,3)=dds5(1,3) &
               +0d0 +( -3.D0*y4**3*y5 -y4*y5**3 +2.D0/3.D0*y5**4*SQRT3 )*FEA124555 &
               +0d0 
            dds3(1,3)=dds2(1,3) +0d0 +( +3d0*y3**2*y4**2 +3d0*y1**2*y5**2 &
               +3d0*y1**2*y4**2 +3d0*y3**2*y5**2 )*FEA233344 &
               +( -3d0*y1**2*y4*y5 -3d0*y1**2*y5**2*SQRT3/3.D0 &
               -3d0*y3**2*y4**2*SQRT3/2.D0 &
               +3d0*y3**2*y5**2*SQRT3/6.D0 )*FEA233345 +0d0 +( y2**3*2d0*y3 &
               +3d0*y1**2*y2**2 +2d0*y1*y2**3 +y2**2*3d0*y3**2 +2d0*y1*y2*3d0*y3**2 &
               +3d0*y1**2*y2*2d0*y3 )*FEA111233 
            dds4(1,3)=dds3(1,3) +0d0 &
               +( -4.D0/3.D0*2d0*y1*y4**2*y5*SQRT3 -2.D0*2d0*y1*y4**3 &
               -2.D0/3.D0*2d0*y1*y4*y5**2 &
               +2d0*y3*y4**3 +2d0*y3*y5**3*SQRT3/3.D0 &
               +5.D0/3.D0*2d0*y3*y4**2*y5*SQRT3 &
               +7.D0/3.D0*2d0*y3*y4*y5**2 )*FEA133444 
            dds5(1,3)=dds4(1,3) +( &
               +3d0*y3**2*y4**2*SQRT3/2.D0 +3d0*y1**2*y4**2*SQRT3/2.D0 &
               +3d0*y1**2*y5**2*SQRT3/6.D0 +3d0*y1**2*y4*y5 &
               +3d0*y3**2*y5**2*SQRT3/6.D0 &
               +3d0*y3**2*y4*y5 )*FEA133345 
            ddv6(1,3)=dds5(1,3) +( &
               +2d0*y1*y4*y5**2*SQRT3/3.D0 +4.D0/3.D0*2d0*y3*y4*y5**2*SQRT3 &
               +2d0*y3*y5**3 +2d0*y3*y4**2*y5 -2d0*y1*y4**3*SQRT3 )*FEA233445 &
               +( -4d0*y1**3*y4*SQRT3 +2.D0*4d0*y3**3*y5 +4d0*y1**3*y5 )*FEA233335 +( &
               +4d0*y1**3*2d0*y3 +2d0*y1*4d0*y3**3 )*FEA222233 
 
            ddv(1,3)= &
               +ddv1(1,3) +ddv2(1,3) +ddv3(1,3) +ddv4(1,3) +ddv5(1,3) +ddv6(1,3)
 
            ddv1(1,4)=0d0 
            ddv2(1,4)=0d0 +0d0 +( +1d0 )*FEA14 +0d0 
            ddv3(1,4)=( y3 &
               +y2 )*FEA124 +0d0 +0d0 +0d0 +( 2d0*y4 )*FEA144 +0d0 +( &
               +2d0*y1 )*FEA114 
            dds2(1,4)=0d0 +( +3d0*y4**2 )*FEA1444 &
               +0d0 
            dds1(1,4)=dds2(1,4) +( +2d0*y1*2d0*y4 )*FEA1144 &
               +( 3d0*y1**2 )*FEA1114 +0d0 +( SQRT3*y3*y5 +y2*2d0*y4 -SQRT3*y2*y5 &
               +y3*2d0*y4 )*FEA1244 
            dds2(1,4)=dds1(1,4) +( -SQRT3*y3*y5 &
               +SQRT3*y2*y5 )*FEA1255 +( -y3**2*1d0/2.D0 +2d0*y1*y3 &
               +2d0*y1*y2 -y2**2*1d0/2.D0 )*FEA1124 +( +SQRT3*y3**2*1d0/2.D0 &
               +SQRT3*y2**2*1d0/2.D0 )*FEA1125 
            ddv4(1,4)=dds2(1,4) +0d0 +0d0 +0d0 +( &
               +y5**2 )*FEA1455 
            dds3(1,4)=0d0 +( &
               +9.D0*2d0*y4*y5**2 -3.D0/2.D0*4d0*y4**3 )*FEA25555 &
               +( -7.D0/2.D0*2d0*y4*y5**2 &
               +5.D0/4.D0*4d0*y4**3 )*FEA24455 
            dds2(1,4)=dds3(1,4) &
               +( -7.D0/12.D0*4d0*y4**3*SQRT3 &
               +3.D0/2.D0*2d0*y4*y5**2*SQRT3 )*FEA24445 &
               +( -5.D0/9.D0*2d0*y1*3d0*y4**2*SQRT3 -2d0*y1*y5**2*SQRT3 )*FEA33445 &
               +( -2d0*y1*3d0*y4**2/3.D0 +2d0*y1*y5**2 )*FEA33455 
            dds1(1,4)=dds2(1,4) &
               +( +3d0*y1**2*2d0*y4*SQRT3/2.D0 )*FEA33345 +( &
               +3d0*y1**2*2d0*y4 )*FEA33344 +( -2.D0*4d0*y1**3 )*FEA33334 +0d0 &
               +( -4.D0/9.D0*y2*3d0*y4**2*SQRT3 &
               +y3*2d0*y4*y5 -y2*2d0*y4*y5 -4.D0/9.D0*y3*3d0*y4**2*SQRT3 )*FEA13445 &
               +( +y2*y5**2 -y2*3d0*y4**2/3.D0 -y3*3d0*y4**2/3.D0 &
               +y3*y5**2 )*FEA13455 
            dds3(1,4)=dds1(1,4) +( +y2**2*2d0*y4 &
               +y3**2*2d0*y4 +2d0*y1*y3*2d0*y4 +2d0*y1*y2*2d0*y4 )*FEA11255 &
               +( 2.D0/3.D0*2d0*y1*y3*2d0*y4*SQRT3 -y2**2*y5 +y3**2*y5 &
               +y3**2*2d0*y4*SQRT3/6.D0 +y2**2*2d0*y4*SQRT3/6.D0 &
               +2.D0/3.D0*2d0*y1*y2*2d0*y4*SQRT3 )*FEA13345 
            dds4(1,4)=dds3(1,4) &
               +( 2d0*y1*y2*y5 +2d0*y1*y3*2d0*y4*SQRT3/3.D0 &
               +2d0*y1*y2*2d0*y4*SQRT3/3.D0 -y2**2*2d0*y4*SQRT3/6.D0 -2d0*y1*y3*y5 &
               -y3**2*2d0*y4*SQRT3/6.D0 )*FEA11245 
            dds2(1,4)=dds4(1,4) &
               +( -y2**3*SQRT3/2.D0 -y3**3*SQRT3/2.D0 )*FEA11135 +( 3d0*y1**2*y3 &
               +3d0*y1**2*y2 -y3**3*1d0/2.D0 &
               -y2**3*1d0/2.D0 )*FEA11134 
            ddv5(1,4)=dds2(1,4) &
               +0d0 +( +2d0*y1*y2**2 +2d0*y1*y3**2 )*FEA11334 +0d0 &
               +( y2*y3*2d0*y4 )*FEA12355 +( -y2*y3**2*1d0/2.D0 -y2**2*y3*1d0/2.D0 &
               +2d0*y1*y2*y3 )*FEA11234 +0d0 &
               +0d0 
            dds3(1,4)=( -8.D0/3.D0*3d0*y1**2*y5**2*SQRT3 )*FEA333555 &
               +( -4d0*y1**3*2d0*y4*SQRT3/6.D0 )*FEA222245 +0d0 +( 4d0*y1**3*y3 &
               +4d0*y1**3*y2 +y3**4 +y2**4 )*FEA133334 +( -y2*y3*3d0*y4**2/3.D0 &
               +y2*y3*y5**2 )*FEA123455 
            dds4(1,4)=dds3(1,4) &
               +( 2.D0/3.D0*SQRT3*y2**2*y3**2 -SQRT3*2d0*y1*y2**2*y3*1d0/3.D0 &
               -SQRT3*2d0*y1*y2*y3**2*1d0/3.D0 )*FEA112335 &
               +( +y2*y3**2*2d0*y4 +y2**2*y3*2d0*y4 &
               +2d0*y1*y2*y3*2d0*y4 )*FEA112355 
            dds2(1,4)=dds4(1,4) +( &
               +3d0*y1**2*y2**2*SQRT3/2.D0 -2d0*y1*y2**3*SQRT3/2.D0 &
               +3d0*y1**2*y3**2*SQRT3/2.D0 -2d0*y1*y3**3*SQRT3/2.D0 )*FEA222335 &
               +( -2d0*y1*y2**2*2d0*y4*SQRT3/6.D0 -2d0*y1*y2**2*y5 &
               +2d0*y1*y3**2*y5 -2d0*y1*y3**2*2d0*y4*SQRT3/6.D0 )*FEA113345 +( &
               +2d0*y1*y3**2*2d0*y4 &
               +2d0*y1*y2**2*2d0*y4 )*FEA223355 
            dds3(1,4)=dds2(1,4) &
               +( y2*y3**2*2d0*y4*SQRT3/6.D0 +y2*y3**2*y5 &
               +2.D0/3.D0*2d0*y1*y2*y3*2d0*y4*SQRT3 -y2**2*y3*y5 &
               +y2**2*y3*2d0*y4*SQRT3/6.D0 )*FEA123345 &
               +( -3d0*y1**2*y2**2*1d0/2.D0 -3d0*y1**2*y3**2*1d0/2.D0 &
               -2d0*y1*y2**3*1d0/2.D0 -2d0*y1*y3**3*1d0/2.D0 )*FEA222334 &
               +( +9.D0*2d0*y1*2d0*y4*y5**2 -3.D0/2.D0*2d0*y1*4d0*y4**3 )*FEA335555 &
               +0d0 
            dds4(1,4)=dds3(1,4) &
               +( -2.D0/5.D0*5d0*y4**4 -2.D0*3d0*y4**2*y5**2 )*FEA244455 &
               +( -2.D0*5d0*y1**4 )*FEA222224 
            dds5(1,4)=dds4(1,4) +( +y5**4 &
               +2.D0*3d0*y4**2*y5**2 -5d0*y4**4/15.D0 )*FEA145555 
            dds1(1,4)=dds5(1,4) &
               +( +3d0*y1**2*y2*y3 -y2**3*y3*1d0/2.D0 -y2*y3**3*1d0/2.D0 )*FEA111234 &
               +( &
               +2.D0/9.D0*5d0*y4**4*SQRT3 -2.D0/3.D0*3d0*y4**2*y5**2*SQRT3 )*FEA244555 &
               +( y2*2d0*y4*y5**2 +y3*2d0*y4*y5**2 -2.D0*y2*3d0*y4**2*y5*SQRT3 &
               +2.D0*y3*3d0*y4**2*y5*SQRT3 )*FEA124455 
            dds3(1,4)=dds1(1,4) +0d0 +0d0 &
               +( 4d0*y1**3*2d0*y4 )*FEA222255 
            dds4(1,4)=dds3(1,4) +( &
               +y3*4d0*y4**3 -4.D0*y3*3d0*y4**2*y5*SQRT3 +y2*4d0*y4**3 &
               +4.D0*y2*3d0*y4**2*y5*SQRT3 )*FEA134444 +( -7.D0/3.D0*2d0*y1*y3*y5**2 &
               +5.D0/3.D0*y2**2*2d0*y4*y5*SQRT3 -7.D0/3.D0*2d0*y1*y2*y5**2 &
               -16.D0/3.D0*y3**2*y5**2 &
               +4.D0/3.D0*2d0*y1*y3*2d0*y4*y5*SQRT3 &
               +3.D0*2d0*y1*y2*3d0*y4**2 -5.D0/3.D0*y3**2*2d0*y4*y5*SQRT3 &
               -4.D0/3.D0*2d0*y1*y2*2d0*y4*y5*SQRT3 &
               +3.D0*2d0*y1*y3*3d0*y4**2 &
               -16.D0/3.D0*y2**2*y5**2 )*FEA233444 
            dds5(1,4)=dds4(1,4) &
               +( +2d0*y1*y3*y5**2*SQRT3 &
               +6.D0*y3**2*2d0*y4*y5 -6.D0*y2**2*2d0*y4*y5 -3.D0*2d0*y1*y3*2d0*y4*y5 &
               +2d0*y1*y2*y5**2*SQRT3 &
               +4.D0*y3**2*y5**2*SQRT3 -3.D0*2d0*y1*y2*3d0*y4**2*SQRT3 &
               +3.D0*2d0*y1*y2*2d0*y4*y5 -3.D0*2d0*y1*y3*3d0*y4**2*SQRT3 &
               +4.D0*y2**2*y5**2*SQRT3 )*FEA113555 
            dds2(1,4)=dds5(1,4) &
               +( -3.D0/2.D0*2d0*y1*2d0*y4*y5**2*SQRT3 &
               +7.D0/12.D0*2d0*y1*4d0*y4**3*SQRT3 )*FEA334445 &
               +( -3.D0*y3*3d0*y4**2*y5 -y3*y5**3 +3.D0*y2*3d0*y4**2*y5 &
               +y2*y5**3 )*FEA124555 +( -7.D0/2.D0*2d0*y1*2d0*y4*y5**2 &
               +5.D0/4.D0*2d0*y1*4d0*y4**3 )*FEA334455 
            dds3(1,4)=dds2(1,4) +0d0 +( &
               +y3**3*2d0*y4 +y2**3*2d0*y4 +3d0*y1**2*y2*2d0*y4 &
               +3d0*y1**2*y3*2d0*y4 )*FEA233344 +( &
               +3d0*y1**2*y2*y5 -3d0*y1**2*y3*y5 -y3**3*2d0*y4*SQRT3/2.D0 &
               -y2**3*2d0*y4*SQRT3/2.D0 )*FEA233345 &
               +( -3.D0*3d0*y1**2*y5**2 +3d0*y1**2*3d0*y4**2 )*FEA111444 &
               +0d0 
            dds4(1,4)=dds3(1,4) +0d0 +( -5.D0/3.D0*y2**2*2d0*y4*y5*SQRT3 &
               +y2**2*3d0*y4**2 -4.D0/3.D0*2d0*y1*y3*2d0*y4*y5*SQRT3 &
               -2.D0*2d0*y1*y2*3d0*y4**2 -2.D0*2d0*y1*y3*3d0*y4**2 &
               +7.D0/3.D0*y2**2*y5**2 -2.D0/3.D0*2d0*y1*y3*y5**2 +y3**2*3d0*y4**2 &
               +4.D0/3.D0*2d0*y1*y2*2d0*y4*y5*SQRT3 &
               +5.D0/3.D0*y3**2*2d0*y4*y5*SQRT3 -2.D0/3.D0*2d0*y1*y2*y5**2 &
               +7.D0/3.D0*y3**2*y5**2 )*FEA133444 
            dds5(1,4)=dds4(1,4) &
               +( -3d0*y1**2*y2*y5 +y3**3*2d0*y4*SQRT3/2.D0 &
               +3d0*y1**2*y3*2d0*y4*SQRT3/2.D0 +3d0*y1**2*y3*y5 &
               +3d0*y1**2*y2*2d0*y4*SQRT3/2.D0 -y2**3*y5 +y2**3*2d0*y4*SQRT3/2.D0 &
               +y3**3*y5 )*FEA133345 
            ddv6(1,4)=dds5(1,4) +( &
               +2d0*y1*y3*y5**2*SQRT3/3.D0 -y2**2*2d0*y4*y5 &
               +4.D0/3.D0*y3**2*y5**2*SQRT3 &
               +2d0*y1*y2*y5**2*SQRT3/3.D0 -2d0*y1*y2*3d0*y4**2*SQRT3 &
               +y3**2*2d0*y4*y5 -2d0*y1*y3*3d0*y4**2*SQRT3 &
               +4.D0/3.D0*y2**2*y5**2*SQRT3 )*FEA233445 &
               +( -4d0*y1**3*y3*SQRT3 -4d0*y1**3*y2*SQRT3 )*FEA233335 &
               +0d0 
 
            ddv(1,4)= +ddv1(1,4) +ddv2(1,4) +ddv3(1,4) +ddv4(1,4) &
               +ddv5(1,4) +ddv6(1,4)
 
            ddv1(1,5)=0d0 
            ddv2(1,5)=0d0 +0d0 +0d0 +0d0 
            ddv3(1,5)=( &
               +SQRT3*y2 -SQRT3*y3 )*FEA124 +( +2d0*y5 )*FEA155 +0d0 +0d0 +0d0 +0d0 &
               +0d0 
            dds2(1,5)=0d0 +0d0 +( &
               +2d0*y1*2d0*y5 )*FEA1155 
            dds1(1,5)=dds2(1,5) +0d0 +0d0 +0d0 &
               +( SQRT3*y3*y4 -SQRT3*y2*y4 )*FEA1244 
            dds2(1,5)=dds1(1,5) +( y3*2d0*y5 &
               +y2*2d0*y5 -SQRT3*y3*y4 +SQRT3*y2*y4 )*FEA1255 &
               +( -SQRT3*y3**2*1d0/2.D0 +SQRT3*y2**2*1d0/2.D0 )*FEA1124 &
               +( 2d0*y1*y2 -y3**2*1d0/2.D0 &
               +y2**2*1d0/2.D0 -2d0*y1*y3 )*FEA1125 
            ddv4(1,5)=dds2(1,5) +0d0 +0d0 &
               +0d0 +( +y4*2d0*y5 )*FEA1455 
            dds3(1,5)=0d0 +( +9.D0*y4**2*2d0*y5 &
               +5.D0/2.D0*4d0*y5**3 )*FEA25555 &
               +( -7.D0/2.D0*y4**2*2d0*y5 &
               -3.D0/4.D0*4d0*y5**3 )*FEA24455 
            dds2(1,5)=dds3(1,5) &
               +( +3.D0/4.D0*4d0*y5**3*SQRT3 +3.D0/2.D0*y4**2*2d0*y5*SQRT3 )*FEA24445 &
               +( -2d0*y1*y4*2d0*y5*SQRT3 )*FEA33445 +( &
               +2d0*y1*y4*2d0*y5 )*FEA33455 
            dds1(1,5)=dds2(1,5) &
               +( -3d0*y1**2*2d0*y5*SQRT3/6.D0 )*FEA33345 +( &
               +3d0*y1**2*2d0*y5 )*FEA33344 +0d0 +0d0 +( -y2*3d0*y5**2 &
               +y3*y4**2 -y2*y4**2 +y3*3d0*y5**2 )*FEA13445 +( +y2*y4*2d0*y5 &
               +y3*y4*2d0*y5 )*FEA13455 
            dds3(1,5)=dds1(1,5) +( 2d0*y1*y3*2d0*y5 &
               +y2**2*2d0*y5 +2d0*y1*y2*2d0*y5 +y3**2*2d0*y5 )*FEA11255 +( &
               +y3**2*2d0*y5*SQRT3/2.D0 +y2**2*2d0*y5*SQRT3/2.D0 -y2**2*y4 &
               +y3**2*y4 )*FEA13345 
            dds4(1,5)=dds3(1,5) +( 2d0*y1*y2*y4 -2d0*y1*y3*y4 &
               +y2**2*2d0*y5*SQRT3/2.D0 &
               +y3**2*2d0*y5*SQRT3/2.D0 )*FEA11245 
            dds2(1,5)=dds4(1,5) &
               +( -3d0*y1**2*y2 +3d0*y1**2*y3 -y2**3*1d0/2.D0 &
               +y3**3*1d0/2.D0 )*FEA11135 +( &
               +y2**3*SQRT3/2.D0 -y3**3*SQRT3/2.D0 )*FEA11134 
            ddv5(1,5)=dds2(1,5) &
               +0d0 +( -SQRT3*2d0*y1*y3**2 +SQRT3*2d0*y1*y2**2 )*FEA11334 +0d0 +( &
               +y2*y3*2d0*y5 )*FEA12355 +( -SQRT3*y2*y3**2*1d0/2.D0 &
               +SQRT3*y2**2*y3*1d0/2.D0 )*FEA11234 +0d0 &
               +0d0 
            dds3(1,5)=( -8.D0/3.D0*3d0*y1**2*y4*2d0*y5*SQRT3 )*FEA333555 &
               +( 4d0*y1**3*2d0*y5*SQRT3/2.D0 )*FEA222245 +0d0 +( +y2**4*SQRT3 &
               +4d0*y1**3*y2*SQRT3 -y3**4*SQRT3 -4d0*y1**3*y3*SQRT3 )*FEA133334 +( &
               +y2*y3*y4*2d0*y5 )*FEA123455 
            dds4(1,5)=dds3(1,5) +( -2d0*y1*y2**2*y3 &
               +2d0*y1*y2*y3**2 )*FEA112335 +( y2**2*y3*2d0*y5 +y2*y3**2*2d0*y5 &
               +2d0*y1*y2*y3*2d0*y5 )*FEA112355 
            dds2(1,5)=dds4(1,5) &
               +( -3d0*y1**2*y2**2*1d0/2.D0 -2d0*y1*y3**3*1d0/2.D0 &
               +3d0*y1**2*y3**2*1d0/2.D0 +2d0*y1*y2**3*1d0/2.D0 )*FEA222335 &
               +( -2d0*y1*y2**2*2d0*y5*SQRT3/2.D0 -2d0*y1*y3**2*2d0*y5*SQRT3/2.D0 &
               -2d0*y1*y2**2*y4 &
               +2d0*y1*y3**2*y4 )*FEA113345 +( +2d0*y1*y2**2*2d0*y5 &
               +2d0*y1*y3**2*2d0*y5 )*FEA223355 
            dds3(1,5)=dds2(1,5) +( +y2*y3**2*y4 &
               +y2*y3**2*2d0*y5*SQRT3/2.D0 -y2**2*y3*y4 &
               +y2**2*y3*2d0*y5*SQRT3/2.D0 )*FEA123345 +( -3d0*y1**2*y2**2*SQRT3/2.D0 &
               +3d0*y1**2*y3**2*SQRT3/2.D0 -2d0*y1*y2**3*SQRT3/2.D0 &
               +2d0*y1*y3**3*SQRT3/2.D0 )*FEA222334 +( +5.D0/2.D0*2d0*y1*4d0*y5**3 &
               +9.D0*2d0*y1*y4**2*2d0*y5 )*FEA335555 +0d0 
            dds4(1,5)=dds3(1,5) &
               +( -2.D0*y4**3*2d0*y5 )*FEA244455 +0d0 
            dds5(1,5)=dds4(1,5) +( &
               +y4*4d0*y5**3 +2.D0*y4**3*2d0*y5 )*FEA145555 
            dds1(1,5)=dds5(1,5) &
               +( -SQRT3*y2*y3**3*1d0/2.D0 +SQRT3*y2**3*y3*1d0/2.D0 )*FEA111234 &
               +( -2.D0/3.D0*y4**3*2d0*y5*SQRT3 )*FEA244555 &
               +( y2*y4**2*2d0*y5 -y2*4d0*y5**3 -y3*4d0*y5**3 &
               +y3*y4**2*2d0*y5 -2.D0*y2*y4**3*SQRT3 &
               +2.D0*y3*y4**3*SQRT3 )*FEA124455 
            dds3(1,5)=dds1(1,5) +0d0 +0d0 +( &
               +4d0*y1**3*2d0*y5 )*FEA222255 
            dds4(1,5)=dds3(1,5) &
               +( 3.D0*y3*4d0*y5**3 -4.D0*y3*y4**3*SQRT3 +4.D0*y2*y4**3*SQRT3 &
               +3.D0*y2*4d0*y5**3 )*FEA134444 &
               +( -y3**2*3d0*y5**2*SQRT3/3.D0 -7.D0/3.D0*2d0*y1*y3*y4*2d0*y5 &
               +5.D0/3.D0*y2**2*y4**2*SQRT3 -7.D0/3.D0*2d0*y1*y2*y4*2d0*y5 &
               -16.D0/3.D0*y3**2*y4*2d0*y5 &
               +4.D0/3.D0*2d0*y1*y3*y4**2*SQRT3 &
               +y2**2*3d0*y5**2*SQRT3/3.D0 -5.D0/3.D0*y3**2*y4**2*SQRT3 &
               -4.D0/3.D0*2d0*y1*y2*y4**2*SQRT3 &
               -16.D0/3.D0*y2**2*y4*2d0*y5 )*FEA233444 
            dds5(1,5)=dds4(1,5) &
               +( 2.D0*y3**2*3d0*y5**2 -2.D0*y2**2*3d0*y5**2 &
               +2d0*y1*y3*y4*2d0*y5*SQRT3 &
               +6.D0*y3**2*y4**2 -6.D0*y2**2*y4**2 -3.D0*2d0*y1*y3*y4**2 &
               +2d0*y1*y2*y4*2d0*y5*SQRT3 +4.D0*y3**2*y4*2d0*y5*SQRT3 &
               +3.D0*2d0*y1*y2*y4**2 -2d0*y1*y2*3d0*y5**2 +2d0*y1*y3*3d0*y5**2 &
               +4.D0*y2**2*y4*2d0*y5*SQRT3 )*FEA113555 
            dds2(1,5)=dds5(1,5) &
               +( -3.D0/2.D0*2d0*y1*y4**2*2d0*y5*SQRT3 &
               -3.D0/4.D0*2d0*y1*4d0*y5**3*SQRT3 )*FEA334445 &
               +( -3.D0*y3*y4**3 +2.D0/3.D0*y2*4d0*y5**3*SQRT3 -y3*y4*3d0*y5**2 &
               +2.D0/3.D0*y3*4d0*y5**3*SQRT3 +3.D0*y2*y4**3 &
               +y2*y4*3d0*y5**2 )*FEA124555 &
               +( -7.D0/2.D0*2d0*y1*y4**2*2d0*y5 &
               -3.D0/4.D0*2d0*y1*4d0*y5**3 )*FEA334455 
            dds3(1,5)=dds2(1,5) &
               +0d0 +( +y2**3*2d0*y5 +3d0*y1**2*y3*2d0*y5 +3d0*y1**2*y2*2d0*y5 &
               +y3**3*2d0*y5 )*FEA233344 +( y2**3*2d0*y5*SQRT3/6.D0 &
               +3d0*y1**2*y2*y4 -3d0*y1**2*y2*2d0*y5*SQRT3/3.D0 -3d0*y1**2*y3*y4 &
               -3d0*y1**2*y3*2d0*y5*SQRT3/3.D0 &
               +y3**3*2d0*y5*SQRT3/6.D0 )*FEA233345 &
               +( -3.D0*3d0*y1**2*y4*2d0*y5 )*FEA111444 +0d0 
            dds4(1,5)=dds3(1,5) +0d0 &
               +( -5.D0/3.D0*y2**2*y4**2*SQRT3 -4.D0/3.D0*2d0*y1*y3*y4**2*SQRT3 &
               -y2**2*3d0*y5**2*SQRT3/3.D0 &
               +7.D0/3.D0*y2**2*y4*2d0*y5 -2.D0/3.D0*2d0*y1*y3*y4*2d0*y5 &
               +y3**2*3d0*y5**2*SQRT3/3.D0 +4.D0/3.D0*2d0*y1*y2*y4**2*SQRT3 &
               +5.D0/3.D0*y3**2*y4**2*SQRT3 -2.D0/3.D0*2d0*y1*y2*y4*2d0*y5 &
               +7.D0/3.D0*y3**2*y4*2d0*y5 )*FEA133444 
            dds5(1,5)=dds4(1,5) &
               +( -3d0*y1**2*y2*y4 +3d0*y1**2*y3*2d0*y5*SQRT3/6.D0 &
               +3d0*y1**2*y2*2d0*y5*SQRT3/6.D0 +3d0*y1**2*y3*y4 &
               +y2**3*2d0*y5*SQRT3/6.D0 -y2**3*y4 +y3**3*2d0*y5*SQRT3/6.D0 &
               +y3**3*y4 )*FEA133345 
            ddv6(1,5)=dds5(1,5) +( &
               +2d0*y1*y3*y4*2d0*y5*SQRT3/3.D0 -y2**2*3d0*y5**2 -y2**2*y4**2 &
               +4.D0/3.D0*y3**2*y4*2d0*y5*SQRT3 +y3**2*3d0*y5**2 &
               +2d0*y1*y2*y4*2d0*y5*SQRT3/3.D0 +y3**2*y4**2 &
               +4.D0/3.D0*y2**2*y4*2d0*y5*SQRT3 )*FEA233445 &
               +( -4d0*y1**3*y2 -2.D0*y2**4 +2.D0*y3**4 +4d0*y1**3*y3 )*FEA233335 &
               +0d0 
 
            ddv(1,5)= +ddv1(1,5) +ddv2(1,5) +ddv3(1,5) +ddv4(1,5) &
               +ddv5(1,5) +ddv6(1,5)
 
            ddv1(1,6)=( +1d0 )*DFEA1 
            ddv2(1,6)=( +y3 +y2 )*DFEA12 +( &
               +2d0*y1 )*DFEA11 +( +y4 )*DFEA14 +0d0 
            ddv3(1,6)=( y3*y4 +y2*y4 &
               +SQRT3*y2*y5 -SQRT3*y3*y5 )*DFEA124 +( +y5**2 )*DFEA155 +( +y3**2 &
               +2d0*y1*y3 +y2**2 +2d0*y1*y2 )*DFEA112 +0d0 +y2*y3*DFEA123 &
               +( y4**2 )*DFEA144 +( +3d0*y1**2 )*DFEA111 +( &
               +2d0*y1*y4 )*DFEA114 
            dds2(1,6)=0d0 +( +y4**3 )*DFEA1444 +( &
               +2d0*y1*y5**2 )*DFEA1155 
            dds1(1,6)=dds2(1,6) +( &
               +2d0*y1*y4**2 )*DFEA1144 +( 3d0*y1**2*y4 )*DFEA1114 +( &
               +4d0*y1**3 )*DFEA1111 +( SQRT3*y3*y4*y5 +y2*y4**2 -SQRT3*y2*y4*y5 &
               +y3*y4**2 )*DFEA1244 
            dds2(1,6)=dds1(1,6) +( y3*y5**2 &
               +y2*y5**2 -SQRT3*y3*y4*y5 +SQRT3*y2*y4*y5 )*DFEA1255 +( -y3**2*y4/2.D0 &
               +2d0*y1*y3*y4 -SQRT3*y3**2*y5/2.D0 +2d0*y1*y2*y4 &
               +SQRT3*y2**2*y5/2.D0 -y2**2*y4/2.D0 )*DFEA1124 +( 2d0*y1*y2*y5 &
               +SQRT3*y3**2*y4/2.D0 +SQRT3*y2**2*y4/2.D0 -y3**2*y5/2.D0 &
               +y2**2*y5/2.D0 -2d0*y1*y3*y5 )*DFEA1125 
            ddv4(1,6)=dds2(1,6) +( &
               +3d0*y1**2*y3 +3d0*y1**2*y2 +y2**3 +y3**3 )*DFEA1112 +( +2d0*y1*y3**2 &
               +2d0*y1*y2**2 )*DFEA1122 +( y2**2*y3 +2d0*y1*y2*y3 &
               +y2*y3**2 )*DFEA1123 +( +y4*y5**2 )*DFEA1455 
            dds3(1,6)=0d0 +( &
               +9.D0*y4**2*y5**2 -3.D0/2.D0*y4**4 +5.D0/2.D0*y5**4 )*DFEA25555 &
               +( -7.D0/2.D0*y4**2*y5**2 -3.D0/4.D0*y5**4 &
               +5.D0/4.D0*y4**4 )*DFEA24455 
            dds2(1,6)=dds3(1,6) +( &
               +3.D0/4.D0*y5**4*SQRT3 -7.D0/12.D0*y4**4*SQRT3 &
               +3.D0/2.D0*y4**2*y5**2*SQRT3 )*DFEA24445 &
               +( -5.D0/9.D0*2d0*y1*y4**3*SQRT3 -2d0*y1*y4*y5**2*SQRT3 )*DFEA33445 &
               +( -2d0*y1*y4**3/3.D0 +2d0*y1*y4*y5**2 )*DFEA33455 
            dds1(1,6)=dds2(1,6) &
               +( +3d0*y1**2*y4**2*SQRT3/2.D0 -3d0*y1**2*y5**2*SQRT3/6.D0 )*DFEA33345 &
               +( +3d0*y1**2*y5**2 +3d0*y1**2*y4**2 )*DFEA33344 &
               +( -2.D0*4d0*y1**3*y4 )*DFEA33334 +( +5d0*y1**4 )*DFEA33333 &
               +( -4.D0/9.D0*y2*y4**3*SQRT3 -y2*y5**3 &
               +y3*y4**2*y5 -y2*y4**2*y5 -4.D0/9.D0*y3*y4**3*SQRT3 &
               +y3*y5**3 )*DFEA13445 +( +y2*y4*y5**2 -y2*y4**3/3.D0 -y3*y4**3/3.D0 &
               +y3*y4*y5**2 )*DFEA13455 
            dds3(1,6)=dds1(1,6) +( 2d0*y1*y3*y5**2 &
               +y2**2*y5**2 +2d0*y1*y2*y5**2 +y2**2*y4**2 +y3**2*y4**2 &
               +2d0*y1*y3*y4**2 +2d0*y1*y2*y4**2 +y3**2*y5**2 )*DFEA11255 &
               +( 2.D0/3.D0*2d0*y1*y3*y4**2*SQRT3 +y3**2*y5**2*SQRT3/2.D0 &
               +y2**2*y5**2*SQRT3/2.D0 -y2**2*y4*y5 +y3**2*y4*y5 &
               +y3**2*y4**2*SQRT3/6.D0 +y2**2*y4**2*SQRT3/6.D0 &
               +2.D0/3.D0*2d0*y1*y2*y4**2*SQRT3 )*DFEA13345 
            dds4(1,6)=dds3(1,6) &
               +( 2d0*y1*y2*y4*y5 +2d0*y1*y3*y4**2*SQRT3/3.D0 &
               +2d0*y1*y2*y4**2*SQRT3/3.D0 -y2**2*y4**2*SQRT3/6.D0 -2d0*y1*y3*y4*y5 &
               +y2**2*y5**2*SQRT3/2.D0 -y3**2*y4**2*SQRT3/6.D0 &
               +y3**2*y5**2*SQRT3/2.D0 )*DFEA11245 
            dds2(1,6)=dds4(1,6) &
               +( -3d0*y1**2*y2*y5 &
               +3d0*y1**2*y3*y5 -y2**3*y4*SQRT3/2.D0 -y2**3*y5/2.D0 &
               +y3**3*y5/2.D0 -y3**3*y4*SQRT3/2.D0 )*DFEA11135 +( 3d0*y1**2*y3*y4 &
               +3d0*y1**2*y2*y4 -y3**3*y4/2.D0 &
               +y2**3*y5*SQRT3/2.D0 -y2**3*y4/2.D0 &
               -y3**3*y5*SQRT3/2.D0 )*DFEA11134 
            ddv5(1,6)=dds2(1,6) &
               +( y2**4 +4d0*y1**3*y3 +4d0*y1**3*y2 +y3**4 )*DFEA23333 +( &
               +2d0*y1*y2**2*y4 -SQRT3*2d0*y1*y3**2*y5 +SQRT3*2d0*y1*y2**2*y5 &
               +2d0*y1*y3**2*y4 )*DFEA11334 +( 2d0*y1*y3**3 +3d0*y1**2*y3**2 &
               +2d0*y1*y2**3 +3d0*y1**2*y2**2 )*DFEA11333 +( y2*y3*y4**2 &
               +y2*y3*y5**2 )*DFEA12355 &
               +( -y2*y3**2*y4/2.D0 -y2**2*y3*y4/2.D0 -SQRT3*y2*y3**2*y5/2.D0 &
               +2d0*y1*y2*y3*y4 +SQRT3*y2**2*y3*y5/2.D0 )*DFEA11234 +( y2**3*y3 &
               +y2*y3**3 +3d0*y1**2*y2*y3 )*DFEA11123 +( 2d0*y1*y2**2*y3 +y2**2*y3**2 &
               +2d0*y1*y2*y3**2 )*DFEA11233 
            dds3(1,6)=( &
               -8.D0/3.D0*3d0*y1**2*y4*y5**2*SQRT3 )*DFEA333555 &
               +( 4d0*y1**3*y5**2*SQRT3/2.D0 -4d0*y1**3*y4**2*SQRT3/6.D0 )*DFEA222245 &
               +( y3**5 +y2**5 +5d0*y1**4*y3 +5d0*y1**4*y2 )*DFEA133333 &
               +( 4d0*y1**3*y3*y4 +4d0*y1**3*y2*y4 +y2**4*y5*SQRT3 +y3**4*y4 &
               +4d0*y1**3*y2*y5*SQRT3 -y3**4*y5*SQRT3 -4d0*y1**3*y3*y5*SQRT3 &
               +y2**4*y4 )*DFEA133334 +( -y2*y3*y4**3/3.D0 &
               +y2*y3*y4*y5**2 )*DFEA123455 
            dds4(1,6)=dds3(1,6) &
               +( 2.D0/3.D0*SQRT3*y2**2*y3**2*y4 -2d0*y1*y2**2*y3*y5 &
               -SQRT3*2d0*y1*y2**2*y3*y4/3.D0 &
               +2d0*y1*y2*y3**2*y5 -SQRT3*2d0*y1*y2*y3**2*y4/3.D0 )*DFEA112335 &
               +( y2**2*y3*y5**2 +y2*y3**2*y5**2 +y2*y3**2*y4**2 +y2**2*y3*y4**2 &
               +2d0*y1*y2*y3*y4**2 &
               +2d0*y1*y2*y3*y5**2 )*DFEA112355 
            dds2(1,6)=dds4(1,6) &
               +( -3d0*y1**2*y2**2*y5/2.D0 -2d0*y1*y3**3*y5/2.D0 &
               +3d0*y1**2*y2**2*y4*SQRT3/2.D0 -2d0*y1*y2**3*y4*SQRT3/2.D0 &
               +3d0*y1**2*y3**2*y5/2.D0 +2d0*y1*y2**3*y5/2.D0 &
               +3d0*y1**2*y3**2*y4*SQRT3/2.D0 -2d0*y1*y3**3*y4*SQRT3/2.D0 )*DFEA222335 &
               +( -2d0*y1*y2**2*y5**2*SQRT3/2.D0 -2d0*y1*y3**2*y5**2*SQRT3/2.D0 &
               -2d0*y1*y2**2*y4**2*SQRT3/6.D0 -2d0*y1*y2**2*y4*y5 &
               +2d0*y1*y3**2*y4*y5 -2d0*y1*y3**2*y4**2*SQRT3/6.D0 )*DFEA113345 +( &
               +2d0*y1*y2**2*y5**2 +2d0*y1*y3**2*y4**2 +2d0*y1*y3**2*y5**2 &
               +2d0*y1*y2**2*y4**2 )*DFEA223355 
            dds3(1,6)=dds2(1,6) &
               +( y2*y3**2*y4**2*SQRT3/6.D0 +y2*y3**2*y4*y5 &
               +y2*y3**2*y5**2*SQRT3/2.D0 &
               +2.D0/3.D0*2d0*y1*y2*y3*y4**2*SQRT3 -y2**2*y3*y4*y5 &
               +y2**2*y3*y4**2*SQRT3/6.D0 +y2**2*y3*y5**2*SQRT3/2.D0 )*DFEA123345 &
               +( -3d0*y1**2*y2**2*y5*SQRT3/2.D0 -3d0*y1**2*y2**2*y4/2.D0 &
               -3d0*y1**2*y3**2*y4/2.D0 -2d0*y1*y2**3*y4/2.D0 &
               +3d0*y1**2*y3**2*y5*SQRT3/2.D0 -2d0*y1*y3**3*y4/2.D0 &
               -2d0*y1*y2**3*y5*SQRT3/2.D0 &
               +2d0*y1*y3**3*y5*SQRT3/2.D0 )*DFEA222334 +( +5.D0/2.D0*2d0*y1*y5**4 &
               +9.D0*2d0*y1*y4**2*y5**2 -3.D0/2.D0*2d0*y1*y4**4 )*DFEA335555 &
               +( 3d0*y1**2*y2**3 +3d0*y1**2*y3**3 )*DFEA222333 
            dds4(1,6)=dds3(1,6) &
               +( -2.D0/5.D0*y4**5 -2.D0*y4**3*y5**2 )*DFEA244455 &
               +( -2.D0*5d0*y1**4*y4 )*DFEA222224 
            dds5(1,6)=dds4(1,6) +( +y4*y5**4 &
               +2.D0*y4**3*y5**2 -y4**5/15.D0 )*DFEA145555 
            dds1(1,6)=dds5(1,6) &
               +( -SQRT3*y2*y3**3*y5/2.D0 +3d0*y1**2*y2*y3*y4 &
               +SQRT3*y2**3*y3*y5/2.D0 -y2**3*y3*y4/2.D0 &
               -y2*y3**3*y4/2.D0 )*DFEA111234 &
               +( +2.D0/9.D0*y4**5*SQRT3 -2.D0/3.D0*y4**3*y5**2*SQRT3 )*DFEA244555 &
               +( y2*y4**2*y5**2 -y2*y5**4 -y3*y5**4 &
               +y3*y4**2*y5**2 -2.D0*y2*y4**3*y5*SQRT3 &
               +2.D0*y3*y4**3*y5*SQRT3 )*DFEA124455 
            dds3(1,6)=dds1(1,6) +( &
               +6d0*y1**5 )*DFEA333333 +( y2**4*y3 +4d0*y1**3*y2*y3 &
               +y2*y3**4 )*DFEA111123 +2d0*y1*y2**2*y3**2*DFEA112233 &
               +( 4d0*y1**3*y4**2 +4d0*y1**3*y5**2 )*DFEA222255 
            dds4(1,6)=dds3(1,6) &
               +( 3.D0*y3*y5**4 +y3*y4**4 -4.D0*y3*y4**3*y5*SQRT3 +y2*y4**4 &
               +4.D0*y2*y4**3*y5*SQRT3 +3.D0*y2*y5**4 )*DFEA134444 &
               +( -y3**2*y5**3*SQRT3/3.D0 -7.D0/3.D0*2d0*y1*y3*y4*y5**2 &
               +5.D0/3.D0*y2**2*y4**2*y5*SQRT3 -7.D0/3.D0*2d0*y1*y2*y4*y5**2 &
               -16.D0/3.D0*y3**2*y4*y5**2 &
               +4.D0/3.D0*2d0*y1*y3*y4**2*y5*SQRT3 +3.D0*2d0*y1*y2*y4**3 &
               +y2**2*y5**3*SQRT3/3.D0 -5.D0/3.D0*y3**2*y4**2*y5*SQRT3 &
               -4.D0/3.D0*2d0*y1*y2*y4**2*y5*SQRT3 &
               +3.D0*2d0*y1*y3*y4**3 &
               -16.D0/3.D0*y2**2*y4*y5**2 )*DFEA233444 
            dds5(1,6)=dds4(1,6) &
               +( 2.D0*y3**2*y5**3 -2.D0*y2**2*y5**3 +2d0*y1*y3*y4*y5**2*SQRT3 &
               +6.D0*y3**2*y4**2*y5 -6.D0*y2**2*y4**2*y5 -3.D0*2d0*y1*y3*y4**2*y5 &
               +2d0*y1*y2*y4*y5**2*SQRT3 &
               +4.D0*y3**2*y4*y5**2*SQRT3 -3.D0*2d0*y1*y2*y4**3*SQRT3 &
               +3.D0*2d0*y1*y2*y4**2*y5 -2d0*y1*y2*y5**3 &
               +2d0*y1*y3*y5**3 -3.D0*2d0*y1*y3*y4**3*SQRT3 &
               +4.D0*y2**2*y4*y5**2*SQRT3 )*DFEA113555 
            dds2(1,6)=dds5(1,6) &
               +( -3.D0/2.D0*2d0*y1*y4**2*y5**2*SQRT3 -3.D0/4.D0*2d0*y1*y5**4*SQRT3 &
               +7.D0/12.D0*2d0*y1*y4**4*SQRT3 )*DFEA334445 +( -3.D0*y3*y4**3*y5 &
               +2.D0/3.D0*y2*y5**4*SQRT3 -y3*y4*y5**3 +2.D0/3.D0*y3*y5**4*SQRT3 &
               +3.D0*y2*y4**3*y5 +y2*y4*y5**3 )*DFEA124555 &
               +( -7.D0/2.D0*2d0*y1*y4**2*y5**2 -3.D0/4.D0*2d0*y1*y5**4 &
               +5.D0/4.D0*2d0*y1*y4**4 )*DFEA334455 
            dds3(1,6)=dds2(1,6) +0d0 +( &
               +y3**3*y4**2 +y2**3*y4**2 +3d0*y1**2*y2*y4**2 +y2**3*y5**2 &
               +3d0*y1**2*y3*y5**2 +3d0*y1**2*y3*y4**2 +3d0*y1**2*y2*y5**2 &
               +y3**3*y5**2 )*DFEA233344 +( y2**3*y5**2*SQRT3/6.D0 &
               +3d0*y1**2*y2*y4*y5 -3d0*y1**2*y2*y5**2*SQRT3/3.D0 -3d0*y1**2*y3*y4*y5 &
               -3d0*y1**2*y3*y5**2*SQRT3/3.D0 -y3**3*y4**2*SQRT3/2.D0 &
               +y3**3*y5**2*SQRT3/6.D0 -y2**3*y4**2*SQRT3/2.D0 )*DFEA233345 &
               +( -3.D0*3d0*y1**2*y4*y5**2 +3d0*y1**2*y4**3 )*DFEA111444 &
               +( y2**3*y3**2 +3d0*y1**2*y2**2*y3 +2d0*y1*y2**3*y3 +y2**2*y3**3 &
               +2d0*y1*y2*y3**3 +3d0*y1**2*y2*y3**2 )*DFEA111233 
            dds4(1,6)=dds3(1,6) &
               +0d0 +( -5.D0/3.D0*y2**2*y4**2*y5*SQRT3 &
               +y2**2*y4**3 -4.D0/3.D0*2d0*y1*y3*y4**2*y5*SQRT3 -2.D0*2d0*y1*y2*y4**3 &
               -y2**2*y5**3*SQRT3/3.D0 -2.D0*2d0*y1*y3*y4**3 &
               +7.D0/3.D0*y2**2*y4*y5**2 -2.D0/3.D0*2d0*y1*y3*y4*y5**2 +y3**2*y4**3 &
               +y3**2*y5**3*SQRT3/3.D0 +4.D0/3.D0*2d0*y1*y2*y4**2*y5*SQRT3 &
               +5.D0/3.D0*y3**2*y4**2*y5*SQRT3 -2.D0/3.D0*2d0*y1*y2*y4*y5**2 &
               +7.D0/3.D0*y3**2*y4*y5**2 )*DFEA133444 
            dds5(1,6)=dds4(1,6) &
               +( -3d0*y1**2*y2*y4*y5 +y3**3*y4**2*SQRT3/2.D0 &
               +3d0*y1**2*y3*y4**2*SQRT3/2.D0 +3d0*y1**2*y3*y5**2*SQRT3/6.D0 &
               +3d0*y1**2*y2*y5**2*SQRT3/6.D0 +3d0*y1**2*y3*y4*y5 &
               +y2**3*y5**2*SQRT3/6.D0 +3d0*y1**2*y2*y4**2*SQRT3/2.D0 -y2**3*y4*y5 &
               +y2**3*y4**2*SQRT3/2.D0 +y3**3*y5**2*SQRT3/6.D0 &
               +y3**3*y4*y5 )*DFEA133345 
            ddv6(1,6)=dds5(1,6) +( &
               +2d0*y1*y3*y4*y5**2*SQRT3/3.D0 -y2**2*y5**3 -y2**2*y4**2*y5 &
               +4.D0/3.D0*y3**2*y4*y5**2*SQRT3 +y3**2*y5**3 &
               +2d0*y1*y2*y4*y5**2*SQRT3/3.D0 -2d0*y1*y2*y4**3*SQRT3 &
               +y3**2*y4**2*y5 -2d0*y1*y3*y4**3*SQRT3 &
               +4.D0/3.D0*y2**2*y4*y5**2*SQRT3 )*DFEA233445 &
               +( -4d0*y1**3*y2*y5 -4d0*y1**3*y3*y4*SQRT3 -2.D0*y2**4*y5 &
               +2.D0*y3**4*y5 -4d0*y1**3*y2*y4*SQRT3 +4d0*y1**3*y3*y5 )*DFEA233335 +( &
               +4d0*y1**3*y3**2 +2d0*y1*y2**4 +2d0*y1*y3**4 &
               +4d0*y1**3*y2**2 )*DFEA222233 
 
            ddv(1,6)= +ddv1(1,6) +ddv2(1,6) &
               +ddv3(1,6) +ddv4(1,6) +ddv5(1,6) +ddv6(1,6)
 
            ddv1(2,2)=0d0 
            ddv2(2,2)=0d0 +( 2d0 )*FEA11 +0d0 +0d0 
            ddv3(2,2)=0d0 &
               +0d0 +( +y1*2d0 +2d0*y3 )*FEA112 +0d0 +0d0 +( +3d0*2d0*y2 )*FEA111 &
               +( -2d0*y4/2.D0 +SQRT3*2d0*y5/2.D0 )*FEA114 
            dds2(2,2)=0d0 +0d0 &
               +( 3.D0/4.D0*2d0*y4**2 +SQRT3*2d0*y4*y5/2.D0 &
               +2d0*y5**2/4.D0 )*FEA1155 
            dds1(2,2)=dds2(2,2) +( &
               +2d0*y4**2/4.D0 -SQRT3*2d0*y4*y5/2.D0 +3.D0/4.D0*2d0*y5**2 )*FEA1144 &
               +( +SQRT3*3d0*2d0*y2*y5/2.D0 -3d0*2d0*y2*y4/2.D0 )*FEA1114 &
               +( 4d0*3d0*y2**2 )*FEA1111 +0d0 
            dds2(2,2)=dds1(2,2) +0d0 +( &
               +SQRT3*2d0*y3*y5/2.D0 -2d0*y3*y4/2.D0 &
               +SQRT3*y1*2d0*y5/2.D0 -y1*2d0*y4/2.D0 )*FEA1124 +( &
               +SQRT3*y1*2d0*y4/2.D0 -SQRT3*2d0*y3*y4/2.D0 -2d0*y3*y5/2.D0 &
               +y1*2d0*y5/2.D0 )*FEA1125 
            ddv4(2,2)=dds2(2,2) +( +y1*3d0*2d0*y2 &
               +3d0*2d0*y2*y3 )*FEA1112 +( 2d0*y3**2 +y1**2*2d0 )*FEA1122 &
               +( y1*2d0*y3 )*FEA1123 +0d0 
            dds3(2,2)=0d0 +0d0 &
               +0d0 
            dds2(2,2)=dds3(2,2) +0d0 +( -2d0*y5**3 &
               +4.D0/9.D0*2d0*y4**3*SQRT3 -2d0*y4**2*y5 )*FEA33445 +( &
               +2d0*y4*y5**2 -2d0*y4**3/3.D0 )*FEA33455 
            dds1(2,2)=dds2(2,2) &
               +( -3d0*2d0*y2*y4*y5 +3d0*2d0*y2*y5**2*SQRT3/3.D0 )*FEA33345 +( &
               +3d0*2d0*y2*y4**2 +3d0*2d0*y2*y5**2 )*FEA33344 +( &
               +4d0*3d0*y2**2*y4 -SQRT3*4d0*3d0*y2**2*y5 )*FEA33334 &
               +( 5d0*4d0*y2**3 )*FEA33333 +0d0 +0d0 
            dds3(2,2)=dds1(2,2) +( &
               +2d0*y3*y4**2 +2d0*y3*y5**2 +y1*2d0*y5**2 +y1*2d0*y4**2 )*FEA11255 +( &
               +y1*2d0*y5**2*SQRT3/2.D0 &
               +2d0*y3*y5**2*SQRT3/2.D0 -y1*2d0*y4*y5 -2d0*y3*y4*y5 &
               +y1*2d0*y4**2*SQRT3/6.D0 &
               +2d0*y3*y4**2*SQRT3/6.D0 )*FEA13345 
            dds4(2,2)=dds3(2,2) &
               +( -y1*2d0*y4**2*SQRT3/6.D0 -2d0*y3*y4*y5 +y1*2d0*y5**2*SQRT3/2.D0 &
               +2d0*y3*y4**2*SQRT3/3.D0 )*FEA11245 
            dds2(2,2)=dds4(2,2) +( &
               +3d0*2d0*y2*y3*y5/2.D0 -y1*3d0*2d0*y2*y4*SQRT3/2.D0 &
               -y1*3d0*2d0*y2*y5/2.D0 &
               +3d0*2d0*y2*y3*y4*SQRT3/2.D0 )*FEA11135 +( -3d0*2d0*y2*y3*y4/2.D0 &
               +y1*3d0*2d0*y2*y5*SQRT3/2.D0 &
               +3d0*2d0*y2*y3*y5*SQRT3/2.D0 &
               -y1*3d0*2d0*y2*y4/2.D0 )*FEA11134 
            ddv5(2,2)=dds2(2,2) &
               +( y1*4d0*3d0*y2**2 +4d0*3d0*y2**2*y3 )*FEA23333 +( -2.D0*2d0*y3**2*y4 &
               +y1**2*2d0*y4 +SQRT3*y1**2*2d0*y5 )*FEA11334 +( +2d0*y3**3 &
               +y1**2*3d0*2d0*y2 +y1**3*2d0 +3d0*2d0*y2*y3**2 )*FEA11333 +0d0 &
               +( -y1*2d0*y3*y4/2.D0 +SQRT3*y1*2d0*y3*y5/2.D0 )*FEA11234 &
               +( y1*3d0*2d0*y2*y3 )*FEA11123 +( y1**2*2d0*y3 &
               +y1*2d0*y3**2 )*FEA11233 
            dds3(2,2)=( 3d0*2d0*y2*y4**3*SQRT3 &
               -3d0*2d0*y2*y4**2*y5 -5.D0/3.D0*3d0*2d0*y2*y4*y5**2*SQRT3 &
               -3d0*2d0*y2*y5**3 )*FEA333555 &
               +( +4d0*3d0*y2**2*y4*y5 +4d0*3d0*y2**2*y4**2*SQRT3/3.D0 )*FEA222245 +( &
               +y1*5d0*4d0*y2**3 +5d0*4d0*y2**3*y3 )*FEA133333 &
               +( -2.D0*4d0*3d0*y2**2*y3*y4 +y1*4d0*3d0*y2**2*y5*SQRT3 &
               +y1*4d0*3d0*y2**2*y4 )*FEA133334 +0d0 
            dds4(2,2)=dds3(2,2) &
               +( 2.D0/3.D0*SQRT3*y1*2d0*y3**2*y4 -y1**2*2d0*y3*y5 &
               -SQRT3*y1**2*2d0*y3*y4/3.D0 )*FEA112335 &
               +( y1*2d0*y3*y5**2 +y1*2d0*y3*y4**2 )*FEA112355 
            dds2(2,2)=dds4(2,2) &
               +( 3d0*2d0*y2*y3**2*y5 -y1**3*2d0*y5/2.D0 -2d0*y3**3*y5 &
               +y1**3*2d0*y4*SQRT3/2.D0 -y1**2*3d0*2d0*y2*y4*SQRT3/2.D0 &
               +y1**2*3d0*2d0*y2*y5/2.D0 )*FEA222335 &
               +( -y1**2*2d0*y5**2*SQRT3/2.D0 -y1**2*2d0*y4**2*SQRT3/6.D0 &
               -y1**2*2d0*y4*y5 -2.D0/3.D0*2d0*y3**2*y4**2*SQRT3 )*FEA113345 &
               +( 2d0*y3**2*y5**2 +2d0*y3**2*y4**2 +y1**2*2d0*y5**2 &
               +y1**2*2d0*y4**2 )*FEA223355 
            dds3(2,2)=dds2(2,2) +( -y1*2d0*y3*y4*y5 &
               +y1*2d0*y3*y4**2*SQRT3/6.D0 +y1*2d0*y3*y5**2*SQRT3/2.D0 )*FEA123345 &
               +( -y1**3*2d0*y5*SQRT3/2.D0 -y1**3*2d0*y4/2.D0 &
               -y1**2*3d0*2d0*y2*y4/2.D0 &
               +3d0*2d0*y2*y3**2*y4 -y1**2*3d0*2d0*y2*y5*SQRT3/2.D0 &
               +2d0*y3**3*y4 )*FEA222334 +( +2d0*y5**4 +3.D0*2d0*y4**4 &
               +4.D0*2d0*y4*y5**3*SQRT3 )*FEA335555 +( y1**3*3d0*2d0*y2 &
               +3d0*2d0*y2*y3**3 )*FEA222333 
            dds4(2,2)=dds3(2,2) +0d0 &
               +( 5d0*4d0*y2**3*y4 &
               -SQRT3*5d0*4d0*y2**3*y5 )*FEA222224 
            dds5(2,2)=dds4(2,2) &
               +0d0 
            dds1(2,2)=dds5(2,2) +( &
               +SQRT3*y1*3d0*2d0*y2*y3*y5/2.D0 -y1*3d0*2d0*y2*y3*y4/2.D0 )*FEA111234 &
               +0d0 +0d0 
            dds3(2,2)=dds1(2,2) +( 6d0*5d0*y2**4 )*FEA333333 &
               +( y1*4d0*3d0*y2**2*y3 )*FEA111123 +y1**2*2d0*y3**2*FEA112233 +( &
               +4d0*3d0*y2**2*y4**2 &
               +4d0*3d0*y2**2*y5**2 )*FEA222255 
            dds4(2,2)=dds3(2,2) +0d0 +( &
               +5.D0/3.D0*y1*2d0*y4**2*y5*SQRT3 -13.D0/3.D0*2d0*y3*y4*y5**2 &
               +4.D0/3.D0*2d0*y3*y5**3*SQRT3 +y1*2d0*y5**3*SQRT3/3.D0 &
               +2d0*y3*y4**3 &
               -16.D0/3.D0*y1*2d0*y4*y5**2 )*FEA233444 
            dds5(2,2)=dds4(2,2) &
               +( &
               +4.D0*2d0*y3*y4*y5**2*SQRT3 -2.D0*y1*2d0*y5**3 -6.D0*y1*2d0*y4**2*y5 &
               -4.D0*2d0*y3*y5**3 &
               +4.D0*y1*2d0*y4*y5**2*SQRT3 )*FEA113555 
            dds2(2,2)=dds5(2,2) &
               +( -2d0*y4**3*y5 -2.D0/3.D0*2d0*y4**4*SQRT3 &
               -3.D0*2d0*y4*y5**3 )*FEA334445 &
               +0d0 +( &
               +2d0*y4**2*y5**2 -2d0*y4**4 &
               -2.D0*2d0*y4*y5**3*SQRT3 )*FEA334455 
            dds3(2,2)=dds2(2,2) &
               +0d0 +( +y1*3d0*2d0*y2*y4**2 +y1*3d0*2d0*y2*y5**2 +3d0*2d0*y2*y3*y4**2 &
               +3d0*2d0*y2*y3*y5**2 )*FEA233344 &
               +( y1*3d0*2d0*y2*y5**2*SQRT3/6.D0 -3d0*2d0*y2*y3*y5**2*SQRT3/3.D0 &
               -3d0*2d0*y2*y3*y4*y5 -y1*3d0*2d0*y2*y4**2*SQRT3/2.D0 )*FEA233345 &
               +( -3.D0*3d0*2d0*y2*y4*y5**2 +3d0*2d0*y2*y4**3 )*FEA111444 &
               +( y1*3d0*2d0*y2*y3**2 +y1**3*2d0*y3 +y1**2*3d0*2d0*y2*y3 &
               +y1*2d0*y3**3 )*FEA111233 
            dds4(2,2)=dds3(2,2) +0d0 &
               +( -5.D0/3.D0*y1*2d0*y4**2*y5*SQRT3 &
               +y1*2d0*y4**3 -y1*2d0*y5**3*SQRT3/3.D0 &
               +4.D0/3.D0*2d0*y3*y4*y5**2 -4.D0/3.D0*2d0*y3*y5**3*SQRT3 &
               +7.D0/3.D0*y1*2d0*y4*y5**2 )*FEA133444 
            dds5(2,2)=dds4(2,2) +( &
               +2.D0/3.D0*3d0*2d0*y2*y3*y5**2*SQRT3 &
               +y1*3d0*2d0*y2*y5**2*SQRT3/6.D0 -y1*3d0*2d0*y2*y4*y5 &
               +y1*3d0*2d0*y2*y4**2*SQRT3/2.D0 )*FEA133345 
            ddv6(2,2)=dds5(2,2) &
               +( -2d0*y3*y4**2*y5 -y1*2d0*y5**3 &
               +4.D0/3.D0*2d0*y3*y4*y5**2*SQRT3 -y1*2d0*y4**2*y5 -2d0*y3*y5**3 &
               +4.D0/3.D0*y1*2d0*y4*y5**2*SQRT3 )*FEA233445 +( &
               +4d0*3d0*y2**2*y3*y4*SQRT3 -2.D0*y1*4d0*3d0*y2**2*y5 &
               -4d0*3d0*y2**2*y3*y5 )*FEA233335 &
               +( 2d0*y3**4 +y1**2*4d0*3d0*y2**2 +4d0*3d0*y2**2*y3**2 &
               +y1**4*2d0 )*FEA222233 
 
            ddv(2,2)= +ddv1(2,2) +ddv2(2,2) &
               +ddv3(2,2) +ddv4(2,2) +ddv5(2,2) +ddv6(2,2)
 
            ddv1(2,3)=0d0 
            ddv2(2,3)=( 1d0 )*FEA12 +0d0 +0d0 &
               +0d0 
            ddv3(2,3)=( -2.D0*y4 )*FEA124 +0d0 +( 2d0*y3 +2d0*y2 )*FEA112 &
               +0d0 +y1*FEA123 +0d0 +0d0 +0d0 
            dds2(2,3)=0d0 +0d0 &
               +0d0 
            dds1(2,3)=dds2(2,3) +0d0 +0d0 +0d0 +( &
               +3.D0/2.D0*y5**2 -y4**2/2.D0 )*FEA1244 
            dds2(2,3)=dds1(2,3) &
               +( -y5**2/2.D0 +3.D0/2.D0*y4**2 )*FEA1255 +( -SQRT3*2d0*y3*y5/2.D0 &
               +SQRT3*2d0*y2*y5/2.D0 -2d0*y2*y4/2.D0 -2d0*y3*y4/2.D0 )*FEA1124 &
               +( -SQRT3*2d0*y3*y4/2.D0 -SQRT3*2d0*y2*y4/2.D0 -2d0*y2*y5/2.D0 &
               +2d0*y3*y5/2.D0 )*FEA1125 
            ddv4(2,3)=dds2(2,3) +( 3d0*y3**2 &
               +3d0*y2**2 )*FEA1112 +( 2d0*y2*2d0*y3 )*FEA1122 +( y1*2d0*y2 +y1**2 &
               +y1*2d0*y3 )*FEA1123 +0d0 
            dds3(2,3)=0d0 +0d0 +0d0 
            dds2(2,3)=dds3(2,3) &
               +0d0 +0d0 +0d0 
            dds1(2,3)=dds2(2,3) +0d0 +0d0 +0d0 +0d0 +( &
               +y4*y5**2*SQRT3 +5.D0/9.D0*y4**3*SQRT3 )*FEA13445 &
               +( y4*y5**2 -y4**3/3.D0 )*FEA13455 
            dds3(2,3)=dds1(2,3) +( &
               +2d0*y2*y4**2 +2d0*y2*y5**2 +2d0*y3*y4**2 +2d0*y3*y5**2 )*FEA11255 +( &
               +2d0*y2*y5**2*SQRT3/2.D0 +2d0*y3*y4*y5 -2d0*y2*y4*y5 &
               +2d0*y3*y4**2*SQRT3/6.D0 +2d0*y3*y5**2*SQRT3/2.D0 &
               +2d0*y2*y4**2*SQRT3/6.D0 )*FEA13345 
            dds4(2,3)=dds3(2,3) +( &
               +2d0*y3*y4*y5 -2d0*y2*y4*y5 +2d0*y3*y4**2*SQRT3/3.D0 &
               +2d0*y2*y4**2*SQRT3/3.D0 )*FEA11245 
            dds2(2,3)=dds4(2,3) +( &
               +3d0*y2**2*y5/2.D0 -3d0*y3**2*y5/2.D0 +3d0*y2**2*y4*SQRT3/2.D0 &
               +3d0*y3**2*y4*SQRT3/2.D0 )*FEA11135 &
               +( -3d0*y2**2*y4/2.D0 -3d0*y3**2*y4/2.D0 &
               +3d0*y2**2*y5*SQRT3/2.D0 &
               -3d0*y3**2*y5*SQRT3/2.D0 )*FEA11134 
            ddv5(2,3)=dds2(2,3) &
               +( +4d0*y2**3 +4d0*y3**3 )*FEA23333 &
               +( -2.D0*2d0*y2*2d0*y3*y4 )*FEA11334 +( +2d0*y2*3d0*y3**2 &
               +3d0*y2**2*2d0*y3 )*FEA11333 +( y1*y4**2 +y1*y5**2 )*FEA12355 &
               +( -y1*2d0*y3*y4/2.D0 -y1*2d0*y2*y4/2.D0 -SQRT3*y1*2d0*y3*y5/2.D0 &
               +y1**2*y4 +SQRT3*y1*2d0*y2*y5/2.D0 )*FEA11234 +( y1*3d0*y2**2 &
               +y1*3d0*y3**2 +y1**3 )*FEA11123 +( y1**2*2d0*y2 +y1*2d0*y2*2d0*y3 &
               +y1**2*2d0*y3 )*FEA11233 
            dds3(2,3)=0d0 +0d0 +( +5d0*y2**4 &
               +5d0*y3**4 )*FEA133333 &
               +( -2.D0*4d0*y2**3*y4 -2.D0*4d0*y3**3*y4 )*FEA133334 +( -y1*y4**3/3.D0 &
               +y1*y4*y5**2 )*FEA123455 
            dds4(2,3)=dds3(2,3) &
               +( 2.D0/3.D0*SQRT3*y1*2d0*y2*2d0*y3*y4 -y1**2*2d0*y2*y5 &
               -SQRT3*y1**2*2d0*y2*y4/3.D0 &
               +y1**2*2d0*y3*y5 -SQRT3*y1**2*2d0*y3*y4/3.D0 )*FEA112335 &
               +( y1*2d0*y2*y5**2 +y1*2d0*y3*y5**2 +y1*2d0*y3*y4**2 +y1*2d0*y2*y4**2 &
               +y1**2*y4**2 +y1**2*y5**2 )*FEA112355 
            dds2(2,3)=dds4(2,3) &
               +( 3d0*y2**2*2d0*y3*y5 -2d0*y2*3d0*y3**2*y5 )*FEA222335 &
               +( -2.D0/3.D0*2d0*y2*2d0*y3*y4**2*SQRT3 )*FEA113345 &
               +( 2d0*y2*2d0*y3*y5**2 &
               +2d0*y2*2d0*y3*y4**2 )*FEA223355 
            dds3(2,3)=dds2(2,3) &
               +( y1*2d0*y3*y4**2*SQRT3/6.D0 +y1*2d0*y3*y4*y5 &
               +y1*2d0*y3*y5**2*SQRT3/2.D0 &
               +2.D0/3.D0*y1**2*y4**2*SQRT3 -y1*2d0*y2*y4*y5 &
               +y1*2d0*y2*y4**2*SQRT3/6.D0 +y1*2d0*y2*y5**2*SQRT3/2.D0 )*FEA123345 +( &
               +3d0*y2**2*2d0*y3*y4 +2d0*y2*3d0*y3**2*y4 )*FEA222334 +0d0 +( &
               +3d0*y2**2*3d0*y3**2 )*FEA222333 
            dds4(2,3)=dds3(2,3) +0d0 &
               +0d0 
            dds5(2,3)=dds4(2,3) +0d0 
            dds1(2,3)=dds5(2,3) &
               +( -SQRT3*y1*3d0*y3**2*y5/2.D0 +y1**3*y4 &
               +SQRT3*y1*3d0*y2**2*y5/2.D0 -y1*3d0*y2**2*y4/2.D0 &
               -y1*3d0*y3**2*y4/2.D0 )*FEA111234 &
               +0d0 +( -3.D0/4.D0*y4**4 &
               +5.D0/4.D0*y5**4 -7.D0/2.D0*y4**2*y5**2 )*FEA124455 
            dds3(2,3)=dds1(2,3) &
               +0d0 +( y1*4d0*y2**3 +y1**4 +y1*4d0*y3**3 )*FEA111123 &
               +y1**2*2d0*y2*2d0*y3*FEA112233 +0d0 
            dds4(2,3)=dds3(2,3) +( &
               +9.D0*y4**2*y5**2 -3.D0/2.D0*y5**4 +5.D0/2.D0*y4**4 )*FEA134444 &
               +( -13.D0/3.D0*2d0*y2*y4*y5**2 -4.D0/3.D0*2d0*y3*y5**3*SQRT3 &
               +4.D0/3.D0*2d0*y2*y5**3*SQRT3 +2d0*y3*y4**3 &
               +2d0*y2*y4**3 &
               -13.D0/3.D0*2d0*y3*y4*y5**2 )*FEA233444 
            dds5(2,3)=dds4(2,3) &
               +( +4.D0*2d0*y3*y5**3 +4.D0*2d0*y2*y4*y5**2*SQRT3 -4.D0*2d0*y2*y5**3 &
               +4.D0*2d0*y3*y4*y5**2*SQRT3 )*FEA113555 
            dds2(2,3)=dds5(2,3) +0d0 &
               +( -7.D0/12.D0*y5**4*SQRT3 +3.D0/2.D0*y4**2*y5**2*SQRT3 &
               +3.D0/4.D0*y4**4*SQRT3 )*FEA124555 +0d0 
            dds3(2,3)=dds2(2,3) +0d0 &
               +( 3d0*y3**2*y4**2 +3d0*y3**2*y5**2 +3d0*y2**2*y4**2 &
               +3d0*y2**2*y5**2 )*FEA233344 &
               +( -3d0*y2**2*y5**2*SQRT3/3.D0 -3d0*y3**2*y5**2*SQRT3/3.D0 &
               -3d0*y2**2*y4*y5 &
               +3d0*y3**2*y4*y5 )*FEA233345 +0d0 +( y1*3d0*y2**2*2d0*y3 +y1**3*2d0*y2 &
               +y1**2*3d0*y2**2 +y1*2d0*y2*3d0*y3**2 +y1**2*3d0*y3**2 &
               +y1**3*2d0*y3 )*FEA111233 
            dds4(2,3)=dds3(2,3) +0d0 +( &
               +4.D0/3.D0*2d0*y2*y4*y5**2 -4.D0/3.D0*2d0*y2*y5**3*SQRT3 &
               +4.D0/3.D0*2d0*y3*y5**3*SQRT3 &
               +4.D0/3.D0*2d0*y3*y4*y5**2 )*FEA133444 
            dds5(2,3)=dds4(2,3) +( &
               +2.D0/3.D0*3d0*y2**2*y5**2*SQRT3 &
               +2.D0/3.D0*3d0*y3**2*y5**2*SQRT3 )*FEA133345 
            ddv6(2,3)=dds5(2,3) &
               +( -2d0*y2*y4**2*y5 +2d0*y3*y4**2*y5 +2d0*y3*y5**3 &
               +4.D0/3.D0*2d0*y2*y4*y5**2*SQRT3 &
               +4.D0/3.D0*2d0*y3*y4*y5**2*SQRT3 -2d0*y2*y5**3 )*FEA233445 &
               +( 4d0*y3**3*y4*SQRT3 +4d0*y2**3*y4*SQRT3 &
               +4d0*y3**3*y5 -4d0*y2**3*y5 )*FEA233335 +( 2d0*y2*4d0*y3**3 &
               +4d0*y2**3*2d0*y3 )*FEA222233 
 
            ddv(2,3)= +ddv1(2,3) +ddv2(2,3) &
               +ddv3(2,3) +ddv4(2,3) +ddv5(2,3) +ddv6(2,3)
 
            ddv1(2,4)=0d0 
            ddv2(2,4)=0d0 +0d0 +( -1d0/2.D0 )*FEA14 +0d0 
            ddv3(2,4)=( &
               +y1 -2.D0*y3 )*FEA124 +( +3.D0/4.D0*2d0*y4 +SQRT3*y5/2.D0 )*FEA155 &
               +0d0 +0d0 +( +2d0*y4/4.D0 -SQRT3*y5/2.D0 )*FEA144 +0d0 &
               +( -2d0*y2*1d0/2.D0 )*FEA114 
            dds2(2,4)=0d0 &
               +( -9.D0/8.D0*y5**2 -3d0*y4**2/8.D0 &
               +3.D0/8.D0*SQRT3*2d0*y4*y5 )*FEA1444 +( 3.D0/4.D0*2d0*y2*2d0*y4 &
               +SQRT3*2d0*y2*y5/2.D0 )*FEA1155 
            dds1(2,4)=dds2(2,4) +( &
               +2d0*y2*2d0*y4/4.D0 -SQRT3*2d0*y2*y5/2.D0 )*FEA1144 &
               +( -3d0*y2**2*1d0/2.D0 )*FEA1114 +0d0 +( -y3*2d0*y4/2.D0 &
               +y1*2d0*y4 -SQRT3*y1*y5 )*FEA1244 
            dds2(2,4)=dds1(2,4) +( &
               +3.D0/2.D0*y3*2d0*y4 +SQRT3*y1*y5 )*FEA1255 +( &
               +y1**2 -2d0*y2*y3*1d0/2.D0 -y3**2*1d0/2.D0 &
               -y1*2d0*y2*1d0/2.D0 )*FEA1124 &
               +( &
               +SQRT3*y1*2d0*y2*1d0/2.D0 -SQRT3*y3**2*1d0/2.D0 &
               -SQRT3*2d0*y2*y3*1d0/2.D0 )*FEA1125 
            ddv4(2,4)=dds2(2,4) &
               +0d0 +0d0 +0d0 +( 5.D0/8.D0*y5**2 &
               +SQRT3*2d0*y4*y5/8.D0 -3.D0/8.D0*3d0*y4**2 )*FEA1455 
            dds3(2,4)=0d0 +( &
               +4.D0*y5**3*SQRT3 +3.D0*4d0*y4**3 )*FEA25555 &
               +( -4d0*y4**3 -2.D0*y5**3*SQRT3 &
               +2d0*y4*y5**2 )*FEA24455 
            dds2(2,4)=dds3(2,4) +( 3d0*y4**2*y5 &
               +3.D0*y5**3 +2.D0/3.D0*4d0*y4**3*SQRT3 )*FEA24445 +( &
               +4.D0/9.D0*2d0*y2*3d0*y4**2*SQRT3 -2d0*y2*2d0*y4*y5 )*FEA33445 +( &
               +2d0*y2*y5**2 -2d0*y2*3d0*y4**2/3.D0 )*FEA33455 
            dds1(2,4)=dds2(2,4) &
               +( -3d0*y2**2*y5 )*FEA33345 +( +3d0*y2**2*2d0*y4 )*FEA33344 +( &
               +4d0*y2**3 )*FEA33334 +0d0 +( -4.D0/9.D0*y1*3d0*y4**2*SQRT3 &
               +y3*y5**2*SQRT3 -y1*2d0*y4*y5 +5.D0/9.D0*y3*3d0*y4**2*SQRT3 )*FEA13445 &
               +( y3*y5**2 &
               +y1*y5**2 -y3*3d0*y4**2/3.D0 &
               -y1*3d0*y4**2/3.D0 )*FEA13455 
            dds3(2,4)=dds1(2,4) &
               +( +2d0*y2*y3*2d0*y4 +y1*2d0*y2*2d0*y4 +y3**2*2d0*y4 &
               +y1**2*2d0*y4 )*FEA11255 +( -y1*2d0*y2*y5 +y3**2*y5 -2d0*y2*y3*y5 &
               +y3**2*2d0*y4*SQRT3/6.D0 +y1*2d0*y2*2d0*y4*SQRT3/6.D0 &
               +2.D0/3.D0*y1**2*2d0*y4*SQRT3 &
               +2d0*y2*y3*2d0*y4*SQRT3/6.D0 )*FEA13345 
            dds4(2,4)=dds3(2,4) &
               +( y1**2*y5 +y1**2*2d0*y4*SQRT3/3.D0 -y1*2d0*y2*2d0*y4*SQRT3/6.D0 &
               +y3**2*y5 -2d0*y2*y3*y5 +y3**2*2d0*y4*SQRT3/3.D0 &
               +2d0*y2*y3*2d0*y4*SQRT3/3.D0 )*FEA11245 
            dds2(2,4)=dds4(2,4) &
               +( -y1*3d0*y2**2*SQRT3/2.D0 +3d0*y2**2*y3*SQRT3/2.D0 &
               +y3**3*SQRT3/2.D0 )*FEA11135 +( -3d0*y2**2*y3*1d0/2.D0 &
               +y1**3 -y3**3*1d0/2.D0 &
               -y1*3d0*y2**2*1d0/2.D0 )*FEA11134 
            ddv5(2,4)=dds2(2,4) &
               +0d0 +( -2.D0*2d0*y2*y3**2 +y1**2*2d0*y2 )*FEA11334 +0d0 &
               +( y1*y3*2d0*y4 )*FEA12355 &
               +( -y1*y3**2*1d0/2.D0 -y1*2d0*y2*y3*1d0/2.D0 +y1**2*y3 )*FEA11234 +0d0 &
               +0d0 
            dds3(2,4)=( 3d0*y2**2*3d0*y4**2*SQRT3 -3d0*y2**2*2d0*y4*y5 &
               -5.D0/3.D0*3d0*y2**2*y5**2*SQRT3 )*FEA333555 &
               +( +4d0*y2**3*y5 +4d0*y2**3*2d0*y4*SQRT3/3.D0 )*FEA222245 +0d0 &
               +( -2.D0*4d0*y2**3*y3 +y1**4 -2.D0*y3**4 +y1*4d0*y2**3 )*FEA133334 &
               +( -y1*y3*3d0*y4**2/3.D0 +y1*y3*y5**2 )*FEA123455 
            dds4(2,4)=dds3(2,4) &
               +( 2.D0/3.D0*SQRT3*y1*2d0*y2*y3**2 -SQRT3*y1**2*2d0*y2*y3*1d0/3.D0 &
               -SQRT3*y1**2*y3**2*1d0/3.D0 )*FEA112335 &
               +( +y1*y3**2*2d0*y4 +y1*2d0*y2*y3*2d0*y4 &
               +y1**2*y3*2d0*y4 )*FEA112355 
            dds2(2,4)=dds4(2,4) +( &
               +y1**3*2d0*y2*SQRT3/2.D0 -y1**2*3d0*y2**2*SQRT3/2.D0 )*FEA222335 &
               +( -y1**2*2d0*y2*2d0*y4*SQRT3/6.D0 -y1**2*2d0*y2*y5 &
               -2.D0/3.D0*2d0*y2*y3**2*2d0*y4*SQRT3 )*FEA113345 &
               +( +2d0*y2*y3**2*2d0*y4 &
               +y1**2*2d0*y2*2d0*y4 )*FEA223355 
            dds3(2,4)=dds2(2,4) &
               +( y1*y3**2*2d0*y4*SQRT3/6.D0 +y1*y3**2*y5 &
               +2.D0/3.D0*y1**2*y3*2d0*y4*SQRT3 -y1*2d0*y2*y3*y5 &
               +y1*2d0*y2*y3*2d0*y4*SQRT3/6.D0 )*FEA123345 &
               +( -y1**3*2d0*y2*1d0/2.D0 -y1**2*3d0*y2**2*1d0/2.D0 +3d0*y2**2*y3**2 &
               +2d0*y2*y3**3 )*FEA222334 +( +3.D0*2d0*y2*4d0*y4**3 &
               +4.D0*2d0*y2*y5**3*SQRT3 )*FEA335555 +0d0 
            dds4(2,4)=dds3(2,4) &
               +( -4d0*y4**3*y5*SQRT3/2.D0 +3d0*y4**2*y5**2 &
               +5d0*y4**4/5.D0 )*FEA244455 &
               +( 5d0*y2**4 )*FEA222224 
            dds5(2,4)=dds4(2,4) +( -7.D0/15.D0*5d0*y4**4 &
               +4d0*y4**3*y5*SQRT3/3.D0 +y5**4 )*FEA145555 
            dds1(2,4)=dds5(2,4) +( &
               +y1**3*y3 -y1*3d0*y2**2*y3*1d0/2.D0 -y1*y3**3*1d0/2.D0 )*FEA111234 &
               +( -4d0*y4**3*y5/3.D0 -y5**4*SQRT3/2.D0 +5d0*y4**4*SQRT3/18.D0 &
               +2d0*y4*y5**3 )*FEA244555 &
               +( y1*2d0*y4*y5**2 -3.D0/4.D0*y3*4d0*y4**3 -7.D0/2.D0*y3*2d0*y4*y5**2 &
               -2.D0*y1*3d0*y4**2*y5*SQRT3 )*FEA124455 
            dds3(2,4)=dds1(2,4) &
               +0d0 +0d0 +( +4d0*y2**3*2d0*y4 )*FEA222255 
            dds4(2,4)=dds3(2,4) +( &
               +9.D0*y3*2d0*y4*y5**2 +y1*4d0*y4**3 +4.D0*y1*3d0*y4**2*y5*SQRT3 &
               +5.D0/2.D0*y3*4d0*y4**3 )*FEA134444 +( &
               +5.D0/3.D0*y1*2d0*y2*2d0*y4*y5*SQRT3 -13.D0/3.D0*2d0*y2*y3*y5**2 &
               -7.D0/3.D0*y1**2*y5**2 &
               +3.D0*y1**2*3d0*y4**2 +y3**2*3d0*y4**2 &
               +2d0*y2*y3*3d0*y4**2 -13.D0/3.D0*y3**2*y5**2 &
               -4.D0/3.D0*y1**2*2d0*y4*y5*SQRT3 &
               -16.D0/3.D0*y1*2d0*y2*y5**2 )*FEA233444 
            dds5(2,4)=dds4(2,4) &
               +( +4.D0*2d0*y2*y3*y5**2*SQRT3 -6.D0*y1*2d0*y2*2d0*y4*y5 &
               +y1**2*y5**2*SQRT3 -3.D0*y1**2*3d0*y4**2*SQRT3 +3.D0*y1**2*2d0*y4*y5 &
               +4.D0*y3**2*y5**2*SQRT3 &
               +4.D0*y1*2d0*y2*y5**2*SQRT3 )*FEA113555 
            dds2(2,4)=dds5(2,4) &
               +( -2d0*y2*3d0*y4**2*y5 -2.D0/3.D0*2d0*y2*4d0*y4**3*SQRT3 &
               -3.D0*2d0*y2*y5**3 )*FEA334445 &
               +( +3.D0*y1*3d0*y4**2*y5 +3.D0/2.D0*y3*2d0*y4*y5**2*SQRT3 +y1*y5**3 &
               +3.D0/4.D0*y3*4d0*y4**3*SQRT3 )*FEA124555 +( &
               +2d0*y2*2d0*y4*y5**2 -2d0*y2*4d0*y4**3 &
               -2.D0*2d0*y2*y5**3*SQRT3 )*FEA334455 
            dds3(2,4)=dds2(2,4) &
               +0d0 +( y3**3*2d0*y4 +y1*3d0*y2**2*2d0*y4 +y1**3*2d0*y4 &
               +3d0*y2**2*y3*2d0*y4 )*FEA233344 +( +y1**3*y5 -3d0*y2**2*y3*y5 &
               +y3**3*y5 -y1*3d0*y2**2*2d0*y4*SQRT3/2.D0 )*FEA233345 &
               +( -3.D0*3d0*y2**2*y5**2 +3d0*y2**2*3d0*y4**2 )*FEA111444 &
               +0d0 
            dds4(2,4)=dds3(2,4) +0d0 +( -5.D0/3.D0*y1*2d0*y2*2d0*y4*y5*SQRT3 &
               +y1*2d0*y2*3d0*y4**2 -2.D0*y1**2*3d0*y4**2 +4.D0/3.D0*2d0*y2*y3*y5**2 &
               +7.D0/3.D0*y1*2d0*y2*y5**2 +4.D0/3.D0*y1**2*2d0*y4*y5*SQRT3 &
               +4.D0/3.D0*y3**2*y5**2 &
               -2.D0/3.D0*y1**2*y5**2 )*FEA133444 
            dds5(2,4)=dds4(2,4) &
               +( -y1**3*y5 +y1**3*2d0*y4*SQRT3/2.D0 -y1*3d0*y2**2*y5 &
               +y1*3d0*y2**2*2d0*y4*SQRT3/2.D0 )*FEA133345 
            ddv6(2,4)=dds5(2,4) &
               +( -2d0*y2*y3*2d0*y4*y5 +y3**2*2d0*y4*y5 &
               +4.D0/3.D0*2d0*y2*y3*y5**2*SQRT3 &
               +4.D0/3.D0*y3**2*y5**2*SQRT3 -y1*2d0*y2*2d0*y4*y5 &
               +y1**2*y5**2*SQRT3/3.D0 -y1**2*3d0*y4**2*SQRT3 &
               +4.D0/3.D0*y1*2d0*y2*y5**2*SQRT3 )*FEA233445 +( y3**4*SQRT3 &
               +4d0*y2**3*y3*SQRT3 -y1**4*SQRT3 )*FEA233335 +0d0 
 
            ddv(2,4)= &
               +ddv1(2,4) +ddv2(2,4) +ddv3(2,4) +ddv4(2,4) +ddv5(2,4) +ddv6(2,4)
 
            ddv1(2,5)=0d0 
            ddv2(2,5)=0d0 +0d0 +( +SQRT3*1d0/2.D0 )*FEA14 &
               +0d0 
            ddv3(2,5)=( +SQRT3*y1 )*FEA124 +( +2d0*y5/4.D0 &
               +SQRT3*y4*1d0/2.D0 )*FEA155 +0d0 +0d0 +( &
               +3.D0/4.D0*2d0*y5 -SQRT3*y4*1d0/2.D0 )*FEA144 +0d0 +( &
               +SQRT3*2d0*y2*1d0/2.D0 )*FEA114 
            dds2(2,5)=0d0 &
               +( 3.D0/8.D0*SQRT3*3d0*y5**2 -9.D0/8.D0*y4*2d0*y5 &
               +3.D0/8.D0*SQRT3*y4**2 )*FEA1444 +( +SQRT3*2d0*y2*y4*1d0/2.D0 &
               +2d0*y2*2d0*y5/4.D0 )*FEA1155 
            dds1(2,5)=dds2(2,5) &
               +( -SQRT3*2d0*y2*y4*1d0/2.D0 +3.D0/4.D0*2d0*y2*2d0*y5 )*FEA1144 +( &
               +SQRT3*3d0*y2**2*1d0/2.D0 )*FEA1114 +0d0 +( &
               +3.D0/2.D0*y3*2d0*y5 -SQRT3*y1*y4 )*FEA1244 
            dds2(2,5)=dds1(2,5) +( &
               +y1*2d0*y5 -y3*2d0*y5/2.D0 +SQRT3*y1*y4 )*FEA1255 &
               +( -SQRT3*y3**2*1d0/2.D0 +SQRT3*2d0*y2*y3*1d0/2.D0 &
               +SQRT3*y1*2d0*y2*1d0/2.D0 )*FEA1124 +( y1**2 -2d0*y2*y3*1d0/2.D0 &
               +y3**2*1d0/2.D0 +y1*2d0*y2*1d0/2.D0 )*FEA1125 
            ddv4(2,5)=dds2(2,5) +0d0 &
               +0d0 +0d0 +( 5.D0/8.D0*y4*2d0*y5 +SQRT3*3d0*y5**2/8.D0 &
               +SQRT3*y4**2*1d0/8.D0 )*FEA1455 
            dds3(2,5)=0d0 +( &
               +4.D0*y4*3d0*y5**2*SQRT3 +4d0*y5**3 )*FEA25555 &
               +( -2.D0*y4*3d0*y5**2*SQRT3 &
               +y4**2*2d0*y5 )*FEA24455 
            dds2(2,5)=dds3(2,5) +( y4**3 &
               +3.D0*y4*3d0*y5**2 )*FEA24445 &
               +( -2d0*y2*3d0*y5**2 -2d0*y2*y4**2 )*FEA33445 +( &
               +2d0*y2*y4*2d0*y5 )*FEA33455 
            dds1(2,5)=dds2(2,5) +( -3d0*y2**2*y4 &
               +3d0*y2**2*2d0*y5*SQRT3/3.D0 )*FEA33345 +( &
               +3d0*y2**2*2d0*y5 )*FEA33344 +( -SQRT3*4d0*y2**3 )*FEA33334 +0d0 &
               +( -y1*3d0*y5**2 +y3*y4*2d0*y5*SQRT3 -y1*y4**2 )*FEA13445 &
               +( y3*y4*2d0*y5 +y1*y4*2d0*y5 )*FEA13455 
            dds3(2,5)=dds1(2,5) +( &
               +2d0*y2*y3*2d0*y5 +y1*2d0*y2*2d0*y5 +y1**2*2d0*y5 &
               +y3**2*2d0*y5 )*FEA11255 +( +y1*2d0*y2*2d0*y5*SQRT3/2.D0 &
               +2d0*y2*y3*2d0*y5*SQRT3/2.D0 -y1*2d0*y2*y4 +y3**2*y4 -2d0*y2*y3*y4 &
               +y3**2*2d0*y5*SQRT3/2.D0 )*FEA13345 
            dds4(2,5)=dds3(2,5) +( y1**2*y4 &
               +y3**2*y4 -2d0*y2*y3*y4 &
               +y1*2d0*y2*2d0*y5*SQRT3/2.D0 )*FEA11245 
            dds2(2,5)=dds4(2,5) +( -y1**3 &
               +3d0*y2**2*y3*1d0/2.D0 -y1*3d0*y2**2*1d0/2.D0 &
               -y3**3*1d0/2.D0 )*FEA11135 &
               +( +y1*3d0*y2**2*SQRT3/2.D0 &
               +3d0*y2**2*y3*SQRT3/2.D0 &
               -y3**3*SQRT3/2.D0 )*FEA11134 
            ddv5(2,5)=dds2(2,5) &
               +0d0 +( +SQRT3*y1**2*2d0*y2 )*FEA11334 +0d0 +( &
               +y1*y3*2d0*y5 )*FEA12355 +( -SQRT3*y1*y3**2*1d0/2.D0 &
               +SQRT3*y1*2d0*y2*y3*1d0/2.D0 )*FEA11234 +0d0 &
               +0d0 
            dds3(2,5)=( -3d0*y2**2*y4**2 -5.D0/3.D0*3d0*y2**2*y4*2d0*y5*SQRT3 &
               -3d0*y2**2*3d0*y5**2 )*FEA333555 &
               +( +4d0*y2**3*y4 )*FEA222245 +0d0 +( +y1*4d0*y2**3*SQRT3 &
               +y1**4*SQRT3 )*FEA133334 +( &
               +y1*y3*y4*2d0*y5 )*FEA123455 
            dds4(2,5)=dds3(2,5) +( -y1**2*2d0*y2*y3 &
               +y1**2*y3**2 )*FEA112335 +( y1*2d0*y2*y3*2d0*y5 +y1*y3**2*2d0*y5 &
               +y1**2*y3*2d0*y5 )*FEA112355 
            dds2(2,5)=dds4(2,5) &
               +( 3d0*y2**2*y3**2 -y1**3*2d0*y2*1d0/2.D0 -2d0*y2*y3**3 &
               +y1**2*3d0*y2**2*1d0/2.D0 )*FEA222335 &
               +( -y1**2*2d0*y2*2d0*y5*SQRT3/2.D0 -y1**2*2d0*y2*y4 )*FEA113345 &
               +( 2d0*y2*y3**2*2d0*y5 &
               +y1**2*2d0*y2*2d0*y5 )*FEA223355 
            dds3(2,5)=dds2(2,5) +( +y1*y3**2*y4 &
               +y1*y3**2*2d0*y5*SQRT3/2.D0 -y1*2d0*y2*y3*y4 &
               +y1*2d0*y2*y3*2d0*y5*SQRT3/2.D0 )*FEA123345 &
               +( -y1**3*2d0*y2*SQRT3/2.D0 -y1**2*3d0*y2**2*SQRT3/2.D0 )*FEA222334 +( &
               +2d0*y2*4d0*y5**3 +4.D0*2d0*y2*y4*3d0*y5**2*SQRT3 )*FEA335555 &
               +0d0 
            dds4(2,5)=dds3(2,5) &
               +( -y4**4*SQRT3/2.D0 -3.D0/10.D0*5d0*y5**4*SQRT3 &
               +y4**3*2d0*y5 )*FEA244455 &
               +( -SQRT3*5d0*y2**4 )*FEA222224 
            dds5(2,5)=dds4(2,5) +( &
               +5d0*y5**4*SQRT3/5.D0 +y4**4*SQRT3/3.D0 &
               +y4*4d0*y5**3 )*FEA145555 
            dds1(2,5)=dds5(2,5) &
               +( -SQRT3*y1*y3**3*1d0/2.D0 &
               +SQRT3*y1*3d0*y2**2*y3*1d0/2.D0 )*FEA111234 &
               +( -y4**4*1d0/3.D0 -y4*4d0*y5**3*SQRT3/2.D0 &
               +y4**2*3d0*y5**2 )*FEA244555 +( y1*y4**2*2d0*y5 -y1*4d0*y5**3 &
               +5.D0/4.D0*y3*4d0*y5**3 -7.D0/2.D0*y3*y4**2*2d0*y5 &
               -2.D0*y1*y4**3*SQRT3 )*FEA124455 
            dds3(2,5)=dds1(2,5) &
               +0d0 +0d0 +( +4d0*y2**3*2d0*y5 )*FEA222255 
            dds4(2,5)=dds3(2,5) +( &
               +9.D0*y3*y4**2*2d0*y5 -3.D0/2.D0*y3*4d0*y5**3 +4.D0*y1*y4**3*SQRT3 &
               +3.D0*y1*4d0*y5**3 )*FEA134444 +( &
               +5.D0/3.D0*y1*2d0*y2*y4**2*SQRT3 -13.D0/3.D0*2d0*y2*y3*y4*2d0*y5 &
               -4.D0/3.D0*y3**2*3d0*y5**2*SQRT3 -7.D0/3.D0*y1**2*y4*2d0*y5 &
               +4.D0/3.D0*2d0*y2*y3*3d0*y5**2*SQRT3 &
               +y1*2d0*y2*3d0*y5**2*SQRT3/3.D0 -13.D0/3.D0*y3**2*y4*2d0*y5 &
               -4.D0/3.D0*y1**2*y4**2*SQRT3 &
               -16.D0/3.D0*y1*2d0*y2*y4*2d0*y5 )*FEA233444 
            dds5(2,5)=dds4(2,5) &
               +( +4.D0*y3**2*3d0*y5**2 &
               +4.D0*2d0*y2*y3*y4*2d0*y5*SQRT3 -2.D0*y1*2d0*y2*3d0*y5**2 &
               -6.D0*y1*2d0*y2*y4**2 &
               +y1**2*y4*2d0*y5*SQRT3 -4.D0*2d0*y2*y3*3d0*y5**2 &
               +3.D0*y1**2*y4**2 -y1**2*3d0*y5**2 +4.D0*y3**2*y4*2d0*y5*SQRT3 &
               +4.D0*y1*2d0*y2*y4*2d0*y5*SQRT3 )*FEA113555 
            dds2(2,5)=dds5(2,5) &
               +( -2d0*y2*y4**3 -3.D0*2d0*y2*y4*3d0*y5**2 )*FEA334445 +( &
               +2.D0/3.D0*y1*4d0*y5**3*SQRT3 &
               +3.D0*y1*y4**3 -7.D0/12.D0*y3*4d0*y5**3*SQRT3 &
               +3.D0/2.D0*y3*y4**2*2d0*y5*SQRT3 +y1*y4*3d0*y5**2 )*FEA124555 +( &
               +2d0*y2*y4**2*2d0*y5 &
               -2.D0*2d0*y2*y4*3d0*y5**2*SQRT3 )*FEA334455 
            dds3(2,5)=dds2(2,5) &
               +0d0 +( +y3**3*2d0*y5 +y1*3d0*y2**2*2d0*y5 +y1**3*2d0*y5 &
               +3d0*y2**2*y3*2d0*y5 )*FEA233344 &
               +( y1*3d0*y2**2*2d0*y5*SQRT3/6.D0 -3d0*y2**2*y3*2d0*y5*SQRT3/3.D0 &
               -y3**3*2d0*y5*SQRT3/3.D0 &
               +y1**3*y4 -y1**3*2d0*y5*SQRT3/3.D0 -3d0*y2**2*y3*y4 &
               +y3**3*y4 )*FEA233345 +( -3.D0*3d0*y2**2*y4*2d0*y5 )*FEA111444 &
               +0d0 
            dds4(2,5)=dds3(2,5) +0d0 &
               +( -5.D0/3.D0*y1*2d0*y2*y4**2*SQRT3 -y1*2d0*y2*3d0*y5**2*SQRT3/3.D0 &
               +4.D0/3.D0*2d0*y2*y3*y4*2d0*y5 -4.D0/3.D0*2d0*y2*y3*3d0*y5**2*SQRT3 &
               +7.D0/3.D0*y1*2d0*y2*y4*2d0*y5 +4.D0/3.D0*y3**2*3d0*y5**2*SQRT3 &
               +4.D0/3.D0*y1**2*y4**2*SQRT3 &
               +4.D0/3.D0*y3**2*y4*2d0*y5 &
               -2.D0/3.D0*y1**2*y4*2d0*y5 )*FEA133444 
            dds5(2,5)=dds4(2,5) &
               +( -y1**3*y4 +2.D0/3.D0*3d0*y2**2*y3*2d0*y5*SQRT3 &
               +y1**3*2d0*y5*SQRT3/6.D0 +y1*3d0*y2**2*2d0*y5*SQRT3/6.D0 &
               +2.D0/3.D0*y3**3*2d0*y5*SQRT3 &
               -y1*3d0*y2**2*y4 )*FEA133345 
            ddv6(2,5)=dds5(2,5) &
               +( -2d0*y2*y3*y4**2 +y3**2*y4**2 +y3**2*3d0*y5**2 -y1*2d0*y2*3d0*y5**2 &
               +4.D0/3.D0*2d0*y2*y3*y4*2d0*y5*SQRT3 &
               +4.D0/3.D0*y3**2*y4*2d0*y5*SQRT3 -y1*2d0*y2*y4**2 -2d0*y2*y3*3d0*y5**2 &
               +y1**2*y4*2d0*y5*SQRT3/3.D0 &
               +4.D0/3.D0*y1*2d0*y2*y4*2d0*y5*SQRT3 )*FEA233445 +( -y1**4 &
               +y3**4 -2.D0*y1*4d0*y2**3 -4d0*y2**3*y3 )*FEA233335 &
               +0d0 
 
            ddv(2,5)= +ddv1(2,5) +ddv2(2,5) +ddv3(2,5) +ddv4(2,5) &
               +ddv5(2,5) +ddv6(2,5)
 
            ddv1(2,6)=( +1d0 )*DFEA1 
            ddv2(2,6)=( y3 +y1 )*DFEA12 &
               +( 2d0*y2 )*DFEA11 +( +SQRT3*y5/2.D0 -y4/2.D0 )*DFEA14 &
               +0d0 
            ddv3(2,6)=( +y1*y4 -2.D0*y3*y4 +SQRT3*y1*y5 )*DFEA124 +( &
               +y5**2/4.D0 +3.D0/4.D0*y4**2 +SQRT3*y4*y5/2.D0 )*DFEA155 +( y3**2 &
               +y1*2d0*y2 +2d0*y2*y3 +y1**2 )*DFEA112 +0d0 +y1*y3*DFEA123 +( &
               +3.D0/4.D0*y5**2 +y4**2/4.D0 -SQRT3*y4*y5/2.D0 )*DFEA144 +( &
               +3d0*y2**2 )*DFEA111 +( -2d0*y2*y4/2.D0 &
               +SQRT3*2d0*y2*y5/2.D0 )*DFEA114 
            dds2(2,6)=0d0 &
               +( 3.D0/8.D0*SQRT3*y5**3 -9.D0/8.D0*y4*y5**2 -y4**3/8.D0 &
               +3.D0/8.D0*SQRT3*y4**2*y5 )*DFEA1444 +( 3.D0/4.D0*2d0*y2*y4**2 &
               +SQRT3*2d0*y2*y4*y5/2.D0 &
               +2d0*y2*y5**2/4.D0 )*DFEA1155 
            dds1(2,6)=dds2(2,6) +( &
               +2d0*y2*y4**2/4.D0 -SQRT3*2d0*y2*y4*y5/2.D0 &
               +3.D0/4.D0*2d0*y2*y5**2 )*DFEA1144 +( &
               +SQRT3*3d0*y2**2*y5/2.D0 -3d0*y2**2*y4/2.D0 )*DFEA1114 &
               +( 4d0*y2**3 )*DFEA1111 +( +3.D0/2.D0*y3*y5**2 -y3*y4**2/2.D0 &
               +y1*y4**2 -SQRT3*y1*y4*y5 )*DFEA1244 
            dds2(2,6)=dds1(2,6) +( &
               +y1*y5**2 -y3*y5**2/2.D0 +3.D0/2.D0*y3*y4**2 &
               +SQRT3*y1*y4*y5 )*DFEA1255 +( -SQRT3*y3**2*y5/2.D0 +y1**2*y4 &
               +SQRT3*2d0*y2*y3*y5/2.D0 -2d0*y2*y3*y4/2.D0 &
               +SQRT3*y1*2d0*y2*y5/2.D0 -y3**2*y4/2.D0 -y1*2d0*y2*y4/2.D0 )*DFEA1124 &
               +( y1**2*y5 &
               +SQRT3*y1*2d0*y2*y4/2.D0 -SQRT3*y3**2*y4/2.D0 -SQRT3*2d0*y2*y3*y4/2.D0 &
               -2d0*y2*y3*y5/2.D0 &
               +y3**2*y5/2.D0 +y1*2d0*y2*y5/2.D0 )*DFEA1125 
            ddv4(2,6)=dds2(2,6) &
               +( y3**3 +y1**3 +y1*3d0*y2**2 +3d0*y2**2*y3 )*DFEA1112 +( 2d0*y2*y3**2 &
               +y1**2*2d0*y2 )*DFEA1122 +( y1*2d0*y2*y3 +y1**2*y3 &
               +y1*y3**2 )*DFEA1123 +( 5.D0/8.D0*y4*y5**2 +SQRT3*y5**3/8.D0 &
               +SQRT3*y4**2*y5/8.D0 -3.D0/8.D0*y4**3 )*DFEA1455 
            dds3(2,6)=0d0 +( &
               +4.D0*y4*y5**3*SQRT3 +3.D0*y4**4 +y5**4 )*DFEA25555 &
               +( -y4**4 -2.D0*y4*y5**3*SQRT3 &
               +y4**2*y5**2 )*DFEA24455 
            dds2(2,6)=dds3(2,6) +( y4**3*y5 &
               +3.D0*y4*y5**3 +2.D0/3.D0*y4**4*SQRT3 )*DFEA24445 +( -2d0*y2*y5**3 &
               +4.D0/9.D0*2d0*y2*y4**3*SQRT3 -2d0*y2*y4**2*y5 )*DFEA33445 +( &
               +2d0*y2*y4*y5**2 -2d0*y2*y4**3/3.D0 )*DFEA33455 
            dds1(2,6)=dds2(2,6) &
               +( -3d0*y2**2*y4*y5 +3d0*y2**2*y5**2*SQRT3/3.D0 )*DFEA33345 +( &
               +3d0*y2**2*y4**2 +3d0*y2**2*y5**2 )*DFEA33344 +( &
               +4d0*y2**3*y4 -SQRT3*4d0*y2**3*y5 )*DFEA33334 +( 5d0*y2**4 )*DFEA33333 &
               +( -4.D0/9.D0*y1*y4**3*SQRT3 -y1*y5**3 +y3*y4*y5**2*SQRT3 -y1*y4**2*y5 &
               +5.D0/9.D0*y3*y4**3*SQRT3 )*DFEA13445 +( y3*y4*y5**2 &
               +y1*y4*y5**2 -y3*y4**3/3.D0 &
               -y1*y4**3/3.D0 )*DFEA13455 
            dds3(2,6)=dds1(2,6) &
               +( +2d0*y2*y3*y4**2 +2d0*y2*y3*y5**2 +y1*2d0*y2*y5**2 +y1**2*y5**2 &
               +y1*2d0*y2*y4**2 +y3**2*y4**2 +y1**2*y4**2 +y3**2*y5**2 )*DFEA11255 +( &
               +y1*2d0*y2*y5**2*SQRT3/2.D0 &
               +2d0*y2*y3*y5**2*SQRT3/2.D0 -y1*2d0*y2*y4*y5 &
               +y3**2*y4*y5 -2d0*y2*y3*y4*y5 +y3**2*y4**2*SQRT3/6.D0 &
               +y1*2d0*y2*y4**2*SQRT3/6.D0 +2.D0/3.D0*y1**2*y4**2*SQRT3 &
               +y3**2*y5**2*SQRT3/2.D0 &
               +2d0*y2*y3*y4**2*SQRT3/6.D0 )*DFEA13345 
            dds4(2,6)=dds3(2,6) &
               +( y1**2*y4*y5 +y1**2*y4**2*SQRT3/3.D0 -y1*2d0*y2*y4**2*SQRT3/6.D0 &
               +y3**2*y4*y5 -2d0*y2*y3*y4*y5 +y3**2*y4**2*SQRT3/3.D0 &
               +y1*2d0*y2*y5**2*SQRT3/2.D0 &
               +2d0*y2*y3*y4**2*SQRT3/3.D0 )*DFEA11245 
            dds2(2,6)=dds4(2,6) &
               +( -y1**3*y5 &
               +3d0*y2**2*y3*y5/2.D0 -y1*3d0*y2**2*y4*SQRT3/2.D0 &
               -y1*3d0*y2**2*y5/2.D0 -y3**3*y5/2.D0 &
               +3d0*y2**2*y3*y4*SQRT3/2.D0 +y3**3*y4*SQRT3/2.D0 )*DFEA11135 &
               +( -3d0*y2**2*y3*y4/2.D0 +y1**3*y4 -y3**3*y4/2.D0 &
               +y1*3d0*y2**2*y5*SQRT3/2.D0 &
               +3d0*y2**2*y3*y5*SQRT3/2.D0 -y3**3*y5*SQRT3/2.D0 &
               -y1*3d0*y2**2*y4/2.D0 )*DFEA11134 
            ddv5(2,6)=dds2(2,6) &
               +( y1*4d0*y2**3 +y1**4 +4d0*y2**3*y3 +y3**4 )*DFEA23333 &
               +( -2.D0*2d0*y2*y3**2*y4 +y1**2*2d0*y2*y4 &
               +SQRT3*y1**2*2d0*y2*y5 )*DFEA11334 +( +2d0*y2*y3**3 +y1**2*3d0*y2**2 &
               +y1**3*2d0*y2 +3d0*y2**2*y3**2 )*DFEA11333 +( y1*y3*y4**2 &
               +y1*y3*y5**2 )*DFEA12355 &
               +( -y1*y3**2*y4/2.D0 -y1*2d0*y2*y3*y4/2.D0 -SQRT3*y1*y3**2*y5/2.D0 &
               +y1**2*y3*y4 +SQRT3*y1*2d0*y2*y3*y5/2.D0 )*DFEA11234 &
               +( y1*3d0*y2**2*y3 +y1*y3**3 +y1**3*y3 )*DFEA11123 +( y1**2*2d0*y2*y3 &
               +y1*2d0*y2*y3**2 &
               +y1**2*y3**2 )*DFEA11233 
            dds3(2,6)=( 3d0*y2**2*y4**3*SQRT3 &
               -3d0*y2**2*y4**2*y5 -5.D0/3.D0*3d0*y2**2*y4*y5**2*SQRT3 &
               -3d0*y2**2*y5**3 )*DFEA333555 &
               +( +4d0*y2**3*y4*y5 +4d0*y2**3*y4**2*SQRT3/3.D0 )*DFEA222245 +( &
               +y1*5d0*y2**4 +5d0*y2**4*y3 +y1**5 +y3**5 )*DFEA133333 &
               +( -2.D0*4d0*y2**3*y3*y4 +y1**4*y4 &
               +y1*4d0*y2**3*y5*SQRT3 -2.D0*y3**4*y4 +y1**4*y5*SQRT3 &
               +y1*4d0*y2**3*y4 )*DFEA133334 +( -y1*y3*y4**3/3.D0 &
               +y1*y3*y4*y5**2 )*DFEA123455 
            dds4(2,6)=dds3(2,6) &
               +( 2.D0/3.D0*SQRT3*y1*2d0*y2*y3**2*y4 -y1**2*2d0*y2*y3*y5 &
               -SQRT3*y1**2*2d0*y2*y3*y4/3.D0 &
               +y1**2*y3**2*y5 -SQRT3*y1**2*y3**2*y4/3.D0 )*DFEA112335 &
               +( y1*2d0*y2*y3*y5**2 +y1*y3**2*y5**2 +y1*y3**2*y4**2 &
               +y1*2d0*y2*y3*y4**2 +y1**2*y3*y4**2 &
               +y1**2*y3*y5**2 )*DFEA112355 
            dds2(2,6)=dds4(2,6) &
               +( 3d0*y2**2*y3**2*y5 -y1**3*2d0*y2*y5/2.D0 -2d0*y2*y3**3*y5 &
               +y1**3*2d0*y2*y4*SQRT3/2.D0 -y1**2*3d0*y2**2*y4*SQRT3/2.D0 &
               +y1**2*3d0*y2**2*y5/2.D0 )*DFEA222335 &
               +( -y1**2*2d0*y2*y5**2*SQRT3/2.D0 -y1**2*2d0*y2*y4**2*SQRT3/6.D0 &
               -y1**2*2d0*y2*y4*y5 -2.D0/3.D0*2d0*y2*y3**2*y4**2*SQRT3 )*DFEA113345 &
               +( 2d0*y2*y3**2*y5**2 +2d0*y2*y3**2*y4**2 +y1**2*2d0*y2*y5**2 &
               +y1**2*2d0*y2*y4**2 )*DFEA223355 
            dds3(2,6)=dds2(2,6) &
               +( y1*y3**2*y4**2*SQRT3/6.D0 +y1*y3**2*y4*y5 &
               +y1*y3**2*y5**2*SQRT3/2.D0 &
               +2.D0/3.D0*y1**2*y3*y4**2*SQRT3 -y1*2d0*y2*y3*y4*y5 &
               +y1*2d0*y2*y3*y4**2*SQRT3/6.D0 &
               +y1*2d0*y2*y3*y5**2*SQRT3/2.D0 )*DFEA123345 &
               +( -y1**3*2d0*y2*y5*SQRT3/2.D0 -y1**3*2d0*y2*y4/2.D0 &
               -y1**2*3d0*y2**2*y4/2.D0 &
               +3d0*y2**2*y3**2*y4 -y1**2*3d0*y2**2*y5*SQRT3/2.D0 &
               +2d0*y2*y3**3*y4 )*DFEA222334 +( +2d0*y2*y5**4 +3.D0*2d0*y2*y4**4 &
               +4.D0*2d0*y2*y4*y5**3*SQRT3 )*DFEA335555 +( y1**3*3d0*y2**2 &
               +3d0*y2**2*y3**3 )*DFEA222333 
            dds4(2,6)=dds3(2,6) &
               +( -y4**4*y5*SQRT3/2.D0 -3.D0/10.D0*y5**5*SQRT3 +y4**3*y5**2 &
               +y4**5/5.D0 )*DFEA244455 &
               +( 5d0*y2**4*y4 -SQRT3*5d0*y2**4*y5 )*DFEA222224 
            dds5(2,6)=dds4(2,6) &
               +( +y5**5*SQRT3/5.D0 -7.D0/15.D0*y4**5 +y4**4*y5*SQRT3/3.D0 &
               +y4*y5**4 )*DFEA145555 
            dds1(2,6)=dds5(2,6) +( -SQRT3*y1*y3**3*y5/2.D0 &
               +y1**3*y3*y4 &
               +SQRT3*y1*3d0*y2**2*y3*y5/2.D0 -y1*3d0*y2**2*y3*y4/2.D0 &
               -y1*y3**3*y4/2.D0 )*DFEA111234 &
               +( -y4**4*y5/3.D0 -y4*y5**4*SQRT3/2.D0 +y4**5*SQRT3/18.D0 &
               +y4**2*y5**3 )*DFEA244555 &
               +( y1*y4**2*y5**2 -3.D0/4.D0*y3*y4**4 -y1*y5**4 &
               +5.D0/4.D0*y3*y5**4 -7.D0/2.D0*y3*y4**2*y5**2 &
               -2.D0*y1*y4**3*y5*SQRT3 )*DFEA124455 
            dds3(2,6)=dds1(2,6) &
               +( 6d0*y2**5 )*DFEA333333 +( y1*4d0*y2**3*y3 +y1**4*y3 &
               +y1*y3**4 )*DFEA111123 +y1**2*2d0*y2*y3**2*DFEA112233 +( &
               +4d0*y2**3*y4**2 +4d0*y2**3*y5**2 )*DFEA222255 
            dds4(2,6)=dds3(2,6) +( &
               +9.D0*y3*y4**2*y5**2 -3.D0/2.D0*y3*y5**4 +y1*y4**4 &
               +4.D0*y1*y4**3*y5*SQRT3 +3.D0*y1*y5**4 &
               +5.D0/2.D0*y3*y4**4 )*DFEA134444 +( &
               +5.D0/3.D0*y1*2d0*y2*y4**2*y5*SQRT3 -13.D0/3.D0*2d0*y2*y3*y4*y5**2 &
               -4.D0/3.D0*y3**2*y5**3*SQRT3 -7.D0/3.D0*y1**2*y4*y5**2 &
               +4.D0/3.D0*2d0*y2*y3*y5**3*SQRT3 +3.D0*y1**2*y4**3 +y3**2*y4**3 &
               +y1*2d0*y2*y5**3*SQRT3/3.D0 &
               +2d0*y2*y3*y4**3 -13.D0/3.D0*y3**2*y4*y5**2 &
               -4.D0/3.D0*y1**2*y4**2*y5*SQRT3 &
               -16.D0/3.D0*y1*2d0*y2*y4*y5**2 )*DFEA233444 
            dds5(2,6)=dds4(2,6) &
               +( +4.D0*y3**2*y5**3 &
               +4.D0*2d0*y2*y3*y4*y5**2*SQRT3 -2.D0*y1*2d0*y2*y5**3 &
               -6.D0*y1*2d0*y2*y4**2*y5 &
               +y1**2*y4*y5**2*SQRT3 -3.D0*y1**2*y4**3*SQRT3 -4.D0*2d0*y2*y3*y5**3 &
               +3.D0*y1**2*y4**2*y5 -y1**2*y5**3 +4.D0*y3**2*y4*y5**2*SQRT3 &
               +4.D0*y1*2d0*y2*y4*y5**2*SQRT3 )*DFEA113555 
            dds2(2,6)=dds5(2,6) &
               +( -2d0*y2*y4**3*y5 -2.D0/3.D0*2d0*y2*y4**4*SQRT3 &
               -3.D0*2d0*y2*y4*y5**3 )*DFEA334445 &
               +( +2.D0/3.D0*y1*y5**4*SQRT3 &
               +3.D0*y1*y4**3*y5 -7.D0/12.D0*y3*y5**4*SQRT3 &
               +3.D0/2.D0*y3*y4**2*y5**2*SQRT3 +y1*y4*y5**3 &
               +3.D0/4.D0*y3*y4**4*SQRT3 )*DFEA124555 +( &
               +2d0*y2*y4**2*y5**2 -2d0*y2*y4**4 &
               -2.D0*2d0*y2*y4*y5**3*SQRT3 )*DFEA334455 
            dds3(2,6)=dds2(2,6) &
               +0d0 +( y3**3*y4**2 +y3**3*y5**2 +y1*3d0*y2**2*y4**2 +y1**3*y4**2 &
               +y1*3d0*y2**2*y5**2 +y1**3*y5**2 +3d0*y2**2*y3*y4**2 &
               +3d0*y2**2*y3*y5**2 )*DFEA233344 &
               +( y1*3d0*y2**2*y5**2*SQRT3/6.D0 -3d0*y2**2*y3*y5**2*SQRT3/3.D0 &
               -y3**3*y5**2*SQRT3/3.D0 &
               +y1**3*y4*y5 -y1**3*y5**2*SQRT3/3.D0 -3d0*y2**2*y3*y4*y5 &
               +y3**3*y4*y5 -y1*3d0*y2**2*y4**2*SQRT3/2.D0 )*DFEA233345 &
               +( -3.D0*3d0*y2**2*y4*y5**2 +3d0*y2**2*y4**3 )*DFEA111444 &
               +( y1*3d0*y2**2*y3**2 +y1**3*2d0*y2*y3 +y1**2*3d0*y2**2*y3 &
               +y1*2d0*y2*y3**3 +y1**2*y3**3 &
               +y1**3*y3**2 )*DFEA111233 
            dds4(2,6)=dds3(2,6) +0d0 &
               +( -5.D0/3.D0*y1*2d0*y2*y4**2*y5*SQRT3 &
               +y1*2d0*y2*y4**3 -2.D0*y1**2*y4**3 -y1*2d0*y2*y5**3*SQRT3/3.D0 &
               +4.D0/3.D0*2d0*y2*y3*y4*y5**2 -4.D0/3.D0*2d0*y2*y3*y5**3*SQRT3 &
               +7.D0/3.D0*y1*2d0*y2*y4*y5**2 +4.D0/3.D0*y3**2*y5**3*SQRT3 &
               +4.D0/3.D0*y1**2*y4**2*y5*SQRT3 &
               +4.D0/3.D0*y3**2*y4*y5**2 &
               -2.D0/3.D0*y1**2*y4*y5**2 )*DFEA133444 
            dds5(2,6)=dds4(2,6) &
               +( -y1**3*y4*y5 +2.D0/3.D0*3d0*y2**2*y3*y5**2*SQRT3 &
               +y1**3*y5**2*SQRT3/6.D0 +y1*3d0*y2**2*y5**2*SQRT3/6.D0 &
               +y1**3*y4**2*SQRT3/2.D0 &
               +2.D0/3.D0*y3**3*y5**2*SQRT3 -y1*3d0*y2**2*y4*y5 &
               +y1*3d0*y2**2*y4**2*SQRT3/2.D0 )*DFEA133345 
            ddv6(2,6)=dds5(2,6) &
               +( -2d0*y2*y3*y4**2*y5 +y3**2*y4**2*y5 +y3**2*y5**3 -y1*2d0*y2*y5**3 &
               +4.D0/3.D0*2d0*y2*y3*y4*y5**2*SQRT3 &
               +4.D0/3.D0*y3**2*y4*y5**2*SQRT3 -y1*2d0*y2*y4**2*y5 -2d0*y2*y3*y5**3 &
               +y1**2*y4*y5**2*SQRT3/3.D0 -y1**2*y4**3*SQRT3 &
               +4.D0/3.D0*y1*2d0*y2*y4*y5**2*SQRT3 )*DFEA233445 &
               +( y3**4*y4*SQRT3 -y1**4*y5 +4d0*y2**3*y3*y4*SQRT3 &
               +y3**4*y5 -2.D0*y1*4d0*y2**3*y5 -y1**4*y4*SQRT3 &
               -4d0*y2**3*y3*y5 )*DFEA233335 &
               +( 2d0*y2*y3**4 +y1**2*4d0*y2**3 +4d0*y2**3*y3**2 &
               +y1**4*2d0*y2 )*DFEA222233 
 
            ddv(2,6)= +ddv1(2,6) +ddv2(2,6) &
               +ddv3(2,6) +ddv4(2,6) +ddv5(2,6) +ddv6(2,6)
 
            ddv1(3,3)=0d0 
            ddv2(3,3)=0d0 +( +2d0 )*FEA11 +0d0 +0d0 
            ddv3(3,3)=0d0 &
               +0d0 +( y2*2d0 +y1*2d0 )*FEA112 +0d0 +0d0 +( 3d0*2d0*y3 )*FEA111 &
               +( -2d0*y4/2.D0 -SQRT3*2d0*y5/2.D0 )*FEA114 
            dds2(3,3)=0d0 +0d0 +( &
               +3.D0/4.D0*2d0*y4**2 &
               +2d0*y5**2/4.D0 -SQRT3*2d0*y4*y5/2.D0 )*FEA1155 
            dds1(3,3)=dds2(3,3) &
               +( 2d0*y4**2/4.D0 +3.D0/4.D0*2d0*y5**2 +SQRT3*2d0*y4*y5/2.D0 )*FEA1144 &
               +( -SQRT3*3d0*2d0*y3*y5/2.D0 -3d0*2d0*y3*y4/2.D0 )*FEA1114 +( &
               +4d0*3d0*y3**2 )*FEA1111 +0d0 
            dds2(3,3)=dds1(3,3) +0d0 &
               +( -y1*2d0*y4/2.D0 -SQRT3*y1*2d0*y5/2.D0 -SQRT3*y2*2d0*y5/2.D0 &
               -y2*2d0*y4/2.D0 )*FEA1124 &
               +( +SQRT3*y1*2d0*y4/2.D0 -SQRT3*y2*2d0*y4/2.D0 &
               +y2*2d0*y5/2.D0 -y1*2d0*y5/2.D0 )*FEA1125 
            ddv4(3,3)=dds2(3,3) &
               +( y2*3d0*2d0*y3 +y1*3d0*2d0*y3 )*FEA1112 +( y2**2*2d0 &
               +y1**2*2d0 )*FEA1122 +( +y1*y2*2d0 )*FEA1123 +0d0 
            dds3(3,3)=0d0 +0d0 &
               +0d0 
            dds2(3,3)=dds3(3,3) +0d0 +( +2d0*y4**2*y5 +2d0*y5**3 &
               +4.D0/9.D0*2d0*y4**3*SQRT3 )*FEA33445 &
               +( 2d0*y4*y5**2 -2d0*y4**3/3.D0 )*FEA33455 
            dds1(3,3)=dds2(3,3) +( &
               +3d0*2d0*y3*y4*y5 +3d0*2d0*y3*y5**2*SQRT3/3.D0 )*FEA33345 &
               +( 3d0*2d0*y3*y4**2 +3d0*2d0*y3*y5**2 )*FEA33344 +( 4d0*3d0*y3**2*y4 &
               +SQRT3*4d0*3d0*y3**2*y5 )*FEA33334 +( +5d0*4d0*y3**3 )*FEA33333 +0d0 &
               +0d0 
            dds3(3,3)=dds1(3,3) +( +y2*2d0*y4**2 +y1*2d0*y4**2 +y1*2d0*y5**2 &
               +y2*2d0*y5**2 )*FEA11255 +( +y1*2d0*y5**2*SQRT3/2.D0 +y2*2d0*y4*y5 &
               +y1*2d0*y4*y5 +y2*2d0*y4**2*SQRT3/6.D0 +y1*2d0*y4**2*SQRT3/6.D0 &
               +y2*2d0*y5**2*SQRT3/2.D0 )*FEA13345 
            dds4(3,3)=dds3(3,3) +( &
               +y2*2d0*y4*y5 +y2*2d0*y4**2*SQRT3/3.D0 -y1*2d0*y4**2*SQRT3/6.D0 &
               +y1*2d0*y5**2*SQRT3/2.D0 )*FEA11245 
            dds2(3,3)=dds4(3,3) &
               +( -y2*3d0*2d0*y3*y5/2.D0 +y1*3d0*2d0*y3*y5/2.D0 &
               +y2*3d0*2d0*y3*y4*SQRT3/2.D0 -y1*3d0*2d0*y3*y4*SQRT3/2.D0 )*FEA11135 &
               +( -y2*3d0*2d0*y3*y4/2.D0 -y1*3d0*2d0*y3*y4/2.D0 &
               -y2*3d0*2d0*y3*y5*SQRT3/2.D0 &
               -y1*3d0*2d0*y3*y5*SQRT3/2.D0 )*FEA11134 
            ddv5(3,3)=dds2(3,3) &
               +( +y2*4d0*3d0*y3**2 +y1*4d0*3d0*y3**2 )*FEA23333 &
               +( -2.D0*y2**2*2d0*y4 -SQRT3*y1**2*2d0*y5 +y1**2*2d0*y4 )*FEA11334 &
               +( y1**2*3d0*2d0*y3 +y1**3*2d0 +y2**2*3d0*2d0*y3 +y2**3*2d0 )*FEA11333 &
               +0d0 +( -y1*y2*2d0*y4/2.D0 -SQRT3*y1*y2*2d0*y5/2.D0 )*FEA11234 +( &
               +y1*y2*3d0*2d0*y3 )*FEA11123 +( +y1*y2**2*2d0 &
               +y1**2*y2*2d0 )*FEA11233 
            dds3(3,3)=( +3d0*2d0*y3*y4**2*y5 &
               +3d0*2d0*y3*y4**3*SQRT3 -5.D0/3.D0*3d0*2d0*y3*y4*y5**2*SQRT3 &
               +3d0*2d0*y3*y5**3 )*FEA333555 +( &
               +4d0*3d0*y3**2*y4**2*SQRT3/3.D0 -4d0*3d0*y3**2*y4*y5 )*FEA222245 &
               +( y1*5d0*4d0*y3**3 +y2*5d0*4d0*y3**3 )*FEA133333 +( &
               +y1*4d0*3d0*y3**2*y4 -2.D0*y2*4d0*3d0*y3**2*y4 &
               -y1*4d0*3d0*y3**2*y5*SQRT3 )*FEA133334 &
               +0d0 
            dds4(3,3)=dds3(3,3) +( 2.D0/3.D0*SQRT3*y1*y2**2*2d0*y4 &
               +y1**2*y2*2d0*y5 -SQRT3*y1**2*y2*2d0*y4/3.D0 )*FEA112335 +( &
               +y1*y2*2d0*y5**2 +y1*y2*2d0*y4**2 )*FEA112355 
            dds2(3,3)=dds4(3,3) &
               +( y2**3*2d0*y5 -y1**2*3d0*2d0*y3*y5/2.D0 -y2**2*3d0*2d0*y3*y5 &
               +y1**3*2d0*y5/2.D0 &
               +y1**3*2d0*y4*SQRT3/2.D0 -y1**2*3d0*2d0*y3*y4*SQRT3/2.D0 )*FEA222335 &
               +( -y1**2*2d0*y5**2*SQRT3/2.D0 -2.D0/3.D0*y2**2*2d0*y4**2*SQRT3 &
               +y1**2*2d0*y4*y5 -y1**2*2d0*y4**2*SQRT3/6.D0 )*FEA113345 &
               +( y2**2*2d0*y5**2 +y2**2*2d0*y4**2 +y1**2*2d0*y4**2 &
               +y1**2*2d0*y5**2 )*FEA223355 
            dds3(3,3)=dds2(3,3) &
               +( y1*y2*2d0*y4**2*SQRT3/6.D0 +y1*y2*2d0*y4*y5 &
               +y1*y2*2d0*y5**2*SQRT3/2.D0 )*FEA123345 +( -y1**3*2d0*y4/2.D0 &
               +y1**3*2d0*y5*SQRT3/2.D0 -y1**2*3d0*2d0*y3*y4/2.D0 +y2**3*2d0*y4 &
               +y2**2*3d0*2d0*y3*y4 +y1**2*3d0*2d0*y3*y5*SQRT3/2.D0 )*FEA222334 &
               +( 3.D0*2d0*y4**4 -4.D0*2d0*y4*y5**3*SQRT3 +2d0*y5**4 )*FEA335555 +( &
               +y1**3*3d0*2d0*y3 +y2**3*3d0*2d0*y3 )*FEA222333 
            dds4(3,3)=dds3(3,3) &
               +0d0 +( +5d0*4d0*y3**3*y4 &
               +SQRT3*5d0*4d0*y3**3*y5 )*FEA222224 
            dds5(3,3)=dds4(3,3) &
               +0d0 
            dds1(3,3)=dds5(3,3) &
               +( -SQRT3*y1*y2*3d0*2d0*y3*y5/2.D0 &
               -y1*y2*3d0*2d0*y3*y4/2.D0 )*FEA111234 &
               +0d0 +0d0 
            dds3(3,3)=dds1(3,3) +( +6d0*5d0*y3**4 )*FEA333333 +( &
               +y1*y2*4d0*3d0*y3**2 )*FEA111123 +y1**2*y2**2*2d0*FEA112233 +( &
               +4d0*3d0*y3**2*y4**2 &
               +4d0*3d0*y3**2*y5**2 )*FEA222255 
            dds4(3,3)=dds3(3,3) +0d0 &
               +( -y1*2d0*y5**3*SQRT3/3.D0 -4.D0/3.D0*y2*2d0*y5**3*SQRT3 &
               -16.D0/3.D0*y1*2d0*y4*y5**2 &
               +y2*2d0*y4**3 -13.D0/3.D0*y2*2d0*y4*y5**2 &
               -5.D0/3.D0*y1*2d0*y4**2*y5*SQRT3 )*FEA233444 
            dds5(3,3)=dds4(3,3) &
               +( 2.D0*y1*2d0*y5**3 +4.D0*y2*2d0*y5**3 +6.D0*y1*2d0*y4**2*y5 &
               +4.D0*y1*2d0*y4*y5**2*SQRT3 &
               +4.D0*y2*2d0*y4*y5**2*SQRT3 )*FEA113555 
            dds2(3,3)=dds5(3,3) &
               +( -2.D0/3.D0*2d0*y4**4*SQRT3 +2d0*y4**3*y5 &
               +3.D0*2d0*y4*y5**3 )*FEA334445 +0d0 &
               +( 2.D0*2d0*y4*y5**3*SQRT3 -2d0*y4**4 &
               +2d0*y4**2*y5**2 )*FEA334455 
            dds3(3,3)=dds2(3,3) +0d0 &
               +( y2*3d0*2d0*y3*y4**2 +y2*3d0*2d0*y3*y5**2 +y1*3d0*2d0*y3*y4**2 &
               +y1*3d0*2d0*y3*y5**2 )*FEA233344 &
               +( -y2*3d0*2d0*y3*y5**2*SQRT3/3.D0 -y1*3d0*2d0*y3*y4**2*SQRT3/2.D0 &
               +y1*3d0*2d0*y3*y5**2*SQRT3/6.D0 +y2*3d0*2d0*y3*y4*y5 )*FEA233345 +( &
               +3d0*2d0*y3*y4**3 -3.D0*3d0*2d0*y3*y4*y5**2 )*FEA111444 &
               +( y1*y2**3*2d0 +y1*y2**2*3d0*2d0*y3 +y1**2*y2*3d0*2d0*y3 &
               +y1**3*y2*2d0 )*FEA111233 
            dds4(3,3)=dds3(3,3) +0d0 +( +y1*2d0*y4**3 &
               +4.D0/3.D0*y2*2d0*y5**3*SQRT3 +y1*2d0*y5**3*SQRT3/3.D0 &
               +4.D0/3.D0*y2*2d0*y4*y5**2 +5.D0/3.D0*y1*2d0*y4**2*y5*SQRT3 &
               +7.D0/3.D0*y1*2d0*y4*y5**2 )*FEA133444 
            dds5(3,3)=dds4(3,3) +( &
               +y1*3d0*2d0*y3*y4**2*SQRT3/2.D0 +2.D0/3.D0*y2*3d0*2d0*y3*y5**2*SQRT3 &
               +y1*3d0*2d0*y3*y5**2*SQRT3/6.D0 &
               +y1*3d0*2d0*y3*y4*y5 )*FEA133345 
            ddv6(3,3)=dds5(3,3) +( &
               +y2*2d0*y4**2*y5 +y2*2d0*y5**3 +4.D0/3.D0*y2*2d0*y4*y5**2*SQRT3 &
               +4.D0/3.D0*y1*2d0*y4*y5**2*SQRT3 +y1*2d0*y5**3 &
               +y1*2d0*y4**2*y5 )*FEA233445 +( y2*4d0*3d0*y3**2*y4*SQRT3 &
               +y2*4d0*3d0*y3**2*y5 +2.D0*y1*4d0*3d0*y3**2*y5 )*FEA233335 &
               +( y2**2*4d0*3d0*y3**2 +y1**4*2d0 +y2**4*2d0 &
               +y1**2*4d0*3d0*y3**2 )*FEA222233 
 
            ddv(3,3)= +ddv1(3,3) +ddv2(3,3) &
               +ddv3(3,3) +ddv4(3,3) +ddv5(3,3) +ddv6(3,3)
 
            ddv1(3,4)=0d0 
            ddv2(3,4)=0d0 +0d0 +( -1d0/2.D0 )*FEA14 &
               +0d0 
            ddv3(3,4)=( y1 -2.D0*y2 )*FEA124 &
               +( 3.D0/4.D0*2d0*y4 -SQRT3*y5/2.D0 )*FEA155 +0d0 +0d0 +( &
               +SQRT3*y5/2.D0 +2d0*y4/4.D0 )*FEA144 +0d0 &
               +( -2d0*y3*1d0/2.D0 )*FEA114 
            dds2(3,4)=0d0 &
               +( -3.D0/8.D0*SQRT3*2d0*y4*y5 -3d0*y4**2/8.D0 &
               -9.D0/8.D0*y5**2 )*FEA1444 &
               +( &
               +3.D0/4.D0*2d0*y3*2d0*y4 &
               -SQRT3*2d0*y3*y5/2.D0 )*FEA1155 
            dds1(3,4)=dds2(3,4) &
               +( 2d0*y3*2d0*y4/4.D0 +SQRT3*2d0*y3*y5/2.D0 )*FEA1144 &
               +( -3d0*y3**2*1d0/2.D0 )*FEA1114 +0d0 +( SQRT3*y1*y5 -y2*2d0*y4/2.D0 &
               +y1*2d0*y4 )*FEA1244 
            dds2(3,4)=dds1(3,4) +( -SQRT3*y1*y5 &
               +3.D0/2.D0*y2*2d0*y4 )*FEA1255 +( -y1*2d0*y3*1d0/2.D0 &
               +y1**2 -y2**2*1d0/2.D0 -y2*2d0*y3*1d0/2.D0 )*FEA1124 +( &
               +SQRT3*y1*2d0*y3*1d0/2.D0 -SQRT3*y2*2d0*y3*1d0/2.D0 &
               -SQRT3*y2**2*1d0/2.D0 )*FEA1125 
            ddv4(3,4)=dds2(3,4) &
               +0d0 +0d0 +0d0 +( -SQRT3*2d0*y4*y5/8.D0 &
               +5.D0/8.D0*y5**2 -3.D0/8.D0*3d0*y4**2 )*FEA1455 
            dds3(3,4)=0d0 &
               +( -4.D0*y5**3*SQRT3 +3.D0*4d0*y4**3 )*FEA25555 +( &
               +2d0*y4*y5**2 -4d0*y4**3 &
               +2.D0*y5**3*SQRT3 )*FEA24455 
            dds2(3,4)=dds3(3,4) +( -3.D0*y5**3 &
               +2.D0/3.D0*4d0*y4**3*SQRT3 -3d0*y4**2*y5 )*FEA24445 +( &
               +2d0*y3*2d0*y4*y5 +4.D0/9.D0*2d0*y3*3d0*y4**2*SQRT3 )*FEA33445 &
               +( 2d0*y3*y5**2 -2d0*y3*3d0*y4**2/3.D0 )*FEA33455 
            dds1(3,4)=dds2(3,4) &
               +( +3d0*y3**2*y5 )*FEA33345 +( 3d0*y3**2*2d0*y4 )*FEA33344 &
               +( 4d0*y3**3 )*FEA33334 +0d0 +( +y1*2d0*y4*y5 +y2*y5**2*SQRT3 &
               +5.D0/9.D0*y2*3d0*y4**2*SQRT3 -4.D0/9.D0*y1*3d0*y4**2*SQRT3 )*FEA13445 &
               +( y2*y5**2 -y2*3d0*y4**2/3.D0 -y1*3d0*y4**2/3.D0 &
               +y1*y5**2 )*FEA13455 
            dds3(3,4)=dds1(3,4) +( +y2**2*2d0*y4 &
               +y2*2d0*y3*2d0*y4 +y1*2d0*y3*2d0*y4 +y1**2*2d0*y4 )*FEA11255 &
               +( 2.D0/3.D0*y1**2*2d0*y4*SQRT3 +y2*2d0*y3*y5 +y1*2d0*y3*y5 -y2**2*y5 &
               +y2*2d0*y3*2d0*y4*SQRT3/6.D0 +y1*2d0*y3*2d0*y4*SQRT3/6.D0 &
               +y2**2*2d0*y4*SQRT3/6.D0 )*FEA13345 
            dds4(3,4)=dds3(3,4) +( &
               +y1**2*2d0*y4*SQRT3/3.D0 +y2*2d0*y3*y5 -y2**2*y5 -y1**2*y5 &
               +y2*2d0*y3*2d0*y4*SQRT3/3.D0 -y1*2d0*y3*2d0*y4*SQRT3/6.D0 &
               +y2**2*2d0*y4*SQRT3/3.D0 )*FEA11245 
            dds2(3,4)=dds4(3,4) +( &
               +y2**3*SQRT3/2.D0 &
               +y2*3d0*y3**2*SQRT3/2.D0 -y1*3d0*y3**2*SQRT3/2.D0 )*FEA11135 &
               +( y1**3 -y2**3*1d0/2.D0 -y2*3d0*y3**2*1d0/2.D0 &
               -y1*3d0*y3**2*1d0/2.D0 )*FEA11134 
            ddv5(3,4)=dds2(3,4) &
               +0d0 +( -2.D0*y2**2*2d0*y3 +y1**2*2d0*y3 )*FEA11334 +0d0 &
               +( y1*y2*2d0*y4 )*FEA12355 &
               +( -y1*y2*2d0*y3*1d0/2.D0 -y1*y2**2*1d0/2.D0 +y1**2*y2 )*FEA11234 +0d0 &
               +0d0 
            dds3(3,4)=( +3d0*y3**2*2d0*y4*y5 &
               +3d0*y3**2*3d0*y4**2*SQRT3 -5.D0/3.D0*3d0*y3**2*y5**2*SQRT3 )*FEA333555 &
               +( +4d0*y3**3*2d0*y4*SQRT3/3.D0 -4d0*y3**3*y5 )*FEA222245 +0d0 &
               +( y1**4 -2.D0*y2**4 +y1*4d0*y3**3 -2.D0*y2*4d0*y3**3 )*FEA133334 &
               +( -y1*y2*3d0*y4**2/3.D0 +y1*y2*y5**2 )*FEA123455 
            dds4(3,4)=dds3(3,4) &
               +( 2.D0/3.D0*SQRT3*y1*y2**2*2d0*y3 -SQRT3*y1**2*y2**2*1d0/3.D0 &
               -SQRT3*y1**2*y2*2d0*y3*1d0/3.D0 )*FEA112335 &
               +( +y1*y2*2d0*y3*2d0*y4 +y1*y2**2*2d0*y4 &
               +y1**2*y2*2d0*y4 )*FEA112355 
            dds2(3,4)=dds4(3,4) +( &
               +y1**3*2d0*y3*SQRT3/2.D0 -y1**2*3d0*y3**2*SQRT3/2.D0 )*FEA222335 &
               +( -2.D0/3.D0*y2**2*2d0*y3*2d0*y4*SQRT3 &
               +y1**2*2d0*y3*y5 -y1**2*2d0*y3*2d0*y4*SQRT3/6.D0 )*FEA113345 +( &
               +y2**2*2d0*y3*2d0*y4 &
               +y1**2*2d0*y3*2d0*y4 )*FEA223355 
            dds3(3,4)=dds2(3,4) &
               +( y1*y2*2d0*y3*2d0*y4*SQRT3/6.D0 +y1*y2*2d0*y3*y5 &
               +2.D0/3.D0*y1**2*y2*2d0*y4*SQRT3 -y1*y2**2*y5 &
               +y1*y2**2*2d0*y4*SQRT3/6.D0 )*FEA123345 &
               +( -y1**3*2d0*y3*1d0/2.D0 -y1**2*3d0*y3**2*1d0/2.D0 +y2**3*2d0*y3 &
               +y2**2*3d0*y3**2 )*FEA222334 &
               +( 3.D0*2d0*y3*4d0*y4**3 -4.D0*2d0*y3*y5**3*SQRT3 )*FEA335555 &
               +0d0 
            dds4(3,4)=dds3(3,4) +( 5d0*y4**4/5.D0 +3d0*y4**2*y5**2 &
               +4d0*y4**3*y5*SQRT3/2.D0 )*FEA244455 +( &
               +5d0*y3**4 )*FEA222224 
            dds5(3,4)=dds4(3,4) +( -4d0*y4**3*y5*SQRT3/3.D0 &
               +y5**4 -7.D0/15.D0*5d0*y4**4 )*FEA145555 
            dds1(3,4)=dds5(3,4) +( &
               +y1**3*y2 -y1*y2**3*1d0/2.D0 -y1*y2*3d0*y3**2*1d0/2.D0 )*FEA111234 &
               +( 4d0*y4**3*y5/3.D0 &
               +5d0*y4**4*SQRT3/18.D0 -2d0*y4*y5**3 -y5**4*SQRT3/2.D0 )*FEA244555 &
               +( -3.D0/4.D0*y2*4d0*y4**3 +y1*2d0*y4*y5**2 -7.D0/2.D0*y2*2d0*y4*y5**2 &
               +2.D0*y1*3d0*y4**2*y5*SQRT3 )*FEA124455 
            dds3(3,4)=dds1(3,4) +0d0 +0d0 &
               +( +4d0*y3**3*2d0*y4 )*FEA222255 
            dds4(3,4)=dds3(3,4) +( +y1*4d0*y4**3 &
               +9.D0*y2*2d0*y4*y5**2 -4.D0*y1*3d0*y4**2*y5*SQRT3 &
               +5.D0/2.D0*y2*4d0*y4**3 )*FEA134444 &
               +( -7.D0/3.D0*y1**2*y5**2 -13.D0/3.D0*y2**2*y5**2 &
               -16.D0/3.D0*y1*2d0*y3*y5**2 &
               +4.D0/3.D0*y1**2*2d0*y4*y5*SQRT3 +y2*2d0*y3*3d0*y4**2 &
               +y2**2*3d0*y4**2 -13.D0/3.D0*y2*2d0*y3*y5**2 &
               -5.D0/3.D0*y1*2d0*y3*2d0*y4*y5*SQRT3 &
               +3.D0*y1**2*3d0*y4**2 )*FEA233444 
            dds5(3,4)=dds4(3,4) +( &
               +4.D0*y2**2*y5**2*SQRT3 +y1**2*y5**2*SQRT3 &
               +6.D0*y1*2d0*y3*2d0*y4*y5 -3.D0*y1**2*2d0*y4*y5 &
               +4.D0*y1*2d0*y3*y5**2*SQRT3 -3.D0*y1**2*3d0*y4**2*SQRT3 &
               +4.D0*y2*2d0*y3*y5**2*SQRT3 )*FEA113555 
            dds2(3,4)=dds5(3,4) &
               +( -2.D0/3.D0*2d0*y3*4d0*y4**3*SQRT3 +2d0*y3*3d0*y4**2*y5 &
               +3.D0*2d0*y3*y5**3 )*FEA334445 +( -3.D0*y1*3d0*y4**2*y5 -y1*y5**3 &
               +3.D0/2.D0*y2*2d0*y4*y5**2*SQRT3 &
               +3.D0/4.D0*y2*4d0*y4**3*SQRT3 )*FEA124555 &
               +( 2.D0*2d0*y3*y5**3*SQRT3 -2d0*y3*4d0*y4**3 &
               +2d0*y3*2d0*y4*y5**2 )*FEA334455 
            dds3(3,4)=dds2(3,4) +0d0 &
               +( y2*3d0*y3**2*2d0*y4 +y1*3d0*y3**2*2d0*y4 +y1**3*2d0*y4 &
               +y2**3*2d0*y4 )*FEA233344 &
               +( -y1**3*y5 -y1*3d0*y3**2*2d0*y4*SQRT3/2.D0 -y2**3*y5 &
               +y2*3d0*y3**2*y5 )*FEA233345 +( &
               +3d0*y3**2*3d0*y4**2 -3.D0*3d0*y3**2*y5**2 )*FEA111444 &
               +0d0 
            dds4(3,4)=dds3(3,4) +0d0 +( -4.D0/3.D0*y1**2*2d0*y4*y5*SQRT3 &
               +4.D0/3.D0*y2**2*y5**2 -2.D0*y1**2*3d0*y4**2 -2.D0/3.D0*y1**2*y5**2 &
               +y1*2d0*y3*3d0*y4**2 +4.D0/3.D0*y2*2d0*y3*y5**2 &
               +5.D0/3.D0*y1*2d0*y3*2d0*y4*y5*SQRT3 &
               +7.D0/3.D0*y1*2d0*y3*y5**2 )*FEA133444 
            dds5(3,4)=dds4(3,4) +( &
               +y1*3d0*y3**2*2d0*y4*SQRT3/2.D0 +y1**3*2d0*y4*SQRT3/2.D0 +y1**3*y5 &
               +y1*3d0*y3**2*y5 )*FEA133345 
            ddv6(3,4)=dds5(3,4) +( -y2**2*2d0*y4*y5 &
               +y1**2*y5**2*SQRT3/3.D0 +y2*2d0*y3*2d0*y4*y5 &
               +4.D0/3.D0*y2**2*y5**2*SQRT3 +4.D0/3.D0*y2*2d0*y3*y5**2*SQRT3 &
               +4.D0/3.D0*y1*2d0*y3*y5**2*SQRT3 &
               +y1*2d0*y3*2d0*y4*y5 -y1**2*3d0*y4**2*SQRT3 )*FEA233445 &
               +( y2*4d0*y3**3*SQRT3 +y2**4*SQRT3 -y1**4*SQRT3 )*FEA233335 &
               +0d0 
 
            ddv(3,4)= +ddv1(3,4) +ddv2(3,4) +ddv3(3,4) +ddv4(3,4) &
               +ddv5(3,4) +ddv6(3,4)
 
            ddv1(3,5)=0d0 
            ddv2(3,5)=0d0 +0d0 +( -SQRT3*1d0/2.D0 )*FEA14 &
               +0d0 
            ddv3(3,5)=( -SQRT3*y1 )*FEA124 +( -SQRT3*y4*1d0/2.D0 &
               +2d0*y5/4.D0 )*FEA155 +0d0 +0d0 +( +3.D0/4.D0*2d0*y5 &
               +SQRT3*y4*1d0/2.D0 )*FEA144 +0d0 &
               +( -SQRT3*2d0*y3*1d0/2.D0 )*FEA114 
            dds2(3,5)=0d0 &
               +( -3.D0/8.D0*SQRT3*y4**2 -3.D0/8.D0*SQRT3*3d0*y5**2 &
               -9.D0/8.D0*y4*2d0*y5 )*FEA1444 &
               +( &
               +2d0*y3*2d0*y5/4.D0 &
               -SQRT3*2d0*y3*y4*1d0/2.D0 )*FEA1155 
            dds1(3,5)=dds2(3,5) &
               +( +3.D0/4.D0*2d0*y3*2d0*y5 +SQRT3*2d0*y3*y4*1d0/2.D0 )*FEA1144 &
               +( -SQRT3*3d0*y3**2*1d0/2.D0 )*FEA1114 +0d0 +( SQRT3*y1*y4 &
               +3.D0/2.D0*y2*2d0*y5 )*FEA1244 
            dds2(3,5)=dds1(3,5) &
               +( y1*2d0*y5 -SQRT3*y1*y4 -y2*2d0*y5/2.D0 )*FEA1255 &
               +( -SQRT3*y1*2d0*y3*1d0/2.D0 -SQRT3*y2*2d0*y3*1d0/2.D0 &
               +SQRT3*y2**2*1d0/2.D0 )*FEA1124 +( -y2**2*1d0/2.D0 &
               +y2*2d0*y3*1d0/2.D0 -y1*2d0*y3*1d0/2.D0 &
               -y1**2 )*FEA1125 
            ddv4(3,5)=dds2(3,5) &
               +0d0 +0d0 +0d0 +( -SQRT3*y4**2*1d0/8.D0 -SQRT3*3d0*y5**2/8.D0 &
               +5.D0/8.D0*y4*2d0*y5 )*FEA1455 
            dds3(3,5)=0d0 &
               +( -4.D0*y4*3d0*y5**2*SQRT3 +4d0*y5**3 )*FEA25555 +( +y4**2*2d0*y5 &
               +2.D0*y4*3d0*y5**2*SQRT3 )*FEA24455 
            dds2(3,5)=dds3(3,5) &
               +( -3.D0*y4*3d0*y5**2 -y4**3 )*FEA24445 +( +2d0*y3*y4**2 &
               +2d0*y3*3d0*y5**2 )*FEA33445 &
               +( 2d0*y3*y4*2d0*y5 )*FEA33455 
            dds1(3,5)=dds2(3,5) +( +3d0*y3**2*y4 &
               +3d0*y3**2*2d0*y5*SQRT3/3.D0 )*FEA33345 +( &
               +3d0*y3**2*2d0*y5 )*FEA33344 +( +SQRT3*4d0*y3**3 )*FEA33334 +0d0 +( &
               +y1*y4**2 +y2*y4*2d0*y5*SQRT3 +y1*3d0*y5**2 )*FEA13445 +( y2*y4*2d0*y5 &
               +y1*y4*2d0*y5 )*FEA13455 
            dds3(3,5)=dds1(3,5) +( y1**2*2d0*y5 &
               +y2**2*2d0*y5 +y1*2d0*y3*2d0*y5 +y2*2d0*y3*2d0*y5 )*FEA11255 +( &
               +y1*2d0*y3*2d0*y5*SQRT3/2.D0 +y2**2*2d0*y5*SQRT3/2.D0 +y2*2d0*y3*y4 &
               +y1*2d0*y3*y4 -y2**2*y4 &
               +y2*2d0*y3*2d0*y5*SQRT3/2.D0 )*FEA13345 
            dds4(3,5)=dds3(3,5) +( &
               +y2*2d0*y3*y4 -y2**2*y4 -y1**2*y4 &
               +y1*2d0*y3*2d0*y5*SQRT3/2.D0 )*FEA11245 
            dds2(3,5)=dds4(3,5) +( +y1**3 &
               +y2**3*1d0/2.D0 -y2*3d0*y3**2*1d0/2.D0 &
               +y1*3d0*y3**2*1d0/2.D0 )*FEA11135 +( &
               +y2**3*SQRT3/2.D0 -y2*3d0*y3**2*SQRT3/2.D0 &
               -y1*3d0*y3**2*SQRT3/2.D0 )*FEA11134 
            ddv5(3,5)=dds2(3,5) &
               +0d0 +( -SQRT3*y1**2*2d0*y3 )*FEA11334 +0d0 +( &
               +y1*y2*2d0*y5 )*FEA12355 +( -SQRT3*y1*y2*2d0*y3*1d0/2.D0 &
               +SQRT3*y1*y2**2*1d0/2.D0 )*FEA11234 +0d0 +0d0 
            dds3(3,5)=( &
               +3d0*y3**2*y4**2 -5.D0/3.D0*3d0*y3**2*y4*2d0*y5*SQRT3 &
               +3d0*y3**2*3d0*y5**2 )*FEA333555 +( -4d0*y3**3*y4 )*FEA222245 +0d0 &
               +( -y1*4d0*y3**3*SQRT3 -y1**4*SQRT3 )*FEA133334 +( &
               +y1*y2*y4*2d0*y5 )*FEA123455 
            dds4(3,5)=dds3(3,5) +( -y1**2*y2**2 &
               +y1**2*y2*2d0*y3 )*FEA112335 +( y1*y2**2*2d0*y5 +y1*y2*2d0*y3*2d0*y5 &
               +y1**2*y2*2d0*y5 )*FEA112355 
            dds2(3,5)=dds4(3,5) &
               +( y2**3*2d0*y3 -y1**2*3d0*y3**2*1d0/2.D0 -y2**2*3d0*y3**2 &
               +y1**3*2d0*y3*1d0/2.D0 )*FEA222335 +( -y1**2*2d0*y3*2d0*y5*SQRT3/2.D0 &
               +y1**2*2d0*y3*y4 )*FEA113345 +( y2**2*2d0*y3*2d0*y5 &
               +y1**2*2d0*y3*2d0*y5 )*FEA223355 
            dds3(3,5)=dds2(3,5) +( &
               +y1*y2*2d0*y3*y4 +y1*y2*2d0*y3*2d0*y5*SQRT3/2.D0 -y1*y2**2*y4 &
               +y1*y2**2*2d0*y5*SQRT3/2.D0 )*FEA123345 +( +y1**3*2d0*y3*SQRT3/2.D0 &
               +y1**2*3d0*y3**2*SQRT3/2.D0 )*FEA222334 &
               +( -4.D0*2d0*y3*y4*3d0*y5**2*SQRT3 +2d0*y3*4d0*y5**3 )*FEA335555 &
               +0d0 
            dds4(3,5)=dds3(3,5) +( +y4**3*2d0*y5 +y4**4*SQRT3/2.D0 &
               +3.D0/10.D0*5d0*y5**4*SQRT3 )*FEA244455 +( &
               +SQRT3*5d0*y3**4 )*FEA222224 
            dds5(3,5)=dds4(3,5) &
               +( -5d0*y5**4*SQRT3/5.D0 -y4**4*SQRT3/3.D0 &
               +y4*4d0*y5**3 )*FEA145555 
            dds1(3,5)=dds5(3,5) &
               +( -SQRT3*y1*y2*3d0*y3**2*1d0/2.D0 &
               +SQRT3*y1*y2**3*1d0/2.D0 )*FEA111234 &
               +( y4**4*1d0/3.D0 -y4**2*3d0*y5**2 -y4*4d0*y5**3*SQRT3/2.D0 )*FEA244555 &
               +( -y1*4d0*y5**3 +5.D0/4.D0*y2*4d0*y5**3 &
               +y1*y4**2*2d0*y5 -7.D0/2.D0*y2*y4**2*2d0*y5 &
               +2.D0*y1*y4**3*SQRT3 )*FEA124455 
            dds3(3,5)=dds1(3,5) +0d0 +0d0 +( &
               +4d0*y3**3*2d0*y5 )*FEA222255 
            dds4(3,5)=dds3(3,5) +( 3.D0*y1*4d0*y5**3 &
               +9.D0*y2*y4**2*2d0*y5 -3.D0/2.D0*y2*4d0*y5**3 &
               -4.D0*y1*y4**3*SQRT3 )*FEA134444 &
               +( -y1*2d0*y3*3d0*y5**2*SQRT3/3.D0 -7.D0/3.D0*y1**2*y4*2d0*y5 &
               -13.D0/3.D0*y2**2*y4*2d0*y5 -4.D0/3.D0*y2*2d0*y3*3d0*y5**2*SQRT3 &
               -16.D0/3.D0*y1*2d0*y3*y4*2d0*y5 &
               +4.D0/3.D0*y1**2*y4**2*SQRT3 &
               +4.D0/3.D0*y2**2*3d0*y5**2*SQRT3 -13.D0/3.D0*y2*2d0*y3*y4*2d0*y5 &
               -5.D0/3.D0*y1*2d0*y3*y4**2*SQRT3 )*FEA233444 
            dds5(3,5)=dds4(3,5) &
               +( 2.D0*y1*2d0*y3*3d0*y5**2 +4.D0*y2*2d0*y3*3d0*y5**2 &
               +4.D0*y2**2*y4*2d0*y5*SQRT3 +y1**2*y4*2d0*y5*SQRT3 &
               +6.D0*y1*2d0*y3*y4**2 -3.D0*y1**2*y4**2 &
               +4.D0*y1*2d0*y3*y4*2d0*y5*SQRT3 -4.D0*y2**2*3d0*y5**2 +y1**2*3d0*y5**2 &
               +4.D0*y2*2d0*y3*y4*2d0*y5*SQRT3 )*FEA113555 
            dds2(3,5)=dds5(3,5) +( &
               +2d0*y3*y4**3 +3.D0*2d0*y3*y4*3d0*y5**2 )*FEA334445 &
               +( -3.D0*y1*y4**3 -y1*y4*3d0*y5**2 &
               +2.D0/3.D0*y1*4d0*y5**3*SQRT3 -7.D0/12.D0*y2*4d0*y5**3*SQRT3 &
               +3.D0/2.D0*y2*y4**2*2d0*y5*SQRT3 )*FEA124555 &
               +( 2.D0*2d0*y3*y4*3d0*y5**2*SQRT3 &
               +2d0*y3*y4**2*2d0*y5 )*FEA334455 
            dds3(3,5)=dds2(3,5) +0d0 +( &
               +y2*3d0*y3**2*2d0*y5 +y1**3*2d0*y5 +y1*3d0*y3**2*2d0*y5 &
               +y2**3*2d0*y5 )*FEA233344 &
               +( -y2**3*2d0*y5*SQRT3/3.D0 -y2*3d0*y3**2*2d0*y5*SQRT3/3.D0 -y1**3*y4 &
               -y1**3*2d0*y5*SQRT3/3.D0 &
               +y1*3d0*y3**2*2d0*y5*SQRT3/6.D0 -y2**3*y4 +y2*3d0*y3**2*y4 )*FEA233345 &
               +( -3.D0*3d0*y3**2*y4*2d0*y5 )*FEA111444 +0d0 
            dds4(3,5)=dds3(3,5) +0d0 &
               +( -4.D0/3.D0*y1**2*y4**2*SQRT3 &
               +4.D0/3.D0*y2**2*y4*2d0*y5 -4.D0/3.D0*y2**2*3d0*y5**2*SQRT3 &
               -2.D0/3.D0*y1**2*y4*2d0*y5 &
               +4.D0/3.D0*y2*2d0*y3*3d0*y5**2*SQRT3 +y1*2d0*y3*3d0*y5**2*SQRT3/3.D0 &
               +4.D0/3.D0*y2*2d0*y3*y4*2d0*y5 +5.D0/3.D0*y1*2d0*y3*y4**2*SQRT3 &
               +7.D0/3.D0*y1*2d0*y3*y4*2d0*y5 )*FEA133444 
            dds5(3,5)=dds4(3,5) +( &
               +2.D0/3.D0*y2**3*2d0*y5*SQRT3 +y1**3*2d0*y5*SQRT3/6.D0 +y1**3*y4 &
               +2.D0/3.D0*y2*3d0*y3**2*2d0*y5*SQRT3 +y1*3d0*y3**2*2d0*y5*SQRT3/6.D0 &
               +y1*3d0*y3**2*y4 )*FEA133345 
            ddv6(3,5)=dds5(3,5) +( -y2**2*y4**2 &
               +y1**2*y4*2d0*y5*SQRT3/3.D0 +y2*2d0*y3*y4**2 +y2*2d0*y3*3d0*y5**2 &
               +4.D0/3.D0*y2**2*y4*2d0*y5*SQRT3 +4.D0/3.D0*y2*2d0*y3*y4*2d0*y5*SQRT3 &
               +4.D0/3.D0*y1*2d0*y3*y4*2d0*y5*SQRT3 -y2**2*3d0*y5**2 &
               +y1*2d0*y3*3d0*y5**2 +y1*2d0*y3*y4**2 )*FEA233445 +( +y2*4d0*y3**3 &
               +2.D0*y1*4d0*y3**3 +y1**4 -y2**4 )*FEA233335 +0d0 
 
            ddv(3,5)= &
               +ddv1(3,5) +ddv2(3,5) +ddv3(3,5) +ddv4(3,5) +ddv5(3,5) +ddv6(3,5)
 
            ddv1(3,6)=( 1d0 )*DFEA1 
            ddv2(3,6)=( y2 +y1 )*DFEA12 +( &
               +2d0*y3 )*DFEA11 +( -SQRT3*y5/2.D0 -y4/2.D0 )*DFEA14 &
               +0d0 
            ddv3(3,6)=( y1*y4 -2.D0*y2*y4 -SQRT3*y1*y5 )*DFEA124 &
               +( 3.D0/4.D0*y4**2 -SQRT3*y4*y5/2.D0 +y5**2/4.D0 )*DFEA155 &
               +( y2*2d0*y3 +y1*2d0*y3 +y1**2 +y2**2 )*DFEA112 +0d0 +y1*y2*DFEA123 +( &
               +3.D0/4.D0*y5**2 +SQRT3*y4*y5/2.D0 +y4**2/4.D0 )*DFEA144 &
               +( 3d0*y3**2 )*DFEA111 &
               +( -2d0*y3*y4/2.D0 -SQRT3*2d0*y3*y5/2.D0 )*DFEA114 
            dds2(3,6)=0d0 &
               +( -3.D0/8.D0*SQRT3*y4**2*y5 -3.D0/8.D0*SQRT3*y5**3 -y4**3/8.D0 &
               -9.D0/8.D0*y4*y5**2 )*DFEA1444 &
               +( +3.D0/4.D0*2d0*y3*y4**2 &
               +2d0*y3*y5**2/4.D0 &
               -SQRT3*2d0*y3*y4*y5/2.D0 )*DFEA1155 
            dds1(3,6)=dds2(3,6) &
               +( 2d0*y3*y4**2/4.D0 +3.D0/4.D0*2d0*y3*y5**2 &
               +SQRT3*2d0*y3*y4*y5/2.D0 )*DFEA1144 &
               +( -SQRT3*3d0*y3**2*y5/2.D0 -3d0*y3**2*y4/2.D0 )*DFEA1114 +( &
               +4d0*y3**3 )*DFEA1111 +( SQRT3*y1*y4*y5 &
               +3.D0/2.D0*y2*y5**2 -y2*y4**2/2.D0 &
               +y1*y4**2 )*DFEA1244 
            dds2(3,6)=dds1(3,6) &
               +( y1*y5**2 -SQRT3*y1*y4*y5 -y2*y5**2/2.D0 &
               +3.D0/2.D0*y2*y4**2 )*DFEA1255 +( -y1*2d0*y3*y4/2.D0 &
               +y1**2*y4 -SQRT3*y1*2d0*y3*y5/2.D0 -SQRT3*y2*2d0*y3*y5/2.D0 &
               +SQRT3*y2**2*y5/2.D0 -y2**2*y4/2.D0 -y2*2d0*y3*y4/2.D0 )*DFEA1124 +( &
               +SQRT3*y1*2d0*y3*y4/2.D0 -SQRT3*y2*2d0*y3*y4/2.D0 -SQRT3*y2**2*y4/2.D0 &
               -y2**2*y5/2.D0 &
               +y2*2d0*y3*y5/2.D0 -y1*2d0*y3*y5/2.D0 &
               -y1**2*y5 )*DFEA1125 
            ddv4(3,6)=dds2(3,6) &
               +( y2*3d0*y3**2 +y1**3 +y1*3d0*y3**2 +y2**3 )*DFEA1112 +( y2**2*2d0*y3 &
               +y1**2*2d0*y3 )*DFEA1122 +( y1*y2**2 +y1**2*y2 &
               +y1*y2*2d0*y3 )*DFEA1123 +( -SQRT3*y4**2*y5/8.D0 -SQRT3*y5**3/8.D0 &
               +5.D0/8.D0*y4*y5**2 -3.D0/8.D0*y4**3 )*DFEA1455 
            dds3(3,6)=0d0 &
               +( -4.D0*y4*y5**3*SQRT3 +3.D0*y4**4 +y5**4 )*DFEA25555 +( &
               +y4**2*y5**2 -y4**4 &
               +2.D0*y4*y5**3*SQRT3 )*DFEA24455 
            dds2(3,6)=dds3(3,6) +( -3.D0*y4*y5**3 &
               +2.D0/3.D0*y4**4*SQRT3 -y4**3*y5 )*DFEA24445 +( +2d0*y3*y4**2*y5 &
               +2d0*y3*y5**3 +4.D0/9.D0*2d0*y3*y4**3*SQRT3 )*DFEA33445 &
               +( 2d0*y3*y4*y5**2 -2d0*y3*y4**3/3.D0 )*DFEA33455 
            dds1(3,6)=dds2(3,6) &
               +( +3d0*y3**2*y4*y5 +3d0*y3**2*y5**2*SQRT3/3.D0 )*DFEA33345 &
               +( 3d0*y3**2*y4**2 +3d0*y3**2*y5**2 )*DFEA33344 +( 4d0*y3**3*y4 &
               +SQRT3*4d0*y3**3*y5 )*DFEA33334 +( +5d0*y3**4 )*DFEA33333 +( &
               +y1*y4**2*y5 +y2*y4*y5**2*SQRT3 &
               +5.D0/9.D0*y2*y4**3*SQRT3 -4.D0/9.D0*y1*y4**3*SQRT3 &
               +y1*y5**3 )*DFEA13445 +( y2*y4*y5**2 -y2*y4**3/3.D0 -y1*y4**3/3.D0 &
               +y1*y4*y5**2 )*DFEA13455 
            dds3(3,6)=dds1(3,6) +( y1**2*y5**2 &
               +y2**2*y4**2 +y2**2*y5**2 +y2*2d0*y3*y4**2 +y1*2d0*y3*y4**2 &
               +y1**2*y4**2 +y1*2d0*y3*y5**2 +y2*2d0*y3*y5**2 )*DFEA11255 &
               +( 2.D0/3.D0*y1**2*y4**2*SQRT3 +y1*2d0*y3*y5**2*SQRT3/2.D0 &
               +y2**2*y5**2*SQRT3/2.D0 +y2*2d0*y3*y4*y5 +y1*2d0*y3*y4*y5 -y2**2*y4*y5 &
               +y2*2d0*y3*y4**2*SQRT3/6.D0 +y1*2d0*y3*y4**2*SQRT3/6.D0 &
               +y2*2d0*y3*y5**2*SQRT3/2.D0 &
               +y2**2*y4**2*SQRT3/6.D0 )*DFEA13345 
            dds4(3,6)=dds3(3,6) +( &
               +y1**2*y4**2*SQRT3/3.D0 +y2*2d0*y3*y4*y5 -y2**2*y4*y5 -y1**2*y4*y5 &
               +y2*2d0*y3*y4**2*SQRT3/3.D0 -y1*2d0*y3*y4**2*SQRT3/6.D0 &
               +y2**2*y4**2*SQRT3/3.D0 &
               +y1*2d0*y3*y5**2*SQRT3/2.D0 )*DFEA11245 
            dds2(3,6)=dds4(3,6) +( &
               +y1**3*y5 +y2**3*y5/2.D0 -y2*3d0*y3**2*y5/2.D0 +y1*3d0*y3**2*y5/2.D0 &
               +y2**3*y4*SQRT3/2.D0 &
               +y2*3d0*y3**2*y4*SQRT3/2.D0 -y1*3d0*y3**2*y4*SQRT3/2.D0 )*DFEA11135 &
               +( y1**3*y4 -y2**3*y4/2.D0 -y2*3d0*y3**2*y4/2.D0 -y1*3d0*y3**2*y4/2.D0 &
               +y2**3*y5*SQRT3/2.D0 -y2*3d0*y3**2*y5*SQRT3/2.D0 &
               -y1*3d0*y3**2*y5*SQRT3/2.D0 )*DFEA11134 
            ddv5(3,6)=dds2(3,6) &
               +( +y1**4 +y2**4 +y2*4d0*y3**3 +y1*4d0*y3**3 )*DFEA23333 &
               +( -2.D0*y2**2*2d0*y3*y4 -SQRT3*y1**2*2d0*y3*y5 &
               +y1**2*2d0*y3*y4 )*DFEA11334 +( y1**2*3d0*y3**2 +y1**3*2d0*y3 &
               +y2**2*3d0*y3**2 +y2**3*2d0*y3 )*DFEA11333 +( y1*y2*y4**2 &
               +y1*y2*y5**2 )*DFEA12355 &
               +( -y1*y2*2d0*y3*y4/2.D0 -y1*y2**2*y4/2.D0 -SQRT3*y1*y2*2d0*y3*y5/2.D0 &
               +y1**2*y2*y4 +SQRT3*y1*y2**2*y5/2.D0 )*DFEA11234 +( y1*y2**3 &
               +y1*y2*3d0*y3**2 +y1**3*y2 )*DFEA11123 +( y1**2*y2**2 +y1*y2**2*2d0*y3 &
               +y1**2*y2*2d0*y3 )*DFEA11233 
            dds3(3,6)=( +3d0*y3**2*y4**2*y5 &
               +3d0*y3**2*y4**3*SQRT3 -5.D0/3.D0*3d0*y3**2*y4*y5**2*SQRT3 &
               +3d0*y3**2*y5**3 )*DFEA333555 +( &
               +4d0*y3**3*y4**2*SQRT3/3.D0 -4d0*y3**3*y4*y5 )*DFEA222245 &
               +( y1*5d0*y3**4 +y2**5 +y1**5 +y2*5d0*y3**4 )*DFEA133333 &
               +( y1**4*y4 -2.D0*y2**4*y4 &
               +y1*4d0*y3**3*y4 -2.D0*y2*4d0*y3**3*y4 -y1*4d0*y3**3*y5*SQRT3 &
               -y1**4*y5*SQRT3 )*DFEA133334 &
               +( -y1*y2*y4**3/3.D0 +y1*y2*y4*y5**2 )*DFEA123455 
            dds4(3,6)=dds3(3,6) &
               +( 2.D0/3.D0*SQRT3*y1*y2**2*2d0*y3*y4 -y1**2*y2**2*y5 &
               -SQRT3*y1**2*y2**2*y4/3.D0 &
               +y1**2*y2*2d0*y3*y5 -SQRT3*y1**2*y2*2d0*y3*y4/3.D0 )*DFEA112335 &
               +( y1*y2**2*y5**2 +y1*y2*2d0*y3*y5**2 +y1*y2*2d0*y3*y4**2 &
               +y1*y2**2*y4**2 +y1**2*y2*y4**2 &
               +y1**2*y2*y5**2 )*DFEA112355 
            dds2(3,6)=dds4(3,6) &
               +( y2**3*2d0*y3*y5 -y1**2*3d0*y3**2*y5/2.D0 -y2**2*3d0*y3**2*y5 &
               +y1**3*2d0*y3*y5/2.D0 &
               +y1**3*2d0*y3*y4*SQRT3/2.D0 -y1**2*3d0*y3**2*y4*SQRT3/2.D0 )*DFEA222335 &
               +( -y1**2*2d0*y3*y5**2*SQRT3/2.D0 -2.D0/3.D0*y2**2*2d0*y3*y4**2*SQRT3 &
               +y1**2*2d0*y3*y4*y5 -y1**2*2d0*y3*y4**2*SQRT3/6.D0 )*DFEA113345 &
               +( y2**2*2d0*y3*y5**2 +y2**2*2d0*y3*y4**2 +y1**2*2d0*y3*y4**2 &
               +y1**2*2d0*y3*y5**2 )*DFEA223355 
            dds3(3,6)=dds2(3,6) &
               +( y1*y2*2d0*y3*y4**2*SQRT3/6.D0 +y1*y2*2d0*y3*y4*y5 &
               +y1*y2*2d0*y3*y5**2*SQRT3/2.D0 &
               +2.D0/3.D0*y1**2*y2*y4**2*SQRT3 -y1*y2**2*y4*y5 &
               +y1*y2**2*y4**2*SQRT3/6.D0 +y1*y2**2*y5**2*SQRT3/2.D0 )*DFEA123345 &
               +( -y1**3*2d0*y3*y4/2.D0 &
               +y1**3*2d0*y3*y5*SQRT3/2.D0 -y1**2*3d0*y3**2*y4/2.D0 +y2**3*2d0*y3*y4 &
               +y2**2*3d0*y3**2*y4 +y1**2*3d0*y3**2*y5*SQRT3/2.D0 )*DFEA222334 &
               +( 3.D0*2d0*y3*y4**4 -4.D0*2d0*y3*y4*y5**3*SQRT3 &
               +2d0*y3*y5**4 )*DFEA335555 +( +y1**3*3d0*y3**2 &
               +y2**3*3d0*y3**2 )*DFEA222333 
            dds4(3,6)=dds3(3,6) +( y4**5/5.D0 &
               +y4**3*y5**2 +y4**4*y5*SQRT3/2.D0 +3.D0/10.D0*y5**5*SQRT3 )*DFEA244455 &
               +( +5d0*y3**4*y4 +SQRT3*5d0*y3**4*y5 )*DFEA222224 
            dds5(3,6)=dds4(3,6) &
               +( -y5**5*SQRT3/5.D0 -y4**4*y5*SQRT3/3.D0 &
               +y4*y5**4 -7.D0/15.D0*y4**5 )*DFEA145555 
            dds1(3,6)=dds5(3,6) &
               +( -SQRT3*y1*y2*3d0*y3**2*y5/2.D0 +y1**3*y2*y4 &
               +SQRT3*y1*y2**3*y5/2.D0 -y1*y2**3*y4/2.D0 &
               -y1*y2*3d0*y3**2*y4/2.D0 )*DFEA111234 &
               +( y4**4*y5/3.D0 &
               +y4**5*SQRT3/18.D0 -y4**2*y5**3 -y4*y5**4*SQRT3/2.D0 )*DFEA244555 &
               +( -3.D0/4.D0*y2*y4**4 -y1*y5**4 +5.D0/4.D0*y2*y5**4 &
               +y1*y4**2*y5**2 -7.D0/2.D0*y2*y4**2*y5**2 &
               +2.D0*y1*y4**3*y5*SQRT3 )*DFEA124455 
            dds3(3,6)=dds1(3,6) +( &
               +6d0*y3**5 )*DFEA333333 +( y1*y2**4 +y1**4*y2 &
               +y1*y2*4d0*y3**3 )*DFEA111123 +y1**2*y2**2*2d0*y3*DFEA112233 +( &
               +4d0*y3**3*y4**2 +4d0*y3**3*y5**2 )*DFEA222255 
            dds4(3,6)=dds3(3,6) &
               +( 3.D0*y1*y5**4 +y1*y4**4 &
               +9.D0*y2*y4**2*y5**2 -3.D0/2.D0*y2*y5**4 -4.D0*y1*y4**3*y5*SQRT3 &
               +5.D0/2.D0*y2*y4**4 )*DFEA134444 &
               +( -y1*2d0*y3*y5**3*SQRT3/3.D0 -7.D0/3.D0*y1**2*y4*y5**2 &
               -13.D0/3.D0*y2**2*y4*y5**2 -4.D0/3.D0*y2*2d0*y3*y5**3*SQRT3 &
               -16.D0/3.D0*y1*2d0*y3*y4*y5**2 &
               +4.D0/3.D0*y1**2*y4**2*y5*SQRT3 +4.D0/3.D0*y2**2*y5**3*SQRT3 &
               +y2*2d0*y3*y4**3 &
               +y2**2*y4**3 -13.D0/3.D0*y2*2d0*y3*y4*y5**2 &
               -5.D0/3.D0*y1*2d0*y3*y4**2*y5*SQRT3 &
               +3.D0*y1**2*y4**3 )*DFEA233444 
            dds5(3,6)=dds4(3,6) &
               +( 2.D0*y1*2d0*y3*y5**3 +4.D0*y2*2d0*y3*y5**3 &
               +4.D0*y2**2*y4*y5**2*SQRT3 +y1**2*y4*y5**2*SQRT3 &
               +6.D0*y1*2d0*y3*y4**2*y5 -3.D0*y1**2*y4**2*y5 &
               +4.D0*y1*2d0*y3*y4*y5**2*SQRT3 -4.D0*y2**2*y5**3 &
               +y1**2*y5**3 -3.D0*y1**2*y4**3*SQRT3 &
               +4.D0*y2*2d0*y3*y4*y5**2*SQRT3 )*DFEA113555 
            dds2(3,6)=dds5(3,6) &
               +( -2.D0/3.D0*2d0*y3*y4**4*SQRT3 +2d0*y3*y4**3*y5 &
               +3.D0*2d0*y3*y4*y5**3 )*DFEA334445 +( -3.D0*y1*y4**3*y5 -y1*y4*y5**3 &
               +2.D0/3.D0*y1*y5**4*SQRT3 -7.D0/12.D0*y2*y5**4*SQRT3 &
               +3.D0/2.D0*y2*y4**2*y5**2*SQRT3 +3.D0/4.D0*y2*y4**4*SQRT3 )*DFEA124555 &
               +( 2.D0*2d0*y3*y4*y5**3*SQRT3 -2d0*y3*y4**4 &
               +2d0*y3*y4**2*y5**2 )*DFEA334455 
            dds3(3,6)=dds2(3,6) +0d0 &
               +( y2*3d0*y3**2*y4**2 +y2*3d0*y3**2*y5**2 +y1*3d0*y3**2*y4**2 &
               +y1**3*y5**2 +y1**3*y4**2 +y2**3*y4**2 +y1*3d0*y3**2*y5**2 &
               +y2**3*y5**2 )*DFEA233344 &
               +( -y2**3*y5**2*SQRT3/3.D0 -y2*3d0*y3**2*y5**2*SQRT3/3.D0 -y1**3*y4*y5 &
               -y1**3*y5**2*SQRT3/3.D0 -y1*3d0*y3**2*y4**2*SQRT3/2.D0 &
               +y1*3d0*y3**2*y5**2*SQRT3/6.D0 -y2**3*y4*y5 &
               +y2*3d0*y3**2*y4*y5 )*DFEA233345 +( &
               +3d0*y3**2*y4**3 -3.D0*3d0*y3**2*y4*y5**2 )*DFEA111444 &
               +( y1*y2**3*2d0*y3 +y1**3*y2**2 +y1**2*y2**3 +y1*y2**2*3d0*y3**2 &
               +y1**2*y2*3d0*y3**2 +y1**3*y2*2d0*y3 )*DFEA111233 
            dds4(3,6)=dds3(3,6) &
               +0d0 +( -4.D0/3.D0*y1**2*y4**2*y5*SQRT3 &
               +4.D0/3.D0*y2**2*y4*y5**2 -4.D0/3.D0*y2**2*y5**3*SQRT3 &
               -2.D0*y1**2*y4**3 -2.D0/3.D0*y1**2*y4*y5**2 &
               +y1*2d0*y3*y4**3 +4.D0/3.D0*y2*2d0*y3*y5**3*SQRT3 &
               +y1*2d0*y3*y5**3*SQRT3/3.D0 +4.D0/3.D0*y2*2d0*y3*y4*y5**2 &
               +5.D0/3.D0*y1*2d0*y3*y4**2*y5*SQRT3 &
               +7.D0/3.D0*y1*2d0*y3*y4*y5**2 )*DFEA133444 
            dds5(3,6)=dds4(3,6) +( &
               +2.D0/3.D0*y2**3*y5**2*SQRT3 +y1*3d0*y3**2*y4**2*SQRT3/2.D0 &
               +y1**3*y4**2*SQRT3/2.D0 +y1**3*y5**2*SQRT3/6.D0 +y1**3*y4*y5 &
               +2.D0/3.D0*y2*3d0*y3**2*y5**2*SQRT3 +y1*3d0*y3**2*y5**2*SQRT3/6.D0 &
               +y1*3d0*y3**2*y4*y5 )*DFEA133345 
            ddv6(3,6)=dds5(3,6) &
               +( -y2**2*y4**2*y5 +y1**2*y4*y5**2*SQRT3/3.D0 +y2*2d0*y3*y4**2*y5 &
               +y2*2d0*y3*y5**3 +4.D0/3.D0*y2**2*y4*y5**2*SQRT3 &
               +4.D0/3.D0*y2*2d0*y3*y4*y5**2*SQRT3 &
               +4.D0/3.D0*y1*2d0*y3*y4*y5**2*SQRT3 -y2**2*y5**3 +y1*2d0*y3*y5**3 &
               +y1*2d0*y3*y4**2*y5 -y1**2*y4**3*SQRT3 )*DFEA233445 &
               +( y2*4d0*y3**3*y4*SQRT3 +y2**4*y4*SQRT3 -y1**4*y4*SQRT3 &
               +y2*4d0*y3**3*y5 +2.D0*y1*4d0*y3**3*y5 &
               +y1**4*y5 -y2**4*y5 )*DFEA233335 +( y2**2*4d0*y3**3 +y1**4*2d0*y3 &
               +y2**4*2d0*y3 +y1**2*4d0*y3**3 )*DFEA222233 
 
            ddv(3,6)= +ddv1(3,6) &
               +ddv2(3,6) +ddv3(3,6) +ddv4(3,6) +ddv5(3,6) +ddv6(3,6)
 
            ddv1(4,4)=0d0 
            ddv2(4,4)=0d0 +0d0 +0d0 +( +2d0 )*FEA44 
            ddv3(4,4)=0d0 &
               +( 3.D0/4.D0*y3*2d0 +3.D0/4.D0*y2*2d0 )*FEA155 +0d0 &
               +( -3d0*2d0*y4/3.D0 )*FEA455 +( y1*2d0 +y2*2d0*1d0/4.D0 &
               +y3*2d0*1d0/4.D0 )*FEA144 +0d0 +0d0 
            dds2(4,4)=( 4d0*3d0*y4**2 &
               +2.D0*2d0*y5**2 )*FEA4444 &
               +( -3.D0/8.D0*SQRT3*y3*2d0*y5 -y3*3d0*2d0*y4/8.D0 -y2*3d0*2d0*y4/8.D0 &
               +y1*3d0*2d0*y4 +3.D0/8.D0*SQRT3*y2*2d0*y5 )*FEA1444 &
               +( 3.D0/4.D0*y2**2*2d0 &
               +3.D0/4.D0*y3**2*2d0 )*FEA1155 
            dds1(4,4)=dds2(4,4) &
               +( y3**2*2d0*1d0/4.D0 +y1**2*2d0 +y2**2*2d0*1d0/4.D0 )*FEA1144 +0d0 &
               +0d0 +( -y2*y3*2d0*1d0/2.D0 +y1*y2*2d0 &
               +y1*y3*2d0 )*FEA1244 
            dds2(4,4)=dds1(4,4) +( &
               +3.D0/2.D0*y2*y3*2d0 )*FEA1255 +0d0 +0d0 
            ddv4(4,4)=dds2(4,4) +0d0 +0d0 &
               +0d0 +( -SQRT3*y3*2d0*y5/8.D0 &
               +SQRT3*y2*2d0*y5/8.D0 -3.D0/8.D0*y2*3d0*2d0*y4 &
               -3.D0/8.D0*y3*3d0*2d0*y4 )*FEA1455 
            dds3(4,4)=( 5d0*4d0*y4**3 &
               -2.D0*3d0*2d0*y4*y5**2 )*FEA44444 &
               +( +9.D0*y1*2d0*y5**2 -3.D0/2.D0*y1*4d0*3d0*y4**2 &
               +3.D0*y2*4d0*3d0*y4**2 +3.D0*y3*4d0*3d0*y4**2 )*FEA25555 &
               +( -y2*4d0*3d0*y4**2 &
               +y3*2d0*y5**2 -y3*4d0*3d0*y4**2 -7.D0/2.D0*y1*2d0*y5**2 +y2*2d0*y5**2 &
               +5.D0/4.D0*y1*4d0*3d0*y4**2 )*FEA24455 
            dds2(4,4)=dds3(4,4) &
               +( y2*3d0*2d0*y4*y5 &
               +2.D0/3.D0*y3*4d0*3d0*y4**2*SQRT3 -7.D0/12.D0*y1*4d0*3d0*y4**2*SQRT3 &
               +3.D0/2.D0*y1*2d0*y5**2*SQRT3 -y3*3d0*2d0*y4*y5 &
               +2.D0/3.D0*y2*4d0*3d0*y4**2*SQRT3 )*FEA24445 +( +y3**2*2d0*y5 &
               +4.D0/9.D0*y2**2*3d0*2d0*y4*SQRT3 -5.D0/9.D0*y1**2*3d0*2d0*y4*SQRT3 &
               +4.D0/9.D0*y3**2*3d0*2d0*y4*SQRT3 -y2**2*2d0*y5 )*FEA33445 &
               +( -y1**2*3d0*2d0*y4/3.D0 -y3**2*3d0*2d0*y4/3.D0 &
               -y2**2*3d0*2d0*y4/3.D0 )*FEA33455 
            dds1(4,4)=dds2(4,4) &
               +( +y1**3*2d0*SQRT3/2.D0 )*FEA33345 +( y3**3*2d0 +y2**3*2d0 &
               +y1**3*2d0 )*FEA33344 +0d0 +0d0 +( -4.D0/9.D0*y1*y2*3d0*2d0*y4*SQRT3 &
               +y1*y3*2d0*y5 -y1*y2*2d0*y5 &
               +5.D0/9.D0*y2*y3*3d0*2d0*y4*SQRT3 &
               -4.D0/9.D0*y1*y3*3d0*2d0*y4*SQRT3 )*FEA13445 &
               +( -y2*y3*3d0*2d0*y4/3.D0 -y1*y2*3d0*2d0*y4/3.D0 &
               -y1*y3*3d0*2d0*y4/3.D0 )*FEA13455 
            dds3(4,4)=dds1(4,4) &
               +( +y2**2*y3*2d0 +y1*y2**2*2d0 +y2*y3**2*2d0 +y1*y3**2*2d0 &
               +y1**2*y3*2d0 +y1**2*y2*2d0 )*FEA11255 +( 2.D0/3.D0*y1**2*y3*2d0*SQRT3 &
               +y2*y3**2*2d0*SQRT3/6.D0 +y1*y3**2*2d0*SQRT3/6.D0 &
               +y1*y2**2*2d0*SQRT3/6.D0 +2.D0/3.D0*y1**2*y2*2d0*SQRT3 &
               +y2**2*y3*2d0*SQRT3/6.D0 )*FEA13345 
            dds4(4,4)=dds3(4,4) +( &
               +y1**2*y3*2d0*SQRT3/3.D0 &
               +y1**2*y2*2d0*SQRT3/3.D0 -y1*y2**2*2d0*SQRT3/6.D0 &
               +y2*y3**2*2d0*SQRT3/3.D0 -y1*y3**2*2d0*SQRT3/6.D0 &
               +y2**2*y3*2d0*SQRT3/3.D0 )*FEA11245 
            dds2(4,4)=dds4(4,4) +0d0 &
               +0d0 
            ddv5(4,4)=dds2(4,4) +0d0 +0d0 +0d0 +( y1*y2*y3*2d0 )*FEA12355 &
               +0d0 +0d0 +0d0 
            dds3(4,4)=( y2**3*3d0*2d0*y4*SQRT3 -y2**3*2d0*y5 &
               +y3**3*2d0*y5 +y3**3*3d0*2d0*y4*SQRT3 )*FEA333555 +( &
               +y2**4*2d0*SQRT3/3.D0 &
               +y3**4*2d0*SQRT3/3.D0 -y1**4*2d0*SQRT3/6.D0 )*FEA222245 +0d0 +0d0 &
               +( -y1*y2*y3*3d0*2d0*y4/3.D0 )*FEA123455 
            dds4(4,4)=dds3(4,4) +0d0 +( &
               +y1*y2*y3**2*2d0 +y1*y2**2*y3*2d0 &
               +y1**2*y2*y3*2d0 )*FEA112355 
            dds2(4,4)=dds4(4,4) +0d0 &
               +( -y1**2*y2**2*2d0*SQRT3/6.D0 -2.D0/3.D0*y2**2*y3**2*2d0*SQRT3 &
               -y1**2*y3**2*2d0*SQRT3/6.D0 )*FEA113345 &
               +( +y2**2*y3**2*2d0 +y1**2*y3**2*2d0 &
               +y1**2*y2**2*2d0 )*FEA223355 
            dds3(4,4)=dds2(4,4) &
               +( y1*y2*y3**2*2d0*SQRT3/6.D0 +2.D0/3.D0*y1**2*y2*y3*2d0*SQRT3 &
               +y1*y2**2*y3*2d0*SQRT3/6.D0 )*FEA123345 +0d0 &
               +( 3.D0*y3**2*4d0*3d0*y4**2 +3.D0*y2**2*4d0*3d0*y4**2 &
               +9.D0*y1**2*2d0*y5**2 -3.D0/2.D0*y1**2*4d0*3d0*y4**2 )*FEA335555 &
               +0d0 
            dds4(4,4)=dds3(4,4) &
               +( y3*5d0*4d0*y4**3/5.D0 -y2*4d0*3d0*y4**2*y5*SQRT3/2.D0 &
               -2.D0/5.D0*y1*5d0*4d0*y4**3 -2.D0*y1*3d0*2d0*y4*y5**2 &
               +y3*3d0*2d0*y4*y5**2 +y3*4d0*3d0*y4**2*y5*SQRT3/2.D0 &
               +y2*3d0*2d0*y4*y5**2 +y2*5d0*4d0*y4**3/5.D0 )*FEA244455 &
               +0d0 
            dds5(4,4)=dds4(4,4) +( -7.D0/15.D0*y2*5d0*4d0*y4**3 &
               +y2*4d0*3d0*y4**2*y5*SQRT3/3.D0 -y3*4d0*3d0*y4**2*y5*SQRT3/3.D0 &
               +2.D0*y1*3d0*2d0*y4*y5**2 -7.D0/15.D0*y3*5d0*4d0*y4**3 &
               -y1*5d0*4d0*y4**3/15.D0 )*FEA145555 
            dds1(4,4)=dds5(4,4) &
               +0d0 +( y3*4d0*3d0*y4**2*y5/3.D0 &
               +y3*5d0*4d0*y4**3*SQRT3/18.D0 -y2*4d0*3d0*y4**2*y5/3.D0 -y3*2d0*y5**3 &
               +2.D0/9.D0*y1*5d0*4d0*y4**3*SQRT3 +y2*5d0*4d0*y4**3*SQRT3/18.D0 &
               +y2*2d0*y5**3 -2.D0/3.D0*y1*3d0*2d0*y4*y5**2*SQRT3 )*FEA244555 &
               +( y1*y2*2d0*y5**2 -3.D0/4.D0*y2*y3*4d0*3d0*y4**2 &
               +y1*y3*2d0*y5**2 -7.D0/2.D0*y2*y3*2d0*y5**2 &
               -2.D0*y1*y2*3d0*2d0*y4*y5*SQRT3 &
               +2.D0*y1*y3*3d0*2d0*y4*y5*SQRT3 )*FEA124455 
            dds3(4,4)=dds1(4,4) +0d0 &
               +0d0 +( y1**4*2d0 +y2**4*2d0 &
               +y3**4*2d0 )*FEA222255 
            dds4(4,4)=dds3(4,4) +( +y1*y3*4d0*3d0*y4**2 &
               +9.D0*y2*y3*2d0*y5**2 -4.D0*y1*y3*3d0*2d0*y4*y5*SQRT3 &
               +y1*y2*4d0*3d0*y4**2 +4.D0*y1*y2*3d0*2d0*y4*y5*SQRT3 &
               +5.D0/2.D0*y2*y3*4d0*3d0*y4**2 )*FEA134444 +( &
               +5.D0/3.D0*y1*y2**2*2d0*y5*SQRT3 +4.D0/3.D0*y1**2*y3*2d0*y5*SQRT3 &
               +3.D0*y1**2*y2*3d0*2d0*y4 +y2*y3**2*3d0*2d0*y4 &
               +y2**2*y3*3d0*2d0*y4 -5.D0/3.D0*y1*y3**2*2d0*y5*SQRT3 &
               -4.D0/3.D0*y1**2*y2*2d0*y5*SQRT3 &
               +3.D0*y1**2*y3*3d0*2d0*y4 )*FEA233444 
            dds5(4,4)=dds4(4,4) +( &
               +6.D0*y1*y3**2*2d0*y5 -6.D0*y1*y2**2*2d0*y5 -3.D0*y1**2*y3*2d0*y5 &
               -3.D0*y1**2*y2*3d0*2d0*y4*SQRT3 &
               +3.D0*y1**2*y2*2d0*y5 &
               -3.D0*y1**2*y3*3d0*2d0*y4*SQRT3 )*FEA113555 
            dds2(4,4)=dds5(4,4) &
               +( -2.D0/3.D0*y3**2*4d0*3d0*y4**2*SQRT3 &
               -3.D0/2.D0*y1**2*2d0*y5**2*SQRT3 -y2**2*3d0*2d0*y4*y5 &
               +7.D0/12.D0*y1**2*4d0*3d0*y4**2*SQRT3 &
               +y3**2*3d0*2d0*y4*y5 -2.D0/3.D0*y2**2*4d0*3d0*y4**2*SQRT3 )*FEA334445 &
               +( -3.D0*y1*y3*3d0*2d0*y4*y5 +3.D0*y1*y2*3d0*2d0*y4*y5 &
               +3.D0/2.D0*y2*y3*2d0*y5**2*SQRT3 &
               +3.D0/4.D0*y2*y3*4d0*3d0*y4**2*SQRT3 )*FEA124555 &
               +( -7.D0/2.D0*y1**2*2d0*y5**2 &
               +y2**2*2d0*y5**2 -y2**2*4d0*3d0*y4**2 -y3**2*4d0*3d0*y4**2 &
               +5.D0/4.D0*y1**2*4d0*3d0*y4**2 &
               +y3**2*2d0*y5**2 )*FEA334455 
            dds3(4,4)=dds2(4,4) +( -6.D0*2d0*y5**4 &
               +9.D0*4d0*3d0*y4**2*y5**2 )*FEA555555 +( y2*y3**3*2d0 +y1*y3**3*2d0 &
               +y1*y2**3*2d0 +y1**3*y2*2d0 +y1**3*y3*2d0 +y2**3*y3*2d0 )*FEA233344 &
               +( -y1*y3**3*2d0*SQRT3/2.D0 -y1*y2**3*2d0*SQRT3/2.D0 )*FEA233345 +( &
               +y3**3*3d0*2d0*y4 +y2**3*3d0*2d0*y4 +y1**3*3d0*2d0*y4 )*FEA111444 &
               +0d0 
            dds4(4,4)=dds3(4,4) +( 9.D0*2d0*y5**4 -6.D0*4d0*3d0*y4**2*y5**2 &
               +6d0*5d0*y4**4 )*FEA444444 +( -5.D0/3.D0*y1*y2**2*2d0*y5*SQRT3 &
               +y1*y2**2*3d0*2d0*y4 -4.D0/3.D0*y1**2*y3*2d0*y5*SQRT3 &
               -2.D0*y1**2*y2*3d0*2d0*y4 -2.D0*y1**2*y3*3d0*2d0*y4 &
               +y1*y3**2*3d0*2d0*y4 +4.D0/3.D0*y1**2*y2*2d0*y5*SQRT3 &
               +5.D0/3.D0*y1*y3**2*2d0*y5*SQRT3 )*FEA133444 
            dds5(4,4)=dds4(4,4) +( &
               +y1*y3**3*2d0*SQRT3/2.D0 +y1**3*y3*2d0*SQRT3/2.D0 &
               +y1**3*y2*2d0*SQRT3/2.D0 &
               +y1*y2**3*2d0*SQRT3/2.D0 )*FEA133345 
            ddv6(4,4)=dds5(4,4) &
               +( -y2**2*y3*2d0*y5 &
               +y2*y3**2*2d0*y5 -y1*y2**2*2d0*y5 -y1**2*y2*3d0*2d0*y4*SQRT3 &
               +y1*y3**2*2d0*y5 -y1**2*y3*3d0*2d0*y4*SQRT3 )*FEA233445 +0d0 &
               +0d0 
 
            ddv(4,4)= +ddv1(4,4) +ddv2(4,4) +ddv3(4,4) +ddv4(4,4) &
               +ddv5(4,4) +ddv6(4,4)
 
            ddv1(4,5)=0d0 
            ddv2(4,5)=0d0 +0d0 +0d0 +0d0 
            ddv3(4,5)=0d0 &
               +( -SQRT3*y3*1d0/2.D0 +SQRT3*y2*1d0/2.D0 )*FEA155 +0d0 +( &
               +2d0*y5 )*FEA455 +( -SQRT3*y2*1d0/2.D0 +SQRT3*y3*1d0/2.D0 )*FEA144 &
               +0d0 +0d0 
            dds2(4,5)=( +2.D0*2d0*y4*2d0*y5 )*FEA4444 &
               +( -3.D0/8.D0*SQRT3*y3*2d0*y4 -9.D0/8.D0*y2*2d0*y5 -9.D0/8.D0*y3*2d0*y5 &
               +3.D0/8.D0*SQRT3*y2*2d0*y4 )*FEA1444 +( -SQRT3*y3**2*1d0/2.D0 &
               +SQRT3*y2**2*1d0/2.D0 )*FEA1155 
            dds1(4,5)=dds2(4,5) +( &
               +SQRT3*y3**2*1d0/2.D0 -SQRT3*y2**2*1d0/2.D0 )*FEA1144 +0d0 +0d0 &
               +( SQRT3*y1*y3 -SQRT3*y1*y2 )*FEA1244 
            dds2(4,5)=dds1(4,5) &
               +( -SQRT3*y1*y3 +SQRT3*y1*y2 )*FEA1255 +0d0 +0d0 
            ddv4(4,5)=dds2(4,5) &
               +0d0 +0d0 +0d0 +( 5.D0/8.D0*y2*2d0*y5 -SQRT3*y3*2d0*y4*1d0/8.D0 &
               +SQRT3*y2*2d0*y4*1d0/8.D0 +y1*2d0*y5 &
               +5.D0/8.D0*y3*2d0*y5 )*FEA1455 
            dds3(4,5)=( -2.D0*3d0*y4**2*2d0*y5 &
               -3.D0*4d0*y5**3 )*FEA44444 &
               +( -4.D0*y3*3d0*y5**2*SQRT3 +9.D0*y1*2d0*y4*2d0*y5 &
               +4.D0*y2*3d0*y5**2*SQRT3 )*FEA25555 +( &
               +y3*2d0*y4*2d0*y5 -2.D0*y2*3d0*y5**2*SQRT3 -7.D0/2.D0*y1*2d0*y4*2d0*y5 &
               +2.D0*y3*3d0*y5**2*SQRT3 &
               +y2*2d0*y4*2d0*y5 )*FEA24455 
            dds2(4,5)=dds3(4,5) &
               +( y2*3d0*y4**2 -3.D0*y3*3d0*y5**2 +3.D0*y2*3d0*y5**2 &
               +3.D0/2.D0*y1*2d0*y4*2d0*y5*SQRT3 -y3*3d0*y4**2 )*FEA24445 +( &
               +y3**2*2d0*y4 -y2**2*2d0*y4 -y1**2*2d0*y5*SQRT3 )*FEA33445 &
               +( y3**2*2d0*y5 +y1**2*2d0*y5 &
               +y2**2*2d0*y5 )*FEA33455 
            dds1(4,5)=dds2(4,5) +( -y2**3 &
               +y3**3 )*FEA33345 +0d0 +0d0 +0d0 +( +y1*y3*2d0*y4 &
               +y2*y3*2d0*y5*SQRT3 -y1*y2*2d0*y4 )*FEA13445 +( y2*y3*2d0*y5 &
               +y1*y2*2d0*y5 +y1*y3*2d0*y5 )*FEA13455 
            dds3(4,5)=dds1(4,5) +0d0 &
               +( -y1*y2**2 +y2*y3**2 &
               +y1*y3**2 -y2**2*y3 )*FEA13345 
            dds4(4,5)=dds3(4,5) +( y1**2*y2 &
               +y2*y3**2 -y2**2*y3 -y1**2*y3 )*FEA11245 
            dds2(4,5)=dds4(4,5) +0d0 &
               +0d0 
            ddv5(4,5)=dds2(4,5) +0d0 +0d0 +0d0 +0d0 +0d0 +0d0 &
               +0d0 
            dds3(4,5)=( -y2**3*2d0*y4 &
               +y3**3*2d0*y4 -5.D0/3.D0*y2**3*2d0*y5*SQRT3 &
               -5.D0/3.D0*y3**3*2d0*y5*SQRT3 &
               -8.D0/3.D0*y1**3*2d0*y5*SQRT3 )*FEA333555 &
               +( +y2**4 -y3**4 )*FEA222245 +0d0 +0d0 +( &
               +y1*y2*y3*2d0*y5 )*FEA123455 
            dds4(4,5)=dds3(4,5) +0d0 &
               +0d0 
            dds2(4,5)=dds4(4,5) +0d0 +( -y1**2*y2**2 +y1**2*y3**2 )*FEA113345 &
               +0d0 
            dds3(4,5)=dds2(4,5) +( +y1*y2*y3**2 -y1*y2**2*y3 )*FEA123345 +0d0 &
               +( -4.D0*y3**2*3d0*y5**2*SQRT3 +9.D0*y1**2*2d0*y4*2d0*y5 &
               +4.D0*y2**2*3d0*y5**2*SQRT3 )*FEA335555 +0d0 
            dds4(4,5)=dds3(4,5) &
               +( -y2*4d0*y4**3*SQRT3/2.D0 -2.D0*y1*3d0*y4**2*2d0*y5 &
               +y3*3d0*y4**2*2d0*y5 +y3*4d0*y4**3*SQRT3/2.D0 &
               +y2*3d0*y4**2*2d0*y5 )*FEA244455 +0d0 
            dds5(4,5)=dds4(4,5) +( &
               +y1*4d0*y5**3 +y2*4d0*y4**3*SQRT3/3.D0 -y3*4d0*y4**3*SQRT3/3.D0 &
               +y3*4d0*y5**3 +y2*4d0*y5**3 &
               +2.D0*y1*3d0*y4**2*2d0*y5 )*FEA145555 
            dds1(4,5)=dds5(4,5) +0d0 &
               +( y3*4d0*y4**3*1d0/3.D0 -y2*4d0*y4**3*1d0/3.D0 &
               -y2*4d0*y5**3*SQRT3/2.D0 -y3*2d0*y4*3d0*y5**2 &
               +y2*2d0*y4*3d0*y5**2 -2.D0/3.D0*y1*3d0*y4**2*2d0*y5*SQRT3 &
               -y3*4d0*y5**3*SQRT3/2.D0 )*FEA244555 &
               +( y1*y2*2d0*y4*2d0*y5 &
               +y1*y3*2d0*y4*2d0*y5 -7.D0/2.D0*y2*y3*2d0*y4*2d0*y5 &
               -2.D0*y1*y2*3d0*y4**2*SQRT3 &
               +2.D0*y1*y3*3d0*y4**2*SQRT3 )*FEA124455 
            dds3(4,5)=dds1(4,5) +0d0 +0d0 &
               +0d0 
            dds4(4,5)=dds3(4,5) +( &
               +9.D0*y2*y3*2d0*y4*2d0*y5 -4.D0*y1*y3*3d0*y4**2*SQRT3 &
               +4.D0*y1*y2*3d0*y4**2*SQRT3 )*FEA134444 +( -7.D0/3.D0*y1**2*y3*2d0*y5 &
               +5.D0/3.D0*y1*y2**2*2d0*y4*SQRT3 -13.D0/3.D0*y2**2*y3*2d0*y5 &
               -7.D0/3.D0*y1**2*y2*2d0*y5 -16.D0/3.D0*y1*y3**2*2d0*y5 &
               +4.D0/3.D0*y1**2*y3*2d0*y4*SQRT3 -13.D0/3.D0*y2*y3**2*2d0*y5 &
               -5.D0/3.D0*y1*y3**2*2d0*y4*SQRT3 -4.D0/3.D0*y1**2*y2*2d0*y4*SQRT3 &
               -16.D0/3.D0*y1*y2**2*2d0*y5 )*FEA233444 
            dds5(4,5)=dds4(4,5) &
               +( +4.D0*y2**2*y3*2d0*y5*SQRT3 +y1**2*y3*2d0*y5*SQRT3 &
               +6.D0*y1*y3**2*2d0*y4 -6.D0*y1*y2**2*2d0*y4 -3.D0*y1**2*y3*2d0*y4 &
               +y1**2*y2*2d0*y5*SQRT3 +4.D0*y1*y3**2*2d0*y5*SQRT3 &
               +3.D0*y1**2*y2*2d0*y4 +4.D0*y2*y3**2*2d0*y5*SQRT3 &
               +4.D0*y1*y2**2*2d0*y5*SQRT3 )*FEA113555 
            dds2(4,5)=dds5(4,5) &
               +( -3.D0/2.D0*y1**2*2d0*y4*2d0*y5*SQRT3 -y2**2*3d0*y4**2 &
               +y3**2*3d0*y4**2 &
               +3.D0*y3**2*3d0*y5**2 -3.D0*y2**2*3d0*y5**2 )*FEA334445 &
               +( -3.D0*y1*y3*3d0*y4**2 -y1*y3*3d0*y5**2 +3.D0*y1*y2*3d0*y4**2 &
               +3.D0/2.D0*y2*y3*2d0*y4*2d0*y5*SQRT3 +y1*y2*3d0*y5**2 )*FEA124555 &
               +( 2.D0*y3**2*3d0*y5**2*SQRT3 -7.D0/2.D0*y1**2*2d0*y4*2d0*y5 &
               +y2**2*2d0*y4*2d0*y5 -2.D0*y2**2*3d0*y5**2*SQRT3 &
               +y3**2*2d0*y4*2d0*y5 )*FEA334455 
            dds3(4,5)=dds2(4,5) &
               +( -6.D0*2d0*y4*4d0*y5**3 +9.D0*4d0*y4**3*2d0*y5 )*FEA555555 +0d0 +( &
               +y1**3*y2 -y1**3*y3 -y2**3*y3 +y2*y3**3 )*FEA233345 &
               +( -3.D0*y2**3*2d0*y5 -3.D0*y3**3*2d0*y5 -3.D0*y1**3*2d0*y5 )*FEA111444 &
               +0d0 
            dds4(4,5)=dds3(4,5) &
               +( 9.D0*2d0*y4*4d0*y5**3 -6.D0*4d0*y4**3*2d0*y5 )*FEA444444 &
               +( -5.D0/3.D0*y1*y2**2*2d0*y4*SQRT3 -4.D0/3.D0*y1**2*y3*2d0*y4*SQRT3 &
               +4.D0/3.D0*y2**2*y3*2d0*y5 &
               +7.D0/3.D0*y1*y2**2*2d0*y5 -2.D0/3.D0*y1**2*y3*2d0*y5 &
               +4.D0/3.D0*y1**2*y2*2d0*y4*SQRT3 +4.D0/3.D0*y2*y3**2*2d0*y5 &
               +5.D0/3.D0*y1*y3**2*2d0*y4*SQRT3 -2.D0/3.D0*y1**2*y2*2d0*y5 &
               +7.D0/3.D0*y1*y3**2*2d0*y5 )*FEA133444 
            dds5(4,5)=dds4(4,5) &
               +( -y1**3*y2 +y1**3*y3 -y1*y2**3 &
               +y1*y3**3 )*FEA133345 
            ddv6(4,5)=dds5(4,5) +( -y2**2*y3*2d0*y4 &
               +y1**2*y3*2d0*y5*SQRT3/3.D0 +y2*y3**2*2d0*y4 &
               +4.D0/3.D0*y2**2*y3*2d0*y5*SQRT3 &
               +4.D0/3.D0*y2*y3**2*2d0*y5*SQRT3 -y1*y2**2*2d0*y4 &
               +4.D0/3.D0*y1*y3**2*2d0*y5*SQRT3 +y1**2*y2*2d0*y5*SQRT3/3.D0 &
               +y1*y3**2*2d0*y4 +4.D0/3.D0*y1*y2**2*2d0*y5*SQRT3 )*FEA233445 +0d0 &
               +0d0 
 
            ddv(4,5)= +ddv1(4,5) +ddv2(4,5) +ddv3(4,5) +ddv4(4,5) &
               +ddv5(4,5) +ddv6(4,5)
 
            ddv1(4,6)=0d0 
            ddv2(4,6)=0d0 +0d0 +( -y3*1d0/2.D0 &
               +y1 -y2*1d0/2.D0 )*DFEA14 +( +2d0*y4 )*DFEA44 
            ddv3(4,6)=( y1*y3 &
               +y1*y2 -2.D0*y2*y3 )*DFEA124 +( 3.D0/4.D0*y3*2d0*y4 -SQRT3*y3*y5/2.D0 &
               +3.D0/4.D0*y2*2d0*y4 +SQRT3*y2*y5/2.D0 )*DFEA155 +0d0 &
               +( -3d0*y4**2/3.D0 +y5**2 )*DFEA455 +( y1*2d0*y4 &
               +y2*2d0*y4/4.D0 -SQRT3*y2*y5/2.D0 +SQRT3*y3*y5/2.D0 &
               +y3*2d0*y4/4.D0 )*DFEA144 +0d0 +( -y2**2*1d0/2.D0 -y3**2*1d0/2.D0 &
               +y1**2 )*DFEA114 
            dds2(4,6)=( 4d0*y4**3 +2.D0*2d0*y4*y5**2 )*DFEA4444 &
               +( -3.D0/8.D0*SQRT3*y3*2d0*y4*y5 -9.D0/8.D0*y2*y5**2 &
               -y3*3d0*y4**2/8.D0 -y2*3d0*y4**2/8.D0 -9.D0/8.D0*y3*y5**2 &
               +y1*3d0*y4**2 +3.D0/8.D0*SQRT3*y2*2d0*y4*y5 )*DFEA1444 &
               +( 3.D0/4.D0*y2**2*2d0*y4 +3.D0/4.D0*y3**2*2d0*y4 -SQRT3*y3**2*y5/2.D0 &
               +SQRT3*y2**2*y5/2.D0 )*DFEA1155 
            dds1(4,6)=dds2(4,6) &
               +( y3**2*2d0*y4/4.D0 +y1**2*2d0*y4 +y2**2*2d0*y4/4.D0 &
               +SQRT3*y3**2*y5/2.D0 -SQRT3*y2**2*y5/2.D0 )*DFEA1144 &
               +( y1**3 -y2**3*1d0/2.D0 -y3**3*1d0/2.D0 )*DFEA1114 +0d0 &
               +( SQRT3*y1*y3*y5 -y2*y3*2d0*y4/2.D0 +y1*y2*2d0*y4 -SQRT3*y1*y2*y5 &
               +y1*y3*2d0*y4 )*DFEA1244 
            dds2(4,6)=dds1(4,6) +( -SQRT3*y1*y3*y5 &
               +3.D0/2.D0*y2*y3*2d0*y4 +SQRT3*y1*y2*y5 )*DFEA1255 &
               +( -y1*y3**2*1d0/2.D0 +y1**2*y3 &
               +y1**2*y2 -y2**2*y3*1d0/2.D0 -y2*y3**2*1d0/2.D0 &
               -y1*y2**2*1d0/2.D0 )*DFEA1124 &
               +( +SQRT3*y1*y3**2*1d0/2.D0 &
               +SQRT3*y1*y2**2*1d0/2.D0 -SQRT3*y2*y3**2*1d0/2.D0 &
               -SQRT3*y2**2*y3*1d0/2.D0 )*DFEA1125 
            ddv4(4,6)=dds2(4,6) &
               +0d0 +0d0 +0d0 +( 5.D0/8.D0*y2*y5**2 -SQRT3*y3*2d0*y4*y5/8.D0 &
               +SQRT3*y2*2d0*y4*y5/8.D0 -3.D0/8.D0*y2*3d0*y4**2 +y1*y5**2 &
               +5.D0/8.D0*y3*y5**2 &
               -3.D0/8.D0*y3*3d0*y4**2 )*DFEA1455 
            dds3(4,6)=( 5d0*y4**4 &
               -2.D0*3d0*y4**2*y5**2 -3.D0*y5**4 )*DFEA44444 &
               +( -4.D0*y3*y5**3*SQRT3 +9.D0*y1*2d0*y4*y5**2 -3.D0/2.D0*y1*4d0*y4**3 &
               +4.D0*y2*y5**3*SQRT3 +3.D0*y2*4d0*y4**3 +3.D0*y3*4d0*y4**3 )*DFEA25555 &
               +( -y2*4d0*y4**3 &
               +y3*2d0*y4*y5**2 -2.D0*y2*y5**3*SQRT3 -y3*4d0*y4**3 &
               -7.D0/2.D0*y1*2d0*y4*y5**2 &
               +2.D0*y3*y5**3*SQRT3 +y2*2d0*y4*y5**2 &
               +5.D0/4.D0*y1*4d0*y4**3 )*DFEA24455 
            dds2(4,6)=dds3(4,6) &
               +( y2*3d0*y4**2*y5 -3.D0*y3*y5**3 +2.D0/3.D0*y3*4d0*y4**3*SQRT3 &
               +3.D0*y2*y5**3 -7.D0/12.D0*y1*4d0*y4**3*SQRT3 &
               +3.D0/2.D0*y1*2d0*y4*y5**2*SQRT3 -y3*3d0*y4**2*y5 &
               +2.D0/3.D0*y2*4d0*y4**3*SQRT3 )*DFEA24445 +( +y3**2*2d0*y4*y5 &
               +4.D0/9.D0*y2**2*3d0*y4**2*SQRT3 -5.D0/9.D0*y1**2*3d0*y4**2*SQRT3 &
               +4.D0/9.D0*y3**2*3d0*y4**2*SQRT3 -y2**2*2d0*y4*y5 &
               -y1**2*y5**2*SQRT3 )*DFEA33445 &
               +( y3**2*y5**2 -y1**2*3d0*y4**2/3.D0 -y3**2*3d0*y4**2/3.D0 &
               +y1**2*y5**2 &
               +y2**2*y5**2 -y2**2*3d0*y4**2/3.D0 )*DFEA33455 
            dds1(4,6)=dds2(4,6) &
               +( -y2**3*y5 +y3**3*y5 +y1**3*2d0*y4*SQRT3/2.D0 )*DFEA33345 &
               +( y3**3*2d0*y4 +y2**3*2d0*y4 +y1**3*2d0*y4 )*DFEA33344 +( y3**4 &
               +y2**4 -2.D0*y1**4 )*DFEA33334 +0d0 &
               +( -4.D0/9.D0*y1*y2*3d0*y4**2*SQRT3 +y1*y3*2d0*y4*y5 &
               +y2*y3*y5**2*SQRT3 -y1*y2*2d0*y4*y5 &
               +5.D0/9.D0*y2*y3*3d0*y4**2*SQRT3 &
               -4.D0/9.D0*y1*y3*3d0*y4**2*SQRT3 )*DFEA13445 &
               +( y2*y3*y5**2 &
               +y1*y2*y5**2 -y2*y3*3d0*y4**2/3.D0 -y1*y2*3d0*y4**2/3.D0 &
               -y1*y3*3d0*y4**2/3.D0 &
               +y1*y3*y5**2 )*DFEA13455 
            dds3(4,6)=dds1(4,6) +( +y2**2*y3*2d0*y4 &
               +y1*y2**2*2d0*y4 +y2*y3**2*2d0*y4 +y1*y3**2*2d0*y4 +y1**2*y3*2d0*y4 &
               +y1**2*y2*2d0*y4 )*DFEA11255 &
               +( 2.D0/3.D0*y1**2*y3*2d0*y4*SQRT3 -y1*y2**2*y5 +y2*y3**2*y5 &
               +y1*y3**2*y5 -y2**2*y3*y5 +y2*y3**2*2d0*y4*SQRT3/6.D0 &
               +y1*y3**2*2d0*y4*SQRT3/6.D0 +y1*y2**2*2d0*y4*SQRT3/6.D0 &
               +2.D0/3.D0*y1**2*y2*2d0*y4*SQRT3 &
               +y2**2*y3*2d0*y4*SQRT3/6.D0 )*DFEA13345 
            dds4(4,6)=dds3(4,6) &
               +( y1**2*y2*y5 +y1**2*y3*2d0*y4*SQRT3/3.D0 &
               +y1**2*y2*2d0*y4*SQRT3/3.D0 -y1*y2**2*2d0*y4*SQRT3/6.D0 &
               +y2*y3**2*y5 -y2**2*y3*y5 -y1**2*y3*y5 &
               +y2*y3**2*2d0*y4*SQRT3/3.D0 -y1*y3**2*2d0*y4*SQRT3/6.D0 &
               +y2**2*y3*2d0*y4*SQRT3/3.D0 )*DFEA11245 
            dds2(4,6)=dds4(4,6) &
               +( -y1*y2**3*SQRT3/2.D0 +y2**3*y3*SQRT3/2.D0 &
               +y2*y3**3*SQRT3/2.D0 -y1*y3**3*SQRT3/2.D0 )*DFEA11135 &
               +( y1**3*y3 -y2**3*y3*1d0/2.D0 &
               +y1**3*y2 -y2*y3**3*1d0/2.D0 -y1*y3**3*1d0/2.D0 &
               -y1*y2**3*1d0/2.D0 )*DFEA11134 
            ddv5(4,6)=dds2(4,6) &
               +0d0 +( -2.D0*y2**2*y3**2 +y1**2*y2**2 +y1**2*y3**2 )*DFEA11334 +0d0 &
               +( y1*y2*y3*2d0*y4 )*DFEA12355 &
               +( -y1*y2*y3**2*1d0/2.D0 -y1*y2**2*y3*1d0/2.D0 &
               +y1**2*y2*y3 )*DFEA11234 +0d0 &
               +0d0 
            dds3(4,6)=( y2**3*3d0*y4**2*SQRT3 -y2**3*2d0*y4*y5 &
               +y3**3*2d0*y4*y5 -5.D0/3.D0*y2**3*y5**2*SQRT3 &
               +y3**3*3d0*y4**2*SQRT3 -5.D0/3.D0*y3**3*y5**2*SQRT3 &
               -8.D0/3.D0*y1**3*y5**2*SQRT3 )*DFEA333555 &
               +( +y2**4*y5 +y2**4*2d0*y4*SQRT3/3.D0 &
               +y3**4*2d0*y4*SQRT3/3.D0 -y3**4*y5 &
               -y1**4*2d0*y4*SQRT3/6.D0 )*DFEA222245 &
               +0d0 +( y1**4*y3 -2.D0*y2**4*y3 +y1**4*y2 +y1*y3**4 -2.D0*y2*y3**4 &
               +y1*y2**4 )*DFEA133334 +( -y1*y2*y3*3d0*y4**2/3.D0 &
               +y1*y2*y3*y5**2 )*DFEA123455 
            dds4(4,6)=dds3(4,6) &
               +( 2.D0/3.D0*SQRT3*y1*y2**2*y3**2 -SQRT3*y1**2*y2**2*y3*1d0/3.D0 &
               -SQRT3*y1**2*y2*y3**2*1d0/3.D0 )*DFEA112335 &
               +( +y1*y2*y3**2*2d0*y4 +y1*y2**2*y3*2d0*y4 &
               +y1**2*y2*y3*2d0*y4 )*DFEA112355 
            dds2(4,6)=dds4(4,6) +( &
               +y1**3*y2**2*SQRT3/2.D0 -y1**2*y2**3*SQRT3/2.D0 &
               +y1**3*y3**2*SQRT3/2.D0 -y1**2*y3**3*SQRT3/2.D0 )*DFEA222335 &
               +( -y1**2*y2**2*2d0*y4*SQRT3/6.D0 -y1**2*y2**2*y5 &
               -2.D0/3.D0*y2**2*y3**2*2d0*y4*SQRT3 &
               +y1**2*y3**2*y5 -y1**2*y3**2*2d0*y4*SQRT3/6.D0 )*DFEA113345 +( &
               +y2**2*y3**2*2d0*y4 +y1**2*y3**2*2d0*y4 &
               +y1**2*y2**2*2d0*y4 )*DFEA223355 
            dds3(4,6)=dds2(4,6) &
               +( y1*y2*y3**2*2d0*y4*SQRT3/6.D0 +y1*y2*y3**2*y5 &
               +2.D0/3.D0*y1**2*y2*y3*2d0*y4*SQRT3 -y1*y2**2*y3*y5 &
               +y1*y2**2*y3*2d0*y4*SQRT3/6.D0 )*DFEA123345 &
               +( -y1**3*y2**2*1d0/2.D0 -y1**3*y3**2*1d0/2.D0 -y1**2*y2**3*1d0/2.D0 &
               -y1**2*y3**3*1d0/2.D0 &
               +y2**3*y3**2 +y2**2*y3**3 )*DFEA222334 +( 3.D0*y3**2*4d0*y4**3 &
               +3.D0*y2**2*4d0*y4**3 -4.D0*y3**2*y5**3*SQRT3 &
               +9.D0*y1**2*2d0*y4*y5**2 -3.D0/2.D0*y1**2*4d0*y4**3 &
               +4.D0*y2**2*y5**3*SQRT3 )*DFEA335555 +0d0 
            dds4(4,6)=dds3(4,6) &
               +( y3*5d0*y4**4/5.D0 -y2*4d0*y4**3*y5*SQRT3/2.D0 &
               -2.D0/5.D0*y1*5d0*y4**4 -2.D0*y1*3d0*y4**2*y5**2 &
               +y3*3d0*y4**2*y5**2 +y3*4d0*y4**3*y5*SQRT3/2.D0 +y2*3d0*y4**2*y5**2 &
               +y2*5d0*y4**4/5.D0 )*DFEA244455 +( y2**5 -2.D0*y1**5 &
               +y3**5 )*DFEA222224 
            dds5(4,6)=dds4(4,6) +( &
               +y1*y5**4 -7.D0/15.D0*y2*5d0*y4**4 &
               +y2*4d0*y4**3*y5*SQRT3/3.D0 -y3*4d0*y4**3*y5*SQRT3/3.D0 +y3*y5**4 &
               +y2*y5**4 &
               +2.D0*y1*3d0*y4**2*y5**2 -7.D0/15.D0*y3*5d0*y4**4 &
               -y1*5d0*y4**4/15.D0 )*DFEA145555 
            dds1(4,6)=dds5(4,6) &
               +( &
               +y1**3*y2*y3 -y1*y2**3*y3*1d0/2.D0 -y1*y2*y3**3*1d0/2.D0 )*DFEA111234 &
               +( y3*4d0*y4**3*y5/3.D0 &
               +y3*5d0*y4**4*SQRT3/18.D0 -y2*4d0*y4**3*y5/3.D0 -y2*y5**4*SQRT3/2.D0 &
               -y3*2d0*y4*y5**3 &
               +2.D0/9.D0*y1*5d0*y4**4*SQRT3 +y2*5d0*y4**4*SQRT3/18.D0 &
               +y2*2d0*y4*y5**3 -2.D0/3.D0*y1*3d0*y4**2*y5**2*SQRT3 &
               -y3*y5**4*SQRT3/2.D0 )*DFEA244555 &
               +( y1*y2*2d0*y4*y5**2 -3.D0/4.D0*y2*y3*4d0*y4**3 &
               +y1*y3*2d0*y4*y5**2 -7.D0/2.D0*y2*y3*2d0*y4*y5**2 &
               -2.D0*y1*y2*3d0*y4**2*y5*SQRT3 &
               +2.D0*y1*y3*3d0*y4**2*y5*SQRT3 )*DFEA124455 
            dds3(4,6)=dds1(4,6) +0d0 &
               +0d0 +( y1**4*2d0*y4 +y2**4*2d0*y4 &
               +y3**4*2d0*y4 )*DFEA222255 
            dds4(4,6)=dds3(4,6) +( +y1*y3*4d0*y4**3 &
               +9.D0*y2*y3*2d0*y4*y5**2 -4.D0*y1*y3*3d0*y4**2*y5*SQRT3 &
               +y1*y2*4d0*y4**3 +4.D0*y1*y2*3d0*y4**2*y5*SQRT3 &
               +5.D0/2.D0*y2*y3*4d0*y4**3 )*DFEA134444 +( -7.D0/3.D0*y1**2*y3*y5**2 &
               +5.D0/3.D0*y1*y2**2*2d0*y4*y5*SQRT3 -13.D0/3.D0*y2**2*y3*y5**2 &
               -7.D0/3.D0*y1**2*y2*y5**2 -16.D0/3.D0*y1*y3**2*y5**2 &
               +4.D0/3.D0*y1**2*y3*2d0*y4*y5*SQRT3 +3.D0*y1**2*y2*3d0*y4**2 &
               +y2*y3**2*3d0*y4**2 &
               +y2**2*y3*3d0*y4**2 -13.D0/3.D0*y2*y3**2*y5**2 &
               -5.D0/3.D0*y1*y3**2*2d0*y4*y5*SQRT3 &
               -4.D0/3.D0*y1**2*y2*2d0*y4*y5*SQRT3 &
               +3.D0*y1**2*y3*3d0*y4**2 &
               -16.D0/3.D0*y1*y2**2*y5**2 )*DFEA233444 
            dds5(4,6)=dds4(4,6) &
               +( +4.D0*y2**2*y3*y5**2*SQRT3 +y1**2*y3*y5**2*SQRT3 &
               +6.D0*y1*y3**2*2d0*y4*y5 -6.D0*y1*y2**2*2d0*y4*y5 &
               -3.D0*y1**2*y3*2d0*y4*y5 &
               +y1**2*y2*y5**2*SQRT3 &
               +4.D0*y1*y3**2*y5**2*SQRT3 -3.D0*y1**2*y2*3d0*y4**2*SQRT3 &
               +3.D0*y1**2*y2*2d0*y4*y5 -3.D0*y1**2*y3*3d0*y4**2*SQRT3 &
               +4.D0*y2*y3**2*y5**2*SQRT3 &
               +4.D0*y1*y2**2*y5**2*SQRT3 )*DFEA113555 
            dds2(4,6)=dds5(4,6) &
               +( -2.D0/3.D0*y3**2*4d0*y4**3*SQRT3 &
               -3.D0/2.D0*y1**2*2d0*y4*y5**2*SQRT3 -y2**2*3d0*y4**2*y5 &
               +7.D0/12.D0*y1**2*4d0*y4**3*SQRT3 +y3**2*3d0*y4**2*y5 &
               +3.D0*y3**2*y5**3 -2.D0/3.D0*y2**2*4d0*y4**3*SQRT3 &
               -3.D0*y2**2*y5**3 )*DFEA334445 &
               +( -3.D0*y1*y3*3d0*y4**2*y5 -y1*y3*y5**3 +3.D0*y1*y2*3d0*y4**2*y5 &
               +3.D0/2.D0*y2*y3*2d0*y4*y5**2*SQRT3 +y1*y2*y5**3 &
               +3.D0/4.D0*y2*y3*4d0*y4**3*SQRT3 )*DFEA124555 &
               +( 2.D0*y3**2*y5**3*SQRT3 -7.D0/2.D0*y1**2*2d0*y4*y5**2 &
               +y2**2*2d0*y4*y5**2 -y2**2*4d0*y4**3 -y3**2*4d0*y4**3 &
               -2.D0*y2**2*y5**3*SQRT3 &
               +5.D0/4.D0*y1**2*4d0*y4**3 &
               +y3**2*2d0*y4*y5**2 )*DFEA334455 
            dds3(4,6)=dds2(4,6) &
               +( -6.D0*2d0*y4*y5**4 +9.D0*4d0*y4**3*y5**2 )*DFEA555555 &
               +( y2*y3**3*2d0*y4 +y1*y3**3*2d0*y4 +y1*y2**3*2d0*y4 +y1**3*y2*2d0*y4 &
               +y1**3*y3*2d0*y4 +y2**3*y3*2d0*y4 )*DFEA233344 +( &
               +y1**3*y2*y5 -y1**3*y3*y5 -y1*y3**3*2d0*y4*SQRT3/2.D0 -y2**3*y3*y5 &
               +y2*y3**3*y5 -y1*y2**3*2d0*y4*SQRT3/2.D0 )*DFEA233345 &
               +( -3.D0*y2**3*y5**2 &
               +y3**3*3d0*y4**2 -3.D0*y3**3*y5**2 -3.D0*y1**3*y5**2 +y2**3*3d0*y4**2 &
               +y1**3*3d0*y4**2 )*DFEA111444 +0d0 
            dds4(4,6)=dds3(4,6) &
               +( 9.D0*2d0*y4*y5**4 -6.D0*4d0*y4**3*y5**2 +6d0*y4**5 )*DFEA444444 &
               +( -5.D0/3.D0*y1*y2**2*2d0*y4*y5*SQRT3 &
               +y1*y2**2*3d0*y4**2 -4.D0/3.D0*y1**2*y3*2d0*y4*y5*SQRT3 &
               -2.D0*y1**2*y2*3d0*y4**2 &
               +4.D0/3.D0*y2**2*y3*y5**2 -2.D0*y1**2*y3*3d0*y4**2 &
               +7.D0/3.D0*y1*y2**2*y5**2 -2.D0/3.D0*y1**2*y3*y5**2 &
               +y1*y3**2*3d0*y4**2 +4.D0/3.D0*y1**2*y2*2d0*y4*y5*SQRT3 &
               +4.D0/3.D0*y2*y3**2*y5**2 &
               +5.D0/3.D0*y1*y3**2*2d0*y4*y5*SQRT3 -2.D0/3.D0*y1**2*y2*y5**2 &
               +7.D0/3.D0*y1*y3**2*y5**2 )*DFEA133444 
            dds5(4,6)=dds4(4,6) &
               +( -y1**3*y2*y5 +y1*y3**3*2d0*y4*SQRT3/2.D0 &
               +y1**3*y3*2d0*y4*SQRT3/2.D0 +y1**3*y3*y5 &
               +y1**3*y2*2d0*y4*SQRT3/2.D0 -y1*y2**3*y5 +y1*y2**3*2d0*y4*SQRT3/2.D0 &
               +y1*y3**3*y5 )*DFEA133345 
            ddv6(4,6)=dds5(4,6) +( -y2**2*y3*2d0*y4*y5 &
               +y1**2*y3*y5**2*SQRT3/3.D0 +y2*y3**2*2d0*y4*y5 &
               +4.D0/3.D0*y2**2*y3*y5**2*SQRT3 &
               +4.D0/3.D0*y2*y3**2*y5**2*SQRT3 -y1*y2**2*2d0*y4*y5 &
               +4.D0/3.D0*y1*y3**2*y5**2*SQRT3 &
               +y1**2*y2*y5**2*SQRT3/3.D0 -y1**2*y2*3d0*y4**2*SQRT3 &
               +y1*y3**2*2d0*y4*y5 -y1**2*y3*3d0*y4**2*SQRT3 &
               +4.D0/3.D0*y1*y2**2*y5**2*SQRT3 )*DFEA233445 +( y2*y3**4*SQRT3 &
               +y2**4*y3*SQRT3 -y1**4*y3*SQRT3 -y1**4*y2*SQRT3 )*DFEA233335 &
               +0d0 
 
            ddv(4,6)= +ddv1(4,6) +ddv2(4,6) +ddv3(4,6) +ddv4(4,6) &
               +ddv5(4,6) +ddv6(4,6)
 
            ddv1(5,5)=0d0 
            ddv2(5,5)=0d0 +0d0 +0d0 +( 2d0 )*FEA44 
            ddv3(5,5)=0d0 +( &
               +y1*2d0 +y2*2d0*1d0/4.D0 +y3*2d0*1d0/4.D0 )*FEA155 +0d0 +( &
               +y4*2d0 )*FEA455 +( +3.D0/4.D0*y3*2d0 +3.D0/4.D0*y2*2d0 )*FEA144 +0d0 &
               +0d0 
            dds2(5,5)=( +4d0*3d0*y5**2 +2.D0*y4**2*2d0 )*FEA4444 &
               +( 3.D0/8.D0*SQRT3*y2*3d0*2d0*y5 -3.D0/8.D0*SQRT3*y3*3d0*2d0*y5 &
               -9.D0/8.D0*y2*y4*2d0 -9.D0/8.D0*y3*y4*2d0 )*FEA1444 &
               +( +y1**2*2d0 +y3**2*2d0*1d0/4.D0 &
               +y2**2*2d0*1d0/4.D0 )*FEA1155 
            dds1(5,5)=dds2(5,5) +( &
               +3.D0/4.D0*y3**2*2d0 +3.D0/4.D0*y2**2*2d0 )*FEA1144 +0d0 +0d0 +( &
               +3.D0/2.D0*y2*y3*2d0 )*FEA1244 
            dds2(5,5)=dds1(5,5) +( y1*y3*2d0 &
               +y1*y2*2d0 -y2*y3*2d0*1d0/2.D0 )*FEA1255 +0d0 +0d0 
            ddv4(5,5)=dds2(5,5) &
               +0d0 +0d0 +0d0 +( 5.D0/8.D0*y2*y4*2d0 +SQRT3*y2*3d0*2d0*y5/8.D0 &
               +y1*y4*2d0 -SQRT3*y3*3d0*2d0*y5/8.D0 &
               +5.D0/8.D0*y3*y4*2d0 )*FEA1455 
            dds3(5,5)=( -2.D0*y4**3*2d0 &
               -3.D0*y4*4d0*3d0*y5**2 )*FEA44444 &
               +( -4.D0*y3*y4*3d0*2d0*y5*SQRT3 +9.D0*y1*y4**2*2d0 &
               +4.D0*y2*y4*3d0*2d0*y5*SQRT3 +5.D0/2.D0*y1*4d0*3d0*y5**2 &
               +y2*4d0*3d0*y5**2 +y3*4d0*3d0*y5**2 )*FEA25555 +( &
               +y3*y4**2*2d0 -2.D0*y2*y4*3d0*2d0*y5*SQRT3 -7.D0/2.D0*y1*y4**2*2d0 &
               -3.D0/4.D0*y1*4d0*3d0*y5**2 &
               +2.D0*y3*y4*3d0*2d0*y5*SQRT3 &
               +y2*y4**2*2d0 )*FEA24455 
            dds2(5,5)=dds3(5,5) +( -3.D0*y3*y4*3d0*2d0*y5 &
               +3.D0/4.D0*y1*4d0*3d0*y5**2*SQRT3 +3.D0*y2*y4*3d0*2d0*y5 &
               +3.D0/2.D0*y1*y4**2*2d0*SQRT3 )*FEA24445 +( -y2**2*3d0*2d0*y5 &
               +y3**2*3d0*2d0*y5 -y1**2*y4*2d0*SQRT3 )*FEA33445 +( y3**2*y4*2d0 &
               +y1**2*y4*2d0 +y2**2*y4*2d0 )*FEA33455 
            dds1(5,5)=dds2(5,5) +( &
               +y2**3*2d0*SQRT3/3.D0 &
               +y3**3*2d0*SQRT3/3.D0 -y1**3*2d0*SQRT3/6.D0 )*FEA33345 +( +y3**3*2d0 &
               +y2**3*2d0 +y1**3*2d0 )*FEA33344 +0d0 +0d0 +( -y1*y2*3d0*2d0*y5 &
               +y2*y3*y4*2d0*SQRT3 +y1*y3*3d0*2d0*y5 )*FEA13445 +( y2*y3*y4*2d0 &
               +y1*y2*y4*2d0 +y1*y3*y4*2d0 )*FEA13455 
            dds3(5,5)=dds1(5,5) &
               +( y1**2*y3*2d0 +y2**2*y3*2d0 +y1*y2**2*2d0 +y1**2*y2*2d0 &
               +y1*y3**2*2d0 +y2*y3**2*2d0 )*FEA11255 +( +y1*y3**2*2d0*SQRT3/2.D0 &
               +y1*y2**2*2d0*SQRT3/2.D0 +y2**2*y3*2d0*SQRT3/2.D0 &
               +y2*y3**2*2d0*SQRT3/2.D0 )*FEA13345 
            dds4(5,5)=dds3(5,5) +( &
               +y1*y2**2*2d0*SQRT3/2.D0 &
               +y1*y3**2*2d0*SQRT3/2.D0 )*FEA11245 
            dds2(5,5)=dds4(5,5) +0d0 &
               +0d0 
            ddv5(5,5)=dds2(5,5) +0d0 +0d0 +0d0 +( +y1*y2*y3*2d0 )*FEA12355 &
               +0d0 +0d0 &
               +0d0 
            dds3(5,5)=( -5.D0/3.D0*y2**3*y4*2d0*SQRT3 &
               -5.D0/3.D0*y3**3*y4*2d0*SQRT3 -y2**3*3d0*2d0*y5 &
               +y3**3*3d0*2d0*y5 -8.D0/3.D0*y1**3*y4*2d0*SQRT3 )*FEA333555 &
               +( y1**4*2d0*SQRT3/2.D0 )*FEA222245 +0d0 +0d0 +( &
               +y1*y2*y3*y4*2d0 )*FEA123455 
            dds4(5,5)=dds3(5,5) +0d0 &
               +( y1*y2**2*y3*2d0 +y1*y2*y3**2*2d0 &
               +y1**2*y2*y3*2d0 )*FEA112355 
            dds2(5,5)=dds4(5,5) +0d0 &
               +( -y1**2*y2**2*2d0*SQRT3/2.D0 -y1**2*y3**2*2d0*SQRT3/2.D0 )*FEA113345 &
               +( y2**2*y3**2*2d0 +y1**2*y2**2*2d0 &
               +y1**2*y3**2*2d0 )*FEA223355 
            dds3(5,5)=dds2(5,5) +( &
               +y1*y2*y3**2*2d0*SQRT3/2.D0 +y1*y2**2*y3*2d0*SQRT3/2.D0 )*FEA123345 &
               +0d0 +( +5.D0/2.D0*y1**2*4d0*3d0*y5**2 &
               +y2**2*4d0*3d0*y5**2 -4.D0*y3**2*y4*3d0*2d0*y5*SQRT3 &
               +y3**2*4d0*3d0*y5**2 +9.D0*y1**2*y4**2*2d0 &
               +4.D0*y2**2*y4*3d0*2d0*y5*SQRT3 )*FEA335555 +0d0 
            dds4(5,5)=dds3(5,5) &
               +( -2.D0*y1*y4**3*2d0 -3.D0/10.D0*y2*5d0*4d0*y5**3*SQRT3 +y3*y4**3*2d0 &
               +y2*y4**3*2d0 +3.D0/10.D0*y3*5d0*4d0*y5**3*SQRT3 )*FEA244455 &
               +0d0 
            dds5(5,5)=dds4(5,5) +( -y3*5d0*4d0*y5**3*SQRT3/5.D0 &
               +y2*5d0*4d0*y5**3*SQRT3/5.D0 +y1*y4*4d0*3d0*y5**2 +y3*y4*4d0*3d0*y5**2 &
               +y2*y4*4d0*3d0*y5**2 &
               +2.D0*y1*y4**3*2d0 )*FEA145555 
            dds1(5,5)=dds5(5,5) +0d0 &
               +( -y2*y4*4d0*3d0*y5**2*SQRT3/2.D0 -y3*y4**2*3d0*2d0*y5 &
               +y2*y4**2*3d0*2d0*y5 -2.D0/3.D0*y1*y4**3*2d0*SQRT3 &
               -y3*y4*4d0*3d0*y5**2*SQRT3/2.D0 )*FEA244555 &
               +( y1*y2*y4**2*2d0 -y1*y2*4d0*3d0*y5**2 -y1*y3*4d0*3d0*y5**2 &
               +5.D0/4.D0*y2*y3*4d0*3d0*y5**2 &
               +y1*y3*y4**2*2d0 &
               -7.D0/2.D0*y2*y3*y4**2*2d0 )*FEA124455 
            dds3(5,5)=dds1(5,5) &
               +0d0 +0d0 +( +y2**4*2d0 +y1**4*2d0 &
               +y3**4*2d0 )*FEA222255 
            dds4(5,5)=dds3(5,5) +( 3.D0*y1*y3*4d0*3d0*y5**2 &
               +9.D0*y2*y3*y4**2*2d0 -3.D0/2.D0*y2*y3*4d0*3d0*y5**2 &
               +3.D0*y1*y2*4d0*3d0*y5**2 )*FEA134444 &
               +( -y1*y3**2*3d0*2d0*y5*SQRT3/3.D0 -7.D0/3.D0*y1**2*y3*y4*2d0 &
               -13.D0/3.D0*y2**2*y3*y4*2d0 -4.D0/3.D0*y2*y3**2*3d0*2d0*y5*SQRT3 &
               -7.D0/3.D0*y1**2*y2*y4*2d0 -16.D0/3.D0*y1*y3**2*y4*2d0 &
               +4.D0/3.D0*y2**2*y3*3d0*2d0*y5*SQRT3 &
               +y1*y2**2*3d0*2d0*y5*SQRT3/3.D0 -13.D0/3.D0*y2*y3**2*y4*2d0 &
               -16.D0/3.D0*y1*y2**2*y4*2d0 )*FEA233444 
            dds5(5,5)=dds4(5,5) &
               +( 2.D0*y1*y3**2*3d0*2d0*y5 +4.D0*y2*y3**2*3d0*2d0*y5 &
               +4.D0*y2**2*y3*y4*2d0*SQRT3 -2.D0*y1*y2**2*3d0*2d0*y5 &
               +y1**2*y3*y4*2d0*SQRT3 +y1**2*y2*y4*2d0*SQRT3 &
               +4.D0*y1*y3**2*y4*2d0*SQRT3 -4.D0*y2**2*y3*3d0*2d0*y5 &
               -y1**2*y2*3d0*2d0*y5 &
               +y1**2*y3*3d0*2d0*y5 +4.D0*y2*y3**2*y4*2d0*SQRT3 &
               +4.D0*y1*y2**2*y4*2d0*SQRT3 )*FEA113555 
            dds2(5,5)=dds5(5,5) &
               +( -3.D0/2.D0*y1**2*y4**2*2d0*SQRT3 &
               -3.D0/4.D0*y1**2*4d0*3d0*y5**2*SQRT3 &
               +3.D0*y3**2*y4*3d0*2d0*y5 -3.D0*y2**2*y4*3d0*2d0*y5 )*FEA334445 +( &
               +2.D0/3.D0*y1*y2*4d0*3d0*y5**2*SQRT3 -y1*y3*y4*3d0*2d0*y5 &
               +2.D0/3.D0*y1*y3*4d0*3d0*y5**2*SQRT3 &
               -7.D0/12.D0*y2*y3*4d0*3d0*y5**2*SQRT3 &
               +3.D0/2.D0*y2*y3*y4**2*2d0*SQRT3 +y1*y2*y4*3d0*2d0*y5 )*FEA124555 &
               +( 2.D0*y3**2*y4*3d0*2d0*y5*SQRT3 -7.D0/2.D0*y1**2*y4**2*2d0 &
               +y2**2*y4**2*2d0 -2.D0*y2**2*y4*3d0*2d0*y5*SQRT3 &
               -3.D0/4.D0*y1**2*4d0*3d0*y5**2 &
               +y3**2*y4**2*2d0 )*FEA334455 
            dds3(5,5)=dds2(5,5) &
               +( -6.D0*y4**2*4d0*3d0*y5**2 +9.D0*y4**4*2d0 &
               +6d0*5d0*y5**4 )*FEA555555 +( +y2*y3**3*2d0 +y1*y2**3*2d0 &
               +y1**3*y3*2d0 +y1**3*y2*2d0 +y1*y3**3*2d0 +y2**3*y3*2d0 )*FEA233344 &
               +( y1*y2**3*2d0*SQRT3/6.D0 -y2**3*y3*2d0*SQRT3/3.D0 &
               -y2*y3**3*2d0*SQRT3/3.D0 -y1**3*y2*2d0*SQRT3/3.D0 &
               -y1**3*y3*2d0*SQRT3/3.D0 &
               +y1*y3**3*2d0*SQRT3/6.D0 )*FEA233345 &
               +( -3.D0*y2**3*y4*2d0 -3.D0*y3**3*y4*2d0 -3.D0*y1**3*y4*2d0 )*FEA111444 &
               +0d0 
            dds4(5,5)=dds3(5,5) &
               +( 9.D0*y4**2*4d0*3d0*y5**2 -6.D0*y4**4*2d0 )*FEA444444 &
               +( -y1*y2**2*3d0*2d0*y5*SQRT3/3.D0 &
               +4.D0/3.D0*y2**2*y3*y4*2d0 -4.D0/3.D0*y2**2*y3*3d0*2d0*y5*SQRT3 &
               +7.D0/3.D0*y1*y2**2*y4*2d0 -2.D0/3.D0*y1**2*y3*y4*2d0 &
               +4.D0/3.D0*y2*y3**2*3d0*2d0*y5*SQRT3 +y1*y3**2*3d0*2d0*y5*SQRT3/3.D0 &
               +4.D0/3.D0*y2*y3**2*y4*2d0 -2.D0/3.D0*y1**2*y2*y4*2d0 &
               +7.D0/3.D0*y1*y3**2*y4*2d0 )*FEA133444 
            dds5(5,5)=dds4(5,5) +( &
               +2.D0/3.D0*y2**3*y3*2d0*SQRT3 +y1**3*y3*2d0*SQRT3/6.D0 &
               +y1**3*y2*2d0*SQRT3/6.D0 +y1*y2**3*2d0*SQRT3/6.D0 &
               +2.D0/3.D0*y2*y3**3*2d0*SQRT3 &
               +y1*y3**3*2d0*SQRT3/6.D0 )*FEA133345 
            ddv6(5,5)=dds5(5,5) +( &
               +y1**2*y3*y4*2d0*SQRT3/3.D0 +y2*y3**2*3d0*2d0*y5 -y1*y2**2*3d0*2d0*y5 &
               +4.D0/3.D0*y2**2*y3*y4*2d0*SQRT3 +4.D0/3.D0*y2*y3**2*y4*2d0*SQRT3 &
               +4.D0/3.D0*y1*y3**2*y4*2d0*SQRT3 -y2**2*y3*3d0*2d0*y5 &
               +y1*y3**2*3d0*2d0*y5 +y1**2*y2*y4*2d0*SQRT3/3.D0 &
               +4.D0/3.D0*y1*y2**2*y4*2d0*SQRT3 )*FEA233445 +0d0 +0d0 
 
            ddv(5,5)= &
               +ddv1(5,5) +ddv2(5,5) +ddv3(5,5) +ddv4(5,5) +ddv5(5,5) +ddv6(5,5)
 
            ddv1(5,6)=0d0 
            ddv2(5,6)=0d0 +0d0 +( -SQRT3*y3*1d0/2.D0 &
               +SQRT3*y2*1d0/2.D0 )*DFEA14 +( 2d0*y5 )*DFEA44 
            ddv3(5,6)=( &
               +SQRT3*y1*y2 -SQRT3*y1*y3 )*DFEA124 +( -SQRT3*y3*y4*1d0/2.D0 &
               +y1*2d0*y5 +y2*2d0*y5/4.D0 +SQRT3*y2*y4*1d0/2.D0 &
               +y3*2d0*y5/4.D0 )*DFEA155 +0d0 +( +y4*2d0*y5 )*DFEA455 +( &
               +3.D0/4.D0*y3*2d0*y5 +3.D0/4.D0*y2*2d0*y5 -SQRT3*y2*y4*1d0/2.D0 &
               +SQRT3*y3*y4*1d0/2.D0 )*DFEA144 +0d0 +( &
               +SQRT3*y2**2*1d0/2.D0 -SQRT3*y3**2*1d0/2.D0 )*DFEA114 
            dds2(5,6)=( &
               +4d0*y5**3 +2.D0*y4**2*2d0*y5 )*DFEA4444 &
               +( 3.D0/8.D0*SQRT3*y2*3d0*y5**2 -3.D0/8.D0*SQRT3*y3*y4**2 &
               -3.D0/8.D0*SQRT3*y3*3d0*y5**2 -9.D0/8.D0*y2*y4*2d0*y5 &
               -9.D0/8.D0*y3*y4*2d0*y5 &
               +3.D0/8.D0*SQRT3*y2*y4**2 )*DFEA1444 +( +y1**2*2d0*y5 &
               +y3**2*2d0*y5/4.D0 -SQRT3*y3**2*y4*1d0/2.D0 +SQRT3*y2**2*y4*1d0/2.D0 &
               +y2**2*2d0*y5/4.D0 )*DFEA1155 
            dds1(5,6)=dds2(5,6) +( &
               +3.D0/4.D0*y3**2*2d0*y5 &
               +SQRT3*y3**2*y4*1d0/2.D0 -SQRT3*y2**2*y4*1d0/2.D0 &
               +3.D0/4.D0*y2**2*2d0*y5 )*DFEA1144 +( &
               +SQRT3*y2**3*1d0/2.D0 -SQRT3*y3**3*1d0/2.D0 )*DFEA1114 +0d0 &
               +( SQRT3*y1*y3*y4 &
               +3.D0/2.D0*y2*y3*2d0*y5 -SQRT3*y1*y2*y4 )*DFEA1244 
            dds2(5,6)=dds1(5,6) &
               +( y1*y3*2d0*y5 +y1*y2*2d0*y5 -SQRT3*y1*y3*y4 -y2*y3*2d0*y5/2.D0 &
               +SQRT3*y1*y2*y4 )*DFEA1255 &
               +( -SQRT3*y1*y3**2*1d0/2.D0 -SQRT3*y2*y3**2*1d0/2.D0 &
               +SQRT3*y2**2*y3*1d0/2.D0 +SQRT3*y1*y2**2*1d0/2.D0 )*DFEA1124 &
               +( y1**2*y2 -y2**2*y3*1d0/2.D0 +y2*y3**2*1d0/2.D0 -y1*y3**2*1d0/2.D0 &
               +y1*y2**2*1d0/2.D0 -y1**2*y3 )*DFEA1125 
            ddv4(5,6)=dds2(5,6) +0d0 +0d0 &
               +0d0 +( 5.D0/8.D0*y2*y4*2d0*y5 &
               +SQRT3*y2*3d0*y5**2/8.D0 -SQRT3*y3*y4**2*1d0/8.D0 &
               +SQRT3*y2*y4**2*1d0/8.D0 +y1*y4*2d0*y5 -SQRT3*y3*3d0*y5**2/8.D0 &
               +5.D0/8.D0*y3*y4*2d0*y5 )*DFEA1455 
            dds3(5,6)=( -2.D0*y4**3*2d0*y5 &
               -3.D0*y4*4d0*y5**3 )*DFEA44444 &
               +( -4.D0*y3*y4*3d0*y5**2*SQRT3 +9.D0*y1*y4**2*2d0*y5 &
               +4.D0*y2*y4*3d0*y5**2*SQRT3 +5.D0/2.D0*y1*4d0*y5**3 +y2*4d0*y5**3 &
               +y3*4d0*y5**3 )*DFEA25555 +( &
               +y3*y4**2*2d0*y5 -2.D0*y2*y4*3d0*y5**2*SQRT3 &
               -7.D0/2.D0*y1*y4**2*2d0*y5 -3.D0/4.D0*y1*4d0*y5**3 &
               +2.D0*y3*y4*3d0*y5**2*SQRT3 &
               +y2*y4**2*2d0*y5 )*DFEA24455 
            dds2(5,6)=dds3(5,6) &
               +( y2*y4**3 -3.D0*y3*y4*3d0*y5**2 +3.D0/4.D0*y1*4d0*y5**3*SQRT3 &
               +3.D0*y2*y4*3d0*y5**2 &
               +3.D0/2.D0*y1*y4**2*2d0*y5*SQRT3 -y3*y4**3 )*DFEA24445 &
               +( -y2**2*3d0*y5**2 +y3**2*y4**2 &
               +y3**2*3d0*y5**2 -y2**2*y4**2 -y1**2*y4*2d0*y5*SQRT3 )*DFEA33445 &
               +( y3**2*y4*2d0*y5 +y1**2*y4*2d0*y5 &
               +y2**2*y4*2d0*y5 )*DFEA33455 
            dds1(5,6)=dds2(5,6) +( -y2**3*y4 &
               +y3**3*y4 +y2**3*2d0*y5*SQRT3/3.D0 &
               +y3**3*2d0*y5*SQRT3/3.D0 -y1**3*2d0*y5*SQRT3/6.D0 )*DFEA33345 +( &
               +y3**3*2d0*y5 +y2**3*2d0*y5 +y1**3*2d0*y5 )*DFEA33344 +( &
               +SQRT3*y3**4 -SQRT3*y2**4 )*DFEA33334 +0d0 +( -y1*y2*3d0*y5**2 &
               +y1*y3*y4**2 +y2*y3*y4*2d0*y5*SQRT3 -y1*y2*y4**2 &
               +y1*y3*3d0*y5**2 )*DFEA13445 +( y2*y3*y4*2d0*y5 +y1*y2*y4*2d0*y5 &
               +y1*y3*y4*2d0*y5 )*DFEA13455 
            dds3(5,6)=dds1(5,6) +( y1**2*y3*2d0*y5 &
               +y2**2*y3*2d0*y5 +y1*y2**2*2d0*y5 +y1**2*y2*2d0*y5 +y1*y3**2*2d0*y5 &
               +y2*y3**2*2d0*y5 )*DFEA11255 +( +y1*y3**2*2d0*y5*SQRT3/2.D0 &
               +y1*y2**2*2d0*y5*SQRT3/2.D0 +y2**2*y3*2d0*y5*SQRT3/2.D0 -y1*y2**2*y4 &
               +y2*y3**2*y4 +y1*y3**2*y4 -y2**2*y3*y4 &
               +y2*y3**2*2d0*y5*SQRT3/2.D0 )*DFEA13345 
            dds4(5,6)=dds3(5,6) &
               +( y1**2*y2*y4 +y2*y3**2*y4 -y2**2*y3*y4 -y1**2*y3*y4 &
               +y1*y2**2*2d0*y5*SQRT3/2.D0 &
               +y1*y3**2*2d0*y5*SQRT3/2.D0 )*DFEA11245 
            dds2(5,6)=dds4(5,6) &
               +( -y1**3*y2 +y1**3*y3 &
               +y2**3*y3*1d0/2.D0 -y1*y2**3*1d0/2.D0 -y2*y3**3*1d0/2.D0 &
               +y1*y3**3*1d0/2.D0 )*DFEA11135 +( +y1*y2**3*SQRT3/2.D0 &
               +y2**3*y3*SQRT3/2.D0 -y2*y3**3*SQRT3/2.D0 &
               -y1*y3**3*SQRT3/2.D0 )*DFEA11134 
            ddv5(5,6)=dds2(5,6) &
               +0d0 +( -SQRT3*y1**2*y3**2 +SQRT3*y1**2*y2**2 )*DFEA11334 +0d0 +( &
               +y1*y2*y3*2d0*y5 )*DFEA12355 +( -SQRT3*y1*y2*y3**2*1d0/2.D0 &
               +SQRT3*y1*y2**2*y3*1d0/2.D0 )*DFEA11234 +0d0 &
               +0d0 
            dds3(5,6)=( -y2**3*y4**2 &
               +y3**3*y4**2 -5.D0/3.D0*y2**3*y4*2d0*y5*SQRT3 &
               -5.D0/3.D0*y3**3*y4*2d0*y5*SQRT3 -y2**3*3d0*y5**2 &
               +y3**3*3d0*y5**2 -8.D0/3.D0*y1**3*y4*2d0*y5*SQRT3 )*DFEA333555 &
               +( y1**4*2d0*y5*SQRT3/2.D0 +y2**4*y4 -y3**4*y4 )*DFEA222245 +0d0 +( &
               +y1*y2**4*SQRT3 &
               +y1**4*y2*SQRT3 -y1*y3**4*SQRT3 -y1**4*y3*SQRT3 )*DFEA133334 +( &
               +y1*y2*y3*y4*2d0*y5 )*DFEA123455 
            dds4(5,6)=dds3(5,6) &
               +( -y1**2*y2**2*y3 +y1**2*y2*y3**2 )*DFEA112335 +( y1*y2**2*y3*2d0*y5 &
               +y1*y2*y3**2*2d0*y5 &
               +y1**2*y2*y3*2d0*y5 )*DFEA112355 
            dds2(5,6)=dds4(5,6) &
               +( y2**3*y3**2 -y1**3*y2**2*1d0/2.D0 -y1**2*y3**3*1d0/2.D0 -y2**2*y3**3 &
               +y1**3*y3**2*1d0/2.D0 +y1**2*y2**3*1d0/2.D0 )*DFEA222335 &
               +( -y1**2*y2**2*2d0*y5*SQRT3/2.D0 -y1**2*y3**2*2d0*y5*SQRT3/2.D0 &
               -y1**2*y2**2*y4 &
               +y1**2*y3**2*y4 )*DFEA113345 +( y2**2*y3**2*2d0*y5 +y1**2*y2**2*2d0*y5 &
               +y1**2*y3**2*2d0*y5 )*DFEA223355 
            dds3(5,6)=dds2(5,6) +( &
               +y1*y2*y3**2*y4 +y1*y2*y3**2*2d0*y5*SQRT3/2.D0 -y1*y2**2*y3*y4 &
               +y1*y2**2*y3*2d0*y5*SQRT3/2.D0 )*DFEA123345 +( -y1**3*y2**2*SQRT3/2.D0 &
               +y1**3*y3**2*SQRT3/2.D0 -y1**2*y2**3*SQRT3/2.D0 &
               +y1**2*y3**3*SQRT3/2.D0 )*DFEA222334 +( +5.D0/2.D0*y1**2*4d0*y5**3 &
               +y2**2*4d0*y5**3 -4.D0*y3**2*y4*3d0*y5**2*SQRT3 +y3**2*4d0*y5**3 &
               +9.D0*y1**2*y4**2*2d0*y5 +4.D0*y2**2*y4*3d0*y5**2*SQRT3 )*DFEA335555 &
               +0d0 
            dds4(5,6)=dds3(5,6) &
               +( -y2*y4**4*SQRT3/2.D0 -2.D0*y1*y4**3*2d0*y5 &
               -3.D0/10.D0*y2*5d0*y5**4*SQRT3 &
               +y3*y4**3*2d0*y5 +y3*y4**4*SQRT3/2.D0 +y2*y4**3*2d0*y5 &
               +3.D0/10.D0*y3*5d0*y5**4*SQRT3 )*DFEA244455 +( -SQRT3*y2**5 &
               +SQRT3*y3**5 )*DFEA222224 
            dds5(5,6)=dds4(5,6) &
               +( -y3*5d0*y5**4*SQRT3/5.D0 +y2*5d0*y5**4*SQRT3/5.D0 +y1*y4*4d0*y5**3 &
               +y2*y4**4*SQRT3/3.D0 -y3*y4**4*SQRT3/3.D0 +y3*y4*4d0*y5**3 &
               +y2*y4*4d0*y5**3 &
               +2.D0*y1*y4**3*2d0*y5 )*DFEA145555 
            dds1(5,6)=dds5(5,6) &
               +( -SQRT3*y1*y2*y3**3*1d0/2.D0 &
               +SQRT3*y1*y2**3*y3*1d0/2.D0 )*DFEA111234 &
               +( y3*y4**4*1d0/3.D0 -y2*y4**4*1d0/3.D0 -y2*y4*4d0*y5**3*SQRT3/2.D0 &
               -y3*y4**2*3d0*y5**2 &
               +y2*y4**2*3d0*y5**2 -2.D0/3.D0*y1*y4**3*2d0*y5*SQRT3 &
               -y3*y4*4d0*y5**3*SQRT3/2.D0 )*DFEA244555 &
               +( y1*y2*y4**2*2d0*y5 -y1*y2*4d0*y5**3 -y1*y3*4d0*y5**3 &
               +5.D0/4.D0*y2*y3*4d0*y5**3 &
               +y1*y3*y4**2*2d0*y5 -7.D0/2.D0*y2*y3*y4**2*2d0*y5 &
               -2.D0*y1*y2*y4**3*SQRT3 &
               +2.D0*y1*y3*y4**3*SQRT3 )*DFEA124455 
            dds3(5,6)=dds1(5,6) +0d0 +0d0 +( &
               +y2**4*2d0*y5 +y1**4*2d0*y5 &
               +y3**4*2d0*y5 )*DFEA222255 
            dds4(5,6)=dds3(5,6) +( 3.D0*y1*y3*4d0*y5**3 &
               +9.D0*y2*y3*y4**2*2d0*y5 -3.D0/2.D0*y2*y3*4d0*y5**3 &
               -4.D0*y1*y3*y4**3*SQRT3 &
               +4.D0*y1*y2*y4**3*SQRT3 +3.D0*y1*y2*4d0*y5**3 )*DFEA134444 &
               +( -y1*y3**2*3d0*y5**2*SQRT3/3.D0 -7.D0/3.D0*y1**2*y3*y4*2d0*y5 &
               +5.D0/3.D0*y1*y2**2*y4**2*SQRT3 -13.D0/3.D0*y2**2*y3*y4*2d0*y5 &
               -4.D0/3.D0*y2*y3**2*3d0*y5**2*SQRT3 -7.D0/3.D0*y1**2*y2*y4*2d0*y5 &
               -16.D0/3.D0*y1*y3**2*y4*2d0*y5 &
               +4.D0/3.D0*y1**2*y3*y4**2*SQRT3 +4.D0/3.D0*y2**2*y3*3d0*y5**2*SQRT3 &
               +y1*y2**2*3d0*y5**2*SQRT3/3.D0 -13.D0/3.D0*y2*y3**2*y4*2d0*y5 &
               -5.D0/3.D0*y1*y3**2*y4**2*SQRT3 -4.D0/3.D0*y1**2*y2*y4**2*SQRT3 &
               -16.D0/3.D0*y1*y2**2*y4*2d0*y5 )*DFEA233444 
            dds5(5,6)=dds4(5,6) &
               +( 2.D0*y1*y3**2*3d0*y5**2 +4.D0*y2*y3**2*3d0*y5**2 &
               +4.D0*y2**2*y3*y4*2d0*y5*SQRT3 -2.D0*y1*y2**2*3d0*y5**2 &
               +y1**2*y3*y4*2d0*y5*SQRT3 &
               +6.D0*y1*y3**2*y4**2 -6.D0*y1*y2**2*y4**2 -3.D0*y1**2*y3*y4**2 &
               +y1**2*y2*y4*2d0*y5*SQRT3 &
               +4.D0*y1*y3**2*y4*2d0*y5*SQRT3 -4.D0*y2**2*y3*3d0*y5**2 &
               +3.D0*y1**2*y2*y4**2 -y1**2*y2*3d0*y5**2 +y1**2*y3*3d0*y5**2 &
               +4.D0*y2*y3**2*y4*2d0*y5*SQRT3 &
               +4.D0*y1*y2**2*y4*2d0*y5*SQRT3 )*DFEA113555 
            dds2(5,6)=dds5(5,6) &
               +( -3.D0/2.D0*y1**2*y4**2*2d0*y5*SQRT3 &
               -3.D0/4.D0*y1**2*4d0*y5**3*SQRT3 -y2**2*y4**3 &
               +y3**2*y4**3 &
               +3.D0*y3**2*y4*3d0*y5**2 -3.D0*y2**2*y4*3d0*y5**2 )*DFEA334445 &
               +( -3.D0*y1*y3*y4**3 &
               +2.D0/3.D0*y1*y2*4d0*y5**3*SQRT3 -y1*y3*y4*3d0*y5**2 &
               +2.D0/3.D0*y1*y3*4d0*y5**3*SQRT3 &
               +3.D0*y1*y2*y4**3 -7.D0/12.D0*y2*y3*4d0*y5**3*SQRT3 &
               +3.D0/2.D0*y2*y3*y4**2*2d0*y5*SQRT3 +y1*y2*y4*3d0*y5**2 )*DFEA124555 &
               +( 2.D0*y3**2*y4*3d0*y5**2*SQRT3 -7.D0/2.D0*y1**2*y4**2*2d0*y5 &
               +y2**2*y4**2*2d0*y5 -2.D0*y2**2*y4*3d0*y5**2*SQRT3 &
               -3.D0/4.D0*y1**2*4d0*y5**3 &
               +y3**2*y4**2*2d0*y5 )*DFEA334455 
            dds3(5,6)=dds2(5,6) &
               +( -6.D0*y4**2*4d0*y5**3 +9.D0*y4**4*2d0*y5 +6d0*y5**5 )*DFEA555555 +( &
               +y2*y3**3*2d0*y5 +y1*y2**3*2d0*y5 +y1**3*y3*2d0*y5 +y1**3*y2*2d0*y5 &
               +y1*y3**3*2d0*y5 +y2**3*y3*2d0*y5 )*DFEA233344 &
               +( y1*y2**3*2d0*y5*SQRT3/6.D0 -y2**3*y3*2d0*y5*SQRT3/3.D0 &
               -y2*y3**3*2d0*y5*SQRT3/3.D0 &
               +y1**3*y2*y4 -y1**3*y2*2d0*y5*SQRT3/3.D0 -y1**3*y3*y4 &
               -y1**3*y3*2d0*y5*SQRT3/3.D0 &
               +y1*y3**3*2d0*y5*SQRT3/6.D0 -y2**3*y3*y4 +y2*y3**3*y4 )*DFEA233345 &
               +( -3.D0*y2**3*y4*2d0*y5 -3.D0*y3**3*y4*2d0*y5 &
               -3.D0*y1**3*y4*2d0*y5 )*DFEA111444 &
               +0d0 
            dds4(5,6)=dds3(5,6) &
               +( 9.D0*y4**2*4d0*y5**3 -6.D0*y4**4*2d0*y5 )*DFEA444444 &
               +( -5.D0/3.D0*y1*y2**2*y4**2*SQRT3 -4.D0/3.D0*y1**2*y3*y4**2*SQRT3 &
               -y1*y2**2*3d0*y5**2*SQRT3/3.D0 &
               +4.D0/3.D0*y2**2*y3*y4*2d0*y5 -4.D0/3.D0*y2**2*y3*3d0*y5**2*SQRT3 &
               +7.D0/3.D0*y1*y2**2*y4*2d0*y5 -2.D0/3.D0*y1**2*y3*y4*2d0*y5 &
               +4.D0/3.D0*y2*y3**2*3d0*y5**2*SQRT3 +y1*y3**2*3d0*y5**2*SQRT3/3.D0 &
               +4.D0/3.D0*y1**2*y2*y4**2*SQRT3 +4.D0/3.D0*y2*y3**2*y4*2d0*y5 &
               +5.D0/3.D0*y1*y3**2*y4**2*SQRT3 -2.D0/3.D0*y1**2*y2*y4*2d0*y5 &
               +7.D0/3.D0*y1*y3**2*y4*2d0*y5 )*DFEA133444 
            dds5(5,6)=dds4(5,6) &
               +( -y1**3*y2*y4 +2.D0/3.D0*y2**3*y3*2d0*y5*SQRT3 &
               +y1**3*y3*2d0*y5*SQRT3/6.D0 +y1**3*y2*2d0*y5*SQRT3/6.D0 +y1**3*y3*y4 &
               +y1*y2**3*2d0*y5*SQRT3/6.D0 &
               +2.D0/3.D0*y2*y3**3*2d0*y5*SQRT3 -y1*y2**3*y4 &
               +y1*y3**3*2d0*y5*SQRT3/6.D0 &
               +y1*y3**3*y4 )*DFEA133345 
            ddv6(5,6)=dds5(5,6) +( -y2**2*y3*y4**2 &
               +y1**2*y3*y4*2d0*y5*SQRT3/3.D0 +y2*y3**2*y4**2 &
               +y2*y3**2*3d0*y5**2 -y1*y2**2*3d0*y5**2 &
               +4.D0/3.D0*y2**2*y3*y4*2d0*y5*SQRT3 &
               +4.D0/3.D0*y2*y3**2*y4*2d0*y5*SQRT3 -y1*y2**2*y4**2 &
               +4.D0/3.D0*y1*y3**2*y4*2d0*y5*SQRT3 -y2**2*y3*3d0*y5**2 &
               +y1*y3**2*3d0*y5**2 +y1**2*y2*y4*2d0*y5*SQRT3/3.D0 +y1*y3**2*y4**2 &
               +4.D0/3.D0*y1*y2**2*y4*2d0*y5*SQRT3 )*DFEA233445 +( -y1**4*y2 &
               +y2*y3**4 -2.D0*y1*y2**4 +2.D0*y1*y3**4 &
               +y1**4*y3 -y2**4*y3 )*DFEA233335 +0d0 
 
            ddv(5,6)= +ddv1(5,6) &
               +ddv2(5,6) +ddv3(5,6) +ddv4(5,6) +ddv5(5,6) +ddv6(5,6)
 
            ddv1(6,6)=( y3 +y2 +y1 )*DDFEA1  
            ddv2(6,6)=( y2*y3 +y1*y3 &
               +y1*y2 )*DDFEA12 +( y2**2 +y3**2 +y1**2 )*DDFEA11 &
               +(  -SQRT3*y3*y5/2.D0 -y3*y4/2.D0 +y1*y4 &
               +SQRT3*y2*y5/2.D0 -y2*y4/2.D0 )*DDFEA14 +( y5**2 &
               +y4**2 )*DDFEA44  
            ddv3(6,6)=( y1*y3*y4 +y1*y2*y4 -2.D0*y2*y3*y4 &
               +SQRT3*y1*y2*y5 -SQRT3*y1*y3*y5 )*DDFEA124 &
               +( 3.D0/4.D0*y3*y4**2 -SQRT3*y3*y4*y5/2.D0 +y1*y5**2 +y2*y5**2/4.D0 &
               +3.D0/4.D0*y2*y4**2 +SQRT3*y2*y4*y5/2.D0 +y3*y5**2/4.D0 )*DDFEA155 &
               +( y2*y3**2 +y1*y3**2 +y1**2*y3 +y1*y2**2 +y2**2*y3 &
               +y1**2*y2 )*DDFEA112 +(  -y4**3/3.D0 +y4*y5**2 )*DDFEA455 &
               +y1*y2*y3*DDFEA123 +( y1*y4**2 +3.D0/4.D0*y3*y5**2 +3.D0/4.D0*y2*y5**2 &
               +y2*y4**2/4.D0 -SQRT3*y2*y4*y5/2.D0 +SQRT3*y3*y4*y5/2.D0 &
               +y3*y4**2/4.D0 )*DDFEA144 +( y3**3 +y2**3 +y1**3 )*DDFEA111 &
               +(  -y2**2*y4/2.D0 -y3**2*y4/2.D0 +SQRT3*y2**2*y5/2.D0 &
               +y1**2*y4 -SQRT3*y3**2*y5/2.D0 )*DDFEA114   
            dds2(6,6)=( y4**4 +y5**4 &
               +2.D0*y4**2*y5**2 )*DDFEA4444 &
               +( 3.D0/8.D0*SQRT3*y2*y5**3 -3.D0/8.D0*SQRT3*y3*y4**2*y5 &
               -3.D0/8.D0*SQRT3*y3*y5**3 -9.D0/8.D0*y2*y4*y5**2 -y3*y4**3/8.D0 &
               -y2*y4**3/8.D0 -9.D0/8.D0*y3*y4*y5**2 &
               +y1*y4**3 +3.D0/8.D0*SQRT3*y2*y4**2*y5 )*DDFEA1444 &
               +( 3.D0/4.D0*y2**2*y4**2 +3.D0/4.D0*y3**2*y4**2 +y1**2*y5**2 &
               +y3**2*y5**2/4.D0 -SQRT3*y3**2*y4*y5/2.D0 +SQRT3*y2**2*y4*y5/2.D0 &
               +y2**2*y5**2/4.D0 )*DDFEA1155 
            dds1(6,6)=dds2(6,6) +( y3**2*y4**2/4.D0 &
               +3.D0/4.D0*y3**2*y5**2 +y1**2*y4**2 +y2**2*y4**2/4.D0 &
               +SQRT3*y3**2*y4*y5/2.D0 -SQRT3*y2**2*y4*y5/2.D0 &
               +3.D0/4.D0*y2**2*y5**2 )*DDFEA1144 +( y1**3*y4 &
               +SQRT3*y2**3*y5/2.D0 -SQRT3*y3**3*y5/2.D0 -y2**3*y4/2.D0 &
               -y3**3*y4/2.D0 )*DDFEA1114 &
               +( y2**4 +y1**4 +y3**4 )*DDFEA1111 +( SQRT3*y1*y3*y4*y5 &
               +3.D0/2.D0*y2*y3*y5**2 -y2*y3*y4**2/2.D0 &
               +y1*y2*y4**2 -SQRT3*y1*y2*y4*y5 &
               +y1*y3*y4**2 )*DDFEA1244 
            dds2(6,6)=dds1(6,6) +( y1*y3*y5**2 &
               +y1*y2*y5**2 -SQRT3*y1*y3*y4*y5 -y2*y3*y5**2/2.D0 &
               +3.D0/2.D0*y2*y3*y4**2 +SQRT3*y1*y2*y4*y5 )*DDFEA1255 &
               +(  -y1*y3**2*y4/2.D0 &
               +y1**2*y3*y4 -SQRT3*y1*y3**2*y5/2.D0 -SQRT3*y2*y3**2*y5/2.D0 &
               +y1**2*y2*y4 +SQRT3*y2**2*y3*y5/2.D0 -y2**2*y3*y4/2.D0 &
               +SQRT3*y1*y2**2*y5/2.D0 -y2*y3**2*y4/2.D0 -y1*y2**2*y4/2.D0 )*DDFEA1124 &
               +( y1**2*y2*y5 +SQRT3*y1*y3**2*y4/2.D0 &
               +SQRT3*y1*y2**2*y4/2.D0 -SQRT3*y2*y3**2*y4/2.D0 &
               -SQRT3*y2**2*y3*y4/2.D0 -y2**2*y3*y5/2.D0 &
               +y2*y3**2*y5/2.D0 -y1*y3**2*y5/2.D0 &
               +y1*y2**2*y5/2.D0 -y1**2*y3*y5 )*DDFEA1125 
            ddv4(6,6)=dds2(6,6) &
               +( y2*y3**3 +y1**3*y3 +y1**3*y2 +y1*y2**3 +y1*y3**3 &
               +y2**3*y3 )*DDFEA1112 +( y2**2*y3**2 +y1**2*y3**2 &
               +y1**2*y2**2 )*DDFEA1122 +( y1*y2**2*y3 +y1**2*y2*y3 &
               +y1*y2*y3**2 )*DDFEA1123 +( 5.D0/8.D0*y2*y4*y5**2 &
               +SQRT3*y2*y5**3/8.D0 -SQRT3*y3*y4**2*y5/8.D0 &
               +SQRT3*y2*y4**2*y5/8.D0 -3.D0/8.D0*y2*y4**3 &
               +y1*y4*y5**2 -SQRT3*y3*y5**3/8.D0 &
               +5.D0/8.D0*y3*y4*y5**2 &
               -3.D0/8.D0*y3*y4**3 )*DDFEA1455   
            dds3(6,6)=( y4**5 -2.D0*y4**3*y5**2 &
               -3.D0*y4*y5**4 )*DDFEA44444 &
               +(  -4.D0*y3*y4*y5**3*SQRT3 +9.D0*y1*y4**2*y5**2 -3.D0/2.D0*y1*y4**4 &
               +4.D0*y2*y4*y5**3*SQRT3 +3.D0*y2*y4**4 +5.D0/2.D0*y1*y5**4 &
               +3.D0*y3*y4**4 +y2*y5**4 +y3*y5**4 )*DDFEA25555 +(  -y2*y4**4 &
               +y3*y4**2*y5**2 -2.D0*y2*y4*y5**3*SQRT3 -y3*y4**4 &
               -7.D0/2.D0*y1*y4**2*y5**2 -3.D0/4.D0*y1*y5**4 &
               +2.D0*y3*y4*y5**3*SQRT3 +y2*y4**2*y5**2 &
               +5.D0/4.D0*y1*y4**4 )*DDFEA24455 
            dds2(6,6)=dds3(6,6) &
               +( y2*y4**3*y5 -3.D0*y3*y4*y5**3 +2.D0/3.D0*y3*y4**4*SQRT3 &
               +3.D0/4.D0*y1*y5**4*SQRT3 +3.D0*y2*y4*y5**3 -7.D0/12.D0*y1*y4**4*SQRT3 &
               +3.D0/2.D0*y1*y4**2*y5**2*SQRT3 -y3*y4**3*y5 &
               +2.D0/3.D0*y2*y4**4*SQRT3 )*DDFEA24445 +(  -y2**2*y5**3 &
               +y3**2*y4**2*y5 +y3**2*y5**3 &
               +4.D0/9.D0*y2**2*y4**3*SQRT3 -5.D0/9.D0*y1**2*y4**3*SQRT3 &
               +4.D0/9.D0*y3**2*y4**3*SQRT3 -y2**2*y4**2*y5 &
               -y1**2*y4*y5**2*SQRT3 )*DDFEA33445 &
               +( y3**2*y4*y5**2 -y1**2*y4**3/3.D0 -y3**2*y4**3/3.D0 +y1**2*y4*y5**2 &
               +y2**2*y4*y5**2 -y2**2*y4**3/3.D0 )*DDFEA33455 
            dds1(6,6)=dds2(6,6) &
               +(  -y2**3*y4*y5 +y3**3*y4*y5 +y2**3*y5**2*SQRT3/3.D0 &
               +y1**3*y4**2*SQRT3/2.D0 &
               +y3**3*y5**2*SQRT3/3.D0 -y1**3*y5**2*SQRT3/6.D0 )*DDFEA33345 &
               +( y3**3*y4**2 +y3**3*y5**2 +y2**3*y4**2 +y2**3*y5**2 +y1**3*y5**2 &
               +y1**3*y4**2 )*DDFEA33344 +( y3**4*y4 +SQRT3*y3**4*y5 &
               +y2**4*y4 -2.D0*y1**4*y4 -SQRT3*y2**4*y5 )*DDFEA33334 +( y2**5 +y3**5 &
               +y1**5 )*DDFEA33333 +(  -4.D0/9.D0*y1*y2*y4**3*SQRT3 -y1*y2*y5**3 &
               +y1*y3*y4**2*y5 +y2*y3*y4*y5**2*SQRT3 -y1*y2*y4**2*y5 &
               +5.D0/9.D0*y2*y3*y4**3*SQRT3 -4.D0/9.D0*y1*y3*y4**3*SQRT3 &
               +y1*y3*y5**3 )*DDFEA13445 +( y2*y3*y4*y5**2 &
               +y1*y2*y4*y5**2 -y2*y3*y4**3/3.D0 -y1*y2*y4**3/3.D0 -y1*y3*y4**3/3.D0 &
               +y1*y3*y4*y5**2 )*DDFEA13455 
            dds3(6,6)=dds1(6,6) +( y1**2*y3*y5**2 &
               +y2**2*y3*y4**2 +y2**2*y3*y5**2 +y1*y2**2*y5**2 +y1**2*y2*y5**2 &
               +y1*y2**2*y4**2 +y2*y3**2*y4**2 +y1*y3**2*y4**2 +y1**2*y3*y4**2 &
               +y1**2*y2*y4**2 +y1*y3**2*y5**2 +y2*y3**2*y5**2 )*DDFEA11255 &
               +( 2.D0/3.D0*y1**2*y3*y4**2*SQRT3 +y1*y3**2*y5**2*SQRT3/2.D0 &
               +y1*y2**2*y5**2*SQRT3/2.D0 +y2**2*y3*y5**2*SQRT3/2.D0 -y1*y2**2*y4*y5 &
               +y2*y3**2*y4*y5 +y1*y3**2*y4*y5 -y2**2*y3*y4*y5 &
               +y2*y3**2*y4**2*SQRT3/6.D0 +y1*y3**2*y4**2*SQRT3/6.D0 &
               +y1*y2**2*y4**2*SQRT3/6.D0 +2.D0/3.D0*y1**2*y2*y4**2*SQRT3 &
               +y2*y3**2*y5**2*SQRT3/2.D0 &
               +y2**2*y3*y4**2*SQRT3/6.D0 )*DDFEA13345 
            dds4(6,6)=dds3(6,6) &
               +( y1**2*y2*y4*y5 +y1**2*y3*y4**2*SQRT3/3.D0 &
               +y1**2*y2*y4**2*SQRT3/3.D0 -y1*y2**2*y4**2*SQRT3/6.D0 &
               +y2*y3**2*y4*y5 -y2**2*y3*y4*y5 -y1**2*y3*y4*y5 &
               +y2*y3**2*y4**2*SQRT3/3.D0 &
               +y1*y2**2*y5**2*SQRT3/2.D0 -y1*y3**2*y4**2*SQRT3/6.D0 &
               +y2**2*y3*y4**2*SQRT3/3.D0 &
               +y1*y3**2*y5**2*SQRT3/2.D0 )*DDFEA11245 
            dds2(6,6)=dds4(6,6) &
               +(  -y1**3*y2*y5 +y1**3*y3*y5 &
               +y2**3*y3*y5/2.D0 -y1*y2**3*y4*SQRT3/2.D0 -y1*y2**3*y5/2.D0 &
               -y2*y3**3*y5/2.D0 &
               +y1*y3**3*y5/2.D0 +y2**3*y3*y4*SQRT3/2.D0 &
               +y2*y3**3*y4*SQRT3/2.D0 -y1*y3**3*y4*SQRT3/2.D0 )*DDFEA11135 &
               +( y1**3*y3*y4 -y2**3*y3*y4/2.D0 &
               +y1**3*y2*y4 -y2*y3**3*y4/2.D0 -y1*y3**3*y4/2.D0 &
               +y1*y2**3*y5*SQRT3/2.D0 &
               +y2**3*y3*y5*SQRT3/2.D0 -y2*y3**3*y5*SQRT3/2.D0 -y1*y2**3*y4/2.D0 &
               -y1*y3**3*y5*SQRT3/2.D0 )*DDFEA11134 
            ddv5(6,6)=dds2(6,6) &
               +( y1*y2**4 +y1**4*y3 +y1**4*y2 +y2**4*y3 +y2*y3**4 &
               +y1*y3**4 )*DDFEA23333 +(  -2.D0*y2**2*y3**2*y4 &
               +y1**2*y2**2*y4 -SQRT3*y1**2*y3**2*y5 +SQRT3*y1**2*y2**2*y5 &
               +y1**2*y3**2*y4 )*DDFEA11334 +( y1**2*y3**3 +y1**3*y3**2 +y2**2*y3**3 &
               +y1**2*y2**3 +y1**3*y2**2 +y2**3*y3**2 )*DDFEA11333 +( y1*y2*y3*y4**2 &
               +y1*y2*y3*y5**2 )*DDFEA12355 &
               +(  -y1*y2*y3**2*y4/2.D0 -y1*y2**2*y3*y4/2.D0 &
               -SQRT3*y1*y2*y3**2*y5/2.D0 &
               +y1**2*y2*y3*y4 +SQRT3*y1*y2**2*y3*y5/2.D0 )*DDFEA11234 +( y1*y2**3*y3 &
               +y1*y2*y3**3 +y1**3*y2*y3 )*DDFEA11123 +( y1**2*y2**2*y3 &
               +y1*y2**2*y3**2 &
               +y1**2*y2*y3**2 )*DDFEA11233  
            dds3(6,6)=( y2**3*y4**3*SQRT3 &
               -y2**3*y4**2*y5 &
               +y3**3*y4**2*y5 -5.D0/3.D0*y2**3*y4*y5**2*SQRT3 &
               +y3**3*y4**3*SQRT3 -5.D0/3.D0*y3**3*y4*y5**2*SQRT3 -y2**3*y5**3 &
               +y3**3*y5**3 -8.D0/3.D0*y1**3*y4*y5**2*SQRT3 )*DDFEA333555 &
               +( y1**4*y5**2*SQRT3/2.D0 +y2**4*y4*y5 +y2**4*y4**2*SQRT3/3.D0 &
               +y3**4*y4**2*SQRT3/3.D0 -y3**4*y4*y5 &
               -y1**4*y4**2*SQRT3/6.D0 )*DDFEA222245 &
               +( y1*y3**5 +y1*y2**5 +y2**5*y3 +y1**5*y3 +y1**5*y2 &
               +y2*y3**5 )*DDFEA133333 +( y1**4*y3*y4 -2.D0*y2**4*y3*y4 +y1**4*y2*y4 &
               +y1*y2**4*y5*SQRT3 +y1*y3**4*y4 -2.D0*y2*y3**4*y4 &
               +y1**4*y2*y5*SQRT3 -y1*y3**4*y5*SQRT3 -y1**4*y3*y5*SQRT3 &
               +y1*y2**4*y4 )*DDFEA133334 +(  -y1*y2*y3*y4**3/3.D0 &
               +y1*y2*y3*y4*y5**2 )*DDFEA123455 
            dds4(6,6)=dds3(6,6) &
               +( 2.D0/3.D0*SQRT3*y1*y2**2*y3**2*y4 -y1**2*y2**2*y3*y5 &
               -SQRT3*y1**2*y2**2*y3*y4/3.D0 &
               +y1**2*y2*y3**2*y5 -SQRT3*y1**2*y2*y3**2*y4/3.D0 )*DDFEA112335 &
               +( y1*y2**2*y3*y5**2 +y1*y2*y3**2*y5**2 +y1*y2*y3**2*y4**2 &
               +y1*y2**2*y3*y4**2 +y1**2*y2*y3*y4**2 &
               +y1**2*y2*y3*y5**2 )*DDFEA112355 
            dds2(6,6)=dds4(6,6) &
               +( y2**3*y3**2*y5 -y1**3*y2**2*y5/2.D0 -y1**2*y3**3*y5/2.D0 &
               -y2**2*y3**3*y5 &
               +y1**3*y2**2*y4*SQRT3/2.D0 -y1**2*y2**3*y4*SQRT3/2.D0 &
               +y1**3*y3**2*y5/2.D0 +y1**2*y2**3*y5/2.D0 &
               +y1**3*y3**2*y4*SQRT3/2.D0 -y1**2*y3**3*y4*SQRT3/2.D0 )*DDFEA222335 &
               +(  -y1**2*y2**2*y5**2*SQRT3/2.D0 -y1**2*y3**2*y5**2*SQRT3/2.D0 &
               -y1**2*y2**2*y4**2*SQRT3/6.D0 -y1**2*y2**2*y4*y5 &
               -2.D0/3.D0*y2**2*y3**2*y4**2*SQRT3 &
               +y1**2*y3**2*y4*y5 -y1**2*y3**2*y4**2*SQRT3/6.D0 )*DDFEA113345 &
               +( y2**2*y3**2*y5**2 +y2**2*y3**2*y4**2 +y1**2*y2**2*y5**2 &
               +y1**2*y3**2*y4**2 +y1**2*y3**2*y5**2 &
               +y1**2*y2**2*y4**2 )*DDFEA223355 
            dds3(6,6)=dds2(6,6) &
               +( y1*y2*y3**2*y4**2*SQRT3/6.D0 +y1*y2*y3**2*y4*y5 &
               +y1*y2*y3**2*y5**2*SQRT3/2.D0 &
               +2.D0/3.D0*y1**2*y2*y3*y4**2*SQRT3 -y1*y2**2*y3*y4*y5 &
               +y1*y2**2*y3*y4**2*SQRT3/6.D0 &
               +y1*y2**2*y3*y5**2*SQRT3/2.D0 )*DDFEA123345 &
               +(  -y1**3*y2**2*y5*SQRT3/2.D0 -y1**3*y2**2*y4/2.D0 &
               -y1**3*y3**2*y4/2.D0 -y1**2*y2**3*y4/2.D0 &
               +y1**3*y3**2*y5*SQRT3/2.D0 -y1**2*y3**3*y4/2.D0 &
               +y2**3*y3**2*y4 -y1**2*y2**3*y5*SQRT3/2.D0 +y2**2*y3**3*y4 &
               +y1**2*y3**3*y5*SQRT3/2.D0 )*DDFEA222334 +( 3.D0*y3**2*y4**4 &
               +5.D0/2.D0*y1**2*y5**4 +y2**2*y5**4 &
               +3.D0*y2**2*y4**4 -4.D0*y3**2*y4*y5**3*SQRT3 +y3**2*y5**4 &
               +9.D0*y1**2*y4**2*y5**2 -3.D0/2.D0*y1**2*y4**4 &
               +4.D0*y2**2*y4*y5**3*SQRT3 )*DDFEA335555 +( y1**3*y2**3 +y1**3*y3**3 &
               +y2**3*y3**3 )*DDFEA222333 
            dds4(6,6)=dds3(6,6) &
               +( y3*y4**5/5.D0 -y2*y4**4*y5*SQRT3/2.D0 -2.D0/5.D0*y1*y4**5 &
               -2.D0*y1*y4**3*y5**2 -3.D0/10.D0*y2*y5**5*SQRT3 &
               +y3*y4**3*y5**2 +y3*y4**4*y5*SQRT3/2.D0 +y2*y4**3*y5**2 &
               +3.D0/10.D0*y3*y5**5*SQRT3 +y2*y4**5/5.D0 )*DDFEA244455 &
               +( y2**5*y4 -2.D0*y1**5*y4 -SQRT3*y2**5*y5 +y3**5*y4 &
               +SQRT3*y3**5*y5 )*DDFEA222224 
            dds5(6,6)=dds4(6,6) &
               +(  -y3*y5**5*SQRT3/5.D0 +y2*y5**5*SQRT3/5.D0 &
               +y1*y4*y5**4 -7.D0/15.D0*y2*y4**5 &
               +y2*y4**4*y5*SQRT3/3.D0 -y3*y4**4*y5*SQRT3/3.D0 +y3*y4*y5**4 &
               +y2*y4*y5**4 &
               +2.D0*y1*y4**3*y5**2 -7.D0/15.D0*y3*y4**5 &
               -y1*y4**5/15.D0 )*DDFEA145555 
            dds1(6,6)=dds5(6,6) &
               +(  -SQRT3*y1*y2*y3**3*y5/2.D0 +y1**3*y2*y3*y4 &
               +SQRT3*y1*y2**3*y3*y5/2.D0 -y1*y2**3*y3*y4/2.D0 &
               -y1*y2*y3**3*y4/2.D0 )*DDFEA111234 &
               +( y3*y4**4*y5/3.D0 &
               +y3*y4**5*SQRT3/18.D0 -y2*y4**4*y5/3.D0 -y2*y4*y5**4*SQRT3/2.D0 &
               -y3*y4**2*y5**3 &
               +2.D0/9.D0*y1*y4**5*SQRT3 +y2*y4**5*SQRT3/18.D0 &
               +y2*y4**2*y5**3 -2.D0/3.D0*y1*y4**3*y5**2*SQRT3 &
               -y3*y4*y5**4*SQRT3/2.D0 )*DDFEA244555 &
               +( y1*y2*y4**2*y5**2 -3.D0/4.D0*y2*y3*y4**4 -y1*y2*y5**4 -y1*y3*y5**4 &
               +5.D0/4.D0*y2*y3*y5**4 &
               +y1*y3*y4**2*y5**2 -7.D0/2.D0*y2*y3*y4**2*y5**2 &
               -2.D0*y1*y2*y4**3*y5*SQRT3 &
               +2.D0*y1*y3*y4**3*y5*SQRT3 )*DDFEA124455 
            dds3(6,6)=dds1(6,6) +( y2**6 &
               +y1**6 +y3**6 )*DDFEA333333 +( y1*y2**4*y3 +y1**4*y2*y3 &
               +y1*y2*y3**4 )*DDFEA111123 +y1**2*y2**2*y3**2*DDFEA112233 &
               +( y1**4*y4**2 +y2**4*y4**2 +y2**4*y5**2 +y3**4*y4**2 +y1**4*y5**2 &
               +y3**4*y5**2 )*DDFEA222255 
            dds4(6,6)=dds3(6,6) +( 3.D0*y1*y3*y5**4 &
               +y1*y3*y4**4 &
               +9.D0*y2*y3*y4**2*y5**2 -3.D0/2.D0*y2*y3*y5**4 &
               -4.D0*y1*y3*y4**3*y5*SQRT3 &
               +y1*y2*y4**4 +4.D0*y1*y2*y4**3*y5*SQRT3 +3.D0*y1*y2*y5**4 &
               +5.D0/2.D0*y2*y3*y4**4 )*DDFEA134444 &
               +(  -y1*y3**2*y5**3*SQRT3/3.D0 -7.D0/3.D0*y1**2*y3*y4*y5**2 &
               +5.D0/3.D0*y1*y2**2*y4**2*y5*SQRT3 -13.D0/3.D0*y2**2*y3*y4*y5**2 &
               -4.D0/3.D0*y2*y3**2*y5**3*SQRT3 -7.D0/3.D0*y1**2*y2*y4*y5**2 &
               -16.D0/3.D0*y1*y3**2*y4*y5**2 &
               +4.D0/3.D0*y1**2*y3*y4**2*y5*SQRT3 +4.D0/3.D0*y2**2*y3*y5**3*SQRT3 &
               +3.D0*y1**2*y2*y4**3 +y2*y3**2*y4**3 +y1*y2**2*y5**3*SQRT3/3.D0 &
               +y2**2*y3*y4**3 -13.D0/3.D0*y2*y3**2*y4*y5**2 &
               -5.D0/3.D0*y1*y3**2*y4**2*y5*SQRT3 -4.D0/3.D0*y1**2*y2*y4**2*y5*SQRT3 &
               +3.D0*y1**2*y3*y4**3 &
               -16.D0/3.D0*y1*y2**2*y4*y5**2 )*DDFEA233444 
            dds5(6,6)=dds4(6,6) &
               +( 2.D0*y1*y3**2*y5**3 +4.D0*y2*y3**2*y5**3 &
               +4.D0*y2**2*y3*y4*y5**2*SQRT3 -2.D0*y1*y2**2*y5**3 &
               +y1**2*y3*y4*y5**2*SQRT3 &
               +6.D0*y1*y3**2*y4**2*y5 -6.D0*y1*y2**2*y4**2*y5 -3.D0*y1**2*y3*y4**2*y5 &
               +y1**2*y2*y4*y5**2*SQRT3 &
               +4.D0*y1*y3**2*y4*y5**2*SQRT3 -3.D0*y1**2*y2*y4**3*SQRT3 &
               -4.D0*y2**2*y3*y5**3 &
               +3.D0*y1**2*y2*y4**2*y5 -y1**2*y2*y5**3 &
               +y1**2*y3*y5**3 -3.D0*y1**2*y3*y4**3*SQRT3 &
               +4.D0*y2*y3**2*y4*y5**2*SQRT3 &
               +4.D0*y1*y2**2*y4*y5**2*SQRT3 )*DDFEA113555 
            dds2(6,6)=dds5(6,6) &
               +(  -2.D0/3.D0*y3**2*y4**4*SQRT3 -3.D0/2.D0*y1**2*y4**2*y5**2*SQRT3 &
               -3.D0/4.D0*y1**2*y5**4*SQRT3 -y2**2*y4**3*y5 &
               +7.D0/12.D0*y1**2*y4**4*SQRT3 +y3**2*y4**3*y5 &
               +3.D0*y3**2*y4*y5**3 -2.D0/3.D0*y2**2*y4**4*SQRT3 &
               -3.D0*y2**2*y4*y5**3 )*DDFEA334445 &
               +(  -3.D0*y1*y3*y4**3*y5 +2.D0/3.D0*y1*y2*y5**4*SQRT3 -y1*y3*y4*y5**3 &
               +2.D0/3.D0*y1*y3*y5**4*SQRT3 &
               +3.D0*y1*y2*y4**3*y5 -7.D0/12.D0*y2*y3*y5**4*SQRT3 &
               +3.D0/2.D0*y2*y3*y4**2*y5**2*SQRT3 +y1*y2*y4*y5**3 &
               +3.D0/4.D0*y2*y3*y4**4*SQRT3 )*DDFEA124555 &
               +( 2.D0*y3**2*y4*y5**3*SQRT3 -7.D0/2.D0*y1**2*y4**2*y5**2 &
               +y2**2*y4**2*y5**2 -y2**2*y4**4 -y3**2*y4**4 &
               -2.D0*y2**2*y4*y5**3*SQRT3 -3.D0/4.D0*y1**2*y5**4 &
               +5.D0/4.D0*y1**2*y4**4 &
               +y3**2*y4**2*y5**2 )*DDFEA334455 
            dds3(6,6)=dds2(6,6) &
               +(  -6.D0*y4**2*y5**4 +9.D0*y4**4*y5**2 +y5**6 )*DDFEA555555 &
               +( y2*y3**3*y4**2 +y2*y3**3*y5**2 +y1*y3**3*y4**2 +y1*y2**3*y4**2 &
               +y1**3*y2*y4**2 +y1*y2**3*y5**2 +y1**3*y3*y5**2 +y1**3*y3*y4**2 &
               +y1**3*y2*y5**2 +y2**3*y3*y4**2 +y1*y3**3*y5**2 &
               +y2**3*y3*y5**2 )*DDFEA233344 &
               +( y1*y2**3*y5**2*SQRT3/6.D0 -y2**3*y3*y5**2*SQRT3/3.D0 &
               -y2*y3**3*y5**2*SQRT3/3.D0 &
               +y1**3*y2*y4*y5 -y1**3*y2*y5**2*SQRT3/3.D0 -y1**3*y3*y4*y5 &
               -y1**3*y3*y5**2*SQRT3/3.D0 -y1*y3**3*y4**2*SQRT3/2.D0 &
               +y1*y3**3*y5**2*SQRT3/6.D0 -y2**3*y3*y4*y5 &
               +y2*y3**3*y4*y5 -y1*y2**3*y4**2*SQRT3/2.D0 )*DDFEA233345 &
               +(  -3.D0*y2**3*y4*y5**2 &
               +y3**3*y4**3 -3.D0*y3**3*y4*y5**2 -3.D0*y1**3*y4*y5**2 +y2**3*y4**3 &
               +y1**3*y4**3 )*DDFEA111444 +( y1*y2**3*y3**2 +y1**3*y2**2*y3 &
               +y1**2*y2**3*y3 +y1*y2**2*y3**3 +y1**2*y2*y3**3 &
               +y1**3*y2*y3**2 )*DDFEA111233 
            dds4(6,6)=dds3(6,6) &
               +( 9.D0*y4**2*y5**4 -6.D0*y4**4*y5**2 +y4**6 )*DDFEA444444 &
               +(  -5.D0/3.D0*y1*y2**2*y4**2*y5*SQRT3 &
               +y1*y2**2*y4**3 -4.D0/3.D0*y1**2*y3*y4**2*y5*SQRT3 &
               -2.D0*y1**2*y2*y4**3 -y1*y2**2*y5**3*SQRT3/3.D0 &
               +4.D0/3.D0*y2**2*y3*y4*y5**2 -4.D0/3.D0*y2**2*y3*y5**3*SQRT3 &
               -2.D0*y1**2*y3*y4**3 &
               +7.D0/3.D0*y1*y2**2*y4*y5**2 -2.D0/3.D0*y1**2*y3*y4*y5**2 &
               +y1*y3**2*y4**3 +4.D0/3.D0*y2*y3**2*y5**3*SQRT3 &
               +y1*y3**2*y5**3*SQRT3/3.D0 +4.D0/3.D0*y1**2*y2*y4**2*y5*SQRT3 &
               +4.D0/3.D0*y2*y3**2*y4*y5**2 &
               +5.D0/3.D0*y1*y3**2*y4**2*y5*SQRT3 -2.D0/3.D0*y1**2*y2*y4*y5**2 &
               +7.D0/3.D0*y1*y3**2*y4*y5**2 )*DDFEA133444 
            dds5(6,6)=dds4(6,6) &
               +(  -y1**3*y2*y4*y5 +2.D0/3.D0*y2**3*y3*y5**2*SQRT3 &
               +y1*y3**3*y4**2*SQRT3/2.D0 +y1**3*y3*y4**2*SQRT3/2.D0 &
               +y1**3*y3*y5**2*SQRT3/6.D0 +y1**3*y2*y5**2*SQRT3/6.D0 +y1**3*y3*y4*y5 &
               +y1*y2**3*y5**2*SQRT3/6.D0 +y1**3*y2*y4**2*SQRT3/2.D0 &
               +2.D0/3.D0*y2*y3**3*y5**2*SQRT3 -y1*y2**3*y4*y5 &
               +y1*y2**3*y4**2*SQRT3/2.D0 +y1*y3**3*y5**2*SQRT3/6.D0 &
               +y1*y3**3*y4*y5 )*DDFEA133345 
            ddv6(6,6)=dds5(6,6) &
               +(  -y2**2*y3*y4**2*y5 +y1**2*y3*y4*y5**2*SQRT3/3.D0 &
               +y2*y3**2*y4**2*y5 +y2*y3**2*y5**3 -y1*y2**2*y5**3 &
               +4.D0/3.D0*y2**2*y3*y4*y5**2*SQRT3 &
               +4.D0/3.D0*y2*y3**2*y4*y5**2*SQRT3 -y1*y2**2*y4**2*y5 &
               +4.D0/3.D0*y1*y3**2*y4*y5**2*SQRT3 -y2**2*y3*y5**3 +y1*y3**2*y5**3 &
               +y1**2*y2*y4*y5**2*SQRT3/3.D0 -y1**2*y2*y4**3*SQRT3 &
               +y1*y3**2*y4**2*y5 -y1**2*y3*y4**3*SQRT3 &
               +4.D0/3.D0*y1*y2**2*y4*y5**2*SQRT3 )*DDFEA233445 &
               +( y2*y3**4*y4*SQRT3 -y1**4*y2*y5 &
               +y2**4*y3*y4*SQRT3 -y1**4*y3*y4*SQRT3 +y2*y3**4*y5 -2.D0*y1*y2**4*y5 &
               +2.D0*y1*y3**4*y5 -y1**4*y2*y4*SQRT3 &
               +y1**4*y3*y5 -y2**4*y3*y5 )*DDFEA233335 +( y2**2*y3**4 +y1**4*y3**2 &
               +y1**2*y2**4 +y2**4*y3**2 +y1**2*y3**4 &
               +y1**4*y2**2 )*DDFEA222233   
 
            ddv0 =f2a*2d0+f3a*6d0*coro+f4a*12d0*coro**2+f5a*20d0*coro**3 &
               +f6a*30d0*coro**4+f7a*42d0*coro**5+f8a*56d0*coro**6
 
            ddv(6,6)=ddv0 +ddv1(6,6) +ddv2(6,6) &
               +ddv3(6,6) +ddv4(6,6) +ddv5(6,6) +ddv6(6,6)
 
            ! ddv(i,j) = ddv(j,i)
            do i = 1, 5
               do j = i+1, 6 ! j != i
                  ddv(j,i) = ddv(i,j)
               end do
            end do
 
         end if ! flag > 1 ( hessian )
      end if ! flag > 0 (gradient and hessian)
 
   end subroutine nh3_pot_work
 
   !> @todo
   subroutine nh3_inter_to_work(work, dwdr, ddwdrdr, r, flag)
      implicit none
      real(8), intent(in)  :: r(6)
      integer, intent(in)  :: flag
      real(8), intent(out) :: work(6)
      real(8), intent(out) :: dwdr(6,6), ddwdrdr(6,6,6)
 
      real(8) ::  alpha1, alpha2, alpha3
      real(8) ::  alpha_half, rhoe, sinrho
 
      real(8) :: yterm(3) 
 
      if (flag < 0) then
         write(*,*) 'flag should be 0/1/2'
         stop
      end if
 
      ! set parameters (only for the first time) 
      call set_param()
 
      ! bond lengths and angles
      ! r14  = r(1) ;  r24    = r(2) ;  r34    = r(3)
      alpha1 = r(4) ;  alpha2 = r(5) ;  alpha3 = r(6)
 
      ! calculate working coordinates
      yterm(1:3) = exp(-aa1*(r(1:3)-re14)) 
      work(1:3) = 1.0d+00 - yterm(1:3)
 
      work(4) = (2.d0*alpha1-alpha2-alpha3)/SQRT6
      work(5) = (alpha2-alpha3)/SQRT2
 
      alpha_half=(alpha1+alpha2+alpha3)/6.d0
      sinrho=2d0*sin(alpha_half)/SQRT3
 
      if ( sinrho .ge. 1.0d0 ) then 
         sinrho=1.d0 
      end if
 
      rhoe=pi*rhoedg/1.8d+02
      work(6)=sin(rhoe)-sinrho
 
      if (flag > 0) then
         ! -- 1st derivatives
         ! y1,y2,y3 are functions of bond lengths
         ! each related to one length
         dwdr(1:3,1:6) = 0d0
         dwdr(1,1) = aa1*yterm(1)
         dwdr(2,2) = aa1*yterm(2)
         dwdr(3,3) = aa1*yterm(3)
 
         ! y4,y5,coro are functions of bond angles
         dwdr(4:6,1:3) = 0d0
 
         ! dy4/dr(4:6)
         dwdr(4,4) = 2d0/SQRT6
         dwdr(4,5:6) = -1d0/SQRT6
 
         ! dy5/dr(4:6)
         dwdr(5,4) = 0d0
         dwdr(5,5) = 1d0/SQRT2
         dwdr(5,6) = -1d0/SQRT2
 
         ! dcoro/dr(4:6)
         dwdr(6,4) = -SQRT3*cos(alpha_half)/9d0 
         dwdr(6,5) = dwdr(6,4)
         dwdr(6,6) = dwdr(6,4)
 
         if (flag > 1) then
            ! -- 2nd derivatives
            ! set all to zero
            ddwdrdr(1:6,1:6,1:6) = 0d0
 
            ! ddy1, ddy2, ddy3
            ddwdrdr(1,1,1) = -aa1*dwdr(1,1) ! -aa1*aa1*yterm(1)
            ddwdrdr(2,2,2) = -aa1*dwdr(2,2) ! -aa1*aa1*yterm(2)
            ddwdrdr(3,3,3) = -aa1*dwdr(3,3) ! -aa1*aa1*yterm(3)
 
            ! ddy4, ddy5 are zero
            ! ddcoro
            ddwdrdr(6,4:6,4:6) = sinrho/36d0
         end if
      end if
 
      return
   end subroutine nh3_inter_to_work
 
   !> set parameters for calculation at the first time calling 
   subroutine set_param()
      implicit none
      integer,parameter :: parmax = 307
      double precision ::  param(parmax)
 
      if (.not. param_saved) then
         ! Edit by wyz in 2022.05.04
         rhoedg    =       112.0966d0
         re14      =      1.0103131d0
         aa1       =           2.15d0
         ve        =  -12419377.330d0
         f1a       =           0.00d0
         f2a       =        324077.d0
         f3a       =       -401792.d0
         f4a       =       1133218.d0
         f5a       =      -2662929.d0
         f6a       =       4561862.d0
         f7a       =           0.00d0
         f8a       =           0.00d0
         f1a1      =       -33671.7d0
         f2a1      =         41558.d0
         f3a1      =       -349715.d0
         f4a1      =       1246471.d0
         f5a1      =      -2439730.d0
         f6a1      =           0.00d0
         f0a11     =        38730.1d0
         f1a11     =        -17346.d0
         f2a11     =         57866.d0
         f3a11     =       -462520.d0
         f4a11     =        979992.d0
         f0a12     =        -403.29d0
         f1a12     =          4673.d0
         f2a12     =         41577.d0
         f3a12     =       -199083.d0
         f4a12     =        606382.d0
         f0a14     =        -3651.0d0
         f1a14     =        -17081.d0
         f2a14     =        -44594.d0
         f3a14     =        156439.d0
         f4a14     =       -427396.d0
         f0a44     =       16833.16d0
         f1a44     =         67883.d0
         f2a44     =       -110704.d0
         f3a44     =        324179.d0
         f4a44     =       -536600.d0
         f0a111    =          277.2d0
         f1a111    =         -9581.d0
         f2a111    =         61987.d0
         f3a111    =       -262799.d0
         f0a112    =         -291.8d0
         f1a112    =          2294.d0
         f2a112    =         26099.d0
         f3a112    =       -103127.d0
         f0a114    =        -2068.4d0
         f1a114    =         -6903.d0
         f2a114    =        -77502.d0
         f3a114    =        233336.d0
         f0a123    =         -186.3d0
         f1a123    =          4547.d0
         f2a123    =         16600.d0
         f3a123    =        -88441.d0
         f0a124    =        1956.59d0
         f1a124    =          5405.d0
         f2a124    =         -2685.d0
         f3a124    =           0.00d0
         f0a144    =        -1648.6d0
         f1a144    =         -6830.d0
         f2a144    =        -19840.d0
         f3a144    =         85887.d0
         f0a155    =        -2938.1d0
         f1a155    =         -9414.d0
         f2a155    =        -14906.d0
         f3a155    =       -119102.d0
         f0a455    =         1551.7d0
         f1a455    =        -57777.d0
         f2a455    =         35451.d0
         f3a455    =        388038.d0
         f0a1111   =         3694.9d0
         f1a1111   =         -6990.d0
         f2a1111   =         19038.d0
         f0a1112   =         -545.7d0
         f1a1112   =           0.00d0
         f2a1112   =         14424.d0
         f0a1114   =         -678.5d0
         f1a1114   =           0.00d0
         f2a1114   =        -59594.d0
         f0a1122   =         -189.9d0
         f1a1122   =           0.00d0
         f2a1122   =         14272.d0
         f0a1123   =         -254.3d0
         f1a1123   =          2729.d0
         f2a1123   =           0.00d0
         f0a1124   =          895.1d0
         f1a1124   =          2070.d0
         f2a1124   =         -7392.d0
         f0a1125   =         1399.4d0
         f1a1125   =          5090.d0
         f2a1125   =           0.00d0
         f0a1144   =        -1150.3d0
         f1a1144   =         -7205.d0
         f2a1144   =           0.00d0
         f0a1155   =        -2615.7d0
         f1a1155   =         -8539.d0
         f2a1155   =           0.00d0
         f0a1244   =          508.2d0
         f1a1244   =           0.00d0
         f2a1244   =        -27657.d0
         f0a1255   =         1321.9d0
         f1a1255   =           921.d0
         f2a1255   =           0.00d0
         f0a1444   =         -765.1d0
         f1a1444   =         -3922.d0
         f2a1444   =           0.00d0
         f0a1455   =        -1581.9d0
         f1a1455   =         -8474.d0
         f2a1455   =        -93060.d0
         f0a4444   =          395.6d0
         f1a4444   =          8421.d0
         f2a4444   =         13895.d0
         f0a44444  =         -127.4d0
         f1a44444  =          1827.d0
         f2a44444  =          9830.d0
         f0a33455  =          -714.d0
         f1a33455  =           0.00d0
         f2a33455  =       -108147.d0
         f0a33445  =          809.8d0
         f1a33445  =          3330.d0
         f2a33445  =           0.00d0
         f0a33345  =          1791.d0
         f1a33345  =          4219.d0
         f2a33345  =           0.00d0
         f0a33344  =        -2410.4d0
         f1a33344  =        -11427.d0
         f2a33344  =         50287.d0
         f0a33334  =         -183.1d0
         f1a33334  =           0.00d0
         f2a33334  =           0.00d0
         f0a33333  =         1749.4d0
         f1a33333  =           0.00d0
         f2a33333  =           0.00d0
         f0a25555  =          -52.8d0
         f1a25555  =         -1133.d0
         f2a25555  =        -17384.d0
         f0a24455  =          584.9d0
         f1a24455  =           0.00d0
         f2a24455  =        -53324.d0
         f0a24445  =         1070.1d0
         f1a24445  =          4655.d0
         f2a24445  =        -37568.d0
         f0a23333  =         -374.9d0
         f1a23333  =           0.00d0
         f2a23333  =           0.00d0
         f0a13455  =          1127.d0
         f1a13455  =          7287.d0
         f2a13455  =         27647.d0
         f0a13445  =         -395.0d0
         f1a13445  =         -2037.d0
         f2a13445  =         21076.d0
         f0a13345  =         -679.2d0
         f1a13345  =           0.00d0
         f2a13345  =           0.00d0
         f0a12355  =           221.d0
         f1a12355  =           0.00d0
         f2a12355  =           0.00d0
         f0a11334  =          624.3d0
         f1a11334  =          2211.d0
         f2a11334  =           0.00d0
         f0a11333  =           0.00d0
         f1a11333  =           0.00d0
         f2a11333  =           0.00d0
         f0a11255  =          807.3d0
         f1a11255  =           0.00d0
         f2a11255  =        -27721.d0
         f0a11245  =           0.00d0
         f1a11245  =         -3018.d0
         f2a11245  =         43059.d0
         f0a11234  =           0.00d0
         f1a11234  =           0.00d0
         f2a11234  =           0.00d0
         f0a11233  =           0.00d0
         f1a11233  =           0.00d0
         f2a11233  =           0.00d0
         f0a11135  =           0.00d0
         f1a11135  =         -7232.d0
         f2a11135  =           0.00d0
         f0a11134  =          245.9d0
         f1a11134  =           0.00d0
         f2a11134  =           0.00d0
         f0a11123  =          -400.d0
         f1a11123  =           0.00d0
         f2a11123  =           0.00d0
         f0a555555 =           0.00d0
         f1a555555 =          2257.d0
         f2a555555 =         24161.d0
         f0a444444 =         127.07d0
         f1a444444 =          4563.d0
         f2a444444 =         42425.d0
         f0a335555 =          121.1d0
         f1a335555 =           0.00d0
         f2a335555 =        -19542.d0
         f0a334455 =           867.d0
         f1a334455 =          4471.d0
         f2a334455 =           0.00d0
         f0a334445 =         -1066.d0
         f1a334445 =         -6945.d0
         f2a334445 =           0.00d0
         f0a333555 =           560.d0
         f1a333555 =         -4979.d0
         f2a333555 =           0.00d0
         f0a333333 =           619.d0
         f1a333333 =          7177.d0
         f2a333333 =           0.00d0
         f0a244555 =          -648.d0
         f1a244555 =        -14993.d0
         f2a244555 =           0.00d0
         f0a244455 =         -493.8d0
         f1a244455 =         -8199.d0
         f2a244455 =         56324.d0
         f0a233445 =          1130.d0
         f1a233445 =           0.00d0
         f2a233445 =           0.00d0
         f0a233444 =          -440.d0
         f1a233444 =         -3746.d0
         f2a233444 =           0.00d0
         f0a233345 =           0.00d0
         f1a233345 =           0.00d0
         f2a233345 =       -144810.d0
         f0a233344 =           0.00d0
         f1a233344 =           0.00d0
         f2a233344 =           0.00d0
         f0a233335 =           0.00d0
         f1a233335 =           0.00d0
         f2a233335 =           0.00d0
         f0a223355 =           0.00d0
         f1a223355 =         11246.d0
         f2a223355 =           0.00d0
         f0a222335 =           0.00d0
         f1a222335 =           0.00d0
         f2a222335 =           0.00d0
         f0a222334 =          -578.d0
         f1a222334 =           0.00d0
         f2a222334 =           0.00d0
         f0a222333 =           0.00d0
         f1a222333 =           0.00d0
         f2a222333 =           0.00d0
         f0a222255 =         -1648.d0
         f1a222255 =           0.00d0
         f2a222255 =           0.00d0
         f0a222245 =         -2032.d0
         f1a222245 =        -34578.d0
         f2a222245 =        211698.d0
         f0a222233 =           0.00d0
         f1a222233 =           0.00d0
         f2a222233 =           0.00d0
         f0a222224 =          -458.d0
         f1a222224 =          5104.d0
         f2a222224 =           0.00d0
         f0a145555 =         -1719.d0
         f1a145555 =        -27172.d0
         f2a145555 =         82494.d0
         f0a134444 =         -336.8d0
         f1a134444 =         -4362.d0
         f2a134444 =           0.00d0
         f0a133444 =           661.d0
         f1a133444 =           0.00d0
         f2a133444 =           0.00d0
         f0a133345 =           0.00d0
         f1a133345 =           0.00d0
         f2a133345 =           0.00d0
         f0a133334 =           0.00d0
         f1a133334 =           0.00d0
         f2a133334 =           0.00d0
         f0a133333 =           0.00d0
         f1a133333 =          5970.d0
         f2a133333 =        -72281.d0
         f0a124555 =          1027.d0
         f1a124555 =         10494.d0
         f2a124555 =           0.00d0
         f0a124455 =           254.d0
         f1a124455 =           0.00d0
         f2a124455 =           0.00d0
         f0a123455 =          -733.d0
         f1a123455 =           0.00d0
         f2a123455 =           0.00d0
         f0a123345 =           0.00d0
         f1a123345 =           0.00d0
         f2a123345 =           0.00d0
         f0a113555 =          -828.d0
         f1a113555 =         -2148.d0
         f2a113555 =           0.00d0
         f0a113345 =           0.00d0
         f1a113345 =           0.00d0
         f2a113345 =           0.00d0
         f0a112355 =           0.00d0
         f1a112355 =           0.00d0
         f2a112355 =           0.00d0
         f0a112335 =           0.00d0
         f1a112335 =           0.00d0
         f2a112335 =           0.00d0
         f0a112233 =           0.00d0
         f1a112233 =           0.00d0
         f2a112233 =           0.00d0
         f0a111444 =          -266.d0
         f1a111444 =           0.00d0
         f2a111444 =           0.00d0
         f0a111234 =           0.00d0
         f1a111234 =           0.00d0
         f2a111234 =           0.00d0
         f0a111233 =           0.00d0
         f1a111233 =           0.00d0
         f2a111233 =           0.00d0
         f0a111123 =           0.00d0
         f1a111123 =           0.00d0
         f2a111123 =           0.00d0
         ! call read_param_impl(param,parmax)
         ! ...
         !-------------------------------
         param_saved = .true.
      end if
   end subroutine set_param

!!!   NEEDN'T IT ANYMORE !!!
!!!   !-----------------------------------------------------------
!!!   !> read parameters from .inp file which is format-fixed
!!!   !! @param [in] n, size of parameters
!!!   !! @param [out] param, parameters read from file
!!!   !! @param [out] parnam, names of parameters
!!!   !! @param [out] ivar, status of parameter, 0 is OK
!!!   subroutine read_param_impl(param, n)
!!!      implicit none
!!!      integer, intent(in)  :: n
!!!      real(8), intent(out) :: param(n)
!!! 
!!!      character(len=10) :: parnam(n)
!!!      integer           :: ivar(n)
!!! 
!!!      ! input file unit
!!!      integer :: f_inp = 10
!!!      character(len=13) :: inpfile
!!!      logical :: exists
!!!      integer :: f_out = 6 ! print to screen
!!! 
!!!      ! input title of the job from the first four input lines 
!!!      character(len=80) :: longlabel 
!!! 
!!!      character(len=10) :: label          ! Temporary label 
!!!      integer           :: ivartmp
!!!      double precision  :: paramtmp
!!!      integer :: i
!!! 
!!!      ! open file for read
!!!      inpfile = 'nh3_param.inp'
!!!      inquire(file=inpfile, exist=exists)
!!!      if (.not. exists) then
!!!         write(f_out,*) 'ERROR: input file ', inpfile, ' NOT FOUND.'
!!!         stop
!!!      end if
!!!      open (f_inp, file=inpfile) 
!!! 
!!!      ! input potential parameters 
!!!      ! skip head lines without parameters
!!!      call skiplines(f_inp,14)  
!!!      ! read parameters for N times
!!!      do i = 1, n
!!!         read (f_inp,"(a80)") longlabel
!!!         read (longlabel,"(a10,i4,d18.8)") label,ivartmp,paramtmp
!!!         if (ivartmp/=-1) then 
!!!            parnam(i) = label 
!!!            ivar(i)   = ivartmp
!!!            param(i)  = paramtmp
!!!         else 
!!!            write(f_out,"('Wrong number of parameters (',I6,'), has to be (',I6,')')") i, n
!!!            stop 'Too few parameters'
!!!         endif 
!!!      enddo 
!!! 
!!!      read (f_inp,"(a80)") longlabel
!!!      !write(6,"(a80)") longlabel ! 'end -1 -1'
!!!      read (longlabel,"(a10,i4,d18.8)") label,ivartmp,paramtmp
!!!      if (ivartmp/=-1) then 
!!!         write(f_out,"('Wrong number of parameters (',I6,'), has to be (',I6,')')") i, n
!!!         stop 'Too many parameters'
!!!      endif 
!!! 
!!!      ! close file
!!!      close(f_inp)
!!! 
!!!   end subroutine read_param_impl
!!! 
!!!   !> Skip n lines  in the input file 
!!!   !! @param [in] inpunit, input file unit
!!!   !! @param [in] n, number of lines to skip
!!!   subroutine skiplines( inpunit,n )
!!!      integer,intent(in) :: n,inpunit
!!!      character(len=80)  :: label
!!!      integer            :: i0
!!! 
!!!      do i0=1,n
!!!         read  (inpunit,"(a80)") label 
!!!      enddo
!!! 
!!!   end subroutine skiplines
 
   !---------------------------------------------------------------------
   !> coordinate conversion: cartesian -> internal (3 lengths, 3 angles)
   !! @param [out] r, internal coordinates (size = 6)
   !! @param [out] drdx, 1st derivatives dr(ir)/dx(ix)
   !! @param [out] ddrdxdx, 2nd derivatives d(dr(ir)/dx(ix))/dx(jx)
   !! @param [in] x, cartesian coordinates (x,y,z) of N,Ha,Hb,Hc, 
   !!                x(1:3) - N, x(4:6) - Ha, x(7:9) - Hb, x(10:12) - Hc
   !! @param [in] flag, control flow according to which value is needed
   !---------------------------------------------------------------------
   subroutine nh3_cart_to_inter(r, drdx, ddrdxdx, x, flag)
      implicit none
      real(8), intent(in) :: x(12)
      integer, intent(in) :: flag
      real(8), intent(out) :: r(6), drdx(6,12), ddrdxdx(6,12,12)
 
      !    Ha   
      !      \  
      !       N--Hb    : 3 bond lengths, 3 bond angles
      !      /
      !     Hc
      ! internal coordinates (bond lengths and bond angles)
      real(8) :: rnha, rnhb, rnhc, theta_a, theta_b, theta_c 
      real(8),dimension(6) :: dradx, drbdx, drcdx
      real(8),dimension(6,6) :: ddradxdx, ddrbdxdx, ddrcdxdx
      !real(8),dimension(9) :: dtadx, dtbdx, dtcdx
      !real(8),dimension(9,9) :: ddtadxdx, ddtbdxdx, ddtcdxdx
 
      ! vectors for bonds
      real(8), dimension(3) :: vnha, vnhb, vnhc  
 
      ! unit vectors
      real(8), dimension(3) :: unha, unhb, unhc
 
      ! cosine and sine of angles 
      real(8) :: costa, sinta, costb, sintb, costc, sintc 
 
      ! derivatives for cosine angles
      real(8), dimension(12) :: dcostadx, dcostbdx, dcostcdx
      real(8), dimension(12,12) :: ddcostadxdx, ddcostbdxdx, ddcostcdxdx
 
      ! index
      integer :: i
 
      if (flag < 0) then
         write(*,*) 'flag should be 0/1/2'
         stop
      end if
 
      ! --- calculate internal coordinates ---
 
      ! label N-1, Ha-2, Hb-3, Hc-4
      vnha(1) = x(4)-x(1)
      vnha(2) = x(5)-x(2)
      vnha(3) = x(6)-x(3)
 
      vnhb(1) = x(7)-x(1)
      vnhb(2) = x(8)-x(2)
      vnhb(3) = x(9)-x(3)
 
      vnhc(1) = x(10)-x(1)
      vnhc(2) = x(11)-x(2)
      vnhc(3) = x(12)-x(3)
 
      ! bond lengths
      rnha = sqrt(dot_product(vnha,vnha))
      rnhb = sqrt(dot_product(vnhb,vnhb))
      rnhc = sqrt(dot_product(vnhc,vnhc))
 
      ! check length result
      if (rnha < epsilon(0d0)) then
         print *, 'ERROR: bond length of N-Ha ',rnha,' is too small!'
         stop
      end if
      if (rnhb < epsilon(0d0)) then
         print *, 'ERROR: bond length of N-Hb ',rnhb,' is too small!'
         stop
      end if
      if (rnhc < epsilon(0d0)) then
         print *, 'ERROR: length of N-Hc ',rnhc,' is too small!'
         stop
      end if
 
      ! unit vector
      unha(1:3) = vnha(1:3)/rnha
      unhb(1:3) = vnhb(1:3)/rnhb
      unhc(1:3) = vnhc(1:3)/rnhc
 
      ! bond angles
      costa = dot_product(unha, unhb)
      costb = dot_product(unhb, unhc)
      costc = dot_product(unhc, unha)
 
      ! check dotproduct result
      ! @todo need check? when shall error occur?
      if (costa > 1d0 ) costa = 1d0
      if (costa < -1d0) costa = -1d0
      if (costb > 1d0 ) costb = 1d0
      if (costb < -1d0) costb = -1d0
      if (costc > 1d0 ) costc = 1d0
      if (costc < -1d0) costc = -1d0
 
      theta_a = acos(costa)
      theta_b = acos(costb)
      theta_c = acos(costc)
 
      ! internal coordinates r(6)
      r(1) = rnha
      r(2) = rnhb
      r(3) = rnhc
      r(4) = theta_a
      r(5) = theta_b
      r(6) = theta_c  
 
      if (flag > 0) then
         ! --- calculate 1st derivatives ---
 
         ! drdx(1,:) -- drnha/dx
         call calc_deriv_r12_x(dradx, ddradxdx, unha, rnha)
         drdx(1,1:3)  = dradx(1:3)
         drdx(1,4:6)  = dradx(4:6)
         drdx(1,7:12) = 0d0
 
         ! drdx(2,:) -- drnhb/dx
         call calc_deriv_r12_x(drbdx, ddrbdxdx, unhb, rnhb)
         drdx(2,1:3)   = drbdx(1:3)
         drdx(2,7:9)   = drbdx(4:6)
         drdx(2,4:6)   = 0d0
         drdx(2,10:12) = 0d0
 
         ! drdx(3,:) -- drnhc/dx
         call calc_deriv_r12_x(drcdx, ddrcdxdx, unhc, rnhc)
         drdx(3,1:3)   = drcdx(1:3)
         drdx(3,10:12) = drcdx(4:6)
         drdx(3,4:9)   = 0d0
 
         ! drdx(4,:) -- dtheta_a/dx = (-1/sinta) * dcosta/dx 
         dcostadx(4:6) = (unhb(1:3) - costa*unha(1:3))/rnha 
         dcostadx(7:9) = (unha(1:3) - costa*unhb(1:3))/rnhb
         dcostadx(1:3) = -dcostadx(4:6) - dcostadx(7:9)
         dcostadx(10:12) = 0d0
 
         sinta = sqrt(1d0 - costa*costa)
         drdx(4, 4:9) = -dcostadx(4:9)/sinta
         drdx(4, 1:3) = -drdx(4,4:6) - drdx(4,7:9)
         drdx(4,10:12) = 0d0
 
         ! drdx(5,:) -- dtheta_b/dx = (-1/sintb) * dcostb/dx 
         dcostbdx(7:9) = (unhc(1:3) - costb*unhb(1:3))/rnhb 
         dcostbdx(10:12) = (unhb(1:3) - costb*unhc(1:3))/rnhc
         dcostbdx(1:3) = -dcostbdx(7:9) - dcostbdx(10:12)
         dcostbdx(4:6) = 0d0
 
         sintb = sqrt(1d0 - costb*costb)
         drdx(5, 7:12) = -dcostbdx(7:12)/sintb
         drdx(5, 1:3) = -drdx(5,7:9) - drdx(5,10:12)
         drdx(5,4:6) = 0d0
 
         ! drdx(6,:) -- dtheta_c/dx = (-1/sintc) * dcostc/dx 
         dcostcdx(10:12) = (unha(1:3) - costc*unhc(1:3))/rnhc 
         dcostcdx(4:6) = (unhc(1:3) - costc*unha(1:3))/rnha
         dcostcdx(1:3) = -dcostcdx(4:6) - dcostcdx(10:12)
         dcostcdx(7:9) = 0d0
 
         sintc = sqrt(1d0 - costc*costc)
         drdx(6, 4:6) = -dcostcdx(4:6)/sintc
         drdx(6, 10:12) = -dcostcdx(10:12)/sintc
         drdx(6, 1:3) = -drdx(6,4:6) - drdx(6,10:12)
         drdx(6,7:9) = 0d0
 
         if (flag > 1) then
            ! --- calculate 2nd derivatives --- @todo
 
            ! set all to zero
            ddrdxdx = 0d0
 
            ddrdxdx(1,1:3,1:3) = ddradxdx(1:3,1:3)
            ddrdxdx(1,1:3,4:6) = ddradxdx(1:3,4:6)
            ddrdxdx(1,4:6,1:3) = ddradxdx(4:6,1:3)
            ddrdxdx(1,4:6,4:6) = ddradxdx(4:6,4:6)
 
            ddrdxdx(2,1:3,1:3) = ddrbdxdx(1:3,1:3)
            ddrdxdx(2,1:3,7:9) = ddrbdxdx(1:3,4:6)
            ddrdxdx(2,7:9,1:3) = ddrbdxdx(4:6,1:3)
            ddrdxdx(2,7:9,7:9) = ddrbdxdx(4:6,4:6)
 
            ddrdxdx(3,1:3,1:3)     = ddrcdxdx(1:3,1:3)
            ddrdxdx(3,1:3,10:12)   = ddrcdxdx(1:3,4:6)
            ddrdxdx(3,10:12,1:3)   = ddrcdxdx(4:6,1:3)
            ddrdxdx(3,10:12,10:12) = ddrcdxdx(4:6,4:6)
 
            ! ddrdxdx(4,:,:) -- d(dta/dx)dx
            ddcostadxdx = 0d0
            do i = 1,9
               ddcostadxdx(4:6,i) = (-ddrdxdx(2,1:3,i) + costa*ddrdxdx(1,1:3,i) &
                  -dcostadx(i)*unha(1:3) - dcostadx(4:6)*drdx(1,i))/rnha
               ddcostadxdx(7:9,i) = (-ddrdxdx(1,1:3,i) + costa*ddrdxdx(2,1:3,i) &
                  -dcostadx(i)*unhb(1:3) - dcostadx(7:9)*drdx(2,i))/rnhb
               ddcostadxdx(1:3,i) = -ddcostadxdx(4:6,i) -ddcostadxdx(7:9,i)
 
               ddrdxdx(4,4:6,i) = -(costa*drdx(4,4:6)*drdx(4,i) + ddcostadxdx(4:6,i))/sinta
               ddrdxdx(4,7:9,i) = -(costa*drdx(4,7:9)*drdx(4,i) + ddcostadxdx(7:9,i))/sinta
               ddrdxdx(4,1:3,i) = -ddrdxdx(4,4:6,i)-ddrdxdx(4,7:9,i)
            end do
 
            ! ddrdxdx(5,:,:) -- d(dtb/dx)dx
            ddcostbdxdx = 0d0
            do i = 1,12
               ddcostbdxdx(7:9,i) = (-ddrdxdx(3,1:3,i) + costb*ddrdxdx(2,1:3,i) &
                  -dcostbdx(i)*unhb(1:3) - dcostbdx(7:9)*drdx(2,i))/rnhb
               ddcostbdxdx(10:12,i) = (-ddrdxdx(2,1:3,i) + costb*ddrdxdx(3,1:3,i) &
                  -dcostbdx(i)*unhc(1:3) - dcostbdx(10:12)*drdx(3,i))/rnhc
               ddcostbdxdx(1:3,i) = -ddcostbdxdx(7:9,i) -ddcostbdxdx(10:12,i)
 
               ddrdxdx(5,7:9,i) = -(costb*drdx(5,7:9)*drdx(5,i) + ddcostbdxdx(7:9,i))/sintb
               ddrdxdx(5,10:12,i) = -(costb*drdx(5,10:12)*drdx(5,i) + ddcostbdxdx(10:12,i))/sintb
               ddrdxdx(5,1:3,i) = -ddrdxdx(5,7:9,i)-ddrdxdx(5,10:12,i)
            end do
 
            ! ddrdxdx(6,:,:) -- d(dtc/dx)dx
            ddcostcdxdx = 0d0
            do i = 1,12
               ddcostcdxdx(10:12,i) = (-ddrdxdx(1,1:3,i) + costc*ddrdxdx(3,1:3,i) &
                  -dcostcdx(i)*unhc(1:3) - dcostcdx(10:12)*drdx(3,i))/rnhc
               ddcostcdxdx(4:6,i) = (-ddrdxdx(3,1:3,i) + costc*ddrdxdx(1,1:3,i) &
                  -dcostcdx(i)*unha(1:3) - dcostcdx(4:6)*drdx(1,i))/rnha
               ddcostcdxdx(1:3,i) = -ddcostcdxdx(4:6,i) -ddcostcdxdx(10:12,i)
 
               ddrdxdx(6,4:6,i) = -(costc*drdx(6,4:6)*drdx(6,i) + ddcostcdxdx(4:6,i))/sintc
               ddrdxdx(6,10:12,i) = -(costc*drdx(6,10:12)*drdx(6,i) + ddcostcdxdx(10:12,i))/sintc
               ddrdxdx(6,1:3,i) = -ddrdxdx(6,4:6,i)-ddrdxdx(6,10:12,i)
            end do
         end if
      end if
 
   end subroutine nh3_cart_to_inter
 
   !---------------------------------------------------------------------
   !> calculate dv/dx from dv/dy, where y is a function of x
   !! @param [out] dvdx, 1st derivatives of v respect to x
   !! @param [out] ddvdxdx, 2nd derivates of v respect to x
   !! @param [in] dvdy, 1st derivatives of v respect to y
   !! @param [in] ddvdydy, 2nd derivates of v respect to y
   !! @param [in] dydx, 1st derivatives of y respect to x
   !! @param [in] ddydxdx, 2nd derivatives of y respect to x 
   !! @param [in] nx, dimension of coordinates x
   !! @param [in] ny, dimension of coordinates y
   !! @todo control flow
   !---------------------------------------------------------------------
   subroutine calc_dvdx_from_dvdy(dvdx, ddvdxdx, dvdy, ddvdydy, dydx, ddydxdx, nx, ny, flag)
      implicit none
      integer, intent(in) :: nx, ny
      real(8), intent(in) :: dvdy(ny), ddvdydy(ny,ny)
      real(8), intent(in) :: dydx(ny,nx), ddydxdx(ny,nx,nx)
      integer, intent(in) :: flag
      real(8), intent(out):: dvdx(nx), ddvdxdx(nx,nx)
 
      integer :: ix, iy, jx, jy
      real(8) :: dvtemp, ddvtemp, ddvterm
 
      if (flag <= 0) then
         write(*,*)'this subroutine need flag > 0, 1 for 1st deriv, 2 for 1st and 2nd deriv'
         return
      end if
 
      ! dv/dx(ix) = sum_{iy}(dvdy(iy) * dydx(iy,ix))
      do ix = 1, nx
         dvtemp = 0d0
         do iy = 1, ny
            dvtemp = dvtemp + dvdy(iy)*dydx(iy,ix)
         end do ! loop iy = 1, ny
         dvdx(ix) = dvtemp
      end do ! loop ix = 1, nx
 
      if (flag > 1) then
         !   d(dv/dx(ix))/dx(jx) 
         ! = d/dx(jx)(sum_{iy}(dvdy(iy) * dydx(iy,ix))
         ! = sum_{iy}(dvdy(iy) * ddydxdx(iy,ix,jx)) 
         !  + sum_{iy}[dydx(iy,ix) * sum_{jy}(dydx(jy,jx)*ddvdydy(iy,jy))]
         do ix = 1, nx
            do jx = ix, nx
               ddvtemp = 0d0
               do iy = 1, ny
                  ddvterm = 0d0
                  do jy = 1, ny
                     ddvterm = ddvterm + ddvdydy(iy,jy)*dydx(jy,jx)
                  end do
                  ddvtemp=ddvtemp + ddvterm*dydx(iy,ix) + dvdy(iy)*ddydxdx(iy,ix,jx)
               end do
               ddvdxdx(ix,jx)=ddvtemp
               ddvdxdx(jx,ix)=ddvtemp
            end do ! loop jx = ix,nx
         end do  ! loop ix = 1, nx
      end if
 
      return
   end subroutine calc_dvdx_from_dvdy
 
   !---------------------------------------------------------------------
   !> calculate derivatives of length 1-2 to point in cartesian coordinates
   !! x(1:3) -- x,y,z of atom 1, x(4:6) -- x,y,z of atom 2
   !! @param [out] dr12, 1st derivatives, size is (6)
   !! @param [out] ddr12, 2nd derivatives, size is (6,6)
   !! @param [in] u12, unit vector of bond, size is (3)
   !! @param [in] r12, bond length
   !---------------------------------------------------------------------
   subroutine calc_deriv_r12_x(dr12, ddr12, u12, r12)
      implicit none
 
      real(8), dimension(3) :: u12 ! unit vector
      real(8), dimension(6) :: dr12 
      real(8), dimension(6,6) :: ddr12
      real(8) :: r12
 
      integer :: i, j
 
      ! 1st derivatives 
      dr12(1:3) = -u12(1:3)
      dr12(4:6) = u12(1:3)
 
      ! 2nd derivatives
      ! ddrdx1dx1
      do i = 1,3
         ddr12(i,i) = (1d0-u12(i)*u12(i))/r12
         do j = i+1,3
            ddr12(i,j) = -u12(i)*u12(j)/r12
            ddr12(j,i) = ddr12(i,j)
         end do
      end do
 
      ! ddrdx1dx2
      ddr12(1:3,4:6) = -ddr12(1:3,1:3)
      ! ddrdx2dx1
      ddr12(4:6,1:3) = -ddr12(1:3,1:3)
      ! ddrdx2dx2
      ddr12(4:6,4:6) = ddr12(1:3,1:3)
 
   end subroutine calc_deriv_r12_x
 
   !------------------------------------------------------------------
   !> convert internal coordinates of NH3 to cartesian coordinates
   !! @param [out] x, cartesian x,y,z of N1,H2,H3,H4, unit == r(1:3) 
   !! @param [in] r, r(1:3) length N-H, r(4:6) angle H-N-H, in radian
   !-------------------------------------------------------------------
   subroutine nh3_inter_to_cart(x,r)
      implicit none
      real(8), dimension(6)  :: r
      real(8), dimension(12) :: x
 
      real(8) :: r1, r2, r3, a1, a2, a3
      real(8) :: rH1H2, rH3H1, rH1E, rEH2, rEF, rH3F, rNE, rOE, rON
      real(8) :: rrc12, rrc31, rrc23, cost, sint
      real(8) :: cosH3H1H2, sinH3H1H2, cosNH1H2, sinNH1H2, cosNH1H3
 
      r1=r(1)
      r2=r(2)
      r3=r(3)
      a1=r(4)
      a2=r(5)
      a3=r(6)
 
      ! temp values for convenience
      rrc12 = r1*r2*cos(a1)
      rrc23 = r2*r3*cos(a2)
      rrc31 = r3*r1*cos(a3)
 
      rH1H2 = sqrt(r1*r1 + r2*r2 - 2d0*rrc12)
      rH3H1 = sqrt(r3*r3 + r1*r1 - 2d0*rrc31)
 
      ! project N to line H1H2 as point E
      rH1E = (r1*r1-rrc12)/rH1H2
      rEH2 = rH1H2 - rH1E 
      cosNH1H2=rH1E/r1
      sinNH1H2=sqrt(1d0-cosNH1H2*cosNH1H2)
      rNE=r1*sinNH1H2
 
      ! project H3 to line H1H2 as point F
      rEF = (rrc23-rrc31)/rH1H2
      cosH3H1H2 = (rH1E+rEF)/rH3H1
      sinH3H1H2 = sqrt(1d0-cosH3H1H2*cosH3H1H2)
      rH3F = rH3H1*sinH3H1H2
 
      ! get the dihedral angle between plane NH1H2 and H3H1H2
      cosNH1H3=(r1*r1-rrc31)/(rH3H1*r1)
      cost=(cosNH1H3 - cosH3H1H2*cosNH1H2)/(sinH3H1H2*sinNH1H2)
 
      ! project N to plane H1H2H3 as point O
      sint=sqrt(1d0-cost*cost)
      rOE=rNE*cost
      rON=rNE*sint
 
      ! set O as origin, H1H2 parallel to the x axis
      ! x,y,z of N -- atom 1
      x(1) = 0d0
      x(2) = 0d0
      x(3) = rON
 
      ! x,y,z of H1 -- atom 2
      x(4) = -rH1E
      x(5) = -rOE
      x(6) = 0
 
      ! x,y,z of H2 -- atom 3
      x(7) = rEH2
      x(8) = -rOE
      x(9) = 0
 
      ! x,y,z of H3 -- atom 4
      x(10) = rEF
      x(11) = rH3F - rOE
      x(12) = 0
 
   end subroutine nh3_inter_to_cart
end module nh3_mod
!
! =================================================================
!
 
! EDIT BY WYZ IN 2022/04/06
subroutine nh3pot(v, dvdx, ddvdxdx, cart, flag)
   USE nh3_mod, ONLY: nh3_pot_cart
   implicit none
   real(8), intent(out) :: v, dvdx(12), ddvdxdx(12,12)
   real(8), intent(inout)  :: cart(12)
   integer, intent(in)  :: flag
 
   CALL nh3_pot_cart(v, dvdx, ddvdxdx, cart, flag)
end subroutine