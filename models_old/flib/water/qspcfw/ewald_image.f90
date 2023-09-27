!subroutine for calculating coulombic forces in a periodic system using ewald's method 
!Fourier space term [energy and force]



subroutine ewald_image(cnum, charges, xxx, yyy, zzz, hc, a0, b0, c0, sigma, fcc_x, fcc_y, fcc_z, energy)
      
    implicit none

    real*8, parameter :: PI = dble(4) * atan(dble(1))
    real*8, parameter :: SQRT_PI = sqrt(PI)
    integer, intent(in) :: cnum
    real*8, intent(in) :: xxx(cnum), yyy(cnum), zzz(cnum), charges(cnum)
    real*8, intent(in) ::  hc, a0, b0, c0, sigma
    real*8, intent(out) :: fcc_x(cnum), fcc_y(cnum), fcc_z(cnum), energy

    integer::numtotal,totalnumfourier_3d
    integer,allocatable::indx_fourier_3d(:),indy_fourier_3d(:),indz_fourier_3d(:)

    real*8 :: forces_c(cnum, 3)
    complex*16, allocatable :: exp_xc(:, :), exp_yc(:, :), exp_zc(:, :)
    complex*16 :: tempc1_c(cnum), tempc2_c(cnum), tempc3_c(cnum)
    real*8 :: h0a, h0b, h0c, h0a_square, h0b_square, hc_square, h0c_square , & 
        & temp_forces_c(cnum, 3), temp_energy
    integer :: mx, my, mz, hnum1, hnum2, hnum3, k
    real*8 :: temp0, temp1, temp2, temp3, temp_para


    h0a = 2*PI / a0
    h0b = 2*PI / b0
    h0c = 2*PI / c0
    hnum1 = int(hc / h0a) + 1
    hnum2 = int(hc / h0b) + 1
    hnum3 = int(hc / h0c) + 1

    numtotal = (hnum1+1)*(hnum2+1)*(hnum3+1)
    allocate(indx_fourier_3d(numtotal), indy_fourier_3d(numtotal), indz_fourier_3d(numtotal))

    totalnumfourier_3d=0

    do mx=0,hnum1
        do my=0,hnum2
            do mz=0,hnum3

                totalnumfourier_3d=totalnumfourier_3d+1

                if(totalnumfourier_3d .gt. numtotal) then
                  write(*,*) "num of Fourier is too big, error"
                  write(*,*) "end subroutine obtain_num_fourier_3d"
                  stop
                endif

                indx_fourier_3d(totalnumfourier_3d)=mx
                indy_fourier_3d(totalnumfourier_3d)=my
                indz_fourier_3d(totalnumfourier_3d)=mz

            end do
        end do
    end do

    h0a_square = h0a**2
    h0b_square = h0b**2
    h0c_square = h0c**2
    hc_square = hc**2

    ! 计算指数e^(-hx)等
    allocate(exp_xc(cnum, 0:hnum1), exp_yc(cnum, 0:hnum2), exp_zc(cnum, 0:hnum3))

    exp_xc(:, 0) = dble(1)
    exp_yc(:, 0) = dble(1)
    exp_zc(:, 0) = dble(1)

    tempc1_c = exp(dcmplx(0, h0a * xxx))
    tempc2_c = exp(dcmplx(0, h0b * yyy))
    tempc3_c = exp(dcmplx(0, h0c * zzz))

    do mx = 1, hnum1
        exp_xc(:, mx) = exp_xc(:, mx-1) * tempc1_c
    end do

    do my = 1, hnum2
        exp_yc(:, my) = exp_yc(:, my-1) * tempc2_c
    end do

    do mz = 1, hnum3
        exp_zc(:, mz) = exp_zc(:, mz-1) * tempc3_c
    end do

    temp_energy = 0
    temp_forces_c = 0

    temp_para = - sigma**2 / dble(4)

    do k=1,totalnumfourier_3d

        mx=indx_fourier_3d(k)
        my=indy_fourier_3d(k)
        mz=indz_fourier_3d(k)

        temp1 = h0a_square * dble(mx**2) + h0b_square * dble(my**2)
        if ( temp1 > hc_square ) cycle
        temp0 = exp(temp1 * temp_para) / temp1 / dble(2) ! exp(-k^2/4a^2)/k^2

            if ( my /= 0 .and. mz .eq. 0 ) then
                ! 计算第一卦限mx > 0, my > 0
                call do_rec_sum_zero_ewald(cnum, temp0, temp_energy, temp_forces_c,  &
                    & charges, exp_xc(:, mx), exp_yc(:, my), &
                    & mx, my, h0a, h0b, .true., .true.)
                ! u = delta_u * mz = 0

                ! mx < 0, my < 0
                call do_rec_sum_zero_ewald(cnum, temp0, temp_energy, temp_forces_c, &
                    & charges, exp_xc(:, mx), exp_yc(:, my), &
                    & mx, my, h0a, h0b, .false., .false.)
                ! u = delta_u * mz = 0
            end if

            if ( mx /= 0 .and. mz .eq. 0 ) then
                ! mx > 0, my < 0
                call do_rec_sum_zero_ewald(cnum, temp0, temp_energy, temp_forces_c, &
                    & charges, exp_xc(:, mx), exp_yc(:, my), &
                    & mx, my, h0a, h0b, .true., .false. )
                ! u = delta_u * mz = 0

                ! mx < 0, my > 0
                call do_rec_sum_zero_ewald(cnum, temp0, temp_energy, temp_forces_c, &
                    & charges, exp_xc(:, mx), exp_yc(:, my), &
                    & mx, my, h0a, h0b, .false., .true. )
                ! u = delta_u * mz = 0
            end if

            temp2 = temp1 + h0c_square * dble(mz**2)
            if ( temp2 > hc_square ) cycle
            temp0 = exp(temp2 * temp_para) / temp2

                if ( my /= 0 .and. mz /= 0 ) then
                    ! mx > 0, my > 0
                    call do_rec_sum_nonzero_ewald(cnum, temp0, temp_energy, temp_forces_c, &
                        & charges, exp_xc(:, mx), exp_yc(:, my), exp_zc(:, mz), mx, my, mz,&
                        & h0a, h0b, h0c, .true., .true.)

                    ! mx < 0, my < 0_ewald
                    call do_rec_sum_nonzero_ewald(cnum, temp0, temp_energy, temp_forces_c, &
                        & charges, exp_xc(:, mx), exp_yc(:, my), exp_zc(:, mz), mx, my, mz,&
                        & h0a, h0b, h0c, .false., .false.)

                end if

                if ( mx /= 0 .and. mz /= 0 ) then
                    ! mx > 0, my < 0
                    call do_rec_sum_nonzero_ewald(cnum, temp0, temp_energy, temp_forces_c, &
                        & charges, exp_xc(:, mx), exp_yc(:, my), exp_zc(:, mz), mx, my, mz,&
                        & h0a, h0b, h0c, .true., .false.)

                    ! mx < 0, my > 0
                    call do_rec_sum_nonzero_ewald(cnum, temp0, temp_energy, temp_forces_c, &
                        & charges, exp_xc(:, mx), exp_yc(:, my), exp_zc(:, mz), mx, my, mz,&
                        & h0a, h0b, h0c, .false., .true.)
                end if

                if( mx.eq.0 .and. my.eq.0 .and. mz /= 0 ) then
                    call do_rec_sum_nonzero_ewald(cnum, temp0, temp_energy, temp_forces_c, &
                        & charges, exp_xc(:, mx), exp_yc(:, my), exp_zc(:, mz), mx, my, mz,&
                        & h0a, h0b, h0c, .true., .true.)
                end if

    end do

    temp0 = PI * dble(8) / a0 / b0 / c0
    energy = temp0 * temp_energy / 2
    forces_c = temp0 * temp_forces_c
    fcc_x=forces_c(:,1)
    fcc_y=forces_c(:,2)
    fcc_z=forces_c(:,3)
end subroutine 

subroutine do_rec_sum_zero_ewald(cnum, half_phi_rec, temp_energy, temp_forces_c,  &
        & charges, exp_mxc, exp_myc, mx, my, h0a, h0b, mx_positive, my_positive)
    ! 进行mz=0时的倒空间操作
    implicit none
    ! 电荷数，偶极数，第(mx, my, 0)个波矢
    integer, intent(in) :: cnum, mx, my
    ! exp(-hx*xi), exp(-hy*yi)
    complex*16, intent(in) :: exp_mxc(cnum), exp_myc(cnum)
    ! \mu_xi*hx, \mu_yi*hy, qi, phi_rec/2, x/y方向的波矢最小单元
    real*8, intent(in) :: charges(cnum), half_phi_rec, h0a, h0b
    ! mx/my取正或负
    logical, intent(in) :: mx_positive, my_positive
    ! 没有乘上1/(L_a*L_b)的能量和没乘2//(L_a*L_b)的力
    real*8, intent(inout) ::  temp_energy, temp_forces_c(cnum, 3)
    ! 结构因子，电荷对应S，偶极对应M
    complex*16 :: strc_fac_i_s(cnum), strc_fac_s, strc_fac_i_s_noq(cnum)

    real*8 :: temp_list_c(cnum), temp_hx, temp_hy, temp2

    if ( mx_positive ) then
        if ( my_positive ) then
            strc_fac_i_s_noq = exp_mxc * exp_myc
            strc_fac_i_s = charges * strc_fac_i_s_noq
            strc_fac_s = sum(strc_fac_i_s)
            temp_hx = h0a * mx
            temp_hy = h0b * my
        else
            strc_fac_i_s_noq = exp_mxc * conjg(exp_myc)
            strc_fac_i_s = charges * strc_fac_i_s_noq
            strc_fac_s = sum(strc_fac_i_s)
            temp_hx = h0a * mx
            temp_hy = - h0b * my
        end if
    else
        if ( my_positive ) then
            strc_fac_i_s_noq = conjg(exp_mxc) * exp_myc
            strc_fac_i_s = charges * strc_fac_i_s_noq
            strc_fac_s = sum(strc_fac_i_s)
            temp_hx = - h0a * mx
            temp_hy = h0b * my
        else
            strc_fac_i_s_noq = conjg(exp_mxc) * conjg(exp_myc)
            strc_fac_i_s = charges * strc_fac_i_s_noq
            strc_fac_s = sum(strc_fac_i_s)
            temp_hx = - h0a * mx
            temp_hy = - h0b * my
        end if
    end if

    ! q-q
    temp2 = (dble(strc_fac_s)**2 + aimag(strc_fac_s)**2) * half_phi_rec
    temp_energy = temp_energy + temp2
    temp_list_c = (aimag(strc_fac_i_s)*dble(strc_fac_s) - aimag(strc_fac_s)*dble(strc_fac_i_s)) * half_phi_rec
    temp_forces_c(:, 1) = temp_forces_c(:, 1) + temp_list_c * temp_hx
    temp_forces_c(:, 2) = temp_forces_c(:, 2) + temp_list_c * temp_hy
end subroutine

subroutine do_rec_sum_nonzero_ewald(cnum, phi_rec, temp_energy, temp_forces_c,&
    & charges, exp_mxc, exp_myc, exp_mzc, mx, my, mz,&
    & h0a, h0b, delta_u, mx_positive, my_positive)

    ! 进行mz/=0时的倒空间操作
    implicit none
    ! 电荷数，偶极数，第(mx, my, mz)个波矢
    integer, intent(in) :: cnum, mx, my, mz
    ! exp(i*hx*xi), exp(i*hy*yi), exp(i*u*zi)
    complex*16, intent(in) :: exp_mxc(cnum), exp_myc(cnum), exp_mzc(cnum)
    ! \mu_xi*hx, \mu_yi*hy, \mu_zi*u, qi, phi_rec, x/y/z方向的波矢最小单元
    real*8, intent(in) ::   charges(cnum), phi_rec, h0a, h0b, delta_u
    ! mx/my取正或负
    logical, intent(in) :: mx_positive, my_positive
    ! 没有乘上1/(L_a*L_b)的能量和没乘2//(L_a*L_b)的力
    real*8, intent(inout) ::  temp_energy,  temp_forces_c(cnum, 3)
    ! 结构因子，电荷对应S，偶极对应M
    complex*16 ::  strc_fac_i_s(cnum), strc_fac_s, strc_fac_i_s_noq(cnum)

    real*8 :: temp_list_c(cnum), temp_hx, temp_hy, temp_uz, temp2


    if ( mx_positive ) then
        if ( my_positive ) then
            strc_fac_i_s_noq= exp_mxc * exp_myc * exp_mzc
            strc_fac_i_s = charges * strc_fac_i_s_noq
            strc_fac_s = sum(strc_fac_i_s)
            temp_hx = h0a * mx
            temp_hy = h0b * my
        else
            strc_fac_i_s_noq = exp_mxc * conjg(exp_myc) * exp_mzc
            strc_fac_i_s = charges * strc_fac_i_s_noq
            strc_fac_s = sum(strc_fac_i_s)
            temp_hx = h0a * mx
            temp_hy = - h0b * my
        end if
    else
        if ( my_positive ) then
            strc_fac_i_s_noq = conjg(exp_mxc) * exp_myc * exp_mzc
            strc_fac_i_s = charges * strc_fac_i_s_noq
            strc_fac_s = sum(strc_fac_i_s)
            temp_hx = - h0a * mx
            temp_hy = h0b * my
        else
            strc_fac_i_s_noq = conjg(exp_mxc) * conjg(exp_myc) * exp_mzc
            strc_fac_i_s = charges * strc_fac_i_s_noq
            strc_fac_s = sum(strc_fac_i_s)
            temp_hx = - h0a * mx
            temp_hy = - h0b * my
        end if
    end if

    temp_uz = delta_u * mz

    ! q-q
    temp2 = (dble(strc_fac_s)**2 + aimag(strc_fac_s)**2) * phi_rec
    temp_energy = temp_energy + temp2
    temp_list_c = (aimag(strc_fac_i_s)*dble(strc_fac_s) - aimag(strc_fac_s)*dble(strc_fac_i_s)) * phi_rec
    temp_forces_c(:, 1) = temp_forces_c(:, 1) + temp_list_c * temp_hx
    temp_forces_c(:, 2) = temp_forces_c(:, 2) + temp_list_c * temp_hy
    temp_forces_c(:, 3) = temp_forces_c(:, 3) + temp_list_c * temp_uz
end subroutine