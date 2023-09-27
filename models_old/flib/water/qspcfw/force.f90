! Writed by Liang
! calculate force for bulk water
! the force field is q-SPC/fw
! input unit:
! pbox,coord,cutoff : A
! energy : kcal/mol
! force : kcal/mol/A
! charge : e


subroutine get_force(natoms, pbox, coord, charge, force, energy, cutoff, ewald_parm) bind(c,name="qspcfw_force")
    implicit none

    integer::ii,jj,nn,ii3,jj3

    integer,intent(in)::natoms
    double precision,intent(in)::pbox(3)
    double precision,intent(in)::coord(3*natoms)
    double precision,intent(in)::charge(natoms)
    double precision,intent(out)::energy
    double precision,intent(out)::force(3*natoms)
    double precision,intent(in)::cutoff
    double precision,intent(in)::ewald_parm

    double precision::xxx(natoms),yyy(natoms),zzz(natoms)
    double precision::ffx(natoms), ffy(natoms), ffz(natoms), fx, fy, fz, fff
    integer::nwaters

    !bond & angle variables
    integer::io,ih1,ih2
    double precision::xoh1, xoh2, yoh1, yoh2, zoh1, zoh2, roh1, roh2 
    double precision::xoh1ri, xoh2ri, yoh1ri, yoh2ri, zoh1ri, zoh2ri 
    double precision::costheta, sintheta, theta

    !vdw & elec variable
    double precision::xab, yab, zab, rrr, rrri, sigmari 
    double precision::xabn, yabn, zabn
    integer::ia, ib, j0
    integer:: num_mx, num_my, num_mz, ndim
    integer:: nx, ny, nz, test, totalnumreal_3d
    double precision::nlx, nly, nlz
    integer,allocatable:: ind_x(:), ind_y(:), ind_z(:)
    integer,allocatable:: indx_real_3d(:), indy_real_3d(:), indz_real_3d(:)
    double precision::elec_rcut,elec_rcutsq
    double precision::hc
    double precision::fcc_x(natoms), fcc_y(natoms), fcc_z(natoms), ecc

    double precision::sum_chargesq, cicj, alphar
    

    double precision,parameter::pi = 3.1415926535898
    !-----------Parameter for qSPC/fw like water-----------------
    double precision,parameter::k_bond = 1059.162d0 ! unit in kcal/mol/A^2
    double precision,parameter::eq_bond = 1.000d0 ! unit in A
    double precision,parameter::k_angle = 75.90d0 ! unit in kcal/mol/rad^2
    double precision,parameter::eq_angle = 112.0d0/180.0d0*pi ! unit in rad
    double precision,parameter::sigmaOO = 3.165719505 ! unit in A
    double precision,parameter::epsilonOO = 0.1554 ! unit in kcal/mol
    double precision,parameter::qO = -0.84 ! unit in e 
    double precision,parameter::qH = 0.42 ! unit in e
    !------------------------------------------------------------
    !-----------Parameter for converge and unit------------------
    double precision,parameter::erfcconv = 4.8d0
    double precision,parameter::expconv=24.0d0
    double precision,parameter::r4pie0 = 332.0637751
    !------------------------------------------------------------

!     if(mod(natoms,3).ne.0) then
!         call error('error, number of atoms should be 3 times')
!     end if

!     if(pbox(1)/cutoff<2.0d0.or.pbox(2)/cutoff<2.0d0.or.pbox(3)/cutoff<2.0d0) then
!         call error('error, cutoff must be less than half of box length')
!     end if

    write(*,*) 'check'

    nwaters = natoms/3

    energy = 0.0d0
    force = 0.0d0
    ffx = 0.0d0
    ffy = 0.0d0
    ffz = 0.0d0

    do ii=1,natoms
        ii3=ii*3
        xxx(ii)=coord(ii3-2)
        yyy(ii)=coord(ii3-1)
        zzz(ii)=coord(ii3)
    end do

    !---------------bond & angle term-----------------------------
    do ii=1,nwaters

        ii3 = ii*3
        io = ii3-2
        ih1 = ii3-1
        ih2 = ii3

        ! O-H1
        xoh1 = xxx(io)-xxx(ih1)
        yoh1 = yyy(io)-yyy(ih1)
        zoh1 = zzz(io)-zzz(ih1)

        call image(pbox,xoh1,yoh1,zoh1)

        roh1 = sqrt(xoh1**2+yoh1**2+zoh1**2)
        xoh1ri = xoh1/roh1
        yoh1ri = yoh1/roh1
        zoh1ri = zoh1/roh1

        energy = energy + 0.5d0*k_bond*(roh1-eq_bond)**2 
        fff = -k_bond*(roh1-eq_bond)
        fx = fff*xoh1ri
        fy = fff*yoh1ri
        fz = fff*zoh1ri

        ffx(io) = ffx(io) + fx
        ffy(io) = ffy(io) + fy
        ffz(io) = ffz(io) + fz
        ffx(ih1) = ffx(ih1) - fx
        ffy(ih1) = ffy(ih1) - fy
        ffz(ih1) = ffz(ih1) - fz

        ! O-H2
        xoh2 = xxx(io)-xxx(ih2)
        yoh2 = yyy(io)-yyy(ih2)
        zoh2 = zzz(io)-zzz(ih2)

        call image(pbox,xoh2,yoh2,zoh2)

        roh2 = sqrt(xoh2**2+yoh2**2+zoh2**2)
        xoh2ri = xoh2/roh2
        yoh2ri = yoh2/roh2
        zoh2ri = zoh2/roh2

        energy = energy + 0.5d0*k_bond*(roh2-eq_bond)**2 
        fff = -k_bond*(roh2-eq_bond)
        fx = fff*xoh2ri
        fy = fff*yoh2ri
        fz = fff*zoh2ri

        ffx(io) = ffx(io) + fx
        ffy(io) = ffy(io) + fy
        ffz(io) = ffz(io) + fz
        ffx(ih2) = ffx(ih2) - fx
        ffy(ih2) = ffy(ih2) - fy
        ffz(ih2) = ffz(ih2) - fz

        !H1-O-H2
        costheta = xoh1ri*xoh2ri+yoh1ri*yoh2ri+zoh1ri*zoh2ri
        theta = acos(costheta)
        sintheta = max(1.0d-8,sqrt(1.0d0-costheta**2))

        energy = energy + 0.5d0*k_angle*(theta-eq_angle)**2
        fff = k_angle*(theta-eq_angle)/sintheta

        fx = -fff*(xoh2ri-xoh1ri*costheta)/roh1
        fy = -fff*(yoh2ri-yoh1ri*costheta)/roh1
        fz = -fff*(zoh2ri-zoh1ri*costheta)/roh1
        ffx(ih1) = ffx(ih1) + fx
        ffy(ih1) = ffy(ih1) + fy
        ffz(ih1) = ffz(ih1) + fz
        ffx(io) = ffx(io) - fx
        ffy(io) = ffy(io) - fy
        ffz(io) = ffz(io) - fz

        fx = -fff*(xoh1ri-xoh2ri*costheta)/roh2
        fy = -fff*(yoh1ri-yoh2ri*costheta)/roh2
        fz = -fff*(zoh1ri-zoh2ri*costheta)/roh2
        ffx(ih2) = ffx(ih2) + fx
        ffy(ih2) = ffy(ih2) + fy
        ffz(ih2) = ffz(ih2) + fz
        ffx(io) = ffx(io) - fx
        ffy(io) = ffy(io) - fy
        ffz(io) = ffz(io) - fz

    end do

!     write(*,*) 'check, finish bond & angle'
!     write(*,*) energy

    energy=0.0d0

    ! -------------------------vdw term -------------------------
    ! not add vdw correction for long range vdw yet
    do ii=1,nwaters

        do jj=ii+1,nwaters

            ia=ii*3-2
            ib=jj*3-2

            xab = xxx(ia)-xxx(ib)
            yab = yyy(ia)-yyy(ib)
            zab = zzz(ia)-zzz(ib)

            call image(pbox,xab,yab,zab)

            rrr = sqrt(xab**2+yab**2+zab**2)
            rrri = 1.0d0/rrr

            if(rrr>cutoff) then
                cycle
            end if

            sigmari = sigmaOO*rrri
            energy = energy + 4*epsilonOO*(sigmari**12-sigmari**6)

            fff = -4*epsilonOO*(12*sigmari**12*rrri-6*sigmari**6*rrri)
            fx = fff*xab*rrri
            fy = fff*yab*rrri
            fz = fff*zab*rrri

            ffx(ia) = ffx(ia) + fx
            ffy(ia) = ffy(ia) + fy
            ffz(ia) = ffz(ia) + fz
            ffx(ib) = ffx(ib) - fx
            ffy(ib) = ffy(ib) - fy
            ffz(ib) = ffz(ib) - fz

        end do
    end do

!     write(*,*) 'check, finish vdw'

    ! -------------- electric real term -------------------------

    num_mx = int((erfcconv/ewald_parm)/pbox(1)+0.5d0)+1
    num_my = int((erfcconv/ewald_parm)/pbox(2)+0.5d0)+1
    num_mz = int((erfcconv/ewald_parm)/pbox(3)+0.5d0)+1

    ndim = 4*num_mx*num_my*num_mz

    allocate(ind_x(1:ndim),ind_y(1:ndim),ind_z(1:ndim),stat = test)
    if(test .ne. 0) then
        call error('allocate ind_x,ind_y,ind_z is wrong')
    endif

    elec_rcut=dble(erfcconv)/ewald_parm

    totalnumreal_3d = 0

    do nx = -num_mx,num_mx
        do ny = -num_my,num_my
            do nz = -num_mz,num_mz

                xab = min(abs(nx)*pbox(1),abs(abs(nx)*pbox(1)-0.5d0*pbox(1)))
                yab = min(abs(ny)*pbox(2),abs(abs(ny)*pbox(2)-0.5d0*pbox(2)))
                zab = min(abs(nz)*pbox(3),abs(abs(nz)*pbox(3)-0.5d0*pbox(3)))

                rrr = sqrt(xab**2+yab**2+zab**2)

                if(rrr .le. elec_rcut) then

                    totalnumreal_3d = totalnumreal_3d+1

                    ind_x(totalnumreal_3d) = nx
                    ind_y(totalnumreal_3d) = ny
                    ind_z(totalnumreal_3d) = nz
                    if ( totalnumreal_3d .gt. 8000000 ) then
                        call error('number of ncount_real exceeds the limit')
                    end if
                endif

            enddo ! ny
        enddo  ! nx
    enddo  ! nz

    allocate(indx_real_3d(1:totalnumreal_3d),indy_real_3d(1:totalnumreal_3d),indz_real_3d(1:totalnumreal_3d),stat = test)
    if(test .ne. 0) then
        write(*,*),"allocate indx_real, indy_real, indz_real is wrong"
        stop
    endif

    indx_real_3d(1:totalnumreal_3d) = ind_x(1:totalnumreal_3d)
    indy_real_3d(1:totalnumreal_3d) = ind_y(1:totalnumreal_3d)
    indz_real_3d(1:totalnumreal_3d) = ind_z(1:totalnumreal_3d)

    totalnumreal_3d= 1
    indz_real_3d(1)=0
    indx_real_3d(1)=0
    indx_real_3d(1)=0

    deallocate(ind_x,ind_y,ind_z)

    do ii=1,natoms

!         j0 = ((ii-1)/3)*3+4
        j0 = ii+1

        do jj=j0,natoms

            xab = xxx(ii)-xxx(jj)
            yab = yyy(ii)-yyy(jj)
            zab = zzz(ii)-zzz(jj)

            do nn = 1,totalnumreal_3d

                nlx = real(indx_real_3d(nn))*pbox(1)
                nly = real(indy_real_3d(nn))*pbox(2)
                nlz = real(indz_real_3d(nn))*pbox(3)

                xabn=xab+nlx
                yabn=yab+nly
                zabn=zab+nlz

                rrr = sqrt(xabn**2+yabn**2+zabn**2)
                rrri = 1.0d0/rrr

                if(rrr>elec_rcut) then
                    cycle
                end if

                alphar=ewald_parm*rrr

                cicj=charge(ii)*charge(jj)

                energy = energy + cicj*erfc(alphar)*rrri*r4pie0

                fff=(2*alphar/sqrt(pi)*exp(-alphar**2)+erfc(alphar))*cicj*rrri**2*r4pie0

                fx = fff*xabn*rrri
                fy = fff*yabn*rrri
                fz = fff*zabn*rrri

                ffx(ii) = ffx(ii) + fx
                ffy(ii) = ffy(ii) + fy
                ffz(ii) = ffz(ii) + fz
                ffx(jj) = ffx(jj) - fx
                ffy(jj) = ffy(jj) - fy
                ffz(jj) = ffz(jj) - fz
            end do
            
        end do

    end do

!     write(*,*) 'check, finish electric real'

    !------------------- electric fourier space term ------------

    hc=sqrt(expconv*dble(4))*ewald_parm

    call ewald_image(natoms, charge, xxx, yyy, zzz, hc, pbox(1), pbox(2), pbox(3), &
        dble(1)/ewald_parm, fcc_x, fcc_y, fcc_z, ecc)

    energy = energy + ecc*r4pie0
    ffx = ffx + fcc_x*r4pie0
    ffy = ffy + fcc_y*r4pie0
    ffz = ffz + fcc_z*r4pie0

!     write(*,*) 'check, finish electric image'

    !------------------- electric self term ---------------------
    sum_chargesq=0.0d0
    do ii=1,natoms
        sum_chargesq = sum_chargesq - charge(ii)**2
    end do

    energy = energy + sum_chargesq*ewald_parm/sqrt(pi)*r4pie0

    !----------------------  end  -------------------------------

!     write(*,*) 'check, finish electric self'
!     write(*,*) energy
!     write(*,*) ffx(1)
!     write(*,*) ffy(1)
!     write(*,*) ffz(1)

    do ii=1,natoms

        ii3 = ii*3
        force(ii3-2)=ffx(ii)
        force(ii3-1)=ffy(ii)
        force(ii3)=ffz(ii)

    end do 
        
end subroutine
