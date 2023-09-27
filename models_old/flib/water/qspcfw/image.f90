! subroutine for calculate distance in PBC
! in all 3 directions


subroutine image(pbox,xxx,yyy,zzz)
	implicit none

	double precision,intent(in)::pbox(3)
	double precision,intent(inout)::xxx,yyy,zzz
    double precision::pboxi(3)

	pboxi=1.d0/pbox

	xxx=xxx-pbox(1)*nint(pboxi(1)*xxx)
    yyy=yyy-pbox(2)*nint(pboxi(2)*yyy)
    zzz=zzz-pbox(3)*nint(pboxi(3)*zzz)

end subroutine