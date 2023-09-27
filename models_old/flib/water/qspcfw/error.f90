!subroutine for error

subroutine error(message)
	implicit none

	character(len=256),intent(in)::message

	write(*,*) trim(message)
	stop

end subroutine