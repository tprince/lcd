	module lcd_mod
	contains
	    subroutine forttime(time) bind(c,name='forttime_')
		use, intrinsic :: iso_c_binding
		real(c_float) :: time
	    end subroutine forttime
	end module lcd_mod

