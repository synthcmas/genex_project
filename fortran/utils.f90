module utils
	implicit none
contains
	function poisson_dist(m)
		implicit none
		integer,intent(in)::m
		double precision::poisson_dist
		integer::i
		double precision::factorial

		factorial=1d0
		do i=1,m
			factorial=factorial*dble(i)
		end do

		poisson_dist=exp(-dble(m))*(dble(m)**dble(m))/factorial
	end function poisson_dist		
end module utils