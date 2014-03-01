C Trapezoidal approximation integrator in Fortran

	program trapFort
	
		implicit none

		double precision :: x(30),y(30)
		integer :: i=1
		double precision :: result 

		result=0.0

		open(10,file ='betas-benchmark.out', status='old')
		do i=1, 30, 1
			read(10,20, advance='no') x(i), y(i)
		end do

20		format(1E3.3E3)

C Integrate like a boss

		do i=1, 30, 1
			if (i < 30) then
				result= result+y(i)
			else
				result= result+(y(i)*.5)
			endif
		end do 	

C Tie up loose ends of integration
		result= result + ((2*y(1)-y(2))*.5)
		result= result* x(1) 

		write(*,*) 'Integration is equal to', result

	end program
