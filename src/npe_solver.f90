module internal_functions
	use npe_functions
	implicit none
	
	contains
	
	type(StateVariable) function getInitialConditions()

		getInitialConditions = StateVariable(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)

	end function getInitialConditions

   !--------------------------------------------------------------------
   !> @author 
   !> Rodrigo Souza
   !
   ! DESCRIPTION: 
   !> Runge-Kutta solver 
   !> @brief
   !> Computes (i+1)-eth value from the i-eth.
   !> @param[in] a      
   !> @param[in] b
   !> @param[in] alpha      
   !> @param[in] f_t_y 
   !--------------------------------------------------------------------
	subroutine runge_kutta_4th_order___(t, w, step, f_t_y)
		implicit none

		real :: t, w, step, f_t_y !Subroutine parameters
		real :: h
		real :: k1, k2, k3, k4

		h = step / 2.0 !mid-point

		k1 = h * f_t_y(t, w)
		k2 = h * f_t_y(t + h/2., w + k1/2.)
		k3 = h * f_t_y(t + h/2., w + k2/2.)
		k4 = h * f_t_y(t + h, w + k3)
		
		w = w + (k1 + 2.*k2 + 2.*k3 + k4)/6.

		return
	end subroutine runge_kutta_4th_order___

	subroutine mass_bar_RK	(t, w, step, f_t_y)
		implicit none

		real	:: t, w, step, f_t_y

	end subroutine mass_bar_RK
	
end module internal_functions

!--------------------------------------------------------------------
!> @author 
!> Rodrigo Souza
!
! DESCRIPTION: 
!> Main program
!> @brief
!> Computes the mass-radius relation for the npe EoS.
!--------------------------------------------------------------------
program npe_solver
	use npe_functions
	use internal_functions
	
	implicit none
	
	integer, parameter :: TOTAL_ITEMS = 100

	real	:: t_param = 0
	real	:: INITIAL_T = 0. , FINAL_T = 100.
	real	:: STEP_RUNGE_KUTTA
	integer :: step_index = 0.

	!Array to store the computations, 0-based of course
	type(StateVariable) :: initialConditions(0 : TOTAL_ITEMS - 1)

	STEP_RUNGE_KUTTA = (FINAL_T - INITIAL_T) / TOTAL_ITEMS
	
	initialConditions(0) = getInitialConditions()

	print *, 'STEP_RUNGE_KUTTA:  ', STEP_RUNGE_KUTTA

	do while (step_index <= TOTAL_ITEMS)
	
		print *, '(step_index, t_param) = ', step_index, t_param
		print *, '(initialConditions%V_hat_2prime_bar) = ', initialConditions(step_index)%V_hat_2prime_bar
		
		
		
		step_index = step_index + 1
		t_param = t_param + STEP_RUNGE_KUTTA

	enddo

	stop
end program npe_solver
