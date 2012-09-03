!------------------------------------------------------------------------------
! IAG USP
!------------------------------------------------------------------------------
!
! MODULE: npe_functions
!
!> @author
!> Rodrigo Souza
!
! DESCRIPTION: 
!> Functions related to npe equation of state.
!
!------------------------------------------------------------------------------

module npe_functions
	implicit none

   !--------------------------------------------------------------------
   !> @author 
   !> Rodrigo Souza
   !
   ! DESCRIPTION: 
   !> Type that keeps track of evolving variables 
   !> @brief
   !--------------------------------------------------------------------
	type StateVariable
		real m_bar
		real r_bar
		real V_hat_bar
		real V_hat_prime_bar
		real V_hat_2prime_bar
		real n_p_bar
		real n_e_bar
		real nu_bar
		real nu_prime_bar
		real lambda_bar
		real lambda_prime_bar
		real P_bar
		real P_em_bar
		real E_bar
		real epsilon_bar
	end type

    integer, public, parameter :: TOTAL_UNKNOWNS = 15

    contains


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
	subroutine runge_kutta_4th_order(a, b, alpha, f_t_y)
		implicit none
		
		real :: a, b, alpha, f_t_y
		real :: t, h, w
		real :: k1, k2, k3, k4
		integer :: N, iterations  !< Total iterations

		N = 1000
		
		h = (b-a)/N
		t = a
		w = alpha
		iterations = 1
		
		print *, '(t, w) = ', t, w

		
		do while (iterations <= N)
			k1 = h * f_t_y(t, w)
			k2 = h * f_t_y(t + h/2., w + k1/2.)
			k3 = h * f_t_y(t + h/2., w + k2/2.)
			k4 = h * f_t_y(t + h, w + k3)
			
			w = w + (k1 + 2.*k2 + 2.*k3 + k4)/6.
			t = a + iterations * h
			print *, '(t, w) = ', t, w

			iterations = iterations + 1
		enddo		
		
		return
	end subroutine runge_kutta_4th_order
	
end module npe_functions
