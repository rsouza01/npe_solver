!-----------------------------------------------------------------------
! IAG USP
!-----------------------------------------------------------------------
!
! MODULE: cgs_constants
!
!> @author
!> Rodrigo Souza
!
! DESCRIPTION: 
!> Math constants not defined in the language
!
!-----------------------------------------------------------------------
module cgs_constants
	implicit none

	real, public, parameter :: pi = 355./113.
	real, public, parameter :: SUN_MASS				= 1.989e33; !g
	real, public, parameter :: SUN_RADIUS 				= 6.959e10;  !cm

	real, public, parameter :: LIGHT_SPEED 			= 2.998e10;  !cm s-1

	real, public, parameter :: GRAVITATIONAL_CONSTANT	= 6.674e-8;  !dyn cm2 g-2

	real, public, parameter :: PLANCK_CONSTANT_H 		= 6.626e-27;  !erg s
	real, public, parameter :: PLANCK_CONSTANT_H_BAR 	= 1.054e-27;  !erg s

	real, public, parameter :: BOLTZMANN_CONSTANT 		= 1.381e-16;  !erg K-1

	real, public, parameter :: AMU 					= 1.661e-24;  !g

	real, public, parameter :: PROTON_MASS 			= 1.6726e-24; !g
	real, public, parameter :: NEUTRON_MASS 			= 1.6749e-24; !g
	real, public, parameter :: ELETRON_MASS 			= 9.1096e-28; !g

	real, public, parameter :: ELEMENTARY_CHARGE		= 4.8032e-10; !statcoulomb 

end module cgs_constants
