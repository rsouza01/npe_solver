!Program to calculate the mass (SM) and radius (r) for a
!            EQ
!            Free format is employed, compile with e.g.,
!            f90 -c -free star.f and link with f90 -o prog.exe star.o
!

!            The constants defined 
!            below have the following units
!            [hbarc]=197.327 MeVfm
!            [g] = m^3Kg^-1s^-2 gravitational const 
!                = 6.67259*10^-11 m^3Kg^-1s^-2      
!                  or 6.67259 MeV^-2\hbarc/c^4
!            [sm_solar]=1.989D+30  kg
!            The solar mass is converted into MeV. Then
!            one has the relation
!                   1 kg=1.D+30/1.78266270D0  1 MeV/c^2
!            All taken from the Particle Data Group 
!            (1994)                                 
!            The program runs with dim-less variables
!            P=P/central_energy_density
!            rho=rho/central_energy_density
!            with central_energy_density the central energy-density of the star
!            r=r/r0 with r0=1.D+19 fm= 10km
!            m=m/m_solar

!          Prints out central density in fm^-3, radius in m,
!          neutron star mass in terms of solar masses
!
!     Global constants used to render equations dim-less

      MODULE constants 
!     gravitational constant in fmMeV^-1
         DOUBLE PRECISION, PUBLIC, PARAMETER :: g=6.67259D-45*197.327
!     4*\pi                  
         DOUBLE PRECISION, PUBLIC, PARAMETER :: fourpi=12.56637061 
!     the differential equations of Hartle, excluding the quadropole terms
	   !INTEGER, PUBLIC, PARAMETER  :: i=30
	   !INTEGER, PUBLIC, PARAMETER  :: j=30
         INTEGER, PUBLIC , PARAMETER :: number_differential_eqs= 2
!     number of runge-kutta iterations
         INTEGER, PUBLIC , PARAMETER :: max_rk4_steps=100000
!     Dimensionless step in RK4 procedure
         DOUBLE PRECISION, PUBLIC, PARAMETER :: diff_eq_step = 0.001
!     1 MeV/c^2=1.78266270D-30 kg
         DOUBLE PRECISION, PUBLIC, PARAMETER ::  xkg=1.E+30/1.78266270
!     solar mass in MeV/c^2
         DOUBLE PRECISION, PUBLIC, PARAMETER :: sm_solar=1.989E+30*xkg    
!     Define Bbag in units MeV/fm^3
         DOUBLE PRECISION, PUBLIC, PARAMETER :: Bbag=0.0d0             
         
      END MODULE constants 
	  

!     This module contains the parametrization of the EOS as
!     a polynomial in density. 

      MODULE eos
         USE constants
         DOUBLE PRECISION, PUBLIC :: central_energy_density
         CONTAINS
!
!     dP/dr, TOV equation, dimensionless
!
             FUNCTION tov(r,e_rho,y)
             IMPLICIT NONE
             DOUBLE PRECISION ::  tov, r, e_rho
             DOUBLE PRECISION, DIMENSION(number_differential_eqs), INTENT(IN)  :: y
!     non-relativistic
!             tov=-e_rho*y(2)/(r*r)
!     relativistic
      tov =  -(e_rho+y(1))*(y(2)+(r**3)*y(1))/(r*r-2*r*y(2))
             END FUNCTION tov 
             
      END MODULE eos 

!
!     Main subroutine starts here
!
      subroutine solvetov(number_central_densities,rho,press,smass,radius)
      USE eos
      USE constants 
      IMPLICIT NONE
      INTEGER :: j, i, number_central_densities
      DOUBLE PRECISION, DIMENSION (number_central_densities) :: rho, press
      DOUBLE PRECISION, DIMENSION (number_central_densities) :: smass
      DOUBLE PRECISION, DIMENSION (number_differential_eqs) :: ynew, y, ders
      DOUBLE PRECISION ::  star_radius, pressure, const_2, const_1
      DOUBLE PRECISION, DIMENSION (number_central_densities) :: radius

!     loop over central densities

      DO i=1, number_central_densities

!Energia del campo magnetico


!     central energy density in units MeV/fm^3        
!         Bbag=57 MeV/fm^3
         central_energy_density=rho(i)+Bbag
	               

!     central pressure, in units of MeV/fm^-3

!         pressure=0.
         pressure=press(i)-Bbag

!     dimless const_1 const_2 c1 en km y c2 en MeV

         const_1=1.E-18/SQRT(fourpi*g*central_energy_density)
	  !
         const_2 = fourpi*central_energy_density/SQRT(fourpi*g*central_energy_density)**3
!     Dimensionless pressure at the centre (function of central dens)

         y(1)=pressure/central_energy_density

!     Dimensionless start radius

         star_radius=diff_eq_step/10.

!     Dimensionless mass of the star at the centre 

         y(2)=0.0d0 !(star_radius**3)/3.              

!     Start of Runge-Kutta solution to differential equations
!     number of RK4 steps given by variable j

         j=0                               
         DO WHILE((pressure > 0.0d0).AND.(j <= max_rk4_steps))
            j=j+1
!     new dimless radius 
            star_radius=star_radius+diff_eq_step
!     start values for derivatives
            CALL derivatives(star_radius,y,ders,number_central_densities,rho,press,i)        
!     RK4 procedure
            CALL runge_kutta_4(star_radius,y,ynew,ders,number_central_densities,rho,press,i)
!     new values for y(1) (pressure), y(2) (mass) and y(3) (metric function)
            y=ynew                   
!     new pressure, now in units of MeV/fm^-3
            pressure=y(1)*central_energy_density			
         ENDDO

!     radius in units of m radio en km

         radius(i)=star_radius*const_1

!     mass in units of solar masses                      
         smass(i)=y(2)*const_2/sm_solar

		 write(*,*) i,j
!     smass(i)=star_radius**3/3*const_2/sm_solar
      ENDDO          

      END  subroutine solvetov
!
!        4th-Runge-Kutta solution of coupled equations     
!        See any textbook on numerical methods for details
!
      SUBROUTINE runge_kutta_4(x,y,yout,dydx,number_central_densities,rho,press,i)
      USE constants
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(number_differential_eqs) :: yt, dyt,dym
      DOUBLE PRECISION, DIMENSION(number_differential_eqs),INTENT(IN) :: y, dydx
      DOUBLE PRECISION, DIMENSION(number_differential_eqs),INTENT(OUT) :: yout
      DOUBLE PRECISION :: hh, h6, xh
      DOUBLE PRECISION, INTENT(IN) :: x 
      integer number_central_densities,i
      real(8) rho(number_central_densities),press(number_central_densities)


      hh=diff_eq_step*0.5; h6=diff_eq_step/6. ; xh=x+hh
!     first rk-step
      yt=y+hh*dydx
      CALL derivatives(xh,yt,dyt,number_central_densities,rho,press,i)
!     second rk-step
      yt=y+hh*dyt
      CALL derivatives(xh,yt,dym,number_central_densities,rho,press,i)
!     third rk-step
      yt=y+diff_eq_step*dym;  dym=dyt+dym
      CALL derivatives(x+diff_eq_step,yt,dyt,number_central_densities,rho,press,i)
!     fourth rk-step
      yout=y+h6*(dydx+dyt+2.*dym)

      END SUBROUTINE runge_kutta_4
!
!     Here the expressions for the derivatives are set up 
! 
      SUBROUTINE derivatives(r,y,ders,number_central_densities,rho,press,i)
      USE eos
      USE constants 
      IMPLICIT NONE
      DOUBLE PRECISION ::  e_rho
      DOUBLE PRECISION, INTENT(IN)  :: r
      DOUBLE PRECISION, DIMENSION(number_differential_eqs),INTENT(OUT)  :: ders
      DOUBLE PRECISION, DIMENSION(number_differential_eqs),INTENT(IN)  :: y
      integer number_central_densities,i
      real(8) rho(number_central_densities),press(number_central_densities),pint,Uint
      real(8) myinterp
!	  OPEN (5, FILE = 'dens.dat')
      IF(y(1) > 0.0d0) THEN
!     Here we need to use the equation of state to calculate E as a function of P
!     This is an important point; otherwise E will be constant through the integration      
!      e_rho=3.0d0*y(1)+4.0d0*Bbag/central_energy_density
      pint=central_energy_density*y(1)+Bbag
      Uint=myinterp(number_central_densities,press,rho,pint)
      
      e_rho=(Uint+Bbag)/central_energy_density
!     TOV equation=dp/dr      
	  ders(1)=tov(r,e_rho,y)  
!     dm/dr=r*r*energydensity, gravitational mass as function of r
      ders(2)=(r**2)*e_rho
      ENDIF
!if(i>200 .and. i<210) then     
!write(5,*) r, Uint
!end if 
!13 FORMAT (30(ES14.7,1x))

!DO j=1, 2000, 1

!        WRITE(1,13) (perfilRho(j,s),s=1,30,1)

!END DO

      END SUBROUTINE derivatives
