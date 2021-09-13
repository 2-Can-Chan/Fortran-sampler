
! *************************************************************************************
! Lab for Msc course 2016
! Programmer V.A. Titarev/P. Tsoutsanis/ J. Chan (implemented HLL, HLLC flux)
! *************************************************************************************

! Finite-volume scheme with TVD Runge-Kutta time stepping
! for one-dimensional compressible Euler equations
! Use third order TVD Runge-Kutta in time

! Declaration of variables

 IMPLICIT NONE

 ! Spatial order
 Integer :: SpatialOrder ! 1 or 2

 ! Flux type
 Integer FluxType

 ! Courant number
 Real  CFL

 ! Type of initial condition, number of spatial cells
 Integer  IvType, N

 ! Vector of conservative variables CSV(1:4,-6+n+5), first index is the Runge-Kutta stage
 Real, ALLOCATABLE::  CSV(:,:,:)
 ! Primitive variables, no RK stage, W = (rho,u,P)
 Real, ALLOCATABLE::  PV(:,:,:)
 ! Intercell fluxes
 Real, ALLOCATABLE::  IntercellFlux(:,:,:)


 ! other variables
 Integer  ::  OutFreq = 25, Rkstage=0
 Real     LB, RB, h,Time
 Real, ALLOCATABLE::  CELL_CENTERS(:)
 Real, parameter :: GM=1.4
 Real DT,T_
 Integer IT

 ! ---------------------------------- START THE PROGRAM----------------------


 ! initialise of variables, grid etc
  CALL INITALL

 ! START TIME CYCLE
 IT=0
 do  !
   ! Compute a stable time step
   CALL ComputeTimeStep(dt)
   ! Use third-order TVD Runge-Kutta
   Call ThirdOrderTVD
   ! advance time and time counter
   T_ = T_  + DT
   IT = IT+1

   If ( MOD(IT,OutFreq) ==0) PRINT*,' it= ',it,'  t= ',T_
   If (mod(it,2) .eq. 0)  Call OutputTecplot

   ! check whether we reached the output time
   If ( ABS(T_ - TIME)/TIME .LE. 1D-8) GOTO 101
 enddo

 101 CONTINUE
 Call Output
 Call OutputTecplot
 close(111) ! close the movie file
 print*,' Job finished.'
 Print*,' Number of time steps : ',it
 PRINT*,' The end'


 ! //////// code's subroutines ///////////////

 Contains


 !%%%%%%%%%% initialization of the run %%%%%%%%%%
 Subroutine InitAll
  Integer I
  Real X,U1,U2,U3

 ! read the input file
  Open(1,file='euler.ini')
   Read(1,*) SpatialOrder
   Read(1,*) IvType
   Read(1,*) N
   Read(1,*) CFL
   Read(1,*) FluxType
  Close(1)


 SELECT Case(IVTYPE)
 Case(1)
  LB = 0.
  RB = +1.
  TIME = 0.2
 Case(2)
  LB = -5.
  RB = +5.
 Time = 5.d0
 Case Default
  print*,' Wrong test problem number. Stop the code.'
  stop
 END SELECT

 ! Spatial Cell size
 h = (RB-LB)/N

 ! Vector of conservative variables QC(1:4,-6+n+5), first index is the Runge-Kutta stage
 ALLOCATE(CSV(1:4,3,-7:n+6),PV(4,1:4,-7:n+6),InterCellFlux(1:4,3,-1:n+1))
 ALLOCATE(CELL_CENTERS(1:N))

 ! calculate cell centers
 do I=1,n
   CELL_CENTERS(I) =  LB+I*H - H/2
 enddo

 ! initialise the vector of conserved quantities
  do I=1,N
    X = CELL_CENTERS(I)
    CALL  U0(x,CSV(1,:,i))
  enddo

  ! calculate primitive variables from conservative
  Rkstage = 1
  do I=1,N
   PV(Rkstage,1,i) = CSV(Rkstage,1,I)
   PV(Rkstage,2,I) = CSV(Rkstage,2,I)/CSV(Rkstage,1,I)
   PV(Rkstage,3,I) = (GM-1)*( CSV(Rkstage,3,I) - 0.5*CSV(Rkstage,2,I)*PV(Rkstage,2,I))
   PV(Rkstage,4,I) = sqrt(GM*PV(Rkstage,3,I)/PV(Rkstage,1,I))
  enddo

  ! set flow time to zero
   T_=0.D0

  OPEN(UNIT = 111, FILE = 'movie.dat', STATUS = 'UNKNOWN')
  WRITE(111,*)'TITLE="Solution" '
  WRITE(111,*)'VARIABLES="X" "rho" "u" "p"'
  Call OutputTecplot

 End subroutine


 !%%%%%%%%%%%%%% Set up boundary conditions for given stage of the Runge Kutta marching %%%%%%%%%%

 Subroutine SetBC(Rkstage)
   Integer k,i,Rkstage

   ! set up ghost cells

   Do i=-5,0
    do k=1,3
    CSV(Rkstage,k,i)  = CSV(Rkstage,k,abs(i)+1)
	enddo
    do k=1,4
     PV(Rkstage,k,i)  = PV(Rkstage,k,abs(i)+1)
	enddo
   Enddo

   Do i=1,6
    do k=1,3
     CSV(Rkstage,k,N+i)  = CSV(Rkstage,k,N-1-I)
	enddo
    do k=1,4
     PV(Rkstage,k,N+i)  = PV(Rkstage,k,N-1-I)
	enddo
   Enddo
 End subroutine


 !%%%%%%%%%%% Compute initial data at t=0 for given spatial position 'x' %%%%%%%%%%%%%
 Subroutine U0(X,Q)
  ! U1 = RHO, U2 = RHOU, U3 = E
   Real X,U1,U2,U3,Q(3)
   Real DL,DR,UL,UR,PL,PR,x0
   Real :: pi= 3.141592653589793
   ! IvType = 1 : Sod' Shock Tube  Problem
   ! IvType = 2    Shock - turbulence interaction

  SELECT Case(IVTYPE)
   Case(1)
    DL=1.0 ;   UL=0.0 ;  PL=1.0
    DR=0.125 ;  UR=0.0 ;  PR=0.1

    If (X .LE. 0.5)  THEN
     U1 = DL
     U2 = DL*UL
     U3 = PL/(GM-1) + 0.5*DL*UL**2
    Else
     U1 = DR
     U2 = DR*UR
     U3 = PR/(GM-1) + 0.5*DR*UR**2
    Endif

 Case(2)
  ! Long time shock/turbulence interaction
  ! Mach number 1.1, S=1.5
  DL=   1.51569506726457
  UL=   0.523345519274197
  PL=   1.80500000000000

  IF (X .LE. -4.5)  THEN
   U1 = DL
   U2 = DL*UL
   U3 = PL/(GM-1) + 0.5*DL*UL**2
  ELSE
   U1 = 1 + 0.1d0*sin(20*pi*x)
   U2 = 0.
   U3 = 1./(GM-1)  ! U = 0
  ENDIF

 End Select

 Q(1) = U1
 Q(2) = U2
 Q(3) = U3
end subroutine


!%%%%%%%%%%%%%%%% Time marching algorithm, which uses third order TVD RK method  %%%%%%%
!%%%%%%  Jiang G.S. and Shu C.W. Efficient Implementation of  weighted ENO schemes //J. Comput. Phys. 1996.  V. 126.  pp.202-212.

  Subroutine  ThirdOrderTVD
   Integer i,k

  ! loop stages from 1 to 3
   do Rkstage=1,3
     ! set up boundary conditions
     Call SetBc(Rkstage)
     ! calculate intercell fluxes
     Call ComputeFlux(Rkstage)
     ! perform the update
     CALL Update(Rkstage)
     Do i=1,n
       PV(Rkstage+1,1,i) = CSV(Rkstage+1,1,i)
       PV(Rkstage+1,2,i) = CSV(Rkstage+1,2,i)/CSV(Rkstage+1,1,i)
       PV(Rkstage+1,3,i) = (gm-1)*( CSV(Rkstage+1,3,i) - 0.5*PV(Rkstage+1,1,i)*PV(Rkstage+1,2,i)**2)
       PV(Rkstage+1,4,I) = sqrt(GM*PV(Rkstage+1,3,I)/PV(Rkstage+1,1,I))
     Enddo
   enddo

   ! re-assign the flow variables to stage 1 of RK method
   Do i=1,n
   do k=1,3
     CSV(1,k,i) = CSV(4,k,i)
   enddo
   do k=1,4
     PV(1,k,i) = PV(4,k,i)
   enddo
   Enddo

  End subroutine

 !%%%%%%%%%%% Solution update for each stage of TVD RK method %%%%%%%%%%
  Subroutine UPDATE(Rkstage)
   Integer I,K,Rkstage

   SELECT Case(Rkstage)
    Case(1)
     do i=1,n
      do K=1,3
       CSV(2,K,i)  =  CSV(1,K,i)  - (Dt/H)*(InterCellFlux(1,K,i) - InterCellFlux(1,K,i-1))
      enddo
     enddo

    Case(2)
     do i=1,n
      do K=1,3
       CSV(3,k,i) =  0.75*Csv(1,k,i) + 0.25*CSV(2,k,i) - (0.25*Dt/H)*(InterCellFlux(2,k,i) - InterCellFlux(2,k,i-1))
      enddo
     enddo

    Case(3)
     do i=1,n
      do k=1,3
       CSV(4,K,i) =  (1./3)*CSV(1,K,i) + (2./3)*CSV(3,K,i) - (2./3)*(Dt/H)*(InterCellFlux(3,K,i) - InterCellFlux(3,K,i-1))
	  enddo
     enddo
    END SELECT
  end subroutine


 !%%%%%%%%% write the output file %%%%%%%%%
  Subroutine    Output
  Integer i
	character(len=11)::zonename
	write(zonename, "(A5,I1,A3,I1,A1)") '"Flux',FluxType,'od_', SpatialOrder,'"'

	 202 format(6(2x,e11.4))
   open(1,file='results.dat')
   WRITE(1,*)'TITLE="Solution" '
  WRITE(1,*)'VARIABLES="X" "rho""u""P" "T"' ! "u" "p"
   WRITE(1,*)'ZONE', ',T=',zonename, ',I=',n, ',F="POINT"'
   do i=1,n
    write(1,202) cell_centers(i),CSV(1,1,i),PV(1,2,I),PV(1,3,I),PV(1,3,I)/PV(1,1,I)
   enddo
   close(1)
  End subroutine


  !%%%%%%%%%% Evaluation of the physical flux function from the conserved vector CDS =(rho,rho*u,E)
  Subroutine FLUEVAL(CDS,Flux)
    Real cds(3),p,u,flux(3)

    u = cds(2)/cds(1)
    p = (gM-1)*(Cds(3) - 0.5*cds(1)*u**2)
    Flux(1) = cds(2)
    Flux(2) = cds(2)*u + p
    Flux(3) = (cds(3)+p)*u

   End subroutine


  !%%%%%%%%%% Calculation of a stable time step %%%%%
  Subroutine ComputeTimeStep(dt)
   Integer i
   Real Umax,dt, a

   umax  = 0.0
   Do i=1,n
    ! compute the sound speed
    a = ComputeSoundSpeed(CSV(1,:,i))
    umax = max(umax, a + abs(PV(1,2,i)))
   Enddo

   ! reduce the time step for first 10 time steps
   If ( IT<10) THEN
    dt = MIN(0.1*H/UMAX, TIME-T_)
   Else
    dt = MIN(CFL*H/UMAX, TIME-T_)
   Endif
  End subroutine


 !%%%%%%%%%%%%%%%% Calculation of the sound speed  on the conserved vector CDS
  Real function ComputeSoundSpeed(cds)
    Real cds(3),p,u
    u = cds(2)/cds(1)
    p = (gm-1)*(cds(3) - 0.5*cds(2)*u)
    ComputeSoundSpeed=sqrt(gm*p/cds(1))
  End function


  !%%%%%%%%%%%%%%%%% minmod slope limiter %%%%%%%%%%%%%%%%%
  Real function minmod(x,y)
   Real x,y
   minmod = 0.
   !minmod = 0.5*(sign(x,1.) + sign(y,1.)) * min(abs(x), abs(y))
   !minmod = (max(0.0d0, sign(0.5d0,x)) + max(0.0d0, sign(0.5d0,y))) * min(abs(x), abs(y))
   minmod = 0.5*(sign(1.0d0,x) + sign(1.0d0,y)) * min(abs(x), abs(y))

   !print*,' minmod limiiter is not implemented'
   !stop
  End function

 ! Compute the numerical flux
 Subroutine ComputeFlux(Rkstage)
  Integer i,k,Rkstage
  Real CDL(3),CDR(3),LocalFlux(3)

  ! Loop over the spatial index i
  do I=0,N

   ! call reconstruction procedure at RK stage Rkstage to compute left CDL and right CDR
   ! values of the conserved vector between cells i and i+1
   CALL Reconstruction(CSV(Rkstage,:,i-2:i+3),CDL,CDR)

   ! calculate the numerical flux using reconstructed conserved vectors CDL, CDR
   ! Left   initial data  for the local Riemann problem is given by CDL
   ! Right  initial data  for the local Riemann problem is given by CDR

   Select Case(FluxType)
   Case(1)
      CALL LxF(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(2)
      CALL Rusanov(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(3)
      CALL HLL(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(4)
      CALL HLLC(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case default
    print*,' the flux is not defined. stop the program'
	read*
	stop
  End select

  Enddo
 end subroutine


  !%%%%%%%%% Lax Friedrich flux %%%%%%%%
   Subroutine LxF(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3)

        CALL FLUEVAL(CDL,FL)
        CALL FLUEVAL(CDR,FR)

        Flux = 0.5*(FL+FR) - 0.5*(h/dt)*(CDR-CDL)
  End subroutine


 !%%%%%%%% Rusanov flux %%%%%%%%%%%%%%%%%%
   Subroutine Rusanov(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3),Speed

      print*,'no Rusanov flux found. stop the code!'
      stop
  End subroutine

 !%%%%%%%%%%%%%% HLL flux %%%%%%%%%%%%%%%%%%%%%
   Subroutine HLL(CDL,CDR,Flux)
      Real  CDL(3),CDR(3),Flux(3)
      Real  FL(3), FR(3)
      Real  S_L, S_R, SOS_L, SOS_R

      !Speed of sound at left and right state
      SOS_L = ComputeSoundSpeed(CDL)
      SOS_R = ComputeSoundSpeed(CDR)

      S_L = MIN((CDL(2)/CDL(1) - SOS_L), (CDR(2)/CDR(1) - SOS_R))
      S_R = MAX((CDL(2)/CDL(1) + SOS_L), (CDR(2)/CDR(1) + SOS_R))

      !Get Fluxes
      Call FLUEVAL(CDL, FL)
      Call FLUEVAL(CDR, FR)

      if (S_L >= 0) then
        Flux = FL
      else if (S_L <= 0 .and. S_R >= 0) then
        Flux = (S_R*FL - S_L*FR + S_L*S_R*(CDR - CDL)) / (S_R - S_L)
      else if (S_R >= 0) then
        Flux = FR
      end if

  End subroutine


 !%%%%%%%%%%%% HLLC flux %%%%%%%%%%%%%%%%%%%%%
   Subroutine HLLC(CDL,CDR,Flux)
      Real  CDL(3),CDR(3),Flux(3)
      Real  FL(3), FR(3), CDL_star(3), CDR_star(3)
      Real  S_L, S_R, S_star, a_L, a_R
      Real  u_L, u_R, p_L, p_R, rho_L, rho_R

      u_R = CDR(2) / CDR(1)
      p_R = (gM-1)*(CDR(3) - 0.5*CDR(1)*u_R**2)
      rho_R = CDR(1)
      u_L = CDL(2) / CDL(1)
      p_L = (gM-1)*(CDL(3) - 0.5*CDL(1)*u_L**2)
      rho_L = CDL(1)

      !Get Fluxes
      Call FLUEVAL(CDL, FL)
      Call FLUEVAL(CDR, FR)

      a_L = ComputeSoundSpeed(CDL)
      a_R = ComputeSoundSpeed(CDR)
      S_L = MIN((u_L - a_L), (u_R - a_R))
      S_R = MAX((u_L + a_L), (u_R + a_R))

      S_star = (p_R - p_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R)) / &
                (rho_L*(S_L - u_L) - rho_R*(S_R - u_R))

      ! Get conserved variable at star states
      CDL_star = HLLC_U_star(CDL, S_L, S_star)
      CDR_star = HLLC_U_star(CDR, S_R, S_star)

      if (S_L >= 0) then
        Flux = FL
      else if (S_L <= 0 .and. S_star >= 0) then
        Flux = FL + S_L*(CDL_star - CDL)
      else if (S_star <= 0 .and. S_R >= 0) then
        Flux = FR + S_R*(CDR_star - CDR)
      else
        Flux = FR
      end if

  End subroutine

  ! function to find U* of HLLC. Requires vector of conserved qty and
  ! speed of sound L / R (denote S_K here) and S_star
  ! returns vector of U*
  Function HLLC_U_star(CDS, S_K, S_star)
    Real  CDS(3), S_K, S_star, p_k, u_k
    Real, dimension(3) :: HLLC_U_star     !! switch type

    u_k = CDS(2) / CDS(1)
    p_k = (gM-1)*(CDS(3) - 0.5*CDS(1)*u_K**2)

    ! build interim U* vector
    HLLC_U_star(1) = 1
    HLLC_U_star(2) = S_star
    HLLC_U_star(3) = CDS(3)/CDS(1) + (S_star - u_k)*(S_star + &
                      (p_k/CDS(1)/(S_K - u_k)))

    HLLC_U_star = CDS(1)*(S_K - u_K)/(S_K - S_star) * HLLC_U_star

  End function



 !%%%%%%%%%%%%%%%% Reconstruction procedure%%%%%%%%%%%%%%
 ! Input: one-dimensional array U1D of flow quantities near cell interface i+1/2
 ! Output: left CDL and right CDR values at interface
 Subroutine Reconstruction(U1D,CDL,CDR)
   Integer I, F
   Real U1d(3,-2:3),CDL(3),CDR(3), alpha, beta
   Real r

   select Case(SpatialOrder)

	 Case(1)

	 ! First order
	  Do f=1,3
	   CDL(f) = U1D(f,0)
	   CDR(f) = U1D(f,1)
	  Enddo

	 ! second order TVD
	 Case(2)
     ! set alpha and beta as backward and forward differencing
	 	do f=1,3
		   ! for CDL
		   alpha = U1D(f,0) - U1D(f,-1)   !backward differencing
		   beta = U1D(f,1) - U1D(f,0)   !forward differencing
		   CDL(f) = U1D(f,0) + 0.5*minmod(alpha, beta)

		   ! for CDR
		   alpha = U1D(f,1) - U1D(f,0)
		   beta = U1D(f,2) - U1D(f,1)
		   CDR(f) = U1D(f,1) - 0.5*minmod(alpha, beta)
		end do

	 Case default
	 print*,' Wrong spatial accuracy. Stop the code'
	 stop
  end select

 end subroutine

 Subroutine OutPutTecplot
    Integer i,j
		real(8) x

   55 Format (4(2x,e11.4))
    WRITE(111,*)'ZONE','I=',n, ',F="POINT"'
    WRITE(111,*) ', SOLUTIONTIME=',T_
    DO  I = 1,n
	  x = lb + i*h-h/2
      WRITE(111,55)X,pv(1,1,i),pv(1,2,i),pv(1,3,i)
    Enddo

 End subroutine


 END
