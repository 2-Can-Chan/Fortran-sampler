Program Parallel_solve_PDE
Use MPI
Implicit none
	!Parameters
	integer:: Xmin= -40
	integer:: Xmax= 40
	real:: U= 1.5
	integer:: Pts= 6400	!!total amount of grid points

	!variables
	!MPI specifics
	integer:: my_id, numprocs, rc, ierr, ID, tag, stat, tag_L2R, tag_R2L
	parameter (tag_L2R = 10)
	parameter (tag_R2L = 20)
	integer:: nextID, prevID
	integer:: part	!number of pts in each processor
	integer:: istart, iend	!start and end index
	real:: inmsgE, inmsgS	!buff for incoming msg, E to be put at end, S to be put at strt
	integer req, req1, req2
	integer, dimension(MPI_STATUS_SIZE):: stats
	real,dimension(:),allocatable:: rbuf, Phirbuf	!receive buffs for all_gatherv
	real, dimension(:,:),allocatable:: Rslt_buff
	integer, dimension(:), allocatable:: rcounts, displs	!arrays for all_gatherv
	real, dimension(:), allocatable:: L0_buff		!norm error result buff
	!=======

	integer:: bdry_case, i, scheme, comp_asoln, w_debug, k, w_numsoln
	integer:: Szp		! Size for Phi mtx
	real:: t_end, CFL, dt, dx, t_now, dt_cfl, timer_start, timer_end
	real:: an_L0, p_L0, an_L1, p_L1, an_L2, p_L2		!values for norm errors
	real:: an_L1_tot, p_L1_tot, an_L2_tot, p_L2_tot		!buff for norm error ans
	real, dimension(:), allocatable:: Xp
	real, dimension(:), allocatable:: ASoln		!Analytical soln
	real, dimension(:,:), allocatable:: Phi
	real, dimension(:), allocatable:: Phi_error_buff	!vector holder for error calc
	!!Detailed time recorders
	integer:: itns, itn_ctr, r, c, w_elapsetime
	double precision :: itn_tstart, comm_end, calc_start, calc_end, itn_tend
	double precision :: comm_total, calc_total, itn_ttotal
	real, dimension(:), allocatable:: comm_time_lst		!vector holder for comm time
	real, dimension(:), allocatable:: calc_time_lst		!vector holder for calc time
	real, dimension(:), allocatable:: itn_time_lst		!vector holder for comm time
	real, dimension(:,:), allocatable:: itn_tsummary

	!!name holders
	character(25)::rslt_filename
	character(27)::Oput_filename
	character(10)::Oput_time
	character(6):: sch_name		!scheme name for printing purpose
	common /mpistuff/ numprocs, nextID, prevID
	common /error_buff/ an_L0, p_L0, an_L1, p_L1, an_L2, p_L2


	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

	bdry_case= 0		!boundary case 0 and 1
	scheme= 5			!schemes: 1=FTCS, 2=Upwind, 3=LaxFried, 4=Lax-Wen, 5=MacCor
	t_end= 5.			!time
	CFL= 1.
	w_debug = 0			!write interim files for debug
	comp_asoln = 0			!compare analytical solution
	w_numsoln = 0		!write numerical soln.
	w_elapsetime = 1 !write elapsed time file
	inmsgE = 0.0D0; inmsgS = 0.0D0
	itn_ctr = 0


	!Parameter calcs
	dx= (abs(Xmax) + abs(Xmin))/real(Pts-1)
	dt_cfl = CFL*dx/U

	!get predicted iterations to initiate time vector
	if (mod(t_end, dt_cfl) /= 0) then
		itns = floor(t_end / dt_cfl) + 1
	else
		itns = t_end / dt_cfl
	end if

	!define processor start and end index
	part = Pts / numprocs
	if (my_id == 0) then
		istart = 0
		iend = part * (my_id + 1)

	else if (my_id == numprocs-1) then
		istart = (part * my_id) + 1
		iend = Pts - 1

	else
		istart = (part * my_id) + 1
		iend = part * (my_id + 1)
	end if

	do scheme = 1, 5
		if (scheme == 5) then
			Szp=0
		else
			Szp=1
		end if

		!Parallel allocate memory
		!allocate(Phi(istart : iend, Szp:2)); Phi= 0.0D0	!1 is cur, 2 is nxt, 0 is previous
		itn_ctr = 0
		comm_total = 0.0D0
		calc_total = 0.0D0
		itn_ttotal = 0.0D0
		allocate(Xp(istart : iend)); Xp= 0.0D0	!With this method, array idx are istart to iend
		allocate(rbuf(0:Pts-1))
		allocate(Phirbuf(0:Pts-1))

		if (my_id==0) then
			allocate(Phi(istart-1 : iend+1, Szp:2)); Phi= 0.0D0	!1=current, 2=next, 0=previous
			allocate(rcounts(0:numprocs-1)); rcounts= 0.0D0		!rcounts=(/51, 50, 50, 49/)
			allocate(displs(0:numprocs-1)); displs = 0.0D0		!displs=(/0, 51, 101, 151/)
			allocate(L0_buff(0:numprocs*2-1)); L0_buff = 0.0D0	!format, aL0_n, pL0_n

			timer_start = MPI_WTIME()		!start timer

		else if (my_id == numprocs-1) then
			allocate(Phi(istart-1 : iend+1, Szp:2)); Phi= 0.0D0

			!t_summary mtx, row # = proc #
			!column # = itn #
			if (w_elapsetime) then
				allocate(comm_time_lst(0:numprocs-1)); comm_time_lst= 0.0D0		!record comm time
				allocate(calc_time_lst(0:numprocs-1)); calc_time_lst= 0.0D0		!record calc time
				allocate(itn_time_lst(0:numprocs-1)); itn_time_lst= 0.0D0		!record calc time
				allocate(itn_tsummary(0:2, 0:numprocs-1)); itn_tsummary= 0.0D0

			end if

		else
			allocate(Phi(istart-1 : iend+1, Szp:2)); Phi= 0.0D0
		end if

		if (comp_asoln) then			!allocate memory for analytical soln
			allocate(ASoln(istart:iend)); ASoln= 0.0D0
			allocate(Phi_error_buff(istart:iend)); Phi_error_buff=0.0D0
		end if

		!Master determine recieve buffs array
		if (my_id == 0) then
			do k=1, numprocs-1
				displs(k) = (part*k) + 1
			end do
			do k=0, numprocs-2
				rcounts(k) = displs(k+1) - displs(k)
			end do
			rcounts(numprocs-1) = Pts-displs(numprocs-1)
		end if


		!Initialize boundary and arrays
		Xp = (/(i, i= istart, iend)/)
		Xp = (Xp*(abs(Xmax) + abs(Xmin))/(Pts-1)) + Xmin	!vectorized way to initialize Xp

		select case(bdry_case)
			case(0)
				do i= istart, iend
					Phi(i, 1) = (sign(1., Xp(i)) + 1)/2
				end do

			case(1)
				Phi(istart:iend, 1) = 0.5*exp(-(Xp)**2)	!vectorized speed up

		end select

		!Indv. Processors find analytical solution
		if (comp_asoln) then
			if (bdry_case == 0) then
				ASoln = 0.5*(sign(1., Xp - 1.5*t_end) + 1.)
			else if (bdry_case == 1) then
				ASoln = 0.5* exp(-(Xp - 1.5*t_end)**2)
			end if

			if (w_debug) then
				write(rslt_filename, "(A5, I1, A4)") "MPIan", my_id, ".dat"
				open(4, file=rslt_filename, action="write", status="replace")
					do i= istart, iend
						write(4, *) Xp(i), ASoln(i)
					end do
				close(4)
			end if
		end if


		!Hello Neighbours
		prevID = my_id - 1
		nextID = my_id + 1
		if (my_id==0) then
			prevID = numprocs-1
		else if (my_id == numprocs-1) then
			nextID = 0
		end if


		!PDE soln finding loop
		sch_name = "XXXX"
		t_now=0.
		do
			!exchange boundary values
			!SEND ==>>>>>> (0 snd last pt, inner procs send pts both way, last proc snd frst pt)
			!==>>>>>> RECV ()
			itn_tstart = MPI_WTIME()

			if (my_id /= 0) then	!interact with left neighbour
				call MPI_ISEND(Phi(istart,1),1, MPI_DOUBLE_PRECISION, prevID, tag_R2L, &
					MPI_COMM_WORLD, req, ierr)

				call MPI_RECV(inmsgS, 1, MPI_DOUBLE_PRECISION, prevID, tag_L2R, &
					MPI_COMM_WORLD, stats, ierr)

				call MPI_WAIT(req, stats, ierr)
				Phi(istart-1,1) = inmsgS

				!Debug
				!print 18, t_now, my_id, Phi(istart,1), inmsgS
				!18 format('Tnow:',F6.4,x,'#',I1,x,"send2L...",F6.3," recv-fr-L...",F6.3)
			end if

			if (my_id /= (numprocs-1)) then	!send data to right neighbour
				call MPI_ISEND(Phi(iend, 1), 1, MPI_DOUBLE_PRECISION, nextID, tag_L2R, &
					MPI_COMM_WORLD, req, ierr)

				call MPI_RECV(inmsgE, 1, MPI_DOUBLE_PRECISION, nextID, tag_R2L, &
					MPI_COMM_WORLD, stats, ierr)

				call MPI_WAIT(req, stats, ierr)

				Phi(iend+1,1) = inmsgE

				!Debug
				!print 19, t_now, my_id, Phi(iend,1), inmsgE
				!19 format('Tnow:',F6.4,x,'#',I1,x,"send2R...",F6.3," recv-fr-R...",F6.3)
			end if

			comm_end = MPI_WTIME()


			dt = min(dt_cfl, (t_end - t_now))

			!call subroutine to update Phi with different scheme
			calc_start = MPI_WTIME()
			if (my_id == 0) then
				call scheme_update(1, iend, scheme, Phi)

			else if (my_id == (numprocs-1)) then
				call scheme_update(istart, iend-1, scheme, Phi)

			else
				call scheme_update(istart, iend, scheme, Phi)
			end if

			calc_end = MPI_WTIME()


			!update boundary conditions and update Phi columns
			if (my_id==0) then
				Phi(0,2) = 0.	!Phi(0,2) = 0.
			else if (my_id == (numprocs-1)) then
				if (bdry_case == 0) then
					Phi(iend,2) = 1.
				else if (bdry_case == 1) then
					Phi(iend,2) = 0.
				end if
			end if

			Phi(:,1) = Phi(:,2)

			!update t_now
			t_now = t_now + dt

			if (w_elapsetime) then
				itn_tend = MPI_WTIME()
				itn_ttotal = itn_ttotal + (itn_tend - itn_tstart)
				comm_total = comm_total + (comm_end - itn_tstart)
				calc_total = calc_total + (calc_end - calc_start)

			end if

			!itn_ctr = itn_ctr + 1
			if (t_now >= t_end) then
				exit
			end if

		end do 	!!end numerical calculation do loop



		!!Each Processor write num. soln file
		if (w_debug) then
			write(rslt_filename, "(A5, I1, A4)") "MPIrs", my_id, ".dat"
			open(2, file=rslt_filename, action="write", status="replace")
				do i= istart, iend
					write(2, *) Xp(i), Phi(i,1)
				end do
			close(2)
		end if

		!!Find errors between anlyt soln and calculated soln
		if (comp_asoln) then
			call calc_errors(Phi, ASoln, Phi_error_buff)
			!print*, "proc ", my_id, "A_L1= ", an_L1, "P_L1= ", p_L1

			!gather (normal) opo to coallate L0 norm
			call MPI_GATHER((/an_L0, p_L0/), 2, MPI_DOUBLE_PRECISION, L0_buff, 2, &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

			!reduce (sum) opo for L1 norm
			call MPI_REDUCE(an_L1, an_L1_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							0, MPI_COMM_WORLD, ierr)
			call MPI_REDUCE(p_L1, p_L1_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							0, MPI_COMM_WORLD, ierr)

			!reduce (sum) opo for L2 norm
			call MPI_REDUCE(an_L2, an_L2_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							0, MPI_COMM_WORLD, ierr)
			call MPI_REDUCE(p_L2, p_L2_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							0, MPI_COMM_WORLD)
		end if

		! all_gatherv Phi to CPU 0
		call MPI_GATHERV(Xp(istart:iend), size(Xp), MPI_DOUBLE_PRECISION, rbuf, &
				rcounts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call MPI_GATHERV(Phi(istart:iend,1), size(Xp), MPI_DOUBLE_PRECISION, Phirbuf, &
				rcounts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

		! all_gather to CPU (-1) for detailed time data
		if (w_elapsetime) then
			call MPI_GATHER(comm_total, 1, MPI_DOUBLE_PRECISION, comm_time_lst, &
					1, MPI_DOUBLE_PRECISION, (numprocs-1), MPI_COMM_WORLD, ierr)
			call MPI_GATHER(calc_total, 1, MPI_DOUBLE_PRECISION, calc_time_lst, &
					1, MPI_DOUBLE_PRECISION, (numprocs-1), MPI_COMM_WORLD, ierr)
			call MPI_GATHER(itn_ttotal, 1, MPI_DOUBLE_PRECISION, itn_time_lst, &
					1, MPI_DOUBLE_PRECISION, (numprocs-1), MPI_COMM_WORLD, ierr)
		end if

		! Write output and result files
		if (my_id == 0 .AND. w_numsoln) then
			print*, "Writing combined MPI result file with scheme: ", sch_name
			write(rslt_filename, "(A8, I1, A1, I2, A4)") "MPIcomb-", scheme, "-", numprocs, ".dat"
			open(5, file=rslt_filename, action="write", status="replace")
				do i=0, (size(rbuf)-1)
					write(5, *) rbuf(i), Phirbuf(i)
				end do
			close(5)
			print*, "......Finished writing data file"

			!!write norm error files if calculated analytcal soln.
			if (comp_asoln) then
				open(6, file="NormStat.dat", action="write", status="replace")
					write(6, 99) "Integral L1 norm ([anlyt]) [num]: ", &
									"(", an_L1_tot, ")", p_L1_tot
					write(6, 99) "Integral L2 norm ([anlyt]) [num]: ", &
									"(", an_L2_tot, ")", p_L2_tot
					!Since L0 array values are interleaved with analyt then numerical soln
					!the step below is a compact way of finding max value by striding the array
					write(6, 99) "L0 norm ([anlyt]) [num]: ", "(", MAXVAL(L0_buff(0: &
								SIZE(L0_buff)-1:2)), ")" , MAXVAL(L0_buff(1:SIZE(L0_buff):2))
					99 format (A, A, F11.5, A, F11.5)		!print format for stat file

					write(6, *) NEW_LINE('a'), "L0 from each processor [analyt, num]"
					do i=0, numprocs*2-1, 2
						write(6, 98) L0_buff(i), L0_buff(i+1)
						98 format(F11.5, x,x, F11.5)
					end do
				close(6)
				print*, "......Finished writing norm stat file"
			end if
		end if

		!Cleanup
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		deallocate(Xp)
		deallocate(Phi)
		deallocate(rbuf)
		deallocate(Phirbuf)
		!deallocate(comm_time_lst, calc_time_lst, itn_time_lst)

		if (my_id==0) then
			timer_end = MPI_WTIME()
			call date_and_time(Time=Oput_time)
			write(Oput_filename, "(A7, I1, A1, A10, A4)") "Output-", scheme, "-", Oput_time, ".dat"
			open(7, file=Oput_filename, action="write", status="replace")
				write(7, *) "Num or procs used: ", numprocs
				write(7, *) "Num of points: ", Pts
				write(7, *) "Time elapsed: ", (timer_end - timer_start)
			close(7)
			print*, "... Finished writing output file"

			deallocate(rcounts)
			deallocate(displs)
			deallocate(L0_buff)
		end if

		! write detailed elapse time file
		if (my_id == (numprocs-1) .AND. w_elapsetime) then
			itn_tsummary(0, :) = comm_time_lst
			itn_tsummary(1, :) = calc_time_lst
			itn_tsummary(2, :) = itn_time_lst

			call date_and_time(Time=Oput_time)
			write(Oput_filename, "(A7, I1, A1, I2, A1, A10, A4)") "ExTime-", scheme, "-", numprocs, "-", Oput_time, ".dat"
			open(8, file=Oput_filename, action="write", status="replace")
				write(8, *) "Num or procs used: ", numprocs
				write(8, *) "End time: ", t_end, "BC: ", bdry_case
				write(8, *) "Num of iterations: ", itns
				write(8, *) "Detailed time breakdown per iteration (s): row1 = comm, row2 = calc, row3 = itn"
				!!write(8, (numprocsF8.6)), ((comm_tsummary(i, k), i=0, itns-1), k=0, numprocs-1)
				do r = 0, 2
					do c = 0, numprocs-1
						write(8, "(E11.4, a)", ADVANCE="NO") (itn_tsummary(r, c), achar(9))
					end do
					write(8, *)
				end do

			close(8)
			print*, "... Finished writing execution time summary"

			deallocate(comm_time_lst, calc_time_lst, itn_time_lst, itn_tsummary)
		end if

		call MPI_BARRIER(MPI_COMM_WORLD, ierr)

		if (comp_asoln) then; deallocate(ASoln, Phi_error_buff); end if
	end do  !finish looping trough schemes


	call MPI_FINALIZE(ierr)


!===========================================================
!		Subroutines
!===========================================================
CONTAINS

!!=================== Update Phi with numerical schemes =================
SUBROUTINE scheme_update(S_idx, E_idx, sch_id, Phi)
	integer S_idx, E_idx, sch_id, bdry_case
	real, dimension(istart-1:iend+1, Szp:2), intent(inout):: Phi
	integer numprocs, nextID, prevID
	common /mpistuff/ numprocs, nextID, prevID
	!Don't do the below commented lines. They seem to reset variables and
		!integer U, dt, dx, Pts
		!common /cseparams/ U, dt, dx, Pts, bdry_case, sch_name

	select case (sch_id)
		case (1)	! FTCS 1st order
			Phi(S_idx:E_idx,2) = Phi(S_idx:E_idx,1) - (U*dt/dx)*(Phi(S_idx+1:E_idx+1,1) &
				 - Phi(S_idx-1:E_idx-1,1))	!vect speed up
			sch_name= "FTCS"

		case(2)		! Upwind
			Phi(S_idx:E_idx,2) = Phi(S_idx:E_idx,1) - (U*dt/dx)* &
								(Phi(S_idx:E_idx,1) - Phi(S_idx-1:E_idx-1,1))
			sch_name= "Upwind"

		case(3)		! Lax-Friedrichs
			Phi(S_idx:E_idx,2) = 0.5*(Phi(S_idx+1:E_idx+1,1) + Phi(S_idx-1:E_idx-1,1)) &
				- (U*dt/2/dx)*(Phi(S_idx+1:E_idx+1,1) - Phi(S_idx-1:E_idx-1,1))
			sch_name= "LaxFrd"

		case(4)		!Lax-Wendroff (Hoffmann Eqn. 6-11)
			Phi(S_idx:E_idx,2) = Phi(S_idx:E_idx,1) - (U*dt/2/dx)*(Phi(S_idx+1:E_idx+1,1) &
				- Phi(S_idx-1:E_idx-1,1)) + 0.5*((U*dt/dx)**2)*(Phi(S_idx+1:E_idx+1,1) &
				- 2*Phi(S_idx:E_idx,1) + Phi(S_idx-1:E_idx-1,1))
			sch_name= "LaxWen"

		case(5)		!MaxCorm
			!col.0 is predictor step, col.1 is current step, col.2 is next step
			!For linear problem, MaxCorm == Lax-Wendroff
			! ==>predictor step
			Phi(S_idx:E_idx,0) = Phi(S_idx:E_idx,1) - (U*dt/dx)*(Phi(S_idx+1:E_idx+1,1) &
								- Phi(S_idx:E_idx,1))

			!transfer boundary values from predictor step
			if (my_id /= 0) then
				call MPI_RECV(inmsgS, 1, MPI_DOUBLE_PRECISION, prevID, tag_L2R, &
				MPI_COMM_WORLD, stats, ierr)

				Phi(S_idx-1,0) = inmsgS
			end if

			if (my_id /= (numprocs-1)) then
				call MPI_ISEND(Phi(iend, 0), 1, MPI_DOUBLE_PRECISION, nextID, tag_L2R, &
				MPI_COMM_WORLD, req, ierr)
				call MPI_WAIT(req, stats, ierr)
			end if

			! ==>corrector step
			Phi(S_idx:E_idx,2) = 0.5*((Phi(S_idx:E_idx,1) + Phi(S_idx:E_idx,0)) - &
								(U*dt/dx) *(Phi(S_idx:E_idx,0) - Phi(S_idx-1:E_idx-1,0)))
			sch_name= "MaxCor"

		case default
			print*, "Select a proper scheme you dummy"

	end select

END SUBROUTINE scheme_update


!!=================== Subroutine to find norms / errors =================
SUBROUTINE calc_errors(phi, asoln, phi_err_buff)
	real, dimension(istart-1:iend+1, Szp:2), intent(in):: phi
	real, dimension(istart:iend), intent(in):: asoln
	real, dimension(istart:iend), intent(inout):: phi_err_buff
	real an_L0, p_L0, an_L1, p_L1, an_L2, p_L2
	common /error_buff/ an_L0, p_L0, an_L1, p_L1, an_L2, p_L2

	phi_err_buff = phi(istart:iend,1)

	an_L0 = MAXVAL(ABS(asoln))
	p_L0 = MAXVAL(ABS(phi_err_buff))
	an_L1 = SUM(ABS(asoln)) / SIZE(asoln)
	p_L1 = SUM(ABS(phi_err_buff)) / SIZE(phi_err_buff)
	an_L2 = SQRT(SUM(asoln**2) / SIZE(asoln))
	p_L2 = SQRT(SUM(phi_err_buff**2) / SIZE(phi_err_buff))

END SUBROUTINE calc_errors

End Program Parallel_solve_PDE
