!--------------------------------------------------------------------------------------------------------
MODULE module_Nbody_AffineTides_solver
!--------------------------------------------------------------------------------------------------------
IMPLICIT NONE
integer												:: n_particles			
integer												:: length_Y_per_body_n = 24 ! = 3+3+(3x3)+(3x3)
integer												:: tot_nr_Y_evol_eqs		! = n_particles*length_Y_per_body_n
real*8, dimension(:), allocatable					:: mass, radius, gas_n
real*8, dimension(:,:), allocatable					:: IC_par_pos, IC_par_vel
real*8, dimension(:), allocatable					:: gas_gamma, Mqp_sph
integer, dimension(:), allocatable					:: evoTides_yesno, RigidSph_yesno
real*8, dimension(3,3)								:: I_MAT
real*8												:: PN_gamma = 2.12370597232e-06	!postnewtonian gamma - units: Rsun, Msun, ....
real*8                                              :: rsun_to_au = 0.004650467 !converts Rsun to AU, for writing output in AU
real*8                                              :: vCM_to_kms = 437.892 !converts CoM velocity to km/s by multipliying by sqrt(GMsun/Rsun)
integer												:: use_12PN, use_1PN, use_2PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps
integer												:: outputinfo_screenfiles, downsample, stepc
real*8												:: scale_dt, max_sim_time, evolvetides_threshold, ENDbinsingle_threshold
real*8												:: max_simtime_sec, IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold
integer												:: pass_tidalthreshold_yesno
!External programs:
!EXTERNAL 											:: c_wrapper_gsl_sf_ellint_rj
!EXTERNAL         									:: DGETRF
!EXTERNAL											:: DGETRI
!EXTERNAL											:: DSYEV
!LAPACK: DGETRF, DGETRI
integer												:: LP_INFO
integer, DIMENSION(3)								:: LP_IPIV
real*8, DIMENSION(3,3)								:: LP_A
real*8, dimension(3)								:: LP_WORK		
real*8, dimension(3)								:: LP_W
real*8, dimension(8)								:: LP_EV_WORK
!counters:
integer												:: i,j,k,l,is,sc,nc,ds


CONTAINS
	
		
	!----------------------------------------------------------------------------------------------------
	FUNCTION len3vec(V)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(3), intent(in)	:: V
		real*8								:: len3vec	
		!---------------------------------------
		!calc length of 3-vec:
		len3vec = SQRT(SUM(V**2))
		!---------------------------------------	
	RETURN
	end function len3vec
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------
	
	
	!----------------------------------------------------------------------------------------------------
	FUNCTION cross(a,b)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(3), intent(in)    :: a, b
		real*8, dimension(3)                :: cross
		!---------------------------------------
		!calc cross product of two 3-vecs:
        cross(1) = a(2)*b(3)-a(3)*b(2)
        cross(2) = a(3)*b(1)-a(1)*b(3)
        cross(3) = a(1)*b(2)-a(2)*b(1)
		!---------------------------------------	
	RETURN
	end function cross
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------



	!----------------------------------------------------------------------------------------------------
	FUNCTION CoM_2body(vec_1, mass_1, vec_2, mass_2)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(3), intent(in)	:: vec_1, vec_2
		real*8, intent(in)	                :: mass_1, mass_2
		real*8, dimension(3)				:: CoM_2body
		!---------------------------------------
		!calc center of mass
        CoM_2body = (mass_1*vec_1 + mass_2*vec_2)/(mass_1 + mass_2)
		!---------------------------------------	
	RETURN
	end function CoM_2body
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------

	
	!----------------------------------------------------------------------------------------------------
	FUNCTION CoM_3body(vec_1, mass_1, vec_2, mass_2, vec_3, mass_3)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(3), intent(in)	:: vec_1, vec_2, vec_3
		real*8, intent(in)	                :: mass_1, mass_2, mass_3
		real*8, dimension(3)				:: CoM_3body
		!---------------------------------------
		!calc center of mass
        CoM_3body = (mass_1*vec_1 + mass_2*vec_2 + mass_3*vec_3)/(mass_1 + mass_2 + mass_3)
		!---------------------------------------	
	RETURN
	end function CoM_3body
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------


	!----------------------------------------------------------------------------------------------------
	FUNCTION CoM_4body(vec_1, mass_1, vec_2, mass_2, vec_3, mass_3, vec_4, mass_4)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(3), intent(in)	:: vec_1, vec_2, vec_3, vec_4
		real*8, intent(in)	                :: mass_1, mass_2, mass_3, mass_4
		real*8, dimension(3)				:: CoM_4body
		!---------------------------------------
		!calc center of mass
        CoM_4body = (mass_1*vec_1 + mass_2*vec_2 + mass_3*vec_3 + mass_4*vec_4) &
                            /(mass_1 + mass_2 + mass_3 + mass_4)
		!---------------------------------------	
	RETURN
	end function CoM_4body
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------


	!----------------------------------------------------------------------------------------------------
	FUNCTION Fthresh(m11, m12, m2, a1, e1, r_12)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8,                  intent(in)	:: m11, m12, m2, a1, e1, r_12
		real*8               				:: Fthresh
		!---------------------------------------
		!calculate ratio between bounding force and tidal force from perturber
        Fthresh = 2d0*(m11+m12)*m2*(a1*(1+e1))**3d0 / (m11*m12*r_12**(3d0))
		!---------------------------------------	
	RETURN
	end function Fthresh
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------


	!----------------------------------------------------------------------------------------------------
	FUNCTION func_det_3x3M(M)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(3,3), intent(in)	:: M
		real*8								:: func_det_3x3M	
		!---------------------------------------
		!calc det of 3x3 Matrix (M):
		func_det_3x3M =	M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) - &
						M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) + &
						M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
		!---------------------------------------	
	RETURN
	end function func_det_3x3M
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------


	!----------------------------------------------------------------------------------------------------
	FUNCTION Trace_3x3M(M)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(3,3), intent(in)	:: M
		real*8								:: Trace_3x3M
		!---------------------------------------
		!calc Trace of 3x3 Matrix (M):
		Trace_3x3M = M(1,1) + M(2,2) + M(3,3)
		!---------------------------------------	
	RETURN
	end function Trace_3x3M
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------


	!----------------------------------------------------------------------------------------------------
	FUNCTION Nbody_dt(Y)	!calc optimal (adaptive) time step
	!----------------------------------------------------------------------------------------------------	
		IMPLICIT NONE
		real*8, dimension(tot_nr_Y_evol_eqs), intent(in)	:: Y
		real*8												:: Nbody_dt
		real*8, dimension(n_particles,length_Y_per_body_n)	:: Y_in_body_all_all
		real*8, dimension(n_particles,3)					:: pos, vel
		real*8, dimension(n_particles,n_particles)			:: r_ij, v_ij, a_ij
		real*8												:: min_deltaT_rov, min_deltaT_roa, min_deltaT	
		!---------------------------------------
		!UNPACK:
		Y_in_body_all_all(:,:) 		= RESHAPE(Y, (/ n_particles,length_Y_per_body_n /), ORDER = (/2,1/))
		pos(:,:)					= Y_in_body_all_all(:,1:3)
		vel(:,:)					= Y_in_body_all_all(:,4:6)
		!initilize:
		r_ij(:,:) = 1e10	!=high number to make use MIN(r/..) is not i,i.
		v_ij(:,:) = 1e-10	!=1 to avoid divergence at i,i when calc dt.
		a_ij(:,:) = 1e-10 !=1 to avoid divergence at i,i when calc dt.
		!i,j loop:
		do i=1, n_particles,1
		do j=1, n_particles,1
		if (i .NE. j) then
			!rel pos, vel, acc between i,j.
			r_ij(i,j) 	= len3vec(pos(i,:)-pos(j,:))						!dist between i,j
			v_ij(i,j) 	= len3vec(vel(i,:)-vel(j,:))						!vel betweeen i,j
			a_ij(i,j)	= mass(j)/(r_ij(i,j)**2)							!acc of i caused by j.
		endif
		enddo
		enddo
		!min delta_t for NBsystem:
		min_deltaT_rov	= (MINVAL(r_ij(:,:)/v_ij(:,:)))				!(r/v)		units: time
		min_deltaT_roa	= (MINVAL(r_ij(:,:)/a_ij(:,:)))**(1./2.)	!sqrt(r/a)	units: time
		min_deltaT		= (MINVAL([min_deltaT_rov, min_deltaT_roa]))
		!return:
		Nbody_dt = min_deltaT
		!---------------------------------------	
	end function Nbody_dt
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------

	
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE Calc_Aint_MAT_selfgrav(S_in, Aint_out)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		!in, out:
		real*8, dimension(3,3), intent(in)	:: S_in
		real*8, dimension(3,3), intent(out)	:: Aint_out
		!params in subroutine:
		real*8, dimension(3)				:: Eigen_vals
		real*8, dimension(3,3)				:: Eigen_vecs, Aint_MAT_D, R_MAT, R_MAT_inv
		real*8								:: s1,s2,s3, si
		real*8								:: CSEI_i
		real*8								:: err_ellint_rj = 0.0001
		real*8								:: Aint_Di	
		!---------------------------------------
		!initialize:
		Aint_MAT_D(:,:) = 0.0
		!calc Eigen vals, vecs:
		LP_A = S_in
		CALL DSYEV( 'V', 'U', 3, LP_A, 3, LP_W, LP_EV_WORK, 8, LP_INFO)
		Eigen_vals	= LP_W
		Eigen_vecs	= LP_A 
		R_MAT = Eigen_vecs
		!calc R^-1:
		CALL Calc_Matrix_Inverse(R_MAT,R_MAT_inv)
		!calc ell integral:
	    s1 = Eigen_vals(1)
	    s2 = Eigen_vals(2)
	    s3 = Eigen_vals(3)
	    do sc=1,3,1
	        si = Eigen_vals(sc)
			CALL c_wrapper_gsl_sf_ellint_rj(CSEI_i, s1, s2, s3, si, err_ellint_rj)
	        Aint_Di = (2d0/3d0)*CSEI_i !(S_i_absdet**(1d0/2d0))*((2d0/3d0)*CSEI_i)
	        Aint_MAT_D(sc,sc) = Aint_Di			
		enddo	
	    !form final Aint:
		Aint_out = MATMUL(MATMUL(R_MAT, Aint_MAT_D), R_MAT_inv)
		!---------------------------------------
	END SUBROUTINE Calc_Aint_MAT_selfgrav
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------	
		
		
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE Calc_Matrix_Inverse(M_in, M_inv_out)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		!in, out:
		real*8, dimension(3,3), intent(in)	:: M_in
		real*8, dimension(3,3), intent(out)	:: M_inv_out
		!---------------------------------------
		!calc:
		LP_A		= M_in 			!LP_A is a dummy for in/out
		CALL DGETRF(3, 3, LP_A, 3, LP_IPIV, LP_INFO)
		CALL DGETRI(3, LP_A, 3, LP_IPIV, LP_WORK, 3, LP_INFO)
		!Final inverse:
		M_inv_out	= LP_A			!LP_A is a dummy for in/out
		!---------------------------------------	
	END SUBROUTINE Calc_Matrix_Inverse
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------	


	!----------------------------------------------------------------------------------------------------
	SUBROUTINE Calc_binary_L(pos_1, vel_1, mass_1, pos_2, vel_2, mass_2, Lvec)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		!in, out:
		real*8, dimension(3), 	intent(in)	:: pos_1, vel_1
		real*8, dimension(3), 	intent(in)	:: pos_2, vel_2
		real*8, 				intent(in)	:: mass_1, mass_2
		real*8, dimension(3),	intent(out)	:: Lvec
		!params in subroutine:
		real*8								:: Mtot, mred
		real*8, dimension(3)				:: pos, vel
		!---------------------------------------
		!define/calc:
		pos			= pos_1 - pos_2
		vel			= vel_1 - vel_2 		
		Mtot		= mass_1+mass_2
		mred		= mass_1*mass_2/Mtot
        Lvec        = mred*cross(pos,vel)
		!---------------------------------------	
	END SUBROUTINE Calc_binary_L
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------


	!----------------------------------------------------------------------------------------------------
	SUBROUTINE Calc_binary_info(pos_1, vel_1, mass_1, pos_2, vel_2, mass_2, out_binary_info_arr)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		!in, out:
		real*8, dimension(3), 	intent(in)	:: pos_1, vel_1
		real*8, dimension(3), 	intent(in)	:: pos_2, vel_2
		real*8, 				intent(in)	:: mass_1, mass_2
		real*8, dimension(5),	intent(out)	:: out_binary_info_arr
		!params in subroutine:
		real*8								:: Mtot, mred
		real*8, dimension(3)				:: pos, vel, Lvec
		real*8								:: len_L, len_r, len_v
		real*8								:: E_kin, E_pot, E_tot, a_bin, e_bin
		!---------------------------------------
		!define/calc:
		pos			= pos_1 - pos_2
		vel			= vel_1 - vel_2 		
		Mtot		= mass_1+mass_2
		mred		= mass_1*mass_2/Mtot
		!Lvec		= mred*([	pos(2)*vel(3)-pos(3)*vel(2), &
		!						pos(3)*vel(1)-pos(1)*vel(3), &
		!						pos(1)*vel(2)-pos(2)*vel(1)])						
        Lvec        = mred*cross(pos,vel)
		len_L		= len3vec(Lvec)			
		len_r 		= len3vec(pos)
		len_v		= len3vec(vel)		
		!calc orbital params wrt bin CM:
		E_kin	= (1./2.)*mred*(len_v**2.)
		E_pot	= - Mtot*mred/len_r
		E_tot	= E_kin + E_pot
		a_bin	= - Mtot*mred/(2.*E_tot)
		e_bin	= SQRT(1. + (2.*E_tot*(len_L**2.))/(mred*((Mtot*mred)**2.)))
		!return:
		out_binary_info_arr(:) = [E_kin, E_pot, E_tot, a_bin, e_bin]
		!---------------------------------------	
	END SUBROUTINE Calc_binary_info
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------


	

	!----------------------------------------------------------------------------------------------------
	SUBROUTINE  check_tidal_threshold_sub(Y, Return_pass_tidalthreshold_yesno)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(tot_nr_Y_evol_eqs),	intent(in)		:: Y
		integer,								intent(out)		:: Return_pass_tidalthreshold_yesno
		real*8, dimension(n_particles,length_Y_per_body_n)		:: body_all_all
		real*8, dimension(n_particles,3)						:: pos
		real*8													:: r_ij, T_ij
		real*8													:: acc_tid_ij, acc_bind_i
		!------------------------------------------------------------
		!UNPACK Y:
		!------------------------------------------------------------
		body_all_all(:,:) 	= RESHAPE(Y, (/ n_particles,length_Y_per_body_n /), 	ORDER = (/2,1/))
		pos(:,:)			= body_all_all(:,1:3)	!we here just use short notation: pos
		!------------------------------------------------------------
		!------------------------------------------------------------
		!Check tidal threshold:
		!------------------------------------------------------------
		!---------------------------------------
		!initialize:
		!---------------------------------------
		Return_pass_tidalthreshold_yesno	= 0											!NO - tidal threshold is not passed.
		!---------------------------------------
		!---------------------------------------
		!loop over obj pairs:
		!---------------------------------------
		do i=1, n_particles,1
		do j=1, n_particles,1
		if (i .NE. j) then
		if (Return_pass_tidalthreshold_yesno .EQ. 0) then
			r_ij		= len3vec(pos(i,:)-pos(j,:))
			acc_tid_ij	= mass(j)*(1./((r_ij-radius(i))**2.) - 1./((r_ij+radius(i))**2.))
			acc_bind_i	= mass(i)/(radius(i)**2.)
			T_ij		= ABS(acc_tid_ij/acc_bind_i)		!F_tide/F_binding of obj_i caused by obj_j
			if (T_ij .GT. evolvetides_threshold .and. evoTides_yesno(i) .EQ. 1) Return_pass_tidalthreshold_yesno = 1	!YES - tidal threshold is now passed.
		endif	!if tidal threshold is not passed (= 0)	 
		endif	!i NE j
		enddo	!loop over j
		enddo	!loop over i
		!---------------------------------------
		!------------------------------------------------------------
		!------------------------------------------------------------
		!Retun info:
		!------------------------------------------------------------
		Return_pass_tidalthreshold_yesno	= Return_pass_tidalthreshold_yesno
		!------------------------------------------------------------
	END SUBROUTINE check_tidal_threshold_sub
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------
	
	
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE  analyze_3body_state_info_sub(Y, Return_3body_Info_REAL, Return_3body_Info_INT, Return_3body_Info_REAL_XTRAINFO)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(tot_nr_Y_evol_eqs),	intent(in)		:: Y
		real*8, dimension(10),					intent(out)		:: Return_3body_Info_REAL
		integer, dimension(10),					intent(out)		:: Return_3body_Info_INT
		real*8, dimension(10),					intent(out)		:: Return_3body_Info_REAL_XTRAINFO
		real*8, dimension(n_particles,length_Y_per_body_n)		:: body_all_all
		real*8, dimension(n_particles,3)						:: pos, vel
		real*8, dimension(n_particles,3,3)						:: body_all_q
		integer													:: end_state_flag
		real*8													:: r_ij, radius_tot_ij
		real*8, dimension(5)									:: binary_info_arr_ij, binary_info_arr_ijk
		real*8, dimension(5)									:: binary_info_arr_12, binary_info_arr_13, binary_info_arr_23
		real*8													:: CM_Mtot_ij
		real*8, dimension(3)									:: CM_ij_pos, CM_ij_vel
		real*8, dimension(3,2)									:: pair_index_ij
		real*8, dimension(3)									:: Force_arr
		integer													:: bin_i, bin_j, sin_k
		integer, dimension(1)									:: maxindex
		integer													:: pair_ij_bound_yesno, sin_k_unbound_pair_ij_yesno
		integer													:: sin_k_awayfrom_CMij_yesno, Forcefrac_ENDbinsin_yesno, Forcefrac_IMSbinsin_yesno
		real*8													:: E_tot_ij, a_ij, e_ij
		real*8													:: E_tot_ij_k
		real*8, dimension(3)									:: pos_sin_k_wrtCMij, vel_sin_k_wrtCMij
		real*8, dimension(3)									:: out_pos_CMij_wrt_sink, out_vel_CMij_wrt_sink
		real*8													:: r_dot_v
		real*8													:: F_bin_ij, F_tid_ij_k, Forcefrac_bin_tid
		integer													:: out_end_state_flag, out_bin_i, out_bin_j, out_sin_k
		real*8, dimension(5)									:: out_binary_info_arr_ij, out_binary_info_arr_ijk
		real*8, dimension(3)									:: Etot_ij, Etot_ijk, CM_pos, CM_vel
		real*8													:: Mtot, rdotv_wrtCM_1, rdotv_wrtCM_2, rdotv_wrtCM_3
		integer													:: allpairs_unbound_yesno, all_awayfrom_CM_yesno
		integer													:: out_IMS_bin_yesno
		real*8													:: out_dist_bin_ij
		integer, dimension(2)									:: bin_ij_arr
		integer													:: x, bin_x, obj_id_TD
		real*8, dimension(3,3)									:: q_x
		real*8, dimension(3)									:: PA_a1a2a3
		real*8													:: max_a1a2a3
		real*8													:: bin_ij_Etot, bin_ij_a, bin_ij_tot_rad
		real*8													:: out_dist_bin_12, out_dist_bin_13, out_dist_bin_23
		
		!This section could be uptimized.		
		
		!------------------------------------------------------------
		!UNPACK Y:
		!------------------------------------------------------------
		body_all_all(:,:) 	= RESHAPE(Y, (/ n_particles,length_Y_per_body_n /), 	ORDER = (/2,1/))
		pos(:,:)			= body_all_all(:,1:3)	!we here just use short notation: pos
		vel(:,:)			= body_all_all(:,4:6)	!we here just use short notation: vel
		body_all_q(:,:,:)	= RESHAPE(body_all_all(:,7:15),	(/ n_particles,3,3 /), 	ORDER = (/1,2,3/))		
		!------------------------------------------------------------
		
		
		!------------------------------------------------------------
		!Initilize:
		!------------------------------------------------------------
		!end state:
		end_state_flag = 0
		!return info (0 until endstate is found):
		out_end_state_flag			= 0
		out_bin_i 					= 0
		out_bin_j 					= 0
		out_sin_k 					= 0
		out_binary_info_arr_ij(:)	= 0d0
		out_binary_info_arr_ijk(:)	= 0d0
		Return_3body_Info_INT(:)	= 0
		Return_3body_Info_REAL(:)	= 0d0
		out_IMS_bin_yesno			= 0
		!------------------------------------------------------------
		
		
		!------------------------------------------------------------
		!Find bin-sin candidate ([i,j],k):
		!------------------------------------------------------------
		nc = 1	!counter
		do i=1,		n_particles,1
		do j=i+1,	n_particles,1
			r_ij				= len3vec(pos(i,:)-pos(j,:))
			Force_arr(nc)		= mass(i)*mass(j)/(r_ij**2)
			pair_index_ij(nc,:)	= [i,j]
			nc = nc+1
		enddo
		enddo
		maxindex	= MAXLOC(Force_arr)	!returns 1D array with arr loc of max element
		!bin-sin candidate indices:
		bin_i 	= pair_index_ij(maxindex(1),1)	!	bin i
		bin_j	= pair_index_ij(maxindex(1),2)	!	bin j
		sin_k	= 6 - (bin_i+bin_j)				!	sin k
		!------------------------------------------------------------
		
		
		!------------------------------------------------------------
		!Relative pos and vel between bin[i,j] and sin[k] (and initialize):
		!------------------------------------------------------------
		!calc CMij properties:
		CM_Mtot_ij	= (mass(bin_i) + mass(bin_j))										!written over later
		CM_ij_pos	= (mass(bin_i)*pos(bin_i,:) + mass(bin_j)*pos(bin_j,:))/CM_Mtot_ij	!written over later
		CM_ij_vel	= (mass(bin_i)*vel(bin_i,:) + mass(bin_j)*vel(bin_j,:))/CM_Mtot_ij	!written over later
		!calc rel pos and vel:	
		out_pos_CMij_wrt_sink	= CM_ij_pos - pos(sin_k,:)
		out_vel_CMij_wrt_sink	= CM_ij_vel - vel(sin_k,:)
		!initialize output:
		CALL	Calc_binary_info(pos(bin_i,:),	vel(bin_i,:),	mass(bin_i),	pos(bin_j,:),	vel(bin_j,:),	mass(bin_j),	binary_info_arr_ij)		! system [i,j]
		CALL	Calc_binary_info(pos(sin_k,:),	vel(sin_k,:),	mass(sin_k),	CM_ij_pos,		CM_ij_vel,		CM_Mtot_ij,		binary_info_arr_ijk)	! system [[i,j],k]
		out_bin_i 				= bin_i
		out_bin_j 				= bin_j
		out_sin_k 				= sin_k
		out_binary_info_arr_ij	= binary_info_arr_ij
		out_binary_info_arr_ijk	= binary_info_arr_ijk
		out_dist_bin_ij			= len3vec(pos(bin_i,:)-pos(bin_j,:))
		!------------------------------------------------------------
		
		
		!------------------------------------------------------------
		!Collisions:
		!------------------------------------------------------------
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.
		!------------------------------------------------------------
		!--------------------------------------------------------
		!loop over obj i,j:
		!--------------------------------------------------------
		do i=1, n_particles,1
		do j=1, n_particles,1
		if (i .NE. j) then
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.
			!----------------------------------------------------
			!Proceed if endstate is not found:
			!----------------------------------------------------			
			!-----------------------------------
			!Check if: obj i,j is overlapping:
			!-----------------------------------		
			!calc:
			r_ij			= len3vec(pos(i,:)-pos(j,:))
			radius_tot_ij	= radius(i) + radius(j)				
			!COLLISION:
			if (r_ij .LT. radius_tot_ij)	end_state_flag = 2	!COLLISION
			!-----------------------------------	
			!-----------------------------------
			!If end-state is found:
			!-----------------------------------
			if (end_state_flag .NE. 0) then			
				!Index of single obj:
				k = 6 - (i+j)			
				!calc CM properties for pair [i,j]:
				CM_Mtot_ij	= (mass(i) + mass(j))
				CM_ij_pos	= (mass(i)*pos(i,:) + mass(j)*pos(j,:))/CM_Mtot_ij
				CM_ij_vel	= (mass(i)*vel(i,:) + mass(j)*vel(j,:))/CM_Mtot_ij
				!Get Binary Info for [i,j] and [[i,j],k]:
				!Format:	Calc_binary_info(pos_1,		vel_1,		mass_1,		pos_2,		vel_2,		mass_2,		binary_info_arr) - binary_info_arr: [E_kin, E_pot, E_tot, a_bin, e_bin]
				CALL		Calc_binary_info(pos(i,:),	vel(i,:),	mass(i),	pos(j,:),	vel(j,:),	mass(j),	binary_info_arr_ij)		! system [i,j]		
				CALL		Calc_binary_info(pos(k,:),	vel(k,:),	mass(k),	CM_ij_pos,	CM_ij_vel,	CM_Mtot_ij,	binary_info_arr_ijk)	! system [[i,j],k]
				!set final end-state output from this module:
				out_end_state_flag		= end_state_flag
				out_bin_i 				= i
				out_bin_j 				= j
				out_sin_k 				= k
				out_binary_info_arr_ij	= binary_info_arr_ij
				out_binary_info_arr_ijk	= binary_info_arr_ijk
			endif
			!-----------------------------------
			!----------------------------------------------------
		ENDIF	!endstate	
		endif	!i NE j
		enddo	!loop over j
		enddo	!loop over i
		!--------------------------------------------------------
		!------------------------------------------------------------
		ENDIF	!end_state_flag .EQ. 0
		!------------------------------------------------------------


		!------------------------------------------------------------
		!Disruption / Tidal Inspiral binary:
		!------------------------------------------------------------
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.
		!------------------------------------------------------------
		!----------------------------------------------------
		!Check for tidal disruptions:
		!----------------------------------------------------
		if (pass_tidalthreshold_yesno .EQ. 1) then
		!----------------------------------------------------
		!Define:
		bin_ij_arr(:)	= [bin_i,bin_j]
		obj_id_TD		= -1
		!Loop over the two bin objs:
		do x=1, 2,1
			!-----------------------------------
			!get obj x id index:
			!-----------------------------------
			bin_x = bin_ij_arr(x)
			!-----------------------------------
			!calc:
			!-----------------------------------
			!Read in q matrix:
			q_x = body_all_q(bin_x,:,:)
			!calc principal axis (a1,a2,a3):
			CALL calc_stellar_paxis_a1a2a3(q_x, PA_a1a2a3)
			max_a1a2a3 = maxval(PA_a1a2a3)
			!-----------------------------------
			!check if tidal threshold is passed:
			!-----------------------------------
			if (max_a1a2a3 .GT. tidaldisrup_threshold) then
				end_state_flag = 1		!TIDAL DISRUPTION
				obj_id_TD = bin_x
			endif
			!-----------------------------------		
		enddo
		!----------------------------------------------------	
		endif	!if we evolve tides
		!----------------------------------------------------
		!----------------------------------------------------
		!Check for Inspiral state:	
		!----------------------------------------------------	
		!calc:
		CALL Calc_binary_info(pos(bin_i,:),	vel(bin_i,:),	mass(bin_i),	pos(bin_j,:),	vel(bin_j,:),	mass(bin_j),	binary_info_arr_ij)		! system [i,j]		
		!Define:
		bin_ij_Etot		= binary_info_arr_ij(3)
		bin_ij_a		= binary_info_arr_ij(4)   	
		bin_ij_tot_rad	= radius(bin_i) + radius(bin_j)	!unperturbed radii (sph radii)
		if (bin_ij_Etot .LT. 0.0 .and. (bin_ij_a/bin_ij_tot_rad) .LT. insp_threshold) then
			end_state_flag = 5			!INSPIRAL STATE	
		endif
		!----------------------------------------------------	
		!----------------------------------------------------
		!If end-state is found:
		!----------------------------------------------------
		if (end_state_flag .NE. 0) then
			!calc CM properties for pair [i,j]:
			CM_Mtot_ij	= (mass(bin_i) + mass(bin_j))
			CM_ij_pos	= (mass(bin_i)*pos(bin_i,:) + mass(bin_j)*pos(bin_j,:))/CM_Mtot_ij
			CM_ij_vel	= (mass(bin_i)*vel(bin_i,:) + mass(bin_j)*vel(bin_j,:))/CM_Mtot_ij
			!Get Binary Info for [i,j] and [[i,j],k]:
			!Format:	Calc_binary_info(pos_1,		vel_1,		mass_1,		pos_2,		vel_2,		mass_2,		binary_info_arr) - binary_info_arr: [E_kin, E_pot, E_tot, a_bin, e_bin]
			CALL Calc_binary_info(pos(bin_i,:),	vel(bin_i,:),	mass(bin_i),	pos(bin_j,:),	vel(bin_j,:),	mass(bin_j),	binary_info_arr_ij)		! system [i,j]		
			CALL Calc_binary_info(pos(sin_k,:),	vel(sin_k,:),	mass(sin_k),	CM_ij_pos,		CM_ij_vel,		CM_Mtot_ij,		binary_info_arr_ijk)	! system [[i,j],k]
			!set final end-state output from this module:
			out_end_state_flag		= end_state_flag
			out_bin_i 				= bin_i
			out_bin_j 				= bin_j
			out_sin_k 				= obj_id_TD	!WE HERE SAVE NOT THE SINGLE OBJ, BUT THE ONE THAT IS TIDALLY DISRUPTED!!! SHOULD BE OK!
			out_binary_info_arr_ij	= binary_info_arr_ij
			out_binary_info_arr_ijk	= binary_info_arr_ijk	
		endif				
		!----------------------------------------------------
		!------------------------------------------------------------
		ENDIF	!end_state_flag .EQ. 0
		!------------------------------------------------------------
				
		
		!------------------------------------------------------------
		!Bin-Sin state:
		!------------------------------------------------------------
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.
		!------------------------------------------------------------
		!----------------------------------------------------
		!calc CMij properties:
		!----------------------------------------------------
		CM_Mtot_ij	= (mass(bin_i) + mass(bin_j))
		CM_ij_pos	= (mass(bin_i)*pos(bin_i,:) + mass(bin_j)*pos(bin_j,:))/CM_Mtot_ij
		CM_ij_vel	= (mass(bin_i)*vel(bin_i,:) + mass(bin_j)*vel(bin_j,:))/CM_Mtot_ij
		!----------------------------------------------------
		!Calc properties of bin-sin candidate:
		!----------------------------------------------------
		!-----------------------------------
		!Initialize check flag params:
		!-----------------------------------
		pair_ij_bound_yesno			= 0
		sin_k_unbound_pair_ij_yesno	= 0
		sin_k_awayfrom_CMij_yesno	= 0
		Forcefrac_ENDbinsin_yesno	= 0
		Forcefrac_IMSbinsin_yesno	= 0
		!-----------------------------------
		!Check if: bin[i,j] are bound to each other
		!-----------------------------------
		CALL	Calc_binary_info(pos(bin_i,:),	vel(bin_i,:),	mass(bin_i),	pos(bin_j,:),	vel(bin_j,:),	mass(bin_j),	binary_info_arr_ij)		! system [i,j]
		E_tot_ij	= binary_info_arr_ij(3)
		a_ij		= binary_info_arr_ij(4)	!used later in this section
		e_ij		= binary_info_arr_ij(5)	!used later in this section
		if (E_tot_ij .LT. 0.0)							pair_ij_bound_yesno			= 1	! - ij bound
		!-----------------------------------
		!Check if: sin[k] is unbound to CM bin[i,j]
		!-----------------------------------
		CALL	Calc_binary_info(pos(sin_k,:),	vel(sin_k,:),	mass(sin_k),	CM_ij_pos,		CM_ij_vel,		CM_Mtot_ij,		binary_info_arr_ijk)	! system [[i,j],k]
		E_tot_ij_k	= binary_info_arr_ijk(3)
		if (E_tot_ij_k .GT. 0.0)						sin_k_unbound_pair_ij_yesno	= 1	! - ij,k unbound
		!-----------------------------------
		!Check if: sin[k] is traveling away from CM bin[i,j]
		!-----------------------------------
		pos_sin_k_wrtCMij	= pos(sin_k,:) - CM_ij_pos
		vel_sin_k_wrtCMij	= vel(sin_k,:) - CM_ij_vel
		r_dot_v				= DOT_PRODUCT(pos_sin_k_wrtCMij,vel_sin_k_wrtCMij)
		if (r_dot_v .GT. 0.0)							sin_k_awayfrom_CMij_yesno	= 1	! - k away from i,j
		!-----------------------------------
		!Check if: Tidal force on bin[i,j] is < 'threshold'
		!-----------------------------------
		F_bin_ij	= ABS(mass(bin_i)*mass(bin_j)/((a_ij*(1d0+e_ij))**2))
		F_tid_ij_k	= ABS((2d0*CM_Mtot_ij*mass(sin_k)/((len3vec(pos(sin_k,:)-CM_ij_pos))**3))*(a_ij*(1d0+e_ij)))
		Forcefrac_bin_tid	= F_tid_ij_k/F_bin_ij
		if (Forcefrac_bin_tid .LT. ENDbinsingle_threshold) 	Forcefrac_ENDbinsin_yesno		= 1	! - ENDstate bin,sin state
		if (Forcefrac_bin_tid .LT. IMSbinsingle_threshold) 	Forcefrac_IMSbinsin_yesno		= 1	! - IMSstate bin,sin state
		!-----------------------------------
		!----------------------------------------------------
		!----------------------------------------------------
		!Check if: bin-sin state is an endstate
		!----------------------------------------------------
		if (pair_ij_bound_yesno			.EQ. 1 .and. &
			sin_k_unbound_pair_ij_yesno	.EQ. 1 .and. &
			sin_k_awayfrom_CMij_yesno	.EQ. 1 .and. &
			Forcefrac_ENDbinsin_yesno	.EQ. 1) then
		!BIN-SINGLE:
		end_state_flag		= 3	
		endif
		!----------------------------------------------------
		!Check if: bin-sin state is a resonance IMS:
		!----------------------------------------------------
		if (pair_ij_bound_yesno			.EQ. 1 .and. &
			sin_k_unbound_pair_ij_yesno	.EQ. 0 .and. & 
			Forcefrac_IMSbinsin_yesno 	.EQ. 1) then
		!resonance IMS:
		out_IMS_bin_yesno	= 1
		endif
		!----------------------------------------------------
		!If end-state is found:
		!----------------------------------------------------
		if (end_state_flag .NE. 0) then
			!set final end-state output from this module:
			out_end_state_flag		= end_state_flag
			out_bin_i 				= bin_i
			out_bin_j 				= bin_j
			out_sin_k 				= sin_k
			out_binary_info_arr_ij	= binary_info_arr_ij
			out_binary_info_arr_ijk	= binary_info_arr_ijk	
		endif
		!----------------------------------------------------
		!------------------------------------------------------------
		ENDIF	!end_state_flag .EQ. 0
		!------------------------------------------------------------


		!------------------------------------------------------------
		!Ionization:
		!------------------------------------------------------------
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.
		!------------------------------------------------------------
		!----------------------------------------------------
		!Initialize check flag params:
		!----------------------------------------------------
		allpairs_unbound_yesno	= 0
		all_awayfrom_CM_yesno	= 0
		!----------------------------------------------------
		!----------------------------------------------------
		!Calc pairwise ([i,j]), ([i,j],k) Etot:
		!----------------------------------------------------
		nc = 1
		do i=1,		n_particles,1
		do j=i+1,	n_particles,1
			k = 6 - (i+j) !single
			CM_Mtot_ij	= mass(i)+mass(j)
			CM_ij_pos	= (mass(i)*pos(i,:)+mass(j)*pos(j,:))/CM_Mtot_ij
			CM_ij_vel	= (mass(i)*vel(i,:)+mass(j)*vel(j,:))/CM_Mtot_ij
			CALL Calc_binary_info(pos(k,:), vel(k,:), mass(k), CM_ij_pos, CM_ij_vel, CM_Mtot_ij, binary_info_arr_ijk)
			Etot_ijk(nc)	= binary_info_arr_ijk(3) 
			CALL Calc_binary_info(pos(i,:), vel(i,:), mass(i), pos(j,:), vel(j,:), mass(j), binary_info_arr_ij)
			Etot_ij(nc)		= binary_info_arr_ij(3)
			nc = nc+1
		enddo
		enddo
		if (Etot_ijk(1) .GT. 0.0 .AND. Etot_ijk(2)	.GT. 0.0 .AND. Etot_ijk(3)	.GT. 0.0 .AND. &
			Etot_ij(1) 	.GT. 0.0 .AND. Etot_ij(2)	.GT. 0.0 .AND. Etot_ij(3)	.GT. 0.0) 		allpairs_unbound_yesno	= 1
		!----------------------------------------------------
		!vel direction wrt CM:
		!----------------------------------------------------
		Mtot	= mass(1) + mass(2) + mass(3)									!mass total NB system
		CM_pos	= (mass(1)*pos(1,:) + mass(2)*pos(2,:) + mass(3)*pos(3,:))/Mtot	!pos CM of total NB system
		CM_vel	= (mass(1)*vel(1,:) + mass(2)*vel(2,:) + mass(3)*vel(3,:))/Mtot	!vel CM of total NB system
		rdotv_wrtCM_1 = DOT_PRODUCT((pos(1,:)-CM_pos),(vel(1,:)-CM_vel))
		rdotv_wrtCM_2 = DOT_PRODUCT((pos(2,:)-CM_pos),(vel(2,:)-CM_vel))
		rdotv_wrtCM_3 = DOT_PRODUCT((pos(3,:)-CM_pos),(vel(3,:)-CM_vel))
		if (rdotv_wrtCM_1 .GT. 0.0 .AND. rdotv_wrtCM_2 .GT. 0.0 .AND. rdotv_wrtCM_3 .GT. 0.0)	all_awayfrom_CM_yesno	= 1
		!----------------------------------------------------
		!Check for end-state:
		!----------------------------------------------------
		if (allpairs_unbound_yesno .EQ. 1 .AND. all_awayfrom_CM_yesno .EQ. 1) then
			!IONIZATION:
			end_state_flag = 4
		endif	
		!----------------------------------------------------
		!If end-state is found:
		!----------------------------------------------------
		if (end_state_flag .NE. 0) then
			!set final end-state output from this module:
			out_end_state_flag		= end_state_flag
			!rest of info not important.
		endif
		!----------------------------------------------------
		!------------------------------------------------------------
		ENDIF	!end_state_flag .EQ. 0
		!------------------------------------------------------------


		!------------------------------------------------------------
		!xtra info:
		!------------------------------------------------------------
		out_dist_bin_12			= len3vec(pos(1,:)-pos(2,:))
		out_dist_bin_13			= len3vec(pos(1,:)-pos(3,:))
		out_dist_bin_23 		= len3vec(pos(2,:)-pos(3,:))		
		!------------------------------------------------------------
			

		!------------------------------------------------------------
		!Return info:
		!------------------------------------------------------------
		Return_3body_Info_INT(1:5)				= [out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno]
		Return_3body_Info_REAL(1:10)			= [out_binary_info_arr_ij(:), out_binary_info_arr_ijk(:)]
		Return_3body_Info_REAL_XTRAINFO(1:4)	= [out_dist_bin_ij, out_dist_bin_12, out_dist_bin_13, out_dist_bin_23]				
		Return_3body_Info_REAL_XTRAINFO(5:10)	= [out_pos_CMij_wrt_sink(:), out_vel_CMij_wrt_sink(:)]
		!------------------------------------------------------------


	END SUBROUTINE analyze_3body_state_info_sub
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------
	
	
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE  Nbody_endstate_sub(Y, Return_Nbody_endstate, Return_endstate_binparams)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(tot_nr_Y_evol_eqs),	intent(in)		:: Y
		integer, dimension(5),					intent(out)		:: Return_Nbody_endstate
		real*8, dimension(10),                   intent(out) 	:: Return_endstate_binparams
		real*8, dimension(n_particles,length_Y_per_body_n)		:: body_all_all
		real*8, dimension(n_particles,3)						:: pos, vel
		real*8, dimension(n_particles,3,3)						:: body_all_q
        real*8, dimension(5)                                    :: binary_info_arr_in, binary_info_arr_out
        real*8, dimension(5)                                    :: binary_info_arr_ijkl
        real*8, dimension(5)                                    :: binary_info_arr_ijk, binary_info_arr_ijl, binary_info_arr_ikl
        real*8, dimension(5)                                    :: binary_info_arr_ij, binary_info_arr_ik, binary_info_arr_il
        real*8, dimension(5)                                    :: binary_info_arr_kl
		integer													:: end_state_flag
		integer													:: out_end_state_flag, out_bin_i, out_bin_j, out_bin_k, out_bin_l
		real*8													:: r_ij, radius_tot_ij, r_12, r_kl, r_ijk, r_ijl, r_ijkl
		real*8, dimension(3,3)									:: q_in
		real*8, dimension(3)									:: PA_a1a2a3
		real*8													:: max_a1a2a3
        real*8                                                  :: KE, PE, eps, Mtot, M_in, v_dir, v_dir_ij, v_dir_kl, v_dir_k
        real*8                                                  :: KE1, KE2, M1, M2, vCM
		real*8, dimension(3)									:: r_CM, v_CM, CM, CM_num, vCM_num
		real*8, dimension(3)									:: CM1, CM2, vCM1, vCM2
		real*8, dimension(3)									:: CM_pos, CM_vel, Lin_vec, Lout_vec, Lvec
        real*8                                                  :: fi, fj, fk, fl
        real*8                                                  :: mass_bin_i, mass_bin_j, mass_bin_k, mass_bin_l
        real*8                                                  :: a_bin, e_bin, a_bin_out, e_bin_out, inc_bin
        real*8                                                  :: a_in, e_in, a_out, e_out, inc, q_out
        real*8                                                  :: E_kin, E_pot, E_tot, E_check, L_check, L_ini, E_ini
        real*8, dimension(3)                                    :: L_check_vec, L_ini_vec
        real*8                                                  :: Etot_ij, Etot_ik, Etot_il, Etot_kl
        real*8                                                  :: Etot_ijkl, Etot_ijk, Etot_ijl, Etot_ikl
        real*8                                                  :: term1, term2, Ft, Ft1, Ft2, delta_F, delta_EL, delta_F_triple
        integer                                                 :: unbound, ub1, ub2, ub3, ub4
        integer                                                 :: one_unbound_yesno, two_unbound_yesno, all_unbound_yesno
        integer                                                 :: stability_crit_yesno
		
		!------------------------------------------------------------
		!UNPACK Y:
		!------------------------------------------------------------
		body_all_all(:,:) 	= RESHAPE(Y, (/ n_particles,length_Y_per_body_n /), 	ORDER = (/2,1/))
		pos(:,:)			= body_all_all(:,1:3)	!we here just use short notation: pos
		vel(:,:)			= body_all_all(:,4:6)	!we here just use short notation: vel
		body_all_q(:,:,:)	= RESHAPE(body_all_all(:,7:15),	(/ n_particles,3,3 /), 	ORDER = (/1,2,3/))		
		!------------------------------------------------------------
		!Initialize:
		!------------------------------------------------------------
		end_state_flag		= 0
		out_bin_i 			= 0
		out_bin_j 			= 0
		out_bin_k 			= 0
		out_bin_l 			= 0
        mass_bin_i          = 0.0
        mass_bin_j          = 0.0
        mass_bin_k          = 0.0
        mass_bin_l          = 0.0
        a_bin               = 0.0
        e_bin               = 0.0
        a_bin_out           = 0.0
        e_bin_out           = 0.0
        inc_bin             = 0.0
        delta_F             = 1e-5                      ! Ftid/Frel threshold
        delta_EL            = 1e-3                      ! E,L threshold
		!------------------------------------------------------------
		!Check for tidal disruptions:
		!------------------------------------------------------------
		if (pass_tidalthreshold_yesno .EQ. 1) then
		!---------------------------------------
		!loop over objs:
		!---------------------------------------
		do i=1, n_particles,1
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.
			!-----------------------------------
			!calc:
			!-----------------------------------
			!Read in q matrix:
			q_in = body_all_q(i,:,:)
			!calc principal axis (a1,a2,a3):
			CALL calc_stellar_paxis_a1a2a3(q_in, PA_a1a2a3)
			max_a1a2a3 = maxval(PA_a1a2a3)
			!-----------------------------------
			!check if tidal threshold is passed:
			!-----------------------------------
			if (max_a1a2a3 .GT. tidaldisrup_threshold) then
                end_state_flag = 1	!TIDAL DISRUPTION			
				out_bin_i 			= i
				out_bin_j 			= -1
            endif
			!-----------------------------------
		ENDIF	!endstate
		enddo	!obj loop
		!---------------------------------------
		endif	!if we evolve tides
		!------------------------------------------------------------


		!------------------------------------------------------------
		!Check for collisions (endstate = 2)
		!------------------------------------------------------------
		!---------------------------------------
		!loop over obj pairs:
		!---------------------------------------
		do i=1, n_particles,1
		do j=i+1, n_particles,1
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.
			!-----------------------------------
			!calc:
			!-----------------------------------
			r_ij			= len3vec(pos(i,:)-pos(j,:))
			radius_tot_ij	= radius(i) + radius(j)	
			!-----------------------------------
			!COLLISION:
			!-----------------------------------
			if (r_ij .LT. radius_tot_ij) then
                end_state_flag = 2	!COLLISION
				out_bin_i 			= i
				out_bin_j 			= j
                mass_bin_i          = mass(i)
                mass_bin_j          = mass(j)
                ! calculate pertinent info for output (E_kin, E_pot, E_tot, a_bin, e_bin)
                CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), pos(j,:), vel(j,:), mass(j), binary_info_arr_ij)		! system [i,j]
                a_bin = binary_info_arr_ij(4)
                e_bin = binary_info_arr_ij(5)
			endif
			!-----------------------------------
		ENDIF	!endstate
		enddo	!loop over j
		enddo	!loop over i


		!------------------------------------------------------------
		!Check for inspirals (endstate = 5)
		!------------------------------------------------------------
		do i=1, n_particles,1
		do j=1, n_particles,1
		if (i .NE. j) then
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.
			!-----------------------------------
			!calc:
			!-----------------------------------
			radius_tot_ij	= radius(i) + radius(j)	
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), pos(j,:), vel(j,:), mass(j), binary_info_arr_ij)		! system [i,j]
            E_tot	    	= binary_info_arr_ij(3)
            a_bin   		= binary_info_arr_ij(4)   	
            e_bin           = binary_info_arr_ij(5)
            if (E_tot .LT. 0 .and. (a_bin/radius_tot_ij) .LT. insp_threshold) then
                end_state_flag = 5			!INSPIRAL STATE	
				out_bin_i 			= i
				out_bin_j 			= j
                mass_bin_i          = mass(i)
                mass_bin_j          = mass(j)
			endif
			!-----------------------------------
		ENDIF	!endstate
		endif	!i NE j
		enddo	!loop over j
		enddo	!loop over i


		!------------------------------------------------------------
        !CHECK FOR AN UNBOUND OBJECT IN 3BODY CASE!
		!Check for bound binary with 1 unbound object (endstate = 4)
		!------------------------------------------------------------

		do i=1, n_particles,1
        do j=i+1, n_particles,1
        k = 6 - (i+j)
		IF (end_state_flag .EQ. 0 .AND. n_particles .EQ. 3) THEN		!if endstate not found

        CM_pos = CoM_3body(pos(i,:), mass(i), pos(j,:), mass(j), pos(k,:), mass(k))
        CM_vel = CoM_3body(vel(i,:), mass(i), vel(j,:), mass(j), vel(k,:), mass(k))
        CM1 = CoM_2body(pos(i,:), mass(i), pos(j,:), mass(j))
        vCM1 = CoM_2body(vel(i,:), mass(i), vel(j,:), mass(j))
        Mtot = mass(i) + mass(j)
        CALL    Calc_binary_info(pos(k,:), vel(k,:), mass(k), CM1, vCM1, Mtot, binary_info_arr_ijk)
        Etot_ijk = binary_info_arr_ijk(3)

        CALL   Calc_binary_info(pos(i,:), vel(i,:), mass(i), pos(j,:), vel(j,:), mass(j), binary_info_arr_ij)
        Etot_ij = binary_info_arr_ij(3)
        a_bin = binary_info_arr_ij(4)
        e_bin = binary_info_arr_ij(5)

        v_dir_ij = DOT_PRODUCT((CM1-CM_pos), (vCM1-CM_vel))
        v_dir_k = DOT_PRODUCT((pos(k,:)-CM_pos), (vel(k,:)-CM_vel))

        if (Etot_ijk .GT. 0.0 .AND. Etot_ij .LT. 0.0 .AND. v_dir_ij .GT. 0 .AND. v_dir_k .GT. 0) then
            
                ! see if binary [i,j] is bound below tidal threshold
                r_ijk = len3vec(CM1-pos(k,:))
                Ft = Fthresh(mass(i),mass(j),mass(k),a_bin,e_bin,r_ijk)
                if (Ft .LT. delta_F) then
                    end_state_flag = 4     ! bound binary with unbound object
                    out_bin_i 			= i
                    out_bin_j 			= j
                    mass_bin_i          = mass(i)
                    mass_bin_j          = mass(j)
                endif
                    

        endif !force check

        ENDIF !endstate

        enddo !loop over j
        enddo !loop over i



		!------------------------------------------------------------
        !CHECK FOR UNBOUND OBJECTS IN 4BODY CASE!
		!Check for stable triple with 1 unbound object (endstate = 3)
		!Check for bound binary with 2 unbound object (endstate = 4)
		!Check for total ionization (endstate = 7)
		!------------------------------------------------------------
		IF (end_state_flag .EQ. 0 .AND. n_particles .EQ. 4) THEN		!if endstate not found

		!Initialize check flag params:
        unbound = 0             ! counter for number of unbound objects

        !See how many objects in the system are unbound
		do i=1, n_particles,1
        do j=1, n_particles,1
        do k=1, n_particles,1
        do l=1, n_particles,1
        if (i .NE. j .AND. i .NE. k .AND. i .NE. l .AND. &
            j .NE. k .AND. j .NE. l .AND. k .NE. l .AND. &
            j .LT. k .AND. k .LT. l) then

            !FIXME how are things unbound then bound again? does this make sense?
            !particle i against system jkl
            CM_pos = CoM_3body(pos(j,:), mass(j), pos(k,:), mass(k), pos(l,:), mass(l))
            CM_vel = CoM_3body(vel(j,:), mass(j), vel(k,:), mass(k), vel(l,:), mass(l))
            Mtot = mass(j) + mass(k) + mass(l)
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), CM_pos, CM_vel, Mtot, binary_info_arr_ijkl)
            Etot_ijkl = binary_info_arr_ijkl(3)

            !particle i against system jk
            CM_pos = CoM_2body(pos(j,:), mass(j), pos(k,:), mass(k))
            CM_vel = CoM_2body(vel(j,:), mass(j), vel(k,:), mass(k))
            Mtot = mass(j) + mass(k)
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), CM_pos, CM_vel, Mtot, binary_info_arr_ijk)
            Etot_ijk = binary_info_arr_ijk(3)

            !particle i against system jl
            CM_pos = CoM_2body(pos(j,:), mass(j), pos(l,:), mass(l))
            CM_vel = CoM_2body(vel(j,:), mass(j), vel(l,:), mass(l))
            Mtot = mass(j) + mass(l)
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), CM_pos, CM_vel, Mtot, binary_info_arr_ijl)
            Etot_ijl = binary_info_arr_ijl(3)

            !particle i against system kl
            CM_pos = CoM_2body(pos(k,:), mass(k), pos(l,:), mass(l))
            CM_vel = CoM_2body(vel(k,:), mass(k), vel(l,:), mass(l))
            Mtot = mass(k) + mass(l)
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), CM_pos, CM_vel, Mtot, binary_info_arr_ikl)
            Etot_ikl = binary_info_arr_ikl(3)

            !particle i against particle j
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), pos(j,:), vel(j,:), mass(j), binary_info_arr_ij)
            Etot_ij = binary_info_arr_ij(3)

            !particle i against particle k
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), pos(k,:), vel(k,:), mass(k), binary_info_arr_ik)
            Etot_ik = binary_info_arr_ik(3)

            !particle i against particle l
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), pos(l,:), vel(l,:), mass(l), binary_info_arr_il)
            Etot_il = binary_info_arr_il(3)

            !see if velocity of particle i is away from CM
            CM_pos = CoM_4body(pos(i,:), mass(i), pos(j,:), mass(j), pos(k,:), mass(k), pos(l,:), mass(l))
            CM_vel = CoM_4body(vel(i,:), mass(i), vel(j,:), mass(j), vel(k,:), mass(k), vel(l,:), mass(l))
            v_dir = DOT_PRODUCT((pos(i,:)-CM_pos), (vel(i,:)-CM_vel))

            if (Etot_ijkl .GT. 0 .AND. Etot_ijk .GT. 0 .AND. Etot_ijl .GT. 0 .AND. Etot_ikl .GT. 0 .AND.  &
                    Etot_ij .GT. 0 .AND. Etot_ik .GT. 0 .AND. Etot_il .GT. 0 .AND. v_dir .GT. 0) then
            !if (Etot_ijkl .GT. 0 .AND. v_dir .GT. 0) then
                unbound = unbound + 1
                if (unbound .EQ. 1)   ub1 = i
                if (unbound .EQ. 2)   ub2 = i
                if (unbound .EQ. 3)   ub3 = i
                if (unbound .EQ. 4)   ub4 = i

            endif
        endif   !if statement for independent particles
        enddo   !loop over l
        enddo   !loop over k
        enddo   !loop over j
        enddo   !loop over i
        ENDIF   !endstate


        
        ! 1 UNBOUND: see if the three remaining are in a stable triple system (only works for 4-body)
        if (unbound .EQ. 1 .AND. n_particles .EQ. 4) then
            do i=1, n_particles,1
            do j=i+1, n_particles,1
            do k=1, n_particles,1
            IF (end_state_flag .EQ. 0 .AND. n_particles .EQ. 4) THEN		!if endstate not found
            if (i .NE. j .AND. i .NE. k .AND. j .NE. k .AND. &
                        i .NE. ub1 .AND. j .NE. ub1 .AND. k .NE. ub1) then
                ! assume i,j makes up the inner binary, calculate binary info and angular momentum
                CALL Calc_binary_L(pos(i,:),vel(i,:),mass(i),pos(j,:),vel(j,:),mass(j), Lin_vec)
                CALL Calc_binary_info(pos(i,:),vel(i,:),mass(i),pos(j,:),vel(j,:),mass(j), binary_info_arr_in)
                a_in = binary_info_arr_in(4)
                e_in = binary_info_arr_in(5)

                ! calculate the angular momentum and binary info of the outer system
                CM_pos = CoM_2body(pos(i,:), mass(i), pos(j,:), mass(j))
                CM_vel = CoM_2body(vel(i,:), mass(i), vel(j,:), mass(j))
                M_in = mass(i) + mass(j)               
                Mtot = mass(i) + mass(j) + mass(k)
                CALL Calc_binary_L(CM_pos,CM_vel,M_in,pos(k,:),vel(k,:),mass(k), Lout_vec)
                CALL Calc_binary_info(CM_pos,CM_vel,M_in,pos(k,:),vel(k,:),mass(k), binary_info_arr_out)
                a_out = binary_info_arr_out(4)
                e_out = binary_info_arr_out(5)
            
                ! make sure a_in is less than a_out and inner system is bound
                if (a_in .GT. 0 .AND. a_in .LT. a_out) then

                    ! calculate inclination from angular momentum vectors
                    inc = ACOS(DOT_PRODUCT(Lin_vec,Lout_vec)/(len3vec(Lin_vec)*len3vec(Lout_vec)))

                    ! now, we see if a stable triple was formed...
                    ! Marding & Aarseth criteria:
                    term1 = a_out/a_in * (1d0-e_out)
                    q_out = mass(k)/M_in
                    term2 = 2.8 * ((1d0 + q_out)*((1+e_out)/SQRT(1-e_out))**(2d0/5d0) &
                                * (1d0 - (0.3/3.1415)*inc))

                    ! Threshold that Ftid/Frel < delta = 1e-5
                    CM_pos = CoM_3body(pos(i,:), mass(i), pos(j,:), mass(j), pos(k,:), mass(k))
                    r_ij  = len3vec(CM_pos-pos(ub1,:))
                    Ft = Fthresh(M_in, mass(k), mass(ub1), a_out, e_out, r_ij)

                    ! check if stable triple criteria is met, and write output variables
                    ! we also note that q_out must be less than or equal to 5
                    if (term1 .GT. term2 .AND. q_out .LE. 5 .AND. Ft .LT. delta_F) then
                        end_state_flag = 3			!STABLE TRIPLE STATE STATE	
                        out_bin_i 			= i
                        out_bin_j 			= j
                        out_bin_k 			= k
                        mass_bin_i          = mass(i)
                        mass_bin_j          = mass(j)
                        mass_bin_k          = mass(k)
                        a_bin               = a_in
                        a_bin_out           = a_out
                        e_bin               = e_in
                        e_bin_out           = e_out
                        inc_bin             = inc
                    endif !if flag is set
                endif !if a_in < a_out

            endif
            ENDIF!endstate
            enddo
            enddo
            enddo
        
        endif ! if statement for single unbound object




        ! 2 UNBOUND: see if the two remaining are bound
        if (unbound .EQ. 2 .AND. n_particles .EQ. 4) then
        IF (end_state_flag .EQ. 0 .AND. n_particles .EQ. 4) THEN		!if endstate not found

            ! write info for the bound binary
            do j=1, n_particles, 1
                if (j .NE. ub1 .AND. j .NE. ub2) then
                    out_bin_j = j
                    mass_bin_j = mass(j)
                endif
            enddo
            do i=1, n_particles, 1
                if (i .NE. ub1 .AND. i .NE. ub2 .AND. i .NE. out_bin_j) then
                    out_bin_i = i
                    mass_bin_i = mass(i)
                endif
            enddo
                
            ! if two objects are unbound, check the force contribution on the remaining binary (Fthresh)
            Mtot = mass(out_bin_i)+mass(out_bin_j)
            CM_pos = CoM_2body(pos(out_bin_i,:), mass(out_bin_i), pos(out_bin_j,:), mass(out_bin_j))
            CALL    Calc_binary_info(pos(out_bin_i,:), vel(out_bin_i,:), mass(out_bin_i), & 
                    pos(out_bin_j,:), vel(out_bin_j,:), mass(out_bin_j), binary_info_arr_ij)		! system [i,j]
            E_tot = binary_info_arr_ij(3)
            a_bin = binary_info_arr_ij(4)
            e_bin = binary_info_arr_ij(5)

            ! first, make sure the system is actually bound (PE > KE)
            if (E_tot .LT. 0) then
                ! first unbound object
                r_ijk  = len3vec(CM_pos-pos(ub1,:))
                Ft1 = Fthresh(mass(out_bin_i),mass(out_bin_j),mass(ub1),a_bin,e_bin,r_ijk)
                ! second unbound object
                r_ijl  = len3vec(CM_pos-pos(ub2,:))
                Ft2 = Fthresh(mass(out_bin_i),mass(out_bin_j),mass(ub2),a_bin,e_bin,r_ijl)
                if (Ft1 .LT. delta_F .AND. Ft2 .LT. delta_F)        end_state_flag = 4  !>(N-2) systems are unbound
            endif

        ENDIF!endstate
        endif ! if statement for two unbound objects

        


        ! >2 UNBOUND: write info and set total ionization flag   # FIXME: flagged removed for now
        if (unbound .GT. 2 .AND. n_particles .EQ. 4) then
        IF (end_state_flag .EQ. 0 .AND. n_particles .EQ. 4) THEN		!if endstate not found
            
            ! write info for the output file
            out_bin_i = ub1
            out_bin_j = ub2
            out_bin_k = ub3
            if (unbound .EQ. 4)   out_bin_l = ub4
            mass_bin_i = mass(ub1)
            mass_bin_j = mass(ub2)
            mass_bin_k = mass(ub3)
            if (unbound .EQ. 4)   mass_bin_l = mass(ub4)

!            end_state_flag = 7  !TOTAL IONIZATION

        ENDIF!endstate
        endif ! if statement for total ionization





		!------------------------------------------------------------
		!Check for 2 unbound binary systems moving away from each other (endstate=6)
		!------------------------------------------------------------

        do i=1, n_particles,1
        do j=i+1, n_particles,1
		IF (end_state_flag .EQ. 0 .AND. n_particles .EQ. 4) THEN		!if no end-state yet.
            CALL    Calc_binary_info(pos(i,:), vel(i,:), mass(i), pos(j,:), vel(j,:), mass(j), binary_info_arr_ij) ! system [i,j]
            Etot_ij = binary_info_arr_ij(3)
            ! see if PE is greater than KE
            if (Etot_ij .LT. 0) then
                a_bin = binary_info_arr_ij(4)
                e_bin = binary_info_arr_ij(5)
                do k=1, n_particles,1
                do l=k+1, n_particles,1
                if (k .NE. i .AND. k .NE. j .AND. l .NE. i .AND. l .NE. j) then
                    CALL    Calc_binary_info(pos(k,:), vel(k,:), mass(k), pos(l,:), vel(l,:), mass(l), binary_info_arr_kl) ! system [k,l]
                    Etot_kl = binary_info_arr_kl(3)
                    ! see if PE is greater than KE
                    if (Etot_kl .LT. 0) then
                        a_bin_out = binary_info_arr_kl(4)
                        e_bin_out = binary_info_arr_kl(5)
                        

                        ! now, we see if the two systems are unbound, and if they are below the tidal threshold
                        CM1 = CoM_2body(pos(i,:), mass(i), pos(j,:), mass(j))
                        vCM1 = CoM_2body(vel(i,:), mass(i), vel(j,:), mass(j))
                        M1 = mass(i)+mass(j)
                        CM2 = CoM_2body(pos(k,:), mass(k), pos(l,:), mass(l))
                        vCM2 = CoM_2body(vel(k,:), mass(k), vel(l,:), mass(l))
                        M2 = mass(k)+mass(l)


                        ! see if KE is greater than PE for the two bound systems
                        CALL    Calc_binary_info(CM1, vCM1, M1, CM2, vCM2, M2, binary_info_arr_ijkl)
                        E_tot = binary_info_arr_ijkl(3)

                        !see if systems are moving away from the total CM
                        CM_pos = CoM_4body(pos(i,:), mass(i), pos(j,:), mass(j), pos(k,:), mass(k), pos(l,:), mass(l))
                        CM_vel = CoM_4body(vel(i,:), mass(i), vel(j,:), mass(j), vel(k,:), mass(k), vel(l,:), mass(l))
                        v_dir_ij = DOT_PRODUCT((CM1-CM_pos), (vCM1-CM_vel))
                        v_dir_kl = DOT_PRODUCT((CM2-CM_pos), (vCM2-CM_vel))

                        if (E_tot .GT. 0 .AND. v_dir_ij .GT. 0 .AND. v_dir_kl .GT. 0) then
                            ! see if binary [i,j] is bound below tidal threshold
                            r_ijkl = len3vec(CM1-CM2)
                            Ft1 = Fthresh(mass(i),mass(j),M2,a_bin,e_bin,r_ijkl)
                            ! see if binary [k,l] is bound below tidal threshold    
                            Ft2 = Fthresh(mass(k),mass(l),M1,a_bin_out,e_bin_out,r_ijkl)
                            if (Ft1 .LT. delta_F .AND. Ft2 .LT. delta_F) then
                                end_state_flag = 6  ! two bound systems unbound from one another!

                                ! set the indices of bound systems [i,j] and [k,l]
                                out_bin_i = i
                                mass_bin_i = mass(i)
                                out_bin_j = j
                                mass_bin_j = mass(j)
                                out_bin_k = k
                                mass_bin_k = mass(k)
                                out_bin_l = l
                                mass_bin_l = mass(l)
                            endif
                        endif    !loop for energy and force threshold checks

                    endif   !loop for [k,l] bound
                endif   !loop fo [i,j] bound

                enddo !loop over l
                enddo !loop over k
            endif   !loop for particles being unique
        ENDIF !endstate
        enddo !loop over j
        enddo !loop over i
                        




		!------------------------------------------------------------
		!Check for for energy and angular momentum conservation (endstate = 8) #FIXME Not working, flag removed
		!------------------------------------------------------------
		IF (end_state_flag .EQ. 0) THEN		!if no end-state yet.

        ! calculate initial energy and angular momentum
        E_kin = 0d0
        do i=1, n_particles,1
            E_kin = E_kin + 0.5*mass(i)*len3vec(IC_par_vel(i,:))**2d0
        enddo
        E_pot = 0d0
        do i=1, n_particles,1
            do j=i+1, n_particles,1
                E_pot = E_pot + mass(i)*mass(j)/len3vec(IC_par_pos(i,:)-IC_par_pos(j,:))
            enddo
        enddo
        E_ini = E_kin - E_pot
        L_ini_vec = [0,0,0]
        do i=1, n_particles,1
            do j=i+1, n_particles,1
                CALL Calc_binary_L(IC_par_pos(i,:),IC_par_vel(i,:),mass(i),IC_par_pos(j,:),IC_par_vel(j,:),mass(j),Lvec)
                L_ini_vec = L_ini_vec + Lvec
            enddo
        enddo 
        L_ini = len3vec(L_ini_vec)


        ! calculate the energy and angular momentum at this step
        E_kin = 0d0
        do i=1, n_particles,1
            E_kin = E_kin + 0.5*mass(i)*len3vec(vel(i,:))**2d0
        enddo
        E_pot = 0d0
        do i=1, n_particles,1
            do j=i+1, n_particles,1
                E_pot = E_pot + mass(i)*mass(j)/len3vec(pos(i,:)-pos(j,:))
            enddo
        enddo
        E_check = E_kin - E_pot
        L_check_vec = [0,0,0]
        do i=1, n_particles,1
            do j=i+1, n_particles,1
                CALL Calc_binary_L(pos(i,:),vel(i,:),mass(i),pos(j,:),vel(j,:),mass(j),Lvec)
                L_check_vec = L_check_vec + Lvec
            enddo
        enddo
        L_check = len3vec(L_check_vec)

        if (L_check/L_ini .GT. 1.0-delta_EL .OR. E_check/E_ini .GT. 1.0-delta_EL) then
!            end_state_flag = 8              ! Energy/Momentum criteria conservation not met!        #FIXME
        endif

		ENDIF	!endstate


		!Retun info:
        Return_Nbody_endstate(:) = [0,0,0,0,0]
        Return_endstate_binparams(:) = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        !-----------------------------------
        !if endstate is found
        !-----------------------------------
        if (end_state_flag .NE. 0) then
            out_end_state_flag = end_state_flag
            ! calculate CoM velocity of the bound binary
            vCM = len3vec(CoM_2body(vel(out_bin_i,:), mass(out_bin_i), vel(out_bin_j,:), mass(out_bin_j)))

            Return_Nbody_endstate(:) = [out_end_state_flag, out_bin_i, out_bin_j, out_bin_k, out_bin_l]
            Return_endstate_binparams(:) = [mass_bin_i, mass_bin_j, mass_bin_k, mass_bin_l, a_bin*rsun_to_au, e_bin, &
                            a_bin_out*rsun_to_au, e_bin_out, inc_bin, vCM*vCM_to_kms]
        endif
        
    !------------------------------------------------------------
	END SUBROUTINE Nbody_endstate_sub
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------
	
	
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE  Nbody_info_module_sub(Y, Return_Nbody_info_REAL_1, Return_Nbody_info_REAL_2, Return_Nbody_info_INT_1)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(tot_nr_Y_evol_eqs),	intent(in)		:: Y
		real*8, dimension(10),					intent(out)		:: Return_Nbody_info_REAL_1, Return_Nbody_info_REAL_2
		integer, dimension(10),					intent(out)		:: Return_Nbody_info_INT_1
		real*8, dimension(n_particles,length_Y_per_body_n)		:: body_all_all
		real*8, dimension(n_particles,3)						:: pos, vel
		real*8, dimension(n_particles,3,3)						:: body_all_q, body_all_qdot
		real*8													:: r_ij
		real*8, dimension(5)									:: binary_info_arr_ij
		real*8, dimension(3,2)									:: pair_index_ij
		real*8, dimension(3)									:: Force_arr
		integer													:: bin_i, bin_j
		integer, dimension(1)									:: maxindex
		real*8, dimension(3,3)									:: q_bin_i, q_bin_j, S_bin_i, S_bin_j
		real*8, dimension(3)									:: Eigen_vals_S, a1a2a3_bin_i, a1a2a3_bin_j
		real*8													:: radtidmax_bin_i, radtidmax_bin_j
		real*8, dimension(3)									:: PA_a1a2a3

		!------------------------------------------------------------
		!UNPACK Y:
		!------------------------------------------------------------
		body_all_all(:,:) 			= RESHAPE(Y, (/ n_particles,length_Y_per_body_n /), 	ORDER = (/2,1/))
		pos(:,:)					= body_all_all(:,1:3)
		vel(:,:)					= body_all_all(:,4:6)
		body_all_q(:,:,:)			= RESHAPE(body_all_all(:,7:15),	(/ n_particles,3,3 /), 	ORDER = (/1,2,3/))
		body_all_qdot(:,:,:)		= RESHAPE(body_all_all(:,16:24),(/ n_particles,3,3 /), 	ORDER = (/1,2,3/))
		!------------------------------------------------------------


		!------------------------------------------------------------
		!2-Body INFO:
		!------------------------------------------------------------
		!----------------------------------------------------
		!Find Binary:
		!----------------------------------------------------
		!initialize:
		nc = 1
		Force_arr(:)		= 0
		pair_index_ij(:,:)	= 0
		do i=1,		n_particles,1
		do j=i+1,	n_particles,1
			r_ij				= len3vec(pos(i,:)-pos(j,:))
			Force_arr(nc)		= mass(i)*mass(j)/(r_ij**2)
			pair_index_ij(nc,:)	= [i,j]
			nc = nc+1
		enddo
		enddo
		maxindex	= MAXLOC(Force_arr)	!returns 1D array with arr loc of max element
		!bin-sin candidate indices:
		bin_i 	= pair_index_ij(maxindex(1),1)	!	i
		bin_j	= pair_index_ij(maxindex(1),2)	!	j
		!----------------------------------------------------
		!Calc properties of [i,j] pair:
		!----------------------------------------------------
		!Orbital parameters:
		!return format: [E_kin, E_pot, E_tot, a_bin, e_bin]
		CALL	Calc_binary_info(pos(bin_i,:),	vel(bin_i,:),	mass(bin_i),	pos(bin_j,:),	vel(bin_j,:),	mass(bin_j),	binary_info_arr_ij)		! system [i,j]
		!Distance between bin_i, bin_j:
		r_ij	= len3vec(pos(bin_i,:)-pos(bin_j,:))
		!----------------------------------------------------
		!----------------------------------------------------
		!calc tidal related factors for (binary) pair [i,j]:
		!----------------------------------------------------
		!In this section we calc different properties of the
		!system which is sent back to the main program. We dont
		!define the endstates/usrstates here, but calc the relevant
		!information.
		!-----------------------------------
		!Tidal bulges:
		!-----------------------------------
		!if we evolve tides:
		if (pass_tidalthreshold_yesno .EQ. 1) then
			!Read in q matrix:
			q_bin_i 	= body_all_q(bin_i,:,:)
			q_bin_j 	= body_all_q(bin_j,:,:)
			!calc principal axis (a1,a2,a3), bin_i:
			CALL calc_stellar_paxis_a1a2a3(q_bin_i, PA_a1a2a3)
			radtidmax_bin_i = radius(bin_i)*maxval(PA_a1a2a3) 
			!calc principal axis (a1,a2,a3), bin_j:
			CALL calc_stellar_paxis_a1a2a3(q_bin_j, PA_a1a2a3)
			radtidmax_bin_j = radius(bin_j)*maxval(PA_a1a2a3) 
		endif	!if we evolve tides
		!if we don not evolve tides:
		if (pass_tidalthreshold_yesno .EQ. 0) then
			radtidmax_bin_i = radius(bin_i)
			radtidmax_bin_j = radius(bin_j)	
		endif	!if we do not evolve tides
		!-----------------------------------
		!----------------------------------------------------
		!------------------------------------------------------------
		
		
		!------------------------------------------------------------
		!Return info:
		!------------------------------------------------------------
		!info set 1:
		Return_Nbody_info_REAL_1(1:5)	= [binary_info_arr_ij(:)]
		Return_Nbody_info_INT_1(1:2)	= [bin_i, bin_j]
		!info set 2:
		Return_Nbody_info_REAL_2(1:3)	= [r_ij, radtidmax_bin_i, radtidmax_bin_j]
		!------------------------------------------------------------


	END SUBROUTINE Nbody_info_module_sub
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------
	
	
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE  calc_stellar_paxis_a1a2a3(q_in, a1a2a3_out)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(3,3),	intent(in)		:: q_in
		real*8, dimension(3),	intent(out)		:: a1a2a3_out
		real*8, dimension(3,3)					:: S_mat
		real*8, dimension(3)					:: Eigen_vals_S
		
		!from q calc S matrix:
		S_mat 	= MATMUL(q_in,	TRANSPOSE(q_in))
		!calc principal axis (a1,a2,a3):
		LP_A			= S_mat
		CALL DSYEV( 'V', 'U', 3, LP_A, 3, LP_W, LP_EV_WORK, 8, LP_INFO)
		Eigen_vals_S	= LP_W
		a1a2a3_out	 	= Eigen_vals_S(:)**(1d0/2d0)
		
	END SUBROUTINE calc_stellar_paxis_a1a2a3
	!----------------------------------------------------------------------------------------------------	
	!----------------------------------------------------------------------------------------------------
	
	
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE  TEST_sub(Y, time_t)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(tot_nr_Y_evol_eqs), intent(in)		:: Y
		real*8, intent(in)										:: time_t
		real*8, dimension(n_particles,length_Y_per_body_n)		:: body_all_all
		real*8, dimension(n_particles,3)						:: pos, vel
		real*8, dimension(n_particles,3,3)						:: body_all_q, body_all_qdot
		real*8, dimension(n_particles,3)						:: body_all_a1a2a3
		real*8, dimension(3,3)									:: body_i_q
		real*8, dimension(3,3)									:: body_i_qdot
		real*8, dimension(3,3)									:: body_i_S, body_i_qdqdT
		real*8													:: body_i_S_absdet, body_i_q_absdet
		real*8													:: body_i_OM_sph
		real*8, dimension(3,3)									:: Aint_MAT, A_MAT, A_MAT_body_i_S
		real*8, dimension(3)									:: Eigen_vals_S
		real*8, dimension(n_particles,n_particles)				:: body_ij_E_pot_pmass, body_ij_E_pot_tides
		real*8													:: body_i_Eself_kinetic, body_i_Eself_selfgrav, body_i_Eself_gas
		real*8, dimension(n_particles,3)						:: body_all_Eselfkin_Eselfgrav_Eselfgas
		real*8, dimension(n_particles)							:: body_all_Eselftot
		real*8, dimension(3)									:: vel_CM_vec
		real*8													:: body_i_vel2_wrt_CM
		real*8, dimension(n_particles)							:: body_all_E_kin_wrtCM
		real*8, dimension(3)									:: pos_ij
		real*8													:: r_ij
		real*8													:: E_pot_pmass_ij, E_pot_tides_ij
		real*8, dimension(1,3)									:: X_MAT
		real*8, dimension(3,3)									:: XX_MAT, C_ij, C_ij_body_i_S
		real*8													:: E_tot_kin, E_tot_pot_pmass, E_tot_pot_tides, E_tot_self
		real*8													:: E_tot_external_terms, E_tot_internal_terms
		real*8													:: E_tot_system
		integer													:: bin_i, bin_j, sin_k
		real*8													:: a_bin_ij, e_bin_ij, rperi_ij
		integer, dimension(1)									:: maxindex
		real*8													:: acc_tid, acc_bind, T_ij, T_ji
		real*8, dimension(n_particles*n_particles)				:: Force_arr
		real*8, dimension(n_particles*n_particles,2)			:: pair_index_ij
		real*8													:: CM_Mtot_ij
		real*8, dimension(3)									:: CM_ij_pos, CM_ij_vel
		real*8, dimension(5)									:: binary_info_arr_ij, binary_info_arr_ijk
		real*8, dimension(3)									:: PA_a1a2a3
        real*8                                                  :: timesamp

        !------------------------------------------------------------
        !Set downsampling conditions
        !------------------------------------------------------------
        timesamp = INT(max_sim_time)/downsample    !this currently downsamples each output file 100 times
		
		!------------------------------------------------------------
		!UNPACK Y:
		!------------------------------------------------------------
		body_all_all(:,:) 			= RESHAPE(Y, (/ n_particles,length_Y_per_body_n /), 	ORDER = (/2,1/))
		pos(:,:)					= body_all_all(:,1:3)
		vel(:,:)					= body_all_all(:,4:6)
		body_all_q(:,:,:)			= RESHAPE(body_all_all(:,7:15),	(/ n_particles,3,3 /), 	ORDER = (/1,2,3/))
		body_all_qdot(:,:,:)		= RESHAPE(body_all_all(:,16:24),(/ n_particles,3,3 /), 	ORDER = (/1,2,3/))
		!------------------------------------------------------------


		!------------------------------------------------------------
		!Initialize:
		!------------------------------------------------------------
		body_ij_E_pot_pmass(:,:)	= 0d0
		body_ij_E_pot_tides(:,:)	= 0d0
		!------------------------------------------------------------


		!------------------------------------------------------------
		!Calc:
		!------------------------------------------------------------
		!CM velocity:
		do k=1, 3,1
		vel_CM_vec(k) = SUM(mass(:)*vel(:,k))/SUM(mass(:))
		enddo
		!------------------------------------------------------------


		!------------------------------------------------------------
		!Loop over objects:
		!------------------------------------------------------------
		do i=1, n_particles,1
	
			!--------------------------------------------------------
			!Define: 
			!--------------------------------------------------------
			body_i_q 		= body_all_q(i,:,:)
			body_i_qdot		= body_all_qdot(i,:,:)
			body_i_S 		= MATMUL(body_i_q,	  TRANSPOSE(body_i_q))
			body_i_qdqdT	= MATMUL(body_i_qdot, TRANSPOSE(body_i_qdot))
			body_i_S_absdet = ABS(func_det_3x3M(body_i_S))
			body_i_q_absdet	= ABS(func_det_3x3M(body_i_q))	
			body_i_OM_sph	= - (3d0/(5d0-gas_n(i)))*((mass(i)**2.)/radius(i))
			!--------------------------------------------------------
	
			!--------------------------------------------------------
			!Calc Principal axis: a1,a2,a3:
			!--------------------------------------------------------
			CALL calc_stellar_paxis_a1a2a3(body_i_q, PA_a1a2a3)
			body_all_a1a2a3(i,:)	= PA_a1a2a3
			!--------------------------------------------------------
	
			!--------------------------------------------------------
			!Calc: Internal terms - Stellar self energies:
			!--------------------------------------------------------
		    !Kinetic energy:
			body_i_Eself_kinetic	= (1d0/2d0)*Mqp_sph(i)*Trace_3x3M(body_i_qdqdT)
		    !Gravitational self energy:
			CALL Calc_Aint_MAT_selfgrav(body_i_S, Aint_MAT)	!in: S (=qqT), out: integral appearing in A (A = sqrt(det(S))*Aint)
			A_MAT					= (body_i_S_absdet**(1d0/2d0))*Aint_MAT
			A_MAT_body_i_S			= MATMUL(A_MAT,body_i_S)
		    body_i_Eself_selfgrav	= (1d0/2d0)*body_i_OM_sph*(body_i_S_absdet**(-1d0/2d0))*Trace_3x3M(A_MAT_body_i_S)
		    !Gas energy
		    body_i_Eself_gas		= - (body_i_OM_sph/(3d0*(gas_gamma(i)-1d0)))*(body_i_q_absdet**(1d0-gas_gamma(i)))
			!save:
			body_all_Eselfkin_Eselfgrav_Eselfgas(i,:)	= [body_i_Eself_kinetic,  body_i_Eself_selfgrav,  body_i_Eself_gas]
			body_all_Eselftot(i)						=  body_i_Eself_kinetic + body_i_Eself_selfgrav + body_i_Eself_gas
			!--------------------------------------------------------
	
			!--------------------------------------------------------
			!Calc: External terms - Orbital energies:
			!--------------------------------------------------------
			!E_kin wrt CM:
			body_i_vel2_wrt_CM		= SUM((vel(i,:) - vel_CM_vec(:))**2)
			body_all_E_kin_wrtCM(i)	= (1d0/2d0)*mass(i)*body_i_vel2_wrt_CM	
			!E_pot terms:
			do j=1, n_particles,1
			if (i .NE. j) then
				!define:
			    pos_ij	= pos(i,:)-pos(j,:)
				r_ij	= len3vec(pos_ij)
				!Epot - Point Particle limit:
				E_pot_pmass_ij	= - mass(i)*mass(j)/r_ij	
				body_ij_E_pot_pmass(i,j)	= E_pot_pmass_ij
				!Epot - Tidal Interaction:
			    X_MAT(1,:)		= pos_ij	!X_MAT is here made into 2d (but its really just a 1d vec).
			    XX_MAT			= MATMUL(TRANSPOSE(X_MAT), X_MAT)				
				C_ij        	= (1d0/(r_ij**3))*(3d0*XX_MAT/(r_ij**2) - I_MAT)   !Tidal Tensor (Cij)
				C_ij_body_i_S 	= MATMUL(C_ij,body_i_S) 
				E_pot_tides_ij	= - (1d0/2d0)*mass(j)*Mqp_sph(i)*Trace_3x3M(C_ij_body_i_S)	
				body_ij_E_pot_tides(i,j)	= E_pot_tides_ij
			endif
			enddo	
			!--------------------------------------------------------	
		
		enddo	!loop obj i
		!------------------------------------------------------------


		!------------------------------------------------------------
		!All calc energies below are the total of the full NBsystem in the CM frame.
		!------------------------------------------------------------
		E_tot_kin		= SUM(body_all_E_kin_wrtCM(:))
		E_tot_pot_pmass	= SUM(body_ij_E_pot_pmass(:,:))/2d0
		E_tot_pot_tides = SUM(body_ij_E_pot_tides(:,:))
		E_tot_self		= SUM(body_all_Eselftot(:))

		E_tot_external_terms	= E_tot_kin + E_tot_pot_pmass + E_tot_pot_tides
		E_tot_internal_terms	= E_tot_self
		E_tot_system			= E_tot_external_terms + E_tot_internal_terms
		!------------------------------------------------------------


		!------------------------------------------------------------
		!Calc orbital params for most bound (or only) (bound/unbound pair) Binary:
		!------------------------------------------------------------
		!----------------------------------------------------
		!Find Binary:
		!----------------------------------------------------
		!initialize:
		nc = 1
		Force_arr(:)		= 0
		pair_index_ij(:,:)	= 0
		do i=1,		n_particles,1
		do j=i+1,	n_particles,1
			r_ij				= len3vec(pos(i,:)-pos(j,:))
			Force_arr(nc)		= mass(i)*mass(j)/(r_ij**2)
			pair_index_ij(nc,:)	= [i,j]
			nc = nc+1
		enddo
		enddo
		maxindex	= MAXLOC(Force_arr)	!returns 1D array with arr loc of max element
		!bin-sin candidate indices:
		bin_i 	= pair_index_ij(maxindex(1),1)	!	i
		bin_j	= pair_index_ij(maxindex(1),2)	!	j
		!----------------------------------------------------
		!Calc orbital params for i,j pair:
		!----------------------------------------------------
		![E_kin, E_pot, E_tot, a_bin, e_bin]
		CALL	Calc_binary_info(pos(bin_i,:),	vel(bin_i,:),	mass(bin_i),	pos(bin_j,:),	vel(bin_j,:),	mass(bin_j),	binary_info_arr_ij)		! system [i,j]
		a_bin_ij	= binary_info_arr_ij(4)
		e_bin_ij	= binary_info_arr_ij(5)
		rperi_ij	= a_bin_ij*(1.-e_bin_ij)
		!----------------------------------------------------
		!Calc tidal forces:
		!----------------------------------------------------
		!dist between bin_i, bin_j:
		r_ij		= len3vec(pos(bin_i,:)-pos(bin_j,:))
		!F_tide/F_binding of obj_i caused by obj_j:
		acc_tid		= mass(bin_j)*(1./((r_ij-radius(bin_i))**2.) - 1./((r_ij+radius(bin_i))**2.))
		acc_bind	= mass(bin_i)/(radius(bin_i)**2.)
		T_ij		= ABS(acc_tid/acc_bind)		
		!F_tide/F_binding of obj_j caused by obj_i:
		acc_tid		= mass(bin_i)*(1./((r_ij-radius(bin_j))**2.) - 1./((r_ij+radius(bin_j))**2.))
		acc_bind	= mass(bin_j)/(radius(bin_j)**2.)
		T_ji		= ABS(acc_tid/acc_bind)		
		!----------------------------------------------------
		!------------------------------------------------------------
		
		
		!------------------------------------------------------------
		!3-body: Calc orbital params for both bin, bin-sin:
		!------------------------------------------------------------
		if (n_particles .EQ. 3) then
			!We use the bin_i,bin_j from above.
			sin_k	= 6 - (bin_i+bin_j)				!	k
			!calc CMij properties:
			CM_Mtot_ij	= (mass(bin_i) + mass(bin_j))
			CM_ij_pos	= (mass(bin_i)*pos(bin_i,:) + mass(bin_j)*pos(bin_j,:))/CM_Mtot_ij
			CM_ij_vel	= (mass(bin_i)*vel(bin_i,:) + mass(bin_j)*vel(bin_j,:))/CM_Mtot_ij
			CALL	Calc_binary_info(pos(bin_i,:),	vel(bin_i,:),	mass(bin_i),	pos(bin_j,:),	vel(bin_j,:),	mass(bin_j),	binary_info_arr_ij)		! system [i,j]
			CALL	Calc_binary_info(pos(sin_k,:),	vel(sin_k,:),	mass(sin_k),	CM_ij_pos,		CM_ij_vel,		CM_Mtot_ij,		binary_info_arr_ijk)	! system [[i,j],k]
			!out_binary_info_arr(:) = [E_kin, E_pot, E_tot, a_bin, e_bin]
		endif
		!------------------------------------------------------------

		
		!------------------------------------------------------------
		!Write info to file:
		!------------------------------------------------------------
		!for 4 particle
		if (n_particles .EQ. 4) then
            if (downsample .EQ. 0) then
			    write(10,*) pos(1,:),  pos(2,:), pos(3,:), pos(4,:)
			    write(18,*) vel(1,:),  vel(2,:), vel(3,:), vel(4,:)
			    write(11,*) time_t
            endif
            if (downsample .NE. 1 .AND. time_t .GE. ds*timesamp) then
			    write(10,*) pos(1,:),  pos(2,:), pos(3,:), pos(4,:)
			    write(18,*) vel(1,:),  vel(2,:), vel(3,:), vel(4,:)
			    write(11,*) time_t
                ds = ds + 1
		    endif
        endif
		!for 3 particles
		if (n_particles .EQ. 3) then
			write(10,*) pos(1,:),  pos(2,:), pos(3,:)
			write(18,*) vel(1,:),  vel(2,:), vel(3,:)
			write(11,*) time_t, body_all_a1a2a3(1,:),  body_all_a1a2a3(2,:), body_all_a1a2a3(3,:)
			write(12,*) time_t, body_all_Eselftot(1),  body_all_Eselftot(2), body_all_Eselftot(3)
			write(13,*) time_t, E_tot_kin, E_tot_pot_pmass, E_tot_pot_tides, &
						E_tot_internal_terms, E_tot_external_terms, E_tot_system
			write(14,*) time_t, bin_i, bin_j, a_bin_ij, e_bin_ij, rperi_ij, T_ij, T_ji
			write(15,*) body_all_q(1,:,:) 
			write(15,*) body_all_q(2,:,:) 
			write(15,*) body_all_q(3,:,:) 
			write(21,*) time_t, binary_info_arr_ij(:), binary_info_arr_ijk(:)
		endif
		!for 2 particles
		if (n_particles .EQ. 2) then
			write(10,*) pos(1,:),  pos(2,:)
			write(18,*) vel(1,:),  vel(2,:)
			write(11,*) time_t, body_all_a1a2a3(1,:),  body_all_a1a2a3(2,:)
			write(12,*) time_t, body_all_Eselftot(1),  body_all_Eselftot(2)
			write(13,*) time_t, E_tot_kin, E_tot_pot_pmass, E_tot_pot_tides, &
			 			E_tot_internal_terms, E_tot_external_terms, E_tot_system
			write(14,*) time_t, bin_i, bin_j, a_bin_ij, e_bin_ij, rperi_ij, T_ij, T_ji
		endif
		!------------------------------------------------------------


	END SUBROUTINE TEST_sub
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------



	!----------------------------------------------------------------------------------------------------
	SUBROUTINE  output_Nbody_sub(Y, time_t)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		real*8, dimension(tot_nr_Y_evol_eqs), intent(in)		:: Y
		real*8, intent(in)										:: time_t
		real*8, dimension(n_particles,length_Y_per_body_n)		:: body_all_all
		real*8, dimension(n_particles,3)						:: pos, vel
		real*8, dimension(n_particles,3,3)						:: body_all_q, body_all_qdot

		!------------------------------------------------------------
		!UNPACK Y:
		!------------------------------------------------------------
		body_all_all(:,:) 			= RESHAPE(Y, (/ n_particles,length_Y_per_body_n /), 	ORDER = (/2,1/))
		pos(:,:)					= body_all_all(:,1:3)
		vel(:,:)					= body_all_all(:,4:6)
		body_all_q(:,:,:)			= RESHAPE(body_all_all(:,7:15),	(/ n_particles,3,3 /), 	ORDER = (/1,2,3/))
		body_all_qdot(:,:,:)		= RESHAPE(body_all_all(:,16:24),(/ n_particles,3,3 /), 	ORDER = (/1,2,3/))
		!------------------------------------------------------------

		
		!------------------------------------------------------------
		!Write info to file:
		!------------------------------------------------------------
		!for 2 particle
		if (n_particles .EQ. 2) then
            !if downsample is 0, don't write any position/velocity data
            !if downsample is 1, write every step
            if (downsample .EQ. 1) then
			    write(10,*) pos(1,:), pos(2,:)
			    write(18,*) vel(1,:), vel(2,:)
			    write(11,*) time_t
            endif
            !if downsample is an integer >1, only write the downsample'th step
            if (downsample .GT. 1 .AND. stepc .GE. ds*downsample) then
			    write(10,*) pos(1,:), pos(2,:)
			    write(18,*) vel(1,:), vel(2,:)
			    write(11,*) time_t
                ds = ds + 1
		    endif
        endif

		!for 3 particle
		if (n_particles .EQ. 3) then
            !if downsample is 0, don't write any position/velocity data
            !if downsample is 1, write every step
            if (downsample .EQ. 1) then
			    write(10,*) pos(1,:), pos(2,:), pos(3,:)
			    write(18,*) vel(1,:), vel(2,:), vel(3,:)
			    write(11,*) time_t
            endif
            !if downsample is an integer >1, only write the downsample'th step
            if (downsample .GT. 1 .AND. stepc .GE. ds*downsample) then
			    write(10,*) pos(1,:), pos(2,:), pos(3,:)
			    write(18,*) vel(1,:), vel(2,:), vel(3,:)
			    write(11,*) time_t
                ds = ds + 1
		    endif
        endif

		!for 4 particle
		if (n_particles .EQ. 4) then
            !if downsample is 0, don't write any position/velocity data
            !if downsample is 1, write every step
            if (downsample .EQ. 1) then
			    write(10,*) pos(1,:), pos(2,:), pos(3,:), pos(4,:)
			    write(18,*) vel(1,:), vel(2,:), vel(3,:), vel(4,:)
			    write(11,*) time_t
            endif
            !if downsample is an integer >1, only write the downsample'th step
            if (downsample .GT. 1 .AND. stepc .GE. ds*downsample) then
			    write(10,*) pos(1,:), pos(2,:), pos(3,:), pos(4,:)
			    write(18,*) vel(1,:), vel(2,:), vel(3,:), vel(4,:)
			    write(11,*) time_t
                ds = ds + 1
		    endif
        endif
		!------------------------------------------------------------


	END SUBROUTINE output_Nbody_sub
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------
	
	
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE  F_Y_Ydot(NEQ, T, Y, YDOT)
	!----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		!in/out params for the subroutine:
		integer												:: NEQ		! = tot_nr_Y_evol_eqs
		real*8												:: T
		real*8, dimension(tot_nr_Y_evol_eqs)				:: Y, YDOT
		!for unpack:
		real*8, dimension(n_particles,length_Y_per_body_n)	:: Y_in_body_all_all, Ydot_out_body_all_all
		real*8, dimension(n_particles,3)					:: Y_in_body_all_pos, Y_in_body_all_vel
		real*8, dimension(n_particles,3,3)					:: Y_in_body_all_q, Y_in_body_all_qdot
		real*8, dimension(3)								:: acc_CM_pointmass_NT_ij, acc_CM_pointmass_1PN_ij, acc_CM_pointmass_2PN_ij, acc_CM_pointmass_25PN_ij
		real*8, dimension(3)								:: acc_CM_tidalfield_ij, body_i_acc_CM_sum_j
		real*8, dimension(3,3,3)							:: D_MAT_xk
		real*8, dimension(3)								:: body_i_pos, body_i_vel
		real*8, dimension(3,3)								:: body_i_q, body_i_qdot, body_i_q_T, body_i_S
		real*8, dimension(3,3)								:: body_i_q_T_inv
		real*8, dimension(3)								:: body_j_pos, body_j_vel
		real*8, dimension(3,3)								:: body_j_q, body_j_q_T, body_j_S, M_MAT
		real*8, dimension(3)								:: pos_ij, vel_ij
		real*8												:: r_ij, v_ij
		real*8, dimension(1,3)								:: X_MAT
		real*8, dimension(3,3)								:: XX_MAT, Cp_MAT_xk, Cp_xk_M_MAT
		real*8												:: xk, Cp_M_sum_xk
		real*8, dimension(3)								:: acc_tot_CM_ij
		real*8, dimension(3,3)								:: Tidal_Tensor_C_body_ij, acc_Tidal_Tensor_C_body_ij
		real*8, dimension(3,3)								:: body_i_acc_tidal_tensor_C_sum_j
		real*8												:: q_i_absdet, S_i_absdet, selfgrav_OM_sph
		real*8, dimension(3,3)								:: selfgrav_OM_MAT
		real*8												:: pressure_int_PI_sph, pressure_int_PI
		real*8, dimension(3,3)								:: body_i_qdotdot
		real*8, dimension(3,3)								:: A_MAT, Aint_MAT
		real*8, dimension(3)								:: n_ij
		real*8												:: n_dot_vi, n_dot_vj, vel_ij2
        real*8, dimension(3)                                :: v_CM, vi_CM, vj_CM
        real*8                                              :: vii_CM, vjj_CM, vij_CM
        real*8                                              :: n_dot_vi_CM, n_dot_vj_CM


		!------------------------------------------------------------
		!UNPACK:
		!------------------------------------------------------------
		Y_in_body_all_all(:,:) 			= RESHAPE(Y, (/ n_particles,length_Y_per_body_n /), ORDER = (/2,1/))
		Y_in_body_all_pos(:,:)			= Y_in_body_all_all(:,1:3)
		Y_in_body_all_vel(:,:)			= Y_in_body_all_all(:,4:6)
		Y_in_body_all_q(:,:,:)			= RESHAPE(Y_in_body_all_all(:,7:15),	(/ n_particles,3,3 /), ORDER = (/1,2,3/))
		Y_in_body_all_qdot(:,:,:)		= RESHAPE(Y_in_body_all_all(:,16:24),	(/ n_particles,3,3 /), ORDER = (/1,2,3/))
		!------------------------------------------------------------


		!------------------------------------------------------------
		!CALC EOM for each particle i (tidal modes and orbit):
		!------------------------------------------------------------
		do i=1, n_particles,1
		!------------------------------------------------------------	

			!--------------------------------------------------------
		    !INFO body i:         
		    !--------------------------------------------------------
		    body_i_pos       	= Y_in_body_all_pos(i,:)
		    body_i_vel       	= Y_in_body_all_vel(i,:)
		    body_i_q     		= Y_in_body_all_q(i,:,:)
		    body_i_qdot  		= Y_in_body_all_qdot(i,:,:)         
		    body_i_q_T   		= TRANSPOSE(body_i_q)
		    body_i_S     		= MATMUL(body_i_q, body_i_q_T)
		    !--------------------------------------------------------

		    !--------------------------------------------------------
		    !EOM for Center-of-Mass (CM) of body i:
		    !--------------------------------------------------------
			!initialize:
		    body_i_acc_CM_sum_j(:)  = 0d0
			!--------------------------------------------------------
			!loop over body j:
			!--------------------------------------------------------
			do j=1, n_particles,1
			if (i .NE. j) then			
			    !-----------------------------------
			    !INFO body j:         
			    !-----------------------------------
			    body_j_pos       	= Y_in_body_all_pos(j,:)  
			    body_j_vel       	= Y_in_body_all_vel(j,:)
			    body_j_q     		= Y_in_body_all_q(j,:,:)
			    body_j_q_T   		= TRANSPOSE(body_j_q)
			    body_j_S     		= MATMUL(body_j_q, body_j_q_T)
			    !-----------------------------------
			    !Define:
			    !-----------------------------------
			    pos_ij  = (body_i_pos - body_j_pos) 
			    vel_ij  = (body_i_vel - body_j_vel)
			    r_ij	= len3vec(pos_ij)
			    v_ij	= len3vec(vel_ij)
			    !-----------------------------------                   
			    !Point mass contributions:
				!-----------------------------------
				!define:
				n_ij		= pos_ij/r_ij					!unit vec pointing from j to i
				n_dot_vi	= dot_product(n_ij,body_i_vel)
				n_dot_vj	= dot_product(n_ij,body_j_vel)
				vel_ij2		= dot_product(vel_ij,vel_ij)	
                v_CM        = ((mass(i)*body_i_vel + mass(j)*body_j_vel) / (mass(i) + mass(j)))     ! CM velocity
                vi_CM       = (v_CM - body_i_vel)            ! vi relative to CM
                vj_CM       = (v_CM - body_j_vel)            ! vj relative to CM
                vii_CM      = dot_product(vi_CM,vi_CM)     ! vi^2
                vjj_CM      = dot_product(vj_CM,vj_CM)     ! vj^2
                vij_CM      = dot_product(vi_CM,vj_CM)     ! vi dot vj
                n_dot_vi_CM = dot_product(n_ij,vi_CM)  ! n_ij dot vi_CM
                n_dot_vj_CM = dot_product(n_ij,vj_CM)  ! n_ij dot vj_CM
				!acc term: Newtonian
			    acc_CM_pointmass_NT_ij(:)						 = - (mass(j)/(r_ij**2.))*n_ij
				!acc term: 1PN, 2PN
                !acc term: 1PN
                if (use_1PN .EQ. 1) acc_CM_pointmass_1PN_ij(:)   =  (PN_gamma)*(mass(j)/(r_ij**2d0)) &
*(n_ij*(-vii_CM-2d0*vjj_CM+4d0*vij_CM+(3d0/2d0)*(n_dot_vj_CM)**(2d0)+5d0*(mass(i)/r_ij)+4d0*(mass(j)/r_ij)) &
+(vi_CM-vj_CM)*(4d0*n_dot_vi_CM-3d0*n_dot_vj_CM))
                if (use_1PN .EQ. 0) acc_CM_pointmass_1PN_ij(:)  =   0d0
                !acc term: 2PN
                if (use_2PN .EQ. 1) acc_CM_pointmass_2PN_ij(:)   =  ((PN_gamma)**(2d0))*(mass(j)/(r_ij**2d0)) &
* (n_ij * (-2d0*vjj_CM**(2d0) + 4d0*vjj_CM*vij_CM - 2d0*vij_CM**(2d0) + (3d0/2d0)*vii_CM*n_dot_vj_CM**(2d0) &
+ (9d0/2d0)*vjj_CM*n_dot_vj_CM**(2d0) - 6d0*vij_CM*n_dot_vj_CM**(2d0) - (15d0/8d0)*n_dot_vj_CM**(4d0) &
+ (mass(i)/r_ij) * (-(15d0/4d0)*vii_CM + (5d0/4d0)*vjj_CM - (5d0/2d0)*vij_CM + (39d0/2d0)*n_dot_vi_CM**(2d0) &
- 39d0*n_dot_vi_CM*n_dot_vj_CM + (17d0/2d0)*n_dot_vj_CM**(2d0)) + (mass(j)/r_ij) * (4d0*vjj_CM - 8d0*vij_CM &
+ 2d0*n_dot_vi_CM**(2d0) - 4d0*n_dot_vi_CM*n_dot_vj_CM - 6d0*n_dot_vj_CM**2)) &
+ (vi_CM - vj_CM) * (vii_CM*n_dot_vj_CM + 4d0*vjj_CM*n_dot_vi_CM - 5d0*vjj_CM*n_dot_vj_CM &
- 4d0*vij_CM*n_dot_vi_CM + 4d0*vij_CM*n_dot_vj_CM - 6d0*n_dot_vi_CM*n_dot_vj_CM**(2d0) &
+ (9d0/2d0)*n_dot_vj_CM**(3d0) + (mass(i)/r_ij) * ((-63d0/4d0)*n_dot_vi_CM + (55d0/4d0)*n_dot_vj_CM) &
+ (mass(j)/r_ij)*(-2d0*n_dot_vi_CM - 2d0*n_dot_vj_CM))) + ((PN_gamma)**(2d0)) * (mass(j)/(r_ij**(4d0))) &
* n_ij * (-(57d0/4d0)*mass(i)**(2d0) - 9d0*mass(j)**2d0 - (69d0/2d0)*mass(i)*mass(j))
                if (use_2PN .EQ. 0) acc_CM_pointmass_2PN_ij(:)   =   0d0
				!acc term: 2.5PN
			    if (use_25PN .EQ. 1) acc_CM_pointmass_25PN_ij(:) =	(PN_gamma**(5d0/2d0))*(4d0/5d0)*(mass(i)*mass(j)/(r_ij**3d0))		&
																	* (vel_ij*(-vel_ij2 + 2d0*(mass(i)/r_ij) - 8d0*(mass(j)/r_ij))		&
																	+ n_ij*(n_dot_vi-n_dot_vj)*(3d0*vel_ij2 - 6d0*(mass(i)/r_ij) + (52d0/3d0)*(mass(j)/r_ij))) 
				if (use_25PN .EQ. 0) acc_CM_pointmass_25PN_ij(:) =	0d0
				!-----------------------------------
				!Tidal field contribution:
				!-----------------------------------
				!If we are not evolving tides:
			    acc_CM_tidalfield_ij(:) = 0d0
				!If we are evolving tides:
			    if (pass_tidalthreshold_yesno .EQ. 1) then
				if (evoTides_yesno(i) .EQ. 1 .or. evoTides_yesno(j) .EQ. 1) then   !we evolve tides for at least one of the two objs i,j
				    !D_MAT_xk = d(xx)/dx_k:
				    D_MAT_xk(1,:,:) = RESHAPE([2*pos_ij(1),pos_ij(2),pos_ij(3),pos_ij(2),0d0,0d0,pos_ij(3),0d0,0d0], (/3,3/))
				    D_MAT_xk(2,:,:) = RESHAPE([0d0,pos_ij(1),0d0,pos_ij(1),2*pos_ij(2),pos_ij(3),0d0,pos_ij(3),0d0], (/3,3/))
				    D_MAT_xk(3,:,:) = RESHAPE([0d0,0d0,pos_ij(1),0d0,0d0,pos_ij(2),pos_ij(1),pos_ij(2),2*pos_ij(3)], (/3,3/))
				    !M_MAT = (m_1M_2S_2 + m_2M_1S_1):
				    M_MAT	= (mass(i)*Mqp_sph(j)*body_j_S + mass(j)*Mqp_sph(i)*body_i_S)                    
				    !XX = (xx^T):
				    X_MAT(1,:)	= pos_ij	!X_MAT is here made into 2d (but its really just a 1d vec).
				    XX_MAT		= MATMUL(TRANSPOSE(X_MAT), X_MAT)
					!sum over each pos coordinate (x,y,z):
					do k=1,3,1
				    	xk          = pos_ij(k)
				    	Cp_MAT_xk   = 3d0*(1d0/(r_ij**5))*(D_MAT_xk(k,:,:) - 5d0*xk*(XX_MAT/(r_ij**2)) + xk*I_MAT)		!Tidal Deviation Tensor (Cijk)
				    	Cp_xk_M_MAT	= MATMUL(Cp_MAT_xk, M_MAT)
						Cp_M_sum_xk	= Trace_3x3M(Cp_xk_M_MAT)
				    	!final CM-acc from tides:  
				    	acc_CM_tidalfield_ij(k) = (1d0/mass(i))*(1d0/2d0)*Cp_M_sum_xk 
					enddo
				endif	!if we follow tides on obj i,j
				endif	!if system has passed tidal threshold
			    !-----------------------------------                               
			    !Calc total CM acc and sum:
			    !-----------------------------------
			    !Total contribution to CM acc:
				acc_tot_CM_ij		= acc_CM_pointmass_NT_ij + acc_CM_pointmass_1PN_ij + acc_CM_pointmass_2PN_ij &
                                    + acc_CM_pointmass_25PN_ij + acc_CM_tidalfield_ij
				!save and sum terms:
			    body_i_acc_CM_sum_j	= body_i_acc_CM_sum_j + acc_tot_CM_ij 
			    !-----------------------------------                                        
			endif !(i .NE. j)
			enddo !end loop over body j.
			!--------------------------------------------------------
			!--------------------------------------------------------

		    !--------------------------------------------------------
		    !EOM for TIDES for 'body i':
		    !--------------------------------------------------------
			!-----------------------------------
			!If we do not evolve tides:
			!-----------------------------------
			body_i_qdotdot(:,:) 	= 0d0
			!-----------------------------------
			!-----------------------------------
			!-----------------------------------
			!If we evolve tides:
		    !-----------------------------------
			if (pass_tidalthreshold_yesno .EQ. 1) then
			if (evoTides_yesno(i) .EQ. 1) then               
		    !-----------------------------------
		    !External terms
		    !-----------------------------------
		    !initialize:
		    body_i_acc_tidal_tensor_C_sum_j(:,:)    = 0d0                
		    !loop over body j:
			do j=1, n_particles,1
			if (i .NE. j) then
		        !INFO body j:         
		        body_j_pos	= Y_in_body_all_pos(j,:)
		        !Define:
		        pos_ij		= (body_i_pos - body_j_pos) 
			    r_ij		= len3vec(pos_ij)
		        !calc:
			    X_MAT(1,:)	= pos_ij	!X_MAT is here made into 2d (but its really just a 1d vec).
			    XX_MAT		= MATMUL(TRANSPOSE(X_MAT), X_MAT)				
		        Tidal_Tensor_C_body_ij                  = (1d0/(r_ij**3))*(3d0*XX_MAT/(r_ij**2) - I_MAT)   !Tidal Tensor (Cij)
		        acc_Tidal_Tensor_C_body_ij              = mass(j)*Tidal_Tensor_C_body_ij
		        body_i_acc_tidal_tensor_C_sum_j(:,:)    = body_i_acc_tidal_tensor_C_sum_j(:,:) + acc_Tidal_Tensor_C_body_ij(:,:) 
			endif
			enddo !end loop over body j.
			!-----------------------------------                    
		    !Internal terms
		    !-----------------------------------         
			!q_T^-1:
			CALL Calc_Matrix_Inverse(body_i_q_T, body_i_q_T_inv)
			!det of q and S:
			q_i_absdet	= ABS(func_det_3x3M(body_i_q))	
			S_i_absdet	= ABS(func_det_3x3M(body_i_S))		
		    !calc self gravity tensor:
		    selfgrav_OM_sph	= - (3d0/(5d0-gas_n(i)))*((mass(i)**2)/radius(i))	
			CALL Calc_Aint_MAT_selfgrav(body_i_S, Aint_MAT)	!in: S (=qqT), out: integral appearing in A (A = sqrt(det(S))*Aint)
			A_MAT			= (S_i_absdet**(1d0/2d0))*Aint_MAT
		    selfgrav_OM_MAT = (1d0/2d0)*selfgrav_OM_sph*(S_i_absdet**(-1./2.))*MATMUL(A_MAT,body_i_S)
		    !calc pressure integral:
		    pressure_int_PI_sph = - selfgrav_OM_sph/3d0 !from virial theorem.
		    pressure_int_PI     = pressure_int_PI_sph*(q_i_absdet**(1d0-gas_gamma(i))) 	
			!-----------------------------------
		    !FINAL qdotdot: EOM for q_i:
		    !-----------------------------------
		    body_i_qdotdot  =	  (1d0/Mqp_sph(i))*(pressure_int_PI)*body_i_q_T_inv			&
								+ (1d0/Mqp_sph(i))*MATMUL(selfgrav_OM_MAT,body_i_q_T_inv)	&
								+ MATMUL(body_i_acc_tidal_tensor_C_sum_j,body_i_q)
			!-----------------------------------           	
			endif	!if we follow tides on obj i
			endif	!if system has passed tidal threshold
			!-----------------------------------
			!-----------------------------------	
			!--------------------------------------------------------

		    !--------------------------------------------------------
		    !PACK and SAVE Ydot OUT:
		    !--------------------------------------------------------
		    Ydot_out_body_all_all(i,1:3)    = body_i_vel                		!xdot
		    Ydot_out_body_all_all(i,4:6)    = body_i_acc_CM_sum_j   			!xdotdot
		    Ydot_out_body_all_all(i,7:15)   = RESHAPE(body_i_qdot, (/9/))		!qdot
		    Ydot_out_body_all_all(i,16:24)  = RESHAPE(body_i_qdotdot, (/9/))	!qdotdot
		    !--------------------------------------------------------			
		
		enddo !end loop over body i. 
		!------------------------------------------------------------


		!------------------------------------------------------------
		!make Ydot_out vec and RETURN:
		!------------------------------------------------------------
		YDOT(:) = RESHAPE(TRANSPOSE(Ydot_out_body_all_all), (/tot_nr_Y_evol_eqs/))
		!------------------------------------------------------------


	END SUBROUTINE F_Y_Ydot
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------
	
	
	!----------------------------------------------------------------------------------------------------
	SUBROUTINE main_evolve_Nbody_system_sub(IC_code_version, IC_simparams_INT, IC_simparams_REAL,						&
											dimlen_IC_nbody, IC_nbody_const_posvel_qqdot_etc_arr,						&
											out_endstate_INT, out_endstate_REAL, out_usrstate_INT, out_usrstate_REAL,	&
											out_xtra_info_INT, out_xtra_info_REAL, out_xtra_2_info_REAL)
	!----------------------------------------------------------------------------------------------------
	IMPLICIT NONE
	!Input:
	integer,								intent(in)	:: IC_code_version
	integer, dimension(10),					intent(in)	:: IC_simparams_INT
	real*8, dimension(10),					intent(in)	:: IC_simparams_REAL
	integer,								intent(in)	:: dimlen_IC_nbody
	real*8, dimension(dimlen_IC_nbody,10),	intent(in)	:: IC_nbody_const_posvel_qqdot_etc_arr
	!ouput:	
	integer, dimension(10),					intent(out)	:: out_endstate_INT
	real*8, dimension(10),					intent(out)	:: out_endstate_REAL
	integer, dimension(10,10),				intent(out)	:: out_usrstate_INT
	real*8, dimension(10,10),				intent(out)	:: out_usrstate_REAL
	integer, dimension(10),					intent(out)	:: out_xtra_info_INT
	real*8, dimension(10),					intent(out)	:: out_xtra_info_REAL
	real*8, dimension(10),					intent(out)	:: out_xtra_2_info_REAL		
	!subroutine params:
	integer												:: index_par_i
	real*8, dimension(:,:), allocatable					:: IC_body_all_const_vec
	!real*8, dimension(:,:), allocatable					:: IC_par_pos, IC_par_vel
	real*8, dimension(:,:,:), allocatable				:: IC_par_q, IC_par_qdot
	real*8, dimension(:), allocatable					:: IC_Yevol_vec
	integer, dimension(6)								:: Nbody_solver_params_1_INT
	real*8, dimension(4)								:: Nbody_solver_params_2_REAL
	real*8												:: t_start, t_now, t_simtime_sec_now 	!fortran time tester
	integer												:: end_simulation
	real*8												:: delta_t_evo
	real*8, dimension(10)								:: Return_3body_Info_REAL, endsim_Return_Info_arr_REAL
	integer, dimension(10)								:: Return_3body_Info_INT, endsim_Return_Info_arr_INT
	real*8, dimension(10)								:: Return_Nbody_info_REAL_1, Return_Nbody_info_REAL_2
	integer, dimension(10)								:: endsim_out_xtra_info_INT
	real*8, dimension(10)								:: endsim_out_xtra_info_REAL
	real*8, dimension(10)								:: endsim_out_xtra_2_info_REAL
	integer, dimension(10)								:: Return_Nbody_info_INT_1
	real*8, dimension(10)								:: Return_3body_Info_REAL_XTRAINFO
	integer, dimension(5)								:: Return_Nbody_endstate
	real*8, dimension(10)								:: Return_endstate_binparams
	integer												:: endsim_end_state_flag, NBsystem_state_flag
    integer                                             :: NBsystem_bin_i, NBsystem_bin_j   ! Added by Mike
	integer												:: nrc_usrstate
	integer												:: Return_pass_tidalthreshold_yesno
	real*8												:: dist_ij
	real*8, dimension(3)								:: dist_ij_3tarr
	integer												:: IMS_rp_counter, IMS_binsin_counter
	real*8												:: a_bin_ij, e_bin_ij, rdist_bin_ij, tidradius_total_bin_ij, Etot_bin_ij
	integer												:: flag_usrstate20, flag_usrstate21, flag_usrstate22
	integer												:: bound_bin_ij, bin_i, bin_j
	real*8												:: sphradius_total_bin_ij
	integer, dimension(10)								:: output_userstate_info_INT
	real*8, dimension(10)								:: output_userstate_info_REAL
	integer												:: IMSbin_yesno, IMSbin_0_yesno, IMSbin_1_yesno
	real*8												:: maxrad_bin_i, maxrad_bin_j
	real*8												:: rmin_12, rmin_13, rmin_23, r_12, r_13, r_23
	real*8												:: MRini_IMS_a, MRini_IMS_e
	real*8, dimension(n_particles)                      :: r_track
	real*8, dimension(6)                                :: rb_track
	real*8, dimension(3)								:: Lvec
	real*8												:: E_ini, L_ini, E_kin, E_pot
	!odepack solver params:
	!EXTERNAL											:: F_Y_Ydot
	integer												:: NEQ
	integer												:: IOPT, IOUT, ISTATE, ITASK, ITOL, MF, LRW, LIW
	real*8												:: JAC = 0.0		!'dummy value' when the JAC is not used.
	real*8												:: T, TOUT, endsim_TOUT, RTOL, ATOL
	real*8,		dimension(:), allocatable				:: Y, endsim_Y
	real*8, 	dimension(:), allocatable				:: RWORK 
	integer, 	dimension(:), allocatable				:: IWORK
	!--------------------------------------------------------
	
		
	!--------------------------------------------------------
	!open all files for info output:
	!--------------------------------------------------------
	!---------------------------------------
	!Code version 1:
	!---------------------------------------
	if (IC_code_version .EQ. 1) then
	!---------------------------------------
	open (unit=10, file='Nbody_positions.dat',            			status='REPLACE', action='write')
	open (unit=11, file='Nbody_times.dat',                   		status='REPLACE', action='write')
	!open (unit=12, file='NbodyTides_dataout_Eself.dat',			status='REPLACE', action='write')
	!open (unit=13, file='NbodyTides_dataout_Etot.dat',			status='REPLACE', action='write')
	!open (unit=14, file='NbodyTides_dataout_binij_info.dat',	status='REPLACE', action='write')
	!open (unit=15, file='NbodyTides_dataout_full_q.dat',		status='REPLACE', action='write')
	!open (unit=16, file='NbodyTides_dataout_resonance.dat',		status='REPLACE', action='write')
	!open (unit=17, file='NbodyTides_dataout_user_states.dat',	status='REPLACE', action='write')
	open (unit=18, file='Nbody_velocities.dat',            			status='REPLACE', action='write')
	open (unit=20, file='Nbody_endsim_info.dat',           			status='REPLACE', action='write')
	!open (unit=21, file='NbodyTides_dataout_3bodybinsin.txt', 	status='REPLACE', action='write')
    open (unit=22, file='Nbody_initial_conditions.dat',               status='REPLACE', action='write')
	!‘REPLACE’ : If the file already exists then it is deleted and a new file created with the same
	!name. If the file does not exist then a new file is created.
	!---------------------------------------
	endif	!code version 1
	!---------------------------------------
	!--------------------------------------------------------
	!--------------------------------------------------------
	
	
	!--------------------------------------------------------
	!Read in data:
	!--------------------------------------------------------
	!---------------------------------------
	!Code version 1:
	!---------------------------------------
	if (IC_code_version .EQ. 1) then
	!---------------------------------------
	!nr particles:
	read(*,*),  n_particles
	!calc:
	tot_nr_Y_evol_eqs	= n_particles*length_Y_per_body_n
	!Allocate:
	allocate(IC_par_pos(n_particles,3))
	allocate(IC_par_vel(n_particles,3))
	allocate(IC_par_q(n_particles,3,3))
	allocate(IC_par_qdot(n_particles,3,3))
	allocate(IC_body_all_const_vec(n_particles,2))
	allocate(mass(n_particles))
	allocate(radius(n_particles))
	allocate(gas_n(n_particles))
	allocate(gas_gamma(n_particles))
	allocate(Mqp_sph(n_particles))
	allocate(evoTides_yesno(n_particles))
	allocate(RigidSph_yesno(n_particles))
	allocate(IC_Yevol_vec(tot_nr_Y_evol_eqs))
	allocate(endsim_Y(tot_nr_Y_evol_eqs))
	!solver params:
	read(*,*), Nbody_solver_params_1_INT(:)		! [use_1PN, use_2PN, use_25PN, outputinfo_screenfiles, Identify_3Body_endstate, max_sim_nrsteps, ...]
	read(*,*), Nbody_solver_params_2_REAL(:)	! [scale_dt, max_sim_time, max_simtime_sec, evolvetides_threshold, ENDbinsingle_threshold, IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, ...]
	!Input ini pos, vel, mass, radius,.. for each body i:
	do i=1, n_particles,1 
		!constants:
		read(*,*), IC_body_all_const_vec(i,:)	! [Mass, Radi, RigidSph_yesno, evoTides_yesno, gas_n, gas_gamma, Mqp, 0,0,0]
		!pos,vel:
		read(*,*), IC_par_pos(i,:)				!pos
		read(*,*), IC_par_vel(i,:)				!vel
        write(22,*) IC_body_all_const_vec(i,1:2), IC_par_pos(i,:), IC_par_vel(i,:) !write intiial conditions
	enddo
	!---------------------------------------
	endif	!code version 1
	!---------------------------------------
	!---------------------------------------	

	
	!--------------------------------------------------------
	!Define:
	!--------------------------------------------------------
	!Nbody_solver_params_1_INT:  [use_1PN, use_2PN, use_25PN, outputinfo_screenfiles, Identify_3Body_endstate, max_sim_nrsteps, ...]
	!Nbody_solver_params_2_REAL: [scale_dt, max_sim_time, max_simtime_sec, evolvetides_threshold, ENDbinsingle_threshold, IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, ...]
	!from Nbody_solver_params_1_INT:
	use_1PN					    	= Nbody_solver_params_1_INT(1)
	use_2PN 					  	= Nbody_solver_params_1_INT(2)
    use_25PN                        = Nbody_solver_params_1_INT(3)
	outputinfo_screenfiles			= Nbody_solver_params_1_INT(4)
	max_sim_nrsteps					= Nbody_solver_params_1_INT(5)
	downsample  					= Nbody_solver_params_1_INT(6)
	!from Nbody_solver_params_2_REAL:
	scale_dt						= Nbody_solver_params_2_REAL(1)
	max_sim_time					= Nbody_solver_params_2_REAL(2)
	max_simtime_sec					= Nbody_solver_params_2_REAL(3)
	insp_threshold					= Nbody_solver_params_2_REAL(4)
	!from IC_body_all_const_vec:
	mass(:)							= IC_body_all_const_vec(:,1)
	radius(:)						= IC_body_all_const_vec(:,2)
	!Identity matrix:
	I_MAT(:,:) = 0d0
	I_MAT(1,1) = 1d0
	I_MAT(2,2) = 1d0
	I_MAT(3,3) = 1d0
	!---------------------------------------
	!--------------------------------------------------------
	!--------------------------------------------------------


	!--------------------------------------------------------
	!Make IC Y evolution vec:
	!--------------------------------------------------------
	!fill IC_Yevol_vec: (order: pos,vel,q,qdot)
	do i=1, n_particles,1
		!index shift (is):
		is = (i-1)*length_Y_per_body_n 
		!pos (3):
		IC_Yevol_vec(is+1:is+3)		= IC_par_pos(i,:)
		!vel (3):
		IC_Yevol_vec(is+4:is+6)		= IC_par_vel(i,:)
		!q (9)
		IC_Yevol_vec(is+7:is+15)	= [IC_par_q(i,1,:), IC_par_q(i,2,:), IC_par_q(i,3,:)]
		!qdot (9)
		IC_Yevol_vec(is+16:is+24)	= [IC_par_qdot(i,1,:), IC_par_qdot(i,2,:), IC_par_qdot(i,3,:)]
	enddo
	!--------------------------------------------------------
	!--------------------------------------------------------


	!--------------------------------------------------------
	!set/initialize odepack params:
	!--------------------------------------------------------
	!for SETTING: MF=10:
	NEQ			= tot_nr_Y_evol_eqs
	MF			= 10			
	LRW			= 20 + 16*NEQ
	LIW			= 20
	allocate(Y(NEQ))
	allocate(RWORK(LRW))
	allocate(IWORK(LIW))
	Y(:)		= 0		!Will be set to 'IC_Yevol_vec' below before run.
	RWORK(:)	= 0.0	!initialize
	IWORK(:)	= 0		!initialize
	ITOL 		= 1
	RTOL 		= 1e-12
	ATOL 		= 1e-12
	ITASK 		= 1
	ISTATE 		= 1
	!opt inputs to solver:
	IOPT 		= 1 !must be =1 when MXSTEP is set.
	IWORK(6)	= 500000000 ! = MXSTEP
	!--------------------------------------------------------
	!--------------------------------------------------------


	!--------------------------------------------------------
	!initialize values:
	!--------------------------------------------------------
	!flags/counters:
	NBsystem_state_flag					= 0
	end_simulation						= 0
	stepc								= 0
	IMSbin_yesno						= 0
	IMSbin_0_yesno						= 1
	IMSbin_1_yesno						= 1
	IMS_binsin_counter					= 0
	pass_tidalthreshold_yesno			= 0
	nrc_usrstate						= 1
	flag_usrstate20						= 0
	flag_usrstate21 					= 0
	flag_usrstate22 					= 0
	IMS_rp_counter						= 0
    ds                                  = 0
	!arrays:
	Return_3body_Info_REAL(:)			= 0.0
	Return_3body_Info_INT(:)			= 0
	Return_3body_Info_REAL_XTRAINFO(:)	= 0.0
	out_endstate_INT(:)					= 0
	out_endstate_REAL(:)				= 0.0
	out_usrstate_INT(:,:)				= 0
	out_usrstate_REAL(:,:)				= 0.0
	out_xtra_info_INT(:)				= 0
	out_xtra_info_REAL(:)				= 0.0
	out_xtra_2_info_REAL(:)				= 0.0
	endsim_out_xtra_info_INT(:)			= 0
	endsim_out_xtra_info_REAL(:)		= 0.0
	endsim_out_xtra_2_info_REAL(:)		= 0.0
	output_userstate_info_INT(:)		= 0
	output_userstate_info_REAL(:)		= 0.0
	rmin_12								= 1e10
	rmin_13								= 1e10
	rmin_23								= 1e10
	MRini_IMS_a							= 0.0
	MRini_IMS_e							= 0.0
	!--------------------------------------------------------
	!--------------------------------------------------------

	
	!-----------------------------------------------------------------------------
	!INTEGRATE NBODY SYSTEM:
	!-----------------------------------------------------------------------------
	!start wall clock time:
	call cpu_time(t_start)
	!Set ICs:
	Y(:) 	= IC_Yevol_vec(:)
	T 		= 0d0
	TOUT	= 0d0	!changed below to 'delta_t_evo'

	!--------------------------------------------------------
	!Evolve system:
	!--------------------------------------------------------
	DO WHILE (end_simulation .EQ. 0)			
		
		!----------------------------------------------------
		!Evolve from TOUT,TOUT+delta_t_evo:
		!----------------------------------------------------
		!delta time evo:
		delta_t_evo = scale_dt*Nbody_dt(Y)	!We always use variable time step (Nbody_dt(Y) is based only on Nbody pos,vel).
		!evole system Y dt:
		TOUT	= TOUT + delta_t_evo
		CALL DLSODE(F_Y_Ydot, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
		!We now have: (Y, TOUT)
		!----------------------------------------------------
	
		!----------------------------------------------------
		!N-body section:
		!----------------------------------------------------
		!-----------------------------------	
		!Check for endstate:
		!-----------------------------------
		CALL Nbody_endstate_sub(Y, Return_Nbody_endstate, Return_endstate_binparams)            
		NBsystem_state_flag = Return_Nbody_endstate(1)		! (.NE. 0 if endstate is found.) ( = 0 if no endstate.)

        !-----------------------------------
        !Write to data files
        !-----------------------------------
		if (IC_code_version .EQ. 1) then
		!if (outputinfo_screenfiles .EQ. 1)	CALL output_Nbody_sub(Y,TOUT)
		if (downsample .NE. 0) CALL output_Nbody_sub(Y,TOUT)
		endif	!code version 1					




!!!THIS SECTION IS NOT USED FOR THE bBBH 4-BODY STUDY!!!
		!----------------------------------------------------
		!Check for tidal threshold:
		!----------------------------------------------------
		!check if tidal threshold is passed (ONLY check if it has NOT been found so far!):
		if (pass_tidalthreshold_yesno .EQ. 0) then
			!call subroutine for tidal threshold check:
			CALL check_tidal_threshold_sub(Y, Return_pass_tidalthreshold_yesno)
			!if tidal threshold is passed set 'pass_tidalthreshold_yesno = 1':
			if (Return_pass_tidalthreshold_yesno .EQ. 1) pass_tidalthreshold_yesno = 1
			!print info to screen:
			if (IC_code_version .EQ. 1) then
			if (pass_tidalthreshold_yesno .EQ. 1) print*, '=========== TIDAL THRESHOLD IS NOW PASSED ==========='
			endif	!code version 1
		endif
		!----------------------------------------------------				
				
		if (Identify_3Body_endstate .EQ. 1) then
		!NOTE: if we also use 3body-part then NBsystem_state_flag from here will be overwritten in 3body-part.
		!-----------------------------------
		!Get Info about most-bound-pair [i,j]
		!-----------------------------------
		CALL Nbody_info_module_sub(Y, Return_Nbody_info_REAL_1, Return_Nbody_info_REAL_2, Return_Nbody_info_INT_1)
		!Define:
		Etot_bin_ij		= Return_Nbody_info_REAL_1(3)
		a_bin_ij		= Return_Nbody_info_REAL_1(4)
		e_bin_ij		= Return_Nbody_info_REAL_1(5)
		rdist_bin_ij	= Return_Nbody_info_REAL_2(1)
		maxrad_bin_i	= Return_Nbody_info_REAL_2(2)	!radius along maximum principal axis (rad(sph)*max(principal a))
		maxrad_bin_j	= Return_Nbody_info_REAL_2(3)	!radius along maximum principal axis (rad(sph)*max(principal a))
		tidradius_total_bin_ij	= maxrad_bin_i + maxrad_bin_j
		bin_i 			= Return_Nbody_info_INT_1(1)
		bin_j 			= Return_Nbody_info_INT_1(2)
		sphradius_total_bin_ij	= radius(bin_i)+radius(bin_j)
		if (Etot_bin_ij .LT. 0.0)	bound_bin_ij = 1
		if (Etot_bin_ij .GT. 0.0)	bound_bin_ij = 0
		!format:
		!Return_Nbody_info_REAL_1(1:5)	= [binary_info_arr_ij(:)]
		!Return_Nbody_info_INT_1(1:2)	= [bin_i, bin_j]
		!Return_Nbody_info_REAL_2(1:3)	= [r_ij, radtidmax_bin_i, radtidmax_bin_j]	

		!-----------------------------------
		!user states:
		!-----------------------------------
		!USER STATE: 20 EXAMPLE! YOU CAN ADD MODE USR STATES!!
		if (flag_usrstate20 .EQ. 0 .and. bound_bin_ij .EQ. 1 .and.	(a_bin_ij/tidradius_total_bin_ij)		.LT. 3.0) then
			!SELECT output info:
			output_userstate_info_INT(1:5) = [20, bin_i, bin_j, IMS_rp_counter, IMS_binsin_counter]	!'IMS_rp_counter, IMS_binsin_counter' are only active if we have 3 objs.
			!SAVE to file (code version 1):
			if (IC_code_version .EQ. 1) then
				write(17,*) output_userstate_info_INT	!write INT info
				write(17,*) output_userstate_info_REAL	!write REAL info
				print*, 'usrstate20'
			endif	!code version 1
			!SAVE in arrays (code version 2):
			if (IC_code_version .EQ. 2) then
				out_usrstate_INT(nrc_usrstate,:)	= output_userstate_info_INT(:)
				out_usrstate_REAL(nrc_usrstate,:)	= output_userstate_info_REAL(:)
				nrc_usrstate						= nrc_usrstate + 1
			endif	!code version 2
			!FLAG usrstate:
			flag_usrstate20 = 1
		endif		
		!-----------------------------------
		!----------------------------------------------------
        endif
	
		!----------------------------------------------------
		!3-body section:
		!----------------------------------------------------	
		!-----------------------------------
		!Get Info: 3Body state [[i,j],k]
		!-----------------------------------
		if (Identify_3Body_endstate .EQ. 1) then
			CALL analyze_3body_state_info_sub(Y, Return_3body_Info_REAL, Return_3body_Info_INT, Return_3body_Info_REAL_XTRAINFO)
			!Output format:
			!Return_3body_Info_REAL(1:10)			= [1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
			!Return_3body_Info_INT(1:5)				= [out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno]
			!Return_3body_Info_REAL_XTRAINFO(1:10)	= [out_dist_bin_ij, out_dist_bin_12, out_dist_bin_13, out_dist_bin_23, out_pos_CMij_wrt_sink(:), out_vel_CMij_wrt_sink(:)]			
			NBsystem_state_flag = Return_3body_Info_INT(1)				! (.NE. 0 if endstate is found.) ( = 0 if no endstate.)
			IMSbin_yesno		= Return_3body_Info_INT(5)				! if IMS yes = 1, otherwise 0.
		endif	
		!-----------------------------------
		!3-Body IMS Resonance analysis:
		!-----------------------------------
		if (NBsystem_state_flag .EQ. 0 .and. Identify_3Body_endstate .EQ. 1) then
			
			!IMS bin-sin config 'now'(1) yes/no:
			IMSbin_1_yesno = IMSbin_yesno	!output from 'analyze_3body_state_info_sub' above.
			
			!IF IMS binary is created (label 1):	
			if (IMSbin_1_yesno .EQ. 1 .and. IMSbin_0_yesno .EQ. 0) then
				!write info to file
				if (IC_code_version .EQ. 1) then
					write(16,*) 1, Return_3body_Info_INT(2:4), Return_3body_Info_REAL(3:5), Return_3body_Info_REAL(8:10)
				endif	!code version 1
				!reset peri-center counter and dist array:
				IMS_rp_counter		= 0
				dist_ij_3tarr(:)	= 0.0
				!counter:
				IMS_binsin_counter	= IMS_binsin_counter+1
			endif
			
			!IF IMS binary is destroyed (label 0):
			if (IMSbin_1_yesno .EQ. 0 .and. IMSbin_0_yesno .EQ. 1) then
				!write info to file
				if (IC_code_version .EQ. 1) then
					write(16,*) 0, Return_3body_Info_INT(2:4), Return_3body_Info_REAL(3:5), Return_3body_Info_REAL(8:10)
				endif	!code version 1				
				!reset peri-center counter and dist array:
				IMS_rp_counter		= 0
				dist_ij_3tarr(:)	= 0.0 
				!counter:
				!IMS_binsin_counter = IMS_binsin_counter+1
			endif

			!Count peri-center passages if we are in IMS state:
			if (IMSbin_0_yesno .EQ. 1 .and. IMSbin_1_yesno .EQ. 1) then
				dist_ij				= Return_3body_Info_REAL_XTRAINFO(1)
				dist_ij_3tarr(1:2)	= dist_ij_3tarr(2:3)
				dist_ij_3tarr(3)	= dist_ij
				!check for peri-center passage:
				if (dist_ij_3tarr(1) .GT. dist_ij_3tarr(2) .and. dist_ij_3tarr(3) .GT. dist_ij_3tarr(2)) IMS_rp_counter = IMS_rp_counter+1				
			endif 
			
			!save most recent initial (MRini) IMS (a,e):
			if (IMSbin_1_yesno .EQ. 1 .and. IMSbin_0_yesno .EQ. 0) then
				MRini_IMS_a	= Return_3body_Info_REAL(4)
				MRini_IMS_e	= Return_3body_Info_REAL(5)
			endif			
						
			!Update IMS 'before'(0):
			IMSbin_0_yesno = IMSbin_1_yesno
			
		!-----------------------------------
		!Track rmin between 12,13,23:
		!-----------------------------------
		!Return_3body_Info_REAL_XTRAINFO(1:10)	= [out_dist_bin_ij, out_dist_bin_12, out_dist_bin_13, out_dist_bin_23, out_pos_CMij_wrt_sink(:), out_vel_CMij_wrt_sink(:)]			
		r_12 = Return_3body_Info_REAL_XTRAINFO(2)
		r_13 = Return_3body_Info_REAL_XTRAINFO(3)
		r_23 = Return_3body_Info_REAL_XTRAINFO(4)
		if (r_12 .LT. rmin_12)	rmin_12 = r_12
		if (r_13 .LT. rmin_13)	rmin_13 = r_13
		if (r_23 .LT. rmin_23)	rmin_23 = r_23
		!-----------------------------------
		!----------------------------------------------------
	
		!----------------------------------------------------
		!Extra info for detailed analysis:
		!----------------------------------------------------
		if (IC_code_version .EQ. 1) then
		if (outputinfo_screenfiles .EQ. 1)	CALL TEST_sub(Y,TOUT)
		endif	!code version 1					
		!----------------------------------------------------
		endif	
!!!THIS HAS BEEN REMOVED FOR THE bBBH 4-BODY STUDY!!!




	
		!----------------------------------------------------
		!Check: Time/step limits, DLSODE status:
		!----------------------------------------------------	
		!Time limit (physical time in sim):
		if (TOUT .GT. max_sim_time)					NBsystem_state_flag = 10
		!nr steps limit:
		if (stepc .GT. max_sim_nrsteps)				NBsystem_state_flag = 11
		!Time limit (wall clock, i.e. time the sim runs):
		call cpu_time(t_now)
		t_simtime_sec_now = t_now-t_start
		if (t_simtime_sec_now .GT. max_simtime_sec)	NBsystem_state_flag = 12
		!output status of LSODE:
		if (ISTATE .LT. 0)							NBsystem_state_flag = 13
		!----------------------------------------------------		
	
		!----------------------------------------------------
		!If end-state is found:
		!----------------------------------------------------
		if (NBsystem_state_flag .NE. 0) then
			!stop sim:
			end_simulation = 1
			!endstate flag:
			endsim_end_state_flag		= NBsystem_state_flag
			!endstate T, Y:
			endsim_TOUT					= TOUT
			endsim_Y					= Y
			!endstate output: INT		
			endsim_Return_Info_arr_INT(:)		= Return_3body_Info_INT
			endsim_Return_Info_arr_INT(1)		= NBsystem_state_flag	!put general endstateflag into first entry: makes it easier to analyze data later
			endsim_Return_Info_arr_INT(6)		= IMS_rp_counter
			endsim_Return_Info_arr_INT(7)		= IMS_binsin_counter
			!endstate output: REAL
			endsim_Return_Info_arr_REAL(:)		= Return_3body_Info_REAL
			!xtra info:
			endsim_out_xtra_info_REAL(1:3)		= [rmin_12, rmin_13, rmin_23]
			endsim_out_xtra_info_REAL(4:5)		= [MRini_IMS_a, MRini_IMS_e]		
			!xtra info 2:
			endsim_out_xtra_2_info_REAL(1:6)	= Return_3body_Info_REAL_XTRAINFO(5:10)
		endif
		!----------------------------------------------------		
	
		!----------------------------------------------------
		!Counters:
		!----------------------------------------------------
		stepc = stepc+1 
		!----------------------------------------------------
	
	END DO	!stop sim if 'end_simulation = 1' at this point.
	!--------------------------------------------------------
	!--------------------------------------------------------
	!stop wall clock time:
	call cpu_time(t_now)	
	!-----------------------------------------------------------------------------
	!-----------------------------------------------------------------------------


	!--------------------------------------------------------
	!Save output info:
	!--------------------------------------------------------
	!---------------------------------------
	!Code version 1:
	!---------------------------------------
	if (IC_code_version .EQ. 1) then
	!---------------------------------------
	!write endstate info to file:
	write(20,*) endsim_end_state_flag, Return_Nbody_endstate(2:5), Return_endstate_binparams(:)
!	write(20,*) endsim_Return_Info_arr_INT
!	write(20,*) endsim_Return_Info_arr_REAL
	!Print info to screen:
	print*, 'Time = seconds:	', t_now-t_start
	if (outputinfo_screenfiles .EQ. 1) then
		print*, endsim_end_state_flag
        print*, Return_Nbody_endstate(2:5)
        print*, Return_endstate_binparams(:)
		print*, mass(1), radius(1), endsim_Y(1:6)
        print*, mass(2), radius(2), endsim_Y(25:30)
        if (n_particles .GT. 2)             print*, mass(3), radius(3), endsim_Y(49:54)
        if (n_particles .GT. 3)             print*, mass(4), radius(4), endsim_Y(73:78)
	endif
	!close all files:
	close(10)
	close(11)
	close(12)
	close(13)
	close(14)
	close(15)
	close(16)
	close(17)
	close(18)
	close(20)
	close(21)
    close(22)
	!---------------------------------------
	endif	!code version 1
	!---------------------------------------
	!---------------------------------------
	!Code version 2:
	!---------------------------------------
	if (IC_code_version .EQ. 2) then
	!---------------------------------------
	out_endstate_INT(:)		= endsim_Return_Info_arr_INT(:)
	out_endstate_REAL(:)	= endsim_Return_Info_arr_REAL(:)
	out_xtra_info_REAL(:)	= endsim_out_xtra_info_REAL(:)
	out_xtra_2_info_REAL(:)	= endsim_out_xtra_2_info_REAL(:)
	!---------------------------------------
	endif	!code version 2
	!---------------------------------------
	!--------------------------------------------------------
	!--------------------------------------------------------


	!--------------------------------------------------------
	!De-allocate allocated arrays:
	!--------------------------------------------------------
	!This is important when running on multiple cores.
	if (allocated (IC_par_pos)) deallocate (IC_par_pos)
	if (allocated (IC_par_vel)) deallocate (IC_par_vel)
	if (allocated (IC_par_q)) deallocate (IC_par_q)
	if (allocated (IC_par_qdot)) deallocate (IC_par_qdot)
	if (allocated (IC_body_all_const_vec)) deallocate (IC_body_all_const_vec)
	if (allocated (radius)) deallocate (radius)
	if (allocated (gas_n)) deallocate (gas_n)
	if (allocated (gas_gamma)) deallocate (gas_gamma)
	if (allocated (Mqp_sph)) deallocate (Mqp_sph)
	if (allocated (evoTides_yesno)) deallocate (evoTides_yesno)
	if (allocated (RigidSph_yesno)) deallocate (RigidSph_yesno)
	if (allocated (IC_Yevol_vec)) deallocate (IC_Yevol_vec)
	if (allocated (endsim_Y)) deallocate (endsim_Y)
	if (allocated (mass)) deallocate (mass)	
	if (allocated (Y)) deallocate (Y)
	if (allocated (RWORK)) deallocate (RWORK)
	if (allocated (IWORK)) deallocate (IWORK)
	!--------------------------------------------------------
	!--------------------------------------------------------
	
	
	END SUBROUTINE main_evolve_Nbody_system_sub
	!----------------------------------------------------------------------------------------------------
	!----------------------------------------------------------------------------------------------------
	
	
END MODULE module_Nbody_AffineTides_solver
!--------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------
