
!--------------------------------------------------------------------------------------------------------
SUBROUTINE interfacesub(IN_IC_code_version, IN_IC_simparams_INT, IN_IC_simparams_REAL,	&
						IN_dimlen_IC_nbody, IN_IC_nbody_const_posvel_qqdot_etc_arr,		&
						OUT_out_endstate_INT, OUT_out_endstate_REAL, 					&
						OUT_out_usrstate_INT, OUT_out_usrstate_REAL,					&
						OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL)
!--------------------------------------------------------------------------------------------------------
USE module_Nbody_AffineTides_solver	!from this mod we use 'main_evolve_Nbody_system_sub'
IMPLICIT NONE
!Input:
integer 									:: IN_IC_code_version
integer, dimension(10)						:: IN_IC_simparams_INT
real*8, dimension(10)						:: IN_IC_simparams_REAL
integer										:: IN_dimlen_IC_nbody
real*8, dimension(IN_dimlen_IC_nbody,10)	:: IN_IC_nbody_const_posvel_qqdot_etc_arr
!ouput:
integer, dimension(10)						:: OUT_out_endstate_INT
real*8, dimension(10)						:: OUT_out_endstate_REAL
integer, dimension(10,10)					:: OUT_out_usrstate_INT
real*8, dimension(10,10)					:: OUT_out_usrstate_REAL
integer, dimension(10)						:: OUT_out_xtra_info_INT
real*8, dimension(10)						:: OUT_out_xtra_info_REAL
real*8, dimension(10)						:: OUT_out_xtra_2_info_REAL

!--------------------------------------------------------
!THE LINES BELOW WITH '!f2py' MUST NOT BE DELETED!!!
!--------------------------------------------------------
!f2py intent(in) IN_IC_code_version
!f2py intent(in) IN_IC_simparams_INT
!f2py intent(in) IN_IC_simparams_REAL
!f2py intent(in) IN_dimlen_IC_nbody
!f2py intent(in) IN_IC_nbody_const_posvel_qqdot_etc_arr
!f2py depend(IN_dimlen_IC_nbody) IN_IC_nbody_const_posvel_qqdot_etc_arr
!f2py intent(out) OUT_out_endstate_INT
!f2py intent(out) OUT_out_endstate_REAL
!f2py intent(out) OUT_out_usrstate_INT
!f2py intent(out) OUT_out_usrstate_REAL
!f2py intent(out) OUT_out_xtra_info_INT
!f2py intent(out) OUT_out_xtra_info_REAL
!f2py intent(out) OUT_out_xtra_2_info_REAL
!--------------------------------------------------------

CALL main_evolve_Nbody_system_sub(	IN_IC_code_version, IN_IC_simparams_INT, IN_IC_simparams_REAL,	&
									IN_dimlen_IC_nbody, IN_IC_nbody_const_posvel_qqdot_etc_arr,		&
									OUT_out_endstate_INT, OUT_out_endstate_REAL, 					&
									OUT_out_usrstate_INT, OUT_out_usrstate_REAL,					&
									OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL)

				
END SUBROUTINE interfacesub
!--------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------




!--------------------------------------------------------------------------------------------------------
PROGRAM Nbody_AffineTides_solver
!--------------------------------------------------------------------------------------------------------
USE module_Nbody_AffineTides_solver !from this mod we use 'main_evolve_Nbody_system_sub'
IMPLICIT NONE
!input:
integer										:: IC_code_version
integer, dimension(10)						:: IC_simparams_INT
real*8, dimension(10)						:: IC_simparams_REAL
integer					 					:: dimlen_IC_nbody
real*8, dimension(27,10)					:: IC_nbody_const_posvel_qqdot_etc_arr
!ouput:	
integer, dimension(10)						:: out_endstate_INT
real*8, dimension(10)						:: out_endstate_REAL
integer, dimension(10,10)					:: out_usrstate_INT
real*8, dimension(10,10)					:: out_usrstate_REAL
integer, dimension(10)						:: out_xtra_info_INT
real*8, dimension(10)						:: out_xtra_info_REAL
real*8, dimension(10)						:: out_xtra_2_info_REAL


!--------------------------------------------------------
!Notes:
!--------------------------------------------------------
!If the code is run as a normal fortran code (./...) then this program part
!will be executed. We usually use this for testing and more detailed
!analysis since we from this mode have the possibility of printing
!out all different kinds of data. All in IC/out data is through .txt files,
!so all the input to the subroutine has no effect, except for the 'IC_code_version'
!which must be set = 1.
!code version = 1:
!This is the standard old version which allows for writing out .txt files.
!The ICs are read from an input file, and all ouput is written to files.
!The input and output arrays in the subroutine 'main_evolve_Nbody_system_sub'
!has no effect in this mode. Eveything is handled through input/output files.
!code version = 2:
!No reading or writing to files at all. All info incl. IC and final output
!is handled through the input arrays in the subroutine 'main_evolve_Nbody_system_sub'.
!This mode is inteded to be used in a parallel version of the program.
!--------------------------------------------------------

!--------------------------------------------------------
!Standard settings:
!--------------------------------------------------------
IC_code_version		= 1
IC_simparams_INT	= 1
IC_simparams_REAL	= 1
dimlen_IC_nbody		= 1
IC_nbody_const_posvel_qqdot_etc_arr	= 1
!--------------------------------------------------------

! !--------------------------------------------------------
! !MIKE'S VERSION SETTINGS
! !--------------------------------------------------------
! IC_code_version     = 4

! !--------------------------------------------------------
! !TEST SETTINGS!!!!!!!!
! !--------------------------------------------------------
! IC_code_version			= 2
! !sim params:
! !IC_simparams_INT		! [use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, FREE, outputinfo_screenfiles, nr_nbody_particles, ...]
! !IC_simparams_REAL		! [scale_dt, max_sim_time, evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec, IMSbinsingle_threshold, ...]
! IC_simparams_INT(:)		= [0,0,1,100000000,0,0,3,0,0,0]
! IC_simparams_REAL(:)	= [0.010000, 88584.216627, 0.100000, 0.010000, 100000000.0, 0.500000, 0.0, 0.0, 0.0, 0.0]
! !dimensions:
! dimlen_IC_nbody			= 27	!=3*9
! !obj 1:
! IC_nbody_const_posvel_qqdot_etc_arr(1,:) = [1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(2,:) = [-0.538861, 93.044669, 64.479373, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(3,:) = [0.004368, 0.117192, -0.027313, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(4,:) = [1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(5,:) = [0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(6,:) = [0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(7,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(8,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(9,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! !obj 2:
! IC_nbody_const_posvel_qqdot_etc_arr(10,:) = [1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(11,:) = [-22.038986, 93.044669, 64.479373, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(12,:) = [0.004368, -0.187805, -0.027313, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(13,:) = [1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(14,:) = [0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(15,:) = [0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(16,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(17,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(18,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! !obj 3:
! IC_nbody_const_posvel_qqdot_etc_arr(19,:) = [1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(20,:) = [22.577848, -186.089337, -128.958746, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(21,:) = [-0.008735, 0.070613, 0.054626, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(22,:) = [1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(23,:) = [0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(24,:) = [0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(25,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(26,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! IC_nbody_const_posvel_qqdot_etc_arr(27,:) = [0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
! !--------------------------------------------------------

CALL main_evolve_Nbody_system_sub(	IC_code_version, IC_simparams_INT, IC_simparams_REAL,						&
									dimlen_IC_nbody, IC_nbody_const_posvel_qqdot_etc_arr,						&
									out_endstate_INT, out_endstate_REAL, out_usrstate_INT, out_usrstate_REAL,	&
									out_xtra_info_INT, out_xtra_info_REAL, out_xtra_2_info_REAL)


END PROGRAM Nbody_AffineTides_solver
!--------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------
