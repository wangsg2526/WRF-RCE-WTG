!WRF:MEDIATION_LAYER:PHYSICS
!

MODULE module_wtg_driver4gass
! This is driver for WTG calculations
! coupling_largescale_omega is the subroutine for WTG calculations
! rce_wtg_flag = 1,2,3,4
! WTG_DGW_P_COUPLING is the legacy code from the WTG-intercomparison project 
! rce_wtg_flag = 11, 12
! Shuguang Wang (wangsg@gmail.com)

CONTAINS

!------------------------------------------------------------------
   SUBROUTINE wtg_driver4gass(                                          &
                  itimestep,dt                     &
                 ,bldt,curr_secs,adapt_step_flag                   &
                 ,rublten,rvblten,rthblten                         &
!SW==
                 ,scm_vert_adv_opt                &
                 ,wtg_dt &
                 ,omega_wtg, the_wtg_forcing, qv_wtg_forcing       &
                 ,wtg_qtend, wtg_qtend_v, wtg_qtend_h, wtg_ttend   &
! diagnostics
                 ,th_avg, tmk_avg, qv_avg                                                   &
                 ,phalf_avg, phi_avg, rho_avg                     &
                 ,ql_avg, qi_avg, qs_avg, qg_avg, qr_avg           &
                 ,uuu_avg, vvv_avg                                 &    
                 ,clf_avg, clh_avg                                 &
                 ,rain_avg, evap_avg, lhf_avg, hfx_avg  &
                 ,mp_tten_avg, mp_qten_avg              &  
                 ,qrad_avg         &
                 ,tot_tten_avg, tot_qten_avg             &
                 ,areacloudupdraft_avg,  masscloudupdraft_avg      & 
                 ,areaupdraftcore_avg, massupdraftcore_avg         &
                 ,areabuoyantcloud_avg, massbuoyantcloud_avg       &
                 ,rainnc, achfx, aclhf, acqfx                             &
! diagnostics
                 ,qv_wtg_hadv_forcing                              &
                 ,qv_avg_flag                                     &                       
                 ,qadv_flag, qadv_tau,qforcing_tau, qadv_plev_up,qadv_plev_dn, qadv_popt        &
                 ,th_rce,  rce_wtg_flag, tao_time, rce_wtg_diag    &
                 ,qv_rce, qv_hadv_flag, qv_lateral_drawing_flag    &
                 ,pblh_opt, relax_umean_opt                        &
                 ,relax_umean_u, relax_umean_tscale                &
                 ,rqvwtgten,rthwtgten                              &
                 ,acqvwtg, acqvhadvwtg                             &
                 ,actwtg                              &
! for wave coupling
                 ,wavecoupling_opt                   &  
                 ,wavecoupling_topbc, wavecoupling_dmpopt                   &  
                 ,nl_bg, theta_bg, qv_bg &
                 , z_bg , prs_bg             & 
                 ,wavenumber, wavedamping            &  
! for vertical mode based spectral wtg
                 ,wtg_opt,wtg_cexp, wtg_nmodes  &
                 ,scm_m_vector, scm_zm, scm_m_value,scm_m_weight &
!SW==
                 ,pblh                         &
                 ,u_phy,v_phy,th_phy,rho                           &
                 ,p_phy,pi_phy,p8w,t_phy,dz8w,z, rw                &
                 ,bl_pbl_physics, ra_lw_physics, dx                &
                 ,stepbl        &
                 ,znu, znw, mut, p_top                             &
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,i_start,i_end, j_start,j_end, kts,kte, num_tiles &
             !  Optional moisture tracers
                 ,qv_curr, qc_curr                        &
                 ,qr_curr, qi_curr                        &
                 ,qs_curr, qg_curr                        &
                 ,rqvblten,rqcblten,  rqvmpten,rthmpten, rthraten  &
             !  Optional moisture tracer flags
                 ,f_qv,f_qc,f_qr,f_qi,f_qs,f_qg                    &
                                                                     )       
!------------------------------------------------------------------
#if (EM_CORE==1)
   USE module_state_description, ONLY :                            &
                   YSUSCHEME
#endif

   USE module_model_constants

! *** add new modules of schemes here


!SW==
    USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
!SW==
   !
!------------------------------------------------------------------
   IMPLICIT NONE
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
! Rho_d      dry density (kg/m^3)
! Theta_m    moist potential temperature (K)
! Qv         water vapor mixing ratio (kg/kg)
! Qc         cloud water mixing ratio (kg/kg)
! Qr         rain water mixing ratio (kg/kg)
! Qi         cloud ice mixing ratio (kg/kg)
! Qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- RUBLTEN       U tendency due to 
!                 PBL parameterization (m/s^2)
!-- RVBLTEN       V tendency due to 
!                 PBL parameterization (m/s^2)
!-- RTHBLTEN      Theta tendency due to 
!                 PBL parameterization (K/s)
!-- RQVBLTEN      Qv tendency due to 
!                 PBL parameterization (kg/kg/s)
!-- RQCBLTEN      Qc tendency due to 
!                 PBL parameterization (kg/kg/s)
!-- RQIBLTEN      Qi tendency due to 
!                 PBL parameterization (kg/kg/s)
!-- itimestep     number of time steps
!-- TSK           surface temperature (K)
!-- PBLH          PBL height (m)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (K)
!-- p_phy         pressure (Pa)
!-- pi_phy        exner function (dimensionless)
!-- p8w           pressure at full levels (Pa)
!-- t_phy         temperature (K)
!-- dz8w          dz between full levels (m)
!-- z             height above sea level (m)
!-- DX            horizontal space interval (m)
!-- DT            time step (second)
!-- n_moist       number of moisture species
!-- PSFC          pressure at the surface (Pa)
!-- ZS
!-- DZS
!-- P_QV          species index for water vapor
!-- P_QC          species index for cloud water
!-- P_QR          species index for rain water
!-- P_QI          species index for cloud ice
!-- P_QS          species index for snow
!-- P_QG          species index for graupel
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!
!******************************************************************
!------------------------------------------------------------------ 
!


   INTEGER,    INTENT(IN   )    ::     bl_pbl_physics, ra_lw_physics

   integer, intent(in) :: scm_vert_adv_opt
   INTEGER, INTENT(IN) :: wtg_dt
   INTEGER,    INTENT(IN   )    ::     ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       kts,kte, num_tiles

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                   &
  &                                    i_start,i_end,j_start,j_end

   INTEGER,    INTENT(IN   )    ::     itimestep,STEPBL
   INTEGER ::  STEPwtg
!

   REAL,       DIMENSION( kms:kme ),                              &
               OPTIONAL, INTENT(IN   )    ::               znu,   &
                                                           znw
!
   REAL,       INTENT(IN   )    ::     DT,DX
   REAL,       INTENT(IN   ),OPTIONAL    ::     bldt
   REAL,       INTENT(IN   ),OPTIONAL    ::     curr_secs
   LOGICAL,    INTENT(IN   ),OPTIONAL    ::     adapt_step_flag

!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(IN   )    ::                         p_phy, &
                                                          pi_phy, &
                                                             p8w, &
                                                             rho, &
                                                           t_phy, &
                                                           u_phy, &
                                                           v_phy, &
                                                            dz8w, &
                                                               z, &
                                                              rw, &
                                                          th_phy
!
!
!
   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               INTENT(INOUT)    ::                           &
                                                            PBLH

!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(INOUT)    ::                       RUBLTEN, &
                                                         RVBLTEN, &
                                                        RTHBLTEN

   REAL,       DIMENSION( ims:ime, jms:jme ),             &
              INTENT(INOUT)    ::        acqvwtg, acqvhadvwtg, actwtg


!
! Optional
!
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, INTENT(IN), OPTIONAL ::                             &
                                                      f_qv, f_qc, f_qi, f_qr, f_qs, f_qg      

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                      ! optional moisture tracers
                      ! 2 time levels; if only one then use CURR
                      qv_curr, qc_curr                 &
                     ,qr_curr, qi_curr                        &
                     ,qs_curr, qg_curr                        &
                     ,rqvblten,rqcblten              

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), intent(in) :: rqvmpten,rthmpten, rthraten              

   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               OPTIONAL                                         , &
               INTENT(IN)    ::                           mut
!
   REAL,       OPTIONAL, INTENT(IN)    ::               p_top


!  LOCAL  VAR


   REAL    :: DTMIN,DTBL,DTwtg



   INTEGER :: i,J,K,NK,jj,ij,its,ite,jts,jte
   LOGICAL :: flag_qv, flag_qc, flag_qr, flag_qi, flag_qs, flag_qg
   CHARACTER*256 :: message
   REAL    :: next_bl_time
   LOGICAL :: run_param
   LOGICAL :: do_adapt

!SW==
   REAL, DIMENSION( kms:kme),                 &
         OPTIONAL, INTENT(INOUT) ::  rqvwtgten, rthwtgten

    real, dimension(kms:kme), INTENT(INOUT) ::                  &
                      the_wtg_forcing, qv_wtg_forcing,  qv_wtg_hadv_forcing
    real ::  no_points, th_temp, uuu_sum, vvv_sum
    real ::  dth_sum, tmk_sum, rho_sum, qv_sum, phalf_sum, pfull_sum, z_sum, mu_sum, mu_avg
    integer :: irnc, jrnc
    integer :: k_abovepbl_first
    real    :: pblh_max, pblh_sum, pblh_domain
    integer, INTENT(IN) :: rce_wtg_flag, rce_wtg_diag, qv_hadv_flag 
    logical,intent(in) :: qv_lateral_drawing_flag
    integer,optional,INTENT(IN) ::  pblh_opt, relax_umean_opt
    real,optional,INTENT(IN) ::   relax_umean_u, relax_umean_tscale
    integer :: wtg_adv_opt
    real, DIMENSION( ims:ime , jms:jme ) :: th_2d,rho_2d, qv_2d, z_2d, pfull_2d, phalf_2d
    real, DIMENSION( ims:ime , jms:jme ) :: thetaint_tmp,qvint_tmp, qvhadvint_tmp 
    real, DIMENSION( kms:kme) :: z_avg

    real, DIMENSION( ims:ime , jms:jme ), intent(in) :: aclhf, rainnc, achfx, acqfx

    real, DIMENSION( kms:kme),intent(inout) :: th_avg,  qv_avg, tmk_avg
    real, DIMENSION( kms:kme),intent(inout) :: phi_avg, rho_avg, ql_avg, qi_avg, qr_avg, qs_avg, qg_avg  
    real, DIMENSION( kms:kme), intent(inout) :: phalf_avg, clf_avg, clh_avg
    real, DIMENSION( kms:kme) :: pfull_avg
    real, DIMENSION( kms:kme),intent(inout) :: uuu_avg, vvv_avg
    real, DIMENSION( kms:kme),intent(inout) :: mp_tten_avg, mp_qten_avg, qrad_avg
    real, DIMENSION( kms:kme),intent(inout) :: tot_tten_avg, tot_qten_avg
    real, DIMENSION( kms:kme),intent(inout) :: areacloudupdraft_avg,  masscloudupdraft_avg
    real, DIMENSION( kms:kme),intent(inout) :: areaupdraftcore_avg, massupdraftcore_avg
    real, DIMENSION( kms:kme),intent(inout) :: areabuoyantcloud_avg, massbuoyantcloud_avg 

    real :: ql_sum, qi_sum, qr_sum, qs_sum, qg_sum, clf_sum, clh_sum, mp_tten_sum, mp_qten_sum, qrad_sum, tot_tten_sum, tot_qten_sum
    real :: areacloudupdraft_sum, masscloudupdraft_sum, areaupdraftcore_sum, massupdraftcore_sum, areabuoyantcloud_sum, massbuoyantcloud_sum
    real :: rain_sum, evap_sum, lhf_sum, hfx_sum

    real, intent(inout) :: rain_avg, evap_avg, lhf_avg, hfx_avg

    real, DIMENSION( kms:kme),INTENT(IN) :: qv_rce,th_rce
    real, DIMENSION( kms:kme),intent(inout) :: omega_wtg
    real, DIMENSION( kms:kme) :: dthdp, dqvdeta, qtend, qtend_qvhadv, qtend_avg, ttend
    real ::  qtend1, qtend2, domgdeta
    real :: theta_v_prime

    real :: omega_p, omega_n
    integer ::  upwind_opt
    real :: tao_time 
    CHARACTER*1024     :: outstring
    CHARACTER*1024     :: wrf_err_message

    ! reference potential temperature profiles 
    integer, intent(in) :: nl_bg
    integer, intent(in) :: wavecoupling_opt
    integer, intent(in) :: wavecoupling_dmpopt ! option for damping: specify or constant or pressure depencance
    real, intent(in) :: wavecoupling_topbc  ! top of the wave coupling layer : meters
    real, intent(in) :: wavedamping, wavenumber
    real, DIMENSION( nl_bg), intent(in) :: theta_bg, qv_bg, z_bg, prs_bg 
    real, DIMENSION( nl_bg) ::  tmk_v_bg 
    real, DIMENSION(kms:kme) :: theta_B, qv_B  ! reference values


    !SW for mode based spectral wtg
    integer,intent(in) :: wtg_opt
    real, intent(in) :: wtg_cexp
    integer,intent(in) :: wtg_nmodes
    real,dimension(:),intent(in) :: scm_m_vector
    !real,dimension(1:num_force_layers,1:kme),intent(in) :: scm_m_vector
    real,dimension(1:kme),intent(in) :: scm_zm
    real,dimension(1:kme),intent(in) :: scm_m_value, scm_m_weight
    real,dimension(wtg_nmodes, 1:kme) :: modes_modellevels
    real, dimension(1:wtg_nmodes) ::  am, am_s ! for Spectral-WTG  
    real,dimension(1:kme) :: th_over_dthdz_int,  th_over_dthdz, rho_int
    real,dimension(1:kme, 1:wtg_nmodes) :: mode_tmp, mode_tmp_deWeight
    integer :: m
    integer :: nnlv


    ! local variables for wave coupling
    !real, DIMENSION( 1:kme-1) :: theta_bg_int, qv_bg_int , prs_bg_int, tmk_v_bg_int 
    real, DIMENSION( 1:kde) :: theta_bg_int, qv_bg_int , prs_bg_int, tmk_v_bg_int 
    real, DIMENSION( 1:kde) :: wavedamping_z, wavedamping_dz ! wave damping
    real, DIMENSION( 1:kde) :: th_v_avg
    integer :: np
    integer :: k2

    ! options for moisture advection
    integer,intent(in) :: qv_avg_flag
    integer,intent(in) :: qadv_flag, qadv_popt
    real,intent(in) :: qadv_tau, qadv_plev_up, qadv_plev_dn
    real,intent(in) :: qforcing_tau
    real, DIMENSION( kms:kme) :: qtend_qadv
    real, DIMENSION( kms:kme), intent(inout) :: wtg_qtend, wtg_qtend_v, wtg_qtend_h, wtg_ttend 
    real, DIMENSION( kms:kme) :: tmk_v 

    real :: cpm, evapor,rho_d,coef_tmp
    !qadv_flag = 1
    !qadv_opt = 0
    !qadv_tau = 1.0

    ! --------------------
    !qv_avg_flag=-1 ! q advection on each point
    !pblh_option = 2
    wtg_adv_opt = 2 ! advective form
    upwind_opt = 1 ! use upwind scheme
    qvint_tmp = 0.0
    qvhadvint_tmp = 0.0
    thetaint_tmp = 0.0
   

   phalf_avg = 0.0; 
   pfull_avg = 0.0
    
   DTwtg = 0.0
!SW==

!------------------------------------------------------------------
!
!
  flag_qv = .FALSE. ; IF ( PRESENT( F_QV ) ) flag_qv = F_QV
  flag_qc = .FALSE. ; IF ( PRESENT( F_QC ) ) flag_qc = F_QC

!SW==
  !if (bl_pbl_physics .eq. 0) return !SW==
!SW==
! RAINBL in mm (Accumulation between PBL calls)
!
! Modified for adaptive time step
!
  STEPwtg = nint(wtg_dt*60./DT)
  STEPwtg = max(STEPwtg,1)

  IF ( (itimestep .EQ. 1) .OR. (MOD(itimestep,STEPwtg) .EQ. 0) ) THEN
    run_param = .TRUE.
  ELSE
    run_param = .FALSE.
  ENDIF




  IF( wtg_dt >0 .and. adapt_step_flag ) then
!    CALL wrf_error_fatal('Error, CANNOT do wtg steps and adaptive time steps the same time!')
  ENDIF

  IF (PRESENT(adapt_step_flag)) THEN
    IF ((adapt_step_flag)) THEN
      IF ( (itimestep .EQ. 1) .OR. (wtg_dt .EQ. 0) .OR. &
           ( CURR_SECS + dt >= ( INT( CURR_SECS / ( wtg_dt * 60 ) + 1 ) * wtg_dt * 60) ) ) THEN
        run_param = .TRUE.
      ELSE
        run_param = .FALSE.
      ENDIF
    ENDIF
  ENDIF

 IF (run_param) THEN

!---- 
! CALCULATE CONSTANT
    
   DTMIN=DT/60.
! PBL schemes need PBL time step for updates

    if (PRESENT(adapt_step_flag)) then
       if (adapt_step_flag) then
          do_adapt = .TRUE.
       else
          do_adapt = .FALSE.
       endif
    else
       do_adapt = .FALSE.
    endif

      if (wtg_dt .eq. 0) then
         DTwtg = dt
      ELSE
         if (do_adapt) then
            call wrf_message("WARNING: When using an adaptive time-step the boundary layer"// &
                             " time-step should be 0 (i.e., equivalent to model time-step).  "// &
                             "In order to proceed, for boundary layer calculations, the "// &
                             "boundary layer time-step"// &
                             " will be rounded to the nearest minute, possibly resulting in"// &
                             " innacurate results.")
            DTwtg=wtg_dt*60
         else
            DTwtg=DT*STEPwtg
         endif
      endif
 ENDIF

if (rce_wtg_diag==1 )then
   print*,'itimestep', itimestep
   print*,'STEPBL',STEPBL
   print*,'STEPwtg',STEpwtg
   print*,'DT,DTwtg',DT,DTwtg
   print*,'run_wtg',run_param
endif

!

! SW==
    IF(relax_umean_opt==1)then
      no_points = float((ide-ids)*(jde-jds))
      DO k=kts,min(kte+1,kde)
         uuu_sum = 0.0
         vvv_sum = 0.0
         DO ij = 1 , num_tiles; 
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               uuu_sum = uuu_sum + u_phy(i,k,j)       
               vvv_sum = vvv_sum + v_phy(i,k,j)       
            ENDDO; ENDDO 
         ENDDO
         
         uuu_avg(k) = wrf_dm_sum_real ( uuu_sum )
         uuu_avg(k) = uuu_avg(k) / no_points
         vvv_avg(k) = wrf_dm_sum_real ( vvv_sum )
         vvv_avg(k) = vvv_avg(k) / no_points

        DO ij = 1 , num_tiles; 
           DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
                 RUBLTEN(i,k,j) = RUBLTEN(i,k,j) - (uuu_avg(k)-relax_umean_u)/(relax_umean_tscale*60.0)
                 RVBLTEN(i,k,j) = RVBLTEN(i,k,j) - (vvv_avg(k))/(relax_umean_tscale*60.0)
                 !RUBLTEN(i,k,j) = - (uuu_avg(k)-relax_umean_u)/(relax_umean_tscale*60.0)
           ENDDO; ENDDO
        ENDDO
      ENDDO         
   ENDIF

! rce well tested using YSU and no pbl
      DO k=kts,kde-1
          if(theta_bg(k) < 0) then
               WRITE(wrf_err_message,*) 'theta_bg(',k,') less than 0'
               CALL wrf_error_fatal(TRIM(wrf_err_message))
          endif
      ENDDO


      no_points = float((ide-ids)*(jde-jds))
      DO k=kts,kde-1
         z_sum = 0.0
         DO ij = 1 , num_tiles; 
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               z_2d(i,j) = z(i,k,j) 
               z_sum = z_sum + z_2d(i,j)
            ENDDO; ENDDO
         ENDDO
         z_avg(k) = wrf_dm_sum_real ( z_sum )
         z_avg(k) = z_avg(k) / no_points
      ENDDO
  
      
      DO k=kts,kde-1
      !DO k=kts,min(kte+1,kde)
         !omega_wtg(k)=0.0
         Dth_sum = 0.0
         rho_sum=0.0
         tmk_sum=0.0
         qv_sum = 0.0
         phalf_sum = 0.0
         pfull_sum = 0.0
         
         ql_sum = 0.0; qi_sum = 0.0; qr_sum=0.0; qs_sum=0.0; qg_sum=0.0; clf_sum=0.0; clh_sum=0.0; 
         mp_tten_sum = 0.0; mp_qten_sum=0.0; tot_tten_sum=0.0; tot_qten_sum=0.0; qrad_sum = 0.0
         areacloudupdraft_sum =0.0;  masscloudupdraft_sum = 0.0; areaupdraftcore_sum = 0.0;
         massupdraftcore_sum = 0.0;  areabuoyantcloud_sum = 0.0; massbuoyantcloud_sum = 0.0;
         rain_sum = 0.0; evap_sum=0.0; lhf_sum=0.0; hfx_sum = 0.0;

         Th_2d = 0.0
         rho_2d = 0.0
         phalf_2d = 0.0
         pfull_2d = 0.0
         DO ij = 1 , num_tiles; 
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               th_2d(i,j)=th_phy(i,k,j)
               dth_sum = dth_sum + th_2d(i,j)
               
               rho_2d(i,j)=rho(i,k,j)
               rho_sum = rho_sum + rho_2d(i,j)
                
               qv_2d(i,j)=qv_curr(i,k,j)
               qv_sum = qv_sum + qv_2d(i,j)
               
               phalf_2d(i,j)=p_phy(i,k,j)
               phalf_sum = phalf_sum + phalf_2d(i,j)

               pfull_2d(i,j)=p8w(i,k,j)
               pfull_sum = pfull_sum + pfull_2d(i,j)
         
               tmk_sum = tmk_sum + t_phy(i,k,j)
               ql_sum = ql_sum + qc_curr(i,k,j)
               qi_sum = qi_sum + qi_curr(i,k,j)
               qg_sum = qg_sum + qg_curr(i,k,j)
               qs_sum = qs_sum + qs_curr(i,k,j)
               qr_sum = qr_sum + qr_curr(i,k,j)
               mp_tten_sum = mp_tten_sum + rthmpten (i,k,j)*86400
               mp_qten_sum = mp_qten_sum + rqvmpten (i,k,j)*86400

              if(k==2)then
               rain_sum = rain_sum + rainnc(i,j)
               evap_sum = evap_sum + acqfx(i,j)
               hfx_sum = hfx_sum + achfx(i,j)
               lhf_sum = lhf_sum + aclhf(i,j)
              endif


              if(k>kms) then 
                 tot_tten_sum = tot_tten_sum + (rthblten (i,k,j) + &
                             (0.5*(rw(i,k,j)+rw(i,k+1,j))/mut(i,j)*th_phy(i,k,j) - 0.5*(rw(i,k,j)+rw(i,k-1,j))/mut(i,j)*th_phy(i,k-1,j))/ dz8w(i,k,j)   &
                              )*86400 !---------------------------
                 tot_qten_sum = tot_qten_sum + rqvblten (i,k,j)+ & !------------------------
                             (0.5*(rw(i,k,j)+rw(i,k+1,j))/mut(i,j)*qv_curr(i,k,j) - 0.5*(rw(i,k,j)+rw(i,k-1,j))/mut(i,j)*qv_curr(i,k-1,j))/ dz8w(i,k-1,j)   
              else
                 tot_tten_sum = tot_tten_sum + (rthblten (i,k,j) + &
                             (0.5*(rw(i,k,j)+rw(i,k+1,j))/mut(i,j)*th_phy(i,k,j) ) / dz8w(i,k,j)   &
                              )*86400 !---------------------------
                 tot_qten_sum = tot_qten_sum + rqvblten (i,k,j)+ & !------------------------
                             (0.5*(rw(i,k,j)+rw(i,k+1,j))/mut(i,j)*qv_curr(i,k,j) )/ dz8w(i,k,j)   

              endif


               qrad_sum = qrad_sum + rthraten (i,k,j)*86400 
 
               if( (qc_curr(i,k,j)+qi_curr(i,k,j)) > 1e-5 ) then
                   clf_sum = clf_sum + 1
                  if( rw(i,k,j)/mut(i,j) > 0) then
                      AreaCloudUpdraft_sum = AreaCloudUpdraft_sum + 1
                      MassCloudUpdraft_sum = MassCloudUpdraft_sum + rho(i,k,j)*rw(i,k,j)/mut(i,j)
                  endif
                  if( rw(i,k,j)/mut(i,j) > 1) then
                      AreaUpdraftCore_sum = AreaUpdraftCore_sum + 1
                      MassUpdraftCore_sum = MassUpdraftCore_sum + rho(i,k,j)*rw(i,k,j)/mut(i,j)
                  endif
               endif
               if( (qc_curr(i,k,j)+qi_curr(i,k,j)+qs_curr(i,k,j)+qg_curr(i,k,j)+qr_curr(i,k,j)) > 1e-5 ) then
                   clh_sum = clh_sum + 1
               endif
            ENDDO; ENDDO 
         ENDDO
         !DO j = jps, min(jpe,jde-1)
         !    DO i = ips, min(ipe,ide-1)
         !      dth_sum = dth_sum + th_2d(i,j)
         !    ENDDO
         !ENDDO
         th_avg(k) = wrf_dm_sum_real ( dth_sum )
         th_avg(k) = th_avg(k) / no_points

         tmk_avg(k) = wrf_dm_sum_real ( tmk_sum )
         tmk_avg(k) = tmk_avg(k) / no_points
         
         rho_avg(k) = wrf_dm_sum_real ( rho_sum )
         rho_avg(k) = rho_avg(k) / no_points

         pfull_avg(k) = wrf_dm_sum_real ( pfull_sum )
         pfull_avg(k) = pfull_avg(k) / no_points

         qv_avg(k) = wrf_dm_sum_real ( qv_sum )
         qv_avg(k) = qv_avg(k) / no_points

         phalf_avg(k) = wrf_dm_sum_real ( phalf_sum )
         phalf_avg(k) = phalf_avg(k) / no_points
         
         if(rho_avg(k)>2.and.k<kde-1) then
            print*,'k,rho',k,rho_avg(k),th_avg(k)
            stop
         endif
         qv_avg(k) = wrf_dm_sum_real ( qv_sum )
         qv_avg(k) = qv_avg(k) / no_points
         th_v_avg(k) = th_avg(k)*(1+0.608*qv_avg(k))
      
         ql_avg(k) = wrf_dm_sum_real ( ql_sum )/ no_points 
         qs_avg(k) = wrf_dm_sum_real ( qs_sum )/ no_points 
         qg_avg(k) = wrf_dm_sum_real ( qg_sum )/ no_points 
         qi_avg(k) = wrf_dm_sum_real ( qi_sum )/ no_points 
         qr_avg(k) = wrf_dm_sum_real ( qr_sum )/ no_points 
         phi_avg(k) = z_avg(k)*9.81
         !tmk_avg(k) = (th_avg(k)*(phalf_avg(k)/1e5)**(2.0/7.0)) 
         clf_avg(k) = wrf_dm_sum_real ( clf_sum )/ no_points 
         clh_avg(k) = wrf_dm_sum_real ( clh_sum )/ no_points 

         mp_tten_avg(k) = wrf_dm_sum_real ( mp_tten_sum )/ no_points 
         mp_qten_avg(k) = wrf_dm_sum_real ( mp_qten_sum )/ no_points 
         qrad_avg(k) = wrf_dm_sum_real ( qrad_sum )/ no_points 
         tot_tten_avg(k) = wrf_dm_sum_real ( tot_tten_sum )/ no_points 
         tot_qten_avg(k) = wrf_dm_sum_real ( tot_qten_sum )/ no_points 

       
         AreaCloudUpdraft_avg(k) = wrf_dm_sum_real ( AreaCloudUpdraft_sum ) /no_points 
         MassCloudUpdraft_avg(k) = wrf_dm_sum_real ( MassCloudUpdraft_sum ) /no_points 
         AreaUpdraftCore_avg(k) = wrf_dm_sum_real ( AreaUpdraftCore_sum ) /no_points 
         MassUpdraftCore_avg(k) = wrf_dm_sum_real ( MassUpdraftCore_sum ) /no_points 

         if(k==2)then
           rain_avg = wrf_dm_sum_real ( rain_sum )/ no_points 
           evap_avg = wrf_dm_sum_real ( evap_sum )/ no_points 
           hfx_avg = wrf_dm_sum_real ( hfx_sum )/ no_points 
           lhf_avg = wrf_dm_sum_real ( lhf_sum )/ no_points 
         endif
     ENDDO

     DO k=kts,kde-1
         areabuoyantcloud_sum = 0.0; massbuoyantcloud_sum = 0.0;
         DO ij = 1 , num_tiles; 
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               if( (qc_curr(i,k,j)+qi_curr(i,k,j)) > 1e-5 ) then
                  if( rw(i,k,j) > 0 .and. (th_phy(i,k,j)*(1+0.608*qv_curr(i,k,j))-th_v_avg(k) )>0 ) then
                      AreaBuoyantCloud_sum = AreaBuoyantCloud_sum + 1
                      MassBuoyantCloud_sum = MassBuoyantCloud_sum + rho(i,k,j)*rw(i,k,j)/mut(i,j)
                  endif
               endif
            ENDDO; ENDDO
         ENDDO 
         AreaBuoyantCloud_avg(k) = wrf_dm_sum_real ( AreaBuoyantCloud_sum ) /no_points 
         MassBuoyantCloud_avg(k) = wrf_dm_sum_real ( MassBuoyantCloud_sum ) /no_points
     ENDDO

   IF( False .and. run_param .and. ( rce_wtg_flag > 0 ) .and.(bl_pbl_physics==YSUSCHEME.or. bl_pbl_physics==0))then
     ! ---------------------------------------------------------------------------
 if( rce_wtg_flag < 10) then
     ! Interpolate theta_bg to coordinate levels:

     do k=1,kme-1
           theta_bg_int(k) =  interp_0(theta_bg, z_bg, z_avg(k), nl_bg)
           qv_bg_int(k) =  interp_0(qv_bg, z_bg, z_avg(k), nl_bg)
           qv_bg_int(k) = qv_bg_int(k)*1e-3 
     enddo

      if(rce_wtg_diag==1)then
        print*, 'nl_bg, kme', nl_bg, kme
        print*, 'z_avg = ', z_avg
        do k = 1, nl_bg
           print*, 'k, z_bg, th,qv', k, z_bg(k), theta_bg(k), qv_bg(k)
        enddo
      endif

      do k = kms, kme-1
         theta_B(k) = theta_bg_int(k)
         qv_B(k) = qv_bg_int(k)
         if(rce_wtg_diag==1.and.mod(itimestep,100)==1)print*,'theta_B,qv_B', k,z_avg(k), theta_B(k), qv_B(k)
      enddo
 
     ! ----------------------------------------------------------------------------

     call  coupling_largescale_omega   (                                  &
                 wtg_ttend, wtg_qtend                                     &
                 ,rce_wtg_flag                                            &
                 ,omega_wtg                                               &
                 ,wavecoupling_topbc                                      &
                 ,tao_time                                                & ! wtg parameters
                 ,real(pblh_opt)                                          & ! wtg parameters
                 ,wtg_opt,wtg_cexp, wtg_nmodes                            & !options for modal based WTG
                 ,scm_m_vector, scm_zm, scm_m_value,scm_m_weight          & ! wtg parameters
                 ,wavenumber, wavedamping                                 &
                 ,z_avg, phalf_avg, pfull_avg, rho_avg, th_avg,  qv_avg          &
                 ,theta_B, qv_B                                           &
                 ,qv_lateral_drawing_flag                                 &
                 ,scm_vert_adv_opt                                        &
                 ,itimestep        &
                 ,kms, kme                          )

     ! ----------------------------------------------------------------------------
 elseif( rce_wtg_flag == 11 .or. rce_wtg_flag == 12 )then

      do k = kms, kme
          prs_bg_int(k) = phalf_avg(k)
      enddo

      call interp_p( prs_bg*1e2, qv_bg, nl_bg, phalf_avg, qv_bg_int, kme)
      qv_bg_int = qv_bg_int*1e-3
      call interp_p( prs_bg*1e2, theta_bg, nl_bg, phalf_avg, theta_bg_int, kme)

      do k = kms, kme
         tmk_v(k) = (th_avg(k)*(phalf_avg(k)/1e5)**(2.0/7.0)) * (1+0.608*qv_avg(k))
         tmk_v_bg_int(k) = (theta_bg_int(k)*(prs_bg_int(k)/1e5)**(2.0/7.0)) * (1.0+0.608*qv_bg_int(k))
      enddo



     call WTG_DGW_P_COUPLING(                                             &
        wtg_dgw_flag = rce_wtg_flag,                                           &
        NP=kme,                                                                &
        PFLD = phalf_avg,    PFLDref = prs_bg_int,                             &
        QFLDref = qv_bg_int,  QFLD = qv_avg,                                   &
        THFLDref = theta_bg_int, THFLD = th_avg,                               &
        TvFLDref = tmk_v_bg_int, TvFLD = tmk_v,                                &
        WFLD = omega_wtg,                                                      &
        SQFLD = wtg_qtend,                                                     &
        SQFLD_v = wtg_qtend_v,                                                 &
        SQFLD_h = wtg_qtend_h,                                                 &
        STHFLD  = wtg_ttend,                                                   &
        THO_WTG = tao_time*60.0,  WTG_BL = real(pblh_opt*100),                         &    
        WAVE_NUM=wavenumber, WAVE_DAMP = wavedamping, DGW_TOP = wavecoupling_topbc       ) 

   
    if (rce_wtg_diag==1 )then
      print*,'wtg tau time (minutes)' , tao_time
      print*,'prs_bg=', prs_bg
      print*,'prs_bg_int=', prs_bg_int
      print*,'theta_bg_int=', theta_bg_int
      !print*,'prs_bg=', log(prs_bg/1e3)
      print*,'qv_bg_int=',qv_bg_int
      !print*,'theta, qv, tmk_v', k,theta_bg_int(k), qv_bg_int(k), tmk_v_bg_int(k), tmk_v(k), th_avg(k), prs_bg_int(k)/100.0
      print*,'omega_wtg=', omega_wtg(1:50)
      print*,'wtg_qtend=', wtg_qtend(1:50)
      print*,'wtg_qtend_h=', wtg_qtend_h(1:50)
      print*,'wtg_ttend=', wtg_ttend(1:50)
    endif
     ! ----------------------------------------------------------------------------
  endif

  ENDIF  ! end of big if block
             

     DO k = kds, kde-1
         DO ij = 1 , num_tiles; 
              DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
                 RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j) + wtg_qtend(k)
                 RTHBLTEN(i,k,j) = RTHBLTEN(i,k,j) + wtg_ttend(k)

                 qvint_tmp(i,j) = qvint_tmp(i,j) + qtend(k)*(znw(k+1)-znw(k))
                 cpm = cp * (1. + 0.8 * qv_curr(i,k,j))
                 evapor = p_phy(i,k,j)*r_v*qv_curr(i,k,j)/(r_d+r_v*qv_curr(i,k,j))
                 rho_d = (p_phy(i,k,j)-evapor)/t_phy(i,k,j)/r_d
            ENDDO; ENDDO
         ENDDO

         rqvwtgten(k) = wtg_qtend(k)
         rthwtgten(k) = wtg_ttend(k)

         if (rce_wtg_diag==1 )then       
           if(k==kts)then 
             WRITE(outstring,'(a10,   a12,   a12,   a9,   a13,   a15,   a10,   a10,   a10,   a10,  a10,  a10,  a10 )') &
               'ww_wtg','th_avg', 'th_v_prime','qv_avg', 'qv_prime',&
               '-ttend','-qtend', 'z_avg','z',  'pres'
             CALL wrf_message ( TRIM(outstring) )
           endif
           WRITE(outstring,'(G10.4, F12.2, F12.2, F12.4, g13.5, g13.5, g15.5, f12.5, f13.5, f13.5, f15.5, f13.5, f13.5, f13.5 )') &
             !omega_wtg(k), th_avg(k),th_rce(k),qv_avg(k)*1e3 &
             omega_wtg(k)/rho_avg(k)/(-9.81)*100.0, th_v_avg(k), th_v_avg(k)-theta_bg_int(k)*(1+0.608*qv_bg_int(k)), qv_avg(k)*1e3, (qv_avg(k)-qv_bg_int(k))*1e3 &
             !omega_wtg(k)/rho_avg(k)/(-9.81), th_avg(k),theta_bg_int(k),qv_avg(k)*1e3 &
             ,wtg_ttend(k)*86400, wtg_qtend(k)*1e6,  z_avg(k)/1e3,  z(1,k,1)/1e3, phalf_avg(k), prs_bg_int(k)

           CALL wrf_message ( TRIM(outstring) )
           if(k==kde-1)then 
             WRITE(outstring,'(a10,   a12,   a12,   a9,   a13,   a15,   a10,   a10,   a10,   a10 )') &
               'ww_wtg','th_avg', 'theta_bg','qv_avg',&
               'dthdp','qforcing','rhs','zh','z',  'pres'
             CALL wrf_message ( TRIM(outstring) )
           endif
         endif 
      ENDDO ! end of the big do loop 

     !diagnose

      DO ij = 1 , num_tiles; 
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
            ACTWTG(i,j) = ACTWTG(i,j)+thetaint_tmp(i,j)*mut(i,j)/g*DT  ! unit: mm, assume rhowater=1kg/m3
            ACQVWTG(i,j) = ACQVWTG(i,j)+qvint_tmp(i,j)*mut(i,j)/9.81*DT  ! unit: mm, assume rhowater=1kg/m3
            ACQVHADVWTG(i,j) = ACQVHADVWTG(i,j)+qvhadvint_tmp(i,j)*mut(i,j)/9.81*DT  ! unit: mm
        ENDDO
        ENDDO
     ENDDO


   
!SW==

!
 


   END SUBROUTINE wtg_driver4gass



!-------------------------------------------------------------------
!-------------------------------------------------------------------
 subroutine interp_p( prs_in, v_in, nz_in, prs_out, v_out, nz_out  )
 implicit none
 integer, intent(in) :: nz_in, nz_out
 real, dimension(1:nz_in),intent(in) ::    v_in(nz_in), prs_in(nz_in)
 real, dimension(1:nz_out),intent(inout) ::    v_out(nz_out), prs_out(nz_out)

 integer k, k2

     DO k = 1, nz_out
        k2  = nz_in
        DO WHILE( log(prs_out(k)/1e5) .gt. log(prs_in(k2)/1e5)  .and.  k2 .gt. 1 )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.nz_in)then
          v_out(k) = v_in(k2) + ( v_in(k2) - v_in(k2-1) )   &
                              * ( log(prs_out(k)/1e5) - log(prs_in(k2)/1e5)   )   &
                              / ( log(prs_in(k2)/1e5) - log(prs_in(k2-1)/1e5) )
        else
          v_out(k) = v_in(k2) + ( v_in(k2+1) - v_in(k2) )   &
                              * (      log(prs_out(k)/1e5) - log(prs_in(k2)/1e5) )   &
                              / ( log(prs_in(k2+1)/1e5) - log(prs_in(k2)/1e5) )
        endif
     ENDDO


 end subroutine interp_p

!-------------------------------------------------------------------
!-------------------------------------------------------------------

 subroutine interp_z( z_in, v_in, nz_in, z_out, v_out, nz_out  )
 implicit none
 integer, intent(in) :: nz_in, nz_out
 real, dimension(1:nz_in),intent(in) ::    v_in(nz_in), z_in(nz_in)
 real, dimension(1:nz_out),intent(inout) ::    v_out(nz_out), z_out(nz_out)

 integer k, k2

     DO k = 1, nz_out
        k2  = nz_in
        DO WHILE( z_out(k) .lt. z_in(k2)  .and.  k2 .gt. 1)
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.nz_in)then
          v_out(k) = v_in(k2) + ( v_in(k2) - v_in(k2-1) )   &
                              * ( z_out(k) - z_in(k2)   ) &
                              / ( z_in(k2) - z_in(k2-1) )
        else
          v_out(k) = v_in(k2) + ( v_in(k2+1) - v_in(k2) )   &
                              * ( z_out(k) - z_in(k2))   &
                              / ( z_in(k2+1) - z_in(k2) )
        endif
     ENDDO

 end subroutine interp_z





!-------------------------------------------------------------------
!-------------------------------------------------------------------

   SUBROUTINE coupling_largescale_omega   (                               &
                 ttend, qtend                                             & ! Tendencies of theta and qv
                 ,plargedynamics_flag                                     & ! Flag for large scale dynamics: 1, WTG; 2, DGW
                 ,omega_wtg                                               & ! Large scale pressure velocity
                 ,coupling_topbc                                          & ! Top level to apply large scale pressure velocity, e.g., 10000 Pa 
                 ,tau_time                                                & ! WTG time scale    (minutes)
                 ,wtg_pbl_height                                          & ! WTG height of PBL (Pa)
                 ,wtg_opt, wtg_cexp, wtg_nmodes                           & ! options for modal based WTG
                 ,scm_m_vector, scm_zm, scm_m_value,scm_m_weight          &
                 ,wavenumber, wavedamping                                 & ! DGW: wavenumber, and momentum damping time scale 
                 ,z_avg, p_avg, pfull_avg, rho_avg, th_avg, qv_avg        & ! Domain avarged pressure, density, theta, qvapor
                 ,theta_B, qv_B                                           & ! Reference values of theta (potential temperature) and qv (water vapor)
                 ,qv_lateral_drawing_flag                                 & ! Switch for lateral drawing by inflow
                 ,scm_vert_adv_opt                                        & ! Options to use different schemes of vertical advection 
                 ,itimestep                                               & 
                 ,kms, kme                          ) 

   ! Fortran subroutine for parameterizing large scale dynamics using two methods: WTG (Weak-temperature-gradient) and DGW (Damped-gravity-wave). 
   ! Prepared for the GASS-WTG project by Shuguang Wang (wangsg@gmail.com), Sep 2013.
   ! See technical specifications "GASS - WTG Inter-comparison Project", 2013, Shuguang Wang, Adam H. Sobel, Steven Woolnough, Sharon Sessions, and Gilles Bellon
 
   IMPLICIT NONE

   REAL    , PARAMETER :: r_d          = 287.                ! constant for WRF-ARW
   REAL    , PARAMETER :: cp           = 7.*r_d/2.           ! constant for WRF-ARW

   INTEGER,    INTENT(IN   )    ::   kms, kme   ! define dimension
   INTEGER,    INTENT(IN   )    ::  itimestep   

   real, dimension (kms:kme), intent(inout) :: ttend, qtend  
   ! ttend: tendency of potential temperature
   ! qtend: tendency of water vapor mixing ratio

   real, DIMENSION( kms:kme), intent(in) :: z_avg, p_avg, rho_avg, th_avg, qv_avg
   ! These are domain averaged quantities at the mass levels: z_avg, p_avg, pressure; rho_avg, density; th_avg: potential temperature; qv_avg: water vapor mixing ratio

   real, DIMENSION( kms:kme), intent(in) :: pfull_avg
   ! pfull_avg: domain averaged pressure at all the levels where W is defined (full pressure level in WRF)

   ! grid structure of p_avg, pfull_avg 
   !
   !   p_avg, omega_wtg, 
   !   th_avg, qv_avg ...         | pfull_avg
   !                              |
   !                              | kme:   xxxxxxxxxxxxxx     (top of model, w=0 )
   !    kme-1: ---------------    |                
   !                  .           | kme-1: xxxxxxxxxxxxxx              
   !                  .           |               .
   !     2:    ---------------    |               .
   !                              |  2:    xxxxxxxxxxxxxx 
   !     1:    ---------------    |
   !                              |  1:    xxxxxxxxxxxxxx     (at surface, w=0)
   !
 
    integer, intent(in) :: plargedynamics_flag ! 1: WTG, 2: DGW (Damped-gravity-wave coupling)
 
    real, intent(in) :: coupling_topbc       ! top of the wave coupling layer, unit: Pa 

    integer,intent(in) :: wtg_opt
    real, intent(in) :: tau_time ! WTG relaxation time scale, units: minutes
    real,intent(in)  :: wtg_pbl_height        ! Pressure level above which WTG applies (e.g., 850 hPa), unit: Pa
    real, DIMENSION( kms:kme),intent(inout) :: omega_wtg ! large scale pressure velocity from either WTG or DGW

    real, DIMENSION( kms:kme), intent(in) :: theta_B, qv_B 
    ! reference profiles of potential temperature and water vapor mixing ratio on the same levels of other domain averaged quantities

    ! variables for vertical modes based WTG
    real, intent(in) :: wtg_cexp
    integer,intent(in) :: wtg_nmodes
    real,dimension(:),intent(in) :: scm_m_vector
    !real,dimension(1:num_force_layers,1:kme),intent(in) :: scm_m_vector
    real,dimension(1:kme),intent(in) :: scm_zm
    real,dimension(1:kme),intent(in) :: scm_m_value, scm_m_weight
    real,dimension(wtg_nmodes, 1:kme) :: modes_modellevels
    real, dimension(1:wtg_nmodes) ::  am, am_s ! for Spectral-WTG  
    real,dimension(1:kme) :: th_over_dthdz_int,  th_over_dthdz, rho_int
    real,dimension(1:kme, 1:wtg_nmodes) :: mode_tmp, mode_tmp_deWeight
    integer :: m
    integer :: nnlv


   REAL, DIMENSION( kms:kme ) :: w_subs_modellevels, dthdz
   REAL, DIMENSION( kms:kme ) :: th_scm_target_modellevels

    real, intent(in) :: wavedamping ! 1.1574e-5 s^{-1} (1 day)
    real, intent(in) :: wavenumber  ! 1e-6  (m^{-1}) 

    logical,intent(in) :: qv_lateral_drawing_flag ! Switch for lateral drawing by inflow 

    ! local variables
    real, dimension(kms:kme) :: dthdp, dqvdp, theta_v_prime
    integer :: K
    integer :: k_abovepbl_first
    real, DIMENSION( kms:kme) :: dp, rhs, aa, bb, cc  ! variables for DGW
    real, DIMENSION( kms:kme) :: wavedamping_z, wavedamping_dz ! wave damping
    integer :: NP  ! number of layers for wave coupling
    real :: th_temp
    
    real, DIMENSION( kms:kme) :: qtend_lateral_drawing 

    real :: omega_p, omega_n, qtend1, qtend2, ttend1, ttend2, domgdp 

    real :: pi
    integer, intent(in) :: scm_vert_adv_opt

    CHARACTER*1024     :: err_message

    ! ------------------------------------------
    ! Default values                           | 
    ! coupling_topbc = 10000 (Pa)              |
    ! tau_time = 180 (min, 3 hours )           |
    ! wtg_pbl_height = 85000 (Pa)              |
    ! wavenumber = 0.25e-6    (m^-1)           |
    ! wavedamping = 1.1574e-5 (s^-1, 1 day)    |
    ! qv_lateral_drawing_flag = .True.         |  
    ! ------------------------------------------
 
    ! initialize the variables 
    DO k=kms, kme-1
        qtend(k) = 0.0
        ttend(k) = 0.0
        qtend_lateral_drawing(k) = 0.0
        omega_wtg(k) = 0.0 ! above the top boundary, omega_wtg is zero
    ENDDO 
    
    pi  =atan(1.0)*4.0

    ! find the top layer WTG/DWG is applied
    if (plargedynamics_flag == 1 .or. plargedynamics_flag == 2 )then
      DO k=kms,kme-1
         if(p_avg(k)<=coupling_topbc)then 
             NP=k-1
             exit 
         endif
      ENDDO 
    elseif (plargedynamics_flag == 3 .or. plargedynamics_flag == 4 )then
      DO k=kms,kme-1
         if(z_avg(k)>=coupling_topbc)then 
             NP=k-1
             exit 
         endif
      ENDDO 
    endif

   if(mod(itimestep,100)==1)print*, 'plargedynamics_flag=',plargedynamics_flag

   ! Solve large scale pressure velocity, omega_wtg, using WTG, vertical mode based WTG or DGW 
   plargedynamics_select: SELECT CASE ( plargedynamics_flag )  
   CASE ( 1 )  ! WTG

      DO k=NP,kms,-1
        omega_wtg(k) = 0.0
        if(p_avg(k) <=  wtg_pbl_height )then ! WTG above the nominal top of the boundary layer
          
           k_abovepbl_first = k

           if(k==kme-1)then
              dthdp(k) = (th_avg(k) - th_avg(k-1))/(p_avg(k) - p_avg(k-1))
           else
              dthdp(k) = 0.5*( (th_avg(k+1) - th_avg(k))/(p_avg(k+1) - p_avg(k)) + &
                              (th_avg(k) - th_avg(k-1))/(p_avg(k) - p_avg(k-1)) )
           endif
           if( dthdp(k) > 0) then
             dthdp(k) = max(dthdp(k), 2e-4)
           elseif( dthdp(k) < 0 ) then
             dthdp(k) = min(dthdp(k), -2e-4)
           endif
           omega_wtg(k) = (  th_avg(k)*(1+0.608*qv_avg(k)) - theta_B(k)*(1+0.608*qv_B(k))  )/(tau_time*60.0) / dthdp(k)
         else ! below boundary: linearly interpolation with respect to p with boundary condition W=0 at surface 
           omega_wtg(k) = (omega_wtg(k_abovepbl_first)-0.0 )/(p_avg(k_abovepbl_first)-pfull_avg(kms)) * (p_avg(k)-pfull_avg(kms))
        endif
      ENDDO
      if(mod(itimestep,100)==1)print*, 'tau-time (seconds)', tau_time*60.0
 
   CASE ( 3 )  ! vertical mode based WTG
       do k = kms, kme
         do m=1,wtg_nmodes
            mode_tmp(k,m) = scm_m_vector(m*kme+(k))
            mode_tmp_deWeight(k,m) = mode_tmp(k,m)/scm_m_weight(k)
         enddo
         if(k>=kms+1)then
                 dthdz(k) = ((th_avg(k)*(1+0.608*qv_avg(k))-th_avg(k-1)*(1+0.608*qv_avg(k-1)))/(z_avg(k)-z_avg(k-1))+ &
                             (th_avg(k+1)*(1+0.608*qv_avg(k+1))-th_avg(k)*(1+0.608*qv_avg(k)))/(z_avg(k+1)-z_avg(k)))*0.5
                 dthdz(k) = max(dthdz(k), 0.001)
          elseif(k==1)then ! kms = 1
                 dthdz(1) = ((th_avg(2)*(1+0.608*qv_avg(2))-th_avg(1)*(1+0.608*qv_avg(1)))/(z_avg(2)-z_avg(1))); dthdz(1) = max(dthdz(1), 0.001)
          endif
          th_scm_target_modellevels(k) = th_avg(k)*(1+0.608*qv_avg(k)) - theta_B(k)*(1+0.608*qv_B(k))  
          th_over_dthdz(k) = th_scm_target_modellevels(k)/dthdz(k)
       enddo

       do k = kms,kme
         if( k>10 .and. scm_zm(k) == 0 )then
               nnlv = k-1
               if(mod(itimestep,100)==1) print*,'nnlv=',nnlv
               exit
         endif
       enddo

      modes_modellevels = 0.0
      do k = kms,kme
         if(z_avg(k) <= coupling_topbc ) then
         do m=1,wtg_nmodes
            modes_modellevels(m,k) =interp_0(mode_tmp_deWeight(1:nnlv,m),scm_zm(1:nnlv), z_avg(k), nnlv)
         enddo
         endif
         th_over_dthdz_int(k) = interp_0(th_over_dthdz(1:kme-1),z_avg(1:kme-1), scm_zm(k), kme-1)
         rho_int(k) = interp_0(rho_avg(1:kme-1), z_avg(1:kme-1), scm_zm(k), kme-1)
         if(mod(itimestep,100)==1)print 10111, k,scm_zm(k), mode_tmp(k,1), mode_tmp(k,2), mode_tmp_deWeight(k,1),mode_tmp_deWeight(k,2), z_avg(k), modes_modellevels(1,k), modes_modellevels(2,k)
          10111 format(i10,2x, f10.4, 2x, f10.4, 2x, f10.4, 2x, f10.4, 2x, f10.4,2x, f10.4, 2x, f10.4,2x, f10.4,2x, f10.4,2x, f10.4 )
      enddo

      am = 0.0; am_s = 0.0;
      do m=1,wtg_nmodes ! first 10 modes
        do k = 1, nnlv
            am(m) = am(m) +  rho_int(k)*th_over_dthdz_int(k) *mode_tmp(k,m)*scm_m_weight(k) ! scm_m_weight: N, mode_tmp: N*Vn
            am_s(m) = am_s(m) + mode_tmp(k,m)*mode_tmp(k,m)
        enddo
      enddo
      if(mod(itimestep,100)==1)  then
         print*,'am=', am(1:wtg_nmodes)
         print*,'am_s=', am_s(1:wtg_nmodes)
      endif
      do k=kms, kme
        do m= 1,wtg_nmodes !first 10 modes
         if( z_avg(k) < coupling_topbc)then
             w_subs_modellevels(k) = w_subs_modellevels(k) + am(m)/am_s(m)*modes_modellevels(m,k) /tau_time/60.0*scm_m_value(m)**wtg_cexp
         endif
        enddo
      enddo
      do k=kms,kme
         omega_wtg(k) = -w_subs_modellevels(k)*9.81 ! use pressure velocity
      enddo

   CASE ( 4 ) ! spectral WTG, interpolate to force_layers first, calculate w_subs
      do k = kms+1,kme-1
         dthdz(k) = ((th_avg(k)*(1+0.608*qv_avg(k))-th_avg(k-1)*(1+0.608*qv_avg(k-1)))/(z_avg(k)-z_avg(k-1))+ &
                     (th_avg(k+1)*(1+0.608*qv_avg(k+1))-th_avg(k)*(1+0.608*qv_avg(k)))/(z_avg(k+1)-z_avg(k)))*0.5
         dthdz(k) = max(dthdz(k), 0.001)
      enddo
      dthdz(1) = ((th_avg(2)*(1+0.608*qv_avg(2))-th_avg(1)*(1+0.608*qv_avg(1)))/(z_avg(2)-z_avg(1))); dthdz(1) = max(dthdz(1), 0.001)

      w_subs_modellevels = 0.0
      do k = kms,kme-1
            th_scm_target_modellevels(k) = th_avg(k)*(1+0.608*qv_avg(k)) - theta_B(k)*(1+0.608*qv_B(k)) ! virtual temeprature
                                          !th_avg(k)-theta_B(k)
      enddo

      am = 0.0; am_s = 0.0;
      do m=1,wtg_nmodes ! first 10 modes
      do k = kms, kme-1
         if( z_avg(k) < coupling_topbc)then
            am(m) = am(m) + ( th_scm_target_modellevels(k)/dthdz(k)*sin(m*pi*z_avg(k)/coupling_topbc) + &
                              th_scm_target_modellevels(k+1)/dthdz(k+1)*sin(m*pi*z_avg(k+1)/coupling_topbc) &
                            )*0.5* (z_avg(k+1)-z_avg(k))
            am_s(m) = am_s(m) + ( sin(m*pi*z_avg(k)/coupling_topbc)**2 + sin(m*pi*z_avg(k+1)/coupling_topbc)**2 ) * 0.5* (z_avg(k+1)-z_avg(k))
         endif
      enddo
      enddo

      do k=kms, NP
        do m= 1,wtg_nmodes !first 10 modes
         if( z_avg(k) < coupling_topbc)then
             w_subs_modellevels(k) = w_subs_modellevels(k) + am(m)/am_s(m)*sin(m*pi*z_avg(k)/coupling_topbc)/tau_time/60.0/real(m)
         endif
        enddo
      enddo

      do k=kms,NP
         omega_wtg(k) = -w_subs_modellevels(k)*rho_avg(k)*9.81 ! use pressure velocity
      enddo

      if(mod(itimestep,100)==1)print*, 'NP=', 'am=',am
      if(mod(itimestep,100)==1)print*, 'ams=',am_s
      
   CASE ( 2 ) ! DGW 

      dp = 0.0;
      rhs = 0.0;
      aa = 0.0;
      bb = 0.0;
      cc = 0.0;

     ! ----------------------------------------------------------------------------
     ! aa, bb, cc are coefficients for the tridiagonal matrix for equation (2) in the technical specifications 
     ! Two swaps are used: the forward and the backward swaps
     ! compute grid spacing between midpoint levels (From Zhiming Kuang)

      do k=kms+1,NP
         dp(k)=p_avg(k)-p_avg(k-1)
      enddo
      dp(kms)=2*(p_avg(kms)-pfull_avg(kms))
      dp(np+1)=2*(pfull_avg(np+1)-p_avg(np))

      DO k = kms, kme 
           wavedamping_z (k) = wavedamping ! momentum damping coefficient 
           wavedamping_dz (k) = 0          ! independant of pressure or height 
      ENDDO

      ! Gauss Elimination
      do k=kms,NP
         theta_v_prime(k) = th_avg(k)*(1+0.608*qv_avg(k)) - theta_B(k)*(1+0.608*qv_B(k)) ! virtual potential temprature
         rhs(k)=wavenumber*wavenumber/wavedamping_z(k)*R_d*theta_v_prime(k)*(p_avg(k)/1.e5)**(R_d/Cp)/p_avg(k)*dp(k)*dp(k+1)*0.5
      enddo
     
      ! set up the tridiagonal matrix
      do k=kms,NP
         aa(k)=dp(k+1)/(dp(k)+dp(k+1))  -dp(k+1)*0.5*wavedamping_dz(k)/wavedamping_z(k)
         bb(k)=-1.  +(dp(k+1)-dp(k))*0.5*wavedamping_dz(k)/wavedamping_z(k)
         cc(k)=dp(k)/(dp(k)+dp(k+1))  +dp(k)*0.5*wavedamping_dz(k)/wavedamping_z(k)
      enddo
      
      ! symmetric lower BC
      aa(kms)=0.;
      bb(kms)=-(2*dp(kms+1)+dp(kms))/(dp(kms)+dp(kms+1))  -dp(kms)*0.5*wavedamping_dz(kms)/wavedamping_z(kms)
      cc(kms) = dp(k)/(dp(k)+dp(k+1))  +dp(kms)*0.5*wavedamping_dz(kms+1) /wavedamping_z(kms)
      ! symmetric upper BC
      aa(np)=dp(np+1)/(dp(np)+dp(np+1))  -dp(np+1)*0.5*wavedamping_dz(np)/wavedamping_z(np)
      bb(np) = -(2*dp(np)+dp(np+1))/(dp(np)+dp(np+1))  +dp(np+1)*0.5*wavedamping_dz(np)/wavedamping_z(np)
      cc(np) = 0.

      ! Gaussian Elimination without pivoting
      do k=1,NP-1
         bb(k+1)=bb(k+1)-aa(k+1)/bb(k)*cc(k)
         rhs(k+1)=rhs(k+1)-aa(k+1)/bb(k)*rhs(k)
      enddo
      ! Backward substitution
      rhs(np)=rhs(np)/bb(np)
      do k= NP-1,kms,-1
         rhs(k)=(rhs(k)-cc(k)*rhs(k+1))/bb(k)
      enddo
     
      do k=kms,NP
          omega_wtg(k)=rhs(k)  ! largescale pressure velocity
      enddo

   !ENDIF ! end of the IF block plargedynamics_flag
   END SELECT plargedynamics_select ! end of the IF block plargedynamics_flag
 
   ! ----------------------------------------------------------------------------

   if(qv_lateral_drawing_flag )then ! lateral drawing of moisture into the domain 

          do k=kms,NP
               domgdp=0
               if(k>kms.and.k<kme-1)then
                  domgdp = ( (omega_wtg(k)-omega_wtg(k-1))/(p_avg(k)-p_avg(k-1))+ &
                                  (omega_wtg(k+1)-omega_wtg(k))/(p_avg(k+1)-p_avg(k))  )*0.5
               elseif(k==kms)then
                  domgdp = (omega_wtg(k+1)-omega_wtg(k))/(p_avg(k+1)-p_avg(k))  
               elseif(k==kme-1)then
                  domgdp = (omega_wtg(k)-omega_wtg(k-1))/(p_avg(k)-p_avg(k-1))  
               endif
               if(domgdp>0)then ! inward  flow
                  qtend_lateral_drawing(k) = (qv_B(k)-qv_avg(k))*domgdp
               endif
          enddo
  
   endif
      
   ! -----------------------------------------------------------------------------------------
   ! Now we do vertical advection of theta and qv by omega_wtg using the 1st order upwind scheme
   DO k=kms,NP
         ! -Omega*dQv/dp, vertical advection of Qvapor by large scale vertical velocity
         ! -Omega*dT/dp, vertical advection of theta by large scale vertical velocity
    
     IF(  scm_vert_adv_opt == 1) then
         if(k>kms.and.k<kme-1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend2 = omega_p*(qv_avg(k+1) - qv_avg(k))/(p_avg(k+1)-p_avg(k)) 
                    qtend1 = omega_n*(qv_avg(k) - qv_avg(k-1)) /(p_avg(k)-p_avg(k-1))
                    qtend(k) = -(qtend1+qtend2)
                    ttend2 = omega_p*(th_avg(k+1) - th_avg(k))/(p_avg(k+1)-p_avg(k)) 
                    ttend1 = omega_n*(th_avg(k) - th_avg(k-1)) /(p_avg(k)-p_avg(k-1))
                    ttend(k) = - (ttend1+ttend2)
         elseif(k==kme-1)then ! at the top model level, redundant becasue np < kme-1
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend1 = omega_n*(qv_avg(k) - qv_avg(k-1)) /(p_avg(k)-p_avg(k-1))
                    qtend(k) = -qtend1
                    ttend1 = omega_n*(th_avg(k) - th_avg(k-1)) /(p_avg(k)-p_avg(k-1))
                    ttend(k) = - ttend1
         elseif(k==kms)then   ! first level
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    qtend2 = omega_p*(qv_avg(k+1) - qv_avg(k))/(p_avg(k+1)-p_avg(k)) 
                    qtend(k) = -qtend2
                    ttend2 = omega_p*(th_avg(k+1) - th_avg(k))/(p_avg(k+1)-p_avg(k)) 
                    ttend(k) = -ttend2
         endif
         ! -Omega * DQv/Dp, plus lateral drawing 
         qtend(k) = qtend(k) + qtend_lateral_drawing(k) 
                         
         if( plargedynamics_flag == 1)then ! for WTG, equavilent to do relaxation above pbl
              if( p_avg(k) <= wtg_pbl_height )then
                    ttend(k) =  -(  th_avg(k)*(1+0.608*qv_avg(k)) - theta_B(k)*(1+0.608*qv_B(k))  ) / (tau_time*60.0)
              endif
         endif

     ELSEIF (scm_vert_adv_opt == 800 .or. scm_vert_adv_opt == 0) then
         qtend(k) = 0.0
         ttend(k) = 0.0
         if(k>1.and.k<kme-1)then
              if( omega_wtg(k) > 0 )then
                    qtend(k) = -omega_wtg(k)*(qv_avg(k+1) - qv_avg(k))/(p_avg(k+1)-p_avg(k))
                    ttend(k) = -omega_wtg(k)*(th_avg(k+1) - th_avg(k))/(p_avg(k+1)-p_avg(k))
              elseif( omega_wtg(k) < 0 )then
                    qtend(k) = -omega_wtg(k)*(qv_avg(k) - qv_avg(k-1))/(p_avg(k)-p_avg(k-1))
                    ttend(k) = -omega_wtg(k)*(th_avg(k) - th_avg(k-1))/(p_avg(k)-p_avg(k-1))
              endif
         elseif(k==1)then
              if( omega_wtg(k) > 0 )then
                    qtend(k) = -omega_wtg(k)*(qv_avg(k+1) - qv_avg(k))/(p_avg(k+1)-p_avg(k))
                    ttend(k) = -omega_wtg(k)*(th_avg(k+1) - th_avg(k))/(p_avg(k+1)-p_avg(k))
              endif
         endif
     ENDIF

     if(qv_lateral_drawing_flag )  qtend(k) = qtend(k) +  qtend_lateral_drawing(k)
   ENDDO 
 
   END SUBROUTINE coupling_largescale_omega

!-------------------------------------------------------------------
!-------------------------------------------------------------------

SUBROUTINE WTG_DGW_P_COUPLING (wtg_dgw_flag,NP,PFLD,PFLDref,QFLDref, & 
                               QFLD,THFLDref,THFLD,TvFLDref,TvFLD,WFLD,SQFLD,SQFLD_v,SQFLD_h,STHFLD, & 
                               THO_WTG, WTG_BL, &
                               WAVE_NUM, WAVE_DAMP, DGW_TOP )

!This subroutine calculates the WTG/DGW-derived large-scale vertical velocities and the associated
!heat and moisture tendencies.
!The WTG/DGW-derived large-scale vertical velocities are at the same levels as qv

! QFLD, THFLD, TvFLD : domain-mean profiles of qv, temperature, potential temperature and
!virtual temperature respectively.
! WFLD: domain-mean profile of large-scale vertical velocity
! SQFLD_h: Lateral drawing (domain-mean) of qv by the diagnosed large-scale vertical velocity
! SQFLD_v: Vertical advection (domain-mean) of qv by the diagnosed large-scale vertical velocity
!STHFLD: Vertical advection (domain-mean) of potential temperature by the diagnosed large-scale vertical velocity (no lateral drwaing of TH)
! PFLD: Pressure levels are the same in the test and reference column for models in pressures coordinates
!wtg_dgw_flag ! flag for large scale dynamics: 1 for WTG and 2 for DGW
!NP: number of model levels in the vertical
! the subscript “ref” refer to the reference column.

 IMPLICIT NONE
 REAL, PARAMETER :: Rd=287.058 ! constant of dry air
 REAL, PARAMETER :: gra=9.81 ! acceleration due to gravity
 REAL, intent(in) :: THO_WTG, WTG_BL, DGW_TOP !WTG relaxation time scale (s)
 REAL, intent(in) :: wave_num, wave_damp ! wave damping
 ! THO_WTG=10800.!WTG relaxation time scale (s)
 ! wave_num=0.000001 ! wave number
 ! wave_damp=1./86400.! wave damping
 ! WTG_BL=85000. !85000 Pa:pressure below which WTG applies (Pa)
 ! DGW_TOP=10000.!Pressure above which DGW applies (Pa): it correspond to the nominate tropopause
 
 INTEGER, INTENT (IN) :: wtg_dgw_flag , NP
 REAL, INTENT(IN), DIMENSION(NP) :: THFLDref,THFLD,TvFLDref,TvFLD,&
 QFLDref,QFLD,PFLD,PFLDref

 REAL, INTENT(OUT), DIMENSION(NP) :: WFLD,STHFLD,SQFLD,SQFLD_v,SQFLD_h
 REAL,  DIMENSION(NP) :: STHFLDv,SQFLDv,SQFLDh
 INTEGER:: ipWTG_first_aboveBL! first model level above the nominal !boundary layer
 INTEGER:: ipDGW_first_aboveTOP ! first model level above the nominal !boundary layer
 REAL,DIMENSION(NP):: DPFLD ! Difference between two consecutives model levels
 REAL::DTHvrefDPFLD ! static stability of the reference column
 REAL:: DWDPFLD           ! Derivative of large-scale vertical velocity which respect to pressure
 REAL, DIMENSION(NP)::RHS,aa,bb,cc ! variables for the DGW formulations
 REAL:: const
 INTEGER:: ip !vertical model levels. 

!Grid structure 
!WFLD is defined at THETA levels (xxxxxxxxx)
!WFLD(ip=1)-0; where ip=1 is the model level just below the surface. 
!The first model level just above the surface is counted as model level 2 and the model level at the top of the domian is the model level NP
!
!                           NP: -----------------(top of the model)
!     NP: xxxxxxxxxxxxxxx 
!                         NP-1: ----------------
!   NP-1: xxxxxxxxxxxxxxx
!
!
!
!
!                            3: ----------------
!     3: xxxxxxxxxxxxxxx
!                            2: ---------------
!     2: xxxxxxxxxxxxxxx (first model level above the surface)
!                            1: ---------------(surface)
!     1: xxxxxxxxxxxxxx (WFLD =0)

 STHFLD = 0.0
 SQFLD = 0.0
 SQFLD_v = 0.0
 SQFLD_h = 0.0
 WFLD = 0.0
 STHFLDv = 0.0
 SQFLDv = 0.0
 RHS = 0.0

 DPFLD(1)=0.
 DO ip=2,NP
      DPFLD(ip)=PFLD(ip)-PFLD(ip-1)
 ENDDO

IF (wtg_dgw_flag .EQ. 11) THEN         !WTG formulation
      DO ip=NP-1,1,-1
         IF (PFLD(ip) .LE. WTG_BL) THEN ! WTG applies above the nominal top of the boundary layer
            ipWTG_first_aboveBL=ip
            IF (ip.EQ.NP-1) THEN  ! top of the domain
               DTHvrefDPFLD=(THFLDref(ip)*(1.0+0.608*QFLDref(ip))- &
               THFLDref(ip-1)*(1.0+0.608*QFLDref(ip-1)))/DPFLD(ip)
            ENDIF
            IF (ip.LT.NP-1) THEN
               DTHvrefDPFLD=0.5*((THFLDref(ip+1)*(1.0+0.608* &
               QFLDref(ip+1))-THFLDref(ip)*(1.0+0.608*QFLDref(ip)))/  &
               DPFLD(ip+1)+(THFLDref(ip)*(1.0+0.608*QFLDref(ip))- &
               THFLDref(ip-1)*(1.0+0.608*QFLDref(ip-1)))/DPFLD(ip))
            ENDIF
            IF (DTHvrefDPFLD.GT.0) THEN
               DTHvrefDPFLD= MAX(DTHvrefDPFLD,2e-4) 
            ENDIF
            IF (DTHvrefDPFLD.LT.0) THEN
               DTHvrefDPFLD= MIN(DTHvrefDPFLD,-2e-4)
            ENDIF
            WFLD(ip)=(THFLD(ip)*(1.0+0.608*QFLD(ip))-THFLDref(ip)* &
            (1.0+0.608*QFLDref(ip)))/(THO_WTG*(DTHvrefDPFLD))
         ELSE                   ! Below the nominal boundary layer: linear interpolation which respect to pressure from the first model level above the boundary layer to the value of 0m/s at the first modle level
            WFLD(ip)=WFLD(ipWTG_first_aboveBL )*(PFLD(ip)-PFLD(1))/&
            (PFLD(ipWTG_first_aboveBL )-PFLD(1))
         ENDIF
      ENDDO
      !print*,'WTG_BL=',wtg_bl, ipWTG_first_aboveBL
      !print*,WFLD
 !ELSE IF (wtg_dgw_flag .EQ. 11) THEN         !WTG formulation



 ELSE IF (wtg_dgw_flag .EQ. 12) THEN !DGW formulation
      DO ip=NP,1,-1
         IF (PFLD(ip) .GE. DGW_TOP) THEN
            ipDGW_first_aboveTOP=ip
            exit 
         ENDIF
      ENDDO

      const=wave_num*wave_num*Rd/wave_damp
!GUASS elimination
      DO ip=1,ipDGW_first_aboveTOP
         RHS(ip)=0.5*const*(TvFLD(ip)-TvFLDref(ip))* &
         DPFLD(ip)*DPFLD(ip+1)/PFLDref(ip)
      ENDDO
!set up the triangular matrix
      DO ip=1,ipDGW_first_aboveTOP
         aa(ip)=DPFLD(ip+1)/(DPFLD(ip+1)+DPFLD(ip))
         bb(ip)=-1.
         cc(ip)=DPFLD(ip)/(DPFLD(ip+1)+DPFLD(ip))
      ENDDO
!symmetric lower boundary
      aa(1)=0.
      bb(1)=-(2.*DPFLD(2)+DPFLD(1))/(DPFLD(1)+DPFLD(2))
      cc(1)=DPFLD(1)/(DPFLD(1)+DPFLD(2))
!symmetric upper boundary
      aa(ipDGW_first_aboveTOP)=DPFLD(ipDGW_first_aboveTOP+1)/ &
      (DPFLD(ipDGW_first_aboveTOP)+DPFLD(ipDGW_first_aboveTOP+1))
      bb(ipDGW_first_aboveTOP)=-(2.*DPFLD(ipDGW_first_aboveTOP)+ &
      DPFLD(ipDGW_first_aboveTOP+1))/(DPFLD(ipDGW_first_aboveTOP)+ &
      DPFLD(ipDGW_first_aboveTOP+1))
      cc(ipDGW_first_aboveTOP)=0.
!GUASSIAN elimination without pivoting
      DO ip=2,ipDGW_first_aboveTOP-1
         bb(ip+1)=bb(ip+1)-aa(ip+1)*cc(ip)/bb(ip)
         RHS(ip+1)=RHS(ip+1)-aa(ip+1)*RHS(ip)/bb(ip)
      ENDDO
!backward substitution
      RHS(ipDGW_first_aboveTOP)=RHS(ipDGW_first_aboveTOP)/ &
      bb(ipDGW_first_aboveTOP)
      DO ip=ipDGW_first_aboveTOP-1,1,-1
         RHS(ip)=(RHS(ip)-cc(ip)*RHS(ip+1))/bb(ip)
      ENDDO
      DO ip=1,ipDGW_first_aboveTOP
         WFLD(ip)=RHS(ip)
      ENDDO
 ENDIF
 DO ip=1,NP-1
      IF (ip.EQ.1) THEN
         DWDPFLD=(WFLD(ip+1)-WFLD(ip))/DPFLD(ip+1)
      ENDIF
      IF (ip.EQ.NP-1) THEN
         DWDPFLD=(WFLD(ip)-WFLD(ip-1))/DPFLD(ip)
      ENDIF
      IF (ip.GT.1.AND. ip.LT.NP-1) THEN
         DWDPFLD=0.5*((WFLD(ip)- WFLD(ip-1))/(DPFLD(ip))+ (WFLD(ip+1)- &
         WFLD(ip))/DPFLD(ip+1))
      ENDIF
      IF (DWDPFLD.GT.0) THEN
         SQFLDh(ip)=DWDPFLD*(QFLDref(ip)-QFLD(ip))
      ELSE
         SQFLDh(ip)=0.
      ENDIF
      SQFLD_h(ip)=         SQFLDh(ip)   
 ENDDO
!print*,'SFLDh', SQFLDh

!Vertical advection: -w*dx/dp with x=TH or qv
      DO ip=1,NP-1
         IF (ip.EQ.1 .AND. WFLD(ip).GT.0.0) THEN
            STHFLDv(ip)=-WFLD(ip)*(THFLD(ip+1)-THFLD(ip))/DPFLD(ip+1) 
            SQFLDv(ip)=-WFLD(ip)*(QFLD(ip+1)-QFLD(ip))/DPFLD(ip+1) 
         ENDIF
         IF (ip .EQ. NP-1 .AND. WFLD(ip).LT. 0.0) THEN
            STHFLDv(ip)=-WFLD(ip)*(THFLD(ip)-THFLD(ip-1))/ DPFLD(ip)
            SQFLDv(ip)=-WFLD(ip)*(QFLD(ip)-QFLD(ip-1))/DPFLD(ip) 
         ENDIF
         IF (ip.GT.1 .AND. ip .LT.NP-1) THEN
            IF (WFLD(ip).GT.0.0) THEN
               STHFLDv(ip)=-WFLD(ip)*(THFLD(ip+1)-THFLD(ip))/DPFLD(ip+1)
               SQFLDv(ip)=-WFLD(ip)*(QFLD(ip+1)-QFLD(ip))/DPFLD(ip+1) 
            ENDIF
            IF (WFLD(ip).LT.0.0) THEN
               STHFLDv(ip)=-WFLD(ip)*(THFLD(ip)-THFLD(ip-1))/DPFLD(ip)
               SQFLDv(ip)=-WFLD(ip)*(QFLD(ip)-QFLD(ip-1))/DPFLD(ip)
            ENDIF
         ENDIF
         STHFLD(ip) =   STHFLDv(ip)
         SQFLD_v(ip) =  SQFLDv(ip)
         SQFLD(ip)=     SQFLD_v(ip)+SQFLD_h(ip)
      ENDDO
       
      IF( wtg_dgw_flag .EQ. 1) THEN !D
            DO ip=ipWTG_first_aboveBL,NP-1
              STHFLD(ip) =   -(  THFLD(ip)*(1+0.608*QFLD(ip)) - THFLDref(ip)*(1+0.608*QFLDref(ip))  ) / THO_WTG 
            ENDDO
      ENDIF

 RETURN
 END SUBROUTINE WTG_DGW_P_COUPLING 

!-------------------------------------------------------------------
!-------------------------------------------------------------------

 real function interp_0( v_in,  &
                         z_in, z_out, nz_in  )
 implicit none
 integer nz_in, nz_out
 real    v_in(nz_in), z_in(nz_in)
 real    z_out

 integer kp, k, im, ip
 logical interp, increasing_z 
 real    height, w1, w2
 logical debug
 parameter ( debug = .false. )

! does vertical coordinate increase or decrease with increasing k?
! set offset appropriately

 height = z_out

 if(debug) write(6,*) ' height in interp_0 ',height

 if (z_in(nz_in) .gt. z_in(1)) then

    if(debug) write(6,*) ' monotonic increase in z in interp_0 '
    IF (height > z_in(nz_in)) then
      if(debug) write(6,*) ' point 1 in interp_0 '
      w2 = (z_in(nz_in)-height)/(z_in(nz_in)-z_in(nz_in-1))
      w1 = 1.-w2
      interp_0 = w1*v_in(nz_in) + w2*v_in(nz_in-1)
    ELSE IF (height < z_in(1)) then
      if(debug) write(6,*) ' point 2 in interp_0 '
      w2 = (z_in(2)-height)/(z_in(2)-z_in(1))
      w1 = 1.-w2
      interp_0 = w1*v_in(2) + w2*v_in(1)
    ELSE
      if(debug) write(6,*) ' point 3 in interp_0 '
      interp = .false.
      kp = nz_in
      DO WHILE ( (interp .eqv. .false.) .and. (kp .ge. 2) )
        IF(   ((z_in(kp)   .ge. height) .and.     &
               (z_in(kp-1) .le. height))        )   THEN
          w2 = (height-z_in(kp))/(z_in(kp-1)-z_in(kp))
          w1 = 1.-w2
          interp_0 = w1*v_in(kp) + w2*v_in(kp-1)
          if(debug) write(6,*) ' interp data, kp, w1, w2 ',kp, w1, w2
          if(debug) write(6,*) ' interp data, v_in(kp), v_in(kp-1), interp_0 ', &
                     v_in(kp), v_in(kp-1), interp_0
          interp = .true.
        END IF
        kp = kp-1
      ENDDO
    ENDIF

 else

    if(debug) write(6,*) ' monotonic decrease in z in interp_0 '

    IF (height < z_in(nz_in)) then
      if(debug) write(6,*) ' point 1 in interp_0 '
      w2 = (z_in(nz_in)-height)/(z_in(nz_in)-z_in(nz_in-1))
      w1 = 1.-w2
      interp_0 = w1*v_in(nz_in) + w2*v_in(nz_in-1)
    ELSE IF (height > z_in(1)) then
      if(debug) write(6,*) ' point 2 in interp_0 '
      w2 = (z_in(2)-height)/(z_in(2)-z_in(1))
      w1 = 1.-w2
      interp_0 = w1*v_in(2) + w2*v_in(1)
    ELSE
      if(debug) write(6,*) ' point 3 in interp_0 '
      interp = .false.
      kp = nz_in
      height = z_out
      DO WHILE ( (interp .eqv. .false.) .and. (kp .ge. 2) )
        IF(   ((z_in(kp)   .le. height) .and.     &
               (z_in(kp-1) .ge. height))             )   THEN
          w2 = (height-z_in(kp))/(z_in(kp-1)-z_in(kp))
          w1 = 1.-w2
          interp_0 = w1*v_in(kp) + w2*v_in(kp-1)
          interp = .true.
        END IF
        kp = kp-1
      ENDDO
    ENDIF

 end if

 return
 END FUNCTION interp_0

!-------------------------------------------------------------------
!-------------------------------------------------------------------

END MODULE module_wtg_driver4gass

