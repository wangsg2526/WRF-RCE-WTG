!WRF:MEDIATION_LAYER:PHYSICS
!

MODULE module_wtg_driver0
! This is legacy code used by Wang and Sobel 2011, 2012
! for WTG computation in eta coordinate, and DGW in pressure coordinate
! rce_wtg_flag = 101, 102, 103
! Shuguang Wang (wangsg@gmail.com)

CONTAINS

!------------------------------------------------------------------
   SUBROUTINE wtg_driver0(                                          &
                  itimestep,dt,u_frame,v_frame                     &
                 ,bldt,curr_secs,adapt_step_flag                   &
                 ,rublten,rvblten,rthblten                         &
                 ,tsk                                 &
!SW==
                 ,omega_wtg, the_wtg_forcing, qv_wtg_forcing       &
                 ,qv_wtg_hadv_forcing                              &
                 , qv_avg_flag                                     &                       
                 , qadv_flag, qadv_tau,qforcing_tau, qadv_plev_up,qadv_plev_dn, qadv_popt        &
                 ,th_rce,  rce_wtg_flag, tao_time, rce_wtg_diag    &
                 ,qv_rce, qv_hadv_flag                             &
                 ,pblh_opt, relax_umean_opt                        &
                 ,relax_umean_u, relax_umean_tscale                &
                 ,rqvwtgten,rthwtgten                              &
                 ,acqvwtg, acqvhadvwtg                             &
                 ,actwtg                              &
                 ,scm_vert_adv_opt                &
! shear variables
                 ,shear_opt, shear_relax_tscale       &
                 ,shear_u_top, shear_u_sfc, shear_depth, shear_bot &
                 ,shear_u_profile &
! for wave coupling
                 ,wavecoupling_opt                   &
                 ,wave_ztop                   &
                 ,nl_bg, theta_bg, qv_bg, z_bg              &
                 ,wavenumber, wavedamping            &
! for vertical mode based spectral wtg
                 ,wtg_opt,wtg_cexp, wtg_nmodes  & 
                 ,scm_m_vector, scm_zm, scm_m_value,scm_m_weight &
!SW==
                 ,pblh                          &
                 ,u_phy,v_phy,th_phy,rho                           &
                 ,p_phy,pi_phy,p8w,t_phy,dz8w,z                    &
                 ,bl_pbl_physics, ra_lw_physics, dx                &
                 ,stepbl                                &
                 ,znu, znw, mut, p_top                             &
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,i_start,i_end, j_start,j_end, kts,kte, num_tiles &
             !  Optional moisture tracers
                 ,qv_curr, qc_curr, qr_curr                        &
                 ,qi_curr, qs_curr, qg_curr                        &
                 ,rqvblten,rqcblten,rqiblten                       &
                 ,rqrblten,rqsblten,rqgblten                       &
             !  Optional moisture tracer flags
                 ,f_qv,f_qc,f_qr                                   &
                 ,f_qi,f_qs,f_qg                                   &
                                                                     )       
!------------------------------------------------------------------
   USE module_state_description, ONLY :      YSUSCHEME

   USE module_model_constants

!SW==
    USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
!SW==
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
!-- PSFC          pressure at the surface (Pa)
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

   INTEGER,    INTENT(IN   )    ::     ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       kts,kte, num_tiles

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                   &
  &                                    i_start,i_end,j_start,j_end

   INTEGER,    INTENT(IN   )    ::     itimestep,STEPBL
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
                                                          th_phy
!
!
   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               INTENT(INOUT)    ::                           TSK, &
                                                            PBLH

!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(INOUT)    ::                       RUBLTEN, &
                                                         RVBLTEN, &
                                                        RTHBLTEN

   REAL,       DIMENSION( ims:ime, jms:jme ),             &
              INTENT(INOUT)    ::        acqvwtg, acqvhadvwtg, actwtg

   REAL ,                             INTENT(IN   )  ::  u_frame, &
                                                         v_frame
   integer, intent(in) :: scm_vert_adv_opt
!
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
                                                      f_qv      &
                                                     ,f_qc      &
                                                     ,f_qr      &
                                                     ,f_qi      &
                                                     ,f_qs      &
                                                     ,f_qg

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                      ! optional moisture tracers
                      ! 2 time levels; if only one then use CURR
                      qv_curr, qc_curr, qr_curr                  &
                     ,qi_curr, qs_curr, qg_curr                  &
                     ,rqvblten,rqcblten,rqrblten                 &
                     ,rqiblten,rqsblten,rqgblten

   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               OPTIONAL                                         , &
               INTENT(IN)    ::                           mut
!
   REAL,       OPTIONAL, INTENT(IN)    ::               p_top

!  LOCAL  VAR

   REAL,       DIMENSION( ims:ime, jms:jme )          ::  PSFC
   REAL    :: DTMIN,DTBL
!

   INTEGER :: initflag
!

   INTEGER :: i,J,K,NK,jj,ij,its,ite,jts,jte
   LOGICAL :: radiation
   LOGICAL :: flag_qv, flag_qc, flag_qr, flag_qi, flag_qs, flag_qg
   CHARACTER*256 :: message
   REAL    :: next_bl_time
   LOGICAL :: run_param
   LOGICAL :: do_adapt

!SW==
! shear variables
  real, intent(in) :: shear_relax_tscale
  real, dimension(kms:kme), intent(out) :: shear_u_profile ! local variable
  real, intent(in) :: shear_u_top, shear_u_sfc, shear_depth, shear_bot
  integer, intent(in) :: shear_opt

   REAL, DIMENSION( kms:kme),                 &
         OPTIONAL, INTENT(INOUT) ::  rqvwtgten, rthwtgten

    real, dimension(kms:kme), INTENT(INOUT) ::                  &
                      the_wtg_forcing, qv_wtg_forcing,  qv_wtg_hadv_forcing
    real ::  no_points, th_temp, uuu_sum, vvv_sum
    real ::  dth_sum, rho_sum, qv_sum, phalf_sum, pfull_sum, z_sum, mu_sum, mu_avg

    integer :: irnc, jrnc
    integer :: k_abovepbl_first
    real    :: pblh_max, pblh_sum, pblh_domain
    integer, INTENT(IN) :: rce_wtg_flag, rce_wtg_diag, qv_hadv_flag 
    integer,optional,INTENT(IN) ::  pblh_opt, relax_umean_opt
    real,optional,INTENT(IN) ::   relax_umean_u, relax_umean_tscale
    integer :: wtg_adv_opt
    real, DIMENSION( ims:ime , jms:jme ) :: th_2d,rho_2d, qv_2d, z_2d, pfull_2d, phalf_2d

    real, DIMENSION( ims:ime , jms:jme ) :: thetaint_tmp,qvint_tmp, qvhadvint_tmp 
    real, DIMENSION( kms:kme) :: th_avg, rho_avg, qv_avg, z_avg
    real, DIMENSION( kms:kme),INTENT(IN) :: qv_rce,th_rce
    real, DIMENSION( kms:kme) :: uuu_avg, vvv_avg
    real, DIMENSION( kms:kme) :: w_wtg, omega_wtg
    real, DIMENSION( kms:kme) :: dthdeta, dqvdeta, qtend, qtend_qvhadv, qtend_avg, ttend
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
    real, intent(in) :: wave_ztop
    real, intent(in) :: wavedamping, wavenumber
    real, DIMENSION( nl_bg), intent(in) :: theta_bg, qv_bg, z_bg

!SW for vertical mode based WTG 
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

   REAL, DIMENSION(kms:kme)     :: dthdz,  dqvdz
   !REAL, DIMENSION( kms:kme ),INTENT (INOUT) :: th_scm_prime
   REAL, DIMENSION( kms:kme ) :: w_subs_modellevels
   REAL, DIMENSION( kms:kme ) :: th_scm_target_modellevels


    ! local variables for wave coupling
    real, DIMENSION( kms:kme) :: dp, rhs, aa, bb, cc
    real, DIMENSION( kms:kme) :: phalf_avg, pfull_avg
    real, DIMENSION( 1:kde) :: theta_bg_int, qv_bg_int
    real, DIMENSION( 1:kde) :: th_v_avg
    integer :: np
    integer :: k2

    ! options for moisture advection
    integer,intent(in) :: qv_avg_flag
    integer,intent(in) :: qadv_flag, qadv_popt
    real,intent(in) :: qadv_tau, qadv_plev_up, qadv_plev_dn
    real,intent(in) :: qforcing_tau
    real, DIMENSION( kms:kme) :: qtend_qadv

    real :: cpm, evapor,rho_d,coef_tmp
    !qadv_flag = 1
    !qadv_opt = 0
    !qadv_tau = 1.0

    ! --------------------
    !qv_avg_flag=-1 ! q advection on each point
    omega_wtg = 0.0
    w_wtg = 0.0 
    !pblh_option = 2
    wtg_adv_opt = 2 ! advective form
    upwind_opt = 1 ! use upwind scheme
    qvint_tmp = 0.0
    qvhadvint_tmp = 0.0
    thetaint_tmp = 0.0
    
!SW==

!------------------------------------------------------------------
!
!
  flag_qv = .FALSE. ; IF ( PRESENT( F_QV ) ) flag_qv = F_QV
  flag_qc = .FALSE. ; IF ( PRESENT( F_QC ) ) flag_qc = F_QC
  flag_qr = .FALSE. ; IF ( PRESENT( F_QR ) ) flag_qr = F_QR
  flag_qi = .FALSE. ; IF ( PRESENT( F_QI ) ) flag_qi = F_QI
  flag_qs = .FALSE. ; IF ( PRESENT( F_QS ) ) flag_qs = F_QS
  flag_qg = .FALSE. ; IF ( PRESENT( F_QG ) ) flag_qg = F_QG

!SW==
  !if (bl_pbl_physics .eq. 0) return !SW==
!SW==
!
! Modified for adaptive time step
!
  IF ( (itimestep .EQ. 1) .OR. (MOD(itimestep,STEPBL) .EQ. 0) ) THEN
    run_param = .TRUE.
  ELSE
    run_param = .FALSE.
  ENDIF

  IF (PRESENT(adapt_step_flag)) THEN
    IF ((adapt_step_flag)) THEN
      IF ( (itimestep .EQ. 1) .OR. (bldt .EQ. 0) .OR. &
           ( CURR_SECS + dt >= ( INT( CURR_SECS / ( bldt * 60 ) + 1 ) * bldt * 60) ) ) THEN
        run_param = .TRUE.
      ELSE
        run_param = .FALSE.
      ENDIF
    ENDIF
  ENDIF

 IF (run_param) THEN
  radiation = .false.
  IF (ra_lw_physics .gt. 0) radiation = .true.

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

   if (PRESENT(BLDT)) then
      if (bldt .eq. 0) then
         DTBL = dt
      ELSE
         if (do_adapt) then
            call wrf_message("WARNING: When using an adaptive time-step the boundary layer"// &
                             " time-step should be 0 (i.e., equivalent to model time-step).  "// &
                             "In order to proceed, for boundary layer calculations, the "// &
                             "boundary layer time-step"// &
                             " will be rounded to the nearest minute, possibly resulting in"// &
                             " innacurate results.")
            DTBL=bldt*60
         else
            DTBL=DT*STEPBL
         endif
      endif
   else
      DTBL=DT*STEPBL
   endif
 ENDIF
! SAVE OLD VALUES

!

 
! SW==
   IF(shear_opt ==1)then
      IF(relax_umean_opt==1)then
        CALL wrf_message("With shear, relax_umean_opt=0")
        stop
      endif

      no_points = float((ide-ids)*(jde-jds))
      DO k=kts,kde-1
         z_sum = 0.0
         z_2d = 0.0
         DO ij = 1 , num_tiles; 
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               z_2d(i,j) = z(i,k,j) 
               z_sum = z_sum + z_2d(i,j)
                 !if (rce_wtg_diag==1.and.i==1.and.j==1)then        
                 !    WRITE(outstring,*) 'i,j,k,z,pblh,=',i,j,k,z(i,k,j),pblh(i,j)
                 !    CALL wrf_message ( TRIM(outstring) )
                 !endif
            ENDDO; ENDDO
         ENDDO
         z_avg(k) = wrf_dm_sum_real ( z_sum )
         z_avg(k) = z_avg(k) / no_points
      ENDDO
     
      DO k=kts,min(kte+1,kde)
          
         if(z_avg(k)<=shear_bot)then
            shear_u_profile (k) = shear_u_sfc
         elseif(z_avg(k)>shear_bot.and.z_avg(k)<=(shear_bot+shear_depth))then
            shear_u_profile (k) = (shear_u_top-shear_u_sfc)/(shear_depth)*(z_avg(k)-shear_bot) + shear_u_sfc
         else
            shear_u_profile (k) = shear_u_top
         endif

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
                 RUBLTEN(i,k,j) = RUBLTEN(i,k,j) - (uuu_avg(k)-shear_u_profile (k) )/(shear_relax_tscale*60.0)
                 RVBLTEN(i,k,j) = RVBLTEN(i,k,j) - (vvv_avg(k)-0.0 )/(shear_relax_tscale*60.0)
           ENDDO; ENDDO
        ENDDO
      ENDDO         
   ENDIF

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
                 RVBLTEN(i,k,j) = RVBLTEN(i,k,j) - (vvv_avg(k)-0.0)/(relax_umean_tscale*60.0)
           ENDDO; ENDDO
        ENDDO
      ENDDO         
   ENDIF

! rce well tested using YSU and no pbl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!start WTG here (Wang and Sobel 2011, 2012) ------------------------------------------
!------------------------------------------
!------------------------------------------

   IF(rce_wtg_flag==101.and.(bl_pbl_physics==YSUSCHEME.or. bl_pbl_physics==1000))then
      DO k=kts,kde-1
          if(th_rce(k) < 0) then
               WRITE(wrf_err_message,*) 'the_rce(',k,') less than 0'
               CALL wrf_error_fatal(TRIM(wrf_err_message))
          endif
      ENDDO
      no_points = float((ide-ids)*(jde-jds))

      if(pblh_opt==0)then
        pblh_sum = 0.0        
        DO ij = 1 , num_tiles; 
           DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
            pblh_sum = pblh_sum + pblh(i,j)
           ENDDO; ENDDO 
        ENDDO
        pblh_domain = wrf_dm_sum_real ( pblh_sum )/no_points
      elseif(pblh_opt==1)then
        pblh_max = 0.0        
        DO ij = 1 , num_tiles; 
           DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
            if(pblh_max<pblh(i,j)) then
               pblh_max =  pblh(i,j)
               irnc=i; jrnc=j;
            endif
           ENDDO; ENDDO 
        ENDDO
        CALL wrf_dm_maxval ( pblh_max,  irnc, jrnc)
        pblh_domain = pblh_max
      elseif(pblh_opt>100)then
        pblh_domain = pblh_opt 
      endif

      !if( pblh_domain < 500.0 .or. itimestep*dt < 36000 )then
      !    pblh_domain = 1000.0 
      !endif
      
      !print*,'time',itimestep*dt , 'pblh',pblh_domain
         
      mu_sum = 0.0
      DO ij = 1 , num_tiles; 
         DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               mu_sum = mu_sum + mut(i,j)
         ENDDO; ENDDO
      ENDDO
      mu_sum = wrf_dm_sum_real ( mu_sum )
      mu_avg = mu_sum/no_points


      DO k=kts,kde-1
         z_sum = 0.0
         DO ij = 1 , num_tiles; 
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               z_2d(i,j) = z(i,k,j) 
               z_sum = z_sum + z_2d(i,j)
                 !if (rce_wtg_diag==1.and.i==1.and.j==1)then        
                 !    WRITE(outstring,*) 'i,j,k,z,pblh,=',i,j,k,z(i,k,j),pblh(i,j)
                 !    CALL wrf_message ( TRIM(outstring) )
                 !endif
            ENDDO; ENDDO
         ENDDO
         z_avg(k) = wrf_dm_sum_real ( z_sum )
         z_avg(k) = z_avg(k) / no_points
      ENDDO
  
      if (rce_wtg_diag==1)then        
          !WRITE(outstring,*) 'pblh=',pblh
          !WRITE(*,*) 'pblh=',pblh
          !CALL wrf_message ( TRIM(outstring) )
          WRITE(outstring,*) 'i,j,pblh_domain,=',irnc,jrnc,pblh_domain
          CALL wrf_message ( TRIM(outstring) )
      endif
      
      DO k=kts,kde-1
      !DO k=kts,min(kte+1,kde)
         omega_wtg(k)=0.0
         Dth_sum = 0.0
         rho_sum=0.0
         qv_sum = 0.0
         Th_2d = 0.0
         rho_2d = 0.0
         DO ij = 1 , num_tiles; 
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               th_2d(i,j)=th_phy(i,k,j)
               dth_sum = dth_sum + th_2d(i,j)
               
               rho_2d(i,j)=rho(i,k,j)
               rho_sum = rho_sum + rho_2d(i,j)
                
               qv_2d(i,j)=qv_curr(i,k,j)
               qv_sum = qv_sum + qv_2d(i,j)
         
            ENDDO; ENDDO 
         ENDDO
         !DO j = jps, min(jpe,jde-1)
         !    DO i = ips, min(ipe,ide-1)
         !      dth_sum = dth_sum + th_2d(i,j)
         !    ENDDO
         !ENDDO
         th_avg(k) = wrf_dm_sum_real ( dth_sum )
         th_avg(k) = th_avg(k) / no_points
         
         rho_avg(k) = wrf_dm_sum_real ( rho_sum )
         rho_avg(k) = rho_avg(k) / no_points
         
         if(rho_avg(k)>2.and.k<kde-1) then
            print*,'k,rho',k,rho_avg(k),th_avg(k)
            stop
         endif
         qv_avg(k) = wrf_dm_sum_real ( qv_sum )
         qv_avg(k) = qv_avg(k) / no_points

         !the_wtg_forcing(k) = 0.0 
         DO ij = 1 , num_tiles; 
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
              if( z_avg(k) .gt. pblh_domain )then ! above the pbl height
              !if( z(i,k,j) .gt. pblh_domain )then ! above the pbl height
              !if( z(i,k,j) .gt. pbl(i,j) )then ! above the pbl height
                 RTHBLTEN(i,k,j)= RTHBLTEN(i,k,j) - (th_avg(k)-th_rce(k))/(tao_time*3600.0)
                 cpm = cp * (1. + 0.8 * qv_curr(i,k,j))
                 evapor = p_phy(i,k,j)*r_v*qv_curr(i,k,j)/(r_d+r_v*qv_curr(i,k,j))
                 rho_d = (p_phy(i,k,j)-evapor)/t_phy(i,k,j)/r_d
                 coef_tmp = cpm*rho(i,k,j)/rho_d* pi_phy(i,k,j)*(-1.0)*(znw(k+1)-znw(k))
                 thetaint_tmp(i,j) = thetaint_tmp(i,j) - coef_tmp*(th_avg(k)-th_rce(k))/(tao_time*3600.0)
              endif
            ENDDO; ENDDO
         ENDDO
         if( z_avg(k) .gt. pblh_domain )then ! above the pbl height
            rthwtgten(k) =  -(th_avg(k)-th_rce(k))/(tao_time*3600.0)
         else
            rthwtgten(k) = 0.0
         endif


         if( z_avg(k) .gt. pblh_domain )then ! above the pbl height
             the_wtg_forcing(k) = the_wtg_forcing(k) - (th_avg(k)-th_rce(k))/(tao_time*3600.0)*DT
         endif

         !WRITE(outstring,*)'k=',k,'th_avg',th_avg(k)
         !CALL wrf_message ( TRIM(outstring) )
         !WRITE(outstring,*) grid_str,'Domain average of dpsdt, dmudt (mb/3h): ', xtime, th_sum/no_points, &
      ENDDO

      DO k=kde-2,kts,-1
      !DO k=kde-1,kts,-1
        omega_wtg(k) = 0.0
        if(z_avg(k)>pblh_domain)then
           k_abovepbl_first = k
           if(k==kte)then
              dthdeta(k) = (th_avg(k+1) - th_avg(k))/(znu(k+1) - znu(k))
              dqvdeta(k) = (qv_avg(k+1) - qv_avg(k))/(znu(k+1) - znu(k))
           elseif(k==kde-1)then
              dthdeta(k) = (th_avg(k) - th_avg(k-1))/(znu(k) - znu(k-1))
              dqvdeta(k) = (qv_avg(k) - qv_avg(k-1))/(znu(k) - znu(k-1))
           else
              dthdeta(k) = 0.5*( (th_avg(k+1) - th_avg(k))/(znu(k+1) - znu(k)) + &
                              (th_avg(k) - th_avg(k-1))/(znu(k) - znu(k-1)) )
              dqvdeta(k) = 0.5*( (qv_avg(k+1) - qv_avg(k))/(znu(k+1) - znu(k)) + &
                              (qv_avg(k) - qv_avg(k-1))/(znu(k) - znu(k-1)) )
           endif
           !the_wtg_forcing(k) = -(th_avg(k)-th_rce(k))/(tao_time*3600.0)
           th_temp = -dthdeta(k)*9.81*rho_avg(k)/mu_avg*1e3  ! ensure dtheta/dz > 1 K/km
           if(th_temp<1) then
               if (rce_wtg_diag==11 )then       
               !if (rce_wtg_diag==1 .and. mod((itimestep-1)*dt/86400,1.0)==0)then       
                  WRITE(outstring,'(a30, i3, F9.3, F9.3 )') &
                  'Warning: dthdeta<1 k/km',k,dthdeta(k), th_temp
                  CALL wrf_message ( TRIM(outstring) )
                  dthdeta(k) = -mu_avg/(9.81*rho_avg(k)*1e3) 
               endif
           endif
              
           omega_wtg(k) = (th_avg(k)-th_rce(k))/(tao_time*3600.0) / dthdeta(k)
           

        else ! below boundary: linearly interpolation with z_avg
           omega_wtg(k) = (omega_wtg(k_abovepbl_first)-0.0 )/(z_avg(k_abovepbl_first)) * (z_avg(k))
           ! omega_wtg(k) = (omega_wtg(k_abovepbl_first)-0.0 )/(z_avg(k_abovepbl_first)-z_avg(1)) * (z_avg(k)-z_avg(1)) 
           dthdeta(k)=0.0
        endif
      ENDDO 
      
      DO k=kts,kde-1
         qv_sum = 0.0
         qtend(k) = 0.0 
                 
         if(qv_avg_flag==1.and.k>1.and.k<kde-1)then
               if(wtg_adv_opt==1)then     ! flux form: wrong 
                    qtend2 = (omega_wtg(k+1)*qv_avg(k+1) - omega_wtg(k)*qv_avg(k))/(znu(k+1)-znu(k))
                    qtend1 = (omega_wtg(k)*qv_avg(k) - omega_wtg(k-1)*qv_avg(k-1))/(znu(k)-znu(k-1))
                    qtend(k) = 0.5*(qtend1+qtend2)
               elseif(wtg_adv_opt==2)then     ! advective form
                 if(upwind_opt==1)then ! upwind wind advection
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend2 = omega_p*(qv_avg(k+1) - qv_avg(k))/(znu(k+1)-znu(k)) 
                    qtend1 = omega_n*(qv_avg(k) - qv_avg(k-1)) /(znu(k)-znu(k-1))
                    qtend(k) = (qtend1+qtend2)
                 elseif(upwind_opt==0)then
                    qtend2 = (omega_wtg(k+1)+omega_wtg(k))*0.5*(qv_avg(k+1) - qv_avg(k))/(znu(k+1)-znu(k)) 
                    qtend1 = (omega_wtg(k)+omega_wtg(k-1))*0.5*(qv_avg(k) - qv_avg(k-1)) /(znu(k)-znu(k-1))
                    qtend(k) = 0.5*(qtend1+qtend2)
                 endif
               endif
         elseif(qv_avg_flag==1.and.k==kde-1)then
               if(wtg_adv_opt==1)then     ! flux form: wrong
                    qtend1 = (omega_wtg(k)*qv_avg(k) - omega_wtg(k-1)*qv_avg(k-1))/(znu(k)-znu(k-1))
                    qtend(k) = qtend1 
                    !WRITE(outstring,*) 'k=',k,'omega_wtg',omega_wtg(k-1),omega_wtg(k), &
                    ! 'qv',qv_avg(k-1),qv_avg(k),'znu',znu(k-1),znu(k),'qtend',qtend(k)
                    !CALL wrf_message ( TRIM(outstring) )
               elseif(wtg_adv_opt==2)then     ! advective form:
                  if(upwind_opt==1)then ! upwind wind advection
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend1 = omega_n*(qv_avg(k) - qv_avg(k-1)) /(znu(k)-znu(k-1))
                    qtend(k) = qtend1
                  else
                    qtend1 = (omega_wtg(k)+omega_wtg(k-1))*0.5*(qv_avg(k) - qv_avg(k-1)) /(znu(k)-znu(k-1))
                    qtend(k) = qtend1 
                  endif
                endif
         elseif(qv_avg_flag==1.and.k==1)then
               if(wtg_adv_opt==1)then     ! flux form: wrong
                    qtend2 = (omega_wtg(k+1)*qv_avg(k+1) - omega_wtg(k)*qv_avg(k))/(znu(k+1)-znu(k))
                    qtend(k) = qtend2
               elseif(wtg_adv_opt==2)then     ! advective form:
                 if(upwind_opt==1)then ! upwind wind advection
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    qtend2 = omega_p*(qv_avg(k+1) - qv_avg(k))/(znu(k+1)-znu(k)) 
                    qtend(k) = qtend2
                 else
                    qtend2 = (omega_wtg(k+1)+omega_wtg(k))*0.5*(qv_avg(k+1) - qv_avg(k))/(znu(k+1)-znu(k)) 
                    qtend(k) = qtend2
                 endif
               endif
         endif
             
         !if( abs(qtend(k)) < 1e-10) qtend(k) = 0.0 
                 
         if(qv_hadv_flag==1)then ! horizontl advection of moisture
               domgdeta=0
               qtend_qvhadv(k) = 0.0
               if(k>1.and.k<kde-1)then
                  domgdeta = ( (omega_wtg(k)-omega_wtg(k-1))/(znu(k)-znu(k-1))+ &
                                  (omega_wtg(k+1)-omega_wtg(k))/(znu(k+1)-znu(k))  )*0.5
               elseif(k==1)then
                  domgdeta =      (omega_wtg(k+1)-omega_wtg(k))/(znu(k+1)-znu(k))  
               elseif(k==kde-1)then
                  domgdeta =      (omega_wtg(k)-omega_wtg(k-1))/(znu(k)-znu(k-1))  
               endif
               if(domgdeta>0)then ! region of inward flow
                      qtend_qvhadv(k) = (qv_rce(k)*1e-3-qv_avg(k))*domgdeta
               endif
               qv_wtg_hadv_forcing(k) = qv_wtg_hadv_forcing(k) -qtend_qvhadv(k)*DT  ! accumulated qv tendency
         endif

         ! dry advection
         qtend_qadv(k)=0.0
         if(qadv_flag==1)then  ! relax to qrce
               if( znu(k) < qadv_plev_dn .and. znu(k) > qadv_plev_up )then
                   qtend_qadv(k) = (qv_rce(k)*1e-3-qv_avg(k))/3600.0/qadv_tau
               endif
               qv_wtg_hadv_forcing(k) = qv_wtg_hadv_forcing(k) -qtend_qadv(k)*DT  ! accumulated qv tendency
         elseif(qadv_flag==2)then ! relax to  dry sounding
               if( znu(k) < qadv_plev_dn .and. znu(k) > qadv_plev_up .and. qadv_tau>0)then
                   qtend_qadv(k) = (-qv_avg(k))/3600.0/qadv_tau
               endif
               qv_wtg_hadv_forcing(k) = qv_wtg_hadv_forcing(k) -qtend_qadv(k)*DT  ! accumulated qv tendency
         elseif(qadv_flag==3)then ! relax to  dry sounding at freqeuncty
               if( znu(k) < qadv_plev_dn .and. znu(k) > qadv_plev_up )then
                   qtend_qadv(k) = (-qv_avg(k))/3600.0/qadv_tau*(1.0+0.5*cos(2*3.1415926/3600.0/qforcing_tau*curr_secs))
                   !qtend_qadv(k) = (-qv_avg(k))/3600.0/qadv_tau*(1.0+0.5*cos(2*3.1415926/86400.0/2*curr_secs))
               endif
               qv_wtg_hadv_forcing(k) = qv_wtg_hadv_forcing(k) -qtend_qadv(k)*DT  ! accumulated qv tendency
         endif                
 

         qv_sum = 0.0
         DO ij = 1 , num_tiles; 
              DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
                 if(qv_avg_flag==1)    then
                     !RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j) - qtend(k)
                 elseif(qv_avg_flag==-1)then
                     if(k>1.and.k<kde-1)then
                        omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                        omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                        qtend2 = omega_p*(qv_curr(i,k+1,j) - qv_curr(i,k,j))/(znu(k+1)-znu(k)) 
                        qtend1 = omega_n*(qv_curr(i,k,j) - qv_curr(i,k-1,j)) /(znu(k)-znu(k-1))
                        qtend(k) = (qtend1+qtend2)
                     elseif(k==kde-1)then
                        omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                        qtend1 = omega_n*(qv_curr(i,k,j) - qv_curr(i,k-1,j)) /(znu(k)-znu(k-1))
                        qtend(k) = qtend1
                     elseif(k==1)then
                        omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                        qtend2 = omega_p*(qv_curr(i,k+1,j) - qv_curr(i,k,j))/(znu(k+1)-znu(k)) 
                        qtend(k) = qtend2
                     endif
                     qv_2d(i,j)=qtend(k)
                     qv_sum = qv_sum + qv_2d(i,j)
                 endif
                   qvint_tmp(i,j) = qvint_tmp(i,j) + qtend(k)*(znw(k+1)-znw(k))
                   RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j) - qtend(k)
                 if(qv_hadv_flag==1)then
                     RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j) + max(qtend_qvhadv(k), (-qv_curr(i,k,j))/DT)
                 endif
                 if(qadv_flag>0)then
                     qvhadvint_tmp(i,j) = qvhadvint_tmp(i,j) + qtend_qadv(k)*(znw(k+1)-znw(k))
                     RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j) + max(qtend_qadv(k), (-qv_curr(i,k,j))/DT)
                 endif
                 if(qadv_flag==10 .and. znu(k) < qadv_plev_dn .and. znu(k) > qadv_plev_up)then
                 ! impose drying at each point
                     if(qadv_tau*3600.0<DT.and.qadv_tau>0)then; 
                            WRITE(wrf_err_message,*) 'qtau = ', qadv_tau, ' is less than DT', DT
                            CALL wrf_error_fatal(TRIM(wrf_err_message))
                     endif
                     if(qadv_tau == -1 ) then
                           ! the following case remove water immediatelly
                           RQVBLTEN(i,k,j)  = RQVBLTEN(i,k,j) + (-qv_curr(i,k,j))/DT
                     else
                           RQVBLTEN(i,k,j)  = RQVBLTEN(i,k,j) + (-qv_curr(i,k,j))/3600.0/qadv_tau
                     endif
                 endif
              ENDDO; ENDDO
         ENDDO


        if(qv_avg_flag==-1)then
           qtend_avg(k) = wrf_dm_sum_real ( qv_sum )
           qtend_avg(k) = qtend_avg(k) / no_points
           qv_wtg_forcing(k) = qv_wtg_forcing(k) -qtend_avg(k)*DT  ! accumulated qv tendency
           rqvwtgten(k) = - qtend_avg(k)                           ! wtg qv tendency 
        elseif(qv_avg_flag==1)then
           qv_wtg_forcing(k) = qv_wtg_forcing(k) -qtend(k)*DT  ! accumulated qv tendency
           rqvwtgten(k) = - qtend(k)                           ! wtg qv tendency 
        endif

         !if (rce_wtg_diag==1 )then       
         if (rce_wtg_diag==1 .and. mod((itimestep-1)*dt/86400,1.0)==0)then       
           if(k==kts)then 
             WRITE(outstring,'(a10,   a8,   a8,   a9,   a13,   a15,   a10,   a10,   a10 )') &
               'ww_wtg','th_avg', 'th_rce','qv_avg',&
               'dthdeta','qforcing','z_avg','zh','z'
             CALL wrf_message ( TRIM(outstring) )
           endif
           WRITE(outstring,'(F10.4, F8.2, F8.2, F9.4, g13.5, g15.5, f10.5, f10.5, f10.5 )') &
             !omega_wtg(k), th_avg(k),th_rce(k),qv_avg(k)*1e3 &
             omega_wtg(k)/rho_avg(k)*98000./(-9.81)*100.0, th_avg(k),th_rce(k),qv_avg(k)*1e3 &
             ,dthdeta(k), qtend(k)*1e6, z_avg(k)/1e3, pblh(1,1)/1e3, z(1,k,1)/1e3

           CALL wrf_message ( TRIM(outstring) )
           if(k==kde-1)then 
             WRITE(outstring,'(a10,   a8,   a8,   a9,   a9,   a15,   a10,   a7,   a10 )') &
               'ww_wtg','th_avg', 'th_rce','qv_avg',&
               'dthdeta','qforcing','z_avg','zh','z'
             CALL wrf_message ( TRIM(outstring) )
           endif
         endif 
      ENDDO 
         
      DO ij = 1 , num_tiles; 
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
            ACTWTG(i,j) = ACTWTG(i,j)+thetaint_tmp(i,j)*mut(i,j)/g*DT  ! unit: mm, assume rhowater=1kg/m3
            ACQVWTG(i,j) = ACQVWTG(i,j)+qvint_tmp(i,j)*mut(i,j)/9.81*DT  ! unit: mm, assume rhowater=1kg/m3
            ACQVHADVWTG(i,j) = ACQVHADVWTG(i,j)+qvhadvint_tmp(i,j)*mut(i,j)/9.81*DT  ! unit: mm
        ENDDO
        ENDDO
     ENDDO
     ! convert to w= cm/s
         
      DO k=kts,kde-1
         if( k<kde-1 .and. k>1)then 
            omega_wtg(k) = -omega_wtg(k)/rho_avg(k)/9.81*mu_avg*100.0
         endif
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
! --  Start wave coupling here: 
   ELSEIF(rce_wtg_flag==102.and.(bl_pbl_physics==YSUSCHEME.or. bl_pbl_physics==1000))then !! wave coupling
      no_points = float((ide-ids)*(jde-jds))
      !eps_wavedamping(k) = 0.0 
      DO k=kts,kde-1
         z_sum = 0.0
         z_2d = 0.0
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
         omega_wtg(k)=0.0
         Dth_sum = 0.0
         rho_sum=0.0
         qv_sum = 0.0
         phalf_sum = 0.0
         pfull_sum = 0.0

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
         
            ENDDO; ENDDO 
         ENDDO
         !DO j = jps, min(jpe,jde-1)
         !    DO i = ips, min(ipe,ide-1)
         !      dth_sum = dth_sum + th_2d(i,j)
         !    ENDDO
         !ENDDO
         th_avg(k) = wrf_dm_sum_real ( dth_sum )
         th_avg(k) = th_avg(k) / no_points
         
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
       
        
     ENDDO
     do k=1,kde-1
         !if(phalf_avg(k)<10000.0)then 
         !if( (z_avg(k)>wave_ztop) .and. ( wave_ztop >0 )  )then 
         !    np=k-1
         !    exit 
         !elseif( (phalf_avg(k) < -wave_ztop ) .and. ( wave_ztop < 0 ) )then
         if( phalf_avg(k) < wave_ztop ) then
             np = k-1
             exit
         endif
     enddo


     ! ----------------------------------------------------------------------------
     !     compute grid spacing between midpoint levels
      !np=kde-1
      !print*,'np=',np,'p at top=', phalf_avg(np), pfull_avg(np), znu(np), znw(np)
      do k=2,np
         dp(k)=phalf_avg(k)-phalf_avg(k-1)
      enddo
      dp(1)=2*(phalf_avg(1)-pfull_avg(1))
      dp(np+1)=2*(pfull_avg(np+1)-phalf_avg(np))
     !print*,'dp=',dp
      !dp(1)=(phalf_avg(1)-pfull_avg(1))
      !dp(np+1)=(phalf_avg(np+1)-phalf_avg(np))
      !dp(np+1)=(pfull_avg(np+1)-phalf_avg(np))
     
      !wavedamping = 1/86400.0 ! damping time scale: 1 day as in Kuang 2008 and 2010 
      !wavenumber = 2*pi/10000.0e3 ! wavenumber

      ! Interpolate theta_bg to coordinate levels:
      DO k = kts, kte
        k2 = kte
        DO WHILE( z_avg(k2) .gt. z_bg(k)  .and.  k2 .gt. 1 )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.kte)then
          theta_bg_int(k) = theta_bg(k2) + ( theta_bg(k2) - theta_bg(k2-1) )   &
                              * (     z_avg(k) - z_bg(k2)   )   &
                              / ( z_bg(k2) - z_bg(k2-1) )
        else
          theta_bg_int(k) = theta_bg(k2) + ( theta_bg(k2+1) - theta_bg(k2) )   &
                              * (       z_avg(k) - z_bg(k2) )   &
                              / ( z_bg(k2+1) - z_bg(k2) )
        endif
       !print*,'z_bg,theta_bg=', z_bg(k), theta_bg(k),'z_avg, theta_bg_int=',theta_bg_int(k), z_avg(k)
      ENDDO
      DO k = kts, kte
        k2 = kte
        DO WHILE( z_avg(k2) .gt. z_bg(k)  .and.  k2 .gt. 1 )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.kte)then
          qv_bg_int(k) = qv_bg(k2) + ( qv_bg(k2) - qv_bg(k2-1) )   &
                              * (     z_avg(k) - z_bg(k2)   )   &
                              / ( z_bg(k2) - z_bg(k2-1) )
        else
          qv_bg_int(k) = qv_bg(k2) + ( qv_bg(k2+1) - qv_bg(k2) )   &
                              * (       z_avg(k) - z_bg(k2) )   &
                              / ( z_bg(k2+1) - z_bg(k2) )
        endif
        qv_bg_int(k) = qv_bg_int(k)*1e-3;
        !qv_bg_int(k) = qv_avg(k);
       !print*,'z_bg,theta_bg=', z_bg(k), theta_bg(k),'z_avg, theta_bg_int=',theta_bg_int(k), z_avg(k)
      ENDDO
 


      !     Gauss Elimination
      rhs=0.
      do k=1,np
         theta_v_prime = th_v_avg(k)-theta_bg_int(k)*(1+0.608*qv_bg_int(k))
         rhs(k)=wavenumber*wavenumber/wavedamping*R_d*theta_v_prime*(phalf_avg(k)/1.e5)**(R_d/Cp)/phalf_avg(k)*dp(k)*dp(k+1)*0.5
         !rhs(k)=wavenumber*wavenumber/wavedamping*(1e5/phalf_avg(k))*R_d*(th_v_avg(k)-theta_bg_int(k))*(phalf_avg(k)/1.e5)**(R_d/Cp)/phalf_avg(k)*dp(k)*dp(k+1)*0.5
         !rhs(k)=wavenumber*wavenumber*R_d*(tv_avg(k)-tvirtual_bg(k))/pfull_avg(k)*dp(k)*dp(k+1)*0.5 ! see equation (A.8) in Blossy 2009
      enddo
     
 
      !     set up the tridiagonal matrix
      do k=1,np
         aa(k)=dp(k+1)/(dp(k)+dp(k+1))
         bb(k)=-1.;
         cc(k)=dp(k)/(dp(k)+dp(k+1))
      enddo
      
      !     symmetric lower BC
      aa(1)=0.;
      bb(1)=-(2*dp(2)+dp(1))/(dp(1)+dp(2))
      
      !     symmetric upper BC
      bb(np) = -(2*dp(np)+dp(np+1))/(dp(np)+dp(np+1))
      cc(np) = 0.

      !     Gaussian Elimination without pivoting
      do k=1,np-1
         bb(k+1)=bb(k+1)-aa(k+1)/bb(k)*cc(k)
         rhs(k+1)=rhs(k+1)-aa(k+1)/bb(k)*rhs(k)
      enddo
      !    Backward substitution
      rhs(np)=rhs(np)/bb(np)
      do k= np-1,1,-1
         rhs(k)=(rhs(k)-cc(k)*rhs(k+1))/bb(k)
      end do
      
      do k=1,np
         omega_wtg(k)=rhs(k)  ! omega, this is pressure velocity!!
      enddo
     ! ----------------------------------------------------------------------------
      !do k=kts, kde-1 !SW==: using model water vapor
      !   qv_bg_int(k) = qv_avg(k)
      !   theta_bg_int(k) = th_avg(k)
      !enddo

      DO k=kts,kde-1
         qv_sum = 0.0
         qtend(k) = 0.0 
         ttend(k) = 0.0 
         if(k>1.and.k<kde-1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend2 = omega_p*(qv_avg(k+1) - qv_avg(k))/(phalf_avg(k+1)-phalf_avg(k)) 
                    qtend1 = omega_n*(qv_avg(k) - qv_avg(k-1)) /(phalf_avg(k)-phalf_avg(k-1))
                    qtend(k) = (qtend1+qtend2)
         elseif(k==kde-1)then
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend1 = omega_n*(qv_avg(k) - qv_avg(k-1)) /(phalf_avg(k)-phalf_avg(k-1))
                    qtend(k) = qtend1
         elseif(k==1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    qtend2 = omega_p*(qv_avg(k+1) - qv_avg(k))/(phalf_avg(k+1)-phalf_avg(k)) 
                    qtend(k) = qtend2
         endif

         if(k>1.and.k<kde-1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend2 = omega_p*(th_avg(k+1) - th_avg(k))/(phalf_avg(k+1)-phalf_avg(k)) 
                    qtend1 = omega_n*(th_avg(k) - th_avg(k-1)) /(phalf_avg(k)-phalf_avg(k-1))
                    ttend(k) = (qtend1+qtend2)
         elseif(k==kde-1)then
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend1 = omega_n*(th_avg(k) - th_avg(k-1)) /(phalf_avg(k)-phalf_avg(k-1))
                    ttend(k) = qtend1
         elseif(k==1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    qtend2 = omega_p*(th_avg(k+1) - th_avg(k))/(phalf_avg(k+1)-phalf_avg(k)) 
                    ttend(k) = qtend2
         endif
             

         DO ij = 1 , num_tiles; 
              DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
                 RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j) - qtend(k)
                 RTHBLTEN(i,k,j) = RTHBLTEN(i,k,j) - ttend(k)

                 qvint_tmp(i,j) = qvint_tmp(i,j) + qtend(k)*(znw(k+1)-znw(k))
                 cpm = cp * (1. + 0.8 * qv_curr(i,k,j))
                 evapor = p_phy(i,k,j)*r_v*qv_curr(i,k,j)/(r_d+r_v*qv_curr(i,k,j))
                 rho_d = (p_phy(i,k,j)-evapor)/t_phy(i,k,j)/r_d
                 coef_tmp = cpm*rho(i,k,j)/rho_d* pi_phy(i,k,j)*(-1.0)*(znw(k+1)-znw(k))
                 thetaint_tmp(i,j) = thetaint_tmp(i,j) - coef_tmp*(th_avg(k)-th_rce(k))/(tao_time*3600.0)
            ENDDO; ENDDO
         ENDDO

         !the_wtg_forcing(k) = 0.0
         qv_wtg_forcing(k) = qv_wtg_forcing(k) -qtend(k)*DT  ! accumulated qv tendency
         rqvwtgten(k) = - qtend(k)                           ! wtg qv tendency 
         rthwtgten(k) = -ttend(k) 


         if (rce_wtg_diag==1.and.i==1.and.j==1)then        
         endif
         !if (rce_wtg_diag==1 )then        
         if (rce_wtg_diag==1 )then       
           if(k==kts)then 
             WRITE(outstring,'(a10,   a12,   a12,   a9,   a13,   a15,   a10,   a10,   a10,   a10,  a10,  a10,  a10 )') &
               'ww_wtg','th_avg', 'th_v_prime','qv_avg', 'qv_prime',&
               '-ttend','-qtend','rhs','z_avg','z',  'pres', 'dp'
             CALL wrf_message ( TRIM(outstring) )
           endif
           WRITE(outstring,'(G10.4, F12.2, F12.2, F12.4, g13.5, g13.5, g15.5, f12.5, f10.5, f10.5, f15.5, f13.5, f13.5, f13.5 )') &
             !omega_wtg(k), th_avg(k),th_rce(k),qv_avg(k)*1e3 &
             omega_wtg(k)/rho_avg(k)/(-9.81)*100.0, th_v_avg(k), th_v_avg(k)-theta_bg_int(k)*(1+0.608*qv_bg_int(k)), qv_avg(k)*1e3, (qv_avg(k)-qv_bg_int(k))*1e3 &
             !omega_wtg(k)/rho_avg(k)/(-9.81), th_avg(k),theta_bg_int(k),qv_avg(k)*1e3 &
             ,-ttend(k)*86400, -qtend(k)*1e6, rhs(k), z_avg(k)/1e3,  z(1,k,1)/1e3, phalf_avg(k), dp(k)

           CALL wrf_message ( TRIM(outstring) )
           if(k==kde-1)then 
             WRITE(outstring,'(a10,   a12,   a12,   a9,   a13,   a15,   a10,   a10,   a10,   a10 )') &
               'ww_wtg','th_avg', 'theta_bg','qv_avg',&
               'dthdeta','qforcing','rhs','zh','z',  'pres'
             CALL wrf_message ( TRIM(outstring) )
           endif
         endif 
      ENDDO ! end of the big do loop 

     !diagnose

!     print*,'k, w_wtg, -qtend, -ttend'
!     do k=1,kde-1
!         print*,k, omega_wtg(k)/rho_avg(k)/(-9.81)*100.0, -qtend(k), -ttend(k)*86400
!     enddo
      DO ij = 1 , num_tiles; 
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
            ACTWTG(i,j) = ACTWTG(i,j)+thetaint_tmp(i,j)*mut(i,j)/g*DT  ! unit: mm, assume rhowater=1kg/m3
            ACQVWTG(i,j) = ACQVWTG(i,j)+qvint_tmp(i,j)*mut(i,j)/9.81*DT  ! unit: mm, assume rhowater=1kg/m3
            ACQVHADVWTG(i,j) = ACQVHADVWTG(i,j)+qvhadvint_tmp(i,j)*mut(i,j)/9.81*DT  ! unit: mm
        ENDDO
        ENDDO
     ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
!------------------------------------------
! --  Start vertical mode based spectral WTG: 
   ELSEIF(rce_wtg_flag==103.and.(bl_pbl_physics==YSUSCHEME.or. bl_pbl_physics==1000))then 
   !ELSEIF(rce_wtg_flag==3.and.(bl_pbl_physics==YSUSCHEME.or. bl_pbl_physics==0))then 
      no_points = float((ide-ids)*(jde-jds))
      !eps_wavedamping(k) = 0.0 
      DO k=kts,kde-1
         z_sum = 0.0
         z_2d = 0.0
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
         omega_wtg(k)=0.0
         Dth_sum = 0.0
         rho_sum=0.0
         qv_sum = 0.0
         phalf_sum = 0.0
         pfull_sum = 0.0

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
         
            ENDDO; ENDDO 
         ENDDO
         !DO j = jps, min(jpe,jde-1)
         !    DO i = ips, min(ipe,ide-1)
         !      dth_sum = dth_sum + th_2d(i,j)
         !    ENDDO
         !ENDDO
         th_avg(k) = wrf_dm_sum_real ( dth_sum )
         th_avg(k) = th_avg(k) / no_points
         
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
       
        
     ENDDO
    
      ! Interpolate theta_bg to coordinate levels:
      DO k = kts, kte
        k2 = kte
        DO WHILE( z_avg(k2) .gt. z_bg(k)  .and.  k2 .gt. 1 )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.kte)then
          theta_bg_int(k) = theta_bg(k2) + ( theta_bg(k2) - theta_bg(k2-1) )   &
                              * (     z_avg(k) - z_bg(k2)   )   &
                              / ( z_bg(k2) - z_bg(k2-1) )
        else
          theta_bg_int(k) = theta_bg(k2) + ( theta_bg(k2+1) - theta_bg(k2) )   &
                              * (       z_avg(k) - z_bg(k2) )   &
                              / ( z_bg(k2+1) - z_bg(k2) )
        endif
       !print*,'z_bg,theta_bg=', z_bg(k), theta_bg(k),'z_avg,
       !theta_bg_int=',theta_bg_int(k), z_avg(k)
      ENDDO
      DO k = kts, kte
        k2 = kte
        DO WHILE( z_avg(k2) .gt. z_bg(k)  .and.  k2 .gt. 1 )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.kte)then
          qv_bg_int(k) = qv_bg(k2) + ( qv_bg(k2) - qv_bg(k2-1) )   &
                              * (     z_avg(k) - z_bg(k2)   )   &
                              / ( z_bg(k2) - z_bg(k2-1) )
        else
          qv_bg_int(k) = qv_bg(k2) + ( qv_bg(k2+1) - qv_bg(k2) )   &
                              * (       z_avg(k) - z_bg(k2) )   &
                              / ( z_bg(k2+1) - z_bg(k2) )
        endif
        qv_bg_int(k) = qv_bg_int(k)*1e-3;
      ENDDO
 

      do k = kts+1,kte
         dthdz(k) = ((th_avg(k)-th_avg(k-1))/(z_avg(k)-z_avg(k-1))+(th_avg(k+1)-th_avg(k))/(z_avg(k+1)-z_avg(k))) *0.5
         dthdz(k) = max(dthdz(k), 0.001)
         dqvdz(k) = ((qv_avg(k)-qv_avg(k-1))/(z_avg(k)-z_avg(k-1)) +(qv_avg(k+1)-qv_avg(k))/(z_avg(k+1)-z_avg(k)))*0.5
      enddo
      dthdz(1) = ( (th_avg(2)-th_avg(1))/(z_avg(2)-z_avg(1))); dthdz(1) =max(dthdz(1), 0.001)
      dqvdz(1) = ( (qv_avg(2)-qv_avg(1))/(z_avg(2)-z_avg(1)))

      th_scm_target_modellevels = 0.0
      w_subs_modellevels = 0.0
      do k = kts,kte
          if(wtg_opt == 200)then
              th_scm_target_modellevels(k) = th_avg(k) - interp_0(theta_bg,z_bg, z_avg(k),nl_bg)
          elseif(wtg_opt == 201)then
              th_scm_target_modellevels(k) = th_avg(k) - theta_bg_int(k) 
          elseif(wtg_opt == 210)then
              if(k>=kts+1)then
                  dthdz(k) = ((th_avg(k)*(1+0.608*qv_avg(k))-th_avg(k-1)*(1+0.608*qv_avg(k-1)))/(z_avg(k)-z_avg(k-1))+ &
                             (th_avg(k+1)*(1+0.608*qv_avg(k+1))-th_avg(k)*(1+0.608*qv_avg(k)))/(z_avg(k+1)-z_avg(k)))*0.5
                  dthdz(k) = max(dthdz(k), 0.001)
              elseif(k==1)then
                  dthdz(1) = ((th_avg(2)*(1+0.608*qv_avg(2))-th_avg(1)*(1+0.608*qv_avg(1)))/(z_avg(2)-z_avg(1))); dthdz(1) = max(dthdz(1), 0.001)
              endif
              th_scm_target_modellevels(k) = th_avg(k)*(1+0.608*qv_avg(k))- &
                        interp_0(theta_bg,z_bg, z_avg(k), nl_bg)*(1+0.608*interp_0(qv_bg*1e-3, z_bg, z_avg(k), nl_bg))
          endif
         !th_scm_prime(k)  = th_scm_target_modellevels(k)
         if(11==1.and.mod(itimestep,100)==1)then
             if(k==kts) print*,'th, qv,th_bg, qv_bg='
             print*,th_avg(k), qv_avg(k), interp_0(theta_bg,z_bg,z_avg(k), nl_bg), interp_0(qv_bg*1e-3, z_bg, z_avg(k), nl_bg)
         endif
         th_over_dthdz(k) = th_scm_target_modellevels(k)/dthdz(k)
         do m=1,wtg_nmodes
            mode_tmp(k,m) = scm_m_vector(m*kme+(k))
            mode_tmp_deWeight(k,m) = mode_tmp(k,m)/scm_m_weight(k)
         enddo
     enddo

     !print*,'scm_zm(k)', scm_zm(1:40)
     !print*,'theta_avg', th_avg

     do k = kts,kte
         if( k>10 .and. scm_zm(k) == 0 )then
               nnlv = k-1
                if(mod(itimestep,100)==1) print*,'nnlv=',nnlv
                !if(i==1.and.j==1.and. mod(itimestep,100)==1) print*,'nnlv=',nnlv
               exit
         endif
      !nnlv = 35
     enddo

      modes_modellevels = 0.0
      do k = kts,kte
         if(z_avg(k)<=wave_ztop ) then
         do m=1,wtg_nmodes
            modes_modellevels(m,k) =interp_0(mode_tmp_deWeight(1:nnlv,m),scm_zm(1:nnlv), z_avg(k), nnlv)
            !modes_modellevels(m,k) = interp_0(scm_m_vector(m,:),scm_zm, z_avg(k), kme)
         enddo
         endif
         th_over_dthdz_int(k) = interp_0(th_over_dthdz(1:kte-1),z_avg(1:kte-1), scm_zm(k), kte-1)
         rho_int(k) = interp_0(rho_avg(1:kte-1), z_avg(1:kte-1), scm_zm(k),kte-1)
         if(mod(itimestep,100)==1)print 10111, k,scm_zm(k), mode_tmp(k,1), mode_tmp(k,2), mode_tmp_deWeight(k,1),mode_tmp_deWeight(k,2), z_avg(k), modes_modellevels(1,k), modes_modellevels(2,k)
         !if(i==1.and.j==1 .and. mod(itimestep,100)==1) print 10111, k,scm_zm(k), mode_tmp(k,1), mode_tmp(k,2), mode_tmp_deWeight(k,1),mode_tmp_deWeight(k,2), z_avg(k), modes_modellevels(1,k), modes_modellevels(2,k)
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

      do k=kts, kte
        do m= 1,wtg_nmodes !first 10 modes
         if( z_avg(k) < wave_ztop)then
             !w_subs_modellevels(k) = w_subs_modellevels(k) + am(m)/am_s(m)*modes_modellevels(m,k) /tao_time*scm_m_value(m)
             w_subs_modellevels(k) = w_subs_modellevels(k) + am(m)/am_s(m)*modes_modellevels(m,k) /tao_time*scm_m_value(m)**wtg_cexp
         endif
        enddo
      enddo

      do k=kts,kte
         omega_wtg(k) = -w_subs_modellevels(k)*9.81 ! use pressure velocity
      enddo
      do k=kts,kte
       if(mod(itimestep,100)==1)print 1001, k, omega_wtg(k), z_avg(k), -w_subs_modellevels(k)*dqvdz(k), &
                     -th_scm_target_modellevels(k)/tao_time,dthdz(k),  th_scm_target_modellevels(k)
       1001 format(i10,2x, f10.4, f10.3, e18.10, 2x, e18.10, 2x, e18.10, 2x,e18.10, 2x, f18.10, 2x, f18.10 )
      enddo

      
      !print*,'w_subs_modellevels',w_subs_modellevels

     ! ----------------------------------------------------------------------------
      do k=kts, kde-1 !SW==: using model water vapor
         qv_bg_int(k) = qv_avg(k)
         theta_bg_int(k) = th_avg(k)
      enddo


      DO k=kts,kde-1
       if(  scm_vert_adv_opt == 1) then
         qv_sum = 0.0
         qtend(k) = 0.0 
         ttend(k) = 0.0 
         if(k>1.and.k<kde-1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend2 = omega_p*(qv_bg_int(k+1) - qv_bg_int(k))/(phalf_avg(k+1)-phalf_avg(k)) 
                    qtend1 = omega_n*(qv_bg_int(k) - qv_bg_int(k-1)) /(phalf_avg(k)-phalf_avg(k-1))
                    qtend(k) = (qtend1+qtend2)
         elseif(k==kde-1)then
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend1 = omega_n*(qv_bg_int(k) - qv_bg_int(k-1)) /(phalf_avg(k)-phalf_avg(k-1))
                    qtend(k) = qtend1
         elseif(k==1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    qtend2 = omega_p*(qv_bg_int(k+1) - qv_bg_int(k))/(phalf_avg(k+1)-phalf_avg(k)) 
                    qtend(k) = qtend2
         endif

         if(k>1.and.k<kde-1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend2 = omega_p*(theta_bg_int(k+1) - theta_bg_int(k))/(phalf_avg(k+1)-phalf_avg(k)) 
                    qtend1 = omega_n*(theta_bg_int(k) - theta_bg_int(k-1)) /(phalf_avg(k)-phalf_avg(k-1))
                    ttend(k) = (qtend1+qtend2)
         elseif(k==kde-1)then
                    omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0) 
                    qtend1 = omega_n*(theta_bg_int(k) - theta_bg_int(k-1)) /(phalf_avg(k)-phalf_avg(k-1))
                    ttend(k) = qtend1
         elseif(k==1)then
                    omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                    qtend2 = omega_p*(theta_bg_int(k+1) - theta_bg_int(k))/(phalf_avg(k+1)-phalf_avg(k)) 
                    ttend(k) = qtend2
         endif
      elseif (scm_vert_adv_opt == 800) then
         qtend(k) = 0.0
         ttend(k) = 0.0
         if(k>1.and.k<kde-1)then
              if( omega_wtg(k) > 0 )then
                    qtend(k) = omega_wtg(k)*(qv_bg_int(k+1) - qv_bg_int(k))/(phalf_avg(k+1)-phalf_avg(k))
                    ttend(k) = omega_wtg(k)*(theta_bg_int(k+1) - theta_bg_int(k))/(phalf_avg(k+1)-phalf_avg(k))
              elseif( omega_wtg(k) < 0 )then
                    qtend(k) = omega_wtg(k)*(qv_bg_int(k) - qv_bg_int(k-1))/(phalf_avg(k)-phalf_avg(k-1))
                    ttend(k) = omega_wtg(k)*(theta_bg_int(k) - theta_bg_int(k-1))/(phalf_avg(k)-phalf_avg(k-1))
              endif
         elseif(k==1)then
              if( omega_wtg(k) > 0 )then
                    qtend(k) = omega_wtg(k)*(qv_bg_int(k+1) - qv_bg_int(k))/(phalf_avg(k+1)-phalf_avg(k))
                    ttend(k) = omega_wtg(k)*(theta_bg_int(k+1) - theta_bg_int(k))/(phalf_avg(k+1)-phalf_avg(k))
              endif
         endif
      endif
             

         DO ij = 1 , num_tiles; 
              DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
                 RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j) - qtend(k)
                 RTHBLTEN(i,k,j) = RTHBLTEN(i,k,j) - ttend(k)
            ENDDO; ENDDO
         ENDDO

         !the_wtg_forcing(k) = 0.0
         qv_wtg_forcing(k) = qv_wtg_forcing(k) -qtend(k)*DT  ! accumulated qv tendency
         rqvwtgten(k) = - qtend(k)                           ! wtg qv tendency 
         rthwtgten(k) = -ttend(k) 

      ENDDO ! end of the big do loop 


! ---------------------------------------------------------------------------------------------------------------

   ENDIF  ! end of big if block
 
  END SUBROUTINE wtg_driver0
   
!SW==

!
 
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





END MODULE module_wtg_driver0
