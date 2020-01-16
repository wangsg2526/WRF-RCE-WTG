!WRF:MEDIATION_LAYER:PHYSICS
!

MODULE module_vdiff_driver
CONTAINS

!------------------------------------------------------------------
   SUBROUTINE vdiff_driver(                                          &
                  itimestep,dt,u_frame,v_frame                     &
                 ,bldt,curr_secs,adapt_step_flag                   &
                 ,bldtacttime                                      & 
                 ,rublten,rvblten,rthblten                         &
                 ,tsk,xland,znt,ht                                 &
                 ,ust,pblh,hfx,qfx                          &
                 ,u_phy,v_phy,th_phy,rho                           &
                 ,p_phy,pi_phy,p8w,t_phy,dz8w,z                    &
                 ,exch_h,exch_m                          &
                 ,uz0,vz0,qsfc,f                          &
                 ,lowlyr,u10,v10,t2                                &
                 ,psim,psih,fm,fhh, wspd,br,chklowq         &
                 ,bl_pbl_physics, ra_lw_physics, dx                &
                 ,stepbl,warm_rain                                 &
                 ,ct,lh,snow,xice                       &
                 ,znu, znw, mut, p_top                             &
                 ,xkmv, xkhv  & !SW==
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,i_start,i_end, j_start,j_end, kts,kte, num_tiles &
             ! Optional
                 ,hol, mol, regime                                 &
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

   USE module_model_constants

! *** add new modules of schemes here

   !  pbl scheme:
   !      1. ysupbl
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
! QNC        cloud Liq number concentration (#/kg) !For CAMUWPBL scheme
! QNI        cloud ice number concentration (#/kg) !For CAMUWPBL scheme
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
!-- RQNIBLTEN     Qni tendency due to 
!                 PBL parameterization (#/kg/s) !For CAMUWPBL scheme
!-- id            WRF grid id  (optional, only needed by turbine drag schemes)
!-- itimestep     number of time steps
!-- GLW           downward long wave flux at ground surface (W/m^2)
!-- GSW           downward short wave flux at ground surface (W/m^2)
!-- EMISS         surface emissivity (between 0 and 1)
!-- TSK           surface temperature (K)
!-- TMN           soil temperature at lower boundary (K)
!-- XLAND         land mask (1 for land, 2 for water)
!-- ZNT           roughness length (m)
!-- MAVAIL        surface moisture availability (between 0 and 1)
!-- UST           u* in similarity theory (m/s)
!-- MOL           T* (similarity theory) (K)
!-- HOL           PBL height over Monin-Obukhov length
!-- PBLH          PBL height (m)
!-- CAPG          heat capacity for soil (J/K/m^3)
!-- THC           thermal inertia (Cal/cm/K/s^0.5)
!-- SNOWC         flag indicating snow coverage (1 for snow cover)
!-- HFX           upward heat flux at the surface (W/m^2)
!-- QFX           upward moisture flux at the surface (kg/m^2/s)
!-- REGIME        flag indicating PBL regime (stable, unstable, etc.)
!-- exch_m        exchange coefficient for momentum, m^2/s
!-- exch_h        exchange coefficient for heat, K m/s 
!-- exch_tke      exchange coeff. for TKE [enhanced], m^2/s (gbmpbl scheme)
!-- rthraten      tendency from radiation, used in GBM PBL scheme
!-- el_pbl        length scale from PBL schemes (m)
!-- wu_tur        turbulent flux of momentum (x) (m^2/s^2)
!-- wv_tur        turbulent flux of momentum (y) (m^2/s^2)
!-- wt_tur        turbulent flux of potential temperature  (K m/s)
!-- wq_tur        turbulent flux of water vapor  (- m/s)
!-- te_temf       Total energy from TEMF BL scheme
!-- km_temf       Exchange coefficient for momentum from TEMF BL scheme
!-- kh_temf       Exchange coefficient for heat from TEMF BL scheme
!-- shf_temf      Sensible heat flux from TEMF BL scheme
!-- qf_temf       Water vapor flux from TEMF BL scheme
!-- uw_temf       Momentum flux in U direction from TEMF BL scheme
!-- vw_temf       Momentum flux in V direction from TEMF BL scheme
!-- wupd_temf     Updraft velocity from TEMF BL scheme
!-- mf_temf       Mass flux from TEMF BL scheme
!-- thup_temf     Updraft thetal from TEMF BL scheme
!-- qtup_temf     Updraft qt from TEMF BL scheme
!-- qlup_temf     Updraft ql from TEMF BL scheme
!-- cf3d_temf     3D cloud fraction from TEMF PBL
!-- cfm_temf      Column cloud fraction from TEMF PBL
!-- exch_temf     Surface exchange coefficient (as for moisture) from TEMF surface layer scheme
!-- flhc          Surface exchange coefficient for heat (for TEMF)
!-- flqc          Surface exchange coefficient for moisture (for TEMF)
!-- thz0          potential temperature at roughness length (K)
!-- uz0           u wind component at roughness length (m/s)
!-- vz0           v wind component at roughness length (m/s)
!-- qsfc          specific humidity at lower boundary (kg/kg)
!-- th2           diagnostic 2-m theta from surface layer and lsm
!-- t2            diagnostic 2-m temperature from surface layer and lsm
!-- q2            diagnostic 2-m mixing ratio from surface layer and lsm
!-- lowlyr        index of lowest model layer above ground
!-- rr            dry air density (kg/m^3)
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
!-- TSLB          
!-- ZS
!-- DZS
!-- num_soil_layers number of soil layer
!-- IFSNOW      ifsnow=1 for snow-cover effects
!-- z_at_w      Height above sea level at layer interfaces (m) 
!-- cldfra      Cloud fraction [unitless]
!-- cldfra_old_mp      Cloud fraction [unitless]
!-- rthratenlw  Tendency for LW ( K/s)
!-- tauresx2d   X-COMP OF RESIDUAL STRESS(m^2/s^2)
!-- tauresy2d   Y-COMP OF RESIDUAL STRESS(m^2/s^2)
!-- tpert2d     Convective temperature excess (K)
!-- qpert2d     Convective humidity excess (kg/kg)
!-- wpert2d     Turbulent velocity excess (m/s)
!-- wsedl3d     Sedimentation velocity of stratiform liquid cloud droplet (m/s)
!-- turbtype3d  Turbulent interface types [ no unit ]  
!-- smaw3d      Normalized Galperin instability function for momentum  ( 0<= <=4.964 and 1 at neutral ) [no units]
!
!-- P_QV          species index for water vapor
!-- P_QC          species index for cloud water
!-- P_QR          species index for rain water
!-- P_QI          species index for cloud ice
!-- P_QNC         species index for cloud liq number concentration !For CAMUWPBL scheme
!-- P_QNI         species index for cloud ice number concentration !For CAMUWPBL scheme
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
   INTEGER,    DIMENSION( ims:ime , jms:jme ),                    &
               INTENT(IN   )    ::                        LOWLYR
!
   LOGICAL,      INTENT(IN   )    ::   warm_rain


   REAL,       DIMENSION( kms:kme ),                              &
               OPTIONAL, INTENT(IN   )    ::               znu,   &
                                                           znw
!
   REAL,       INTENT(IN   )    ::     DT,DX
   REAL,       INTENT(IN   ),OPTIONAL    ::     bldt
   REAL,       INTENT(IN   ),OPTIONAL    ::     curr_secs
   LOGICAL,    INTENT(IN   ),OPTIONAL    ::     adapt_step_flag
   REAL,       INTENT(INOUT),OPTIONAL    ::     bldtacttime  

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
   REAL,       DIMENSION( ims:ime , jms:jme ),                    &
               INTENT(IN   )    ::                         XLAND, &
                                                              HT, &
                                                            PSIM, &
                                                            PSIH, &
                                                              FM, &
                                                             FHH, &
                                                              BR, &
                                                               F, &
                                                         CHKLOWQ
!
   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               INTENT(INOUT)    ::                           TSK, &
                                                             UST, &
                                                            PBLH, &
                                                             HFX, &
                                                             QFX, &
                                                             ZNT, &
                                                            QSFC, &
                                                             UZ0, &
                                                             VZ0, &
                                                              CT, &
                                                             U10, &
                                                             V10, &
                                                              T2, &
                                                            WSPD

!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(INOUT)    ::                       RUBLTEN, &
                                                         RVBLTEN, &
                                                        RTHBLTEN, &
                                           EXCH_H,EXCH_M


!

   REAL ,                             INTENT(IN   )  ::  u_frame, &
                                                         v_frame
!


   REAL,       DIMENSION( ims:ime , jms:jme ),                    &
               INTENT(IN)    :: XICE, SNOW, LH


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
                     ,rqiblten,rqsblten,rqgblten        !rqniblten  used in CAMUWPBL


   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               OPTIONAL, INTENT(INOUT)  ::   HOL, MOL, REGIME


   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               OPTIONAL                                         , &
               INTENT(IN)    ::                           mut
!
   REAL,       OPTIONAL, INTENT(IN)    ::               p_top
!
!
!

   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &   !SW==
             intent(in   )   ::                                          xkmv,xkhv !SW==  



!  LOCAL  VAR

   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ) ::v_phytmp
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ) ::u_phytmp

   REAL,       DIMENSION( ims:ime, jms:jme )          ::  TSKOLD, &
                                                          USTOLD, &
                                                          ZNTOLD, &
                                                             ZOL, &
                                                            PSFC

! make these allocatable depending on the setting of idiff
! Typically, we try to avoide allocating and deallocating local storage like this
! so as not to fragment the stack. But at this point, the idiff = 1 case is disabled
! (set to 0 for all cases) and has to be set manually by users who want to work with
! it.  When it becomes a more standard option, this should be redone, either defining
! these as state with package clauses to turn them on and off and passing them in,
! or pass in an integer flag that can be used to dimension the arrays to 1:1:1 as
! local variables.  JM 20100316

   REAL, ALLOCATABLE, DIMENSION( :, :, : )::a_u        ! Implicit component for the momemtum in X-direction
   REAL, ALLOCATABLE, DIMENSION( :, :, : )::a_v        ! Implicit component for the momemtum in Y-direction
   REAL, ALLOCATABLE, DIMENSION( :, :, : )::a_t        ! Implicit component for the Pot. Temp.
   REAL, ALLOCATABLE, DIMENSION( :, :, : )::a_q        ! Implicit component for the water vapor

   REAL, ALLOCATABLE, DIMENSION( :, :, : )::b_u        ! Explicit component for the momemtum in X-direction
   REAL, ALLOCATABLE, DIMENSION( :, :, : )::b_v        ! Explicit component for the momemtum in Y-direction
   REAL, ALLOCATABLE, DIMENSION( :, :, : )::b_t        ! Explicit component for the Pot. Temp.
   REAL, ALLOCATABLE, DIMENSION( :, :, : )::b_q        ! Explicit component for the water vapor

   REAL, ALLOCATABLE, DIMENSION( :, :, : )::sf         ! surfaces
   REAL, ALLOCATABLE, DIMENSION( :, :, : )::vl         ! volumes

   REAL    :: DTMIN,DTBL
!
   INTEGER :: initflag
!
   INTEGER :: i,J,K,NK,jj,ij,its,ite,jts,jte
   LOGICAL :: radiation
   LOGICAL :: flag_bep
   LOGICAL :: flag_myjsfc
   LOGICAL :: flag_qv, flag_qc, flag_qr, flag_qi, flag_qs, flag_qg, flag_qnc,flag_qni !flag_qnc,flag_qnc are used in camuwpbl scheme
   CHARACTER*256 :: message
   REAL    :: next_bl_time
   LOGICAL :: run_param , doing_adapt_dt , decided
   LOGICAL :: do_adapt
   integer iu_bep,iurb,idiff
   real seamask,thsk,zzz,unew,vnew,tnew,qnew,umom,vmom
   REAL :: z0,z1,z2,w1,w2

!------------------------------------------------------------------
!
!!!!!!!if using BEP set flag_bep to true

!
  flag_qv = .FALSE. ; IF ( PRESENT( F_QV ) ) flag_qv = F_QV
  flag_qc = .FALSE. ; IF ( PRESENT( F_QC ) ) flag_qc = F_QC
  flag_qr = .FALSE. ; IF ( PRESENT( F_QR ) ) flag_qr = F_QR
  flag_qi = .FALSE. ; IF ( PRESENT( F_QI ) ) flag_qi = F_QI
  flag_qs = .FALSE. ; IF ( PRESENT( F_QS ) ) flag_qs = F_QS
  flag_qg = .FALSE. ; IF ( PRESENT( F_QG ) ) flag_qg = F_QG

!SW==
  !if (bl_pbl_physics .eq. 0) return
!SW==

! RAINBL in mm (Accumulation between PBL calls)
!
   doing_adapt_dt = .FALSE.
   IF ( PRESENT(adapt_step_flag) ) THEN
      IF ( adapt_step_flag ) THEN
         doing_adapt_dt = .TRUE.
         IF ( bldtacttime .eq. 0. ) THEN
            bldtacttime = CURR_SECS + bldt*60.
         END IF
      END IF
   END IF

!  Do we run through this scheme or not?

!    Test 1:  If this is the initial model time, then yes.
!                ITIMESTEP=1
!    Test 2:  If the user asked for the pbl to be run every time step, then yes.
!                BLDT=0 or STEPBL=1
!    Test 3:  If not adaptive dt, and this is on the requested pbl frequency, then yes.
!                MOD(ITIMESTEP,STEPBL)=0
!    Test 4:  If using adaptive dt and the current time is past the last requested activate pbl time, then yes.
!                CURR_SECS >= BLDTACTTIME

!  If we do run through the scheme, we set the flag run_param to TRUE and we set the decided flag
!  to TRUE.  The decided flag says that one of these tests was able to say "yes", run the scheme.
!  We only proceed to other tests if the previous tests all have left decided as FALSE.

!  If we set run_param to TRUE and this is adaptive time stepping, we set the time to the next
!  pbl run.

   run_param = .FALSE.
   decided = .FALSE.
   IF ( ( .NOT. decided ) .AND. &
        ( itimestep .EQ. 1 ) ) THEN
      run_param   = .TRUE.
      decided     = .TRUE.
   END IF

   IF ( PRESENT(bldt) )THEN
      IF ( ( .NOT. decided ) .AND. &
           ( ( bldt .EQ. 0. ) .OR. ( stepbl .EQ. 1 ) ) ) THEN
         run_param   = .TRUE.
         decided     = .TRUE.
      END IF
   ELSE
      IF ( ( .NOT. decided ) .AND. &
                                   ( stepbl .EQ. 1 )   ) THEN
         run_param   = .TRUE.
         decided     = .TRUE.
      END IF
   END IF

   IF ( ( .NOT. decided ) .AND. &
        ( .NOT. doing_adapt_dt ) .AND. &
        ( MOD(itimestep,stepbl) .EQ. 0 ) ) THEN
      run_param   = .TRUE.
      decided     = .TRUE.
   END IF

   IF ( ( .NOT. decided ) .AND. &
        ( doing_adapt_dt ) .AND. &
        ( curr_secs .GE. bldtacttime ) ) THEN
      run_param   = .TRUE.
      decided     = .TRUE.
      bldtacttime = curr_secs + bldt*60
   END IF


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
            IF ( curr_secs .LT. 2. * dt ) THEN
               call wrf_message("WARNING: When using an adaptive time-step the boundary layer"// &
                                " time-step should be 0 (i.e., equivalent to model time-step).  ")
               call wrf_message("In order to proceed, for boundary layer calculations, the "// &
                                "boundary layer time-step"// &
                                 " will be rounded to the nearest minute," )
                call wrf_message("possibly resulting in innacurate results.")
            END IF
            DTBL=bldt*60
         else
            DTBL=DT*STEPBL
         endif
      endif
   else
      DTBL=DT*STEPBL
   endif

! Added this else clause so that idiff is always initialized regardless of which core we are using
       idiff=0

   IF ( idiff .EQ. 1 ) THEN
     ALLOCATE (a_u(ims:ime,kms:kme,jms:jme))       ! Implicit component for the momemtum in X-direction
     ALLOCATE (a_v(ims:ime,kms:kme,jms:jme))       ! Implicit component for the momemtum in Y-direction
     ALLOCATE (a_t(ims:ime,kms:kme,jms:jme))       ! Implicit component for the Pot. Temp.
     ALLOCATE (a_q(ims:ime,kms:kme,jms:jme))       ! Implicit component for the water vapor
     ALLOCATE (b_u(ims:ime,kms:kme,jms:jme))       ! Explicit component for the momemtum in X-direction
     ALLOCATE (b_v(ims:ime,kms:kme,jms:jme))       ! Explicit component for the momemtum in Y-direction
     ALLOCATE (b_t(ims:ime,kms:kme,jms:jme))       ! Explicit component for the Pot. Temp.
     ALLOCATE (b_q(ims:ime,kms:kme,jms:jme))       ! Explicit component for the water vapor
     ALLOCATE (sf(ims:ime,kms:kme,jms:jme) )       ! surfaces
     ALLOCATE (vl(ims:ime,kms:kme,jms:jme) )       ! volumes
   ENDIF
   
! SAVE OLD VALUES

   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij,i,j,k )

   DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         TSKOLD(i,j)=TSK(i,j)
         USTOLD(i,j)=UST(i,j)
         ZNTOLD(i,j)=ZNT(i,j)

! REVERSE ORDER IN THE VERTICAL DIRECTION

! testing change later

         DO k=kts,kte
            v_phytmp(i,k,j)=v_phy(i,k,j)+v_frame
            u_phytmp(i,k,j)=u_phy(i,k,j)+u_frame
         ENDDO

! PSFC : in Pa

         PSFC(I,J)=p8w(I,kms,J)

         DO k=kts,min(kte+1,kde)
            RTHBLTEN(I,K,J)=0.
            RUBLTEN(I,K,J)=0.
            RVBLTEN(I,K,J)=0.
            IF ( PRESENT( RQCBLTEN )) RQCBLTEN(I,K,J)=0.
            IF ( PRESENT( RQVBLTEN )) RQVBLTEN(I,K,J)=0.
         ENDDO

         IF (flag_QI .AND. PRESENT(RQIBLTEN) ) THEN
            DO k=kts,min(kte+1,kde)
               RQIBLTEN(I,K,J)=0.
            ENDDO
         ENDIF
         !Following if condition is added for CAMUWPBL scheme
      ENDDO
      ENDDO


      
   ENDDO
   !$OMP END PARALLEL DO
!
  !$OMP PARALLEL DO   &
  !$OMP PRIVATE ( ij, i,j,k, its, ite, jts, jte, z0, z1, z2, w1, w2, message, initflag )
  DO ij = 1 , num_tiles

   its = i_start(ij)
   ite = i_end(ij)
   jts = j_start(ij)
   jte = j_end(ij)

        CALL wrf_debug(100,'in YSU PBL')
           IF ( PRESENT( qv_curr )  .AND. PRESENT( qc_curr )  .AND. &
                PRESENT( qi_curr )                            .AND. &
                PRESENT( rqvblten ) .AND. PRESENT( rqcblten ) .AND. &
                PRESENT( rqiblten )                           .AND. &
                PRESENT( hol      ) ) THEN
!
             CALL vdiff(                                              &
               U3D=u_phytmp,V3D=v_phytmp,TH3D=th_phy,T3D=t_phy      &
              ,QV3D=qv_curr,QC3D=qc_curr,QI3D=qi_curr               &
              ,P3D=p_phy,P3DI=p8w,PI3D=pi_phy                       &
              ,RUBLTEN=rublten,RVBLTEN=rvblten                      &
              ,RTHBLTEN=rthblten,RQVBLTEN=rqvblten                  &
              ,RQCBLTEN=rqcblten,RQIBLTEN=rqiblten                  &
              ,FLAG_QI=flag_qi                                      &
              ,XKMV=xkmv, XKHV=xkhv  & !SW==
              ,CP=cp,G=g,ROVCP=rcp,RD=r_D,ROVG=rovg                 &
              ,DZ8W=dz8w,XLV=XLV,RV=r_v,PSFC=PSFC                   &
              ,ZNU=znu,ZNW=znw,MUT=mut,P_TOP=p_top                  &
              ,ZNT=znt,UST=ust,HPBL=pblh                            &
              ,PSIM=fm,PSIH=fhh,XLAND=xland                         &
              ,HFX=hfx,QFX=qfx                                      &
              ,U10=u10,V10=v10                                      &
              ,WSPD=wspd,BR=br,DT=dtbl               &
              ,EP1=ep_1,EP2=ep_2,KARMAN=karman                      &
              ,EXCH_H=exch_h,REGIME=regime                          &
              ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde      &
              ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme      &
              ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte      &
                                                                    )
           ELSE
               WRITE ( message , FMT = '(A,7(L1,1X))' )             &
                 'present: '//                                      &
                 'qv_curr, '//                                      &
                 'qc_curr, '//                                      &
                 'qi_curr, '//                                      &
                 'rqvblten, '//                                     &
                 'rqcblten, '//                                     &
                 'rqiblten, '//                                     &
                 'hol = ' ,                                         &
                  PRESENT( qv_curr ) ,                              &
                  PRESENT( qc_curr ) ,                              &
                  PRESENT( qi_curr ) ,                              &
                  PRESENT( rqvblten ) ,                             &
                  PRESENT( rqcblten ) ,                             &
                  PRESENT( rqiblten ) ,                             &
                  PRESENT( hol      )
               CALL wrf_debug(0,message)
               CALL wrf_error_fatal('Lack arguments to call YSU pbl')
           ENDIF





      
   ENDDO
   !$OMP END PARALLEL DO

   ENDIF

   END SUBROUTINE vdiff_driver

!=============================================================================


   subroutine vdiff(u3d,v3d,th3d,t3d,qv3d,qc3d,qi3d,p3d,p3di,pi3d,               &
                  rublten,rvblten,rthblten,                                    &
                  rqvblten,rqcblten,rqiblten,flag_qi,                          &
                  xkmv, xkhv, & !SW==
                  cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,                    &
                  dz8w,psfc,                                                   &
                  znu,znw,mut,p_top,                                           &
                  znt,ust,hpbl,psim,psih,                                      &
                  xland,hfx,qfx,wspd,br,                                       &
                  dt,                                                  &
                  exch_h,                                                      &
                  u10,v10,                                                     &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                !optional
                  regime                                           )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- th3d        3d potential temperature (k)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- qc3d        3d cloud mixing ratio (kg/kg)
!-- qi3d        3d ice mixing ratio (kg/kg)
!               (note: if P_QI<PARAM_FIRST_SCALAR this should be zero filled)
!-- p3d         3d pressure (pa)
!-- p3di        3d pressure (pa) at interface level
!-- pi3d        3d exner function (dimensionless)
!-- rr3d        3d dry air density (kg/m^3)
!-- rublten     u tendency due to
!               pbl parameterization (m/s/s)
!-- rvblten     v tendency due to
!               pbl parameterization (m/s/s)
!-- rthblten    theta tendency due to
!               pbl parameterization (K/s)
!-- rqvblten    qv tendency due to
!               pbl parameterization (kg/kg/s)
!-- rqcblten    qc tendency due to
!               pbl parameterization (kg/kg/s)
!-- rqiblten    qi tendency due to
!               pbl parameterization (kg/kg/s)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g           acceleration due to gravity (m/s^2)
!-- rovcp       r/cp
!-- rd          gas constant for dry air (j/kg/k)
!-- rovg        r/g
!-- dz8w        dz between full levels (m)
!-- xlv         latent heat of vaporization (j/kg)
!-- rv          gas constant for water vapor (j/kg/k)
!-- psfc        pressure at the surface (pa)
!-- znt         roughness length (m)
!-- ust         u* in similarity theory (m/s)
!-- hpbl        pbl height (m)
!-- psim        similarity stability function for momentum
!-- psih        similarity stability function for heat
!-- xland       land mask (1 for land, 2 for water)
!-- hfx         upward heat flux at the surface (w/m^2)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- wspd        wind speed at lowest model level (m/s)
!-- u10         u-wind speed at 10 m (m/s)
!-- v10         v-wind speed at 10 m (m/s)
!-- br          bulk richardson number in surface layer
!-- dt          time step (s)
!-- rvovrd      r_v divided by r_d (dimensionless)
!-- ep1         constant for virtual temperature (r_v/r_d - 1)
!-- ep2         constant for specific humidity calculation
!-- karman      von karman constant
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------------------
!
   integer,parameter ::  ndiff = 3
   real,parameter    ::  rcl = 1.0
!
   integer,  intent(in   )   ::      ids,ide, jds,jde, kds,kde,                &
                                     ims,ime, jms,jme, kms,kme,                &
                                     its,ite, jts,jte, kts,kte
!
   real,     intent(in   )   ::      dt,cp,g,rovcp,rovg,rd,xlv,rv
!
   real,     intent(in )     ::      ep1,ep2,karman
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in   )   ::                                          qv3d, &
                                                                         qc3d, &
                                                                         qi3d, &
                                                                          p3d, &
                                                                         pi3d, &
                                                                         th3d, &
                                                                          t3d, &
                                                                         dz8w
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in   )   ::                                          p3di
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &   !SW==
             intent(in   )   ::                                          xkmv,xkhv !SW==  
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(inout)   ::                                       rublten, &
                                                                      rvblten, &
                                                                     rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(inout)   ::                                        exch_h
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(inout)   ::                                           u10, &
                                                                          v10
!
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(in   )   ::                                         xland, &
                                                                          hfx, &
                                                                          qfx, &
                                                                           br, &
                                                                         psfc
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(in   )   ::                                                &
                                                                         psim, &
                                                                         psih
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(inout)   ::                                           znt, &
                                                                          ust, &
                                                                         hpbl, &
                                                                          wspd
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in   )   ::                                           u3d, &
                                                                          v3d
!
   logical,  intent(in)      ::                                       flag_qi
!
!optional
!
   real,     dimension( ims:ime, jms:jme )                                   , &
             optional                                                        , &
             intent(inout)   ::                                        regime
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             optional                                                        , &
             intent(inout)   ::                                       rqiblten
!
   real,     dimension( kms:kme )                                            , &
             optional                                                        , &
             intent(in   )   ::                                           znu, &
                                                                          znw
!
   real,     dimension( ims:ime, jms:jme )                                   , &
             optional                                                        , &
             intent(in   )   ::                                           mut
!
   real,     optional, intent(in   )   ::                               p_top
!
!local
   integer ::  i,j,k
   real,     dimension( its:ite, kts:kte*ndiff )  ::                 rqvbl2dt, &
                                                                         qv2d
   real,     dimension( its:ite, kts:kte )  ::                            pdh
   real,     dimension( its:ite, kts:kte )  ::                            xkmv2d, xkhv2d !SW==
 
   real,     dimension( its:ite, kts:kte+1 )  ::                         pdhi
   real,     dimension( its:ite )  ::                                          &
                                                                        dusfc, &
                                                                        dvsfc, &
                                                                        dtsfc, &
                                                                        dqsfc
!
   qv2d(its:ite,:) = 0.0
!
   do j = jts,jte
     if(present(mut))then
!
! For ARW we will replace p and p8w with dry hydrostatic pressure
!
        do k = kts,kte+1
          do i = its,ite
             if(k.le.kte)pdh(i,k) = mut(i,j)*znu(k) + p_top
             pdhi(i,k) = mut(i,j)*znw(k) + p_top
          enddo
        enddo
      else
        do k = kts,kte+1
          do i = its,ite
            if(k.le.kte)pdh(i,k) = p3d(i,k,j)
            pdhi(i,k) = p3di(i,k,j)
          enddo
        enddo
      endif
      do k = kts,kte
        do i = its,ite
          qv2d(i,k) = qv3d(i,k,j)
          qv2d(i,k+kte) = qc3d(i,k,j)
          !xkmv2d(i,k) = xkmv(i,k,j) !SW==
          !xkhv2d(i,k) = xkhv(i,k,j) !SW==
          if(present(rqiblten)) qv2d(i,k+kte+kte) = qi3d(i,k,j)
        enddo
      enddo
!
      call vdiff2d(J=j,ux=u3d(ims,kms,j),vx=v3d(ims,kms,j)                       &
              ,tx=t3d(ims,kms,j)                                               &
              ,qx=qv2d(its,kts)                                                &
              ,p2d=pdh(its,kts),p2di=pdhi(its,kts)                             &
              ,pi2d=pi3d(ims,kms,j)                                            &
              ,xkmv2d=xkmv(ims,kms,j)                                            & !SW==
              ,xkhv2d=xkhv(ims,kms,j)                                            & !SW==
              ,utnp=rublten(ims,kms,j),vtnp=rvblten(ims,kms,j)                 &
              ,ttnp=rthblten(ims,kms,j),qtnp=rqvbl2dt(its,kts),ndiff=ndiff     &
              ,cp=cp,g=g,rovcp=rovcp,rd=rd,rovg=rovg                           &    
              ,xlv=xlv,rv=rv                                                   &
              ,ep1=ep1,ep2=ep2,karman=karman                                   &
              ,dz8w2d=dz8w(ims,kms,j)                                          &
              ,psfcpa=psfc(ims,j),znt=znt(ims,j),ust=ust(ims,j)                &
              ,hpbl=hpbl(ims,j)                                                &
              ,regime=regime(ims,j),psim=psim(ims,j)                           &
              ,psih=psih(ims,j),xland=xland(ims,j)                             &
              ,hfx=hfx(ims,j),qfx=qfx(ims,j)                                   &
              ,wspd=wspd(ims,j),br=br(ims,j)                                   &
              ,dusfc=dusfc,dvsfc=dvsfc,dtsfc=dtsfc,dqsfc=dqsfc                 &
              ,dt=dt,rcl=1.0                              &
              ,exch_hx=exch_h(ims,kms,j)                                       &
              ,u10=u10(ims,j),v10=v10(ims,j)                                   &
              ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde               &
              ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme               &
              ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte   )
!
     do k = kts,kte
       do i = its,ite
         rthblten(i,k,j) = rthblten(i,k,j)/pi3d(i,k,j)
         rqvblten(i,k,j) = rqvbl2dt(i,k)
         rqcblten(i,k,j) = rqvbl2dt(i,k+kte)
         if(present(rqiblten)) rqiblten(i,k,j) = rqvbl2dt(i,k+kte+kte)
       enddo
     enddo
!
   enddo
!
   end subroutine vdiff 
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
   subroutine vdiff2d(j,ux,vx,tx,qx,p2d,p2di,pi2d,                               &
                  xkmv2d, xkhv2d,  & !SW==
                  utnp,vtnp,ttnp,qtnp,ndiff,                                   &
                  cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,                    &
                  dz8w2d,psfcpa,                                               &
                  znt,ust,hpbl,psim,psih,                                      &
                  xland,hfx,qfx,wspd,br,                                       &
                  dusfc,dvsfc,dtsfc,dqsfc,                                     &
                  dt,rcl,                                               &
                  exch_hx,                                                     &
                  u10,v10,                                                     &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                !optional
                  regime )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!     this code is a revised vertical diffusion package ("ysupbl")
!     with a nonlocal turbulent mixing in the pbl after "mrfpbl".
!     the ysupbl (hong et al. 2006) is based on the study of noh
!     et al.(2003) and accumulated realism of the behavior of the
!     troen and mahrt (1986) concept implemented by hong and pan(1996).
!     the major ingredient of the ysupbl is the inclusion of an explicit
!     treatment of the entrainment processes at the entrainment layer.
!     this routine uses an implicit approach for vertical flux
!     divergence and does not require "miter" timesteps.
!     it includes vertical diffusion in the stable atmosphere
!     and moist vertical diffusion in clouds.
!
!     mrfpbl:
!     coded by song-you hong (ncep), implemented by jimy dudhia (ncar)
!              fall 1996
!
!     ysupbl:
!     coded by song-you hong (yonsei university) and implemented by
!              song-you hong (yonsei university) and jimy dudhia (ncar)
!              summer 2002
!
!     further modifications :
!              an enhanced stable layer mixing, april 2008
!              ==> increase pbl height when sfc is stable (hong 2010)
!              pressure-level diffusion, april 2009
!              ==> negligible differences
!              implicit forcing for momentum with clean up, july 2009
!              ==> prevents model blowup when sfc layer is too low
!              incresea of lamda, maximum (30, 0.1 x del z) feb 2010
!              ==> prevents model blowup when delz is extremely large
!              revised prandtl number at surface, peggy lemone, feb 2010
!              ==> increase kh, decrease mixing due to counter-gradient term
!              revised thermal, shin et al. mon. wea. rev. , songyou hong, aug 2011
!              ==> reduce the thermal strength when z1 < 0.1 h
!              revised prandtl number for free convection, dudhia, mar 2012
!              ==> pr0 = 1 + bke (=0.272) when newtral, kh is reduced
!              minimum kzo = 0.01, lo = min (30m,delz), hong, mar 2012
!              ==> weaker mixing when stable, and les resolution in vertical
!              gz1oz0 is removed, and phim phih are ln(z1/z0)-phim,h, hong, mar 2012
!              ==> consider thermal z0 when differs from mechanical z0
!              a bug fix in wscale computation in stable bl, sukanta basu, jun 2012
!              ==> wscale becomes small with height, and less mixing in stable bl
!
!     references:
!
!        hong (2010) quart. j. roy. met. soc
!        hong, noh, and dudhia (2006), mon. wea. rev.
!        hong and pan (1996), mon. wea. rev.
!        noh, chun, hong, and raasch (2003), boundary layer met.
!        troen and mahrt (1986), boundary layer met.
!
!-------------------------------------------------------------------------------
!
   real,parameter    ::  xkzmin = 0.01,xkzmax = 1000.,rimin = -100.
   real,parameter    ::  rlam = 30.,prmin = 0.25,prmax = 4.
   real,parameter    ::  brcr_ub = 0.0,brcr_sb = 0.25,cori = 1.e-4
   real,parameter    ::  afac = 6.8,bfac = 6.8,pfac = 2.0,pfac_q = 2.0
   real,parameter    ::  phifac = 8.,sfcfrac = 0.1
   real,parameter    ::  d1 = 0.02, d2 = 0.05, d3 = 0.001
   real,parameter    ::  h1 = 0.33333335, h2 = 0.6666667
   real,parameter    ::  ckz = 0.001,zfmin = 1.e-8,aphi5 = 5.,aphi16 = 16.
   real,parameter    ::  tmin=1.e-2
   real,parameter    ::  gamcrt = 3.,gamcrq = 2.e-3
   real,parameter    ::  xka = 2.4e-5
   integer,parameter ::  imvdif = 1
!
   integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,                 &
                                    ims,ime, jms,jme, kms,kme,                 &
                                    its,ite, jts,jte, kts,kte,                 &
                                    j,ndiff
!
   real,     intent(in   )   ::     dt,rcl,cp,g,rovcp,rovg,rd,xlv,rv
!
   real,     intent(in )     ::     ep1,ep2,karman
!
   real,     dimension( ims:ime, kms:kme ),                                    &
             intent(in)      ::                                        dz8w2d, &
                                                                         pi2d

   real,     dimension( ims:ime, kms:kme ),                                    &  !SW==
             intent(in)      ::                                        xkmv2d, xkhv2d !SW== 
!
   real,     dimension( ims:ime, kms:kme )                                   , &
             intent(in   )   ::                                            tx
   real,     dimension( its:ite, kts:kte*ndiff )                             , &
             intent(in   )   ::                                            qx
!
   real,     dimension( ims:ime, kms:kme )                                   , &
             intent(inout)   ::                                          utnp, &
                                                                         vtnp, &
                                                                         ttnp
   real,     dimension( its:ite, kts:kte*ndiff )                             , &
             intent(inout)   ::                                          qtnp
!
   real,     dimension( its:ite, kts:kte+1 )                                 , &
             intent(in   )   ::                                          p2di
!
   real,     dimension( its:ite, kts:kte )                                   , &
             intent(in   )   ::                                           p2d
!
!
   real,     dimension( ims:ime )                                            , &
             intent(inout)   ::                                           ust, &
                                                                         hpbl, &
                                                                          znt
   real,     dimension( ims:ime )                                            , &
             intent(in   )   ::                                         xland, &
                                                                          hfx, &
                                                                          qfx
!
   real,     dimension( ims:ime ), intent(inout)   ::                    wspd
   real,     dimension( ims:ime ), intent(in  )    ::                      br
!
   real,     dimension( ims:ime ), intent(in   )   ::                    psim, &
                                                                         psih
!
   real,     dimension( ims:ime ), intent(in   )   ::                  psfcpa
!
   real,     dimension( ims:ime, kms:kme )                                   , &
             intent(in   )   ::                                            ux, &
                                                                           vx
   real,     dimension( ims:ime )                                            , &
             optional                                                        , &
             intent(inout)   ::                                        regime
!
! local vars
!
   real,     dimension( its:ite )            ::                           hol
   real,     dimension( its:ite, kts:kte+1 ) ::                            zq
!
   real,     dimension( its:ite, kts:kte )   ::                                &
                                                                     thx,thvx, &
                                                                          del, &
                                                                          dza, &
                                                                          dzq, &
                                                                         xkzo, &
                                                                           za
!
   real,    dimension( its:ite )             ::                                &
                                                                         rhox, &
                                                                       govrth, &
                                                                  zl1,thermal, &
                                                                       wscale, &
                                                                    brdn,brup, &
                                                                    phim,phih, &
                                                                  dusfc,dvsfc, &
                                                                  dtsfc,dqsfc, &
                                                                        prpbl, &
                                                                        wspd1
!
   real,    dimension( its:ite, kts:kte )    ::                     xkzm,xkzh, &
                                                                        f1,f2, &
                                                                        r1,r2, &
                                                                        ad,au, &
                                                                           cu, &
                                                                           al, &
                                                                         xkzq, &
                                                                         zfac
!
!jdf added exch_hx
!
   real,    dimension( ims:ime, kms:kme )                                    , &
            intent(inout)   ::                                        exch_hx
!
   real,    dimension( ims:ime )                                             , &
            intent(inout)    ::                                           u10, &
                                                                          v10
   real,    dimension( its:ite )    ::                                         &
                                                                         brcr, &
                                                                        sflux, &
                                                                         zol1, &
                                                                    brcr_sbro
!
   real,    dimension( its:ite, kts:kte, ndiff)  ::                     r3,f3
!
   logical, dimension( its:ite )             ::                        pblflg, &
                                                                       sfcflg, &
                                                                       stable
!
   integer ::  n,i,k,l,ic,is
   integer ::  klpbl, ktrace1, ktrace2, ktrace3
!
!
   real    ::  dt2,rdt,spdk2,fm,fh,hol1,gamfac,vpert,prnum,prnum0
   real    ::  ss,ri,qmean,tmean,alph,chi,zk,rl2,dk,sri
   real    ::  brint,dtodsd,dtodsu,rdz,dsdzt,dsdzq,dsdz2,rlamdz
   real    ::  utend,vtend,ttend,qtend
   real    ::  dtstep,govrthv
   real    ::  cont, conq, conw, conwrc
!

   real, dimension( its:ite, kts:kte )     ::                         wscalek
   real, dimension( its:ite, kts:kte )     ::                     xkzml,xkzhl, &
                                                               zfacent,entfac
   real, dimension( its:ite )              ::                            ust3, &
                                                                       wstar3, &
                                                                      wm2, we, &
                                                                       bfxpbl, &
                                                                hfxpbl,qfxpbl, &
                                                                ufxpbl,vfxpbl, &
                                                                        dthvx
   real    ::  prnumfac,bfx0,hfx0,qfx0,delb,dux,dvx,                           &
               dsdzu,dsdzv,wm3,dthx,dqx,wspd10,ross,tem1,dsig,tvcon,conpr,     &
               prfac,prfac2,phim8z
!
!-------------------------------------------------------------------------------
!
   klpbl = kte
!
   cont=cp/g
   conq=xlv/g
   conw=1./g
   conwrc = conw*sqrt(rcl)
   conpr = bfac*karman*sfcfrac
!
!  k-start index for tracer diffusion
!
   ktrace1 = 0
   ktrace2 = 0 + kte
   ktrace3 = 0 + kte*2
!
   do k = kts,kte
     do i = its,ite
       thx(i,k) = tx(i,k)/pi2d(i,k)
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       tvcon = (1.+ep1*qx(i,k))
       thvx(i,k) = thx(i,k)*tvcon
     enddo
   enddo
!
   do i = its,ite
     tvcon = (1.+ep1*qx(i,1))
     rhox(i) = psfcpa(i)/(rd*tx(i,1)*tvcon)
     govrth(i) = g/thx(i,1)
   enddo
!
!-----compute the height of full- and half-sigma levels above ground
!     level, and the layer thicknesses.
!
   do i = its,ite
     zq(i,1) = 0.
   enddo
!
   do k = kts,kte
     do i = its,ite
       zq(i,k+1) = dz8w2d(i,k)+zq(i,k)
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       za(i,k) = 0.5*(zq(i,k)+zq(i,k+1))
       dzq(i,k) = zq(i,k+1)-zq(i,k)
       del(i,k) = p2di(i,k)-p2di(i,k+1)
     enddo
   enddo
!
   do i = its,ite
     dza(i,1) = za(i,1)
   enddo
!
   do k = kts+1,kte
     do i = its,ite
       dza(i,k) = za(i,k)-za(i,k-1)
     enddo
   enddo
!
!
!-----initialize vertical tendencies and
!
   utnp(its:ite,:) = 0.
   vtnp(its:ite,:) = 0.
   ttnp(its:ite,:) = 0.
   qtnp(its:ite,:) = 0.
!
   do i = its,ite
     wspd1(i) = sqrt(ux(i,1)*ux(i,1)+vx(i,1)*vx(i,1))+1.e-9
   enddo
!
!-------------------------------------------------------------------------------------------
!---- compute vertical diffusion
!-------------------------------------------------------------------------------------------
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     compute preliminary variables
!
   dtstep = dt
   dt2 = 2.*dtstep
   rdt = 1./dt2
!
!
   do k = kts,klpbl-1
     do i = its,ite
       xkzo(i,k) = ckz*dza(i,k+1)
     enddo
   enddo
!
   do i = its,ite
     dusfc(i) = 0.
     dvsfc(i) = 0.
     dtsfc(i) = 0.
     dqsfc(i) = 0.
   enddo
!
!
!     compute the first guess of pbl height
!
!     compute the surface variables for pbl height estimation
!     under unstable conditions
!
!
!     enhance the pbl height by considering the thermal
!

!SW==
   do k = kts,kte
     do i = its,ite
          xkzm(i,k) = xkmv2d(i,k)
          xkzh(i,k) = xkhv2d(i,k)
          xkzq(i,k) = xkzh(i,k)
          xkzml(i,k) = xkzm(i,k)
          xkzhl(i,k) = xkzh(i,k)
     enddo
   enddo
!SW==
!
!-------------------------------------------------------------------------------------------
!     compute tridiagonal matrix elements for heat
!-------------------------------------------------------------------------------------------
!
   do k = kts,kte
     do i = its,ite
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.
     f1(i,1) = thx(i,1)-300.+hfx(i)/cont/del(i,1)*dt2
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzh(i,k)*rdz
       
       f1(i,k+1) = thx(i,k+1)-300.
       tem1   = dsig*xkzh(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
       exch_hx(i,k+1) = xkzh(i,k)
     enddo
   enddo
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
     enddo
   enddo
!
   call tridin_ysu(al,ad,cu,r1,au,f1,its,ite,kts,kte,1)
!
!     recover tendencies of heat
!
   do k = kte,kts,-1
     do i = its,ite
       ttend = (f1(i,k)-thx(i,k)+300.)*rdt*pi2d(i,k)
       ttnp(i,k) = ttnp(i,k)+ttend
       dtsfc(i) = dtsfc(i)+ttend*cont*del(i,k)/pi2d(i,k)
     enddo
   enddo
!-------------------------------------------------------------------------------------------
!     compute tridiagonal matrix elements for moisture, clouds, and gases
!-------------------------------------------------------------------------------------------
!
   do k = kts,kte
     do i = its,ite
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
     enddo
   enddo
!
   do ic = 1,ndiff
     do i = its,ite
       do k = kts,kte
         f3(i,k,ic) = 0.
       enddo
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.
     f3(i,1,1) = qx(i,1)+qfx(i)*g/del(i,1)*dt2
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do i = its,ite
         f3(i,1,ic) = qx(i,1+is)
       enddo
     enddo
   endif


   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzq(i,k)*rdz

       f3(i,k+1,1) = qx(i,k+1)

       tem1   = dsig*xkzq(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
!      exch_hx(i,k+1) = xkzh(i,k)
     enddo
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do k = kts,kte-1
         do i = its,ite
           f3(i,k+1,ic) = qx(i,k+1+is)
         enddo
       enddo
     enddo
   endif
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
     enddo
   enddo
!
   do ic = 1,ndiff
     do k = kts,kte
       do i = its,ite
         r3(i,k,ic) = f3(i,k,ic)
       enddo
     enddo
   enddo
!
!     solve tridiagonal problem for moisture, clouds, and gases
!
   call tridin_ysu(al,ad,cu,r3,au,f3,its,ite,kts,kte,ndiff)
!
!     recover tendencies of heat and moisture
!
   do k = kte,kts,-1
     do i = its,ite
       qtend = (f3(i,k,1)-qx(i,k))*rdt
       qtnp(i,k) = qtnp(i,k)+qtend
       dqsfc(i) = dqsfc(i)+qtend*conq*del(i,k)
     enddo
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do k = kte,kts,-1
         do i = its,ite
           qtend = (f3(i,k,ic)-qx(i,k+is))*rdt
           qtnp(i,k+is) = qtnp(i,k+is)+qtend
         enddo
       enddo
     enddo
   endif
!
!-------------------------------------------------------------------------------------------
!     compute tridiagonal matrix elements for momentum
!-------------------------------------------------------------------------------------------
!
   do i = its,ite
     do k = kts,kte
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
       f2(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
! paj: ctopo=1 if topo_wind=0 (default)
! mchen add this line to make sure NMM can still work with YSU PBL
       ad(i,1) = 1.+ust(i)**2/wspd1(i)*rhox(i)*g/del(i,1)*dt2                  &
        *(wspd1(i)/wspd(i))**2
     f1(i,1) = ux(i,1)
     f2(i,1) = vx(i,1)
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzm(i,k)*rdz
       f1(i,k+1) = ux(i,k+1)
       f2(i,k+1) = vx(i,k+1)

       tem1   = dsig*xkzm(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
     enddo
   enddo
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
       r2(i,k) = f2(i,k)
     enddo
   enddo
!
!     solve tridiagonal problem for momentum
!
   call tridi1n(al,ad,cu,r1,r2,au,f1,f2,its,ite,kts,kte,1)
!
!     recover tendencies of momentum
!
   do k = kte,kts,-1
     do i = its,ite
       utend = (f1(i,k)-ux(i,k))*rdt
       vtend = (f2(i,k)-vx(i,k))*rdt
       utnp(i,k) = utnp(i,k)+utend
       vtnp(i,k) = vtnp(i,k)+vtend
       dusfc(i) = dusfc(i) + utend*conwrc*del(i,k)
       dvsfc(i) = dvsfc(i) + vtend*conwrc*del(i,k)
     enddo
   enddo
!
!
!---- end of vertical diffusion
!
!
!
   end subroutine vdiff2d
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridi1n(cl,cm,cu,r1,r2,au,f1,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, kts,kte, nt
!
   real, dimension( its:ite, kts+1:kte+1 )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(in   )  ::                                                 cm, &
                                                                           r1
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(in   )  ::                                                 r2
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(inout)  ::                                                 au, &
                                                                           cu, &
                                                                           f1
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(inout)  ::                                                 f2
!
   real    :: fk
   integer :: i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
   l = ite
   n = kte
!
   do i = its,l
     fk = 1./cm(i,1)
     au(i,1) = fk*cu(i,1)
     f1(i,1) = fk*r1(i,1)
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
!
   do k = kts+1,n-1
     do i = its,l
       fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
       au(i,k) = fk*cu(i,k)
       f1(i,k) = fk*(r1(i,k)-cl(i,k)*f1(i,k-1))
     enddo
   enddo
!
   do it = 1,nt
     do k = kts+1,n-1
       do i = its,l
         fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
         f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
       enddo
     enddo
   enddo
!
   do i = its,l
     fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
     f1(i,n) = fk*(r1(i,n)-cl(i,n)*f1(i,n-1))
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
       f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
     enddo
   enddo
!
   do k = n-1,kts,-1
     do i = its,l
       f1(i,k) = f1(i,k)-au(i,k)*f1(i,k+1)
     enddo
   enddo
!
   do it = 1,nt
     do k = n-1,kts,-1
       do i = its,l
         f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
       enddo
     enddo
   enddo
!
   end subroutine tridi1n
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridin_ysu(cl,cm,cu,r2,au,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, kts,kte, nt
!
   real, dimension( its:ite, kts+1:kte+1 )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(in   )  ::                                                 cm
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(in   )  ::                                                 r2
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(inout)  ::                                                 au, &
                                                                           cu
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(inout)  ::                                                 f2
!
   real    :: fk
   integer :: i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
   l = ite
   n = kte
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       au(i,1) = fk*cu(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
!
   do it = 1,nt
     do k = kts+1,n-1
       do i = its,l
         fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
         au(i,k) = fk*cu(i,k)
         f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
       enddo
     enddo
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
       f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
     enddo
   enddo
!
   do it = 1,nt
     do k = n-1,kts,-1
       do i = its,l
         f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
       enddo
     enddo
   enddo
!
   end subroutine tridin_ysu
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!===============================================================================
!===============================================================================
!===============================================================================
END MODULE module_vdiff_driver
