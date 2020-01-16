MODULE module_force_scm

! Modified by Shuguang Wang for simulaitons with parameterized large-scale
! dynamcs: WTG (weak temperature gradient), spectral WTG, and DGW (damped gravity
! wave)

CONTAINS

   SUBROUTINE force_scm(itimestep, dt, scm_force, dx, num_force_layers       &
                             , scm_th_adv, scm_qv_adv                        &
                             , scm_ql_adv                                    &
                             , scm_wind_adv, scm_vert_adv,  scm_vert_adv_opt                &
                             , scm_th_relax_largescale                       &
                             , scm_th_relax_height                       &
                             , scm_qv_relax_largescale                       &
                             , scm_qv_relax_height                       &
                             , scm_th_hadv_largescale                       &
                             , scm_qv_hadv_largescale                       &
                             , scm_force_ql_largescale                       &
                             , scm_force_wind_largescale                     &
!SW==
           ,scm_dwdt         &
           ,scm_qv_hadv_tau         &
           ,scm_force_radiation         &
           ,scm_force_wave        &
           ,scm_force_wtg        &
                  ,qradiation_scm       &
                  ,th_scm_target        &
                  ,th_scm_prime       &
                  ,tao_time             &
                 ,pblh_opt                   &
                 ,wtg_opt                   &
                 ,wtg_cexp                  &
                 ,wtg_nmodes  &
                 ,scm_m_vector &
                 ,scm_zm &
                 ,scm_m_value &
                 ,scm_m_weight &
! for wave coupling
                 ,wavecoupling_opt                   &
                 ,wavecoupling_topbc, wavecoupling_dmpopt                   &
                 ,nl_bg, theta_bg, qv_bg, w_bg, z_bg              &
                 ,wavenumber, wavedamping            &
!SW==
                             , u_base, v_base, z_base                        &
                             , z_force, z_force_tend                         &
                             , w_subs, w_subs_tend                           &
                             ,th_largescale                                  &
                             ,th_largescale_tend                             &
                             ,qv_largescale                                  &
                             ,qv_largescale_tend                             &
                             ,qv_largescale_tend_deriv                             &
                             ,ql_largescale                                  &
                             ,ql_largescale_tend                             &
                             ,u_largescale                                   &
                             ,u_largescale_tend                              &
                             ,v_largescale                                   &
                             ,v_largescale_tend                              &
                             ,tau_largescale                                 &
                             ,tau_largescale_tend                            &
                             , z, z_at_w, th_phy, qv, ql, u_phy, v_phy                   &
                             ,  p_phy, pi_phy, p8w, t_phy,  rho             & !SW==
                             , thten, qvten, qlten, uten, vten               &
                             , scmthten, scmqvten                & !SW==
                             , mut, znu, znw, omega_wtg   &
                             , ids, ide, jds, jde, kds, kde                  &
                             , ims, ime, jms, jme, kms, kme                  &
                             , ips, ipe, jps, jpe, kps, kpe                  &
                             , i_start,i_end,j_start,j_end,kts,kte,num_tiles       &
                            )

     USE module_model_constants
     USE module_init_utilities, ONLY : interp_0
!SW==
    USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
!SW==
 
   IMPLICIT NONE


   INTEGER,    INTENT(IN   )                 :: itimestep
   INTEGER,    INTENT(IN   )                 :: num_force_layers, scm_force
   REAL,       INTENT(IN   )                 :: dt,dx
   LOGICAL,    INTENT(IN   )                 :: scm_th_adv, &
                                                scm_qv_adv, &
                                                scm_ql_adv, &
                                                scm_wind_adv, &
                                                scm_vert_adv, &
                                                scm_force_ql_largescale, &
                                                scm_th_relax_largescale,     &
                                                scm_qv_relax_largescale,     &
                                                scm_th_hadv_largescale ,     &
                                                scm_qv_hadv_largescale ,     &
                                                scm_force_radiation  ,       &
                                                scm_force_wave,       &
                                                scm_force_wtg,       &
                                                scm_force_wind_largescale, &
                                                scm_qv_hadv_tau,  &
                                                scm_dwdt

   REAL,       INTENT(IN   )                 ::  scm_th_relax_height     ,     &
                                                 scm_qv_relax_height     

   REAL,       DIMENSION( kms:kme ),                              &
                         INTENT(IN   )    ::               znu,   &
                                                           znw

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN   ) :: z, qv, ql, z_at_w

   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(IN   )    ::                         p_phy, &
                                                           pi_phy, &
                                                             p8w, &
                                                             rho, &
                                                           t_phy, &
                                                           u_phy, &
                                                           v_phy, &
                                                          th_phy

   REAL , DIMENSION( ims:ime , jms:jme ), INTENT(  IN) :: mut

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: thten, qvten
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: scmthten, scmqvten
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: qlten
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: uten, vten
   REAL, DIMENSION( kms:kme ), INTENT(INOUT)               :: u_base, v_base
   REAL, DIMENSION( kms:kme ), INTENT(INOUT)               :: z_base
   REAL, DIMENSION(num_force_layers), INTENT (INOUT)       :: z_force

   REAL, DIMENSION(num_force_layers), INTENT (IN) :: z_force_tend
   REAL, DIMENSION(num_force_layers), INTENT (IN) :: w_subs_tend

   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: w_subs
   REAL, DIMENSION( kms:kme ) :: w_subs_modellevels

   real, DIMENSION( kms:kme) :: omega_wtg

   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: th_largescale
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: th_largescale_tend
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: u_largescale
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: u_largescale_tend
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: v_largescale
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: v_largescale_tend
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: qv_largescale
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: qv_largescale_tend
   REAL, DIMENSION(kms:kme), INTENT (INOUT) :: qv_largescale_tend_deriv
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: ql_largescale
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: ql_largescale_tend
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: tau_largescale
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) :: tau_largescale_tend
!SW==
   REAL, DIMENSION(num_force_layers), INTENT (INOUT) ::   qradiation_scm       &
                                                         ,th_scm_target      
   REAL, DIMENSION( kms:kme ) :: th_scm_target_modellevels
   REAL, DIMENSION( kms:kme ),INTENT (INOUT) :: th_scm_prime
!SW==

!SW==
   INTEGER, INTENT(IN) ::        kts,kte,num_tiles

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                       &
     &           i_start,i_end,j_start,j_end
   REAL, INTENT(In) :: tao_time 
   REAL, DIMENSION( kms:kme ) :: tao_ptime

!SW==

   INTEGER,    INTENT(IN   )    ::     ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       ips,ipe, jps,jpe, kps,kpe

   
! Local
   REAL, DIMENSION(num_force_layers) :: dthdz_flayers
   INTEGER                      :: i,j,k
   LOGICAL                      :: debug = .false.
   REAL                         :: t_x, t_y, qv_x, qv_y, ql_x, ql_y
   REAL                         :: u_x, u_y, v_x, v_y
   REAL, DIMENSION(kms:kme)     :: th_adv_tend, qv_adv_tend, ql_adv_tend
   REAL, DIMENSION(kms:kme)     :: u_adv_tend, v_adv_tend
   REAL, DIMENSION(kms:kme)     :: dthdz, dudz, dvdz, dqvdz, dqldz
   REAL                         :: wtmp, wtmp_1, dthdz_1, dthdz_2, dqvdz_1,dqvdz_2
   REAL, DIMENSION(kms:kme)     :: w_dthdz, w_dudz, w_dvdz, w_dqvdz, w_dqldz, w_dthdz_vert, w_dqvdz_vert
   REAL, DIMENSION(kms:kme)     :: adv_timescale_x, adv_timescale_y
   CHARACTER*256                :: message
! Large-scale forcing WA 1/8/10
   REAL                         :: t_ls, qv_ls, ql_ls, qrad_ls, qv_ls_tmp2
   REAL                         :: u_ls, v_ls
   REAL, DIMENSION(kms:kme)     :: th_ls_tend, qv_ls_tend, ql_ls_tend
   REAL, DIMENSION(kms:kme)     :: u_ls_tend, v_ls_tend
   REAL, DIMENSION(kms:kme)     :: qrad_tend
   REAL, DIMENSION(kms:kme)     :: ls_timescale

   integer :: ij    
   real ::  no_points, uuu_sum, vvv_sum, zzz_sum, z8w_sum, dth_sum, qv_sum,  phalf_sum, pfull_sum, rho_sum, mu_sum, mu_avg
   real, DIMENSION( kms:kme) :: uuu_avg, vvv_avg, zzz_avg, z8w_avg,qv_avg, th_avg, rho_avg, phalf_avg, pfull_avg , th_v_avg
   real :: qtend1, qtend2,  omega_p, omega_n

   integer, intent(in) :: scm_vert_adv_opt

   integer :: scm_qv_hadv_opt
   real, dimension(kms:kme) :: div_wtg
!SW== wave coupling variables
    integer, intent(in) :: nl_bg
    integer, intent(in) :: pblh_opt 
    integer, intent(in) :: wtg_opt
    real, intent(in) :: wtg_cexp
    integer, intent(in) :: wavecoupling_opt
    integer, intent(in) :: wavecoupling_dmpopt ! option for damping: specify or constant or pressure depencance
    real, intent(in) :: wavecoupling_topbc  ! top of the wave coupling layer : meters
    real, intent(in) :: wavedamping, wavenumber
    real, DIMENSION( nl_bg), intent(in) :: theta_bg, qv_bg, z_bg, w_bg
                                                           
!SW for modal-wtg
    integer,intent(in) :: wtg_nmodes
    real,dimension(:),intent(in) :: scm_m_vector
    !real,dimension(1:num_force_layers,1:kme),intent(in) :: scm_m_vector
    real,dimension(1:kme),intent(in) :: scm_zm
    real,dimension(1:kme),intent(in) :: scm_m_value, scm_m_weight
    real,dimension(1:kme) :: th_over_dthdz_int,  th_over_dthdz, rho_int 
    real,dimension(wtg_nmodes, 1:kme) :: modes_modellevels
    real, dimension(1:wtg_nmodes) ::  am, am_s ! for Spectral-WTG  
    real,dimension(1:kme, 1:wtg_nmodes) :: mode_tmp, mode_tmp_deWeight
    integer :: m                                          
    integer :: nnlv 
    real :: pi    
!SW==
  logical ::  ltogadiag
  integer :: k_pbl = 0
  integer :: k_wtg_top = 0
  logical :: scm_force_ls_relax, scm_force_ls_tend
    integer :: ktmp
    integer,dimension(kms:kme) :: kneg_ind
    real, dimension(kms:kme) :: w_subs_useful, z_useful
    real,DIMENSION(num_force_layers) ::  dthdz_tmp
    real,DIMENSIOn(kms:kme) :: zf_tmp
    real :: ztmp
    CHARACTER*1024     :: outstring


   pi  =atan(1.0)*4.0

  modes_modellevels = 0.0

  dthdz_tmp = 0.0
  zf_tmp = 0.0
  ltogadiag = .false.
  scm_force_ls_relax = .true.
  scm_force_ls_tend = .true.

  scm_force_ls_relax = .false.
  scm_force_ls_tend = .true.

!SW==

   IF ( scm_force .EQ. 0 ) return
 
! NOTES
! z is kts:kte
! z_at_w is kms:kme



! 0 everything in case we don't set it later
    scmthten = 0.0
    scmqvten = 0.0
     th_adv_tend = 0.0
     qv_adv_tend = 0.0
     ql_adv_tend = 0.0
     u_adv_tend  = 0.0
     v_adv_tend  = 0.0
     th_ls_tend = 0.0
     qv_ls_tend = 0.0
     ql_ls_tend = 0.0
     u_ls_tend  = 0.0
     v_ls_tend  = 0.0
     w_dthdz     = 0.0
     w_dqvdz     = 0.0
     w_dqldz     = 0.0
     th_scm_target_modellevels = 0.0
     w_subs_modellevels = 0.0
     adv_timescale_x = 0.0
     adv_timescale_y = 0.0

    th_over_dthdz_int = 0.0
    th_over_dthdz = 0.0
    rho_int = 0.0

     

     CALL wrf_debug(100,'k z_base  u_base  v_base')
     do k = kms,kme-1
       WRITE( message, '(i4,3f12.4)' ) k,z_base(k),u_base(k),v_base(k)
       CALL wrf_debug ( 100, message )
    enddo

    DO ij = 1 , num_tiles;
     CALL wrf_debug(100,'its,ite,jts,jte=')
     WRITE( message, '(4i5)' )j_start(ij),j_end(ij),i_start(ij),i_end(ij)
     CALL wrf_debug ( 100, message )
    ENDDO
     CALL wrf_debug(100,'kts,kte,kms,kme=')
     WRITE( message, '(4i5)' )kts,kte,kms,kme
     CALL wrf_debug ( 100, message )


    IF(scm_force==1)then !calculate domain mean quantities:
      no_points = float((ide-ids)*(jde-jds))
      DO k=kts,min(kte+1,kde)
         uuu_sum = 0.0
         vvv_sum = 0.0
         zzz_sum = 0.0
         z8w_sum = 0.0
         qv_sum = 0.0
         dth_sum = 0.0
         phalf_sum = 0.0
         pfull_sum = 0.0
         rho_sum = 0.0

         DO ij = 1 , num_tiles;
            DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
               uuu_sum = uuu_sum + u_phy(i,k,j)
               vvv_sum = vvv_sum + v_phy(i,k,j)
               zzz_sum = zzz_sum + z(i,k,j)
               z8w_sum = z8w_sum + z_at_w(i,k,j)
               dth_sum = dth_sum + th_phy(i,k,j)
               qv_sum = qv_sum + qv(i,k,j)

               rho_sum = rho_sum + rho(i,k,j)
               phalf_sum = phalf_sum + p_phy(i,k,j)
               pfull_sum = pfull_sum + p8w(i,k,j)

               if(k==kts)  mu_sum = mu_sum + mut(i,j)

            ENDDO; ENDDO
         ENDDO

         uuu_avg(k) = wrf_dm_sum_real ( uuu_sum )
         uuu_avg(k) = uuu_avg(k) / no_points

         vvv_avg(k) = wrf_dm_sum_real ( vvv_sum )
         vvv_avg(k) = vvv_avg(k) / no_points

         zzz_avg(k) = wrf_dm_sum_real ( zzz_sum )
         zzz_avg(k) = zzz_avg(k) / no_points

         z8w_avg(k) = wrf_dm_sum_real ( z8w_sum )
         z8w_avg(k) = z8w_avg(k) / no_points

         th_avg(k) = wrf_dm_sum_real ( dth_sum )
         th_avg(k) = th_avg(k) / no_points

         qv_avg(k) = wrf_dm_sum_real ( qv_sum )
         qv_avg(k) = qv_avg(k) / no_points

         rho_avg(k) = wrf_dm_sum_real ( rho_sum )
         rho_avg(k) = rho_avg(k) / no_points

         pfull_avg(k) = wrf_dm_sum_real ( pfull_sum )
         pfull_avg(k) = pfull_avg(k) / no_points

         phalf_avg(k) = wrf_dm_sum_real ( phalf_sum )
         phalf_avg(k) = phalf_avg(k) / no_points

         if(rho_avg(k)>2.and.k<kde-1) then
            print*,'k,rho',k,rho_avg(k),th_avg(k)
            stop
         endif
         th_v_avg(k) = th_avg(k)*(1+0.608*qv_avg(k))

         mu_avg = wrf_dm_sum_real ( mu_sum )/no_points


      ENDDO
     CALL wrf_debug(100,'k zzz_avg  u_avg  v_avg')
     do k = kms,kme-1
       WRITE( message, '(i4,3f12.4)' ) k,zzz_avg(k),uuu_avg(k),vvv_avg(k)
       CALL wrf_debug ( 100, message )
    enddo
   ENDIF

  IF(scm_force_wave)THEN
    if(wavecoupling_opt==22)then  !calculate temperature anamalies using  toga obs as background state, and calculate omega using wave coupling approach
        do k=kts,kte
          th_scm_target_modellevels(k) = interp_0(th_scm_target,z_force,zzz_avg(k),num_force_layers)  !should be th_v_prime
        enddo
      call wave_coupling (                                            &
                  qvten, thten                                        &
                 , omega_wtg, znu, znw                                &
                 , wavecoupling_opt                                   &          ! on model variables
                 , wavecoupling_topbc, wavecoupling_dmpopt            &
                 , phalf_avg, pfull_avg, rho_avg, zzz_avg, th_avg, th_v_avg, qv_avg              &
                 , num_force_layers, th_largescale, qv_largescale*1e3, z_force,th_scm_target_modellevels                        &
                 ,th_scm_prime                                     &
                 , wavenumber, wavedamping                            &
                 , qv                                                 &
                 , dt                                                 &
                 , scm_dwdt                              &
                 , scm_vert_adv_opt  &
                 ,ids,ide, jds,jde, kds,kde                           &
                 ,ims,ime, jms,jme, kms,kme                           &
                 ,i_start,i_end, j_start,j_end, kts,kte, num_tiles                    )
    else
        CALL wrf_error_fatal('please set wavecoupling_opt==22');
    endif
  ENDIF !end of wavecoupling


 !$OMP PARALLEL DO   &
 !$OMP PRIVATE ( ij,i,j,k )
 DO ij = 1 , num_tiles;
   DO j=j_start(ij),j_end(ij);   DO i=i_start(ij),i_end(ij)
      do k = kts,kte
          dthdz(k) = 0.0
          dqvdz(k) = 0.0
          w_dthdz(k) = 0.0
          w_dthdz_vert(k) = 0.0
          w_dqvdz(k) = 0.0
          w_dqvdz_vert(k) = 0.0
          th_ls_tend(k) = 0.0
          qv_ls_tend(k) = 0.0
          u_ls_tend(k) = 0.0
          v_ls_tend(k) = 0.0
          qrad_tend(k) = 0.0
      enddo


       do k = kts,kte
          ls_timescale(k) = interp_0(tau_largescale,z_force,z(i,k,j),num_force_layers)
       enddo

    if ( scm_th_relax_largescale ) then
       !if (  scm_force_ls_relax) then
        do k = kts,kte-1
             ztmp = zzz_avg(k)-scm_th_relax_height
             if( ztmp >= 1000.0 .and.zzz_avg(k)<=z_force(num_force_layers))then
                t_ls = interp_0(th_largescale,z_force,zzz_avg(k),num_force_layers)
                if(t_ls >= 200.0 .and. t_ls < 800.0)&
                th_ls_tend(k) = th_ls_tend(k) + (t_ls-th_avg(k))/3600.0/1.0
                !th_ls_tend(k) = th_ls_tend(k) + (t_ls-th_avg(k))/3600.0/6.0
             elseif( ztmp>=0 .and. (ztmp<1000.0 ).and.zzz_avg(k)<=z_force(num_force_layers))then
                t_ls = interp_0(th_largescale,z_force,zzz_avg(k),num_force_layers)
                if(t_ls >= 200.0 .and. t_ls < 800.0)&
                th_ls_tend(k) = th_ls_tend(k) + (t_ls-th_avg(k))/3600.0*(1.0/1.0*ztmp/1000.0)
                !th_ls_tend(k) = th_ls_tend(k) + (t_ls-th_avg(k))/3600.0*(1.0/6.0*ztmp/1000.0)
             endif
        enddo
    endif

    if( scm_th_hadv_largescale) then ! WA if upstream is empty, use tendency only not value+tend
          do k = kts,kte-1
             if(z(i,k,j)<= z_force(num_force_layers))then
                t_ls = interp_0(th_largescale_tend,z_force,z(i,k,j),num_force_layers) !SW==
                th_ls_tend(k) =    th_ls_tend(k) +  t_ls
             endif
          enddo
    endif

    
    if ( scm_qv_relax_largescale  ) then
    !if ( scm_qv_relax_largescale .and. (.not. scm_qv_hadv_largescale) ) then
          do k = kms,kme-1
             ztmp = zzz_avg(k)-scm_qv_relax_height
             if( ztmp >= 1000.0 .and.zzz_avg(k)<=z_force(num_force_layers))then
                qv_ls = interp_0(qv_largescale,z_force,zzz_avg(k),num_force_layers)
                qv_ls_tend(k) = qv_ls_tend(k) + (qv_ls-qv_avg(k))/3600.0/3.0
             elseif( ztmp>=0 .and. (ztmp<1000.0 ).and.zzz_avg(k)<=z_force(num_force_layers))then
                qv_ls = interp_0(qv_largescale,z_force,zzz_avg(k),num_force_layers)
                qv_ls_tend(k) = qv_ls_tend(k) + (qv_ls-qv_avg(k))/3600.0*(1.0/3.0*ztmp/1000.0)
             endif
          enddo
     endif

     if( scm_qv_hadv_largescale  ) then ! WA if upstream is empty, use tendency only not value+tend
          do k = kts,kte
             if(z(i,k,j)< z_force(num_force_layers))then
               qv_ls = interp_0(qv_largescale_tend,z_force,z(i,k,j),num_force_layers)
               qv_ls_tend(k) =  qv_ls_tend(k) + qv_ls
             endif
          enddo

    endif



    if (scm_force_radiation ) then
          do k = kts,kte
             qrad_ls = interp_0(qradiation_scm,z_force,z(i,k,j),num_force_layers) !SW==
             qrad_tend(k) =  qrad_ls
          enddo
    endif


    if ( scm_force_wind_largescale ) then
       if ( u_largescale(1) > -900.0) then
          do k = kts,kte
           if(zzz_avg(k) < z_force(num_force_layers)  )then
             u_ls = interp_0(u_largescale,z_force,zzz_avg(k),num_force_layers)
             v_ls = interp_0(v_largescale,z_force,zzz_avg(k),num_force_layers)
             u_ls_tend(k) = (u_ls-uuu_avg(k))/ls_timescale(k)
             v_ls_tend(k) = (v_ls-vvv_avg(k))/ls_timescale(k)
           endif
          enddo
       endif
    endif


    CALL wrf_debug ( 200, 'j, k, th_adv_ten, qv_adv_ten, ql_adv_ten, u_adv_ten, v_adv_ten')

  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  IF ( scm_force_wtg )then
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
    if( wtg_opt==1 .or. wtg_opt==2 .or. wtg_opt==3) then  ! interpolate to force_layers first, calculate w_subs, and finally advect moisture 
      do k = kts,kte
         if(k>kts)then
            dqvdz(k) = ((qv_avg(k)-qv_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1)) + (qv_avg(k+1)-qv_avg(k))/(zzz_avg(k+1)-zzz_avg(k)))*0.5
         elseif(k==1)then
            dqvdz(1) = ( (qv_avg(2)-qv_avg(1))/(zzz_avg(2)-zzz_avg(1)))
         endif
         w_subs_modellevels(k) = 0.0
         omega_wtg(k) = 0.0
      enddo

      if(wtg_opt ==1)then ! for  temperature
          do k = kts,kte
             if(k>1)then
                  dthdz(k) = ((th_avg(k)-th_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1))+ (th_avg(k+1)-th_avg(k))/(zzz_avg(k+1)-zzz_avg(k))) *0.5
                  dthdz(k) = max(dthdz(k), 0.001)
             elseif(k==1)then
                  dthdz(1) = ( (th_avg(2)-th_avg(1))/(zzz_avg(2)-zzz_avg(1))); dthdz(1) = max(dthdz(1), 0.001)
             endif
             th_scm_target_modellevels(k) = th_avg(k)-interp_0(th_largescale,z_force,zzz_avg(k),num_force_layers)
             th_scm_prime(k)  = th_scm_target_modellevels(k) 
          enddo
      elseif(wtg_opt ==2)then ! for virtual temperature
          do k = kts,kte
              if(k>=kts+1)then
                  dthdz(k) = ((th_avg(k)*(1+0.608*qv_avg(k))-th_avg(k-1)*(1+0.608*qv_avg(k-1)))/(zzz_avg(k)-zzz_avg(k-1))+ &
                             (th_avg(k+1)*(1+0.608*qv_avg(k+1))-th_avg(k)*(1+0.608*qv_avg(k)))/(zzz_avg(k+1)-zzz_avg(k))) *0.5
                  dthdz(k) = max(dthdz(k), 0.001)
              elseif(k==1)then
                  dthdz(1) = ( (th_avg(2)*(1+0.608*qv_avg(2))-th_avg(1)*(1+0.608*qv_avg(1)))/(zzz_avg(2)-zzz_avg(1))); dthdz(1) = max(dthdz(1), 0.001)
              endif
              th_scm_target_modellevels(k) = th_avg(k)*(1+0.608*qv_avg(k))- &
                        interp_0(th_largescale,z_force,zzz_avg(k),num_force_layers)*(1+0.608*interp_0(qv_largescale,z_force,zzz_avg(k),num_force_layers))
              th_scm_prime(k)  = th_scm_target_modellevels(k) 
          enddo
      endif

      do k = kte, kts, -1
         if(zzz_avg(k) < wavecoupling_topbc  )then
             k_wtg_top = k
             exit
         endif
      enddo

      do k=k_wtg_top,1,-1
        if(zzz_avg(k).gt.real(pblh_opt))then
             w_subs_modellevels(k) = th_scm_target_modellevels(k)/tao_time/dthdz(k)
             k_pbl = k
        endif
      enddo
      do k=1,k_pbl-1
         w_subs_modellevels(k) = zzz_avg(k)*w_subs_modellevels(k_pbl)/zzz_avg(k_pbl)
      enddo
    
      do k=1,k_wtg_top
         omega_wtg(k) = -w_subs_modellevels(k)*rho_avg(k)*9.81 ! use pressure velocity
      enddo



      if(scm_vert_adv_opt == 800)then ! use W at the middle point

         call vertical_advection_omega_800(                            &
                 w_dqvdz, w_dthdz, omega_wtg, znu, znw                              &
                 ,phalf_avg, pfull_avg, rho_avg, zzz_avg, th_avg, qv_avg    &
                 ,kms,kme, kts,kte,kds,kde                   )
         w_dqvdz = -w_dqvdz
         w_dthdz = -w_dthdz

      endif
  ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------
  elseif( wtg_opt >= 200 .and. wtg_opt < 300)then ! vertical mode based wtg, wtg based on mode decomposition, correct mode projection
      do k = kts+1,kte
         dthdz(k) = ((th_avg(k)-th_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1))+ (th_avg(k+1)-th_avg(k))/(zzz_avg(k+1)-zzz_avg(k))) *0.5
         dthdz(k) = max(dthdz(k), 0.001)
         dqvdz(k) = ((qv_avg(k)-qv_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1)) + (qv_avg(k+1)-qv_avg(k))/(zzz_avg(k+1)-zzz_avg(k)))*0.5
      enddo
      dthdz(1) = ( (th_avg(2)-th_avg(1))/(zzz_avg(2)-zzz_avg(1))); dthdz(1) = max(dthdz(1), 0.001) 
      dqvdz(1) = ( (qv_avg(2)-qv_avg(1))/(zzz_avg(2)-zzz_avg(1))) 


      w_subs_modellevels = 0.0
      do k = kts,kte
           if(wtg_opt ==200)then
              th_scm_target_modellevels(k) = th_avg(k)-interp_0(th_largescale,z_force,zzz_avg(k),num_force_layers)
           elseif(wtg_opt ==210)then ! for virtual temperature
              if(k>=kts+1)then 
                  dthdz(k) = ((th_avg(k)*(1+0.608*qv_avg(k))-th_avg(k-1)*(1+0.608*qv_avg(k-1)))/(zzz_avg(k)-zzz_avg(k-1))+ &
                             (th_avg(k+1)*(1+0.608*qv_avg(k+1))-th_avg(k)*(1+0.608*qv_avg(k)))/(zzz_avg(k+1)-zzz_avg(k))) *0.5
                  dthdz(k) = max(dthdz(k), 0.001)
              elseif(k==1)then
                  dthdz(1) = ( (th_avg(2)*(1+0.608*qv_avg(2))-th_avg(1)*(1+0.608*qv_avg(1)))/(zzz_avg(2)-zzz_avg(1))); dthdz(1) = max(dthdz(1), 0.001)
              endif
              th_scm_target_modellevels(k) = th_avg(k)*(1+0.608*qv_avg(k))- &
                        interp_0(th_largescale,z_force,zzz_avg(k),num_force_layers)*(1+0.608*interp_0(qv_largescale,z_force,zzz_avg(k),num_force_layers))

           elseif(wtg_opt ==220)then
              th_scm_target_modellevels(k) = th_avg(k)-interp_0(th_largescale,z_force,zzz_avg(k),num_force_layers)
              dthdz(k) = 1.0 ! assume 4 K/km 
           endif
           th_scm_prime(k)  = th_scm_target_modellevels(k) 
         th_over_dthdz(k) = th_scm_target_modellevels(k)/dthdz(k)

         do m=1,wtg_nmodes
            mode_tmp(k,m) = scm_m_vector(m*kme+(k))
            mode_tmp_deWeight(k,m) = mode_tmp(k,m)/scm_m_weight(k) 
         enddo
      enddo
      do k = kts,kte
            if( k>10 .and. scm_zm(k) == 0 )then
               nnlv = k-1
                if(i==1.and.j==1.and. mod(itimestep,100)==1) print*,'nnlv=', nnlv
               exit
            endif
      !nnlv = 35
      enddo

      do k = kts,kte
         if(zzz_avg(k)<=wavecoupling_topbc ) then
         do m=1,wtg_nmodes
            modes_modellevels(m,k) = interp_0(mode_tmp_deWeight(1:nnlv,m),scm_zm(1:nnlv), zzz_avg(k), nnlv)
         enddo
         endif
         th_over_dthdz_int(k) = interp_0(th_over_dthdz(1:kte-1), zzz_avg(1:kte-1), scm_zm(k), kte-1) 
         rho_int(k) = interp_0(rho_avg(1:kte-1), zzz_avg(1:kte-1), scm_zm(k), kte-1) 
      enddo

      am = 0.0; am_s = 0.0;
      do m=1,wtg_nmodes ! first 10 modes
      do k = 1, nnlv
            am(m) = am(m) +  rho_int(k)*th_over_dthdz_int(k) *mode_tmp(k,m)*scm_m_weight(k) ! scm_m_weight: N, mode_tmp: N*Vn 
            am_s(m) = am_s(m) + mode_tmp(k,m)*mode_tmp(k,m) 
      enddo
      enddo
      if(i==1.and.j==1 .and. mod(itimestep,100)==1)  then
         print*,'am=', am(1:wtg_nmodes) 
         print*,'am_s=', am_s(1:wtg_nmodes) 
      endif

      do k=kts, kte
        do m= 1,wtg_nmodes !first 10 modes
         if( zzz_avg(k) < wavecoupling_topbc)then
             w_subs_modellevels(k) = w_subs_modellevels(k) + am(m)/am_s(m)* modes_modellevels(m,k) /tao_time*scm_m_value(m)**wtg_cexp
         endif
        enddo
      enddo
   
 
      do k=kts,kte
         omega_wtg(k) = -w_subs_modellevels(k)*9.81 ! use pressure velocity
      enddo

      if(scm_vert_adv_opt == 1)then 
         call vertical_advection_omega(                            &
                 w_dqvdz, w_dthdz, omega_wtg, znu, znw                              &
                 ,phalf_avg, pfull_avg, rho_avg, zzz_avg, th_avg, qv_avg    &
                 ,kms,kme, kts,kte,kds,kde                   )
      elseif(scm_vert_adv_opt == 800)then 
         call vertical_advection_omega_800(                            &
                 w_dqvdz, w_dthdz, omega_wtg, znu, znw                              &
                 ,phalf_avg, pfull_avg, rho_avg, zzz_avg, th_avg, qv_avg    &
                 ,kms,kme, kts,kte,kds,kde                   )
      endif
         w_dqvdz = -w_dqvdz
         w_dthdz = -w_dthdz


  ! -----------------------------------------------------------------------------------------------------------------------------------------------------

  elseif( wtg_opt==10) then  ! spectral WTG, interpolate to force_layers first, calculate w_subs, and finally advect moisture 
      do k = kts+1,kte
         dthdz(k) = ((th_avg(k)-th_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1))+ (th_avg(k+1)-th_avg(k))/(zzz_avg(k+1)-zzz_avg(k))) *0.5
         dthdz(k) = max(dthdz(k), 0.001)
         dqvdz(k) = ((qv_avg(k)-qv_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1)) + (qv_avg(k+1)-qv_avg(k))/(zzz_avg(k+1)-zzz_avg(k)))*0.5
      enddo
      dthdz(1) = ( (th_avg(2)-th_avg(1))/(zzz_avg(2)-zzz_avg(1))); dthdz(1) = max(dthdz(1), 0.001) 
      dqvdz(1) = ( (qv_avg(2)-qv_avg(1))/(zzz_avg(2)-zzz_avg(1))) 


      w_subs_modellevels = 0.0
      do k = kts,kte
            th_scm_target_modellevels(k) = th_avg(k)-interp_0(th_largescale,z_force,zzz_avg(k),num_force_layers)
            th_scm_prime(k)  = th_scm_target_modellevels(k) 
      enddo

      am = 0.0; am_s = 0.0;
      do m=1,10 ! first 10 modes
      do k = kts, kte
         if( zzz_avg(k) < wavecoupling_topbc)then
            am(m) = am(m) + ( th_scm_target_modellevels(k)/dthdz(k)*sin(m*pi*zzz_avg(k)/wavecoupling_topbc) + &
                              th_scm_target_modellevels(k+1)/dthdz(k+1)*sin(m*pi*zzz_avg(k+1)/wavecoupling_topbc) &
                            )*0.5* (zzz_avg(k+1)-zzz_avg(k)) 
            am_s(m) = am_s(m) + ( sin(m*pi*zzz_avg(k)/wavecoupling_topbc)**2 + sin(m*pi*zzz_avg(k+1)/wavecoupling_topbc)**2 ) * 0.5* (zzz_avg(k+1)-zzz_avg(k))
         endif
      enddo
      enddo
      
      do k=kts, kte
        do m= 1,10 !first 10 modes
         if( zzz_avg(k) < wavecoupling_topbc)then
             w_subs_modellevels(k) = w_subs_modellevels(k) + am(m)/am_s(m)*sin(m*pi*zzz_avg(k)/wavecoupling_topbc)/tao_time/real(m)
         endif
        enddo
      enddo
   
 
      do k=kts,kte
         omega_wtg(k) = -w_subs_modellevels(k)*rho_avg(k)*9.81 ! use pressure velocity
      enddo

      if(1==1)then  
         call vertical_advection_omega(                            &
                 w_dqvdz, w_dthdz, omega_wtg, znu, znw                              &
                 ,phalf_avg, pfull_avg, rho_avg, zzz_avg, th_avg, qv_avg    &
                 ,kms,kme, kts,kte,kds,kde                   )
         w_dqvdz = -w_dqvdz
         w_dthdz = -w_dthdz

      else
        do k = kts,kte 
          if(zzz_avg(k) < wavecoupling_topbc )then
              !w_dthdz(k) = -w_subs_modellevels(k)*dthdz(k)
              w_dthdz(k) = -th_scm_target_modellevels(k)/tao_time
              w_dqvdz(k) = -w_subs_modellevels(k)*dqvdz(k)
          endif
        enddo
      endif

  ! -----------------------------------------------------------------------------------------------------------------------------------------------------

  endif 
  ENDIF ! END of WTG
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------

! Impose W
  IF ( scm_vert_adv ) then
     if(scm_vert_adv_opt==1 .or. scm_vert_adv_opt==11)then ! perform vertical advection on forcing levels 
        do k = kts,kte
           if(k==kts)then
             dthdz(k) =  (th_avg(k+1)-th_avg(k))/(zzz_avg(k+1)-zzz_avg(k))  
             dqvdz(k) =  ((qv_avg(k+1)-qv_avg(k))/(zzz_avg(k+1)-zzz_avg(k))) 
           else
             dthdz(k) = 0.5*((th_avg(k)-th_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1)) + (th_avg(k+1)-th_avg(k))/(zzz_avg(k+1)-zzz_avg(k))  )
             dqvdz(k) = 0.5*((qv_avg(k)-qv_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1)) + ((qv_avg(k+1)-qv_avg(k))/(zzz_avg(k+1)-zzz_avg(k))) )
           endif
        enddo
        do k = kts,kte
          if( scm_vert_adv_opt == 11 )then  !turn on vertical advection by the background time-independant large scale vertical velocity
            if(zzz_avg(k) < z_bg(nl_bg)) then
               wtmp = interp_0(w_bg,z_bg,z8w_avg(k),nl_bg)
            else
               wtmp=0.0
            endif
          else
            if(k==kts)then
              wtmp = 0.0 
            else
              wtmp = interp_0(w_subs,z_force,z8w_avg(k),num_force_layers)
            endif
            wtmp_1 = interp_0(w_subs,z_force,z8w_avg(k+1),num_force_layers)
          endif

          if(k==kts)then
             dthdz_1 = 0.0
             dqvdz_1 = 0.0
          else          
             dthdz_1=(th_avg(k)-th_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1))
             dqvdz_1=(qv_avg(k)-qv_avg(k-1))/(zzz_avg(k)-zzz_avg(k-1))
          endif

          dthdz_2=(th_avg(k+1)-th_avg(k))/(zzz_avg(k+1)-zzz_avg(k))
          dqvdz_2=(qv_avg(k+1)-qv_avg(k))/(zzz_avg(k+1)-zzz_avg(k))

          w_dthdz_vert(k) = -max(wtmp,0.0)*dthdz_1 - min(wtmp_1,0.0)*dthdz_2
          w_dqvdz_vert(k) = -max(wtmp,0.0)*dqvdz_1 - min(wtmp_1,0.0)*dqvdz_2
        enddo

     elseif(scm_vert_adv_opt==2)then ! perform vertical advection on eta model levels
        do k = kts,kte
              if(z_force(num_force_layers) > zzz_avg(k))then
                   omega_wtg(k) = interp_0 (w_subs, z_force,zzz_avg(k),num_force_layers)
                   omega_wtg(k) = -omega_wtg(k)*rho_avg(k)*9.81/mu_avg
              else
                   omega_wtg(k)  = 0.0
              endif
        enddo
        do k = kts,kte
                 if(k>1.and.k<kde-1)then
                        omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                        omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0)
                        qtend2 = omega_p*(qv_avg(k+1) - qv_avg(k))/(znu(k+1)-znu(k))
                        qtend1 = omega_n*(qv_avg(k) - qv_avg(k-1)) /(znu(k)-znu(k-1))
                        w_dqvdz_vert(k) = -(qtend1+qtend2)
                 elseif(k==kde-1)then
                        omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0)
                        qtend1 = omega_n*(qv_avg(k) - qv_avg(k-1)) /(znu(k)-znu(k-1))
                        w_dqvdz_vert(k) = -qtend1
                 elseif(k==1)then
                        omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                        qtend2 = omega_p*(qv_avg(k+1) - qv_avg(k))/(znu(k+1)-znu(k))
                        w_dqvdz_vert(k) = -qtend2
                 endif
                 if(k>1.and.k<kde-1)then
                        omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                        omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0)
                        qtend2 = omega_p*(th_avg(k+1) - th_avg(k))/(znu(k+1)-znu(k))
                        qtend1 = omega_n*(th_avg(k) - th_avg(k-1)) /(znu(k)-znu(k-1))
                        w_dthdz_vert(k) = -(qtend1+qtend2)
                     elseif(k==kde-1)then
                        omega_n = min((omega_wtg(k-1)+omega_wtg(k))*0.5, 0.0)
                        qtend1 = omega_n*(th_avg(k) - th_avg(k-1)) /(znu(k)-znu(k-1))
                        w_dthdz_vert(k) = -qtend1
                     elseif(k==1)then
                        omega_p = max((omega_wtg(k+1)+omega_wtg(k))*0.5, 0.0)
                        qtend2 = omega_p*(th_avg(k+1) - th_avg(k))/(znu(k+1)-znu(k))
                        w_dthdz_vert(k) = -qtend2
                 endif
      enddo

     endif
    ENDIF ! END of Impose-W
 

    do k = kts,kte-1
       thten(i,k,j) = thten(i,k,j) + th_adv_tend(k) +              &
                      w_dthdz(k)               &                               !WTG or wavecoupling
                      + w_dthdz_vert(k)        &
                      + th_ls_tend(k) + qrad_tend(k)
       scmthten(i,k,j) =  w_dthdz(k) + w_dthdz_vert(k) + th_ls_tend(k)


       qvten(i,k,j) = qvten(i,k,j) + qv_adv_tend(k) +              &
                         w_dqvdz(k)               &
                      + w_dqvdz_vert(k)        &
                      + qv_ls_tend(k)
       scmqvten(i,k,j) = w_dqvdz(k) + w_dqvdz_vert(k) + qv_ls_tend(k)
       uten(i,k,j)  = uten(i,k,j) + u_adv_tend(k) +                &
                      + u_ls_tend(k)
       vten(i,k,j)  = vten(i,k,j) + v_adv_tend(k) +                &
                      + v_ls_tend(k)
    enddo


    ENDDO; ENDDO ! end of i%j
 ENDDO !end of num_tiles
 !$OMP END PARALLEL DO

      if(ltogadiag) then
        print*,'k_pbl=',k_pbl
        print*,'w_subs=',w_subs
      endif

    RETURN

   END SUBROUTINE force_scm

   SUBROUTINE wave_coupling (                                      &
                  qvten, thten                                     &
                 ,omega_wtg, znu, znw                              &
                 ,wavecoupling_opt                                 &  
                 ,wavecoupling_topbc, wavecoupling_dmpopt          & 
                 , phalf_avg, pfull_avg, rho_avg, z_avg, th_avg, th_v_avg, qv_avg                    & 
                 ,nl_bg, theta_bg, qv_bg, z_bg, th_scm_target_modellevels                            & 
                 ,th_scm_prime                                     &
                 ,wavenumber, wavedamping                          &  
                 ,qv_curr                                          &
                 , dt                                              &
                 , scm_dwdt                              &
                 ,scm_vert_adv_opt  &
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,i_start,i_end, j_start,j_end, kts,kte, num_tiles                    )
 
    USE module_model_constants

    USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
   IMPLICIT NONE

   INTEGER,    INTENT(IN   )    ::     ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       kts,kte
   REAL,       DIMENSION( kms:kme ),                              &
                          INTENT(IN   )    ::               znu,   &
                                                           znw
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         INTENT(IN) ::        qv_curr 
                      ! optional moisture tracers
                      ! 2 time levels; if only one then use CURR


   INTEGER :: i,J,K,NK,jj,ij,its,ite,jts,jte

   INTEGER, INTENT(IN) ::    num_tiles

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                       &
     &           i_start,i_end,j_start,j_end

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: thten, qvten

   real, DIMENSION( kms:kme), intent(in) :: z_avg, phalf_avg, pfull_avg, rho_avg, th_v_avg, th_avg, qv_avg
   REAL :: dt
    real ::  no_points, th_temp, uuu_sum
    real ::  dth_sum, rho_sum, qv_sum, phalf_sum, pfull_sum, z_sum, mu_sum,  mu_avg
    integer :: irnc, jrnc
    integer :: k_abovepbl_first
    real, DIMENSION( ims:ime , jms:jme ) :: thetaint_tmp,qvint_tmp, qvhadvint_tmp
    real, DIMENSION( kms:kme) :: omega_wtg
    real, DIMENSION( kms:kme) :: dthdeta, dqvdeta, qtend, qtend_qvhadv, qtend_avg, ttend
    real ::  qtend1, qtend2, domgdeta
    real, dimension(kms:kme) :: theta_v_prime

    real :: omega_p, omega_n
    real :: tao_time
    CHARACTER*1024     :: outstring
    CHARACTER*1024     :: wrf_err_message

    ! reference potential temperature profiles 
    integer, intent(in) :: nl_bg
    integer, intent(in) :: wavecoupling_opt
    integer, intent(in) :: wavecoupling_dmpopt ! option for damping: specify or constant or pressure depencance
    real, intent(in) :: wavecoupling_topbc  ! top of the wave coupling layer : meters
    real, intent(in) :: wavedamping, wavenumber
    real, DIMENSION( nl_bg), intent(in) :: theta_bg, qv_bg, z_bg
    real, DIMENSION( kms:kme) :: th_scm_target_modellevels, th_scm_prime


    ! local variables for wave coupling
    real, DIMENSION( kms:kme) :: dp, rhs, aa, bb, cc
    real, DIMENSION( 1:kme) :: theta_bg_int, qv_bg_int
    real, DIMENSION( 1:kme) :: wavedamping_z, wavedamping_dz ! wave damping
    integer :: np  !number of layers for wave coupling
    integer :: k2

   real :: domgdt_tmp

    real :: cpm, evapor,rho_d,coef_tmp
   CHARACTER*512                :: message

    logical scm_dwdt
    integer, intent(in) :: scm_vert_adv_opt

 
    !scm_dwdt = .True.
    qvint_tmp = 0.0
    qvhadvint_tmp = 0.0
    thetaint_tmp = 0.0


    dp = 0.0;
    rhs = 0.0;
    aa = 0.0;
    bb = 0.0;
    cc = 0.0;


      
     k=1
     do k=1,kde-1
         if(z_avg(k)>wavecoupling_topbc)then 
             np=k-1
             exit 
         endif
     enddo
     do k=np+1, kde
        omega_wtg(k) = 0.0 ! above troposphere  omega_wtg is zero
     enddo

     ! ----------------------------------------------------------------------------
     !     compute grid spacing between midpoint levels
      do k=2,np
         dp(k)=phalf_avg(k)-phalf_avg(k-1)
      enddo
      dp(1)=2*(phalf_avg(1)-pfull_avg(1))
      dp(np+1)=2*(pfull_avg(np+1)-phalf_avg(np))
     
      !wavedamping = 1/86400.0 ! damping time scale: 1 day as in Kuang 2008 and 2010 
      !wavenumber = 2*pi/10000.0e3 ! wavenumber

      ! Interpolate theta_bg to coordinate levels:
      DO k = kts, kte
        k2 = nl_bg 
        DO WHILE( z_bg(k2) .gt. z_avg(k)  .and.  k2 .gt. 1 )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt. nl_bg)then
          theta_bg_int(k) = theta_bg(k2) + ( theta_bg(k2) - theta_bg(k2-1) )   &
                              * (     z_avg(k) - z_bg(k2)   )   &
                              / ( z_bg(k2) - z_bg(k2-1) )
        else
          theta_bg_int(k) = theta_bg(k2) + ( theta_bg(k2+1) - theta_bg(k2) )   &
                              * (       z_avg(k) - z_bg(k2) )   &
                              / ( z_bg(k2+1) - z_bg(k2) )
        endif
      ENDDO


      DO k = kts, kte
        k2  = nl_bg
        DO WHILE( z_bg(k2) .gt. z_avg(k)  .and.  k2 .gt. 1 )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.nl_bg)then
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

      DO k = kts, kte 
      if(wavecoupling_dmpopt==1)then
           wavedamping_z (k) = wavedamping
           wavedamping_dz (k) = 0
      elseif(wavecoupling_dmpopt==2)then
           wavedamping_z (k) = wavedamping * phalf_avg(k)/1.0e5
           wavedamping_dz (k) = wavedamping /1.0e5
      else
           WRITE(wrf_err_message,*) ' wrong option for wavecoupling_dmpopt: 1 or 2'
           CALL wrf_error_fatal(TRIM(wrf_err_message))
      endif
      ENDDO

      !     Gauss Elimination
      rhs=0.
      do k=1,np
         if(wavecoupling_opt==22)then
            theta_v_prime(k) = th_v_avg(k)-theta_bg_int(k)*(1+0.608*qv_bg_int(k))
            th_scm_prime(k) =  theta_v_prime(k) 
         endif
         rhs(k)=wavenumber*wavenumber/wavedamping_z(k)*R_d*theta_v_prime(k)*(phalf_avg(k)/1.e5)**(R_d/Cp)/phalf_avg(k)*dp(k)*dp(k+1)*0.5
         !print*,'k,prime', th_v_avg(k), theta_bg_int(k)*(1+0.608*qv_bg_int(k)), theta_v_prime
      enddo
     

    if(wavecoupling_dmpopt==1)then
       !     set up the tridiagonal matrix
       do k=1,np
         aa(k)=dp(k+1)/(dp(k)+dp(k+1))
         bb(k)=-1.;
         cc(k)=dp(k)/(dp(k)+dp(k+1))
       enddo
      !     symmetric lower BC
      aa(1)=0.;
      bb(1)=-(2*dp(2)+dp(1))/(dp(1)+dp(2))
      cc(1) = dp(k)/(dp(k)+dp(k+1))
      !     symmetric upper BC
      aa(np)=dp(np+1)/(dp(np)+dp(np+1))
      bb(np) = -(2*dp(np)+dp(np+1))/(dp(np)+dp(np+1))
      cc(np) = 0.

    elseif(wavecoupling_dmpopt==2)then
       do k=1,np
         aa(k)=dp(k+1)/(dp(k)+dp(k+1))  -dp(k+1)*0.5*wavedamping_dz(k)/wavedamping_z(k)
         bb(k)=-1.  +(dp(k+1)-dp(k))*0.5*wavedamping_dz(k)/wavedamping_z(k)
         cc(k)=dp(k)/(dp(k)+dp(k+1))  +dp(k)*0.5*wavedamping_dz(k)/wavedamping_z(k)
       enddo
      !     symmetric lower BC
       aa(1)=0.;
       bb(1)=-(2*dp(2)+dp(1))/(dp(1)+dp(2))  -dp(1)*0.5*wavedamping_dz(1)/wavedamping_z(1)
       cc(1) = dp(k)/(dp(k)+dp(k+1))  +dp(1)*0.5*wavedamping_dz(2) /wavedamping_z(1)
       !     symmetric upper BC
       aa(np)=dp(np+1)/(dp(np)+dp(np+1))  -dp(np+1)*0.5*wavedamping_dz(np)/wavedamping_z(np)
       bb(np) = -(2*dp(np)+dp(np+1))/(dp(np)+dp(np+1))  +dp(np+1)*0.5*wavedamping_dz(np)/wavedamping_z(np)
       !bb(np) = -1 
       cc(np) = 0.
    else
       stop 'wrong wavecoupling_dmpopt'
    endif

      
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
     
      if(scm_dwdt)then
         do k=1,np
           domgdt_tmp = wavedamping_z(k)*(-omega_wtg(k)+rhs(k)) 
           omega_wtg(k) = omega_wtg(k) + domgdt_tmp*dt  
           omega_wtg(k) = omega_wtg(k) + 4800.0*domgdt_tmp  
           !omega_wtg(k) = omega_wtg(k) + wavedamping_z(k)*(-omega_wtg(k)+rhs(k))*dt ! old version 
         enddo
         

      else 
        do k=1,np
           omega_wtg(k)=rhs(k)  ! omega, this is pressure velocity!!
        enddo
      endif

     ! ----------------------------------------------------------------------------
      do k=kts,kde ! on model state
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
                 qvten(i,k,j) = qvten (i,k,j) - qtend(k)
                 thten(i,k,j) = thten(i,k,j) - ttend(k)
            ENDDO; ENDDO
         ENDDO
     
     ENDDO

   END SUBROUTINE wave_coupling


   ! ----------------------------------------------------------------------------
   ! ----------------------------------------------------------------------------


   SUBROUTINE vertical_advection_omega(                            &   
                 qtend, ttend, omega_wtg, znu, znw                              &
                 ,phalf_avg, pfull_avg, rho_avg, z_avg, th_avg, qv_avg    & 
                 ,kms,kme, kts,kte, kds,kde                  )

   INTEGER,    INTENT(IN   )    ::  kms,kme, kts,kte, kds,kde
   REAL,       DIMENSION( kms:kme ),INTENT(IN   )    ::      znu,   znw

   INTEGER :: i,J,K,NK
   real, DIMENSION( kms:kme), intent(in) :: z_avg, phalf_avg, pfull_avg, rho_avg, th_avg, qv_avg
   real, DIMENSION( kms:kme) :: omega_wtg
   real :: omega_p, omega_n
   real, DIMENSION( kms:kme),INTENT(INOUT   )  :: qtend,  ttend
   real :: qtend1, qtend2
   real, DIMENSION( 1:kme) :: theta_bg_int, qv_bg_int

      do k=kts,kde ! on model state
        qv_bg_int(k) = qv_avg(k)
        theta_bg_int(k) = th_avg(k)
      enddo


      DO k=kts,kde-1
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
     
     ENDDO
   END SUBROUTINE vertical_advection_omega
   ! ----------------------------------------------------------------------------
   ! ----------------------------------------------------------------------------


   SUBROUTINE vertical_advection_omega_800(                            &   
                 qtend, ttend, omega_wtg, znu, znw                              &
                 ,phalf_avg, pfull_avg, rho_avg, z_avg, th_avg, qv_avg    & 
                 ,kms,kme, kts,kte, kds,kde                  )

   INTEGER,    INTENT(IN   )    ::  kms,kme, kts,kte, kds,kde
   REAL,       DIMENSION( kms:kme ),INTENT(IN   )    ::      znu,   znw

   INTEGER :: i,J,K,NK
   real, DIMENSION( kms:kme), intent(in) :: z_avg, phalf_avg, pfull_avg, rho_avg, th_avg, qv_avg
   real, DIMENSION( kms:kme) :: omega_wtg
   real :: omega_p, omega_n
   real, DIMENSION( kms:kme),INTENT(INOUT   )  :: qtend,  ttend
   real :: qtend1, qtend2
   real, DIMENSION( 1:kme) :: theta_bg_int, qv_bg_int

      do k=kts,kde ! on model state
        qv_bg_int(k) = qv_avg(k)
        theta_bg_int(k) = th_avg(k)
      enddo


      DO k=kts,kde-2
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

     ENDDO
   END SUBROUTINE vertical_advection_omega_800


   ! ----------------------------------------------------------------------------



END MODULE module_force_scm
