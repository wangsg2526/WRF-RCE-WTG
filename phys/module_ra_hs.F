!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_hs

CONTAINS

!------------------------------------------------------------------
   SUBROUTINE HSRAD(RTHRATEN,p8w,p_phy,pi_phy,dz8w,t_phy,          &
                     t8w, rho_phy, R_d,G,CP,dt,xlat,degrad,        &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte, hs_lw_type                     )



   IMPLICIT NONE

!------------------------------------------------------------------
!------------------------------------------------------------------
   INTEGER,    INTENT(IN   ) ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte  

   REAL, INTENT(IN    )      ::        DEGRAD

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT)  ::                              RTHRATEN

   REAL, INTENT(IN   )   ::                   R_d,CP,G,dt

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN ) ::                                     dz8w, &
                                                             p8w, &
                                                           p_phy, &
                                                          pi_phy, &
                                                           t_phy, &
                                                             t8w, &
                                                         rho_phy  
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(IN ) ::                                     xlat

   INTEGER :: i,j,K,NK
   real :: delty,delthez,p0,sec_p_d,sigb,kka,kks,kkf,rcp
   real :: ttmp,teq,sig,sigterm,kkt,t_tend
   REAL, DIMENSION(  ims:ime, kms:kme ):: zght
   integer,intent(in) :: hs_lw_type

!------------------------------------------------------------------
! Newtonian relaxation scheme from Held and Suarez, Bull. Amer. Met.
! Soc., Vol. 75, No. 10., p1825-1830, 1994.  (box on page 1826)
! CEN and MIR  31-JUL-04

   delty   = 60.0
   delthez = 10.0
   p0      = 100000.0
   sec_p_d = 86400.
   sigb    = 0.7
   kka     = 1.0/40.0   ! units of per day
   kks     = 0.25
   kkf     = 1.0
   rcp     = R_d/CP
   !zght    = 0.0
   
   j_loop: DO J=jts,MIN(jte,jde-1)
   k_loop: DO K=kts,MIN(kte,kde-1)
   i_loop: DO I=its,MIN(ite,ide-1)
        

      if(1==2)then
      ttmp = 315.0 - delty*(sin(xlat(i,j)*degrad))**2.0- &
               delthez*alog(p_phy(i,k,j)/p0)*(cos(xlat(i,j)*degrad))**2.0
               
      teq=max(200.0,ttmp*(p_phy(i,k,j)/p0)**rcp)

      sig=p_phy(i,k,j)/p8w(i,1,j)
      sigterm=max(0.0,(sig-sigb)/(1.0-sigb))

      kkt=kka+(kks-kka)*sigterm*(cos(xlat(i,j)*degrad))**4.0

      t_tend=-kkt*(t_phy(i,k,j)-teq)/sec_p_d  ! t_tend in kinetic K/s

      RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
      endif

!      if(i==3.and.j==4.and.k==5)then
!        print*,'hs_lw_type=',hs_lw_type
!      endif
      
      if(hs_lw_type==1)then
        if(t_phy(i,k,j).gt.207.5)then
            t_tend = -1.5/sec_p_d
        !else
        !    t_tend = (200-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      elseif(hs_lw_type==1000)then ! for GASS-WTG 

        if(p_phy(i,k,j) >= 200e2 )then ! troposphere
            t_tend = -1.5/sec_p_d
        elseif(p_phy(i,k,j) <= 100e2) then ! stratasphere
            t_tend = (200-t_phy(i,k,j))/1.0/sec_p_d ! relax at 1 day time scale
        else ! TTL transition
            t_tend = -1.5/sec_p_d  * (p_phy(i,k,j)-100e2)/100e2  + (200-t_phy(i,k,j))/1.0/sec_p_d  * (200e2-p_phy(i,k,j))/100e2
        endif

        RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)

      elseif(hs_lw_type==203)then
        if(t_phy(i,k,j).gt.(195.0+3.0*5))then
            t_tend = -3.0/sec_p_d
        else
            t_tend = (195.0-t_phy(i,k,j))/5.0/sec_p_d
        endif
        RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
      elseif(hs_lw_type==206)then
        if(t_phy(i,k,j).gt.(195.0+6.0*5))then
            t_tend = -6.0/sec_p_d
        else
            t_tend = (195.0-t_phy(i,k,j))/5.0/sec_p_d
        endif
        RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
      elseif(hs_lw_type==209)then
        if(t_phy(i,k,j).gt.(195.0+9.0*5))then
            t_tend = -9.0/sec_p_d
        else
            t_tend = (195.0-t_phy(i,k,j))/5.0/sec_p_d
        endif
        RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
      elseif(hs_lw_type==220)then
        if(t_phy(i,k,j).gt.(195.0+20.0*5))then
            t_tend = -20.0/sec_p_d
        else
            t_tend = (195.0-t_phy(i,k,j))/5.0/sec_p_d
        endif
        RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
      elseif(hs_lw_type==2)then
            if(k==1)then
                zght(i,k)=0.0
            else 
                zght(i,k) = zght(i,k-1)+dz8w(i,k,j)
            endif

            if(zght(i,k)<12e3)then
               t_tend = -2.0/sec_p_d
            elseif(zght(i,k)<=15e3)then
               t_tend = 2.0/3.0e3*(zght(i,k)-12.0e3)-2.0
               t_tend = t_tend/sec_p_d
            else
               t_tend = 0.0 
            endif
            !if(i==1.and.j==1)then
            !       print*,'t_tend',i,j,k,t_tend*sec_p_d, zght(k)
            !endif
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
      elseif(hs_lw_type==3)then
            if(k==1)then
                zght(i,k)=0.0
            else 
                zght(i,k) = zght(i,k-1)+dz8w(i,k,j)
            endif
            
            if(zght(i,k)<12e3)then
               t_tend = -2.0/sec_p_d
            elseif(zght(i,k)<=15e3)then
               t_tend = 2.0/3.0e3*(zght(i,k)-12.0e3)-2.0
               t_tend = t_tend/sec_p_d
            else
               t_tend = (190.0-t_phy(i,k,j))/5.0/sec_p_d
            endif
            !if(i==1.and.j==1)then
            !       print*,'t_tend',i,j,k,t_tend*sec_p_d, zght(k)
            !endif
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
      elseif(hs_lw_type==4)then
            if(k==1)then
                zght(i,k)=0.0
            else 
                zght(i,k) = zght(i,k-1)+dz8w(i,k,j)
            endif
            
            if(zght(i,k)<12e3)then
               t_tend = -2.0/sec_p_d
            elseif(zght(i,k)<=17e3)then
               t_tend = 2.0/5.0e3*(zght(i,k)-12.0e3)-2.0
               t_tend = t_tend/sec_p_d
            else
               t_tend = (190.0-t_phy(i,k,j))/5.0/sec_p_d
            endif
            !if(i==1.and.j==1)then
            !       print*,'t_tend',i,j,k,t_tend*sec_p_d, zght(i,k)
            !endif
      elseif(hs_lw_type==10)then
        if(t_phy(i,k,j).gt.187.5)then
            t_tend = -1.5/sec_p_d
        !else
        !    t_tend = (200-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      elseif(hs_lw_type==11)then
        if(t_phy(i,k,j).gt.167.5)then
            t_tend = -1.5/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      elseif(hs_lw_type==12)then
        if(t_phy(i,k,j).gt.177.5)then
            t_tend = -1.5/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      elseif(hs_lw_type==101)then
        if(t_phy(i,k,j).gt.189.15)then
            t_tend = -1.2/sec_p_d
        !else
        !    t_tend = (200-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
     elseif(hs_lw_type==103)then
        if(t_phy(i,k,j).gt.209.15)then
            t_tend = -1.2/sec_p_d
        !else
        !    t_tend = (200-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      else
           stop 'Wrong lw_type!'
      endif
     !if(i==20.and.j==20.and.(k==10.or.k==40))then
     !    print*,'t_tend',i,j,k,RTHRATEN(I,K,J)
     !endif

   ENDDO i_loop
   ENDDO k_loop
   ENDDO j_loop                                          

   END SUBROUTINE HSRAD

!====================================================================
   SUBROUTINE hsinit(RTHRATEN,restart,                              &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte                   )
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           :: restart
   INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,  &
                                     ims, ime, jms, jme, kms, kme,  &
                                     its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::        &
                                                          RTHRATEN
   INTEGER :: i, j, k, itf, jtf, ktf

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   IF(.not.restart)THEN
     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RTHRATEN(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO
   ENDIF

   END SUBROUTINE hsinit








!====================================================================
   SUBROUTINE HSRAD_SW(RTHRATEN,p8w,p_phy,pi_phy,dz8w,t_phy,          &
                     t8w, rho_phy, R_d,G,CP,dt,xlat,degrad,        &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte, hs_sw_type                   )



   IMPLICIT NONE

!------------------------------------------------------------------
!------------------------------------------------------------------
   INTEGER,    INTENT(IN   ) ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte  

   REAL, INTENT(IN    )      ::        DEGRAD

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT)  ::                              RTHRATEN

   REAL, INTENT(IN   )   ::                   R_d,CP,G,dt

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN ) ::                                     dz8w, &
                                                             p8w, &
                                                           p_phy, &
                                                          pi_phy, &
                                                           t_phy, &
                                                             t8w, &
                                                         rho_phy  
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(IN ) ::                                     xlat

   INTEGER :: i,j,K,NK
   real :: delty,delthez,p0,sec_p_d,sigb,kka,kks,kkf,rcp
   real :: ttmp,teq,sig,sigterm,kkt,t_tend
   
 integer,intent(in) :: hs_sw_type

!------------------------------------------------------------------
! Newtonian relaxation scheme from Held and Suarez, Bull. Amer. Met.
! Soc., Vol. 75, No. 10., p1825-1830, 1994.  (box on page 1826)
! CEN and MIR  31-JUL-04

   delty   = 60.0
   delthez = 10.0
   p0      = 100000.0
   sec_p_d = 86400.
   sigb    = 0.7
   kka     = 1.0/40.0   ! units of per day
   kks     = 0.25
   kkf     = 1.0
   rcp     = R_d/CP

   
   j_loop: DO J=jts,MIN(jte,jde-1)
   k_loop: DO K=kts,MIN(kte,kde-1)
   i_loop: DO I=its,MIN(ite,ide-1)
        

      if(1==2)then
      ttmp = 315.0 - delty*(sin(xlat(i,j)*degrad))**2.0- &
               delthez*alog(p_phy(i,k,j)/p0)*(cos(xlat(i,j)*degrad))**2.0
               
      teq=max(200.0,ttmp*(p_phy(i,k,j)/p0)**rcp)

      sig=p_phy(i,k,j)/p8w(i,1,j)
      sigterm=max(0.0,(sig-sigb)/(1.0-sigb))

      kkt=kka+(kks-kka)*sigterm*(cos(xlat(i,j)*degrad))**4.0

      t_tend=-kkt*(t_phy(i,k,j)-teq)/sec_p_d  ! t_tend in kinetic K/s

      RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
      endif
!      if(i==3.and.j==4.and.k==5)then
!        print*,'hs_sw_type=',hs_sw_type
!      endif
      if(hs_sw_type==1)then
        if(t_phy(i,k,j).le.207.5)then
        !    t_tend = -1.5/sec_p_d
        !else
            t_tend = (200-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      elseif(hs_sw_type==1000)then ! for GASS-WTG
         ! do nothing here!
      elseif(hs_sw_type==2)then
      elseif(hs_sw_type==3)then
        ! do nothing
        if(t_phy(i,k,j).le.190.0)then
        !    t_tend = -1.5/sec_p_d
        !else
            t_tend = (190.0-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      
      elseif(hs_sw_type==10)then
        if(t_phy(i,k,j).le.187.5)then
            t_tend = (180.0-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      elseif(hs_sw_type==11)then
        if(t_phy(i,k,j).le.167.5)then
            t_tend = (160.0-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
      elseif(hs_sw_type==12)then
        if(t_phy(i,k,j).le.177.5)then
            t_tend = (170.0-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif

     elseif(hs_sw_type==101)then
        if(t_phy(i,k,j).le.189.15)then
        !    t_tend = -1/sec_p_d
        !else
            t_tend = (183.15-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
     elseif(hs_sw_type==103)then
        if(t_phy(i,k,j).le.209.15)then
        !    t_tend = -1/sec_p_d
        !else
            t_tend = (203.15-t_phy(i,k,j))/5.0/sec_p_d
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
        endif
     elseif(hs_sw_type > 200)then
     else
       stop 'Wrong hs_sw_type !'
     endif
     !if(i==20.and.j==20.and.(k==10.or.k==40))then
     !    print*,'t_tend',i,j,k,RTHRATEN(I,K,J)
     !endif

   ENDDO i_loop
   ENDDO k_loop
   ENDDO j_loop                                          

   END SUBROUTINE HSRAD_SW

!====================================================================
   SUBROUTINE hsinit_sw(RTHRATEN,restart,                              &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte                   )
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           :: restart
   INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,  &
                                     ims, ime, jms, jme, kms, kme,  &
                                     its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::        &
                                                          RTHRATEN
   INTEGER :: i, j, k, itf, jtf, ktf

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   IF(.not.restart)THEN
     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RTHRATEN(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO
   ENDIF

   END SUBROUTINE hsinit_sw

!====================================================================


END MODULE module_ra_hs
