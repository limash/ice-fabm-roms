      SUBROUTINE biology (ng,tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group         Phil Wallhead   !
!                                                   Andre Staalstrom   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!  This routine computes the  biological sources and sinks for the     !
!  any model in the FABM family. Then, it adds those terms to the      !
!  global biological fields.                                           !
!  The code below is adapted from npzd_Franks.h                        !
!                                                                      !
!***********************************************************************

      USE mod_param     !PWA: Provides access to BOUNDS and FMODELS
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean     !PWA: Provides access to OCEAN (hence "t")
      USE mod_stepping
!!!Inserted PWA 08/03/2017
      USE mod_forces
      USE fabm
      USE fabm_config
!!!End insertion PWA 08/03/2017

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
!!!PWA Inserted 08/03/2017
     &                   FORCES(ng) % srflx,                            &
     &                   OCEAN(ng) % state_sf,                          &
     &                   OCEAN(ng) % state_bt,                          &
!!!End insertion PWA 08/03/2017
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
!!!PWA Modified 08/03/2017
!     &                         Hz, z_r, z_w,                            &
!     &                         t)
     &                         Hz, z_r, z_w, srflx,                     &
     &                         state_sf, state_bt, t)
!!!End modification PWA 08/03/2017
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology !PWA: This provides access to idbio
      USE mod_ncparam
      USE mod_scalars
!!!PWA Inserted 08/03/2017
      USE fabm
      USE fabm_config
!!!End insertion PWA 08/03/2017
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
!!!PWA Inserted 08/03/2017
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: state_sf(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: state_bt(LBi:,LBj:,:,:)
!!!End insertion PWA 08/03/2017
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
!!!PWA Inserted 08/03/2017
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: state_sf(LBi:UBi,LBj:UBj,3,NBS(ng))
      real(r8), intent(inout) :: state_bt(LBi:UBi,LBj:UBj,3,NBB(ng))
!!!End insertion PWA 08/03/2017
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!!!PWA: Note that by declaring the variables with index ranges, the index labels are also passed to the subroutine.
!!!     Here we pass chunks defined by tile boundaries including ghost points (LBi, UBi, etc.)
!!!     even though this subroutine only acts on ranges Istr:Iend etc. which exclude ghost points.

!
!  Local variable declarations.
!
!!!PWA Modified 08/03/2017
!      integer, parameter :: Nsink = 1
!      ...
!      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
      integer :: Iter, i, ibio, itrc, ivar, j, k, cnt, ks, nx, j0 !, isink
!      integer :: icheckmax
!!!PWA: Rather than prescribe here the number of sinking variables,
!!!easier to just loop over NBT and take action for i where max(w(:,k,i))>0

! local variables
      real(r8)  yday, hour, dBdt1 !, dBdt1max
      integer :: iday, month, year
!      real(r8) :: dtdays
      real(r8) :: cff, cffL, cffR, cu, dltL, dltR

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: w !ROMS vertical sinking velocity [m/s] positive downward
      integer, dimension(IminS:ImaxS,N(ng)) :: ksource
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: dBio_dt    !FABM SMS term (units [C/s])
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ext_bio           !FABM PAR attenuation due to biology [/m]
      real(r8), dimension(IminS:ImaxS,NBT) :: flux_sf             !FABM surface fluxes (units [Matter/m2/s])
      real(r8), dimension(IminS:ImaxS,NBS(ng)) :: sms_sf          !FABM surface-attached SMS terms (units [Matter/m2/s])
      real(r8), dimension(IminS:ImaxS,NBT) :: flux_bt             !FABM bottom fluxes (units [Matter/m2/s])
      real(r8), dimension(IminS:ImaxS,NBB(ng)) :: sms_bt          !FABM bottom-attached SMS terms (units [Matter/m2/s])
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

      logical,parameter :: repair = .false.
      logical :: valid_int,valid_sf,valid_bt
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),100) :: diagi
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,100) :: diagh
      integer :: idiag
!!!End modification PWA 08/03/2017

#include "set_bounds.h" !PWA: Sets Istr,Iend etc. using BOUNDS (note that LBi,UBi etc. are passed through subroutine args)



!!!PWA Inserted 08/03/2017
!Link environmental data to FABM (for variables that must be relinked for each time step since nstp varies)
      CALL FMODELS(ng)%fmodel(tile)%link_bulk_data(                     &
     &                  FMODELS(ng)%id_temp(tile),                      &
     &                  t(Istr:Iend,Jstr:Jend,1:N(ng),nstp,1))
      CALL FMODELS(ng)%fmodel(tile)%link_bulk_data(                     &
     &                  FMODELS(ng)%id_salt(tile),                      &
     &                  t(Istr:Iend,Jstr:Jend,1:N(ng),nstp,2))

!Link biogeochemical state variables to FABM (must be redone for each time step since nstp varies)
      DO itrc=1,NBT
        ibio=idbio(itrc)
        CALL fabm_link_bulk_state_data(FMODELS(ng)%fmodel(tile),        &
     &   itrc,t(Istr:Iend,Jstr:Jend,1:N(ng),nstp,ibio))
      END DO

      DO itrc=1,NBS(ng)
        CALL fabm_link_surface_state_data(FMODELS(ng)%fmodel(tile),     &
     &   itrc,state_sf(Istr:Iend,Jstr:Jend,nstp,itrc))
      END DO

      DO itrc=1,NBB(ng)
        CALL fabm_link_bottom_state_data(FMODELS(ng)%fmodel(tile),      &
     &   itrc,state_bt(Istr:Iend,Jstr:Jend,nstp,itrc))
      END DO
!!!End insertion PWA 08/03/2017


!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------

!!!PWA Inserted 08/03/2017
!      dtdays=dt(ng)*sec2day
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
      ydayc(ng) = yday-1.0_r8

      IF (iic(ng).le.icheckmax(ng)) THEN
!        write(*,*) "PWA: dtdays = ", dtdays
        write(*,*) "tdays = ", tdays(ng)
        write(*,*) "yday, ydayc = ", yday, ydayc(ng)
      END IF
      nx = Iend-Istr+1
!      dBdt1max = 1e3   !Maximum absolute rate of change from FABM model (def = 1e3)
                     !Set to zero to examine initial FABM output
!      icheckmax = 0  !Maximum time step counter value (iig) for which MIN/MAX dBio_dt is checked
                     !Set to zero to switch off all checks
                     !Set to large value to check all time steps
!!!End insertion PWA 08/03/2017

!!!Deleted PWA 08/03/2017
!
!  Avoid computing source/sink terms if no biological iterations.
!
!      IF (BioIter(ng).le.0) RETURN
!      ...
!      Wbio(1)=wDet(ng)                ! Small detritus
!!!End deletion PWA 08/03/2017

!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
!!!Inserted PWA 08/03/2017
        j0 = j-Jstr+1 !FABM chunks are always indexed 1:chunksize
                      !Therefore, all calls to FABM APIs should use index j0 instead of j
!!!End insertion PWA 08/03/2017
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!

!!!Deleted PWA 08/03/2017
! ...
!  Determine Correction for negativity.
! ...
!  If correction needed, determine the largest pool to debit.
! ...
!  Update new values.
! ...
!!!PWA: This could be a useful way to conserve mass, but for now we stick to simple capping
!!!     following the fennel.h approach (see below)
!!!Deleted PWA 08/03/2017

!!!Inserted PWA 08/03/2017
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              Bio(i,k,ibio)=Bio_old(i,k,ibio)
            END DO
          END DO
        END DO
!
!  Extract potential temperature and salinity imposing caps.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8) !PWA: Should not be necessary
!            Bio(i,k,itemp)=t(i,j,k,nstp,itemp)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
          END DO
        END DO
!!!End insertion PWA 08/03/2017


!!!Deleted PWA 08/03/2017
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  ...
!
!          END DO
!!!Deleted PWA 08/03/2017



!!!Inserted PWA 08/03/2017

!-----------------------------------------------------------------------
!  Check the max/min surface PAR for this j if req'd
!-----------------------------------------------------------------------
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Istr,Iend,Jstr,Jend = ",Istr,Iend,Jstr,Jend
          write(*,*) "iic,icheckmax,j,j0 = ",iic(ng),icheckmax(ng),j,j0
          write(*,*) "MIN,MAX(rmask(Istr:Iend,j)) = ",                  &
     &      MINVAL(rmask(Istr:Iend,j)), MAXVAL(rmask(Istr:Iend,j))
          write(*,*) "MIN,MAX(srflx(Istr:Iend,j)) = ",                  &
     &      MINVAL(srflx(Istr:Iend,j)), MAXVAL(srflx(Istr:Iend,j))
        END IF

!-----------------------------------------------------------------------
!  Generate the biological contributions to PAR attenuation
!-----------------------------------------------------------------------
! (Vectorized over i, looping over k)
        DO k=1,N(ng)
          CALL fabm_get_light_extinction(FMODELS(ng)%fmodel(tile),      &
     &        1,nx,j0,k,ext_bio(Istr:Iend,k))
        END DO
! Check the output max/min if req'd
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "MIN,MAX(ext_bio(Istr:Iend,1:N(ng))) = ",          &
     &      MINVAL(ext_bio(Istr:Iend,:)), MAXVAL(ext_bio(Istr:Iend,:))
          IF (MAXVAL(rmask(Istr:Iend,j)).eq.1) THEN
!            write(*,*) "srflx(Istr:Iend,j) = ",srflx(Istr:Iend,j)
          END IF
        END IF

!-----------------------------------------------------------------------
!  Generate the light field (one column at a time, looping over i)
!-----------------------------------------------------------------------
        DO i=1,nx
          IF (rmask(Istr+i-1,Jstr+j-1).eq.1) THEN
!            write(*,*) ext_bio(Istr+i-1,:)
            CALL fabm_get_light(FMODELS(ng)%fmodel(tile),1,N(ng),i,j0)
          END IF
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_get_light"
        END IF

!-----------------------------------------------------------------------
!  Generate surface fluxes and SMS for surface-attached variables
!-----------------------------------------------------------------------
        flux_sf = 0.0_r8
        sms_sf = 0.0_r8
        CALL fabm_do_surface(FMODELS(ng)%fmodel(tile),1,nx,j0,          &
     &    flux_sf(Istr:Iend,1:NBT),sms_sf(Istr:Iend,1:NBS(ng)))
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do_surface"
          write(*,*) "MIN,MAX(flux_sf(Istr:Iend,:)) = ",                &
     &      MINVAL(flux_sf(Istr:Iend,:)), MAXVAL(flux_sf(Istr:Iend,:))
          write(*,*) "MIN,MAX(sms_sf(Istr:Iend,:)) = ",                 &
     &      MINVAL(sms_sf(Istr:Iend,:)), MAXVAL(sms_bt(Istr:Iend,:))
        END IF

!-----------------------------------------------------------------------
!  Generate bottom fluxes and SMS terms for bottom-attached variables
!-----------------------------------------------------------------------
        flux_bt = 0.0_r8
        sms_bt = 0.0_r8
        CALL fabm_do_bottom(FMODELS(ng)%fmodel(tile),1,nx,j0,           &
     &    flux_bt(Istr:Iend,1:NBT),sms_bt(Istr:Iend,1:NBB(ng)))
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do_bottom"
          write(*,*) "MIN,MAX(flux_bt(Istr:Iend,:)) = ",                &
     &      MINVAL(flux_bt(Istr:Iend,:)), MAXVAL(flux_bt(Istr:Iend,:))
          write(*,*) "MIN,MAX(sms_bt(Istr:Iend,:)) = ",                 &
     &      MINVAL(sms_bt(Istr:Iend,:)), MAXVAL(sms_bt(Istr:Iend,:))
        END IF

!-----------------------------------------------------------------------
!  Generate SMS terms for pelagic variables
!-----------------------------------------------------------------------
        DO k=1,N(ng)
          dBio_dt(Istr:Iend,k,1:NT(ng))=0.0_r8
          CALL fabm_do(FMODELS(ng)%fmodel(tile),1,nx,j0,k,              &
     &              dBio_dt(Istr:Iend,k,NT(ng)-NBT+1:NT(ng)))
!Note: the biological tracer indices idbio range over (NAT+NPT+NCS+NNS)+1:NT = NT-NBT+1:NT, see rfabm_mod.h
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do"
          write(*,*) "MIN,MAX(dBio_dt(Istr:Iend,:,NT-NBT+1:NT)) = ",    &
     &      MINVAL(dBio_dt(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng))),     &
     &      MAXVAL(dBio_dt(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng)))
        END IF

!-----------------------------------------------------------------------
!  Add contributions from surface and bottom fluxes
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO i=Istr,Iend
            dBio_dt(i,N(ng),ibio) = dBio_dt(i,N(ng),ibio) +             &
     &            flux_sf(i,itrc)*Hz_inv(i,N(ng))
            dBio_dt(i,1,ibio) = dBio_dt(i,1,ibio) +                     &
     &            flux_bt(i,itrc)*Hz_inv(i,1)
          !E.g surface/bottom attached variables may have [mass/m2] while pelagic variables may have [mass/m3]
          END DO
        END DO

!-----------------------------------------------------------------------
!  Check the model state and FABM computations in detail if req'd
!-----------------------------------------------------------------------
        IF (iic(ng).le.icheckmax(ng)) THEN

          ! Check the model state
          DO k=1,N(ng)
            CALL fabm_check_state(FMODELS(ng)%fmodel(tile),1,nx,j0,k,   &
     &                            repair,valid_int)
            IF (.not.(valid_int.or.repair)) THEN
              write(*,*) "fabm_check_state returns invalid"
              STOP
            END IF
          END DO
          CALL fabm_check_surface_state(FMODELS(ng)%fmodel(tile),1,nx,  &
     &                                  j0,repair,valid_sf)
          CALL fabm_check_bottom_state(FMODELS(ng)%fmodel(tile),1,nx,   &
     &                                  j0,repair,valid_bt)
          IF (.not.(valid_sf.and.valid_bt).and..not.repair) STOP

          dBdt1 = MAXVAL(ABS(dBio_dt(Istr:Iend,:,NT(ng)-NBT+1:NT(ng))))

          IF (dBdt1.ge.dBdt1max(ng)) THEN
            write(*,*) "dBdt1 = MAX(ABS(dBio_dt1)) = ", dBdt1

            cnt = 0
            DO idiag=1,ndiagi(ng)
              IF (FMODELS(ng)%fmodel(tile)%diagnostic_variables(idiag)% &
     &                        save) THEN
                cnt = cnt + 1
                diagi(Istr:Iend,Jstr:Jend,1:N(ng),cnt) =                &
     &  fabm_get_bulk_diagnostic_data(FMODELS(ng)%fmodel(tile),idiag)
              END IF
            END DO

            cnt = 0
            DO idiag=1,ndiagh(ng)
              IF (FMODELS(ng)%fmodel(tile)%                             &
     &            horizontal_diagnostic_variables(idiag)%save) THEN
                cnt = cnt + 1
                diagh(Istr:Iend,Jstr:Jend,cnt) =                        &
     &  fabm_get_horizontal_diagnostic_data(FMODELS(ng)%fmodel(tile),   &
     &  idiag)
                write(*,*) idiag, FMODELS(ng)%fmodel(tile)%             &
     &          horizontal_diagnostic_variables(idiag)%name
                write(*,*) diagh(Istr:Iend,j,cnt)
              END IF
            END DO

            DO k=1,N(ng)
              dBdt1 = MAXVAL(ABS(dBio_dt(Istr:Iend,k,                   &
     &                   NT(ng)-NBT+1:NT(ng))))
              write(*,*) "MAX(ABS(dBio_dt1)) for k = ", k, " = ",dBdt1

              !Output temperature and salinity for this layer
              write(*,*) "k = ", k, " temperature:"
              write(*,*) t(Istr:Iend,j,k,nstp,1)
              write(*,*) "k = ", k, " salinity:"
              write(*,*) t(Istr:Iend,j,k,nstp,2)

              !Loop over NBT, outputting B and dBdt for each variable
              DO itrc=1,NBT
                ibio=idbio(itrc)
                write(*,*) "k, itrc, ibio = ", k, itrc, ibio
                write(*,*) FMODELS(ng)%fmodel(tile)%                    &
     &                     state_variables(itrc)%name,                  &
     &                     FMODELS(ng)%fmodel(tile)%                    &
     &                     state_variables(itrc)%long_name,             &
     &                     FMODELS(ng)%fmodel(tile)%                    &
     &                     state_variables(itrc)%units
                write(*,*) "t1 = ", t(Istr:Iend,j,k,nstp,ibio)
                write(*,*) "dBio_dt1 = ", dBio_dt(Istr:Iend,k,ibio)
                IF (k.eq.1) THEN
                  write(*,*) "flux_bt for itrc = ", itrc, " = "
                  write(*,*) flux_bt(Istr:Iend,itrc)
                END IF
                IF (k.eq.N(ng)) THEN
                  write(*,*) "flux_sf for itrc = ", itrc, " = "
                  write(*,*) flux_sf(Istr:Iend,itrc)
                END IF
              END DO

              !Output the cell thickness and biological attenuation for this layer
              write(*,*) "k = ", k
              write(*,*) "Hz1 = ", Hz(Istr:Iend,j,k)
              write(*,*) "ext_bio1 = ", ext_bio(Istr:Iend,k)
              !NOTE: ext_bio is the total extinction [/m] due to biological material (phyto + detritus + DOM)

              cnt = 0
              DO idiag=1,ndiagi(ng)
                IF (FMODELS(ng)%fmodel(tile)%                           &
     &                          diagnostic_variables(idiag)%save) THEN
                  cnt = cnt + 1
                  write(*,*) idiag, FMODELS(ng)%fmodel(tile)%           &
     &                       diagnostic_variables(idiag)%name
                  write(*,*) diagi(Istr:Iend,j,k,cnt)
                END IF
              END DO
            END DO

            !Finally output the surface shortwave and max/min cell thickness
            write(*,*) "srflx1 = ", srflx(Istr:Iend,j)
            write(*,*) "MIN,MAX(Hz) = ", MINVAL(Hz), MAXVAL(Hz)
!            STOP
          END IF
        END IF

!-----------------------------------------------------------------------
!  Update tracers with rates of change from FABM (Euler step)
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,ibio)=Bio(i,k,ibio)+dBio_dt(i,k,ibio)*dt(ng)
            END DO
          END DO
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Updated water column with FABM dBio_dt"
          write(*,*) "MIN,MAX(Bio(Istr:Iend,1:N(ng),NT-NBT+1:NT)) = ",  &
            MINVAL(Bio(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng))),         &
     &      MAXVAL(Bio(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng)))
        END IF

!-----------------------------------------------------------------------
!  Update surface states with rate of change from FABM
!-----------------------------------------------------------------------
        DO itrc=1,NBS(ng)
          DO i=Istr,Iend
            state_sf(i,j,nnew,itrc) = MAX(state_sf(i,j,nstp,itrc) +     &
     &                                 sms_sf(i,itrc)*dt(ng), 0.0_r8)
          !Note: we DO impose a zero lower bound on non-tracer state variables (cf. "t" below)
          END DO
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Updated surface states with FABM sms_sf"
          write(*,*) "MIN,MAX(state_sf(Istr:Iend,j,nnew,1:NBS(ng))) = ",&
            MINVAL(state_sf(Istr:Iend,j,nnew,1:NBS(ng))),               &
     &      MAXVAL(state_sf(Istr:Iend,j,nnew,1:NBS(ng)))
        END IF

!-----------------------------------------------------------------------
!  Update bottom states with rate of change from FABM
!-----------------------------------------------------------------------
        DO itrc=1,NBB(ng)
          DO i=Istr,Iend
            state_bt(i,j,nnew,itrc) = MAX(state_bt(i,j,nstp,itrc) +     &
     &                                 sms_bt(i,itrc)*dt(ng), 0.0_r8)
          !Note: we DO impose a zero lower bound on non-tracer state variables (cf. "t" below)
          END DO
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Updated bottom states with FABM sms_bt"
          write(*,*) "MIN,MAX(state_bf(Istr:Iend,j,nnew,1:NBB(ng))) = ",&
            MINVAL(state_bt(Istr:Iend,j,nnew,1:NBB(ng))),               &
     &      MAXVAL(state_bt(Istr:Iend,j,nnew,1:NBB(ng)))
        END IF

!-----------------------------------------------------------------------
!  Get vertical sinking velocities from FABM
!-----------------------------------------------------------------------
        DO k=1,N(ng)
          CALL fabm_get_vertical_movement(FMODELS(ng)%fmodel(tile),     &
     &              1,nx,j0,k,w(Istr:Iend,k,NT(ng)-NBT+1:NT(ng)))
        END DO
        w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng)) = -1*                  &
     &      w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng))
        !FABM outputs vertical velocities in [m/s] positive upward
        !ROMS code below expects w in [m/s] positive downward
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_get_vertical_movement"
          write(*,*) "MIN,MAX(w(Istr:Iend,1:N(ng),NT-NBT+1:NT)) = ",    &
            MINVAL(w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng))),           &
     &      MAXVAL(w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng)))
        END IF
!!!End insertion PWA 08/03/2017



!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
!!!Modified PWA 08/03/2017
!        SINK_LOOP: DO isink=1,Nsink !!!PWA: This is not convenient for arbitrary fabm model
!                                    !!!PWA: Instead we loop over all biol tracers and use IF statement (leaves indentation unchanged since no BioIter loop)
!        ibio=idsink(isink)
        SINK_LOOP: DO itrc=1,NBT
          ibio=idbio(itrc)
          IF (MAXVAL(ABS(w(Istr:Iend,1:N(ng),ibio))).gt.0) THEN
!!!End modification PWA 08/03/2017
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,ibio)
              END DO
            END DO
  !
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
!!!Modified PWA 08/03/2017
!            cff=dtdays*ABS(Wbio(isink))
!!!End modification PWA 08/03/2017
            DO k=1,N(ng)
              DO i=Istr,Iend
!!!Inserted PWA 08/03/2017
                cff=dt(ng)*ABS(w(i,k,ibio))
!!!End insertion PWA 08/03/2017
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO
!!!Modified PWA 08/03/2017
!          END DO SINK_LOOP
!        END DO ITER_LOOP
          END IF
        END DO SINK_LOOP
!!!End modification PWA 08/03/2017
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of nutrients
!  when advection causes tracer concentration to go negative.
!-----------------------------------------------------------------------
!!!NOTE (PWA): this is slightly different to ROMS 3.5 approach
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO

      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile
