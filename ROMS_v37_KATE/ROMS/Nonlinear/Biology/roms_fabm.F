#include "cppdefs.h"
      MODULE roms_fabm_mod
!
! This module houses all the subroutines needed to link ROMS to FABM
!
!!Phil Wallhead 2015

#if defined RFABM && defined BIOLOGY

      implicit none

      PRIVATE
      PUBLIC  :: roms_fabm_init


      CONTAINS



      SUBROUTINE roms_fabm_init

      USE fabm
      USE fabm_config

      USE mod_param   !Provides access to BOUNDS, FMODELS, ndiagi etc.
      USE mod_parallel!Provides access to MyRank
      USE mod_forces  !Provides access to FORCES
      USE mod_grid    !Provides access to GRID
      USE mod_ocean   !Provides access to OCEAN (hence "t")
      USE mod_biology !Provides access to idbio via the python-generated file rfabm_mod.h
      USE mod_scalars !Provides access to ydayc = No. days since start of year

      integer :: ng, tile, Istr, Iend, Jstr, Jend, nx, ny, j, k
      integer :: itrc, ibio, idiag

# ifdef DISTRIBUTE
      tile = MyRank
# else
      tile = 0  !Zero to access the full domain (as for "t", see mod_arrays -> allocate ocean)
# endif
      DO ng=1,Ngrids
        !  Read FABM's runtime config from fabm.yaml in current working directory
        CALL fabm_create_model_from_yaml_file(FMODELS(ng)%fmodel(tile))
        write(*,*) "Done fabm_create_model_from_yaml_file -PWA"

        ndiagi(ng) = size(FMODELS(ng)%fmodel(tile)%diagnostic_variables)
        idiagimin(ng) = 126 !First diagnostic variable (see numbered list from fabm_describe_model.py)
        ndiagimax(ng) = 10  !Maximum number of diagnostic variables to display if dBdt1max exceeded
                        !(must not exceed 100)

        ndiagh(ng) = size(FMODELS(ng)%fmodel(tile)%                     &
     &                   horizontal_diagnostic_variables)
        idiaghmin(ng) = 1 !First horizontal diagnostic variable (see numbered list from fabm_describe_model.py)
        ndiaghmax(ng) = 10  !Maximum number of diagnostic variables to display if dBdt1max exceeded
                        !(must not exceed 100)

        ! Set the bulk diagnostic calculator vector (logical: 1 to calculate, 0 to not calculate)
        DO idiag=1,ndiagi(ng)
          FMODELS(ng)%fmodel(tile)%diagnostic_variables(idiag)%         &
     &                save = .FALSE.
        END DO
        DO idiag=1,10 !Include this loop to output the first 10 as well as the chosen chunk (total must be <= 100 diagnostics)
          FMODELS(ng)%fmodel(tile)%diagnostic_variables(idiag)%         &
     &                save = .TRUE.
        END DO
        DO idiag=idiagimin(ng),MIN(ndiagi(ng),idiagimin(ng)+            &
     &                             ndiagimax(ng)-1)
          FMODELS(ng)%fmodel(tile)%diagnostic_variables(idiag)%         &
     &                save = .TRUE.
        END DO

        ! Set the horizontal diagnostic calculator vector (logical: 1 to calculate, 0 to not calculate)
        DO idiag=1,ndiagh(ng)
          FMODELS(ng)%fmodel(tile)%                                     &
     &        horizontal_diagnostic_variables(idiag)%save = .FALSE.
        END DO
        DO idiag=idiaghmin(ng),MIN(ndiagh(ng),idiaghmin(ng)+            &
     &                             ndiaghmax(ng)-1)
          FMODELS(ng)%fmodel(tile)%                                     &
     &        horizontal_diagnostic_variables(idiag)%save = .TRUE.
        END DO


        !  Provide FABM with the extents of the spatial domain
        Istr = BOUNDS(ng)%Istr(tile) !Note: this is the lower bound index NOT INCLUDING GHOST POINTS (cf. LBi)
        Iend = BOUNDS(ng)%Iend(tile)
        Jstr = BOUNDS(ng)%Jstr(tile)
        Jend = BOUNDS(ng)%Jend(tile)

        nx = Iend-Istr+1
        ny = Jend-Jstr+1

        CALL fabm_set_domain(FMODELS(ng)%fmodel(tile),nx,ny,N(ng))
        write(*,*) "Done fabm_set_domain -PWA"

        CALL FMODELS(ng)%fmodel(tile)%set_surface_index(N(ng))
        write(*,*) "Done fabm set_surface_index -PWA"

        CALL FMODELS(ng)%fmodel(tile)%set_bottom_index(1)
        write(*,*) "Done fabm set_bottom_index -PWA"

        ! write(*,*) "Istr,Iend,Jstr,Jend = ",Istr,Iend,Jstr,Jend
        ! write(*,*) "nx,ny = ",nx,ny
        ! write(*,*) ng,N(ng),Istr,Iend,Jstr,Jend
        ! write(*,*) OCEAN(ng)%t(Istr,Jstr,1,1,1)
        ! write(*,*) NBT
        ! write(*,*) size(FMODELS(ng)%fmodel(tile)%state_variables)
        ! write(*,*) idbio
        ! write(*,*) OCEAN(ng)%t(Istr,Jstr,1,1,NBT)
        ! write(*,*) GRID(ng)%rmask(1,1)

        CALL fabm_set_mask(FMODELS(ng)%fmodel(tile),                    &
     &                  GRID(ng)%rmask(Istr:Iend,Jstr:Jend))
        ! This passes the masking grid at rho points (size[Istr:Iend,Jstr:Jend]) to FABM, see mod_grid.F
        write(*,*) "Done fabm_set_mask -PWA"


        !! Link the FABM state variables
        ! First link the bulk state variables to FABM (do this for nstp=1, cf. ana_initial.h)
        DO itrc=1,NBT
          ibio=idbio(itrc)
          CALL fabm_link_bulk_state_data(FMODELS(ng)%fmodel(tile),      &
     &     itrc,OCEAN(ng)%t(Istr:Iend,Jstr:Jend,1:N(ng),1,ibio))
        END DO
        write(*,*) "Done fabm_link_bulk_state_data calls -PWA"

        ! Next link the surface-attached state variables (if any)
        DO itrc=1,NBS(ng)
          CALL fabm_link_surface_state_data(FMODELS(ng)%fmodel(tile),   &
     &     itrc,OCEAN(ng)%state_sf(Istr:Iend,Jstr:Jend,1,itrc))
        END DO
        write(*,*) "Done fabm_link_surface_state_data calls -PWA"

        ! Next link the bottom-attached (benthic) state variables (if any)
        DO itrc=1,NBB(ng)
          CALL fabm_link_bottom_state_data(FMODELS(ng)%fmodel(tile),    &
     &     itrc,OCEAN(ng)%state_bt(Istr:Iend,Jstr:Jend,1,itrc))
        END DO
        write(*,*) "Done fabm_link_bottom_state_data calls -PWA"





        !  Set the pointers for cell thickness and bgc forcings
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%cell_thickness)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_bulk_data(                 & !cell thickness (m)
     &     standard_variables%cell_thickness,                           &
     &     GRID(ng)%Hz(Istr:Iend,Jstr:Jend,:))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%temperature)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_bulk_data(                 & !temperature [degC]
     &     standard_variables%temperature,                              &
     &     OCEAN(ng)%t(Istr:Iend,Jstr:Jend,1:N(ng),1,1))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%practical_salinity)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_bulk_data(                 & !salinity [psu]
     &     standard_variables%practical_salinity,                       &
     &     OCEAN(ng)%t(Istr:Iend,Jstr:Jend,1:N(ng),1,2))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%density)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_bulk_data(                 & !density [kg/m3]
     &     standard_variables%density,                                  &
     &     OCEAN(ng)%rhof(Istr:Iend,Jstr:Jend,:))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%pressure)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_bulk_data(                 & !pressure [dbar]
     &     standard_variables%pressure,                                 &
     &     OCEAN(ng)%pres(Istr:Iend,Jstr:Jend,:))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%surface_downwelling_shortwave_flux)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_horizontal_data(           & !surface downwelling shortwave flux [W/m2]
     &     standard_variables%surface_downwelling_shortwave_flux,       &
     &     FORCES(ng)%srflxWm2(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%wind_speed)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_horizontal_data(           & !wind speed [m/s]
     &     standard_variables%wind_speed,                               &
     &     FORCES(ng)%Swind(Istr:Iend,Jstr:Jend))
        END IF
!! Absorption of silt (needed by light_iop) is a custom variable for ERSEM, not a standard variable
!! Name may be misleading, it is actually the PAR attenuation due to silt
!! A simple fix is to modify the fabm.yaml file by adding these lines:
! silt_abs:
!   model: bulk_constant
!   parameters:
!     standard_name: absorption_of_silt
!     value: 0.1
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &  standard_variables%mole_fraction_of_carbon_dioxide_in_air)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_horizontal_data(           & !partial pressure of atmospheric CO2 [uatm]
     &  standard_variables%mole_fraction_of_carbon_dioxide_in_air,      & !note: the name in FABM is not quite correct (only equivalent if total atmospheric pressure = 1 atm)
     &     FORCES(ng)%pCO2atm(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%bottom_stress)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_horizontal_data(           & !bottom stress [Pa]
     &     standard_variables%bottom_stress,                            &
     &     FORCES(ng)%bstrPa(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%bottom_depth_below_geoid)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_horizontal_data(           & !bottom depth below geoid [m]
     &     standard_variables%bottom_depth_below_geoid,                 &
     &     GRID(ng)%h(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%longitude)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_horizontal_data(           & !longitude [degrees]
     &     standard_variables%longitude,                                &
     &     GRID(ng)%lonr(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &   standard_variables%latitude)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_horizontal_data(           & !latitude [degrees]
     &     standard_variables%latitude,                                 &
     &     GRID(ng)%latr(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &  standard_variables%number_of_days_since_start_of_the_year)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_scalar(                    & !year day (0-365) [days]
     &     standard_variables%number_of_days_since_start_of_the_year,   &
     &     ydayc(ng))
        END IF



        write(*,*) "Done fabm_link_data calls -PWA"


        ! Set variable ids for pointers that will change depending on time step (via nstp)
        FMODELS(ng)%id_temp(tile) =                                     &
     &    FMODELS(ng)%fmodel(tile)%get_bulk_variable_id(                &
     &     standard_variables%temperature)
        FMODELS(ng)%id_salt(tile) =                                     &
     &    FMODELS(ng)%fmodel(tile)%get_bulk_variable_id(                &
     &     standard_variables%practical_salinity)
        write(*,*) "Done setting runtime variable identifiers -PWA"


        ! Check this fabm model is ready
        CALL fabm_check_ready(FMODELS(ng)%fmodel(tile))
        write(*,*) "Done fabm_check_ready -PWA"


# ifdef FABM_INITIAL
        ! Call fabm_initialize_state for each layer and j coordinate
        ! NOTE: This will OVERRIDE any ROMS tracer initialization via .nc input
        !       --- call to ROMS_initialise occurs before call to roms_fabm_init in ocean.h
        DO k=1,N(ng)
          DO j=1,ny
           CALL fabm_initialize_state(FMODELS(ng)%fmodel(tile),1,nx,j,k)
          END DO
        END DO
        write(*,*) "Done fabm_initialize_state calls -PWA"
        IF (NBS(ng).gt.0) THEN
          DO j=1,ny
           CALL fabm_initialize_surface_state(FMODELS(ng)%fmodel(tile), &
     &                                        1,nx,j)
          END DO
          write(*,*) "Done fabm_surface_state calls -PWA"
        END IF
        IF (NBB(ng).gt.0) THEN
          DO j=1,ny
            CALL fabm_initialize_bottom_state(FMODELS(ng)%fmodel(tile), &
     &                                        1,nx,j)
          END DO
          write(*,*) "Done fabm_bottom_state calls -PWA"
        END IF
# endif

      END DO

      RETURN
      END SUBROUTINE roms_fabm_init

#endif

      END MODULE roms_fabm_mod
