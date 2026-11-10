! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Subroutine to calculate cosmogenic 14CO production and
!          returns values to UKCA_EMISSION_CTL.
!          Based on routine HCOX_Cosmogenic_Mod in GEOS-Chem
!          provided by Lee Murray.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_EMISSION_CTL.F90
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
MODULE ukca_cosmogenic_14CO_mod

USE ereport_mod,          ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE ukca_config_constants_mod, ONLY: g
USE ukca_constants,        ONLY: pi_over_180, recip_pi_over_180

IMPLICIT NONE

PRIVATE

! Public subroutines
PUBLIC :: ukca_cosmogenic_14CO_init
PUBLIC :: ukca_cosmogenic_14CO_ctl
   
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_COSMOGENIC_14CO_MOD'

! Should be parameters in suite later
! Data from https://cosmicrays.oulu.fi/phi/Phi_Table_2017.txt
CHARACTER(LEN=*),PARAMETER :: file_solar_modulation = '/home/vagrant/umdir/vn13.9/ctldata/UKCA/cosmo14CO/Phi_Table_1952-2025.txt'
! Data from https://www.ngdc.noaa.gov/geomag/data/poles/NP.xy
CHARACTER(LEN=*),PARAMETER :: file_magnetic_pole = '/home/vagrant/umdir/vn13.9/ctldata/UKCA/cosmo14CO/magnetic_pole_1700-2025.txt'
! Data provided by Stepan Poluianov and Ilya Usoskin
CHARACTER(LEN=*),PARAMETER :: file_c14_prod = '/home/vagrant/umdir/vn13.9/ctldata/UKCA/cosmo14CO/Q_14C_20261109.nc'

INTEGER :: first_year_phi
INTEGER :: last_year_phi

INTEGER :: first_year_mpole
INTEGER :: last_year_mpole

REAL, ALLOCATABLE :: phi_data(:,:)      ! (year_index, month)
REAL, ALLOCATABLE :: mpole_lon_data(:)  ! (year_index)
REAL, ALLOCATABLE :: mpole_lat_data(:)  ! (year_index)

! 14C cosmogenic production table
REAL, ALLOCATABLE :: C14_prod_table(:,:,:)  ! (ndepth, ngcr, nphi)
REAL, ALLOCATABLE :: C14_prod_phi(:)    !
REAL, ALLOCATABLE :: C14_prod_gcr(:)    !
REAL, ALLOCATABLE :: C14_prod_depth(:)  !

CONTAINS

SUBROUTINE ukca_cosmogenic_14CO_init()
  IMPLICIT NONE
  ! Initialise data arrays
  CALL read_phi_data(file_solar_modulation)
  CALL read_magnetic_pole_data(file_magnetic_pole)
  CALL read_c14_prod_table(file_c14_prod)
END SUBROUTINE ukca_cosmogenic_14CO_init

SUBROUTINE ukca_cosmogenic_14CO_ctl(                                                     &
      rows, row_length, p_levels, iyear, imonth, lat, lon, pressure, emiss14CO)

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

! No of latitudes
INTEGER, INTENT(IN) :: rows
! No of longtitudes
INTEGER, INTENT(IN) :: row_length
! No of levels
INTEGER, INTENT(IN) :: p_levels
! Year
INTEGER, INTENT(IN) :: iyear
! Month
INTEGER, INTENT(IN) :: imonth

REAL, INTENT(IN) :: lat(row_length, rows)
REAL, INTENT(IN) :: lon(row_length, rows)
REAL, INTENT(IN) :: pressure(row_length, rows, p_levels)

! Array of cosmogenic 14CO emissions with unit kg(CO?)/m-2/s
REAL, INTENT(OUT) :: emiss14CO(row_length, rows, p_levels)

! Local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ::i,j,k ! loop indices

REAL :: phi_now
REAL :: mpole_lon_now
REAL :: mpole_lat_now
REAL :: gcr
REAL :: depth

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_COSMOGENIC_14CO_CTL'

!       ...Initialize variables and arrays
emiss14CO(:,:,:) = 0.0

! 1) Read in solar modulation potentential phi
WRITE(*,*) 'iyear', iyear
WRITE(*,*) 'imonth', imonth

phi_now = get_phi(iyear, imonth)
WRITE(*,*) 'phi_now', phi_now

! 2) Read in magnetic pole construction

CALL get_magnetic_pole(iyear, mpole_lon_now, mpole_lat_now)
WRITE(*,*) 'mpole_lon_now', mpole_lon_now
WRITE(*,*) 'mpole_lat_now', mpole_lat_now


DO i = 1,rows
  DO j = 1,row_length
    WRITE(*,*) 'longitude, i , j', i, j, lon(j,i)
    WRITE(*,*) 'latitude, i, j',   i, j, lat(j,i)
    
    ! 3) Determine geomagnetic cutoff rigidity in [GV]
    gcr = geomagnetic_cutoff_rigidity( lon(j,i), lat(j,i), mpole_lon_now, mpole_lat_now)
    WRITE(*,*) 'geomagnetic_cutoff_rigidity', gcr
    
    DO k = 1, p_levels 
      !4) Determine atmospheric depth from TOA in [g cm-2]
      ! |dz| * rho = dp / g = (p - 0)/g
      ! Conversion from [kg/m2] to [g cm-2], factor 10-1
      depth = 1.0e-1 * pressure(j,i,k) / g
      !WRITE(*,*) 'depth', depth
      
! 5) Get 14CO production [atoms/g/s]
! Assuming 93% yield of 14C + O2 -> 14CO [Mak et al., 1994]

! 6) Unit conversion      

    END DO
  END DO
END DO



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_cosmogenic_14CO_ctl


SUBROUTINE read_phi_data(filename)

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: filename

  INTEGER :: iu, ios
  INTEGER :: year
  INTEGER :: values(12)
  CHARACTER(LEN=256) :: line

  !----------------------------------------
  ! Determine year range
  !----------------------------------------

  first_year_phi = HUGE(0)
  last_year_phi  = -HUGE(0)

  OPEN(newunit=iu, file=filename, status='OLD', action='READ')

  DO
     READ(iu,'(A)',iostat=ios) line
     IF (ios /= 0) EXIT

     IF (.NOT. is_data_line(line)) CYCLE

     READ(line,*,iostat=ios) year
     IF (ios /= 0) CYCLE

     first_year_phi = MIN(first_year_phi,year)
     last_year_phi  = MAX(last_year_phi ,year)

  END DO

  CLOSE(IU)

  !----------------------------------------
  ! Allocate with actual year bounds
  !----------------------------------------

  ALLOCATE(phi_data(first_year_phi:last_year_phi,12))
  phi_data = -999.0

  !----------------------------------------
  ! Read data
  !----------------------------------------

  OPEN(newunit=iu, file=filename, status='OLD', action='READ')

  DO
     READ(iu,'(A)',iostat=ios) line
     IF (ios /= 0) EXIT

     IF (.NOT. is_data_line(line)) CYCLE

     READ(line,*,iostat=ios) year, values
     IF (ios /= 0) CYCLE

     phi_data(YEAR,:) = REAL(values)

  END DO

  CLOSE(iu)

END SUBROUTINE read_phi_data


FUNCTION get_phi(year,month) RESULT(value_phi)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: year
  INTEGER, INTENT(IN) :: month

  ! Return value
  REAL :: value_phi

  ! Local
  CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message
  INTEGER                            :: errcode    ! Error code for ereport

  IF (year < first_year_phi .OR. year > last_year_phi) THEN
     WRITE(cmessage,'(A,I0)') 'Year outside of solar modulation potential dataset: ', year
     errcode = year
     CALL ereport (ModuleName, errcode, cmessage)
  END IF

  IF (month < 1 .OR. month > 12) THEN
     WRITE(cmessage,'(A,I0)') 'Invalid month for solar modulation potential: ', month
     errcode = month
     CALL ereport (ModuleName, errcode, cmessage)
  END IF

  value_phi = phi_data(year,month)

END FUNCTION get_phi


LOGICAL FUNCTION is_data_line(line)

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: line

  is_data_line = .FALSE.

  ! Skip headers and separator lines
  IF (LEN_TRIM(line) == 0) RETURN
  IF (INDEX(line,'Year') > 0) RETURN
  IF (INDEX(line,'+++') > 0) RETURN

  IF (line(1:1) >= '0' .AND. line(1:1) <= '9') THEN
     is_data_line = .TRUE.
  END IF

END FUNCTION is_data_line

SUBROUTINE read_magnetic_pole_data(filename)

   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN) :: filename

   INTEGER :: iu
   INTEGER :: year
   INTEGER :: ios

   REAL :: lon
   REAL :: lat
   REAL :: year_real

   ! First pass: determine year range

   first_year_mpole = HUGE(0)
   last_year_mpole  = -HUGE(0)

   OPEN(newunit=iu, file=filename, status='OLD', action='READ')

   DO
      READ(iu,*,iostat=ios) lon, lat, year_real
      IF (ios /= 0) EXIT

      year = NINT(year_real)

      first_year_mpole = MIN(first_year_mpole,year)
      last_year_mpole  = MAX(last_year_mpole,year)
   END DO

   CLOSE(iu)

   ALLOCATE(mpole_lon_data(first_year_mpole:last_year_mpole))
   ALLOCATE(mpole_lat_data(first_year_mpole:last_year_mpole))

   ! Second pass: load data

   OPEN(newunit=iu, file=filename, status='OLD', action='READ')

   DO
      READ(iu,*,iostat=ios) lon, lat, year_real
      IF (ios /= 0) EXIT

      year = NINT(year_real)

      mpole_lon_data(year) = lon
      mpole_lat_data(year) = lat
   END DO

   CLOSE(iu)

END SUBROUTINE read_magnetic_pole_data

SUBROUTINE get_magnetic_pole(year, lon, lat)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: year
  REAL,    INTENT(OUT) :: lon
  REAL,    INTENT(OUT) :: lat

  ! Local
  CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message
  INTEGER                            :: errcode    ! Error code for ereport
  
  IF (year < first_year_mpole .OR. year > last_year_mpole) THEN
    WRITE(cmessage,'(A,I0)') 'Year outside of magnetic pole dataset: ', year
    errcode = year
    CALL ereport (ModuleName, errcode, cmessage)
  END IF

  lon = mpole_lon_data(year)
  lat = mpole_lat_data(year)

END SUBROUTINE get_magnetic_pole


FUNCTION geomagnetic_cutoff_rigidity( lon, lat, pole_lon, pole_lat)

  ! Arguments
  REAL, INTENT(IN) :: lat, lon           ! Input, geograph coords [Degrees]
  REAL, INTENT(IN) :: pole_lat, pole_lon ! Input, magnetic pole location [Degrees]
  ! Return value
  REAL:: geomagnetic_cutoff_rigidity     ! Output [GV]

  ! Local
  REAL :: gmlat

  ! Get geomagnetic latitude
  gmlat = geomag_lat( lon, lat, pole_lon, pole_lat )

  ! Use dipole approximation [Smart and Shea, 2005]
  geomagnetic_cutoff_rigidity = ( 14.5 * ( cos( gmlat * pi_over_180 ) )**4 )

END FUNCTION geomagnetic_cutoff_rigidity

FUNCTION geomag_lat( lon, lat, pole_lon, pole_lat )
  ! Returns the geomagnetic latitude for a given lat, lon and year
  ! i.e., the latitude of the sphere with the central axis
  !       being the Earth's geomagnetic dipole
    
  ! Arguments
  REAL, INTENT(IN) :: lat, lon   ! Input, geographic coords  [Degrees]
  REAL, INTENT(IN) :: pole_lat, pole_lon ! Input, magnectic north pole [Degrees]
  ! Return value
  REAL :: geomag_lat ! Output, geomagnetic lat   [Degrees]
    
  ! Determine the Geomagnetic Latitude relative to mag. North [Radians]
  geomag_lat = asin( sin(lat*pi_over_180) * sin(pole_lat*pi_over_180) + &
        cos(lat*pi_over_180) * cos(pole_lat*pi_over_180) * cos(lon*pi_over_180-pole_lon*pi_over_180) )
    
  ! Convert to Degrees
  geomag_lat = geomag_lat * recip_pi_over_180
    
END FUNCTION geomag_lat


SUBROUTINE read_c14_prod_table(filename)

  USE netcdf
  USE ISO_FORTRAN_ENV, ONLY : INT32

   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN) :: filename

   INTEGER(INT32) :: ncid
   INTEGER(INT32) :: status

   INTEGER(INT32) :: dimid_phi
   INTEGER(INT32) :: dimid_rc
   INTEGER(INT32) :: dimid_depth

   INTEGER(INT32) :: nphi
   INTEGER(INT32) :: nrc
   INTEGER(INT32) :: ndepth

   INTEGER(INT32) :: varid_q
   INTEGER(INT32) :: varid_phi
   INTEGER(INT32) :: varid_rc
   INTEGER(INT32) :: varid_depth

   !----------------------------------------
   ! Open file
   !----------------------------------------

   status = NF90_OPEN(TRIM(filename), NF90_NOWRITE, ncid)
   IF (status /= NF90_NOERR) THEN
      WRITE(*,*) NF90_STRERROR(status)
      STOP
   END IF

   !----------------------------------------
   ! Read dimensions
   !----------------------------------------

   status = NF90_INQ_DIMID(ncid, 'phi', dimid_phi)
   status = NF90_INQUIRE_DIMENSION(ncid, dimid_phi, len=nphi)

   status = NF90_INQ_DIMID(ncid, 'Rc', dimid_rc)
   status = NF90_INQUIRE_DIMENSION(ncid, dimid_rc, len=nrc)

   status = NF90_INQ_DIMID(ncid, 'depthi', dimid_depth)
   status = NF90_INQUIRE_DIMENSION(ncid, dimid_depth, len=ndepth)

   !----------------------------------------
   ! Allocate arrays
   !----------------------------------------

   ALLOCATE(c14_prod_table(ndepth,nrc,nphi))
   ALLOCATE(c14_prod_phi(nphi))
   ALLOCATE(c14_prod_gcr(nrc))
   ALLOCATE(c14_prod_depth(ndepth))

   !----------------------------------------
   ! Get variable IDs
   !----------------------------------------

   status = NF90_INQ_VARID(ncid, 'Q', varid_q)
   status = NF90_INQ_VARID(ncid, 'phi', varid_phi)
   status = NF90_INQ_VARID(ncid, 'Rc', varid_rc)
   status = NF90_INQ_VARID(ncid, 'depthi', varid_depth)

   !----------------------------------------
   ! Read variables
   !----------------------------------------

   status = NF90_GET_VAR(ncid, varid_q,      c14_prod_table)
   status = NF90_GET_VAR(ncid, varid_phi,    c14_prod_phi)
   status = NF90_GET_VAR(ncid, varid_rc,     c14_prod_gcr)
   status = NF90_GET_VAR(ncid, varid_depth,  c14_prod_depth)

   !----------------------------------------
   ! Close file
   !----------------------------------------

   status = NF90_CLOSE(ncid)

END SUBROUTINE read_c14_prod_table
 
END MODULE ukca_cosmogenic_14CO_mod
