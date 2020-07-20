!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module forcing_sst

!BOP
! !MODULE: forcing_sst
!
! !DESCRIPTION:
!  Contains routines and variables used for determining when to send
!  the coupler observational SST instead of POP SST (or a weighted
!  average of the two).
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES

   use kinds_mod
   use domain
   use constants
   use broadcast
   use io
   use forcing_tools
   use time_management
!   use prognostic
!   use grid
   use exit_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_sst,      &
             get_sst_data

! !PUBLIC DATA MEMBERS:

   real (r8), public :: &! public for use in restart
      sst_interp_last    ! time last interpolation was done

   real (r8), public, dimension(:,:,:,:,:), allocatable :: &
      OBS_SST_DATA  ! data to possibly send to coupler
                    ! One dimension is depth, which should be 1
                    ! for surface temp

   logical (log_kind), public :: &
      sst_use_obs

   character (char_len), public :: &
      sst_data_type                ! keyword for period of forcing data

   real (r8), public, dimension(nx_block,ny_block,max_blocks_clinic) :: &
                                                              sst_alpha

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  internal module variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(12) :: sst_data_time

   real (r8), dimension(20) :: &
      sst_data_renorm   ! factors to convert data to model units

   real (r8) ::        &
      sst_data_inc,    &! time increment between values of forcing data
      sst_data_next,   &! time to be used for the next value of forcing data that is needed
      sst_data_update, &! time when new forcing value to be added to interpolation set
      sst_interp_inc,  &! time increment between interpolation
      sst_interp_next   ! time next interpolation will be done

   integer (int_kind) ::    &
      sst_interp_order,     &! order of temporal interpolation
      sst_data_time_min_loc  ! index of the third dimension of sst_data_time containing the minimum forcing time

   character (char_len) :: &
      sst_filename,        &! name of file conainting forcing data
      sst_file_fmt,        &! format (bin or netcdf) of forcing file
      sst_mask_filename,   &! name of file conainting forcing data
      sst_mask_file_fmt,   &! format (bin or netcdf) of forcing file
      sst_interp_freq,     &! keyword for period of temporal interpolation
      sst_interp_type,     &
      sst_data_label,      &
      sst_formulation

   character (char_len), dimension(:), allocatable :: &
      sst_data_names    ! names for required input data fields

   integer (int_kind), dimension(:), allocatable :: &
      sst_bndy_loc,    &! location and field type for
      sst_bndy_type     !   ghost cell updates

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_sst
! !INTERFACE:

 subroutine init_sst

! !DESCRIPTION:
!  Initializes observational SST forcing by reading in the data.
!
! !REVISION HISTORY:
!  same as module

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! dummy loop index
      nml_error           ! namelist i/o error flag

   character (char_len) :: &
      forcing_filename,    &! full filename of forcing data file
      long_name             ! long name for input data field

   type (datafile) :: &
      sst_data_file,  &! data file descriptor for interior pot temp data
      mask_data_file   ! data file descriptor for interior pot temp data

   type (io_field_desc) :: &
      sst_data_in          ! io field descriptor for input pot temp data

   type (io_dim) :: &
      i_dim, j_dim   ! dimension descriptors for horiz dims

   namelist /forcing_sst_nml/ sst_data_type,           &
        sst_interp_type,      sst_interp_freq,         &
        sst_data_inc,         sst_interp_inc,          &
        sst_filename,         sst_file_fmt,            &
        sst_mask_filename,    sst_mask_file_fmt,       &
        sst_use_obs

!-----------------------------------------------------------------------
!
!  read observational SST namelist input
!    after setting default values.
!
!-----------------------------------------------------------------------

   sst_data_type     = 'none'
   sst_interp_type   = 'nearest'
   sst_interp_freq   = 'never'
   sst_filename      = 'unknown-sst'
   sst_file_fmt      = 'nc'
   sst_mask_filename = 'unknown-mask'
   sst_mask_file_fmt = 'nc'
   sst_use_obs       = .false.
   sst_formulation   = 'restoring'
   sst_data_renorm   = c1

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
         read(nml_in, nml=forcing_sst_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR: reading forcing_sst_nml')
   endif

   call broadcast_scalar(sst_data_type,         master_task)
   call broadcast_scalar(sst_interp_type,       master_task)
   call broadcast_scalar(sst_interp_freq,       master_task)
   call broadcast_scalar(sst_filename,          master_task)
   call broadcast_scalar(sst_file_fmt,          master_task)
   call broadcast_scalar(sst_mask_filename,     master_task)
   call broadcast_scalar(sst_mask_file_fmt,     master_task)
   call broadcast_array (sst_data_renorm,       master_task)
   call broadcast_scalar(sst_use_obs,  master_task)

!-----------------------------------------------------------------------
!
!  convert data_type to 'monthly-calendar' if input is 'monthly'
!
!-----------------------------------------------------------------------

   if (sst_data_type == 'monthly') &
       sst_data_type = 'monthly-calendar'

!-----------------------------------------------------------------------
!
!  convert interp_type to corresponding integer value.
!
!-----------------------------------------------------------------------

   select case (sst_interp_type)
   case ('nearest')
     sst_interp_order = 1

   case ('linear')
     sst_interp_order = 2

   case ('4point')
     sst_interp_order = 4

   case default
     call exit_POP(sigAbort, &
       'init_sst: Unknown value for sst_interp_type')

   end select

!-----------------------------------------------------------------------
!
!    set values of the SST array (OBS_SST_DATA)
!    depending on the type of the SST data.
!
!-----------------------------------------------------------------------

   select case (sst_data_type)

   case ('none')

      !*** no interior forcing, therefore no interpolation in time
      !*** needed, nor are any new values to be used

      sst_data_next = never
      sst_data_update = never
      sst_interp_freq = 'never'

   case ('annual')

      !*** annual mean climatological interior potential temperature
      !*** (read in from a file) that is constant in time, therefore
      !*** no new values will be needed.

      allocate(OBS_SST_DATA(nx_block,ny_block,max_blocks_clinic,1,1))

      allocate(sst_data_names(1), &
               sst_bndy_loc  (1), &
               sst_bndy_type (1))

      OBS_SST_DATA = c0
      sst_data_names(1) = 'SST'
      sst_bndy_loc  (1) = field_loc_center
      sst_bndy_type (1) = field_type_scalar

      forcing_filename = sst_filename

      sst_data_file = construct_file(sst_file_fmt,                    &
                                   full_name=trim(forcing_filename),  &
                                   record_length=rec_type_dbl,        &
                                   recl_words=nx_global*ny_global)

      call data_set(sst_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)

      sst_data_in = construct_io_field(trim(sst_data_names(1)), &
                             dim1=i_dim, dim2=j_dim,            &
                             field_loc = sst_bndy_loc(1),       &
                             field_type = sst_bndy_type(1),     &
                             d2d_array = OBS_SST_DATA(:,:,:,1,1))

      call data_set (sst_data_file, 'define', sst_data_in)
      call data_set (sst_data_file, 'read',   sst_data_in)
      call data_set (sst_data_file, 'close')
      call destroy_io_field(sst_data_in)
      call destroy_file(sst_data_file)

      if (sst_data_renorm(1) /= c1) &
         OBS_SST_DATA = OBS_SST_DATA*sst_data_renorm(1)

      sst_data_next = never
      sst_data_update = never
      sst_interp_freq = 'never'

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a31,a)') ' Observational SST Annual file read: ', &
                                 trim(forcing_filename)
     endif

   case ('monthly-equal','monthly-calendar')

      !*** monthly mean climatological interior potential temperature.
      !*** All 12 months are read in from a file. interpolation order
      !*** may be specified with namelist input.

      allocate(OBS_SST_DATA(nx_block,ny_block,max_blocks_clinic,1,0:12))

      allocate(sst_data_names(12), &
               sst_bndy_loc  (12), &
               sst_bndy_type (12))

      OBS_SST_DATA = c0
      call find_forcing_times(  sst_data_time,         &
               sst_data_inc,    sst_interp_type,       &
               sst_data_next,   sst_data_time_min_loc, &
               sst_data_update, sst_data_type)

      forcing_filename = sst_filename
      sst_data_file = construct_file(sst_file_fmt,                     &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

      call data_set(sst_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)

      do n=1,12
         write(sst_data_names(n),'(a3,i2.2)') 'SST',n
         sst_bndy_loc (n) = field_loc_center
         sst_bndy_type(n) = field_type_scalar

         sst_data_in = construct_io_field(                      &
                             trim(sst_data_names(n)),           &
                             dim1=i_dim, dim2=j_dim,            &
                             field_loc = sst_bndy_loc(n),       &
                             field_type = sst_bndy_type(n),     &
                             d2d_array = OBS_SST_DATA(:,:,:,1,n))

         call data_set (sst_data_file, 'define', sst_data_in)
         call data_set (sst_data_file, 'read',   sst_data_in)
         call destroy_io_field(sst_data_in)
      enddo

      call data_set (sst_data_file, 'close')
      call destroy_file(sst_data_file)

      if (sst_data_renorm(1) /= c1) &
         OBS_SST_DATA = OBS_SST_DATA*sst_data_renorm(1)

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a32,a)') ' Observational SST Monthly file read: ', &
                                 trim(forcing_filename)
      endif

   case ('n-hour')

      !*** interior potential temperature specified every n-hours,
      !*** where the n-hour increment (sst_data_inc) should
      !*** be specified with namelist input. only as many times as
      !*** are necessary based on the order of the temporal
      !*** interpolation scheme reside in memory at any given time.

      allocate(OBS_SST_DATA(nx_block,ny_block,max_blocks_clinic,1,0:sst_interp_order))

      allocate(sst_data_names(1), &
               sst_bndy_loc  (1), &
               sst_bndy_type (1))

      OBS_SST_DATA = c0
      sst_data_names(1) = 'SST'
      sst_bndy_loc  (1) = field_loc_center
      sst_bndy_type (1) = field_type_scalar

      call find_forcing_times(  sst_data_time,         &
               sst_data_inc,    sst_interp_type,       &
               sst_data_next,   sst_data_time_min_loc, &
               sst_data_update, sst_data_type)

      do n = 1, sst_interp_order
         call get_forcing_filename(forcing_filename, &
                                   sst_filename,     &
                                   sst_data_time(n), &
                                   sst_data_inc)

         sst_data_file = construct_file(sst_file_fmt,                  &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

         call data_set(sst_data_file, 'open_read')

         i_dim = construct_io_dim('i',nx_global)
         j_dim = construct_io_dim('j',ny_global)

         sst_data_in = construct_io_field(                      &
                             trim(sst_data_names(1)),           &
                             dim1=i_dim, dim2=j_dim,            &
                             field_loc = sst_bndy_loc(1),       &
                             field_type = sst_bndy_type(1),     &
                             d2d_array = OBS_SST_DATA(:,:,:,1,n))

         call data_set (sst_data_file, 'define', sst_data_in)
         call data_set (sst_data_file, 'read',   sst_data_in)
         call data_set (sst_data_file, 'close')
         call destroy_io_field(sst_data_in)
         call destroy_file(sst_data_file)

         if (my_task == master_task) then
            write(stdout,blank_fmt)
            write(stdout,'(a31,a)') ' Observational SST n-hour file read: ', &
                                    trim(forcing_filename)
         endif
      enddo

      if (sst_data_renorm(1) /= c1) &
         OBS_SST_DATA = OBS_SST_DATA*sst_data_renorm(1)

   case default

     call exit_POP(sigAbort, &
       'init_sst: Unknown value for sst_data_type')

   end select

   ! Set sst_alpha. For SST sent to coupler will be
   ! sst_alpha*POP_SST + (1-sst_alpha)*OBS_SST_DATA
   ! POP SST is the top level of the first TRACER field.
   sst_alpha = c1
   if (trim(sst_data_type).ne.'none') then
     mask_data_file = construct_file(sst_mask_file_fmt,                   &
                                     full_name = trim(sst_mask_filename), &
                                     record_length = rec_type_dbl,        &
                                     recl_words = nx_global*ny_global)

     call data_set(mask_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)

      sst_data_in = construct_io_field('alpha',                 &
                             dim1=i_dim, dim2=j_dim,            &
                             field_loc = field_loc_center,      &
                             field_type = field_type_scalar,    &
                             d2d_array = sst_alpha)

      call data_set (mask_data_file, 'define', sst_data_in)
      call data_set (mask_data_file, 'read',   sst_data_in)
      call data_set (mask_data_file, 'close')
      call destroy_io_field(sst_data_in)
      call destroy_file(mask_data_file)

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a32,a)') ' Observational SST Mask (alpha) read: ', &
                                 trim(sst_mask_filename)
      endif
      
   end if

!-----------------------------------------------------------------------
!
!  now check interpolation period (sst_interp_freq) to set
!    the time for the next temporal interpolation
!    (sst_interp_next).
!
!  if no interpolation is to be done, set next interpolation time
!    to a large number so the interior PT update test in routine
!    set_surface_forcing will always be false.
!
!  if interpolation is to be done every n-hours, find the first
!    interpolation time greater than the current time.
!
!  if interpolation is to be done every timestep, set next interpolation
!    time to a large negative number so the interior PT update
!    test in routine set_surface_forcing will always be true.
!
!-----------------------------------------------------------------------

   select case (sst_interp_freq)

   case ('never')

     sst_interp_next = never
     sst_interp_last = never
     sst_interp_inc  = c0

   case ('n-hour')
     call find_interp_time(sst_interp_inc, &
                           sst_interp_next)

   case ('every-timestep')

     sst_interp_next = always
     sst_interp_inc  = c0

   case default

     call exit_POP(sigAbort, &
        'init_sst: Unknown value for sst_interp_freq')

   end select

   if(nsteps_total == 0) sst_interp_last = thour00

!-----------------------------------------------------------------------
!
!  echo forcing options to stdout.
!
!-----------------------------------------------------------------------

   sst_data_label = 'Observational SST Forcing'
   if (sst_use_obs .and. my_task == master_task) &
      write(stdout,'(a57)') &
      'Observational SST Restoring enabled'
   call echo_forcing_options(         sst_data_type,   &
                     sst_formulation, sst_data_inc,    &
                     sst_interp_freq, sst_interp_type, &
                     sst_interp_inc,  sst_data_label)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_sst

!***********************************************************************
!BOP
! !IROUTINE: get_sst_data
! !INTERFACE:

 subroutine get_sst_data

! !DESCRIPTION:
!  Determines whether new interior temperature forcing data is required
!  and reads the data if necessary.  Also interpolates data to current
!  time if required.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  check if new data is necessary for interpolation.  if yes, then
!    shuffle time indices in sst_data_time arrays
!    and read in new data if necessary ('n-hour' case).  also
!    increment values of sst_data_time_min_loc,
!    sst_data_next and sst_data_update. note that no new
!    data is necessary for 'analytic' and 'annual' cases.
!
!-----------------------------------------------------------------------

   select case(sst_data_type)

   case ('monthly-equal','monthly-calendar')

      sst_data_label = 'SST Monthly'
      if (thour00 >= sst_data_update) then
         call update_forcing_data(            sst_data_time,   &
               sst_data_time_min_loc, sst_interp_type,         &
               sst_data_next,         sst_data_update,         &
               sst_data_type,         sst_data_inc,            &
               OBS_SST_DATA(:,:,:,:,1:12),sst_data_renorm,     &
               sst_data_label,        sst_data_names,          &
               sst_bndy_loc,          sst_bndy_type,           &
               sst_filename,          sst_file_fmt)
      endif

      if (thour00 >= sst_interp_next .or. nsteps_run==0) then
         call interpolate_forcing(OBS_SST_DATA(:,:,:,:,0),       &
                                  OBS_SST_DATA(:,:,:,:,1:12),    &
                   sst_data_time, sst_interp_type,   &
                   sst_data_time_min_loc,                    &
                   sst_interp_freq, sst_interp_inc,  &
                   sst_interp_next, sst_interp_last, &
                   nsteps_run)

         if (nsteps_run /= 0) sst_interp_next = &
                              sst_interp_next + &
                              sst_interp_inc
      endif

   case('n-hour')

      sst_data_label = 'SST n-hour'
      if (thour00 >= sst_data_update) then
         call update_forcing_data(            sst_data_time,   &
               sst_data_time_min_loc, sst_interp_type, &
               sst_data_next,         sst_data_update, &
               sst_data_type,         sst_data_inc,    &
               OBS_SST_DATA(:,:,:,:,1:sst_interp_order),   &
               sst_data_renorm,                                &
               sst_data_label,        sst_data_names,  &
               sst_bndy_loc,          sst_bndy_type,   &
               sst_filename,          sst_file_fmt)
      endif

      if (thour00 >= sst_interp_next .or. nsteps_run==0) then
         call interpolate_forcing(OBS_SST_DATA(:,:,:,:,0),         &
               OBS_SST_DATA(:,:,:,:,1:sst_interp_order),   &
               sst_data_time,         sst_interp_type, &
               sst_data_time_min_loc, sst_interp_freq, &
               sst_interp_inc,        sst_interp_next, &
               sst_interp_last,       nsteps_run)

         if (nsteps_run /= 0) sst_interp_next = &
                              sst_interp_next + &
                              sst_interp_inc
      endif

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine get_sst_data

 end module forcing_sst

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
