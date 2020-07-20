!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module forcing_uv_interior

!   Time-stamp: <2020-07-20 09:20:47 ivan>

!   Contains routines and variables used for U and V velocity restoring.
!   For use with 1D version of POP (pencil)
!   Code based on forcing_pt_interior.F90 and forcing_s_interior.F90
!   Created by Ivan Lima on Mon Sep 16 2019 13:55:08 -0400

    use kinds_mod
    use domain
    use constants
    use broadcast
    use io
    use forcing_tools
    use time_management
    use prognostic
    use grid
    use exit_mod
    use tavg

    implicit none
    private
    save

    public :: init_uv_interior, &
        get_uv_interior_data,   &
        set_uv_interior

    real (r8), public ::       & ! public for use in restart
        uvel_interior_interp_last, &  ! time last interpolation was done
        vvel_interior_interp_last     ! time last interpolation was done

!-----------------------------------------------------------------------
!
!   internal module variables
!
!-----------------------------------------------------------------------

    integer (int_kind) ::     &
        tavg_UVEL_INTERIOR,   & ! tavg_id for UVEL restoring term
        tavg_VVEL_INTERIOR,   & ! tavg_id for VVEL restoring term
        tavg_UVEL_INTERIOR_2, & ! tavg_id for UVEL restoring term
        tavg_VVEL_INTERIOR_2    ! tavg_id for VVEL restoring term

    real (r8), dimension(:,:,:,:,:), allocatable :: &
        UVEL_INTERIOR_DATA, & ! data to restore UVEL towards
        VVEL_INTERIOR_DATA    ! data to restore VVEL towards

    real (r8), dimension(:,:,:), allocatable :: &
        UV_RESTORE_RTAU  ! inverse restoring timescale for UVEL & VVEL restoring

    integer (int_kind), dimension(:,:,:), allocatable :: &
        UV_RESTORE_MAX_LEVEL ! maximum level for applying restoring

    real (r8), dimension(12) :: &
        uv_interior_data_time    !

    real (r8), dimension(20) :: &
        uv_interior_data_renorm  ! factors to convert to model units

    real (r8) ::                 &
        uv_interior_data_inc,    &! time increment between values of forcing data
        uv_interior_data_next,   &! time to be used for the next value of forcing data that is needed
        uv_interior_data_update, &! time when new forcing value to be added to interpolation set
        uv_interior_interp_inc,  &! time increment between interpolation
        uvel_interior_interp_next, &! time next interpolation will be done
        vvel_interior_interp_next, &! time next interpolation will be done
        uv_interior_restore_tau, &! restoring timescale (non-variable)
        uv_interior_restore_rtau  ! reciprocal of restoring timescale

    integer (int_kind) ::              &
        uv_interior_interp_order,      &! order of temporal interpolation
        uv_interior_data_time_min_loc, &! index of the third dimension of uv_interior_data_time containing the minimum forcing time
        uv_interior_restore_max_level

    character (char_len) ::           &
        uv_interior_data_type,        &! keyword for period of forcing data
        uv_interior_filename,         &! name of file conainting forcing data
        uv_interior_file_fmt,         &! format (bin or netcdf) of forcing file
        uv_interior_interp_freq,      &! keyword for period of temporal interpolation
        uv_interior_interp_type,      &!
        uv_interior_data_label,       &!
        uv_interior_formulation,      &!
        uv_interior_restore_filename, &!
        uv_interior_restore_file_fmt

    character (char_len), dimension(:), allocatable :: &
        uvel_interior_data_names, & ! names for required input data fields
        vvel_interior_data_names

    integer (int_kind), dimension(:), allocatable :: &
        uvel_interior_bndy_loc,  & ! location and field type for
        uvel_interior_bndy_type, & ! ghost cell updates
        vvel_interior_bndy_loc,  &
        vvel_interior_bndy_type

    logical (log_kind) :: &
        uv_interior_variable_restore, &
        uv_interior_surface_restore,  &! Flag to include surface layer when restoring
        uv_interior_nudge              ! Flag to determine if we nudge U & V or just
                                       ! apply forcing directly [default = nudge]

!***********************************************************************

contains

!***********************************************************************

subroutine init_uv_interior

!   Initializes UVEL and VVEL forcing by reading in the 3D fields. Also performs
!   initial book-keeping concerning when new data is needed for the temporal
!   interpolation and when the forcing will need to be updated.

!-----------------------------------------------------------------------
!
!   local variables
!
!-----------------------------------------------------------------------

    integer (int_kind) :: &
        n,                & ! dummy loop index
        nml_error          ! namelist i/o error flag

    character (char_len) :: &
        forcing_filename,   & ! full filename of forcing data file
        long_name            ! long name for input data field

    type (datafile) :: &
        uv_int_data_file  ! data file descriptor for U and V data

    type (io_field_desc) :: &
        uv_data_in  ! io field descriptor for U and V data

    type (io_dim) ::   &
        i_dim, j_dim, & ! dimension descriptors for horiz dims
        k_dim           ! dimension descriptor  for depth

    namelist /forcing_uv_interior_nml/ uv_interior_data_type,           &
        uv_interior_data_inc,         uv_interior_interp_type,         &
        uv_interior_interp_freq,      uv_interior_interp_inc,          &
        uv_interior_restore_tau,      uv_interior_filename,            &
        uv_interior_file_fmt,         uv_interior_restore_max_level,   &
        uv_interior_data_renorm,      uv_interior_formulation,         &
        uv_interior_variable_restore, uv_interior_restore_filename,    &
        uv_interior_restore_file_fmt, uv_interior_surface_restore,     &
        uv_interior_nudge

!-----------------------------------------------------------------------
!
!   Set default values and read input namelist
!
!-----------------------------------------------------------------------

    uv_interior_formulation       = 'restoring'
    uv_interior_data_type         = 'none'
    uv_interior_data_inc          = 1.e20_r8
    uv_interior_interp_type       = 'nearest'
    uv_interior_interp_freq       = 'never'
    uv_interior_interp_inc        = 1.e20_r8
    uv_interior_restore_tau       = 1.e20_r8
    uv_interior_filename          = 'unknown-uv_interior'
    uv_interior_file_fmt          = 'bin'
    uv_interior_restore_max_level = 0
    uv_interior_data_renorm       = c1
    uv_interior_variable_restore  = .false.
    uv_interior_restore_filename  = 'unknown-uv_interior_restore'
    uv_interior_restore_file_fmt  = 'bin'
    uv_interior_surface_restore   = .false.
    uv_interior_nudge             = .true.

    if (my_task == master_task) then
        open (nml_in, file=nml_filename, status='old', iostat=nml_error)
        if (nml_error /= 0) then
            nml_error = -1
        else
            nml_error =  1
        endif
        ! keep reading until find right namelist
        do while (nml_error > 0)
            read(nml_in, nml=forcing_uv_interior_nml,iostat=nml_error)
        end do
        if (nml_error == 0) close(nml_in)
    end if

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
        call exit_POP(sigAbort,'ERROR: reading forcing_uv_interior_nml')
    endif

    call broadcast_scalar(uv_interior_formulation,       master_task)
    call broadcast_scalar(uv_interior_data_type,         master_task)
    call broadcast_scalar(uv_interior_data_inc,          master_task)
    call broadcast_scalar(uv_interior_interp_type,       master_task)
    call broadcast_scalar(uv_interior_interp_freq,       master_task)
    call broadcast_scalar(uv_interior_interp_inc,        master_task)
    call broadcast_scalar(uv_interior_restore_tau,       master_task)
    call broadcast_scalar(uv_interior_filename,          master_task)
    call broadcast_scalar(uv_interior_file_fmt,          master_task)
    call broadcast_scalar(uv_interior_restore_max_level, master_task)
    call broadcast_scalar(uv_interior_variable_restore,  master_task)
    call broadcast_scalar(uv_interior_restore_filename,  master_task)
    call broadcast_scalar(uv_interior_restore_file_fmt,  master_task)
    call broadcast_scalar(uv_interior_surface_restore,   master_task)
    call broadcast_scalar(uv_interior_nudge,             master_task)
    call broadcast_array (uv_interior_data_renorm,       master_task)

!-----------------------------------------------------------------------
!
!   convert data_type to 'monthly-calendar' if input is 'monthly'
!
!-----------------------------------------------------------------------

    if (uv_interior_data_type == 'monthly') &
        uv_interior_data_type = 'monthly-calendar'

!-----------------------------------------------------------------------
!
!   calculate inverse of restoring time scale and convert to seconds.
!
!-----------------------------------------------------------------------

    uv_interior_restore_rtau = c1 / (seconds_in_day * uv_interior_restore_tau)

!-----------------------------------------------------------------------
!
!   convert interp_type to corresponding integer value.
!
!-----------------------------------------------------------------------

    select case (uv_interior_interp_type)

    case ('nearest')
        uv_interior_interp_order = 1

    case ('linear')
        uv_interior_interp_order = 2

    case ('4point')
        uv_interior_interp_order = 4

    case default
        call exit_POP(sigAbort, &
        'init_uv_interior: Unknown value for uv_interior_interp_type')

    end select

!-----------------------------------------------------------------------
!
!   set values of the UVEL and VVEL array (UVEL_INTERIOR_DATA,
!   VVEL_INTERIOR_DATA) depending on the type of the U and V data.
!
!-----------------------------------------------------------------------

    select case (uv_interior_data_type)

    case ('none')

        !*** no interior forcing, therefore no interpolation in time
        !*** needed, nor are any new values to be used

        uv_interior_data_next = never
        uv_interior_data_update = never
        uv_interior_interp_freq = 'never'

    case ('annual')

        !*** annual mean climatological UVEl and VVEL
        !*** (read in from a file) that is constant in time, therefore
        !*** no new values will be needed.

        allocate(UVEL_INTERIOR_DATA(nx_block,ny_block,km,max_blocks_clinic,1), &
            VVEL_INTERIOR_DATA(nx_block,ny_block,km,max_blocks_clinic,1))

        allocate(uvel_interior_data_names(1),   &
            uvel_interior_bndy_loc  (1),        &
            uvel_interior_bndy_type (1),        &
            vvel_interior_data_names(1),        &
            vvel_interior_bndy_loc  (1),        &
            vvel_interior_bndy_type (1))

        UVEL_INTERIOR_DATA = c0
        VVEL_INTERIOR_DATA = c0
        if (uv_interior_nudge) then
            uvel_interior_data_names(1) = 'UVEL'
            vvel_interior_data_names(1) = 'VVEL'
        else
            uvel_interior_data_names(1) = 'UVEL_INTERIOR'
            vvel_interior_data_names(1) = 'VVEL_INTERIOR'
        end if
        uvel_interior_bndy_loc  (1) = field_loc_center
        vvel_interior_bndy_loc  (1) = field_loc_center
        uvel_interior_bndy_type (1) = field_type_scalar
        vvel_interior_bndy_type (1) = field_type_scalar

        forcing_filename = uv_interior_filename

        uv_int_data_file = construct_file(uv_interior_file_fmt,   &
            full_name=trim(forcing_filename),                       &
            record_length=rec_type_dbl,                             &
            recl_words=nx_global*ny_global)

        call data_set(uv_int_data_file, 'open_read')

        i_dim = construct_io_dim('i',nx_global)
        j_dim = construct_io_dim('j',ny_global)
        k_dim = construct_io_dim('k',km)

        uv_data_in = construct_io_field(trim(uvel_interior_data_names(1)),    &
            dim1=i_dim, dim2=j_dim, dim3=k_dim,                                 &
            field_loc = uvel_interior_bndy_loc(1),                              &
            field_type = uvel_interior_bndy_type(1),                            &
            d3d_array = UVEL_INTERIOR_DATA(:,:,:,:,1))

        call data_set (uv_int_data_file, 'define', uv_data_in)
        call data_set (uv_int_data_file, 'read',   uv_data_in)
        call data_set (uv_int_data_file, 'close')
        call destroy_io_field(uv_data_in)
        call destroy_file(uv_int_data_file)

        uv_data_in = construct_io_field(trim(vvel_interior_data_names(1)),    &
            dim1=i_dim, dim2=j_dim, dim3=k_dim,                                 &
            field_loc = vvel_interior_bndy_loc(1),                              &
            field_type = vvel_interior_bndy_type(1),                            &
            d3d_array = VVEL_INTERIOR_DATA(:,:,:,:,1))

        call data_set (uv_int_data_file, 'define', uv_data_in)
        call data_set (uv_int_data_file, 'read',   uv_data_in)
        call data_set (uv_int_data_file, 'close')
        call destroy_io_field(uv_data_in)
        call destroy_file(uv_int_data_file)

        if (uv_interior_data_renorm(1) /= c1) &
            UVEL_INTERIOR_DATA = UVEL_INTERIOR_DATA*uv_interior_data_renorm(1)
            VVEL_INTERIOR_DATA = VVEL_INTERIOR_DATA*uv_interior_data_renorm(1)

        uv_interior_data_next = never
        uv_interior_data_update = never
        uv_interior_interp_freq = 'never'

        if (my_task == master_task) then
            write(stdout,blank_fmt)
            write(stdout,'(a33,a)') ' UVEL and VVEL Annual file read: ', &
                trim(forcing_filename)
        endif

    case ('monthly-equal','monthly-calendar')

        !*** monthly mean climatological UVEL and VVEL.
        !*** All 12 months are read in from a file. interpolation order
        !*** may be specified with namelist input.

        allocate(UVEL_INTERIOR_DATA(nx_block,ny_block,km,max_blocks_clinic,0:12), &
            VVEL_INTERIOR_DATA(nx_block,ny_block,km,max_blocks_clinic,0:12))

        allocate(uvel_interior_data_names(12),  &
            uvel_interior_bndy_loc  (12),       &
            uvel_interior_bndy_type (12),       &
            vvel_interior_data_names(12),       &
            vvel_interior_bndy_loc  (12),       &
            vvel_interior_bndy_type (12))

        UVEL_INTERIOR_DATA = c0
        VVEL_INTERIOR_DATA = c0
        call find_forcing_times(     uv_interior_data_time,         &
            uv_interior_data_inc,    uv_interior_interp_type,       &
            uv_interior_data_next,   uv_interior_data_time_min_loc, &
            uv_interior_data_update, uv_interior_data_type)

        forcing_filename = uv_interior_filename
        uv_int_data_file = construct_file(uv_interior_file_fmt, &
            full_name=trim(forcing_filename),                   &
            record_length=rec_type_dbl,                         &
            recl_words=nx_global*ny_global)

        call data_set(uv_int_data_file, 'open_read')

        i_dim = construct_io_dim('i',nx_global)
        j_dim = construct_io_dim('j',ny_global)
        k_dim = construct_io_dim('k',km)

        do n=1,12
            if (uv_interior_nudge) then
                write(uvel_interior_data_names(n),'(a4,i2.2)') 'UVEL',n
                write(vvel_interior_data_names(n),'(a4,i2.2)') 'VVEL',n
            else
                write(uvel_interior_data_names(n),'(a13,i2.2)') 'UVEL_INTERIOR',n
                write(vvel_interior_data_names(n),'(a13,i2.2)') 'VVEL_INTERIOR',n
            end if
            uvel_interior_bndy_loc (n) = field_loc_center
            vvel_interior_bndy_loc (n) = field_loc_center
            uvel_interior_bndy_type(n) = field_type_scalar
            vvel_interior_bndy_type(n) = field_type_scalar

            uv_data_in = construct_io_field(              &
                trim(uvel_interior_data_names(n)),          &
                dim1=i_dim, dim2=j_dim, dim3=k_dim,         &
                field_loc = uvel_interior_bndy_loc(n),      &
                field_type = uvel_interior_bndy_type(n),    &
                d3d_array = UVEL_INTERIOR_DATA(:,:,:,:,n))

            call data_set (uv_int_data_file, 'define', uv_data_in)
            call data_set (uv_int_data_file, 'read',   uv_data_in)
            call destroy_io_field(uv_data_in)

            uv_data_in = construct_io_field(              &
                trim(vvel_interior_data_names(n)),          &
                dim1=i_dim, dim2=j_dim, dim3=k_dim,         &
                field_loc = vvel_interior_bndy_loc(n),      &
                field_type = vvel_interior_bndy_type(n),    &
                d3d_array = VVEL_INTERIOR_DATA(:,:,:,:,n))

            call data_set (uv_int_data_file, 'define', uv_data_in)
            call data_set (uv_int_data_file, 'read',   uv_data_in)
            call destroy_io_field(uv_data_in)
        enddo

        call data_set (uv_int_data_file, 'close')
        call destroy_file(uv_int_data_file)

        if (uv_interior_data_renorm(1) /= c1) &
            UVEL_INTERIOR_DATA = UVEL_INTERIOR_DATA*uv_interior_data_renorm(1)
            VVEL_INTERIOR_DATA = VVEL_INTERIOR_DATA*uv_interior_data_renorm(1)

        if (my_task == master_task) then
            write(stdout,blank_fmt)
            write(stdout,'(a43,a)') ' Interior UVEL and VVEL Monthly file read: ', &
                trim(forcing_filename)
        endif

    case ('n-hour')

        !*** velocites specified every n-hours,
        !*** where the n-hour increment (uv_interior_data_inc) should
        !*** be specified with namelist input. only as many times as
        !*** are necessary based on the order of the temporal
        !*** interpolation scheme reside in memory at any given time.

        allocate(UVEL_INTERIOR_DATA(nx_block,ny_block,km,max_blocks_clinic,0:uv_interior_interp_order),&
            VVEL_INTERIOR_DATA(nx_block,ny_block,km,max_blocks_clinic,0:uv_interior_interp_order))

        allocate(uvel_interior_data_names(1),   &
            uvel_interior_bndy_loc  (1),        &
            uvel_interior_bndy_type (1),        &
            vvel_interior_data_names(1),        &
            vvel_interior_bndy_loc  (1),        &
            vvel_interior_bndy_type (1))

        UVEL_INTERIOR_DATA = c0
        VVEL_INTERIOR_DATA = c0
        if (uv_interior_nudge) then
            uvel_interior_data_names(1) = 'UVEL'
            vvel_interior_data_names(1) = 'VVEL'
        else
            uvel_interior_data_names(1) = 'UVEL_INTERIOR'
            vvel_interior_data_names(1) = 'VVEL_INTERIOR'
        end if
        uvel_interior_bndy_loc  (1) = field_loc_center
        vvel_interior_bndy_loc  (1) = field_loc_center
        uvel_interior_bndy_type (1) = field_type_scalar
        vvel_interior_bndy_type (1) = field_type_scalar

        call find_forcing_times(     uv_interior_data_time,         &
            uv_interior_data_inc,    uv_interior_interp_type,       &
            uv_interior_data_next,   uv_interior_data_time_min_loc, &
            uv_interior_data_update, uv_interior_data_type)

        do n = 1, uv_interior_interp_order
            call get_forcing_filename(forcing_filename, &
                uv_interior_filename,                   &
                uv_interior_data_time(n),               &
                uv_interior_data_inc)

            uv_int_data_file = construct_file(uv_interior_file_fmt, &
                full_name=trim(forcing_filename),                   &
                record_length=rec_type_dbl,                         &
                recl_words=nx_global*ny_global)

            call data_set(uv_int_data_file, 'open_read')

            i_dim = construct_io_dim('i',nx_global)
            j_dim = construct_io_dim('j',ny_global)
            k_dim = construct_io_dim('k',km)

            uv_data_in = construct_io_field(              &
                trim(uvel_interior_data_names(1)),          &
                dim1=i_dim, dim2=j_dim, dim3=k_dim,         &
                field_loc = uvel_interior_bndy_loc(1),      &
                field_type = uvel_interior_bndy_type(1),    &
                d3d_array = UVEL_INTERIOR_DATA(:,:,:,:,n))

            call data_set (uv_int_data_file, 'define', uv_data_in)
            call data_set (uv_int_data_file, 'read',   uv_data_in)
            call data_set (uv_int_data_file, 'close')
            call destroy_io_field(uv_data_in)
            call destroy_file(uv_int_data_file)

            uv_data_in = construct_io_field(              &
                trim(vvel_interior_data_names(1)),          &
                dim1=i_dim, dim2=j_dim, dim3=k_dim,         &
                field_loc = vvel_interior_bndy_loc(1),      &
                field_type = vvel_interior_bndy_type(1),    &
                d3d_array = VVEL_INTERIOR_DATA(:,:,:,:,n))

            call data_set (uv_int_data_file, 'define', uv_data_in)
            call data_set (uv_int_data_file, 'read',   uv_data_in)
            call data_set (uv_int_data_file, 'close')
            call destroy_io_field(uv_data_in)
            call destroy_file(uv_int_data_file)

            if (my_task == master_task) then
                write(stdout,blank_fmt)
                write(stdout,'(a31,a)') ' Interior UVEL and VVEL n-hour file read: ', &
                    trim(forcing_filename)
            endif
        enddo

        if (uv_interior_data_renorm(1) /= c1) &
            UVEL_INTERIOR_DATA = UVEL_INTERIOR_DATA*uv_interior_data_renorm(1)
            VVEL_INTERIOR_DATA = VVEL_INTERIOR_DATA*uv_interior_data_renorm(1)

    case default

        call exit_POP(sigAbort, &
        'init_uv_interior: Unknown value for uv_interior_data_type')

    end select

!-----------------------------------------------------------------------
!
!   now check interpolation period (uv_interior_interp_freq) to set
!    the time for the next temporal interpolation
!    (uvel_interior_interp_next).
!
!   if no interpolation is to be done, set next interpolation time
!    to a large number so the interior update test in routine
!    set_surface_forcing will always be false.
!
!   if interpolation is to be done every n-hours, find the first
!    interpolation time greater than the current time.
!
!   if interpolation is to be done every timestep, set next interpolation
!    time to a large negative number so the interior update
!    test in routine set_surface_forcing will always be true.
!
!-----------------------------------------------------------------------

    select case (uv_interior_interp_freq)

    case ('never')

        uvel_interior_interp_next = never
        uvel_interior_interp_last = never
        vvel_interior_interp_next = never
        vvel_interior_interp_last = never
        uv_interior_interp_inc  = c0

    case ('n-hour')
        call find_interp_time(uv_interior_interp_inc, &
            uvel_interior_interp_next)
        call find_interp_time(uv_interior_interp_inc, &
            vvel_interior_interp_next)

    case ('every-timestep')

        uvel_interior_interp_next = always
        vvel_interior_interp_next = always
        uv_interior_interp_inc  = c0

    case default

        call exit_POP(sigAbort, &
            'init_uv_interior: Unknown value for uv_interior_interp_freq')

    end select

    if(nsteps_total == 0) uvel_interior_interp_last = thour00
    if(nsteps_total == 0) vvel_interior_interp_last = thour00

!-----------------------------------------------------------------------
!
!   allocate and read in arrays used for variable interior restoring
!   if necessary.
!
!-----------------------------------------------------------------------

    if (uv_interior_variable_restore) then

        allocate(UV_RESTORE_MAX_LEVEL(nx_block,ny_block,max_blocks_clinic), &
            UV_RESTORE_RTAU(nx_block,ny_block,max_blocks_clinic))

        forcing_filename = uv_interior_restore_filename

        uv_int_data_file = construct_file(uv_interior_restore_file_fmt, &
            full_name=trim(forcing_filename),                           &
            record_length=rec_type_dbl,                                 &
            recl_words=nx_global*ny_global)

        call data_set(uv_int_data_file, 'open_read')

        i_dim = construct_io_dim('i',nx_global)
        j_dim = construct_io_dim('j',ny_global)

        uv_data_in = construct_io_field('UV_RESTORE_MAX_LEVEL', &
            dim1=i_dim, dim2=j_dim,                             &
            field_loc = field_loc_center,                       &
            field_type = field_type_scalar,                     &
            d2d_array = UV_RESTORE_RTAU)

        call data_set (uv_int_data_file, 'define', uv_data_in)
        call data_set (uv_int_data_file, 'read',   uv_data_in)
        UV_RESTORE_MAX_LEVEL = nint(UV_RESTORE_RTAU)
        call destroy_io_field(uv_data_in)

        uv_data_in = construct_io_field('UV_RESTORE_RTAU', dim1=i_dim, dim2=j_dim,  &
            field_loc = field_loc_center,                                           &
            field_type = field_type_scalar,                                         &
            d2d_array = UV_RESTORE_RTAU)

        call data_set (uv_int_data_file, 'define', uv_data_in)
        call data_set (uv_int_data_file, 'read',   uv_data_in)
        UV_RESTORE_RTAU = UV_RESTORE_RTAU/seconds_in_day ! convert days to secs
        call destroy_io_field(uv_data_in)
        call data_set (uv_int_data_file, 'close')
        call destroy_file(uv_int_data_file)

    endif

!-----------------------------------------------------------------------
!
!   echo forcing options to stdout.
!
!-----------------------------------------------------------------------

    uv_interior_data_label = 'Interior UVEL and VVEL forcing'
    if (uv_interior_variable_restore .and. my_task == master_task) &
    write(stdout,'(a52)') &
    'Variable interior U and V velocity restoring enabled'
    call echo_forcing_options(uv_interior_data_type,        &
        uv_interior_formulation, uv_interior_data_inc,      &
        uv_interior_interp_freq, uv_interior_interp_type,   &
        uv_interior_interp_inc,  uv_interior_data_label)

!-----------------------------------------------------------------------
!
!   define tavg diagnostic field
!
!-----------------------------------------------------------------------

    call define_tavg_field(tavg_UVEL_INTERIOR, 'UVEL_INTERIOR', 3,  &
        long_name='UVEL Restoring Values',                          &
        units='cm/sec^2', grid_loc='3111',                          &
        coordinates='ULONG ULAT z_t time')

    call define_tavg_field(tavg_VVEL_INTERIOR, 'VVEL_INTERIOR', 3,  &
        long_name='VVEL Restoring Values',                          &
        units='cm/sec^2', grid_loc='3111',                          &
        coordinates='ULONG ULAT z_t time')

    call define_tavg_field(tavg_UVEL_INTERIOR_2, 'UVEL_INTERIOR_2', 3,  &
        long_name='UVEL Restoring Values',                              &
        units='cm/sec^2', grid_loc='3111',                              &
        coordinates='ULONG ULAT z_t time')

    call define_tavg_field(tavg_VVEL_INTERIOR_2, 'VVEL_INTERIOR_2', 3,  &
        long_name='VVEL Restoring Values',                              &
        units='cm/sec^2', grid_loc='3111',                              &
        coordinates='ULONG ULAT z_t time')

!-----------------------------------------------------------------------

end subroutine init_uv_interior

!***********************************************************************

subroutine get_uv_interior_data

!   Determines whether new forcing data is required and reads the data.
!   Interpolates data to current time.

!-----------------------------------------------------------------------

    select case(uv_interior_data_type)

    case ('monthly-equal','monthly-calendar')

        uv_interior_data_label = 'UV_INTERIOR Monthly'
        if (thour00 >= uv_interior_data_update) then
            call update_forcing_data(uv_interior_data_time,                     &
                uv_interior_data_time_min_loc,      uv_interior_interp_type,    &
                uv_interior_data_next,              uv_interior_data_update,    &
                uv_interior_data_type,              uv_interior_data_inc,       &
                UVEL_INTERIOR_DATA(:,:,:,:,1:12),   uv_interior_data_renorm,    &
                uv_interior_data_label,             uvel_interior_data_names,   &
                uvel_interior_bndy_loc,             uvel_interior_bndy_type,    &
                uv_interior_filename,               uv_interior_file_fmt)
        endif
        ! write(stdout,'(a30,12(f12.6,1X))') 'IVAN uvel_interior_data_time ', uv_interior_data_time
        ! write(stdout,'(a32,a20)') 'IVAN uvel_interior_interp_type ', uv_interior_interp_type
        ! write(stdout,*) 'IVAN uvel_interior_data_time_min_loc ', uv_interior_data_time_min_loc
        ! write(stdout,'(a32,a20)') 'IVAN uvel_interior_interp_freq ', uv_interior_interp_freq
        ! write(stdout,*) 'IVAN uvel_interior_interp_inc ', uv_interior_interp_inc
        ! write(stdout,*) 'IVAN uvel_interior_interp_next ', uvel_interior_interp_next
        ! write(stdout,*) 'IVAN uvel_interior_interp_last ', uvel_interior_interp_last
        ! write(stdout,*) 'IVAN uvel_nsteps_run ', nsteps_run
        ! write(stdout,*) 'IVAN -------------------------------------------------------------------'
        if (thour00 >= uvel_interior_interp_next .or. nsteps_run==0) then
            call interpolate_forcing(UVEL_INTERIOR_DATA(:,:,:,:,0), &
                UVEL_INTERIOR_DATA(:,:,:,:,1:12),                   &
                uv_interior_data_time, uv_interior_interp_type,     &
                uv_interior_data_time_min_loc,                      &
                uv_interior_interp_freq, uv_interior_interp_inc,    &
                uvel_interior_interp_next, uvel_interior_interp_last,   &
                nsteps_run)
            if (nsteps_run /= 0) uvel_interior_interp_next =  &
                uvel_interior_interp_next + uv_interior_interp_inc
        endif

        if (thour00 >= uv_interior_data_update) then
            call update_forcing_data(uv_interior_data_time,                     &
                uv_interior_data_time_min_loc,      uv_interior_interp_type,    &
                uv_interior_data_next,              uv_interior_data_update,    &
                uv_interior_data_type,              uv_interior_data_inc,       &
                VVEL_INTERIOR_DATA(:,:,:,:,1:12),   uv_interior_data_renorm,    &
                uv_interior_data_label,             vvel_interior_data_names,   &
                vvel_interior_bndy_loc,             vvel_interior_bndy_type,    &
                uv_interior_filename,               uv_interior_file_fmt)
        endif
        ! write(stdout,'(a30,12(f12.6,1X))') 'IVAN vvel_interior_data_time ', uv_interior_data_time
        ! write(stdout,'(a32,a20)') 'IVAN vvel_interior_interp_type ', uv_interior_interp_type
        ! write(stdout,*) 'IVAN vvel_interior_data_time_min_loc ', uv_interior_data_time_min_loc
        ! write(stdout,'(a32,a20)') 'IVAN vvel_interior_interp_freq ', uv_interior_interp_freq
        ! write(stdout,*) 'IVAN vvel_interior_interp_inc ', uv_interior_interp_inc
        ! write(stdout,*) 'IVAN vvel_interior_interp_next ', vvel_interior_interp_next
        ! write(stdout,*) 'IVAN vvel_interior_interp_last ', vvel_interior_interp_last
        ! write(stdout,*) 'IVAN vvel_nsteps_run ', nsteps_run
        ! write(stdout,*) 'IVAN ###################################################################'
        if (thour00 >= vvel_interior_interp_next .or. nsteps_run==0) then
            call interpolate_forcing(VVEL_INTERIOR_DATA(:,:,:,:,0), &
                VVEL_INTERIOR_DATA(:,:,:,:,1:12),                   &
                uv_interior_data_time, uv_interior_interp_type,     &
                uv_interior_data_time_min_loc,                      &
                uv_interior_interp_freq, uv_interior_interp_inc,    &
                vvel_interior_interp_next, vvel_interior_interp_last,   &
                nsteps_run)
            if (nsteps_run /= 0) vvel_interior_interp_next =  &
                vvel_interior_interp_next + uv_interior_interp_inc
        endif

    case('n-hour')

        uv_interior_data_label = 'UV_INTERIOR n-hour'
        if (thour00 >= uv_interior_data_update) then
            call update_forcing_data(uv_interior_data_time,             &
                uv_interior_data_time_min_loc, uv_interior_interp_type, &
                uv_interior_data_next,         uv_interior_data_update, &
                uv_interior_data_type,         uv_interior_data_inc,    &
                UVEL_INTERIOR_DATA(:,:,:,:,1:uv_interior_interp_order), &
                uv_interior_data_renorm,                                &
                uv_interior_data_label,        uvel_interior_data_names,&
                uvel_interior_bndy_loc,        uvel_interior_bndy_type, &
                uv_interior_filename,          uv_interior_file_fmt)
        endif
        if (thour00 >= uvel_interior_interp_next .or. nsteps_run==0) then
            call interpolate_forcing(UVEL_INTERIOR_DATA(:,:,:,:,0),     &
                UVEL_INTERIOR_DATA(:,:,:,:,1:uv_interior_interp_order), &
                uv_interior_data_time,         uv_interior_interp_type, &
                uv_interior_data_time_min_loc, uv_interior_interp_freq, &
                uv_interior_interp_inc,        uvel_interior_interp_next, &
                uvel_interior_interp_last,     nsteps_run)
            if (nsteps_run /= 0) uvel_interior_interp_next =    &
                uvel_interior_interp_next + uv_interior_interp_inc
        endif


        if (thour00 >= uv_interior_data_update) then
            call update_forcing_data(uv_interior_data_time,             &
                uv_interior_data_time_min_loc, uv_interior_interp_type, &
                uv_interior_data_next,         uv_interior_data_update, &
                uv_interior_data_type,         uv_interior_data_inc,    &
                VVEL_INTERIOR_DATA(:,:,:,:,1:uv_interior_interp_order), &
                uv_interior_data_renorm,                                &
                uv_interior_data_label,        vvel_interior_data_names,&
                vvel_interior_bndy_loc,        vvel_interior_bndy_type, &
                uv_interior_filename,          uv_interior_file_fmt)
        endif
        if (thour00 >= uvel_interior_interp_next .or. nsteps_run==0) then
            call interpolate_forcing(VVEL_INTERIOR_DATA(:,:,:,:,0),     &
                VVEL_INTERIOR_DATA(:,:,:,:,1:uv_interior_interp_order), &
                uv_interior_data_time,         uv_interior_interp_type, &
                uv_interior_data_time_min_loc, uv_interior_interp_freq, &
                uv_interior_interp_inc,        vvel_interior_interp_next, &
                vvel_interior_interp_last,     nsteps_run)
            if (nsteps_run /= 0) vvel_interior_interp_next =    &
                vvel_interior_interp_next + uv_interior_interp_inc
        endif

    end select

!-----------------------------------------------------------------------

end subroutine get_uv_interior_data

!***********************************************************************

subroutine set_uv_interior(k,this_block,DU_INTERIOR,DV_INTERIOR)

!   Computes UVEL and VVEL restoring terms using updated & interpolated data.

!   INPUT PARAMETERS:

    integer (int_kind), intent(in) :: &
        k     ! vertical level index

    type (block), intent(in) :: &
        this_block   ! block information for this block

!   OUTPUT PARAMETERS:

!    real (r8), dimension(nx_block,ny_block), intent(inout) :: &
!        U_SOURCE, &  ! source term for UVEL
!        V_SOURCE     ! source term for VVEL
!                     ! add restoring terms to these arrays

    real (r8), dimension(nx_block,ny_block), intent(out) :: &
        DU_INTERIOR, &  ! UVEL restoring term for this block and level
        DV_INTERIOR     ! VVEL restoring term for this block and level

!-----------------------------------------------------------------------
!
!   local variables
!
!-----------------------------------------------------------------------

    integer (int_kind) ::  &
        bid, & ! local block address for this block
        now    ! index for interpolated data

!-----------------------------------------------------------------------
!
!   do restoring if required
!
!-----------------------------------------------------------------------

    if ((k > 1 .or. uv_interior_surface_restore).and.                     &
        uv_interior_data_type.ne.'none') then

!   set index for location of interpolated data.

        select case(uv_interior_data_type)

        case('annual')
            now = 1

        case ('monthly-equal','monthly-calendar')
            now = 0

        case('n-hour')
            now = 0

        end select

!     compute restoring

        bid = this_block%local_id

        if (uv_interior_nudge) then
            if (uv_interior_variable_restore) then
                DU_INTERIOR = UV_RESTORE_RTAU(:,:,bid) *        &
                    merge((UVEL_INTERIOR_DATA(:,:,k,bid,now) -  &
                    UVEL(:,:,k,newtime,bid)),                   &
                    c0, k <= UV_RESTORE_MAX_LEVEL(:,:,bid))

                DV_INTERIOR = UV_RESTORE_RTAU(:,:,bid) *        &
                    merge((VVEL_INTERIOR_DATA(:,:,k,bid,now) -  &
                    VVEL(:,:,k,newtime,bid)),                   &
                    c0, k <= UV_RESTORE_MAX_LEVEL(:,:,bid))
            else
                if (k <= uv_interior_restore_max_level) then
                    DU_INTERIOR = uv_interior_restore_rtau *    &
                        (UVEL_INTERIOR_DATA(:,:,k,bid,now) -    &
                        UVEL(:,:,k,newtime,bid))

                    DV_INTERIOR = uv_interior_restore_rtau *    &
                        (VVEL_INTERIOR_DATA(:,:,k,bid,now) -    &
                        VVEL(:,:,k,newtime,bid))
                else
                    DU_INTERIOR = c0
                    DV_INTERIOR = c0
                endif
            endif
        else
            DU_INTERIOR = UVEL_INTERIOR_DATA(:,:,k,bid,now)
            DV_INTERIOR = VVEL_INTERIOR_DATA(:,:,k,bid,now)
        end if

        call accumulate_tavg_field(DU_INTERIOR, tavg_UVEL_INTERIOR, bid, k)
        call accumulate_tavg_field(DV_INTERIOR, tavg_VVEL_INTERIOR, bid, k)
        call accumulate_tavg_field(DU_INTERIOR, tavg_UVEL_INTERIOR_2, bid, k)
        call accumulate_tavg_field(DV_INTERIOR, tavg_VVEL_INTERIOR_2, bid, k)

!        ! to check if forcing is read correctly
!        call accumulate_tavg_field(UVEL_INTERIOR_DATA(:,:,k,bid,12), tavg_UVEL_INTERIOR, bid, k)
!        call accumulate_tavg_field(VVEL_INTERIOR_DATA(:,:,k,bid,12), tavg_VVEL_INTERIOR, bid, k)
!        ! to check interpolated reference field
!        call accumulate_tavg_field(UVEL_INTERIOR_DATA(:,:,k,bid,now), tavg_UVEL_INTERIOR, bid, k)
!        call accumulate_tavg_field(VVEL_INTERIOR_DATA(:,:,k,bid,now), tavg_VVEL_INTERIOR, bid, k)

!       add restoring term to other source terms

!        U_SOURCE = U_SOURCE + DU_INTERIOR
!        V_SOURCE = V_SOURCE + DV_INTERIOR

    endif ! k=1

!-----------------------------------------------------------------------

end subroutine set_uv_interior

!***********************************************************************

! subroutine uv_restore(UVEL, VVEL, this_block)

! !-----------------------------------------------------------------------
! !
! !   INPUT PARAMETERS
! !
! !-----------------------------------------------------------------------

!     type (block), intent(in) :: this_block ! block info for the current block

! !-----------------------------------------------------------------------
! !
! !   INPUT/OUTPUT PARAMETERS
! !
! !-----------------------------------------------------------------------

!     real(r8), dimension(nx_block, ny_block, km), intent(inout) :: UVEL, VVEL

! !-----------------------------------------------------------------------
! !
! !  local variables
! !
! !-----------------------------------------------------------------------

!     integer (int_kind) :: k     ! depth level index

!     real(r8), dimension(nx_block, ny_block) :: URESTORE, VRESTORE

! !-----------------------------------------------------------------------

!     do k = 1,km
!         call set_uv_interior(k,this_block,URESTORE,VRESTORE)
!     enddo

! !-----------------------------------------------------------------------

! end subroutine uv_restore

!***********************************************************************

end module forcing_uv_interior

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
