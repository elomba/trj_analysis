subroutine gpu_and_header(startEvent,stopEvent)
    use cudafor
    use mod_input, only : nthread
    use mod_common, only : shmsize, maxthread
    implicit None
    type(cudaEvent), intent(inout) :: startEvent, stopEvent
    integer :: istat
    ! Get CUDA properties from device 0 (can be set from environmente variable CUDA_VISIBLE_DEVICES)

    istat = cudaSetDevice(0)
    if (istat == 0) then
        istat = cudaGetDeviceProperties(gpu_properties, 0)
    else
        write(*,"('*** Unrecoverable error: no GPU available !!')")
        stop
    end if
    shmsize = gpu_properties%sharedMemPerBlock
    maxthread = gpu_properties%maxThreadsPerBlock
    ! Initialize CUDA timing
    istat = cudaEventCreate(startEvent)
    istat = cudaEventCreate(stopEvent)
    ! Print program header
    write(*,"(/80('*')/'*',78(' '),'*')")
    write(*,"('*    Program trj_analysis: analyzing LAMMPS trajectory in NETCDF format',t80,'*')")
    write(*,"('*',t80,'*')")
    write(*,"('*    Using GPU with CUDA nvfortran/nvcc >= 11.6',t80,'*')")
    write(*,"('*',t80,'*')")
    write(*,"('*    Version 0.2.30 October 2024',,t80,'*')")
    write(*,"('*',78(' '),'*'/80('*')/)")
    call printDevPropShort(gpu_propertoies, 0)
     ! Check that maximum number of threads is not surpassed
    if (nthread > maxthread/8) then
        nthread = maxthread/8
        write(*,'("** Warning: number of threads reset to",I3)')nthread
    endif
end subroutine gpu_and_header


subroutine basic_init(use_cell,run_clusters,run_dyn,nmol)
    implicit none
    logical, intent(in) :: use_cell,run_clusters,run_dyn
    integer, intent(in) :: nmol
    if (use_cell) call cells_init_pre_nc_read(nmol)
    if (run_clusters) call clusters_init(nmol)
    if (run_dyn) call dyn_init()
end subroutine basic_init

subroutine init_modules(use_cell,run_rdf,run_sq,run_clusters)
    implicit none
    logical, intent(IN) :: use_cell,run_rdf,run_sq,run_clusters 
     if (use_cell) call cells_init_post_nc_read() 
     if (run_rdf) call RDF_init(nsp)
     if (run_sq) call sq_init(nmol, nsp, nbcuda)
     if (run_clusters) call clusters_sq_init()
end subroutine init_modules

subroutine clean_memory(run_sq,run_rdf,run_clusters,use_cell,run_dyn,confined)
    implicit none
    logical, intent(IN) :: run_sq,run_rdf,run_clusters,use_cell,run_dyn,confined 
    if (run_sq) call sq_clear()
    if (run_rdf) call rdf_clear()
    if (run_clusters) call clusters_clear()
    if (use_cell) call cells_clear()
    if (run_dyn) call dyn_clear()
    if (confined) call prof_clear()
    if (run_thermo) call thermo_clear()
    call common_clear()
    call log_clear()
    call input_clear()
end subroutine clean_memory

subroutine print_results(run_sq,  run_rdf, run_dyn, run_clusters, run_thermo, nsp, lsmax, nmol, nqmin)
    implicit none
    integer, intent(in) :: nsp, lsmax, nmol, nqmin
    logical, intent(in) :: run_sq,  run_rdf, run_dyn, run_clusters, run_thermo 

    ! Print out S(Q)'s
    if (run_sq) then
        call sq_transfer_gpu_cpu()
        call printSQ(Nmol)
    end if

    do i = 1, nsp
        write (*, '(" ** ",i6," atoms of type ",i2)') ntype(i), i
    end do

    ! Print partial pair distribution functions
    if (run_rdf) then
         call printrdf(rcl, lsmax)
    end if
    if (run_dyn) then
        call print_rtcor()
    endif
    if (run_clusters) then
        call print_clusinfo(nqmin, Nmol)
        if (run_thermo) then
            call printPotEngCl()
            call printPotEngClCl()
        end if
    end if
    
end subroutine print_results

subroutine reformat_input_conf(io,final_conf,current_conf,ntypes,nsp)
    implicit none
    integer, intent(inout) :: final_conf
    integer, intent(in) :: io, ntypes, nsp, current_conf
    if (io<0) then
            final_conf=current_conf-1
    else
        ! Pre configuration analysis data transformations, corrections and format
        if (ntypes == nsp) then
            ! All species selected
            call trans_ncdfinput()
        else
            ! Only some species selected from the configuration
            call select_ncdfinput()
        endif
    endif
end subroutine reformat_input_conf