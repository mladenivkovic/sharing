program eval_galaxies

  use,intrinsic :: iso_c_binding
  use omp_lib
  implicit none
  include 'fftw3.f03'

  !-------------------------------------------
  ! Global variables and parameters
  !-------------------------------------------

  integer, parameter  :: dp = kind(1.0d0) 

  real(dp), parameter :: pi = 3.14159265359
  real(dp), parameter :: twopicube = (2.d0*pi)**3
  real(dp)            :: nccube
  real(dp), parameter :: Mpc = 3.086e+24
  real(dp), parameter :: boxlen = 1.d0
  real(dp)            :: aexp             ! expansion factor
  real(dp)            :: H0               ! Hubble factor [km s-1 Mpc-1]
  real(dp)            :: unit_l           ! length of box [Mpc]
  real(dp)            :: dk, dx           ! physical values of cell size in real/k space
  integer, parameter  :: nsamples = 100   ! number of samples for output; Number of bins for histograms
  integer             :: nhalos, ngalaxies, ncpu
  character(len=12)   :: srcdir           ! output_XXXXX 

  ! read-in data arrays
  real(dp), dimension(:), allocatable :: x, y, z, stellar_mass
  integer, dimension(:), allocatable  :: halo_id, clmp_id

  ! computed global arrays
  integer  :: nc
  real(dp) :: cellvolume
  real(dp), dimension(:,:,:), allocatable    :: density_field
  complex(dp), dimension(:,:,:), allocatable :: Pk_field
  real(dp), dimension(:), allocatable        :: Pk, correlation
  real(dp), dimension(:), allocatable        :: Pk_counts, correlation_counts
  real(dp), dimension(:), allocatable        :: distances, distances_k

  ! FFTW stuff
  integer*8                                  :: plan_forward, plan_backward


  !-----------------------------
  ! Other variables
  !-----------------------------
  integer :: i,j,k, void, ix, iy, iz, ik, ig, c
  real(dp):: mean_density, kmax, lmax, d
  character(len=80) :: outfname
  character(len=4)  :: add_name
  character(len=3), dimension(1:2) :: cases=(/'sub','all'/)






  write(*,*) "============================"
  write(*,*) "get_correlation.f90 started."
  write(*,*) "============================"


  !----------------------------------
  write(*,*) "Reading data."
  !----------------------------------

  nc = 128
  call read_cmdlineargs() ! read command line arguments
  call read_info()        ! read data from info files
  call read_data()        ! read actual data and nhalos, ngalaxies


  ! Some definitions
  lmax       = 0.5d0*unit_l
  dx         = unit_l/nc
  dk         = 1.d0/unit_l
  kmax       = (nc/2) * dk  ! choose max k to be at half "box size" so you have full spheric samples only
  cellvolume = dx**3
  nccube     = nc**3


  write(*,*)                "Working with parameters:"
  write(*,'(A16,A25)')      "  srcdir:    ", srcdir
  write(*,'(A16,I25)')      "  ncpu:      ", ncpu
  write(*,'(A16,I25)')      "  nc:        ", nc
  write(*,'(A16,I25)')      "  nhaloes:   ", nhalos
  write(*,'(A16,I25)')      "  ngalaxies: ", ngalaxies
  write(*,'(A16,E25.8,A9)') "  lmax:      ", lmax, "[Mpc]"
  write(*,'(A16,E25.8,A9)') "  kmax:      ", kmax, "[Mpc-1]"




  !-----------------------------------------------
  write(*,*) "Calculating Stellar Mass Function"
  !-----------------------------------------------
  call smf()




  allocate(distances(0:nsamples));                distances=0
  allocate(distances_k(0:nsamples));              distances_k=0;
  allocate(Pk_counts(0:nsamples));                Pk_counts=0;
  allocate(correlation_counts(0:nsamples));       correlation_counts=0;
  allocate(Pk(1:nsamples));                       Pk=0;
  allocate(correlation(1:nsamples));              correlation=0;

  allocate(density_field(1:nc, 1:nc, 1:nc))
  allocate(Pk_field(1:nc/2+1, 1:nc, 1:nc))




  !----------------------------------
  ! Initiate fftw stuff
  !----------------------------------
  void = fftw_init_threads()
  if (void==0) then
    write(*,*) "Error in fftw_init_threads, quitting"
    stop
  endif
  call fftw_plan_with_nthreads(omp_get_max_threads())
  call dfftw_plan_dft_r2c_3d(plan_forward,  nc, nc, nc, density_field, Pk_field, FFTW_ESTIMATE)
  call dfftw_plan_dft_c2r_3d(plan_backward, nc, nc, nc, Pk_field, density_field, FFTW_ESTIMATE)





  !-------------------------------------------
  ! Do everything twice:
  ! If c = 1, include orphans
  ! If c = 2, don't include orphan galaxies
  !-------------------------------------------

  do c = 1, 2
    write(*,*)
    add_name = cases(c)

    !---------------------------------------------------------
    write(*,*) "Computing density field for case ", add_name
    !---------------------------------------------------------

    call get_density_field(c)
    ! call write_galaxy_images() ! call before changing density

    mean_density = 0
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,k)
      !$OMP DO COLLAPSE(3) REDUCTION (+:mean_density)
        do i = 1, nc
          do j = 1, nc
            do k = 1, nc
              mean_density = mean_density + density_field(i,j,k)
            enddo
          enddo
        enddo
      !$OMP END DO

      !$OMP SINGLE
        mean_density = mean_density/nccube
      !$OMP END SINGLE

      !$OMP FLUSH(mean_density)

      !$OMP DO COLLAPSE(3)
        do i=1, nc
          do j = 1, nc
            do k = 1, nc
              density_field(i,j,k) = density_field(i,j,k)/mean_density - 1
            enddo
          enddo
        enddo
      !$OMP END DO

    !$OMP END PARALLEL




    !-----------------------------------------------------------------------
    write(*,*) "Computing and histogramming P(k) FFT for case ", add_name
    !-----------------------------------------------------------------------

    call dfftw_execute_dft_r2c(plan_forward, density_field, Pk_field)

    !$OMP PARALLEL PRIVATE(i,j,k,ix,iy,iz,ik,ig,d)
      !$OMP DO COLLAPSE(3)
        do i=1, nc/2+1
          do j=1, nc
            do k=1, nc
              Pk_field(i,j,k) = Pk_field(i,j,k)*conjg(Pk_field(i,j,k))/twopicube
            enddo
          enddo
        enddo
      !$OMP ENDDO

      !$OMP DO
        ! Get bin distances, reset values. Use logarithmic bins
        do i = 1, nsamples
          distances_k(i) = dk * (kmax/dk)**(real(i,dp)/real(nsamples,dp))
          Pk_counts(i)=0;           
          Pk(i) = 0;
        enddo
      !$OMP ENDDO

      !$OMP DO COLLAPSE(3)
        do i = 1, nc/2+1
          do j = 1, nc
            do k = 1, nc
              ix = i-1
              if (j-1<nc/2+1) then; iy = j-1; else; iy = -nc+j-1; endif;
              if (k-1<nc/2+1) then; iz = k-1; else; iz = -nc+k-1; endif;
              d = sqrt(real(ix**2+iy**2+iz**2,kind=dp))*dk
              if (d<=kmax .and. d>0.d0) then
                ig = int( log(d/dk)*nsamples/log(kmax/dk) ) -1
                ig = max(ig, 1)
                do ik=ig, nsamples
                  if (d<=distances_k(ik)) exit
                enddo
                !$OMP ATOMIC
                Pk(ik) = Pk(ik)+real(Pk_field(i,j,k))
                !$OMP ATOMIC
                Pk_counts(ik) = Pk_counts(ik) + 1
              endif
            enddo
          enddo
        enddo
      !$OMP END DO

    !$OMP END PARALLEL





    !-----------------------------------------------------------------------
    write(*,*) "Computing and histogramming xi(r) FFT for case ", add_name
    !-----------------------------------------------------------------------

    call dfftw_execute_dft_c2r(plan_backward, Pk_field, density_field)
    
    !$OMP PARALLEL PRIVATE(i,j,k,ix,iy,iz,ik,ig,d)

      !$OMP DO COLLAPSE(3)
        do i=1, nc
          do j=1, nc
            do k=1, nc
              density_field(i,j,k) = density_field(i,j,k)/nccube
            enddo
          enddo
        enddo
      !$OMP END DO

      !$OMP DO
        ! Get bin distances, reset values
        do i = 1, nsamples
          distances(i) = i*lmax/nsamples
          correlation_counts(i)=0;  
          correlation(i)=0;
        enddo
      !$OMP END DO

      !$OMP DO COLLAPSE(3)
        do i = 1, nc
          do j = 1, nc
            do k = 1,nc
              if (i-1<nc/2+1) then; ix = i-1; else; ix = -nc+i-1; endif
              if (j-1<nc/2+1) then; iy = j-1; else; iy = -nc+j-1; endif
              if (k-1<nc/2+1) then; iz = k-1; else; iz = -nc+k-1; endif
              d = sqrt(real(ix**2+iy**2+iz**2,kind=dp))*dx
              if (d<=lmax) then
                ig = int(d*nsamples/lmax) - 1
                ig = max(ig, 1)
                do ik=ig, nsamples
                  if (d<=distances(ik)) exit
                enddo
                !$OMP ATOMIC
                correlation(ik) = correlation(ik)+density_field(i,j,k)
                !$OMP ATOMIC
                correlation_counts(ik) = correlation_counts(ik) + 1
              endif
            enddo
          enddo
        enddo
      !$OMP END DO




      !---------------------------------
      ! Write results to file
      !---------------------------------

      !$OMP SINGLE
        outfname=TRIM(TRIM(srcdir)//"/Pk_"//TRIM(add_name)//".txt")
        open(unit=666, file=outfname, form='formatted')
        write(666,*) "Results generated by eval_galaxies.f90"
        write(666,'(2A22,x)') "k", "P(k)"
        do i = 1, nsamples
          if (Pk_counts(i)>0) then
            ! write center of bins, multiply by 2*pi to get values for FFT convention without 2pi in
            ! exponent. => *0.5 for center value and *2 from *2*pi cancel out
            write(666, '(2E22.10,x)') (distances_k(i-1)+distances_k(i))*pi, Pk(i)/Pk_counts(i)
          endif
        enddo
        write(*,*) "Finished P(k), written results to "//TRIM(outfname)
        close(666)
      !$OMP END SINGLE

      !$OMP SINGLE
        outfname=TRIM(TRIM(srcdir)//"/correlation_"//TRIM(add_name)//".txt")
        open(unit=667, file=outfname, form='formatted')
        write(667,*) "Results generated by eval_galaxies.f90"
        write(667,'(2A22,x)') "r", "correlation(r)"
        do i = 1, nsamples
          if (correlation_counts(i) > 0 ) then
            write(667, '(2E22.10,x)') 0.5*(distances(i-1)+distances(i)), correlation(i)/correlation_counts(i)
          endif
        enddo
        write(*,*) "Finished xi(r), written results to "//TRIM(outfname)
        close(667)
      !$OMP END SINGLE
    !$OMP END PARALLEL

  enddo

  


  write(*,*) "eval_galaxies.f90 finished."

contains

  !====================================
  subroutine get_density_field(which)
  !====================================
    !------------------------------------------------------------
    ! Computes the density field and overdensity field delta
    ! which = 1: main + satellites
    ! which = 2: main + satellites + orphans
    !------------------------------------------------------------

    implicit none
    integer, intent(in)             :: which
    integer                         :: i, j, k, gal
    logical                         :: found_it

    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,k,gal, found_it)

      !$OMP DO COLLAPSE(3)
        ! reset density
        do i=1, nc
          do j=1, nc
            do k=1, nc
              density_field(i,j,k) = 0.d0
            enddo
          enddo
        enddo
      !$OMP END DO

      !$OMP DO
        do gal=1, ngalaxies
          if (which==2 .or. clmp_id(gal)>0) then
            i = int(x(gal)/boxlen*nc)+1
            j = int(y(gal)/boxlen*nc)+1
            k = int(z(gal)/boxlen*nc)+1
            !$OMP ATOMIC
            density_field(i,j,k) = density_field(i,j,k)+stellar_mass(gal)
          endif
        enddo
      !$OMP END DO

      !$OMP DO COLLAPSE(3)
        do i=1, nc
          do j=1, nc
            do k=1, nc
              density_field(i,j,k) = density_field(i,j,k)/cellvolume
            enddo
          enddo
        enddo
      !$OMP END DO
    !$OMP END PARALLEL
  end subroutine get_density_field





  !========================================
  subroutine read_cmdlineargs()
  !========================================
    !-----------------------
    ! reads cmd line args
    !-----------------------
    character(len=12) :: arg
    integer :: ncells_temp

    do i = 1, iargc()
      call getarg(i, arg)
      if (i==1) then
        read(arg, '(A12)') srcdir
      else if (i==2) then
        read(arg, *) ncells_temp
        if (ncells_temp>=0) nc = ncells_temp
      endif
    end do
  end subroutine read_cmdlineargs



  
  !==================================
  subroutine read_data()
  !==================================
    !------------------------------------------
    ! reads in file sizes so you can allocate
    ! arrays of appropriate sizes.
    ! Essentially, get nhaloes and ngalaxies
    !------------------------------------------

    use omp_lib
    character(len=5)            :: outputnr, cpunr
    character(len=80)           :: fname
    character(len=125)          :: junk
    integer, dimension(0:ncpu)  :: counts_g, counts_h
    integer                     :: i, ind, id, io, ifile
    logical                     :: file_exists

    counts_g = 0; counts_g(0)=1;
    counts_h = 0; counts_h(0)=1;
    outputnr = srcdir(8:12)

    !$OMP PARALLEL PRIVATE(i, id, io, ifile, file_exists, fname, cpunr, ind, junk)
      id = omp_get_thread_num()

      !----------------------------------
      ! Read file sizes
      !----------------------------------

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/galaxies_'//outputnr//'.txt'//cpunr)
          inquire(file=fname, exist=file_exists)
          if (file_exists) then
            open(unit=600+id, file=fname)
            read(600+id, *) ! skip header
            do 
              read(600+id,*, iostat=io)
              if (io/=0) exit
              counts_g(ifile) = counts_g(ifile)+1
            enddo
            close(600+id)
          else
            write(*,*) "Didn't find file ", fname
          endif
        enddo
      !$OMP ENDDO

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/halo_'//outputnr//'.txt'//cpunr)
          inquire(file=fname, exist=file_exists)
          if (file_exists) then
            open(unit=600+id, file=fname)
            read(600+id,*) ! skip header
            do 
              read(600+id,*, iostat=io)
              if (io/=0) exit
              counts_h(ifile) = counts_h(ifile)+1
            enddo
            close(600+id)
          else
            write(*,*) "Didn't find file ", fname
          endif
        enddo
      !$OMP ENDDO

    
      !$OMP DO ORDERED
        do i=0, ncpu-1
          !$OMP ORDERED
            counts_g(i+1) = counts_g(i+1)+counts_g(i)
            counts_h(i+1) = counts_h(i+1)+counts_h(i)
          !$OMP END ORDERED
        enddo
      !$OMP ENDDO



      !----------------------------------
      ! Allocate Arrays
      !----------------------------------

      !$OMP SINGLE
        nhalos = counts_h(ncpu)-1
        ngalaxies = counts_g(ncpu)-1
        allocate(x(1:ngalaxies));             x=0;
        allocate(y(1:ngalaxies));             y=0;
        allocate(z(1:ngalaxies));             z=0;
        allocate(stellar_mass(1:ngalaxies));  stellar_mass=0;
        allocate(clmp_id(1:ngalaxies));       clmp_id=0;
        allocate(halo_id(1:nhalos));          halo_id=0;
      !$OMP END SINGLE
      !$OMP BARRIER



      !----------------------------------
      ! Read galaxy data
      !----------------------------------

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/galaxies_'//outputnr//'.txt'//cpunr)
          open(unit=600+id, file=fname)
          read(600+id,*) ! skip header
          do ind=counts_g(ifile-1), counts_g(ifile)-1
            read(600+id,'(I20,x,4(E20.12,x))') clmp_id(ind), stellar_mass(ind), x(ind), y(ind), z(ind)
          enddo
          close(600+id)
        enddo
      !$OMP END DO


      !----------------------------------
      ! Read halo data
      !----------------------------------

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/halo_'//outputnr//'.txt'//cpunr)
          open(unit=600+id, file=fname)
          read(600+id,*) ! skip header
          do ind=counts_h(ifile-1), counts_h(ifile)-1
            read(600+id,'(I10,x,A125)') halo_id(ind), junk
          enddo
          close(600+id)
        enddo
      !$OMP END DO

    !$OMP END PARALLEL
  end subroutine read_data





  !========================================
  subroutine read_info()
  !========================================
    !-------------------------------
    ! Read in galaxy data.
    !-------------------------------
    implicit none
    character(len=5)  :: outputnr
    character(len=13) :: junk
    character(len=80) :: fname
    logical :: file_exists
    integer :: i

    outputnr = srcdir(8:12)
    fname=TRIM(srcdir//'/info_'//outputnr//'.txt')
    
    inquire(file=fname, exist=file_exists)
    if (file_exists) then

      open(unit=666, file=fname)
      read(666, '(A13,I11)') junk, ncpu
      do i = 1, 8
        ! skip lines
        read(666,*)
      enddo

      read(666, '(A13,E23.15)') junk, aexp
      read(666, '(A13,E23.15)') junk, H0

      do i = 1, 4
        ! skip lines
        read(666,*)
      enddo

      read(666, '(A13,E23.15)') junk, unit_l
      unit_l = unit_l / Mpc

      close(666)

    else
      write(*,*) "Didn't find file", fname
      stop
    endif
  end subroutine read_info



  !========================================
  subroutine smf()
  !========================================
    !-----------------------------------------------
    ! Computes and writes stellar mass function
    !-----------------------------------------------

    ! log10(M_galaxy) boundaries:
    real(dp), parameter                 :: minmass = 5.5d0
    real(dp), parameter                 :: maxmass = 13.0d0
    real(dp), allocatable, dimension(:) :: smf_all, smf_sub, smf_main, logmass
    real(dp)                            :: dm, lm, vol, f
    integer                             :: i,j,guess
    character(len=80)                   :: fname


    allocate(smf_all(1:nsamples));  smf_all =0;
    allocate(smf_sub(1:nsamples));  smf_sub =0;
    allocate(smf_main(1:nsamples)); smf_main=0;
    allocate(logmass(0:nsamples));  logmass =0;


    dm = (maxmass-minmass)/nsamples
    vol = unit_l**3


    !$OMP PARALLEL PRIVATE(i,j,guess,lm,f)
      !$OMP DO
        do i = 0, nsamples 
          logmass(i) = minmass+i*dm
        enddo
      !$OMP ENDDO

      !$OMP DO
        do i = 1, ngalaxies
          ! determine index
          lm = log10(stellar_mass(i))
          guess = int((lm-minmass)/dm) - 1
          guess = min(guess, 1)
          do
            if (logmass(guess)>lm) then
              smf_all(guess) = smf_all(guess) + 1
              exit
            endif
            guess = guess + 1
          enddo
          
          ! check if subhalo or main
          !$OMP ATOMIC
          smf_all(guess) = smf_all(guess) + 1
          if (clmp_id(i) > 0) then
            ! found non-orphan
            !$OMP ATOMIC
            smf_sub(guess) = smf_sub(guess) + 1
            do j=1, nhalos
              if (clmp_id(i) == halo_id(j)) then
                ! found main
                !$OMP ATOMIC
                smf_main(guess) = smf_main(guess) + 1
              endif
            enddo
          endif
        enddo
      !$OMP ENDDO

      !$OMP DO
        do i=1, nsamples
          f = (logmass(i)-logmass(i-1))*vol
          smf_all(i) = smf_all(i)/f
          smf_sub(i) = smf_sub(i)/f
          smf_main(i) = smf_main(i)/f
        enddo
      !$OMP ENDDO

    !$OMP END PARALLEL


    fname=TRIM(TRIM(srcdir)//"/smf.txt")
    open(unit=666, file=fname, form='formatted')
    write(666, '(4A20)') "Mass bins", "SMF subhalos+orphans", "SMF subhalos", "SMF centrals only"
    do i=1, nsamples
      write(666, '(F20.6, 3E20.12)') (logmass(i-1)+logmass(i))*0.5, smf_all(i), smf_sub(i), smf_main(i)
    enddo
    close(666)

    write(*,*) "SMF stuff finished, results written to ", fname


    deallocate(smf_main, smf_all, smf_sub, logmass)
  end subroutine smf





  !==================================
  subroutine title(n,nchar)
  !==================================
    implicit none
    integer::n
    character(LEN=5)::nchar

    character(LEN=1)::nchar1
    character(LEN=2)::nchar2
    character(LEN=3)::nchar3
    character(LEN=4)::nchar4
    character(LEN=5)::nchar5

    if(n.ge.10000)then
       write(nchar5,'(i5)') n
       nchar = nchar5
    elseif(n.ge.1000)then
       write(nchar4,'(i4)') n
       nchar = '0'//nchar4
    elseif(n.ge.100)then
       write(nchar3,'(i3)') n
       nchar = '00'//nchar3
    elseif(n.ge.10)then
       write(nchar2,'(i2)') n
       nchar = '000'//nchar2
    else
       write(nchar1,'(i1)') n
       nchar = '0000'//nchar1
    endif
  end subroutine title





  !========================================
  subroutine write_galaxy_images()
  !========================================
    !--------------------------------------------------------
    ! write images to file that can be
    ! plotted directly with plot_fortran_galaxies.py
    !--------------------------------------------------------

    real(dp), dimension(1:nc, 1:nc) :: density_image
    real(dp) :: pixel
    integer :: i, j, k

    character(len=80) :: fname

    write(*,*) "Creating density field image for plotting."

    !$OMP DO
    do i=1, nc
      do j=1, nc
        pixel = 0.d0
        do k = 1, nc
          pixel = pixel + density_field(i, j, k) 
        enddo
        density_image(i,j) = pixel
      enddo
    enddo
    !$OMP END DO


    write(*,*) "Writing density field image to file."

    fname = TRIM('density_image.txt')
    open(unit=666, file=fname, form='formatted')
    ! note switched j/i for pyplot.imshow
    do j=1, nc
      do i=1, nc
        write(666, '(E14.7)', advance='no') density_image(i,j)
      enddo
      write(666,*) ! write newline
    enddo
    close(666)

    write(*,*) fname, " written."
 
    fname = TRIM('galaxy_field.txt')
    open(unit=666, file=fname, form='formatted')
    do i=1, ngalaxies
      write(666,*) x(i), y(i)
    enddo
    close(666)

    write(*,*) fname, " written."
  end subroutine write_galaxy_images 

end program
