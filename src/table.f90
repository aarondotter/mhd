program table
  use MHDeos

  implicit none

  call table_for_Bill

contains


  subroutine table_for_Bill
    integer, parameter :: nmax=1000 !set by MHD commons
    double precision, allocatable :: logTs(:), logRhos(:) ! (nt)
    double precision, allocatable :: res(:,:) !(nt,nres)
    double precision :: logT, logRho_min, logRho_max, dlogRho
    integer :: io, j, n
    character(len=128) :: abunfile, datafile, outfile


    logRho_min = -15.0  
    logRho_max = 3.2
    logT = 2.9
    n = 3


    ! abun_z_0.0, abun_z_0.02, abun_z_0.2, abun_z_0.4, abun_z_0.5, abun_z_0.6

    ! logT = 7.7, logRho_min = -15, logRho_max = 3.2
    ! logT = 6.0, logRho_min = -15, logRho_max = 3.2

    ! logT = 5.9, logRho_min = -15, logRho_max = 1.2
    ! logT = 3.6, logRho_min = -15, logRho_max = 1.2

    ! logT = 3.5, logRho_min = -15, logRho_max = -1.2
    ! logT = 2.9, logRho_min = -15, logRho_max = -1.2


    ! abun_z_0.7, abun_z_1.0

    ! logT = 7.7, logRho_min = -15, logRho_max = 3.2
    ! logT = 6.0, logRho_min = -15, logRho_max = 3.2

    ! logT = 5.9, logRho_min = -15, logRho_max = 0.1
    ! logT = 3.6, logRho_min = -15, logRho_max = 0.1

    ! logT = 3.5, logRho_min = -15, logRho_max = -1.2
    ! logT = 2.9, logRho_min = -15, logRho_max = -1.2


    !names of data files for MHD
    datafile='eosdat07'
    abunfile='abun_z_0.02.dat'
    outfile = 'eostab2.dat'

    dlogRho = (logRho_max - logRho_min)/dble(n - 1)
    write(*,*) 'logRho_min', logRho_min
    write(*,*) 'logRho_max', logRho_max
    write(*,*) 'dlogRho', dlogRho

    allocate(logTs(n), logRhos(n))
    allocate(res(n, nres))

    logTs(1:n) = logT
    do j=1,n
       logRhos(j) = logRho_min + dble(j-1)*dlogRho
       write(*,*) 'logRhos(j)', j, logRhos(j)
    enddo

    write(*,*) 'read in data files'
    !read in data files
    call mhd_init(datafile,abunfile)

    write(*,*) 'process T,Rho arrays through MHD'
    !process T,Rho arrays through MHD
    call eosDT_get( logTs, logRhos, res)

    write(*,*) 'write results', trim(outfile)
    io = 22 !unit for output table
    open( unit=io , file=trim(outfile))
    call write_result(io,logTs,res)      
    close(io)

    write(*,*)
    do j=1,n
       write(*,*) 'logRhos(j)', j, logRhos(j)
    enddo
    write(*,*) 'logT', logT
    write(*,*)

  end subroutine table_for_Bill


end program table
