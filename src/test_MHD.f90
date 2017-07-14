program test_MHD
  use MHDeos

  implicit none

  logical, parameter :: do_test = .true.

  if(do_test) then
     call test
  else
     call table_for_Bill
  endif
  
contains

  subroutine test
    integer :: io, nt
    integer, parameter :: nmax=1000 !set by MHD commons

    double precision, allocatable :: tl(:) ! (nt)
    double precision, allocatable :: rhol(:) ! (nt)
    double precision, allocatable :: res(:,:) !(nt,nres)
    character(len=128) :: abunfile, datafile

    !names of data files for MHD
    datafile='eosdat07'
    abunfile='abun.dat'

    !allocate arrays
    nt = 8
    allocate(tl(nt), rhol(nt), res(nt,nres))

    !fill T, Rho arrays 
    tl = [ 3.7617403d0, 4.2769411d0, 4.8326198d0, 5.5349764d0, &
         6.2063285d0, 6.6277368d0, 6.9877448d0, 7.1531246d0 ]

    rhol = [ -6.6830237d0, -4.9882721d0, -3.1878478d0, -1.9824848d0, &
         -0.9629616d0,  0.2434761d0,  1.5861873d0,  2.0437896d0 ]

    !read in data files
    call mhd_init(datafile,abunfile)

    !process T,Rho arrays through MHD
    call eosDT_get( tl, rhol, res)

    io = 22 !unit for output table
    open( unit=io , file='eostab.dat')
    call write_result(io,tl,res)      
    close(io)

  end subroutine test

  
  subroutine table_for_Bill
    integer, parameter :: nmax=1000 !set by MHD commons
    double precision, allocatable :: logT(:), logRho(:), logQ(:) ! (nt)
    double precision, allocatable :: res(:,:) !(nt,nres)
    double precision :: logT_min, logT_max, logQ_min, logQ_max, dlogT, dlogQ
    integer :: i, io, j, nlogT, nlogQ, nmid_logT
    character(len=128) :: abunfile, datafile, outfile

    if(command_argument_count()>0) then
       call get_command_argument(1,outfile)
    else
       outfile = 'eostab2.dat'
    endif
    
    !names of data files for MHD
    datafile='eosdat07'
    abunfile='abun2.dat'
      
    logT_min = 3.1
    logT_max = 7.7
    logQ_min = -5.0
    logQ_max = 3.0
      
    nlogT = 6
    nlogQ = 3
    nmid_logT = nlogT/2

    allocate(logT(nlogT*nlogQ), logQ(nlogT*nlogQ), logRho(nlogT*nlogQ))
    allocate(res(nlogT*nlogQ, nres))
      
    dlogT = (logT_max - logT_min)/dble(nlogT - 1)
    dlogQ = (logQ_max - logQ_min)/dble(nlogQ - 1)

    do i=1,nlogT
       do j=1,nlogQ
          logT((i-1)*nlogQ + j) = logT_min + dlogT*dble(i-1)
          logQ((i-1)*nlogQ + j) = logQ_min + dlogQ*dble(j-1)
       enddo
    enddo
    
    logRho = logQ + 2d0*logT - 12.0d0
        
    !read in data files
    call mhd_init(datafile,abunfile)

    !process T,Rho arrays through MHD
    call eosDT_get( logT, logRho, res)

    io = 22 !unit for output table
    open( unit=io , file=trim(outfile))
    call write_result(io,logT,res)      
    close(io)
 
  end subroutine table_for_Bill

  
  subroutine write_result(io,tl,res)
    integer, intent(in) :: io
    double precision, intent(in) :: tl(:), res(:,:)
    integer :: i,n
    write(io,'(99a16)') 'logT','logRho', 'entropy', 'logE', &
         'chiRho', 'chiT', 'logE_logRho', 'logE_logT', 'grad_ad', &
         'Cp', 'Cv', 'Gamma_1', 'Gamma_2', 'Gamma_3', 'f_H+', 'f_He+', &
         'f_He++', 'f_H2', 'eta', 'logPrad', 'logPgas'
    n=size(tl)
    do i=1,n
       write(io,'(1p99e16.8)') tl(i), res(i,:)
    enddo
  end subroutine write_result

end program test_MHD
