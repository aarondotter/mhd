program table
  use MHDeos
  
  implicit none
  
  call table_for_Bill
  
contains
  
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
    abunfile='abun.dat'
      
    logT_min = 6.0d0
    logT_max = 6.0d0
    logQ_min = -5.4666668574015302
    logQ_max = -5.4666668574015302
    
    nlogT = 1
    nlogQ = 1
    nmid_logT = nlogT/2

    dlogT = 0.0d0
    dlogQ = 0.0d0
    
    allocate(logT(nlogT*nlogQ), logQ(nlogT*nlogQ), logRho(nlogT*nlogQ))
    allocate(res(nlogT*nlogQ, nres))

    if(nlogT > 1)then
       dlogT = (logT_max - logT_min)/dble(nlogT - 1)
    else
       dlogT = 0.0d0
    endif
    if(nlogQ > 1)then
       dlogQ = (logQ_max - logQ_min)/dble(nlogQ - 1)
    else
       dlogT = 0.0d0
    endif

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

end program table
