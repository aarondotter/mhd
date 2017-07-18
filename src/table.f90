program table
  use MHDeos

  implicit none

  call table_for_Bill

contains

     subroutine write_MESA_table(io,tl,res)
      integer, intent(in) :: io
      double precision, intent(in) :: tl(:), res(:,:)
      integer :: i, j, n
      double precision, parameter :: crad = 7.5657673816464059d-15
      double precision :: &
         logT, logRho, logE, &
         chiRho, chiT, dlgE_dlgRho, dlgE_dlgT, &
         f_H_plus1, f_He_plus1, &
         f_He_plus2, f_H2, eta, logPgas, logS, &
         dlnPgas_dlnT, dlnPgas_dlnRho, d2lnPgas_dlnT_dlnRho, &
         dlnE_dlnT, dlnE_dlnRho, d2lnE_dlnT_dlnRho, &
         dlnS_dlnT, dlnS_dlnRho, d2lnS_dlnT_dlnRho, &
         mu, log_free_e, dse, dpe, dsp, dPrad_dlnT, dP_dlnT, dPgas_dlnT, &
         T, rho, energy, Pgas, Prad, P, dlnP_dlnRho, dlnP_dlnT
      logT = tl(1)
      write(io,'(a22)') 'logT'
      write(io,'(1p99e22.14)') logT
      write(io,*)
      write(io,'(99a22)') &
         'logRho', 'logPgas', 'logE', 'logS', 'dlnPgas_dlnT', 'dlnPgas_dlnRho', &
         'd2lnPgas_dlnT_dlnRho', 'dlnE_dlnT', 'dlnE_dlnRho', 'd2lnE_dlnT_dlnRho', &
         'dlnS_dlnT', 'dlnS_dlnRho', 'd2lnS_dlnT_dlnRho', 'mu', 'log_free_e', &
         'eta', 'f_H+', 'f_He+', 'f_He++', 'f_H2', 'dse', 'dpe', 'dsp'
      n=size(tl)
      do i=1,n
      
         logT = tl(i)
         j = 1
         
         logRho               = res(i, 1)
         logPgas              = res(i, 2)
         logE                 = res(i, 3)
         logS                 = res(i, 4)
         chiT                 = res(i, 5)
         chiRho               = res(i, 6)
         d2lnPgas_dlnT_dlnRho = res(i, 7)
         dlgE_dlgT            = res(i, 8)
         dlgE_dlgRho          = res(i, 9)
         d2lnE_dlnT_dlnRho    = res(i,10)
         dlnS_dlnT            = res(i,11)
         dlnS_dlnRho          = res(i,12)
         d2lnS_dlnT_dlnRho    = res(i,13)
         mu                   = res(i,14)
         log_free_e           = res(i,15)
         eta                  = res(i,16)
         f_H_plus1            = res(i,17)
         f_He_plus1           = res(i,18)
         f_He_plus2           = res(i,19)
         f_H2                 = res(i,20)
         dse                  = res(i,21)
         dpe                  = res(i,22)
         dsp                  = res(i,23)

         T = 10d0**logT
         Prad = crad*T*T*T*T/3d0
         rho = 10d0**logRho
         energy = 10d0**logE
         Pgas = 10d0**logPgas
         P = Pgas + Prad


         
         dlnP_dlnT = chiT
         dlnP_dlnRho = chiRho
         dPrad_dlnT = 4d0*Prad
         dP_dlnT = P*dlnP_dlnT
         dPgas_dlnT = dP_dlnT - dPrad_dlnT
         dlnPgas_dlnT = dPgas_dlnT/Pgas
         
         dlnPgas_dlnRho = dlnP_dlnRho
         dlnE_dlnT = dlgE_dlgT
         dlnE_dlnRho = dlgE_dlgRho

         write(io,'(1p99e22.14)') &
            logRho, logPgas, logE, logS, dlnPgas_dlnT, dlnPgas_dlnRho, &
            d2lnPgas_dlnT_dlnRho, dlnE_dlnT, dlnE_dlnRho, d2lnE_dlnT_dlnRho, &
            dlnS_dlnT, dlnS_dlnRho, d2lnS_dlnT_dlnRho, mu, log_free_e, &
            eta, f_H_plus1, f_He_plus1, f_He_plus2, f_H2, dse, dpe, dsp
         
      enddo
   end subroutine write_MESA_table


  subroutine table_for_Bill
    integer, parameter :: nmax=1000 !set by MHD commons
    double precision, allocatable :: logTs(:), logRhos(:) ! (nt)
    double precision, allocatable :: res(:,:) !(nt,nres)
    double precision :: logT, logRho_min, logRho_max, dlogRho
    integer :: io, j, n
    character(len=128) :: abunfile, datafile, outfile


    logRho_min = -15.0d0
    logRho_max = -1.0d0
    logT = 7.0d0
    n = 15


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
    abunfile='abun.dat'
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
    call write_MESA_table(io,logTs,res)      
    close(io)

    write(*,*)
    do j=1,n
       write(*,*) 'logRhos(j)', j, logRhos(j)
    enddo
    write(*,*) 'logT', logT
    write(*,*)

  end subroutine table_for_Bill


end program table
