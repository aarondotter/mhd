program test_MHD
  use MHDeos

  implicit none

  double precision :: logPgas_res, logE_res, logS_res, logEgas_res, logSgas_res

  double precision :: dlnPgas_dlnT, dlnPgas_dlnd, d2lnPgas_dlnd_dlnT
  double precision :: dlnE_dlnT, dlnE_dlnd, d2lnE_dlnd_dlnT
  double precision :: dlnS_dlnT, dlnS_dlnd, d2lnS_dlnd_dlnT

  double precision :: X, Z

  integer :: ierr
  logical :: doing_1st_call = .true.


  logical :: doing_d_dlnd
  double precision :: lnd, lnT, logT, logRho, logRho_min, logRho_max, logT_min, logT_max, dlogT, dlogRho
  integer :: i_var ! 1 = logPgas, 2 = logE, 3 = logS
  character (len=132) :: results_filename, abund_filename

  double precision, parameter :: ln10 = 2.3025850929940455D0 ! = log(10d0)

  Z=0.02d0
  X=0.6d0

  call write_abun_file(99, X, Z, abund_filename, ierr)
  if(ierr/=0) stop 'problem in write_abun_file'

  logT_max = 8.3d0
  logT_min = 3.3d0
  dlogT = 0.05d0
  dlogRho = 1.0d0
  logRho_min = -4.0d0
  logRho_max = -4.0d0
  results_filename='test_logRho-4_varyT.data'

  call get_results( &
       logT_Min, logT_max, dlogT, logRho_min, logRho_max, dlogRho, &
       results_filename, abund_filename)

contains

  subroutine write_abun_file(io,X,Z,abunfile,ierr)
    integer, intent(in) :: io
    double precision, intent(in) :: X, Z
    character(len=132), intent(out) :: abunfile
    integer, intent(out) :: ierr
    integer, parameter :: n_elem = 6 ! H - Ne, eosdat06
    integer :: i, j, num_levels, levels(n_elem), atomic_number(n_elem)
    double precision :: Zmet(n_elem-2), num_frac(n_elem), mass_frac(n_elem)
    double precision :: weight(n_elem), Y
    character(len=4) :: id(n_elem)

    ierr=0        
    num_frac = 0d0
    mass_frac = 0d0

    !check valid range of X, Z
    if(X+Z > 1d0 .or. X+Z < 0d0) then
       ierr=-1
       write(*,*) ' X+Z out of range! '
       write(*,*) ' X = ', X
       write(*,*) ' Z = ', Z
       write(*,*) ' X + Z = ', X+Z
       return
    elseif(X>1d0.or.X<0d0)then
       ierr=-1
       write(*,*) 'X out of range!'
       write(*,*) ' X = ', X
    elseif(Z>1d0.or.Z<0d0)then
       ierr=-1
       write(*,*) 'Z out of range!'
       write(*,*) ' Z = ', Z
    endif

    ! H He  C  N  O  Ne  Fe
    atomic_number = [ 1, 2, 6, 7, 8, 10]
    levels        = [ 3, 2, 6, 7, 8, 10]
    id = ['hydr', 'heli', 'carb', 'nitr', 'oxyg', 'neon']
    weight = [ 1.0079d0, 4.0026d0, 12.011d0, 14.0067d0, 15.9994d0, 20.179d0]

    !                  C           N            O           Ne
    Zmet = [ 1.906614d-1, 5.58489d-2, 5.429784d-1, 2.105114d-1]


    !calculate mass fractions then number fractions
    Y = 1d0 - X - Z
    mass_frac(1) = X
    mass_frac(2) = Y
    mass_frac(3:n_elem) = Z*Zmet

    mass_frac = mass_frac/sum(mass_frac)

    num_frac = mass_frac/weight

    num_frac = num_frac/sum(num_frac)

    !write abun file
    abunfile = 'abun.dat'
    open(io,file=trim(abunfile),status='unknown')
    write(io,'(a)') 'full partition functions: h-fe'
    write(io,'(i5)') n_elem

    !have to treat H separately
    write(io,'(2x,a4,2i5,f10.4,1pe15.8)') id(1), atomic_number(1), levels(1)+2, &
         weight(1), num_frac(1)
    !other elements normal
    do i =2, n_elem
       write(io,'(2x,a4,2i5,f10.4,1pe15.8)') id(i), atomic_number(i), levels(i)+1, &
            weight(i), num_frac(i)
    enddo

    num_levels = sum(levels)
    write(io,'(i5)') num_levels


    do i = 1, n_elem
       if(i==1)then
          write(io,'(2i5)') i, 1
          write(io,'(2i5)') i, 2
          write(io,'(2i5)') i, 4
       else
          do j=1,levels(i)
             write(io,'(2i5)') atomic_number(i), j
          enddo
       endif
    enddo

    close(io)

  end subroutine write_abun_file

  subroutine get_results( &
       logT_min, logT_max, dlogT, logRho_min, logRho_max, dlogRho, &
       results_file, abund_file)
    double precision, intent(in) :: logT_min, logT_max, logRho_min, logRho_max, dlogRho, dlogT
    character (len=*) , intent(in):: results_file, abund_file
    integer, parameter :: io = 22
    integer :: i, j, nT, nR

    results_filename = results_file
    abund_filename = abund_file 

    ! for bicubic splines must have at least 4 points in logRho and logT

    ! abun_z_0.0, abun_z_0.02, abun_z_0.2, abun_z_0.4, abun_z_0.5, abun_z_0.6

    ! logT = 7.700, logRho_min = -10, logRho_max = 3.2
    ! logT = 7.275, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.850, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.425, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.000, logRho_min = -10, logRho_max = 3.2

    ! logT = 6.00, logRho_min = -15, logRho_max = 1.2
    ! logT = 5.52, logRho_min = -15, logRho_max = 1.2
    ! logT = 5.04, logRho_min = -15, logRho_max = 1.2
    ! logT = 4.56, logRho_min = -15, logRho_max = 1.2
    ! logT = 4.08, logRho_min = -15, logRho_max = 1.2
    ! logT = 3.60, logRho_min = -15, logRho_max = 1.2

    ! logT = 3.600, logRho_min = -15, logRho_max = -1.2
    ! logT = 3.425, logRho_min = -15, logRho_max = -1.2
    ! logT = 3.250, logRho_min = -15, logRho_max = -1.2
    ! logT = 3.075, logRho_min = -15, logRho_max = -1.2
    ! logT = 2.900, logRho_min = -15, logRho_max = -1.2


    ! abun_z_0.7, abun_z_1.0

    ! logT = 7.7, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.8, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.0, logRho_min = -10, logRho_max = 3.2

    ! logT = 6.0, logRho_min = -15, logRho_max = 0.1
    ! logT = 4.8, logRho_min = -15, logRho_max = 0.1
    ! logT = 3.6, logRho_min = -15, logRho_max = 0.1

    ! logT = 3.6, logRho_min = -15, logRho_max = -1.2
    ! logT = 3.2, logRho_min = -15, logRho_max = -1.2
    ! logT = 2.9, logRho_min = -15, logRho_max = -1.2

    nR = floor((logRho_max - logRho_min)/dlogRho + 1d-6) + 1
    nT = floor((logT_max - logT_min)/dlogT + 1d-6) + 1
    write(*,*) 'logT_min', logT_min
    write(*,*) 'logT_max', logT_max
    write(*,*) 'dlogT', dlogT
    write(*,*) 'nT', nT
    write(*,*)
    write(*,*) 'logRho_min', logRho_min
    write(*,*) 'logRho_max', logRho_max
    write(*,*) 'dlogRho', dlogRho
    write(*,*) 'nR', nR
    write(*,*)

    open(unit=io, file=trim(results_filename))
    call write_logT_header(io)

    do i=1,nT
       logT = logT_min + dble(i-1)*dlogT
       do j=1,nR
          logRho = logRho_min + dble(j-1)*dlogRho
          write(*,*) 'logT logRho', logT, logRho
          write(*,*) 'logPgas partials'
          call get_partials(1, dlnPgas_dlnT, dlnPgas_dlnd, d2lnPgas_dlnd_dlnT)     
          write(*,*) 'logE partials'
          call get_partials(2, dlnE_dlnT, dlnE_dlnd, d2lnE_dlnd_dlnT)     
          write(*,*) 'logS partials'
          call get_partials(3, dlnS_dlnT, dlnS_dlnd, d2lnS_dlnd_dlnT)     
          write(*,*) 'other values'
          call eval1(logRho, logT, io)
       end do
    enddo

    close(io)

  end subroutine get_results


  ! uses global logT, logRho
  subroutine get_partials(which, d_dlnT, d_dlnd, d2_dlnd_dlnT)
    integer, intent(in) :: which
    double precision, intent(out) :: d_dlnT, d_dlnd, d2_dlnd_dlnT
    double precision :: hx, hy, err_d_dlnd, err_d_dlnT, err_d2_dlnd_dlnT, err_tol
    logical, parameter :: do_2nd_partial = .false.
    logical, parameter :: do_1st_partial = .true.

    lnT = logT*ln10
    lnd = logRho*ln10

    ! may need to adjust these as we get experience
    !hx = max(1d-7, abs(lnd*1d-5))   
    !hy = max(1d-7, abs(lnT*1d-5))
    !hx = 1d-5  
    !hy = 1d-5 
    hx = 1d-6  
    hy = 1d-6 

    ! pick which variable
    i_var = which ! 1 = logPgas, 2 = logE, 3 = logS
    write(*,*) 'get_partials: i_var logT logRho', i_var, logT, logRho
    write(*,*)

    if(do_1st_partial)then
       ! 1st wrt lnT
       doing_d_dlnd = .false.
       d_dlnT = dfridr(hy,err_d_dlnT)
       write(*,*) 'd_dlnT, err_d_dlnT', d_dlnT, err_d_dlnT
       write(*,*)

       ! 1st wrt lnd
       doing_d_dlnd = .true.
       d_dlnd = dfridr(hx,err_d_dlnd)
       write(*,*) 'd_dlnd, err_d_dlnd', d_dlnd, err_d_dlnd
       write(*,*)
    else
       d_dlnT = 0d0
       d_dlnd = 0d0
       err_d_dlnT = 0d0
       err_d_dlnd = 0d0
    endif

    if(do_2nd_partial)then
       ! mixed 2nd
       hx = 1d-3  
       hy = 1d-3 
       d2_dlnd_dlnT = dfridr2(hx,hy,err_d2_dlnd_dlnT)
       write(*,*) 'd2_dlnd_dlnT, err_d_dlnd', d2_dlnd_dlnT, err_d2_dlnd_dlnT
       write(*,*)
    else
       d2_dlnd_dlnT = 0d0
       err_d2_dlnd_dlnT = 0d0
    endif

    err_tol = 1d-6
    if (err_d_dlnT > err_tol) write(*,*) 'BAD ERROR for d_dlnT', err_d_dlnT
    if (err_d_dlnd > err_tol) write(*,*) 'BAD ERROR for d_dlnd', err_d_dlnd
    if (err_d2_dlnd_dlnT > err_tol) write(*,*) 'BAD ERROR for d2_dlnd_dlnT', err_d2_dlnd_dlnT

  end subroutine get_partials


  double precision function get1_val(logRho_current, logT_current) result(val)
    double precision, intent(in) :: logRho_current, logT_current
    integer, parameter :: io = -1 ! no file output
    call eval1(logRho_current, logT_current, io)
    doing_1st_call = .false.
    if (i_var == 1) then
       val = logPgas_res*ln10
    else if (i_var == 2) then
       val = logE_res*ln10
    else if (i_var == 3) then
       val = logS_res*ln10
    else
       stop 'bad i_var'
    end if
    !write(*,*) 'logRho, logT, val/ln10', logRho_current, logT_current, val/ln10
  end function get1_val


  double precision function dfridr_func(delta_x) result(val)
    double precision, intent(in) :: delta_x
    double precision :: logRho_current, logT_current
    if (doing_d_dlnd) then
       logRho_current = (lnd + delta_x)/ln10
       logT_current = logT
    else
       logT_current = (lnT + delta_x)/ln10
       logRho_current = logRho
    end if
    val = get1_val(logRho_current, logT_current)
  end function dfridr_func


  double precision function dfridr2_func(delta_x, delta_y) result(val)
    double precision, intent(in) :: delta_x, delta_y
    double precision :: logRho_current, logT_current
    logRho_current = (lnd + delta_x)/ln10
    logT_current = (lnT + delta_y)/ln10
    val = get1_val(logRho_current, logT_current)
  end function dfridr2_func


  double precision function dfridr(hx,err) ! from Frank
    double precision, intent(in) :: hx
    double precision, intent(out) :: err
    !  this routine returns the first derivative of a function func(x)
    !  at the point x, by ridders method of polynomial extrapolation.
    !  value hx is the initial step size;
    !  it should be an increment for which func changes substantially.
    !  an estimate of the error in the first derivative is returned in err.
    integer, parameter :: ntab = 20
    integer :: i,j
    double precision :: errt,fac,hh,a(ntab,ntab)
    double precision, parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
    dfridr = 0d0
    hh = hx
    ! 2nd order central difference
    a(1,1) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
    write(*,*) '       dfdx hh', 1, a(1,1), hh
    err = big
    ! succesive columns in the neville tableu will go to smaller stepsizes
    ! and higher orders of extrapolation
    do i=2,ntab
       hh = hh/con
       a(1,i) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
       ! compute extrapolations of various orders; the error stratagy is to compare
       ! each new extrapolation to one order lower but both at the same stepsize
       ! and at the previous stepsize
       fac = con2
       do j=2,i
          a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
          fac = con2*fac
          errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt <= err) then
             err = errt
             dfridr = a(j,i)
             write(*,*) '    dfridr err', i, j, dfridr, err
          end if
       end do
       ! if higher order is worse by a significant factor safe, then bail
       if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
          write(*,*) '    higher order is worse', err, a(i,i), a(i-1,i-1)
          return
       end if
    end do
  end function dfridr


  double precision function dfridr2(hx,hy,err)

    !  this routine returns the second derivative of a dfridr2_function dfridr2_func(x,y)
    !  at the point x, y, by ridders method of polynomial extrapolation.
    !  values hx and hy are input as a guess of the initial step size;
    !  they should not be small, but be an increment which dfridr2_func changes substantially.
    !  an estimate of the error in the second derivative is returned in err.

    !  second order difference expression for the derivative

    double precision, intent(in) :: hx, hy
    double precision, intent(out) :: err

    integer :: i,j
    integer, parameter :: ntab = 20
    double precision :: hh,gg,fac
    double precision, parameter :: con=1.4d0, con2=con*con, big=1.0d50, safe=2.0d0
    double precision :: a(ntab,ntab),errt

    dfridr2 = 0d0
    err = 0d0

    if (hx == 0.0 .or. hy == 0.0) stop 'bad hx,hy to routine dfridr2'
    hh = hx
    gg = hy

    a(1,1) = (dfridr2_func(hh,gg) - dfridr2_func(hh,-gg) &
         - dfridr2_func(-hh,gg) + dfridr2_func(-hh,-gg)) &
         /(4.0d0*gg*hh)
    write(*,*) '        dfdx hh', 1, a(1,1), hh, gg
    err = big

    !  successive columns in the neville tableu will go to smaller stepsizes
    !  and higher orders of extrapolation

    do i=2,ntab
       hh = hh/con
       gg = gg/con
       a(1,i) = (dfridr2_func(hh,gg) - dfridr2_func(hh,-gg) &
            - dfridr2_func(-hh,gg) + dfridr2_func(-hh,-gg)) &
            /(4.0d0*gg*hh)
       fac = con2

       !  compute extrapolations of various orders; the error strategy is to compare
       !  each new extrapolation to one order lower but both at the same stepsize
       !  and at the previous stepsize

       do j=2,i
          a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1.0d0)
          fac = con2 * fac
          errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt <= err) then
             err = errt
             dfridr2 = a(j,i)
             write(*,*) '    dfridr2 err', i, j, dfridr2, err
          end if
       enddo

       !  if higher order is worse by a significant factor safe, then bail

       if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
          write(*,*) '     higher order is worse', err, a(i,i), a(i-1,i-1)
          return
       end if
    enddo

  end function dfridr2


  subroutine eval1(logRho, logT, io)
    double precision, intent(in) :: logRho, logT
    integer, intent(in) :: io
    integer :: nt
    double precision, allocatable :: tl(:) ! (nt)
    double precision, allocatable :: rhol(:) ! (nt)
    double precision, allocatable :: res(:,:) !(nt,nres)
    character(len=128) :: datafile
    datafile='eosdat06'
    nt = 1
    allocate(tl(nt), rhol(nt), res(nt,nres))
    tl = [ logT ]
    rhol = [ logRho ]
    if (io > 0) write(*,*) 'logRho, logT', logRho, logT
    if (doing_1st_call) then !read in data files
       call mhd_init(datafile,abund_filename)
       doing_1st_call = .false.
    end if
    call eosDT_get(tl, rhol, res)
    call write_result(io,tl,res)      
  end subroutine eval1

  subroutine write1(logRho, logT, outfile)
    double precision, intent(in) :: logRho, logT
    character (len=*), intent(in) :: outfile
    integer, parameter :: io = 22
    open(unit=io, file=trim(outfile))
    call eval1(logRho, logT, io)
    close(io)
  end subroutine write1


  subroutine write2(logRho1, logT1, outfile1, logRho2, logT2, outfile2)
    double precision, intent(in) :: logRho1, logT1, logRho2, logT2
    character (len=*), intent(in) :: outfile1, outfile2
    call write1(logRho1, logT1, outfile1)
    write(*,*)
    write(*,*) 'done ' // trim(outfile1)
    write(*,*)
    write(*,*)
    call write1(logRho2, logT2, outfile2)
    write(*,*)
    write(*,*) 'done ' // trim(outfile2)
    write(*,*)
    write(*,*)
  end subroutine write2


  subroutine test_write2
    double precision :: logRho1, logT1, logRho2, logT2
    character (len=132) :: outfile1, outfile2
    logRho1 = -4d0
    logT1 = 5d0
    outfile1 = 'eosMHD_test1.data'
    logRho2 = -2d0
    logT2 = 6d0
    outfile2 =  'eosMHD_test2.data'
    !call write2(logRho1, logT1, outfile1, logRho2, logT2, outfile2)
    ! switch order to test that get same results
    call write2(logRho2, logT2, outfile2, logRho1, logT1, outfile1)
  end subroutine test_write2


  subroutine table_for_Bill
!!$    integer, parameter :: nmax=1000 !set by MHD commons
!!$    double precision, allocatable :: logTs(:), logRhos(:) ! (nt)
!!$    double precision, allocatable :: res(:,:) !(nt,nres)
!!$    double precision :: logT, logRho_min, logRho_max, dlogRho
!!$    integer :: io, j, n
!!$    character(len=128) :: datafile, outfile
!!$
!!$    logRho_min = -9.4d0  
!!$    logRho_max = 2.6d0
!!$    logT = 7.275d0

    ! for bicubic splines must have at least 4 points in logRho and logT

    ! abun_z_0.0, abun_z_0.02, abun_z_0.2, abun_z_0.4, abun_z_0.5, abun_z_0.6

    ! logT = 7.700, logRho_min = -10, logRho_max = 3.2
    ! logT = 7.275, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.850, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.425, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.000, logRho_min = -10, logRho_max = 3.2

    ! logT = 6.00, logRho_min = -15, logRho_max = 1.2
    ! logT = 5.52, logRho_min = -15, logRho_max = 1.2
    ! logT = 5.04, logRho_min = -15, logRho_max = 1.2
    ! logT = 4.56, logRho_min = -15, logRho_max = 1.2
    ! logT = 4.08, logRho_min = -15, logRho_max = 1.2
    ! logT = 3.60, logRho_min = -15, logRho_max = 1.2

    ! logT = 3.600, logRho_min = -15, logRho_max = -1.2
    ! logT = 3.425, logRho_min = -15, logRho_max = -1.2
    ! logT = 3.250, logRho_min = -15, logRho_max = -1.2
    ! logT = 3.075, logRho_min = -15, logRho_max = -1.2
    ! logT = 2.900, logRho_min = -15, logRho_max = -1.2


    ! abun_z_0.7, abun_z_1.0

    ! logT = 7.7, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.8, logRho_min = -10, logRho_max = 3.2
    ! logT = 6.0, logRho_min = -10, logRho_max = 3.2

    ! logT = 6.0, logRho_min = -15, logRho_max = 0.1
    ! logT = 4.8, logRho_min = -15, logRho_max = 0.1
    ! logT = 3.6, logRho_min = -15, logRho_max = 0.1

    ! logT = 3.6, logRho_min = -15, logRho_max = -1.2
    ! logT = 3.2, logRho_min = -15, logRho_max = -1.2
    ! logT = 2.9, logRho_min = -15, logRho_max = -1.2


!!$    !names of data files for MHD
!!$    datafile='eosdat07'
!!$    outfile = 'eosMHD.data'
!!$
!!$    dlogRho = 0.8d0
!!$    n = floor((logRho_max - logRho_min)/dlogRho + 1d-6) + 1
!!$    write(*,*) 'logRho_min', logRho_min
!!$    write(*,*) 'logRho_max', logRho_max
!!$    write(*,*) 'dlogRho', dlogRho
!!$    write(*,*) 'n', n
!!$    write(*,*)
!!$
!!$    allocate(logTs(n), logRhos(n))
!!$    allocate(res(n, nres))
!!$
!!$    logTs(1:n) = logT
!!$    do j=1,n
!!$       logRhos(j) = logRho_min + dble(j-1)*dlogRho
!!$       write(*,*) 'logRhos(j)', j, logRhos(j)
!!$    enddo
!!$
!!$    write(*,*) 'read in data files'
!!$    !read in data files
!!$    call mhd_init(datafile,abund_filename)
!!$
!!$    write(*,*) 'process T,Rho arrays through MHD'
!!$    !process T,Rho arrays through MHD
!!$    call eosDT_get( logTs, logRhos, res)
!!$
!!$    write(*,*) 'write results ', trim(outfile)
!!$    io = 22 !unit for output table
!!$    open( unit=io , file=trim(outfile))
!!$    call write_result(io,logTs,res)      
!!$    close(io)
!!$
!!$    write(*,*)
!!$    do j=1,n
!!$       write(*,*) 'logRhos(j)', j, logRhos(j)
!!$    enddo
!!$    write(*,*) 'logT', logT
!!$    write(*,*)

  end subroutine table_for_Bill


  subroutine write_logT_header( io)
    integer, intent(in) :: io
    write(io,'(99a22)') 'logT', &
         'logRho', 'logPgas', 'logE', 'logS', 'dlnPgas_dlnT', 'dlnPgas_dlnd', &
         'd2lnPgas_dlnd_dlnT', 'dlnE_dlnT', 'dlnE_dlnd', 'd2lnS_dlnd_dlnT', &
         'dlnS_dlnT', 'dlnS_dlnd', 'd2lnS_dlnd_dlnT', 'mu', 'log_free_e', &
         'eta', 'f_H+', 'f_He+', 'f_He++', 'f_H2', 'dse', 'dpe', 'dsp', &
         'rho', 'T', 'Pgas', 'Egas', 'gas_gamma'
  end subroutine write_logT_header


  subroutine write_result(io,tl,res)
    integer, intent(in) :: io
    double precision, intent(in) :: tl(:), res(:,:)
    integer :: j
    double precision :: &
         logT, logRho, chiRho, chiT, grad_ad, &
         Cp, Cv, Gamma_1, Gamma_2, Gamma_3, f_H_plus1, f_He_plus1, &
         f_He_plus2, f_H2, eta, Prad, T, rho, Pgas, entropy, energy, P, &
         dPrad_dT, dS_dRho, dS_dT, dE_dT, dE_dRho, dP_dT, mu, &
         log_free_e, dse, dpe, dsp, Egas, gas_gamma

    if (size(tl) /= 1) stop 'write_result expects size(tl) = 1'

    logT = tl(1)
    j = 1
    logRho = res(1,j); j=j+1
    logS_res = res(1,j) - logRho; j=j+1
    logE_res = res(1,j) - logRho; j=j+1
    chiRho = res(1,j); j=j+1
    chiT = res(1,j); j=j+1
    logEgas_res = res(1,j) - logRho; j=j+1
    logSgas_res = res(1,j); j=j+1
    grad_ad = res(1,j); j=j+1
    Cp = res(1,j); j=j+1
    Cv = res(1,j); j=j+1
    Gamma_1 = res(1,j); j=j+1
    Gamma_3 = res(1,j); j=j+1
    mu = res(1,j); j=j+1
    f_H_plus1 = res(1,j); j=j+1
    f_He_plus1 = res(1,j); j=j+1
    f_He_plus2 = res(1,j); j=j+1
    f_H2 = res(1,j); j=j+1
    eta = res(1,j); j=j+1
    Prad = res(1,j); j=j+1
    logPgas_res = res(1,j); j=j+1

    if (io <= 0) return

    Gamma_2 = 0d0
    log_free_e = 0
    T = 10d0**logT
    rho = 10d0**logRho
    Prad = Prad/rho
    Pgas = 10d0**logPgas_res
    Egas = 10d0**logEgas_res
    gas_gamma = rho*Egas/Pgas + 1d0
    entropy = 10d0**logS_res
    energy = 10d0**logE_res
    P = Pgas + Prad
    dPrad_dT = 4d0*Prad/T
    dS_dT = dlnS_dlnT*entropy/T
    dS_dRho = dlnS_dlnd*entropy/rho
    dE_dT = dlnE_dlnT*energy/T
    dE_dRho = dlnE_dlnd*energy/rho
    dP_dT = dlnPgas_dlnT*Pgas/T + dPrad_dT

    ! dse = T ∂S/∂T|_⍴ / ∂E/∂T|_⍴  - 1.0d0
    ! dpe = (⍴^2 ∂E/∂⍴|_T + T ∂P/∂T_⍴) / P - 1.0d0
    ! dsp = -(∂S/∂⍴|_T * ⍴^2) / ∂P/∂T|_⍴ - 1.0d0
    dse = T*dS_dT/dE_dT - 1d0
    dpe = (rho*rho*dE_dRho + T*dP_dT)/P - 1d0
    dsp = -rho*rho*dS_dRho/dP_dT - 1d0

    write(io,'(1p99e22.14)') &
         logT, logRho, logPgas_res, logE_res, logS_res, dlnPgas_dlnT, dlnPgas_dlnd, &
         d2lnPgas_dlnd_dlnT, dlnE_dlnT, dlnE_dlnd, d2lnS_dlnd_dlnT, &
         dlnS_dlnT, dlnS_dlnd, d2lnS_dlnd_dlnT, mu, log_free_e, &
         eta, f_H_plus1, f_He_plus1, f_He_plus2, f_H2, dse, dpe, dsp, &
         rho, T, Pgas, Egas, gas_gamma

  end subroutine write_result

end program test_MHD
