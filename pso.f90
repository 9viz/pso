program pso
  implicit none

  integer, parameter :: DIMENSIONS=3*10
  integer, parameter :: NINDIVIDUALS=4
  real, parameter :: CLENGTH=8.0
  real, parameter :: CCUTOFF=1.6
  real, parameter :: CCUTOFFSQ=CCUTOFF**2
  character(len=3) :: LOGFILE="OUT"

  type individual
     real :: X(DIMENSIONS)      ! The coordinates.
     real :: V(DIMENSIONS)      ! The velocity.
     real :: p_best             ! The best score.
     real :: X_best(DIMENSIONS) ! The local best coordinates.
     integer :: g_best          ! The index of global best coordinates.
  end type individual

  type(individual) :: inds(NINDIVIDUALS)

  integer :: i, ngen, gbest, log_unit
  real :: score

  open(newunit=log_unit, file=LOGFILE)
  inds = fresh_individuals()
  call log_inds(inds, log_unit, 0)

  ngen = 1
  ! It would be nice to do each of these individuals in parallel, and
  ! terminate them one by one.
  do
     if(ngen==10) exit

!$OMP do
     do i=1,NINDIVIDUALS
        ! In first iteration, the second term will be zero since they
        ! are both set to the same value.
        gbest = inds(i)%g_best
        inds(i)%V = inds(i)%V + rand_vec(0.0, 1.0)*(inds(i)%X_best - inds(i)%X) &
                              + rand_vec(0.0, 1.0)*(inds(gbest)%X_best - inds(i)%X)
        inds(i)%X = inds(i)%X + inds(i)%V

        score = score_gen(inds(i)%X, i, ngen)

        if(score < inds(i)%p_best) then
           inds(i)%X_best = inds(i)%X
           inds(i)%p_best = score
        end if
     end do
!$OMP end do

     call update_g_best(inds)
     call log_inds(inds, log_unit, ngen)
     ngen = ngen + 1
  end do

  close(log_unit)

contains
  ! Return the new g_best value for I-th individual among INDS.
  function new_g_best(inds, i) result(g_best)
    implicit none
    integer, intent(in) :: i
    type(individual), intent(in) :: inds(NINDIVIDUALS)
    integer :: j, g_best, start, end

    if(i==1) then
       g_best=NINDIVIDUALS
       start=i
       end=i+1
    else if(i==NINDIVIDUALS) then
       g_best=1
       start=i-1
       end=i
    else
       g_best=i-1
       start=i
       end=i+1
    end if

    do j=start,end,1
       if(inds(j)%p_best < inds(g_best)%p_best) g_best=j
    end do
  end function new_g_best

  ! Update the g_best value for all individuals in INDS.
  subroutine update_g_best(inds)
    implicit none
    type(individual), intent(inout) :: inds(NINDIVIDUALS)
    integer :: i

    do i=1,NINDIVIDUALS
       inds(i)%g_best = new_g_best(inds, i)
    end do
  end subroutine update_g_best

  ! Return a random number in the range [MIN, MAX).
  function rand_range(min, max) result(r)
    implicit none
    real, intent(in) :: min, max
    real :: r

    call random_seed()
    call random_number(r)

    r = (max-min) * r + min
  end function rand_range

  ! Return a vector of random numbers of dimension DIMENSIONS.
  ! The range of random numbers is [MIN, MAX).
  function rand_vec(min, max) result(vec)
    implicit none
    real, intent(in) :: min, max
    real :: vec(DIMENSIONS)
    integer :: i

    do i=1,DIMENSIONS
       vec(i) = rand_range(min, max)
    end do
  end function rand_vec

  ! Return a vector of three random numbers in range [MIN,MAX).
  function rand_vec3(min, max) result(vec)
    implicit none
    real, intent(in) :: min, max
    real :: vec(3)
    integer :: i

    do i=1,3
       vec(i) = rand_range(min,max)
    end do
  end function rand_vec3

  ! Return the Nth atom's coordinate in X.
  function ncoords(n, X) result(coords)
    implicit none
    real, intent(in) :: X(DIMENSIONS)
    integer, intent(in) :: n
    real :: coords(3)

    coords = X((n-1)*3+1:n*3)
  end function ncoords

  ! Return .true. if Ith atom does not collide with (I-1) atoms in X.
  function check_conflict(i, X) result(conflictp)
    implicit none
    integer, intent(in) :: i
    real, intent(in) :: X(DIMENSIONS)
    logical :: conflictp
    real :: icoords(3), jcoords(3), diff(3)
    integer :: j

    conflictp = .false.

    icoords = ncoords(i, X)

    do j=1,i-1
       if(conflictp) exit
       jcoords = ncoords(j, X)
       diff = jcoords - icoords
       if(dot_product(diff, diff) < CCUTOFFSQ) conflictp = .true.
    end do
  end function check_conflict

  function new_X() result(nX)
    implicit none
    real :: nX(DIMENSIONS)
    integer :: natom, l, u
    logical :: conflictp
    real :: icoords(3)

    nX = 0.0
    natom = DIMENSIONS / 3

    do i=1,natom
       conflictp = .false.
       l = (i-1)*3
       u = i*3

       icoords = rand_vec3(-CLENGTH/2, CLENGTH/2)
       nX(l:u) = icoords
       conflictp = check_conflict(i, nX)

       if(conflictp) then
          do
             if(.not. conflictp) exit

             icoords = icoords + 0.003*icoords
             if(any(icoords > CLENGTH/2) .or. &
                  any(icoords < -CLENGTH/2)) &
                  icoords = rand_vec3(-CLENGTH/2, CLENGTH/2)

             nX(l:u) = icoords

             conflictp = check_conflict(i, nX)
          end do
       end if
    end do
  end function new_X

  ! Return a batch of fresh individuals for further use.
  function fresh_individuals() result(inds)
    implicit none
    type(individual) :: inds(NINDIVIDUALS)
    integer :: i
    real :: randv(DIMENSIONS)

!$OMP do
    do i=1,NINDIVIDUALS
       randv = rand_vec(0.0, 1.0)

       inds(i)%X = new_X()
       inds(i)%V = randv
       inds(i)%p_best = score_gen(inds(i)%X, i, 0)
       inds(i)%X_best = inds(i)%X
    end do
!$OMP end do

    call update_g_best(inds)
  end function fresh_individuals

  ! Write the g16 input file for coord vec X to FILENAME.
  subroutine write_g16_input(X, filename)
    implicit none
    real, intent(in) :: X(DIMENSIONS)
    character(len=*), intent(in) :: filename
    integer :: i, stat, unit

    open(newunit=unit, file=filename, iostat=stat)
    if(stat < 0) then
       print *, "COULD NOT CREATE FILE", filename, " FOR G16 INPUT"
       stop 1
    end if

    ! *PREAMBLE*
    write(unit, "(A)") &
         "%NProcShared=4", &
         "%Chk="//filename(:len(filename)-3)//"chk", &
         "#n UHF/6-31G(d,p) SP", &
         "", &
         filename, &
         "", &
         "0 1"

    do i=1,DIMENSIONS/3
       write(unit, "(A)", advance="no") "C"
       write(unit, *) ncoords(i, X)
    end do

    write(unit, "(A)") "", ""

    close(unit)
  end subroutine write_g16_input

  ! Return the score for the coordinates X.
  ! NIND is the index of the individual, NGEN is the index of the
  ! generation.
  ! TODO: Do the actual thing.
  function score_gen(X, nind, ngen) result(score)
    implicit none
    real, intent(in) :: X(DIMENSIONS)
    integer, intent(in) :: nind, ngen
    real :: score
    character(len=50), allocatable :: filename
    integer :: stat, unit

    filename = ""
    write(filename, "(i0,a,i0,a)") nind, "_", ngen, ".com"
    call write_g16_input(X, trim(filename))

    call execute_command_line("./calc-energy "//filename, &
         wait=.true., exitstat=stat)
    ! This can happen if atoms are too close.
    if(stat==0) then
       score = 10000000.0
       return
    end if

    open(newunit=unit, file=filename//"_energy", status="old", iostat=stat)
    if(stat < 0) then
       score = 10000000.0
    else
       read(unit,*) score
    end if
  end function score_gen

  ! Log the individuals INDS' status to file with unit UNIT.
  ! NGEN is the current generation no.
  subroutine log_inds(inds, unit, ngen)
    implicit none
    type(individual), intent(in) :: inds(NINDIVIDUALS)
    integer, intent(in) :: unit, ngen
    integer :: i

    write(unit, "(A,I0,A)") "**** GENERATION ", ngen, " ****"
    do i=1,NINDIVIDUALS
       write(unit, "(A,I0,A)") "**** INDIVIDUAL ", i, " ****"

       write(unit, *) "p_best=", inds(i)%p_best
       write(unit, *) "X_best=", inds(i)%X_best
       write(unit, *) "X=", inds(i)%X
       write(unit, *) "V=", inds(i)%V
       write(unit, *) "g_best=", inds(i)%g_best

       write(unit, "(A)") "**** END INDIVIDUAL ****"
    end do

    write(unit, "(A)") "**** END GENERATION *****"
  end subroutine log_inds
end program pso
