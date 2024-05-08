program pso
  implicit none

  integer, parameter :: DIMENSIONS=2
  integer, parameter :: NINDIVIDUALS=4
  real, parameter :: K_D = 4.6
  real, parameter :: D_0 = 1.113
  real, parameter :: K_A = 0.36
  real, parameter :: A_0 = 1.919 ! rad ≡ 110°
  real, parameter :: PI=3.14159265359
  character(len=11), parameter :: LOGFILE="OUT.methane"

  ! The coordinates are the bond angle and the bond distance.
  type individual
     real :: X(DIMENSIONS)      ! The coordinates.
     real :: V(DIMENSIONS)      ! The velocity.
     real :: best_score         ! The best score.
     real :: X_best(DIMENSIONS) ! The local best coordinates.
  end type individual

  type(individual) :: inds(NINDIVIDUALS)

  integer :: i, ngen, log_unit
  real :: score, g_best, g_X_best(DIMENSIONS)

  g_best = 1000

  open(newunit=log_unit, file=LOGFILE)
  inds = fresh_individuals()
  call log_inds(inds, log_unit, 0)

  ngen = 1
  ! It would be nice to do each of these individuals in parallel, and
  ! terminate them one by one.
  do
     if(ngen==100) exit

     !$OMP do
     do i=1,NINDIVIDUALS
        ! In first iteration, the second term will be zero since they
        ! are both set to the same value.
        inds(i)%V = inds(i)%V + rand_vec(0.0, 1.0)*(inds(i)%X_best - inds(i)%X) &
                              + rand_vec(0.0, 1.0)*(g_X_best - inds(i)%X)
        inds(i)%X = inds(i)%X + inds(i)%V

        score = score_gen(inds(i)%X)

        if(score < inds(i)%best_score) then
           inds(i)%X_best = inds(i)%X
           inds(i)%best_score = score
        end if
     end do
     !$OMP end do

     call update_g_best(inds)
     call log_inds(inds, log_unit, ngen)
     ngen = ngen + 1
  end do

  close(log_unit)

contains
  ! Update the global best energy.
  subroutine update_g_best(inds)
    implicit none
    type(individual), intent(inout) :: inds(NINDIVIDUALS)
    integer :: i

    do i=1,NINDIVIDUALS
       if(inds(i)%best_score < g_best) then
          g_best = inds(i)%best_score
          g_X_best = inds(i)%X
       end if
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

  function new_X() result(nX)
    implicit none
    ! (Bond angle, Bond distance)
    real :: nX(DIMENSIONS)

    nX(1) = rand_range(0.0, 2*PI)
    nX(2) = rand_range(0.5, 3.0)
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
       inds(i)%best_score = score_gen(inds(i)%X)
       inds(i)%X_best = inds(i)%X
    end do
    !$OMP end do

    call update_g_best(inds)
  end function fresh_individuals

  ! Return the score for the coordinates X.
  function score_gen(X) result(score)
    implicit none
    real, intent(in) :: X(DIMENSIONS)
    real :: score, diff_d, diff_a

    ! if(X(2) < 0 .or. X(2) > 5.0 .or. X(1) > 2*pi .or. X(1) < 0) then
       ! score = 100000000.0
    ! else
       ! See allinger1977.pdf:4, table II.
       ! The unit is in kcal/mol.
       diff_d = X(2) - D_0
       diff_a = X(1) - A_0
       score = 4 * (K_D * diff_d**2) / 2 &
            + 6 * (K_A * diff_a**2) / 2
    ! end if
  end function score_gen

  ! Log the individuals INDS' status to file with unit UNIT.
  ! NGEN is the current generation no.
  subroutine log_inds(inds, unit, ngen)
    implicit none
    type(individual), intent(in) :: inds(NINDIVIDUALS)
    integer, intent(in) :: unit, ngen
    integer :: i

    do i=1,NINDIVIDUALS
       write(unit, "(f12.4,f12.4)", advance="no") 360*inds(i)%X_best(1)/2/PI, inds(i)%X(2)
    end do
    write(unit, *) 360*g_X_best(1)/2/PI, g_X_best(2:)

    ! write(unit, "(A,I0,A)") "**** GENERATION ", ngen, " ****"
    ! write(unit, *) "g_best=", g_best
    ! write(unit, *) "g_X_best=", 360*g_X_best(1)/2/PI, g_X_best(2:)
    ! do i=1,NINDIVIDUALS
    !    write(unit, "(A,I0,A)") "**** INDIVIDUAL ", i, " ****"

    !    write(unit, *) "best_score=", inds(i)%best_score
    !    write(unit, *) "X_best=", 360*inds(i)%X_best(1)/2/PI, inds(i)%X_best(2:)
    !    write(unit, *) "X=", 360*inds(i)%X(1)/2/PI, inds(i)%X(2:)
    !    write(unit, *) "V=", inds(i)%V

    !    write(unit, "(A)") "**** END INDIVIDUAL ****"
    ! end do

    ! write(unit, "(A)") "**** END GENERATION *****"
  end subroutine log_inds
end program pso
