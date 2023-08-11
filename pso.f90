program pso
  implicit none

  integer, parameter :: DIMENSIONS=3*3
  integer, parameter :: NINDIVIDUALS=1

  type individual
     real :: X(DIMENSIONS)      ! The coordinates.
     real :: V(DIMENSIONS)      ! The velocity.
     real :: p_best             ! The best score.
     real :: X_best(DIMENSIONS) ! The local best coordinates.
     integer :: g_best          ! The index of global best coordinates.
  end type individual

  type(individual) :: inds(NINDIVIDUALS)

  integer :: i, ngen
  real :: score

  inds = fresh_individuals()

  ngen = 1
  ! It would be nice to do each of these individuals in parallel, and
  ! terminate them one by one.
  do
     if(ngen==1000) donep=.true.

     do i=1,NINDIVIDUALS
        ! In first iteration, the second term will be zero since they
        ! are both set to the same value.
        inds(i)%V = inds(i)%V + rand_vec()*(inds(i)%X_best - inds(i)%X) &
                  + rand_vec()*(inds(inds(i)%g_best) - inds(i)%X)
        inds(i)%X = inds(i)%X + inds(i)%V

        score = score_gen(inds(i)%X)

        if(score < inds(i)%p_best) then
           inds(i)%X_best = inds(i)%X
           inds(i)%p_best = score
        end if
     end do

     call update_g_best(inds)
  end do

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

  ! Return a batch of fresh individuals for further use.
  ! TODO: Random X must make sense!
  function fresh_individuals() result(inds)
    implicit none
    type(individual) :: inds(NINDIVIDUALS)
    integer :: i
    real :: randv(DIMENSIONS), randx(DIMENSIONS)

    do i=1,NINDIVIDUALS
       randv = rand_vec()
       randx = rand_vec()

       inds(i)%X = randx
       inds(i)%V = randv
       inds(i)%p_best = score_gen(randx)
       inds(i)%X_best = randx
    end do

    call update_g_best(inds)
  end function fresh_individuals

  ! Return the score for the coordinates X.
  ! TODO: Do the actual thing.
  function score_gen(X) result(score)
    implicit none
    real, intent(in) :: X(DIMENSIONS)
    real :: score

    score = rand(1)
  end function score_gen

  ! Return a vector of random numbers of dimension DIMENSIONS.
  ! TODO: Range my friend!
  function rand_vec() result(vec)
    implicit none
    real :: vec(DIMENSIONS)
    integer :: i

    do i=1,DIMENSIONS
       vec(i) = rand(1)
    end do
  end function rand_vec
end program pso
