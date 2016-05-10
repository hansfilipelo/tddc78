program laplsolv
    use omp_lib

    !-----------------------------------------------------------------------
    ! Serial program for solving the heat conduction problem
    ! on a square using the Jacobi method.
    ! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
    ! Modified by Berkant Savas (besav@math.liu.se) April 2006
    !-----------------------------------------------------------------------
    integer, parameter                      :: n=1000, maxiter=1000, nr_threads=16
    double precision,parameter              :: tol=1.0E-3
    double precision,dimension(0:n+1,0:n+1) :: T
    double precision,dimension(n)           :: tmp
    double precision                        :: error,x
    double precision                        :: t1,t0
    integer                                 :: i,j,k
    integer                                 :: quote,my_id
    integer,dimension(nr_threads)           :: start_it,stop_it
    character(len=20)                       :: str
    double precision,dimension(0:n-1,0:nr_threads-1)   :: padding_before,padding_after

    ! Set boundary conditions and initial values for the unknowns
    T=0.0D0
    T(0:n+1 , 0)     = 1.0D0
    T(0:n+1 , n+1)   = 1.0D0
    T(n+1   , 0:n+1) = 2.0D0

    quote = n/nr_threads

    do i=0,nr_threads-2
        start_it(i) = (i*quote)+1
        stop_it(i) = (i+1)*quote
    end do
    start_it(nr_threads-1) = ((nr_threads-1)*quote)+1
    stop_it(nr_threads-1) = n

    call omp_set_num_threads(nr_threads)

    ! Solve the linear system of equations using the Jacobi method
    t0 = omp_get_wtime()

    do k=1,maxiter

        error=0.0D0

        ! Calculate start and stop criteria
        do i=0,nr_threads-2
            padding_before(0:n-1, i) = T(1:n, start_it(i)-1)
            padding_after(0:n-1, i) = T(1:n, stop_it(i)+1)
        end do
        padding_before(0:n-1, nr_threads-1) = T(1:n,start_it(nr_threads-1)-1)
        padding_after(0:n-1, nr_threads-1) = T(1:n,stop_it(nr_threads-1)+1)

        !$omp parallel private(j,tmp,my_id) shared(T,padding_before,padding_after) reduction(max: error)
        my_id = OMP_GET_THREAD_NUM()

        do j=start_it(my_id),stop_it(my_id)-1
            tmp=T(1:n,j)
            T(1:n,j) =(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+padding_before(0:n-1,my_id))/4.0D0
            error=max(error,maxval(abs(tmp-T(1:n,j))))
            padding_before(0:n-1,my_id)=tmp
        end do

        tmp=T(1:n,stop_it(my_id))
        T(1:n,stop_it(my_id)) = &
            (T(0:n-1,stop_it(my_id))+T(2:n+1,stop_it(my_id))+padding_after(0:n-1,my_id)+padding_before(0:n-1,my_id))/4.0D0
        error=max(error,maxval(abs(tmp-T(1:n,stop_it(my_id)))))
        !$omp end parallel

        if (error<tol) then
            exit
        end if

    end do

    t1 = omp_get_wtime()

    write(unit=*,fmt=*) 'Time:',t1-t0,'Number of Iterations:',k
    write(unit=*,fmt=*) 'Temperature of element T(1,1)  =',T(1,1)

    ! Uncomment the next part if you want to write the whole solution
    ! to a file. Useful for plotting.

    !open(unit=7,action='write',file='result.dat',status='unknown')
    !write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
    !do i=0,n+1
    !   write (unit=7,fmt=str) T(i,0:n+1)
    !end do
    !close(unit=7)

end program laplsolv
