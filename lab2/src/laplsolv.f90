program laplsolv
    use omp_lib

    !-----------------------------------------------------------------------
    ! Serial program for solving the heat conduction problem
    ! on a square using the Jacobi method.
    ! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
    ! Modified by Berkant Savas (besav@math.liu.se) April 2006
    !-----------------------------------------------------------------------
    integer, parameter                      :: n=1000, maxiter=1000, nr_threads=4
    double precision,parameter              :: tol=1.0E-3
    double precision,dimension(0:n+1,0:n+1) :: T
    double precision,dimension(n)           :: tmp1,tmp2
    double precision                        :: error,x
    real                                    :: t1,t0
    integer                                 :: i,j,k,chunk_size
    character(len=20)                       :: str
    double precision,dimension(0:n-1,0:nr_threads-1) :: padding

    ! Set boundary conditions and initial values for the unknowns
    T=0.0D0
    T(0:n+1 , 0)     = 1.0D0
    T(0:n+1 , n+1)   = 1.0D0
    T(n+1   , 0:n+1) = 2.0D0

    chunk_size = n/nr_threads

    ! Solve the linear system of equations using the Jacobi method
    call cpu_time(t0)

    do k=1,maxiter

        tmp1=T(1:n,0)
        error=0.0D0

        do i=0,nr_threads-1
            padding(0:n-1,i) = T(1:n,(i+1)*chunk_size+1)
        end do

        !$omp parallel do schedule(STATIC,chunk_size) private(j,tmp1,tmp2) firstprivate(padding) shared(T) reduction(max: error)
        do j=1,n
            tmp2=T(1:n,j) !! FEEEEEEL
            if((.NOT. j==chunk_size*(OMP_GET_THREAD_NUM()+1)) .OR. j==n) then
                T(1:n,j) =(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+tmp1)/4.0D0
            else
                T(1:n,j) =(T(0:n-1,j)+T(2:n+1,j)+padding(0:n-1,OMP_GET_THREAD_NUM())+tmp1)/4.0D0
            end if
            error=max(error,maxval(abs(tmp2-T(1:n,j))))
            tmp1=tmp2
        end do
        !$omp end parallel do

        if (error<tol) then
            write(unit=*,fmt=*) 'k = ',k,' , threadnr = ',OMP_GET_THREAD_NUM()
            write(unit=*,fmt=*) 'error = ',error
            exit
        end if

    end do

    call cpu_time(t1)

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
