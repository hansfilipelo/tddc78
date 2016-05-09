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
    double precision,dimension(n)           :: tmp, padding_before, padding_after
    double precision                        :: error,x
    double precision                        :: t1,t0
    integer                                 :: i,j,k,chunk_size
    character(len=20)                       :: str

    ! Set boundary conditions and initial values for the unknowns
    T=0.0D0
    T(0:n+1 , 0)     = 1.0D0
    T(0:n+1 , n+1)   = 1.0D0
    T(n+1   , 0:n+1) = 2.0D0

    chunk_size = n/nr_threads

    ! Solve the linear system of equations using the Jacobi method
    t0 = omp_get_wtime()

    do k=1,maxiter

        error=0.0D0

        !$omp parallel private(j,tmp1,tmp,padding_before,padding_after) shared(T) reduction(max: error)
        padding_before = T(1:n,OMP_GET_THREAD_NUM()*chunk_size)
        padding_after = T(1:n,(OMP_GET_THREAD_NUM()+1)*chunk_size+1)

        !$omp do schedule(STATIC,chunk_size)
        do j=1,n
            tmp=T(1:n,j)
            if((.NOT. j==chunk_size*(OMP_GET_THREAD_NUM()+1))) then
                T(1:n,j) =(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+padding_before)/4.0D0
            else
                T(1:n,j) =(T(0:n-1,j)+T(2:n+1,j)+padding_after+padding_before)/4.0D0
            end if

            error=max(error,maxval(abs(tmp-T(1:n,j))))
            padding_before=tmp
        end do
        !$omp end do
        !$omp end parallel

        if (error<tol) then
            write(unit=*,fmt=*) 'k = ',k,' , threadnr = ',OMP_GET_THREAD_NUM()
            write(unit=*,fmt=*) 'error = ',error
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
