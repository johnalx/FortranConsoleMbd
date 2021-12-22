
    program FortranConsoleMbd
    use mod_rigidbody
    implicit none

    ! Variables
    type(world) :: sim
    integer :: i,n, steps, iter
    real(wp) :: theta, h, t_end,fps, time
    integer(li) :: tic,toc,rate

    ! Body of FortranConsoleMbd
    n = 6
    sim = world(n, body_cylinder(m=0.05_wp, r=0.02_wp, h=0.05_wp))
    do i=1, n
        theta = (pi*(i-1))/(n-1)
        sim%current(i)%ori = rot(i_, pi/2) .o. rot(k_, theta)
        sim%current(i)%vee = 10.0_wp * j_
        sim%current(i)%omg = 50*pi*i_ - 3*pi*k_
    end do
    
    t_end = 1.0_wp
    steps = 2**21
    h = (t_end - sim%time)/steps
    iter = 0
    call SYSTEM_CLOCK(tic, rate)
    call show_momentum_head()
    call show_momentum(sim ,1)
    do while(sim%time < t_end)
        iter = iter + 1
        call sim%integrate(h)        
        if( mod(iter, 2**16)==0) then
            call show_momentum(sim ,1)
        end if
    end do
    call SYSTEM_CLOCK(toc,rate)
    time = dble(toc-tic)/rate
    fps = (iter*n)/time
    print '(a,i0,a,f0.1,a,f0.3)', 'iter=', iter, ', kfps=', fps/1000, ', time=', time
    ! iter=2097152, kfps=1330.5, time=9.457         <seq> - Native
    ! iter=2097152, kfps=675.0, time=18.641         <omp> - OpenMP parallelism
    ! iter=2097152, kfps=1842.3, time=6.830         <Qipo> - Interprocedureal Optimization
    contains
    
    subroutine show_momentum_head()
        print '(a7,1x,3(a13,1x),1x,3(a13,1x))', "TIME", "L_X", "L_Y", "L_Z", "H_X", "H_Y", "H_Z"
    end subroutine
    subroutine show_momentum(sim, k)
    type(world), intent(in) :: sim
    integer, intent(in) :: k
    type(loading), allocatable :: mom(:)
    type(loading) :: H
        mom = sim%momentum(sim%current)
        H = mom(k)
        print '(f7.4,1x,3(g13.3,1x),1x,3(g13.3,1x))', sim%time, H%force, H%torque
    end subroutine
    
    subroutine show_pos_head()
        print '(a7,1x,3(a13,1x),1x,3(a13,1x))', "TIME", "P_X", "P_Y", "P_Z", "V_X", "V_Y", "V_Z"
    end subroutine
    subroutine show_pos(sim,k)
    type(world), intent(in) :: sim
    integer, intent(in) :: k    
        print '(f7.4,1x,3(f13.3,1x),1x,3(f13.3,1x))', sim%time, sim%current(k)%pos, sim%current(k)%vee
    end subroutine
    subroutine show_ori_head()
        print '(a7,1x,4(a13,1x),1x,3(a13,1x))', "TIME", "Q_X", "Q_Y", "Q_Z", "Q_W", "W_X", "W_Y", "W_Z"
    end subroutine
    subroutine show_ori(sim, k)
    type(world), intent(in) :: sim
    integer, intent(in) :: k    
        print '(f7.4,1x,4(f13.3,1x),1x,3(f13.3,1x))', sim%time, sim%current(k)%ori, sim%current(k)%omg
    end subroutine
    
    subroutine show_mag_head()
        print '(a7,1x,*(a13,1x))', "TIME", "POS", "ORI", "VEE", "OMG"
    end subroutine
    subroutine show_mag(sim, k)
    type(world), intent(in) :: sim
    integer, intent(in) :: k    
    type(state) :: st
        st = sim%current(k)
        print '(f7.4,1x,*(f13.3,1x))', sim%time, norm2(st%pos), norm2(st%ori), norm2(st%vee), norm2(st%omg)
    end subroutine

    subroutine helloworld_omp()
    use omp_lib
    integer :: id, threads
    !$omp parallel
    threads = omp_get_num_threads()
    id = omp_get_thread_num()

    print *, 'Hello World', id, '/', threads

    !$omp end parallel 
    end subroutine
    
    function matvec(A,x) result(y)
    use omp_lib
    real(wp), intent(in) :: A(:,:), x(:)
    real(wp), allocatable :: y(:)
    integer :: i,k, n,m
    integer :: id, threads
    reaL(wp) :: tmp
        n = size(A,1)
        m = size(A,2)
        allocate(y(n))
        !$omp parallel private(id)
        threads = omp_get_num_threads()
        id = omp_get_thread_num()
        print *, 'Start ', id, '/', threads
        !$omp do private(i,k) reduction(+:tmp)
        do i=1, n
            tmp = 0.0_wp
            do k=1, m
                tmp = tmp + A(i,k) * x(k)
            end do
            !$omp critical(dosum)
            y(i) = tmp
            !$omp end critical(dosum)
        end do
        !$omp end do
        print *, 'End ', id, '/', threads
        !$omp end parallel
    end function

    end program FortranConsoleMbd

