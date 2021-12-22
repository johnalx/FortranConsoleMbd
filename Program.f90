
    program FortranConsoleMbd
    use mod_rigidbody
    implicit none

    ! Variables
    type(world) :: sim
    integer :: i,n, steps, iter
    real(wp) :: theta, h, t_end,fps, time
    integer(li) :: tic,toc,rate
    type(rigidbody):: rb

    ! Body of FortranConsoleMbd
    n = 6
    rb = body_cylinder(m=0.05_wp, r=0.02_wp, h=0.05_wp)
    sim = world(n, rb)
    do i=1, n
        theta = (pi*(i-1))/(n-1)
        sim%current(i)%ori = rot(i_, pi/2) .o. rot(k_, theta)
        call rb%set_motion(sim%current(i), motion(10.0_wp * j_, 50*pi*i_ - 3*pi*k_))
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
    ! iter=2097152, kfps=1574.6, time=7.991         <mom> - Momentum state
    ! iter=2097152, kfps=1885.4, time=6.674         <Qipo> - Interprocedureal Optimization
    contains
    
    subroutine show_momentum_head()
        print '(a7,1x,3(a13,1x),1x,3(a13,1x))', "TIME", "L_X", "L_Y", "L_Z", "H_X", "H_Y", "H_Z"
    end subroutine
    subroutine show_momentum(sim, k)
    type(world), intent(in) :: sim
    integer, intent(in) :: k
    type(loading) :: H
        H = loading(sim%current(k)%mom, sim%current(k)%agl)
        print '(f7.4,1x,3(g13.3,1x),1x,3(g13.3,1x))', sim%time, H%force, H%torque
    end subroutine
    
    subroutine show_pos_head()
        print '(a7,1x,3(a13,1x),1x,3(a13,1x))', "TIME", "P_X", "P_Y", "P_Z", "V_X", "V_Y", "V_Z"
    end subroutine
    subroutine show_pos(sim,k)
    type(world), intent(in) :: sim
    integer, intent(in) :: k    
    type(state) :: st
    type(motion) :: v
        st = sim%current(k)
        v = sim%bodies(k)%motion(st)
        print '(f7.4,1x,3(f13.3,1x),1x,3(f13.3,1x))', sim%time, st%pos, v%vee
    end subroutine
    subroutine show_ori_head()
        print '(a7,1x,4(a13,1x),1x,3(a13,1x))', "TIME", "Q_X", "Q_Y", "Q_Z", "Q_W", "W_X", "W_Y", "W_Z"
    end subroutine
    subroutine show_ori(sim, k)
    type(world), intent(in) :: sim
    integer, intent(in) :: k    
    type(state) :: st
    type(motion) :: v
        st = sim%current(k)
        v = sim%bodies(k)%motion(st)
        print '(f7.4,1x,4(f13.3,1x),1x,3(f13.3,1x))', sim%time, st%ori, v%omg
    end subroutine
    
    subroutine show_mag_head()
        print '(a7,1x,*(a13,1x))', "TIME", "POS", "ORI", "VEE", "OMG"
    end subroutine
    subroutine show_mag(sim, k)
    type(world), intent(in) :: sim
    integer, intent(in) :: k    
    type(state) :: st
    type(motion) :: v
        st = sim%current(k)
        v = sim%bodies(k)%motion(st)
        print '(f7.4,1x,*(f13.3,1x))', sim%time, norm2(st%pos), norm2(st%ori), norm2(v%vee), norm2(v%omg)
    end subroutine
    

    end program FortranConsoleMbd

