    module mod_rigidbody
    use mod_quaternions
    implicit none (type, external)
    
    real(wp), parameter :: gravity(3) = [0d0,-10d0,0d0]  
    
    type rigidbody
        real(wp) :: mass, mmoi(3)        
    contains
        procedure :: weight => body_get_weight
        procedure :: inertia => body_get_inertia
        procedure :: motion => body_get_motion
        procedure :: set_motion => body_set_motion
    end type rigidbody
    
    type state
        real(wp) :: pos(3), ori(4)
        real(wp) :: mom(3), agl(3)
    end type state
    
    type motion
        real(wp) :: vee(3), omg(3)
    end type
        
    type loading
        real(wp) :: force(3), torque(3)
    end type
    
    type contact
        real(wp) :: direction(3), impulse
        real(wp) :: pos(3)
    end type
    
    type world
        real(wp) :: time
        type(rigidbody), allocatable :: bodies(:)
        type(state), allocatable :: current(:)
    contains
        procedure :: integrate => world_rk_step
        procedure :: calc_rate => world_calc_rate
        procedure :: motion => world_motion_state
    end type world
    
    interface world
        procedure :: new_world
    end interface
    
    interface operator (+)
        procedure :: state_add
    end interface
    interface operator (-)
        procedure :: state_sub, state_neg
    end interface
    interface operator (*)
        procedure :: state_scale_left, state_scale_right
    end interface

    contains

    pure function body_sphere(m, r) result(rb)
        real(wp), intent(in) :: m, r
        type(rigidbody) :: rb
        rb = body_ellipsoid(m,r,r,r)
    end function
    pure function body_ellipsoid(m, a,b,c) result(rb)
        real(wp), intent(in) :: m, a,b,c
        type(rigidbody) :: rb
        rb = rigidbody(m, m*[2*a**5/2,2*b**2/5,2*c**2/5])
    end function
    
    pure function body_rod(m, l) result(rb)
        real(wp), intent(in) :: m, l
        type(rigidbody) :: rb
        rb = body_cylinder(m,0d0,l)
    end function
    pure function body_disk(m, r) result(rb)
        real(wp), intent(in) :: m, r
        type(rigidbody) :: rb
        rb =  body_cylinder(m,r,0d0)
    end function
    pure function body_cylinder(m, r, h) result(rb)
        real(wp), intent(in) :: m, r,h
        type(rigidbody) :: rb
        rb = rigidbody(m, m*[r**2/4+h**2/12,r**2/4+h**2/12,r**2/2])
    end function
    pure function body_prism(m, a, b, c) result(rb)
        real(wp), intent(in) :: m, a,b,c
        type(rigidbody) :: rb
        rb = rigidbody(m, m*[a**2/12,b**2/12,c**2/12])
    end function
    pure function body_cone(m, r, h) result(rb)
        real(wp), intent(in) :: m, r, h
        type(rigidbody) :: rb
        rb = rigidbody(m, m*[3*h**2/80+3*r**2/20, 3*h**2/80+3*r**2/20, 3*r**2/10])
    end function
    
    
    elemental function state_neg(s) result(r)
    type(state), intent(in) :: s
    type(state) :: r
        r%pos = -s%pos
        r%ori = -s%ori
        r%mom = -s%mom
        r%agl = -s%agl
    end function
    elemental function state_add(g,s) result(r)
    type(state), intent(in) :: g, s
    type(state) :: r
        r%pos = g%pos + s%pos
        r%ori = g%ori + s%ori
        r%mom = g%mom + s%mom
        r%agl = g%agl + s%agl
    end function
    elemental function state_sub(g,s) result(r)
    type(state), intent(in) :: g, s
    type(state) :: r
        r%pos = g%pos - s%pos
        r%ori = g%ori - s%ori
        r%mom = g%mom - s%mom
        r%agl = g%agl - s%agl
    end function
    elemental function state_scale_left(f,s) result(r)
    real(wp), intent(in) :: f
    type(state), intent(in) :: s
    type(state) :: r
        r%pos = f * s%pos
        r%ori = f * s%ori
        r%mom = f * s%mom
        r%agl = f * s%agl
    end function    
    
    elemental function state_scale_right(s,f) result(r)
    real(wp), intent(in) :: f
    type(state), intent(in) :: s
    type(state) :: r
        r%pos = f * s%pos
        r%ori = f * s%ori
        r%mom = f * s%mom
        r%agl = f * s%agl
    end function
    
    pure function new_world(n,rb) result(w)
    integer, intent(in) :: n
    type(rigidbody), intent(in) :: rb
    type(world) :: w
    integer :: i
        w%time= 0.0_wp
        allocate(w%current(n))
        allocate(w%bodies(n))
        do i=1,n
            w%bodies(i) = rb
            w%current(i)%pos = o_
            w%current(i)%ori = q_eye
            w%current(i)%mom = o_
            w%current(i)%agl = o_
        end do
    end function
    
    pure function body_get_inertia(rb,rot,inverse) result(I)
    class(rigidbody), intent(in) :: rb
    real(wp), intent(in) :: rot(3,3)
    logical, optional, intent(in) :: inverse
    real(wp) :: I(3,3)
        if( present(inverse) .and. inverse) then
            I = reshape( [1/rb%mmoi(1), 0.0_wp, 0.0_wp, &
                0.0_wp, 1/rb%mmoi(2), 0.0_wp, &
                0.0_wp, 0.0_wp, 1/rb%mmoi(3)], [3,3] )            
        else
            I = reshape( [rb%mmoi(1), 0.0_wp, 0.0_wp, &
                0.0_wp, rb%mmoi(2), 0.0_wp, &
                0.0_wp, 0.0_wp, rb%mmoi(3)], [3,3] )            
        end if
        I = matmul(rot, matmul(I, transpose(rot)))
    end function
    
    pure function body_get_weight(rb, g) result(w)
    class(rigidbody), intent(in) :: rb
    real(wp), intent(in) :: g(3)
    type(loading) :: w
        w%force = rb%mass * g
        w%torque = o_
    end function
        
    elemental function body_get_motion(rb, current) result(v)
    class(rigidbody), intent(in) :: rb
    type(state), intent(in) :: current
    type(motion) :: v
    real(wp) :: R(3,3), M(3,3)
        R = rot(current%ori)
        M = rb%inertia(R, .true.)
        v%vee = current%mom / rb%mass
        v%omg = matmul(M, current%agl )
    end function
    
    elemental subroutine body_set_motion(rb, current, v)
    class(rigidbody), intent(in) :: rb
    type(state), intent(inout) :: current
    type(motion),intent(in) :: v
    real(wp) :: R(3,3), I(3,3)
        R = rot(current%ori)
        I = rb%inertia(R, .false.)
        current%mom = rb%mass * v%vee
        current%agl = matmul(I, v%omg)
    end subroutine
    
    pure function world_calc_rate(self, current) result(rate)
    class(world), intent(in) :: self
    type(state), intent(in), allocatable :: current(:)
    type(state), allocatable :: rate(:)
    type(loading) :: f
    type(motion) :: v
    type(rigidbody) :: rb
    integer :: k, n     
    
        n = size(self%bodies)        
        allocate(rate(n))        
        do k=1, n
            rb = self%bodies(k)
            f = rb%weight(gravity)            
            v = rb%motion(current(k))
                                    
            rate(k)%pos = v%vee
            rate(k)%ori = q_der(current(k)%ori, v%omg)
            
            ! d(p)/dt = F
            rate(k)%mom = f%force
            ! d(H)/dt = τ
            rate(k)%agl = f%torque
        end do        
    end function
    
    pure subroutine world_rk_step(self, h)
    class(world), intent(inout) :: self
    real(wp), intent(in) :: h
    type(state), allocatable :: next(:), K0(:), K1(:), K2(:), K3(:)
    integer :: n 
        n = size(self%bodies)
        allocate(next(n))
        next = self%current
        K0 = self%calc_rate(next)
        
        next = self%current + h/2 * K0
        K1 = self%calc_rate(next)
        
        next = self%current + h/2 * K1
        K2 = self%calc_rate(next)
        
        next = self%current + h * K2
        K3 = self%calc_rate(next)
        
        self%time = self%time + h
        self%current = self%current + (h/6)*(K0 + 2.0_wp*K1 + 2.0_wp*K2 + K3)
        
    end subroutine
        
    !pure function world_momentum_state(self, current) result(H)
    !class(world), intent(in) :: self
    !type(state), intent(in), allocatable :: current(:)    
    !type(loading), allocatable :: H(:)
    !integer :: n
    !    n = size(self%bodies)
    !    H = self%bodies%momentum(current)
    !end function
    
    pure function world_motion_state(self, current) result(v)
    class(world), intent(in) :: self
    type(state), intent(in), allocatable :: current(:)
    type(motion), allocatable :: v(:)
    integer :: n
        n = size(self%bodies)
        v = self%bodies%motion(current)
    end function
    
    end module