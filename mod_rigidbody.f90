    module mod_rigidbody
    use mod_quaternions
    implicit none (type, external)
    
    real(wp), parameter :: gravity(3) = [0d0,-10d0,0d0]  
    
    type rigidbody
        real(wp) :: mass, mmoi(3)        
    contains
        procedure :: inertia => body_calc_inertia
        procedure :: momentum => world_momentum_body
    end type rigidbody
    
    type state
        real(wp) :: pos(3), vee(3)
        real(wp) :: ori(4), omg(3)
    end type state
    
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
        procedure :: momentum => world_momentum_state
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
    
    pure function body_calc_inertia(rb,rot,inverse) result(I)
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
    
    elemental function world_momentum_body(rb, current) result(H)
    class(rigidbody), intent(in) :: rb
    type(state), intent(in) :: current
    type(loading) :: H
    real(wp) :: R(3,3), I(3,3)
        R = rot( current%ori )
        I = rb%inertia(R)
        H%force = rb%mass * current%vee
        H%torque = matmul(I, current%omg )        
    end function
    
    elemental function state_neg(s) result(r)
    type(state), intent(in) :: s
    type(state) :: r
        r%pos = -s%pos
        r%vee = -s%vee
        r%ori = -s%ori
        r%omg = -s%omg
    end function
    elemental function state_add(g,s) result(r)
    type(state), intent(in) :: g, s
    type(state) :: r
        r%pos = g%pos + s%pos
        r%vee = g%vee + s%vee
        r%ori = g%ori + s%ori
        r%omg = g%omg + s%omg
    end function
    elemental function state_sub(g,s) result(r)
    type(state), intent(in) :: g, s
    type(state) :: r
        r%pos = g%pos - s%pos
        r%vee = g%vee - s%vee
        r%ori = g%ori - s%ori
        r%omg = g%omg - s%omg
    end function
    elemental function state_scale_left(f,s) result(r)
    real(wp), intent(in) :: f
    type(state), intent(in) :: s
    type(state) :: r
        r%pos = f * s%pos
        r%vee = f * s%vee
        r%ori = f * s%ori
        r%omg = f * s%omg
    end function
    
    elemental function state_scale_right(s,f) result(r)
    real(wp), intent(in) :: f
    type(state), intent(in) :: s
    type(state) :: r
        r%pos = f * s%pos
        r%vee = f * s%vee
        r%ori = f * s%ori
        r%omg = f * s%omg
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
            w%current(i)%vee = o_
            w%current(i)%ori = q_eye
            w%current(i)%omg = o_
        end do
    end function
    
    function world_calc_rate(self, current) result(rate)
    class(world), intent(in) :: self
    type(state), intent(in), allocatable :: current(:)
    type(state), allocatable :: rate(:)
    real(wp) :: R(3,3), I(3,3), M(3,3), H(3), F(3), T(3)
    integer :: k,n     
    
        n = size(self%bodies)        
        allocate(rate(n))        
        !!$omp parallel do private(k,R,I,M,H,F,T)
        do k=1, n
            F = self%bodies(k)%mass * gravity
            T = o_
            rate(k)%pos = current(k)%vee
            rate(k)%ori = q_der(current(k)%ori, current(k)%omg)
            R = rot( current(k)%ori )
            ! I = R*I_body*tr(R)
            I = self%bodies(k)%inertia(R)
            M = self%bodies(k)%inertia(R, .true.)
            ! H = I*ω
            H = matmul(I, current(k)%omg)
            ! a = F/m
            rate(k)%vee = F / self%bodies(k)%mass
            ! α = (τ - ω×H)/I
            rate(k)%omg = matmul(M, T - (current(k)%omg .x. H))
        end do
        !!$omp end parallel do
    end function
    
    subroutine world_rk_step(self, h)
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
        
    pure function world_momentum_state(self, current) result(H)
    class(world), intent(in) :: self
    type(state), intent(in), allocatable :: current(:)
    type(loading), allocatable :: H(:)
    real(wp) :: R(3,3), I(3,3)
    integer :: k, n
        n = size(self%bodies)
        H = self%bodies%momentum(current)
    end function
    
    end module