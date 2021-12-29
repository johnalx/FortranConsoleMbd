    module mod_rigidbody
    use mod_quaternions
    implicit none (type, external)
    
    type rigidbody
        real(wp) :: mass, mmoi(3)        
    contains
        procedure :: weight => body_get_weight
        procedure :: inertia => body_get_inertia
        procedure :: motion => body_get_motion
        procedure :: set_motion => body_set_motion
        procedure :: rate => body_calc_rate
        procedure :: ke => body_calc_kinetic_energy
        procedure, pass :: write => rb_write
        procedure, pass :: read => rb_read
        generic, public :: write(formatted) => write
        generic, public :: read(formatted)  => read
    end type rigidbody
    
    type state
        type(vector) :: pos
        type(quaternion) :: ori
        type(vector) :: mom, agl
    end type state
    
    type motion
        type(vector) :: vee, omg
    end type
        
    type loading
        type(vector) :: force, torque
    end type
    
    type contact
        type(vector) :: direction
        real(wp) :: impulse
        type(vector) :: pos
    end type
    
    type world
        real(wp) :: time
        type(vector) :: gee
        type(rigidbody), allocatable :: bodies(:)
        type(state), allocatable :: current(:)        
    contains
        procedure :: integrate => world_rk_step
        procedure :: calc_rate => world_calc_rate
        procedure :: motion => world_motion_state
        procedure :: ke => world_calc_kinetic_energy
        procedure :: pe => world_calc_potential_energy
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
    
    interface show
        procedure :: rb_show
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
    
    pure function new_world(n,rb, gee) result(w)
    integer, intent(in) :: n
    type(rigidbody), intent(in) :: rb
    type(vector), optional, intent(in) :: gee
    type(world) :: w
    integer :: i
        w%time= 0.0_wp
        if( present(gee) ) then
            w%gee = gee
        else
            w%gee = o_
        end if
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
    real(wp) :: I(3,3), A(3,3), d(3)
        if( present(inverse) .and. inverse) then
            d = 1/rb%mmoi
            !I = reshape( [1/rb%mmoi(1), 0.0_wp, 0.0_wp, &
            !    0.0_wp, 1/rb%mmoi(2), 0.0_wp, &
            !    0.0_wp, 0.0_wp, 1/rb%mmoi(3)], [3,3] )            
        else
            d = rb%mmoi
            !I = reshape( [rb%mmoi(1), 0.0_wp, 0.0_wp, &
            !    0.0_wp, rb%mmoi(2), 0.0_wp, &
            !    0.0_wp, 0.0_wp, rb%mmoi(3)], [3,3] )            
        end if
        ! A = I*tr(rot)
        A(:, 1) = d * rot(1, :)
        A(:, 2) = d * rot(2, :)
        A(:, 3) = d * rot(3, :)
        ! I = matmul(rot, matmul(I, transpose(rot)))
        I = matmul(rot, A)
    end function
    
    pure function body_get_weight(rb, gee) result(w)
    class(rigidbody), intent(in) :: rb
    type(vector), intent(in) :: gee
    type(loading) :: w
        w%force = rb%mass * gee
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
    
    elemental function body_calc_rate(rb, current, f) result(rate)
    class(rigidbody), intent(in) :: rb
    type(state), intent(in) :: current
    type(loading), intent(in) :: f
    type(state) :: rate
    type(motion) :: v
        v = rb%motion(current)
        ! d(r)/dt = v
        rate%pos = v%vee
        ! d(q)/dt = 1/2*ω*q
        rate%ori = q_der(current%ori, v%omg)            
        ! d(p)/dt = F
        rate%mom = f%force
        ! d(H)/dt = τ
        rate%agl = f%torque        
    end function
    
    elemental function body_calc_kinetic_energy(rb, current) result(ke)
    class(rigidbody), intent(in) :: rb
    type(state), intent(in) :: current
    real(wp) :: ke
    type(motion) :: v
        v = rb%motion(current)
        ke = 0.5_wp * ( dot_product(v%vee, current%mom) + dot_product(v%omg, current%agl) )
    end function
    
    pure function world_calc_kinetic_energy(self) result(ke)
    class(world), intent(in) :: self
    real(wp) :: ke
    !type(motion) :: v
    integer :: k, n
        n = size(self%bodies)
        ke = 0.0_wp
        do k=1, n
            !v = self%bodies(k)%motion(self%current(k))
            !ke = ke + 0.5_wp * ( dot_product(v%vee, self%current(k)%mom) + dot_product(v%omg, self%current(k)%agl) )
            ke = ke +  self%bodies(k)%ke( self%current(k) )
        end do
        ! ke = sum( self%bodies%ke( self%current ) )
    end function
    
    pure function world_calc_potential_energy(self) result(pe)
    class(world), intent(in) :: self
    real(wp) :: pe
    type(loading) :: fa
    integer :: k, n
        n = size(self%bodies)
        pe = 0.0_wp
        do k=1, n
            fa = self%bodies(k)%weight(self%gee)
            pe = pe + dot_product( self%current(k)%pos, fa%force)
        end do
    end function
    
    pure function world_calc_rate(self, current) result(rate)
    class(world), intent(in) :: self
    type(state), intent(in), allocatable :: current(:)
    type(state), allocatable :: rate(:)
    type(rigidbody) :: rb
    type(loading), allocatable :: fa(:)
    integer :: k, n     
    
        n = size(self%bodies)        
        allocate(rate(n))
        allocate(fa(n))
        
        do concurrent (k=1:n)
            rb = self%bodies(k)
            fa(k) = rb%weight(self%gee)                   
            rate(k) = rb%rate(current(k), fa(k))
        end do  
        
        !rate = self%bodies%rate(current, fa)        
            
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
            
    pure function world_motion_state(self, current) result(v)
    class(world), intent(in) :: self
    type(state), intent(in), allocatable :: current(:)
    type(motion), allocatable :: v(:)
    integer :: n
        n = size(self%bodies)
        v = self%bodies%motion(current)
    end function
    
    subroutine rb_write (rb, unit, iotype, v_list, iostat, iomsg)
    class(rigidbody), intent(in) :: rb
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in)  :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    character(len=:), allocatable :: fmt
        if( iotype == 'LISTDIRECTED' ) then
            write (unit, *, iostat=iostat) rb%mass, rb%mmoi
        else
            fmt = '(a,' // iotype(3:) // ',a,' // iotype(3:) // ',a,' // iotype(3:) // ',a,' // iotype(3:) // ',a)'
            write (unit, fmt, iostat=iostat) "[m=",rb%mass,", I=(",rb%mmoi(1),", ",rb%mmoi(2),", ",rb%mmoi(3),")]"
        end if
    end subroutine

    subroutine rb_read (rb, unit, iotype, v_list, iostat, iomsg)
    class(rigidbody), intent(inout) :: rb
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in)  :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    character(len=:), allocatable :: fmt
        read (unit, *, iostat=iostat) rb%mass, rb%mmoi
    end subroutine
    
    subroutine rb_show(label, rb, fmt)
    character(len=*), intent(in) :: label
    type(rigidbody), intent(in) :: rb
    character(len=*), optional, intent(in) :: fmt
        if(present(fmt)) then
            print '(a,DT "' // fmt // '")', label, rb
        else
            print '(a,DT "g0")', label, rb
        end if
    end subroutine
    
        
    end module