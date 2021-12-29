    module mod_quaternions
    use mod_vectors
    implicit none (type, external)
    
    enum, bind(c)
        enumerator :: vector_scalar_layout
        enumerator :: scalar_vector_layout
    end enum
    
    type quaternion
        real(wp) :: data(4)
    contains
        procedure :: layout => q_layout
        procedure :: vector => q_vector
        procedure :: scalar => q_scalar
        procedure :: unit => q_unit
        procedure :: sum_sq => q_sum_sq
        procedure :: mag => q_mag
        procedure :: axis => q_axis
        procedure :: angle => q_angle
        procedure :: rot => q_rot_matrix
        procedure, pass :: write => q_write
        procedure, pass :: read => q_read
        generic, public :: write(formatted) => write
        generic, public :: read(formatted)  => read        
    end type

    type(quaternion), parameter :: q_zero = quaternion([0d0,0d0,0d0,0d0])
    type(quaternion), parameter :: q_eye = quaternion([0d0,0d0,0d0,1d0])

    interface rot
        procedure :: q_rot_matrix, q_axis_angle, q_rot_vector
    end interface
        
    interface inv
        procedure :: q_inv
    end interface
    
    interface operator (.o.)
        procedure :: q_mul
    end interface
    
    interface dot_product
        procedure :: q_dot
    end interface
    
    interface norm2
        procedure :: q_mag
    end interface
    
    interface matmul
        procedure :: q_mul_matrix_left !, q_mul_matrix_right        
    end interface
    
    interface assignment (=)
        procedure asgn_quat_from_array, asgn_array_from_quat
    end interface    
    
    interface operator (+)
        procedure :: q_add
    end interface
    interface operator (-)
        procedure :: q_neg, q_sub
    end interface
    interface operator (*)
        procedure :: q_scale_left, q_scale_right
    end interface
    interface operator (/)
        procedure :: q_div
    end interface
    
    interface show
        procedure :: q_show
    end interface

    contains

    !-- QUATERNIONS
    !
    pure subroutine asgn_quat_from_array(q, a)
    type(quaternion), intent(out) :: q
    real(wp), intent(in) :: a(4)
        q = quaternion( a )
    end subroutine
    
    pure subroutine asgn_array_from_quat(a, q)
    real(wp), intent(out) :: a(4)
    type(quaternion), intent(in) :: q
        a = q%data
    end subroutine    
    
    pure function q_layout(q) result(u)
    class(quaternion), intent(in) :: q
    integer :: u
        u = vector_scalar_layout
    end function
    
    pure function q_angle(q) result(u)
    class(quaternion), intent(in) :: q
    real(wp) :: u, c, s
        c = q%data(4)
        s = norm2(q%data(1:3))
        u = 2d0*atan(s/c)
    end function

    pure function q_axis(q) result(k)
    class(quaternion), intent(in) :: q
    real(wp) :: k(3), m, s2
        k = q%data(1:3)
        s2 = q%data(4)**2
        if(s2<1.0_wp-1.0e-9_wp) then
            m = sqrt(1.0_wp - s2)
        else
            m = norm2(k)
        end if
        k = k/m
    end function
    
    elemental function q_vector(q) result(v)
    class(quaternion), intent(in) :: q
    type(vector) :: v
        v = vector(q%data(1:3))
    end function
    
    elemental function q_scalar(q) result(s)
    class(quaternion), intent(in) :: q
    real(wp) :: s
        s = q%data(4)
    end function
    
    pure function q_axis_angle(axis,angle,layout) result(q)
    type(vector), intent(in) :: axis
    real(wp), intent(in) :: angle
    integer, optional, intent(in) :: layout
    type(quaternion) :: q
        if( present(layout) .and. layout == scalar_vector_layout) then
            q%data = [ cos(angle/2), axis%data*sin(angle/2) ]
        else
            q%data = [ axis%data*sin(angle/2), cos(angle/2) ]
        end if
    end function
    
    pure function q_rot_matrix(q) result(R)
    class(quaternion), intent(in) :: q
    real(wp) :: R(3,3)
    real(wp) :: v(3), s, vx(3,3), vxvx(3,3)
        v = q%data(1:3)
        s = q%data(4)
        vx = a_cross_op(v)
        vxvx = -a_mmoi(v)
        R = eye3 + 2*s*vx + 2*vxvx
    end function
    
    pure function q_rot_vector(q, t) result(u)
    type(quaternion), intent(in) :: q
    real(wp), intent(in) :: t(3)
    real(wp) :: u(3)
    real(wp) :: v(3), s, vxt(3), vxvxt(3)
        v = q%data(1:3)
        s = q%data(4)
        vxt = v .x. t
        vxvxt = v .x. vxt
        u = t + 2*s*vxt + 2*vxvxt
    end function

    pure function q_inv(q) result(qinv)
    type(quaternion), intent(in) :: q
    real(wp) :: qi(4), qm2
    type(quaternion) :: qinv
        qm2 = q%sum_sq()
        qinv%data = [-q%data(1:3)/qm2,q%data(4)/qm2]
    end function
    !
    pure function q_mul(q_1, q_2) result(q_3)
    type(quaternion), intent(in) :: q_1, q_2
    type(quaternion) :: q_3
    real(wp) :: qop(4,4)
        qop(:,1) = [ q_1%data(4), q_1%data(3),-q_1%data(2),-q_1%data(1)]
        qop(:,2) = [-q_1%data(3), q_1%data(4), q_1%data(1),-q_1%data(2)]
        qop(:,3) = [ q_1%data(2),-q_1%data(1), q_1%data(4),-q_1%data(3)]
        qop(:,4) = [ q_1%data(1), q_1%data(2), q_1%data(3), q_1%data(4)]
        q_3 %data = matmul(qop, q_2%data)
    end function
        
    pure function q_der(q, omg) result(qp)
    type(quaternion), intent(in) :: q
    type(vector), intent(in) :: omg
    type(quaternion) :: q_w, qp
        q_w%data = [ omg%data/2, 0d0]
        qp = q_mul(q_w, q)
    end function

    pure function q_omg(q,qp) result(omg)
    type(quaternion), intent(in) :: q, qp
    type(vector) :: omg
    real(wp) :: v(3), vp(3), s, sp
        s = q%data(4)
        v = q%data(1:3)
        sp = qp%data(4)
        vp = qp%data(1:3)
        
        omg%data = 2*(s*vp-sp*v+(v .x. vp))
    end function

    pure function q_step_omg(q,omg,h) result(q_next)
    type(quaternion), intent(in) :: q
    type(vector), intent(in) :: omg
    real(wp), intent(in) :: h
    real(wp) :: u
    type(vector) :: k
    type(quaternion) :: q_w, qp, q_next
    real(wp), parameter :: rot_tol = 1d-4*deg
        u = h*norm2(omg%data)
        if( u>rot_tol) then
            k = (h/u)*omg
            q_w = q_axis_angle(k, u)
            q_next = q_mul(q_w, q)
        else
            qp = q_der(q, omg)
            q_next = q + h*qp
        end if
    end function

    pure function q_step_qp(q,qp,h) result(q_next)
    type(quaternion), intent(in) :: q, qp
    real(wp), intent(in) :: h
    real(wp) :: u, eps
    type(vector) :: k, omg
    type(quaternion) :: q_w, q_next
    real(wp), parameter :: rot_tol = 1d-4*deg
        eps = q%sum_sq()-1d0
        if( eps>rot_tol ) then
            omg = q_omg(q, qp)
            u = norm2(omg%data)
            k = omg/u
            q_w = q_axis_angle(k, u)
            q_next = q_mul(q_w, q)
        else
            q_next = q + h*qp
        end if
    end function
        
    elemental function q_dot(q_1,q_2) result(m)
    type(quaternion), intent(in) :: q_1, q_2
    real(wp) :: m
        m = dot_product(q_1%data, q_2%data)
    end function
    elemental function q_sum_sq(q) result(m2)
    class(quaternion), intent(in) :: q
    real(wp) :: m2
        m2 = q%data(1)**2 + q%data(2)**2 + q%data(3)**2 + q%data(4)**2
    end function
    elemental function q_mag(q) result(m)
    class(quaternion), intent(in) :: q
    real(wp) :: m
        m = sqrt(q_sum_sq(q))
    end function
    elemental function q_unit(q) result(q_u)
    class(quaternion), intent(in) :: q
    type(quaternion) :: q_u
    real(wp) :: m2
        m2 = q_sum_sq(q)
        if( m2 >= 1.0e-11_wp ) then
            q_u = q/sqrt(m2)
        else
            q_u = q
        end if
    end function
    elemental function q_add(q_1,q_2) result(q_3)
    type(quaternion), intent(in) :: q_1, q_2
    type(quaternion) :: q_3
        q_3 = quaternion( q_1%data + q_2%data )
    end function
    elemental function q_sub(q_1,q_2) result(q_3)
    type(quaternion), intent(in) :: q_1, q_2
    type(quaternion) :: q_3
        q_3 = quaternion( q_1%data - q_2%data )
    end function
    elemental function q_neg(q_1) result(q_3)
    type(quaternion), intent(in) :: q_1
    type(quaternion) :: q_3
        q_3 = quaternion( -q_1%data )
    end function
    elemental function q_scale_left(f, q_1) result(q_3)
    real(wp), intent(in) :: f
    type(quaternion), intent(in) :: q_1
    type(quaternion) :: q_3
        q_3 = quaternion( f*q_1%data )
    end function
    elemental function q_scale_right(q_1, f) result(q_3)
    type(quaternion), intent(in) :: q_1
    real(wp), intent(in) :: f
    type(quaternion) :: q_3
        q_3 = quaternion( q_1%data*f )
    end function
    elemental function q_div(q_1, f) result(q_3)
    type(quaternion), intent(in) :: q_1
    real(wp), intent(in) :: f
    type(quaternion) :: q_3
        q_3 = quaternion( q_1%data/f )
    end function
    
    pure function q_mul_matrix_left(A, v) result(w)
    real(wp), intent(in) :: A(4,4)
    type(quaternion), intent(in) :: v
    type(quaternion) :: w
        w = quaternion( matmul(A, v%data ) )
    end function
        
    subroutine q_write (q, unit, iotype, v_list, iostat, iomsg)
    class(quaternion), intent(in) :: q
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in)  :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    character(len=:), allocatable :: fmt
        if( iotype == 'LISTDIRECTED' ) then
            write (unit, *, iostat=iostat) q%data
        else
            fmt = '(a,' // iotype(3:) // ',a,' // iotype(3:) // ',a,' // iotype(3:) // ',a,' // iotype(3:) // ',a)'
            write (unit, fmt, iostat=iostat) "<",q%data(1),", ",q%data(2),", ",q%data(3),": ",q%data(4),">"
        end if
    end subroutine

    subroutine q_read (q, unit, iotype, v_list, iostat, iomsg)
    class(quaternion), intent(inout) :: q
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in)  :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    character(len=:), allocatable :: fmt
        read (unit, *, iostat=iostat) q%data
    end subroutine
    
    subroutine q_show(label,q,fmt)
    character(len=*), intent(in) :: label
    type(quaternion), intent(in) :: q
    character(len=*), optional, intent(in) :: fmt
        if(present(fmt)) then
            print '(a,DT "' // fmt // '")', label, q
        else
            print '(a,DT "g0")', label, q
        end if
    end subroutine
    
    end module