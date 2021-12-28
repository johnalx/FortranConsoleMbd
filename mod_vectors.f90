    module mod_vectors
    use, intrinsic :: iso_fortran_env, only : sp=>real32, wp => real64, li => int64
    implicit none ( type, external )

    integer, parameter :: dims = 3

    enum, bind(c)
        enumerator :: x_axis
        enumerator :: y_axis
        enumerator :: z_axis
    end enum

    type vector
        real(wp) :: data(dims)
    contains
        procedure :: mag => v_mag
        procedure, pass :: write => v_write
        procedure, pass :: read => v_read
        generic, public :: write(formatted) => write
        generic, public :: read(formatted)  => read
    end type

    real(wp), parameter :: nan64 =  transfer(-2251799813685248_li, 1._wp)
    real(wp), parameter :: pi = 4d0*atan(1d0), deg = pi/180d0
    real(wp), parameter :: eye3(3,3) = reshape( &
        [1.0_wp, 0.0_wp, 0.0_wp, &
        0.0_wp, 1.0_wp, 0.0_wp, &
        0.0_wp, 0.0_wp, 1.0_wp],[3,3])
    type(vector), parameter :: o_ = vector([0.0_wp, 0.0_wp, 0.0_wp])
    type(vector), parameter :: i_ = vector([1.0_wp, 0.0_wp, 0.0_wp])
    type(vector), parameter :: j_ = vector([0.0_wp, 1.0_wp, 0.0_wp])
    type(vector), parameter :: k_ = vector([0.0_wp, 0.0_wp, 1.0_wp])

    interface eye
    procedure m_identity, m_elemental, a_elemental
    end interface
    interface null
    procedure m_null
    end interface

    interface norm2
    procedure v_mag
    end interface

    interface dot_product
    procedure v_dot
    end interface

    interface dot
    procedure v_dot
    end interface

    interface cross
    procedure a_cross_op, a_cross, v_cross_op, v_cross
    end interface
    interface operator (.x.)
    procedure a_cross_op, a_cross, v_cross_op, v_cross
    end interface

    interface matmul
    procedure :: v_mul_matrix_left !, v_mul_matrix_right
    end interface

    interface diag
    procedure a_diag, D_diag
    end interface
    interface inv
    procedure m_inv_3
    end interface
    interface solve
    procedure a_solve, a_solve2
    end interface

    interface vector
    procedure xyz_to_vector, axis_to_vector
    end interface

    interface assignment (=)
    procedure asgn_vector_from_array, asgn_array_from_vector
    end interface

    interface operator (+)
    procedure :: v_add
    end interface
    interface operator (-)
    procedure :: v_neg, v_sub
    end interface
    interface operator (*)
    procedure :: v_scale_left, v_scale_right
    end interface
    interface operator (/)
    procedure :: v_div
    end interface
    
    interface show
        procedure :: v_show
    end interface

    contains

    pure function i_null(n, index) result(non)
    integer, intent(in) :: n, index(:)
    integer, allocatable :: non(:)
    integer i, k, nk
    nk = size(index)
    k = 0
    allocate(non(n-nk))
    do i=1, n
        if ( any(index==i) ) then
        else
            k = k + 1
            non(k) = i
        end if
    end do
    end function

    pure function m_identity(n) result(A)
    integer, intent(in) :: n
    real(wp) :: A(n,n)
    integer :: i
    A = 0d0
    forall(i=1:n)
        A(i,i) = 1d0
    end forall
    end function

    pure function a_elemental(n,i) result(v)
    integer, intent(in) :: n, i
    real(wp) :: v(n)
    v = 0d0
    v(i) = 1d0
    end function

    pure function m_elemental(n,index) result(A)
    integer, intent(in) :: n, index(:)
    real(wp), allocatable :: A(:,:)
    integer :: i, nk
    nk = size(index)
    allocate(A(n, nk), source=0d0)
    do i=1, nk
        A(index(i),i) = 1d0
    end do
    end function

    pure function m_null(n, index) result(A)
    integer, intent(in) :: n, index(:)
    real(wp), allocatable:: A(:,:)
    integer, allocatable :: non(:)
    integer :: i, nk
    nk = size(index)
    non = i_null(n, index)
    allocate(A(n,n-nk), source=0d0)
    do i=1,n-nk
        A(non(i),i) = 1d0
    end do
    end function

    pure function a_cross(v_1, v_2) result(v_3)
    real(wp), intent(in) :: v_1(3), v_2(3)
    real(wp) :: v_3(3)
    v_3 = [ v_1(2)*v_2(3)-v_1(3)*v_2(2), &
        v_1(3)*v_2(1)-v_1(1)*v_2(3), &
        v_1(1)*v_2(2)-v_1(2)*v_2(1) ]
    end function

    pure function a_cross_op(v) result(vx)
    real(wp), intent(in) :: v(3)
    real(wp) :: vx(3,3)
    vx(1,:) = [0d0, -v(3), v(2)]
    vx(2,:) = [v(3), 0d0, -v(1)]
    vx(3,:) = [-v(2), v(1), 0d0]
    end function

    pure function a_mmoi(v) result(mm)
    real(wp), intent(in) :: v(3)
    real(wp) :: mm(3,3)
    mm(1,:) = [v(2)**2+v(3)**2, -v(1)*v(2), -v(1)*v(3)]
    mm(2,:) = [-v(1)*v(2), v(1)**2+v(3)**2, -v(2)*v(3)]
    mm(3,:) = [-v(1)*v(3), -v(2)*v(3), v(1)**2+v(2)**2]
    end function

    pure function a_mag(v) result(m)
    real(wp), intent(in) :: v(:)
    real(wp) :: m
    m = sqrt( sum( v**2 ))
    end function

    pure function a_dot(v, w) result(m)
    real(wp), intent(in) :: v(:), w(:)
    real(wp) :: m
    m = sum(v*w)
    end function

    pure function a_diag(v) result(D)
    real(wp), intent(in) :: v(:)
    real(wp),allocatable :: D(:,:)
    integer :: n, i
    n = size(v)
    allocate(D(n,n))
    D = 0d0
    forall(i=1:n)
        D(i,i) = v(i)
    end forall
    end function

    pure function D_diag(D) result(v)
    real(wp), intent(in) :: D(:,:)
    real(wp), allocatable :: v(:)
    integer :: n,i
    n = min(size(D,1),size(D,2))
    allocate(v(n))
    do i=1,n
        v(i) = D(i,i)
    end do
    end function

    pure function a_matmul(A, v) result(w)
    ! evaluate w = A*v
    real(wp), intent(in) :: A(:,:), v(:)
    real(wp) :: w(size(A,1))
    integer :: i
    do i=1,size(A,1)
        w(i) = sum(A(i,:)*v(:))
    end do
    end function

    pure function a_matmul2(v, A) result(w)
    ! evaluate w = v*A
    real(wp), intent(in) :: A(:,:), v(:)
    real(wp) :: w(size(A,2))
    integer :: i
    do i=1,size(A,2)
        w(i) = sum(v(:)*A(:,i))
    end do
    end function

    pure function m_inv_3(A) result(B)
    ! invert 3x3 matrix
    real(wp), intent(in) :: A(3,3)
    real(wp) :: B(3,3), d, t2, t3, t4, t7, t8, t9, t6
    t2 = A(1,1)*A(2,2)*A(3,3)
    t3 = A(1,2)*A(2,3)*A(3,1)
    t4 = A(1,3)*A(2,1)*A(3,2)
    t7 = A(1,1)*A(2,3)*A(3,2)
    t8 = A(1,2)*A(2,1)*A(3,3)
    t9 = A(1,3)*A(2,2)*A(3,1)
    d = t2+t3+t4-t7-t8-t9
    t6 = 1.0D0/d
    B(1,1) = t6*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    B(1,2) = -t6*(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    B(1,3) = t6*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    B(2,1) = -t6*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    B(2,2) = t6*(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    B(2,3) = -t6*(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    B(3,1) = t6*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    B(3,2) = -t6*(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    B(3,3) = t6*(A(1,1)*A(2,2)-A(1,2)*A(2,1))
    end function

    pure function a_solve(A, b) result(x)
    ! solve b = A*x for x
    real(wp), intent(in) :: A(3,3), b(3)
    real(wp) :: x(3), d, w(3)

    d = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    w = b/d
    x = [ &
        w(1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+w(2)*(A(1,3)*A(3,2)-A(1,2)*A(3,3))+w(3)*(A(1,2)*A(2,3)-A(1,3)*A(2,2)), &
        w(1)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+w(2)*(A(1,1)*A(3,3)-A(1,3)*A(3,1))+w(3)*(A(1,3)*A(2,1)-A(1,1)*A(2,3)), &
        w(1)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))+w(2)*(A(1,2)*A(3,1)-A(1,1)*A(3,2))+w(3)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))]
    end function

    pure function a_solve2(A, B) result(x)
    ! solve b = A*x for x
    real(wp), intent(in) :: A(3,3), B(3,3)
    real(wp) :: x(3,3), d, w(3,3)
    integer :: i

    d = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

    forall (i=1:3)
        w(:,i) = B(:,i)/d
        x(:,i) = [ &
            w(1,i)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+w(2,i)*(A(1,3)*A(3,2)-A(1,2)*A(3,3))+w(3,i)*(A(1,2)*A(2,3)-A(1,3)*A(2,2)), &
            w(1,i)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+w(2,i)*(A(1,1)*A(3,3)-A(1,3)*A(3,1))+w(3,i)*(A(1,3)*A(2,1)-A(1,1)*A(2,3)), &
            w(1,i)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))+w(2,i)*(A(1,2)*A(3,1)-A(1,1)*A(3,2))+w(3,i)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))]
    end forall

    end function

    function a_rotate_x(v,q) result(u)
    real(wp), intent(in) :: v(3), q
    real(wp) :: u(3), cq, sq
    cq = cos(q)
    sq = sin(q)
    u = [ v(1), &
        cq*v(2)-sq*v(3), &
        sq*v(2)+cq*v(3)]
    end function
    function a_rotate_y(v,q) result(u)
    real(wp), intent(in) :: v(3), q
    real(wp) :: u(3), cq, sq
    cq = cos(q)
    sq = sin(q)
    u = [ cq*v(1)+sq*v(3), &
        v(2), &
        -sq*v(1)+cq*v(3)]
    end function
    function a_rotate_z(v,q) result(u)
    real(wp), intent(in) :: v(3), q
    real(wp) :: u(3), cq, sq
    cq = cos(q)
    sq = sin(q)
    u = [ cq*v(1)-sq*v(1), &
        sq*v(1)+cq*v(2), &
        v(3)]
    end function

    function m_mmul(A,x,b) result(y)
    !tex: Matrix vector calculation
    ! $$ y = A\,x + b $$
    ! or
    ! $$ y = A\,x  $$
    real(wp), intent(in) :: A(:,:), x(:)
    real(wp), intent(in), optional :: b(:)
    real(wp), allocatable :: y(:)
    if( present(b) ) then
        y = matmul(A,x)  + b
    else
        y = matmul(A,x)
    end if
    end function

    pure function axis_to_vector(axis) result(v)
    integer, intent(in) :: axis
    type(vector) :: v
    select case(axis)
    case(x_axis)
        v = i_
    case(y_axis)
        v = j_
    case(z_axis)
        v = k_
        case default
        v = o_
    end select
    end function

    pure function xyz_to_vector(x,y,z) result(v)
    real(wp), intent(in) :: x,y,z
    type(vector) :: v
    v%data = [x,y,z]
    end function

    pure subroutine asgn_vector_from_array(v, a)
    type(vector), intent(out) :: v
    real(wp), intent(in) :: a(dims)
    v%data = a
    end subroutine

    pure subroutine asgn_array_from_vector(a, v)
    real(wp), intent(out) :: a(dims)
    type(vector), intent(in) :: v
    a = v%data
    end subroutine

    elemental function v_mag(v) result(m)
    class(vector), intent(in) :: v
    real(wp) :: m
    m = norm2(v%data)
    end function

    elemental function v_dot(u,v) result(m)
    class(vector), intent(in) :: u,v
    real(wp) :: m
    m = dot_product(u%data, v%data)
    end function

    elemental function v_cross(u, v) result(w)
    class(vector), intent(in) :: u, v
    type(vector) :: w
    w%data = [ &
        u%data(2)*v%data(3) - u%data(3)*v%data(2), &
        u%data(3)*v%data(1) - u%data(1)*v%data(2), &
        u%data(1)*v%data(2) - u%data(2)*v%data(1) ]
    end function

    pure function v_cross_op(v) result(w)
    class(vector), intent(in) :: v
    real(wp) :: w(3,3)
    !     |  0 -z  y |
    ! w = |  z  0 -x |
    !     | -y  x  0 |
    w = reshape([ &
        0.0_wp, v%data(3), -v%data(2), &
        -v%data(1), 0.0_wp, v%data(1), &
        v%data(2), -v%data(1), 0.0_wp], [3,3])
    end function

    elemental function v_scale_left(f, v) result(w)
    real(wp), intent(in) :: f
    class(vector), intent(in) :: v
    type(vector) :: w
    w%data = f*v%data
    end function

    elemental function v_scale_right(v, f) result(w)
    class(vector), intent(in) :: v
    real(wp), intent(in) :: f
    type(vector) :: w
    w%data = f*v%data
    end function

    elemental function v_neg(v) result(w)
    class(vector), intent(in) :: v
    type(vector) :: w
    w%data = -v%data
    end function

    elemental function v_add(v, u) result(w)
    class(vector), intent(in) :: v, u
    type(vector) :: w
    w%data = v%data + u%data
    end function

    elemental function v_sub(v, u) result(w)
    class(vector), intent(in) :: v, u
    type(vector) :: w
    w%data = v%data + u%data
    end function

    elemental function v_div(v, f) result(w)
    class(vector), intent(in) :: v
    real(wp), intent(in) :: f
    type(vector) :: w
    w%data = v%data/f
    end function

    pure function v_mul_matrix_left(A, v) result(w)
    real(wp), intent(in) :: A(dims,dims)
    type(vector), intent(in) :: v
    type(vector) :: w
    w%data = matmul(A, v%data)
    end function

    subroutine v_write (v, unit, iotype, v_list, iostat, iomsg)
    class(vector), intent(in) :: v
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in)  :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    character(len=:), allocatable :: fmt
        if( iotype == 'LISTDIRECTED' ) then
            write (unit, *, iostat=iostat) v%data
        else
            fmt = '(a,' // iotype(3:) // ',a,' // iotype(3:) // ',a,' // iotype(3:) // ',a)'
            write (unit, fmt, iostat=iostat) "(",v%data(1),", ",v%data(2),", ",v%data(3),")"
        end if
    end subroutine

    subroutine v_read (v, unit, iotype, v_list, iostat, iomsg)
    class(vector), intent(inout) :: v
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in)  :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    character(len=:), allocatable :: fmt
        read (unit, *, iostat=iostat) v%data
    end subroutine
    
    subroutine v_show(label,vec,fmt)
    character(len=*), intent(in) :: label
    type(vector), intent(in) :: vec
    character(len=*), optional, intent(in) :: fmt
        if(present(fmt)) then
            print '(a,DT "' // fmt // '")', label, vec
        else
            print '(a,DT "g0")', label, vec
        end if
    end subroutine
    

    end module