    module mod_vectors
    use, intrinsic :: iso_fortran_env, only : sp=>real32, wp => real64, li => int64
    implicit none ( type, external )
    
    real(wp), parameter :: nan64 =  transfer(-2251799813685248_li, 1._wp)
    real(wp), parameter :: pi = 4d0*atan(1d0), deg = pi/180d0
    reaL(wp), parameter :: o_(3) = [0d0,0d0,0d0]
    reaL(wp), parameter :: i_(3) = [1d0,0d0,0d0]
    reaL(wp), parameter :: j_(3) = [0d0,1d0,0d0]
    reaL(wp), parameter :: k_(3) = [0d0,0d0,1d0]
    real(wp), parameter :: eye3(3,3) = reshape([i_,j_,k_],[3,3])
    
    interface eye
        procedure a_identity, a_elemental, v_elemental
    end interface    
    interface null
        procedure a_null
    end interface
    
    interface mag
        procedure v_mag
    end interface
    !interface dot
    !    procedure v_dot
    !end interface
    interface cross
        procedure v_cross_op, v_cross
    end interface
    interface operator (.x.)
        procedure v_cross_op, v_cross
    end interface
    interface diag
        procedure v_diag, D_diag
    end interface
    interface inv
        procedure mat3_inv
    end interface
    interface solve
        procedure v_solve, v_solve2
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
    
    pure function a_identity(n) result(A)
    integer, intent(in) :: n
    real(wp) :: A(n,n)
    integer :: i
        A = 0d0
        forall(i=1:n)
            A(i,i) = 1d0
        end forall
    end function
    
    pure function v_elemental(n,i) result(v)
    integer, intent(in) :: n, i
    real(wp) :: v(n)
        v = 0d0
        v(i) = 1d0        
    end function
    
    pure function a_elemental(n,index) result(A)
    integer, intent(in) :: n, index(:)
    real(wp), allocatable :: A(:,:)
    integer :: i, nk
        nk = size(index)
        allocate(A(n, nk), source=0d0)
        do i=1, nk
            A(index(i),i) = 1d0
        end do
    end function
    
    pure function a_null(n, index) result(A)
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
    
    pure function v_cross(v_1,v_2) result(v_3)
    real(wp), intent(in) :: v_1(3), v_2(3)
    real(wp) :: v_3(3)
        v_3 = [ v_1(2)*v_2(3)-v_1(3)*v_2(2), &
                v_1(3)*v_2(1)-v_1(1)*v_2(3), &
                v_1(1)*v_2(2)-v_1(2)*v_2(1) ]
    end function
    
    pure function v_cross_op(v) result(vx)
    real(wp), intent(in) :: v(3)
    real(wp) :: vx(3,3)
        vx(1,:) = [0d0, -v(3), v(2)]
        vx(2,:) = [v(3), 0d0, -v(1)]
        vx(3,:) = [-v(2), v(1), 0d0]
    end function
    
    pure function v_mmoi(v) result(mm)
    real(wp), intent(in) :: v(3)
    real(wp) :: mm(3,3)
        mm(1,:) = [v(2)**2+v(3)**2, -v(1)*v(2), -v(1)*v(3)]
        mm(2,:) = [-v(1)*v(2), v(1)**2+v(3)**2, -v(2)*v(3)]
        mm(3,:) = [-v(1)*v(3), -v(2)*v(3), v(1)**2+v(2)**2]
    end function
    
    pure function v_mag(v) result(m)
    !$omp declare simd(v_mag)
    real(wp), intent(in) :: v(:)
    real(wp) :: m
        m = sqrt( sum( v**2 ))
    end function
    
    pure function v_dot(v, w) result(m)
    !$omp declare simd(v_dot)
    real(wp), intent(in) :: v(:), w(:)
    real(wp) :: m    
        m = sum(v*w)
	end function
    
    pure function v_diag(v) result(D)
    !$omp declare simd(v_diag)
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
    !$omp declare simd(D_diag)
    real(wp), intent(in) :: D(:,:)
    real(wp), allocatable :: v(:)
    integer :: n,i
        n = min(size(D,1),size(D,2))
        allocate(v(n))        
        do i=1,n
            v(i) = D(i,i)
        end do
    end function
    
    pure function v_matmul(A, v) result(w)
    ! evaluate w = A*v
    real(wp), intent(in) :: A(:,:), v(:)
	real(wp) :: w(size(A,1))
    integer :: i        
        do i=1,size(A,1)
            w(i) = sum(A(i,:)*v(:))
        end do
    end function
    
    pure function v_matmul2(v, A) result(w)
    ! evaluate w = v*A
    real(wp), intent(in) :: A(:,:), v(:)
	real(wp) :: w(size(A,2))
    integer :: i
        do i=1,size(A,2)
            w(i) = sum(v(:)*A(:,i)) 
        end do
    end function
    
    pure function mat3_inv(A) result(B)
    !$omp declare simd(mat3_inv)
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

    pure function v_solve(A, b) result(x)
    !$omp declare simd(v_solve)
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
    
    pure function v_solve2(A, B) result(x)
    !$omp declare simd(v_solve2)
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
        
        !do concurrent (i=1:3)
        !    w = B(:,i)/d
        !    x(:,i) = [ &
        !        w(1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+w(2)*(A(1,3)*A(3,2)-A(1,2)*A(3,3))+w(3)*(A(1,2)*A(2,3)-A(1,3)*A(2,2)), &
        !        w(1)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+w(2)*(A(1,1)*A(3,3)-A(1,3)*A(3,1))+w(3)*(A(1,3)*A(2,1)-A(1,1)*A(2,3)), &
        !        w(1)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))+w(2)*(A(1,2)*A(3,1)-A(1,1)*A(3,2))+w(3)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))]
        !end do
    end function
    
    function v_rotate_x(v,q) result(u)
    real(wp), intent(in) :: v(3), q
    real(wp) :: u(3), cq, sq
        cq = cos(q)
        sq = sin(q)
        u = [ v(1), &
            cq*v(2)-sq*v(3), &
            sq*v(2)+cq*v(3)]
    end function
    function v_rotate_y(v,q) result(u)
    real(wp), intent(in) :: v(3), q
    real(wp) :: u(3), cq, sq
        cq = cos(q)
        sq = sin(q)
        u = [ cq*v(1)+sq*v(3), &
            v(2), &
            -sq*v(1)+cq*v(3)]
    end function
    function v_rotate_z(v,q) result(u)
    real(wp), intent(in) :: v(3), q
    real(wp) :: u(3), cq, sq
        cq = cos(q)
        sq = sin(q)
        u = [ cq*v(1)-sq*v(1), &
            sq*v(1)+cq*v(2), &
            v(3)]
    end function
    
    function a_mmul(A,x,b) result(y)
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
    
    end module