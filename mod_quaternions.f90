    module mod_quaternions
    use mod_vectors
    implicit none (type, external)

    reaL(wp), parameter :: q_zero(4) = [0d0,0d0,0d0,0d0]
    reaL(wp), parameter :: q_eye(4) = [0d0,0d0,0d0,1d0]

    interface rot
        procedure q_rot_matrix, q_axis_angle, q_rot_vector
    end interface
    
    interface inv
        procedure q_inv
    end interface
    
    interface operator (.o.)
        procedure q_mul, q_mul_op
    end interface

    contains

    !-- QUATERNIONS
    pure function q_axis_angle(k,u) result(q)
    real(wp), intent(in) :: k(3), u
    real(wp) :: q(4)
        q = [ k*sin(u/2), cos(u/2) ]
    end function
    !
    pure function q_angle(q) result(u)
    real(wp), intent(in) :: q(4)
    real(wp) :: u, c, s
        c = q(4)
        s = v_mag(q(1:3))
        u = 2d0*atan(s/c)
    end function

    pure function q_axis(q) result(k)
    real(wp), intent(in) :: q(4)
    real(wp) :: k(3)
        k = q(1:3)
        k = k/v_mag(k)
    end function
    
    pure function q_rot_vector(q,t) result(u)
    real(wp), intent(in) :: t(3), q(4)
    real(wp) :: u(3)
    real(wp) :: v(3), s, vxt(3), vxvxt(3)
        v = q(1:3)
        s = q(4)
        vxt = v .x. t
        vxvxt = v .x. vxt
        u = t + 2*s*vxt + 2*vxvxt
    end function

    pure function q_rot_matrix(q) result(R)
    real(wp), intent(in) :: q(4)
    real(wp) :: R(3,3)
    real(wp) :: v(3), s, vx(3,3), vxvx(3,3)
        v = q(1:3)
        s = q(4)
        vx = v_cross_op(v)
        vxvx = -v_mmoi(v)
        R = eye3 + 2*s*vx + 2*vxvx
    end function

    pure function q_inv(q) result(qi)
    real(wp), intent(in) :: q(4)
    real(wp) :: qi(4), qm2, qc(4)
        qm2 = dot_product(q,q)
        qc = [-q(1:3),q(4)]
        qi = qc/qm2
    end function
    !
    pure function q_mul(q_1, q_2) result(q_3)
    real(wp), intent(in) :: q_1(4), q_2(4)
    real(wp) :: q_3(4)
    !real(wp) :: v_1(3), v_2(3), s_1, s_2
        !v_1 = q_1(1:3)
        !s_1 = q_1(4)
        !v_2 = q_2(1:3)
        !s_2 = q_2(4)
        !
        !q_3(1:3) = s_1*v_2 + s_2*v_1 + v_cross(v_1, v_2)
        !q_3(4) = s_1*s_2 - dot_product(v_1, v_2)
    
        q_3 = matmul(q_mul_op(q_1), q_2)
    end function
    
    pure function q_mul_op(q) result(op)
    real(wp), intent(in) :: q(4)
    real(wp) :: op(4,4)
        op(:,1) = [q(4),q(3),-q(2),-q(1)]
        op(:,2) = [-q(3),q(4),q(1),-q(2)]
        op(:,3) = [q(2),-q(1),q(4),-q(3)]
        op(:,4) = q
    end function
    
    pure function q_der(q,omg) result(qp)
    real(wp), intent(in) :: q(4), omg(3)
    real(wp) :: qp(4), w(4)
        w = [ omg, 0d0]
        qp = 0.5d0*q_mul(w, q)
    end function

    pure function q_omg(q,qp) result(omg)
    real(wp), intent(in) :: q(4), qp(4)
    real(wp) :: omg(3), v(3), vp(3)
        v=q(1:3)
        vp=qp(1:3)
        omg=2*(q(4)*vp-qp(4)*v+(v .x. vp))
    end function

    pure function q_step_omg(q,omg,h) result(q_next)
    real(wp), intent(in) :: q(4), omg(3), h
    real(wp) :: k(3), qw(4), qp(4), u, q_next(4)
    real(wp), parameter :: rot_tol = 1d-4*deg
        u = h*mag(omg)
        if( u>rot_tol) then
            k = h*omg/u
            qw = q_axis_angle(k, u)
            q_next = q_mul(qw, q)
        else
            qp = q_der(q,omg)
            q_next = q + h*qp
        end if
    end function

    pure function q_step_qp(q,qp,h) result(q_next)
    real(wp), intent(in) :: q(4), qp(4), h
    real(wp) :: omg(3), q_next(4), u, qw(4), k(3), eps
    real(wp), parameter :: rot_tol = 1d-4*deg
        eps = dot_product(q, q)-1d0
        if( eps>rot_tol ) then
            omg = q_omg(q, qp)
            u = mag(omg)
            k = omg/u
            qw = q_axis_angle(k, u)
            q_next = q_mul(qw, q)
        else
            q_next = q + h*qp
        end if
    end function

    end module