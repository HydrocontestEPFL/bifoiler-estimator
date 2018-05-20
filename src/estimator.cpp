#include <Eigen/Dense>

namespace bifoiler {

class MEKF<typename Scalar> {
public:
    using Scalar = typename Scalar;
    using State = Eigen::Matrix<Scalar, 18, 1>; // [v, w, r, a, bg, ba];
    using SysState = Eigen::Matrix<Scalar, 13, 1>; // [v, w, r, q_BI];
    using Control = typename Dynamics::Control;
    using Measurement = typename Observation::Measurement;
    using StateCov = Eigen::Matrix<Scalar, State::RowsAtCompileTime, State::RowsAtCompileTime>;
    using Quaternion = Eigen::Matrix<Scalar, 4, 1>; // [w, x, y, z];
    enum {
        nx = State::RowsAtCompileTime,
        nz = Measurement::RowsAtCompileTime,
    };

private:
    Dynamics f;
    Observation h;
    State x; // state vector
    Quaternion qref;
    StateCov P; // state covariance
    Eigen::Matrix<Scalar, nz, nz> S; // innovation covariance
    Eigen::Matrix<Scalar, nx, nz> K; // Kalman gain
    Eigen::Matrix<Scalar, nx, nx> I; // identity

    static SysState _to_system_state(const State &x)
    {
        SysState x_sys;
        x_sys.block<9,1>(0,0) = x.block<9,1>(0,0);
        Matrix<Scalar, 3, 1> a;
        Matrix<Scalar, 4, 1> q, qref;
        a = x.block<3,1>(9,0);
        // qref = x.block<4,1>(9,0)
        q = _quat_error_mult(a, qref);
        q = quatnormalize(q);
        x_sys.block<4,1>(9,0) = q;
    }

public:
    void update(const Control &u, const Measurement &z);

};

/*

function [q, eta] = quat_error_mult(a,qr)


% Factor
f             = 1;

% Scalar Part
% eta         = 1;
eta           = 1 - norm(a)^2/8;

%Vector Part
% nu          = a;
nu            = a/2;

% Error Parametrisation
dq          = f*[eta; nu];

%Erroneous Quaternion
q           = quatmul(dq,qr);

q           = q/norm(q);
*/

/*
%% (M)EKF dynamics/ Composition of total Jacobian

% Model Jacobian to get linearized model
A_sys       = dyn_jac(xs,u);

% df/da = df/dq*dq/da; Jq := dq/da
Jq          = q.jacobian(a);

% Matrix valued Chainrule
Aa          = A_sys(1:9,10:13)*Jq;

% Attitude Dynamics (from Nokland)
f_att_nl    = [1/2*(eye(3)*eta + skew_mat(a))*(W-bg); zeros(3,1)];

% Attitude Dynamics (from Markley)
f_att_mk    = [-bg-skew_mat(W)*a; zeros(3,1)];


Jatt        = f_att_mk.jacobian(xe);
%Jatt(1:3,4:6) = zeros(3);

%Jatt        = [zeros(6,3)

A           = [A_sys(1:9,1:9) Aa zeros(9,6); Jatt; zeros(3,length(xe))];
A_func      = Function('A_func',{xe,u,qr},{A});

%% Discretization

% Euler
F           = eye(size(A)) + t_samp*A; % + 1/2*(t_samp^2*A*A); 

% Jordan Normal Form; jordan() and expm() don't work with casadi
% [V,J]       = jordan(F);
% F           = expm(F*t_samp);

F_func      = Function('F_func',{xe,u,qr},{F});
*/

void MEKF::update(const Control &u, const Measurement &z)
{
    /*
    qref            = sys_state_est(10:13);

    %% Nonlinear Prediction

    x0              = sys_state_est;

    % Numerical integration
    % out = CVODES_INT('x0',x0, 'p',u0);
    out             = CVODES_INT('x0',x0,'p',u);

    % x_p = [v_pred w_pred r_pred q_pred];
    x_p             = full(out.xf);

    % Evaluate symbolic Statepropag.  at current time step
    F               = F_func(state_est,u, qref);

    % xe_p = [v_pred w_pred r_pred a_pred bg_pred ba_pred];
    % Propagation of the latter three is trivial, derivatives are zero.
    % Hardcode is not nice but 'end' is spooky...
    xe_p            = [x_p(1:9); F(10:18,:)*state_est];
    %xe_p            = F*state_est;


    %% Covariance Prediction

    P_last          = state_var;

    % For now assume discrete dynamics. If assume model to be continuous, have
    % to use Jacobians and integrate ricatti equation here.
    P_p             = F*P_last*F' + Q;
    */

    // prediction using numerical integration
    SysState x_sys, x_sys_pred
    x_sys = _to_system_state(x);
    x_sys_p = Dynamics.integrate(x_sys, u);


    // TODO: how to get F_func? why does it need qref?
    F = F_func(x, u, qref);

    State xe_p;
    xe_p.block<9,1>(0,0) = x_sys_p.block<9,1>(0,0);
    xe_p.block<9,1>(9,0) = F.block<9,18>(9,0)*x;


    /*

    %% Measurement Update

    H               = H_func(state_est,u, qref,r_ant);
    % full() because Matlab functions like size() not applicable to casadi-sparse
    K               = P_p*H'*inv(full(H*P_p*H') + R);
    I               = eye(length(state_est));
    state_var       = (I - K*H)*P_last*(I-K*H)' + K*R*K';

    y               = h_func(xe_p,u,qref,r_ant);
    %y               = H*xe_p;

    state_est       = xe_p + K*(z_m - y);

    %% Transforming Estimator State xe to System State xs

    a               = state_est(10:12);

    q               = quat_error_mult(a,qref);

    % Enforcing Unit norm constraint
    q               = quatnormalize(full(q'))';

    xem             = state_est;
    Pm              = state_var;
    sys_state_est   = [xem(1:9); q];
    xs              = sys_state_est;

    % Reset MEKF
    state_est(10:12)= zeros(3,1);
    */
}

} // namespace bifoiler
