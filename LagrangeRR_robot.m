function tau = LagrangeRR_robot(l, m, c, I, g)
% =========================================== %
% This function derives Lagrange Equations    %
% of the RR-robot to get torques of each link %
% =========================================== %

syms theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2

% ======================== %
% Acquire nessecary values %
% ======================== %

s_1 = { % Link 0->1 position vector
    @(theta1) l(1) * cos(theta1)                    % s_1{1}(t)
    @(theta1) l(1) * sin(theta1)                    % s_1{2}(t)
    0
};
s_2 = { % Link 1->2 position vector
    @(theta1, theta2) l(2) * cos(theta1 + theta2)   % s_2{1}(t)
    @(theta1, theta2) l(2) * sin(theta1 + theta2)   % s_2{2}(t)
    0
};

% Angular velocity as a 3D-vector
omega_1 = [ 0; 0; dtheta1 ];
omega_2 = [ 0; 0; dtheta1 + dtheta2 ];

% Recursive Forward Velocity Analysis
v_1 = cross(omega_1, s_1);
v_c1 = cross(omega_1, s_1)/2;
v_c2 = v_1 + cross(omega_2, s_2)/2;

% =============================== %
% Calculate translational and     %
% potential energies of each link %
% to get Lagrange expression L.   %
% =============================== %

T1 = simplify( 1/2 * m(1) * transpose(v_c1) * v_c1 + 1/2 * I(1) * omega_1(3)^2 );
T2 = simplify( 1/2 * m(2) * transpose(v_c2) * v_c2 + 1/2 * I(2) * omega_2(3)^2 );

V1 = m(1) * g * c(1) * sin(theta1);
V2 = m(2) * g * (l(1)*sin(theta1) + c(2)*sin(theta1 + theta2));

L = T1 + T2 - V1 - V2;

% ===================================== %
% Acquire partial derivatives nessecary %
% and calculate joint torques.          %
% ===================================== %
pd_T1 = diff(L, dtheta1);
pd_V1 = diff(L, theta1);
pd_T2 = diff(L, dtheta2);
pd_V2 = diff(L, theta2);

% Get tau1
tau1 = diff(pd_T1, theta1) * dtheta1 + diff(pd_T1, theta2) * dtheta2 + ...
    diff(pd_T1, dtheta1) * ddtheta1 + diff(pd_T1, dtheta2) * ddtheta2 - ...
    pd_V1;

% Get tau2
tau2 = diff(pd_T2, theta1) * dtheta1 + diff(pd_T2, theta2) * dtheta2 + ...
    diff(pd_T2, dtheta1) * ddtheta1 + diff(pd_T2, dtheta2) * ddtheta2 - ...
    pd_V2;

tau = [tau1, tau2];