% ========================================= %
% ROBOT DYNAMICS ASSIGNMENT 1 - Inverse dynamics
% Provided a 2-dof RR planar manipulator with lengths:
%   - Link lengths: l_1 = 1.0m, l_2 0 0.5m, 
%   - Link masses: m_1 = 19.5kg, m_2 = 9.75kg
% Simulate robot dynamics using:
%
%   (a) known joint motion profile
%           - theta_1(t) = 3 * sin(pi * t)
%           - theta_2(t) = 0.5 * sin(3 * pi * t + 45*pi/180)
% Plot:
%   (1) Joint angles, velocities, and acceleration
%       theta(t), dtheta(t), ddtheta(t)
%   (2) EE position, velocity, and acceleration
%   (3) joint torques
% ========================================= %

clear; close; clc

syms t % t is only used for differentiating theta(t)
syms theta1 dtheta1 ddtheta1 theta2 dtheta2 ddtheta2

% ============================== %
% Initialize provided parameters %
% ============================== %
l = [ 1.0, 0.5 ];                                   % Link length
m = [ 19.5, 9.75 ];                                 % Link mass
c = [ l(1)/2, l(2)/2 ];                             % mass center
I = [ 1/3 * m(1) * l(1)^2, 1/3 * m(2) * l(2)^2 ];   % moment of inertia
g = 9.801;                                          % Gravity constant

% ================================================= %
% Initialize a bunch of arrays used for plotting    %
% theta, dtheta, ddtheta, pos, vel, acc, tau        %
% ================================================= %
fps = 60;                                           % Your PC most likely has a 60Hz screen
animation_time = 2;                                 % Duration of plot animation [s]
num_frames = animation_time * fps; 

joint = {zeros(num_frames), zeros(num_frames)};
djoint = {zeros(num_frames), zeros(num_frames)};
ddjoint = {zeros(num_frames), zeros(num_frames)};
P1 = {zeros(num_frames), zeros(num_frames)};
P2 = {zeros(num_frames), zeros(num_frames)};
dP1 = {zeros(num_frames), zeros(num_frames)};
dP2 = {zeros(num_frames), zeros(num_frames)};
ddP1 = {zeros(num_frames), zeros(num_frames)};
ddP2 = {zeros(num_frames), zeros(num_frames)};
tau = {zeros(num_frames), zeros(num_frames)};
time = zeros(num_frames);
frames = repmat(struct(getframe), num_frames, 2); % the "getframe"-command provides sufficient struct

% ================================================ %
% Set motion planner and initial state of RR-robot %
% ================================================ %

theta = { % Get theta(t)
    @(t) 3*sin(pi*t);                               % theta{1}(t)
    @(t) 0.5*sin(3*pi*t + 45*pi/180)                % theta{2}(t)
};
dtheta = { % Get angular velocity by differentiation
    str2func([ '@(t)' char(diff(theta{1}, t)) ])    % dtheta{1}(t)
    str2func([ '@(t)' char(diff(theta{2}, t)) ])    % dtheta{2}(t)
};
ddtheta = { % Get angular acceleration by differentiation
    str2func([ '@(t)' char(diff(dtheta{1}, t)) ])   % ddtheta{1}(t)
    str2func([ '@(t)' char(diff(dtheta{2}, t)) ])   % ddtheta{2}(t)
};

% ========================================================= %
% Get link position, velocity, and acceleration             %
% using recursive method.                                   %
% ========================================================= %
% The alternative is to get Jacobian J and its time         %
% derivative dJ, yet this requires more processes compared  %
% to using recursive method.                                %
% ========================================================= %

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

[ omega_1, domega_1 ] = deal([ 0; 0; dtheta1 ], [ 0; 0; ddtheta1 ]);
[ omega_2, domega_2 ] = deal(omega_1 + [ 0; 0; dtheta2 ], domega_1 + [ 0; 0; ddtheta2 ]);

% Get link 1 velocity and acceleration function
v_1 = cross(omega_1, s_1);
a_1 = cross(domega_1, s_1) + cross(omega_1, cross(omega_1, s_1));

ds_1 = { % Link 0->1 velocity vector
    str2func([ '@(theta1, dtheta1)' char(v_1(1)) ]) % ds_1{1}(theta1)
    str2func([ '@(theta1, dtheta1)' char(v_1(2)) ]) % ds_1{2}(theta1)
};
dds_1 = { % Link 0->1 acceleration vector
    str2func([ '@(theta1, dtheta1, ddtheta1)' char(a_1(1)) ]) % dds_1{1}(theta1)
    str2func([ '@(theta1, dtheta1, ddtheta1)' char(a_1(2)) ]) % dds_1{2}(theta1)
};

% Get link 2 velocity and acceleration
v_2 = v_1 + cross(omega_2, s_2);
a_2 = a_1 + cross(domega_2, s_2) + cross(omega_2, cross(omega_2, s_2));

ds_2 = {
    str2func([ '@(theta1, dtheta1, theta2, dtheta2)' char(v_2(1)) ]) % ds_2{1}(theta2)
    str2func([ '@(theta1, dtheta1, theta2, dtheta2)' char(v_2(2)) ]) % ds_2{2}(theta2)
};
dds_2 = {
    str2func([ '@(theta1, dtheta1, ddtheta1, theta2, dtheta2, ddtheta2)' char(a_2(1)) ]) % dds_2{1}(theta2)
    str2func([ '@(theta1, dtheta1, ddtheta1, theta2, dtheta2, ddtheta2)' char(a_2(2)) ]) % dds_2{2}(theta2)
};

% ========================================== %
% Calculate desired values frame by frame    %
% ========================================== %

i = 0;
for t = linspace(0, animation_time, num_frames)
    i = i + 1; time(i) = t;

    % theta(t), dtheta(t), & ddtheta(t)
    [ joint{1}(i), joint{2}(i) ] = deal(theta{1}(t), theta{2}(t));
    [ djoint{1}(i), djoint{2}(i) ] = deal(dtheta{1}(t), dtheta{2}(t));
    [ ddjoint{1}(i), ddjoint{2}(i) ] = deal(ddtheta{1}(t), ddtheta{2}(t));
    
    % Link 1 position, velocity, and acceleration
    %[P1{1}(i), P1{2}(i)] = deal(s_1{1}(theta{1}(t)), s_1{2}(theta{1}(t)));
    P1{1}(i) = s_1{1}(theta{1}(t)); % x-coordinate of P1 -> X1
    P1{2}(i) = s_1{2}(theta{1}(t)); % y-coordinate of P1 -> Y1
    dP1{1}(i) = ds_1{1}(theta{1}(t), dtheta{1}(t)); % dx/dt
    dP1{2}(i) = ds_1{2}(theta{1}(t), dtheta{1}(t)); % dy/dt
    ddP1{1}(i) = dds_1{1}(theta{1}(t), dtheta{1}(t), ddtheta{1}(t)); % d^2x/dt^2
    ddP1{2}(i) = dds_1{2}(theta{1}(t), dtheta{1}(t), ddtheta{1}(t)); % d^2y/dt^2
    
    % Link 2 position, velocity, and acceleration
    %[P2{1}(i), P2{2}(i)] = deal(P1{1}(i) + s_2{1}(theta{1}(t), theta{2}(t)), P1{2}(i) + s_2{2}(theta{1}(t), theta{2}(t)));
    P2{1}(i) = P1{1}(i) + s_2{1}(theta{1}(t), theta{2}(t)); % x-coordinate of P2 -> X2
    P2{2}(i) = P1{2}(i) + s_2{2}(theta{1}(t), theta{2}(t)); % y-coordinate of P2 -> Y2
    dP2{1}(i) = ds_2{1}(theta{2}(t), dtheta{2}(t), theta{2}(t), dtheta{2}(t)); % dx/dt
    dP2{2}(i) = ds_2{2}(theta{2}(t), dtheta{2}(t), theta{2}(t), dtheta{2}(t)); % dy/dt
    ddP2{1}(i) = dds_2{1}(theta{2}(t), dtheta{2}(t), ddtheta{2}(t), theta{2}(t), dtheta{2}(t), ddtheta{2}(t)); % d^2x/dt^2
    ddP2{2}(i) = dds_2{2}(theta{2}(t), dtheta{2}(t), ddtheta{2}(t), theta{2}(t), dtheta{2}(t), ddtheta{2}(t)); % d^2y/dt^2
    
    % Coefficients of dynamic equation
    H12 = m(2)*l(1)*c(2)*cos(theta{2}(t)) + I(2); 
    H11 = I(1) + m(2)*(l(1)^2 + l(1)*c(2)*cos(theta{2}(t))) + H12;
    H22 = I(2); 

    h = m(2)*c(2)*sin(theta{2}(t));
    G1 = m(1)*g*c(1)*cos(theta{1}(t)) + m(2)*g*( c(2)*cos(theta{1}(t) + theta{2}(t)) + l(1) * cos(theta{1}(t) ) ); 
    G2 = m(2)*g*c(2)*cos(theta{1}(t) + theta{2}(t));
    
    % Get actuator torques
    tau{1}(i) = [H11, H12] * [ddtheta{1}(t); ddtheta{2}(t)] - h*(dtheta{2}(t))^2 ...
        -2*h*dtheta{1}(t)*dtheta{2}(t) + G1;
    tau{2}(i) = [H12, H22] * [ddtheta{1}(t); ddtheta{2}(t)] + h*(dtheta{1}(t))^2 + G2;  
end

%% plot the motion profiles
figure(1)
subplot(3, 1, 1)
hold on
plot(time, joint{1});
plot(time, joint{2});
hold off
legend('joint 1', 'joint 2')
grid on; 
xlabel('time [sec]'); ylabel('theta [rad]'); 
subplot(3, 1, 2)
hold on
plot(time, djoint{1})
plot(time, djoint{2})
hold off
grid on; 
xlabel('time [sec]'); ylabel('dtheta [rad/s]'); 
subplot(3, 1, 3)
hold on
plot(time, ddjoint{1})
plot(time, ddjoint{2})
hold off
grid on; 
xlabel('time [sec]'); ylabel('ddtheta [rad/s^2]');

%% plot the link motion profiles
figure(4)
subplot(3, 1, 1)
hold on
plot(time, P1{1});
plot(time, P1{2});
plot(time, P2{1});
plot(time, P2{2});
hold off
legend('Link 1 x-coordinate', 'Link 1 y-coordinate', ...
    'Link 2 x-coordinate', 'Link 2 y-coordinate')
grid on; 
xlabel('time [sec]'); ylabel('Pos [m]'); 
subplot(3, 1, 2)
hold on
plot(time, dP1{1});
plot(time, dP1{2});
plot(time, dP2{1});
plot(time, dP2{2});
hold off
grid on; 
xlabel('time [sec]'); ylabel('Vel [m/s]'); 
subplot(3, 1, 3)
hold on
plot(time, ddP1{1});
plot(time, ddP1{2});
plot(time, ddP2{1});
plot(time, ddP2{2});
hold off
grid on; 
xlabel('time [sec]'); ylabel('Acc [m/s^2]');

%% plot the motor torques
figure(2)
hold on
plot(time, tau{1})
plot(time, tau{2})
hold off
legend('joint 1', 'joint 2')
grid on; 
xlabel('time [sec]'); ylabel('torques [Nm/rad]'); 


%%  plot links
figure(3);
hold on 
axis equal; grid on; axis([-2 2 -2 2]); 
xlabel('x [m]'); ylabel('y [m]'); 

for j = 1 : num_frames  
    set(gca,'NextPlot','replaceChildren'); % Do this every iteration of loop to keep all vectors & dots
    % Origo
    plot(0, 0, 'k*', 'MarkerSize', 6, 'Linewidth', 2); hold on
    

    % Plot end-points of link 1 & 2 (aka. EE)
    plot(P1{1}(j), P1{2}(j), 'b*', 'MarkerSize', 3, 'Linewidth', 2); hold on
    plot(P2{1}(j), P2{2}(j), 'r*', 'MarkerSize', 3, 'Linewidth', 2); hold on
    
    % Plot link 1
    line([0 P1{1}(j)], [0 P1{2}(j)], 'Color', '#AAAAAA', 'Linewidth', 2); % Link
    %quiver(P1{1}(j), P1{2}(j), dP1{1}(j), dP1{2}(j), 'Color', 'b', 'Linewidth', 1); hold on
    %quiver(P1{1}(j), P1{2}(j), ddP1{1}(j), ddP1{2}(j), 'Color', 'm', 'Linewidth', 1); hold on

    % Plot link 2
    line([P1{1}(j) P2{1}(j)], [P1{2}(j) P2{2}(j)], 'Color', '#888888', 'Linewidth', 2); hold on
    %quiver(P2{1}(j), P2{2}(j), dP2{1}(j), dP2{2}(j), 'Color', 'm', 'Linewidth', 1); hold on
    %quiver(P2{1}(j), P2{2}(j), ddP2{1}(j), ddP2{2}(j), 'Color', 'r', 'Linewidth', 1); hold on
    hold off

    % Store frame F
    frames(j) = getframe;
end
hold off
%% Redo link animation
f3 = figure(3);
hold on 
axis equal; grid on; axis([-2 2 -2 2]); 
xlabel('x [m]'); ylabel('y [m]'); 
set(gca,'NextPlot','replaceChildren');
movie(frames, 1, fps);