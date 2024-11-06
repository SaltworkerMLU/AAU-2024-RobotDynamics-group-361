% ========================================= %
% ROBOT DYNAMICS ASSIGNMENT 1
% Provided a 2-dof RR planar manipulator with lengths:
%   - Link lengths: l_1 = 1.0m, l_2 0 0.5m, 
%   - Link masses: m_1 = 19.5kg, m_2 = 9.75kg
% Simulate robot dynamics using:
%
% (b) Constant EE velocity: v_2 = [0, 0.5]
% The trajectory is given in Caresian space, i.i. we know the 
% end-effector motion profile:
% ============================================================= %
%
% The EE moves at a constant velocity, v=[0, 0.5], and the 
% initial angles are theta_1 = 10*pi/180, theta_2 = 90*pi/180.
%
% Run simulation for this case
% ============================================================= %
%   - EE velocity: v_2 = [0, 0.5]
%   - EE acceleration: a_2 = [0, 0] % Because constant velocity
%   - Initial angles: 
%       * theta_0{1} = 10*pi/180
%       * theta_0{2} = 90*pi/180
% ============================================================= %

clear; close; clc
syms theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2

% ============================== %
% Initialize provided parameters %
% ============================== %
l = [ 1.0, 0.5 ];                                   % Link length
m = [ 19.5, 9.75 ];                                 % Link mass
c = [ l(1)/2, l(2)/2 ];                             % mass center
I = [ 1/12 * m(1) * l(1)^2, 1/12 * m(2) * l(2)^2 ]; % moment of inertia
g = 9.801;                                          % Gravity constant

% ================================================= %
% Initialize a bunch of arrays used for plotting    %
% theta, dtheta, ddtheta, pos, vel, acc, tau        %
% ================================================= %
fps = 60;                                           % Your PC most likely has a 60Hz screen
animation_time = 1;                                 % Duration of plot animation [s]
num_frames = animation_time * fps;

theta = {zeros(num_frames), zeros(num_frames)};
dtheta = {zeros(num_frames), zeros(num_frames)};
ddtheta = {zeros(num_frames), zeros(num_frames)};
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

% Yep, this is the motion planner -> moveL()
v_2 = [0, 0.5]; % m/s
a_2 = [0, 0]; % m/s^2

theta0 = [10*pi/180, 90*pi/180];

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

P1_0{1} = s_1{1}(theta0(1));
P1_0{2} = s_1{2}(theta0(1));
P2_0{1} = P1_0{1} + s_2{1}(theta0(1), theta0(2));
P2_0{2} = P1_0{2} + s_2{2}(theta0(1), theta0(2));

% ========================================== %
% Calculate desired values frame by frame    %
% ========================================== %

i = 0;
for t = linspace(0, animation_time, num_frames)
    i = i + 1; time(i) = t;

    % Link 2 position, velocity, and acceleration
    P2{1}(i) = P2_0{1} + v_2(1) * t;
    P2{2}(i) = P2_0{2} + v_2(2) * t;
    dP2{1}(i) = v_2(1);
    dP2{2}(i) = v_2(2);
    ddP2{1}(i) = a_2(1);
    ddP2{2}(i) = a_2(2);

    % Inverse position analysis
    x = P2_0{1} + v_2(1) * t;
    y = P2_0{2} + v_2(2) * t;
    theta{2}(i) = acos( (x^2 + y^2) - (l(1)^2 + l(2)^2) ) / ( 2 * l(1) * l(2) ); % try this: % -acos(...) 
    t1a = l(2) * sin(theta{2}(i));
    t1b = l(1) + l(2)*cos(theta{2}(i));
    theta{1}(i) = atan( (-t1a*x + t1b*y) / (t1a*y + t1b*x) );
    
    % Link 1 position, velocity, and acceleration
    P1{1}(i) = s_1{1}(theta{1}(i)); % x-coordinate of P1 -> X1
    P1{2}(i) = s_1{2}(theta{1}(i)); % y-coordinate of P1 -> Y1
end

% Get angular velocity and acceleration of P1 and P2
dtheta{1} = gradient(theta{1}(1:i)) ./ gradient(time(1:i));
dtheta{2} = gradient(theta{2}(1:i)) ./ gradient(time(1:i));
ddtheta{1} = gradient(dtheta{1}(1:i)) ./ gradient(time(1:i));
ddtheta{2} = gradient(dtheta{2}(1:i)) ./ gradient(time(1:i));

% Get cartesian velocity and acceleration of P1 and P2
dP1{1} = gradient(P1{1}(1:i)) ./ gradient(time(1:i));
dP1{2} = gradient(P1{2}(1:i)) ./ gradient(time(1:i));
ddP1{1} = gradient(dP1{1}(1:i)) ./ gradient(time(1:i));
ddP1{2} = gradient(dP1{2}(1:i)) ./ gradient(time(1:i));

% Get torques using Lagrange Equations
eq_tau = LagrangeRR_robot(l, m, c, I, g);

i = 0;
for t = linspace(0, animation_time, num_frames)
    i = i + 1; time(i) = t;
    tau{1}(i) = subs(eq_tau(1), {theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2}, ...
    {theta{1}(i) theta{2}(i) dtheta{1}(i) dtheta{2}(i) ddtheta{1}(i) ddtheta{2}(i)});

    tau{2}(i) = subs(eq_tau(2), {theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2}, ...
    {theta{1}(i) theta{2}(i) dtheta{1}(i) dtheta{2}(i) ddtheta{1}(i) ddtheta{2}(i)});  
end

%% plot the motion profiles
figure(1)
subplot(3, 1, 1)
hold on
plot(time, theta{1});
plot(time, theta{2});
hold off
legend('joint 1', 'joint 2')
grid on; 
xlabel('time [sec]'); ylabel('theta [rad]'); 

subplot(3, 1, 2)
hold on
plot(time, dtheta{1})
plot(time, dtheta{2})
hold off
grid on; 
xlabel('time [sec]'); ylabel('dtheta [rad/s]'); 
subplot(3, 1, 3)
hold on
plot(time, ddtheta{1})
plot(time, ddtheta{2})
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
legend(' joint 1', 'joint 2')
grid on; 
xlabel('time [sec]'); ylabel('torques [Nm/rad]'); 

%%  plot links
figure(3);
hold on 
axis equal; grid on; axis([-2 2 -2 2]); 
title('Robot start position')
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
    quiver(P1{1}(j), P1{2}(j), dP1{1}(j), dP1{2}(j), 'Color', 'm', 'Linewidth', 1); hold on
    quiver(P1{1}(j), P1{2}(j), ddP1{1}(j), ddP1{2}(j), 'Color', 'r', 'Linewidth', 1); hold on

    % Plot link 2
    line([P1{1}(j) P2{1}(j)], [P1{2}(j) P2{2}(j)], 'Color', '#888888', 'Linewidth', 2); hold on
    quiver(P2{1}(j), P2{2}(j), dP2{1}(j), dP2{2}(j), 'Color', 'm', 'Linewidth', 1); hold on
    quiver(P2{1}(j), P2{2}(j), ddP2{1}(j), ddP2{2}(j), 'Color', 'r', 'Linewidth', 1); hold on

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