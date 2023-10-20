clear
clc
% Time span for EKF Simulation
t_span_EKF = [0, 30];
%Inertia Tensor
%J = diag([1.016, 1.022, 1.795]);
%Stowed Inertia Tensor
J = diag([0.831, 0.838, 1.379]); 
q_e = [1, 0, 0, 0]';
w_e = [0.5 0.5 0.5]'; 
x_EKF_0 = [q_e; w_e];
[t_EKF, x_EKF] = ode45(@(t_EKF, x_EKF) xdot(x_EKF,J,0,0), t_span_EKF, x_EKF_0);

%% EKF Simulation 
syms q [4 1]
syms w [3 1]
syms T_rw [3 1] 
syms T_d [3 1]
w_dot = J \ (T_rw + T_d - cross(w, J*w));
OMEGA = [0, w(3), -w(2), w(1); -w(3), 0, w(1), w(2); w(2), -w(1), 0, w(3); -w(1), -w(2), -w(3), 0];
q_dot = 0.5 * OMEGA * q;
R_bi = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2)+q(4)*q(3)), 2*(q(1)*q(3)-q(4)*q(2)) ;
2*(q(1)*q(2)-q(4)*q(3)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)+q(4)*q(1) );
2*(q(1)*q(3)+q(4)*q(2)), 2*(q(2)*q(3)-q(4)*q(1)), -q(1)^2-q(2)^2+q(3)^2+q(4)^2];
q_EKF_e(1,1) = 0.5*rand + 0.5; 
q_EKF_e(2,1) = 0.1*rand;
q_EKF_e(3,1) = 0.1*rand;
q_EKF_e(4,1) = 0.1*rand;
w_EKF_e(1,1) = 0.5*rand+0.5;
w_EKF_e(2,1) = 0.5*rand+0.5;
w_EKF_e(3,1) = 0.5*rand+0.5;
T_rw_e = [0, 0, 0]';
s_i = approxECISunPosition(juliandate('20-June-2023','dd-mmm-yyyy'));
for i = 1:length(x_EKF)
    s_b(:,i) = [x_EKF(i,1)^2-x_EKF(i,2)^2-x_EKF(i,3)^2+x_EKF(i,4)^2, 2*(x_EKF(i,1)*x_EKF(i,2)+x_EKF(i,4)*x_EKF(i,3)), 2*(x_EKF(i,1)*x_EKF(i,3)-x_EKF(i,4)*x_EKF(i,2)) ;
        2*(x_EKF(i,1)*x_EKF(i,2)-x_EKF(i,4)*x_EKF(i,3)), -x_EKF(i,1)^2+x_EKF(i,2)^2-x_EKF(i,3)^2+x_EKF(i,4)^2, 2*(x_EKF(i,2)*x_EKF(i,3)+x_EKF(i,4)*x_EKF(i,1) );
        2*(x_EKF(i,1)*x_EKF(i,3)+x_EKF(i,4)*x_EKF(i,2)), 2*(x_EKF(i,2)*x_EKF(i,3)-x_EKF(i,4)*x_EKF(i,1)), -x_EKF(i,1)^2-x_EKF(i,2)^2+x_EKF(i,3)^2+x_EKF(i,4)^2]*s_i';
end
s_noise = s_b + 0.1*randn(3,length(t_EKF));
% Defining noise covariance matricies
R = diag([0.1^2 0.1^2 0.1^2]);
P = diag([(q_EKF_e(1) - q_e(1))^2, (q_EKF_e(2) - q_e(2))^2, (q_EKF_e(3) - q_e(3))^2, (q_EKF_e(4) - q_e(4))^2, (w_EKF_e(1) - w_e(1))^2, (w_EKF_e(2) - w_e (2))^2, (w_EKF_e(3) - w_e(3))^2]);
Q = diag([0.05^2, 0.05^2, 0.05^2, 0.05^2, 0.05^2, 0.05^2, 0.05^2]); % Creating state space model
f = [q_dot; w_dot];
% Simulating EKF
for i = 2:length(t_EKF)
t_EKF_prev = t_EKF(i-1);
t_EKF_current = t_EKF(i);
dt = (t_EKF_current - t_EKF_prev)/2;
[x_EKF_post(:,i), P_post] = EKF(q, w, T_rw, f, q_EKF_e, w_EKF_e, T_rw_e, dt, P, Q, R, R_bi, s_noise(:,i),s_i');
cov_err(:,i) = diag(P_post);
q_eq_EKF = x_EKF_post(1:4,i); 
w_eq_EKF = x_EKF_post(5:7,i); 
P = P_post;
end

figure
plot(t_EKF, x_EKF(:,1)-x_EKF_post(1,:)', 'b-', 'LineWidth', 1.5)
hold on
plot(t_EKF, x_EKF(:,2)-x_EKF_post(2,:)', 'r-', 'LineWidth', 1.5)
plot(t_EKF, x_EKF(:,3)-x_EKF_post(3,:)', 'k-', 'LineWidth', 1.5)
plot(t_EKF, x_EKF(:,4)-x_EKF_post(4,:)', 'g-', 'LineWidth', 1.5)
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Time [$\mathrm{s}$]','Interpreter','latex', 'FontSize',18) 
ylabel('Quaternions Error','Interpreter','latex', 'FontSize',18)
legend({'${q_1}$','${q_2}$','${q_3}$','${q_4}$'},'Interpreter','latex','FontSize',16,Location='northeast')
grid on
hold off
saveas(gcf,'q_kalman_error.png')

figure
plot(t_EKF, x_EKF(:,5)-x_EKF_post(5,:)', 'b-', 'LineWidth', 1.5)
hold on
plot(t_EKF, x_EKF(:,6)-x_EKF_post(6,:)', 'r-', 'LineWidth', 1.5)
plot(t_EKF, x_EKF(:,7)-x_EKF_post(7,:)', 'k-', 'LineWidth', 1.5)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Time [$\mathrm{s}$]','Interpreter','latex', 'FontSize',18) 
ylabel('Angular Velocity Error [rad/s]','Interpreter','latex', 'FontSize',18)
legend({'${\omega_x}$','${\omega_y}$','${\omega_z}$'},'Interpreter','latex','FontSize',16)

figure
plot(t_EKF, x_EKF(:,1), 'b-', 'LineWidth', 1.5)
hold on
plot(t_EKF, x_EKF_post(1,:), 'bo', 'LineWidth', 0.7)
plot(t_EKF, x_EKF(:,2), 'r-', 'LineWidth', 1.5)
plot(t_EKF, x_EKF_post(2,:), 'ro', 'LineWidth', 0.7)
plot(t_EKF, x_EKF(:,3), 'k-', 'LineWidth', 1.5)
plot(t_EKF, x_EKF_post(3,:), 'ko', 'LineWidth', 0.7)
plot(t_EKF, x_EKF(:,4), 'g-', 'LineWidth', 1.5)
plot(t_EKF, x_EKF_post(4,:), 'go', 'LineWidth', 0.7)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Time [$\mathrm{s}$]','Interpreter','latex', 'FontSize',18) 
ylabel('Quaternions','Interpreter','latex', 'FontSize',18)
legend({'${q_1}$ (Actual)','${q_1}$ (Estimated)','${q_2}$ (Actual)','${q_2}$ (Estimated)','${q_3}$ (Actual)','${q_3}$ (Estimated)','${q_4}$ (Actual)','${q_4}$ (Estimated)'},'Interpreter','latex','FontSize',16)

figure
plot(t_EKF, x_EKF(:,5), 'b-', 'LineWidth', 1.5)
hold on
plot(t_EKF, x_EKF_post(5,:), 'bo', 'LineWidth', 0.7)
plot(t_EKF, x_EKF(:,6), 'r-', 'LineWidth', 1.5)
plot(t_EKF, x_EKF_post(6,:), 'ro', 'LineWidth', 0.7)
plot(t_EKF, x_EKF(:,7), 'k-', 'LineWidth', 1.5)
plot(t_EKF, x_EKF_post(7,:), 'ko', 'LineWidth', 0.7)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
legend({'${\omega_x}$ (Actual)','${\omega_x}$ (Estimated)','${\omega_y}$ (Actual)','${\omega_y}$ (Estimated)','${\omega_z}$ (Actual)','${\omega_z}$ (Estimated)'},'Interpreter','latex', 'FontSize',16)
xlabel('Time [s]','Interpreter','latex', 'FontSize',18) 
ylabel('$\omega$ [rad/s]','Interpreter','latex', 'FontSize',18)

%% Controller Simulation 
% Load data from mission analysis team
load('mission_altitude.mat');
load('mission_time.mat');
% Function for mission profile data
[xData1, yData1] = prepareCurveData(mission_time, mission_altitude);
altitude_function = fit(xData1, yData1,'pchip');
dt = 0.1; % simulation time step
t = 0:dt:50; % simulation time span
r_init = 250000; % initial altitude
% Reaction wheel properties
theta = 0; % roll angle
beta = 30; % tilt angle
T_max = 0.02; % max torque
% Controller gains
K1 = 0.005*eye(3); 
K2 = 1.5*eye(3);
% Desired conditions
q_d = eul2quat([deg2rad(0), deg2rad(0), deg2rad(0)]); 
eta_d = q_d(1);
eps_d = [q_d(2), q_d(3), q_d(4)]'; 
w_d = [0, 0, 0]';
w_d_dot = [0, 0, 0]';
% Initial conditions
q_eq1 = eul2quat([deg2rad(0), deg2rad(0), deg2rad(0)]); 
q_eq = [q_eq1(2), q_eq1(3), q_eq1(4), q_eq1(1)]';
w_eq = [deg2rad(0), deg2rad(0), deg2rad(0)]'; 
x_0 = [q_eq; w_eq];
% Run simulation
x_actual = x_0;
start_point = find(mission_altitude >= r_init);
t0 = mission_time(start_point(1));
for i = 1:length(t)
    r(i) = r_init + altitude_function(t0 + t(i));
    [T_c(:,i), T_g(:,i), eta_err, eps_err, w_err, u(:,i), T_m(:,i)] = IBS(x_actual(4,i), x_actual(1:3,i), x_actual(5:7,i), eta_d, eps_d, w_d, w_d_dot, K1, K2, theta, beta, r(i), T_max, J);
    if t(i) < 1
        T_g(:,i) = [0.1,-0.1,0.05]';
    end
    x_actual(:,i+1) = xdot(x_actual(:,i), J, T_c(:,i), T_g(:,i))*dt + x_actual(:,i);
end
x_actual = x_actual(:,1:end-1);
eul_actual = 180/pi*quat2eul([x_actual(4,:)', x_actual(1,:)', x_actual(2,:)', x_actual(3,:)']);
stop = find(eul_actual(:,3) <= 1);
settling_time = t(stop(1));

figure
plot(t, eul_actual(:,3), 'b-', 'LineWidth', 1.5)
hold on
plot(t, eul_actual(:,2), 'r-', 'LineWidth', 1.5)
plot(t, eul_actual(:,1), 'k-', 'LineWidth', 1.5)
load('eul_actual_OWF.mat')
plot(t, eul_actual(:,3), 'b:', 'LineWidth', 2)
plot(t, eul_actual(:,2), 'r:', 'LineWidth', 2)
plot(t, eul_actual(:,1), 'k:', 'LineWidth', 2)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
legend({'X','Y','Z','X (OWF)','Y (OWF)','Z (OWF)'},'Interpreter','latex','FontSize',16) 
xlabel('Time [$\mathrm{s}$]','Interpreter','latex','FontSize',18)
ylabel('Euler Angles [${^{\circ}}$]','Interpreter','latex','FontSize',18)

figure
plot(t, x_actual(5,:), 'b-', 'LineWidth', 1.5)
hold on
plot(t, x_actual(6,:), 'r-', 'LineWidth', 1.5)
plot(t, x_actual(7,:), 'k-', 'LineWidth', 1.5)
load('x_actual_OWF.mat')
plot(t, x_actual(5,:), 'b:', 'LineWidth', 2)
plot(t, x_actual(6,:), 'r:', 'LineWidth', 2)
plot(t, x_actual(7,:), 'k:', 'LineWidth', 2)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
legend({'${\omega_x}$','${\omega_y}$','${\omega_z}$','${\omega_x}$ (OWF)','${\omega_y}$ (OWF)','${\omega_z}$ (OWF)'},'Interpreter','latex','FontSize',16) 
xlabel('Time [s]','Interpreter','latex','FontSize',18)
ylabel('Angular Velocities [rad/s] ','Interpreter','latex','FontSize',18)

figure
plot(t, x_actual(1,:), 'b-', 'LineWidth', 1.5)
hold on
plot(t, x_actual(2,:), 'r-', 'LineWidth', 1.5)
plot(t, x_actual(3,:), 'k-', 'LineWidth', 1.5)
plot(t, x_actual(4,:), 'g-', 'LineWidth', 1.5)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
legend({'q1','q2','q3','q4'},'Interpreter','latex','FontSize',16) 
xlabel('Time [s]','Interpreter','latex','FontSize',18)
ylabel('Quaternions','Interpreter','latex','FontSize',18)

figure
plot(t, T_c(1,:),'b-','LineWidth',1.5)
hold on
plot(t, T_c(2,:),'r-','LineWidth',1.5)
plot(t, T_c(3,:),'k-','LineWidth',1.5)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
legend({'${T_x}$','${T_y}$','${T_z}$'},'Interpreter','latex','FontSize',16)
xlabel('Time [s]','Interpreter','latex','FontSize',18)
ylabel('$T_{rw}$ [Nm]','Interpreter','latex','FontSize',18)

figure
plot(t, u(1,:), 'b-', 'LineWidth', 1.5)
hold on
plot(t, u(2,:), 'r-', 'LineWidth', 1.5)
plot(t, u(3,:), 'k-', 'LineWidth', 1.5)
plot(t, u(4,:), 'g-', 'LineWidth', 1.5)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
legend({'Wheel 1','Wheel 2','Wheel 3','Wheel 4'},'Interpreter','latex','FontSize',16)
xlabel('Time [s]','Interpreter','latex','FontSize',18) 
ylabel('RWs Torques [Nm]','Interpreter','latex','FontSize',18) 

%J_w = 1.59*10^-3*eye(4); %RWT150
J_w = 6.52*10^-4*eye(4); %RW100
%J_w = 1.9*10^-4*eye(4); %RW35
%J_w = 0.9*10^-4*eye(4); %RW25
for i = 2:length(t)
w_wheel(1,i) = trapz(t(1:i),u(1,1:i)/J_w(1,1)); 
w_wheel(2,i) = trapz(t(1:i),u(2,1:i)/J_w(2,2)); 
w_wheel(3,i) = trapz(t(1:i),u(3,1:i)/J_w(3,3)); 
w_wheel(4,i) = trapz(t(1:i),u(4,1:i)/J_w(4,4));
end
w_rpm = w_wheel *9.5493;

figure
plot(t, w_rpm(1,:), 'b-', 'LineWidth', 1.5)
hold on
plot(t, w_rpm(2,:), 'r-', 'LineWidth', 1.5)
plot(t, w_rpm(3,:), 'k-', 'LineWidth', 1.5)
plot(t, w_rpm(4,:), 'g-', 'LineWidth', 1.5)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
legend({'Wheel 1','Wheel 2','Wheel 3','Wheel 4'},'Interpreter','latex','FontSize',16,Location='northeast')
xlabel('Time [s]','Interpreter','latex','FontSize',18) 
ylabel('RWs Rotational Speed [RPM]','Interpreter','latex','FontSize',18)

for i = 2:length(t)
Energy(:,i-1) = abs(0.5*J_w*(w_rpm(:,i-1) ./ 9.5493).^2);
Power(:,i-1) = abs(u(:,i-1).*(w_rpm(:,i-1) ./ 9.5493));
end
Energy_tot = (trapz(t(1:end-1),Power(1,:)) + trapz(t(1:end-1),Power(2,:)) + trapz(t(1:end-1),Power(3,:)) + trapz(t(1:end-1),Power(4,:)));

figure
plot(t(1:end-1), Power(1,:), 'b-', 'LineWidth', 1.5)
hold on
plot(t(1:end-1), Power(2,:), 'r-', 'LineWidth', 1.5)
plot(t(1:end-1), Power(3,:), 'k-', 'LineWidth', 1.5)
plot(t(1:end-1), Power(4,:), 'g-', 'LineWidth', 1.5)
hold off
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
legend({'Wheel 1','Wheel 2','Wheel 3','Wheel 4'},'Interpreter','latex','FontSize',16) 
xlabel('Time [s]','Interpreter','latex','FontSize',18) 
ylabel('Power [W]','Interpreter','latex','FontSize',18)

% figure
% plot(settling_time(6), Energy_tot(6), 'bo', 'LineWidth', 1.5)
% hold on
% plot(settling_time, Energy_tot, 'b-', 'LineWidth', 1.5)
% hold off
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex') 
% xlabel('Settling Time [s]','Interpreter','latex','FontSize',18) 
% ylabel('Energy [J]','Interpreter','latex','FontSize',18) 
% legend('Design Point','Interpreter','latex','FontSize',16) 
% grid on
% saveas(gcf,'Energy_tuning.png')
% 
% figure
% plot(settling_time(6), Peak_Power(6), 'bo', 'LineWidth', 1.5)
% hold on
% plot(settling_time, Peak_Power, 'b-', 'LineWidth', 1.5)
% hold off
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex') 
% xlabel('Settling Time [s]','Interpreter','latex','FontSize',18) 
% ylabel('Peak Power [W]','Interpreter','latex','FontSize',18) 
% legend('Design Point','Interpreter','latex','FontSize',16) 
% grid on
% saveas(gcf,'Power_tuning.png')

%% 3D plots

% m = min(Energy_tot);
% M = min(m);
% [a,b] = find(Energy_tot == M);
% figure
% plot3(r_init/1000,110*ones(7,1),Energy_tot(:,8)','m-','LineWidth',1.5)
% hold on
% surf(r_init/1000,settling_time',Energy_tot')
% shading interp
% hold off
% grid on
% box on
% colormap(jet)
% C=colorbar;
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('$h_0$ [km]','Interpreter','latex','FontSize',18)
% ylabel('$t_s$ [s]','Interpreter', 'latex', 'FontSize',18)
% zlabel('Energy [J]','Interpreter', 'latex', 'FontSize',18)
% legend({'Minimum point'},'Interpreter','latex','FontSize',16) 
% view(40,40);
% saveas(gcf,'Energy_3D.png')
% 
% % n = min(Peak_Power);
% % N = min(n);
% % [c,d] = find(Peak_Power == N);
% figure
% % plot3(beta(d),theta(c),Peak_Power(c,d),'md',MarkerFaceColor='m')
% % hold on
% surf(r_c_init/1000,settling_time',Peak_Power')
% shading interp
% % hold off
% grid on
% box on
% colormap(jet)
% C=colorbar;
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('$h_0$ [km]','Interpreter','latex','FontSize',18)
% ylabel('$t_s$ [s]','Interpreter', 'latex', 'FontSize',18)
% zlabel('Peak Power [W]','Interpreter', 'latex', 'FontSize',18)
% %legend({'Minimum point'},'Interpreter','latex','FontSize',16,Location='northeast')
% % ylabel(C,'Number of Ribs','Interpreter', 'latex', 'FontSize',18)
% view(40,40);
% saveas(gcf,'Peak_Power_3D.png')
