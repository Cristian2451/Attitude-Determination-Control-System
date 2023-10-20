function [T, T_gg, eta_e, eps_e, w_e, u, T_m] = IBS(eta, eps, w, eta_d, eps_d, w_d, w_d_dot, K1, K2, theta, beta, r, T_max, J)
% Magnetic torque
[b0] = wrldmagm(200000, 35, -115.28, decyear(2023,6,20))*10^(-9);
A = quat2dcm([eta,eps(1),eps(2),eps(3)]);
B = A*b0;
m = [0.1,0.1,0.1]';
T_m = m.*B;
% Gravity gradient torque
mu = 3.986e14;
w_o = sqrt(mu/((r+6371000)^3));
eul = 180/pi*quat2eul([eta, eps(1), eps(2), eps(3)]);
T_gg = zeros(3,1);
T_gg(1,1) = 3*w_o^2*(J(3,3)-J(2,2))*cosd(eul(2))^2*sind(2*eul(3));
T_gg(2,1) = 3*w_o^2*(J(3,3)-J(1,1))*cosd(eul(3))*sind(2*eul(2));
T_gg(3,1) = 3*w_o^2*(J(2,2)-J(1,1))*sind(eul(3))*sind(2*eul(2));
% Control law
eta_e = eta_d*eta + eps_d'*eps;
eps_e = eta*eps_d - eta_d*eps - cp_op(eps_d)*eps; 
w_e = w_d - w;
G = [sign(eta_e)*eps_e'; eta_e*eye(length(eps_e)) + cp_op(eps_e)]'; 
z1 = [1 - abs(eta_e); eps_e];
alpha1 = -K1*G*z1;
z2 = w_e - alpha1;
z1_dot = 0.5*G'*(alpha1 + z2); 
eta_e_dot = -0.5 * eps_e' * w_e;
eps_e_dot = 0.5 * (eta_e * eye(length(eps_e)) + cp_op(eps_e))*w_e;
G_dot = [sign(eta_e)*eps_e_dot'; eta_e_dot*eye(length(eps_e_dot)) + cp_op(eps_e_dot)]';
alpha1_dot = -K1*(G_dot*z1 + G*z1_dot);
% Control torque
T = (G*z1 + K2*z2 - J*w_d_dot - cross(w, J*w) + T_gg + J*alpha1_dot);
% Reaction wheel placement
L = [cosd(beta)*cosd(theta), -cosd(beta)*cosd(theta), -cosd(beta)*cosd(theta), cosd( beta)*cosd(theta);
cosd(beta)*cosd(theta), cosd(beta)*cosd(theta), -cosd(beta)*cosd(theta), -cosd( beta)*cosd(theta);
sind(beta), sind(beta), sind(beta), sind(beta)];
T_p = L'/(L*L');
%Reaction wheel torques
u = eye(4)\T_p*T;
% Saturation
if abs(u(1)) > T_max
u(1) = sign(u(1))*T_max;
end
if abs(u(2)) > T_max
u(2) = sign(u(2))*T_max;
end
if abs(u(3)) > T_max
u(3) = sign(u(3))*T_max;
end
if abs(u(4)) > T_max
u(4) = sign(u(4))*T_max;
end 
% Wheel failures
%u(1) = 0;

% Update control torque
T = T_p\eye(4)*u;
end