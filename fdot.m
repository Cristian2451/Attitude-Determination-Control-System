function f_dot = fdot(t, x, J, K1, K2, eta_d, eps_d, w_d, w_d_dot)
% Error Dynamics & Disturbances
eta = x(4); 
eps = x(1:3); 
w = x(5:7);
[T_rw, T_gg, eta_err, eps_err, w_err,u] = IBS(J, eta, eps, w, eta_d, eps_d, w_d, w_d_dot, K1, K2);
% if t >= 100 && t <= 110 
%     T_gg = T_gg + [2;2;2];
% end
f_dot(1:3,:) = 0.5 * (eta_err * eye(length(eps_err)) + cp_op(eps_err))*w_err;
f_dot(4,:) = -0.5 * eps_err' * w_err;
f_dot(5:7,:) = J \ (J*w_d_dot + cross(w, J*w) - T_rw - T_gg);
end