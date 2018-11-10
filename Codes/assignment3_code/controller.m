function [F, M] = controller(t, state, des_state, params)
%CONTROLLER  Controller for the quadrotor
%
%   state: The current state of the robot with the following fields:
%   state.pos = [x; y; z], state.vel = [x_dot; y_dot; z_dot],
%   state.rot = [phi; theta; psi], state.omega = [p; q; r]
%
%   des_state: The desired states are:
%   des_state.pos = [x; y; z], des_state.vel = [x_dot; y_dot; z_dot],
%   des_state.acc = [x_ddot; y_ddot; z_ddot], des_state.yaw,
%   des_state.yawdot
%
%   params: robot parameters

%   Using these current and desired states, you have to compute the desired
%   controls


% =================== Your code goes here ===================

% Thrust
F = 0;

% Moment
M = zeros(3,1);

Kpz = 10;
Kpx = 5;
Kpy= 5;
Kdz = 2;
Kdx = 2;
Kdy= 2;

% des_state.yaw= des_state.yaw+3;
% des_state.yawdot = des_state.yawdot+1;

z_acc = des_state.acc(3)+ Kpz*(des_state.pos(3)-state.pos(3)) + Kdz*(des_state.vel(3)-state.vel(3));
x_acc = des_state.acc(1)+ Kpx*(des_state.pos(1)-state.pos(1)) + Kdx*(des_state.vel(1)-state.vel(1));
y_acc = des_state.acc(2)+ Kpy*(des_state.pos(2)-state.pos(2)) + Kdy*(des_state.vel(2)-state.vel(2));

F = params.mass*(z_acc+params.gravity);
eul_des=1/params.gravity*([x_acc, -y_acc;...
                           y_acc x_acc]...
                              *[sin(des_state.yaw) cos(des_state.yaw)].');
phi_des = eul_des(1);
theta_des= eul_des(2);
                          
  K= [200 10;200 10;200 10];

r_des = [sin(state.rot(2)) 0 cos(state.rot(1))*cos(state.rot(2))]*[0 0 des_state.yawdot].';

E = [phi_des-state.rot(1) 0-state.omega(1);...
     theta_des-state.rot(2) 0-state.omega(2);...
     des_state.yaw-state.rot(3)  r_des-state.omega(3)];

 M = sum(K.*E, 2);
% =================== Your code ends here ===================

end
