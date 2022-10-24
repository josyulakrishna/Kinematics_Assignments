function [coefftheta1,coefftheta2] = trajGen(theta1_bvp,theta2_bvp)
%TRAJGEN takes the bvp {(positions, time),(velocities,time),(accelerations,time)}
% %   Detailed explanation goes here
% example inputs: 
% theta1_bvp = struct; 
% theta1_bvp.positions = [-pi/4,0; pi/4,2];
% theta1_bvp.velocities = [0,0; 0,2];
% theta1_bvp.accels = [0,0; 0,2];
% theta2_bvp = struct; 
% theta2_bvp.positions = [0,0; pi/2,2];
% theta2_bvp.velocities = [0,0; 0,2];
% theta2_bvp.accels = [0,0; 0,2];

syms a1 b1 c1 d1 e1 f1 
syms a2 b2 c2 d2 e2 f2 
% for positions 
eqn_p={};
for i=1:size(theta1_bvp.positions,1)
    eqn_p{end+1}=a1*theta1_bvp.positions(i,2)^5+b1*theta1_bvp.positions(i,2)^4+c1*theta1_bvp.positions(i,2)^3+d1*theta1_bvp.positions(i,2)^2+e1*theta1_bvp.positions(i,2)^1+f1==theta1_bvp.positions(i,1)

end

