% Assignment 2 - Question 2, Forward Kinematics
clear all;
%M 
M = [-1, 0, 0, 0; 
    0, 0, 1, 3; 
    0, 1, 0, 2; 
    0, 0, 0, 1];
%Screw Axes w,v
S1 = [0,0,1, 0,0,0];
S2 = [1,0,0, 0,2,0];
S3 = [0,0,0, 0,1,0];

Slist = [S1',S2',S3'];
%Body Axes 
B1 = [0,1,0,3,0,0];
B2 = [-1,0,0,0,3,0];
B3 = [0,0,0,0,0,1];

Blist = [B1',B2',B3'];
%thetalist 
thetalist = [pi/2;pi/2;1];
%Configuration Space frame 
Ts = FKinSpace(M,Slist,thetalist);
%Configuration Body Frame 
Tb = FKinBody(M,Blist,thetalist);

% OUTPUT: 
% 
% 
% Tb =
% 
%    -0.0000    1.0000   -0.0000   -0.0000
%    -1.0000   -0.0000    0.0000    0.0000
%          0    0.0000    1.0000    6.0000
%          0         0         0    1.0000
% 
% 
% Ts =
% 
%    -0.0000    1.0000   -0.0000   -0.0000
%    -1.0000   -0.0000    0.0000    0.0000
%          0    0.0000    1.0000    6.0000
%          0         0         0    1.0000
