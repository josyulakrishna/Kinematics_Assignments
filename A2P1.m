% Assignment 2 - Question 1, Forward Kinematics
clear all; 
%lengths
l0 = 1;
l1 = 1;
l2 = 1;
%M 
M = [1, 0, 0, 0; 
    0, 1, 0, l1+l2; 
    0, 0, 1, l0; 
    0, 0, 0, 1];
%Screw Axes 
S1 = [0,0,1,0,0,0];
S2 = [0,0,1,l1,0,0];
S3 = [0,0,1,l1+l2,0,0];
S4 = [0,0,0,0,0,1];
Slist = [S1', S2', S3', S4'];
%Body Axes 
B1 = [0,0,1,-(l1+l2),0,0];
B2 = [0,0,1,-l2,0,0];
B3 = [0,0,1,0,0,0];
B4 = [0,0,0,0,0,1];
Blist = [B1', B2', B3', B4'];
%thetalist 
thetalist = [0;pi/2;-pi/2;1];
%Configuration Space frame 
Ts = FKinSpace(M,Slist,thetalist);
%Configuration Body Frame 
Tb = FKinBody(M,Blist,thetalist);

% OUTPUT
% Ts =
% 
%      1     0     0    -1
%      0     1     0     1
%      0     0     1     2
%      0     0     0     1
% 
% 
% Tb =
% 
%      1     0     0    -1
%      0     1     0     1
%      0     0     1     2
%      0     0     0     1

