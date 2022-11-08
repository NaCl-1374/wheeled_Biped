clear all;
R=0.25/2;
L =0.5;  %摆杆重心到驱动轮轴距离
LM =0.5; %摆杆重心到机体转轴距离
l =0.2;  %机体重心到机体转轴距离
mw =1;   %转子质量
mp =2;   %摆杆质量
M=20;    %机体质量
Iw =0.5*mw*(R^2+(R-0.03)^2);%转子惯量
Ip =1/4*mp*0.08^2+1/12*mp*L^2;%摆杆惯量
IM =1/4*M*0.3^2+1/12*M*0.7^2;%机体惯量
g =9.8;
% 来自calculate_robot_model.m
B=[ 0,0;
    -1/(Ip + L^2*M + LM^2*M + L^2*mp + 2*L*LM*M), 1/(Ip + L^2*M + LM^2*M + L^2*mp + 2*L*LM*M);
    0,0;
    R/(Iw + M*R^2 + R^2*mp + R^2*mw),0;
    0,0;
    0,1/(M*l^2 + IM)];

A=[ 0,1,0,0,0,0;
    (g*(L*mp + L*M + LM*M))/(Ip + L^2*M + LM^2*M + L^2*mp + 2*L*LM*M), 0, 0, 0,0, 0;
    0,0,0,1,0,0;
    0, 0, 0, 0,0, 0;
    0,0,0,0,0,1;
    0, 0, 0, 0, (M*g*l)/(M*l^2 + IM), 0];
dt=0.001;
C=eye(6);
D=zeros(6,2);
sys_c = ss(A, B, C, D);            
sys_d = c2d(sys_c, dt);
[A_d, B_d, C_d, D_d] = ssdata(sys_d);
Q=diag([10,1,1,1,10,8]);
R=diag([0.1,0.1]);
%离散
[K,P]=dlqr(A_d,B_d,Q,R);
%连续
[Kc,Pc]=lqr(A,B,Q,R);
