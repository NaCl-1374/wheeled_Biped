clear;
clear all;



%% %仿真参数
runtime=10;%仿真时长
delta_t = 0.0005;%仿真时间间隔
ControlHz=1000;
IterControlRate=fix(1/ControlHz/delta_t);%mpc频率控制
t=0:delta_t:runtime-delta_t;

%% 
%仿真循环
State_X=zeros(6,1);
State_X(1)=0;
State_X(5)=0.2;

log_x=[];
log_U=[];


%%

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
%%
Q=diag([1,1,1000,100,5000,1]);
R=diag([50,10]);
[K,P]=dlqr(A_d,B_d,Q,R);
[Kc,Pc]=lqr(A,B,Q,R);
U=[0,0]';
%%
for iter = 1:runtime/delta_t
    %控制器频率控制
    
    if rem(iter,IterControlRate)==0
    U = control_LQR(State_X,Kc);
    end
    log_x=[log_x,State_X];
    log_U=[log_U,U];

    State_X=State_X+sfuntmpl(State_X,U)*delta_t;
    
    
end
subplot(2,2,1)
plot(t,log_x(1,:),'r')
hold on
plot(t,log_x(2,:),'g')
legend('th','dth')
title('theta')
subplot(2,2,2)
plot(t,log_x(3,:),'g')
hold on
plot(t,log_x(4,:),'r')
legend('x','dx')
title('X')
subplot(2,2,3)
plot(t,log_x(5,:),'b')
hold on
plot(t,log_x(6,:),'r')
legend('fai','dfai')
title('fai')
subplot(2,2,4)
plot(log_U(1,:),'b')
hold on
plot(log_U(2,:),'r')
legend('T轮子','T大腿')
title('U')


function DX=sfuntmpl(X,U)

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

th=X(1);
dth=X(2);
x=X(3);
dx=X(4);
fai=X(5);
dfai=X(6);
T=U(1);
Tp=U(2);

ddfai =(Tp + M*g*l*sin(fai) - L*M*dth^2*l*sin(fai + th) - LM*M*dth^2*l*sin(fai + th))/(M*l^2 + IM);
ddth =(Tp - T + L*M*g*sin(th) + LM*M*g*sin(th) + L*g*mp*sin(th) - L*M*dfai^2*l*sin(fai + th) - LM*M*dfai^2*l*sin(fai + th))/(Ip + L^2*M + LM^2*M + L^2*mp + 2*L*LM*M);
ddx =(T + R*(M*(- l*sin(fai)*dfai^2 + sin(th)*(L + LM)*dth^2) + L*dth^2*mp*sin(th)))/(((R*(M + mp))/(R*mw + Iw/R) + 1)*(R*mw + Iw/R));

DX=[dth,ddth,dx,ddx,dfai,ddfai]';
end

function U = control_LQR(X,K)

% K=[-518.630952756606,-163.865462396828,-0.918278733001755,-3.02651334530555,83.2139675667216,19.0724987540786;-107.009642498362,-33.8105176313671,-0.395934550794602,-1.16875605176547,95.7004883416066,22.1728131144323];
U=-(K*X);
end
