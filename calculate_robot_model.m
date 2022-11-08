clear all;
syms t T Tp N P Nm Pm Nf R L LM l mw mp Iw Ip IM g M
syms ddot_x ddot_theta ddot_fai dot_fai dot_x dot_theta
syms x(t) theta(t) fai(t) %T(t) Tp(t)
syms th xx ffai
%% 基本方程
eq01=mw*diff(x,t,2)==Nf-N;
eqNf=rhs(isolate(eq01,Nf));
eq02=Iw*diff(x,t,2)/R==T-Nf*R;
eq=isolate(subs(eq02,Nf,eqNf),diff(x,t,t));
% eq=diff(x,t,2)==(T-N*R)/(Iw/R+mw*R);
eq1= N-Nm==mp*diff(x+L*sin(theta),t,2);
eq2= P-Pm-mp*g==mp*diff(L*cos(theta),t,2);
eq21= Ip*diff(theta,t,2)==(P*L+Pm*LM)*sin(theta)-(N*L+Nm*LM)*cos(theta)-T+Tp;

eq3= Nm==M*diff((L+LM)*sin(theta)-l*sin(fai)+x,t,2);
eq4= Pm-M*g==M*diff((L+LM)*cos(theta)+l*cos(fai),t,2) ;
eq41= IM*diff(fai,t,2)==Tp+Nm*l*cos(fai)+Pm*l*sin(fai);
%%
eqN=rhs(isolate(eq1+eq3,N));%求N
eqP=rhs(isolate(eq2+eq4,P));%求P
eqNm=rhs(isolate(eq3,Nm));%求Nm
eqPm=rhs(isolate(eq4,Pm));%求Pm
%%
% pretty(eqN);
% pretty(eqP);
% pretty(eqNm);
% pretty(eqPm);
%%
eq=subs(eq,[P,N,Pm,Nm],[eqP,eqN,eqPm,eqNm]);
%%
eq21=subs(eq21,[P,N,Pm,Nm],[eqP,eqN,eqPm,eqNm]);
%%
eq41=subs(eq41,[Nm,Pm],[eqNm,eqPm]);
%%
DDX=isolate(eq,diff(x(t), t, t));
DDth=isolate(eq21,diff(theta(t), t, t));
DDFAI=isolate(eq41,diff(fai(t), t, t));

ddx=rhs(isolate(eq,diff(x(t), t, t)));
ddth=rhs(isolate(eq21,diff(theta(t), t, t)));
ddfai=rhs(isolate(eq41,diff(fai(t), t, t)));

ddx=simplify(collect(ddx));
ddth=simplify(collect(ddth));
ddfai=simplify(collect(ddfai));


% pretty(simplify(collect(ddth)))

%% 输出

ddx=subs(ddx,[theta(t),diff(theta(t), t),x(t),diff(x(t), t),fai(t),diff(fai(t), t)]...
    ,[th,dot_theta,xx,dot_x,ffai,dot_fai]);
ddth=subs(ddth,[theta(t),diff(theta(t), t),x(t),diff(x(t), t),fai(t),diff(fai(t), t)]...
    ,[th,dot_theta,xx,dot_x,ffai,dot_fai]);
ddfai=subs(ddfai,[theta(t),diff(theta(t), t),x(t),diff(x(t), t),fai(t),diff(fai(t), t)]...
    ,[th,dot_theta,xx,dot_x,ffai,dot_fai]);

jac=jacobian([ddth;ddx;ddfai],[th,dot_theta,xx,dot_x,ffai,dot_fai]);
simplify(subs(jac,[th,dot_theta,xx,dot_x,ffai,dot_fai,T,Tp],[0,0,xx,0,0,0,0,0]))

jacU=jacobian([ddth;ddx;ddfai],[T,Tp]);
simplify(subs(jacU,[th,dot_theta,xx,dot_x,ffai,dot_fai,T,Tp],[0,0,xx,0,0,0,0,0]))
