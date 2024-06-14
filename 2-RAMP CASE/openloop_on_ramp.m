%Traffic Modelling and Control of a STM corridor

clear all; close all

stm_alpha = 0.0167; 

stm_ks   = 0.031;  % k(s-2) 
stm_ks2  = 0.010;  % k(s-1) 
stm_ks3  = 0.033;  % k(s)  
stm_ks4  = 0.024; % k(s+1)
stm_ks5  = 0.025; % k(s+2)

stm_ftd=[7 100 120 60   1100;
         7 95 120 60   1000;
         7 95  100 53.8   1000;
         7 100 120 60   3600;
         7 100 120 63   1600;
         7 100 120 63   1600;];

%initialisation of input to simulink model
dt=1; dts=10;       % sampling rate 10 secs
SIM=183;            % simulation length (minutes)
SIM=SIM*60/dts;     % simulation length (samples)
t_in=[0:1:SIM-1]';  % time (samples)
z0=0; 

r6 = [10*ones(300,1);  140*ones(500,1);  10*ones(298,1)];
p6 = [5 100 110 63 1000]; v_in=ftdcurve(stm_ftd(6,:),r6,1);

q_in=[10*ones(150,1); linspace(10,70,150)'; 70*ones(500,1); linspace(70,10,150)'; 10*ones(148,1)]; 
q_in=0.80*q_in;  

%Ramp Input
q_ramp=[3*ones(200,1); linspace(3,30,200)'; 30*ones(400,1); linspace(30,3,150)'; 3*ones(148,1);];

sim('openloop_on_ramp_mdl')

%Collect all values
q1=qq(:,1); q2=qq(:,2); q3=qq(:,3); q4=qq(:,4); q5=qq(:,5); q6=qq(:,6);
r1=dd(:,1); r2=dd(:,2); r3=dd(:,3); r4=dd(:,4); r5=dd(:,5); 
v1=vv(:,1); v2=vv(:,2); v3=vv(:,3); v4=vv(:,4); v5=vv(:,5); 

figure(1);
subplot(321); plot(r1,q1,'.'); hold; plot(stm_ftd(1,4),linspace(0,120),'k.'); plot(stm_ftd(1,4),max(q1),'r.');title('s-2');
subplot(322); plot(r2,q2,'.'); hold; plot(stm_ftd(2,4),linspace(0,120),'k.'); plot(stm_ftd(2,4),max(q2),'r.');title('s-1 {load}');
subplot(323); plot(r3,q3,'.'); hold; plot(stm_ftd(3,4),linspace(0,120),'k.'); plot(stm_ftd(3,4),max(q3),'r.');title('s {junction}');
subplot(324); plot(r4,q4,'.'); hold; plot(stm_ftd(4,4),linspace(0,120),'k.'); plot(stm_ftd(4,4),max(q4),'r.');title('s+1');
subplot(325); plot(r5,q5,'.'); hold; plot(stm_ftd(5,4),linspace(0,120),'k.'); plot(stm_ftd(5,4),max(q5),'r.');title('s+2');

figure(2);
x=[0:1:500]';y_1=ftdcurve(stm_ftd(1,:),x,1);subplot(231);plot(x,y_1); hold; axis([0 500 0 150]);title('s-2');
x=[0:1:500]';y_2=ftdcurve(stm_ftd(2,:),x,1);subplot(232);plot(x,y_2); hold; axis([0 500 0 150]) ;title('s-1 {load}');
x=[0:1:500]';y_3=ftdcurve(stm_ftd(3,:),x,1);subplot(233);plot(x,y_3); hold; axis([0 500 0 150]) ;title('s {junction}');
x=[0:1:500]';y_4=ftdcurve(stm_ftd(4,:),x,1);subplot(234);plot(x,y_4); hold; axis([0 500 0 150]);title('s+1'); 
x=[0:1:500]';y_5=ftdcurve(stm_ftd(5,:),x,1);subplot(235);plot(x,y_5); hold; axis([0 500 0 150]);title('s+2');
subplot(231);plot(r1,v1,'.'); plot(stm_ftd(1,4),stm_ftd(1,2),'r.');title('s-2')
subplot(232);plot(r2,v2,'.'); plot(stm_ftd(2,4),stm_ftd(2,2),'r.');title('s-1')
subplot(233);plot(r3,v3,'.'); plot(stm_ftd(3,4),stm_ftd(3,2),'r.');title('s')
subplot(234);plot(r4,v4,'.'); plot(stm_ftd(4,4),stm_ftd(4,2),'r.');title('s+1')
subplot(235);plot(r5,v5,'.'); plot(stm_ftd(5,4),stm_ftd(5,2),'r.');title('s+2')

figure(11);plot([r3 r2]); hold; ylabel('density @ junction r(s,t) [veh/km]'); xlabel('Simulation Period [min]'); 

UnmLoad=r2;
flow_ups=q2;
flow_jun=q3;
old_v3=v3;
ramp_load=sum(q_ramp(450:900));


figure(14); plot(q3);hold; plot(q_ramp,'r.'); plot(q_in,'k.')
xlabel('time {min}'); ylabel('Flow {veh/min}');

disp('when cap occurs;'); index=find(q3==max(q3))
disp('load'); max(r2)
disp('junction'); max(r3)


%time delays
r1_1d = lag(r1,1,mean(r1)); r2_1d = lag(r2,1,mean(r1)); r2_2d = lag(r2,2,mean(r2)); r3_1d = lag(r3,1,mean(r3)); r3_2d = lag(r3,2,mean(r3)); r4_1d = lag(r4,1,mean(r4)); r4_2d = lag(r4,2,mean(r4)); r5_1d = lag(r5,1,mean(r5)); r5_2d = lag(r5,2,mean(r5));
v2_2d = lag(v2,2,mean(v2)); v3_1d = lag(v3,1,mean(v3)); v3_2d = lag(v3,2,mean(v3)); v4_1d = lag(v4,1,mean(v4)); v4_2d = lag(v4,2,mean(v4)); v5_2d = lag(v5,2,mean(v5));




% ------------------------------- LLM ---------------------------------
%Define an F matrix as follows:
Z1=[r2 r1_1d r3_1d];         %for s-1 
Z2=[r3 r2_1d r4_1d q_ramp] ; %for s
Z3=[r4 r3_1d ];              %for s+1 

Z1 = Z1(3:end,1:end); Z2 = Z2(3:end,1:end); Z3 = Z3(3:end,1:end);  dtc=1;  
[z1, m1, s1]=prepz(Z1, [], 1, 10, 0, dtc);  % 1 minute sampling
[z2, m2, s2]=prepz(Z2, [], 1, 10, 0, dtc);  % 1 minute sampling
[z3, m3, s3]=prepz(Z3, [], 1, 10, 0, dtc);  % 1 minute sampling
z1=Z1; z2=Z2; z3=Z3;

%============== link s-1
nn=[1 1 1 1 1 0];
[th, stats]=riv(z1, nn);%, [1 0 0 0 0]);  % least squares estimation
RT2=stats(3); YIC=stats(2); [at, bt]=getpar(th);
RT2
YIC

%============== link s
nn= [1 1 1 1 1 1 1 0];
[th2, stats2]=riv(z2, nn);%, [1 0 0 0 0]);  % least squares estimation
RT2_2=stats2(3); YIC_2=stats2(2); [at_2, bt_2]=getpar(th2);
RT2_2
YIC_2

%=============== link s+1
nn=[1 1 1 0];
[th3, stats3]=riv(z3, nn);%, [1 0 0 0 0]);
RT2_3=stats3(3); YIC_3=stats3(2); [at_3, bt_3]=getpar(th3);
RT2_3
YIC_3

%================ Extraction of model parameters =======================
aa=[-at(2) -at_2(2) -at_3(2)]
bb=[bt(1, 2) bt_2(1, 2) bt_3(1, 2)]
cc=[0 bt_2(2, 2) 0]
dd=[0 bt_2(3,2) 0]

%Formulation of the NMSS matrix
F = [aa(1)  bb(1)      0  0;...
        bb(2)  aa(2)   cc(2) 0;...
        0   bb(3)   aa(3) 0;...
        -bb(2) -aa(2)  -cc(2) 1;];
g = [0 dd(2) 0 -dd(2)]';
% 
% F = [aa(2)  bb(2)  cc(2) 0;...
%      cc(1)  aa(1)   0    0;...
%      bb(3)   0     aa(3) 0;...
%     -aa(2) -bb(2) -cc(2) 1];
% 
% g = [dd(2) 0  0 -dd(2)]';

Q=eye(size(F)); Q(end,end)=1; r=1; k_llm = dlqr(F,g,Q,r); k_llm(end)=-k_llm(end)

