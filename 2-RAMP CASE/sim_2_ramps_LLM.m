close all
clear all

global F g r


%--------------------------------------------------------------------
% STM simulation
k1  = 0.0182;    % k(s-2) 
k2  = 0.026;   % k(s-1) 
k3  = 0.024;   % k(s)  
k4  = 0.022;   % k(s+1)
k5  = 0.028;   % k(s+2)
k6  = 0.025;   % k(s+3)
k7  = 0.026;   % k(s+4)
k8  = 0.016;   % k(s+5)

alpha1=0.0167; alpha2=alpha1; alpha3=alpha1; alpha4=alpha1; alpha5=alpha1; alpha6=alpha1; alpha7=alpha1; alpha8=alpha1;

p1=[7 100 110 63 360]; p2=p1; p3=p1; p4=p1; p5=p1; p6=p1; p7=p1; p8=p1; p9=p1;

% Set boundary conditions
kappa=linspace(3,70,50); kappa2=linspace(60,10,40); kappa3=linspace(70,60,50); %kappa4=linspace(5,15,40); kappa5=linspace(15,5,40);
q_in=[3*ones(400,1); kappa'; 70*ones(350,1); kappa3'; 60*ones(250,1); kappa2'; 10*ones(300,1)];  
v_in=[55*ones(1440,1);]; v9=v_in; r9=ftdcurve(p9,v9,2); 

% initialise other variables
r1=zeros(size(q_in)); r2=zeros(size(q_in)); r3=zeros(size(q_in)); r4=zeros(size(q_in)); r5=zeros(size(q_in)); r6=zeros(size(q_in)); r7=zeros(size(q_in));
v1=zeros(size(q_in)); v2=zeros(size(q_in)); v3=zeros(size(q_in)); v4=zeros(size(q_in)); v5=zeros(size(q_in)); v6=zeros(size(q_in)); v7=zeros(size(q_in));
q1=zeros(size(q_in)); q2=zeros(size(q1)); q3=zeros(size(q1)); q4=zeros(size(q1)); q5=zeros(size(q1)); q6=zeros(size(q1)); q7=zeros(size(q1));

% === Open Loop Simulation ===
r2(1)=0.1; r3(1)=0.1; r4(1)=0.1; r5(1)=0.1; r6(1)=0.1; r7(1)=0.1; r8(1)=0.1; r9(1)=0.1; q2(1)=0.1; q3(1)=0.1; q4(1)=0.1; q5(1)=0.1; q6(1)=0.1; q7(1)=0.1; q8(1)=0.1; q9(1)=0.1; 
v2(1)=ftdcurve(p2,r2(1),1); v3(1)=ftdcurve(p3,r3(1),1); v4(1)=ftdcurve(p4,r4(1),1); v5(1)=ftdcurve(p5,r5(1),1); v6(1)=ftdcurve(p6,r6(1),1); v7(1)=ftdcurve(p7,r7(1),1); v8(1)=ftdcurve(p8,r8(1),1); v9(1)=ftdcurve(p9,r9(1),1);

q_ramp=13*ones(size(q1)); q_ramp2=13*ones(size(q1));

h=waitbar(0,'Open Loop Simulation  Please wait... ');

for n=2:length(q1);
    
    %link 2 [s-2]
    r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
    v1(n)=ftdcurve(p1,r1(n),1);
    
    %link 2 [s-1]
    r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
    v2(n)=ftdcurve(p2,r2(n),1);
        
    %Ramp [1]    
    r_ramp(n)=r_ramp(n-1) + ((q_urban(n-1)-q_ramp(n-1))* kramp);    
    q_ramp(n)=alpha1*min(v_ramp(n-1),v3(n-1))*r_ramp(n-1);
    v_ramp(n)=ftdcurve(p7,r_ramp(n),1);
    
    %link 3 [s]    ------- JUNCTION 1
    r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
    v3(n)=ftdcurve(p3,r3(n),1);
    
    %link 4 [s+1]
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
    v4(n)=ftdcurve(p4,r4(n),1);
    
    %link 5 [s+2]   
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    v5(n)=ftdcurve(p5,r5(n),1);
    
    %Ramp [2]    
    r_ramp2(n)=r_ramp2(n-1) + ((q_urban2(n-1)-q_ramp2(n-1))* kramp2);    
    q_ramp2(n)=alpha1*min(v_ramp2(n-1),v6(n-1))*r_ramp2(n-1);
    v_ramp2(n)=ftdcurve(p8,r_ramp2(n),1);
    
    %link 6 [s+3]  ------- JUNCTION 2
    r6(n)=r6(n-1)+((q_ramp2(n-1)+q5(n-1)-q6(n-1))*k6);
    q6(n)=alpha6*min(v6(n-1),v7(n-1))*r6(n-1);
    v6(n)=ftdcurve(p6,r6(n),1);

    %link 7 [s+4]
    r7(n)=r7(n-1)+((q6(n-1)-q7(n-1))*k7);
    q7(n)=alpha7*min(v7(n-1),v8(n-1))*r7(n-1);
    v7(n)=ftdcurve(p7,r7(n),1);

    %link 8 [s+5]
    r8(n)=r8(n-1)+((q7(n-1)-q8(n-1))*k8);
    q8(n)=alpha8*min(v8(n-1),v9(n-1))*r8(n-1);
    v8(n)=ftdcurve(p8,r8(n),1);

    waitbar(n/length(q1),h);
end
close(h)

fig('Slacked Plots'); plot([r2 r3 r4 r5 r6 r7 r8']); legend('r2','r3','r4','r5','r6','r7','r8');

 figure(5); 
o1=occ2den(r1,4,2,3,0);     subplot(421); plot(q1,o1,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');
o2=occ2den(r2,4,2,3,0);     subplot(422); plot(q2,o2,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample'); 
o3=occ2den(r3,4,2,3,0);     subplot(423); plot(q3,o3,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');
o4=occ2den(r4,4,2,3,0);     subplot(424); plot(q4,o4,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');
o5=occ2den(r5,4,2,3,0);     subplot(425); plot(q5,o5,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');
o6=occ2den(r6,4,2,3,0);     subplot(426); plot(q6,o6,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');
o7=occ2den(r7,4,2,3,0);     subplot(427); plot(q7,o7,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');
o8=occ2den(r8,4,2,3,0);     subplot(428); plot(q8,o8,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');

break
oramp1=occ2den(r_ramp,4,2,3,0); subplot(426); plot(q_ramp,oramp1,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');










figure(2); subplot(211);hold; plot(r3);title('Junction1'); subplot(212);hold; plot(r6);title('Junction2');

r1_1d = lag(r1,1,mean(r1)); r1_2d = lag(r1,2,mean(r1)); r2_1d = lag(r2,1,mean(r2)); r2_2d = lag(r2,2,mean(r2)); r3_1d = lag(r3,1,mean(r3)); r3_2d = lag(r3,2,mean(r3)); r4_1d = lag(r4,1,mean(r4)); r4_2d = lag(r4,2,mean(r4));
r5_1d = lag(r5,1,mean(r5)); r5_2d = lag(r5,2,mean(r5)); r6_1d = lag(r6,1,mean(r6)); r6_2d = lag(r6,2,mean(r6)); r7_1d = lag(r7,1,mean(r7)); r7_2d = lag(r7,2,mean(r7)); r8_1d = lag(r8',1,mean(r8));r8_2d = lag(r8',2,mean(r8));
v1_1d = lag(v1,1,mean(v1)); v1_2d = lag(v1,2,mean(v1)); v2_1d = lag(v2,1,mean(v2)); v2_2d = lag(v2,2,mean(v2)); v3_1d = lag(v3,1,mean(v3)); v3_2d = lag(v3,2,mean(v3)); v4_1d = lag(v4,1,mean(v1)); v4_2d = lag(v4,2,mean(v4));
v5_1d = lag(v5,1,mean(v5)); v5_2d = lag(v5,2,mean(v5)); v6_1d = lag(v3,1,mean(v6)); v6_2d = lag(v6,2,mean(v6)); v7_1d = lag(v7,1,mean(v7)); v7_2d = lag(v7,2,mean(v7)); v8_2d = lag(v8',2,mean(v8));

% ============= Implement a Constant gain controller ============
%Define an F matrix as follows:
Z1=[r2 r1_1d r3_1d];         %for s-1 
Z2=[r3 r2_1d r4_1d q_ramp];  %for s
Z3=[r4 r3_1d r5_1d];         %for s+1 
Z4=[r5 r4_1d r6_1d];         %for s+2
Z5=[r6 r5_1d r7_1d q_ramp2]; %for s+3 
Z6=[r7 r6_1d r8_1d];         %for s+4

% truncate data  see whether this makes any difference its not in F matrix anyway
Z1 = Z1(3:end,1:end); Z2 = Z2(3:end,1:end); Z3 = Z3(3:end,1:end);
Z4 = Z4(3:end,1:end); Z5 = Z5(3:end,1:end); Z6 = Z6(3:end,1:end);  dtc=1;  

[z1, m1, s1]=prepz(Z1, [], 1, 10, 0, dtc);  % 1 minute sampling
[z2, m2, s2]=prepz(Z2, [], 1, 10, 0, dtc);  % 1 minute sampling
[z3, m3, s3]=prepz(Z3, [], 1, 10, 0, dtc);  % 1 minute sampling
[z4, m4, s4]=prepz(Z4, [], 1, 10, 0, dtc);  % 1 minute sampling
[z5, m5, s5]=prepz(Z5, [], 1, 10, 0, dtc);  % 1 minute sampling
[z6, m6, s6]=prepz(Z6, [], 1, 10, 0, dtc);  % 1 minute sampling

z1=Z1; z2=Z2; z3=Z3; z4=Z4; z5=Z5; z6=Z6;

%============== link s-1
nn=[1 1 1 1 1 0];
[th, stats]=riv(z1, nn, [1 0 0 0 0]);  % least squares estimation
RT2=stats(3); YIC=stats(2); [at, bt]=getpar(th);

%============== link s
nn=[1 1 1 1 1 1 1 0];
[th2, stats2]=riv(z2, nn, [1 0 0 0 0]);  % least squares estimation
RT2_2=stats2(3); YIC_2=stats2(2); [at_2, bt_2]=getpar(th2);

%=============== link s+1
nn=[1 1 1 1 1 0];
[th3, stats3]=riv(z3, nn, [1 0 0 0 0]);  % least squares estimation
RT2_3=stats3(3); YIC_3=stats3(2); [at_3, bt_3]=getpar(th3);

%=============== link s+2
nn=[1 1 1 1 1 0];
[th4, stats4]=riv(z4, nn, [1 0 0 0 0]);  % least squares estimation
RT2_4=stats(3); YIC_4=stats(2); [at_4, bt_4]=getpar(th4);

%=============== link s+3
nn=[1 1 1 1 1 1 1 0];
[th5, stats5]=riv(z5, nn, [1 0 0 0 0]);  % least squares estimation
RT2_5=stats(3); YIC_5=stats(2); [at_5, bt_5]=getpar(th5);

%=============== link s+4
nn=[1 1 1 1 1 0];
[th6, stats6]=riv(z6, nn, [1 0 0 0 0]);  % least squares estimation
RT2_6=stats(3); YIC_6=stats(2); [at_6, bt_6]=getpar(th6);

aa=[-at(2) -at_2(2) -at_3(2) -at_4(2) -at_5(2) -at_6(2)]
bb=[bt(1, 2) bt_2(1, 2) bt_3(1, 2) bt_4(1, 2) bt_5(1, 2) bt_6(1, 2)]
cc=[bt(2, 2) bt_2(2, 2) bt_3(2, 2) bt_4(2, 2) bt_5(2, 2) bt_6(2, 2)]
dd=[0 bt_2(3, 2) 0 0 bt_5(3, 2) 0]

% F=[aa(1)  0    cc(1)  0     0    0    0    0;...
%     0    aa(4)  0    cc(4)  0    0    0    0;...
%    bb(2)  0    aa(2)  0    cc(2) 0    0    0;...
%     0    bb(5)  0    aa(5)  0   cc(5) 0    0;...
%     0    cc(3) bb(3)  0    aa(3) 0    0    0;...
%     0     0     0    bb(6)  0    0   aa(6) 0;...
%   -bb(2)  0   -aa(2)  0   -cc(2) 0    1    0;...
%     0   -bb(5)  0   -aa(5)  0  -cc(5) 0    1;];
% 
% g=[0 0 dd(2)  0     0   0 -dd(2) 0;... 
%    0 0   0  dd(5)   0   0    0  -dd(5);]'; 

F = [aa(1) cc(1)   0    0     0     0     0   0;...
     bb(2) aa(2) cc(2)  0     0     0     0   0;...
       0   bb(3) aa(3) cc(3)  0     0     0   0;...
       0     0   bb(4) aa(4) cc(4)  0     0   0;...
       0     0     0   bb(5) aa(5) cc(5)  0   0;...
       0     0     0    0    bb(6) aa(6)  0   0;...
   -bb(2) -aa(2) -cc(2) 0     0     0     1   0;...
      0     0     0  -bb(5) -aa(5) -cc(5) 0   1;];

g=[0 dd(2) 0  0     0   0 -dd(2) 0;... 
   0  0    0  0 dd(5) 0   0 -dd(5);]'; 

Q=eye(size(F)); r=eye(2); Q(end,end)=0.01; Q(end-1,end-1)=0.01;  r(end,end)=10; r(1,1)=10;

k_llm = dlqri(F,g,Q,r); k_llm(1,end)=-k_llm(1,end); k_llm(2,end)=-k_llm(2,end);





% ------------------ CONTROL IMPLEMENTATION ---------------------------
e=zeros(size(r2)); e2=zeros(size(r2));

uc=50*ones(size(r2)); uc2=70*ones(size(r2));

h=waitbar(0,'Please wait.... Control in Progress');

for n=3:length(r2);
    
 %link 2 [s-2]
    r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
    v1(n)=ftdcurve(p1,r1(n),1);
    
 %link 2 [s-1]
    r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
    v2(n)=ftdcurve(p2,r2(n),1);
    
 %link 3 [s]    ------- JUNCTION 1
    r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
    v3(n)=ftdcurve(p3,r3(n),1);
    
 %link 4 [s+1]
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
    v4(n)=ftdcurve(p4,r4(n),1);
    
 %link 5 [s+2]   
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    v5(n)=ftdcurve(p5,r5(n),1);
    
 %link 6 [s+3]  ------- JUNCTION 2
    r6(n)=r6(n-1)+((q_ramp2(n-1)+q5(n-1)-q6(n-1))*k6);
    q6(n)=alpha6*min(v6(n-1),v7(n-1))*r6(n-1);
    v6(n)=ftdcurve(p6,r6(n),1);

 %link 7 [s+4]
    r7(n)=r7(n-1)+((q6(n-1)-q7(n-1))*k7);
    q7(n)=alpha7*min(v7(n-1),v8(n-1))*r7(n-1);
    v7(n)=ftdcurve(p7,r7(n),1);

 %link 8 [s+5]
    r8(n)=r8(n-1)+((q7(n-1)-q8(n-1))*k8);
    q8(n)=alpha8*min(v8(n-1),v9(n-1))*r8(n-1);
    v8(n)=ftdcurve(p8,r8(n),1); 
    
    if  (n>=250) & (n<1400) 
        
     e(n) = uc(n) - r3(n); 
     
     e2(n)= uc2(n) - r6(n);

q_ramp(n) =q_ramp(n-1)  + k_llm(1,8)*e2(n) + k_llm(1,7)*e(n) - k_llm(1,1)*(r2(n)-r2(n-1)) - k_llm(1,2)*(r3(n)-r3(n-1)) - ...
    k_llm(1,3)*(r4(n)-r4(n-1)) - k_llm(1,4)*(r5(n)-r5(n-1)) - k_llm(1,5)*(r6(n)-r6(n-1)) - k_llm(1,6)*(r7(n)-r7(n-1));
 
q_ramp2(n)=q_ramp2(n-1) + k_llm(2,8)*e2(n) + k_llm(2,7)*e(n) - k_llm(2,1)*(r2(n)-r2(n-1)) - k_llm(2,2)*(r3(n)-r3(n-1)) - ...
    k_llm(2,3)*(r4(n)-r4(n-1)) - k_llm(2,4)*(r5(n)-r5(n-1)) - k_llm(2,5)*(r6(n)-r6(n-1)) - k_llm(2,6)*(r7(n)-r7(n-1));
          
     if q_ramp(n)<0;  q_ramp(n)=0;  elseif q_ramp(n)>45; q_ramp(n)=45; end
   
     if q_ramp2(n)<0; q_ramp2(n)=0; elseif q_ramp2(n)>45; q_ramp2(n)=45; end  
     
    end
    
    waitbar(n/length(q1),h);
end

close(h)

IAE  = abs(sum(uc(250:1400)-r3(250:1400)))
IAE2 = abs(sum(uc2(250:1400)-r6(250:1400)))
J_initial=IAE+IAE2

figure(2); 
subplot(211); plot(r3,'g'); plot(uc,'r');  plot(q_ramp,'m.');axis([0 1500 0 70]);
subplot(212); plot(r6,'g'); plot(uc2,'r'); plot(q_ramp2,'m.');axis([0 1500 0 120]);


break
%Optimisation !!!

% Initilisation of Q matrices
par1=1;    par2=0;     par3=0;     par4=0;     par5=0;     par6=0;     par7=0;      par8=0;
par9=1;    par10=0;    par11=0;    par12=0;    par13=0;    par14=0;    par15=0;     par16=1;
par17=0;   par18=0;    par19=0;    par20=0;    par21=0;    par22=1;    par23=0;     par24=0;     
par25=0;   par26=0;    par27=1;    par28=0;    par29=0;    par30=0;    par31=1;     par32=0; 
par33=0;   par34=1;    par35=0;    par36=1;    %par37=1;    par38=0;    par39=0;     par40=0;
% par41=0;   par42=0;    par43=0;    par44=0;    par45=0;    par46=1;    par47=0;     par48=0;
% par49=0;   par50=0;    par51=0;    par52=0;    par53=0;    par54=0;    par55=0.01;  par56=0;
% par57=0;   par58=0;    par59=0;    par60=0;    par61=0;    par62=0;    par63=0;     par64=0.02;

% 
% Q_vx=[par1  par2  par3  par4  par5  par6  par7  par8;...
%         par2  par9 par10 par11 par12 par13 par14 par15;...     
%         par3 par10 par16 par17 par18 par19 par20 par21;...     
%         par4 par11 par17 par22 par23 par24 par25 par26;...     
%         par5 par12 par18 par23 par27 par28 par29 par30;...     
%         par6 par13 par19 par24 par28 par31 par32 par33;...     
%         par7 par14 par20 par25 par29 par32 par34 par35;...     
%         par8 par15 par21 par26 par30 par33 par35 par36;];

QQ=[par1;    par2;     par3;     par4;     par5;     par6;     par7;      par8;...
    par9;   par10;    par11;    par12;    par13;    par14;    par15;     par16;...
   par17;   par18;    par19;    par20;    par21;    par22;    par23;     par24;...   
   par25;   par26;    par27;    par28;    par29;    par30;    par31;     par32;...
   par33;   par34;    par35;    par36; ]; 

Q_vx_opt=LSQNONLIN('LLM_fun_2I20',QQ,zeros(size(QQ)),1000*ones(size(QQ)));


 






