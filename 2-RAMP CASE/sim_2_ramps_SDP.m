close all
clear all

global F g r

%--------------------------------------------------------------------
% STM simulation
k1  = 0.0282;    % k(s-2) 
k2  = 0.02556;   % k(s-1) 
k3  = 0.02854;   % k(s)  
k4  = 0.02582;   % k(s+1)
k5  = 0.02482;   % k(s+2)
k6  = 0.02582;   % k(s+3)
k7  = 0.02682;   % k(s+4)
k8  = 0.01682;   % k(s+5)

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

    waitbar(n/length(q1),h);
end
close(h)

fig('Slacked Plots'); plot([r2 r3 r4 r5 r6 r7 r8']); legend('r2','r3','r4','r5','r6','r7','r8');

figure(2); subplot(211);hold; plot(r3);title('Junction1'); subplot(212);hold; plot(r6);title('Junction2');

r1_1d = lag(r1,1,mean(r1)); r1_2d = lag(r1,2,mean(r1)); r2_1d = lag(r2,1,mean(r2)); r2_2d = lag(r2,2,mean(r2)); r3_1d = lag(r3,1,mean(r3)); r3_2d = lag(r3,2,mean(r3)); r4_1d = lag(r4,1,mean(r4)); r4_2d = lag(r4,2,mean(r4));
r5_1d = lag(r5,1,mean(r5)); r5_2d = lag(r5,2,mean(r5)); r6_1d = lag(r6,1,mean(r6)); r6_2d = lag(r6,2,mean(r6)); r7_1d = lag(r7,1,mean(r7)); r7_2d = lag(r7,2,mean(r7)); r8_1d = lag(r8',1,mean(r8));r8_2d = lag(r8',2,mean(r8));
v1_1d = lag(v1,1,mean(v1)); v1_2d = lag(v1,2,mean(v1)); v2_1d = lag(v2,1,mean(v2)); v2_2d = lag(v2,2,mean(v2)); v3_1d = lag(v3,1,mean(v3)); v3_2d = lag(v3,2,mean(v3)); v4_1d = lag(v4,1,mean(v1)); v4_2d = lag(v4,2,mean(v4));
v5_1d = lag(v5,1,mean(v5)); v5_2d = lag(v5,2,mean(v5)); v6_1d = lag(v3,1,mean(v6)); v6_2d = lag(v6,2,mean(v6)); v7_1d = lag(v7,1,mean(v7)); v7_2d = lag(v7,2,mean(v7)); v8_2d = lag(v8',2,mean(v8));

% Regressors
z1=[r3_2d r2_2d r4_2d q_ramp];
z2=[r6_2d r5_2d r7_2d q_ramp2];

% Dependent States
inter_1=min(v3_2d,v4_2d); inter_2=min(v2_2d,v3_2d); inter_3=min(v4_2d,v5_2d);
inter_4=min(v6_2d,v7_2d); inter_5=min(v5_2d,v6_2d); inter_6=min(v7_2d,v8_2d);
x1=[inter_1 inter_2 inter_3 ones(size(q_ramp))]; x1=x1(3:end,:); z1=z1(3:end,:);
x2=[inter_4 inter_5 inter_6 ones(size(q_ramp2))]; x2=x2(3:end,:); z2=z2(3:end,:);

% SDP Analysis
TVP=0;nvr=[-1 -1 -1 0];
LHS=r3(3:end); LHS2=r6(3:end);
    [fit,fitse,par,parse,zs,pars,parses,rsq,nvre]=sdp(LHS,z1,x1,TVP,nvr);
figure(3); plot([ fit LHS ]);title(['Fit: ', num2str(rsq),]); legend('SDP estimated','Simulated Data' )
    [fit2,fitse2,par2,parse2,zs2,pars2,parses2,rsq2,nvre2]=sdp(LHS2,z2,x2,TVP,nvr);
figure(4); plot([ fit2 LHS2 ]);title(['Fit: ', num2str(rsq),]); legend('SDP estimated','Simulated Data' )



% ------------------ CONTROL IMPLEMENTATION ---------------------------
e=zeros(size(r2)); e2=zeros(size(r2)); z=zeros(size(r2)); z2=zeros(size(r2));
state_A=zeros(size(r2)); state_B=zeros(size(r2)); state_C=zeros(size(r2)); state_D=zeros(size(r2)); state_E=zeros(size(r2)); state_F=zeros(size(r2)); par_A=zeros(size(r2));
par_B=zeros(size(r2)); par_C=zeros(size(r2)); par_D=zeros(size(r2)); par_E=zeros(size(r2)); par_F=zeros(size(r2));

x1_irw1=irwsm(sort(x1(:,1)),1,1);
x1_irw2=irwsm(sort(x1(:,2)),1,1);
x1_irw3=irwsm(sort(x1(:,3)),1,1);
x2_irw1=irwsm(sort(x2(:,1)),1,1);
x2_irw2=irwsm(sort(x2(:,2)),1,1);
x2_irw3=irwsm(sort(x2(:,3)),1,1);



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
    
 if (n>=250) & (n<1400) 

     state_A(n-2)=min(v3(n-2),v4(n-2));if state_A(n-2)<=min(x1(:,1));state_A(n-2)=min(x1(:,1));elseif state_A(n-2)>max(x1(:,1)); state_A(n-2)=max(x1(:,1));end 
     state_B(n-2)=min(v2(n-2),v3(n-2));if state_B(n-2)<=min(x1(:,2));state_B(n-2)=min(x1(:,2));elseif state_B(n-2)>max(x1(:,2)); state_B(n-2)=max(x1(:,2));end
     state_C(n-2)=min(v4(n-2),v5(n-2));if state_C(n-2)<=min(x1(:,3));state_C(n-2)=min(x1(:,3));elseif state_C(n-2)>max(x1(:,3)); state_C(n-2)=max(x1(:,3));end    
     state_D(n-2)=min(v6(n-2),v7(n-2));if state_D(n-2)<=min(x2(:,1));state_D(n-2)=min(x2(:,1));elseif state_D(n-2)>max(x2(:,1)); state_D(n-2)=max(x2(:,1));end 
     state_E(n-2)=min(v5(n-2),v6(n-2));if state_E(n-2)<=min(x2(:,2));state_E(n-2)=min(x2(:,2));elseif state_E(n-2)>max(x2(:,2)); state_E(n-2)=max(x2(:,2));end
     state_F(n-2)=min(v7(n-2),v8(n-2));if state_F(n-2)<=min(x2(:,3));state_F(n-2)=min(x2(:,3));elseif state_F(n-2)>max(x2(:,3)); state_F(n-2)=max(x2(:,3));end
    
     par_A(n) = interp1(x1_irw1,pars(:,1), state_A(n-2));
     par_B(n) = interp1(x1_irw2,pars(:,2), state_B(n-2));
     par_C(n) = interp1(x1_irw3,pars(:,3), state_C(n-2));
     par_D(n) = interp1(x2_irw1,pars2(:,1),state_D(n-2));
     par_E(n) = interp1(x2_irw2,pars2(:,2),state_E(n-2));
     par_F(n) = interp1(x2_irw3,pars2(:,3),state_F(n-2));

   
 F = [par_A(n)    0     par_B(n)   0      par_C(n)   0     0   0;...
         0     par_D(n)    0    par_E(n)     0    par_F(n) 0   0;...
         0        0        0       0         0       0     0   0;...
         0        0        0       0         0       0     0   0;...
         0        0        0       0         0       0     0   0;...
         0        0        0       0         0       0     0   0;...
     -par_A(n)    0    -par_B(n)   0     -par_C(n)   0     1   0;...
         0    -par_D(n)    0   -par_E(n)     0   -par_F(n) 0   1]; 
 
 g = [par(1,4)  0      0  0  0  0  -par(1,4)  0;...
      0     par2(1,4)  0  0  0  0   0     -par2(1,4);]';

     Q = eye(size(F)); r=eye(2); 
         
     Q(end,end)=0.01; Q(end-1,end-1)=0.01;  r(end,end)=10; r(1,1)=10;
     
     k=dlqri(F,g,Q,r); 
        
     k_1(n)=k(1,1); k_2(n) = k(1,3);  k_3(n) = k(1,5);  k_4(n) = -k(1,7); 
     k_5(n)=k(2,2); k_6(n) = k(2,4);  k_7(n) = k(2,6);  k_8(n) = -k(2,8); 
   
     e(n) = uc(n) - r3(n); e2(n) = uc2(n) - r6(n);   
     
     %z(n) = z(n-1) + e(n); z2(n) = z2(n-1) + e2(n); 
        
     %q_ramp(n)  = k_4(n)*z(n)  - k_1(n)*r3(n) - k_2(n)*r2(n-1) - k_3(n)*r4(n-1); 
    
     %q_ramp2(n) = k_8(n)*z2(n) - k_5(n)*r6(n) - k_6(n)*r5(n-1) - k_7(n)*r7(n-1); 
   
     q_ramp(n) = q_ramp(n-1) + k_4(n)*e(n) - k_1(n)*(r3(n)-r3(n-1)) - k_2(n)*(r2(n)-r2(n-1)) - k_3(n)*(r4(n)-r4(n-1));
     
     q_ramp2(n) = q_ramp2(n-1) + k_8(n)*e2(n) - k_5(n)*(r6(n)-r6(n-1)) - k_6(n)*(r5(n)-r5(n-1)) - k_7(n)*(r7(n)-r7(n-1));

     if q_ramp(n)<0;  q_ramp(n)=0;  elseif q_ramp(n)>45; q_ramp(n)=45; end
   
     if q_ramp2(n)<0; q_ramp2(n)=0; elseif q_ramp2(n)>45; q_ramp2(n)=45; end  

   end
    
    waitbar(n/length(q1),h);
end

 close(h)

IAE  = abs(sum(uc(250:1400)-r3(250:1400)))
IAE2 = abs(sum(uc2(250:1400)-r6(250:1400)))
J_initial=IAE+IAE2

posit=r3(250:1400)-uc(250:1400);
index=find(posit>0); sum(posit(index))
posit2=r6(250:1400)-uc2(250:1400);
index2=find(posit2>0); sum(posit2(index2))

figure(2); 
subplot(211); plot(r3,'k'); plot(uc,'r');  plot(q_ramp,'m.');axis([0 1500 0 70]);
subplot(212); plot(r6,'k'); plot(uc2,'r'); plot(q_ramp2,'m.');axis([0 1500 0 120]);
  