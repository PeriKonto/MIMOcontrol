function J=LLM_fun_2I20(QQ)

global F g r


Q = [QQ(1)  QQ(2)  QQ(3)  QQ(4)   QQ(5) QQ(6)  QQ(7)  QQ(8);...
     QQ(2)  QQ(9) QQ(10) QQ(11) QQ(12) QQ(13) QQ(14) QQ(15);...     
     QQ(3) QQ(10) QQ(16) QQ(17) QQ(18) QQ(19) QQ(20) QQ(21);...     
     QQ(4) QQ(11) QQ(17) QQ(22) QQ(23) QQ(24) QQ(25) QQ(26);...     
     QQ(5) QQ(12) QQ(18) QQ(23) QQ(27) QQ(28) QQ(29) QQ(30);...     
     QQ(6) QQ(13) QQ(19) QQ(24) QQ(28) QQ(31) QQ(32) QQ(33);...     
     QQ(7) QQ(14) QQ(20) QQ(25) QQ(29) QQ(32) QQ(34) QQ(35);...     
     QQ(8) QQ(15) QQ(21) QQ(26) QQ(30) QQ(33) QQ(35) QQ(36);];
  
% INDEX=find(Q<0); Q(INDEX)=0;

k_llm = dlqri(F,g,Q,r); 

k_llm(1,end)=-k_llm(1,end); k_llm(2,end)=-k_llm(2,end);

% STM simulation
k1=0.0282; k2=0.02556; k3=0.02854; k4=0.02582; k5=0.02482; k6=0.02582; k7=0.02682; k8=0.01682; 
alpha1=0.0167; alpha2=alpha1; alpha3=alpha1; alpha4=alpha1; alpha5=alpha1; alpha6=alpha1; alpha7=alpha1; alpha8=alpha1;
p1=[7 100 110 63 360]; p2=p1; p3=p1; p4=p1; p5=p1; p6=p1; p7=p1; p8=p1; p9=p1;

% Set boundary conditions
kappa=linspace(3,70,50); kappa2=linspace(60,10,40); kappa3=linspace(70,60,50); q_in=[3*ones(400,1); kappa'; 70*ones(350,1); kappa3'; 60*ones(250,1); kappa2'; 10*ones(300,1)];  
v_in=[55*ones(1440,1);]; v9=v_in; r9=ftdcurve(p9,v9,2); 

% initialise other variables
r1=zeros(size(q_in)); r2=zeros(size(q_in)); r3=zeros(size(q_in)); r4=zeros(size(q_in)); r5=zeros(size(q_in)); r6=zeros(size(q_in)); r7=zeros(size(q_in)); r8=zeros(size(q_in)); v1=zeros(size(q_in)); v2=zeros(size(q_in)); v3=zeros(size(q_in)); v4=zeros(size(q_in)); v5=zeros(size(q_in));
v6=zeros(size(q_in)); v7=zeros(size(q_in)); v8=zeros(size(q_in));
q1=zeros(size(q_in)); q2=zeros(size(q1)); q3=zeros(size(q1)); q4=zeros(size(q1)); q5=zeros(size(q1)); q6=zeros(size(q1));
q7=zeros(size(q1));   q8=zeros(size(q1));
r1(1)=0.1; r2(1)=0.1; r3(1)=0.1; r4(1)=0.1; r5(1)=0.1; r6(1)=0.1; r7(1)=0.1; r8(1)=0.1; r9(1)=0.1;
q1(1)=0.1; q2(1)=0.1; q3(1)=0.1; q4(1)=0.1; q5(1)=0.1; q6(1)=0.1; q7(1)=0.1; q8(1)=0.1; q9(1)=0.1; 
v1(1)=ftdcurve(p1,r1(1),1); v2(1)=ftdcurve(p2,r2(1),1); v3(1)=ftdcurve(p3,r3(1),1); v4(1)=ftdcurve(p4,r4(1),1); v5(1)=ftdcurve(p5,r5(1),1); v6(1)=ftdcurve(p6,r6(1),1); v7(1)=ftdcurve(p7,r7(1),1); v8(1)=ftdcurve(p8,r8(1),1); v9(1)=ftdcurve(p9,r9(1),1);
q_ramp=13*ones(size(q1)); q_ramp2=13*ones(size(q1));

% ------------------ CONTROL IMPLEMENTATION ---------------------------
e=zeros(size(r2)); e2=zeros(size(r2)); uc=50*ones(size(r2)); uc2=70*ones(size(r2));

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

IAE  = abs(sum(uc(250:1400)-r3(250:1400)));
IAE2 = abs(sum(uc2(250:1400)-r6(250:1400)));

J=IAE+IAE2

if ~finite(J)
  J=1e+23;
end 



