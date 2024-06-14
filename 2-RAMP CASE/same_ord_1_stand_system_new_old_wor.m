close all; clear all; %clc

% STM simulation
k1  = 0.031;  % k(s-2) 
k2  = 0.022;  % k(s-1) 
k3  = 0.033;  % k(s)  
k4  = 0.024;  % k(s+1)
k5  = 0.025;  % k(s+2)

alpha1=0.0167; alpha2=alpha1; alpha3=alpha1; alpha4=alpha1; alpha5=alpha1; 
p1=[7 100 120  33.5  1100]; 
p2=[7 100 120  35.1  1000]; 
p3=[7 95  100  53.8  1100]; 
p4=[7 95  120  60    3600]; 
p5=[7 110 120  63    5600];
p6=[7 100 120  63    3600];

%v_in=[30*ones(1440,1); ];  v6=v_in; r6=ftdcurve(p6,v6,2);
r6  =[10*ones(300,1);  140*ones(500,1);  10*ones(300,1)];
p6=[5 100 110 63 1000]; v6=ftdcurve(p6,r6,1);

q_in=[10*ones(150,1); linspace(10,70,150)'; 70*ones(500,1); linspace(70,10,150)'; 10*ones(150,1)]; 
q_in=0.80*q_in;  

%Ramp Input
q_ramp=[3*ones(200,1); linspace(3,30,200)'; 30*ones(400,1); linspace(30,3,150)'; 3*ones(150,1);];
% q_ramp=zeros(size(q_ramp));
q_ramp_load=sum(q_ramp(250:1050));
% initialise other variables
r1=NaN*ones(size(q_in)); r2=NaN*ones(size(q_in)); r3=NaN*ones(size(q_in)); r4=NaN*ones(size(q_in)); r5=NaN*ones(size(q_in)); 
v1=NaN*ones(size(q_in)); v2=NaN*ones(size(q_in)); v3=NaN*ones(size(q_in)); v4=NaN*ones(size(q_in)); v5=NaN*ones(size(q_in)); 
q1=NaN*ones(size(q_in)); q2=NaN*ones(size(q_in)); q3=NaN*ones(size(q_in)); q4=NaN*ones(size(q_in)); q5=NaN*ones(size(q_in));
r1(1:2)=2; r2(1:2)=2; r3(1:2)=2; r4(1:2)=2; r5(1:2)=2; q1(1:2)=2; q2(1:2)=2; q3(1:2)=2; q4(1:2)=2; q5(1:2)=2; 
v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1); v5(1:2)=ftdcurve(p5,r5(1:2),1);

dt=1;
% === Open Loop Simulation ===
h=waitbar(0,'Open Loop Simulation ');

for n=3:length(q_in);
    
    %link 1 [s-2]
    r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1)*dt;
    v1(n)=ftdcurve(p1,r1(n),1); 
   
    %link 2 [s-1]
    r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1)*dt;
    v2(n)=ftdcurve(p2,r2(n),1);
    B(n)=k3*alpha2*min(v2(n-2),v3(n-2));
    
    %link 3 [s]
    r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1)*dt;
    v3(n)=ftdcurve(p3,r3(n),1);
    A(n)=-(k3-k4)*alpha3*min(v3(n-2),v4(n-2));
    
    %link 4 [s+1]
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1)*dt;
    v4(n)=ftdcurve(p4,r4(n),1);
    C(n)=-k4*alpha4*min(v4(n-2),v5(n-2));
    
    %link 5 [s+2]
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1)*dt;
    v5(n)=ftdcurve(p5,r5(n),1);
    
    r3sdp(n)=- r4(n) + r3(n-1) + r4(n-1) + A(n)*r3(n-2) + B(n)*r2(n-2) + C(n)*r4(n-2);
    r3sdp_ramp(n)=- r4(n) + r3(n-1) + r4(n-1) + A(n)*r3(n-2) + B(n)*r2(n-2) + C(n)*r4(n-2) + k3*q_ramp(n-1);
    r3sdp_red(n)= r3(n-1) + A(n)*r3(n-2) + B(n)*r2(n-2) + C(n)*r4(n-2);
    
    waitbar(n/length(r2),h);
end
close(h)


figure(12);
subplot(321); plot(r1,q1,'.'); hold; plot(p1(4),linspace(0,120),'k.'); plot(p1(4),max(q1),'r.');title('s-2');
subplot(322); plot(r2,q2,'.'); hold; plot(p2(4),linspace(0,120),'k.'); plot(p2(4),max(q2),'r.');title('s-1 {load}');
subplot(323); plot(r3,q3,'.'); hold; plot(p3(4),linspace(0,120),'k.'); plot(p3(4),max(q3),'r.');title('s {junction}');
subplot(324); plot(r4,q4,'.'); hold; plot(p4(4),linspace(0,120),'k.'); plot(p4(4),max(q4),'r.');title('s+1');
subplot(325); plot(r5,q5,'.'); hold; plot(p5(4),linspace(0,120),'k.'); plot(p5(4),max(q5),'r.');title('s+2');

figure(11);plot([r3 r2]); hold; ylabel('density @ junction r(s,t) [veh/km]'); xlabel('Simulation Period [min]'); plot(r3sdp_ramp,'m');

UnmLoad=r2; flow_ups=q2; flow_jun=q3; old_v3=v3; ramp_load=sum(q_ramp(250:900));


figure(13);
x=[0:1:500]';y_1=ftdcurve(p1,x,1);subplot(231);plot(x,y_1); hold; axis([0 500 0 150]);title('s-2');
x=[0:1:500]';y_2=ftdcurve(p2,x,1);subplot(232);plot(x,y_2); hold; axis([0 500 0 150]) ;title('s-1 {load}');
x=[0:1:500]';y_3=ftdcurve(p3,x,1);subplot(233);plot(x,y_3); hold; axis([0 500 0 150]) ;title('s {junction}');
x=[0:1:500]';y_4=ftdcurve(p4,x,1);subplot(234);plot(x,y_4); hold; axis([0 500 0 150]);title('s+1'); 
x=[0:1:500]';y_5=ftdcurve(p5,x,1);subplot(235);plot(x,y_5); hold; axis([0 500 0 150]);title('s+2');
% x=[0:1:500]';y_6=ftdcurve(pramp,x,1);subplot(236);plot(x,y_6); hold; axis([0 500 0 150]);title('ramp');
subplot(231);plot(r1,v1,'.'); plot(p1(4),p1(2),'r.');title('s-2')
subplot(232);plot(r2,v2,'.'); plot(p2(4),p2(2),'r.');title('s-1')
subplot(233);plot(r3,v3,'.'); plot(p3(4),p3(2),'r.');title('s')
subplot(234);plot(r4,v4,'.'); plot(p4(4),p4(2),'r.');title('s+1')
subplot(235);plot(r5,v5,'.'); plot(p5(4),p5(2),'r.');title('s+2')
% subplot(236);plot(r_ramp,v_ramp,'.'); plot(pramp(4),pramp(2),'r.');title('ramp')

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
% [z1, m1, s1]=prepz(Z1, [], 1, 10, 0, dtc);  % 1 minute sampling
% [z2, m2, s2]=prepz(Z2, [], 1, 10, 0, dtc);  % 1 minute sampling
% [z3, m3, s3]=prepz(Z3, [], 1, 10, 0, dtc);  % 1 minute sampling
z1=Z1; z2=Z2; z3=Z3;

%============== link s-1
nn=[1 1 1 1 1 0];
[th, stats]=riv(z1, nn);
RT2=stats(3); YIC=stats(2); [at, bt]=getpar(th);
RT2
YIC

%============== link s
nn= [1 1 1 1 1 1 1 0];
[th2, stats2]=riv(z2, nn);
RT2_2=stats2(3); YIC_2=stats2(2); [at_2, bt_2]=getpar(th2);
RT2_2
YIC_2

%=============== link s+1
nn=[1 1 1 0];
[th3, stats3]=riv(z3, nn);
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

Q=eye(size(F)); Q(end,end)=1; r=0.01; k_llm = dlqr(F,g,Q,r); k_llm(end)=-k_llm(end);

[th,stats]=rivid(z1,[1 1 1 1 1 0; 2 1 1 2 2 0]);                 [at_1, bt_1]=getpar(th);
[th2,stats2]=rivid(z2,[1 1 1 1 1 1 1 0; 2 1 1 1 2 2 1 0]); [at_2, bt_2]=getpar(th2);
[th3,stats3]=rivid(z3,[1 1 1 0; 2 1 2 0]);                 [at_3, bt_3]=getpar(th3);       


%location s+1
% example fixing parameter values
%fix parameters at a1 =-1
nn=[2 1 2 0];[th3, stats3]=riv(z3, nn);  RT2_3_ini=stats3(3); YIC_3_ini=stats3(2); [at_3_ini, bt_3_ini]=getpar(th3);
a0=[-1 0 0 ]; P0=[0 1000 1000]; flags=[];
[th3, stats3]=riv(z3, nn,flags,a0,P0);
RT2_3_fin=stats3(3); YIC_3_fin=stats3(2); [at_3_fin, bt_3_fin]=getpar(th3);

%location s-1
%nn=[2 1 1 2 2 0]; [th, stats]=riv(z1, nn); RT2=stats(3); YIC=stats(2); [at_ini, bt_in]=getpar(th);



% %--------------------- TVP----------------------------------

% %============== link s-1
% nn=[1 1 1 1 1];
% nvr=dtfmopt(z1(:,1),z1(:,2:3),nn,[]); 
% [tfs,fit,fitse,par_1_tvp,parse,e]=dtfm(z1(:,1),z1(:,2:3),nn,[],nvr);
% 
% %============== link s
% nn=[1 1 1 1 1 1 1];
% nvr=dtfmopt(z2(:,1),z2(:,2:4),nn,[]); 
% [tfs2,fit2,fitse2,par_2_tvp,parse2,e2]=dtfm(z2(:,1),z2(:,2:4),nn,[],[nvr(1) nvr(2) nvr(3) 0]);
% 
% %============== link s+1
% nn=[1 1 1];
% nvr=dtfmopt(z3(:,1),z3(:,2),nn,[]); 
% [tfs3,fit3,fitse3,par_3_tvp,parse3,e3]=dtfm(z3(:,1),z3(:,2),nn,[],nvr);
% 
% save tvp_pars  par_1_tvp par_2_tvp par_3_tvp

%-----------------------------------------------------------------


%-----------------------------  SDP Modelling--------------------------
z=[r3_1d r2_1d r4_1d q_ramp]; 
inter_1=max(r3_2d,r4_2d); 
inter_2=max(r2_2d,r3_2d);
inter_3=max(r4_2d,r5_2d);
LHS = r3(3:end); TVP=[0 0 0 0]; nvr=[-1 -1 -1 0]; 
x=[inter_1 inter_2 inter_3 zeros(size(inter_1))];  x=x(3:end,:); z=z(3:end,:);

[fit,fitse,par,parse,zs,pars,parses,rsq,nvre]=sdp(LHS,z,x,TVP,nvr); figure(3); plot([fit LHS]);title(['Fit: ', num2str(rsq),]); legend('SDP estimated','Simulated Data' )
[x1_sort I1]=sort(x(:,1)); A_i=A(I1); [x2_sort I2]=sort(x(:,2)); B_i=B(I1); [x3_sort I3]=sort(x(:,3)); C_i=C(I1); 

figure(4);
subplot(311); plot(x1_sort, pars(:,1),'.' , x1_sort, pars(:,1) + parses(:,1),'k-.' , x1_sort,  pars(:,1)-parses(:,1),'k-.');axis([0 90 0.95 0.97]);
subplot(312); plot(x2_sort, pars(:,2),'.' , x2_sort, pars(:,2) + parses(:,2),'k-.' , x2_sort,  pars(:,2)-parses(:,2),'k-.');axis([0 90 0.032 0.033]);
subplot(313); plot(x3_sort, pars(:,3),'.' , x3_sort, pars(:,3) + parses(:,3),'k-.' , x3_sort,  pars(:,3)-parses(:,3),'k-.');axis([0 120 -0.001 0.01]);

% figure(5);
% subplot(311); plot(x1_sort,pars(:,1),'.' , x1_sort,pars(:,1) + parses(:,1),'k-.' , x1_sort, pars(:,1)-parses(:,1),'k-.');
% subplot(312); plot(x2_sort, pars(:,2),'.', x2_sort, pars(:,2)+ parses(:,2),'k-.' , x2_sort, pars(:,2)-parses(:,2),'k-.');
% subplot(313); plot(x3_sort,pars(:,3),'.' , x3_sort,pars(:,3) + parses(:,3),'k-.' , x3_sort, pars(:,3)-parses(:,3),'k-.');


% % Parameterisation of relationships 
% figure(4)
% c =polyfit(x1_sort,  pars(:,1),4); Y1=polyval(c,x1_sort);   subplot(311); hold; plot(x1_sort,Y1,'g'); 
% c2=polyfit(x2_sort,  pars(:,2),4); Y2=polyval(c2,x2_sort);  subplot(312); hold; plot(x2_sort,Y2,'g');  
% c3=polyfit(x3_sort,  pars(:,3),4); Y3=polyval(c3,x3_sort);  subplot(313); hold; plot(x3_sort,Y3,'g'); 
% 
% %anfis
% epoch_n = 20;
% in_fis  = genfis1([x1_sort  pars(:,1)],5,'gaussmf');   out_fis1 = anfis([x1_sort  pars(:,1)],in_fis,epoch_n);
% in_fis  = genfis1([x2_sort  pars(:,2)],5, 'gaussmf');  out_fis2 = anfis([x2_sort  pars(:,2)],in_fis,epoch_n);
% in_fis  = genfis1([x3_sort  pars(:,3)],5,'gaussmf');   out_fis3 = anfis([x3_sort  pars(:,3)],in_fis,epoch_n);
% subplot(311); plot(x1_sort,evalfis(x1_sort, out_fis1),'r');  legend('SDP Estimate', 'Error 1', 'Error 2', 'Poly','Anfis'); 
% subplot(312); plot(x2_sort,evalfis(x2_sort, out_fis2),'r');  legend('SDP Estimate', 'Error 1', 'Error 2', 'Poly','Anfis'); 
% subplot(313); plot(x3_sort,evalfis(x3_sort, out_fis3),'r');  legend('SDP Estimate', 'Error 1', 'Error 2', 'Poly','Anfis'); 
% 
% % PROPER USAGE OF ANFIS
% % ====first graph========================= 
% %mftypes {1:dsigmf; 2:gauss2mf, 3:gaussmf, 4:gbellmf, 5:pimf, 6:psigmf, 7:trapmf, 8:sigmf, 9:trimf, 10:smf} 
% numMFs=[2:14]; MFtype=3; iout=1;
% [x0, Y0, Yhat, modsel, AIC, BIC, Rt2, SSer, Vred, intmftype, FS]=ANFIS1_id(x1_sort,pars(:,1),[],numMFs,MFtype,iout); 
% [in,out,rules]=fismat2InOut(FS{modsel},numMFs(modsel),1);%how to extract the parameters
% vx = FIS2Vec(in,out); %how to put everything in a vector
% recons_Yhat=FISnon_SISO(x1_sort,1,numMFs(modsel),FS{modsel},vx,rules,{'gaussmf'});
% OPTIONS = OPTIMSET('Display','on','MaxIter',100); 
% VX=fminsearch('myfunction',vx, OPTIONS, []  ,[],  x1_sort, pars(:,1), rules, MFtype, numMFs(modsel),1,[]);
% optmised_Yhat=FISnon_SISO(x1_sort,1,numMFs(modsel),[],VX,rules,{'gaussmf'});
% 
% % ====second graph========================= 
% numMFs=[2:14]; MFtype=3; iout=1;
% [x02, Y02, Yhat2, modsel2, AIC, BIC, Rt2, SSer, Vred, intmftype, FS]=ANFIS1_id(x2_sort,pars(:,2),[],numMFs,MFtype,iout); 
% [in,out,rules2]=fismat2InOut(FS{modsel2},numMFs(modsel2),1); vx2 = FIS2Vec(in,out);
% recons_Yhat2=FISnon_SISO(x2_sort,1,numMFs(modsel2),FS{modsel2},vx2,rules2,{'gaussmf'});
% VX2=fminsearch('myfunction',vx2, OPTIONS, []  ,[],  x2_sort, pars(:,2), rules2, MFtype, numMFs(modsel2),1,[]);
% optmised_Yhat2=FISnon_SISO(x2_sort,1,numMFs(modsel2),[],VX2,rules2,{'gaussmf'});
% 
% % =====third graph===========================
% numMFs=[2:14]; MFtype=3; iout=1;
% [x03, Y03, Yhat3, modsel3, AIC, BIC, Rt2, SSer, Vred, intmftype, FS]=ANFIS1_id(x3_sort,pars(:,3),[],numMFs,MFtype,iout); 
% [in,out,rules3]=fismat2InOut(FS{modsel3},numMFs(modsel3),1); vx3 = FIS2Vec(in,out);
% recons_Yhat3=FISnon_SISO(x3_sort,1,numMFs(modsel3),FS{modsel3},vx3,rules3,{'gaussmf'});
% VX3=fminsearch('myfunction',vx3, OPTIONS, []  ,[],  x3_sort, pars(:,3), rules3, MFtype, numMFs(modsel3),1,[]);
% optmised_Yhat3=FISnon_SISO(x3_sort,1,numMFs(modsel3),[],VX3,rules3,{'gaussmf'});
% %--------------------------------------------------------------------------------------------------------------
% 
% figure(4);
% subplot(311);plot(x1_sort, optmised_Yhat,'k');   legend('SDP Estimate', 'Error 1', 'Error 2', 'Poly','Anfis', 'ANFIS opt'); 
% subplot(312);plot(x2_sort, optmised_Yhat2,'k'); % legend('SDP Estimate', 'Error 1', 'Error 2', 'Poly','Anfis', 'ANFIS opt'); 
% subplot(313);plot(x3_sort, optmised_Yhat3,'k'); % legend('SDP Estimate', 'Error 1', 'Error 2', 'Poly','Anfis', 'ANFIS opt'); 


% % OPEN LOOP SIMULATION
% par_A=zeros(size(r2)); par_B=par_A; par_C=par_A; par_Ax=par_A; par_Bx=par_A; par_Cx=par_A;
% h=waitbar(0,'Open Loop Simulation Validation of F structure');
% r1(1:2)=2; r2(1:2)=2; r3(1:2)=2; r4(1:2)=2; r5(1:2)=2; q1(1:2)=2; q2(1:2)=2; q3(1:2)=2; q4(1:2)=2; q5(1:2)=2; 
% v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1); v5(1:2)=ftdcurve(p5,r5(1:2),1);
% 
% for n=3:length(r2);
%     
%     %------------------------------------------------------------------------------------------------------------------------------------   
%     %sdp model   
%     r3sdp_ramp(n)=- r4(n) + r3(n-1) + r4(n-1) + A(n)*r3(n-2) + B(n)*r2(n-2) + C(n)*r4(n-2) + k3*q_ramp(n-1);
%     %------------------------------------------------------------------------------------------------------------------------------------   
%     %system 
%     %link 1 [s-2]
%     r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
%     q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
%     v1(n)=ftdcurve(p1,r1(n),1); 
%     %link 2 [s-1]
%     r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
%     q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
%     v2(n)=ftdcurve(p2,r2(n),1);
%     %link 3 [s]
%     r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
%     q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
%     v3(n)=ftdcurve(p3,r3(n),1);
%     %link 4 [s+1]
%     r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
%     q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
%     v4(n)=ftdcurve(p4,r4(n),1);
%     %link 5 [s+2]
%     r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
%     q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
%     v5(n)=ftdcurve(p5,r5(n),1);
%     
%     %-----------------------------------------------------------------------------------------------------------------------------------
%     state_A(n-2) = max(r3(n-2),r4(n-2));    % if state_A(n-2)<=min(x(:,1));state_A(n-2)=min(x(:,1));elseif state_A(n-2)>max(x(:,1));state_A(n-2)=max(x(:,1));end 
%     state_B(n-2) = max(r2(n-2),r3(n-2));    % if state_B(n-2)<=min(x(:,2));state_B(n-2)=min(x(:,2));elseif state_B(n-2)>max(x(:,2));state_B(n-2)=max(x(:,2));end
%     state_C(n-2) = max(r4(n-2),r5(n-2));    % if state_C(n-2)<=min(x(:,3));state_C(n-2)=min(x(:,3));elseif state_C(n-2)>max(x(:,3));state_C(n-2)=max(x(:,3));end
%     
%     %non-parametric   
%     par_A_h(n-2) =  interp1(x1_irw,  par(:,1), state_A(n-2));
%     par_B_h(n-2) =  interp1(x2_irw,  par(:,2), state_B(n-2));
%     par_C_h(n-2) =  interp1(x3_irw,  par(:,3), state_C(n-2));
%     
%     %parametric
%     %  par_A(n-2)=  FISnon_SISO(state_A(n-2), 1, numMFs(modsel),  [], VX,  rules,  {'gaussmf'});
%     %  par_B(n-2)=  FISnon_SISO(state_B(n-2), 1, numMFs(modsel2), [], VX2, rules2, {'gaussmf'});
%     %  par_C(n-2)=  FISnon_SISO(state_C(n-2), 1, numMFs(modsel3), [], VX3, rules3, {'gaussmf'});
%     
%     %using true nonlinearities 
% %     par_Ax(n-2) =  interp1(x1_sort,  A_i', state_A(n-2)); 
% %     par_Bx(n-2) =  interp1(x2_sort,  B_i', state_B(n-2));
% %     par_Cx(n-2) =  interp1(x3_sort,  C_i', state_C(n-2));
%     
%     FF(n)     =  par_A_h(n-2)*r3(n-1) +  par_B_h(n-2)*r2(n-1) +  par_C_h(n-2)*r4(n-1) +   pars(1,4)*q_ramp(n-1);
%     %  FF_par(n) =  -par_A(n-2)*r3(n-1)   +  par_B(n-2)*r2(n-1)   +  par_C(n-2)*r4(n-1)   +   pars(1,4)*q_ramp(n-1);
%     llm(n)    =  aa(2)*r3(n-1)        +  bb(2)*r2(n-1)        +  cc(2)*r4(n-1)        +   dd(2)*q_ramp(n-1);
%     %  Rx(n)     =  1*r3(n-1)-par_Ax(n-2)*r3(n-2)+par_Bx(n-2)*r2(n-2)+par_Cx(n-2)*r4(n-2)+k3*q_ramp(n-1); %Deterministic
%     r3sdp_ramp_red(n)= r3(n-1) + A(n)*r3(n-2) + B(n)*r2(n-2) + C(n)*r4(n-2) + k3*q_ramp(n-1); 
%     waitbar(n/length(r2),h);
% end
% close(h)
% disp('SDP non-par');
% sum(r3'-FF)
% % sum(r3'-FF_par)
% disp('LLM');
% sum(r3'-llm)
% %sum(r3'-Rx)
% disp('Deterministic');
% sum(r3'-r3sdp_ramp)
% sum(r3'-r3sdp_ramp_red)
% 
% figure(11); plot(FF,'r'); plot(llm,'g'); %plot(FF_par,'c')
% legend('r3','r2', 'SDP_open','LLM_open','SDP_param');
% 
% figure(5);
% subplot(311);hold; plot(state_A,  par_A_h,'go'); xlabel('Dependent State min(r(s,t),r(s-1,t-1))'); ylabel('Parameter a');
% subplot(312);hold; plot(state_B,  par_B_h,'go'); xlabel('Dependent State min(r(s-1,t),r(s,t-1))'); ylabel('Parameter b');
% subplot(313);hold; plot(state_C,  par_C_h,'go'); xlabel('Dependent State min(r(s+1,t),r(s+2,t-1))'); ylabel('Parameter c');
% 

x1_irw=irwsm(x1_sort,1,1);
x2_irw=irwsm(x2_sort,1,1);
x3_irw=irwsm(x3_sort,1,1);


%  %%%%%%%%%%%%%%%%%%%%%%%%%%% control implementations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(1) NON-PARAMETRIC
%initialization of variables
e=zeros(length(q1),1); z=e; state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1));
par_A_h=zeros(size(q1)); par_B_h=par_A_h; par_C_h=par_A_h; 
uc=50*ones(length(q1),1); 

h=waitbar(0,'Please wait.... Control in Progress (non-parametric)');

for n=3:length(r2);
    
    %link 1 [s-2]
    r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1)*dt;
    v1(n)=ftdcurve(p1,r1(n),1); 
    %link 2 [s-1]
    r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1)*dt;
    v2(n)=ftdcurve(p2,r2(n),1);
    B(n)=k3*alpha2*min(v2(n-2),v3(n-2));
    %link 3 [s]
    r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1)*dt;
    v3(n)=ftdcurve(p3,r3(n),1);
    A(n)=-(k3-k4)*alpha3*min(v3(n-2),v4(n-2));
    %link 4 [s+1]
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4); 
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1)*dt;
    v4(n)=ftdcurve(p4,r4(n),1);
    C(n)=-k4*alpha4*min(v4(n-2),v5(n-2));
    %link 5 [s+2]
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1)*dt;
    v5(n)=ftdcurve(p5,r5(n),1);
    
    if  (n>=250) & (n<1050) 
         
        % read-off States 
        state_A(n-2) = max(r3(n-2),r4(n-2)); if state_A(n-2) <= min(x(:,1)); state_A(n-2)=min(x(:,1)); elseif state_A(n-2)>max(x(:,1)); state_A(n-2)=max(x(:,1));end 
        state_B(n-2) = max(r2(n-2),r3(n-2)); if state_B(n-2) <= min(x(:,2)); state_B(n-2)=min(x(:,2)); elseif state_B(n-2)>max(x(:,2)); state_B(n-2)=max(x(:,2));end
        state_C(n-2) = max(r4(n-2),r5(n-2)); if state_C(n-2) <= min(x(:,3)); state_C(n-2)=min(x(:,3)); elseif state_C(n-2)>max(x(:,3)); state_C(n-2)=max(x(:,3));end
        
        % estimated non-parametric parameters 
        par_A_h(n-2) =  interp1(x1_irw,pars(:,1), state_A(n-2));
        par_B_h(n-2) =  interp1(x2_irw,pars(:,2), state_B(n-2));
        par_C_h(n-2) =  interp1(x3_irw,pars(:,3), state_C(n-2));
        
        %      F = [par_A_h(n-2) par_B_h(n-2) par_C_h(n-2) 0;
        %               0 0 0 0;
        %               0 0 0 0;
        %              -par_A_h(n-2) -par_ B_h(n-2) -par_C_h(n-2) 1];
        %         g =[pars(1,4) 0 0 -pars(1,4)];
        
        F = [ 0 0 0 0;
            par_B_h(n-2) par_A_h(n-2) par_C_h(n-2) 0;
            0 0 0 0;
            -par_B_h(n-2) -par_A_h(n-2) -par_C_h(n-2) 1];
        
        g =[0 pars(1,4) 0 -pars(1,4)];
        
        Q = eye(size(F)); Q(end,end)=1; r=0.01;  k=dlqri(F,g',Q,r); 
        
        k_1(n)=k(1,1); k_2(n) = k(1,2);  k_3(n) = k(1,3);  k_4(n) = -k(1,4); 
        
        e(n) = uc(n) - r3(n);  
        
        % z(n) = z(n-1) + e(n);  q_ramp(n) = k_4(n)*z(n) - k_1(n)*r3(n) - k_2(n)*r2(n) - k_3(n)*r4(n); 
        
        q_ramp(n)=q_ramp(n-1) + k_4(n)*e(n) - k_2(n)*(r3(n)-r3(n-1)) - k_1(n)*(r2(n)-r2(n-1)) - k_3(n)*(r4(n)-r4(n-1));  

        if q_ramp(n)<0; q_ramp(n)=0;  elseif q_ramp(n)>=30; q_ramp(n)=30; end  
        
     end       
        waitbar(n/length(q1),h);    
end
    close(h)

    figure(11); plot(uc,'r'); plot(r3,'g')
    disp('SDP-NON-parametric');
    IAE_sdp=sum(abs(r3(250:1050)-uc(250:1050)));
    resid=r3(250:1050)-uc(250:1050); index=find(resid>0); IAE2_sdp=sum(resid(index));
    Control_effort_sdp=sum(abs(q_ramp(250:1050)));
    %load1=sum((UnmLoad(250:1050)-r2(250:1050)))./sum(r2(250:1050))%;load1=load1*100 
    % max_flow=sum(q2(250:1050)-flow_ups(250:1050))./sum(flow_ups(250:1050))%; max_flow=max_flow*100
    % max_flow2=sum(q3(250:1050)-flow_jun(250:1050))./sum(flow_jun(250:1050))%;max_flow2=max_flow2*100
     chang_V3_sdp=sum(v3(250:1050)-old_v3(250:1050))./sum(old_v3(250:1050))
    
    figure(4);
    subplot(311);hold; plot(state_A,  par_A_h,'go'); xlabel('Dependent State min(r(s,t),r(s-1,t-1))');   ylabel('Parameter a');
    subplot(312);hold; plot(state_B,  par_B_h,'go'); xlabel('Dependent State min(r(s-1,t),r(s,t-1))');   ylabel('Parameter b');
    subplot(313);hold; plot(state_C,  par_C_h,'go'); xlabel('Dependent State min(r(s+1,t),r(s+2,t-1))'); ylabel('Parameter c');
    
    disp('-----Residual of intepolation----')
    disp('***Par_A***');
    sum(par(250:1040,1) - par_A_h(250:1040) )
       disp('Acceptable level +/-');
       abs(sum( par(250:1040,1)-par(250:1040,1)-parse(250:1040,1) ))
    disp('***Par_B***');   
    sum(par(250:1040,2) - par_B_h(250:1040) )
       disp('Acceptable level +/-');
       abs(sum( par(250:1040,2)-par(250:1040,2)-parse(250:1040,2) ))
    disp('***Par_C***');    
    sum(par(250:1040,3) - par_C_h(250:1040) )
       disp('Acceptable level +/-');
       abs(sum( par(250:1040,3)-par(250:1040,3)-parse(250:1040,3) ))
              
       
%   Delta_T=(max(flow_jun)-q3(250:1040)) ./ ( q2(250:1040)+q_ramp(250:1040)-q3(250:1040) ); 

    % %  %%%%%%%%%%%%%%%%%%%%%%%%%%% control implementations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %(2) NON-PARAMETRIC using values from LLM
    % %initialization of variables
    % par_A=zeros(size(q1)); par_B=par_A; par_C=par_A; par_D=par_A; e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); state_D=zeros(size(q1)); state_D2=zeros(size(q1));
    % par_A_h=zeros(size(q1)); par_B_h=par_A; par_C_h=par_A; par_D_h=par_A; par_D2_h=par_A;
    % uc=53*ones(length(q1),1);
    % %mean_pars1=mean(pars(:,1)); 
    % pars1=aa(2)*ones(size(pars(:,1)));
    % %mean_pars2=mean(pars(:,2));
    % pars2=bb(2)*ones(size(pars(:,2)));
    % %mean_pars3=mean(pars(:,3)); 
    % pars3=cc(2)*ones(size(pars(:,3)));
    % 
    % h=waitbar(0,'Please wait.... Control in Progress (non-parametric)');
    % 
    % for n=3:length(r2);
    %     
    %   %link 1 [s-2]
    %     r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    %     q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
    %     v1(n)=ftdcurve(p1,r1(n),1); 
    %   %link 2 [s-1]
    %     r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    %     q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
    %     v2(n)=ftdcurve(p2,r2(n),1);
    %   %link 3 [s]
    %     r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    %     q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
    %     v3(n)=ftdcurve(p3,r3(n),1);
    %   %link 4 [s+1]
    %     r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    %     q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
    %     v4(n)=ftdcurve(p4,r4(n),1);
    %   %link 5 [s+2]
    %     r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
    %     q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    %     v5(n)=ftdcurve(p5,r5(n),1);
    %     
    %     if  (n>=250) & (n<1050) 
    %         
    %         % read-off States 
    % %         state_A(n-2) = max(r3(n-2),r4(n-2)); %if state_A(n-2) <= min(x(:,1)); state_A(n-2)=min(x(:,1)); elseif state_A(n-2)>max(x(:,1)); state_A(n-2)=max(x(:,1));end 
    % %         state_B(n-2) = max(r2(n-2),r3(n-2)); %if state_B(n-2) <= min(x(:,2)); state_B(n-2)=min(x(:,2)); elseif state_B(n-2)>max(x(:,2)); state_B(n-2)=max(x(:,2));end
    % %         state_C(n-2) = max(r4(n-2),r5(n-2)); %if state_C(n-2) <= min(x(:,3)); state_C(n-2)=min(x(:,3)); elseif state_C(n-2)>max(x(:,3)); state_C(n-2)=max(x(:,3));end
    %         
    %         % estimated non-parametric parameters 
    %         par_A_h = pars1(n);
    %         par_B_h = pars2(n);
    %         par_C_h = pars3(n);
    %         
    %         F = [-par_A_h par_B_h par_C_h 0;
    %               0 0 0 0;
    %               0 0 0 0;
    %              par_A_h -par_B_h -par_C_h 1];
    % 
    % %         g =[pars(1,4) 0 0 -pars(1,4)];
    % 
    %         g =[dd(2) 0 0 -dd(2)];
    %  
    %         Q = eye(size(F)); Q(end,end)=0.1; r=1;  k=dlqri(F,g',Q,r); 
    %         
    %         k_1(n)=k(1,1); k_2(n) = k(1,2);  k_3(n) = k(1,3);  k_4(n) = -k(1,4); 
    %         
    %         e(n) = uc(n) - r3(n);  z(n) = z(n-1) + e(n); 
    %         
    %         q_ramp(n) = k_4(n)*z(n) - k_1(n)*r3(n) - k_2(n)*r2(n) - k_3(n)*r4(n); 
    % 
    %         if q_ramp(n)<0; q_ramp(n)=0;  elseif q_ramp(n)>=30; q_ramp(n)=30; end  
    %         
    %     end
    %     waitbar(n/length(q1),h);    
    % end
    % close(h)
    % 
    % figure(11);hold; plot(uc,'r'); plot(r3,'k'); plot(e)
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    % %(2) PARAMETRIC
    % par_A=zeros(size(q1)); par_B=par_A; par_C=par_A; par_D=par_A; e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); state_D=zeros(size(q1)); state_D2=zeros(size(q1));
    % par_A_h=zeros(size(q1)); par_B_h=par_A; par_C_h=par_A; par_D_h=par_A; par_D2_h=par_A; r3=zeros(size(q1)); r3(1:2)=0; 
    % % q_ramp=[5*ones(500,1); kappa4'; 20*ones(360,1); kappa5'; 5*ones(500,1);];
    % uc=53.7*ones(length(q1),1);
    % 
    % h=waitbar(0,'Please wait.... Control in Progress (parametric)');
    % 
    % for n=3:length(r2);
    %     
    %   %link 1 [s-2]
    %     r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    %     q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
    %     v1(n)=ftdcurve(p1,r1(n),1); 
    %   %link 2 [s-1]
    %     r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    %     q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
    %     v2(n)=ftdcurve(p2,r2(n),1);
    %   %link 3 [s]
    %     r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    %     q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
    %     v3(n)=ftdcurve(p3,r3(n),1);
    %   %link 4 [s+1]
    %     r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    %     q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
    %     v4(n)=ftdcurve(p4,r4(n),1);
    %   %link 5 [s+2]
    %     r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
    %     q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    %     v5(n)=ftdcurve(p5,r5(n),1);
    %     
    %     if  (n>=250) & (n<1050) 
    %         
    %         % read-off States 
    %         state_A(n-2) = max(r3(n-2),r4(n-2)); %if state_A(n-2) <= min(x(:,1)); state_A(n-2)=min(x(:,1)); elseif state_A(n-2)>max(x(:,1)); state_A(n-2)=max(x(:,1));end 
    %         state_B(n-2) = max(r2(n-2),r3(n-2)); %if state_B(n-2) <= min(x(:,2)); state_B(n-2)=min(x(:,2)); elseif state_B(n-2)>max(x(:,2)); state_B(n-2)=max(x(:,2));end
    %         state_C(n-2) = max(r4(n-2),r5(n-2)); %if state_C(n-2) <= min(x(:,3)); state_C(n-2)=min(x(:,3)); elseif state_C(n-2)>max(x(:,3)); state_C(n-2)=max(x(:,3));end
    %         
    %         % estimated parametric parameters 
    %          par_A(n-2)=-FISnon_SISO(state_A(n-2), 1, numMFs(modsel),  [], VX,  rules,  {'gaussmf'});
    %          par_B(n-2)= FISnon_SISO(state_B(n-2), 1, numMFs(modsel2), [], VX2, rules2, {'gaussmf'});
    %          par_C(n-2)= FISnon_SISO(state_C(n-2), 1, numMFs(modsel3), [], VX3, rules3, {'gaussmf'});
    % 
    % %         F = [-par_A(n-2) par_B(n-2) par_C(n-2) 0;
    % %               0 0 0 0;
    % %               0 0 0 0;
    % %               par_A(n-2) -par_B(n-2) -par_C(n-2) 1];
    % 
    %           F = [0 0 0 0;
    %               -par_A(n-2) par_B(n-2) par_C(n-2) 0;
    %               0 0 0 0;
    %               par_A(n-2) -par_B(n-2) -par_C(n-2) 1];
    % 
    %             g =[0 pars(1,4) 0 -pars(1,4)];
    % 
    % %         g =[pars(1,4) 0 0 -pars(1,4)];
    %         
    %         Q = eye(size(F)); Q(end,end)=0.1; r=0.1;  k=dlqri(F,g',Q,r); 
    %         
    %         k_1(n)=k(1,1); k_2(n) = k(1,2);  k_3(n) = k(1,3);  k_4(n) = -k(1,4); 
    %         
    %         e(n) = uc(n) - r3(n);  z(n) = z(n-1) + e(n); 
    %         
    %         q_ramp(n) = k_4(n)*z(n) - k_2(n)*r3(n) - k_1(n)*r2(n) - k_3(n)*r4(n); 
    % 
    %          if q_ramp(n)<0; q_ramp(n)=0;  elseif q_ramp(n)>=30; q_ramp(n)=30; end 
    %         
    %     end
    %     waitbar(n/length(q1),h);    
    % end
    % close(h)
    % 
    % figure(11); plot(uc,'r'); plot(r3)
    % 
    % IAE=sum(abs(r3(250:1050)-uc(250:1050)))
    % resid=r3(250:1050)-uc(250:1050); index=find(resid>0); IAE2=sum(resid(index)) 
    % Control_effort=sum(abs(q_ramp(250:1050)))
    % load1=sum((UnmLoad(250:1050)-r2(250:1050)))./sum(r2(250:1050))%;load1=load1*100 
    % max_flow=sum(q2(250:1050)-flow_ups(250:1050))./sum(flow_ups(250:1050))%; max_flow=max_flow*100
    % max_flow2=sum(q3(250:1050)-flow_jun(250:1050))./sum(flow_jun(250:1050))%;max_flow2=max_flow2*100
    % chang_V3=sum(v3(250:1050)-old_v3(250:1050))./sum(old_v3(250:1050))
    
    
    
    
    
    
    % 
    % %(3) PURE DETERMINISTIC
    % par_A=zeros(size(q1)); par_B=par_A; par_C=par_A; par_D=par_A; e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); state_D=zeros(size(q1)); state_D2=zeros(size(q1));
    % par_A_h=zeros(size(q1)); par_B_h=par_A; par_C_h=par_A; par_D_h=par_A; par_D2_h=par_A; r3=zeros(size(q1)); r3(1:2)=0; 
    % q_ramp=[5*ones(500,1); kappa4'; 20*ones(360,1); kappa5'; 5*ones(500,1);]; uc=63*ones(length(q1),1);
    h=waitbar(0,'Please wait.... Control in Progress (pure deterministic)');
    
    for n=3:length(r2);
    %     
      %link 1 [s-2]
        r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
        q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1)*dt;
        v1(n)=ftdcurve(p1,r1(n),1); 
      %link 2 [s-1]
        r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
        q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1)*dt;
        v2(n)=ftdcurve(p2,r2(n),1);
        B(n)=k3*alpha2*min(v2(n-2),v3(n-2));
      %link 3 [s]
        r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
        q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1)*dt;
        v3(n)=ftdcurve(p3,r3(n),1);
        A(n)=-(k3-k4)*alpha3*min(v3(n-2),v4(n-2));
      %link 4 [s+1]
        r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
        q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1)*dt;
        v4(n)=ftdcurve(p4,r4(n),1);
        C(n)=-k4*alpha4*min(v4(n-2),v5(n-2));
      %link 5 [s+2]
        r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
        q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1)*dt;
        v5(n)=ftdcurve(p5,r5(n),1);
        
        if  (n>=250) & (n<1050) 
            
            % read-off States 
            state_A(n-2) = max(r3(n-2),r4(n-2)); if state_A(n-2) <= min(x(:,1)); state_A(n-2)=min(x(:,1)); elseif state_A(n-2)>max(x(:,1)); state_A(n-2)=max(x(:,1));end 
            state_B(n-2) = max(r2(n-2),r3(n-2)); if state_B(n-2) <= min(x(:,2)); state_B(n-2)=min(x(:,2)); elseif state_B(n-2)>max(x(:,2)); state_B(n-2)=max(x(:,2));end
            state_C(n-2) = max(r4(n-2),r5(n-2)); if state_C(n-2) <= min(x(:,3)); state_C(n-2)=min(x(:,3)); elseif state_C(n-2)>max(x(:,3)); state_C(n-2)=max(x(:,3));end
            
            % estimated non-parametric parameters 
            par_A_h(n-2) = -interp1(x1_irw,  A_i',  state_A(n-2),'spline'); 
            par_B_h(n-2) =  interp1(x2_irw,  B_i', state_B(n-2),'spline');
            par_C_h(n-2) = -interp1(x3_irw,  C_i',  state_C(n-2),'spline');
            
%             F = [par_A_h(n-2) par_B_h(n-2) par_C_h(n-2) 0;
%                   0 0 0 0;
%                   0 0 0 0;
%                  -par_A_h(n-2) -par_B_h(n-2) -par_C_h(n-2) 1];
%     
%             g =[k3 0 0 -k3];
                       
            F = [ 0 0 0 0;
                par_B_h(n-2) par_A_h(n-2) par_C_h(n-2) 0;
                0 0 0 0;
                -par_B_h(n-2) -par_A_h(n-2) -par_C_h(n-2) 1];
            
            g =[0 k3 0 -k3];
            
            Q = eye(size(F)); Q(end,end)=1; r=0.01;  k=dlqri(F,g',Q,r); 

            k_1(n)=k(1,1); k_2(n) = k(1,2);  k_3(n) = k(1,3);  k_4(n) = -k(1,4); 
                        
             e(n) = uc(n) - r3(n);  
        
           % z(n) = z(n-1) + e(n);  q_ramp(n) = k_4(n)*z(n) - k_1(n)*r3(n) - k_2(n)*r2(n) - k_3(n)*r4(n); 
           % q_ramp(n) = k_4(n)*z(n) - k_1(n)*r3(n) - k_2(n)*r2(n) - k_3(n)*r4(n); 
 
        q_ramp(n)=q_ramp(n-1) + k_4(n)*e(n) - k_2(n)*(r3(n)-r3(n-1)) - k_1(n)*(r2(n)-r2(n-1)) - k_3(n)*(r4(n)-r4(n-1));  

             if q_ramp(n)<0; q_ramp(n)=0;  elseif q_ramp(n)>=30; q_ramp(n)=30; end 
            
        end
        waitbar(n/length(q1),h);    
    end
    close(h)
    
    figure(11); plot(uc,'r'); plot(r3,'k')
    
    IAE_det=sum(abs(r3(250:1050)-uc(250:1050)));
    resid=r3(250:1050)-uc(250:1050); index=find(resid>0); IAE2_det=sum(resid(index));
    Control_effort_det=sum(abs(q_ramp(250:1050)));
    % load1=sum((UnmLoad(250:1050)-r2(250:1050)))./sum(r2(250:1050))%;load1=load1*100 
    % max_flow=sum(q2(250:1050)-flow_ups(250:1050))./sum(flow_ups(250:1050))%; max_flow=max_flow*100
    % max_flow2=sum(q3(250:1050)-flow_jun(250:1050))./sum(flow_jun(250:1050))%;max_flow2=max_flow2*100
    % chang_V3=sum(v3(250:1050)-old_v3(250:1050))./sum(old_v3(250:1050))
    % 
    
    
    % % *** TVP Controller ***
     uc=50*ones(length(q1),1); 
    load tvp_pars
    h=waitbar(0,'Please wait.... Control in Progress (TVP)');
    for n=3:length(r2);   
        
      %link 1 [s-2]
        r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
        q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1)*dt;
        v1(n)=ftdcurve(p1,r1(n),1); 
      %link 2 [s-1]
        r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
        q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1)*dt;
        v2(n)=ftdcurve(p2,r2(n),1);
        B(n)=k3*alpha2*min(v2(n-2),v3(n-2));
      %link 3 [s]
        r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
        q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1)*dt;
        v3(n)=ftdcurve(p3,r3(n),1);
        A(n)=-(k3-k4)*alpha3*min(v3(n-2),v4(n-2));
      %link 4 [s+1]
        r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
        q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1)*dt;
        v4(n)=ftdcurve(p4,r4(n),1);
        C(n)=-k4*alpha4*min(v4(n-2),v5(n-2));
      %link 5 [s+2]
        r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
        q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1)*dt;
        v5(n)=ftdcurve(p5,r5(n),1);
        
        if  (n>=250) & (n<1050)         
    
           %Formulation of the NMSS matrix
            F = [-par_1_tvp(n,1)    par_1_tvp(n,2)      0             0;...
                  par_2_tvp(n,2)   -par_2_tvp(n,1)   par_2_tvp(n,3)   0;...
                        0           par_3_tvp(n,2)  -par_3_tvp(n,1)   0;...
                 -par_2_tvp(n,2)    par_2_tvp(n,1)  -par_2_tvp(n,3)   1;];
        
    %        F = [aa(1)  bb(1)      0  0;...  for reference
    %             bb(2)  aa(2)   cc(2) 0;...
    %                0   bb(3)   aa(3) 0;...
    %            -bb(2) -aa(2)  -cc(2) 1;];
    
    %        g = [0 k3 0 -k3]';
    
         g = [0 par_2_tvp(1,4) 0 -par_2_tvp(1,4)]';
          
            Q=eye(size(F)); r=0.01; Q(end, end)=1;
            
            k = dlqr(F,g,Q,r); k(end)=-k(end);
            
            e(n) = uc(n) - r3(n);  
    
            q_ramp(n) = q_ramp(n-1) + k(4)*e(n) - k(1)*(r2(n)-r2(n-1)) - k(2)*(r3(n)-r3(n-1)) - k(3)*(r4(n)-r4(n-1));

             if q_ramp(n)<0; q_ramp(n)=0;  elseif q_ramp(n)>=30; q_ramp(n)=30; end       
    
        end
    waitbar(n/length(q1),h);
    end
    close(h)
    figure(11); plot(r3,'m')
    disp('TVP controller');
    IAE_tvp=sum(abs(r3(250:1050)-uc(250:1050)));
    resid=r3(250:1050)-uc(250:1050); index=find(resid>0); IAE2_tvp=sum(resid(index)); 
    Control_effort_tvp=sum(abs(q_ramp(250:1050)));
    % load2=sum((UnmLoad(250:1050)-r2(250:1050)))./sum(r2(250:1050))%load2=load2*100 
    % 
    % max_flow=sum(q2(250:1050)-flow_ups(250:1050))./sum(flow_ups(250:1050))%; max_flow=max_flow*100
    % max_flow2=sum(q3(250:1050)-flow_jun(250:1050))./sum(flow_jun(250:1050))%;max_flow2=max_flow2*100
     chang_V3_tvp=sum(v3(250:1050)-old_v3(250:1050))./sum(old_v3(250:1050));
    
    
    
    
    % *** LLM Controller on SDP model ***
    %q_ramp=[5*ones(500,1); kappa4'; 20*ones(360,1); kappa5'; 5*ones(500,1);];
     uc=50*ones(length(q1),1);
    
    h=waitbar(0,'Please wait.... Control in Progress (LLM)');
    for n=3:length(r2);   
        
    %link 1 [s-2]
    r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1)*dt;
    v1(n)=ftdcurve(p1,r1(n),1); 
    %link 2 [s-1]
    r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1)*dt;
    v2(n)=ftdcurve(p2,r2(n),1);
    B(n)=k3*alpha2*min(v2(n-2),v3(n-2));
    %link 3 [s]
    r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1)*dt;
    v3(n)=ftdcurve(p3,r3(n),1);
    A(n)=-(k3-k4)*alpha3*min(v3(n-2),v4(n-2));
    %link 4 [s+1]
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4); 
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1)*dt;
    v4(n)=ftdcurve(p4,r4(n),1);
    C(n)=-k4*alpha4*min(v4(n-2),v5(n-2));
    %link 5 [s+2]
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1)*dt;
    v5(n)=ftdcurve(p5,r5(n),1);
    
        if  (n>=250) & (n<1050)      
         
            e(n) = uc(n) - r3(n);
            
            q_ramp(n) = q_ramp(n-1) + k_llm(4)*e(n) - k_llm(1)*(r2(n)-r2(n-1)) - k_llm(2)*(r3(n)-r3(n-1)) - k_llm(3)*(r4(n)-r4(n-1));
            
            %           q_ramp(n) = q_ramp(n-1) + k_llm(4)*e(n) - k_llm(1)*(r3(n)-r3(n-1)) - k_llm(2)*(r2(n)-r2(n-1)) - k_llm(3)*(r4(n)-r4(n-1)); 
            
            if q_ramp(n)<0; q_ramp(n)=0;  elseif q_ramp(n)>=30; q_ramp(n)=30; end       
            
        end
        
 waitbar(n/length(q1),h);
 
    end
 
    close(h)
    
    figure(11); plot(r3,'r')
    disp('LLM');
    IAE_llm=sum(abs(r3(250:1050)-uc(250:1050)));
    resid=r3(250:1050)-uc(250:1050); index=find(resid>0); IAE2_llm=sum(resid(index)); 
    Control_effort_llm=sum(abs(q_ramp(250:1050))); 
    % load2=sum((UnmLoad(250:1050)-r2(250:1050)))./sum(r2(250:1050))%load2=load2*100 
    % max_flow=sum(q2(250:1050)-flow_ups(250:1050))./sum(flow_ups(250:1050)); max_flow=max_flow*100
    % max_flow2=sum(q3(250:1050)-flow_jun(250:1050))./sum(flow_jun(250:1050));max_flow2=max_flow2*100
      chang_V3_llm=sum(v3(250:1050)-old_v3(250:1050))./sum(old_v3(250:1050));
    
    % figure(14);plot(q3,'r')
    
%     figure(12);
%     subplot(321); plot(r1,q1,'m.'); plot(p1(4),max(q1),'y.');title('s-2');
%     subplot(322); plot(r2,q2,'m.'); plot(p2(4),max(q2),'y.');title('s-1 {load}');
%     subplot(323); plot(r3,q3,'m.'); plot(p3(4),max(q3),'y.');title('s {junction}');
%     subplot(324); plot(r4,q4,'m.'); plot(p4(4),max(q4),'y.');title('s+1');
%     subplot(325); plot(r5,q5,'m.'); plot(p5(4),max(q5),'y.');title('s+2');
    
    
    
    
    
    
    
    % 
    % % *** Demand Capacity Controller on SDP model ***
    % q_ramp=[5*ones(500,1); kappa4'; 15*ones(360,1); kappa5'; 5*ones(500,1);];
    % o_crit=occ2den(63,4,2,3,0); q_cap=60.75; o4=zeros(size(r2));
    % q_ramp=[5*ones(500,1); kappa4'; 15*ones(360,1); kappa5'; 5*ones(500,1);];
    % 
    % h=waitbar(0,'Please wait.... Control in Progress (Demand-Capacity)');
    % for n=3:length(r2);   
    %     
    %     %link 1 [s-2]
    %     r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    %     q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
    %     v1(n)=ftdcurve(p1,r1(n),1); 
    %     %link 2 [s-1]
    %     r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    %     q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
    %     v2(n)=ftdcurve(p2,r2(n),1);
    %     %link 3 [s]
    %     r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    %     q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
    %     v3(n)=ftdcurve(p3,r3(n),1);
    %     %link 4 [s+1]
    %     r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    %     q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
    %     v4(n)=ftdcurve(p4,r4(n),1);
    %     o4(n)=occ2den(r4(n),4,2,3,0);
    %     %link 5 [s+2]
    %     r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
    %     q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    %     v5(n)=ftdcurve(p5,r5(n),1);
    %     
    %     if  (n>=250) & (n<1050)         
    % 
    %         if o4(n)<= o_crit
    %             
    %             q_ramp(n)=q_cap-q2(n-1);
    %             
    %         else
    %             q_ramp(n)=5; %minimal size
    %             
    %         end
    %         waitbar(n/length(q1),h);
    %     end
    %    
    % end
    %  close(h)
    % figure(1); plot(r3,'y')
    % 
    % IAE=sum(abs(r3(250:1050)-uc(250:1050)))
    % resid=r3(250:1050)-uc(250:1050); index=find(resid>0); IAE2=sum(resid(index)) 
    % Control_effort=sum(abs(q_ramp(250:1050)))
    % load2=sum((UnmLoad(250:1050)-r2(250:1050)))./sum(r2(250:1050))%load2=load2*100 
    % max_flow=sum(q2(250:1050)-flow_ups(250:1050))./sum(q2(250:1050))%; max_flow=max_flow*100
    % max_flow2=sum(q3(250:1050)-flow_jun(250:1050))./sum(q3(250:1050))%;max_flow2=max_flow2*100
    %        
    %         
    
    
    
    
    
    % 
    % 
    % *** ALINEA Controller on SDP model ***
%     q_ramp=[5*ones(500,1); kappa4'; 15*ones(360,1); kappa5'; 5*ones(500,1);];
    %o_crit=occ2den(63,4,2,3,0);
    q_cap=max(flow_jun); o4=zeros(size(r2));
    %o_crit=8.4464; q_cap=70.4851;
    %q_ramp=[5*ones(500,1); kappa4'; 15*ones(360,1); kappa5'; 5*ones(500,1);];
    
    h=waitbar(0,'Please wait.... Control in Progress (ALINEA)');
    for n=3:length(r2);   
        
        %link 1 [s-2]
        r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
        q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
        v1(n)=ftdcurve(p1,r1(n),1); 
        %link 2 [s-1]
        r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
        q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
        v2(n)=ftdcurve(p2,r2(n),1);
        %link 3 [s]
        r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
        q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
        v3(n)=ftdcurve(p3,r3(n),1);
        o3(n)=occ2den(r3(n),4,2,3,0);
        %link 4 [s+1]
        r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
        q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
        v4(n)=ftdcurve(p4,r4(n),1);
        %link 5 [s+2]
        r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);
        q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
        v5(n)=ftdcurve(p5,r5(n),1);
        
    %     if r3(n)>80; break; end
        
        if  (n>=250) & (n<1050)         
           
            q_ramp(n)= q_ramp(n-1) + 70*(50-r3(n));
                
             if q_ramp(n)<0; q_ramp(n)=0;  elseif q_ramp(n)>=30; q_ramp(n)=30; end 
    
            waitbar(n/length(q1),h);
        end
       
    end
     close(h)
     
    figure(11); plot(r3,'y')
    disp('Alinea');
    IAE_alina=sum(abs(r3(250:1050)-uc(250:1050)));
    resid=r3(250:1050)-uc(250:1050); index=find(resid>0); IAE2_alinea=sum(resid(index)); 
    Control_effort_alinea=sum(abs(q_ramp(250:1050)));
%     load2=sum((UnmLoad(250:1050)-r2(250:1050)))./sum(r2(250:1050))%load2=load2*100 
%     max_flow=sum(q2(250:1050)-flow_ups(250:1050))./sum(flow_ups(250:1050))%; max_flow=max_flow*100
%     max_flow2=sum(q3(250:1050)-flow_jun(250:1050))./sum(flow_jun(250:1050))%;max_flow2=max_flow2*100
    chang_V3_alinea=sum(v3(250:1050)-old_v3(250:1050))./sum(old_v3(250:1050));
    
    
ALINEA=['    ALINEA     '  num2str(IAE_alina)  '  ' num2str(IAE2_alinea)  '    '  num2str(Control_effort_alinea)];
LLM   =['      LLM        '  num2str(IAE_llm)  '   ' num2str(IAE2_llm)  '    '  num2str(Control_effort_llm)];
TVP   =['      TVP        '  num2str(IAE_tvp)  '   ' num2str(IAE2_tvp)  '    '  num2str(Control_effort_tvp)];
SDP   =['      SDP        '  num2str(IAE_sdp)  '   ' num2str(IAE2_sdp)  '    '  num2str(Control_effort_sdp)];
SDPdet=['    SDP-Det      '  num2str(IAE_det)  '   ' num2str(IAE2_det)  '    '  num2str(Control_effort_det)];
fprintf('    Controller      IAE      CARS      Ramp Effort   \n  ');
disp(ALINEA)
disp(LLM)
disp(TVP)
disp(SDP)
disp(SDPdet)
