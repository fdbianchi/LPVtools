
% LPVtools: Two-link Flexible Manipulator Example (control using griding)
% 
% Ref: Apkarian & Adams. Advanced gain-scheduling techniques for  uncertain 
% systems. IEEE Trans. on Control Systems Technology, 6(1), 21–32. 1998
%
%   Plant:          PWA LPV (ppss)
%                   1 parameter,
%                   2 inputs, 2 outputs
%   Constraint:     Hinf
%   Lyapunov Fcn:   Constant/Parameter dependent
%   Ctrller Fcn:    No specified

% fbianchi - 2020-07-05

% cleaning
clearvars
close all
clc

% ========================================================================
% Modelling

% parameter set
nv     = 7;
range  = [0 pi];        % rad/sec
rate   = [-100 100];
points = pgrid(range,nv);
pv     = pset.Gral(points,rate,'theta');
% plot(pv)

% Nonlinear model
% inertia matrix
M_pi2 = [34.7077  9.7246 23.6398  5.9114;
          9.7246  9.8783  9.7246  5.9114;
         23.6398  9.7246 17.5711  5.9114;
          5.9114  5.9114  5.9114  3.7233];

M_pi  = [17.0296  0.8856  9.7776  0.8430;
          0.8856  9.8783  4.7016  5.9114;
          9.7776  4.7016  7.5242  3.0311;
          0.8430  5.9114  3.0311  3.7233];
M =@(p) M_pi2 + cos(p)*(M_pi2-M_pi);
% damping
Dp = diag([0 0 0.09 0.05]);
% sitfness
K = diag([0 0 89.1473 45.6434]);
% torque matrix
F = [eye(2);zeros(2)];

A =@(p) [zeros(4),  eye(4);
        -M(p)\K,   -M(p)\Dp];
B =@(p) [zeros(4,2);
         M(p)\F];
C = eye(2,8);

a = zeros(8,8,nv);
b = zeros(8,2,nv);
for ii = 1:nv
    p = points(ii);
    a(:,:,ii) = A(p);
    b(:,:,ii) = B(p);
end
pdG = ppss(a,b,C,0,pv,'InputName','u','OutputName','y');

% frequency response
nf = 250;
w  = logspace(-2,2,nf);
figure
subPoints = [0 pi/2 pi];
lcolor = lines;
for ii = 1:length(subPoints)
    p = subPoints(ii);
    a = A(p);
    b = B(p);
    c = C;
    d = zeros(2);
    G = ss(a,b,c,d);
    sv = sigma(G,w);
    hl = loglog(w,sv,'Color',lcolor(ii,:));
    Hl(ii) = hl(1);
    hold on
end
legend(Hl,'G(0)','G(pi/2)','G(pi)');
ylabel('\sigma(G)')
xlabel('Freq(rad/sec)')
title('Plant frequency response')

% ------------------------------------------------------------------------
% controller design
s = tf('s');

% weithing functions 
Wf = balreal(minreal(4*(s+0.1)^2/(s+100)^2*eye(2)));
Wf.u = 'u';      Wf.y = {'z(3)' 'z(4)'};
% Wp = balreal(minreal(0.1075*(s+1.066)/(s+0.03)*eye(2)));
Wp = balreal(minreal(0.1075*(s+1.066)/(s+0.03)*eye(2)));
Wp.u = 'e';     Wp.y = {'z(1)' 'z(2)'};

% to deal with B2(p)
ctrlfcn.pdOut = false;
Wb = ss(eye(2));
% ctrlfcn.pdOut = true;
% Wb = balreal(minreal(1/(0.001*s+1)*eye(2)));
Wb.u = 'u';     Wb.y = 'v';
pdG.u = 'v';

% augmented plant
sb1 = sumblk('e = y + w',2);
s1 = sumblk('ec(1) = e(1) + w(3)');
s2 = sumblk('ec(2) = e(2) + w(4)');
sb2 = append(s1,s2);
pdGau = connect(pdG,sb1,sb2,Wb,Wp,Wf,{'w','u'},{'z','ec'});


% -------------------------------------------------------------------------
% Controller synthesis: H-infinity
ys = sprintfc('ec(%d)',1:2);     % y = 5:6;
us = sprintfc('u(%d)',1:2);      % u = 5:6;
zs = sprintfc('z(%d)',1:4);      % z = 1:4;
ws = sprintfc('w(%d)',1:4);      % w = 1:4;

% Hinf controller @ theta = pi/2
Gau = ss(pdGau,pi/2);
[Khinf,CL,ghinf,INFO] = hinfsyn(Gau,2,2,'method','lmi');

% LPV designs:
opts = lpvsettings('solver','mosek');
const(1) = synConst.Gain(ws,zs);
const(2) = synConst.Poles('MaxFreq',1e3);

% case 0: constant Lyapunov functions
[pdK0,constR] = lpvsyn(pdGau,ys,us,const,[],ctrlfcn,opts);
g0 = constR(1).bound;

% case 1: parameter dependent Lyapunov functions, dX=0
lyapSet.parfcnX  =@(p) cos(p);
lyapSet.dparfcnX = 0;
lyapSet.parfcnY  =@(p) cos(p);
lyapSet.dparfcnY = 0;
[pdK1,constR] = lpvsyn(pdGau,ys,us,const,lyapSet,[],opts);
g1 = constR(1).bound;

% case 2: PWA parameter dependent Lyapunov functions, dX=0
[pdK2,constR] = lpvsyn(pdGau,ys,us,const,'pwa',ctrlfcn,opts);
g2 = constR(1).bound;

fprintf('\n--------------------------------------\n')
fprintf('Comparison with no scalings and dX/dt=0\n')
fprintf('hinfsyn:                  %5.4f\n',ghinf)
fprintf('X & Y cte:                %5.4f\n',g0)
fprintf('X(p) & Y(p):              %5.4f\n',g1)
fprintf('X(p) & Y(p) (PWA):        %5.4f\n',g2)
fprintf('--------------------------------------\n')


%% =======================================================================
% Using scalings
% 
% scaling
S0l = blkdiag(0.064*eye(2),eye(4));
S0r = inv(S0l);
pdGauSc = S0l*pdGau*S0r;
pdGauSc.u = [ws, us];
pdGauSc.y = [zs, ys];

% Hinf controller @ theta = pi/2
Gau = ss(pdGauSc,pi/2);
[Khinf2,CL,ghinf,INFO] = hinfsyn(Gau,2,2,'method','lmi');

% LPV designs:
% case 0: constant Lyapunov functions
opts = lpvsettings('solver','mosek');%,'XYpenalty',1,'penalty',1e-4);
[pdK02,constR] = lpvsyn(pdGauSc,ys,us,const,[],ctrlfcn,opts);
g0 = constR(1).bound;

% case 1: parameter dependent Lyapunov functions, dX=0
% X=cte, Y(p)
lyapSet.parfcnX  = 0;
lyapSet.dparfcnX = 0;
lyapSet.parfcnY  =@(p) cos(p);
lyapSet.dparfcnY = 0;
[pdK22,constR] = lpvsyn(pdGauSc,ys,us,const,lyapSet,[],opts);
g22 = constR(1).bound;
 
% X(p), Y=cte
lyapSet.parfcnX  =@(p) cos(p);
lyapSet.dparfcnX = 0;
lyapSet.parfcnY  = 0;
lyapSet.dparfcnY = 0;
[pdK23,constR] = lpvsyn(pdGauSc,ys,us,const,lyapSet,[],opts);
g23 = constR(1).bound;
 
% X(p), Y(p)
lyapSet.parfcnX  =@(p) cos(p);
lyapSet.dparfcnX = 0;
lyapSet.parfcnY  =@(p) cos(p);
lyapSet.dparfcnY = 0;
[pdK24,constR] = lpvsyn(pdGauSc,ys,us,const,lyapSet,[],opts);
g24 = constR(1).bound;

fprintf('Comparison scaling (dX=0, dY=0)\n')
fprintf('hinfsyn:                  %5.4f\n',ghinf)
fprintf('X=cte & Y=cte:            %5.4f\n',g0)
fprintf('X=cte & Y(p):             %5.4f\n',g22)
fprintf('X(p)  & Y=cte:            %5.4f\n',g23)
fprintf('X(p)  & Y(p):             %5.4f\n',g24)
fprintf('--------------------------------------\n')

% controller frequency response
nf = 250;
wc = logspace(-2,2,nf);
figure
subPoints = [0 pi/2 pi];
for ii = 1:length(subPoints)
    p = subPoints(ii);
    Ki = subs(pdK23,p);
    sv = sigma(Ki,wc);
    hl = loglog(wc,sv,'Color',lcolor(ii,:));
    Hl(ii) = hl(1);
    hold on
end
legend(Hl,'K(0)','K(pi/2)','K(pi)');
ylabel('\sigma(G)')
xlabel('Freq(rad/sec)')
title('Controller frequency response (X(p)  & Y=cte)')

%% =======================================================================
% Checking different parameter variation rates
% 

lyapSet.parfcnX  =@(p) cos(p);
lyapSet.dparfcnX =@(p) -sin(p);
lyapSet.parfcnY  = 0;
lyapSet.dparfcnY = 0;
rateSet = [0 10 100 1000 10000];
opts = lpvsettings('solver','mosek','XYpenalty',1,'penalty',1e-4);
for ii = 1:length(rateSet) 
    pdGauSc.parset.rate = rateSet(ii)*[-1 1];
    [pdK3,constR] = lpvsyn(pdGauSc,ys,us,const,lyapSet,[],opts);
    g3(ii) = constR(1).bound;
end
figure
semilogx(rateSet,g3,'rs')
xlabel('rate of parameters (deg/sec)')
ylabel('gamma')
title('Performance vs parameter rate')


%% =======================================================================
% Non linear simulations

p0 = pi/2;

simFile = 'TwoLink_mdl.slx';
load_system(simFile);

% LTI controller no scaling
pdK = ppss(Khinf2);

set_param([simFile(1:end-4) '/LPV Controller'],'Commented','off')
set_param([simFile(1:end-4) '/LPV Controller X(p)'],'Commented','on')
set_param([simFile(1:end-4) '/sw'],'sw','1')

sim(simFile);
t = SimData.time;
ref1 = SimData.signals(1).values(:,1);
theta1 = SimData.signals(1).values(:,2);
ref2 = SimData.signals(2).values(:,1);
theta2 = SimData.signals(2).values(:,2);
uc = SimData.signals(3).values;
p = SimData.signals(4).values;

% color lines
clines = lines(7);
figure('Position', [680   10   800   780])

ha1 = subplot(4,1,1);
plot(t,ref1,'Color',clines(1,:),'LineWidth',1); hold on
plot(t,theta1,'Color',clines(2,:),'LineWidth',1);
ylabel('\theta_1'); xlabel('time (s)')
title(sprintf('LTI controller, with scaling, @ p=%5.2f',p0))
ha2 = subplot(4,1,2);
plot(t,ref2,'Color',clines(1,:),'LineWidth',1); hold on
plot(t,theta2,'Color',clines(2,:),'LineWidth',1);
ylabel('\theta_2'); xlabel('time (s)')
ha3 = subplot(4,1,3);
plot(t,uc(:,1),'Color',clines(1,:),'LineWidth',1); hold on
plot(t,uc(:,2),'Color',clines(2,:),'LineWidth',1);
ylabel('u'); xlabel('time (s)')
ha4 = subplot(4,1,4);
plot(t,p,'Color',clines(1,:),'LineWidth',1); hold on
ylabel('p'); xlabel('time (s)')

% LPV controller no scaling (X and Y cte)
pdK = pdK02;

set_param([simFile(1:end-4) '/LPV Controller'],'Commented','off')
set_param([simFile(1:end-4) '/LPV Controller X(p)'],'Commented','on')
set_param([simFile(1:end-4) '/sw'],'sw','1')

sim(simFile);
t = SimData.time;
ref1 = SimData.signals(1).values(:,1);
theta1 = SimData.signals(1).values(:,2);
ref2 = SimData.signals(2).values(:,1);
theta2 = SimData.signals(2).values(:,2);
uc = SimData.signals(3).values;
p = SimData.signals(4).values;

% color lines
clines = lines(7);
figure('Position', [680   10   800   780])

ha1 = subplot(4,1,1);
plot(t,ref1,'Color',clines(1,:),'LineWidth',1); hold on
plot(t,theta1,'Color',clines(2,:),'LineWidth',1);
ylabel('\theta_1'); xlabel('time (s)')
title(sprintf('LPV controller, no scaling, X and Y cte, @ p=%5.2f',p0))
ha2 = subplot(4,1,2);
plot(t,ref2,'Color',clines(1,:),'LineWidth',1); hold on
plot(t,theta2,'Color',clines(2,:),'LineWidth',1);
ylabel('\theta_2'); xlabel('time (s)')
ha3 = subplot(4,1,3);
plot(t,uc(:,1),'Color',clines(1,:),'LineWidth',1); hold on
plot(t,uc(:,2),'Color',clines(2,:),'LineWidth',1);
ylabel('u'); xlabel('time (s)')
ha4 = subplot(4,1,4);
plot(t,p,'Color',clines(1,:),'LineWidth',1); hold on
ylabel('p'); xlabel('time (s)')

% LPV controller with scaling (X(p) and Y cte)
pdK = pdK23;

set_param([simFile(1:end-4) '/LPV Controller'],'Commented','on')
set_param([simFile(1:end-4) '/LPV Controller X(p)'],'Commented','off')
set_param([simFile(1:end-4) '/sw'],'sw','0')

sim(simFile);
t = SimData.time;
ref1 = SimData.signals(1).values(:,1);
theta1 = SimData.signals(1).values(:,2);
ref2 = SimData.signals(2).values(:,1);
theta2 = SimData.signals(2).values(:,2);
uc = SimData.signals(3).values;
p = SimData.signals(4).values;

% color lines
clines = lines(7);
figure('Position', [680   10   800   780])

ha1 = subplot(4,1,1);
plot(t,ref1,'Color',clines(1,:),'LineWidth',1); hold on
plot(t,theta1,'Color',clines(2,:),'LineWidth',1);
ylabel('\theta_1'); xlabel('time (s)')
title(sprintf('LPV controller, with scaling, X(p) and Y=cte, @ p=%5.2f',p0))
ha2 = subplot(4,1,2);
plot(t,ref2,'Color',clines(1,:),'LineWidth',1); hold on
plot(t,theta2,'Color',clines(2,:),'LineWidth',1);
ylabel('\theta_2'); xlabel('time (s)')
ha3 = subplot(4,1,3);
plot(t,uc(:,1),'Color',clines(1,:),'LineWidth',1); hold on
plot(t,uc(:,2),'Color',clines(2,:),'LineWidth',1);
ylabel('u'); xlabel('time (s)')
ha4 = subplot(4,1,4);
plot(t,p,'Color',clines(1,:),'LineWidth',1); hold on
ylabel('p'); xlabel('time (s)')
