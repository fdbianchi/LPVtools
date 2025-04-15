
% LPVtools: Two-link Flexible Manipulator Example (control using griding)
% 
% Ref: Apkarian & Adams. Advanced gain-scheduling techniques for  uncertain 
% systems. IEEE Trans. on Control Systems Technology, 6(1), 21–32. 1998
%
%   Plant:          general LPV (pgss)
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
% damping
Dp = diag([0 0 0.09 0.05]);
% sitfness
K = diag([0 0 89.1473 45.6434]);
% torque matrix
F = [eye(2);zeros(2)];

% LFT description
M0 = M_pi2;
Mp = M_pi2 - M_pi;
[u,s,v] = svds(Mp,2);
L = u*sqrt(s);  W = sqrt(s)*v'; 

A = [zeros(4), eye(4);...
    -M0\[K Dp]];
B1 = [zeros(4,2); -M0\F];
B2 = [zeros(4,2); -M0\L];
C1 = eye(2,8);
C2 = -W*(M0\[K Dp]);
D21 = -W*(M0\F);
D22 = -W*(M0\L);
   
% expressing the nonlinear model as sum of terms
V1 = [1 0; 0 0];
V2 = [0 1;-1 0];
V3 = [0 0; 0 1];
Ap = cat(3,A,B2*V1*C2,B2*V2*C2,B2*V3*C2);
Bp = cat(3,B1,B2*V1*D21,B2*V2*D21,B2*V3*D21);

% parameter set
nv     = 7;
range  = [0 pi];        % rad/sec
rate   = [-100 100];
points = pgrid(range,nv);
pv     = pset.Grid(range,nv,rate,'theta');

% parameter dependence
E = eye(4); E(3,:) = [];
phi = -W*(M0\L);
fcn =@(p) E*reshape(inv(eye(2)-phi*cos(p))*cos(p),4,1);

% lpv plant
pdG = pgss(Ap,Bp,C1,0,pv,fcn,'InputName','u','OutputName','y');

% checking results
M =@(p) M_pi2 + cos(p)*(M_pi2-M_pi);
Af =@(p) [zeros(4),  eye(4);
        -M(p)\K,   -M(p)\Dp];
Bf =@(p) [zeros(4,2);
        -M(p)\F];
for ii = 1:nv
    p = points(ii);
    Gtest(:,:,ii) = ss(Af(p),Bf(p),C1,0);
end
Gtest.u = 'u'; Gtest.y = 'y';
% sigma(pdG,Gtest)


% frequency response
nf = 250;
w  = logspace(-2,2,nf);
figure
subPoints = [0 pi/2 pi];
lcolor = lines;
% Glpv = ss(pdG,subPoints);
Glpv = Gtest(:,:,1:3:7);
for ii = 1:length(subPoints)
    sv = sigma(Glpv(:,:,ii),w);
    hl = loglog(w,sv,'Color',lcolor(ii,:));
    Hl(ii) = hl(1);
    hold on
end
legend(Hl,'$G(0)$','$G(\pi/2)$','$G(\pi)$','interpreter','latex');
ylabel('\sigma(G)')
xlabel('Freq(rad/sec)')


% ------------------------------------------------------------------------
% controller design

% augmented plant
sb1 = sumblk('e = y + w',2);
s1 = sumblk('ec(1) = e(1) + w(3)');
s2 = sumblk('ec(2) = e(2) + w(4)');
sb2 = append(s1,s2);
% for designing
pdGau = connect(pdG,sb1,sb2,{'w','u'},{'e','u','ec'});
% for testing
pdGaut = connect(pdG,sb1,{'w','u'},{'y','e'});

% design #1: H-infitity for plant @ pi/2
%
% weigths:
s = tf('s');
Wf = 4*(s+0.1)^2/(s+100)^2*eye(2);
Wp = 0.1075*(s+1.066)/(s+0.03)*eye(2);
Wout = balreal(append(Wp,Wf,eye(2)));
pdGaw = Wout*pdGau;
% scaling
S0l = blkdiag(0.064*eye(2),eye(4));
S0r = inv(S0l);
Gau = S0l*ss(pdGaw,pi/2)*S0r;
% design
[Khinf,CL,ghinf,INFO] = hinfsyn(Gau,2,2,'method','lmi');
% test
figure
sv = sigma(Khinf,w);
hl = loglog(w,sv,'Color',lcolor(1,:));
ylabel('$\sigma(G)$','interpreter','latex')
xlabel('Freq(rad/sec)')
title('Singular values H-infinity controller (design #1)')

Gaut = ss(pdGaut,subPoints);
t = linspace(0,20,200);
y11 = zeros(200,3); y12 = y11;
for ii = 1:length(subPoints)
    Gcl = lft(Gaut(:,:,ii),Khinf);
    yaux = step(Gcl,t);
    y11(:,ii) = -yaux(:,1,1);
    y12(:,ii) = -yaux(:,2,2);
end
figure
subplot(2,1,1)
plot(t,y11)
ylabel('$\theta_1$','interpreter','latex')
title('Step response H-infinity controller (with scaling)')
subplot(2,1,2)
plot(t,y12)
ylabel('$\theta_2$','interpreter','latex')
xlabel('Time (s)')
legend('$G(0)$','$G(\pi/2)$','$G(\pi)$','interpreter','latex',...
    'location','southeast','Orientation','horizontal');
legend boxoff

% design #2: LPV
%
% weigths:
pdGaw = Wout*pdGau;
% scaling
pdGawSc = S0l*pdGaw*S0r;
% specifications
opts = lpvsettings('solver','mosek');
const(1) = synConst.Gain(1:4,1:4);
const(2) = synConst.Poles('MaxFreq',1e3);
% const(2) = synConst.Poles('MinDecay',0);

% case 21: constant Lyapunov functions
[pdK0,constR] = lpvsyn(pdGawSc,5:6,5:6,const,[],[],opts);
g0 = constR(1).bound;

% test
figure
sv = sigma(ss(pdK0,pi/2),w);
hl = loglog(w,sv,'Color',lcolor(1,:));
ylabel('$\sigma(G)$','interpreter','latex')
xlabel('Freq(rad/sec)')
title('Singular values LPV controller (with scaling, X & Y cte)')

pdGcl = lft(pdGaut,pdK0);
Gcl = ss(pdGcl,subPoints);
y11 = zeros(200,3); y12 = y11;
for ii = 1:length(subPoints)
    yaux = step(Gcl(:,:,ii),t);
    y11(:,ii) = -yaux(:,1,1);
    y12(:,ii) = -yaux(:,2,2);
end
figure
subplot(2,1,1)
plot(t,y11)
ylabel('$\theta_1$','interpreter','latex')
title('Step response LPV controller (with scaling, X & Y cte)')
subplot(2,1,2)
plot(t,y12)
ylabel('$\theta_2$','interpreter','latex')
xlabel('Time (s)')
legend('$G(0)$','$G(\pi/2)$','$G(\pi)$','interpreter','latex',...
    'location','southeast','Orientation','horizontal');
legend boxoff

% case 1: parameter dependent Lyapunov functions, dX=0
% lyapSet.parfcnX  =@(p) cos(p);
lyapSet.parfcnX  = fcn;
lyapSet.dparfcnX = 0;
% lyapSet.parfcnY  =@(p) cos(p);
lyapSet.parfcnY  = fcn;
lyapSet.dparfcnY = 0;
[pdK1,constR] = lpvsyn(pdGawSc,5:6,5:6,const,lyapSet,[],opts);
g1 = constR(1).bound;

% test
figure
sv = sigma(ss(pdK1,pi/2),w);
hl = loglog(w,sv,'Color',lcolor(1,:));
ylabel('$\sigma(G)$','interpreter','latex')
xlabel('Freq(rad/sec)')
title('Singular values LPV controller (with scaling, X(p) & Y(p))')

Gcl = lft(pdGaut,pdK1,subPoints);
y11 = zeros(200,3); y12 = y11;
for ii = 1:length(subPoints)
    yaux = step(Gcl(:,:,ii),t);
    y11(:,ii) = -yaux(:,1,1);
    y12(:,ii) = -yaux(:,2,2);
end
figure
subplot(2,1,1)
plot(t,y11)
ylabel('$\theta_1$','interpreter','latex')
title('Step response LPV controller (with scaling, X(p) & Y(p))')
subplot(2,1,2)
plot(t,y12)
ylabel('$\theta_2$','interpreter','latex')
xlabel('Time (s)')
legend('$G(0)$','$G(\pi/2)$','$G(\pi)$','interpreter','latex',...
    'location','southeast','Orientation','horizontal');
legend boxoff


fprintf('\n--------------------------------------\n')
fprintf('Comparison with scalings and dX/dt=0\n')
fprintf('hinfsyn:                  %5.4f\n',ghinf)
fprintf('X & Y cte:                %5.4f\n',g0)
fprintf('X(p) & Y(p):              %5.4f\n',g1)
fprintf('--------------------------------------\n')




