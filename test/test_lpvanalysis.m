% =========================================================================
%
% Testing basic use of lpvanalysis
%
% fbianchi - 2024-02-02
%
% =========================================================================

% cleaning
clearvars
close all
clc

fprintf('\nTesting lpvanalysis\n')

% ========================================================================
% Model from 
% F. Wu, X. H. Yang, A. Packard, and G. S. Becker, ‘Induced L2-norm control 
% for LPV systems with bounded parameter variation rates’, Int. J. of Robust 
% and Nonlinear Control, vol. 6, no. 9–10, pp. 983–998, Nov. 1996

% grid of points
nv    = 15;
range = [-pi,pi];
rate  = 2*[-1 1];
grd   = pgrid(range,nv);
% parameter set
pv = pset.Grid(range,nv,rate,'rho');

% LPV model: described as affine function on functions of the parameters
%
A11 = [0.75   2.00;
       0      0.50];
tau = 10;

% Independent term:
A0 = [A11 zeros(2);
      zeros(2) -tau*eye(2)];
B  = [zeros(2); tau*eye(2)];
C  = [eye(2) zeros(2)];
D  = zeros(2,2);
% term depending on sen(p)
A1 = zeros(4); A1(1,4) = 1; A1(2,3) = -1;
% term depending on cos(p)
A2 = zeros(4); A2(1,3) = 1; A2(2,4) = 1;
A  = cat(3,A0,A1,A2);

% state feedback
dk = zeros(2,4,3);
dk(:,:,3) = [-eye(2)*(eye(2)+A11) zeros(2)];
dk(:,:,2) = [-[0 -1;1 0]*(eye(2)+A11) zeros(2)];

Acl = zeros(4,4,3);
for ii = 1:3
    Acl(:,:,ii) = A(:,:,ii) + B*dk(:,:,ii);
end

% parameter function
func  = @(p) [sin(p);cos(p)];

% lpv model
pdG = pgss(Acl,B,C,D,pv,func,...
    'InputName','u',...
    'OutputName','y');

fprintf('\nPlant\n')
disp(pdG)

% Eigenvalues
eig(pdG);

% graphic functions
figure
sigma(pdG);

% parameter dependant Lyapunov functions,
lyapSetA = createLyapFcn('cl',@(p) [sin(p);cos(p)],@(p) [cos(p);-sin(p)]);
lyapSetA.inv = false;
const(1) = synConst.Gain(1:2,1:2);
% const(2) = synConst.Poles('MinDecay',0.19,'MaxFreq',1500);

opts = lpvsettings('solver','mosek','verb',true);
% opts = lpvsettings('solver','sedumi');
rateSet = linspace(0,15,11);
for ii = 1:11
    rate = rateSet(ii)*[-1 1];
    pdG.parset = pset.Grid(range,nv,rate,'rho');
    constOutA = lpvanalysis(pdG,const,lyapSetA,[],opts);
    ghinf(ii) = constOutA(1).bound;
end

figure
idx = isinf(ghinf);
plot(rateSet(idx),0*rateSet(idx),'bx')
hold on
plot(rateSet(~idx),ghinf(~idx),'rs')
Ha = gca;
Ha.XTick = rateSet;
xlabel('Parameter variation rate (rad/sec)')
ylabel('gamma')
title('Performance vs parameter rate')

