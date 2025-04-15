
% LPVtools: Basic example with b2 parameter dependent
%
%   Plant:          general LPV (pgss)
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   Constant

% fbianchi - 2025-02-11


% cleaning
clearvars; 
clc; 
close all


% ========================================================================
% Modelling

% parameter set
range = [0 1];
pv = pset.Grid(range,5);

% system matrix          
A(:,:,1) = [-2 1; -1 -2];
A(:,:,2) = [-1.3 0.1; -0.1 -1.3];
A(:,:,3) = [-0.5 0; 0 -0.5];
B(:,:,1) = [1;1];
B(:,:,2) = [0.1;0];
B(:,:,3) = [0;0.2];
% the rest of matrices are parameter independent
C = [10 10];
D = 0;

% functions
fcn =@(p) [sin(p); cos(p)];

% LPV model
pdG = pgss(A,B,C,D,pv,fcn);
pdG.u = 'u';    pdG.y = 'y';

% ========================================================================
% Control design

% Augmented plant
sb    = sumblk('e = r - y');
pdGau = connect(pdG,sb,{'r','u'},{'e','u','e'});
%
% weigths
W1 = tf(10,[1 0.01]);
W2 = tf([0.04 0.1],[0.004 1]);
Wout = append(W1,W2,1);
% augmented plant + weigths
pdGaw = Wout*pdGau;
% filtering the control for classical design
F = tf(1,[0.01 1]);
Win = append(1,F);
pdGawF = pdGaw*Win;

% constraints
const(1) = synConst.Gain(1,1:2);
const(2) = synConst.Poles('MaxFreq',1000);

% synthesis
opt = lpvsettings('solver','mosek');
% synthesis
ctrlfcn.pdIn = 1;
ctrlfcn.pdOut = 0;
[pdK,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn,opt);
glpv = constOut(1).bound;

% classical design
[pdKf,constOut] = lpvsyn(pdGawF,3,2,const,[],[],opt);
glpvf = constOut(1).bound;


% ========================================================================
% Evaluation/Validation

p = linspace(0,1,5);

K   = ss(pdK,p);
Gau = ss(pdGau,p);
Gcl_lpv = lft(Gau,K);

Kf = ss(pdKf)*F;
Gcl_lpvf = lft(Gau,Kf);

% checking the design at frozen parameter values
[~,nv] = size(p);
Gaw = ss(pdGaw,p);
for ii = 1:nv
   [Klti(:,:,ii),constOut] = lpvsyn(Gaw(:,:,ii),3,2,const);
   glti(ii) = constOut(1).bound;
   [Khinf(:,:,ii),~,ghinf(ii),INFO] = hinfsyn(Gaw(:,:,ii),1,1);
end
Gcl_lti = lft(Gau,Klti);
Gcl_hinf = lft(Gau,Khinf);

step(Gcl_lpv,Gcl_lpvf,Gcl_lti,Gcl_hinf,1);
legend('LPV controller', 'LPV controller filtered', 'LTI controllers', 'Hinf controllers')

fprintf('\n')
fprintf('---------------------------------------------------------------------\n')
fprintf('Comparison\n')
fprintf('LPV:                   %5.4f\n',glpv)
fprintf('LPV (filtered):        %5.4f\n',glpvf)
for ii = 1:nv
    fprintf('LTI @ p=[%3.2f]:   %5.4f (hinfsyn -> %5.4f)\n',...
        [p(:,ii);glti(ii);ghinf(ii)]);
end
fprintf('---------------------------------------------------------------------\n')

