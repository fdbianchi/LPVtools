
% LPVtools: Basic example b2 and c3 parameter dependent
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
C(:,:,1) = [10 10];
C(:,:,2) = [0.2 0];
C(:,:,3) = [0 0.3];
% the rest of matrices are parameter independent
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
% filtering control and measure for classical design
F = tf(1,[0.01 1]);
Woutf = append(1,1,F);
Win = append(1,F);
pdGawF = Woutf*pdGaw*Win;

% constraints
const(1) = synConst.Gain(1,1:2);
const(2) = synConst.Poles('MaxFreq',1000);

% synthesis
opt = lpvsettings('solver','mosek');
% synthesis
ctrlfcn.pdIn = 0;
ctrlfcn.pdOut = 0;
ctrlfcn.dk = 1;
[pdK1,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn,opt);
glpv1 = constOut(1).bound;

ctrlfcn.pdIn = 1;
ctrlfcn.pdOut = 0;
[pdK2,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn,opt);
glpv2 = constOut(1).bound;

ctrlfcn.pdIn = 0;
ctrlfcn.pdOut = 1;
[pdK3,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn,opt);
glpv3 = constOut(1).bound;

ctrlfcn.pdIn = 1;
ctrlfcn.pdOut = 1;
ctrlfcn.dk = 0;
[pdK4,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn,opt);
glpv4 = constOut(1).bound;

% classical design
[pdKf,constOut] = lpvsyn(pdGawF,3,2,const,[],[],opt);
glpvf = constOut(1).bound;


% ========================================================================
% Evaluation/Validation

p = linspace(0,1,5);

Gau = ss(pdGau,p);
K1 = ss(pdK1,p);
Gcl_lpv1 = ss(lft(pdGau,pdK1),p);
K2 = ss(pdK2,p);
Gcl_lpv2 = ss(lft(pdGau,pdK2),p);
K3 = ss(pdK3,p);
Gcl_lpv3 = ss(lft(pdGau,pdK3),p);
K4 = ss(pdK4,p);
Gcl_lpv4 = ss(lft(pdGau,pdK4),p);

Kf = F*ss(pdKf)*F;
Gcl_lpvf = ss(lft(pdGau,F*pdKf*F),p);

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

step(Gcl_lpv1,Gcl_lpv2,Gcl_lpv3,Gcl_lpv4,Gcl_lpvf,Gcl_lti,Gcl_hinf,1);
legend('LPV ctrller (pdIn=false & pdOut=false)',...
       'LPV ctrller (pdIn=true & pdOut=false)',...
       'LPV ctrller (pdIn=false & pdOut=true)',...
       'LPV ctrller (pdIn=true & pdOut=true, dk=0)',...
       'LPV controller filtered', 'LTI controllers', 'Hinf controllers')

fprintf('\n')
fprintf('---------------------------------------------------------------------\n')
fprintf('Comparison\n')
fprintf('LPV (pdIn=false & pdOut=false): %5.4f\n',glpv1)
fprintf('LPV (pdIn=true & pdOut=false):  %5.4f\n',glpv2)
fprintf('LPV (pdIn=false & pdOut=true):  %5.4f\n',glpv3)
fprintf('LPV (pdIn=true & pdOut=true):   %5.4f (dk=0)\n',glpv4)
fprintf('LPV (filtered):                 %5.4f\n',glpvf)
for ii = 1:nv
    fprintf('LTI @ p=[%3.2f]:   %5.4f (hinfsyn -> %5.4f)\n',...
        [p(:,ii);glti(ii);ghinf(ii)]);
end
fprintf('---------------------------------------------------------------------\n')

