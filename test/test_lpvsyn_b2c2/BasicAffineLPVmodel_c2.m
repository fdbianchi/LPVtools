
% LPVtools: Basic example with c2 parameter dependent
%
%   Plant:          affine LPV (pass)
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
range = [ 0 1;
          0 1];
pv = pset.Box(range);

% system matrix          
A(:,:,1) = [-2 1; -1 -2];
A(:,:,2) = [-1.3 0.1; -0.1 -1.3];
A(:,:,3) = [-0.5 0; 0 -0.5];

B(:,:,1) = [1;1];
B(:,:,2) = [0;0];
B(:,:,3) = [0;0];

C(:,:,1) = [10 10];
C(:,:,2) = [0.2 0];
C(:,:,3) = [0 0.3];

D = zeros(1,1,3);

% LPV model
pdG = pass(A,B,C,D,pv);
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
% filtering for classical design
F = tf(1,[0.01 1]);
Woutf = append(1,1,F);
pdGawF = Woutf*pdGaw;

% constraints
const(1) = synConst.Gain(1,1:2);
const(2) = synConst.Poles('MaxFreq',1000);

opt = lpvsettings('solver','mosek');
% synthesis
ctrlfcn.pdIn = 0;
ctrlfcn.pdOut = 1;
[pdK,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn,opt);
glpv = constOut(1).bound;

% classical design
[pdKf,constOut] = lpvsyn(pdGawF,3,2,const,[],[],opt);
glpvf = constOut(1).bound;

% ========================================================================
% Evaluation/Validation

K   = ss(pdK);
Gau = ss(pdGau);
% Gcl_lpv = lft(Gau,K);
Gcl_lpv = lft(pdGau,pdK);

Kf = ss(pdKf)*F;
% Gcl_lpvf = lft(Gau,Kf);
Gcl_lpvf = lft(pdGau,pdKf*F);

% checking the design at frozen parameter values
[~,nv] = size(pv);
Gaw = ss(pdGaw);
for ii = 1:nv
   [Klti(:,:,ii),constOut] = lpvsyn(Gaw(:,:,ii),3,2,const,[],[],opt);
   glti(ii) = constOut(1).bound;
   [Khinf(:,:,ii),~,ghinf(ii),INFO] = hinfsyn(Gaw(:,:,ii),1,1);
end
Gcl_lti  = lft(Gau,Klti);
Gcl_hinf = lft(Gau,Khinf);

step(Gcl_lpv,Gcl_lpvf,Gcl_lti,Gcl_hinf,1);
legend('LPV controller', 'LPV controller filtered', 'LTI controllers', 'Hinf controllers')

fprintf('\n')
fprintf('-------------------------------------------------------\n')
fprintf('Comparison with local Hinf synthesis\n')
fprintf('LPV:                   %5.4f\n',glpv)
fprintf('LPV (filtered):        %5.4f\n',glpvf)
for ii = 1:nv
    fprintf('LTI @ p=[%3.2f,%3.2f]:   %5.4f (hinfsyn -> %5.4f)\n',...
        [pv.points(:,ii);glti(ii);ghinf(ii)]);
end
fprintf('---------------------------------------------------------------------\n')

