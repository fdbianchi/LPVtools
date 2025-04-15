
% LPVtools: Basic example with c2 parameter dependent
%
%   Plant:          Polytopic LPV (ppss)
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
vert = [0 1 0;
        0 0 1];
pv = pset.Gral(vert);

% system matrix          
A(:,:,1) = [-2.0 1.0; -1.0 -2.0];
A(:,:,2) = [-1.3 0.5; -0.5 -1.3];
A(:,:,3) = [-3.0 1.0;  1.0 -3.0];
C(:,:,1) = [10 10];
C(:,:,2) = [10.2 10];
C(:,:,3) = [10 10.3];

% all these matrices are the same at eache vertices
B = [1.0; 1.0];
D = 0;

% LPV model
pdG = ppss(A,B,C,D,pv);
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
n = 5;
theta = linspace(0,pi/2,n);
p = [1-sin(theta);1-cos(theta)];
plot(pv)
plot(p(1,:),p(2,:),'r.')

K = ss(pdK,p);
Gau = ss(pdGau,p);
Gcl_lpv1 = lft(Gau,K);
Gcl_lpv = ss(lft(pdGau,pdK),p);

Kf = ss(pdKf,p)*F;
Gcl_lpvf1 = lft(Gau,Kf);
Gcl_lpvf = ss(lft(pdGau,pdKf*F),p);

% checking the design at frozen parameter values
Gaw = ss(pdGaw,p);
for ii = 1:n
   [Klti(:,:,ii),constOut] = lpvsyn(Gaw(:,:,ii),3,2,const);
   glti(ii) = constOut(1).bound;
   [Khinf(:,:,ii),~,ghinf(ii),INFO] = hinfsyn(Gaw(:,:,ii),1,1);
end
Gcl_lti  = lft(Gau,Klti);
Gcl_hinf = lft(Gau,Khinf);

figure
step(Gcl_lpv,Gcl_lpvf,Gcl_lti,Gcl_hinf,1);
legend('LPV controller', 'LPV controller filtered', 'LTI controllers', 'Hinf controllers')

fprintf('\n')
fprintf('-------------------------------------------------------\n')
fprintf('Comparison with local Hinf synthesis\n')
fprintf('LPV:                   %5.4f\n',glpv)
fprintf('LPV (filtered):        %5.4f\n',glpvf)
for ii = 1:n
    fprintf('LTI @ p=[%3.2f,%3.2f]:   %5.4f (hinfsyn -> %5.4f)\n',...
        [p(:,ii);glti(ii);ghinf(ii)]);
end
fprintf('---------------------------------------------------------------------\n')

