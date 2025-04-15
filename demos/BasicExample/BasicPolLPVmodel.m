
% LPVtools: Basic example
%
%   Plant:          Polytopic LPV (ppss)
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   Constant
%   Ctrller Fcn:    No specified

% fbianchi - 2021-04-02


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

% all these matrices are the same at eache vertices
B = [1.0; 1.0];
C = [10 10];
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

% constraints
const(1) = synConst.Gain(1,1:2);
const(2) = synConst.Poles('MaxFreq',1000);

% synthesis
[pdK,constOut] = lpvsyn(pdGaw,3,2,const);
glpv = constOut(1).bound;

% checking with analysis function
pdGcl = lft(pdGaw,pdK);
[constA,objA] = lpvanalysis(pdGcl,const);
glpvA = constA(1).bound;

% affine Lyapunov function
[constA2,objA2,X] = lpvanalysis(pdGcl,const,'pwa');
glpvA2 = constA2(1).bound;

% ========================================================================
% Evaluation/Validation
n = 5;
theta = linspace(0,pi/2,n);
p = [1-sin(theta);1-cos(theta)];
plot(pv)
plot(p(1,:),p(2,:),'r.')

K   = ss(pdK,p);
Gau = ss(pdGau,p);
Gaw = ss(pdGaw,p);
Gcl_lpv = lft(Gau,K);

% checking the design at frozen parameter values
for ii = 1:n
   [Klti(:,:,ii),constOut] = lpvsyn(Gaw(:,:,ii),3,2,const);
   glti(ii) = constOut(1).bound;
   [Khinf(:,:,ii),~,ghinf(ii),INFO] = hinfsyn(Gaw(:,:,ii),1,1);
end
Gcl_lti  = lft(Gau,Klti);
Gcl_hinf = lft(Gau,Khinf);

figure
step(Gcl_lpv,Gcl_lti,Gcl_hinf,1);
legend('LPV controller', 'LTI controllers', 'Hinf controllers')

fprintf('\n')
fprintf('---------------------------------------------------------------------\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('LPV:                   %5.4f (analysis: %5.4f (X cte),  %5.4f (X PWA))\n',glpv,glpvA,glpvA2)
for ii = 1:n
    fprintf('LTI @ p=[%3.2f,%3.2f]:   %5.4f (hinfsyn -> %5.4f)\n',...
        [p(:,ii);glti(ii);ghinf(ii)]);
end
fprintf('---------------------------------------------------------------------\n')

