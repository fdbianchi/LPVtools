
% LPVtools: Basic example - state feedback
%
%   Plant:          affine LPV (pass)
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   Constant
%   Ctrller Fcn:    No specified

% fbianchi - 2021-03-27


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
C(:,:,2) = [0 0];
C(:,:,3) = [0 0];

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

% constraints
const(1) = synConst.Gain(1,1:2);
const(2) = synConst.Poles('MaxFreq',1000);

% synthesis: output feedback
[pdK_of,constOut] = lpvsyn(pdGaw,3,2,const);
g_of = constOut(1).bound;

% synthesis: state feedback
[pdK_sf,constOut] = lpvsyn(pdGaw,0,2,const);
g_sf = constOut(1).bound;

% checking with analysis function
pdGcl_sf = lft(pdGaw,pdK_sf);
[constA,objA] = lpvanalysis(pdGcl_sf,const);
g_sf_A = constA(1).bound;

% synthesis: state feedback robust
[K_rob,constOut] = lpvsyn(pdGaw,0,2,const,[],0);
g_sf_rob = constOut(1).bound;

% checking with analysis function
pdGcl_rob = lft(pdGaw,K_rob);
[constA,objA] = lpvanalysis(pdGcl_rob,const);
g_sf_rob_A = constA(1).bound;

% ========================================================================
% Evaluation/Validation

Kof = ss(pdK_of);
Gau = ss(pdGau);
Gcl_of = lft(Gau,Kof);

% state feedback interconnection
A = pdGaw.A;
B = pdGaw.B;
C = [-10 -10 0 0;
       0   0 0 0; eye(4)];
D = [  1  0;
       0  1;zeros(4,2)];
pdGaw_sf = pass(A,B,C,D,pv); 
Gaw_sf   = ss(pdGaw_sf);
K_sf     = ss(pdK_sf);
Gcl_sf   = lft(Gaw_sf,K_sf);
K_rob    = repmat(K_rob,[1 1 4]);
Gcl_rob  = lft(Gaw_sf,K_rob);

figure
step(Gcl_of,Gcl_sf,Gcl_rob,1);
legend('output feedback controller', 'state feedback controllers', 'robust state feedback controllers')

% checking the design at frozen parameter values
[~,nv] = size(pv);
Gaw = ss(pdGaw);
for ii = 1:nv
   [Klti(:,:,ii),constOut] = lpvsyn(Gaw(:,:,ii),0,2,const);
   glti(ii) = constOut(1).bound;
   [~,~,~,INFO] = hinfsyn(Gaw(:,:,ii),1,1);
   Khinf(:,:,ii) = INFO.KFI;
   ghinf(ii) = INFO.GAMFI;
end
Gcl_lti  = lft(Gau,Klti);
Gcl_hinf = lft(Gau,Khinf);

fprintf('\n')
fprintf('-------------------------------------------------\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('LPV-Output feedback:    %5.4f\n',g_of)
fprintf('LPV-State feedback:     %5.4f (analysis: %5.4f)\n',g_sf,g_sf_A)
fprintf('Robust-State feedback:  %5.4f (analysis: %5.4f)\n',g_sf_rob,g_sf_rob_A)
fprintf('-------------------------------------------------\n')
fprintf('Comparison at frozen parameter values\n')
for ii = 1:nv
    fprintf('LTI @ p=[%1.0f,%1.0f]:   %5.4f (hinfsyn -> %5.4f)\n',...
        [pv.points(:,ii);glti(ii);ghinf(ii)]);
end
fprintf('-------------------------------------------------\n')

Gcl_lti  = lft(Gaw_sf,Klti);
Gcl_hinf = lft(Gaw_sf,Khinf(1,1:4,:));

figure
step(Gcl_sf,Gcl_lti,Gcl_hinf,1);
legend('LPV controller', 'LTI controllers', 'Hinf controllers')
