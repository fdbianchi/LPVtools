
% LPVtools: Basic example
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

% synthesis
[pdK,constOut] = lpvsyn(pdGaw,3,2,const);
glpv = constOut(1).bound;

% ========================================================================
% Evaluation/Validation

K   = ss(pdK);
Gau = ss(pdGau);
Gcl_lpv = lft(Gau,K);

% checking the design at frozen parameter values
[~,nv] = size(pv);
Gaw = ss(pdGaw);
for ii = 1:nv
   [Klti(:,:,ii),constOut] = lpvsyn(Gaw(:,:,ii),3,2,const);
   glti(ii) = constOut(1).bound;
   [Khinf(:,:,ii),~,ghinf(ii),INFO] = hinfsyn(Gaw(:,:,ii),1,1);
end
Gcl_lti  = lft(Gau,Klti);
Gcl_hinf = lft(Gau,Khinf);

step(Gcl_lpv,Gcl_lti,Gcl_hinf,1);
legend('LPV controller', 'LTI controllers', 'Hinf controllers')

% checking with analysis function
pdGcl = lft(pdGaw, pdK);
[constA,objA] = lpvanalysis(pdGcl,const);
glpvA = constA(1).bound;

% ========================================================================
% parameter dependent Lyapunov function X(p) = X0 + p(1)X1 + p(2)X2


% synthesis
% - using multiconvexity
[pdKp1,constOut] = lpvsyn(pdGaw,3,2,const,'affX');
glpvp1 = constOut(1).bound;

% - using gridding
lyapSet = createLyapFcn(@(p) p,0);
[pdKp2,constOut] = lpvsyn(pdGaw,3,2,const,lyapSet);
glpvp2 = constOut(1).bound;

% analysis 
% denser grid
points = pgrid(range,4);
pvTest = pset.Gral(points);
% updated plant
pdGawp = pdGaw;
pdGawp.parset = pvTest;
% affine Lyapunov function
lyapSetA = createLyapFcn('cl',@(p) p,0);
pdGclp1 = ppss(lft(pdGawp, pdKp1),pvTest);
constA1 = lpvanalysis(pdGclp1, const, lyapSetA);
glpvp1_A = constA1(1).bound;

pdGclp2 = lft(pdGawp, pdKp2);
pdGclp2 = ppss(lft(pdGawp, pdKp2),pvTest);
constA2 = lpvanalysis(pdGclp2, const, lyapSetA);
glpvp2_A = constA2(1).bound;


fprintf('\n')
fprintf('-------------------------------------------------------\n')
fprintf('Comparison with local Hinf synthesis\n')
fprintf('LPV:             %5.4f (analysis: %5.4f)\n',glpv,glpvA)
for ii = 1:nv
    fprintf('LTI @ p=[%1.0f,%1.0f]:   %5.4f (hinfsyn -> %5.4f)\n',...
        [pv.points(:,ii);glti(ii);ghinf(ii)]);
end
fprintf('\n')
fprintf('-------------------------------------------------------\n')
fprintf('Comparison using parameter-dependent Lyapunov functions\n')
fprintf('LPV (cte)        %5.4f (analysis: %5.4f)\n',glpv,glpvA)
fprintf('LPV (X aff)      %5.4f (analysis: %5.4f)\n',glpvp1,glpvp1_A)
fprintf('LPV (fcnX(p)=p)  %5.4f (analysis: %5.4f)\n',glpvp2,glpvp2_A)
fprintf('-------------------------------------------------------\n')


% ========================================================================
% Testing gridding
% 
% % constraints
% const(1) = synConst.Gain(1,1:2);
% const(2) = synConst.Poles('MaxFreq',1000);
% 
% [pdK,constOut,obj,synSet] = lpvsyn(pdGaw,3,2,const);
% pvDenser = pset.Grid(range,2*nv);
% [bool,constraint,obj] = lpvsynCheck(synSet,pvDenser);
% 

