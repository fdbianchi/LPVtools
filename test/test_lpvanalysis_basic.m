
% =========================================================================
%
% Testing basic use of lpvanalysis - Basic
%
% fbianchi - 2024-12-02
%
% =========================================================================

% cleaning
clearvars
close all
clc

fprintf('\nTesting lpvanalysis\n')

%% -----------------------------------------------------------------------
% LTI

fprintf('\nCase LTI\n')

rng(23656);
G = rss(5,2,3);
G.D = zeros(2,3);

fprintf('\nEigenvalues:\n')
E1 = eig(G);
disp(E1)

const = synConst.Gain(1:3,1:2);
constOutA = lpvanalysis(G,const);
ghinf = constOutA(1).bound;

fprintf('\nH-infinity Performance LPVANALYSIS:\n')
disp(ghinf)

fprintf('Norm infinity:\n')
gnorm_inf = norm(G,inf);
disp(gnorm_inf)

const = synConst.GainH2(1:3,1:2);
constOutA = lpvanalysis(G,const);
gh2 = constOutA(1).bound;

fprintf('\nH2 Performance LPVANALYSIS:\n')
disp(gh2)

fprintf('Norm-2:\n')
gnorm2 = norm(G,2);
disp(gnorm2)


%% -----------------------------------------------------------------------
% LPV: 

% LPV description
A(:,:,1) = [ 0  1; 0  0];
A(:,:,2) = [-1  0; 0  0];
A(:,:,3) = [ 0  0;-1  0];
A(:,:,4) = [ 0  0;-1  0];
B(:,:,1) = [-1/5; 0];
B(:,:,2) = [0;    1];
C = eye(2);
D = zeros(2,1);

range = [3.0, 4.5;
         5.0, 6.5;
         9.0, 20.25];
pv = pset.Box(range);
pdG = pass(A,B,C,D,pv);

% Comparison models
figure
sigma(pdG)

fprintf('\nCase LPV\n')
disp(pdG)

% fprintf('\nEigenvalues:\n')
eig(pdG);

const = synConst.Gain(1,1:2);
constOutA = lpvanalysis(pdG,const);
ghinf = constOutA(1).bound;

% lmitool
psG = psys(pdG);
glmi = quadperf(psG);

fprintf('\nH-infinity Performance LPVANALYSIS:\n')
disp(ghinf)

fprintf('\nH-infinity Performance QUADPERF:\n')
disp(ghinf)


const = synConst.Gain(1,1:2);
constOutA = lpvanalysis(pdG,const,'aff');
ghinf2 = constOutA(1).bound;

fprintf('\nH-infinity Performance LPVANALYSIS,\n\twith parameter dependend Lyapunov function, dp/dt=0:\n')
disp(ghinf2)

fprintf('Norm infinity at each vertex:\n')
gnorm_inf = norm(ss(pdG),inf);
disp(gnorm_inf)





