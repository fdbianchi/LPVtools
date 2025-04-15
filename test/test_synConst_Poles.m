
% -------------------------------------------------------------------------
% How to use and Examples for synConst.Poles
%
% -------------------------------------------------------------------------

% fbianchi - 2020-03-31

clc
clear all

fprintf('\nTesting synConst.Poles\n\n')

% empty
g0 = synConst.Poles;
fprintf('Empty object\n')
disp(g0)

% basic syntaxis
g1 = synConst.Poles('MinDecay',1);
fprintf('\nMinimun decay rate\n')
disp(g1)
g2 = synConst.Poles('MaxFreq',10);
fprintf('\nMaximun frequency\n')
disp(g2)
g3 = synConst.Poles('MaxFreq',10,'MinDamping',0.5);
fprintf('\nMaximun frequency and minimum damping\n')
disp(g3)

% concatenation of constraints
g4(1) = synConst.Gain([1 3],[2 3]);
g4(3) = g3;
fprintf('\nConstraint concatenation\n')
disp(g4)




