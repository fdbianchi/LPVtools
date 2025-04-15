
% -------------------------------------------------------------------------
% How to use and Examples for synConst.Gain
%
% -------------------------------------------------------------------------

% fbianchi - 2020-03-31

clc
clear all

fprintf('\nTesting synConst.Gain\n\n')

% empty constraint
g0 = synConst.Gain;
fprintf('Empty object\n')
disp(g0)

% minimum syntaxis
g1a = synConst.Gain([1 3],[2 3]);
fprintf('\nInput-output map with indices\n')
disp(g1a)
g1b = synConst.Gain({'e','y'},{'u','w'});
fprintf('\nInput-output map with names\n')
disp(g1b)

% copy
g2 = g1b;
fprintf('\nCopying object\n')
disp(g2)

% more complex syntax
g3a = synConst.Gain([1 3],[2 3],'factor',0.2);
fprintf('\nInput-output map with indices, and factor\n')
disp(g3a)
g3b = synConst.Gain([1 3],[2 3],'bound',3);
fprintf('\nInput-output map with indices, and bound\n')
disp(g3b)


% g3c = synConst.Gain([1 3],[2 3],'Win',2) 
% g3d = synConst.Gain([1 3],[2 3],'Win',tf(1,[1 1])) 
% win = append(tf(1,[2 2]),2);
% win = {tf(1,[2 2]),2};
% % win = append(tf(1,[2 2]),2,3);
% g3e = synConst.Gain([1 3],[2 3],'Win',win) 
% wout = append(tf(1,[2 2]),2);
% % wout = {tf(1,[2 2]),2};
% % wout = append(tf(1,[2 2]),2,3);
% g3f = synConst.Gain([1 3],[2 3],'Wout',wout) 


fprintf('\nTesting synConst.GainH2\n\n')

% empty constraint
g0 = synConst.GainH2;
fprintf('Empty object\n')
disp(g0)

% minimum syntaxis
g1a = synConst.GainH2([1 3],[2 3]);
fprintf('\nInput-output map with indices\n')
disp(g1a)
g1b = synConst.GainH2({'e','y'},{'u','w'});
fprintf('\nInput-output map with names\n')
disp(g1b)

% copy
g2 = g1b;
fprintf('\nCopying object\n')
disp(g2)

% more complex syntax
g3a = synConst.GainH2([1 3],[2 3],'factor',0.2);
fprintf('\nInput-output map with indices, and factor\n')
disp(g3a)
g3b = synConst.GainH2([1 3],[2 3],'bound',3);
fprintf('\nInput-output map with indices, and bound\n')
disp(g3b)


fprintf('\nTesting synConst.GainH2g\n\n')

% empty constraint
g0 = synConst.GainH2g;
fprintf('Empty object\n')
disp(g0)

% minimum syntaxis
g1a = synConst.GainH2g([1 3],[2 3]);
fprintf('\nInput-output map with indices\n')
disp(g1a)
g1b = synConst.GainH2g({'e','y'},{'u','w'});
fprintf('\nInput-output map with names\n')
disp(g1b)

% copy
g2 = g1b;
fprintf('\nCopying object\n')
disp(g2)

% more complex syntax
g3a = synConst.GainH2g([1 3],[2 3],'factor',0.2);
fprintf('\nInput-output map with indices, and factor\n')
disp(g3a)
g3b = synConst.GainH2g([1 3],[2 3],'bound',3);
fprintf('\nInput-output map with indices, and bound\n')
disp(g3b)

