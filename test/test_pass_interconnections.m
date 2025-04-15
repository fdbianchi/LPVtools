% -------------------------------------------------------------------------
% Testing interconnections for pass
%
% -------------------------------------------------------------------------

% fbianchi - 2024-12-09

% cleaning
clc
clear all
close all

%% 
fprintf('\n------------------------------------------------------------\n')
fprintf('PASS: Checking interconnections\n\n')


pv = pset.Box([0 1;2 6]);

% model with input and output channel parameter dependent
a0 = [ 0 1; 0 0]; a1 = [-1 0; 0 0]; a2 = [ 0 0;-1 0];
b0 = [ 0;1];      b1 = [ 0;0];      b2 = [ 1;0];
c0 = [-1 0;0 1];  c1 = [ 1 0; 0 0]; c2 = [ 0 0; 0 0];
d0 = [ 0;1];      d1 = [ 0;0];      d2 = [ 0;0];

a = cat(3,a0,a1,a2);
b = cat(3,b0,b1,b2);
c = cat(3,c0,c1,c2);
d = cat(3,d0,d1,d2);

pdG1 = pass(a,b,c,d,pv,'InputName','u1','OutputName','y1');
G1 = ss(pdG1);

dgn1 = ispd(pdG1);
str_dgn1 = {'false','false','false'}; str_dgn1(dgn1) = {'true'};
fprintf('model pdG1: pd A: %s, pd in: %s, pd out: %s\n',str_dgn1{:});

% model with input channel parameter dependent
a0 = [ 0 1; 0 0]; a1 = [-1 0; 0 0]; a2 = [ 0 0;-1 0];
b0 = [ 1 0;0 1];  b1 = [ 0 0;0 0];  b2 = [ 0 0; 0 0];
c0 = [-1 0;0 1];  c1 = [ 1 0; 0 0]; c2 = [ 0 0; 1 0];
d0 = [ 0 0;1 0];  d1 = [ 0 0; 0 0]; d2 = [ 0 0; 0 0];

a = cat(3,a0,a1,a2);
b = cat(3,b0,b1,b2);
c = cat(3,c0,c1,c2);
d = cat(3,d0,d1,d2);

pdG2 = pass(a,b,c,d,pv,'InputName','u2','OutputName','y2');
G2 = ss(pdG2);

dgn2 = ispd(pdG2);
str_dgn2 = {'false','false','false'}; str_dgn2(dgn2) = {'true'};
fprintf('model pdG2: pd A: %s, pd in: %s, pd out: %s\n',str_dgn2{:});

% model with output channel parameter dependent
pdG3 = pdG1(2,1);
pdG3.u = 'u3'; pdG3.y = 'y';
G3 = ss(pdG3);

dgn3 = ispd(pdG3);
str_dgn3 = {'false','false','false'}; str_dgn3(dgn3) = {'true'};
fprintf('model pdG3: pd A: %s, pd in: %s, pd out: %s\n',str_dgn3{:});

fprintf('\n')

%% ====================================================================
% 1) Series

%    LTI & LPV

G4 = rss(3,3,2);
G5 = rss(2,1,2);
G6 = diag([5, 2]);

pdGs1 = G4*pdG1;
s1 = G4*G1;

figure
bodemag(pdGs1,'r',s1,'kx')
legend('Using LPV tools','Using LTI tools')
title('Series interconnection, G4*pdG1')

pdGs2 = pdG1*G5;
s2 = G1*G5;

figure
bodemag(pdGs2,'r',s2,'kx')
legend('Using LPV tools','Using LTI tools')
title('Series interconnection, pdG1*G5')

pdGs3 = G6*pdG1;
s3 = G6*G1;

figure
bodemag(pdGs3,'r',s3,'kx')
legend('Using LPV tools','Using LTI tools')
title('Series interconnection, diag([5, 2])*pdG1')


%    LPV & LPV

pdGp1 = pdG2*pdG1;
p1 = G2*G1;

figure
bodemag(pdGp1,'r',p1,'kx')
legend('Using LPV tools','Using LTI tools')
title('Series interconnection, pdG2*pdG1')

pdGp2 = pdG1*pdG3;
p2 = G1*G3;

figure
bodemag(pdGp2,'r',p2,'kx')
legend('Using LPV tools','Using LTI tools')
title('Series interconnection, pdG1*pdG3')

pdGp3 = series(pdG1,pdG2,1,2);
p3 = series(G1,G2,1,2);

figure
bodemag(pdGp3,'r',p3,'kx')
legend('Using LPV tools','Using LTI tools')
title('Series interconnection, pdG2(:,2)*pdG1(1,:)')

% ====================================================================
% 2) Parallel

%    LTI & LPV

G5 = rss(3,2,1);
G6 = [3;5];

pdGs1 = G5 + pdG1;
s1 = G5 + G1;

figure
bodemag(pdGs1,'r',s1,'kx')
legend('Using LPV tools','Using LTI tools')
title('Parallel interconnection, G4 + pdG1')

pdGs2 = pdG1 + G6;
s2 = G1 + G6;

figure
bodemag(pdGs2,'r',s2,'kx')
legend('Using LPV tools','Using LTI tools')
title('Parallel interconnection, pdG1 + [3;5]')


%    LPV & LPV
pdGp1 = parallel(pdG1,pdG2,1,2,2,1);
p1 = parallel(G1,G2,1,2,2,1);

figure
bodemag(pdGp1,'r',p1,'kx')
legend('Using LPV tools','Using LTI tools')
title('Parallel interconnection, parallel(pdG1,pdG2,1,2,2,1)')

pdGp2 = pdG1(1,1) + pdG3;
p2 = G1(1,1,:) + G3;

figure
bodemag(pdGp2,'r',p2,'kx')
legend('Using LPV tools','Using LTI tools')
title('Parallel interconnection, pdG1(1,1) + pdG3')


% ====================================================================
% 3) Append

%    LTI & LPV

G5 = rss(3,2,1);
G6 = [3;5];

pdGs1 = append(G5,pdG1);
s1 = append(G5,G1);

figure
bodemag(pdGs1,'r',s1,'kx')
legend('Using LPV tools','Using LTI tools')
title('Append interconnection, append(G4,pdG1)')

pdGs2 = append(pdG1,G6);
s2 = append(G1,G6);

figure
bodemag(pdGs2,'r',s2,'kx')
legend('Using LPV tools','Using LTI tools')
title('Append interconnection, append(pdG1,[3;5])')


%    LPV & LPV

pdGp1 = append(pdG2,pdG1);
p1 = append(G2,G1);

figure
bodemag(pdGp1,'r',p1,'kx')
legend('Using LPV tools','Using LTI tools')
title('Append interconnection, append(pdG2,pdG1)')

% ====================================================================
% 4) cat

% 4.1) horzcat

%    LTI & LPV

pdGs1 = horzcat(G5,pdG1);
s1 = horzcat(G5,G1);

figure
bodemag(pdGs1,'r',s1,'kx')
legend('Using LPV tools','Using LTI tools')
title('horizontally concatenation, horzcat(G5,pdG1)')

pdGs2 = horzcat(pdG1,G6);
s2 = horzcat(G1,G6);

figure
bodemag(pdGs2,'r',s2,'kx')
legend('Using LPV tools','Using LTI tools')
title('horizontally concatenation, horzcat(pdG1,[3;5])')

%    LPV & LPV

pdGp1 = horzcat(pdG2,pdG1);
p1 = horzcat(G2,G1);

figure
bodemag(pdGp1,'r',p1,'kx')
legend('Using LPV tools','Using LTI tools')
title('horizontally concatenation, horzcat(pdG2,pdG1)')

% 4.2) vertcat

%    LTI & LPV

pdGs1 = vertcat(G5,pdG1);
s1 = vertcat(G5,G1);

figure
bodemag(pdGs1,'r',s1,'kx')
legend('Using LPV tools','Using LTI tools')
title('Vertically concatenation, vertcat(G5,pdG1)')

pdGs2 = vertcat(pdG1,G6);
s2 = vertcat(G1,G6);

figure
bodemag(pdGs2,'r',s2,'kx')
legend('Using LPV tools','Using LTI tools')
title('Vertically concatenation, vertcat(pdG1,[3;5])')

%    LPV & LPV

pdGp1 = vertcat(pdG2(:,1),pdG1);
p1 = vertcat(G2(:,1),G1);

figure
bodemag(pdGp1,'r',p1,'kx')
legend('Using LPV tools','Using LTI tools')
title('Vertically concatenation, vertcat(pdG2(:,1),pdG1)')
