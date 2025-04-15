
% LPVtools: Example of multi-objective synthesis (LTI)
%
% Scherer, C. W., Gahinet, P., & Chilali, M. (1997). Multiobjective output-
% feedback control via LMI optimization. IEEE Transactions on Automatic 
% Control, 42(7), 896–911. doi:10.1109/9.599969 (Example 7)
%
%   Plant:          LTI (ppss)
%   Constraint:     Hinf/H2
%   Lyapunov Fcn:   Constant
%   Ctrller Fcn:    No specified

% fbianchi - 2021-05-29

% cleanig 
clc
clearvars
close all

% ========================================================================
% Modelling

% augmented plant
A = [0  10   2; 
    -1   1   0; 
     0   2  -5];
B = [1   0; 
     0   1;
     1   0];
C = [1   0   0; 
     0   0   0; 
     0   1   0; 
     0   0   1; 
     0   0   0; 
     0   1   0];
D = [0   0; 
     0   1; 
     0   0; 
     0   0; 
     0   1; 
     2   0];

G = ss(A,B,C,D);
G.u = {'w','u'};
G.y = {'zi1','zi2','z21','z22','z23','y'};
pdG = ppss(G);


% ========================================================================
% Control design

% 1) designing only with H2 constraints
G2 = G(3:6,:);
[Kcs_h2,CL,gcs_h2,INFO] = h2syn(G2,1,1);
Gcl = lft(G,Kcs_h2);
gcs_hinf = norm(Gcl,inf);

% using LPVtools 
const(1) = synConst.GainH2(1,[3 4 5]);
[Klpv_h2,constOut] = lpvsyn(G,6,2,const);
glpv_h2 = constOut(1).bound;

Gcl = lft(G,Klpv_h2);
glpv_hinf = norm(Gcl,inf);

[constA,objA] = lpvanalysis(Gcl,const);
glpv_h2A = constA(1).bound;

% using LPVtools - H2 generalized
const(1) = synConst.GainH2g(1,[3 4 5]);
[Klpv_h2g,constOut] = lpvsyn(G,6,2,const);
glpv_h2g = constOut(1).bound;

Gclg = lft(G,Klpv_h2g);
glpv_hinfg = norm(Gclg,inf);

[constA,objA] = lpvanalysis(Gclg,const);
glpv_h2gA = constA(1).bound;

fprintf('\n----------------------------------------------------\n')
fprintf('Comparison: min H2\n')
fprintf('h2hinfsyn:  %6.4f (norm_hinf= %6.4f)\n',gcs_h2,gcs_hinf)
fprintf('lpvsyn:     %6.4f (norm_hinf= %6.4f)\n',glpv_h2,glpv_hinf)
fprintf('lpvsyn:     %6.4f               (analysis)\n',glpv_h2A)
fprintf('lpvsyn:     %6.4f - generalized (norm_hinf= %6.4f)\n',glpv_h2g,glpv_hinfg)
fprintf('lpvsyn:     %6.4f               (analysis)\n',glpv_h2gA)
fprintf('\n----------------------------------------------------\n')


% 2) minimizing H2 s.t. to Hinf constraints

hinfmax = 23.6;

% multiobjective from control system toolbox
[Kcs_h2,~,gcs_h2] = h2hinfsyn(G,1,1,3,[0 1],'HINFMAX',hinfmax);

% using LPVtools
const(1) = synConst.Gain(1,[1 2],'bound',hinfmax);
const(2) = synConst.GainH2(1,[3 4 5]);
[Klpv_h2,constOut] = lpvsyn(G,6,2,const);
glpv_h2 = constOut(2).bound;

% using LPVtools - H2 generalized
const(1) = synConst.Gain(1,[1 2],'bound',hinfmax);
const(2) = synConst.GainH2g(1,[3 4 5]);
[Klpv_h2g,constOut] = lpvsyn(G,6,2,const);
glpv_h2g = constOut(2).bound;

fprintf('\n-----------------------------------\n')
fprintf('Comparison: min H2, s.t Hinf=23.6\n')
fprintf('h2hinfsyn:  %6.4f\n',gcs_h2(2))
fprintf('lpvsyn:     %6.4f\n',glpv_h2)
fprintf('lpvsyn:     %6.4f - generalized\n',glpv_h2g)
fprintf('-----------------------------------\n')



% ========================================================================
% simulations

% H2 generalized: ||w||_2 -> peak
% pulse with ||w||_2 = 1
h = 0.01;
A = 1/sqrt(h);

sim('MultiObj_LTI_01_sim.slx');

t = output.time;
u = output.signals(1).values;
y = output.signals(2).values;
% inf-norm
n2g = max(sqrt(dot(y',y')));

figure
subplot(2,1,1)
plot(t,u)
ylabel('u')
xlabel('time (s)')
title(sprintf('Example of Generalized H2 norm; ||y||_\\infty = %5.2f',n2g))

subplot(2,1,2)
plot(t,y)
ylabel('y')
xlabel('time (s)')


% H2: impulse -> ||z||_2
[y,t] = impulse(Gcl(3:5),10);
Ts = mean(diff(t));
% 2-norm
n2g = sqrt(sum(dot(y',y')*Ts));

figure
plot(t,y)
ylabel('y')
xlabel('time (s)')
title(sprintf('Example of H2 norm; ||y||_2 = %5.2f',n2g))



