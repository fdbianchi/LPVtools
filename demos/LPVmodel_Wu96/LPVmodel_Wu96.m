
% LPVtools: Nonaffine LPV system, control using griding
% 
% Ref: Wu, F. et al "Induced L2-Norm Control for LPV Systems with Bounded
% Parameter Variation Rates", Int. Journal on robust and Nonlinear Control,
% 1996
%
%   Plant:          general LPV (pgss)
%                   1 parameter,
%                   2 inputs, 2 outputs
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   Constant/Parameter dependent
%   Ctrller Fcn:    No specified

% fbianchi - 2021-04-03

% cleaning
clearvars
close all
clc

% ========================================================================
% Modelling

% grid of points
nv    = 7;
range = [-pi,pi];
rate  = 0*[-1 1];
grd   = pgrid(range,nv);
% parameter set
pv = pset.Grid(range,nv,rate,'rho');

% LPV model: described as affine function on functions of the parameters
%
% Independent term:
A0 = [0.75   2.00    0      0;
      0      0.50    0      0;
      0      0      -10.00  0;
      0      0       0    -10.00];
B  = [0   0    0;
      3   0    0;
      0  10    0;
      0   0   10];
C  = [1   0    0    0;
      0   1    0    0];
D  = zeros(2,3);
% term depending on sen(p)
A1 = zeros(4); A1(1,4) = 1; A1(2,3) = -1;
% term depending on cos(p)
A2 = zeros(4); A2(1,3) = 1; A2(2,4) = 1;
A  = cat(3,A0,A1,A2);

% parameter function
func  = @(p) [sin(p);cos(p)];

% lpv model
pdG = pgss(A,B,C,D,pv,func,...
    'InputName',{'f','u(1)','u(2)'},...
    'OutputName','v');

% graphic functions
% figure
% bodemag(pdG);


% ========================================================================
% Control design

opt = lpvsettings('solver','sedumi');

% weighting functions
sb1 = sumblk('e = v - r',2);
Wp = ss(eye(2));                    Wp.u = 'e';  Wp.y = 'ep';
Wn = tf(10*[1 10],[1 1000])*eye(2); Wn.u = 'dn'; Wn.y = 'n';
sb2 = sumblk('y = v - r + n',2);
Wf = ss(1);                         Wf.u = 'df'; Wf.y = 'f';
Wu = ss(1/280*eye(2));              Wu.u = 'u';  Wu.y = 'eu';
Wr = tf(20,[1 0.2])*eye(2);         Wr.u = 'dr'; Wr.y = 'r';

% augmented plant
pdGau = connect(pdG,sb1,sb2,{'r','n','f','u'},{'e','u','y'});

Win  = append(Wr,Wn,Wf,1,1);    
Win.u = {'dr(1)','dr(2)','dn(1)','dn(2)','df','u(1)','u(2)'};
Wout = append(Wp,Wu,1,1);
Wout.y = {'ep(1)' 'ep(2)' 'eu(1)' 'eu(2)' 'y(1)' 'y(2)'};
pdGwe = Wout*pdGau*Win;

% H-infinity
y = sprintfc('y(%d)',1:2);
u = sprintfc('u(%d)',1:2);
z = {'dr(1)' 'dr(2)' 'dn(1)' 'dn(2)'};
w = {'ep(1)' 'ep(2)' 'eu(1)' 'eu(2)'};
const1 = synConst.Gain(z,w);

% Constant Lyapunov functions
[pdK,constOut] = lpvsyn(pdGwe,y,u,const1,[],[],opt);
ghinf_1 = constOut.bound;

% Parameter dependant Lyapunov function: X = X0; Y(p) = [sin(p);cos(p)]; 
%
% dY = 0
lyapSet = createLyapFcn(0,0,@(p) [sin(p);cos(p)],0);
[pdK0,constOut,~,synSet] = lpvsyn(pdGwe,y,u,const1,lyapSet,[],opt);
ghinf_2 = constOut.bound;
% checking in a denser grid
nvd = 2*nv;
pvDenser = pset.Grid(range,nvd,rate,'rho');
[bool,constOut,obj] = lpvsynCheck(synSet,pvDenser,opt);
ghinf_2C = constOut.bound;

% dY~=0
lyapSet = createLyapFcn(0,0,@(p) [sin(p);cos(p)],@(p) [cos(p);-sin(p)]);
ghinf_3 = zeros(1,11);
for ii = 1:11
    rate   = (ii-1)*[-1 1];
    pdGwe.parset = pset.Grid(range,nv,rate,'rho');
    [pdKi,constOut] = lpvsyn(pdGwe,y,u,const1,lyapSet,[],opt);
    if (ii == 6)
        pdK5 = pdKi;
    end
    ghinf_3(ii) = constOut.bound;
end

%%
% parameter dependant Lyapunov functions, dY~=0, and pole constraints
const2(1) = synConst.Gain(z,w);
const2(2) = synConst.Poles('MinDecay',0.19,'MaxFreq',1500);
ratef = 3*[-1 1];
pdGwe.parset = pset.Grid(range,nv,ratef,'rho');
[pdK4,constOut,~,synSet] = lpvsyn(pdGwe,y,u,const2,lyapSet,[],opt);
ghinf_4 = constOut(1).bound;
%
% checking in a denser grid
pvDenser = pset.Grid(range,nvd,ratef,'rho');
[chk1,constChk] = lpvsynCheck(synSet,pvDenser,opt);
ghinf_4C = constChk(1).bound;
%
% G = ss(pdGwe,pvDenser.points);
% K = ss(pdK4,pvDenser.points);
% Gcl = lft(G,K);
% E = eig(Gcl);
% for ii=1:length(pvDenser.points)
%     [wn,ep,poles] = damp(E(:,:,ii));
%     ReMax = max(real(poles));
%     fprintf('\tMinDecay = %6.4f\t MinDamping = %6.4f \tmaxFreq = %6.4f\n',...
%         -ReMax,min(ep),max(wn));
% end
%
% redesign with denser grid
pdGwe_new = pdGwe;
pdGwe_new.parset = pvDenser;
[pdK4n,constOut,~,synSet] = lpvsyn(pdGwe_new,y,u,const2,lyapSet,[],opt);
ghinf_4n = constOut(1).bound;
%
% checking in a denser grid
nvdd = 3*nv;
pvDenser = pset.Grid(range,nvdd,ratef,'rho');
[chk2,constChk] = lpvsynCheck(synSet,pvDenser,opt);
ghinf_4Cn = constChk(1).bound;
%
% G = ss(pdGwe,pvDenser.points);
% K = ss(pdK4n,pvDenser.points);
% Gcl = lft(G,K);
% E = eig(Gcl);
% for ii=1:length(pvDenser.points)
%     [wn,ep,poles] = damp(E(:,:,ii));
%     ReMax = max(real(poles));
%     fprintf('\tMinDecay = %6.4f\t MinDamping = %6.4f \tmaxFreq = %6.4f\n',...
%         -ReMax,min(ep),max(wn));
% end


% =========================================================================
% Evaluation/Validation

fprintf('\n')
fprintf('-------------------------------------------------------------------------------------\n')
fprintf('Comparison:\n')
fprintf('  constraint Hinf:\n')
fprintf('    X=X0, dX=0, Y=Y0,   dY=0:   %6.4f\n',ghinf_1)
fprintf('    X=X0, dX=0, Y=Y(p), dY=0:   %6.4f (%6.4f in a denser grid)\n',...
        ghinf_2,ghinf_2C)
for ii = 1:11
    fprintf('    X=X0, dX=0, Y=Y(p), dY~=0:  %6.4f (parameter rate %2.0f)\n',...
        ghinf_3(ii),ii-1);
end
fprintf('  constraint Hinf + Poles:\n')
fprintf('    X=X0, dX=0, Y=Y(p), dY~=0:   %6.4f (desing w/ a grid of %1.0f points)\n',...
        ghinf_4,nv);
if chk1, ppStr = 'ok'; else ppStr = 'no'; end
fprintf('    (parameter rate %1.0f)           %6.4f, const = %s (checking on a grid of %1.0f points)\n',...
        3,ghinf_4C,ppStr,2*nv);
fprintf('                                 %6.4f (redesing w/ a grid of %1.0f points)\n',...
        ghinf_4n,2*nv);
if chk2, ppStr = 'ok'; else ppStr = 'no'; end
fprintf('                                 %6.4f, const = %s (checking on a grid of %1.0f points)\n',...
        ghinf_4Cn,ppStr,3*nv);
fprintf('-------------------------------------------------------------------------------------\n')

figure
plot(0:10,ghinf_3,'rs')
xlabel('nu (rad/sec)')
ylabel('gamma')
title('Performance vs parameter rate')

return

