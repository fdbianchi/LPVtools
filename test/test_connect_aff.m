
% -------------------------------------------------------------------------
% test for connect in case of affine and general
%
% -------------------------------------------------------------------------

% fbianchi - 2023-07-07

% cleanig 
clc
clearvars
close all

%% ========================================================================
% affine example

% Operating range
Zmin =  0.5;   Zmax =   4;
Mmin =  0.0;   Mmax = 106;  
% system matrices
a0 = [ 0 1; 0 0]; a1 = [-1 0; 0 0]; a2 = [ 0 0;-1 0];
b0 = [0;1];       b1 = [0;0];       b2 = [0;0];
c0 = [-1 0;0 1];  c1 = [ 0 0; 0 0]; c2 = [ 0 0; 0 0];
d0 = [0;0];       d1 = [0;0];       d2 = [0;0];

% -------------------------------------------------------------------------
% Using lmitool
% parameter set
pv = pvec('box',[Zmin Zmax; Mmin Mmax]);
Vx = polydec(pv);
% Affine model:
s0 = ltisys(a0,b0,c0,d0);
s1 = ltisys(a1,b1,c1,d1,0); % Z_al 
s2 = ltisys(a2,b2,c2,d2,0); % M_al 
pdGlmi = psys(pv,[s0 s1 s2]);

% Weight on S
nf1 = 2.0101;  df1 = [1.0000e+00   2.0101e-01];
% Weight on KS
nf2 = [9.6785e+00   2.9035e-02   0            0];
df2 = [1.0000e+00   1.2064e+04   1.1360e+07   1.0661e+10];
W1lmi = ltisys('tf',nf1,df1);
W2lmi = ltisys('tf',nf2,df2);

% augmented plant
inputs  = 'r';
outputs = 'e=r-G(1);K';
K_in    = 'K:e;G(2)';
G_in    = 'G:K';
[pdGauLmi,r] = sconnect(inputs,outputs,K_in,G_in,pdGlmi);

% augmented plant + weight
pdGawLmi = smult(pdGauLmi,sdiag(W1lmi,W2lmi,eye(2)));


% -------------------------------------------------------------------------
% Using LPVtool

% system matrices
A = cat(3,a0,a1,a2);
B = cat(3,b0,b1,b2);
C = cat(3,c0,c1,c2);
D = cat(3,d0,d1,d2);
% parameter set
Pset = pset.Box([Zmin Zmax; Mmin Mmax],[-10 10],{'Z','M'});

% affine lpv model
pdG = pass(A,B,C,D,Pset,'InputName','u','OutputName','y');

% Weights
W1 = mat2lti(W1lmi);
W2 = mat2lti(W2lmi);

% augmented plant
er  = sumblk('e = r - y(1)');
pdGau = connect(pdG,er,{'r','u'},{'e','u','e','y(2)'});

% augmented plant + weight
Wout  = append(W1,W2,1,1);   
Wout.y = {'et','ut','e','y(2)'};
pdGaw = Wout*pdGau;
pdGaw.y = {'et','ut','e','y(2)'};
pdGaw.u = {'r','u'};

% -------------------------------------------------------------------------
% comparison
pdGaw_lmi = ppss(pdGawLmi);

figure
bodemag(pdGaw, pdGaw_lmi)
legend('Lmitool','LPVtool')




%% ========================================================================
% general example

% grid of points
nv    = 7;
range = [-pi,pi];
rate  = 2*[-1 1];
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

% weighting functions
sb1 = sumblk('e = v - r',2);
Wp = ss(eye(2));                    Wp.u = 'e';  Wp.y = 'ep';
Wn = tf(10*[1 10],[1 1000])*eye(2); Wn.u = 'dn'; Wn.y = 'n';
sb2 = sumblk('y = v - r + n',2);
Wf = ss(1);                         Wf.u = 'df'; Wf.y = 'f';
Wu = ss(1/280*eye(2));              Wu.u = 'u';  Wu.y = 'eu';
Wr = tf(20,[1 0.2])*eye(2);         Wr.u = 'dr'; Wr.y = 'r';

% augmented plant
pdGwe = connect(pdG,sb1,sb2,Wp,Wn,Wf,Wu,Wr,{'dr','dn','df','u'},{'ep','eu','y'});


% comparison with polytopic
G = ss(pdG);
Gwe = connect(G,sb1,sb2,Wp,Wn,Wf,Wu,Wr,{'dr','dn','df','u'},{'ep','eu','y'});

figure
bodemag(pdGwe, Gwe)
legend('LPVtool','LTI @ grid points')






