
% =========================================================================
%
% Testing basic use of lpvsyn
%
% fbianchi - 2020-04-24
%
% =========================================================================

close all
clear all
clc

% Operating range
%
Zmin =  0.5;   Zmax =   4;
Mmin =  0.0;   Mmax = 106;  

% system matrices
a0 = [ 0 1; 0 0];  a1 = [-1 0; 0 0];   a2 = [ 0 0;-1 0];
b0 = [ 0;1];       b1 = [ 0; 0];       b2 = [ 0; 0];
c0 = [-1 0; 0 1];  c1 = [ 0 0; 0 0];   c2 = [ 0 0; 0 0];
d0 = [ 0; 0];      d1 = [ 0; 0];       d2 = [ 0; 0];

%% ===============================================================
% Design with lmitool
%
% parameter set
pv = pvec('box',[Zmin Zmax ; Mmin Mmax]);
Vx = polydec(pv);
%
% Affine model:
s0 = ltisys(a0,b0,c0,d0);
s1 = ltisys(a1,b1,c1,d1,0); % Z_al 
s2 = ltisys(a2,b2,c2,d2,0); % M_al 

pdG0 = psys(pv,[s0 s1 s2]);

% Controller design
%
nf1 = 2.0101;  df1 = [1.0000e+00   2.0101e-01];
w1  = ltisys('tf',nf1,df1);
%  Filter w2 to shape KS 
nf2 = [9.6785e+00   2.9035e-02   0            0];
df2 = [1.0000e+00   1.2064e+04   1.1360e+07   1.0661e+10];
w2  = ltisys('tf',nf2,df2);
%

%  Specify the loop-shaping control structure with SCONNECT
inputs  = 'r';
outputs = 'e=r-G(1);K';
K_in    = 'K:e;G(2)';
G_in    = 'G:K';
[pdGau0,r] = sconnect(inputs,outputs,K_in,G_in,pdG0);

%  Augment with the shaping filters
pdGwe0 = smult(pdGau0,sdiag(w1,w2,eye(2)));

% Controller design with Hinf
[gLMI,pdK0] = hinfgs(pdGwe0,r);
pdK0 = addpv(pdK0,pv);


%% ===============================================================
% Design with new tools
% 

% Using with lmitool objects
[gLMIb,pdK0b] = lpvsyn(pdGwe0,r);

fprintf('\n')
fprintf('Performance hinfgs:   %6.4f\n',gLMI);
fprintf('Performance lpvsyn:   %6.4f (old syntaxix)\n',gLMIb);
fprintf('\n')

% ------------------------------------------------------------------------
% New objects
% Construction as affine models
%
%
% Affine model:
clear sys
sys(:,:,1) = ss(a0,b0,c0,d0);
sys(:,:,2) = ss(a1,b1,c1,d1); % Z_al 
sys(:,:,3) = ss(a2,b2,c2,d2); % M_al 
%
% lpv model
pdG1 = pass(sys,[Zmin Zmax ; Mmin Mmax]);
% input names
pdG1.u = 'u';
pdG1.y = {'y1','y2'};

% Augmented plant construction
%
% Filter w1 to shape S
% n1 = 2.0101;  d1 = [1.0000e+00   2.0101e-01];
w1 = tf(nf1,df1); w1.u = 'e'; w1.y = 'et';
%
%  Filter w2 to shape KS 
% n2 = [9.6785e+00   2.9035e-02   0 ];
% d2 = [1.0000e+00   1.2064e+04   1.1360e+07   1.0661e+10];
w2 = tf(nf2,df2); w2.u = 'u'; w2.y = 'ut';
%
Wout = append(w1,w2,1,1);
Wout.y = {'et','ut','e','y2'};
%
% Specify the loop-shaping control structure
er = sumblk('e = r - y1');
pdGau1 = connect(pdG1,er,{'r','u'},{'e','u','e','y2'});

pdGau1 = Wout*pdGau1;


% ------------------> Affine case <----------------------
%
% Hinf default
[pdKaff1,const] = lpvsyn(pdGau1,[2,1]);
gAff1 = const.bound;

% Hinf explicit
const1 = synConst.Gain({'r'},{'et','ut'});
[pdKaff2,const] = lpvsyn(pdGau1,{'e','y2'},'u',const1);
gAff2 = const.bound;

% Hinf (two channels)
const2(1) = synConst.Gain({'r'},{'et'});
const2(2) = synConst.Gain({'r'},{'ut'},'bound',0.05);
[pdKaff3,const] = lpvsyn(pdGau1,{'e','y2'},'u',const2);
gAff3 = const.bound;

% Hinf + Poles
const2(1) = synConst.Gain({'r'},{'et','ut'});
const2(2) = synConst.Poles('MinDecay',0.2);
[pdKaff4,const] = lpvsyn(pdGau1,{'e','y2'},'u',const2);
gAff4 = const(1).bound;

% Hinf - X(p), dX=0
lyapSet.parfcnX  = @(p) p;
lyapSet.dparfcnX = 0;
lyapSet.parfcnY  = 0;
lyapSet.dparfcnY = 0;
[pdKaff5,const] = lpvsyn(pdGau1,{'e','y2'},'u',const1,lyapSet);
gAff5 = const.bound;

% Hinf - affX, dX=0
[pdKaff6,const] = lpvsyn(pdGau1,{'e','y2'},'u',const1,'affX');
gAff6 = const.bound;
 
% Hinf - affX, dX=0 - gral controller
ctrlFcn = [1 2 3];
[pdKaff7,const] = lpvsyn(pdGau1,{'e','y2'},'u',const1,'affX',ctrlFcn);
gAff7 = const.bound;
 
% Hinf+Poles - affX, dX=0 - gral controller
const3(1) = synConst.Gain({'r'},{'et','ut'});
const3(2) = synConst.Poles('MinDecay',0.2);
[pdKaff8,const] = lpvsyn(pdGau1,{'e','y2'},'u',const3,'affX',ctrlFcn);
gAff8 = const(1).bound;

% Hinf - X cte - gral controller
[pdKaff9,const] = lpvsyn(pdGau1,{'e','y2'},'u',const1,[],ctrlFcn);
gAff9 = const.bound;
 
% Hinf+Poles - X cte - gral controller
[pdKaff10,const] = lpvsyn(pdGau1,{'e','y2'},'u',const3,[],ctrlFcn);
gAff10 = const(1).bound;

% % pdG affine, Y affine parameter dependent + gral controller + poles
% % constraint
% const3(1) = synConst.Gain({'r'},{'et','ut'});
% const3(2) = synConst.Poles('MinDecay',0.2);
% [pdKaff9,~,const] = lpvsyn(pdGau1,{'e','y2'},'u',const3,lyapSet,ctrlFcn);
% gAff9 = const(1).bound;

fprintf('--------------------------------------------------------------------------\n')
fprintf('Performance hinfgs:   %6.4f\n',gLMI);
fprintf('--------------------------------------------------------------------------\n')
fprintf('Affine pdG (pass)\n')
fprintf('Performance lpvsyn:   %6.4f (Hinf, nargin=2)           - Ctrl class: %s\n',...
    gAff1,class(pdKaff1));
fprintf('Performance lpvsyn:   %6.4f (Hinf, nargin=4)           - Ctrl class: %s\n',...
    gAff2,class(pdKaff2));
fprintf('Performance lpvsyn:   %6.4f (2xHinf)                   - Ctrl class: %s\n',...
    gAff3,class(pdKaff3));
fprintf('Performance lpvsyn:   %6.4f (Hinf+PP)                  - Ctrl class: %s\n',...
    gAff4,class(pdKaff4));
fprintf('Performance lpvsyn:   %6.4f (Hinf, X(p))               - Ctrl class: %s\n',...
    gAff5,class(pdKaff5));
fprintf('Performance lpvsyn:   %6.4f (Hinf, affX)               - Ctrl class: %s\n',...
    gAff6,class(pdKaff6));
fprintf('Performance lpvsyn:   %6.4f (Hinf, affX, pdK(p))       - Ctrl class: %s\n',...
    gAff7,class(pdKaff7));
fprintf('Performance lpvsyn:   %6.4f (Hinf+PP, affX, pdK(p))    - Ctrl class: %s\n',...
    gAff8,class(pdKaff8));
fprintf('Performance lpvsyn:   %6.4f (Hinf, X cte, pdK(p))      - Ctrl class: %s\n',...
    gAff9,class(pdKaff9));
fprintf('Performance lpvsyn:   %6.4f (Hinf+PP, X cte, pdK(p))   - Ctrl class: %s\n',...
    gAff10,class(pdKaff10));
fprintf('--------------------------------------------------------------------------\n')
fprintf('\n')


%%
% -----------------> Polytopic case <---------------------
% Hinf explicit
[pdKpol1,const] = lpvsyn(ppss(pdGau1),{'e','y2'},'u',const1);
gpol1 = const.bound;

% Hinf (two channels)
const2(1) = synConst.Gain({'r'},{'et'});
const2(2) = synConst.Gain({'r'},{'ut'},'bound',0.05);
[pdKpol2,const] = lpvsyn(ppss(pdGau1),{'e','y2'},'u',const2);
gpol2 = const.bound;

% Hinf - X(p), dX=0
lyapSet.parfcnX  = @(p) p;
lyapSet.dparfcnX = 0;
lyapSet.parfcnY  = 0;
lyapSet.dparfcnY = 0;
[pdKpol3,const] = lpvsyn(ppss(pdGau1),{'e','y2'},'u',const1,lyapSet);
gpol3 = const.bound;

% Hinf+Poles -  X(p), dX=0
[pdKpol4,const] = lpvsyn(ppss(pdGau1),{'e','y2'},'u',const3,lyapSet);
gpol4 = const(1).bound;

% Hinf - X(p), dX=0 - gral controller
% ctrlFcn =@(p) p;
% [pdKpol5,~,const] = lpvsyn(ppss(pdGau1),{'e','y2'},'u',const1,lyapSet,ctrlFcn);
% gpol5 = const(1).bound;

% Hinf+Poles - X(p), dX=0 - gral controller
% ctrlFcn =@(p) p;
% [pdKpol6,~,const] = lpvsyn(ppss(pdGau1),{'e','y2'},'u',const3,lyapSet,ctrlFcn);
% gpol6 = const(1).bound;


fprintf('--------------------------------------------------------------------------\n')
fprintf('Performance hinfgs:   %6.4f\n',gLMI);
fprintf('--------------------------------------------------------------------------\n')
fprintf('Polytopic pdG (ppss)\n')
fprintf('Performance lpvsyn:   %6.4f (Hinf, nargin=4)           - Ctrl class: %s\n',...
    gpol1,class(pdKpol1));
fprintf('Performance lpvsyn:   %6.4f (2xHinf, nargin=4)         - Ctrl class: %s\n',...
    gpol2,class(pdKpol2));
fprintf('Performance lpvsyn:   %6.4f (Hinf, X(p))               - Ctrl class: %s\n',...
    gpol3,class(pdKpol3));
fprintf('Performance lpvsyn:   %6.4f (Hinf+PP, X(p))            - Ctrl class: %s\n',...
    gpol2,class(pdKpol4));
% fprintf('Performance lpvsyn:   %6.4f (Hinf, X(p), pdK(p))       - Ctrl class: %s\n',...
%     gpol5,class(pdKpol5));
% fprintf('Performance lpvsyn:   %6.4f (Hinf+PP, X(p), pdK(p))    - Ctrl class: %s\n',...
%     gpol6,class(pdKpol6));
fprintf('--------------------------------------------------------------------------\n')
fprintf('\n')

return

