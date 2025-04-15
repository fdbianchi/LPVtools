% =========================================================================
%
% Testing basic use of pcss
%
% fbianchi - 2020-07-02
%
% =========================================================================

close all
clear all
clc

% Operating range
%
Zmin =  0.5;   Zmax =   4;
Mmin =  0.0;   Mmax = 106;  
vert = pgrid([Zmin, Zmax; Mmin, Mmax]);

% system matrices
a0 = [ 0 1; 0 0];  a1 = [-1 0; 0 0];   a2 = [ 0 0;-1 0];
b0 = [ 0;1];       b1 = [ 0; 0];       b2 = [ 0; 0];
c0 = [-1 0; 0 1];  c1 = [ 0 0; 0 0];   c2 = [ 0 0; 0 0];
d0 = [ 0; 0];      d1 = [ 0; 0];       d2 = [ 0; 0];

% Affine model:
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
nf1 = 2.0101;  df1 = [1.0000e+00   2.0101e-01];
w1 = tf(nf1,df1); w1.u = 'e'; w1.y = 'et';
%
%  Filter w2 to shape KS 
nf2 = [9.6785e+00   2.9035e-02   0 ];
df2 = [1.0000e+00   1.2064e+04   1.1360e+07   1.0661e+10];
w2 = tf(nf2,df2); w2.u = 'u'; w2.y = 'ut';
%
Wout = append(w1,w2,1,1);
Wout.y = {'et','ut','e','y2'};
%
% Specify the loop-shaping control structure
er = sumblk('e = r - y1');
pdGau1 = connect(pdG1,er,{'r','u'},{'e','u','e','y2'});

pdGau1 = Wout*pdGau1;

% controller design
const1 = synConst.Gain({'r'},{'et','ut'});
lyapSet.parfcnX  = @(p) p;
lyapSet.dparfcnX = 0;
lyapSet.parfcnY  = 0;
lyapSet.dparfcnY = 0;
[pdK,const] = lpvsyn(pdGau1,{'e','y2'},'u',const1,lyapSet);
g = const.bound;


%% ================================================================
% dimension functions

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking dimensions\n')

[ny1,nu1] = iosize(pdK);
ns1 = order(pdK);
np1 = npar(pdK);
[ny2,nu2,ns2,np2] = size(pdK);

if (nu1 ~= nu2)
    error('NU in IOSIZE and SIZE are not coincident')
end
if (ny1 ~= ny2)
    error('NY in IOSIZE and SIZE are not coincident')
end
fprintf('Sys dimensions using: size    iosize\n')
fprintf('Number of inputs:     %2.0f      %2.0f\n',nu1,nu2)
fprintf('Number of outputs:    %2.0f      %2.0f\n',ny1,ny2)

if (ns1 ~= ns2)
    error('NS in ORDER and SIZE are not coincident')
end
fprintf('Sys dimensions using: size    order\n')
fprintf('Number of states:     %2.0f      %2.0f\n',ns1,ns2)

if (np1 ~= np2)
    error('NP in NPAR and SIZE are not coincident')
end
fprintf('Sys dimensions using: size    npar\n')
fprintf('Number of parameters: %2.0f      %2.0f\n',np1,np2)



%% ================================================================
% Sub refering (NOT IMPLEMENTED YET)

% % this returns an lpvsys corresponding to the 1st output and the 1st input
% pdK1 = pdK(1,1)
% pdK1 = pdK('u','e')
% 
% % this returns an lpvsys corresponding to the 1st output and the 2nd input
% pdK2 = pdK(1,2)
% pdK2 = pdK('u','y2')
% 
% % all elements
% pdK3 = pdK(1,:)

%% ===============================================================
% Characteristics

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking characteristics\n\n')

% checking if is empty
bool = isempty(pdK);
if (bool == true)
    error('Error in ISEMPTY')
else
    disp('EMPTY working correctly')
end

% evaluation at frozen parameter values
K = subs(pdK,[Zmin;Mmin]);
Gau = subs(pdGau1,[Zmin;Mmin]);
Gcl = lft(Gau,K);
if (norm(Gcl,inf) <= g) && (max(real(eig(Gcl))) < 0)
    disp('SUBS working correctly')
else
    error('Error in SUBS')
end

% eigenvalues
fprintf('\nEigenvalues\n')
eig(pdK);
E = eig(pdK);
P = pole(pdK);
Et = zeros(ns1,1,4);
Gs = ss(pdK);
for ii = 1:4
    Et(:,ii) = eig(Gs.A(:,:,ii));
    if all(E(:,:,ii) == P(:,:,ii)) && all(Et(:,:,ii) == E(:,:,ii))
        fprintf('p = [%6.2f, %6.2f]: \t',vert(:,ii))
        fprintf('eigenvalues: [%g, %g, %g, %g, %g, %g]:\n',E(:,:,ii))
    else
        error('Error in EIG or POLE')
    end
end

% % transmission zeros
% fprintf('\nTransmission zeros\n')
% zz = tzero(pdG1)
% tzero(pdG1)

% DC gain
fprintf('\nDC gains\n')
dcgain(pdK);
gg = dcgain(pdK)


%% ================================================================
% Conversions
%
% to ss:
sys1a = ss(pdK);
gr    = pgrid([Zmin Zmax; Mmin Mmax]);
sys1b = ss(pdK,gr);
% to tf
sys1c = tf(pdK);
% to zpk
sys1d = zpk(pdK);

% checking Bode
bodemag(sys1a,sys1b,sys1c,sys1d);
legend('ss','ss(grid)','tf','zpk')


%% ========================================================================
% Graphics functions
figure
bodemag(pdK)

figure
bode(pdK)

figure
sigma(pdK)

figure
nyquist(pdK)

figure
step(pdK)

figure
initial(pdK,rand(order(pdK),1))

figure
impulse(pdK)

figure
pzmap(pdK)


