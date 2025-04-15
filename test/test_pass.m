
% -------------------------------------------------------------------------
% How to use and Examples for PASS
%
% -------------------------------------------------------------------------

% fbianchi - 2020-06-29

% *****************
clc
clear all
close all
% *****************

% =========================================================================
% Model data
%
% Operating range
Zmin =  0.5;   Zmax =   4;
Mmin =  0.0;   Mmax = 106;  
% names
pnames = {'p','q'};

% parameter set
setP = pset.Box([Zmin Zmax ; Mmin Mmax],[],pnames);

% system matrices
a0 = [ 0 1; 0 0]; a1 = [-1 0; 0 0]; a2 = [ 0 0;-1 0];
b0 = [ 0;1];      b1 = [ 0;0];      b2 = [ 0;0];
c0 = [-1 0;0 1];  c1 = [ 1 0; 0 0]; c2 = [ 0 0; 0 0];
d0 = [ 0;1];      d1 = [ 0;0];      d2 = [ 0;0];


% =========================================================================
% Model definitions
%
% ------------------------------------------------------------------------
% Using matrices
a = cat(3,a0,a1,a2);
b = cat(3,b0,b1,b2);
c = cat(3,c0,c1,c2);
d = cat(3,d0,d1,d2);

pdG1a = pass(a,b,c,d,[Zmin Zmax; Mmin Mmax]);
pdG1b = pass(a,b,c,d,setP);
pdG1  = pass(a,b,c,d,setP,...
            'InputName','u','OutputName','y');

% ------------------------------------------------------------------------
% Using ss objects
sys(:,:,1) = ss(a0,b0,c0,d0);
sys(:,:,2) = ss(a1,b1,c1,d1); % Z_al 
sys(:,:,3) = ss(a2,b2,c2,d2); % M_al 
%
% sintaxis 1
pdG2a = pass(sys,[Zmin Zmax; Mmin Mmax]);
pdG2a.u = 'u';           % input names
pdG2a.y = 'y';

% sintaxis 2
pdG2 = pass(sys,setP);
pdG2.u = 'u';           % input names
pdG2.y = 'y';

% ------------------------------------------------------------------------
% from uss object
p1 = ureal('p1',1,'Range',[Zmin Zmax]);
p2 = ureal('p2',1,'Range',[Mmin Mmax]);
A = a0 + p1*a1 + p2*a2;
B = b0 + p1*b1 + p2*b2;
C = c0 + p1*c1 + p2*c2;
D = d0 + p1*d1 + p2*d2;
sys = ss(A,B,C,D);
pdG3 = pass(sys);
pdG3.u = 'u';           % input names
pdG3.y = {'y1','y2'};

% ------------------------------------------------------------------------
% from psys object
pv = pvec('box',[Zmin Zmax ; Mmin Mmax]);
s0 = ltisys(a0,b0,c0,d0);
s1 = ltisys(a1,b1,c1,d1,0); % Z_al 
s2 = ltisys(a2,b2,c2,d2,0); % M_al 
pdG4 = pass(psys(pv,[s0 s1 s2]));
pdG4.u = 'u';           % input names
pdG4.y = {'y1','y2'};

%% ===============================================================

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking plots\n')

figure
bodemag(pdG1,'r.',pdG2,'kx',pdG3,'gs',pdG4)
legend('pdG1','pdG2','pdG3','pdG4')
title('Bodemag plot')

figure
bode(pdG1,'r.',pdG2,'kx',pdG3,'gs',pdG4)
legend('pdG1','pdG2','pdG3','pdG4')
title('Bode plot')

figure
sigma(pdG1,'r.',pdG2,'kx',pdG3,'gs',pdG4)
legend('pdG1','pdG2','pdG3','pdG4')
title('Sigma plot')

figure
nyquist(pdG1,'r.',pdG2,'kx',pdG3,'gs',pdG4)
legend('pdG1','pdG2','pdG3','pdG4')
title('Nyquist plot')

figure
step(pdG1,'r.',pdG2,'kx',pdG3,'gs',pdG4)
legend('pdG1','pdG2','pdG3','pdG4')
title('Step response')

figure
initial(pdG1,'r.',pdG2,'kx',pdG3,'gs',pdG4,[10;5])
legend('pdG1','pdG2','pdG3','pdG4')
title('Initial conditions response')

figure
impulse(pdG1,'r.',pdG2,'kx',pdG3,'gs',pdG4)
legend('pdG1','pdG2','pdG3','pdG4')
title('Impulse response')

figure
pzmap(pdG1,'r.',pdG2,'kx',pdG3,'gs',pdG4)
legend('pdG1','pdG2','pdG3','pdG4')
title('Pole-zero map')


%% ========================================================================
% Dimensions

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking dimensions\n')

[ny1,nu1] = iosize(pdG1);
ns1 = order(pdG1);
np1 = npar(pdG1);
nv1 = nsys(pdG1);
[ny2,nu2,ns2,np2,nv2] = size(pdG1);

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

if (nv1 ~= nv2)
    error('NV in NPAR and SIZE are not coincident')
end
fprintf('Sys dimensions using: size    nsys\n')
fprintf('Number of systems:    %2.0f      %2.0f\n',nv1,nv2)


%% ================================================================
% Characteristics

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking characteristics\n\n')

% checking if is empty
bool = isempty(pdG1);
if (bool == true)
    error('Error in ISEMPTY')
else
    disp('EMPTY working correctly')
end

% checking if parameter dependent
dx = ispd(pdG1);
disp('')
if (all(dx == [1 0 1]))
    disp('ISPD working correctly')
else
    error('Error in ISPD')
end

% evaluation at frozen parameter values
syse = subs(pdG1,[Zmin;Mmin]);
disp('')
if (all(syse.a == (a0 + Zmin*a1 + Mmin*a2)))
    disp('SUBS working correctly')
else
    error('Error in SUBS')
end

% eigenvalues
fprintf('\nEigenvalues\n')
eig(pdG1);
E = eig(pdG1);
P = pole(pdG1);
Ap = zeros(2,2,4); Et = zeros(2,1,4);
for ii = 1:4
    p = setP.points(:,ii);
    Ap(:,:,ii) = a0 + p(1)*a1 + p(2)*a2;
    Et(:,ii) = eig(Ap(:,:,ii));
    if all(E(:,:,ii) == P(:,:,ii)) && all(Et(:,:,ii) == E(:,:,ii))
        fprintf('p = [%6.2f, %6.2f]: \t',p)
        fprintf('eigenvalues: [%g, %g]:\n',E(:,:,ii))
    else
        error('Error in EIG or POLE')
    end
end

% transmission zeros
fprintf('\nTransmission zeros\n')
zz = tzero(pdG1)
tzero(pdG1)

% DC gain
fprintf('\nDC gains\n')
dcgain(pdG1);
gg = dcgain(pdG1)


%% ================================================================
% Sub refering

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking Sub refering\n')

% this returns the independent term
fprintf('\nIndependent term\n')
s1 = pdG1(1);
if ~isequal(a0,s1.A) || ~isequal(b0,s1.B) || ~isequal(c0,s1.C) || ~isequal(d0,s1.D)
    error('Error in SUBSREf with 1 argument')
else
    disp('SUBSREf with 1 argument: OK')
end

% this returns a pass corresponding to the 2nd output and the 1st input
fprintf('\nLPV Submodel for output 2, input 1\n') 
pdG1s1 = pdG1(2,1);
pdG1s2 = pdG1('y(2)','u');
if ~isequal(a,pdG1s1.A) || ~isequal(b(:,1,:),pdG1s1.B) ||...
        ~isequal(c(2,:,:),pdG1s1.C) || ~isequal(d(2,1,:),pdG1s1.D)
    error('Error in SUBSREf with 2 argument')
elseif ~isequal(a,pdG1s2.A) || ~isequal(b(:,1,:),pdG1s2.B) || ...
        ~isequal(c(2,:,:),pdG1s2.C) || ~isequal(d(2,1,:),pdG1s2.D)
    error('Error in SUBSREf with 2 argument')
else
    disp('SUBSREf with 2 argument: OK')
end

% this returns the independent term corresponding to the 2nd output and 
% the 1st input
fprintf('\nIndependent term for output 2, input 1\n') 
idx1 = 1; idx2 = 1;
s21 = pdG1(idx1,idx2,1);
s22 = pdG1('y(1)','u',1);
if ~isequal(a(:,:,1),s21.A) || ~isequal(b(:,idx2,1),s21.B) ||...
        ~isequal(c(idx1,:,1),s21.C) || ~isequal(d(idx1,idx2,1),s21.D)
    error('Error in SUBSREf with 3 argument')
elseif ~isequal(a(:,:,1),s22.A) || ~isequal(b(:,idx2,1),s22.B) ||...
        ~isequal(c(idx1,:,1),s22.C) || ~isequal(d(idx1,idx2,1),s22.D)
    error('Error in SUBSREf with 3 argument')
else
    disp('SUBSREf with 3 argument: OK')
end


%% ================================================================
% Conversions
%

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking conversions\n')

% to ss:
sys1a = ss(pdG1);
gr    = pgrid([Zmin Zmax; Mmin Mmax]);
sys1b = ss(pdG1,gr);
% to tf
sys1c = tf(pdG1);
% to zpk
sys1d = zpk(pdG1);

% to uss:
sys2 = uss(pdG1);
pval = pdG1.parset.points;
sys1e = usubs(sys2,'p',pval(1,:),'q',pval(2,:));

% checking Bode
bodemag(sys1a,sys1b,sys1c,sys1d,sys1e);
legend('ss','ss(grid)','tf','zpk','uss')
title('Model conversions')

% to psys:
sys3 = psys(pdG1);

% normalization
pdG1n = normalize(pdG1);
pdG1d = denormalize(pdG1n, pdG1.parset);
figure
bodemag(pdG1,'-r',pdG1n,'kx',pdG1d,'gs');
legend('Original','Normalized','Denormalized')
title('Model normalizations')


%% ===============================================================
% interconnections

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking interconnections\n')

% series connection
G1 = tf(1,[1 2]); G1.u = 'w'; G1.y = 'z';
G2 = ss(-eye(2),eye(2),eye(2),zeros(2));  G2.u = 'p'; G2.y = 'q';
G3 = ss(eye(2));  G2.u = 'p'; G2.y = 'q';
s1 = ss(G2*pdG1);
s2 = G2*ss(pdG1);
s3 = G3*ss(pdG1);
figure
bodemag(s1,'r',s2,'kx')
legend('Using LPV','Using LTI')
title('Series interconnection, LTI and LPV')

pdG2 = pdG1(1,1); pdG2.u = 'w'; pdG2.y = 'z';
pdG12 = pdG1*pdG2;
s12 = ss(pdG1)*ss(pdG2);
figure
bodemag(pdG12,'r',s12,'kx')
legend('Using LPV','Using LTI')
title('Series interconnection, LPV and LPV')


% parallel connection
G3 = ss(ones(2,1));
G3 = ss(-eye(2),[1;1],eye(2),zeros(2,1));
s1 = ss(parallel(pdG1,G3));
s2 = ss(pdG1) + G3;
figure
bodemag(s1,'r',s2,'kx')
legend('Using LPV','Using LTI')
title('Parallel interconnection')

% uminus and minus functions
s1 = ss(parallel(pdG1,-G3));
s2 = ss(pdG1) - G3;
figure
bodemag(s1,'r',s2,'kx')
legend('Using LPV','Using LTI')
title('Minus interconnection')

% append
s1 = ss(append(pdG1,G3));
s2 = append(ss(pdG1),G3);
figure
bodemag(s1,'r',s2,'kx')
legend('Using LPV','Using LTI')
title('Append interconnection')


%% ===============================================================
% extreme cases

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking extreme cases\n\n')

% Constant A
disp('Extreme cases: A constant')
a = rand(2);
b = rand(2,2,3);
c = rand(1,2,2);
d = 0;
pv = pset.Box([1 3;2 5]);
pdG = pass(a,b,c,d,pv)

% Constant A, B
disp('Extreme cases: A & B constant')
a = rand(2);
b = rand(2);
c = rand(1,2,3);
d = 0;
pdG = pass(a,b,c,d,pv)

% Constant A, B, C and D
disp('Extreme cases: A, B, C and D constant')
sys = rss(2,3,2);
pdG = pass(sys);
G = ss(pdG);
if ~isequal(pdG.A,G.A) || ~isequal(pdG.B,G.B) ||...
        ~isequal(pdG.C,G.C) || ~isequal(pdG.D,G.D)
    error('Error ')
else
    disp('OK')
end

% submodels not parameter dependent
a = diag([-2,-3]);
b = cat(3,eye(2),diag([0,1]));
c = eye(2);
d = zeros(2);
pdG = pass(a,b,c,d,[2 3]);

pdG_1 = pdG(:,1);
pdG_2 = pdG(:,2);






