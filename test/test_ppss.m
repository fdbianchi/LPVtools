
% ------------------------------------------------------------------------
% How to use and Examples for PPSS
%
% ------------------------------------------------------------------------

% fbianchi - 2020-06-29

% *****************
clc
clear all
close all
% *****************

% ========================================================================
% Model data
%
% Operating range
Zmin = 0.5;   Zmax =   4;
Mmin = 0.0;   Mmax = 106;  
vert = pgrid([Zmin, Zmax; Mmin, Mmax]);
% names
pnames = {'p','q'};

% parameter set (two options, same result)
set1 = pset.Hull(vert);
set2 = pset.Box([Zmin Zmax ; Mmin Mmax]);

% system matrices
a0 = [ 0 1; 0 0]; a1 = [-1 0; 0 0]; a2 = [ 0 0;-1 0];
b0 = [0;1];
c0 = [-1 0;0 1];
d0 = [0;1];

% ========================================================================
% Model constructions
%
% ------------------------------------------------------------------------
% Using matrices
a = zeros(2,2,4);
b = zeros(2,1,4);
c = zeros(2,2,4);
d = zeros(2,1,4);
for ii = 1:4
    a(:,:,ii) = a0 + vert(1,ii)*a1 + vert(2,ii)*a2;
    b(:,:,ii) = b0;
    c(:,:,ii) = c0;
    d(:,:,ii) = d0;
end

pdG1a = ppss(a,b,c,d,vert);
pdG1b = ppss(a,b,c,d,set1);
pdG1c = ppss(a,b,c,d,set2,...
            'InputName','u','OutputName','y');

        
% ------------------------------------------------------------------------
% From an affine model
sysAff(:,:,1) = ss(a0,b0,c0,d0);
sysAff(:,:,2) = ss(a1,0*b0,0*c0,0*d0);
sysAff(:,:,3) = ss(a2,0*b0,0*c0,0*d0);
pdGa = pass(sysAff,[Zmin Zmax ; Mmin Mmax]);
%
% politopic model:
for ii = 1:4
    sysPol(:,:,ii) = ss(a0 + a1*vert(1,ii) + a2*vert(2,ii), b0, c0, d0);
end
%
% ------------------------------------------------------------------------
% sintaxis 1: in this case the parset created by ppss is gral and subs
% function cannot by used
pdG1 = ppss(sysPol,vert);
pdG1.u = 'u';               
pdG1.y = {'y1','y2'};

% sintaxis 2
pdG2 = ppss(sysPol,set1);
pdG2.u = 'u';               
pdG2.y = {'y1','y2'};

% sintaxis 3: from psys object
pv = pvec('pol',vert);
sysPol = [];
for ii=1:4
    sysPol = [sysPol ltisys(a0 + a1*vert(1,ii) + a2*vert(2,ii), b0, c0, d0)];
end
pds = psys(sysPol);
pds = addpv(pds,pv);
pdG3 = ppss(pds);
pdG3.u = 'u';               
pdG3.y = {'y1','y2'};

% sintaxis 4: from affine model
pdG4 = ppss(pdGa);
pdG4.u = 'u';               
pdG4.y = {'y1','y2'};


%% ===============================================================
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


%% ===============================================================
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
if (all(dx == [1 0 0]))
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
Et = zeros(2,1,4);
for ii = 1:4
    Et(:,ii) = eig(a(:,:,ii));
    if all(E(:,:,ii) == P(:,:,ii)) && all(Et(:,:,ii) == E(:,:,ii))
        fprintf('p = [%6.2f, %6.2f]: \t',vert(:,ii))
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



%% ===============================================================
% Sub refering

idx1 = 2; idx2 = 1; idx3 = 1;

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking Sub refering\n')

% this returns the LTI model of the first vertice
fprintf('\nFirst term in the polytopic representation\n')
s1 = pdG1(1);
if ~isequal(a(:,:,idx3),s1.A) || ~isequal(b(:,:,idx3),s1.B) ||...
        ~isequal(c(:,:,idx3),s1.C) || ~isequal(d(:,:,idx3),s1.D)
    error('Error in SUBSREf with 1 argument')
else
    disp('SUBSREf with 1 argument: OK')
end

% this returns a pass corresponding to the 2nd output and the 1st input
fprintf('\nLPV Submodel for output 2, input 1\n') 
pdG1s1 = pdG1(idx1,idx2);
pdG1s2 = pdG1('y2','u');
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
fprintf('\nFirst term for output 2, input 1\n') 
idx1 = 1; idx2 = 1;
s21 = pdG1(idx1,idx2,idx3);
s22 = pdG1('y1','u',idx3);
if ~isequal(a(:,:,idx3),s21.A) || ~isequal(b(:,idx2,idx3),s21.B) ||...
        ~isequal(c(idx1,:,idx3),s21.C) || ~isequal(d(idx1,idx2,idx3),s21.D)
    error('Error in SUBSREf with 3 argument')
elseif ~isequal(a(:,:,idx3),s22.A) || ~isequal(b(:,idx2,idx3),s22.B) ||...
        ~isequal(c(idx1,:,idx3),s22.C) || ~isequal(d(idx1,idx2,idx3),s22.D)
    error('Error in SUBSREf with 3 argument')
else
    disp('SUBSREf with 3 argument: OK')
end



%% ===============================================================
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

% checking Bode
bodemag(sys1a,sys1b,sys1c,sys1d);
legend('ss','ss(grid)','tf','zpk')
title('Model conversions')

% to psys:
sys3 = psys(pdG1);



%% ===============================================================
% interconnections

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking interconnections\n')

% series connection
G1 = tf(1,[1 2]); G1.u = 'w'; G1.y = 'z';
G2 = ss(-eye(2),eye(2),eye(2),zeros(2));  G2.u = 'p'; G2.y = 'q';
G2 = ss(eye(2));  G2.u = 'p'; G2.y = 'q';
s1 = ss(G2*pdG1);
s2 = G2*ss(pdG1);
figure
bodemag(s1,'r',s2,'kx')
legend('Using LPV','Using LTI')
title('Series interconnection')


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
% Graphics functions

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

%% ===============================================================
% extreme cases

fprintf('\n------------------------------------------------------------\n')
fprintf('Checking extreme cases\n\n')

% Constant A
disp('Extreme cases: A constant')
a = rand(2);
b = rand(2,2,3);
c = rand(1,2,3);
d = 0;
pv = pset.Gral([1 3 5]);
pdG = ppss(a,b,c,d,pv)

% Constant A, B
disp('Extreme cases: A & B constant')
a = rand(2);
b = rand(2);
c = rand(1,2);
d = 0;
pdG = ppss(a,b,c,d)

% Constant A, B, C and D
disp('Extreme cases: A, B, C and D constant')
sys = rss(2,3,2);
pdG = ppss(sys);
G = ss(pdG);
if ~isequal(pdG.A,G.A) || ~isequal(pdG.B,G.B) ||...
        ~isequal(pdG.C,G.C) || ~isequal(pdG.D,G.D)
    error('Error ')
else
    disp('OK')
end

