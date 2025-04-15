
% -------------------------------------------------------------------------
% test for parsysdata
%
% -------------------------------------------------------------------------

% fbianchi - 2020-07-01
clearvars; clc

% LTI objects
ni = 5; no = 3; ns = 5;
sys = rss(ns,no,ni);

% 1) ios numeric
ny = 1; nz = no - ny;
nu = 2; nw = ni - nu;
[a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(sys,[ny nu]);

sys11 = ss(a,b1,c1,d11); norm11 = norm(sys11 - sys(1:nz,1:nw),inf);
sys12 = ss(a,b2,c1,d12); norm12 = norm(sys12 - sys(1:nz,nw+1:ni),inf);
sys21 = ss(a,b1,c2,d21); norm21 = norm(sys21 - sys(nz+1:no,1:nw),inf);
sys22 = ss(a,b2,c2,d22); norm22 = norm(sys22 - sys(nz+1:no,nw+1:ni),inf);

fprintf('Checking vector input\n')
fprintf(' ||G11 - G(1:nz,1:nw)||       = %d\n',norm11)
fprintf(' ||G12 - G(1:nz,nw+1:ni)||    = %d\n',norm12)
fprintf(' ||G21 - G(nz+1:no,1:nw)||    = %d\n',norm21)
fprintf(' ||G22 - G(nz+1:no,nw+1:ni)|| = %d\n',norm22)
fprintf('\n')

% 2) ios struct (numeric)
ios.dist = 1:nw;
ios.ctrl = nw+1:ni;
ios.perf = 1:nz;
ios.meas = nz+1:no;
[a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(sys,ios);

sys11 = ss(a,b1,c1,d11); norm(sys11 - sys(1:nz,1:nw),inf); 
sys12 = ss(a,b2,c1,d12); norm(sys12 - sys(1:nz,nw+1:ni),inf);
sys21 = ss(a,b1,c2,d21); norm(sys21 - sys(nz+1:no,1:nw),inf);
sys22 = ss(a,b2,c2,d22); norm(sys22 - sys(nz+1:no,nw+1:ni),inf);

fprintf('Checking struct input\n')
fprintf(' ||G11 - G(1:nz,1:nw)||       = %d\n',norm11)
fprintf(' ||G12 - G(1:nz,nw+1:ni)||    = %d\n',norm12)
fprintf(' ||G21 - G(nz+1:no,1:nw)||    = %d\n',norm21)
fprintf(' ||G22 - G(nz+1:no,nw+1:ni)|| = %d\n',norm22)
fprintf('\n')

% 3) ios struct (numeric) + sys with groups
clear ios
ios.dist = 1:nw;
ios.perf = 1:nz;
sys.InputGroup.ctrl = nw+1:ni;
sys.OutputGroup.meas = nz+1:no;
[a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(sys,ios);

sys11 = ss(a,b1,c1,d11); norm(sys11 - sys(1:nz,1:nw),inf);
sys12 = ss(a,b2,c1,d12); norm(sys12 - sys(1:nz,nw+1:ni),inf);
sys21 = ss(a,b1,c2,d21); norm(sys21 - sys(nz+1:no,1:nw),inf);
sys22 = ss(a,b2,c2,d22); norm(sys22 - sys(nz+1:no,nw+1:ni),inf);

fprintf('Checking struct + sys groups input\n')
fprintf(' ||G11 - G(1:nz,1:nw)||       = %d\n',norm11)
fprintf(' ||G12 - G(1:nz,nw+1:ni)||    = %d\n',norm12)
fprintf(' ||G21 - G(nz+1:no,1:nw)||    = %d\n',norm21)
fprintf(' ||G22 - G(nz+1:no,nw+1:ni)|| = %d\n',norm22)
fprintf('\n')

% 4) ios struct (names)
sys.u = 'u';
sys.y = 'y';
ios.dist = sprintfc('u(%d)',1:nw);
ios.ctrl = sprintfc('u(%d)',nw+1:ni);
ios.perf = sprintfc('y(%d)',1:nz);
ios.meas = sprintfc('y(%d)',nz+1:no);
[a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(sys,ios);

sys11 = ss(a,b1,c1,d11); norm(sys11 - sys(1:nz,1:nw),inf);
sys12 = ss(a,b2,c1,d12); norm(sys12 - sys(1:nz,nw+1:ni),inf);
sys21 = ss(a,b1,c2,d21); norm(sys21 - sys(nz+1:no,1:nw),inf);
sys22 = ss(a,b2,c2,d22); norm(sys22 - sys(nz+1:no,nw+1:ni),inf);

fprintf('Checking struct with names\n')
fprintf(' ||G11 - G(1:nz,1:nw)||       = %d\n',norm11)
fprintf(' ||G12 - G(1:nz,nw+1:ni)||    = %d\n',norm12)
fprintf(' ||G21 - G(nz+1:no,1:nw)||    = %d\n',norm21)
fprintf(' ||G22 - G(nz+1:no,nw+1:ni)|| = %d\n',norm22)
fprintf('\n')

% 5) ios struct (names) + sys with groups
clear ios
ios.dist = sprintfc('u(%d)',1:nw);
ios.perf = sprintfc('y(%d)',1:nz);
sys.InputGroup.ctrl = nw+1:ni;
sys.OutputGroup.meas = nz+1:no;
[a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(sys,ios);

sys11 = ss(a,b1,c1,d11); norm(sys11 - sys(1:nz,1:nw),inf); 
sys12 = ss(a,b2,c1,d12); norm(sys12 - sys(1:nz,nw+1:ni),inf);
sys21 = ss(a,b1,c2,d21); norm(sys21 - sys(nz+1:no,1:nw),inf);
sys22 = ss(a,b2,c2,d22); norm(sys22 - sys(nz+1:no,nw+1:ni),inf);

fprintf('Checking struct (names) + sys groups input\n')
fprintf(' ||G11 - G(1:nz,1:nw)||       = %d\n',norm11)
fprintf(' ||G12 - G(1:nz,nw+1:ni)||    = %d\n',norm12)
fprintf(' ||G21 - G(nz+1:no,1:nw)||    = %d\n',norm21)
fprintf(' ||G22 - G(nz+1:no,nw+1:ni)|| = %d\n',norm22)
fprintf('\n')

% 6) MAT objects & ios numeric
sysM = lti2mat(sys);
ny = 1; nz = no - ny;
nu = 2; nw = ni - nu;
[a,b1,b2,c1,c2,d11,d12,d21,d22] = parsysdata(sys,[ny nu]);

sys11 = ss(a,b1,c1,d11); norm(sys11 - sys(1:nz,1:nw),inf); 
sys12 = ss(a,b2,c1,d12); norm(sys12 - sys(1:nz,nw+1:ni),inf);
sys21 = ss(a,b1,c2,d21); norm(sys21 - sys(nz+1:no,1:nw),inf);
sys22 = ss(a,b2,c2,d22); norm(sys22 - sys(nz+1:no,nw+1:ni),inf);

fprintf('Checking with lmitool object\n')
fprintf(' ||G11 - G(1:nz,1:nw)||       = %d\n',norm11)
fprintf(' ||G12 - G(1:nz,nw+1:ni)||    = %d\n',norm12)
fprintf(' ||G21 - G(nz+1:no,1:nw)||    = %d\n',norm21)
fprintf(' ||G22 - G(nz+1:no,nw+1:ni)|| = %d\n',norm22)
fprintf('\n')

% 7) flags
fprintf('Checking flags\n')

fprintf('  Dimensions\n')
[ns2,nz2,nw2,ny2,nu2] = parsysdata(sys,ios,'dims');
fprintf('    ns: expected = %2.0f, obtained = %2.0f\n',ns,ns2)
fprintf('    nz: expected = %2.0f, obtained = %2.0f\n',nz,nz2)
fprintf('    nw: expected = %2.0f, obtained = %2.0f\n',nw,nw2)
fprintf('    ny: expected = %2.0f, obtained = %2.0f\n',ny,ny2)
fprintf('    nu: expected = %2.0f, obtained = %2.0f\n',nu,nu2)

fprintf('\n')
fprintf('  Matrices\n')
a   = parsysdata(sys,ios,'a');      nr = norm(a - sys.A);
fprintf('   ||A - sys.A||                    = %d\n',nr)
b1  = parsysdata(sys,ios,'b1');     nr = norm(b1 - sys.B(:,1:nw));
fprintf('   ||B1 - sys.B(:,1:nw)||           = %d\n',nr)
b2  = parsysdata(sys,ios,'b2');     nr = norm(b2 - sys.B(:,nw+1:end));
fprintf('   ||B2 - sys.B(:,nw+1:end)||       = %d\n',nr)
c1  = parsysdata(sys,ios,'c1');     nr = norm(c1 - sys.C(1:nz,:));
fprintf('   ||C1 - sys.C(1:nz,:)||           = %d\n',nr)
c2  = parsysdata(sys,ios,'c2');     nr = norm(c2 - sys.C(nz+1:end,:));
fprintf('   ||C1 - sys.C(nz+1:end,:)||       = %d\n',nr)
d11 = parsysdata(sys,ios,'d11');    nr = norm(d11 - sys.D(1:nz,1:nw));
fprintf('   ||D11 - sys.D(1:nz,1:nw)||       = %d\n',nr)
d12 = parsysdata(sys,ios,'d12');    nr = norm(d12 - sys.D(1:nz,nw+1:ni));
fprintf('   ||D12 - sys.D(1:nz,nw+1:ni)||    = %d\n',nr)
d21 = parsysdata(sys,ios,'d21');    nr = norm(d21 - sys.D(nz+1:no,1:nw));  
fprintf('   ||D21 - sys.D(nz+1:no,1:nw)||    = %d\n',nr)
d22 = parsysdata(sys,ios,'d22');    nr = norm(d22 - sys.D(nz+1:no,nw+1:ni));
fprintf('   ||D22 - sys.D(nz+1:no,nw+1:ni)|| = %d\n',nr)


