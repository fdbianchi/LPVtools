
% -------------------------------------------------------------------------
% How to use and Examples for PSET.BOX
%
% -------------------------------------------------------------------------

% fbianchi - 2021-03-29

% *****************
clc; 
clear all; 
close all
% *****************

%% ------------------------------------------------------------------------
% 1d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 1 parameter\n\n')

% Parameter range
prange = [0.5 4];
% Parameter rate
prate  = [-1 2];
% Parameter names
pnames = {'p'};
% Parameter set definition:
% 1) given only the parameter range
setBox1a = pset.Box(prange);
disp('Set defition with range only')
disp(setBox1a)
% 2) given the parameter range and the rate limits
setBox1b = pset.Box(prange,prate);
disp('Set defition with range & rate')
disp(setBox1b)
% 3) idem 2) and given parameter names
setBox1c = pset.Box(prange,prate,pnames);
disp('Set defition with range, rate and names')
disp(setBox1c)
% 4) idem 3) without rate limits for all parameters
setBox1d = pset.Box(prange,[],pnames);
disp('Set defition with range and names only')
disp(setBox1d)

% check if the set is empty
t = isempty(setBox1a);
if ~t
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setBox1a);
if (np == 1) && (nv == 2)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 =  1;    % in the set
[bool1,msg1] = checkval(setBox1a,p1);
p2 = -1;    % out the set
[bool2,msg2] = checkval(setBox1a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setBox1a)
plot(p1(1),0,'r*')
plot(p2(1),0,'m*')

% convex decomposition
[alpha,idx] = cvxdec(setBox1a,p1);
p3 = setBox1a.points(:,idx)*alpha;
plot(p3(1),0,'ks')
if all(abs(p1 - p3) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


%% ------------------------------------------------------------------------
% 2d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 2 parameters\n\n')

% Parameter range
prange = [0.5 4; -5 6];
% Parameter rate
prate  = [-1 1; -2 2];
% Parameter names
pnames = {'p','q'};

% Parameter set definition:
setBox2a = pset.Box(prange);
disp('Set defition with range only')
disp(setBox2a)
% 2) given the parameter range and the rate limits
setBox2b = pset.Box(prange,prate);
disp('Set defition with range & rate')
disp(setBox2b)
% 3) idem 2) and given parameter names
setBox2c = pset.Box(prange,prate,pnames);
disp('Set defition with range, rate and names')
disp(setBox2c)
% 4) idem 3) without rate limits for all parameters
setBox2d = pset.Box(prange,[],pnames);
disp('Set defition with range and names only')
disp(setBox2d)

% check if the set is empty
t = isempty(setBox2a);
if ~t
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setBox2a);
if (np == 2) && (nv == 4)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = [ 1; 3];    % in the set
[bool1,msg1] = checkval(setBox2a,p1);
p2 = [-1; 3];    % out the set
[bool2,msg2] = checkval(setBox2a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setBox2a)
plot(p1(1),p2(2),'r*')
plot(p2(1),p2(2),'m*')

% convex decomposition
[alpha,idx] = cvxdec(setBox2a,p1);
p3 = setBox2a.points(:,idx)*alpha;
plot(p3(1),p3(2),'ks')
if all(abs(p1 - p3) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


%% ------------------------------------------------------------------------
% 3d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 3 parameters\n')

% Parameter set definition:
setBox3a = pset.Box([1 3;1 6;-2 2]);
disp('Set defition with range only')
disp(setBox2a)

% check if a value is in the set
p1 = [2 3 0]';    % in the set
[bool1,msg1] = checkval(setBox3a,p1);
p2 = [2 7 0]';    % out the set
[bool2,msg2] = checkval(setBox3a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setBox3a)
plot3(p1(1),p1(2),p1(3),'r*')
plot3(p2(1),p2(2),p2(3),'m*')

[alpha,idx] = cvxdec(setBox3a,p1);
p3 = setBox3a.points(:,idx)*alpha;
plot(p3(1),p3(2),'ks')
if all(abs(p1 - p3) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


%% ------------------------------------------------------------------------
% conversion from lmitool object

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Conversion from lmitool object\n')

% from box lmitool obj to pset
setLMI1 = pvec('box',prate, prate);
setBox4 = pset.Box('pvec',setLMI1);
disp('from PVEC to PSET.BOX')
disp(setBox4)

setLMI2 = pvec(setBox4);
disp('from PSET.BOX to PVEC')
pvinfo(setLMI2)


%% ------------------------------------------------------------------------
% subset

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Subsets\n')

disp('Original set')
disp(setBox3a)

disp('Subset with parameters 1 and 3')
setBox5 = setBox3a([1 3]);
disp(setBox5)

%% ------------------------------------------------------------------------
% extreme cases

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Extreme cases\n')

fprintf('\nEmpty sets\n')
setBox6 = pset.Box();
bool = isempty(setBox6);
if (bool == 1)
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end
[np,nv] = size(setBox6);
if (np == 0) && (nv == 0)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end
[bool1,msg1] = checkval(setBox6,[]);
[bool2,msg2] = checkval(setBox6,1);
if (bool1 == 1) && (bool2 < 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
try 
    [alpha,idx] = cvxdec(setBox6,1);
catch ME
    disp(ME.message)
    disp('CVXDEC working correctly')
end


fprintf('\nSet with the same vertices\n')
setBox7 = pset.Box([1 1]);
[bool1,msg1] = checkval(setBox7,1);
[bool2,msg2] = checkval(setBox7,1.5);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
p1 = 1;
[alpha,idx] = cvxdec(setBox7,p1);
p1r = setBox7.points(:,idx)*alpha;
if all(abs(p1 - p1r) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

fprintf('\nSet with the some coincident vertices\n')
setBox8 = pset.Box([2 5;1 1;2 9]);
p1 = [3;1;4]; 
[bool1,msg1] = checkval(setBox8,p1);
p2 = [3;1.5;4]; 
[bool2,msg2] = checkval(setBox8,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
[alpha,idx] = cvxdec(setBox8,p1);
p1r = setBox8.points(:,idx)*alpha;
if all(abs(p1 - p1r) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


fprintf('\nChecking plots\n')

setBox8 = pset.Box([1 1]);
figure; plot(setBox8); title('Set 1D with coincident vertices')

setBox8 = pset.Box([2 1;2 2]);
figure; plot(setBox8); title('Set 2D with 2 coincident vertices')
setBox8 = pset.Box([1 1;2 2]);
figure; plot(setBox8); title('Set 2D with 4 coincident vertices')

setBox8 = pset.Box([0 1;2 2;3 5]);
figure; plot(setBox8); title('Set 3D with 4 coincident vertices')
setBox8 = pset.Box([0 1;2 2;5 5]);
figure; plot(setBox8); title('Set 2D with 12 coincident vertices')
setBox8 = pset.Box([1 1;2 2;5 5]);
figure; plot(setBox8); title('Set 2D with 16 coincident vertices')
