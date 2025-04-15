% -------------------------------------------------------------------------
% How to use and Examples for PSET.GRID
%
% -------------------------------------------------------------------------

% fbianchi - 2021-03-26

% *****************
clc; 
clear all; 
close all
% *****************

%% ------------------------------------------------------------------------
% 1d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 1 parameter\n')

% range
p1min  = 0.5;   p1max =   4;
prange = [p1min p1max];
prate  = [-1 1];
% Names
pnames = {'p'};

% 1) given only the parameter range
setGrid1a = pset.Grid(prange,3,prate);
disp('Set defition with range, # points and rate')
disp(setGrid1a)
% 2) given a subset for each parameter
subsets{1} = linspace(p1min,p1max,3);
setGrid1b = pset.Grid(subsets,prate);
disp('Set defition with grid and rate')
disp(setGrid1b)
% 3) idem 2) and given parameter names
setGrid1c = pset.Grid(subsets,prate,pnames);
disp('Set defition with grid, rate and names')
disp(setGrid1c)

% check if the set is empty
bool = isempty(setGrid1a);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setGrid1a);
if (np == 1) && (nv == 3)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 =  1;    % in the set
[bool1,msg1] = checkval(setGrid1b,p1);
p2 = -1;    % out the set
[bool2,msg2] = checkval(setGrid1b,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setGrid1b)
plot(p1(1),0,'r*')
plot(p2(1),0,'m*')

% convex decomposition
[alpha,idx] = cvxdec(setGrid1a,p1);
p3 = setGrid1a.points(:,idx)*alpha;
plot(p3(1),0,'ks')
p3s = setGrid1a.points(:,idx([1 2]));
plot(p3s(1,:),0,'k-')
if all(abs(p1 - p3) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end



%% ------------------------------------------------------------------------
% 2d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 2 parameters\n')

% range
p1min  = 0.5;   p1max =   4;
p2min  = 0.0;   p2max = 100;
prange = [p1min p1max; p2min p2max];
prate  = [-1 1];
% Names
pnames = {'p','q'};

% 1) given only the parameter range
setGrid2a = pset.Grid(prange,[3 4],prate);
disp('Set defition with range, # points and rate')
disp(setGrid2a)
% 2) given a subset for each parameter
subsets{1} = linspace(p1min,p1max,3);
subsets{2} = linspace(p2min,p2max,4);
setGrid2b = pset.Grid(subsets,prate);
disp('Set defition with grid and rate')
disp(setGrid2b)
% 3) idem 2) and given parameter names
setGrid2c = pset.Grid(subsets,prate,pnames);
disp('Set defition with grid, rate and names')
disp(setGrid2c)

% check if the set is empty
bool = isempty(setGrid2a);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setGrid2a);
if (np == 2) && (nv == 12)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = [1;30];     % in the set
[bool1,msg1] = checkval(setGrid2a,p1);
p2 = [-1;30];    % out the set
[bool2,msg2] = checkval(setGrid2a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setGrid2b)
plot(p1(1),p2(2),'r*')
plot(p2(1),p2(2),'m*')

% convex decomposition
[alpha,idx] = cvxdec(setGrid2b,p1);
p3 = setGrid2b.points(:,idx)*alpha;
plot(p3(1),p3(2),'ks')
p3s = setGrid2b.points(:,idx([1 2 4 3 1]));
plot(p3s(1,:),p3s(2,:),'k-')
if all(abs(p1 - p3) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


%% ------------------------------------------------------------------------
% 3d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 3 parameters\n')

% given a subset for each parameter
subsets{1} = linspace(1,3,3);
subsets{2} = linspace(1,6,4);
subsets{3} = linspace(-2,2,4);
setGrid3a = pset.Grid(subsets);
disp('Set defition with grid')
disp(setGrid3a)

% check if the set is empty
bool = isempty(setGrid3a);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setGrid3a);
if (np == 3) && (nv == 3*4*4)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = [2 3 0]';     % in the set
[bool1,msg1] = checkval(setGrid3a,p1);
p2 = [2 7 0]';    % out the set
[bool2,msg2] = checkval(setGrid3a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setGrid3a)
plot3(p1(1),p1(2),p1(3),'r*')
plot3(p2(1),p2(2),p2(3),'m*')

[alpha,idx] = cvxdec(setGrid3a,p1);
p3 = setGrid3a.points(:,idx)*alpha;
plot(p3(1),p3(2),'ks')
if all(abs(p1 - p3) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

% ------------------------------------------------------------------------
% subset

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Subsets\n')

disp('Original set')
disp(setGrid3a)

disp('Subset with parameters 1 and 3')
setGrid4 = setGrid3a([1 3]);
disp(setGrid4)


%% ------------------------------------------------------------------------
% extreme cases

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Extreme cases\n')

fprintf('\nEmpty sets\n')
setGrid5 = pset.Grid();
bool = isempty(setGrid5);
if (bool == 1)
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end
[np,nv] = size(setGrid5);
if (np == 0) && (nv == 0)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end
[bool1,msg1] = checkval(setGrid5,[]);
[bool2,msg2] = checkval(setGrid5,1);
if (bool1 == 1) && (bool2 < 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
try 
    [alpha,idx] = cvxdec(setGrid5,1);
catch ME
    disp(ME.message)
    disp('CVXDEC working correctly')
end

fprintf('\nSet with the same vertices\n')
setGrid6 = pset.Grid([1 1]);
[bool1,~] = checkval(setGrid6,1);
[bool2,~] = checkval(setGrid6,1.5);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
p1 = 1;
[alpha,idx] = cvxdec(setGrid6,p1);
p1r = setGrid6.points(:,idx)*alpha;
if all(abs(p1 - p1r) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

fprintf('\nSet with the some coincident vertices\n')
setGrid8 = pset.Grid([2 5;1 1;2 9]);
p1 = [3;1;4];  
[bool1,msg1] = checkval(setGrid8,p1);
p2 = [3;1.5;4]; 
[bool2,msg2] = checkval(setGrid8,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
[alpha,idx] = cvxdec(setGrid8,p1);
p1r = setGrid8.points(:,idx)*alpha;
if all(abs(p1 - p1r) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


fprintf('\nChecking plots\n')

setGrid9 = pset.Grid([1 1]);
figure; plot(setGrid9); title('Set 1D with coincident vertices')

setGrid9 = pset.Grid([2 2]);
figure; plot(setGrid9); title('Set 2D with 1 vertex')
setGrid9 = pset.Grid([2 2;2 2]);
figure; plot(setGrid9); title('Set 2D with 2 coincident vertices')
setGrid9 = pset.Grid([1 2;2 2]);
figure; plot(setGrid9); title('Set 2D with 2 vertices, with 1 coord coincident')

setGrid9 = pset.Grid([1 1;2 2;5 5]);
figure; plot(setGrid9); title('Set 3D with 1 vertex')
setGrid9 = pset.Grid([0 1;2 2;3 5]);
figure; plot(setGrid9); title('Set 3D with 4 coincident vertices')
setGrid9 = pset.Grid([0 1;2 2;5 5]);
figure; plot(setGrid9); title('Set 2D with 12 coincident vertices')
setGrid9 = pset.Grid([1 1;2 2;5 5]);
figure; plot(setGrid9); title('Set 2D with 16 coincident vertices')


