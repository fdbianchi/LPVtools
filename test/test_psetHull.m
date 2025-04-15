
% -------------------------------------------------------------------------
% How to use and Examples for PSET.HULL
%
% -------------------------------------------------------------------------

% fbianchi - 2021-03-26

% *****************
clc; 
clear all; 
close all
% *****************

rng(23656);

%% -------------------------------------------------------------------------
% 1d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 1 parameter\n\n')


points = rand(1,13);

% only points
setHull1a = pset.Hull(points);
disp('Set defition with points only')
disp(setHull1a)
% points + rate
setHull1b = pset.Hull(points,[-1 1]);
disp('Set defition with points and rates')
disp(setHull1b)
% point + rate + names
setHull1c = pset.Hull(points,[-1 1],'q');
disp('Set defition with points, rates and names')
disp(setHull1c)

% check if the set is empty
bool = isempty(setHull1a);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setHull1a);
if (np == 1) && (nv == 2)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = 0.8;    % in the set
[bool1,~] = checkval(setHull1a,p1);
p2 = -1.2;    % out the set
[bool2,~] = checkval(setHull1a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setHull1b)

% testing cvxdec
p1 = 0.5;
[alpha,idx]=cvxdec(setHull1b,p1);
ps = setHull1b.points(:,idx);
pc = ps*alpha;
plot(ps(1,:),[0 0],'bs')
plot(p1(1,:),0,'rx')
plot(pc(1,:),0,'ms')
if all(abs(p1 - pc) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


%% -------------------------------------------------------------------------
% 2d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 2 parameters\n\n')

points = rand(2,13);

% only points
setHull2a = pset.Hull(points);
disp('Set defition with range only')
disp(setHull2a)
% points + rate
setHull2b = pset.Hull(points,[-1 1;-1 1]);
disp('Set defition with ranges and rates')
disp(setHull2b)
% point + rate + names
setHull2c = pset.Hull(points,[-1 1],{'q','tr'});
disp('Set defition with ranges, rates and names')
disp(setHull2c)

% check if the set is empty
bool = isempty(setHull2a);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setHull2a);
if (np == 2) && (nv == 7)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = [0.4;0.5];    % in the set
[bool1,msg1] = checkval(setHull2a,p1);
p2 = [0.4;1];    % out the set
[bool2,msg2] = checkval(setHull2a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setHull2b)

% testing cvxdec
p1 = [0.5;0.4];
[alpha,idx] = cvxdec(setHull2b,p1);
ps = setHull2b.points(:,idx);
pc = ps*alpha;
plot(ps(1,:),ps(2,:),'bs')
plot(p1(1,:),p1(2,:),'rx')
plot(pc(1,:),pc(2,:),'ms')
if all(abs(p1 - pc) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

%% -------------------------------------------------------------------------
% 3d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.BOX with 3 parameters\n\n')

points = rand(3,13);

% definition with points
setHull3 = pset.Hull(points);
disp('Set defition with range only')
disp(setHull3)

% check if the set is empty
bool = isempty(setHull3);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setHull3);
if (np == 3) && (nv == 10)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = [0.5;0.4;0.6];    % in the set
[bool1,msg1] = checkval(setHull3,p1);
p2 = [1.0;0.4;0.6];    % out the set
[bool2,msg2] = checkval(setHull3,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% testing cvxdec
p = [0.5;0.4;0.6];
[alpha,idx] = cvxdec(setHull3,p);
ps = setHull3.points(:,idx);
pc = ps*alpha;
if all(abs(p1 - pc) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

% plot
figure
plot(setHull3)
plot3(ps(1,:),ps(2,:),ps(3,:),'bs')
plot3(p(1,:),p(2,:),p(3,:),'rx')
plot3(pc(1,:),pc(2,:),pc(3,:),'ms')


% ------------------------------------------------------------------------
% subset

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Subsets\n')

disp('Original set')
P = [1 1 4 3;2 2 3 4;2 3 7 8];
setHull8 = pset.Hull(P);
disp(setHull8)

disp('Subset with parameters 1 and 3')
setHull8a = setHull8([1 2]);
disp(setHull8a)

figure
plot(setHull8)
hold on
xa = setHull8a.points(1,:);
ya = setHull8a.points(2,:);
patch(xa,ya,P(3,1)*[1 1 1],0.5*[1 1 1],'EdgeColor',0.5*[1 1 1],'FaceAlpha',0.3)


%% ------------------------------------------------------------------------
% extreme cases

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Extreme cases\n')

fprintf('\nEmpty sets\n')
setHull4 = pset.Box();
bool = isempty(setHull4);
if (bool == 1)
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end
[np,nv] = size(setHull4);
if (np == 0) && (nv == 0)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end
[bool1,msg1] = checkval(setHull4,[]);
[bool2,msg2] = checkval(setHull4,1);
if (bool1 == 1) && (bool2 < 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
try 
    [alpha,idx] = cvxdec(setHull4,1);
catch ME
    disp(ME.message)
    disp('CVXDEC working correctly')
end

fprintf('\nSet with the same vertices\n')
setHull5 = pset.Hull(5);
[bool1,msg1] = checkval(setHull5,5);
[bool2,msg2] = checkval(setHull5,1.5);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
p1 = 5;
[alpha,idx] = cvxdec(setHull5,p1);
p1r = setHull5.points(:,idx)*alpha;
if all(abs(p1 - p1r) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

fprintf('\nSet with the some coincident vertices\n')
setHull6 = pset.Hull([1;3;5]);
p1 = [1;3;5];  
[bool1,msg1] = checkval(setHull6,p1);
p2 = [3;1.5;4]; 
[bool2,msg2] = checkval(setHull6,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
[alpha,idx] = cvxdec(setHull6,p1);
p1r = setHull6.points(:,idx)*alpha;
if all(abs(p1 - p1r) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


fprintf('\nChecking plots\n')

setHull7 = pset.Hull([1 1]);
figure; plot(setHull7); title('Set 1D with coincident vertices')

setHull7 = pset.Hull([1;2]);
figure; plot(setHull7); title('Set 2D with 1 vertex')
setHull7 = pset.Hull([2 1;2 2]);
figure; plot(setHull7); title('Set 2D with 2 vertices and 2 coincident coord')
setHull7 = pset.Hull([1 1;2 2]);
figure; plot(setHull7); title('Set 2D with 2 coincident vertices')

setHull7 = pset.Hull([1;2;5]);
figure; plot(setHull7); title('Set 3D with 1 vertex')
setHull7 = pset.Hull([0 1;2 2;3 5]);
figure; plot(setHull7); title('Set 3D with 2 vertices')
setHull7 = pset.Hull([0 1;2 2;5 5]);
figure; plot(setHull7); title('Set 2D with 2 vertices and 2 coincident coord')
setHull7 = pset.Hull([1 1;2 2;5 5]);
figure; plot(setHull7); title('Set 2D with 2 coincident vertices')
