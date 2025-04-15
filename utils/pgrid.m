function grd = pgrid(range,np)

% PGRID(range,np) returns a grid of points given the range and the number 
%   of points
%
% Use:
%   grd = PGRID(range,np)
%
% where:
%   - range:    (n x 2) numeric matrix with n number of parameters
%   - np:       scalar or vector of dimension n with the number of points in each
%               parameter direction
%   - grd:      (n x np) numeric matrix, each column corresponds to one
%               point in the grid
%

% fbianchi - 2020-02-19

% error checking
if (~isnumeric(range) || size(range,2)~=2)
    error('PGRID:inputError','RANGE must be a numeric matrix with two columns')
end

% number of rows
nr = size(range,1);

% default values
if (nargin < 2)
    np = 2;
elseif (~isvector(np))
    error('PGRID:inputError','NP must be a numeric vector')
elseif (~isscalar(np) && length(np) > nr)
    error('PGRID:inputError',...
        'NP must be a scalar or a numeric vector of lenght %3.0f',nr)
end
    
% sorting the range just in case
range = sort(range,2);

% in case np is a scalar, the grid will have the same number of points in 
% each direction
if length(np) == 1
    np = np*(diff(range,1,2) > 0) + (diff(range,1,2) == 0);
end

% pre-allocate memory
grd(nr,prod(np)) = 0;

% creates the grid
row = linspace(range(1,1),range(1,2),np(1));
grd(1,:) = repmat(row,1,prod(np(2:end)));
for ii = 2:nr
    row = linspace(range(ii,1),range(ii,2),np(ii));
    tmp = row(ones(1,prod(np(1:ii-1))),:); 
    row = tmp(:)';
    grd(ii,:) = repmat(row,1,prod(np(ii+1:end)));
end

