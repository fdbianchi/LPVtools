function y = cat(varargin)
% CAT (overloaded)

switch varargin{1}
    case 1
        y = vertcat(varargin{2:end});
    case 2
        y = horzcat(varargin{2:end});
    otherwise
        error('PGSS:cat:inputError','dim > 2 not implemented')
end