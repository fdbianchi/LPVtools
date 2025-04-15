function [pdG,ny,nu] = setIOctrl(pdG, varargin)

% *** for internal use ***
%
% SETIOCTRL defines the InputGroup and OutputGroup properties according to
% the measure and control signals
%
% [pdG,ny,nu] = SETIOCTRL(pdG, r) returns a new pdG model where:
%               - pdG.OutputGroup.meas = nz+1:no    (no = nz + r(1))
%               - pdG.InputGroup.ctrl  = nw+1:ni    (ni = nw + r(2)) 
%               when r is a scalar (r=nu) or r=[0 nu], state-feedback is
%               assumed.
%
% [pdG,ny,nu] = SETIOCTRL(pdG, meas, ctrl) returns a new pdG model where
%               - pdG.OutputGroup.meas = meas
%               - pdG.InputGroup.ctrl  = ctrl
%               meas and ctrl are numeric indices or signal names 
%               When meas is empty or zero, state-feedback is assumed.

% fbianchi - 2020-06-30

% getting signal indices
if (nargin == 2)
    % vector input
    ioCtrl = varargin{1};

    if (isnumeric(ioCtrl) && isscalar(ioCtrl))
        % case ioCtrl is a scalar with # of controls (state feedback)
        
        nu = ioCtrl;
        ny = 0;
        [no,ni] = iosize(pdG);
        nw = ni - nu;

        meas_ind = 1:no;
        ctrl_ind = nw+1:ni;

    elseif (isnumeric(ioCtrl) && all(size(ioCtrl) == [1 2]))
        % case ioCtrl is a vector with # of measures and # of controls

        nu = ioCtrl(2);
        ny = ioCtrl(1);
        [no,ni] = iosize(pdG);
        nz = no - ny;
        nw = ni - nu;

        if (ny == 0)
            meas_ind = 1:no;
        else            
            meas_ind = nz+1:no;
        end
        
        ctrl_ind = nw+1:ni;
        
    else
        error('SETIOCTRL:inputError',...
            'R must be a numeric vector of dimension 1 or 2')
        
    end
    
elseif (nargin == 3)
    % vectors with meas & ctrl indices or names
    
    [no,ni] = iosize(pdG);
    
    % measured outputs
    meas = varargin{1};
    if isempty(meas)
        % state-feedback case
        meas_ind = 1:no;
        ny = 0;
        
    elseif isnumeric(meas)
        % numeric specifications
        measUn = unique(meas,'stable');
        if length(measUn) ~= length(meas)
            warning('SETIOCTRL:inputWarn',...
                'The repeated measure signals will be ignored')
        end
        ny = length(measUn);
        if (max(measUn) > no) && (min(measUn) < 0)
            error('SETIOCTRL:inputError',...
                'Measure indices exceed the number of pdG inputs')
        elseif (measUn == 0)
            meas_ind = 1:no;
            ny = 0;
        else
            meas_ind = measUn;
        end
        
    elseif iscellstr(meas) || ischar(meas)
        % list of names
        if iscellstr(meas)
            measUn = unique(meas,'stable');
            ny = length(measUn);
        else
            measUn = meas;
            ny = 1;
        end            
        if length(measUn) ~= length(meas)
            warning('SETIOCTRL:inputWarn',...
                'The repeated measure signals will be ignored')
        end
        i_names = unique(pdG.y,'stable');
        meas_ind = find(ismember(i_names,measUn))';
        if isempty(meas_ind) || (length(meas_ind) < ny)
            error('SETIOCTRL:inputError',...
                'At least one measure name is not an input of pdG')
        end
        
    else
        error('SETIOCTRL:inputError',...
            'MEAS must be vector with signal indices or names')
    end

    % control inputs
    ctrl = varargin{2};
    if isnumeric(ctrl) && isvector(ctrl)
        ctrlUn = unique(ctrl);
        if length(ctrlUn) ~= length(ctrl)
            warning('SETIOCTRL:inputWarn',...
                'The repeated control signals will be ignored')
        end
        nu = length(ctrlUn);
        if (max(ctrlUn) > ni) && (min(ctrlUn) <= 0)
            error('SETIOCTRL:inputError',...
                'Control indices exceed the number of pdG inputs')
        else
            ctrl_ind = ctrlUn;
        end
        
    elseif iscellstr(ctrl) || ischar(ctrl)
        % list of names
        if iscellstr(ctrl)
            ctrlUn = unique(ctrl,'stable');
            nu = length(ctrlUn);
        else
            ctrlUn = ctrl;
            nu = 1;
        end
        if length(ctrlUn) ~= length(ctrl)
            warning('SETIOCTRL:inputWarn',...
                'The repeated control signals will be ignored')
        end
        i_names = unique(pdG.u,'stable');
        ctrl_ind = find(ismember(i_names,ctrlUn))';
        if isempty(ctrl_ind) || (length(ctrl_ind) < nu)
            error('SETIOCTRL:inputError',...
                'At least one control name is not an input of pdG')
        end
        
    else
        error('SETIOCTRL:inputError',...
            'CTRL must be vector with signal indices or names')
    end
    
else
    error('SETIOCTRL:inputError',...
        'Invalid number of input argument')
    
end

% new model
pdG.InputGroup.ctrl = ctrl_ind;
pdG.OutputGroup.meas = meas_ind;

