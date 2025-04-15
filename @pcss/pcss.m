classdef pcss
    
    % PCSS implements an online computed controller
    %
    % PCSS is used by lpvsyn to store all the information needed to
    % undo the change of controller variables:
    %
    %   ak = N\(Ak - X*(A-B2*Dk*C2)*Y - Bk*C2*Y - X*B2*Ck)/Mt;
    %   bk = N\(Bk - X*B2*Dk);
    %   ck = (Ck - Dk*C2*Y)/Mt;
    %   dk = Dk;
    %
    % when the Lyapunov functions are parameter dependent
    %
    % To create a PCSS object, use:
    %
    %   pdK = pcss(pdKhat, pdG22, pdY, pdX);
    %
    % where 
    %   - pdKhat: a pass/ppss/pcss object with the auxiliary controller
    %           variables
    %   - pdG22: a pass/ppss/pcss with the plant pdG(meas,ctrl)
    %   - xFcn: a pass/ppss/pcss with X Lyapunov function (only matrix D)
    %   - yFcn: a pass/ppss/pcss with Y Lyapunov function (only matrix D)
    %
    %
    % PCSS Methods:
    %   pcss       - class constructor
    %   iosize     - returns the number of input and outputs
    %   order      - returns the number of states
    %   size       - returns the model dimensions
    %   isempty    - checks if a empty object
    %   subs       - evaluates an LPV model at a frozen parameter values
    %   ss         - converts to ss-object (tf & zpk also available)
    %
    % other overloaded functions: 
    %   eig, pole, tzero, dcgain
    %   bode, bodemag, step, etc.
    %
    % See also ss, pass, ppss, pgss
    
    % fbianchi - 2021-07-02
    
    properties
        ctrller     % controller matrices before variable change 
                    %   [pass/ppss/pgss object]
        plant       % plant pdG(meas,ctrl), [pass/ppss/pgss object]
        xFcn        % X Lyapunov function, [pass/ppss/pgss object]
        yFcn        % y Lyapunov function, [pass/ppss/pgss object]
        
    end
    
    methods
        
        function obj = pcss(varargin)
            
            % PCSS implements an online computed controller
            %
            % PCSS is used by lpvsyn to store all the information needed to
            % undo the change of controller variables:
            %
            %   ak = N\(Ak - X*(A-B2*Dk*C2)*Y - Bk*C2*Y - X*B2*Ck)/Mt;
            %   bk = N\(Bk - X*B2*Dk);
            %   ck = (Ck - Dk*C2*Y)/Mt;
            %   dk = Dk;
            %
            % when the Lyapunov functions are parameter dependent
            %
            % To create a PCSS object, use:
            %
            %   pdK = pcss(pdKhat, pdG22, pdY, pdX);
            %
            % where 
            %   - pdKhat: a pass/ppss/pcss object with the auxiliary controller
            %           variables
            %   - pdG22: a pass/ppss/pcss with the plant pdG(meas,ctrl)
            %   - xFcn: a pass/ppss/pcss with X Lyapunov function (only matrix D)
            %   - yFcn: a pass/ppss/pcss with Y Lyapunov function (only matrix D)

            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return

            % ------------------------------------------------------------
            % to manage the copy of objects
            elseif isa(varargin{1},'pcss')
                obj = varargin{1};            
            
           
            % ------------------------------------------------------------
            % pdK = pcss(pdKhat,pdG22,pdY,pdX)
            else
                
                % controller matrices before change of variable
                if isa(varargin{1},'p_ss')
                    obj.ctrller = varargin{1};
                    nk = order(obj.ctrller);
                    [nu,ny] = iosize(obj.ctrller);
                    if (nk == 0)
                        nk = ny;    % for state-feedback case
                    end

                else
                    error('PCSS:PCSS:inputError',....
                        'ctrller must be an LPV model (pass, ppss or pgss)')

                end
                
                % plant model channel (2,2)
                if ~isempty(varargin{2})
                    if isa(varargin{2},'p_ss')
                        
                        obj.plant = varargin{2};
                        
                        ns = order(obj.plant);
                        if (ns ~= nk)
                            error('PCSS:PCSS:inputError',....
                                'The order of PLANT must be %d',nk)
                        end
                        [no,ni] = iosize(obj.plant);
                        if (no ~= ny)
                            error('PCSS:PCSS:inputError',....
                                'The number of outputs of PLANT must be %d',ny)
                        end
                        if (ni ~= nu)
                            error('PCSS:PCSS:inputError',....
                                'The number of inputs of PLANT must be %d',nu)
                        end
                        
                    else
                        error('PCSS:PCSS:inputError',....
                            'plant must be an LPV model (pass, ppss or pgss)')
                        
                    end
                end
            
                % Y Lyapunov function
                if (isa(varargin{3},'pgss') || isa(varargin{3},'ppss') ...
                         || isnumeric(varargin{3}))
                    obj.yFcn = varargin{3};

                    if (size(obj.yFcn,1) ~= nk)
                        error('PCSS:PCSS:inputError',....
                            'Y must be a matrix of %d x %d',nk,nk)
                    end

                else
                    error('PCSS:PCSS:inputError',....
                        'Y must be a 2D matrix or an LPV model (pass or ppss)')

                end
                
                % X Lyapunov function
                if (nargin > 3)
                    if (isa(varargin{4},'pgss') || isa(varargin{4},'ppss') ...
                            || isnumeric(varargin{4}))
                        
                        obj.xFcn = varargin{4};
                        if (size(obj.xFcn,1) ~= nk)
                            error('PCSS:PCSS:inputError',....
                                'X must be a matrix of %d x %d',nk,nk)
                        end
                        
                    else
                        error('PCSS:PCSS:inputError',....
                            'X must be a 2D matrix or an LPV model (pass or ppss)')
                        
                    end
                end
            end
        end
     
        function pv = parset(obj)
            
            % PARSET(obj) returns the parameter set 
            
            pv = obj.ctrller.parset;
        end
        
    end
    
end

