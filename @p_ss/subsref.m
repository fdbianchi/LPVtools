function G = subsref(pdG,S)

% SUBSREF method for p_ss class
%
% G = pdG(idx) returns idx-th LTI model in the description
%
% G = pdG(ny,nu) returns the LPV model corresponding to the input nu and
%               output ny
%
% G = pdG(ny,nu,idx) returns idx-th LTI model corresponding to the 
%               input nu and output ny
%
% See also pass, ppss, pgss, pcss

% fbianchi - 2021-03-31
% fbianchi - 2024-12-13 - rev.

if (nargin == 1)
    G = pdG;
    return
end

switch S(1).type
    case '()'
        % use the third subscript to select an element of the system set    
        if (length(S(1).subs) >= 3)
            ind1 = S(1).subs{1};
            ind2 = S(1).subs{2};
            ind3 = S(1).subs{3};
            
        elseif (length(S(1).subs) == 2)
            ind1 = S(1).subs{1};
            ind2 = S(1).subs{2};
            ind3 = ':';
            
        elseif (length(S(1).subs) == 1)
            ind1 = ':';
            ind2 = ':';
            ind3 = S(1).subs{1};
            
        end
        
        if ~(isnumeric(ind3) || strcmp(ind3,':'))
            error('P_SS:SUBSREF:inputError',...
                'The number of system must be numeric')
            
        elseif (~strcmp(ind3,':') && (ind3 > nsys(pdG)))
            error('P_SS:SUBSREF:inputError',...
                'Index exceeds matrix dimensions.')
            
        end
        
        % ToDo: implement a more efficient way
%         % subsystem
%         warning('off','Control:ltiobject:RepeatedChannelNames')
%         G = ss(pdG.A(:,:,ind3),pdG.B(:,:,ind3),...
%                 pdG.C(:,:,ind3),pdG.D(:,:,ind3),...
%                 'StateName',pdG.StateName,...
%                 'InputName',pdG.InputName,...
%                 'OutputName',pdG.OutputName,...
%                 'InputGroup',pdG.InputGroup,...
%                 'OutputGroup',pdG.OutputGroup);
%         warning('on','Control:ltiobject:RepeatedChannelNames')
%         % input/output selection
%         G = G(ind1,ind2);
%         % if colon then the result is p_ss
%         if strcmp(ind3,':')
%             if isa(pdG,'pass')
%                 G = pass(G,pdG.parset);
%             elseif isa(pdG,'ppss')
%                 G = ppss(G,pdG.parset);
%             elseif isa(pdG,'pgss')
%                 G = pgss(G,pdG.parset,pdG.parfcn);
%             end                
%             
%             dgn = ispd(G);
%             if sum(dgn) == 0
%                 warning('The resulting model is not parameter dependent')
%             end
%         end
        
        % #### new implementation ####
        if strcmp(ind3,':')
            % if colon then the result is p_ss
            
            [no,ni] = size(pdG.D(:,:,1));
            % outputs
            if strcmp(ind1,':')
                ind1 = 1:no;
                
            elseif iscellstr(ind1) || ischar(ind1)
                [aux1,aux2] = ismember(ind1,pdG.OutputName);
                if any(~aux1)
                    % no in the outputs then search in the OutputGroup
                    aux1 = ismember(ind1,fields(pdG.OutputGroup));
                    if any(~aux1)
                        error('P_SS:SUBSREF:inputError',...
                            '%s is not a model output',ind1(~aux1))
                    else
                        ind1 = pdG.OutputGroup.(ind1);
                    end 
                else
                    ind1 = aux2;
                end
            elseif isnumeric(ind1)
                if any(ind1 <= 0) || any(ind1 > no)
                    error('P_SS:SUBSREF:inputError','First is out of range')
                else
                    lo = length(ind1); ind1 = reshape(ind1,1,lo);
                end
            else
                error('P_SS:SUBSREF:inputError',...
                    'First index invalid')
            end
            % inputs
            if strcmp(ind2,':')
                ind2 = 1:ni;
                
            elseif iscellstr(ind2) || ischar(ind2)
                [aux1,aux2] = ismember(ind2,pdG.InputName);
                if any(~aux1)
                    % no in the outputs then search in the InputGroup
                    aux1 = ismember(ind2,fields(pdG.InputGroup));
                    if any(~aux1)
                        error('P_SS:SUBSREF:inputError',...
                            '%s is not a model output',ind2(~aux1))
                    else
                        ind2 = pdG.InputGroup.(ind2);
                    end
                else
                    ind2 = aux2;
                end
            elseif isnumeric(ind2)
                if any(ind2 <= 0) || any(ind2 > ni)
                    error('P_SS:SUBSREF:inputError','First is out of range')
                else
                    lo = length(ind2); ind2 = reshape(ind2,1,lo);
                end
            else
                error('P_SS:SUBSREF:inputError',...
                    'Second index invalid')
            end
            
            A = pdG.A(:,:,ind3);
            B = pdG.B(:,ind2,ind3);
            C = pdG.C(ind1,:,ind3);
            D = pdG.D(ind1,ind2,ind3);
            
            if isa(pdG,'pass')
                G = pass(A,B,C,D,pdG.parset,...
                    'StateName',pdG.StateName,...
                    'InputName',pdG.InputName(ind2),...
                    'OutputName',pdG.OutputName(ind1));
            elseif isa(pdG,'ppss')
                G = ppss(A,B,C,D,pdG.parset,...
                    'StateName',pdG.StateName,...
                    'InputName',pdG.InputName(ind2),...
                    'OutputName',pdG.OutputName(ind1));
            elseif isa(pdG,'pgss')
                G = pgss(A,B,C,D,pdG.parset,pdG.parfcn,...
                    'StateName',pdG.StateName,...
                    'InputName',pdG.InputName(ind2),...
                    'OutputName',pdG.OutputName(ind1));
            end
            
            gSet = fields(pdG.InputGroup);
            for ii = 1:length(gSet)
                ind4 = pdG.InputGroup.(gSet{ii});
                [~,gIdx] = ismember(ind4,ind2);
                if any(gIdx)
                    G.InputGroup.(gSet{ii}) = gIdx(gIdx > 0);
                end
            end
            gSet = fields(pdG.OutputGroup);
            for ii = 1:length(gSet)
                ind4 = pdG.OutputGroup.(gSet{ii});
                [~,gIdx] = ismember(ind4,ind1);
                if any(gIdx)
                    G.OutputGroup.(gSet{ii}) = gIdx(gIdx > 0);
                end
            end

            dgn = ispd(G);
            if sum(dgn) == 0
                warning('The resulting model is not parameter dependent')
            end
        else
            
            % resulting model is LTI
            warning('off','Control:ltiobject:RepeatedChannelNames')
            G = ss(pdG.A(:,:,ind3),pdG.B(:,:,ind3),...
                pdG.C(:,:,ind3),pdG.D(:,:,ind3),...
                'StateName',pdG.StateName,...
                'InputName',pdG.InputName,...
                'OutputName',pdG.OutputName,...
                'InputGroup',pdG.InputGroup,...
                'OutputGroup',pdG.OutputGroup);
            warning('on','Control:ltiobject:RepeatedChannelNames')
            % input/output selection
            G = G(ind1,ind2);
        end
        
        % rest of subindices
        if (length(S) > 1)
            G = builtin('subsref',G,S(2:end));
        end
        
    case '.'
        % leave other subref without definition
        G = builtin('subsref',pdG,S);
        
    otherwise
        error('P_SS:SUBSREF:inputError','Indexing type not supported')

end

    
    

