function varargout = ParseVarargin(params,defparams,legalvalues,var_arg,char2logic)
%PARSEVARARGIN(PARAMS,DEFPARAMS,LEGALPARAMS,VAR_ARG)parses VAR_ARG and 
% returns a variable number of outputs (equal to the number of PARAMS) with
% the specified or default values (DEFPARMS).
% VAR_ARG must be constructed so that each property name is followed by its 
% property value (ie:,'param1',value1,'param2',value2).
%
%Performs several checks: 1) returns error if property names are not valid 
% (ie, are not among PARAMS). 2) [optinal] checks validity of property values
% (requires LEGALVALUES). 
% All input must be cells (see examples). 
%
% LEGALVALUES can be:
%   - empty -> no check on property value
%   - a cell array of strings to be compared to the parameter value
%   - a two cell element. First element is a function handle which has to return 
%          true or false depending on the parameter value. The second element
%          is the error message in case the function returned false.
%
%NB: varargout are forced to be lowercase if and only if the LEGALVALUE for
%   the parameters are provided (i.e., it is not empty). 
%
%If CHAR2LOGIC is provided (and ~= 0), the function performs automatic convertion of 
% chars to logical in case of: 'on','go','true' -> 1; 'off','stop','false' -> 0.
%
%Usage example:
%
%   varargin = {'tail','right',...
%               'zscore','true',...
%               'path_file', '/home/LAB_G1'};
%               'delta',42}
% 
%   params   = {'tail','zscore',   'path_file',  'delta'};
%   defparms = {'both',   'off', '/home/PIPPO',      3, };                          
%
%   legalvalues{1} = {'both','right','left'};     %forced to be lowercase
%   legalvalues{2} = {'true','false','off','on'}; %forced to be lowercase
%   legalvalues{3} = [];                          %NOT forced to be lowercase  
%   legalvalues{4} = {@(x) (sum(mod(x,1))==0),'Only integer values are allowed.'};
%
%   [tail,zscore_flag,delta] = parse_varargin(params,defparms,legalvalues,varargin)

%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com

debug_mode = 0; %use this flag to print complete error msg, including error lines

if nargin == 0
    help parse_varargin
    return
end

if nargin < 5
    char2logic = 0;
else
    if char2logic  % Automatic convertion of the following strings to logical values
        char2logical_true  = {'on','go','true'}; %these values will be converted to logical 1
        char2logical_false = {'off','stop','false'}; %these values will be converted to logical 0
        tobeconverted = [char2logical_true,char2logical_false];
    end
end

n_params = length(params);

%check valid property names
l = 1;
while l <= size(var_arg,2)
    check_property_name = sum(strcmpi(params, var_arg{l}));
    if ~check_property_name 
        %print valid property names
        str = [];
        for m = 1:n_params
            str = [str,' "',params{m},'"'];
        end
        msg = sprintf(['There is no "',var_arg{l},'" property.\nValid property names are:',str,'.\n']);
        throw_error(msg,debug_mode)
    end
    l = l +2;
end

for l = 1:n_params 
    inputExist = find(strcmpi(var_arg, params{l}));
    if inputExist
        %-------------check validity if desired-----------
        if ~isempty(legalvalues) && ~isempty(legalvalues{l}) 
            
            if iscell(legalvalues{l}) && isa(legalvalues{l}{1},'function_handle')
                %this is a function handle
                if ~legalvalues{l}{1}(var_arg{inputExist+1}) %evaluate the parameter
                    % result is false, print attached error message
                    msg = sprintf(['Invalid value for parameter: "',var_arg{inputExist},'"\n',legalvalues{l}{2},'\n']);
                    throw_error(msg,debug_mode)
                end
            else
                if iscell(legalvalues{l}) 
                    valid_check = sum(strcmpi(legalvalues{l}, var_arg{inputExist+1}));
                else % numbers must not be in cells
                    valid_check = sum(legalvalues{l} == var_arg{inputExist+1});
                end
                if valid_check == 0
                    if iscell(legalvalues{l}) 
                        str = ['\nLegal values are:'];
                        for m = 1:length(legalvalues{l})
                            str = [str,' "',legalvalues{l}{m},'"'];
                        end
                        str = [str,'.\n'];
                    else
                        str = ['\nLegal values are:'];
                        for m = 1:length(legalvalues{l})
                            str = [str,' ',num2str(legalvalues{l}(m))];
                        end
                        str = [str,'.\n'];
                    end
                    msg = sprintf(['Invalid value for parameter: "',var_arg{inputExist},'"',str]);
                    throw_error(msg,debug_mode)
                else % they are valid, force to be lowercase
                    var_arg{inputExist+1} = lower(var_arg{inputExist+1});
                end
            end 
        end
        %-------------------------------------------------
        %input exists and is valid
        if char2logic && sum(strcmpi(tobeconverted,var_arg{inputExist+1})) %convert to logical if appropriate
            varargout{l} = char2logical(var_arg{inputExist+1},char2logical_true,char2logical_false);
        else
            varargout{l} = var_arg{inputExist+1}; %they are trasnformed to lowercase if and only if legalvalues{l} was not empty
        end
    else %assign default value
        if char2logic && ischar(defparams{l}) && sum(strcmpi(tobeconverted,defparams{l})) %convert to logical if appropriate
            varargout{l} = char2logical(defparams{l},char2logical_true,char2logical_false);
        else        
            varargout{l} = defparams{l}; %they are returned as they are (no lowercase forcing)
        end
    end
end

return
end

function throw_error(msg,debug)
if ~debug
    stack = dbstack(2);
    stack(1).line = 0;
else
    stack = dbstack('-completenames');
end
errorStruct.message = msg;
errorStruct.stack = stack;
error(errorStruct);
return
end

function out = char2logical(inp,truecase,falsecase)
%convert char varargin to logical
switch inp
    case truecase
        out = 1;
    case falsecase
        out = 0;
end
out = logical(out);
return
end