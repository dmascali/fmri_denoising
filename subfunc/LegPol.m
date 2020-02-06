function X = LegPol(N,order,exclude_constant,varargin)
%return leg pol for regression up to order 2
% -N     = number of time points
% -order =  desired polynomial order (up to 2), i.e., :
%     -1 : empty matrix
%      0 : constant term
%      1 : constant + linear terms
%      2 : constant + linear + quadratic terms
%
% Use the exclude_costant flag (optional) to exlude the constant term (NB
% you normally don't want to do that!)
%
%Additional options can be specified using the following parameters (each 
%parameter must be followed by its value ie,'param1',value1,'param2',value2):
%
%  'concat'    : An array of integer values for specifing the starting index
%                of each run (index starts from 1). E.g., [1 240 480].
%                This option should be always used when you want to
%                concatenate multiple runs. In this way the basis functions
%                are constructed separately for each run.{default = []}.
%
%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com

if order == -1
    X = [];
    return
end

%--------------VARARGIN----------------------------------------------------
params   = {'concat'}; 
defparms = {      []};
legalvalues{1} = [];
[concat_index] = ParseVarargin(params,defparms,legalvalues,varargin,1);
% -------------------------------------------------------------------------

if nargin < 3
    exclude_constant = 0;
end

%---- find out how many runs are present-----------------------------------
if ~isempty(concat_index)
    run_number = length(concat_index);
    if concat_index(1) ~= 1
        error('concat_index(1) must be 1.')
    end
    if run_number < 2 % no concat case
        n = N;
        run_number = 1;
    else
        n = zeros(run_number,1);
        for l = 2:run_number
            n(l-1) = concat_index(l) - concat_index(l-1);  
        end
        n(end) = N - concat_index(end) + 1;
    end
else
    n = N;
    run_number = 1;
end
%--------------------------------------------------------------------------

X_dim = zeros(run_number,1);
X_diag{1} = [];  

for r = 1:run_number

    if order == -1
        X_diag{r} = [];
        return
    elseif order == 0
        if ~exclude_constant
            X_diag{r} = ones(n(r),1);
        else
            X_diag{r} = [];
        end
    elseif order == 1
        if ~exclude_constant
            X_diag{r} = [ones(n(r),1), linspace(-1,1,n(r))'];    
        else
            X_diag{r} = linspace(-1,1,n(r))';    
        end
    elseif order == 2
        if ~exclude_constant
            % the 2nd order is: P2(x) = 0.5*(3*x.^2 -1)
            t = linspace(-1,1,n(r))';
            X_diag{r} = [ones(n(r),1), t, (0.5*(3*t.^2-1))];    
        else
            X_diag{r} = [linspace(-1,1,n(r))',(linspace(-1,1,n(r)).^2)'];   
        end  
    end
    
    %get size
    X_dim(r) = size(X_diag{r},2);
    
end

% assemble pieces  
X = [];
for l_row = 1:run_number
    X_row = [];
    for l_col = 1:run_number
        if l_row == l_col
            X_row = [X_row,X_diag{l_row}];
        else
            X_row = [X_row,zeros(n(l_row),X_dim(l_col))];
        end
    end
    X = [X;X_row];
end

return
end