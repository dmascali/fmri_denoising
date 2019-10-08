function X = LegPol(N,order,exclude_constant)
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
%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com

if nargin < 3
    exclude_constant = 0;
end

if order == -1
    X = [];
    return
elseif order == 0
    if ~exclude_constant
        X = ones(N,1);
    else
        X = [];
    end
elseif order == 1
    if ~exclude_constant
        X = [ones(N,1), linspace(-1,1,N)'];    
    else
        X = linspace(-1,1,N)';    
    end
elseif order == 2
    if ~exclude_constant
        X = [ones(N,1), linspace(-1,1,N)',(linspace(-1,1,N).^2)'];    
    else
        X = [linspace(-1,1,N)',(linspace(-1,1,N).^2)'];   
    end  
end

return
end