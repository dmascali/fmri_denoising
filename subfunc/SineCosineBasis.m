function [X] = SineCosineBasis(N,TR,F1,F2,invert,varargin)
% Similar to AFNI function 1dBport creates a set of columns of sines and
% cosines for the purpose of bandpassing via regression. 
%  -N  = number of time points
%  -TR = repetition time
%  -F1 = lower  frequency edge 
%  -F2 = higher frequency edge (use inf for the nyquist frequency)
%  Waves (f) that satisfy the equations F1 <= f <= F2 are created.
%  -invert [0/1]. Use invert = 1 to create waves outside the band of
%                 interest (F1>f & f>F2). This is the option you should use 
%                 if you want bandpassing data via linear regression. 
%                 Otherwise (invert = 0) the desired frequencyies
%                 are created (F1<=f<=F2). 
%
%Additional options can be specified using the following parameters (each 
%parameter must be followed by its value ie,'param1',value1,'param2',value2):
%
%  'concat'    : An array of integer values for specifing the starting index
%                of each run (index starts from 1). E.g., [1 240 480].
%                This option should be always used when you want to
%                concatenate multiple runs. In this way the basis functions
%                are constructed separately for each run.{default = []}.
%  'polort'    : is an integer for including polynomials up to and including
%                degree "polort". {default = -1; i.e., no polort}
%
% NB: By default, the constant component is not produced.
%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com

%--------------VARARGIN----------------------------------------------------
params   = {'concat','polort'}; 
defparms = {      [],      -1};
legalvalues{1} = [];
legalvalues{2} = [-1 0 1 2];
[concat_index,polort] = ParseVarargin(params,defparms,legalvalues,varargin,1);
% -------------------------------------------------------------------------

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
    %add Legendre pol if requested
    X_diag{r} = LegPol(n(r),polort);
    %create sines and cosines
    X_diag{r} = [X_diag{r}, basis_constructor(n(r),TR,F1,F2,invert)];
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

function Y = basis_constructor(N,TR,F1,F2,invert)

%define frequency grid
deltaf = 1/(N*TR);
nyquist = 1/(2*TR);
freq = deltaf:deltaf:nyquist;

index_selector = ( F1 <= freq & freq <= F2);

if invert
    index_selector = not(index_selector);
end

Nreg = sum(index_selector);

if Nreg == 0
    Y = [];
    disp('No frequency found.');
    return
end

index = find(index_selector);

x = 0:1:(N-1);
t = x.*TR;

Y = [];
for l = 1:Nreg
    Y(:,end+1) = cos(2*pi*freq(index(l))*t)';
    Y(:,end+1) = sin(2*pi*freq(index(l))*t)';
end

% Replace the nyquist frequency if present
if freq(index(end)) == nyquist
    Y(:,end-1:end) = [];
    %in this case N is even
    tmp = ones(N,1);
    tmp(1:2:N,1) = -1;
    Y = [Y,tmp];
end

%figure; plot(Y);


return
end
