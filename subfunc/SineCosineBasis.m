function [Y] = SineCosineBasis(N,TR,F1,F2,invert)
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
% NB: the constant component is not produced.
%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com

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
    Y = [Y,cos(2*pi*freq(index(l))*t)'];
    Y = [Y,sin(2*pi*freq(index(l))*t)'];
end

% The sine of the nyquist frequency must be removed 
if freq(index(end)) == nyquist
    Y = Y(:,1:end-1);
end

%figure; plot(Y);
return
end
