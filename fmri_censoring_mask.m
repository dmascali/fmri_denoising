function [tmask,n_cens] = fmri_censoring_mask(vect,thr,varargin)
%Given a vector VECT the function creates a binary temporal mask which is 
%   0 where VECT  > THR  (to be censored)
%   1 where VECT <= THR  (to be kept)
% this is the same convention of ANFI 3dTproject (censored points are 
% indicated with 0s)
% If an array of THR is provided, multiple masks will be created.
%
%Additional options can be specified using the following parameters (each 
% parameter must be followed by its value ie,'param1',value1,'param2',value2)
%
%Besides the standard mode (based on absolute THR values) you can create
% masks based on different algorithm:
%   - 'mode' = ['standard'/'random'/'top'/'bottom'] {default = 'standard'}
%       -standard- works as described above.
%       -random  - create random masks. In this modality, the THR values
%                  specify the number of volumes to be removed (or the
%                  percentage of volumes, depending on the'Units' property).
%       -top     - remove the top N (or percent, depending on the'Units' 
%                  property) points based on the VECT values. 
%       -bottom  - remove the bottom N (or percent, depending on the'Units' 
%                  property) points based on the VECT values. 
%   - 'Units' = ['N'/'percent']  {default = 'N'}
%      specify the units (number of points or percentage) for modality 
%      'random','top' and 'bottom'  
%
% In the standard mode only, you can specify the following parameters:
%   - 'preTR' = [M (integer value)],
%      for censoring M previous TR volumes {default = []}
%   - 'postTR' = [M (integer value)],
%      for censoring M post TR volumes {default = []}
%
% Use 'verbose' property (1/0) to show output messages {defalut = 1}
%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com


%--------------VARARGIN----------------------
params  =  {'preTR','postTR',    'mode',  'units','verbose'};
defParms = {     [],      [],'standard',      'N',        1};
legalValues{1} = [];
legalValues{2} = [];
legalValues{3} = {'standard','random','top','bottom'};
legalValues{4} = {'percent','N'};
legalValues{5} = [0 1];
[pre_TR,post_TR,Mode,Units,verbose] = ParseVarargin(params,defParms,legalValues,varargin);
% %------------------------------------------

%vect must be a row vector
if not(isrow(vect))
    vect = vect';
end

N_masks = length(thr);
if strcmpi(Mode,'standard'); Units = 'vector unit';end
if verbose; fprintf('\nCreating %d temporal mask(s) in modality ''%s'' (units: ''%s'')...',N_masks,Mode,upper(Units)); end 
N = length(vect);
tmask = zeros(N_masks,N);

switch Mode
    case {'standard'}
        for l = 1:N_masks
            tmask(l,:) = cens_standard(vect,thr(l),pre_TR,post_TR,N);
        end
    case {'random','top','bottom'}
        %------------------------------------------------------------------
        switch Units %this section is common to all
            case {'n'}
                if ~isempty(find(thr> N,1))
                    error('In the current modality THR represents the number of volumes to be randomly censored, thus, thr must be <= N.');
                end                 
            case {'percent'}
                if ~isempty(find(thr> 100 | thr < 0,1))
                    error('In the current modality THR represents the percentage of volumes to be randomly censored, thus, thr must be <= 100.');
                end 
                % convert thr from percent to N
                thr = round(thr.*N/100);
        end
        %------------------------------------------------------------------
        switch Mode
            case {'random'}
                for l = 1:N_masks
                    indx = randsample(N,thr(l),'false'); %random sample without replacements
                    tmask(l,indx) = 1;
                end
            case {'top','bottom'}
                if strcmpi(Mode,'top')
                    sort_str = 'descend';
                else
                    sort_str = 'ascend';
                end
                for l = 1:N_masks
                    [~,indx] = sort(vect,2,sort_str);
                    indx = indx(1:thr(l));
                    tmask(l,indx) = 1;
                end
        end
        %------------------------------------------------------------------
end

n_cens = sum(tmask,2)';

% from logic to double plus invert 0s with 1s (to follow ANFI's convention)
tmask = uint8(~tmask);

if verbose; fprintf('done!\n'); end

return
end

function tmask = cens_standard(vect,thr,pre_TR,post_TR,N)

tmask = (vect > thr);
indx = find(tmask);

if ~isempty(pre_TR) && pre_TR > 0
   for l = 1:pre_TR
       indxtemp = indx; indxtemp(indxtemp <= l) = [];
       tmask(indxtemp-l) = 1;
   end
end

if ~isempty(post_TR) && post_TR > 0
   for l = 1:post_TR
       indxtemp = indx; indxtemp(indxtemp >= (N-l+1)) = [];
       tmask(indxtemp+l) = 1;
   end
end

return
end

