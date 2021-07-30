function [res,model_info] = fmri_cleaning(data,polort,ort,varargin)
%FMRI_CLEANING mimics (most of) AFNI 3dTproject:
%  "This program projects (detrends) out various 'nuisance' time series from
%   each voxel in the input dataset.  Note that all the projections are 
%   done via linear regression, including the frequency filtering.  In this
%   way, you can bandpass time-censored data, and at the same time, remove 
%   other time series of no interest (e.g., physiological estimates, motion
%   parameters)." (from the help of 3dTproject)
%
%   Despite FMRI_CLEANING was designed specifically to clean fMRI data, it 
%   can be use to clean (i.e., detrend) any vector or matrix.
%
%Input Options:
%
% -DATA can be a matrix or the path to a nifti file, if DATA is a matrix, the 
%   last dimension must be time (e.g., [XxYxZxTIME] or [VOXELSxTIME]).
% -POLORT is an integer for removing polynomials up to and including degree
%   "polort". Polort up to order 5 are supported. Use POLORT = -1 to avoid
%   polynomial regression (NOT advised, unless data is mean centred).
%   For concatenated datasets, each run gets a separate set.
% -ORT are the confounding timeseries. ORT must be a matrix [TIMExVARIABLES].
%   Data will be orthogonalised with respect to ORT. For concatenated datasets
%   you can construct ORT so to have one regressor across all runs or 
%   separate each run regressors (e.g., you can concatenate realignment 
%   parameters (RP) across runs. You can use 6 RP regressors or have 
%   6*n_run regressors). 
%   Use ORT = [] to avoid regression of confounders. 
%   The mean will be removed from ort.
%
%Additional options can be specified using the following parameters (each 
%parameter must be followed by its value ie,'param1',value1,'param2',value2):
%
%  'passband'  : [TR,F1,F2] Set this option to regress out undesired 
%                frequencies using a basis of sines and cosines. TR is the 
%                repetition time, F1 and F2 are the frequency edges of the 
%                bandpass filter.
%                Use F1= 0, to perform a low-pass filter
%                Use F2=inf to perform a high-pass filter
%                {default: []}. 
%                For concatenated datasets, each run gets a separate set.
%  'cens'      : a binary vector (lenght of data) indicating volumes to be
%                considered. That is, zero-indexed volumes will be censored.
%                You can use fmri_censoring_mask to construct such vectors.
%                {default: []}. 
%  'censmode'   : ['zero' / 'kill'] Specifies how censored points are treated: 
%               'Zero', put zero values in their place {default: 'zero'}.
%               'kill', remove those time points NB: output dataset will be
%                       shorter than the input DATA.
%  'concat'    : An array of integer values for specifing the starting index
%                of each run (index starts from 1). E.g., [1 240 480].
%                This option should be always used when you want to
%                concatenate multiple runs. In this way POLORT and PASSBAND
%                regressors are constructed for each run separately.
%                Conversly, ORT run on the entire timeseries. If you want
%                them to be split accordingly to the run, you must do it
%                manually. {default = []}. 
%  'mask'      : [matrix or nifti path] Only non nonzero elements in mask
%                will be considered for cleaning.
%                {default = []}.
%  'writeNii'  : ['on'/'off'] Write the residuals in Nifti format. This 
%                option works only if the input data is provided as a nifti 
%                file. The file will be created in the working directory
%                with the suffix "_cleaned".
%  'restoremean':['on'/'off'] restore mean value after regression. 
%                {default = 'off'}
%  'removePol0': ['on'/'off'], in case you want POLOROT > 0 but don't want the 
%                constant term this option remove Pol 0 (ie., the costant)
%                NB: dangerous option. {default = 'off'} 
%
%Output variables
% -RES are the residuals of the linear model, i.e., the cleaned timeseries
% -MODEL_INFO is a structure containing information about the model as well
%   as the model itself (i.e., the explenatory variables)
%   You can plot the model:
%       figure; imagesc(model_info.X); xlabel('explenatory variables');
%       ylabel('time');
%
%Requirements:
% SPM (https://www.fil.ion.ucl.ac.uk/spm/) is required if DATA is passed as
% a Nifti file.

%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com


%--------------VARARGIN----------------------------------------------------
params   = {'censmode','concat','ortdemean','removePol0','restoremean','writeNii','passband' 'cens', 'mask'}; 
defparms = {'zero',         [],          1,       'off',        'off',     'off',        [],     [],     []};
legalvalues{1} = {'zero','kill'};
legalvalues{2} = {@(x) (isempty(x) || (~ischar(x) && sum(mod(x,1))==0 && sum((x < 0)) == 0)),'Only positive integers are allowed, which represent the starting indexes of the runs.'};
legalvalues{3} = [0 1]; %NB: In case of multiple sessions with one separate set of ort for each session, deaming all X or separately each ort session is the same. 
legalvalues{4} = {'on','off'};
legalvalues{5} = {'on','off'};
legalvalues{6} = {'on','off'};
legalvalues{7} = {@(x) (isempty(x) || (~ischar(x) && numel(x)==3) && x(2) <= x(3)),'Passband requires a 3-element vector: [TR,F1,F2], with F1 <= F2.'};
legalvalues{8} = {@(x) (isempty(x) || (~ischar(x) && length(unique(x)) <= 2) && max(x) <= 1),'Cens requires a binary vector, where 1-s indicate volumes to be preserved.'};
legalvalues{9} = [];
[cenmode,concat_index,ortdemean,removePol0,restoremean,writeNii,passband,cens,mask] = ParseVarargin(params,defparms,legalvalues,varargin,1);
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
%------LOADING DATA and reshape--------------------------------------------
if ischar(data)  %in case data is a path to a nifti file
    [~,name] = fileparts(data); name = remove_nii_ext(name);
    [~,hdr] = evalc('spm_vol(data);'); % to avoid an annoying messange in case of .gz
    data = spm_read_vols(hdr);
    s = size(data);
    n_dimension = length(s);
    if restoremean; MEAN=mean(data,n_dimension); end
    data = reshape(data,[prod(s(1:end-1)),s(end)])';
else
    % In this case the last dimension must be time!
    s = size(data);
    n_dimension = length(s);
    if n_dimension > 2
        if restoremean; MEAN=mean(data,n_dimension); end
        data = reshape(data,[prod(s(1:end-1)),s(end)])';
    elseif n_dimension == 2 && ~isvector(data)
        if restoremean; MEAN=mean(data,n_dimension); end
        data = data';
    elseif isvector(data) %only in this case, whatever orientation is fine
        if isrow(data)
            vector_flipped = 1;
            data = data';
        else
            vector_flipped = 0;
        end
        if restoremean; MEAN=mean(data); end
    end
    %override writeNii
    if writeNii; warning('You have to provide DATA as Nifti if you want the ouptut to be written to disk as Nifti.'); end
    writeNii = 0;
end
colons = repmat({':'},1,n_dimension-1);   %might be used for selecting non temporal dimensions
%NB: data must be provided with the last dimension as time, but at the end
%data will be reshaped with the first dimension as time (this allows easy
%handling of preloaded data)
%--------------------------------------------------------------------------

if polort == -1 || removePol0
    warning('The intercept is not modeled. Be sure data have been demeaned.');
end

%------------------- Masking if required ----------------------------------
if ischar(mask)  %in case the mask is a path to a nifti file
    [~,maskhdr] = evalc('spm_vol(mask);'); % to avoid an annoying messange in case of .gz
    mask = spm_read_vols(maskhdr);
    smask = size(mask);
    if ~logical(smask == s(1:3))
        error('The input mask does not have the same dimension of data');
    end
    n_dimension = length(s);
    if restoremean %we need to mask also the mean; if it's not binary, we 
        %get an error later on
        MEAN(mask == 0) = 0;
    end
    mask = reshape(mask,[1,smask(1)*smask(2)*smask(3)]);
elseif ~isempty(mask)
    smask = size(mask);
    if  length(smask) > 2  %it's a 3d volume
        if ~isequal(s(1:end-1),sr) 
            error('The input mask does not have the same dimension of data');
        end
        if restoremean %we need to mask also the mean; if it's not binary, we 
        %get an error later on
            MEAN(mask == 0) = 0;
        end
        % ok, reshape
        ROI = reshape(mask,[1,sr(1)*sr(2)*sr(3)]);
    elseif isvector(mask)  %in this case the data is assumed to be 2D
        if numel(mask)~= s(1)
            error('The input mask does not have the same dimension of data');
        end
        if restoremean %we need to mask also the mean; if it's not binary, we 
            %get an error later on
            MEAN(mask == 0) = 0;
        end
        % no need to reshape
    else
        error('The input mask does not have the same dimension of data');
    end
end
% extract indexes, and mask data
if ~isempty(mask)
%     un = unique(mask(:));
%     if length(un) > 2 || sum(uint8(un)) > 1
%         %the mask is not binary, make it binary
%         mask(mask>0) = 1;
%     end
    mask_index = find(mask);
    data = data(:,mask_index);
else
    mask_index = [];
end
%--------------------------------------------------------------------------

% %------------------- Automask if required ---------------------------------
% automask = 0; % In most of the cases, this option is useless
% if automask
%     stdv = std(data);
%     badIndex = zeros(1,length(stdv));
%     badIndex(stdv == 0) = 1;
%     badIndex(isnan(stdv)) = 1;
%     data(:,logical(badIndex)) = [];
% end
% %--------------------------------------------------------------------------

%---- find out how many runs are present-----------------------------------
if ~isempty(concat_index)
    run_number = length(concat_index);
    if concat_index(1) ~= 1
        error('concat_index(1) must be 1.')
    end
    if run_number < 2 % no concat case
        N = size(data,1);
        run_number = 1;
    else
        N = zeros(run_number,1);
        for l = 2:run_number
            N(l-1) = concat_index(l) - concat_index(l-1);  
        end
        N(end) = size(data,1) - concat_index(end) + 1;
    end
else
    N = size(data,1);
    run_number = 1;
end
%--------------------------------------------------------------------------

% ---Construct POLORT and PASSBAND regressors------------------------------
X_dim = zeros(run_number,1);
for r = 1:run_number
    %add Legendre pol
    X_diag{r} = LegPol(N(r),polort,removePol0);
    if ~isempty(passband)
        % sines and cosines to apply a band_pass filter
        X_diag{r} = [X_diag{r},SineCosineBasis(N(r),passband(1),passband(2),passband(3),1)];
    end
    X_dim(r) = size(X_diag{r},2);
end

X = [];
for l_row = 1:run_number
    X_row = [];
    for l_col = 1:run_number
        if l_row == l_col
            X_row = [X_row,X_diag{l_row}];
        else
            X_row = [X_row,zeros(N(l_row),X_dim(l_col))];
        end
    end
    X = [X;X_row];
end
%--------------------------------------------------------------------------

%---------------add ort regressors-----------------------------------------
if ~isempty(ort)
    %demean ort
    if ortdemean
        ort = ort -mean(ort,1);
    end
    M_ort = size(ort,2);
    try
        X = [X,ort];
    catch
        error('ort regressors don''t have the correct size');
    end
else
    M_ort = 0;
end
%--------------------------------------------------------------------------

%get now some basic info
N_tot = size(data,1);
M = size(X,2);
M_band = size(X,2) - M_ort -(polort+1)*run_number;

%censoring initialization
if ~isempty(cens)
    if length(cens) ~= N_tot
        error('Cens length doesn''t match the time points in data.');
    end
    N_cens = sum(cens == 0);
else 
    N_cens = 0;
end

%calculate dof
dof = N_tot - M - N_cens;

% Print Model info
fprintf('-> fMRI denosing. Model info:');
fprintf('\nNumber of time points \t\t%d',N_tot);
fprintf('\nNumber of total regressors \t%d',M);
fprintf('\nNumber of Legendre pol \t\t%d',(polort+1)*run_number);
fprintf('\nNumber of passband regressors \t%d',M_band);
fprintf('\nNumber of ort regressors \t%d',M_ort);
fprintf('\nNumber of censored points \t%d',N_cens);
fprintf('\nNumber of Residual tDoF \t%d\n',dof);

if dof <= 0
    error('Not enough tDoF to proceed.');
elseif dof < 10
    warning('Low number of tDoF. Results might be not reliable.');
end

%do censoring
if ~isempty(cens)
    if strcmpi(cenmode,'kill')
        data(cens == 0,:) = [];
        X(cens == 0,:) = [];
    elseif strcmpi(cenmode,'zero')
        data(cens == 0,:) = 0;
        X(cens == 0,:) = 0;
    end
end

%let's do the glm
res = data - X*(X\data);    
Nfinal = size(res,1);  %temporal points, they may have changed depending on censoring

% if masking was applied we have to restore the original shape
if ~isempty(mask)
    RES = zeros(Nfinal,length(mask));
    RES(:,mask_index) = res;
    res = RES; 
end

% % if automask was applied we have to restore the original shape
% if automask
%     RES = zeros(Nfinal,length(badIndex));
%     RES(:,not(logical(badIndex))) = res;
%     res = RES; clear RES;
% end

%reshape res to original data size
if n_dimension > 2 && ~isvector(data)
    res = reshape(res',[s(1:end-1),Nfinal]);
else % last dim must return to be time
    res = res';
end

%restore mean if requested
if restoremean 
    if isempty(cens) || strcmpi(cenmode,'kill')
        res = res + MEAN;
    else % zero volumes must be preserved 
        MEAN = repmat(MEAN,[ones(1,n_dimension-1),Nfinal]);        
        MEAN(colons{:},~cens) = 0;
        res = res + MEAN;
    end
end

if writeNii
    output_name = [name,'_cleaned.nii'];
    for l = 1:Nfinal
        hdr(l).fname = output_name;
        hdr(l).private.dat.fname = output_name;
        if ~isempty(cens) && strcmpi(cenmode,'kill')
            hdr(l).private.dat.dim(end) =  Nfinal;
        end
        spm_write_vol(hdr(l),res(colons{:},l));
    end
end

% construct model_info variable
model_info.N = N_tot;
model_info.Mtot = M;
model_info.Mpol = (polort+1)*run_number;
model_info.Mband = M_band;
model_info.Mort = M_ort;
model_info.Mcens = N_cens;
model_info.dof = dof;
model_info.X = X;
if ~isempty(cens)
    model_info.cenmode = cenmode;
else
    model_info.cenmode = [];
end

return
end

function s = remove_nii_ext(s)
% in case you pass a .gz, fileparts remove the last ext and not the .nii
indx =  strfind(s,'.nii');
s(indx:end) = [];
return
end
