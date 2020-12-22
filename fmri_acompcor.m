function X = fmri_acompcor(data,rois,dime,varargin)
%FMRI_ACOMPCOR(DATA,ROIS,DIME) extracts signals from DATA using ROIS as masks.
%Inputs:
% -DATA can be a matrix or the path to a nifti file. If DATA is a matrix, the 
%   last dimension must be time. (e.g., [XxYxZxTIME] or [VOXELSxTIME])
% -ROIS is a cell array containg either matrices or paths to nifti files.
%   (e.g., {ROI1, ROI2, ROI3} ). ROIS must be binary.
% -DIME is a vector that specifies for each ROI the number/type of signals: 
%  - DIME = 0            Only the mean signal will be extracted.
%  - DIME > 0 & integer  The first n=DIME principal components will be 
%                        extracted (following the aCompCor approach). 
%                        E.g., DIME = [0 5 5], extracts the mean from the 
%                        first ROI and the first 5 principal components from
%                        each of the remaining ROIs.
%  - 0<DIME<1;           In this case DIME specifies the desired percentage
%                        of variance that the components must explain. The 
%                        first n principal components that satisfy the
%                        variance requirement are extracted. 
%                        E.g., [0.5 0.5], extracts enough principal components
%                        to explain at least 50% of the variance in each ROI.  
%                        This method was proposed by Muschelli et al. (2014)
%
%Additional options can be specified using the following parameters (each 
% parameter must be followed by its value ie,'param1',value1,'param2',value2):
%
%   'confounds'   : If a matrix of confound variables (time x variables) is 
%                   provided, DATA will be orthogonalized with respect to 
%                   the confounds before extracting any signal. In this
%                   way, an fMRI denoising model containing both the
%                   confounds and the aCompCor signals will be more
%                   predictive. {default = []}
%   'filter'      : As above, but regress out undesired frequencies using
%                   a basis of sines and cosines. Use a 3 element vector,
%                   i.e., ...'filter',[TR,F1,F2]. Where TR is the
%                   repetition time, F1 and F2 are the frequency edges of
%                   the bandpass filter. {default= []}
%   'PolOrder'    : Legendre Polynomial Order, as bove, for regressing
%                   constant, linear or quadratic terms. PolOrder can be: 
%                   -1 : skip (not possilbe when DIME > 0)
%                    0 : constant term 
%                    1 : constant + linear terms {default}
%                    2 : constant + linear + quadratic terms
%   'Concat'      : An array of integer values for specifing the starting index
%                   of each run in case of concatenated runs (index starts 
%                   from 1). E.g., [1 240 480]. {default = []}.
%   'derivatives' : is a vector of roi length that specifies the degree of
%                   derivatives to be computed on the extracted signals
%                   {default=[],which is the same as zeros(1,length(rois))
%   'squares'     : is a vector of roi length that specifies (flag 1/0)
%                   whether to compute the squares of the extracted signals
%                   (as well as the squares of derivatives, if present)
%                   {default=[],which is the same as zeros(1,length(rois))
%
% If DIME = 0, you can either extract the mean (def) or the median,
% by using the following parameter:
%
%   'type'        : ['mean','median'] By defaul the function extracts the
%                   mean. This parameter works only with DIME = 0. 
%
% Advanced/experimantal options (usually you don't need to change the 
% default values) are:
%
%   'firstmean'   : ['on'/'off'] If 'on', the first extracted component is 
%                   the mean signal, then PCA is performed on data 
%                   ortogonalised with respect to the mean signal.
%                   {default='off'}
%   'FullOrt'     : ['on'/'off'] If 'on', for each ROI (excluding the first
%                   one) data is ortogonalised with respect to the
%                   signals/components extracted up to that point (including
%                   derivatives or square terms if present). In this way,
%                   the full set of extracted signals are orthogonal. 
%                   {default='off'}
%   'DatNormalise': ['on'/'off'], if set to 'on' the data is normalised by
%                   its temporal variance before performing PCA {default='off'}
%   'SigNormalise': ['on'/'off'], if set to 'on' the extracted signals are
%                   normalised to unit variance {defaul='on'}
%   'DataCleared' : ['true'/'false'] if 'true' the function avoids to seek
%                   for NaNs or voxels with zero signals (to save time)
%                   {default='false'}
%
% NB1: By default, data is detrended (costant and linear trends are removed)
%      before any computation (unless 'PolOrder' is specified)
% NB2: All extracted signals are normalised to unit variance (unless 'SigNormalise'
%      is set to 'off')
% NB3: When DIME = 0, ortogonalizing the data with respect to 'confounds', 
%     'filter' or 'PolOrder' has no effect on the final denoising model.
%
%Requirements:
% SPM (https://www.fil.ion.ucl.ac.uk/spm/) is required if DATA and ROIs are 
% passed as Nifti files.
%
%References:
% - Behzadi et al. (2007) Neuroimage 37, 90-101   
% - Whitfield-Gabrieli and Nieto-Castanon (2012) Brain Connect. 2(3), 154-41
% - Muschelli et al. (2014) Neuroimage 96, 22-35
%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com


if nargin < 3
    error('Not enough input.');
end

%--------------VARARGIN----------------------
params  =  {'confounds','firstmean','derivatives','squares','DatNormalise','DataCleared','filter','PolOrder','FullOrt', 'SigNormalise', 'concat', 'type'};
defParms = {         [],      'off',           [],       [],          'off',      'false',      [],        1     'off',           'on',       [], 'mean'};
legalValues{1} = [];
legalValues{2} = {'on','off'};
legalValues{3} = [];
legalValues{4} = [];
legalValues{5} = {'on','off'};
legalValues{6} = {'true','false'};
legalValues{7} = [];
legalValues{8} = [-1 0 1 2];
legalValues{9} = {'on','off'};
legalValues{10} ={'on','off'};
legalValues{11} = [];
legalValues{12} = {'mean','median'};
[confounds,firstmean,deri,squares,DatNormalise,DataCleared,freq,PolOrder,FullOrt,SigNormalise,ConCat,MetricType] = ParseVarargin(params,defParms,legalValues,varargin,1);
% --------------------------------------------

if ~iscell(rois)
    error('Please provide rois as cell, i.e., rois = {''path1'',''path2.nii''} or rois = {matrix1,matrix2}');
end
n_rois = length(rois);
if length(dime)~=n_rois
    error('Please specify a dime for each rois e.g., dime = [5 5]');
end
if ~isempty(deri)
    %check if there is one value for each roi
    if length(deri) ~= n_rois
        error('Please specify a derivative value for each roi e.g., ...,''derivatives'',[1 0]');
    end
else
    deri = zeros(1,n_rois);
end

if ~isempty(squares)
    %check if there is one value for each roi
    if length(squares) ~= n_rois
        error('Please specify a square flag for each roi e.g., ...,''squares'',[1 0]');
    end
else
    squares= zeros(1,n_rois);
end

%--------------------------------------------------------------------------
%------LOADING DATA and reshape--------------------------------------------
if ischar(data)  %in case data is a path to a nifti file
    data = spm_read_vols(spm_vol(data));
    s = size(data);
    data = reshape(data,[s(1)*s(2)*s(3),s(4)])';
else
    % In this case the last dimension must be time!
    s = size(data);
    n_dimension = length(s);
    if n_dimension > 2
        data = reshape(data,[prod(s(1:end-1)),s(end)])';
    else
        data = data';
    end
end
%NB: data must be provided with the last dimension as time, but at the end
%data will be reshaped with the first dimension as time (this allows easy
%handling of preloaded data)
%--------------------------------------------------------------------------

X = []; %output variable

for r = 1:n_rois
    %----------------------------------------------------------------------
    %------LOADING ROI, Checking compatibility with data and reshape-------
    if ischar(rois{r})
        ROI = spm_read_vols(spm_vol(rois{r}));
        sr = size(ROI);
        %check if s and sr are identical in the first 3 dimensions
        if ~logical(sr == s(1:3))
            error(sprintf('ROI %d does not have the same dimension of data',r));
        end
        ROI = reshape(ROI,[1,sr(1)*sr(2)*sr(3)]);
    else
        ROI = rois{r};
        sr = size(ROI);
        if length(sr) > 2  %it's a 3d volume
            if ~isequal(s(1:end-1),sr) 
                error(sprintf('ROI %d does not have the same dimension of data',r));
            end
            % ok, reshape
            ROI = reshape(ROI,[1,sr(1)*sr(2)*sr(3)]);
        elseif isvector(ROI)  %in this case the data is assumed to be 2D
            if numel(ROI)~= s(1)
                error(sprintf('ROI %d does not have the same dimension of data',r));
            end
            % no need to reshape
        else
            error(sprintf('ROI %d does not have the same dimension of data',r));
        end
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %check if ROI is binary
    un = unique(ROI(:));
    if length(un) > 2 || sum(uint8(un)) > 1
        error(sprintf('ROI %d is not binary',r));
    end
    ROI = uint8(ROI);
    %----------------------------------------------------------------------
    % check if PolOrder is compatible with dime
    if dime(r) > 0 && PolOrder == -1
        warning('PCA should run on data demeaned. Changing PolOrder from -1 to 1');
        PolOrder = 1;
    end
    %----------------------------------------------------------------------
    % data extraction
    indx = find(ROI);
    V = data(:,indx);
    %----------------------------------------------------------------------
    % remove voxels whose variance is equal zero (no signal in those voxels)
    % and Nan values
    if ~DataCleared
        stdv = std(V);
        a = zeros(1,length(stdv));
        a(stdv == 0) = 1;
        a(isnan(stdv)) = 1;
        V(:,logical(a)) = [];
    end
    %------------------Orthogonalise V-------------------------------------
    COV = [];
    if firstmean && dime(r) > 0 % as done in CONN: first extract the mean signal (mS), then compute PCA over data ortogonalised with respect to mS. 
        % to get a "clean" mean signal, I have to remove trends from V.
        % Here PolOrder can't be -1
        Xtrends = LegPol(size(V,1),PolOrder,0,'concat',ConCat);
        V = V -Xtrends*(Xtrends\V); 
        mS = mean(V,2);
        % add the mean to COV
        COV = [COV,mS];
    end
    if PolOrder ~= -1  %regress trends
        COV = [COV,LegPol(size(V,1),PolOrder,0,'concat',ConCat)];
    end
    if ~isempty(freq) 
        COV = [COV,SineCosineBasis(size(V,1),freq(1),freq(2),freq(3),1,'concat',ConCat)];
    end
    if ~isempty(confounds)  
        COV = [COV,confounds];
    end
    if FullOrt  %include also already extracted signals (if present)
        COV = [COV,X];
    end
    if ~isempty(COV)
        V = V-COV*(COV\V);
    end
    %----------------------------------------------------------------------
    if dime(r) > 0
        % force the mean to be zero (the distribution of mean values may be
        % slightly shifted if cofounds have been regressed). Also, again
        % remove trends to avoid SVD failure
        Xtrends = LegPol(size(V,1),PolOrder,0,'concat',ConCat); V = V -Xtrends*(Xtrends\V); 
        if DatNormalise
            %tvariance normalization
            V = bsxfun(@rdivide,V,std(V));
        end
        try % PCA may still fail to converge. In this case normalising the data usually solves the problem 
            [U,P] = svd(V);
        catch ME 
            if DatNormalise
                % nothing to do!
                rethrow(ME)
            else % let's try with DataNormalise
                warning('PCA failed to converge. Second attempt using temporal variance normalisation...');
                V = bsxfun(@rdivide,V,std(V));
                [U,P] = svd(V);
            end
        end
        % How many components to extract?
        if mod(dime(r),1) == 0  % dime is an integer, and specifies the number of components
            D = dime(r);
            if firstmean; D = D-1; end %remove one dimension
        else  % dime specifies the percentage of variance to extract
            latent = diag(P.^2./(size(V,1)-1));
            latent = latent./sum(latent); %normalise
            indexes = find(cumsum(latent)>dime(r));
            D = indexes(1);
        end
        % Select the components
        comp = U(:,1:D)*diag(diag(P(1:D,1:D)));      
        if firstmean %add the mean signal as first component 
            comp = [mS,comp];
        end
    else  %if dime == 0 simply compute the straight average or median
        switch MetricType
            case {'mean'}
                comp = mean(V,2);
            case {'median'}
                comp = median(V,2);
        end
    end
    %----------------------------------------------------------------------
    % derivatives computation, if requested
    if deri(r) > 0
        d = [];
        for l = 1:deri(r)
            d = [d, [zeros(l,size(comp,2));diff(comp,l)] ]; % I have to add l zeros as first rows...
        end
    else 
        d = [];
    end
    Xtmp = [comp,d];
    %----------------------------------------------------------------------
    % squares computation, if requested
    if squares(r) > 0
        Xtmp = [Xtmp,Xtmp.^2];
    end
   
    X = [X,Xtmp];
end

if SigNormalise
    %variance normalise the extracted components
    X = X./std(X,0,1);
end

return
end