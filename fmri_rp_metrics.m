function RPmetrics = fmri_rp_metrics(rp,type,ReferenceImage)
%FMRI_RP_METRICS computes various Realignment Parameters (RP) based metrics
%
%Inputs:
% -RP can be the path to an RP file or the pre-loaded RP matrix. 
% -TYPE is a string specifying with which software was used to compute RP
%       {'fsl','spm','afni','hcp'}
% -REFERENCEIMAGE (is optional and is used for FD Jenkinson computation.  
%       Provide the path of the nifti volume used as base for realignment.
%
%Ouptut:
% -RPmetrics is a structure containing:
%           .FDpower (L1 norm)
%           .tra_abs (absolute traslations)
%           .rot_abs (absolute rotations)
%           .FDafni  (L2 norm)
%           .FDjenk  (rms relative to previous volume)
%           .rms     (rms relative to the reference volume)
% 
%Requirements:
% SPM (https://www.fil.ion.ucl.ac.uk/spm/) is required if FDjenk is
% desired.
%
%Acknowledgements:
% For FDjenk computation, this function adopts a code written by:
% YAN Chao-Gan 120930 (ycg.yan@gmail.com) that is part of the toolbox:
% DPARSF_V2.2_130309
%
%References:
% - Power, Jonathan D., et al. "Spurious but systematic correlations in 
%   functional connectivity MRI networks arise from subject motion." 
%   Neuroimage 59.3 (2012): 2142-2154.
% - Yan, Chao-Gan, et al. "A comprehensive assessment of regional variation 
%   in the impact of head micromovements on functional connectomics." 
%   Neuroimage 76 (2013): 183-201.
%__________________________________________________________________________
% Daniele Mascali
% Enrico Fermi Center, MARBILab, Rome
% danielemascali@gmail.com

switch lower(type)
    case {'fsl'} %ok
        rp_indx.rot = [1 2 3];
        rp_indx.tra = [4 5 6];  
        rp_rot_unit = 'rad';
    case {'spm'} %should be ok. Indeed they are plotted by SPM in degrees but they are saved in radians
        rp_indx.rot = [4 5 6];
        rp_indx.tra = [1 2 3];         
        rp_rot_unit = 'rad';
    case {'afni'} %ok
        rp_indx.rot = [1 2 3];
        rp_indx.tra = [4 5 6];         
        rp_rot_unit = 'deg';
    case {'hcp'} %verified
        rp_indx.rot = [4 5 6];
        rp_indx.tra = [1 2 3];  
        rp_rot_unit = 'deg';
end
rp_average_brain_radius = 50; %mm; 50 according to Power,

if ischar(rp)
    rp = load(rp);
end

%FD power  %L1 norm
[RPmetrics.FDpower,RPmetrics.tra_abs,RPmetrics.rot_abs] = FD_Power(rp,rp_indx,rp_rot_unit,rp_average_brain_radius);

%FD afni   %L2 norm + different rotational weight
RPmetrics.FDafni = FD_afni(rp,rp_indx,rp_rot_unit);

%FD Jenkinson
if nargin > 2
    [rel_rms, abs_rms] = FD_Jenkinson(rp,ReferenceImage,rp_indx,rp_rot_unit);
    RPmetrics.FDjenk = rel_rms; %rms relative to the previous volume
    RPmetrics.rms = abs_rms; %RMS relative to the reference volume
end

return
end

function [fd,tra_abs,rot_abs] = FD_Power(rp,rp_indx,rp_rot_unit,rp_average_brain_radius)

tra = rp(:,rp_indx.tra);
rot = rp(:,rp_indx.rot);

fd_tra = sum( abs( diff(tra)  ) ,2);
%rotaions must be converted from rad to mm. (assuming a radius of 50 mm)
switch lower(rp_rot_unit)
    case {'rad'}
        rot_mm = rot*rp_average_brain_radius; 
    case {'mm'} 
        %do nothing
        rot_mm = rot;
    case {'deg'}
        rot_mm = rot*rp_average_brain_radius*2*pi/360;         
end
fd_rot = sum( abs( diff(rot_mm)  ) ,2);
fd = sum([fd_tra,fd_rot],2); 
fd = [0; fd]; % due to the differenziation I need to add a 0.

%extract also abs rot and tra
% let's subtract the value at time 0, so in case the ref volume is not the
% first one the position still refers to the first scan. 
rp_mm = [tra,rot_mm];
rp_mm_dm = bsxfun(@minus, rp_mm, rp_mm(1,:));
tra_abs = sum(abs(rp_mm_dm(:,1:3)),2);
rot_abs = sum(abs(rp_mm_dm(:,4:6)),2);


return
end

function fd = FD_afni(rp,rp_indx,rp_rot_unit)
% AFNI use L2 instead of L1 (as it used in FDpower), moreover the rotations
% are scaled in a different way (rot in deg equals mm)

tra = rp(:,rp_indx.tra);
rot = rp(:,rp_indx.rot);

%rotaions must be converted to deg
switch lower(rp_rot_unit)
    case {'rad'}
        rot = 180*rot./pi; 
    case {'deg'}
        %do nothing       
end
fd = sqrt(sum(diff([tra,rot]).^2,2));
fd = [0; fd]; % due to the differenziation I need to add a 0.


return
end

function [rel_rms, abs_rms] = FD_Jenkinson(RP,ReferenceImage,rp_indx,rp_rot_unit)
% function [rel_rms, abs_rms] = y_FD_Jenkinson(RealignmentParameterFile,ReferenceImage)
% Calculate FD Jenkinson (relative RMS) and absolute RMS based on SPM's realignment parameters
% Reference: Jenkinson, M., Bannister, P., Brady, M., Smith, S., 2002. Improved optimization for the robust and accurate linear registration and motion correction of brain images. Neuroimage 17, 825-841.
%            Jenkinson, M. 1999. Measuring transformation error by RMS deviation. Internal Technical Report TR99MJ1, FMRIB Centre, University of Oxford. Available at www.fmrib.ox.ac.uk/analysis/techrep for downloading.
% Input:
% 	RealignmentParameterFile  -   The realignment parameter file for a given participant generated by SPM. E.g., rp***.txt
%   ReferenceImage            -   The reference image for realignment (usually the first time point (one-pass) or the mean image after an initial motion correction (two-pass))
% Output:
%	rel_rms      -   relative RMS (FD Jenkinson)
%	abs_rms      -   absolute RMS
%-----------------------------------------------------------
% Written by YAN Chao-Gan 120930.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com



%RP=load(RealignmentParameterFile);

% if necessary convert RP to SPM format
if sum(rp_indx.tra) ~= 6  
    RPtmp = nan(size(RP));
    RPtmp(:,1:3) = RP(:,4:6);
    RPtmp(:,4:6) = RP(:,1:3);
    RP = RPtmp; clear RPtmp;
end

switch lower(rp_rot_unit)
    case {'rad'}
        %do nothing
    case {'deg'}
        %converting to rad
        RP(:,4:6) = RP(:,4:6)*pi/180;         
end
%-------------------------------------

rmax = 80.0; %The default radius (as in FSL) of a sphere represents the brain


nTimePoint=size(RP,1);
sinq1=sin(RP(:,4));
sinq2=sin(RP(:,5));
sinq3=sin(RP(:,6));
cosq1=cos(RP(:,4));
cosq2=cos(RP(:,5));
cosq3=cos(RP(:,6));

%[RefData, RefHead] = rest_ReadNiftiImage(ReferenceImage,1);
RefHead = spm_vol(ReferenceImage);
center = RefHead.mat*([0.5*(RefHead.dim(1));0.5*(RefHead.dim(2));0.5*(RefHead.dim(3));1]);
center = center(1:3); %Get the coordinate for the center

abs_rms = zeros(nTimePoint,1);
for t=1:nTimePoint

    M1=[1       0        0     0;...
        0    cosq1(t)  sinq1(t)  0;...
        0    -sinq1(t) cosq1(t)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t)  0    sinq2(t)     0;...
        0        1       0        0;...
        -sinq2(t) 0    cosq2(t)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t)   sinq3(t)   0     0;...
        -sinq3(t)  cosq3(t)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t,1);...
        0    1     0     RP(t,2);...
        0    0     1     RP(t,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform=MT*M1*M2*M3;
    
    MA1=eye(4);
    MA2=(M_RigidBodyTransform);
    
    M = MA1*inv(MA2) - eye(4);
    
    A = M(1:3,1:3);
    
    T = M(1:3,4);
    
    abs_rms(t) = sqrt(rmax*rmax/5*trace(A'*A) + (T+A*center)'*(T+A*center));
end


rel_rms = zeros(nTimePoint-1,1);
for t=2:nTimePoint
    M1=[1       0        0     0;...
        0    cosq1(t)  sinq1(t)  0;...
        0    -sinq1(t) cosq1(t)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t)  0    sinq2(t)     0;...
        0        1       0        0;...
        -sinq2(t) 0    cosq2(t)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t)   sinq3(t)   0     0;...
        -sinq3(t)  cosq3(t)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t,1);...
        0    1     0     RP(t,2);...
        0    0     1     RP(t,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform=MT*M1*M2*M3;
    
    
    M1=[1       0        0     0;...
        0    cosq1(t-1)  sinq1(t-1)  0;...
        0    -sinq1(t-1) cosq1(t-1)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t-1)  0    sinq2(t-1)     0;...
        0        1       0        0;...
        -sinq2(t-1) 0    cosq2(t-1)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t-1)   sinq3(t-1)   0     0;...
        -sinq3(t-1)  cosq3(t-1)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t-1,1);...
        0    1     0     RP(t-1,2);...
        0    0     1     RP(t-1,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform_1=MT*M1*M2*M3;
    
    MA1=(M_RigidBodyTransform_1);
    MA2=(M_RigidBodyTransform);
    
    M = MA1*inv(MA2) - eye(4);
    
    A = M(1:3,1:3);
    
    T = M(1:3,4);
    
    rel_rms(t-1) = sqrt(rmax*rmax/5*trace(A'*A) + (T+A*center)'*(T+A*center));
    
end

rel_rms=[0;rel_rms]; %The FD_Jenkinson at time point t means the movement from time point t-1 to time point t. (Put the FD_Jenkinson for the first time point to "0".)



return
end