% Sua Bae
% 2017-12-21
%  import stTrans from stRFInfo
%
% 2017-11-28
%   view angle: 30 deg for xz plane, 23.3584 deg for yz plane
%
% 2017-11-09
%
%  using clBFGrid3D
%
function stBFpm = SetIQBFParam(stRFInfo)
    
    stTrans = stRFInfo.stTrans;
    
    stBFpm.sVersion = '2.0';
    stBFpm.bPlot = 1;
    stBFpm.bSaveEachAngle = 0;
    
    %% 1. Demodulation (IQ Data!)
    stBFpm.nDeciFactor      = 2; % factor of decimation
    stBFpm.nDemodFreq       = stTrans.nFc; % demodulation frequency
    stBFpm.nDemodLPFFreq    = min(stRFInfo.nFc/2, stRFInfo.nFs/2/stBFpm.nDeciFactor); % demodulation LPF frequency % "nFc/2" (100% band limit) and "nFs/2/nDeciFactor" (band limit for decimation)
    stBFpm.nDemodLPFTap     = 127; % demodulation LPF tap
    
    %% 2. Beamforming grid
    %%%%--  Dim, interval, offset
    nMultiBeam = 4;
    % num of BF points
    nRdim = round(stRFInfo.nSample/stBFpm.nDeciFactor);
    nTdim = stTrans.nNumEle*nMultiBeam;
    % interval of BF points
    dr = 1/(stRFInfo.nFs)*(stRFInfo.nSoundSpeed)/2*stBFpm.nDeciFactor; % [m] dr = sample interval of RF data 
    dt = stTrans.nPitchAngle_x/nMultiBeam; % [deg]
    % offset of BF points
    nRos = stTrans.nRadius; % [m]
    nTos = 0; % [deg] 

    %%%%-- BF grid gen using clBFGrid_rt (for beamforming)
    stBFpm.stG_rt = clBFGrid_rt(nRdim, nTdim, dr, dt, nRos, nTos);
    display(['ROI size: R0: ' num2str(stBFpm.stG_rt.aR0(1)*1e3) '~' num2str(stBFpm.stG_rt.aR0(end)*1e3) 'mm, ' ...
                       'T: ' num2str(stBFpm.stG_rt.aT(1)) '~' num2str(stBFpm.stG_rt.aT(end)) 'deg']);
    display(['          dr: ' num2str(dr*1e3) 'mm (' num2str(dr/(stRFInfo.nSoundSpeed/stRFInfo.nFc)) ' lambda), ' ...
                       'dt: ' num2str(dt) 'deg']);  
                    
   %% 3. DSC grid
    %%%%--  Dim, interval
    % interval of points
    dx = stTrans.nPitch_x/nMultiBeam; % [m] 
    dz = 1/(stRFInfo.nFs)*(stRFInfo.nSoundSpeed)/2*stBFpm.nDeciFactor; % [m]
    % num of points
    nROISize_z = 160e-3 + stTrans.nRadius*(1-cosd(stTrans.nMaxTheta));
    nROISize_x = (160e-3 + stTrans.nRadius)*sind(stTrans.nMaxTheta)*2;
    nXdim = round(nROISize_x/dx);
    nZdim = round(nROISize_z/dz);
    
    %%%%-- Dsc grid gen using clDscGrid_convex (for dsc)
    stBFpm.stG = clDscGrid_convex(nXdim, nZdim, dx, dz, stTrans.nMaxTheta, stTrans.nRadius);
    display(['ROI size: X: ' num2str(stBFpm.stG.aX(1)*1e3) '~' num2str(stBFpm.stG.aX(end)*1e3) 'mm, ' ...
                       'Z0: ' num2str(stBFpm.stG.aZ0(1)*1e3) '~' num2str(stBFpm.stG.aZ0(end)*1e3) 'mm']);
    display(['          dx: ' num2str(dx*1e3) 'mm (' num2str(dx/(stRFInfo.nSoundSpeed/stRFInfo.nFc)) ' lambda), ' ...
                       'dz: ' num2str(dz*1e3) 'mm (' num2str(dz/(stRFInfo.nSoundSpeed/stRFInfo.nFc)) ' lambda)']);  


    
    %% 3. RX
    stBFpm.nRxFnum          = 1.0;
    stBFpm.sRxApodWindow    = 'tukey50';
    
    %% 4. Tx offset delay
    stBFpm.nTxOffsetDelay_sec   = 0;        
         
    
    
end