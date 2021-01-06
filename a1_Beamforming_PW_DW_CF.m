% Sua Bae
% For convex PWSTF journal paper
%
% 2017-12-21
%   Add CF
%
% 2017-12-18 
%   Beamforming for PWSTF, DWSTF in convex array imaging
%
% 
%
clc;
clear;
close all;
addpath('src');

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transducer spec :: Convex
load('stTrans_C5-2');

for phantom_idx = 9:11

    for acqType_idx = 1:4
        % 1. PWSTF 
        % 2. DWSTF centerRes
        % 3. CF Fnum2p4
        % 4. CF Fnum5
        % 5. PWSTF EA
        switch acqType_idx        
            case {1,2,5};   nDN = 32;
            case {3,4};     nDN = 256;
        end
        for dataNumber_idx = nDN
            
            %%% get path
            stD = clDir(phantom_idx, acqType_idx, dataNumber_idx, '..\Data\');   
            sDataDir = stD.sDataDir
            
            %%% Beamformer parameter setup
            load([sDataDir 'RcvData\stRFInfo.mat']);        
            stBFpm = SetIQBFParam(stRFInfo);

            %%% Create folders and save stBFpm
            sRfDir = [sDataDir 'RcvData\'];
            sBfDir = [sDataDir 'BfData\' ];
            if(~isdir(sBfDir)); mkdir(sBfDir); end  
            save([sBfDir 'stBFpm.mat'], 'stBFpm');
            
            %%% Beamforming
            IQBeamformer(sRfDir, sBfDir, stRFInfo, stBFpm);
        end
    end
end

                 
                 