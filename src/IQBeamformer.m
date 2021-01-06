% Sua Bae
% 
% 2017-12-21
%   convex, conventional focusing
% 
% 2017-12-19 (LINEAR & CONVEX, PLANE WAVE & DIVERGING WAVE)
%   update for using stG_rt
%   include DSC for convex
%   include diverging wave beamforming
%
% 2017-08-28
%  IQ beamformer for plane wave imaging
%  lienar & convex
% 
function IQBeamformer(sRfDir, sBfDir, stRFInfo, stBFpm)
%vIQData, stRFInfo, stBFpm

    if ~strcmp(stBFpm.sVersion,'2.0')
        error('version mismatch');
    end
    
    %% 1. BF Parameters Download
    display('    1. BF Parameters Download');
    
    % Transducer
    stTrans             = stRFInfo.stTrans;
    nTransEleNum_x      = stTrans.nNumEle_x; % num of elements in lateral
    nTransElePitch_x    = stTrans.nPitch_x; % [meter]
    nTransRadius        = stTrans.nRadius; % [meter]
    
    % stTx
    stTx            = stRFInfo.stTx;    
        
    % sampling freq
    nFs             = stRFInfo.nFs; % sampling frequency before demodulation [Hz]     
    
    % Decimation
    nDeciFactor     = stBFpm.nDeciFactor; % factor of demcimation 
    
    % Demodulation
    nDemodFreq      = stBFpm.nDemodFreq; % demodualtion frequency
    nDemodLPFFreq   = stBFpm.nDemodLPFFreq; % demodulation LPF frequency
    nDemodLPFTap    = stBFpm.nDemodLPFTap; % demodulation LPF tap
    
    % IQ beamformer
    nChannel_x      = stTrans.nNumEle_x; % num. of channels
    nSoundSpeed     = stRFInfo.nSoundSpeed; % [m/s]
    nSec2Rad_iq     = 2*pi*stBFpm.nDemodFreq; % 2*pi*f0: for phase rotation (phase = 2*pi*f0*delay)
    
    % BF information 
    nRxFnum         = stBFpm.nRxFnum;
    sRxApodWindow   = stBFpm.sRxApodWindow;  
    nTxOffsetDelay_sec = stBFpm.nTxOffsetDelay_sec; % Tx offset delay [sec]

    
    sTransType = stTrans.sType;
    if strcmp(sTransType,'Linear')
        [mP_z, mP_x] = ndgrid(stBFpm.stG.aZ, stBFpm.stG.aX);
    elseif strcmp(sTransType,'Convex')
        [mR, mT] = ndgrid(stBFpm.stG_rt.aR, stBFpm.stG_rt.aT);
        [mP_x, mP_z] = rt2xz(stBFpm.stG_rt, mR, mT);
    else
        error('undefined');
    end
    

    %%   2. Calculate the Elements Coordinates & Position of origin of scanline 
    display('    2. Calculate the Elements Coordinates & Position of origin of scanline ');
    % the center of the transducer = the origin (0,0)
    
    %%%         1) Tx element position
    aEleIdx_x = -(nTransEleNum_x-1)/2 : 1 : (nTransEleNum_x-1)/2;        
    switch sTransType
        case 'Linear'
            aEl_x = aEleIdx_x*nTransElePitch_x; % x-pos of elements
            aEl_z = zeros(size(aEl_x)); % z-pos of elements
        case 'Convex'
            aEl_a = aEleIdx_x*nTransElePitch_x/nTransRadius/pi*180; % Transducer element theta position [degree] % size = [nChannel]
            aEl_x = nTransRadius*sind(aEl_a);
            aEl_z = nTransRadius*cosd(aEl_a);   
    end      

    %%%         2) Rx channel position
    aChanIdx_x = -(nChannel_x-1)/2 : 1 : (nChannel_x-1)/2;        
    switch sTransType
        case 'Linear'
            aCh_x = aChanIdx_x*nTransElePitch_x; % x-pos of elements
            aCh_z = zeros(size(aCh_x)); % z-pos of elements
        case 'Convex'
            aCh_a = aChanIdx_x*nTransElePitch_x/nTransRadius/pi*180; % Transducer element theta position [degree] % size = [nChannel]
            aCh_x = nTransRadius*sind(aCh_a);
            aCh_z = nTransRadius*cosd(aCh_a);   
    end

    %%   3. Beamforming & compounding
    display('    3. Beamforming');    
    
    if strcmp(stTx.sType,'CF')
        mDiff = bsxfun(@minus, stBFpm.stG_rt.aT, stTx.mFocalPos(:,5)'); % theta of BF - theta of TX, size: (stG_rt.nTdim x stTx.nNum)
        [~, aTxIdx] = min(abs(mDiff),[],2); % data index (txidx) for reconstruction of each BF scanline, size: (stG_rt.nTdim x 1)
    end
    
    mIQBFOut = zeros(size(mP_x));
    mCompNum = zeros(size(mP_x));
    if(stBFpm.bPlot); figure; end;
    for txidx = 1:stTx.nNum
                
        % Tx info 
        if strcmp(stTx.sType,'PW')
            nTxAngle = stTx.aPwAngle_deg(txidx); % [degree]
            display(['    ' num2str(txidx) '-th angle is being processed: ' num2str(nTxAngle) 'deg']);        
        elseif strcmp(stTx.sType,'DW')
            aVSPos   = stTx.mVSPos(txidx,:); % (x,y,z) of virtual source position  [meter]  
            display(['    ' num2str(txidx) '-th virtual source data is being processed: x = ' num2str(aVSPos(1)*1e3) ' mm, z = ' num2str(aVSPos(3)*1e3) ' mm']);   
        elseif strcmp(stTx.sType,'CF')
            aFocalPos   = stTx.mFocalPos(txidx,:); % (x,z,z0,r,tht) of virtual source position [meter, meter, meter, meter, deg]
            aThtIdx_bf = find(aTxIdx == txidx); % theta indices for reconstruction using txidx-th data
            display(['    ' num2str(txidx) '-th scanline data (' num2str(aFocalPos(5)) ' deg) is being processed for ' num2str(aThtIdx_bf') '-th BF scanline (' num2str(stBFpm.stG_rt.aT(aThtIdx_bf)') ' deg)']);   
        end
        
        % load data
        load([sRfDir 'mRcvData_tx_' num2str(txidx) '.mat']);
        
        %%%     0) Demodulation   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mIQData = QDM(mRcvData, stRFInfo.nFs, nDemodFreq, nDemodLPFFreq, nDemodLPFTap, nDeciFactor);     
        nIQSample = size(mIQData,1); % num. of samples

    
        %%%     1) Select BF points reconstructed by using the tidx-th data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
        if strcmp(stTx.sType, 'CF')
            % Select Beamforming Points
            mLogic_PR = zeros(size(mP_x,1),size(mP_x,2));
            mLogic_PR(:,aThtIdx_bf) = 1;
            mLogic_PR = logical(mLogic_PR);
            
        elseif strcmp(stTx.sType, 'PW')||strcmp(stTx.sType, 'DW')
            
            % leftmost and rightmost active elements
            aApod = TxApodGen(stTrans, stTx, txidx);
            lftmidx = find(aApod,1,'first'); 
            rftmidx = find(aApod,1,'last');
            nLftm_x = aEl_x(lftmidx);
            nRtm_x  = aEl_x(rftmidx);
            nLftm_z = aEl_z(lftmidx);
            nRtm_z  = aEl_z(rftmidx);       
            if strcmp(stTx.sType,'PW')
                % plane wave propagation region
                mLogic_PR = (mP_x >= (nLftm_x + (mP_z-nLftm_z)*tand(nTxAngle))) & (mP_x <= (nRtm_x + (mP_z-nRtm_z)*tand(nTxAngle)));
            elseif strcmp(stTx.sType,'DW')
                nAng_L  = atand( (nLftm_x-aVSPos(1))./(nLftm_z-aVSPos(3)) );% left diverging angle
                nAng_R  = atand( (nRtm_x-aVSPos(1))./(nRtm_z-aVSPos(3)) );% left diverging angle
                % diverging wave propagation region
                mAng_vs2bf = atand( (mP_x-aVSPos(1))./(mP_z-aVSPos(3)) ); % angle between virtual source and BF point   
                mLogic_PR = (mAng_vs2bf >= nAng_L) & (mAng_vs2bf <= nAng_R);
            end

        else
            error('undefined');
        end

        % Select Beamforming Points
        aBF_x = mP_x(mLogic_PR);
        aBF_z = mP_z(mLogic_PR);
        
        %%%     2) Beamforming  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Tx delay
        if strcmp(stTx.sType,'PW')
            aTxDist = aBF_x*sind(nTxAngle) + aBF_z*cosd(nTxAngle); % [meter] Distance from the PW at origin to target [meter] 
        elseif strcmp(stTx.sType,'DW')
            aTxDist = sqrt((aVSPos(1)-aBF_x).*(aVSPos(1)-aBF_x) + (aVSPos(3)-aBF_z).*(aVSPos(3)-aBF_z)); % [meter] Distance from the virtual source to target [meter] 
        elseif strcmp(stTx.sType,'CF')
            % multibeam delay type: CASE1 (the same tx delay for all multibeams)
            switch sTransType
                case 'Linear'
                    aTxDist = aBF_z; % [meter] Distance from the x-axis to target [meter] 
                case 'Convex'
                    aTxDist = sqrt(aBF_x.^2 + aBF_z.^2); % [meter] Distance from the origin to target [meter] 
            end
        end
        aTxDelay = aTxDist/nSoundSpeed; % [sec]
        if strcmp(stRFInfo.sStartAcqTime,'placed at origin')
            aTxDelay = aTxDelay; % do  nothing
        elseif strcmp(stRFInfo.sStartAcqTime,'first firing time')
%             switch sTransType
%                 case 'Linear'
%                     if strcmp(stTx.sType,'PW')
%                         aTxDelay = aTxDelay + max(abs(stRFInfo.stTxDelay(txidx).aDelay))/2;
%                     elseif strcmp(stTx.sType,'DW')
%                         aTxDelay = aTxDelay - abs(aVSPos(3));
%                     end
%                 case 'Convex'
                    aTxDelay = aTxDelay - stRFInfo.stTx.aDelayOffset(txidx);
%             end
        else
            error('undefined');
        end
        aTxDelay = aTxDelay + stBFpm.nTxOffsetDelay_sec;   
        
        % Rx Aperture size
        switch sTransType
            case 'Linear'
                aAptSize_m = aBF_z/nRxFnum + nTransElePitch_x; % plus nTransElePitch_x to prevent nApod being NaN
            case 'Convex'
                aRadial_distance_from_Txdcr = sqrt(aBF_x.^2 + aBF_z.^2) - nTransRadius; % aRadial_distance_from_Txdcr to each BF point
                aAptSize_m = aRadial_distance_from_Txdcr/nRxFnum + nTransElePitch_x; % plus nTransElePitch_x to prevent nApod being NaN
                % nMaxAptSize_m = (30/180*pi)*nTransRadius*2; % maximum aperture size is limited by acceptance angle
                nMaxAptSize_m = nTransRadius*sind(30)*2; % maximum aperture size is limited by acceptance angle
                aAptSize_m = min(aAptSize_m, nMaxAptSize_m);
        end
        
        % Interpolation & summation
        aInph = zeros(size(aBF_x)); % beamformed data (only within LWR)
        aQuad = zeros(size(aBF_x)); % beamformed data (only within LWR)
        for cidx = 1 : nChannel_x

            % 1) Apodization window
            if(nRxFnum~=0)
                switch sTransType
                    case 'Linear'
                        aDistance_m = abs(aBF_x - aCh_x(cidx));  % from x-pos of BF point to each element      
                    case 'Phased'
                        aDistance_m = abs(aCh_x(cidx)*ones(numel(aBF_data),1));  % from center to each element       
                    case 'Convex'
                        aBF_a = atand(aBF_x./aBF_z);
                        aDistance_m = abs(nTransRadius*sind(aCh_a(cidx) - aBF_a));  % from specific point (intersection of scanline and array) to each element              
                end
                
                aApod = ApodGen(aDistance_m,aAptSize_m,sRxApodWindow);
%                 idx_a_certain_BFPoint = 154;
%                 aApod_for_a_certain_BFpoint(cidx) = aApod(idx_a_certain_BFPoint);
            else
                aApod = ones(size(aAptSize_m));
            end            

            % 2) Rx delay, Round trip delay
            aRxDelay = sqrt( (aBF_x - aCh_x(cidx)).^2 + (aBF_z - aCh_z(cidx)).^2 ) /nSoundSpeed; % [sec]
            aTimeDelay = aTxDelay + aRxDelay + nTxOffsetDelay_sec; % [sec]

            % 3) calc address
            aAddress = aTimeDelay*(nFs/nDeciFactor);     
            aAddress = max(min(aAddress, nIQSample-1),0);% saturation, if Add=0, summation X
            aLogic   = ( aAddress > 0 ).* (aAddress < nIQSample-1);

            % 4) interpolate data
            sInterpMode = 'Spline'; % 'Spline' or 'InterpFilter';
            switch sInterpMode
                case 'Spline'
                    aInph_intp = interpn(0:1:(nIQSample-1), real(mIQData(:,cidx)), aAddress,'spline');
                    aQuad_intp = interpn(0:1:(nIQSample-1), imag(mIQData(:,cidx)), aAddress,'spline');
                case 'InterpFilter'
                    if (nIntpFactor == 1)
                        aInphChData_intp = real(mIQData(:,cidx));
                        aQuadChData_intp = imag(mIQData(:,cidx));
                    else
                        aInphChData_intp = InterpFunc(nIntpFactor, 31, real(mIQData(:,cidx)));
                        aQuadChData_intp = InterpFunc(nIntpFactor, 31, imag(mIQData(:,cidx)));
                        % aInphChData_intp = InterpFunc(nIntpFactor, 129, real(vIQData(:,cidx,scidx)));
                        % aQuadChData_intp = InterpFunc(nIntpFactor, 129, imag(vIQData(:,cidx,scidx)));
                    end
                    aInph_intp = aInphChData_intp( round(aAddress*nIntpFactor) )';
                    aQuad_intp = aQuadChData_intp( round(aAddress*nIntpFactor) )';
            end

            % 5) compansate phase and envelop
%             aCompen_phase = (aTimeDelay-2*aBF_z/nSoundSpeed)*nSec2Rad_iq; %[sec]->[rad]
            aCompen_phase = (aTimeDelay)*nSec2Rad_iq; %[sec]->[rad]
            aInph_cps = aInph_intp.*cos(aCompen_phase) - aQuad_intp.*sin(aCompen_phase);
            aQuad_cps = aInph_intp.*sin(aCompen_phase) + aQuad_intp.*cos(aCompen_phase);
            
            % 6) accumulate
            aInph = aInph + aInph_cps.*aLogic.*aApod;  
            aQuad = aQuad + aQuad_cps.*aLogic.*aApod;   

        end
        
        %%%     3) Compounding  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mIQBFOut(mLogic_PR) = mIQBFOut(mLogic_PR) + (aInph + 1j*aQuad);
        mCompNum(mLogic_PR) = mCompNum(mLogic_PR) + 1;
        if(stBFpm.bSaveEachAngle)
            mIQBFOut_eachAng                = zeros(size(mP_x));
            mIQBFOut_eachAng(mLogic_PR)    = (aInph + 1j*aQuad); 
            stBfData_b4c.mBfData            = mIQBFOut_eachAng;
            stBfData_b4c.stRFInfo            = stRFInfo;
            stBfData_b4c.stBFpm              = stBFpm;
            save([sBfDir, 'stBfData_b4c_tx_' num2str(txidx)],'stBfData_b4c');
        end
        if(stBFpm.bPlot); imagesc(mP_x(1,:),mP_z(:,1),db(abs(mIQBFOut)/max(max(abs(mIQBFOut)))));axis equal; axis tight; caxis([-100 0]);pause(0.1); 
        end;

    
    end
    
    %%%     4) Normalization   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    mIQBFOut = mIQBFOut./mCompNum;
    mIQBFOut(isnan(mIQBFOut)) = 0;
    
    
    %%   4. DSC to catesian grid (only for Convex data)
    
    switch sTransType
        case 'Linear'
            mBfData = abs(mIQBFOut);
            
        case 'Convex'
            %%%  FROM RT GRID TO XZ GRID
            
            mBfData_rt = abs(mIQBFOut);     % (rt grid)    

            stG_rt = stBFpm.stG_rt; % BF grid (rt grid)
            stG     = stBFpm.stG; % DSC grid (xz grid)

            mBfData = DSC_rt2xz(mBfData_rt, stG_rt, stG);
            
%             % (x,z) of Cartesian grid (xz grid)
%             [mZ_cts,mX_cts] = ndgrid(stG.aZ,stG.aX);
% 
%             % (r,t) of Cartesian grid
%             [aR_cts,aT_cts] = xz2rt(stG_rt, mX_cts(:), mZ_cts(:));
% 
%             % ndgrid of RT grid
%             [mR_rt, mT_rt] = ndgrid(stG_rt.aR,stG_rt.aT);
% 
%             % interpolate data to Cartesian grid
%             aBfData_cts = interpn(mR_rt, mT_rt, mBfData_rt, aR_cts, aT_cts);
%             aBfData_cts(isnan(aBfData_cts)) = 0;
% 
%             % reshape and export!
%             mBfData = reshape(aBfData_cts,stG.nZdim,stG.nXdim); %  (xz grid)

            % plot
            if(stBFpm.bPlot);figure;imagesc(stG.aX,stG.aZ0,db(abs(mBfData)/max(abs(mBfData(:)))));axis equal; axis tight; caxis([-100 0]);title('after dsc'); end;
        
            % save before DSC
            stBfData.mBfData_rt  = mBfData_rt; %  (xz grid)
    end
    
    %%   4. Save
    
    stBfData.mBfData     = mBfData;
    stBfData.mCompNum    = mCompNum; % when 'Convex', it is on RT grid
    stBfData.stRFInfo    = stRFInfo;
    stBfData.stBFpm      = stBFpm;
    save([sBfDir 'stBfData.mat'], 'stBfData');
    display(sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock));

end