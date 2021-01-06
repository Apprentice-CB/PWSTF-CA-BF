% Sua Bae
%
% 2017-12-26
%  update variable names, inputs
%
% 2015-12-17
% calculate Contrast of cyst
%
% -------IN/OUT-----------------------------------
% mData         : image data (envelop data)
% mNoise         : noise data (envelop data)
% stG           : grid of volume data, type = [clBFGrid_rt]
% sUnit         : 'pixel', 'mm' : unit of position
% sShape_dark   : 'rect', 'circ'
% aPos          : [x of center of rectangle, z of center of rectangle, width, hight]
% sDataName     : data name
% bPlot         : plot
%
% stCon         : .nCon (contrast in dB), .nE_dark, .E_bright
% ------------------------------------------------
%
% if     'rect' aPos: 
% elseif 'circ' aPos: [x of center, z of center, radius]
% elseif 'rectwh' (rect with a hole) aPos: [x of center of rect, z of center of rect, 
%                                           width, hight, 
%                                           x of center of hole, z of center of hole, radius of hole]
%
function stSNR = EvalSNR_convex(mData, mNoise, stG, sUnit, aPos, sDataName, bPlot)

    % use envelop
    mData = abs(mData);
    
    % to ensure mData has the size of 'nZdim x nXdim'
    if (size(mData,1)==stG.nZdim)&&(size(mData,2)==stG.nXdim)
        % do nothing
    elseif (size(mData,1)==stG.nXdim)&&(size(mData,2)==stG.nZdim)
        mData = mData';
    else
        error('size mismatch!!');
    end  
    
    % use envelop
    mNoise = abs(mNoise);
    
    % to ensure mData has the size of 'nZdim x nXdim'
    if (size(mNoise,1)==stG.nZdim)&&(size(mNoise,2)==stG.nXdim)
        % do nothing
    elseif (size(mNoise,1)==stG.nXdim)&&(size(mNoise,2)==stG.nZdim)
        mNoise = mNoise';
    else
        error('size mismatch!!');
    end  
    
    % grid 
    aXaxis = stG.aX;
    aZaxis = stG.aZ0;

    if bPlot
        figure('position',[100 100 1500 400]);
        subplot(1,3,1); % total image region
            mLogOut_80dB = LogCompression(mData, 80);
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mLogOut_80dB);
    end
    
    % calc indices
    if strcmp(sUnit,'mm')
        x = aPos(1); % x pos of center of rect
        z = aPos(2); % z pos of center of rect
        w = aPos(3); % width
        h = aPos(4); % hight 
        if bPlot; subplot(1,3,1); rectangle('Position', [x-w/2, z-h/2, w, h]*1e3, 'EdgeColor', [1 1 0]); end; % plot

        % convert 'mm' to 'pixel'            
        xs_p = round(interp1(aXaxis,1:size(mData,2),x-w/2)); % x pixel of upper-left corner
        xe_p = round(interp1(aXaxis,1:size(mData,2),x+w/2)); % x pixel of lower-right corner
        zs_p = round(interp1(aZaxis,1:size(mData,1),z-h/2)); % z pixel of upper-left corner
        ze_p = round(interp1(aZaxis,1:size(mData,1),z+h/2)); % z pixel of lower-right corner   

    elseif strcmp(sUnit,'pixel')
        x_p = aPos(1); % x pixel of top-left corner
        z_p = aPos(2); % z pixel of top-left corner
        w_p = aPos(3); % width pixel
        h_p = aPos(4); % hight pixel
    
        xs_p = round(x_p-w_p/2); % x pixel of upper-left corner
        xe_p = round(x_p+w_p/2); % x pixel of lower-right corner
        zs_p = round(z_p-h_p/2); % z pixel of upper-left corner
        ze_p = round(z_p+h_p/2); % z pixel of lower-right corner
    else
        error('The string sUnit should be mm or pixel');
    end

    % widnowing
    mSignal_w = mData(zs_p:ze_p,xs_p:xe_p); % 2-D windowed sample
    aSignal   = mSignal_w(:); % rearrange into 1-D
    
    mNoise_w = mNoise(zs_p:ze_p,xs_p:xe_p); % 2-D windowed sample
    aNoise   = mNoise_w(:); % rearrange into 1-D
    

    if bPlot     
        subplot(1,3,2);
            mPlot_sig = zeros(size(mData,1),size(mData,2));
            mPlot_sig(zs_p:ze_p,xs_p:xe_p) = mLogOut_80dB(zs_p:ze_p,xs_p:xe_p); % only dark region
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mPlot_sig);    
        subplot(1,3,3);
            mPlot_noise = zeros(size(mNoise,1),size(mNoise,2));
            mNoise_dB = db(mNoise/max(mData(:)));
            mPlot_noise(zs_p:ze_p,xs_p:xe_p) = mNoise_dB(zs_p:ze_p,xs_p:xe_p); % only dark region
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mPlot_noise);  caxis([-150 0]); 
    end

    
    if bPlot
        subplot(1,3,1);
            set(gca,'FontName','Times New Roman', 'FontSize', 12);  
            axis equal; axis([stG.aX(1) stG.aX(end) stG.aZ0(1) stG.aZ0(end)]*1e3); %axis([-30 30 -30 30 0 70]);
            xlabel('{\itx} (mm)'); ylabel('{\itz} (mm)'); colormap(gray);
            title([sDataName]);
        subplot(1,3,2);
            set(gca,'FontName','Times New Roman', 'FontSize', 12);  
            axis equal; axis([stG.aX(1) stG.aX(end) stG.aZ0(1) stG.aZ0(end)]*1e3); %axis([-30 30 -30 30 0 70]);
            xlabel('{\itx} (mm)'); ylabel('{\itz} (mm)'); colormap(gray);
            title(['signal (80dB)']);
        subplot(1,3,3);
            set(gca,'FontName','Times New Roman', 'FontSize', 12);  
            axis equal; axis([stG.aX(1) stG.aX(end) stG.aZ0(1) stG.aZ0(end)]*1e3); %axis([-30 30 -30 30 0 70]);
            xlabel('{\itx} (mm)'); ylabel('{\itz} (mm)'); colormap(gray);
            title(['noise (150dB)']);
        drawnow;  
    end    
    
    %% SNR
    
    nE_signal  = sum(aSignal.^2);
    nE_noise   = sum(aNoise.^2);
    nSNR        = 10*log10(nE_signal/nE_noise);
    stSNR.sDataName  = sDataName;
    stSNR.aPos       = aPos;     
    stSNR.nE_signal  = nE_signal;    % measured energy of signal
    stSNR.nE_noise   = nE_noise;     % measured energy of noise
    stSNR.nSNR       = nSNR;
    
    


end
