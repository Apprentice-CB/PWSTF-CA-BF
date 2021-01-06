% Sua Bae
%   2017-11-24
%      Measure contrast of PSF!!
%   
% -------IN/OUT-----------------------------------
% vData         : volume data (envelop data)
% stG           : grid of volume data, type = [clBFGrid3D] or [clGrid3D]
% aPos          : x,y,z positions of a point target, dim = (1 x 3)
% nRadius_sph   : radius of sphere
% sDataName     : data name
% bPlot         : plot
%
% stCon         : .nCon (contrast in dB)
% ------------------------------------------------
function stCon = EvalConPSF_convex(mData, stG, aPos, nSphereRadius, sDataName, bPlot)
    
    % grid    
    [mX, mZ0] = ndgrid(stG.aX,stG.aZ0);
        
    % to ensure mData has the size of 'nXdim x nZdim'
    if (size(mData,1)==stG.nXdim)&&(size(mData,2)==stG.nZdim)
        % do nothing
    elseif (size(mData,1)==stG.nZdim)&&(size(mData,2)==stG.nXdim)
        mData = mData';
    else
        error('size mismatch!!');
    end  
    
    % set total area
    mWithinWindow   = true(stG.nXdim, stG.nZdim);
%     nXlen_w         = 30e-3;
%     nZlen_w         = 10e-3;
%     mWithinWindow   = ( (abs(vX-aPos(1)) < nXlen_w/2)  & (abs(vZ-aPos(3)) < nZlen_w/2) );
    
    % set spherical region
    mWithinSph      = sqrt( (mX-aPos(1)).^2 + (mZ0-aPos(3)).^2 ) < nSphereRadius;

    % contrast
    nPowerTotal  = sum(sum(mData(mWithinWindow).^2));
    nPowerSphere = sum(mData(mWithinSph).^2);
    nCon = 10*log10( (nPowerTotal - nPowerSphere) / nPowerTotal );  % proportion of outer region to total
    
    if bPlot
        figure('position',[100 100 1500 400]);
        subplot(1,3,1);
            %%% image recon   
            mLogOut_80dB = LogCompression(abs(mData), 80);
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mLogOut_80dB');       
            set(gca,'FontName','Times New Roman', 'FontSize', 12);  
            axis equal; axis([stG.aX(1) stG.aX(end) stG.aZ0(1) stG.aZ0(end)]*1e3); %axis([-30 30 -30 30 0 70]);
            xlabel('{\itx} (mm)'); ylabel('{\itz} (mm)');
            % sphere 
            hold on;
            x = nSphereRadius*1e3*cos(-pi:0.01:pi) + aPos(1)*1e3;
            z = nSphereRadius*1e3*sin(-pi:0.01:pi) + aPos(3)*1e3;
            plot(x,z,'Color','w','LineWidth',1);
            % window
            if ~all(mWithinWindow)
                x1 = (aPos(1) - nXlen_w/2)*1e3; x2 = (aPos(1) + nXlen_w/2)*1e3;
                z1 = (aPos(3) - nZlen_w/2)*1e3; z2 = (aPos(3) + nZlen_w/2)*1e3;
                line([x1,x2,x2,x1,x1],[z1,z1,z2,z2,z1]);
            end
            hold off;
            title([sDataName]);drawnow;
        subplot(1,3,2);
            %%% image recon  
            mWindowRegion = mData;
            mWindowRegion(~mWithinWindow) = 0; % only window region
            mWindowRegion(1) = max(mData(:)); % for dynamic range Xmax
            mLogOut_80dB = LogCompression(mWindowRegion, 80);            
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mLogOut_80dB');    
            set(gca,'FontName','Times New Roman', 'FontSize', 12);  
            axis equal; axis([stG.aX(1) stG.aX(end) stG.aZ0(1) stG.aZ0(end)]*1e3); %axis([-30 30 -30 30 0 70]);
            xlabel('{\itx} (mm)'); ylabel('{\itz} (mm)');
            title(['total region']);drawnow;
            clear mWindowRegion
        subplot(1,3,3);
            %%% 3D image recon              
            mSphereRegion = mData;
            mSphereRegion(~mWithinWindow) = 0;
            mSphereRegion(mWithinSph) = 0;
            mSphereRegion(1) = max(mData(:)); % for dynamic range Xmax
            mLogOut_80dB = LogCompression(mSphereRegion, 80);            
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mLogOut_80dB');    
            set(gca,'FontName','Times New Roman', 'FontSize', 12);  
            axis equal; axis([stG.aX(1) stG.aX(end) stG.aZ0(1) stG.aZ0(end)]*1e3); %axis([-30 30 -30 30 0 70]);
            xlabel('{\itx} (mm)'); ylabel('{\itz} (mm)');
            title(['background region']);drawnow;
            clear mSphereRegion        
    end    
    nCon
    stCon.sDataName     = sDataName;
    stCon.aPos          = aPos;
    stCon.nSphereRadius = nSphereRadius;
    if ~all(mWithinWindow); stCon.aWindowSize   = [nXlen_w, nZlen_w];
    else stCon.aWindowSize   = [0, 0, 0]; end;
    stCon.nPowerTotal   = nPowerTotal;
    stCon.nPowerSphere  = nPowerSphere;
    stCon.nCon          = nCon;

end
