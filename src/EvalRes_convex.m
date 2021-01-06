% Sua Bae
%   2017-11-21
%      using 'GetFWHM2D', 'InterpSlice'
%   
% -------IN/OUT-----------------------------------
% mData         : envelop data
% stG           : grid of volume data, type = [clBFGrid3D] or [clGrid3D]
% mPointPos     : x,y,z position of point targets for measurement, dim = (num of points x 3)
% sDataName     : for saving
% bPlot         : plot
%
% stRes             : lateral, elevational, axial resolution
% stRes.sResType    : 'lateral resolution', 'elevational resolution', 'axial resolution'
% stRes.sDataName   : sDataName (INPUT)
% stRes.mPointPos   : mPointPos;
% stRes.aFWHM       : mFWHM(:,1);
% ------------------------------------------------
function stRes = EvalRes_convex(mData, stG, aPos_scat, sDataName, bPlot)
    
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
    
    
    if bPlot; figure('position',[50 50 1000 300]); end;
    

    % get an interpolated image centered at the pidx-th scatterer position
    stG_i = clGrid(1200, 800, ...
                   0.005e-3, 0.005e-3, ...
                   aPos_scat(1), aPos_scat(3));
    [mX_itp, mZ_itp] = ndgrid(stG_i.aX, stG_i.aZ);
    mData_itp = interpn(mX, mZ0, mData, mX_itp, mZ_itp);

    % find maximum position (it will be close to the scatterer position)
    [~, maxidx] = max(mData_itp(:));
    maxidx_x = idx2xidx(stG_i,maxidx);
    maxidx_z = idx2zidx(stG_i,maxidx);
    aPos_max(1) = stG_i.aX(maxidx_x);
    aPos_max(3) = stG_i.aZ(maxidx_z);
    
    % log comp
    nDR = 80;
    mLog_xz = LogCompression(mData_itp, nDR);

    % plot if needed
    if bPlot
        subplot(1,2,1); 
            imagesc(stG_i.aX*1e3, stG_i.aZ*1e3, mLog_xz');
            axis equal; axis tight; colormap(gray); title(['lateral res::' sDataName]); xlabel('x (mm)');ylabel('z (mm)');
        subplot(1,2,2); 
            imagesc(stG_i.aX*1e3, stG_i.aZ*1e3, mLog_xz');
            axis equal; axis tight; colormap(gray); title(['axial res::' sDataName]); xlabel('x (mm)');ylabel('z (mm)');
    end 

    % lateral resolution
    [aFWHM(1), x1, z1, x2, z2] = GetFWHM2D(stG_i.aX, stG_i.aZ, mLog_xz);
    % axial resolution
    [aFWHM(2), z3, x3, z4, x4] = GetFWHM2D(stG_i.aZ, stG_i.aX, mLog_xz');

    % plot if needed
    if bPlot
        subplot(1,2,1); 
            hold on; contour(stG_i.aX*1e3, stG_i.aZ*1e3, mLog_xz', [-6 -6], 'edgecolor','r');  
            scatter([x1,x2]*1e3,[z1,z2]*1e3,'x');hold off; text(stG_i.aX(1)*1e3+.5,stG_i.aZ(end)*1e3-.5,['FWHM=' num2str(round(aFWHM(1)*1e5)*1e-2) 'mm'],'Color','w');
        subplot(1,2,2); 
            hold on; contour(stG_i.aX*1e3, stG_i.aZ*1e3, mLog_xz', [-6 -6], 'edgecolor','r');  
            scatter([x3,x4]*1e3,[z3,z4]*1e3,'x');hold off; text(stG_i.aX(1)*1e3+.5,stG_i.aZ(end)*1e3-.5,['FWHM=' num2str(round(aFWHM(2)*1e5)*1e-2) 'mm'],'Color','w');                
        drawnow;
    end 

    stRes(1).sResType   = 'lateral resolution';
    stRes(1).sDataName  = sDataName;
    stRes(1).aPos_scat  = aPos_scat;
    stRes(1).nFWHM      = aFWHM(1);

    stRes(2).sResType   = 'axial resolution';
    stRes(2).sDataName  = sDataName;
    stRes(2).aPos_scat  = aPos_scat;
    stRes(2).nFWHM      = aFWHM(2);


end
