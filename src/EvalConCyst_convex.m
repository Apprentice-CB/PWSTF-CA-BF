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
% stG           : grid of volume data, type = [clBFGrid_rt]
% sUnit         : 'pixel', 'mm' : unit of position
% aPos_dark     : position of cyst region 
% sShape_dark   : 'rect', 'circ'
% aPos_bright   : position of bkg region 
% sShape_bright : 'rect', 'circ', 'rectwh'
% sDataName     : data name
% bPlot         : plot
%
% stCon         : .nCon (contrast in dB), .nE_dark, .E_bright
% ------------------------------------------------
%
% if     'rect' aPos: [x of top-left corner, z of top-left corner, width, hight]
% elseif 'circ' aPos: [x of center, z of center, radius]
% elseif 'rectwh' (rect with a hole) aPos: [x of center of rect, z of center of rect, 
%                                           width, hight, 
%                                           x of center of hole, z of center of hole, radius of hole]
%
function stCon = EvalConCyst_convex(mData, stG, sUnit, aPos_dark, sShape_dark, aPos_bright, sShape_bright, sDataName, bPlot)

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
    
    % grid 
    aXaxis = stG.aX;
    aZaxis = stG.aZ0;

    if bPlot
        figure('position',[100 100 1500 400]);
        subplot(1,3,1); % total image region
            mLogOut_80dB = LogCompression(mData, 80);
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mLogOut_80dB);
    end
    %% Dark region
    if strcmp(sShape_dark,'rect')        
        
        if strcmp(sUnit,'mm')
            x = aPos_dark(1); % x pos of top-left corner
            z = aPos_dark(2); % z pos of top-left corner
            w = aPos_dark(3); % width
            h = aPos_dark(4); % hight            
            if bPlot; subplot(1,3,1); rectangle('Position', [x y w h]*1e3); end; % plot
            % convert 'mm' to 'pixel'            
            xs_p = round(interp1(aXaxis,1:size(mData,2),x-w/2)); % x pixel of upper-left corner
            xe_p = round(interp1(aXaxis,1:size(mData,2),x+w/2)); % x pixel of lower-right corner
            zs_p = round(interp1(aZaxis,1:size(mData,1),z-h/2)); % z pixel of upper-left corner
            ze_p = round(interp1(aZaxis,1:size(mData,1),z+h/2)); % z pixel of lower-right corner
            
        elseif strcmp(sUnit,'pixel')
            x_p = aPos_dark(1); % x pixel of top-left corner
            z_p = aPos_dark(2); % z pixel of top-left corner
            w_p = aPos_dark(3); % width pixel
            h_p = aPos_dark(4); % hight pixel
            if bPlot; subplot(1,3,1); rectangle('Position', [x_p z_p w_p h_p]*1e3); end; % plot
            
            xs_p = round(x_p-w_p/2); % x pixel of upper-left corner
            xe_p = round(x_p+w_p/2); % x pixel of lower-right corner
            zs_p = round(z_p-h_p/2); % z pixel of upper-left corner
            ze_p = round(z_p+h_p/2); % z pixel of lower-right corner
        else
            error('The string sUnit should be mm or pixel');
        end
        
        mDark = mData(zs_p:ze_p,xs_p:xe_p); % 2-D windowed sample
        aDark = reshape(mDark,1,numel(mDark)); % rearrange into 1-D
        
        if bPlot     
            subplot(1,3,2);
            mPlot_dark = zeros(size(mData,1),size(mData,2));
            mPlot_dark(zs_p:ze_p,xs_p:xe_p) = mLogOut_80dB(zs_p:ze_p,xs_p:xe_p); % only dark region
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mPlot_dark);    
        end

        
    elseif strcmp(sShape_dark,'circ')
        
        if strcmp(sUnit,'mm')        
            c_x = aPos_dark(1); % x pos of center of circle
            c_z = aPos_dark(2); % z pos of center of circle
            r = aPos_dark(3);   % radius of circle
            if bPlot; subplot(1,3,1); rectangle('Position', [c_x-r c_z-r 2*r 2*r]*1e3, 'Curvature', [1 1], 'EdgeColor', [1 1 1]); end; % plot
            [mZ, mX] = ndgrid(aZaxis,aXaxis);
            aIndex_w = find((mX-c_x).^2 + (mZ-c_z).^2 < r^2); % find indices within circle
            aDark = mData(aIndex_w);    
            
%             aTest = mImage;
%             aTest(aIndex) = 1000000;
            
        elseif strcmp(sUnit,'pixel')
            error('position for circle should be in mm (because of radius) ');  
        else
            error('The string sUnit should be mm or pixel');
        end       
        
        if bPlot     
            subplot(1,3,2);
            mPlot_dark = zeros(size(mData,1),size(mData,2));
            mPlot_dark(aIndex_w) = mLogOut_80dB(aIndex_w); % only dark region
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mPlot_dark);  caxis([-80 0]);  
        end
    
    else
        error('sShape_dark should be circ or rect');
    end
    

    %% Bright region
    if strcmp(sShape_bright,'rect')        
        
        if strcmp(sUnit,'mm')
            x = aPos_bright(1); % x pos of top-left corner
            z = aPos_bright(2); % z pos of top-left corner
            w = aPos_bright(3); % width
            h = aPos_bright(4); % hight            
            if bPlot; subplot(1,3,1); rectangle('Position', [x z w h]*1e3); end; % plot
            % convert 'mm' to 'pixel'            
            xs_p = round(interp1(aXaxis,1:size(mData,2),x-w/2)); % x pixel of upper-left corner
            xe_p = round(interp1(aXaxis,1:size(mData,2),x+w/2)); % x pixel of lower-right corner
            zs_p = round(interp1(aZaxis,1:size(mData,1),z-h/2)); % z pixel of upper-left corner
            ze_p = round(interp1(aZaxis,1:size(mData,1),z+h/2)); % z pixel of lower-right corner
            
        elseif strcmp(sUnit,'pixel')
            x_p = aPos_bright(1); % x pixel of top-left corner
            z_p = aPos_bright(2); % z pixel of top-left corner
            w_p = aPos_bright(3); % width pixel
            h_p = aPos_bright(4); % hight pixel
            if bPlot; subplot(1,3,1); rectangle('Position', [x_p z_p w_p h_p]*1e3); end; % plot
            
            xs_p = round(x_p-w_p/2); % x pixel of upper-left corner
            xe_p = round(x_p+w_p/2); % x pixel of lower-right corner
            zs_p = round(z_p-h_p/2); % z pixel of upper-left corner
            ze_p = round(z_p+h_p/2); % z pixel of lower-right corner
        else
            error('CalcCNR:: The string sUnit should be mm or pixel');
        end
        
        mBright = mData(zs_p:ze_p,xs_p:xe_p); % 2-D windowed sample
        aBright = reshape(mBright,1,numel(mBright)); % rearrange into 1-D
        
        
        
    elseif strcmp(sShape_bright,'circ')
        
        if strcmp(sUnit,'mm')        
            c_x = aPos_bright(1); % x pos of center of circle
            c_z = aPos_bright(2); % z pos of center of circle
            r = aPos_bright(3);   % radius of circle
            if bPlot; subplot(1,3,1); rectangle('Position', [c_x-r c_z-r 2*r 2*r]*1e3, 'Curvature', [1 1], 'EdgeColor', [1 1 1]); end; % plot
            
            %-- inside of circle
            [mZ, mX] = ndgrid(aZaxis,aXaxis);
            aIndex_w = find((mX-c_x).^2 + (mZ-c_z).^2 < r^2); % find indices within circle
            aBright = mData(aIndex_w);     
            
        elseif strcmp(sUnit,'pixel')
            error('CalcCNR:: position for circle should be in mm (because of radius) ');  
        else
            error('CalcCNR:: The string sUnit should be mm or pixel');
        end   
        
        if bPlot     
            subplot(1,3,2);
            mPlot_bright = zeros(size(mData,1),size(mData,2));
            mPlot_bright(aIndex_w) = mLogOut_80dB(aIndex_w); % only dark region
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mPlot_bright);    
        end
    
    elseif strcmp(sShape_bright, 'rectwh') % rectangular with a hole
        if strcmp(sUnit,'mm')        
            x = aPos_bright(1); % x pos of center of rect
            z = aPos_bright(2); % z pos of center of rect
            w = aPos_bright(3); % width
            h = aPos_bright(4); % hight 
            c_x = aPos_bright(5); % x pos of center of hole
            c_z = aPos_bright(6); % z pos of center of hole
            r = aPos_bright(7);   % radius of hole
            if bPlot; subplot(1,3,1); rectangle('Position', [x-w/2, z-h/2, w, h]*1e3, 'EdgeColor', [1 1 0]); end; % plot
            if bPlot; subplot(1,3,1); rectangle('Position', [c_x-r c_z-r 2*r 2*r]*1e3, 'Curvature', [1 1], 'EdgeColor', [1 1 0]); end; % plot hole
            
            %-- rectangular region
            % convert 'mm' to 'pixel'            
            xs_p = round(interp1(aXaxis,1:size(mData,2),x-w/2)); % x pixel of upper-left corner
            xe_p = round(interp1(aXaxis,1:size(mData,2),x+w/2)); % x pixel of lower-right corner
            zs_p = round(interp1(aZaxis,1:size(mData,1),z-h/2)); % z pixel of upper-left corner
            ze_p = round(interp1(aZaxis,1:size(mData,1),z+h/2)); % z pixel of lower-right corner            
            mLogic = zeros(size(mData,1),size(mData,2));
            mLogic(zs_p:ze_p,xs_p:xe_p) = 1;            
            
            %-- hole
            [mZ, mX] = ndgrid(aZaxis,aXaxis);
            aIndex_hole = find((mX-c_x).^2 + (mZ-c_z).^2 < r^2);
            mLogic(aIndex_hole) = 0; % leave out inside of hole
            
            %-- rect with hole
            aIndex_recwh = find(mLogic == 1);
            aBright = mData(aIndex_recwh);  
            
%             [mZ, mX] = ndgrid(aZaxis,aXaxis);
%             % windowing rect region
%             mZ_w = mZ(zs_p:ze_p,xs_p:xe_p); % 2-D windowed Z grid
%             mX_w = mX(zs_p:ze_p,xs_p:xe_p); % 2-D windowed X grid
%             mImage_w = mData(zs_p:ze_p,xs_p:xe_p); % 2-D windowed image
            
%             aIndex_w = find((mX_w-c_x).^2 + (mZ_w-c_z).^2 > r^2); % find indices outside of circle
%             aBright = mImage_w(aIndex_w);     
            
        elseif strcmp(sUnit,'pixel')
            error('position for circle should be in mm (because of radius) ');  
        else
            error('The string sUnit should be mm or pixel');
        end 
        
        if bPlot     
            subplot(1,3,3);
            mPlot_bright = zeros(size(mData,1),size(mData,2));
            mPlot_bright(aIndex_recwh) = mLogOut_80dB(aIndex_recwh); % only dark region
            imagesc(stG.aX*1e3,stG.aZ0*1e3,mPlot_bright);    
        end
    else
        error('sShape_bright should be circ, rect, rectwh');
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
            title(['dark region']);
        subplot(1,3,3);
            set(gca,'FontName','Times New Roman', 'FontSize', 12);  
            axis equal; axis([stG.aX(1) stG.aX(end) stG.aZ0(1) stG.aZ0(end)]*1e3); %axis([-30 30 -30 30 0 70]);
            xlabel('{\itx} (mm)'); ylabel('{\itz} (mm)'); colormap(gray);
            title(['bright region']);
        drawnow;  
    end    
    
    %% Contrast
    
    nE_dark   = sum(aDark.^2);
    nE_bright = sum(aBright.^2);
    nCon = 10*log10(nE_dark/nE_bright);
    stCon.sDataName     = sDataName;
    stCon.sShape_dark   = sShape_dark;    % shape of dark region
    stCon.aPos_dark     = aPos_dark;      % position and size of dark region
    stCon.sShape_bright = sShape_bright;  % shape of bright region
    stCon.aPos_bright   = aPos_bright;    % position and size of bright region
    stCon.nE_dark     = nE_dark;       % measured energy of dark region
    stCon.nE_bright   = nE_bright;     % measured energy of bright region
    stCon.nCon          = nCon;
    
    %% CNR 
%     
%     B_mean = mean(aBright)
%     D_mean = mean(aDark)
%     B_std = std(aBright)
%     D_std = std(aDark)    
%     
%     nCNR = (B_mean-D_mean)/sqrt(B_std^2 + D_std^2)


%     nCNR = sum(aBright)/sum(aDark);

    


end
