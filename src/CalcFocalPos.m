% Sua Bae
% 2017-12-21
%
%  stTrans      : transducer structure (.nPitch_x, .nNumEle_x)
%  nFocalDepth  : focal depth (from transducer)
%
%  mFocalPos: dim = [ num of foci x 5(x,z,z0,r,tht)]
%
function [mFocalPos] = CalcFocalPos(stTrans, nNumTx, nFocalDepth, bPlot)

    da = stTrans.nPitchAngle_x*stTrans.nNumEle_x/nNumTx; % angle interval between two scanlines
    
    aTheta = ( -(nNumTx-1)/2:1:(nNumTx-1)/2 )' * da;
    aRadius = (nFocalDepth + stTrans.nRadius)*ones(nNumTx,1);
    
    aX = aRadius.*sind(aTheta);
    aZ = aRadius.*cosd(aTheta);
    aZ0 = aZ - stTrans.nRadius;
    
    if(bPlot)
        figure;
        scatter(aX*1e3,aZ0*1e3); 
        hold on; 
        line([0,200*sind(-stTrans.nMaxTheta)], [0,200*cosd(-stTrans.nMaxTheta)] - stTrans.nRadius*1e3);
        line([0,200*sind(stTrans.nMaxTheta)],  [0,200*cosd(stTrans.nMaxTheta)] - stTrans.nRadius*1e3);
        set(gca,'YDir','reverse');
        axis([-100 100 0 120]); box on; grid on; title(['N = ' num2str(nNumTx) ' at d = ' num2str(nFocalDepth*1e3) ' mm']);
    end
    
    
    mFocalPos = [aX(:), aZ(:), aZ0(:), aRadius(:), aTheta(:)];
    
end