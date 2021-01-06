% Sua Bae
% 2017-11-30
%
function stScat = ScatGen(phantom_idx, stTrans)

    % 1. points center
    if phantom_idx == 1
        nNumScat = 11;
        stScat.mScatXYZPos = [zeros(nNumScat,1), zeros(nNumScat,1), (10e-3:10e-3:110e-3)'];
        stScat.aScatMag    = 1e26*ones(nNumScat,1);
        
    elseif phantom_idx == 2
        nNumScat = 1;
        stScat.mScatXYZPos = [0, 0, 20e-3'];
        stScat.aScatMag    = 1e26;
        
    elseif phantom_idx == 3
        nNumScat = 1;
        stScat.mScatXYZPos = [0, 0, 40e-3'];
        stScat.aScatMag    = 1e26;
        
    elseif phantom_idx == 4
        nNumScat = 1;
        stScat.mScatXYZPos = [0, 0, 60e-3'];
        stScat.aScatMag    = 1e26;
        
    elseif phantom_idx == 5
        nNumScat = 1;
        stScat.mScatXYZPos = [0, 0, 80e-3'];
        stScat.aScatMag    = 1e26;        
    
    elseif phantom_idx == 6
        nNumScat = 1;
        stScat.mScatXYZPos = [0, 0, 100e-3'];
        stScat.aScatMag    = 1e26;
        
    elseif phantom_idx == 7
        nNumScat = 11;
        aPos_r = (10e-3:10e-3:110e-3)' + stTrans.nRadius;
        aPos_x = aPos_r*sind(30);
        aPos_z = aPos_r*cosd(30) - stTrans.nRadius;
        stScat.mScatXYZPos = [aPos_x(:), zeros(nNumScat,1), aPos_z(:)];
        stScat.aScatMag    = 1e26*ones(nNumScat,1);
        
    elseif phantom_idx == 8
        nNumScat = 1;
        nPos_r = 20e-3 + stTrans.nRadius;
        nPos_x = nPos_r*sind(30);
        nPos_z = nPos_r*cosd(30) - stTrans.nRadius;
        stScat.mScatXYZPos = [nPos_x, 0, nPos_z];
        stScat.aScatMag    = 1e26;
        
    elseif phantom_idx == 9
        nNumScat = 1;
        nPos_r = 40e-3 + stTrans.nRadius;
        nPos_x = nPos_r*sind(30);
        nPos_z = nPos_r*cosd(30) - stTrans.nRadius;
        stScat.mScatXYZPos = [nPos_x, 0, nPos_z];
        stScat.aScatMag    = 1e26;
        
    elseif phantom_idx == 10
        nNumScat = 1;
        nPos_r = 60e-3 + stTrans.nRadius;
        nPos_x = nPos_r*sind(30);
        nPos_z = nPos_r*cosd(30) - stTrans.nRadius;
        stScat.mScatXYZPos = [nPos_x, 0, nPos_z];
        stScat.aScatMag    = 1e26;
        
    elseif phantom_idx == 11
        nNumScat = 1;
        nPos_r = 80e-3 + stTrans.nRadius;
        nPos_x = nPos_r*sind(30);
        nPos_z = nPos_r*cosd(30) - stTrans.nRadius;
        stScat.mScatXYZPos = [nPos_x, 0, nPos_z];
        stScat.aScatMag    = 1e26;
        
    elseif phantom_idx == 12
        nNumScat = 1;
        nPos_r = 100e-3 + stTrans.nRadius;
        nPos_x = nPos_r*sind(30);
        nPos_z = nPos_r*cosd(30) - stTrans.nRadius;
        stScat.mScatXYZPos = [nPos_x, 0, nPos_z];
        stScat.aScatMag    = 1e26;
    end
            
    figure; scatter(stScat.mScatXYZPos(:,1)*1e3,stScat.mScatXYZPos(:,3)*1e3,'filled');  axis equal;
    hold on; scatter(stTrans.mElePos(:,1)*1e3, stTrans.mElePos(:,3)*1e3 - stTrans.nRadius*1e3, 5*ones(stTrans.nNumEle,1),'filled');  axis equal;
    axis([-80 80 -10 120]); set(gca,'YDir','reverse'); grid on;
    
end