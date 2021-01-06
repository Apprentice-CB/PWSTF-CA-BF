% 2017-12-21
%
% Sua Bae
%
%
% <Linear>
%   - CF: Tx delay at each element assuming the focused wave locates at x-axis (t=0 when the last element fired)
% <Convex>
%   Tx delay at each element assuming the wave (focused, DW, PW) locates at the origin when t=0
function  [aTxDelay] = TxDelayGen(stTrans, stTx, txidx, nSoundSpeed)

    switch stTx.sType
        case 'CF'        
            aFocalPos = stTx.mFocalPos(txidx,:); % (x,z,z0,r,t) [m]
            
            aFocus2Ele = sqrt((stTrans.mElePos(:,1) - aFocalPos(1)).^2 + (stTrans.mElePos(:,3) - aFocalPos(2)).^2); % distance from element to focus [meter]         
            aTxDelay = (aFocalPos(4) - aFocus2Ele)/nSoundSpeed;% full aperture delay [sec] <- [meter]
            
        case 'PW'

            switch stTrans.sType
            case 'Linear'
                error('undefined');
            case 'Convex'
                nPwAngle = stTx.aPwAngle_deg(txidx);
                aPw2Ele = stTrans.nRadius.*cosd(stTrans.mElePos(:,4) - nPwAngle); % Distance from the PW at origin to element [meter] % aEleToPW: dim = 1x(nEleNum)
                aTxDelay =  aPw2Ele/nSoundSpeed; % full aperture delay [sec] <- [meter]
            otherwise
                error('undefined');
            end

            
        case 'DW'

            switch stTrans.sType
            case 'Linear'
                error('undefined');
            case 'Convex'
                % only when virtual source is behind of transducer
                aVSPos = stTx.mVSPos(txidx,:);
                aVs2Ele = sqrt((stTrans.mElePos(:,1)-aVSPos(1)).^2 + (stTrans.mElePos(:,3)-aVSPos(3)).^2); % Distance from the PW at origin to element [meter] % aEleToPW: dim = 1x(nEleNum)
                aTxDelay =  aVs2Ele/nSoundSpeed; % full aperture delay [sec] <- [meter]
            otherwise
                error('undefined');
            end


    end

            
end