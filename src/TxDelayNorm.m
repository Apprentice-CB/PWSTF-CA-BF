% 2017-12-22
%
% Sua Bae
%
function  [aDelay_norm, nDelayOffset] = TxDelayNorm(aDelay, aApod, sStartAcqTime)
    
    aDelay_norm = zeros(1,numel(aApod));
    if strcmp(sStartAcqTime,'first firing time')
        aDelay_apod = aDelay(logical(aApod));
        nDelayOffset = min(aDelay_apod);
        aDelay_apod_norm = aDelay_apod - min(aDelay_apod);
        aDelay_norm(logical(aApod)) = aDelay_apod_norm;
    else
        error('undefined');
    end
    
            
end