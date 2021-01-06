% 2013-03-14
% Sua Bae
%
% fInphase : Inphase components
% fQuadrature : Quadrature components
%
% fData : RF signal
% nFs : Sampling frequency
% nDemodFreq : Demodulation frequency
% nDemodLPFFreq : LPF cutoff frequency
% nDeci : Decimation ratio
%
function [mOutData] = QDM(mData, nFs, nDemodFreq, nDemodLPFFreq, nDemodLPFTap, nDeci)

    %% Demodulation    
    n =0:size(mData,1)-1; 
    cos_n = cos(2*pi*nDemodFreq/nFs*n);
    sin_n = sin(2*pi*nDemodFreq/nFs*n);
    
    mInphase_temp = zeros(size(mData));
    mQuadrature_temp = zeros(size(mData));
    
    for idx = 1:size(mData,2)
        mInphase_temp (:,idx) = mData(:,idx).* cos_n';
        mQuadrature_temp (:,idx) = mData(:,idx).* -sin_n';
    end
    
    %% LPF    
    aLPFCoefs = fir1(nDemodLPFTap-1, nDemodLPFFreq/nFs*2, 'low');
%     aLPFCoefs = fir1(nDemodLPFTap-1, nDemodLPFFreq/nFs*2/4, 'low');
    mInphase = convn(mInphase_temp, aLPFCoefs', 'same');
    mQuadrature = convn(mQuadrature_temp, aLPFCoefs', 'same');
    
    %% Decimation
%     fOutData = fInphase + 1i*fQuadrature;
%     fOutData = fInphase(1:nDeci:end, :) - 1i*fQuadrature(1:nDeci:end, :); 
    mOutData = mInphase(1:nDeci:end, :) + 1i*mQuadrature(1:nDeci:end, :);  % demod시 sin에 -적용!
    
    
%     aAx = linspace(-30,30,2048);
%     figure;
%     subplot(3,1,1);
%         plot(aAx,FFTDB(mData(:,128))); title('data'); grid on; ylim([-150 0]);
%     subplot(3,1,2);
%         plot(aAx,FFTDB(mInphase_temp(:,128))); title('*cos'); grid on; ylim([-150 0]);
% %     subplot(3,1,3);
% %         plot(aAx,FFTDB(mInphase(:,128))); title('after LPF'); grid on; ylim([-150 0]);
% %         hold on; plot(aAx,FFTDB(aLPFCoefs));
%     subplot(3,1,3);
%         plot(aAx,FFTDB(mInphase(1:nDeci:end, 128))); title('after LPF & decimation'); grid on; ylim([-150 0]);
%         hold on; plot(aAx,FFTDB(aLPFCoefs(1:nDeci:end)));      
        
%     figure;
%     subplot(3,1,1);
%         plot(aAx,FFTDB(mData(:,128)));  title('data'); grid on; ylim([-150 0]);
%     subplot(3,1,2);
%         plot(aAx,FFTDB(mQuadrature_temp(:,128))); title('*sin'); grid on; ylim([-150 0]);
% %     subplot(3,1,3);
% %         plot(aAx,FFTDB(mQuadrature(:,128))); title('after LPF'); grid on; ylim([-150 0]);
% %         hold on; plot(aAx,FFTDB(aLPFCoefs));
%     subplot(3,1,3);
%         plot(aAx,FFTDB(mQuadrature(1:nDeci:end, 128))); title('after LPF & decimation'); grid on; ylim([-150 0]);
%         hold on; plot(aAx,FFTDB(aLPFCoefs(1:nDeci:end)));     
    

%     aAx = linspace(-30,30,2048);
%     figure('Position',[100,100,400,700]);
%     subplot(4,1,1);
%         nDeci = 1;
%         plot(aAx,FFTDB(mInphase(1:nDeci:end, 128))); title(['after LPF & decimation (x' num2str(nDeci) ')']); grid on; ylim([-150 0]);
%         hold on; plot(aAx,FFTDB(aLPFCoefs(1:nDeci:end)));  
%     subplot(4,1,2);
%         nDeci = 4;
%         plot(aAx,FFTDB(mInphase(1:nDeci:end, 128))); title(['after LPF & decimation (x' num2str(nDeci) ')']); grid on; ylim([-150 0]);
%         hold on; plot(aAx,FFTDB(aLPFCoefs(1:nDeci:end)));  
%     subplot(4,1,3);
%         nDeci = 8;
%         plot(aAx,FFTDB(mInphase(1:nDeci:end, 128))); title(['after LPF & decimation (x' num2str(nDeci) ')']); grid on; ylim([-150 0]);
%         hold on; plot(aAx,FFTDB(aLPFCoefs(1:nDeci:end)));  
%     subplot(4,1,4);
%         nDeci = 16;
%         plot(aAx,FFTDB(mInphase(1:nDeci:end, 128))); title(['after LPF & decimation (x' num2str(nDeci) ')']); grid on; ylim([-150 0]);
%         hold on; plot(aAx,FFTDB(aLPFCoefs(1:nDeci:end)));  
        
end