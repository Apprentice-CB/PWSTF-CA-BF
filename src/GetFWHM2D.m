% Sua Bae
% 2017-11-21
%
% mData: log compressed data
function [nFWHM, nLEdge_x, nLEdge_z, nREdge_x, nREdge_z] = GetFWHM2D(aX, aZ, mData)

        nShow_dB = -6;
        
        % contour
        figure(101);
        [C, h]=contour(aX, aZ, mData', [nShow_dB nShow_dB], 'edgecolor','k');  
        
        
        % find left and right edge position
        C(C==nShow_dB)=1e100; % Some elements at first row are nShow_dB.
        [nLEdge_x, llidx] = min(C(1,2:end));
        nLEdge_z = C(2,llidx+1);

        C(C==1e100)=-1e10; % Some elements at first row are nShow_dB.
        [nREdge_x, rridx] = max(C(1,2:end));
        nREdge_z = C(2,rridx+1);

        
        nFWHM = sqrt( (nLEdge_x-nREdge_x).^2 + (nLEdge_z-nREdge_z).^2 );
        close(101);
end