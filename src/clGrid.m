% Sua Bae
% 2017-12-19
%
classdef clGrid
   
    properties (SetAccess = public)
        % common attributes
        nXdim    % x dimension
        nZdim    % z dimension
        
        dx       % x interval of grid
        dz       % z interval of grid
        
        nXos     % x offset
        nZos     % z offset
        
        aX       % x axis
        aZ       % z axis
        
    end
    
    %-- Constructor
    methods (Access = public)
        
        function h = clGrid(nXdim, nZdim, dx, dz, nXos, nZos)           
%             h.creation_date=sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);

            h.nXdim = nXdim; h.nZdim = nZdim; 
            h.dx    = dx;    h.dz    = dz;  
            h.nXos  = nXos;  h.nZos  = nZos; 
            
            h.aX  = (-(nXdim-1)/2:1:(nXdim-1)/2)' * dx + nXos;
            h.aZ  = (-(nZdim-1)/2:1:(nZdim-1)/2)' * dz + nZos;
            
        end
    end
    
    
    
    methods (Access = public)   
        
        function xidx = idx2xidx(h, idx)
            xidx = mod(idx-1, h.nXdim)+1;
        end
        
        function zidx = idx2zidx(h, idx)
            zidx = floor((idx-1)/(h.nXdim))+1;
        end
        
        function idx = xzidx2idx(h, xidx, zidx)
            idx = xidx + (zidx-1)*h.nXdim;
        end
        
        function [x,z] = rt2xz(h, r, t)  
            x = r.*sind(t);
            z = r.*cosd(t);
        end
        
        function [r,t] = xz2rt(h, x, z)  
            r = sqrt(x.*x + z.*z);
            t = tand(x./z);
        end
    end
    
    
    
    
end