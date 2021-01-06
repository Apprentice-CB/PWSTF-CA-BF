% Sua Bae
% 2017-12-18
%
classdef clDscGrid_convex
   
    properties (SetAccess = public)
        % common attributes
        nXdim    % x dimension
        nZdim    % z dimension
        
        dx       % x interval of grid
        dz       % z interval of grid
        
        nMaxTheta % maximum angle of transducer element (deg)
        nRadius  % radius of transducer element
        
        aX       % 1D R axis
        aZ       % 1D Z axis  (0 at the origin of convex array)
        aZ0      % 1D Z0 axis (0 at transducer surface)
        
    end
    
    %-- Constructor
    methods (Access = public)
        
        function h = clDscGrid_convex(nXdim, nZdim, dx, dz, nMaxTheta, nRadius)           
%             h.creation_date=sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);

            h.nXdim = nXdim; h.nZdim = nZdim; 
            h.dx    = dx;    h.dz    = dz;  
            h.nMaxTheta  = nMaxTheta;  h.nRadius  = nRadius; 
            
            h.aX  = (-(nXdim-1)/2:1:(nXdim-1)/2)' * dx;
            h.aZ  = (0:1:(nZdim-1))' * dz + nRadius*cosd(nMaxTheta);
            h.aZ0 = (0:1:(nZdim-1))' * dz - nRadius*(1-cosd(nMaxTheta));
            
        end
    end
    
    
    
    methods (Access = public)        
        function [x,z] = rt2xz(h, r, t)  
            x = r.*sind(t);
            z = r.*cosd(t);
        end
    end

    methods (Access = public)        
        function [r,t] = xz2rt(h, x, z)  
            r = sqrt(x.*x + z.*z);
            t = tand(x./z);
        end
    end
    
    
    
    
end