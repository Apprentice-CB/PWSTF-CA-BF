% Sua Bae
% 2017-12-18
%
classdef clBFGrid_rt
   
    properties (SetAccess = public)
        % common attributes
        nRdim    % number of BF points per each scanline (R dimension)
        nTdim    % number of BF scanlines (Theta dimension)
        
        dr       % r interval of grid
        dt       % tht interval of grid
        
        nRos     % r offset of grid  % = stTrans.nRadius
        nTos     % tht offset of grid
        
        aR       % 1D R axis
        aR0      % 1D R0 axis : 0 at transducer surface (= aR - nRos)
        aT       % 1D theta axis
        
    end
    
    %-- Constructor
    methods (Access = public)
        
        function h = clBFGrid_rt(nRdim, nTdim, dr, dt, nRos, nTos)           
%             h.creation_date=sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);

            h.nRdim = nRdim; h.nTdim = nTdim; 
            h.dr    = dr;    h.dt    = dt;  
            h.nRos  = nRos;  h.nTos  = nTos; 
            
            h.aR  = (0:1:(nRdim-1))' * dr + nRos;
            h.aR0 = (0:1:(nRdim-1))' * dr ;
            h.aT  = (-(nTdim-1)/2:1:(nTdim-1)/2)' * dt + nTos;
            
        end
    end
    
    
    
    methods (Access = public)        
        function [x,z] = rt2xz(h, r, t)  
            x = r.*sind(t);
            z = r.*cosd(t);
        end
    end
    
    methods (Access = public)        
        function [x,z] = r0t2xz(h, r0, t)  
            x = (r0+h.nRos).*sind(t);
            z = (r0+h.nRos).*cosd(t);
        end
    end
    
    methods (Access = public)        
        function [r,t] = xz2rt(h, x, z)  
            r = sqrt(x.*x + z.*z);
            t = atand(x./z);
        end
    end
    
    
    methods (Access = public)        
        function [r0,t] = xz2r0t(h, x, z)  
            r0 = sqrt(x.*x + z.*z) - h.nRos;
            t = atand(x./z);
        end
    end
    
    
    
    
end