% Sua Bae
%  2017-11-13
%
% Xdata : input
% Ydata : output
% nDR : Dynamic Range [dB]
% nYmin : min of output (default: 0)
% nYmax : max of output (default: 255)   
% nXmax : max of input (default: max(Xdata(:))
function  Ydata = LogCompression(Xdata, nDR, nYmin, nYmax, nXmax)

    %% ymin, ymax 
    if ~exist('nYmin','var')
        ymin = -nDR;
    else
        ymin = nYmin;
    end
    if ~exist('nYmax','var')
        ymax = 0;
    else
        ymax = nYmax;
    end
    %% xmax 
    if ~exist('nXmax','var')
        xmax = max(Xdata(:));
    else
        xmax = nXmax;
    end
    
    %% xmin 
    xmin = xmax / ( 10^(nDR/20) ) ; 

    %% y = A*log10(x/xmax)+B 
    A = (ymax-ymin)/(log10(xmax/xmin));
    B = ymax;
    
    %% Log compression
    Xdata = min(Xdata,xmax) ; 
    Xdata = max(Xdata,xmin) ;      
    Ydata = A*log10(Xdata/xmax)+B;
   
    %%
    Ydata = min(Ydata,ymax) ; 
    Ydata = max(Ydata,ymin) ;  
    
end