% Sua Bae
% 2017-12-20
%
% mBfData_rt :data on BF grid
% stG_rt     : BF grid (rt grid)
% stG        : DSC grid (xz grid)
function mBfData = DSC_rt2xz(mBfData_rt, stG_rt, stG)

    % (x,z) of Cartesian grid (xz grid)
    [mZ_cts,mX_cts] = ndgrid(stG.aZ,stG.aX);

    % (r,t) of Cartesian grid
    [aR_cts,aT_cts] = xz2rt(stG_rt, mX_cts(:), mZ_cts(:));

    % ndgrid of RT grid
    [mR_rt, mT_rt] = ndgrid(stG_rt.aR,stG_rt.aT);

    % interpolate data to Cartesian grid
    aBfData_cts = interpn(mR_rt, mT_rt, mBfData_rt, aR_cts, aT_cts);
    aBfData_cts(isnan(aBfData_cts)) = 0;

    % reshape and export!
    mBfData = reshape(aBfData_cts,stG.nZdim,stG.nXdim); %  (xz grid)
    
end