% Sua Bae
% 2017-12-18
%
classdef clDir
   
    properties (SetAccess = public)
        % common attributes
        phantom_idx     % phantom index
        acqType_idx    % acquisition type index
        dataNumber_idx   % data number index 
        
        sPhantom        % phantom
        sAcqType        % acquisition type
        
        sPhtmDir        % phantom directory
        sAcqDir         % acquisition directory
        sDataDir        % data directory
    end
    
    %-- Constructor
    methods (Access = public)
        
        function h = clDir(phantom_idx, acqType_idx, dataNumber_idx, sTopDir)           
%             h.creation_date=sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
            h.phantom_idx       = phantom_idx;
            h.acqType_idx       = acqType_idx;
            h.dataNumber_idx    = dataNumber_idx;
            
            %%--Phantom
            if     phantom_idx == 1;  h.sPhantom = 'p';
            elseif phantom_idx == 2;  h.sPhantom = 'p';
            elseif phantom_idx == 3;  h.sPhantom = 'p';
            elseif phantom_idx == 4;  h.sPhantom = 'p';
            elseif phantom_idx == 5;  h.sPhantom = 'p';
            elseif phantom_idx == 6;  h.sPhantom = 'p';
            elseif phantom_idx == 7;  h.sPhantom = 'p';
            elseif phantom_idx == 8;  h.sPhantom = 'p';
            elseif phantom_idx == 9;  h.sPhantom = 'p';
            elseif phantom_idx == 10; h.sPhantom = 'p';
            elseif phantom_idx == 11; h.sPhantom = 'p';
            elseif phantom_idx == 12; h.sPhantom = 'p';
            elseif phantom_idx == 100; h.sPhantom = 'noise without tx';
            elseif phantom_idx == 101; h.sPhantom = 'noise with tx';
            else   error('undefined');
            end
            h.sPhtmDir = [sTopDir num2str(phantom_idx) '. ' h.sPhantom '\'];
            
            %%--Acquisition Type
            if     acqType_idx == 1; h.sAcqType = 'PWSTF';
            elseif acqType_idx == 2; h.sAcqType = 'DWSTF centerRes';
            elseif acqType_idx == 3; h.sAcqType = 'CF Fnum2p4';
            elseif acqType_idx == 4; h.sAcqType = 'CF Fnum5';
            elseif acqType_idx == 5; h.sAcqType = 'PWSTF EA';
            else   error('undefined');
            end
            h.sAcqDir = [h.sPhtmDir num2str(acqType_idx) '. ' h.sAcqType '\'];
            
            %%--Data number
            h.sDataDir = [h.sAcqDir 'Data_' num2str(h.dataNumber_idx) '\'];
            
        end
    end
    
    
    methods (Access = public)
        function CreateFolders(h)

            % %%%%% Phantom Select %%%%% %
            if ~exist(h.sPhtmDir,'dir')
                display(['creating folder:: ' h.sPhtmDir]);
                mkdir(h.sPhtmDir); 
            end

            % %%%%% Acquisition Type Select %%%%% %
            if ~exist(h.sAcqDir,'dir')
                display(['creating folder:: ' h.sAcqDir]);
                mkdir(h.sAcqDir); 
            end

            % %%%%% Data Number Select %%%%% %            
            if ~exist(h.sDataDir,'dir')
                display(['creating folder:: ' h.sDataDir]);
                mkdir(h.sDataDir); 
            end

            % %%%%% Final Data Directory %%%%% %
            display(['creating RcvData folder at ' h.sDataDir]);

            % %%%%% Create folders %%%% %    
            sRFDir = [h.sDataDir 'RcvData\'];
            if(~isdir(sRFDir)); 
                mkdir(sRFDir); 
            end  

        end
    end

    
    
end