classdef aoSystem < handle

    properties (SetObservable=true)        
        
        % --------------------- PARAMETERS STRUCTURE --------------------%
        parms;
        % ------------------------- COMPONENTS -------------------------%        
        sci;
        tel;
        atm;                
        lGs;
        nGs;
        tGs;
        wfs;
        lowfs;
        truthwfs;
        dm;        
        ttdm;
        cam;        
        % ------------------------- CALIBRATION -------------------------%            
        wfs2dm;
        truth2dm;
        matrices;
        ncpa;
        % --------------------------- AO LOOP ---------------------------%
        rtc;
        loopStatus;
        % ------------------------- TELEMETRY ---------------------------%
        loopData;
        psf;
        % --------------------------- FLAGs -----------------------------%
        flagCalibration = false;
        flagSimulation  = false;
    end
    
    
    
    
    methods
        
        function obj = aoSystem(parms,varargin)
            
            %Check inputs
            inputs = inputParser;
            inputs.addRequired('parms', @isstruct);            
            inputs.addParameter('runSimulation', false,@islogical);
            inputs.parse(parms,varargin{:});
            
            obj.parms = parms;
            flagSimu  = inputs.Results.runSimulation;
            
            % 1/ Check inputs            
            parms = obj.checkParms(parms);
            
            % 2/ Instantiation
            obj = obj.initComponents(parms);
                                    
            if flagSimu   
                % 3/ Calibration
                obj = obj.calibrateSystem();
            
                % 4/ Simulation
                obj = obj.runSimulation();                               
            end
        end
                                
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        %                       INSTANTIATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = initComponents(obj,parms)
                                                        
            
            %1/ Science
            obj.sci = obj.initSource(parms.sci,parms.atm);
            
            %2/ Telescope
            if ~isfield(parms,'seg')
                parms.seg = [];            
            end            
            obj.tel = obj.initTelescope(parms.tel,parms.seg);
            %3/ Calibration
            telCalib       = utilities.duplicateTelescope(obj.tel);
            telCalib.pupil = double(logical(real(telCalib.pupil)));                        
            srcCalib = source();
                       
            %4/ Atmosphere
            obj.atm = obj.initAtmosphere(parms.atm);
            
            %5/ Guide stars                       
            obj.nGs = obj.initSource(parms.nGs,parms.atm); % NGSs constellation
            if isfield(parms,'lGs')
                 parms.lGs.apertureSize     = parms.tel.D/parms.wfs.nLenslet;
                 parms.lGs.apertureDistance = parms.tel.D;
                obj.lGs = obj.initSource(parms.lGs,parms.atm); % LGSs constellation
            end
            if isfield(parms,'tGs')
                obj.tGs = obj.initSource(parms.tGs,parms.atm); % NGS for the truth sensor
            end
            
            %6/ WFS
            if isfield(parms,'lGs') && ~isempty(obj.lGs)
                srcCalib.photometry = obj.lGs(1).photometry;
                obj.wfs             = obj.initWfs(parms.wfs,telCalib,srcCalib );        
                srcCalib.photometry = obj.nGs(1).photometry;
                obj.lowfs           = obj.initWfs(parms.lowfs,telCalib,srcCalib );     
            else
                srcCalib.photometry = obj.nGs(1).photometry;
                obj.wfs             = obj.initWfs(parms.wfs,telCalib,srcCalib );
            end
            
            if isfield(parms,'lGs') && ~isempty(obj.tGs)
                srcCalib.photometry = obj.tGs.photometry;
                obj.truthwfs   = obj.initWfs(parms.truwfs,telCalib,srcCalib );
            end
            
            %7/ DM
            obj.dm  = obj.initDm(parms.dm,obj.tel.resolution,obj.wfs.validActuator);
            %obj.ttdm= obj.initDm(parms.ttdm,obj.tel.resolution,true(2));
            
            %8/ IMAGERS
            srcCalib.photometry = obj.sci(1).photometry;
            obj.cam = obj.initImager(parms.cam,telCalib,srcCalib);                        
        end
                       
        function obj = initLoop(obj,parms)
                                    
            %1/ WFS configuration            
            obj.wfs.camera.exposureTime = parms.wfs.exposureTime;
            obj.wfs.camera.clockRate    = parms.wfs.frameRate;
            obj.wfs.camera.photonNoise  = parms.wfs.photonNoise;
            obj.wfs.camera.readOutNoise = parms.wfs.ron;
            obj.wfs.framePixelThreshold = parms.wfs.pixelThreshold;
            obj.wfs.lenslets.throughput = parms.wfs.throughput;
            
            if ~isempty(obj.lowfs)               
                obj.lowfs.camera.exposureTime = parms.lowfs.exposureTime;
                obj.lowfs.camera.clockRate    = parms.lowfs.frameRate;
                obj.lowfs.camera.photonNoise  = parms.lowfs.photonNoise;
                obj.lowfs.camera.readOutNoise = parms.lowfs.ron;
                obj.lowfs.lenslets.throughput = parms.lowfs.throughput;
            end
            
            if ~isempty(obj.truthwfs)                               
                obj.truthwfs.camera.exposureTime = parms.truwfs.exposureTime;
                obj.truthwfs.camera.clockRate    = parms.truthwfs.frameRate;
                obj.truthwfs.camera.photonNoise  = parms.truwfs.photonNoise;
                obj.truthwfs.camera.readOutNoise = parms.truwfs.ron;
                obj.truthwfs.framePixelThreshold = parms.truwfs.pixelThreshold;
                obj.truthwfs.lenslets.throughput = parms.truwfs.throughput;
            end
            
            
            %2/ Main loop status           
            obj.loopStatus               = parms.loopStatus;     
            obj.loopStatus.ho.frameRate  = parms.wfs.frameRate;
            if ~isempty(obj.lowfs)
                obj.loopStatus.lo.frameRate  = parms.lowfs.frameRate;
            end
            if ~isempty(obj.truthwfs)
                obj.loopStatus.tru.frameRate = parms.truwfs.frameRate;
            end

            obj.loopStatus.nIteration = parms.cam.exposureTime + parms.cam.startDelay;
            obj.loopStatus.startDelay = parms.cam.startDelay;
            
            %3/ Grab external information
            if isfield(parms,'ncpa')
                obj.ncpa.map = parms.ncpa;
            else
                obj.ncpa.map = [];
            end
                                    
            %4/ Reset
            flush(obj.cam);
            obj.tel = obj.tel + obj.atm;
                                                          
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        %                         CALIBRATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = calibrateSystem(obj)
            
            
            %1/ Calibrate the interaction matrix + reconstructor     
            telCalib       = utilities.duplicateTelescope(obj.tel);
            telCalib.pupil = double(logical(real(telCalib.pupil)));                        
            srcCalib       = source();
            
            % MAIN WFS                                    
            if isempty(obj.lGs)
                srcCalib.photometry = obj.nGs(1).photometry;
            else
                srcCalib.photometry = obj.lGs(1).photometry;
            end
            obj.wfs2dm = obj.calibrateInteractionMatrix(obj.parms.loopStatus.ho,telCalib,obj.dm,obj.wfs,srcCalib);                                             
            
            % TRUTH SENSOR
            if ~isempty(obj.truthwfs)                
                srcCalib.photometry = obj.tGs.photometry;
                obj.truth2dm = obj.calibrateInteractionMatrix(obj.parms.loopStatus.tru,telCalib,obj.dm,obj.truthwfs,srcCalib);              
            end

            
            %2/ Calibrate the dm filter
            dmLow                 = obj.initDm(obj.parms.dm,obj.dm.nActuator,obj.dm.validActuator);
            FmapLow               = 2*dmLow.modes.modes(obj.dm.validActuator,:);
            dmFilter              = pinv(full(FmapLow),1e-1);
            obj.matrices.dmFilter = dmFilter;
                
            %3/ Calibrate tip-tilt filters    
            
            % Tip-tilt slopes filters
            t  = ones(obj.wfs.nValidLenslet,1);
            tt = [t 0*t; 0*t,t];
            obj.matrices.SlopeTTRem = eye(2*obj.wfs.nValidLenslet) - tt*pinv(tt);
            obj.matrices.slopes2Tilt = pinv(tt);
                                                               
            % Tip-tilt commands filters
            zern = zernike(2:3,'resolution',obj.dm.nActuator,'pupil',obj.dm.validActuator);
            TT   = dmFilter*zern.modes(obj.dm.validActuator,:);
            obj.matrices.DMTTRem       = eye(obj.dm.nValidActuator) - TT*pinv(TT);
            obj.matrices.tilt2Commands = TT;
            obj.matrices.commands2Tilt = pinv(TT);
                
            %4/ Calibrate tomographic reconstructor  
            if numel(obj.nGs) > 1 || numel(obj.lGs) > 1    
                if ~isempty(obj.lowfs)
                    ast      = obj.lGs;
                    mmseStar = source('wavelength',obj.lGs(1).photometry);
                    obj.matrices.offSlopes2Command= slopesLinearMMSE(obj.wfs,obj.tel,obj.atm,ast,...
                        'mmseStar',mmseStar,'isTTRM',true,'wavefrontMask', obj.wfs.validActuator);
                else
                    ast = obj.nGs;
                    mmseStar = source('wavelength',obj.nGs(1).photometry);
                    obj.matrices.offSlopes2Command= slopesLinearMMSE(obj.wfs,obj.tel,obj.atm,ast,...
                        'mmseStar',mmseStar);
                end
            end
            
            if numel(obj.nGs) > 1 && ~isempty(obj.lowfs)  
                mmseStar = source('wavelength',obj.nGs(1).photometry);
                obj.matrices.offTipTilt2Command = slopesLinearMMSE(obj.lowfs,obj.tel,obj.atm,...
                    obj.nGs,'mmseStar',mmseStar);
            end
            
            obj.flagCalibration = true;
        end
                        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        %                         AO LOOP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
        
        function obj = runSimulation(obj)
            
            if ~obj.flagCalibration
                obj = obj.calibrateSystem();
            end
                
            %1/ Initialize the loop
            obj   = obj.initLoop(obj.parms);            
            hwait = waitbar(0,'Loop is being closed...');                                    
            dt    = 0;
            count = 0;
            nIter = obj.loopStatus.nIteration;
            % Loop status
            gho   = obj.loopStatus.ho.gain;
            dho   = obj.loopStatus.ho.latency/(obj.tel.samplingTime*obj.loopStatus.ho.frameRate);
            gpolc = obj.loopStatus.ho.gainPolc;
            
            if ~isempty(obj.lowfs)
                gtt      = obj.loopStatus.lo.gain;
                gpolc_tt = obj.loopStatus.lo.gainPolc;
                dtt      = obj.loopStatus.lo.latency/(obj.tel.samplingTime*obj.loopStatus.lo.frameRate);
            end
            if ~isempty(obj.truthwfs)
                gtru  = obj.loopStatus.tru.gain;
                dtru  = obj.loopStatus.tru.latency/(obj.tel.samplingTime*obj.loopStatus.tru.frameRate);
            end
            
            %2/ Initialize the buffers for recording the telemetry
            [nCom,nSl]= size(obj.wfs2dm.M);          
            if ~isempty(obj.lGs)
                nS    = numel(obj.lGs);
            else
                nS    = numel(obj.nGs);
            end                        
            wfsSl = zeros(nSl,nS,nIter);
            
            nDM   = numel(obj.dm);
            dmCom = zeros(nCom,nDM,nIter);
            
            if ~isempty(obj.lowfs)
               ttSl  = zeros(2,numel(obj.nGs),nIter); 
               ttCom = zeros(2,nIter);
            end
            
            if ~isempty(obj.truthwfs)
               truthSl = zeros(2*sum(obj.truthwfs.validSlopes(:)),numel(obj.tGs),nIter); 
               truthCom= zeros(nCom,1,nIter);
            end
            
            
            %3/ Initialize matrices
            
            %3.1./ Reconstructors
            if ~isempty(obj.lowfs)
                % Slopes to commands
                if isempty(obj.parms.loopStatus.lo.tomography) || strcmp(obj.parms.loopStatus.lo.tomography,'mean')
                    Rtt    = obj.nGs(1).wavelength/8;                    
                else
                    spolc_tt = 0*ttSl;
                    Rtt      = obj.matrices.offTipTilt2Command; % Tip-tilt removed automatically if necessary
                end
                            
                T2U    = obj.matrices.tilt2Commands;
                % Tip-tilt filtering
                S2SnoT = obj.matrices.SlopeTTRem;
                U2UnoT = obj.matrices.DMTTRem;
            end
            
            
            if isempty(obj.parms.loopStatus.ho.tomography) || strcmp(obj.parms.loopStatus.ho.tomography,'mean')
                Rho = obj.wfs2dm.M;
                if ~isempty(obj.lowfs)
                    Rho = Rho*S2SnoT;
                end
            else
                spolc = 0*wfsSl;
                Rho   = obj.matrices.offSlopes2Command; % Tip-tilt removed automatically if necessary                
            end
            
            
            if ~isempty(obj.truthwfs)
                Rtru = obj.truth2dm.M;
            end
            
            flush(obj.cam)
            %4/ run the loop
            
            for kIter = 1:nIter          
                tic;
                
                %% 4.1\ Updating phase screens
                +obj.tel;
                                                
                %% 4.2\ Propagating sources to WFSs
                if ~isempty(obj.lowfs)
                    obj.lGs = obj.lGs.*obj.tel*obj.dm*obj.wfs;
                    obj.nGs = obj.nGs.*obj.tel*obj.dm*obj.lowfs;
                else
                    obj.nGs = obj.nGs.*obj.tel*obj.dm*obj.wfs;
                end                
                if ~isempty(obj.tGs)
                    obj.tGs = obj.tGs.*obj.tel*obj.dm*obj.truthwfs;
                end
                                                
                %% 4.3\ Buffering telemetry
                dmCom(:,:,kIter) = obj.dm.coefs;
                
                %4.3.1. Ho WFS
                if mod(kIter,(obj.wfs.camera.exposureTime/obj.tel.samplingTime)*obj.wfs.camera.clockRate) == 0
                    wfsSl(:,:,kIter) = obj.wfs.slopes;
                else
                    if kIter >1
                        wfsSl(:,:,kIter) = wfsSl(:,:,kIter-1);
                    else
                        wfsSl(:,:,kIter) = wfsSl(:,:,1);
                    end
                end
                
                
                %4.3.2. TT WFS
                if ~isempty(obj.lowfs)                    
                    if mod(kIter,(obj.lowfs.camera.exposureTime/obj.tel.samplingTime)*obj.lowfs.camera.clockRate) == 0
                        ttSl(:,:,kIter) = obj.lowfs.slopes;
                    else
                        if kIter >1
                            ttSl(:,:,kIter) = ttSl(:,:,kIter-1);
                        else
                            ttSl(:,:,kIter) = ttSl(:,:,1);
                        end
                    end                    
                end
                
                %4.3.3. Truth WFS
                if ~isempty(obj.tGs)
                    if mod(kIter,(obj.truthwfs.camera.exposureTime/obj.tel.samplingTime)*obj.truthwfs.camera.clockRate) == 0
                        truthSl(:,:,kIter) = obj.truthwfs.slopes;
                    else
                        if kIter >1
                            truthSl(:,:,kIter) = truthSl(:,:,kIter-1);
                        else
                            truthSl(:,:,kIter) = truthSl(:,:,1);
                        end
                    end
                end
                
                %% ------------------------------------------------------ %
                %4.4\ Close the main AO loop                 
                %4.4.1\ Introduce a temporal delay
                sho = obj.introduceDelay(kIter,wfsSl,dho);
                %4.4.2.\ Reconstruct the DM commands                 
                if ~isempty(obj.loopStatus.ho.tomography)
                    if strcmp(obj.loopStatus.ho.tomography,'mean')
                        %GLAO
                        du = Rho*mean(sho,2); 
                    else
                        %LTAO/MCAO
                        spolc(:,:,kIter) = obj.polc(sho,dmCom(:,:,kIter),obj.wfs2dm.D);
                        du = Rho*-spolc(:,:,kIter);
                        du = obj.matrices.dmFilter*du;
                    end
                else
                    du  = Rho*sho;
                end
                
                
                %4.4.3\ Update the DM command
                if kIter > 1
                    dmCom(:,:,kIter) = obj.updateCommands(dmCom(:,:,kIter-1),du,gho,gpolc);
                else
                    dmCom(:,:,kIter) = dmCom(:,:,1);
                end
                                       
                %% ------------------------------------------------------ %
                %4.5\ Close the TT loop if different from the main one
                
                if ~isempty(obj.lowfs)
                    %4.5.1
                    stt = obj.introduceDelay(kIter,ttSl,dtt);
                    %4.5.2
                    if ~isempty(obj.loopStatus.lo.tomography)
                        if strcmp(obj.loopStatus.lo.tomography,'mean')
                            %GLAO
                            du = Rtt*mean(stt,2); 
                        else
                            %LTAO/MCAO
                            spolc_tt(:,:,kIter) = obj.polc(stt,ttCom(:,kIter),8/obj.nGs(1).wavelength);
                            du = Rtt*-spolc_tt(:,:,kIter);
                        end
                    else
                        du  = Rtt*stt;
                    end
                    %4.5.3                    
                    if kIter > 1
                        ttCom(:,kIter) = obj.updateCommands(ttCom(:,kIter-1),du,gtt,gpolc_tt);
                    else
                        ttCom(:,kIter) = ttCom(:,1);
                    end
                    
                    dmCom(:,:,kIter) = U2UnoT*dmCom(:,:,kIter) - T2U*ttCom(:,kIter);
                end                                
                
                %% ------------------------------------------------------ %
                %4.6\ Close the Truth sensor loop if any
                if ~isempty(obj.truthwfs)
                    %4.6.1
                    stru = obj.introduceDelay(kIter,truthSl,dtru);
                    %4.6.2
                    du   = Rtru*stru;
                     if kIter > 1
                        truthCom(:,kIter) = obj.updateCommands(truthCom(:,kIter-1),du,gtru,1);
                    else
                        truthCom(:,kIter) = truthCom(:,1);
                     end
                    %4.6.3
                    dmCom(:,:,kIter) = dmCom(:,:,kIter) - truthCom(:,kIter);
                end
                   
                %% ------------------------------------------------------ %
                
                                                
                %4.7\ Propagates electric field to the imaging camera
                obj.dm.coefs = squeeze(dmCom(:,:,kIter));
                obj.sci = obj.sci.*obj.tel;
                if ~isempty(obj.ncpa.map)
                    obj.sci.phase  = +obj.sci(1).waveNumber*obj.ncpa.map;
                end
                obj.sci = obj.sci*obj.dm*obj.cam;

                % Updating the waiting bar and remaining time
                waitbar(kIter/nIter);                
                dt  = dt + toc();                
                if kIter>1 && kIter<nIter
                    fprintf(1, repmat('\b',1,count)); %delete line before
                    count = fprintf('Remaining simulation time: %0.5g s',dt*(nIter-kIter)/kIter);
                elseif kIter==nIter
                    fprintf(1, repmat('\b',1,count)); %delete line before
                    count = fprintf('Remaining simulation time: %0.5g s\n',dt*(nIter-kIter)/kIter);
                end
                
            end
            
            
            %% ------------------------------------------------------ %
            close(hwait);  
            
            %5\ Save telemetry            
            %5.1\ Get the control loop data
            obj.loopData.slopes = wfsSl(:,:,obj.loopStatus.startDelay+1:end);
            obj.loopData.dmcom  = dmCom(:,:,obj.loopStatus.startDelay+1:end);
            if ~isempty(obj.lowfs)
                obj.loopData.tiptilt = Rtt*ttSl(:,:,obj.loopStatus.startDelay+1:end);
                obj.loopData.tiltCom = ttCom(:,obj.loopStatus.startDelay+1:end);
            else
                obj.loopData.tiptilt = obj.matrices.slopes2Tilt*squeeze(obj.loopData.slopes);
                obj.loopData.tiltCom = obj.matrices.commands2Tilt*squeeze(obj.loopData.dmcom);
            end
            if ~isempty(obj.truthwfs)
                obj.loopData.truthslopes = truthSl(:,obj.loopStatus.startDelay+1:end);
            end
            
            if ~isempty(obj.parms.loopStatus.ho.tomography) && strcmp(obj.parms.loopStatus.ho.tomography,'mmse')
                obj.loopData.slopesPolc = spolc(:,obj.loopStatus.startDelay+1:end);
            end
            
            %5.2\ Get the psf
            npix   = obj.cam.resolution(1);
            ps     = obj.cam.pixelScale*constants.radian2mas;
            tExp   = obj.cam.exposureTime*obj.tel.samplingTime;
            tmp    = kaspPSFStats.empty(0,numel(obj.sci));
            for iSrc=1:numel(obj.sci)
                xi      = 1 + (iSrc-1)*npix;
                xf      = iSrc*npix;
                psfi    = obj.cam.frame(:,xi:xf);
                % define the psf class
                zen     = obj.sci(iSrc).zenith;
                azi     = obj.sci(iSrc).azimuth;
                wvl     = obj.sci(iSrc).wavelength;
                tmp(iSrc) = kaspPSFStats(psfi,double(obj.tel.pupilLogical),obj.tel.D,zen,azi,wvl,ps,tExp,'flagMoffat',true);
            end
            obj.psf = tmp;
            
            obj.flagSimulation = true;
        end
                                       
    end
    
    
    methods (Static)
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        
        function parms = checkParms(parms)
            
            %% 1\ Check Telescope
            D = parms.tel.D;
            %1.1\ Field of view
            sepMax = 2*max(horzcat(hypot(parms.sci.x,parms.sci.y),hypot(parms.nGs.x,parms.nGs.y)));
            if isfield(parms,'lGs') && ~isempty(parms.lGs.x)
                sepMax = max(horzcat(2*hypot(parms.lGs.x,parms.lGs.y),sepMax));
            end
            parms.tel.fieldOfViewInArcsec = max(parms.tel.fieldOfViewInArcsec,sepMax);
            
            %1.2.\ Resolution
            wvl    = [parms.sci.photometry.wavelength,parms.nGs.photometry.wavelength];
            
            if isfield(parms,'lGs') && ~isempty(parms.lGs.x)
                wvl = [wvl,parms.lGs.photometry.wavelength];
            end
            if isfield(parms,'tGs') && ~isempty(parms.tGs.x)
                 wvl = [wvl,parms.lGs.photometry.wavelength];
            end
                
                
            wvl_min= min(wvl);           
            r0_min = parms.atm.r0*(wvl_min/parms.atm.photometry.wavelength)^(6/5);
            sampR0 = parms.tel.resolution*r0_min/D;            
%             if sampR0 < 2
%                 parms.tel.resolution = round(3*D/r0_min);
%                 fprintf("Resizing the pupil resolution with %d pixels to sample the r0 with 3 pixels at %.3d nm \n",parms.tel.resolution,wvl_min*1e9);
%             end
            
            %% 2\ WFS undersampling
            if isfield(parms,'lGs') && ~isempty(parms.lGs.x)
                samp  = parms.wfs.nLenslet*constants.radian2mas*parms.lGs.photometry.wavelength/D/parms.wfs.pixelScale;             
                sampTT= constants.radian2mas*parms.nGs.photometry.wavelength/D/parms.lowfs.pixelScale;
            else
                samp  = parms.wfs.nLenslet*constants.radian2mas*parms.nGs.photometry.wavelength/D/parms.wfs.pixelScale;             
                sampTT= 1;
            end
                                                            
            if samp < 0.5
                % Get the pupil resolution we should have regarding the WFS sampling
                resHR = round(parms.wfs.nLenslet*parms.wfs.nPx/samp);
                if isfield(parms,'lGs') && ~isempty(parms.lGs.x)
                    resHR = max(resHR,round(parms.lowfs.nPx/sampTT));
                end
                %
                resHR = max(parms.tel.resolution,resHR);                
                % Make sure that there's a integer number of pixel per subap
                resHR = parms.wfs.nLenslet*parms.wfs.nPx*round(resHR/(parms.wfs.nLenslet*parms.wfs.nPx));                
                % Update the telescope resolution
                parms.tel.resolution = resHR;
                fprintf("Resizing the pupil at high resolution with %d pixels to take into account the WFS undersampling  \n",parms.tel.resolution);
            end
            
           %% 3\ Loop status
           
           % Check the POLC gain -> set to 1 if no POLC
           if ~strcmp(parms.loopStatus.ho.control,'polc') || ~strcmp(parms.loopStatus.ho.tomography,'mmse')
               parms.loopStatus.ho.gainPolc = 1;               
           end
           
           if isfield(parms.loopStatus,'lo')
               if ~strcmp(parms.loopStatus.lo.control,'polc') || ~strcmp(parms.loopStatus.lo.tomography,'mmse')
               parms.loopStatus.lo.gainPolc = 1;
               end
           end
               
           % Check that MMSE -> POLC if multiple stars
           if strcmp(parms.loopStatus.ho.tomography,'mmse') && (~numel(parms.nGs.x)>1 || ~numel(parms.lGs.x)>1)
               parms.loopStatus.ho.control = 'polc';
           else
               parms.loopStatus.ho.control = 'integrator';
               parms.loopStatus.ho.gainPolc = 1;
           end
        
           
           if isfield(parms.loopStatus,'lo')
               if strcmp(parms.loopStatus.lo.tomography,'mmse') && numel(parms.nGs.x)>1
                   parms.loopStatus.lo.control = 'polc';
               else
                   parms.loopStatus.lo.control = 'integrator';
                   parms.loopStatus.lo.gainPolc = 1;
               end
           end
            
            % Check that GLAO -> integrator controller
            if strcmp(parms.loopStatus.ho.tomography,'mean')
                parms.loopStatus.ho.control = 'integrator';
                parms.loopStatus.ho.gainPolc = 1;
            end
            
             if isfield(parms.loopStatus,'lo') 
                 if strcmp(parms.loopStatus.lo.tomography,'mean')
                     parms.loopStatus.lo.control = 'integrator';
                     parms.loopStatus.lo.gainPolc = 1;
                 end
             end
            
        end
        
        function sci = initSource(parms,parmsAtm)
            
            if ~isempty(parms.x) && ~isempty(parms.y) && isscalar(parms.magnitude)
                
                [zenith,azimuth] = utilities.arcsec2polar(parms.x,parms.y);
                
                if isfield(parms,'height') && ~isempty(parms.height) && parms.height~=Inf
                    airmass = 1/parmsAtm.zenithAngle;
                    meanH   = mean(parms.height) * airmass;
                    sci = laserGuideStar(parms.apertureSize,parms.apertureDistance,...
                        meanH,parms.spotFwhm,parms.nPhoton,parms.naProfile,'zenith',zenith,'azimuth',azimuth,'height',...
                        parms.height*airmass,'wavelength',parms.photometry);
                else
                    sci  = source('zenith',zenith,'azimuth',azimuth,...
                        'wavelength',parms.photometry,'magnitude',parms.magnitude);
                end
            else
                sci = [];
            end
        end
        
        function tel = initTelescope(parmsTel,parmsSeg)
            
            % Create the segmented pupil
            
            if isfield(parmsSeg,'usePupilClass') && ~isempty(parmsSeg.usePupilClass) && parmsSeg.usePupilClass
                
                seg = segment(parmsSeg.nSides,parmsSeg.radius,parmsSeg.resolution);
                
                if ~isfield(parmsSeg,'pupilUnit')
                    pupilUnit='m';
                end
                
                
                if ~isfield(parmsSeg,'segmentCoordinates')
                    if strcmp(pupilUnit,'m')
                        segmentCoordinates=[0 0];
                    elseif strcmp(pupilUnit,'px')
                        [sx,sy]=size(seg);
                        segmentCoordinates=[round(sx/2) round(sy/2)];
                    end
                else
                    segmentCoordinates = parmsSeg.segmentCoordinates;
                end
                
                pupil = pupilGenerator('segRef',seg,...
                    'segCoord',segmentCoordinates,...
                    'D',parmsTel.D,'obstructionRatio',parmsTel.obstructionRatio,...
                    'unit',parmsSeg.pupilUnit,'flagNoGap',1);
                
                if isfield(parmsSeg,'pupilSpider')
                    pupil.spiders = parmsSeg.pupilSpider;
                    pupil.mkSpiders(parmsSeg.pupilSpider);
                end
                if isfield(parmsSeg,'pupilReflexivity')
                    pupil.coeffReflexion = parmsSeg.pupilReflexivity;
                    pupil.applyReflexivity(1:pupil.nSegments,parmsSeg.pupilReflexivity);
                end
                if isfield(parmsSeg,'pupilCoeffPhaseModes')
                    pupil.coeffPhaseModes = parmsSeg.pupilCoeffPhaseModes;
                    pupil.applyPhaseError(1:pupil.nSegments,parmsSeg.pupilCoeffPhaseModes);
                end
                
                if exist('parmsSegoptionalPupilCommand','var')
                    for k=1:length(parmsSeg.optionalPupilCommand)
                        eval(char(parmsSeg.optionalPupilCommand(k)));
                    end
                end
                pupilClassMatrix = pupil.resize(parmsTel.resolution,0,'nearest');
            end
            
            
            % Create the telescope
            
            tel = telescope(parmsTel.D,'resolution',parmsTel.resolution,...
                'obstructionRatio',parmsTel.obstructionRatio,'fieldOfViewInArcsec',...
                parmsTel.fieldOfViewInArcsec ,'samplingTime',parmsTel.samplingTime);
            
            if isfield(parmsTel,'pupilShape') && ~isempty(parmsTel.pupilShape)
                tel.pupil = parmsTel.pupilShape;
                disp('___ PUPIL___');
                disp('User-defined pupil');
                disp('----------------------------------------------------');
            elseif exist('pupilClassMatrix','var')
                tel.pupil = abs(pupilClassMatrix);
                disp('___ PUPIL___');
                disp('Set telescope pupil from pupil Class');
                disp('----------------------------------------------------');
            else
                disp('___ PUPIL___');
                disp('Default built-in telescope pupil');
                disp('----------------------------------------------------');
            end
        end
        
        function atm = initAtmosphere(parms)
            
            airmass = 1/cos(parms.zenithAngle); % zenithAngleInRad
            atm = atmosphere(parms.photometry,parms.r0*airmass^(-3/5),mean(parms.L0),'layeredL0',parms.L0,...
                'fractionnalR0',parms.fractionalR0,'altitude',parms.altitude*airmass,...
                'windSpeed',parms.windSpeed,'windDirection',parms.windDirection);
            
        end
        
        function wfs = initWfs(parmsWfs,tel,src)
            
            %1/ Define the WFS class
            lond       = parmsWfs.nLenslet*constants.radian2mas*src.photometry.wavelength/tel.D;
            samp       = lond/parmsWfs.pixelScale; % samp = 1 -> Nyquist sampling
            resolution = parmsWfs.nLenslet*parmsWfs.nPx;
            fov        = parmsWfs.pixelScale*parmsWfs.nLenslet/lond; % fov in lambda/D units
            
            if samp >= 0.5
                wfs = shackHartmann(parmsWfs.nLenslet,resolution,parmsWfs.minLightRatio);                
            else
                wfs                          = shackHartmann(parmsWfs.nLenslet,2*tel.resolution,parmsWfs.minLightRatio);
                wfs.lenslets.nyquistSampling = 0.5;
                wfs.camera.resolution        = [resolution resolution];
            end
            
            %2/ Initialize the WFS
            src = src.*tel*wfs;
            wfs.INIT;
            +wfs;
            
            %3/ Calibrate the WFS optical gain
            pixelScale      = wfs.lenslets.nyquistSampling*wfs.lenslets.fieldStopSize*lond/parmsWfs.nPx;
            wfs.slopesUnits = calibrateWfsPixelScale(wfs,src,pixelScale,parmsWfs.nPx);
            fprintf('WFS pixel scale set to %4.2f arcsec/pixel\n',pixelScale);
            
        end
        
        function dm = initDm(parms,resolution,validMap)
            
            switch parms.influenceType
                
                case 'gaussian'
                    bif = gaussianInfluenceFunction(parms.crossCoupling,parms.pitch);
                    
                case 'monotonic'
                    bif = influenceFunction('monotonic',parms.crossCoupling,parms.influenceCentre);
                    
                case 'overshoot'
                    bif = influenceFunction('overshoot',parms.crossCoupling,parms.influenceCentre);
                    
                otherwise
                    bif = influenceFunction(parms.influencePoints,parms.crossCoupling,parms.influenceCentre);
            end
            
            
            if isfield(parms,'validActuMap')
                validActuatorMap = logical(parms.validActuMap);
            elseif ~isfield(parms,'validActuMap') && ~isempty(validMap)
                validActuatorMap = validMap;
            else
                validActuatorMap = true(parms.nActuators);
            end
            
            dm      = deformableMirror(parms.nActuators,'modes',bif,...
                'resolution',resolution,'validActuator',validActuatorMap);
            
        end
        
        function cam = initImager(parms,tel,src)
            
            %1/ Define the imager class
            loD  = constants.radian2mas*src.wavelength/tel.D;
            samp = loD/parms.pixelScale/2;
            
            if samp >= 0.5
                cam = imager('nyquistSampling',samp,'fieldStopSize',parms.resolution/2/samp);
                cam.resolution = round(cam.resolution);
            else
                cam = imager('nyquistSampling',0.5,'fieldStopSize',parms.resolution/samp/2);
                cam.resolution = round(parms.resolution);
            end
            
            %2/ Get the reference frame
            src               = src.*tel*cam;
            cam.referenceFrame= repmat(cam.frame,[1,numel(src)]);
            flush(cam,numel(src));
            
            %3/ Initialize the detector paraaeters
            cam.clockRate   = parms.clockRate;
            cam.exposureTime= parms.exposureTime;
            cam.startDelay  = parms.startDelay;
            cam.readOutNoise= parms.ron;
            cam.photonNoise = parms.photonNoise;
            cam.pixelScale  = cam.imgLens.fieldStopSize*loD/cam.resolution(1)/constants.radian2mas;
            
        end
        
        function dmCalib = calibrateInteractionMatrix(parm,tel,dm,wfs,src)
            
            %1/ Calibrate the interaction matrix
            if strcmp(parm.controlBasis,'modal')
                nMode = dm.nValidActuator;
                zer   = zernike(2:nMode+1,'resolution',tel.resolution);
                dmCal = deformableMirror(nMode,'modes',zer,'resolution',...
                    tel.resolution,'validActuator',true(1,nMode));
            else
                dmCal = dm;
            end
            
            src     = src.*tel;
            dmCalib = calibration(dmCal,wfs,src,src.wavelength/40,1);
            
            %2/ Calibrate the reconstructor
            
            
            if strcmp(parm.controlBasis,'modal')
                
                % Manage the model gain
                if isfield(parm,'modalGain')
                    if length(parm.modalGain) >= nMode
                        modalGain = parm.modalGainmodalGain(1:nMode);
                    elseif length(parm.modalGainmodalGain) < nMode
                        modalGain = ones(nMode);
                        modalGain(1:length(modalGain)) = parm.modalGain;
                    end
                else
                    modalGain = 1;
                end
                
                Fmat  = 2*dm.modes.modes(tel.pupilLogical,:);
                iFmat = pinv(full(Fmat),1e-1);
                M2U   = iFmat*zer.modes(tel.pupilLogical,:);
                dmCalib.nThresholded = parm.nModeFiltered;
                dmCalib.M            = M2U*bsxfun(@times,dmCalib.M,modalGain);
            else
                dmCalib.nThresholded = parm.nModeFiltered;
            end
        end
        
        function s2 = introduceDelay(k,s,lat)
                        
            
            if lat ==0
                s2 = s(:,:,k);
            elseif lat ==1 && k>1
                s2 = s(:,:,k-1);
            elseif lat == 2 && k>2
                s2 = s(:,:,k-2);
            elseif lat > 1 && lat < 2 && k>2
                s2 = lat*s(:,:,k-2) + (1-lat)*s(:,:,k-1);
            else
                s2 = s(:,:,1);
            end
        end
        
        function u2 = updateCommands(u,du,gloop,gpolc)
            u2 = gpolc*u - gloop*du;
        end
        
        function s2 = polc(s,u,D)
            s2  = bsxfun(@minus,s,D*u);
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        %                             DEMO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function demoAOsystem()
            
            % NGS SCAO
            devel.parFileTemplate()
            sys = aoSystem(parm);
            sys = sys.runSimulation();
            imagesc(sys.psf.image);
        end
        
        
    end
end