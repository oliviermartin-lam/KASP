classdef errorBreakdownSimulation < handle
    
    properties (SetObservable=true)
        % ------------------------- SUB-CLASSES ------------------------- %
        sys;                    % SYSTEM CLASS, see systemSimulator.m
        fao;                    % Spatial Fourier AO class
        wfs;
        coefs;
        coefsTT;
        P;
        telSq;       
        ast;
        ngs;
        astTT;
        astSrc;
        ttgsOn;
        paraFlag;
        psf;
        % ------------------------- PARAMETERS -------------------------- %
        nGs;
        nTT;
        nTimes;
        resolution;
        pitch;
        nSize;
        nActu;
        nIter;
        pixelScale;
        nSrc;
        noiseLgs;
        noiseTT;
        nFrames;                        
        startDelay;
        zerTT;
        zerHO;
        % --------------------------- OUTPUTS --------------------------- %
        window;
        psdFit;
        nPixelPerPitch;
        % LGS
        slopesGs;        
        slopesServo;
        slopesAlias;
        slopesNoise;        
        slopesAliasOL;                      
        %TIP-TILT
        slopesTT;        
        slopesServoTT;
        slopesAliasTT;
        slopesNoiseTT;                 
        % ----------------------------- WFE ------------------------------%
        wfeFourier;
        wfeSimu;
        phiAliasOL;
        sAliasOL;
    end
   
    
    methods 
        
        function obj = errorBreakdownSimulation(sys,varargin)
            
            inputs = inputParser;
            inputs.addRequired('sys',@(x) isa(x,'aoSystem') );
            inputs.addParameter('nPixelPerPitch',2*sys.wfs.lenslets.nLenslet,@isnumeric);
            inputs.parse(sys,varargin{:});
            %% Inputs from system
            obj.sys        = sys;
            obj.nPixelPerPitch = inputs.Results.nPixelPerPitch;           
            %% GEOMETRY           
            
            obj.nTimes     = floor(sys.src(1).wavelength/sys.tel.D/sys.cam.pixelScale);
            obj.resolution = sys.tel.resolution;%obj.nPixelPerPitch * (obj.nActu - 1);
            obj.nSize      = obj.nTimes * obj.resolution;
             
            %% Define a telescope class without atmosphere
            obj.P   = telescope(sys.tel.D,'resolution',obj.resolution,'obstructionRatio',...
                sys.tel.obstructionRatio,'fieldOfViewInArcsec',sys.tel.fieldOfView*constants.radian2arcsec...
                ,'samplingTime',sys.tel.samplingTime);
            
            obj.P.pupil = sys.tel.pupil;
            
            %% Define a telescope class with a square pupil
            obj.telSq = telescope(sys.tel.D,'resolution',obj.resolution,'fieldOfViewInArcsec',...
                sys.tel.fieldOfView*constants.radian2arcsec,'samplingTime',sys.tel.samplingTime);
            
            obj.telSq.pupil = ones(obj.resolution);            
                      
            %% Define new asterisms
            
            obj.astSrc = source('zenith',[sys.src.zenith],'azimuth',[sys.src.azimuth]...
                ,'wavelength',[sys.ast.photometry],'magnitude',[sys.ast.magnitude]);
            
             obj.ngs = source('zenith',[sys.ast.zenith],'azimuth',[sys.ast.azimuth]...
                ,'wavelength',[sys.ast.photometry],'magnitude',[sys.ast.magnitude]);
            if obj.nTT
                obj.ngs.photometry = obj.sys.ast_tt.photometry;
            end
            %% Defining DM classes
                                           
            obj.coefs   = zeros(obj.sys.dm.nValidActuator,4,obj.nIter);
            obj.coefsTT = zeros(obj.sys.dm.nValidActuator,4,obj.nIter);
                                                                      
            if ~any(obj.sys.FiF)
                bif      = gaussianInfluenceFunction(obj.sys.bif.mechCoupling,obj.sys.tel.D/(obj.sys.dm.nActuator-1));
                obj.sys.dmSq = deformableMirror(obj.sys.dm.nActuator,'modes',bif,'resolution',obj.sys.tel.resolution...
                    ,'validActuator',logical(ones(obj.sys.dm.nActuator)));
                F            = bif.modes;
                obj.sys.FiF  = pinv(full(F),1e-1);
            end
                
            %% Instantiating ouputs
            obj.window      = mywindow('hann',obj.resolution);          
            obj.psdFit      = zeros(obj.resolution);
         
            %% Create Zernike modes
            obj.zerTT = zernike(2:3,'D',obj.sys.tel.D,'resolution',...
                obj.resolution,'pupil',obj.P.pupil, ...
                'fieldOfViewInArcmin',obj.sys.tel.fieldOfView/60.);
            
            obj.zerHO = zernike(4:100,'D',obj.sys.tel.D,'resolution',...
                obj.resolution,'pupil',obj.P.pupil, ...
                'fieldOfViewInArcmin',obj.sys.tel.fieldOfView/60.);
            
        end
        
        
        function obj = runSplitSimulation(obj)
            
            %Allocation in memory: outputs
            varRes      = zeros(obj.sys.nIter,obj.sys.nSrc);         
            varAtm      = zeros(obj.sys.nIter,1);
            varPara     = zeros(obj.sys.nIter,1);
            varPerp     = zeros(obj.sys.nIter,1);
            varAliasOL  = zeros(obj.sys.nIter,1);
            varAlias    = zeros(obj.sys.nIter,1);
            varServo    = zeros(obj.sys.nIter,1);
            varNoise    = zeros(obj.sys.nIter,1);
            varResWfe   = zeros(obj.sys.nIter,1);
            varAliasTT  = zeros(obj.sys.nIter,1);            
            varServoTT  = zeros(obj.sys.nIter,1);
            varNoiseTT  = zeros(obj.sys.nIter,1);
            varResWfeTT = zeros(obj.sys.nIter,1);
            varAniso    = zeros(obj.sys.nIter,obj.sys.nSrc);
            varAnisoTT  = zeros(obj.sys.nIter,obj.sys.nSrc);
            %Allocation in memory: LGS slopes
            sGs         = zeros(obj.sys.wfs.nSlope,obj.sys.nGs,obj.sys.nIter);
            sAliasOL    = zeros(obj.sys.wfs.nSlope,obj.sys.nGs,obj.sys.nIter);
            sAlias      = zeros(obj.sys.wfs.nSlope,obj.sys.nGs,obj.sys.nIter);
            sServo      = zeros(obj.sys.wfs.nSlope,obj.sys.nGs,obj.sys.nIter);
            sNoise      = zeros(obj.sys.wfs.nSlope,obj.sys.nGs,obj.sys.nIter);
            phNoise     = obj.sys.wfs.camera.photonNoise;
            roNoise     = obj.sys.wfs.camera.readOutNoise;
            if obj.sys.nTT
                sTT       = zeros(2,obj.sys.nTT,obj.sys.nIter);
                sAliasTT  = zeros(2,obj.sys.nTT,obj.sys.nIter);
                sServoTT  = zeros(2,obj.sys.nTT,obj.sys.nIter);
                sNoiseTT  = zeros(2,obj.sys.nTT,obj.sys.nIter);
                phNoiseTT = obj.sys.wfs_tt.camera.photonNoise;
                roNoiseTT = obj.sys.wfs_tt.camera.readOutNoise;
            else
                phNoiseTT = 0;
                roNoiseTT = 0;
            end
            
            % System instantiation
            obj.sys.dm.coefs   = zeros(obj.sys.dm.nValidActuator,1);            
            obj.sys.wfs.slopes = zeros(obj.sys.wfs.nSlope,1);
            obj.coefs          = zeros(obj.sys.dm.nValidActuator,4,obj.sys.nIter);           

            if obj.nTT
                obj.coefsTT           = zeros(obj.sys.dm.nValidActuator,4,obj.sys.nIter);
                obj.sys.wfs_tt.slopes = zeros(2,1);
            end
            % Temporary tel and atmosphere
            obj.telSq       = psfrTools.duplicateTelescope(obj.telSq);
            obj.telSq.pupil = ones(obj.telSq.resolution);
            atmSq           = psfrTools.duplicateAtmosphere(obj.sys.atm);
            obj.telSq       = obj.telSq + atmSq;               
            tmpTel          = psfrTools.duplicateTelescope(obj.sys.tel);            
            tmpAtm          = psfrTools.duplicateAtmosphere(obj.sys.atm);
            tmpTel          = tmpTel + tmpAtm;
            pupLog          = obj.P.pupilLogical;                                 
            
            dt    = 0;
            count = 0;
            hwait = waitbar(0,'Loop is being closed...'); 
            for k=1:obj.sys.nIter
                %% 1. MAIN AO LOOP : grab system telmetry
                tic;
                % Updating phase screens
                +tmpTel;
                obj.sys.rtc.kIteration = k;
                % Propagating sources through telescope+atmosphere
                obj.sys.src       = obj.sys.src.*tmpTel;                
                obj.sys.ast       = obj.sys.ast.*tmpTel;                
                
                if obj.sys.nTT
                    obj.sys.ast_tt      = obj.sys.ast_tt.*tmpTel;
                    obj.sys.rtc.ast_tt  = obj.sys.ast_tt;
                end
                
                % Buffering telemetry              
                if k > obj.sys.rtc.delay
                    obj.sys.rtc.compensatorData(:,k-1) = obj.coefs(:,1,k-1);
                    obj.sys.rtc.compensator.coefs      = obj.coefs(:,1,k-1);
                    obj.sys.rtc.sensorData(:,:,k-1)    = sGs(:,:,k-1);
                    if obj.sys.nTT
                        obj.sys.rtc.tipTiltCompensatorData(:,k-1) = obj.coefsTT(:,1,k-1); 
                        obj.sys.rtc.tipTiltSensorData(:,:,k-1)    = sTT(:,:,k-1);
                    end
                end
                                
                % Putting back noise if required
                obj.sys.wfs.camera.photonNoise  = phNoise;
                obj.sys.wfs.camera.readOutNoise = roNoise;
                obj.sys.wfs.framePixelThreshold = roNoise;
                
                if obj.nTT
                    obj.sys.wfs_tt.camera.photonNoise  = phNoiseTT;
                end
                    
                % Closing the loop                            
                obj.sys.rtc.closedAOLoop();
                obj.coefs(:,1,k) = obj.sys.rtc.compensatorData(:,k);                
                sGs(:,:,k)       = obj.sys.rtc.sensorData(:,:,k);
                
                % Propagates science target through the system
                obj.sys.src = obj.sys.src*obj.sys.dm*obj.sys.cam;
                for i=1:obj.sys.nSrc
                    varRes(k,i) = var(obj.sys.src(i).phase(tmpTel.pupilLogical));
                end
                obj.sys.ast = obj.sys.ast.*tmpTel*obj.sys.dm;
                phiRes      = obj.sys.ast.phase;
                
                if obj.sys.nTT
                    obj.coefsTT(:,1,k) = obj.sys.rtc.tipTiltCompensatorData(:,k);
                    sTT(:,:,k)         = obj.sys.rtc.tipTiltSensorData(:,:,k);
                    obj.ngs            = obj.ngs.*tmpTel*obj.sys.dm;
                    phiResTT           = obj.ngs.phase;
                end
                
                %% 2. SPLIT PHASE : get parallel/perpendicular phases
                +obj.telSq;
                % Nullifying noise
                obj.sys.wfs.camera.photonNoise  = 0;
                obj.sys.wfs.camera.readOutNoise = 0;
                obj.sys.wfs.framePixelThreshold = 0;
                
                if obj.sys.nTT
                    obj.sys.wfs_tt.camera.photonNoise  = 0;
                    obj.sys.wfs_tt.camera.readOutNoise = 0;
                    obj.sys.rtc.ast_tt                 = obj.ngs;
                end
                
                % Getting phases in sensing directions
                obj.sys.ast = obj.sys.ast.*obj.telSq;
                phiAtmo     = obj.sys.ast.phase;
                phiPara     = psfrTools.dmFilter(phiAtmo,obj.sys.dmSq,obj.sys.FiF);                
                phiPerp     = phiAtmo - phiPara;
                
                if obj.sys.nTT
                    obj.ngs   = obj.ngs.*obj.telSq;
                    phiAtmoTT = obj.ngs.phase;
                    phiPara   = psfrTools.dmFilter(phiAtmoTT,obj.sys.dmSq,obj.sys.FiF);                            
                    phiPerpTT = phiAtmoTT - phiParaTT;
                end
                
                % Computing variance/covariance of the phase
                varAtm(k)   = var(phiAtmo(pupLog));
                varPara(k)  = var(phiPara(pupLog));
                varPerp(k)  = var(phiPerp(pupLog));

                % Grabbing the fitting PSD
                tmp_phi     = obj.P.pupil.*(phiPerp - mean(phiPerp(pupLog)));
                obj.psdFit  = obj.psdFit + abs(fft2(obj.window.*tmp_phi)).^2;
                
                %% 3. ANISOPLANATISM : get anisoplanatism errors
                
                for iSrc=1:obj.sys.nSrc
                    if obj.sys.nTT
                        % Focal + Angular anisoplanatism
                        obj.sys.ast      = obj.sys.ast.*obj.telSq;
                        z                = obj.zerHO\obj.sys.ast;
                        phiHOoff         = z.modes*z.c;
                        obj.astSrc(iSrc) = obj.astSrc(iSrc).*obj.telSq;
                        zOn              = obj.zerHO\obj.astSrc;
                        phiHOon          = zOn.modes*zOn.c*obj.sys.ast.waveNumber;
                        
                        phiAniso = reshape(phiHOon-phiHOoff,z.resolution,z.resolution);
                        varAniso(k,iSrc) = var(phiAniso(pupLog));
                        
                        % Tip-tilt anisoplanatism
                        obj.sys.ast_tt = obj.sys.ast_tt.*obj.telSq;
                        z              = obj.zerTT\obj.sys.ast_tt;
                        phiTToff       = z.modes*z.c;
                        zOn            = obj.zerTT\obj.astSrc(iSrc);
                        phiTTon        = zOn.modes*zOn.c*obj.sys.ast_tt.waveNumber;
                                                    
                        phiTTAniso = reshape(phiTTon-phiTToff,z.resolution,z.resolution)*obj.sys.ast_tt.waveNumber;
                        varAnisoTT(k,iSrc)  = var(phiTTAniso(pupLog));
                    else
                        obj.astSrc(iSrc) = obj.astSrc(iSrc).*obj.telSq;
                        phiAtmoOn        = obj.astSrc(iSrc).phase;                        
                        phiAniso         = phiAtmo - phiAtmoOn;
                        varAniso(k,iSrc) = var(phiAniso(pupLog));
                    end
                end
                
                %% 4. SERVO-LAG ERROR : pure temporal servo-lag error
                
                % Putting the parallel phase in front of WFSs
                obj.sys.ast       = obj.sys.ast.*obj.P;
                obj.sys.ast.phase = phiPara+obj.sys.ast.phase;
                if obj.sys.nTT
                    obj.ngs       = obj.ngs.*obj.P;
                    obj.ngs.phase = phiParaTT+obj.sys.ngs.phase;
                end
                
                % Buffering telemetry
                if k > obj.sys.rtc.delay
                    obj.sys.rtc.compensatorData(:,k-1) = obj.coefs(:,2,k-1);                    
                    obj.sys.rtc.compensator.coefs      = obj.coefs(:,2,k-1);  
                    obj.sys.rtc.sensorData(:,:,k-1)    = sServo(:,:,k-1);
                    if obj.nTT
                        obj.sys.rtc.tipTiltCompensatorData(:,k-1) = obj.coefsTT(:,2,k-1);
                        obj.sys.rtc.tipTiltSensorData(:,:,k-1)    = sServoTT(:,:,k-1);
                    end
                end
                % Closing the loop                
                obj.sys.rtc.closedAOLoop();
                obj.coefs(:,2,k) = obj.sys.rtc.compensatorData(:,k);                
                sServo(:,:,k)    = obj.sys.rtc.sensorData(:,:,k);
                
                % Residual High-order modes phase
                phiServo = phiPara - 2.*obj.sys.rtc.compensator.surface.*obj.sys.ast.waveNumber;
                
                if obj.sys.nTT
                    obj.sys.ast       = obj.sys.ast.*obj.P;
                    obj.sys.ast.phase = phiServo;
                    z                 = obj.zerHO\obj.sys.ast;
                    phiServo          = reshape(z.modes*z.c,z.resolution,z.resolution)*obj.sys.ast.waveNumber;
                    % Residual tip-tilt modes phase
                    obj.coefsTT(:,2,k)= obj.sys.rtc.tipTiltCompensatorData(:,k);
                    sServoTT(:,:,k)   = obj.sys.rtc.tipTiltSensorData(:,:,k);
                    phiServoTT        = phiParaTT - 2.*obj.sys.dm.surface.*obj.sys.ast_tt.waveNumber;                    
                    obj.ngs           = obj.ngs.*obj.P;
                    obj.ngs.phase     = phiServoTT;
                    z                 = obj.zerTT\obj.ngs;
                    phiServoTT        = reshape(z.modes*z.c,z.resolution,z.resolution)*obj.ngs.waveNumber;
                    varServoTT(k)     = var(phiServoTT(pupLog));
                end
                
                varServo(k) = var(phiServo(pupLog));
                                              
                
                %% 5. ALIASING : propagation of aliased phased from phase sampling by WFSs
                
                % ------- 5.1. Open-loop aliasing : loop gain set to zero
                %Phase propagation
                obj.sys.ast       = obj.sys.ast.*obj.P;
                obj.sys.ast.phase = phiPerp;
                obj.sys.ast       = obj.sys.ast*obj.sys.wfs;
                
                if obj.sys.SHflag
                    sAliasOL(:,k) = obj.sys.wfs.slopes;
                else
                    sAliasOL(:,k) = obj.sys.G*obj.sys.ast.phase(pupLog);
                end
                obj.sAliasOL=sAliasOL;
                % Phase reconstruction
                obj.sys.dm.coefs = -obj.sys.R*sAliasOL(:,k);
                if obj.sys.nTT
                    obj.sys.dm.coefs = obj.sys.DMTTRem*obj.sys.dm.coefs;
                end
                
                phiAliasOL       = 2.*obj.sys.dm.surface*obj.sys.ast.waveNumber;
                obj.phiAliasOL       = phiAliasOL;
                
                varAliasOL(k)    = var(phiAliasOL(pupLog));                                
                
                % ---------- 5.2 CLOSED-LOOP ALIASING
                % Putting the perpendicular phase in front of WFSs
                obj.sys.ast       = obj.sys.ast.*obj.P;
                obj.sys.ast.phase = phiPerp;
                if obj.sys.nTT
                    obj.ngs       = obj.ngs.*obj.P;
                    obj.ngs.phase = phiPerpTT;
                end
                
                % Buffering telemetry
                if k > obj.sys.rtc.delay
                    obj.sys.rtc.compensatorData(:,k-1) = obj.coefs(:,3,k-1);
                    obj.sys.rtc.compensator.coefs      = obj.coefs(:,3,k-1);
                    obj.sys.rtc.sensorData(:,:,k-1)    = sAlias(:,:,k-1);
                    if obj.nTT
                        obj.sys.rtc.tipTiltCompensatorData(:,k-1) = obj.coefsTT(:,3,k-1);
                        obj.sys.rtc.tipTiltSensorData(:,:,k-1)    = sAliasTT(:,:,k-1);
                    end
                end
                % Closing the loop
                obj.sys.rtc.closedAOLoop();
                obj.coefs(:,3,k) = obj.sys.rtc.compensatorData(:,k);
                sAlias(:,:,k)    = obj.sys.rtc.sensorData(:,:,k);
                
                % Residual High-order modes phase
                phiAlias  = - 2.*obj.sys.rtc.compensator.surface.*obj.sys.ast.waveNumber;
                               
                if obj.sys.nTT
                    obj.sys.ast       = obj.sys.ast.*obj.P;
                    obj.sys.ast.phase = phiAlias; 
                    z                 = obj.zerHO\obj.sys.ast;
                    phiAlias          = reshape(z.modes*z.c,z.resolution,z.resolution)*obj.sys.ast.waveNumber;
                    % Residual tip-tilt modes phase
                    obj.coefsTT(:,3,k)= obj.sys.rtc.tipTiltCompensatorData(:,k);
                    sAliasTT(:,:,k)   = obj.sys.rtc.tipTiltSensorData(:,:,k);
                    phiAliasTT        = -2.*obj.sys.rtc.compensator.surface.*obj.sys.ast_tt.waveNumber;
                    obj.ngs           = obj.ngs.*obj.P;
                    obj.ngs.phase     = phiAliasTT;
                    z                 = obj.zerTT\obj.ngs;
                    phiAliasTT        = reshape(z.modes*z.c,z.resolution,z.resolution)*obj.ngs.waveNumber;
                    varAliasTT(k)     = var(phiAliasTT(pupLog));
                end
                
                varAlias(k)       = var(phiAlias(pupLog));
                                                
                %% 6. NOISE ERROR : WFS noise propagation
                phiNoise   = 0;
                phiNoiseTT = 0;
                if phNoise || roNoise || phNoiseTT || roNoiseTT
                    % Putting back noise
                    obj.sys.wfs.camera.photonNoise  = phNoise;
                    obj.sys.wfs.camera.readOutNoise = roNoise;
                    obj.sys.wfs.framePixelThreshold = roNoise/2;
                    
                    if obj.sys.nTT
                        obj.sys.wfs_tt.camera.photonNoise  = phNoiseTT;
                        obj.sys.wfs_tt.camera.readOutNoise = roNoiseTT;
                    end
                    
                    % Propagating
                    obj.sys.ast       = obj.sys.ast.*obj.P;
                    obj.sys.ast.phase = obj.sys.ast.phase+phiPara;
                    if obj.sys.nTT
                        obj.ngs       = obj.ngs.*obj.P;
                        obj.ngs.phase = obj.sys.ngs.phase+phiParaTT;
                    end
                    % Buffering telemetry
                    if k > obj.sys.rtc.delay
                        obj.sys.rtc.compensatorData(:,k-1) = obj.coefs(:,4,k-1);
                        obj.sys.rtc.compensator.coefs      = obj.coefs(:,4,k-1);
                        obj.sys.rtc.sensorData(:,:,k-1)    = sNoise(:,:,k-1);
                        if obj.sys.nTT
                            obj.sys.rtc.tipTiltCompensatorData(:,k-1) = obj.coefsTT(:,4,k-1);   
                            obj.sys.rtc.tipTiltSensorData(:,:,k-1)    = sNoiseTT(:,:,k-1);
                        end
                    end
                    % Closing the loop                    
                    obj.sys.rtc.closedAOLoop();
                    obj.coefs(:,4,k) = obj.sys.rtc.compensatorData(:,k);
                    sNoise(:,:,k)    = obj.sys.rtc.sensorData(:,:,k);
                    
                    % Noise high-order modes phase
                    phiNoise = phiPara - 2.*obj.sys.rtc.compensator.surface.*obj.sys.ast.waveNumber- phiServo;

                    if obj.sys.nTT                        
                        obj.sys.ast        = obj.sys.ast.*obj.P;
                        obj.sys.ast.phase  = phiNoise;
                        z                  = obj.zerHO\obj.sys.ast;
                        phiNoise           = reshape(z.modes*z.c,z.resolution,z.resolution)*obj.sys.ast.waveNumber;
                        
                        % Noise tip-tilt modes phase
                        obj.coefsTT(:,4,k) = obj.sys.rtc.tipTiltCompensatorData(:,k);
                        sNoiseTT(:,:,k)    = obj.sys.rtc.tipTiltSensorData(:,:,k);
                        phiNoiseTT         = phiParaTT - 2.*obj.sys.dm.surface*obj.sys.ast_tt.waveNumber - phiServoTT;
                        obj.ngs            = obj.ngs.*obj.P;
                        obj.ngs.phase      = phiNoiseTT;
                        z                  = obj.zerTT\obj.ngs;
                        phiNoiseTT         = reshape(z.modes*z.c,z.resolution,z.resolution)*obj.ngs.waveNumber;
                        varNoiseTT(k)      = var(phiNoiseTT(pupLog));
                    end
                    varNoise(k) = var(phiNoise(pupLog));
                end                
             
                 %% 6. RESIDUAL WFE ERROR : Simulated residual minus calculated one                
                phiResWfe         = phiRes - phiServo - phiPerp - phiNoise - phiAlias;               
                
                if obj.sys.nTT  
                    obj.sys.ast       = obj.sys.ast.*obj.P;
                    obj.sys.ast.phase = phiResWfe;
                    z                 = obj.zerHO\obj.sys.ast;
                    phiResWfe         = reshape(z.modes*z.c,z.resolution,z.resolution)*obj.sys.ast.waveNumber;
                    
                    % Aliasing tip-tilt modes phase
                    phiResWfeTT       = phiResTT - phiServoTT - phiPerpTT - phiNoiseTT - phiAliasTT;
                    obj.ngs           = obj.ngs.*obj.P;
                    obj.ngs.phase     = phiResWfeTT;
                    z                 = obj.zerTT\obj.ngs;
                    phiResWfeTT       = reshape(z.modes*z.c,z.resolution,z.resolution)*obj.ngs.waveNumber;
                    varResWfeTT(k)    = var(phiResWfeTT(pupLog));
                end
                varResWfe(k)    = var(phiResWfe(pupLog));
                
                
             waitbar(k/obj.sys.nIter);
                                       
             dt  = dt + toc();
             
             if k>1 & k<obj.sys.nIter
                 fprintf(1, repmat('\b',1,count)); %delete line before
                 count = fprintf('Remaining simulation time: %0.5g s',dt*(obj.sys.nIter-k)/k);
             elseif k==obj.sys.nIter
                 fprintf(1, repmat('\b',1,count)); %delete line before
                 count = fprintf('Remaining simulation time: %0.5g s\n',dt*(obj.sys.nIter-k)/k);
             end
            end
            close(hwait);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%---------- STORING TELEMETRY -------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Slopes            
            obj.slopesAliasOL= sAliasOL;
            obj.slopesAlias  = sAlias;
            obj.slopesNoise  = sNoise;               
            obj.slopesServo  = sServo;
            obj.slopesGs     = sGs;                        
                        
            % TT slopes
            if obj.sys.nTT
                obj.slopesTT      = sTT;                
                obj.slopesServoTT = sServoTT;
                obj.slopesNoiseTT = sNoiseTT;
                obj.slopesAliasTT = sAliasTT;
            end
                    
            % FITTING
            ps          = 1./obj.sys.tel.D;
            s2          = sum(obj.window(:).^2);
            obj.psdFit  = fftshift(obj.psdFit.*(1./obj.resolution^2/ps^2/obj.sys.nIter/s2) ...
                .*(obj.sys.ast(1).wavelength/obj.sys.src(1).wavelength)^2);
            
            % FILLING THE SYSTEM CLASS
            obj.sys.slopesGs = obj.slopesGs;
            obj.sys.coefsHO  = obj.coefs(:,1,:);
            if obj.sys.nTT
                obj.sys.coefsTT = obj.coefsTT(:,1,:);
                obj.slopesTT    = obj.slopesTT;
            end
             
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%---------- STORING PSF AND GETTING STATISTICS -------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            tmppsf = obj.sys.cam.frame;
            npix   = obj.sys.cam.resolution(1);
            ps     = obj.sys.cam.pixelScale*constants.radian2mas;
            t      = obj.sys.cam.exposureTime*obj.sys.tel.samplingTime;
            
            for i=1:obj.sys.nSrc
                xi        = 1 + (i-1)*npix;
                xf        = i*npix;
                psfi      = tmppsf(:,xi:xf);         
                SR        = 100*obj.sys.cam.strehl(i);
                tmp(i)    = psfStats(obj.sys.src(i),obj.sys.tel,psfi,ps,t,'flagMoffat',obj.sys.flagMoffat,'strehl',100*obj.sys.cam.strehl(i));
            end
            
            obj.psf     = tmp; 
            obj.sys.psf = tmp;

     
            %% %%%%%%%%%%%%%%%%%%%%%%%
            %%----------- WFE ERROR BREAKDOWN  -----------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            rmsNano    = @(x,l) 1e9*sqrt(x).*l/2/pi;
            init       = obj.sys.startDelay;
            
            wvl        = obj.sys.ast.wavelength;            
            wfeFit     = rmsNano(mean(varPerp(init:end)),wvl);            
            wfeAliasOL = rmsNano(mean(varAliasOL(init:end)),wvl);            
            wfeAlias   = rmsNano(mean(varAlias(init:end)),wvl);
            wfeServo   = rmsNano(mean(varServo(init:end)),wvl);
            wfeNoise   = rmsNano(mean(varNoise(init:end)),wvl);
            wfeAniso   = rmsNano(mean(varAniso(init:end,:),1),wvl.*ones(1,obj.sys.nSrc));
            wfeRes     = rmsNano(mean(varResWfe(init:end),1),wvl);
            if obj.sys.nTT
                wvlTT      = obj.ngs.wavelength; 
                wfeServoTT = rmsNano(mean(varServoTT(init:end)),wvlTT);
                wfeAliasTT = rmsNano(mean(varAliasTT(init:end)),wvlTT);
                wfeNoiseTT = rmsNano(mean(varNoiseTT(init:end)),wvlTT);
                wfeAnisoTT = rmsNano(mean(varAnisoTT(init:end),1),wvlTT.*ones(1,obj.nSrc));
                wfeResTT   = rmsNano(mean(varResWfeTT(init:end),1),wvlTT);
            else
                wfeServoTT = 0;
                wfeAliasTT = 0;
                wfeNoiseTT = 0;
                wfeAnisoTT = 0;
                wfeResTT   = 0;
            end
            
            % Residual wfe
            wfeTotal   = sqrt(wfeFit^2 + wfeServo^2 + wfeAlias^2 + ...
                wfeNoise^2 + wfeAniso.^2 + wfeServoTT^2 + wfeAliasTT^2 ... 
                + wfeNoiseTT^2 + wfeAnisoTT.^2);
                       
            wfeAOres  = rmsNano(mean(varRes(init:end,:),1),[obj.sys.src.wavelength]);
            obj.sys.varRes = wfeAOres;
            % Strehl ratios
            wfeSRth    = 100.*exp(-(2*pi.*wfeTotal.*1e-9./[obj.sys.src.wavelength]).^2);
            wfeSRres   = 100.*exp(-(2*pi.*wfeAOres.*1e-9./[obj.sys.src.wavelength]).^2);

            % Storaging
            obj.wfeSimu       = zeros(17,obj.sys.nSrc);            
            obj.wfeSimu(1,:)  = obj.psf.Strehl;
            obj.wfeSimu(2,:)  = wfeSRres;
            obj.wfeSimu(3,:)  = wfeSRth;
            
            obj.wfeSimu(4,:)  = wfeAOres;
            obj.wfeSimu(5,:)  = wfeTotal;
            
            obj.wfeSimu(6,:)  = wfeFit;
            obj.wfeSimu(7,:)  = wfeServo;       
            obj.wfeSimu(8,:)  = wfeAlias;
            obj.wfeSimu(9,:)  = wfeAliasOL;            
            obj.wfeSimu(10,:) = wfeNoise;    
            obj.wfeSimu(11,:) = wfeAniso;
            obj.wfeSimu(12,:) = wfeRes;
            
            obj.wfeSimu(13,:) = wfeServoTT;            
            obj.wfeSimu(14,:) = wfeNoiseTT;     
            obj.wfeSimu(15,:) = wfeAnisoTT;
            obj.wfeSimu(16,:) = wfeAliasTT;
            obj.wfeSimu(17,:) = wfeResTT;
            
        end

        
        %% ----------------- SIMULATION ERROR BREAKDOWN ----------------- %
        function displayWfeTable(obj,iSrc)
            if nargin < 2
                iSrc =1;
            end
            fprintf('-------------------------------\n');
            fprintf('Error Term\t\t[nm]\t\n');
            fprintf('-------------------------------\n');
            fprintf('Atmospheric Fitting\t%.4g\n',obj.wfeSimu(6,iSrc));
            fprintf('Telescope Fitting\t%.4g\n',0);
            fprintf('NCPA calibration\t%.4g\n',0);
            fprintf('-------------------------------\n');
            fprintf('DM bandwidth\t\t%.4g\n',obj.wfeSimu(7,iSrc));
            fprintf('DM measurements\t\t%.4g\n',hypot(obj.wfeSimu(8,iSrc),obj.wfeSimu(10,iSrc)));
            fprintf('Ang./Focal aniso.\t%.4g\n',obj.wfeSimu(11,iSrc));
            fprintf('LGS focus error\t\t%.4g\n',0);
            fprintf('LGS high order error\t%.4g\n',0);
            fprintf('-------------------------------\n');
            if obj.nTT
                fprintf('Tip-tilt bandwidth\t%.4g\n',obj.wfeSimu(13,iSrc));
                fprintf('Tip-tilt measurements\t%.4g\n',hypot(obj.wfeSimu(14,iSrc),obj.wfeSimu(16,iSrc)));
                fprintf('Tip-tilt anisoplanatism\t%g\n',obj.wfeSimu(15,iSrc));
                fprintf('-------------------------------\n');     
            end                   
            fprintf('Calibration errors\t%.4g\n',0);
            fprintf('Miscellaneous\t\t%.4g\n',0);
            fprintf('-------------------------------\n');    
            fprintf('Wavefront errors\t%s%.4g/%s%.4g\n','wfs:',obj.wfeSimu(4,iSrc),'budget:',obj.wfeSimu(5,iSrc));
            fprintf('%s-band Strehl\t\t%s%.3g/%s%.3g\n',obj.sys.src(iSrc).photometry.char,'cam:',obj.psf(iSrc).Strehl,'budget:',obj.wfeSimu(3,iSrc));
            fprintf('-------------------------------\n');    
            fprintf('GS star %s-band mag.\t%g\n',obj.sys.ast.photometry.char,obj.sys.ast.magnitude);
            fprintf('GS star separation\t%g %s\n',obj.sys.ast.zenith*constants.radian2arcsec,'["]');
            if obj.nTT
                fprintf('TT star %s-band mag.\t%g\n',obj.sys.ast_tt.photometry.char,obj.sys.ast_tt.magnitude);
                fprintf('TT star separation\t%g %s\n',obj.sys.ast_tt.zenith*constants.radian2arcsec,'["]');
            end
            fprintf('Zenith angle\t\t%g %s\n',obj.sys.tel.elevation*180/pi,'[deg]');
            fprintf('-------------------------------\n');    
            
        end
        function [wfe,coef] = slopesNgs2wfe(obj,s,wvl)
            obj.sys.wfs.slopes = s;
            z    = obj.zerNGS\obj.sys.wfs;
            coef = z.c;
            coef = coef - repmat(mean(coef,2),1,size(coef,2));
            cov  = coef*coef'./size(coef,2);                        
            wfe  = 1e9*sqrt(trace(cov))*wvl/2/pi;
        end
        
        function [wfe,coef] = slopesHO2wfe(obj,s,wvl)
            obj.sys.wfs.slopes = s;
            z    = obj.zerHO\obj.sys.wfs;
            coef = z.c;
            coef = coef - repmat(mean(coef,2),1,size(coef,2));
            cov  = coef*coef'./size(coef,2);                        
            wfe  = 1e9*sqrt(trace(cov))*wvl/2/pi;
        end
        
        function [wfe,coef] = slopesTT2wfe(obj,s,wvl)
            obj.sys.wfs_tt.slopes = s;
            z    = obj.zerTT\obj.sys.wfs_tt;
            coef = z.c;
            coef = coef - repmat(mean(coef,2),1,size(coef,2));
            cov  = coef*coef'./size(coef,2);                        
            wfe  = 1e9*sqrt(trace(cov))*wvl/2/pi;
        end
        
        function out = getNgsAoError(obj)
            
            % -------------------- LGS MEASUREMENTS --------------------- %
            % Total
            s           = obj.slopesGs;
            [wfeNGSmeas,cHO]  = obj.slopesNgs2wfe(s,obj.sys.ast.wavelength);
            % Bandwidth error
            s           = obj.slopesServo;
            wfeNGSservo = obj.slopesNgs2wfe(s,obj.sys.ast.wavelength);
            % Aliasing error
            s           = obj.slopesAlias;
            wfeNGSalias = obj.slopesNgs2wfe(s,obj.sys.ast.wavelength);
            % Noise error
            s           = obj.slopesNoise;
            wfeNGSnoise = obj.slopesNgs2wfe(s,obj.sys.ast.wavelength);
            % Calibration error
            wfeNGScalib = sqrt(wfeNGSmeas^2 - wfeNGSservo^2 - wfeNGSalias^2-...
                wfeNGSnoise^2);
            
            % -------------- FOCAL ANISOPLANATISM ERROR ----------------- %
            if any(obj.sys.ast.directionVector ~= obj.sys.src.directionVector)
                %  Zernike Covariance matrix without cone effect
                obj.sys.wfs.slopes = obj.slopesNGS;
                z     = obj.zerNGS\obj.sys.wfs;
                con   = z.c;
                con   = con - repmat(mean(con,2),1,size(con,2));
                COn   = con*con'./size(con,2);
                CHO   = cHO*cHO'/size(cHO,2);
                %  Zernike cross-covariance matrix
                Ccross = con*cHO'./size(con,2);
                
                % Aniso error
                CfocAni     = COn + CHO - 2.*Ccross;
                wfeAniso = 1e9*sqrt(trace(CfocAni))*obj.sys.ast.wavelength/2/pi;
            else
                wfeAniso = 0;
            end
            
            out = [wfeNGSmeas,wfeNGSservo,wfeNGSalias,wfeNGSnoise,wfeNGScalib,wfeAniso];
        end
        
        function out = getHighOrderError(obj)
            
            % -------------------- LGS MEASUREMENTS --------------------- %
            % Total
            s           = obj.slopesLGS;
            [wfeLGSmeas,cHO]  = obj.slopesHO2wfe(s,obj.sys.ast.wavelength);
            % Bandwidth error
            s           = obj.slopesLGSservo;
            wfeLGSservo = obj.slopesHO2wfe(s,obj.sys.ast.wavelength);
            % Aliasing error
            s           = obj.slopesLGSalias;
            wfeLGSalias = obj.slopesHO2wfe(s,obj.sys.ast.wavelength);
            % Noise error
            s           = obj.slopesLGSnoise;
            wfeLGSnoise = obj.slopesHO2wfe(s,obj.sys.ast.wavelength);
            % Calibration error
            wfeLGScalib = sqrt(wfeLGSmeas^2 - wfeLGSservo^2 - wfeLGSalias^2-...
                wfeLGSnoise^2);
            
            % -------------- FOCAL ANISOPLANATISM ERROR ----------------- %
            if any(obj.sys.ast.directionVector ~= obj.sys.src.directionVector)
                %  Zernike Covariance matrix without cone effect
                obj.sys.wfs.slopes = obj.slopesNGS;
                z     = obj.zerHO\obj.sys.wfs;
                con   = z.c;
                con   = con - repmat(mean(con,2),1,size(con,2));
                COn   = con*con'./size(con,2);
                CHO   = cHO*cHO'/size(cHO,2);
                %  Zernike cross-covariance matrix
                Ccross = con*cHO'./size(con,2);
                
                % Aniso error
                CfocAni     = COn + CHO - 2.*Ccross;
                wfeFocAniso = 1e9*sqrt(trace(CfocAni))*obj.sys.ast.wavelength/2/pi;
            else
                wfeFocAniso = 0;
            end
            
            out = [wfeLGSmeas,wfeLGSservo,wfeLGSalias,wfeLGSnoise,wfeLGScalib,wfeFocAniso];
        end
        
        function out = getTipTiltError(obj)
            
            
            % ------------------ TIP-TILT MEASUREMENTS ------------------ %
            % Total
            s           = obj.slopesTT;
            [wfeTTmeas,coff]  = obj.slopesTT2wfe(s,obj.sys.ast_tt.wavelength);
            % Bandwidth error
            s           = obj.slopesTTservo;
            wfeTTbw    = obj.slopesTT2wfe(s,obj.sys.ast_tt.wavelength);
            % Aliasing error
            s           = obj.slopesTTalias;
            wfeTTalias = obj.slopesTT2wfe(s,obj.sys.ast_tt.wavelength);
            % Noise error
            s           = obj.slopesTTnoise;
            wfeTTnoise = obj.slopesTT2wfe(s,obj.sys.ast_tt.wavelength);
            % Calibration error
            wfeTTcalib = sqrt(wfeTTmeas^2 - wfeTTbw^2 - wfeTTalias^2-...
                wfeTTnoise^2);
            
            % ------------------ ANISOKINETISM ERROR -------------------- %
            
            %  Zernike Covariance matrix on-axis
            obj.sys.wfs_tt.slopes = obj.slopesTTOn;
            z     = obj.zerTT\obj.sys.wfs_tt;
            con   = z.c;
            con   = con - repmat(mean(con,2),1,size(con,2));
            CttOn = con*con'./size(con,2);
            
            %  Zernike cross-covariance matrix
            CttOff   = coff*coff'/size(coff,2);
            CttCross = con*coff'./size(con,2);
            
            % Aniso error
            obj.CaniTT = CttOn + CttOff - 2.*CttCross;
            wfeTTAniso = 1e9*sqrt(trace(obj.CaniTT))*obj.sys.ast_tt.wavelength/2/pi;
            
            out = [wfeTTmeas,wfeTTbw,wfeTTalias,wfeTTnoise,wfeTTcalib,wfeTTAniso];
        end
                
        %% ------------------ FOURIER ERROR BREAKDOWN ------------------- %
        function out = getWFEbreakdown(obj,wvl)
            
            %% Fourier domain instantiation
            obj.fao = fourierModel(obj.sys,'nTimes',10,'resolution',201);                 
            fx            = obj.fao.fx;
            fy            = obj.fao.fy;
                      
            obj.fao.sys.atm.wavelength = wvl;
            rad2nm                 = wvl/2./pi*1e9;
            
            % Fitting error
            wfeFit = sqrt(trapz(obj.fao.fyExt(:,1),trapz(obj.fao.fxExt(1,:),fittingPSD(obj.fao,fx,fy))))...
                *rad2nm;
            % Aliasing error
            wfeAlias  = sqrt(trapz(fy(:,1),trapz(fx(1,:),aliasingPSD(obj.fao,fx,fy))))...
                *rad2nm;
            % Noise error
            wfeNoise  = sqrt(trapz(fy(:,1),trapz(fx(1,:),noisePSD(obj.fao,fx,fy))))...
                *rad2nm;
            % AnisoServo-lag error
            wfeAS = zeros(obj.sys.nSrc,1);
            for i=1:obj.nSrc
                psd      = obj.fao.anisoServoLagPSD(fx,fy,i);
                wfeAS(i) = sqrt(trapz(fy(:,1),trapz(fx(1,:),psd)))*rad2nm;
            end
            
            % Total
            wfeTotal  = sqrt(wfeFit^2 + wfeNoise^2 + wfeAlias^2 + wfeAS^2);
            %SR
            wfeSR     = 100.*exp(-(2*pi*wfeTotal*1e-9/wvl)^2);
            
            %Open-loop aliasing
            wfeAliasOL = sqrt(0.074*(obj.sys.pitch/obj.fao.sys.atm.r0)^(5/3))*rad2nm;
            
            % Servo-lag error
            psd      = obj.fao.servoLagPSD(fx,fy);
            wfeServo = sqrt(trapz(fy(:,1),trapz(fx(1,:),psd)))*rad2nm;
            
            %Concatenation of values
            out    = zeros(8,1);
            out(1) = wfeSR;
            out(2) = wfeTotal;
            out(3) = wfeFit;
            out(4) = wfeAS;
            out(5) = wfeAlias;
            out(6) = wfeNoise;
            out(7) = wfeAliasOL;
            out(8) = wfeServo;
        end
        
        %% ----------------------- DISPLAY TOOLS ------------------------ %
        
%         function displayWFEbreakdown(obj)
%             %Display properties
%             thick  = 0.2;
%             %Grabbing wfe breakdown
%             budgF = obj.wfeFourier(2:end);
%             budgS = obj.wfeSimu(2:end-1);
%             
%             figure
%             hold on
%             n = size(budgF,1);
%             x = linspace(1,2*n,n);
%             bar(x-0.25,budgS,thick,'blue');
%             bar(x+0.25,budgF,thick,'red');function out = zernikeResidualVariance(N,atm,tel)
%             hold off
%             
%             ylabel('WF error [nm]');
%             labs = {'sigma_tot','sigma_Fit','sigma_AnisoServo','sigma_Alias','sigma_Noise'...
%                 'sigma_AliasOL','sigma_Servo'};
%             for i=1:length(labs)
%                 m = max(budgF(i),budgS(i));
%                 text(x(i)-2*thick,m*1.01+10,texlabel(labs(i)));
%             end
%             title('WFE decomposition');
%             legend('Simulation','Fourier');
%             
%             figure
%             hold on
%             sr = [obj(iSrc).SR obj(iSrc).wfeSimu(1) obj(iSrc).wfeSimu(2) ...
%                 obj(iSrc).wfeFourier(1)];
%             n = size(sr,2);
%             x = linspace(1,n,n);
%             ylabel('Strehl ratio [%]');
%             hold on
%             bar(1,sr(1),thick,'r');
%             bar(2,sr(2),thick,'b');
%             bar(3,sr(3),thick,'g');
%             bar(4,sr(4),thick,'y');
%             
%             legend('Image','Residual phase','WFE budget','Fourier');
%             hold off;
            
%             fprintf('-----------------------------------------------------------\n');
%             fprintf('Source location x/y        ["] : %g/%g \n',obj(iSrc).x,obj(iSrc).y);
%             fprintf('Image SR                   [%s]: %g \n','%',obj(iSrc).SR);
%             fprintf('SR from the residual phase [%s]: %g \n','%',obj(iSrc).wfeSimu(1));
%             fprintf('SR from the sum of error   [%s]: %g \n','%',obj(iSrc).wfeSimu(2));
%             fprintf('SR from Fourier            [%s]: %g \n','%',obj(iSrc).wfeFourier(1));
%             fprintf('-----------------------------------------------------------\n');
            
%        end
        
    end
end

