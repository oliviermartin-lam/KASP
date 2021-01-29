classdef psfStats < handle
    
    properties (SetObservable=true)       
        % system
        pupil;
        D;
        % Image
        image;
        procImg;
        xSrc;
        ySrc;
        psfResolution;
        pixelScaleInMas;
        nyquistSampling;
        psfFieldOfView;
        cropResolution;
        wavelength;
        focalGrid;
        % Noise
        exposureTime;
        telArea;
        bg;
        ron;
        % PSF Estimates
        Strehl;
        StrehlAccuracy;
        FWHM;
        FWHMAccuracy;
        FluxInFWHM;
        Ellipticity;
        EnsquaredEnergy;                
        %Astrophysics estimates
        flux;
        nPhoton;
        magnitude;
        zeroPoint;        
        %Moffat
        MoffatImage;
        MoffatResidual;
        MoffatParam;
        MoffatStats;
        MoffatStd;
        flagMoffat;
        %Gaussian
        GaussianImage;
        GaussianResidual;
        GaussianParam;
        GaussianStats;
        GaussianStd;
        flagGaussian;
        fitDL;
        psfDL;
    end
    
    
    methods
        
        function obj = psfStats(image,pupil,D,zenith,azimuth,wavelength,psInMas,expTime,varargin)
            
            % ---------  Checking inputs
            inputs = inputParser;
            inputs.addRequired('image', @isnumeric );
            inputs.addRequired('pupil', @isnumeric );
            inputs.addRequired('D', @isnumeric );
            inputs.addRequired('zenith', @isnumeric );
            inputs.addRequired('azimuth', @isnumeric );
            inputs.addRequired('wavelength', @isnumeric );
            inputs.addRequired('psInMas', @isnumeric );      
            inputs.addRequired('expTime', @isnumeric );   
            inputs.addParameter('zeroPoint', 1e11, @isnumeric );
            inputs.addParameter('cropResolution', [], @isnumeric );                        
            inputs.addParameter('flagMoffat',false,@(x) isnumeric(x) || islogical(x));
            inputs.addParameter('flagGaussian',false,@(x) isnumeric(x)|| islogical(x));
            inputs.addParameter('strehl',0,@isnumeric);
            inputs.addParameter('fitDL',false,@islogical);
            inputs.parse(image,pupil,D,zenith,azimuth,wavelength,psInMas,expTime,varargin{:});


            % ------------ Instantiating the src class
            obj.pupil      = inputs.Results.pupil;
            obj.D          = inputs.Results.D;
            obj.image      = inputs.Results.image;
            obj.wavelength = wavelength;
            obj.xSrc       = tan(zenith)*cos(azimuth).*constants.radian2arcsec;
            obj.ySrc       = tan(zenith)*sin(azimuth).*constants.radian2arcsec;
            obj.zeroPoint  = inputs.Results.zeroPoint;
            % ------------ Instantiating the cam class
            
            obj.pixelScaleInMas = inputs.Results.psInMas;
            loverD              = constants.radian2mas*wavelength/D/2;
            obj.nyquistSampling = loverD/obj.pixelScaleInMas;
            obj.exposureTime    = inputs.Results.expTime;
            obj.telArea         = psfrTools.getAreaFromPupil(pupil,D);
            % ------------ Manage the image
            obj.flagMoffat      = inputs.Results.flagMoffat;
            obj.flagGaussian    = inputs.Results.flagGaussian;
            obj.psfResolution   = size(obj.image,1);
            obj.cropResolution  = inputs.Results.cropResolution;
                                           
            if obj.cropResolution
                n = obj.psfResolution/obj.cropResolution;
                if n > 1
                    tmp   = lamTools.cropImage(obj.image,size(obj.image)/n);
                    [a,b] = find(tmp == max(tmp(:)));
                    da    = floor(obj.cropResolution/2+1-a);
                    db    = floor(obj.cropResolution/2+1-b);
                    if da~=0
                        tmp = circshift(tmp,[da,0]);
                    end
                    if db~=0
                        tmp = circshift(tmp,[0,db]);
                    end
                    obj.image         = tmp;
                    obj.psfResolution = size(obj.image,1);
                end
            end
            
            obj.psfDL = psfrTools.telescopePsf(obj.pupil,obj.nyquistSampling*2);
            obj.psfDL = lamTools.crop(obj.psfDL,obj.psfResolution);
            obj.psfDL = obj.psfDL/sum(obj.psfDL(:));
            
            obj.psfFieldOfView = obj.psfResolution * obj.pixelScaleInMas;
            
            % --------------- Computing image statistics            
            obj.focalGrid  = psfrTools.getFocalGrid(obj.psfResolution,obj.pixelScaleInMas);
        
            % ---------- Strehl ratio
            obj.Strehl    = inputs.Results.strehl;
            if ~obj.Strehl
                [obj.Strehl] = psfrTools.getStrehl(obj.image,obj.pupil,obj.nyquistSampling);
                obj.Strehl = 100*obj.Strehl;
            end
            
            % Flux, background value and read-out noise
            [obj.flux,obj.bg,obj.ron,msk] = psfrTools.getFlux(obj.image);
            
            % Magnitude
            obj.zeroPoint = obj.zeroPoint; % In ph/s
            obj.nPhoton   = obj.flux/obj.telArea/obj.exposureTime;
            obj.magnitude = -2.5*log10(obj.nPhoton/obj.zeroPoint);
            
            mx                 = max(obj.image(:));
            n                  = sum(sum(~msk));
            obj.StrehlAccuracy = obj.Strehl*(obj.ron/mx + n*obj.ron/obj.flux);

            
            % ---------- Ensquared Energy
            nEE                 = floor(obj.psfFieldOfView/loverD/2); 
            obj.EnsquaredEnergy = zeros(1,nEE);
            for k=1:nEE
                obj.EnsquaredEnergy(k) = psfrTools.getEnsquaredEnergy(k,obj.image,obj.D,obj.nyquistSampling);
            end
            
            % --------------- Moffat fitting is desired
            if obj.flagMoffat
                % Focal plane grid
                obj.focalGrid  = psfrTools.getFocalGrid(obj.psfResolution,obj.pixelScaleInMas);            
                obj = obj.getMoffat();
                obj.Ellipticity = obj.MoffatStats(4);
                obj.FWHM        = max(obj.MoffatStats(2:3));       
                be              = obj.MoffatParam(4);
                da              = 2*sqrt(2^(1/be)-1)*max(obj.MoffatStd(2),obj.MoffatStd(3));
                db              = obj.MoffatStd(4)*obj.FWHM*2^(1/be-1)*log(2)/(be^2*(2^(1/be)-1));
                obj.FWHMAccuracy= hypot(da,db);
            elseif obj.flagGaussian
                % --------------- Gaussian fitting is desired    
                % Focal plane grid
                obj = obj.getGaussian();
                obj.Ellipticity = obj.GaussianStats(4);
                obj.FWHM        = max(obj.GaussianStats(2:3));
                obj.FWHMAccuracy= 2*sqrt(2*log(2))*hypot(obj.GaussianStd(2),obj.GaussianStd(3));
            else
                % ---------- FWHM
                [Fx,Fy,obj.Ellipticity] = psfrTools.getFWHM(obj.image,obj.pixelScaleInMas,8,'contour');
                obj.FWHM = max(Fx,Fy);
            end
                                   
        end
                     
        function obj = getMoffat(obj)
            
            psf                = obj.image/sum(obj.image(:));
            xdata              = obj.focalGrid;                   
            if obj.fitDL
                x0 = [1,5,5,2,0,0,0,0];     
                f  = @(x,xdata) max(psf(:))*psfrTools.convolve(psfrTools.moffatModel(x(1:end-1),xdata),obj.psfDL)+x(end);
                % how do we get the FWHM then ?
            else
                x0 = [max(psf(:)),2,2,1,0,0,0,0];     
                f  = @(x,xdata) psfrTools.moffatModel(x(1:end-1),xdata)+x(end);
            end
            
            %[beta,R,J]         = nlinfit(xdata,psf(:),f,x0);
            [beta,~,R,~,~,~,J] = lsqcurvefit(f,x0,xdata,psf);
            obj.MoffatStd      = diff(nlparci(beta,R,'jacobian',J),1,2);
            obj.MoffatParam    = beta;
            obj.MoffatImage    = f(obj.MoffatParam,obj.focalGrid);
            obj.MoffatResidual = psf - obj.MoffatImage;
            % Getting FWHM and Ellipticity
            obj.MoffatStats    = zeros(4,1);
            ax                 = abs(obj.MoffatParam(2));    % X Spreading factor
            ay                 = abs(obj.MoffatParam(3));    % Y Spreading factor
            b                  = obj.MoffatParam(4);    % slope or Ellipticity
            obj.MoffatStats(1) = 100*obj.MoffatParam(1);            
            obj.MoffatStats(2) = 2*ax*sqrt(2^(1./b)-1);
            obj.MoffatStats(3) = 2*ay*sqrt(2^(1./b)-1);
            obj.MoffatStats(4) = max([ax./ay ay./ax]);                       
            idx                = obj.MoffatImage(:) >= (max(obj.MoffatImage(:))/2.);
            obj.MoffatStats(5) = 100*sum(obj.MoffatImage(idx))/sum(obj.MoffatImage(:)); %% Ensquared energy in FWHM diameter
        end
        
        function obj = getGaussian(obj)
            
            psf                  = obj.image;
            xdata                = obj.focalGrid;            
            if obj.fitDL
                x0 = [1,10,10,0,0,0,0];    
                f = @(x,xdata) max(psf(:))*psfrTools.convolve(psfrTools.gaussianModel(x,xdata),obj.psfDL)+x(end);
            else
                x0 = [max(psf(:)),5,5,0,0,0,0];    
                f  = @(x,xdata) psfrTools.gaussianModel(x,xdata)+x(end);
            end
            
            %[beta,R,~,C]         = nlinfit(xdata,psf(:),f,x0);
            [beta,~,R,~,~,~,J]   = lsqcurvefit(f,x0,xdata,psf);
            obj.GaussianStd      = diff(nlparci(beta,R,'jacobian',J),1,2);
            obj.GaussianParam    = beta;
            obj.GaussianImage    = f(obj.GaussianParam,obj.focalGrid);
            obj.GaussianResidual = psf - obj.GaussianImage;
            % Getting FWHM and Ellipticity
            obj.GaussianStats    = zeros(4,1);
            ax                   = abs(obj.GaussianParam(2));    % X Spreading factor
            ay                   = abs(obj.GaussianParam(3));    % Y Spreading factor
            obj.GaussianStats(1) = 100*obj.GaussianParam(1);
            obj.GaussianStats(2) = 2*ax*sqrt(2*log(2));
            obj.GaussianStats(3) = 2*ay*sqrt(2*log(2));
            obj.GaussianStats(4) = max([ax./ay ay./ax]);
            
            idx                  = obj.GaussianImage(:) >= (max(obj.GaussianImage(:))/2.);
            obj.GaussianStats(5) = 100*sum(obj.GaussianImage(idx))/sum(obj.GaussianImage(:));
            obj.Ellipticity      = max([ax./ay ay./ax]);
        end
        
        
    end
    
end
