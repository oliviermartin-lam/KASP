parm = [];

%% ATM
parm.atm.zenithAngle            = pi/6;                         % rad
airmass                         = 1/cos(parm.atm.zenithAngle);
parm.atm.photometry             = photometry.V0;                % photometric band
parm.atm.windSpeed              = [10,25,15];                   % wind speed vector in m/s
parm.atm.windDirection          = [0,pi/3,pi/2];                % wind direction vector in rad
parm.atm.r0                     = 0.12/airmass^(3/5);           % coherence lenght in meters at 0.5microns
parm.atm.L0                     = 25;                           % Outer scale in meters
parm.atm.altitude               = [0,5,10]*1e3*airmass;         % altitude vector in meters
parm.atm.fractionalR0           = [.5,.3,.2];                   % fractional Cn2 profile (sum=1)

%% SCIENCE                                              
parm.sci.magnitude              = 0;                            % stellar magnitude
parm.sci.photometry             = photometry.K;                 % photometric band
parm.sci.height                 = Inf;                          % height (default Inf)
parm.sci.x                      = 0;                            % Sources x cartesian position in arcsec
parm.sci.y                      = 0;                            % Sources x cartesian position in arcsec

%% TEL
parm.tel.D                      = 11.25;                        % telescope primary mirror diameter
parm.tel.obstructionRatio       = 0.2356;                       % central obscuration ratio
parm.tel.fieldOfViewInArcsec    = 10;                           % fieldOfViewInArcsec
parm.tel.resolution             = 200;                          % resolution 
parm.tel.samplingTime           = 1e-3;                         % sampling time in sec
parm.tel.pupilShape             = [];                           % User-defined pupil-shape

%% SEGMENTED PUPIL
parm.seg.usePupilClass          = true;
parm.seg.nSides                 = 6;                            % Number of segment sides
parm.seg.radius                 = 0.9;                          % Segment radius
parm.seg.resolution             = parm.tel.resolution;         % Segment resolution in pixel
xSeg                            = parm.seg.radius*[1.5 0 -1.5 -1.5 0 1.5 3 1.5 0 -1.5 -3 -3 -3 -1.5 0  1.5 3 3 4.5 3 1.5 0 -1.5 -3 -4.5 -4.5 -4.5 -4.5 -3 -1.5 0 1.5 3 4.5 4.5 4.5];
ySeg                            = sqrt(3)/2*parm.seg.radius*[1 2 1 -1 -2 -1  2  3  4   3   2  0 -2  -3 -4   -3 -2 0 3 4 5 6 5 4 3 1 -1 -3 -4 -5 -6 -5 -4 -3 -1 1];
parm.seg.segmentCoordinates     = [xSeg' ySeg'];                % Segment coordinates
parm.seg.pupilUnit              = 'm';                          % Units for coordinates: pixel | meter
parm.seg.pupilSpider.n          = 3;                            % Number of spiders
parm.seg.pupilSpider.angle      = [0 pi/3 2*pi/3];              % Spiders angular positions
parm.seg.pupilSpider.width      = 0.1;                        % Spider width in pupilUnit
parm.seg.pupilSpider.usrDefined = 0;
parm.seg.pupilSpider.usrMask    = [];

diffPiston                      = 0*100e-9;                               % Statistical diff piston in m;
rampPiston                      = 0*80e-9;                                % Piston ramp max value in m                  
piston                          = diffPiston*randn(1,length(parm.seg.segmentCoordinates));
piston                          = piston-mean(piston);
ramp                            = keckTools.keck_pistonRamp(rampPiston);
parm.seg.pupilCoeffPhaseModes   = piston' + ramp;               % Vector of piston errors
parm.seg.pupilReflexivity       = [];                           % Vector of transmission disparities

%% NGS SOURCES
parm.nGs.magnitude              = 11;                           % magnitude                         
parm.nGs.photometry             = photometry.R;                 % photometric band                    
parm.nGs.x                      = 0;                            % NGSs x cartesian position in arcsec
parm.nGs.y                      = 0;                            % NGSs y cartesian position in arcsec                   

%% WAVE-FRONT SENSOR
parm.wfs.nLenslet               = 20;                           % number of lenslets
parm.wfs.nPx                    = 10;                           % number of pixels per lenslet
parm.wfs.pixelScale             = parm.wfs.nLenslet*constants.radian2mas*parm.nGs.photometry.wavelength/parm.tel.D/2;          % Pixel scale in mas          
parm.wfs.minLightRatio          = 0.5;                          % [0-1] partial illumination ratio to consider valid a sub-aperture
parm.wfs.exposureTime           = 1e-3;                         % WFS exposure time
parm.wfs.frameRate              = 1;                            % slaved to tel.samplingTime
parm.wfs.ron                    = 0;                            % ron in e-
parm.wfs.photonNoise            = false;                         % true/false statement
parm.wfs.pixelThreshold         = 0.5;                          % Pixel threshold in the cog algorithm
parm.wfs.throughput             = .6;                           % [0-1] throughput ratio 

%% DEFORMABLE MIRROR
parm.dm.pitch                   = 0.5625;                          % dm pitch
parm.dm.nActuators              = 21;                           % Number of actuators par row
parm.dm.crossCoupling           = 0.2;                          % actuator cross-coupling
parm.dm.influenceType           = 'gaussian';                   % type of influence function: gaussian | monotonic | overshoot | other

%% IMAGING CAMERA
parm.cam.pixelScale             = 9.94;                         % Pixel scale in mas                         
parm.cam.resolution             = 200;                          % camera resolution                               
parm.cam.clockRate              = 1;                            % clock-rate (default 1)                   
parm.cam.exposureTime           = 1000;                         % exposure time in time-steps x clockRate                   
parm.cam.startDelay             = 50;                           % transient period for loop to converge in time-steps x clockRate                   
parm.cam.ron                    = 0;                            % ron in e-                  
parm.cam.photonNoise            = false;                        % true/false photon noise   

%% LOOPS SETTINGS
parm.loopStatus.ho.nModeFiltered= 5;                            % # modes to truncate during SVD
parm.loopStatus.ho.controlBasis = [];                           % controlBasis: classic (default) | modal
parm.loopStatus.ho.gain         = 0.5;                          % integrator control loop gain
parm.loopStatus.ho.latency      = 0e-3;                         % loop latency (default tel.samplingTime)
parm.loopStatus.ho.tomography   = [];                           % tomographic algorithm
parm.loopStatus.ho.control      = 'integrator';                 % controller
parm.loopStatus.ho.gainPolc     = 1;                            % POLC control loop gain
%% NCPA
parm.ncpa                       = [];                           % NCPA map in meter




