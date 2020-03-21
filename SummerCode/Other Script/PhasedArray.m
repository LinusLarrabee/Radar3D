%% Design Frequency and Array parameters
% Choose the design frequency to be 1.8 GHz.
freq = 1.8e9;
c = physconst('lightspeed');
lambda = c/freq;
N = 9;
dx= 0.49*lambda;

%% Create Resonant Dipole
dipole_L = lambda/2;
dipole_W = lambda/200;
mydipole = dipole;
mydipole.Length = dipole_L;
mydipole.Width = dipole_W;
mydipole.TiltAxis = 'Z';
mydipole.Tilt = 90;
fmin = freq - .05*freq;
fmax = freq + .05*freq;
minX = 0.0001;          % Minimum value of reactance to achieve
trim = 0.0005;          % The amount to shorten the length at each iteration
resonant_dipole = dipole_tuner(mydipole,freq,fmin,fmax,minX,trim);
Z_resonant_dipole = impedance(resonant_dipole,freq)


%% Create Linear Array
dipole_array = linearArray;
dipole_array.Element = resonant_dipole;
dipole_array.NumElements = N;
dipole_array.ElementSpacing = dx;
%% Plot 2D Radiation Pattern
patternazfig1 = figure;
az_angle = 1:0.25:180;
pattern(dipole_array,freq,az_angle,0,'CoordinateSystem','rectangular')
axis([0 180 -25 15])

%% Beam Scanning
scanangle = [135 0];
ps = phaseshift(freq,dx,scanangle,N);
dipole_array.PhaseShift = ps;
patternazfig2 = figure;
pattern(dipole_array,freq,az_angle,0,'CoordinateSystem','rectangular')
axis([0 180 -25 15])
