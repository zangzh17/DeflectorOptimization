function Gradient = TopoGradient1D(x,OptParm)

%Extract common values for easier use
Wavelength = OptParm.Input.Wavelength;
Period = OptParm.Geometry.Period;
nBot = OptParm.Geometry.Substrate;
nTop = OptParm.Geometry.Top;
nDevice = OptParm.Geometry.Device;

%Compute incident k-vector
kParallelForward = nBot*sind(OptParm.Input.Theta);

% Compute total Fourier orders
NFourier = OptParm.Simulation.Fourier;

% Define polarization values
if strcmp(OptParm.Input.Polarization,'TE') 
    Polarizations = 1;
elseif strcmp(OptParm.Input.Polarization,'TM')
    Polarizations = -1;
elseif strcmp(OptParm.Input.Polarization,'Both')
    Polarizations = [1, -1];
else
    error('Invalid polarization');
end
NumPol = length(Polarizations); 

% If no starting point is given, generate a random starting point
if isempty(OptParm.Optimization.Start)
    DevicePattern = OptParm.Optimization.RandomStart;
else
    DevicePattern = RandomStart1D(OptParm.Geometry.PixelNum, OptParm.Geometry.Level, OptParm.Optimization.RandomStart);
end
StartPattern = DevicePattern;


