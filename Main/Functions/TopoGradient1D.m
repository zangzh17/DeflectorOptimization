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

% Define grid for the device
[xGrid, ~, GridScale] = DefineGrid(OptParm.Simulation.Grid, Period, Wavelength);
N = length(xGrid); %Number of x grid points

% If no starting point is given, generate a random starting point
if isempty(OptParm.Optimization.Start)
    DevicePattern = OptParm.Optimization.RandomStart;
else
    DevicePattern = RandomStart(N,1,Period,...
        OptParm.Optimization.RandomStart,0,0);
end
StartPattern = DevicePattern * (OptParm.Geometry.Level -1);

% Define full device stack
Nz = ceil(OptParm.Geometry.Thickness/OptParm.Simulation.ZGrid);
DeviceProfile = {[0, ones(1,Nz)*OptParm.Simulation.ZGrid, 0],...
    [1, Nz, Nz+1]};

% Generate binarization parameter B
BVector = GenerateBVector(MaxIterations, OptParm.Optimization.Binarize);
% Generate thresholding parameters for robustness
[ThresholdVectors, NRobustness] = GenerateThreshVectors(OptParm);

% Compute blur radii in grid units
BlurGridLarge = OptParm.Optimization.Filter.BlurRadiusLarge/GridScale;
BlurGrid = OptParm.Optimization.Filter.BlurRadius/GridScale;
Figs = [];
% Initializing plot for geometries
if OptParm.Display.PlotGeometry 
    Figs.FigGeo = OptParm.Display.FigGeo;
end

% Initializing plot for geometries
if OptParm.Display.PlotEfficiency
    Figs.FigEff = OptParm.Display.FigEff;
end


% Store efficiency at each iteration
AbsoluteEfficiency = zeros(MaxIterations,NRobustness,NumPol);
RelativeEfficiency = zeros(MaxIterations,NRobustness,NumPol);

%Main optimization loop
for iter = iterStart:MaxIterations
    tic;
    % First filter to enforce binarization
    
    
    % Begin robustness loop
    % Can be changed to parfor as necessary
    for robustIter = 1:NRobustness


