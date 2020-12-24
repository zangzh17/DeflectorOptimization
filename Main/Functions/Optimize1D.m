function Gradient = TopoGradient1D(x,OptParm)

%Extract common values for easier use
Wavelengths = OptParm.Input.Wavelength;
Period = OptParm.Geometry.Period;
nBot = OptParm.Geometry.Substrate;
nTop = OptParm.Geometry.Top;
nDevice = OptParm.Geometry.Device;

NumWave = length(Wavelengths); 
Wavelength0 = mean(Wavelengths);

%Compute incident k-vector
kParallelForward = nBot*sind(OptParm.Input.Theta);

% Compute total Fourier orders
NFourier = OptParm.Simulation.Fourier;
NLevel = OptParm.Geometry.Level;

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
% Grid points per wavelength 
[xGrid, ~, GridScale] = DefineGrid(OptParm.Simulation.Grid, Period, Wavelength0);
N = length(xGrid); %Number of x grid points

% If no starting point is given, generate a random starting point
if isempty(OptParm.Optimization.Start)
    DeviceIn = OptParm.Optimization.Start;
    DevicePattern = FineGrid(DeviceIn,Period,N/length(DeviceIn),0,1);
else
    DevicePattern = RandomStart(N,1,Period,...
        OptParm.Optimization.RandomStart,0,0);
end
DevicePattern = DevicePattern * (NLevel -1);
StartPattern = DevicePattern;

% Define device stack of Nz layers for simulation
Nz = ceil(OptParm.Geometry.Thickness/OptParm.Simulation.ZGrid);
DeviceProfile = {[0, ones(1,Nz)*OptParm.Simulation.ZGrid, 0],...
    [1, 2+(1:Nz), 2]};

% Generate binarization parameter B
BVector = GenerateBVector(MaxIterations, OptParm.Optimization.Binarize);
% Generate thresholding parameters for robustness
[ThresholdVectors, NRobustness] = GenerateThreshVectors(OptParm);

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
AbsoluteEfficiency = zeros(MaxIterations,NRobustness,NumPol,NumWave);
RelativeEfficiency = zeros(MaxIterations,NRobustness,NumPol,NumWave);
% Compute blur radii [in grid units]
BlurGridLarge = OptParm.Optimization.Filter.BlurRadiusLarge/GridScale;
BlurGrid = OptParm.Optimization.Filter.BlurRadius/GridScale;

%Main optimization loop
for iter = iterStart:MaxIterations
    tic;
    % initial blurring of radius R=BlurGridLarge
    FilteredPattern = DensityFilter2D(DevicePattern,BlurGridLarge);
    % enforce binarization
    BinaryPattern = LevelFilter(FilteredPattern,BVector(iter),0.5);
    
    GradientsAll = cell([NRobustness, 1]);
    DispPattern = cell([NRobustness, 1]);
    
    % Begin robustness loop
    % Can be changed to parfor as necessary
    parfor robustIter = 1:NRobustness
        % Second filter to model physical edge deviations
        FilteredPattern2 = GaussFilter2D(BinaryPattern,BlurGrid);
        FinalPattern = ThreshFilter(FilteredPattern2,BVector(iter),ThresholdVectors(robustIter, iter));
        DispPattern{robustIter} = FinalPattern;
        
        % Define textures for each layer
        LayerTextures = TexturesGen1D(PatternIn,NLevel,GridScale,Nz,nTop,nDevice);
        % Initialize empty field matrix
        FieldProductWeighted = zeros(NumPol,NumWave,OptParm.Simulation.ZGrid,N);
        
        % Begin polarization loop
        % Can be changed to parfor as necessary
        for polIter = 1:NumPol
            % Set simulation parameters in Reticolo
            ReticoloParm = SetReticoloParm(OptParm, Polarizations, polIter);
            for wIter = 1:NumWave
                % res1 computes the scattering matrices of each layer
                LayerResults = res1(Wavelengths(wIter),Period,LayerTextures,NFourier,kParallelForward,0,ReticoloParm);
                % res2 computes the scattering matrix of the full device stack
                DeviceResults = res2(LayerResults,DeviceProfile);
                if OptParm.Optimization.Target(polIter) == 0 && OptParm.Input.Theta == 0
                    FieldConvention = -1; % For normal output, opposite sign convention
                else
                    FieldConvention = 1;
                end
                if (Polarizations(polIter)==1) %For TE polarization
                    % Extract simulation results
                    TransmittedLight = DeviceResults.TEinc_bottom_transmitted;
                    
                    % Find appropriate target
                    TargetIndex = find((TransmittedLight.order(:,1)==OptParm.Optimization.Target(polIter))&(TransmittedLight.order(:,2)==0));
                    
                    % Store efficiencies
                    AbsEff = TransmittedLight.efficiency_TE(TargetIndex);
                    RelativeEfficiency(iter,robustIter,polIter,wIter) = TransmittedLight.efficiency_TE(TargetIndex)/sum(sum(TransmittedLight.efficiency));
                    
                    % Compute input field for adjoint simulation
                    AdjointIncidence = [0,FieldConvention*exp(1i*angle(conj(TransmittedLight.amplitude_TE(TargetIndex))))];
                    normalization = sqrt(2/3); % Normalize field for polarizations
                    
                    % res3 computes internal fields for each layer
                    [ForwardField,~,~] = res3(xGrid,yGrid,LayerResults,DeviceProfile,[0,1],ReticoloParm);
                    
                elseif (Polarizations(polIter)==-1) %For TM polarization
                    % Extract simulation results
                    TransmittedLight = DeviceResults.TMinc_bottom_transmitted;
                    
                    % Find appropriate target and store efficiencies
                    TargetIndex = find((TransmittedLight.order(:,1)==OptParm.Optimization.Target(polIter))&(TransmittedLight.order(:,2)==0));
                    AbsEff = TransmittedLight.efficiency_TM(TargetIndex);
                    RelativeEfficiency(iter,robustIter,polIter,wIter) = TransmittedLight.efficiency_TM(TargetIndex)/sum(sum(TransmittedLight.efficiency));
                    
                    % Compute input field for adjoint simulation
                    AdjointIncidence = [-FieldConvention*exp(1i*angle(conj(TransmittedLight.amplitude_TM(TargetIndex)))),0];
                    normalization = (3/2)*sqrt(2/3); % Normalize field for polarizations
                    
                    % res3 computes internal fields for each layer
                    [ForwardField,~,~] = res3(xGrid,yGrid,LayerResults,DeviceProfile,[1,0],ReticoloParm);
                end
                
                kParallelAdjoint = -nTop*TransmittedLight.K(TargetIndex,1); % Get appropriate adjoint k vector
                ReticoloParm.res3.sens = 1; % Reverse illumination direction
                
                % Recompute layer scattering matrices for reverse direction
                LayerResults = res1(Wavelengths(wIter),Period,LayerTextures,NFourier,kParallelAdjoint,0,ReticoloParm);
                % Compute adjoint internal field
                [AdjointField,~,RefractiveIndex] = res3(xGrid,yGrid,LayerResults,DeviceProfile,AdjointIncidence,ReticoloParm);
                
                % Begin to compute 3D gradient by overlap of forward and adjoint fields
                FieldProduct = ForwardField(:,:,:,1).*AdjointField(:,:,:,1) + ...
                    ForwardField(:,:,:,2).*AdjointField(:,:,:,2) + ForwardField(:,:,:,3).*AdjointField(:,:,:,3);
                
                % Weight field overlap by refractive index and efficiency
                FieldProductWeighted(polIter,wIter,:,:) = 0.5*normalization*RefractiveIndex.*FieldProduct.*(1-AbsEff);
                AbsoluteEfficiency(iter,robustIter,polIter,wIter) = AbsEff;
                
            end
        end
        
        % Compute raw gradient for each pattern averaged over polarization
        FieldAll = 2*squeeze(mean(mean(sum(FieldProductWeighted,1),2),));
        Gradient = real(-1i*FieldAll);
            

        
        
        
        
        
        
        
        
        
        
        
        
        


