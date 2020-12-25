function Gradient = Optimize1D(OptParm)

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
[xGrid, ~, xGridScale] = DefineGrid(OptParm.Simulation.Grid, Period, Wavelength0);
Nx = length(xGrid); %Number of x grid points

% If no starting point is given, generate a random starting point
if isempty(OptParm.Optimization.Start)
    DeviceIn = OptParm.Optimization.Start;
    DevicePattern = FineGrid(DeviceIn,Period,Nx/length(DeviceIn),0,1);
else
    DevicePattern = RandomStart(Nx,1,Period,...
        OptParm.Optimization.RandomStart,0,0);
end
% pattern height in unit of levels
% max:(NLevel -1)
% min: 0
% allowing NLevel-1 layers, NLevel levels
DevicePattern = DevicePattern * (NLevel -1);
StartPattern = DevicePattern;

% Define device stack of Nz layers for simulation
% OptParm.Simulation.ZGrid: number of layers per level
Nz = ceil(NLevel*OptParm.Simulation.ZGrid);
% thickness of one layer in nm
zGridScale = OptParm.Geometry.Thickness/Nz;
% thickness of one layer in level
zLevelScale = (NLevel -1)/Nz;
DeviceProfile = {[0, ones(1,Nz)*zGridScale, 0],...
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
BlurGridLarge = OptParm.Optimization.Filter.BlurRadiusLarge/xGridScale;
BlurGrid = OptParm.Optimization.Filter.BlurRadius/xGridScale;

%Main optimization loop
for iter = iterStart:MaxIterations
    tic;
    % [First filter] to enforce discretization
    FilteredPattern = DensityFilter2D(DevicePattern,BlurGridLarge);
    BinaryPattern = LevelFilter(FilteredPattern,BVector(iter),0.5);
    
    GradientsAll = cell([NRobustness, 1]);
    DispPattern = cell([NRobustness, 1]);
    
    % Begin robustness loop
    % Can be changed to parfor as necessary
    parfor robustIter = 1:NRobustness
        % [Second filter] to model physical edge deviations
        FilteredPattern2 = GaussFilter2D(BinaryPattern,BlurGrid);
        FinalPattern = LevelFilter(FilteredPattern2,BVector(iter),ThresholdVectors(robustIter, iter));
        DispPattern{robustIter} = FinalPattern;
        
        % Define textures for each layer
        LayerTextures = TexturesGen1D(FinalPattern,NLevel,xGridScale,Nz,nTop,nDevice);
        % Initialize empty field matrix
        FieldProductWeighted = zeros(NumPol,NumWave,Nz,Nx);
        
        % Set simulation parameters in Reticolo
        ReticoloParm = res0;
        ReticoloParm.res3.npts = [0, ones(1,Nz), 0]; %Number of points to sample field in each layer
        ReticoloParm.res3.sens = -1; % Default to illumination from below
        ReticoloParm.res1.champ = 1; % Accurate fields
        
        % Begin polarization loop
        % Can be changed to parfor as necessary
        for polIter = 1:NumPol
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
                    % the 1st index of fields refer to the z-coordinate
                    % the 2nd index of fields refer to the x-coordinate and the 3rd to the ycoordinate.
                    % the 4th index (1-6) of fields refer to [Ex,Ey,Ez,Hx,Hy,Hz]
                    [ForwardField,~,~] = res3(xGrid,yGrid,LayerResults,DeviceProfile,[1,0],ReticoloParm);
                end
                
                kParallelAdjoint = -nTop*TransmittedLight.K(TargetIndex,1); % Get appropriate adjoint k vector
                ReticoloParm.res3.sens = 1; % Reverse illumination direction
                
                % Recompute layer scattering matrices for reverse direction
                LayerResults = res1(Wavelengths(wIter),Period,LayerTextures,NFourier,kParallelAdjoint,0,ReticoloParm);
                % Compute adjoint internal field
                [AdjointField,~,RefractiveIndex] = res3(xGrid,yGrid,LayerResults,DeviceProfile,AdjointIncidence,ReticoloParm);
                
                % Begin to compute 3D gradient by overlap of forward and adjoint fields
                % Prod = Ex.*Eax + Ey.*Eay + Ez.*Eaz
                FieldProduct = ForwardField(:,:,:,1).*AdjointField(:,:,:,1) + ...
                    ForwardField(:,:,:,2).*AdjointField(:,:,:,2) + ForwardField(:,:,:,3).*AdjointField(:,:,:,3);
                
                % Weight field overlap by refractive index and efficiency
                FieldProductWeighted(polIter,wIter,:,:) = 0.5*normalization*RefractiveIndex.*FieldProduct.*(1-AbsEff);
                AbsoluteEfficiency(iter,robustIter,polIter,wIter) = AbsEff;
                
            end
        end
        
        % Compute raw gradient for each pattern averaged over pol/wave
        FieldAll = 2*squeeze(mean(mean(FieldProductWeighted,1),2));
        
        % height in number of layers ([0,Nz])
        HeightLayers = FinalPattern/zLevelScale;
        % index of layers ([1,Nz])
        LayerSelect = ceil(HeightLayers);
        LayerIndex = sub2ind(size(FieldAll),LayerSelect,1:Nx);
        % choose surface field
        FieldAll = FieldAll(LayerIndex);
        Gradient = real(-1i*FieldAll);
        
        % Back propagate gradient through robustness filters
        % Second deviations filter
        Gradient = FilteredGrad2D(Gradient,FilteredPattern2,BVector(iter),ThresholdVectors(robustIter, iter),BlurGrid);
        % First blur & binarize filter
        Gradient = FilteredGrad2D(Gradient,FilteredPattern,BVector(iter),0.5,BlurGridLarge);
        GradientsAll{robustIter} = Gradient;
    end
    
    % Sum gradient over all robustness variants
    Gradients=zeros(size(DevicePattern));
    for robustIter = 1:NRobustness
        Gradients= Gradients + GradientsAll{robustIter};
    end
    Gradient = Gradients;
    
    % Normalize gradient to step size
    CurrStepSize = OptParm.Optimization.Gradient.StepSize*OptParm.Optimization.Gradient.StepDecline^iter;
    Gradient = CurrStepSize*Gradient/max(max(abs(Gradient))); 
    
    % Remove unusable terms from normalization
    Gradient((Gradient+DevicePattern)>NLevel-1) = NLevel-1 -DevicePattern((Gradient+DevicePattern)>NLevel-1);
    Gradient((DevicePattern+Gradient)<0) = -DevicePattern((DevicePattern+Gradient)<0);
    Gradient = CurrStepSize*Gradient/max(max(abs(Gradient))); % Re-normalize gradient
    
    % Add gradient to device
    DevicePattern = DevicePattern + Gradient;
    
    % Ensure valid geometry       
    DevicePattern(DevicePattern<0) = 0;                                              
    DevicePattern(DevicePattern>NLevel-1) = NLevel-1;   
    
    % Apply blur if necessary
    DevicePattern = BlurGeomPostGrad(DevicePattern, iter, OptParm, xGridScale);
    
    % Show Progress in Figs
    [~,RobustInd] = min(abs(OptParm.Optimization.Robustness.EndDeviation));
    ShowProgress(OptParm, {xGrid,yGrid}, DispPattern{RobustInd}, iter, MaxIterations, AbsoluteEfficiency, RelativeEfficiency, Figs)
    toc;
end
% Compute different pattern varients
FilteredPattern = DensityFilter2D(DevicePattern,BlurGridLarge);
BinaryPattern = LevelFilter(FilteredPattern,BVector(iter),0.5);
FilteredPattern2 = GaussFilter2D(BinaryPattern,BlurGrid);
FinalPattern = LevelFilter(FilteredPattern2,BVector(iter),0.5);

% Save outputs
OptOut.BasePattern = DevicePattern;
OptOut.BinaryPattern = FilteredPattern;
OptOut.FinalPattern = FinalPattern;
OptOut.StartPattern = StartPattern;
OptOut.AbsoluteEfficiency = AbsoluteEfficiency;
OptOut.RelativeEfficiency = RelativeEfficiency;

end
        
        
        
        
        


