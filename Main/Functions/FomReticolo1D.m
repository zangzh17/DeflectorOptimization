function FOM = FomReticolo1D(height,OptParm)
% range of height: [0,OptParm.Geometry.Level-1], integers
nBot = OptParm.Geometry.Substrate;
nTop = OptParm.Geometry.Top;
nDevice = OptParm.Geometry.Device;
Period = OptParm.Geometry.Period;
NFourier = OptParm.Simulation.Fourier;
kParallelForward = nBot*sind(OptParm.Input.Theta);
Direction = OptParm.Input.Direction;
% 偏振
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
% 波长
Wavelengths = OptParm.Input.Wavelength;
NumWave = length(Wavelengths); 
% 容差
Start = OptParm.Optimization.Robustness.StartDeviation;
NRobustness = length(Start);
% 结构texture初始化
LayerTextures = cell(1,OptParm.Geometry.Level + 1);
LayerTextures{1} = {nTop};
LayerTextures{end} = {nBot};
for ii=1:OptParm.Geometry.Level-1
    index = (height(:)'>ii-1);
    index = [index,~index(end)];
    index = diff(index);
    lowIndexPos = find(index==1)*OptParm.Geometry.Pixel;
    highIndexPos = find(index==-1)*OptParm.Geometry.Pixel;
    index = [ones(size(lowIndexPos))*nTop, ones(size(highIndexPos))*nDevice];
    pos = [lowIndexPos,highIndexPos];
    [pos, I] = sort(pos); index = index(I);
    LayerTextures{ii+1} = {pos, index};
end

% begin loop
[RobustnessInd,PolInd,WaveInd] = ndgrid(1:NRobustness,1:NumPol,1:NumWave);
robustIter = RobustnessInd(:);
polIter = PolInd(:);
wIter = WaveInd(:);
AbsEff = zeros(1,length(robustIter));
for ii = 1:length(robustIter)
    % 修改厚度
    DeviceProfile = {[0, ones(1,OptParm.Geometry.Level-1)*(OptParm.Geometry.Thickness*(1+Start(robustIter(ii))))/(OptParm.Geometry.Level-1), 0],...
        [1, 2:OptParm.Geometry.Level, OptParm.Geometry.Level+1]};
    % prepare for res0
    ReticoloParm = res0(Polarizations(polIter(ii)));
    if Direction == 1
        % illumination from top
        ReticoloParm.res3.sens=1;
    elseif Direction == -1
        % illumination from bottom
        ReticoloParm.res3.sens=-1;
    end
    % prepare for res1
    LayerResults = res1(Wavelengths(wIter(ii)),Period,LayerTextures,NFourier,kParallelForward,ReticoloParm);
    % prepare for res2
    DeviceResults = res2(LayerResults,DeviceProfile);
    % get DE of the required order
    TransmittedLight = DeviceResults.inc_bottom_transmitted;
    AbsEff(ii) = TransmittedLight.efficiency(TransmittedLight.order==OptParm.Optimization.Target);
end
AbsEff = reshape(AbsEff,size(RobustnessInd));
Weights = OptParm.Optimization.Robustness.Weights/sum(OptParm.Optimization.Robustness.Weights);
FOM = Weights*mean(mean(AbsEff,3),2);
end
