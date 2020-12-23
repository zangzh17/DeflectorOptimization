function FOM = FigureOfMerit1D(level,OptParm)
% range of level: [0,OptParm.Geometry.Level-1], integers
nBot = OptParm.Geometry.Substrate;
nTop = OptParm.Geometry.Top;
nDevice = OptParm.Geometry.Device;
Period = OptParm.Geometry.Period;
NFourier = OptParm.Simulation.Fourier;
kParallelForward = nBot*sind(OptParm.Input.Theta);
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
AbsEff = zeros(NRobustness,NumPol,NumWave);
% 结构texture初始化
LayerTextures = cell(1,OptParm.Geometry.Level + 1);
LayerTextures{1} = {nTop};
LayerTextures{end} = {nBot};
for ii=1:OptParm.Geometry.Level-1
    index = (level(:)'>ii-1);
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
for robustIter = 1:NRobustness
    % 修改厚度
    DeviceProfile = {[0, ones(1,OptParm.Geometry.Level-1)*(OptParm.Geometry.Thickness+Start(robustIter))/(OptParm.Geometry.Level-1), 0],...
    [1, 2:OptParm.Geometry.Level, OptParm.Geometry.Level+1]};
    for polIter = 1:NumPol
        % prepare for res0
        ReticoloParm = res0(Polarizations(polIter));
        % % illumination from bottom
        % ReticoloParm.res3.sens=-1;
        % illumination from top
        ReticoloParm.res3.sens=1;
        for wIter = 1:NumWave
            % prepare for res1
            LayerResults = res1(Wavelengths(wIter),Period,LayerTextures,NFourier,kParallelForward,ReticoloParm);
            % prepare for res2
            DeviceResults = res2(LayerResults,DeviceProfile);
            % get DE of the required order
            TransmittedLight = DeviceResults.inc_bottom_transmitted;
            AbsEff(robustIter,polIter,wIter) = TransmittedLight.efficiency(TransmittedLight.order==OptParm.Optimization.Target);
        end
    end
end
FOM = OptParm.Optimization.Robustness.Weights*mean(mean(AbsEff,3),2);
end
