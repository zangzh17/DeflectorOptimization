function FOM = FomFdtd1D(h,InPattern,OptParm,fig_handle)
% h: fdtd_handle
% range of height: [0,OptParm.Geometry.Level-1], integers
%% basic setting
Period = OptParm.Geometry.Period;
NLevel = OptParm.Geometry.Level;
Theta = OptParm.Input.Theta;
% 入射角
appputvar(h,'Theta',Theta);
code=strcat('switchtolayout;',...
    'setnamed("source","angle theta",Theta);');
if Theta~=0
    code=strcat(code,'setnamed("source","plane wave type","BFAST");');
else
    code=strcat(code,'setnamed("source","plane wave type","Bloch");');
end
appevalscript(h,code);
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
% 容器
AbsEff = zeros(NRobustness,NumPol,NumWave);
% 细网格初始化
[xGrid,~, xGridScale] = DefineGrid(OptParm.Simulation.Grid, Period, Wavelength0);
Nx = length(xGrid); %Number of x grid points
FinePattern = FineGrid(InPattern,Period,Nx/length(InPattern),0);
% number of layers
Nz = NLevel -1;
% thickness of one layer in nm
zGridScale = OptParm.Geometry.Thickness/Nz;
% robustness parameters
Start = OptParm.Optimization.Robustness.StartDeviation;
BlurGrid = OptParm.Optimization.Filter.BlurRadius/xGridScale;
ThresholdVector = 0.5*(1-erf(Start/OptParm.Optimization.Filter.BlurRadius));
BValue = OptParm.Optimization.Binarize.Max;
% begin loop
% filter to model physical edge deviations
FilteredPattern = GaussFilter2D(FinePattern,BlurGrid);
for robustIter = 1:NRobustness
    FinalPattern = LevelFilter(FilteredPattern,BValue,ThresholdVector(robustIter));
    % plot
    if nargin>3
        set(groot,'CurrentFigure',fig_handle);
        plot(xGrid,FinalPattern); hold on
    end
    % FDTD setting & structure
    appputvar(h,'Height',FinalPattern*1e-9);
    code=strcat('switchtolayout;',...
        'select("model");',...
        'set("Height",Height);',...
        'set("Height0",0);',...
        'set("Direction",1);');
    appevalscript(h,code);
    for polIter = 1:NumPol
        % set polarization
        if (Polarizations(polIter)==1) %For TE polarization
            code='setnamed("source","polarization angle",90);';
        elseif (Polarizations(polIter)==-1) %For TM polarization
            code='setnamed("source","polarization angle",0);';
        end
        % set wavelength
        appputvar(h,'F_points',length(Wavelengths));
        appputvar(h,'F1',Wavelengths(1));
        appputvar(h,'F2',Wavelengths(end));
        code=strcat(code,'setglobalmonitor("frequency points",F_points);');
        code=strcat(code,'setnamed("source","wavelength start",F1);');
        code=strcat(code,'setnamed("source","wavelength stop",F2);');
        % run simulation
        code=strcat(code,'run;');
        % Extract simulation results
        code=strcat(code,'select("grating_T");',...
            'runanalysis("grating_T");',...
            'Ord=getdata("grating_T","n");',...
            'Eff=getdata("grating_T","T_grating");');
        appevalscript(h,code);
        Orders = appgetvar(h,'Ord');
        Eff = appgetvar(h,'Eff');
        TargetIndex = (Orders==OptParm.Optimization.Target);
        for waveIter = 1:NumWave
            AbsEff(robustIter,polIter,waveIter) = Eff(TargetIndex,waveIter);
        end
    end
end
if nargin>3
    set(groot,'CurrentFigure',fig_handle);
    xlabel('Position/nm')
    ylabel('height/level')
end
FOM = OptParm.Optimization.Robustness.Weights*mean(mean(AbsEff,3),2);
end
