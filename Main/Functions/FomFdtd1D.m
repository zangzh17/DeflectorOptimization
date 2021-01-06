function FOM = FomFdtd1D(InPattern,OptParm,sim_file_name)
% h: fdtd_handle
% range of height: [0,OptParm.Geometry.Level-1], integers
%% Open FDTD session

h=appopen('fdtd');
sim_file_path = getAttachedFilesFolder(sim_file_name);

%Pass the path variables to FDTD
appputvar(h,'sim_file_path',sim_file_path);
appputvar(h,'sim_file_name',sim_file_name);
code='cd(sim_file_path);load(sim_file_name);';
appevalscript(h,code);
%% basic setting
Period = OptParm.Geometry.Period;
NLevel = OptParm.Geometry.Level;
Theta = OptParm.Input.Theta;
Direction = OptParm.Input.Direction;
Thickness = OptParm.Geometry.Thickness;
% 仿真设置
appputvar(h,'Theta',Theta);
appputvar(h,'Direction',Direction);
appputvar(h,'Period',OptParm.Geometry.Period(1)*1e-9);
appputvar(h,'dx',OptParm.Simulation.GridScale*1e-9);
appputvar(h,'dy',OptParm.Simulation.ZGridScale*1e-9);
code=strcat('switchtolayout;',...
    'setnamed("::model","Period",Period);',...
    'setnamed("mesh","dx",dx);',...
    'setnamed("mesh","dy",dy);',...
    'setnamed("source","angle theta",Theta);',...
    'setnamed("::model","Direction",Direction);');
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
Wavelength0 = mean(Wavelengths);
NumWave = length(Wavelengths);
% 细网格初始化
[xGrid,~, xGridScale] = DefineGrid(OptParm.Simulation.Grid, Period, Wavelength0);
Nx = length(xGrid); %Number of x grid points
FinePattern = FineGrid(InPattern,Period,Nx/length(InPattern),0);
% number of layers
Nz = NLevel -1;
% thickness of one layer in nm
zGridScale = OptParm.Geometry.Thickness/Nz;
% 容差
Deviation = OptParm.Optimization.Robustness.StartDeviation;
NRobustness = length(Deviation);
% 容器
AbsEff = zeros(NRobustness,NumPol,NumWave);

% begin loop
for robustIter = 1:NRobustness
    % filter to model physical edge deviations
    BlurGrid = OptParm.Geometry.Pixel/xGridScale/2*Deviation(robustIter);
    if BlurGrid>0
        FinalPattern = GaussFilter2D(FinePattern,BlurGrid);
    else
        FinalPattern = FinePattern;
    end
    % FDTD setting & structure
    Height = FinalPattern/(NLevel -1)*Thickness;
    appputvar(h,'Height',Height*1e-9);
    code=strcat('switchtolayout;',...
        'select("::model");',...
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
        appputvar(h,'F1',Wavelengths(1)*1e-9);
        appputvar(h,'F2',Wavelengths(end)*1e-9);
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
FOM = OptParm.Optimization.Robustness.Weights*mean(mean(AbsEff,3),2);
appclose(h);
end
