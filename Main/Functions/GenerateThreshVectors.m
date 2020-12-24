% Interpolate between the start and end edge deviations to create
% and edge deviation value for each iteration
function [ThresholdVectors, NRobustness] = GenerateThreshVectors(OptParm)
    % Extract necessary parameters
    MaxIterations = OptParm.Optimization.Iterations;
    Start = OptParm.Optimization.Robustness.StartDeviation;
    End = OptParm.Optimization.Robustness.EndDeviation;
    Ramp = OptParm.Optimization.Robustness.Ramp;
    NRobustness = length(Start);
    
    % Error if the dimensions don't match
    if NRobustness ~= length(End)
       error('Robustness vectors are not the same length!');
    end
    
    % Generate edge deviation vector
    DeviationVectors = zeros(NRobustness, MaxIterations);
    for ii = 1:NRobustness
       DeviationVectors(ii, 1:Ramp) = linspace(Start(ii), End(ii), Ramp);
       DeviationVectors(ii, (Ramp+1):end) = End(ii);
    end
    
    % Generate corresponding binarize values
	% 0.5*(1-erf(x))
    % For very small deviations, Threshold = 0.5
    % For deviations = 0.2*BlurRadius, Threshold is about 0.4
    % For deviations = 0.5*BlurRadius, Threshold is about 0.24
    % For deviations = BlurRadius, Threshold is about 0.1
    % For deviations = 2*BlurRadius, Threshold is about 0.002
    ThresholdVectors = 0.5*(1-erf(DeviationVectors/OptParm.Optimization.Filter.BlurRadius));
end