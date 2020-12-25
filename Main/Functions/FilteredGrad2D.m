% Computes the base pattern gradient from the density and level
% filtered gradient
function GradientOut = FilteredGrad2D(GradientIn,PatternIn,Bin,Midpoint,Radius)

% Compute the derivative of the level filter
% level filter is a piecewise function for each step
% centered around the Midpoint of each step
LowLevels = floor(PatternIn);
PatternIn = PatternIn - LowLevels;
if Bin~=0
    PatternLow = Bin*(exp(-Bin*(1-PatternIn/Midpoint)))+exp(-Bin);
    PatternHigh = exp(-Bin) + Bin*exp(-Bin*(PatternIn-Midpoint)/(1-Midpoint));
else
    PatternLow = ones(size(PatternIn))/(2*Midpoint);
    PatternHigh = ones(size(PatternIn))/(2*(1-Midpoint));
end

% Combine the two pieces of the threshold derivative
LowIndex = PatternIn<=Midpoint;
HighIndex = PatternIn>Midpoint;
PatternDeriv = zeros(size(PatternIn));
PatternDeriv(LowIndex) = PatternLow(LowIndex);
PatternDeriv(HighIndex) = PatternHigh(HighIndex);

% Apply chain rule
Gradient = PatternDeriv.*GradientIn;

% Chain rule of the DensityFilter is another DensityFilter
GradientOut = DensityFilter2D(Gradient,Radius);
end