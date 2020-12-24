% continuous -> discrete as B increases
% Midpoint within (0,1) sets the threshold between levels
function PatternOut = LevelFilter(PatternIn,Bin,Midpoint)
if nargin<3
    Midpoint = 0.5;
end
LowLevels = floor(PatternIn);
PatternOut = LowLevels + ThreshFilter(PatternIn - LowLevels,Bin,Midpoint);
end