
% Breaks the pattern down into cells of size Pitch x Pitch
% Each recieves a random refractive index
% Then the entire pattern is blurred by a moving filter
function PatternOut = RandomStart1D(N,Level,RandParm)
    CellPitch = RandParm.Pitch;
    RandAverage = RandParm.Average;
    RandSigma = RandParm.Sigma;
    NCells = ceil(N/CellPitch);

    % Normally distributed index of refractions
    RandomIndices = RandAverage*ones(1,NCells) + RandSigma*randn(1,NCells);
    RandomPattern = zeros(1,N);
    for ii = 1:N
        RandomPattern = RandomIndices(ceil(ii/CellPitch));
    end
    % Smooth and blur the pattern
    
    RandomPattern = smooth(RandomPattern, 1.5*CellPitch);
    
    RandomPattern = round(RandomPattern*(Level-1));
    % Ensure the pattern is proper
    RandomPattern(RandomPattern<0) = 0;
    RandomPattern(RandomPattern>Level-1) = Level-1;

    PatternOut = RandomPattern;
end