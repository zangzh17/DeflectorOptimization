%Gaussian weighted averaging within radius Sigma
function PatternOut = GaussFilter1D(PatternIn,Sigma)
    PatternOut = imgaussfilt(PatternIn, Sigma, 'padding', 'circular');
end