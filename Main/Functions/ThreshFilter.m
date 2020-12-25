% Linear -> step function as Bin increases
% Midpoint within [0,1] is the midpoint of the function
% PatternIn is within [0,1], PatternOut is within [0,1]
function PatternOut = ThreshFilter(PatternIn,Bin,Midpoint)
    % The threshold filter is a piecewise function centered around Midpoint
    if Bin~=0
        % [0,Midpoint] to [0, Midpoint]
        PattNormLow = 1-PatternIn/Midpoint;
        PatternLow = Midpoint*(exp(-Bin*PattNormLow)-PattNormLow*exp(-Bin));
        % [Midpoint,1] to [Midpoint,1]
        PattNormHigh = (PatternIn-Midpoint)/(1-Midpoint);
        PatternHigh = Midpoint + (1-Midpoint)*(1-exp(-Bin*PattNormHigh)+PattNormHigh*exp(-Bin));
        
    % For Bin=0, maps Midpoint to .5; and linearly on either side
    elseif Bin==0
        % [0,Midpoint] to [0, Midpoint]
        PatternLow = PatternIn/(2*Midpoint);
        % [Midpoint,1] to [Midpoint,1]
        PatternHigh = (PatternIn-1)/(2-2*Midpoint)+1;   
    end
    LowIndex = PatternIn<=Midpoint;
    HighIndex = PatternIn>Midpoint;
    PatternOut = zeros(size(PatternIn));
    PatternOut(LowIndex) = PatternLow(LowIndex);
    PatternOut(HighIndex) = PatternHigh(HighIndex);
end