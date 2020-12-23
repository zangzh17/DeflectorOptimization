
% Breaks the pattern down into cells of size Pitch x Pitch
% Each recieves a random array within [0,1]
% Then the entire pattern is blurred by a gaussian filter
function PatternOut = RandomStart(Nx,Ny,Period,RandParm,SymX,SymY)
    Pitch = RandParm.Pitch;
    RandAverage = RandParm.Average;
    RandSigma = RandParm.Sigma;
    
    if length(Period)>1
        GridSize = Period(1)/Nx;
        NCellsX = ceil(2*Period(1)/Pitch);
        NCellsY = ceil(2*Period(2)/Pitch);
    else
        if Nx==1
            GridSize = Period/Ny;
            NCellsX = 1;
            NCellsY = ceil(2*Period/Pitch);
            Period = [Pitch,Period];
        elseif Ny==1
            GridSize = Period/Nx;
            NCellsY = 1;
            NCellsX = ceil(2*Period/Pitch);
            Period = [Period,Pitch];
        end
    end
    
    % Normally distributed index of refractions
    RandomIndices = RandAverage*ones(NCellsX,NCellsY) + RandSigma*randn(NCellsX,NCellsY);
    RandomIndices = EnforceSymmetry(RandomIndices, SymX, SymY);
    
    % Smooth and blur the pattern
    [RandomPattern,~,~] = FineGrid(RandomIndices,Period,[Nx/NCellsX, Ny/NCellsY],0);
    DiskFilter = fspecial('disk',1.1*Pitch/GridSize);
    RandomPattern = imfilter(RandomPattern,DiskFilter,'circular');

    % Ensure the pattern is proper
    RandomPattern = EnforceSymmetry(RandomPattern, SymX, SymY);
    RandomPattern(RandomPattern<0) = 0;
    RandomPattern(RandomPattern>1) = 1;

    PatternOut = RandomPattern;
end