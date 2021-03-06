% Transfers the pattern onto a finer/coarser grid specifid by the scaling
% factor Scale
function [PatternOut,XOut,YOut] = FineGrid(PatternIn,Period,Scale,SmoothGeom)

Pattern = PatternIn;
[Nx, Ny] = size(PatternIn); % Size of input pattern
% Compute size of output pattern
if length(Scale)==1
    if Nx==1
        NxScaled = 1;
        NyScaled = round(Scale.*Ny);
        if length(Period)==1, Period = [0,Period]; end
    elseif Ny==1
        NxScaled = round(Scale.*Nx);
        NyScaled = 1;
        if length(Period)==1, Period = [Period,0]; end
    end
else
    NxScaled = round(Scale(1).*Nx);
    NyScaled = round(Scale(2).*Ny);
end


% Blur pattern if option is specfied
if SmoothGeom == 1
    Pattern = imgaussfilt(Pattern,max(Nx,Ny)/50,'Padding','circular');
end

% Add extra pixels in order to allow for periodic interpolation
if length(Scale)==1
    if Nx==1
        Pattern = [Pattern(end),Pattern];
    elseif Ny==1
        Pattern = [Pattern(end);Pattern];
    end
else
    Pattern = [Pattern(:,end),Pattern];
    Pattern = [Pattern(end,:);Pattern];
end

% Input grid
if length(Scale)==1
    if Nx==1
        [X,Y] = meshgrid(0,0:Ny);
    elseif Ny==1
        [X,Y] = meshgrid(0:Nx,0);
    end
else
    [X,Y] = meshgrid(0:Nx,0:Ny);
end

[XScaled,YScaled] = meshgrid(linspace(Nx/NxScaled,Nx,NxScaled),...
        linspace(Ny/NyScaled,Ny,NyScaled));
% Output coordinates
% Scaled grid
[XOut,YOut] = meshgrid(linspace(0,Period(1),NxScaled),...
        linspace(0,Period(2),NyScaled));

if SmoothGeom == 1
    % Use interpolation if smooth geometry is desired
    if length(Scale)==1
        if Nx==1
            PatternOut = interp1(Y,Pattern',YScaled,'PCHIP')';
        elseif Ny==1
            PatternOut = interp1(X,Pattern',XScaled,'PCHIP')';
        end
    else
        PatternOut = interp2(X,Y,Pattern',XScaled,YScaled,'cubic')';
    end
else
    % Place pixels on a square grid otherwise
    PatternOut = zeros([NxScaled NyScaled]);
    for ii = 1:NxScaled
        for jj = 1:NyScaled
            PatternOut(ii,jj) = Pattern(ceil(ii*Nx/NxScaled),...
                ceil(jj*Ny/NyScaled));
        end
    end
end

end