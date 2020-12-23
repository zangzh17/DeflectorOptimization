% Compute the simulation grid for given geometry
 % Grid: Grid points per wavelength 
function [xGrid, yGrid, dr] = DefineGrid(Grid, Period, Wavelength)
    
    %Number of grid points
    NGrid = ceil(Grid*Period/Wavelength);
    if length(Period)==1
        Nx = NGrid;
        Ny = 1;
        %Device period
        Px = Period;
        Py = 0;
    else
        Nx = NGrid(1);
        Ny = NGrid(2);
        %Device period
        Px = Period(1);
        Py = Period(2);
    end
    %Compute external grid coordinates
    xBounds = linspace(0,Px,Nx+1); 
    yBounds = linspace(0,Py,Ny+1);
    
    %Compute size of each grid box
    dx = xBounds(2) - xBounds(1);
    dy = yBounds(2) - yBounds(1);
    
    %Compute coordinates of center of each box
    xGrid = xBounds(2:end)- 0.5*dx;
    yGrid = yBounds(2:end)- 0.5*dy;
    
    %Compute average grid size
    if length(Period)==1
        dr = dx;
    else
        dr = mean([dx dy]);
    end
    
end