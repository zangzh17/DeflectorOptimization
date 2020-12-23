function gif(varargin)

%% Using Information 
% creat first frame:
%
% gif('filename.gif')
%
%
% record frame:
%
% gif('filename') (use gca)
% or
% gif('filename','frame',gcf) (use full current figure rather than current axis area)
% (frame must be a figure handle or axis handle)


% % Define persistent variables: 
% persistent frame
%% demo
if nargin==0
    % Some sample data:
    t = sin(linspace(0,2*pi,30));
    [X,Y,Z] = peaks(500);
    
    % Plot the first frame:
    h = surf(X,Y,Z*t(1));
    shading interp
    axis([-3 3 -3 3 -9 9])
    
    % Make it fancy:-
    camlight
    set(gca,'color','k')
    set(gcf,'color','k')
    caxis([min(Z(:)) max(Z(:))])
    % first frame
    gif_init('myfile')
    gif('myfile')
    % rest frame
    for k = 2:30
        set(h,'Zdata',Z*t(k))
        disp(k)
        gif('myfile')
    end
    close
    web('myfile.gif')
    return
end
%% Parse Inputs

gif_filename = [varargin{1},'.gif'];


%% Perform work:

% Set defaults frame
frame = gcf;
% set user defined frame
tmp = strcmpi(varargin,'frame');
if any(tmp)
    frame = varargin{find(tmp)+1};
    assert(ishandle(frame)==1,'Error: frame must be a figure handle or axis handle.')
end

% gif setting
DelayTime = 1/30;
DitherOption = 'dither';
LoopCount = Inf;

% Get frame: 
f = getframe(frame); 

% Convert the frame to a colormap and corresponding indices: 
[imind,cmap] = rgb2ind(f.cdata,256,DitherOption);    

% Write the file:     
if exist(gif_filename,'file')==0
    % first frame
   imwrite(imind,cmap,gif_filename,'gif','LoopCount',LoopCount,'DelayTime',DelayTime)
else
    % rest frame
   imwrite(imind,cmap,gif_filename,'gif','WriteMode','append','DelayTime',DelayTime)
end

end