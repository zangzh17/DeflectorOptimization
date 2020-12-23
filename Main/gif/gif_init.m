function gif_init(gif_filename)
gif_filename = [gif_filename,'.gif'];
% Check for an existing .gif file by the same name:
if exist(gif_filename,'file')==2
    % Ask the user if (s)he wants to overwrite the existing file:
    choice = questdlg(['The file  ',gif_filename,' already exists. Overwrite it?'], ...
        'The file already exists.','Overwrite','Cancel','Cancel');
    % Overwriting basically means deleting and starting from scratch:
    if strcmp(choice,'Overwrite')
        delete(gif_filename)
    else
        clear gif_filename firstframe DelayTime DitherOption LoopCount frame
        error('The giffing has been canceled.')
    end
end
end