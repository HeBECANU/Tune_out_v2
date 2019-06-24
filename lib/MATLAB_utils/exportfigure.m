function exportfigure(figfile,savedir,caption)

%% Slick save function for figure exporting

%%

% timenow = posixtime(datetime);
% timestamp = sprintf( '%.0f',timenow);
timestamp = [];
png_title = fullfile(savedir,strcat(caption,'_',timestamp,'.png'));
fig_title = fullfile(savedir,strcat(caption,'_',timestamp,'.fig'));

if exist(fig_title,'file') || exist(png_title,'file')
    if exist(png_title,'file')
      warning('Filename (PNG) already exists!')
    end
    if exist(fig_title,'file')
      warning('Filename (FIG) already exists!')  
    end
    
    % 28 leftarrow
    % 29 rightarrow
    % 30 uparrow
    % 31 downarrow
    warning('Press BSPACE to overwrite or ENTER to skip');
    waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    if value == 8
        savefig(figfile,fig_title);
        saveas(figfile,png_title);
    elseif value == 13
        return
    end
end
    
end