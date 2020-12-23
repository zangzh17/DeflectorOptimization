clear
%% gen periodic line patterns
prefix = '0912';
period_list = [37,60,100];
period_num = containers.Map([37,60,100],[6,5,4]); 
offset_list = [0,0.2,0.4];
I = [];
cnt = 0;
file_num = 12;% maximum waveforms per file page
sep_size = 100;
for kk = 1:length(period_list)
    period = period_list(kk);
    num = period_num(period);
    for nn = 1:length(offset_list)
        offset = offset_list(nn);
        period_sampling_num = round(period/0.5);
        tmp = uint8(linspace(0,1-offset,period_sampling_num)*255);
        tmp = 255 - tmp;
        tmp = repmat(tmp,1,num);
        if kk~=length(period_list) || nn~=length(offset_list)
            tmp = [tmp,zeros(1,round(sep_size/0.5))];
        end
        if mod(cnt,file_num)==0
            I = tmp;            
        else
            I = [I,tmp];
        end
        if mod(cnt,file_num)==file_num-1 || (nn==length(offset_list) && kk==length(period_list))
            save(['gt_',prefix,num2str(ceil(cnt/file_num),'%02d'),'.mat'],'I');
        end
        cnt = cnt+1;
    end
end

%% gen line bmp
bmp_path = './';
prefix = '0912';
y_size = 300; % x/y size in um

I_tot = [];
for kk = 1:1
    load(['gt_',prefix,num2str(kk,'%02d'),'.mat']);
    rep_y = round(y_size/0.5);
    % gen 1D map
    I = repmat(I,rep_y,1);
    I = uint8(round(I));
    if size(I,2)>size(I_tot,2)
        I_tot = [I_tot,zeros(size(I_tot,1),size(I,2)-size(I_tot,2));I];
    else
        I_tot = [I_tot;I,zeros(size(I,1),size(I_tot,2)-size(I,2))];
    end

    imwrite(I, [bmp_path,prefix,num2str(kk,'%02d'),'.bmp']);
end
