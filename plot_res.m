fold=1;
rand_v=1;
run=1;
res_dir=strcat(output,'/',sub,'_Rand_',num2str(rand_v),'_Run_',num2str(run),'_ResMS_fold-',num2str(fold));
s_time=-5;
e_time=5;
bin_num=101;
tile=linspace(s_time,e_time,bin_num);
volume_avg_res=[];
    for i=1:length(tile)
        tile_str{i}=sprintf('%g', tile(i));
        volume_avg_res(i)=mean(niftiread([res_dir,'/ResMS',tile_str{i},'.nii']),"all","omitnan");
    end
plot(tile,volume_avg_res);

minind=find(volume_avg_res==min(volume_avg_res));
filename=tile_str{minind}