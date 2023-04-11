fold=2;
res_dir=strcat(output,'/',sub,'_ResMS_fold-',num2str(fold));
s_time=2.7;
e_time=3.7;
bin_num=101;
tile=linspace(s_time,e_time,bin_num);
volume_avg_res=[];
    for i=1:length(tile)
        tile_str{i}=sprintf('%g', tile(i));
        volume_avg_res(i)=mean(niftiread([res_dir,'/ResMS',tile_str{i},'.nii']),"all","omitnan");
    end
plot(tile,volume_avg_res);