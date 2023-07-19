%% Dumb for-loop for looping over runs and random versions
load('C:\Users\haozi\OneDrive\Desktop\Postdoc\Wagner\random_offsets_for_validation\test.mat');

sub_list={'sub-029'};
rand_list=[1];
run_list=[1,2,3,4,5,6,7,8];

%save output
vname={'sub','rand_v','run','estimated shift (x in onset-x)'};
vtype={'string','double','double','double'};
minoffset=table('Size',[0,4],'VariableNames',vname,'VariableTypes',vtype);

for i=1:length(sub_list)
    sub=sub_list{i};
    for j=1:length(rand_list)
        rand_v=rand_list(j);
        for k=1:length(run_list)
            run=run_list(k);
            x=zoom_wrapper(fmriprep_dir,derivative_dir,behav_dir,sub,run,rand_v,output,TR,mask,'start_time',-5,'end_time',5);
            minoffset=[minoffset;{sub,rand_v,run,x}];
        end
    end
end