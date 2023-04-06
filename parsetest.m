function par= parsetest(fmriprep_dir,behav_dir,sub,output,TR,age,varargin)

%% input parser
p=inputParser;

default.bin_num=101; %number of endpoints between the starting and the end time points (bin number + 1)

default.start_time=0;
default.end_time=10;

addParameter(p,'bin_num',default.bin_num,@isint);
addParameter(p,'start_time',default.start_time,@isnumeric);
addParameter(p,'end_time',default.end_time,@isnumeric);

parse(p,varargin{:});

disp(p.Results.start_time);
disp(p.Results.end_time);
disp(p.Results.bin_num);

par= p;
end