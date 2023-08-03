%% compile results and generate figures across subjects
function results(ss_list,data_dir,run_list,output_dir)

%read in ss list and run list
ss_open=fopen(ss_list,'r');
SSID=textscan(ss_open,'%s', 'Delimiter', '\n');
SSID=SSID{1};%this is where all lines are stored
SSID(cellfun('isempty',SSID))=[];%get rid of empty cells/lines

run_open=fopen(run_list,'r');
runs=textscan(run_open,'%s','Delimiter','\n');
runs=runs{1};
runs(cellfun('isempty',runs))=[];

% Loop over subject, run and fold to generate figures and accuracy measure

    %the plots are always based on fold-1 to give a general idea on how the
    %GLM did over the entire search range

    %the accuracy is pulled out from the highest fold to have the highest
    %potential temporal resolution

% generate some plots

    %for each subject

    %across subjects

end