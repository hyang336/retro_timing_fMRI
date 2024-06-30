# Find mapping between 0.05 Hz V1 BOLD phase with task onset
import nibabel as nib
import nilearn as nl
import nilearn.image
import nilearn.masking
import numpy as np
from scipy.signal import filtfilt
from scipy.signal import butter
import scipy.signal as signal
import matplotlib.pyplot as plt
import os
import pandas as pd
import argparse
import re

#parse arguments
parser = argparse.ArgumentParser(description='Extract narrow-frequncy phase from mean ROI BOLD signal')
parser.add_argument('--targetF', type=float, help='which frequency to extract phase from', default=0.05)
parser.add_argument('--ROI', type=str, help='which region to extract signal from', default='V1')

args = parser.parse_args()
targetF=args.targetF
ROI=args.ROI

# parameters
yaams_dir = '/oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/'
haams_dir = '/oak/stanford/groups/awagner/yaams-haams/fmri/haams/bids_trimmed/derivatives/fmriprep-23.0.1/'
output_dir = '/oak/stanford/groups/awagner/HY/retro_timing_fMRI/'

if ROI=='V1':
    mask = 'V1_exvivo-bilat.nii.gz'
elif ROI=='lFusiform':
    mask = 'fusiform-lh.nii.gz'

fs=0.5 # sampling frequency in Hz
highcut=targetF+0.01
lowcut=targetF-0.01
sub_list=['sub-015','sub-022','sub-024','sub-026','sub-027','sub-028','sub-029',
           'sub-030','sub-031','sub-033','sub-034','sub-035','sub-036','sub-037',
            'sub-040','sub-042','sub-043','sub-044','sub-045','sub-047','sub-048',
            'sub-049','sub-050','sub-051','sub-052','sub-053','sub-055','sub-056',
            'sub-058','sub-059','sub-061','sub-063','sub-064','sub-065','sub-066',
            'sub-067','sub-068','sub-070','sub-071','sub-072','sub-073','sub-1021',
            'sub-1024','sub-1028','sub-1034','sub-1040','sub-1060','sub-1071',
           'sub-1093','sub-1111','sub-1126','sub-1130','sub-1142','sub-1150','sub-1153',
           'sub-1166','sub-1177','sub-1182','sub-1200','sub-1208','sub-1211','sub-1216',
           'sub-1240','sub-1261','sub-1269','sub-1302','sub-1322','sub-1331','sub-1338',
           'sub-1342','sub-1344','sub-1352','sub-1354','sub-1355','sub-1358','sub-1362',
           'sub-1400','sub-1407','sub-1479','sub-1506']

# initialize pd dataframe to store results
results=pd.DataFrame(columns=['subject','run','phase'])

# regex for subject ID
sub_old_re=re.compile('sub-\d{4}')

# loop over yaam subjects
for sub in sub_list:
    if sub_old_re.match(sub):
        subdir=haams_dir
    else:
        subdir=yaams_dir
    #make subject directory if it doesn't exist
    if not os.path.exists(os.path.join(output_dir,sub)):
        os.makedirs(os.path.join(output_dir,sub),exist_ok=True)
    print('Processing subject:',sub)
    sub_v1=nib.load(os.path.join(subdir,sub,'masks',mask))  
    # detect all preprocessed runs in native space in func using regex
    func_dir=os.path.join(subdir,sub,'func')
    func_list=os.listdir(func_dir)
    func_list=[f for f in func_list if 'space-T1w_desc-preproc_bold.nii.gz' in f]
    # resample V1 mask to functional resolution, assuming all runs are corrected registed to native T1w so we can just use an arbitray run
    run_nii=nib.load(os.path.join(func_dir,func_list[0]))
    ROI_resampled=nilearn.image.resample_to_img(sub_v1,run_nii,interpolation='nearest')
    # replace non-zero values with 1
    ROI_resampled_data=ROI_resampled.get_fdata()
    ROI_resampled_data[ROI_resampled_data>0]=1
    ROI_resampled=nib.Nifti1Image(ROI_resampled_data,affine=ROI_resampled.affine,header=ROI_resampled.header)    
    #save the resampled mask for diagnostic purposes
    nib.save(ROI_resampled,os.path.join(output_dir,sub,ROI+'_resampled_mask.nii.gz'))

    #initialize figures
    plt.clf()
    fig1, axs1 = plt.subplots(len(func_list), 1, figsize=(10, 5 * len(func_list)))
    fig2, axs2 = plt.subplots(len(func_list), 1, figsize=(10, 5 * len(func_list)))

    # loop over runs
    for i, run in enumerate(func_list):
        # Extract BOLD time series from V1
        print('Processing run:',run)
        run_nii=nib.load(os.path.join(func_dir,run))
        ROI_data=nilearn.masking.apply_mask(run_nii,ROI_resampled)
        #average over voxels in V1 in time domain before filtereing
        ROI_data_mean=np.mean(ROI_data,axis=1)
        #pad the data with values at beginning and end to avoid edge effects
        ROI_data_mean=np.pad(ROI_data_mean,(20,20),'edge')
        # Extract 0.05 Hz phase from each voxel, plot as a function of time
        # First, bandpass filter the data from 0.04 to 0.06 Hz (target frequency is 0.05 Hz)
        b,a=butter(8,[lowcut,highcut],'bandpass',fs=fs)        
        ROI_data_mean_filt=filtfilt(b,a,ROI_data_mean)        
        
        #extract only the run-xx part of run name
        runXX=run.split('_')[3]               

        # Compute phase of the analytic signal
        analytic_signal = signal.hilbert(ROI_data_mean_filt)
        phase=np.angle(analytic_signal)

        # save the filtered signal as a text file
        np.savetxt(os.path.join(output_dir,sub,runXX+'_'+ROI+'_mean_BOLD_signal_bandpass.txt'),ROI_data_mean_filt)
        #save phase as a text file
        np.savetxt(os.path.join(output_dir,sub,runXX+'_'+ROI+'_phase_bandpass.txt'),phase)
        #store results in dataframe
        results=pd.concat([results,pd.DataFrame.from_dict({'subject':sub,'run':runXX,'phase':[phase[0]]})])

        #plot phase
        axs1[i].plot(phase)
        
        #plot ROI mean data
        axs2[i].plot(ROI_data_mean_filt)

    if i >=0:
        axs1[-1].set_xlabel('TR')
        axs1[-1].set_ylabel('phase')
        axs2[-1].set_xlabel('TR')
        axs2[-1].set_ylabel('Mean BOLD signal')

    fig1.savefig(os.path.join(output_dir,sub,ROI+'_phase_bandpass.png'))
    fig2.savefig(os.path.join(output_dir,sub,ROI+'_mean_BOLD_signal_bandpass.png'))

# convert phase to seconds and account for the fact that phase go from [-pi,pi]
results['phase_in_sec']=(1/targetF)*(results['phase']+np.pi)/(2*np.pi)
#save results
results.to_csv(os.path.join(output_dir,ROI+'_phase_results.csv'))
# Note that this is in trimmed data, so the beginning phase is with the onset of the first ITI, we probably need to look at some later period to find consistent phase since the first ITI is visually indistinguishible from the dummy period





