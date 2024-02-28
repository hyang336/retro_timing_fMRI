#script to apply mask on nii, and convert nii to h5, for residual images from GLM in MNI space
import h5py
import numpy as np
import nibabel as nib
import os
import sys
import glob
import re
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert nii to h5')
    parser.add_argument('input', help='input directory')
    parser.add_argument('output', help='output directory')
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    files = glob.glob(os.path.join(input_dir, 'ResMS*.nii'))

    #load in mask, which is the same for all images
    mask = nib.load('mask.nii')
    mask_data = mask.get_data()

    #load in files
    for f in files:
        img = nib.load(f)
        data = img.get_data()
        data = data * mask_data
        data = data.astype(np.float32)
        data = np.expand_dims(data, axis=0)
        #save as h5
        filename = os.path.basename(f)
        filename = os.path.splitext(filename)[0]
        filename = filename + '.h5'
        with h5py.File(os.path.join(output_dir, filename), 'w') as hf:
            hf.create_dataset('data', data=data)
            
if __name__ == '__main__':
    main()







