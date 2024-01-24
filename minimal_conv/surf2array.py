#Convert surface data (cifti) based on corresponding geometry file (gifti) to 2D array
import nibabel as nib
import numpy as np

#file with correspondance between surface vertices and 3D spatial coordinates
L_surf=nib.load('/home/users/hyang336/retro_timing_fMRI/ROIs/MMP_left/Q1-Q6_RelatedParcellation210.L.flat.32k_fs_LR.surf.gii')
L_surf_array=[x.data for x in L_surf.darrays]
L_surf_3Dcoord=L_surf_array[0]
L_surf_triangles=L_surf_array[1]

#file with parcellation labels and vertices
MMP_label=nib.load('/home/users/hyang336/retro_timing_fMRI/ROIs/MMP_left/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Areas_Group_Colors_Atlas_ROIs2.32k_fs_LR.dlabel.nii')
#label for each voxel/vertex
label_codes=MMP_label.dataobj[0]

#MMP_data=MMP_label.get_fdata(dtype=np.float32)#this thing is 1 by 91282 (grayordinates-surface and subcortical)
MMP_header_axes=[MMP_label.header.get_axis(i) for i in range(MMP_label.ndim)] #using cifti2_axes API from nibabel to extract data
#need to figure out how to extract ROI from the two axes, and the correspondence between the BrainModel axis and the surface GIFTI (maybe it is redundant)
MMP_label_axis=MMP_header_axes[0]
#label code to name dictionary
label_dict=MMP_label_axis.label[0]
#Left motor key
L_M1_key=[k for k,v in label_dict.items() if v[0]=='L_4_ROI'][0]
#left and right V1


#find the index (based on the 91k grayordinate) on the left cortical surface that has the label L_4_ROI (i.e. primary motor)
# MMP_brainmodel_axis=MMP_header_axes[1]
L_M1_grayord=np.where(label_codes.astype(int)==L_M1_key)[0]#according to Glasser et al. (2016) fig.1, these indices on the L cortex do not need to be modified.

#


##!! would be too memory-intensive to check which triangle is detached, need to come up with a general solution that can handle any triangles
# #detect triangles with only 1 (:.:) or 0 (:. :.) shared vertices
# triangles, _ = L_surf_triangles.shape
# results_rows_1shared = []
# results_rows_0shared = []

# #use numpy array broadcasting to do pairwise comparisons to figure out which elements are the same
# #if you don't understand this, run:
# #arr=np.array([[1,2,3],[1,4,5],[1,2,7],[10,9,6]])
# #arr.shape
# #arr[None,:,:]
# #arr[:,None,:]
# #arr[:,:,None]
# #arr[:,None,:]==arr[None,:,:]

# #symmetric matrix counting the number of overlapping elements across row pairs, this matrix is huge and you probably need tenth of GB memory to fit it in.
# vert_1shared=np.sum(L_surf_triangles[:,None,:]==L_surf_triangles[None,:,:],axis=-1)#each element (i,j) is the number of overlapping element between row i and row j in the original array


# for i in range(triangles):
#     share_vert_count=0
#     for j in range(triangles):
#         if i!=j: #don't compare a triangle with itself
#             share_vert_count += len(np.intersect1d(L_surf_triangles[i],L_surf_triangles[j]))
#     if share_vert_count == 0:
#         results_rows_0shared.append(i)
#     elif share_vert_count ==1:
#         results_rows_1shared.append(i)

#