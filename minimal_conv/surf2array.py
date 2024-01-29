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

#find all triangles on surface mesh that contain any one of the ROI vertices
triangle_mask_vert_logic=np.isin(L_surf_triangles,L_M1_grayord)
triangle_mask_3vert=np.all(triangle_mask_vert_logic,axis=1)
triangle_L_M1=L_surf_triangles[triangle_mask_3vert]

#check if all triangles have at least one other triangle that shares an edge with it (i.e. the ROI mesh is continuous)
share_edge=np.full(triangle_L_M1.shape[0],False)
triangles, _ = triangle_L_M1.shape
for i in range(triangles):
    for j in range(triangles):
        if i!=j:
            if len(np.intersect1d(triangle_L_M1[i],triangle_L_M1[j]))==2:
                share_edge[i]=True

#if there are triangles that doesn't share edge with other triangles, check if they at least share vertex
dangle_triangles=np.where(share_edge==False)
share_vertex=np.full(triangle_L_M1.shape[0],True)
if len(dangle_triangles[0])!=0:
    for i in range(len(dangle_triangles[0])):
        dangle_triangle_idx=dangle_triangles[0][i]
        for j in range(triangles):
            if dangle_triangle_idx!=j:#dangle_triangle is a tuple with its 1st element being a numpy array of indices of triangles
                if len(np.intersect1d(triangle_L_M1[dangle_triangle_idx],triangle_L_M1[j]))==1:
                    share_vertex[dangle_triangle_idx]=False

#if the list of triangles that don't share edges is not equal to the list of triangles sharing vertex, then these triangles are orphans, remove these triangles from the ROI
orphan_triangles=np.where(share_edge!=share_vertex)
triangle_L_M1_cont=np.delete(triangle_L_M1,orphan_triangles[0],axis=0)

#since the ROI is defined on the group (template) surface space, we do not need to account for different real-world position of the vertices

### Converting triangular mesh to matrix seems to necessitate interpolation or rasterization, both will cause info loss. So the alternative strategy is to apply conv directly on the mesh, via something like graph neuralnet