function create_atlasnifti (sch_vector,name,symm)
% script to create Schaefer nifti file for rendering
%  sch_vector is a 1x1000 vector with a value for each node in the Schaefer
%  1000 parcellation
%  name is the name of the output file
%  symm is a flag that if >0 means that ordering is symmetrical
addpath('/path/to/nifti-tools'); %path for tools for nifti and analyze image
atlasnii=load_nii('/path/to/parcellations/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii');

% create new image
newnii=atlasnii;
% set datatype to correct single type (rather than uint8)
newnii.img=zeros(size(atlasnii.img),'single');
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
newnii.hdr.dime.datatype=16;
newnii.hdr.dime.bitpix=32;

if (symm>0) 
    % convert symmetrical to standard
    schvals(1:2:1000)=sch_vector(1:500);
    schvals(2:2:1000)=sch_vector(1000:-1:501);
else
    sch_vals=sch_vector;
end

for i=1:1000
    % find the indices of voxels for a given Schaefer region
    vox=find(atlasnii.img==i);
    % set value for each voxel in Schaefer region 
    newnii.img(vox)=sch_vals(i);
end


%save new nifti file
save_nii (newnii, name);

