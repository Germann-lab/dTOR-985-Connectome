% Script to seed the structural connectome with an ROI input, calculating a weight for each fiber streamline based on intersection with the ROI
% Suresh Joel, General Electric Global Research, April 2018

% The following lines select the script input and require the user to specify this input as instructed
load('/PATH_TO_FIBER_VOX_FILE_OF_CHOICE/FIBER_VOX_FILE_OF_CHOICE.mat'); % Provide path and filename for voxel-resampled connectome file (this determines the resolution in which the streamline sampling will be conducted)
vta_path='/PATH_TO_INPUT_ROI_FILE'; % Insert path to ROI file here
vta_file = fullfile(vta_path,'INPUT_ROI_FILE.nii'); % Insert ROI file here (must be in unzipped nifti format and matching the resolution of the voxel-resampled file being used)
v=spm_vol(vta_file);
im=spm_read_vols(v);
im=im./max(im(:));
im(im<0)=0;

% The following lines initialize fiber weights to zero and then assign weights based on fiber-ROI intersection
fibers_wt=zeros(length(fibers_vox),1);

for i=1:length(fibers_vox)
    for j=1:size(fibers_vox{i},1)
        if (any(fibers_vox{i}(j,:)<=0) || fibers_vox{i}(j,1)>size(im,1) || fibers_vox{i}(j,2)>size(im,2) || fibers_vox{i}(j,3)>size(im,3))
            continue
        else
            fibers_wt(i)=fibers_wt(i) + im(fibers_vox{i}(j,1),fibers_vox{i}(j,2),fibers_vox{i}(j,3));
        end
    end
    if rem(i,10000)==0
        disp(['Fiber number:',num2str(i)])
    end
end

save([vta_file(1:end-4),'_fiber_wts.mat'],'fibers_wt');