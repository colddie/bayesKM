
%info=dicom_folder_info('X:\storage0\home\sedata\PIG_TPP_BAM15_20180614\PET\TPP DYN 2.5HRS VPFX\A')
%info=dicom_folder_info('X:\storage1\tsun\20180702_CardiacPhantom\ref')
%info=dicom_folder_info('X:\storage1\tsun\ecv\pig_e\pet')
% info=dicom_folder_info('X:\GCMI_DATA\TPP\TPP_PIG_20181219\CTs\Precontrast')
% info=dicom_folder_info('X:\storage1\tsun\cbct\phantomplanCT')
% info=dicom_folder_info('/home/tsun/work/MK6240_datasets_for_Tao/SMK005/Dicom_MPRAGE')
info=dicom_folder_info('/home/tsun/work/MK6240_datasets_for_Tao/SMK005/Dicom_PET_dyn_1')

for nseries = 1: numel(info)
    % mkdir(num2str(nseries))  
    dir =  [pwd,'/',info(nseries).SeriesTime,info(nseries).SeriesDescription];
    mkdir(dir);
  for i = 1: numel(info(1,nseries).Filenames) 
    copyfile(num2str(info(1,nseries).Filenames{i,1}), [dir,'/',num2str(i)]);     
  end
end


% then read dicoms into matlab, as niread_dicom2 fails to read thw whole dataset
% hdr = dicom_read_header();
% img = dicom_read_volume(hdr);
% writeraw('post.raw',img,'int16')