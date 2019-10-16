% displays slices of 4d brain phantom
% parameters: 
% @phantompath     output path on file system for series of 3d volumes
% @slicenum        slice number to be displayed   
% author: michael manhart | michael.manhart@cs.fau.de
%         pattern recognition lab, university of erlangen-nuremberg
% last change: 16.01.2012

function show_phantom(phantompath, slicenum)
if nargin < 2
    error('usage: show_phantom(phantompath, slicenumber)');
end
if slicenum < 1 || slicenum > 256
    error('usage: wrong slice number');
end

if ~exist(phantompath,'dir')
    error(strcat('directory not found ', phantompath));
end
maxIdx = 1;
filename = fullfile(phantompath, int2str(maxIdx)); 
while exist(filename,'file')
    maxIdx = maxIdx + 1;
    filename = fullfile(phantompath, int2str(maxIdx)); 
end
maxIdx = maxIdx - 1;

if maxIdx < 1
    error('phantom files not found');
end

for i=1:maxIdx
    i
    filename = fullfile(phantompath, int2str(i)); 
    fid = fopen(filename, 'rb');
    if fid == -1
        error(strcat('could not open ', filename));
    end
    fseek(fid,256*256*4*(slicenum-1),'bof');
    brain = fread(fid,256*256, 'float32');
    brain = reshape(brain,256,256);
    fclose(fid);
    if i==1
        firstimage = brain;
    end

    figure(1);   imshow(brain-firstimage,[0 80],'InitialMagnification',200);
    if i == 20
        keyboard
        M=tissueclasses();
        colormap(M)
    end
    colormap gray
    pause;
end


end