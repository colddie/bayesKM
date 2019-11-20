
% get AIF curve
to = 1:2:49;
aifs = [0 0 0 0 25 105 220 350 440 485 430 300 180 110 104 108 115 125 115 108 98 90 98 108 112];
ts = 0:0.1:49;
aif = interp1(to,aifs,ts,'spline');

% read phantom images in one voxel











% compute patlak for each voxel












% simulate projections







% reconstruct frames and save them







% derive AIF from images (optional)






% compute patlak for each voxel