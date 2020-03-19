% creates colormap for displaying segmented brain
% author: michael manhart | michael.manhart@informatik.uni-erlangen.de
%         pattern recognition lab, university of erlangen-nuremberg
% last change: 16.01.2012

function M = tissueclasses()
% classes of tissue
M = [0      0       0;
     127    127     127; % GM
     255    255     255; % WM
     220    220     0  ; % GMR
     220    0       0  ; % GMSR
     255    255     0  ; % WMR
     255    0       0  ; % WMSR
     255    1       255; % AI
     0      0       125; % VO
     0      255      0; % CSF
     0      0      220; % HEWM
     0      0      255]/255; % HEGM