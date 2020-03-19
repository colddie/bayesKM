% creates 4D perfusion phantom brain volume
% parameters: 
% @phantompath     output path on file system for series of 3d volumes
% @dt              time step between 3d volumes in s 
% @addskull        indicates if skull is included in 4D volume
% author: michael manhart | michael.manhart@cs.fau.de
%         pattern recognition lab, university of erlangen-nuremberg
% last change: 17.04.2013

function create_phantom(phantompath, dt, addskull)
% define perfusion parameters for different matters
% gm -> gray matter wm -> white matter
% gm_cbv -> healthy gray matter blood flow
% gmr_cbv -> reduced gray matter blood flow
% gmsr_cbv -> severely reduced gray matter blood flow etc...
% cbv in ml/100ml  cbf in ml/100 ml/min  mtt in seconds
gm_cbv = 3.3;     gm_cbf = 53;    gm_mtt = (gm_cbv/gm_cbf)*60;   
gmr_cbv = 3;      gmr_cbf = 16;   gmr_mtt = (gmr_cbv/gmr_cbf)*60;  
gmsr_cbv = 0.71;  gmsr_cbf = 5.3; gmsr_mtt = (gmsr_cbv/gmsr_cbf)*60;   

wm_cbv = 1.9;     wm_cbf = 25;    wm_mtt = (wm_cbv/wm_cbf)*60;   
wmr_cbv = 1.7;    wmr_cbf = 7.5;  wmr_mtt = (wmr_cbv/wmr_cbf)*60;    
wmsr_cbv = 0.42;  wmsr_cbf = 2.5; wmsr_mtt = (wmsr_cbv/wmsr_cbf)*60;  
% define perfusion parameter variation w.r.t. mr data
                  cbf_var = 14;   mtt_var = 0.7;  
                  cbfr_var = 4.25;mttr_var = 0.75;
                  cbfsr_var = 1.4;mttsr_var = 1;                    
% tissue attenuation offsets
HU_GM = 35; HU_WM = 29; HU_V = 40; HU_CSF = 12;
HU_GM_VAR = 3; HU_WM_VAR = 3; HU_CSF_VAR = 3;




% classes of tissue
GM = 1; WM = 2; GMR = 3; GMSR = 4; WMR = 5; WMSR = 6; AI = 7; VO = 8; CSF = 9; HEGM = 10; HEWM = 11;

if nargin < 1
    error('usage: create_phantom(phantompath, dt = 1, addskull = 1)');
end
if nargin < 2
    dt = 1;
end
if nargin < 3
    addskull = 1;
end

if ~exist(phantompath,'dir')
    error(['directory not found ' phantompath]);
end
% load brain maps
if ~exist('brain.mat','file')
    error('brain.mat not found! create stroke annotation with strokecreator before running create_phantom!');
end
brainstruct = load('brain.mat','brain');
brain = uint8(brainstruct.brain);
brain(brain == HEWM) = WM;
brain(brain == HEGM) = GM;
clear brainstruct;
if ~exist('mrbrain.raw','file')
    error('mrbrain.raw not found!');
end
% load brain mr data for perfusion parameter variation
fid = fopen('mrbrain.raw','rb');
mrbrain = zeros(256,256,256);
for i=1:256
    mrbrain(:,:,i) = fread(fid,[256 256],'uint16');
end
fclose(fid);
% load skull data
skull = zeros(256,256,256);
if addskull
    fid = fopen('skull.raw','rb');
    for i=1:256
       skull(:,:,i) = fread(fid,[256 256],'int16');
    end
    fclose(fid);
end

disp('welcome to create_phantom');
disp(['output path is: ' phantompath ' , time step is: ' num2str(dt) ' , addskull is: ' num2str(addskull)]);
% check if phantom files to create already exist
t = 0:dt:49;
fileexists = 0;

if exist(fullfile(phantompath, 'baseline'),'file') || exist(fullfile(phantompath, 'cbf'),'file') ...
    || exist(fullfile(phantompath, 'cbv'),'file') || exist(fullfile(phantompath, 'mtt'),'file') ...
    || exist(fullfile(phantompath, 'ttp'),'file')
    fileexists = 1;
else
    for i=1:size(t,2)
      filename = fullfile(phantompath, int2str(i));  
      if exist(filename,'file')
         fileexists = 1;
         break;
      end
    end
end
if fileexists
    answer = input('baseline file exists. clean up directory?','s');
    if strcmpi(answer,'y') || strcmpi(answer,'yes')
        filename = fullfile(phantompath, 'baseline');
        if exist(filename,'file') delete(filename); end
        filename = fullfile(phantompath, 'cbf');
        if exist(filename,'file') delete(filename); end
        filename = fullfile(phantompath, 'cbv');
        if exist(filename,'file') delete(filename); end
        filename = fullfile(phantompath, 'mtt');
        if exist(filename,'file') delete(filename); end
        filename = fullfile(phantompath, 'ttp');
        if exist(filename,'file') delete(filename); end
        for i=1:size(t,2)
            filename = fullfile(phantompath, int2str(i));  
            if exist(filename,'file') delete(filename); end
        end
    else
        error('ouput path already contains phantom files!');
    end
end
% create aif 0..49 s in 0.1 sec sampling interval
to = 1:2:49;
aifs = [0 0 0 0 25 105 220 350 440 485 430 300 180 110 104 108 115 125 115 108 98 90 98 108 112];
ts = 0:0.1:49;
aif = interp1(to,aifs,ts,'spline');
 
for z=1:size(brain,3)
    disp(['processing slice ' num2str(z)]);
    baselineslice = single(zeros(size(brain,1),size(brain,2)));
    tacslice = single(zeros(size(brain,1),size(brain,2),size(t,2)));
    cbfslice = single(zeros(size(brain,1),size(brain,2)));
    cbvslice = single(zeros(size(brain,1),size(brain,2))); 
    mttslice = single(zeros(size(brain,1),size(brain,2)));
    ttpslice = single(zeros(size(brain,1),size(brain,2))); 
    mrslice = mrbrain(:,:,z);
    brainslice = brain(:,:,z);
    % process white matter  
    wm_idx = find(brainslice==WM | brainslice==WMR | brainslice==WMSR);
    if ~isempty(wm_idx)
        mrvalues = mrslice(wm_idx);
        mrvalues = mrvalues - mean(mrvalues);
        mrvalues_bound = 2*sqrt(var(mrvalues));
        mrvalues(mrvalues < -mrvalues_bound) = -mrvalues_bound;
        mrvalues(mrvalues >  mrvalues_bound) =  mrvalues_bound;
        mrvalues = (mrvalues/mrvalues_bound);
        baselineslice(wm_idx) = HU_WM - (mrvalues*HU_WM_VAR);
          
        for i=1:size(wm_idx,1)
            idx = wm_idx(i);
            switch(brainslice(idx))
                case WM
                    cbf = wm_cbf - mrvalues(i)*cbf_var;
                    mtt = wm_mtt + mrvalues(i)*mtt_var;
                case WMR
                    cbf = wmr_cbf - mrvalues(i)*cbfr_var;
                    mtt = wmr_mtt + mrvalues(i)*mttr_var;
                case WMSR
                    cbf = wmsr_cbf - mrvalues(i)*cbfsr_var;
                    mtt = wmsr_mtt + mrvalues(i)*mttsr_var; 
            end
            tac = GetTissueTac(aif,ts,0.1,mtt,cbf);
            for ti=1:size(t,2)
                tac_index = round(t(ti)/0.1) + 1;
                [slice_x slice_y] = ind2sub(size(baselineslice),idx);
                tacslice(slice_x,slice_y,ti) = tac(tac_index)+baselineslice(idx);
            end
            cbfslice(idx) = cbf;
            cbvslice(idx) = cbf*mtt/60;
            mttslice(idx) = mtt;
            [~, ttp] = max(tac); 
            ttpslice(idx) = ttp*0.1;
        end
    end
    % process gray matter   
    gm_idx = find(brainslice==GM | brainslice==GMR | brainslice==GMSR);
    if ~isempty(gm_idx)
        mrvalues = mrslice(gm_idx);
        mrvalues = mrvalues - mean(mrvalues);
        mrvalues_bound = 2*sqrt(var(mrvalues));
        mrvalues(mrvalues < -mrvalues_bound) = -mrvalues_bound;
        mrvalues(mrvalues >  mrvalues_bound) =  mrvalues_bound;
        mrvalues = (mrvalues/mrvalues_bound);
        baselineslice(gm_idx) = HU_GM - (mrvalues*HU_GM_VAR);
        for i=1:size(gm_idx,1)
            idx = gm_idx(i);
            switch(brainslice(idx))
                case GM
                    cbf = gm_cbf - mrvalues(i)*cbf_var;
                    mtt = gm_mtt + mrvalues(i)*mtt_var;
                case GMR
                    cbf = gmr_cbf - mrvalues(i)*cbfr_var;
                    mtt = gmr_mtt + mrvalues(i)*mttr_var;                 
                case GMSR
                    cbf = gmsr_cbf - mrvalues(i)*cbfsr_var;
                    mtt = gmsr_mtt + mrvalues(i)*mttsr_var; 
            end
            tac = GetTissueTac(aif,ts,0.1,mtt,cbf);
            for ti=1:size(t,2)
                tac_index = round(t(ti)/0.1) + 1;
                [slice_x slice_y] = ind2sub(size(baselineslice),idx);
                tacslice(slice_x,slice_y,ti) = tac(tac_index)+baselineslice(idx);
            end
            cbfslice(idx) = cbf;
            cbvslice(idx) = cbf*mtt/60;
            mttslice(idx) = mtt;
            [~, ttp] = max(tac); 
            ttpslice(idx) = ttp*0.1;
            
            if cbf > 0
                keyboard; end
        end
    end

   % process arteries
    ai_idx = find(brainslice==AI);
    if ~isempty(ai_idx)
        baselineslice(ai_idx) = HU_V;
        for i=1:size(ai_idx,1)
            idx = ai_idx(i);
            for ti=1:size(t,2)
                tac_index = round(t(ti)/0.1) + 1;
                [slice_x slice_y] = ind2sub(size(baselineslice),idx);
                tacslice(slice_x,slice_y,ti) = aif(tac_index)+HU_V;
            end
        end
    end
   % process csf
    csf_idx = find(brainslice==CSF);
    if ~isempty(csf_idx)
        mrvalues = mrslice(csf_idx);
        mrvalues = mrvalues - mean(mrvalues);
        mrvalues_bound = 2*sqrt(var(mrvalues));
        mrvalues(mrvalues < -mrvalues_bound) = -mrvalues_bound;
        mrvalues(mrvalues >  mrvalues_bound) =  mrvalues_bound;
        mrvalues = (mrvalues/mrvalues_bound);
        baselineslice(csf_idx) = HU_CSF - (mrvalues*HU_CSF_VAR);
        for i=1:size(csf_idx,1)
            idx = csf_idx(i);
            for ti=1:size(t,2)
                [slice_x slice_y] = ind2sub(size(baselineslice),idx);
                tacslice(slice_x,slice_y,ti) = baselineslice(idx);
            end
        end
    end    
    % save baseline slice
    skullslice = skull(:,:,z);
    count = size(brain,1)*size(brain,2);
    filename = fullfile(phantompath, 'baseline');
    fid = fopen(filename,'ab');
    if fid == -1
        error(['error opening ' filename]);
    end
    skullslice(baselineslice ~= 0) = 0;
    baselineslice = baselineslice + skullslice;
    if fwrite(fid,baselineslice, 'float32') ~= count
      fclose(fid);
      error(['coud not write file ' filename]);
    end
    fclose(fid);
    % save temporal slices
    for i=1:size(t,2)
      filename = fullfile(phantompath, int2str(i));
      fid = fopen(filename,'ab');
      if fid == -1
          error(['error opening ' filename]);
      end
      tacslice(:,:,i) = tacslice(:,:,i) + skullslice;
      if fwrite(fid,tacslice(:,:,i), 'float32') ~= count
          fclose(fid);
          error(['coud not write file ' filename]);
      end
      fclose(fid);
    end
    
    % save reference perfusion maps
    filename = fullfile(phantompath, 'cbf');
    fid = fopen(filename,'ab');
    if fid == -1
        error(['error opening ' filename]);
    end
    if fwrite(fid,cbfslice, 'float32') ~= count
      fclose(fid);
      error(['coud not write file ' filename]);
    end
    fclose(fid);
    filename = fullfile(phantompath, 'cbv');
    fid = fopen(filename,'ab');
    if fid == -1
        error(['error opening ' filename]);
    end
    if fwrite(fid,cbvslice, 'float32') ~= count
      fclose(fid);
      error(['coud not write file ' filename]);
    end
    fclose(fid); 
    filename = fullfile(phantompath, 'mtt');
    fid = fopen(filename,'ab');
    if fid == -1
        error(['error opening ' filename]);
    end
    if fwrite(fid,mttslice, 'float32') ~= count
      fclose(fid);
      error(['coud not write file ' filename]);
    end
    fclose(fid);
    filename = fullfile(phantompath, 'ttp');
    fid = fopen(filename,'ab');
    if fid == -1
        error(['error opening ' filename]);
    end
    if fwrite(fid,ttpslice, 'float32') ~= count
      fclose(fid);
      error(['coud not write file ' filename]);
    end
    fclose(fid);
end
end



%--------------------------------------------------------------------------
% Compute the tissue time-attenuation curve using the indicator
% dilution theory.
%--------------------------------------------------------------------------
function tac = GetTissueTac(aif,t,dt,mtt,cbf)
cbf         = cbf/6000;              % ml/100ml/min = 1/(100*60s)
tmin        = mtt; 
r           = cbf *exp( -(t-mtt)); 
r(t<tmin)   = cbf;                   %   ...residue function
tac         = dt * conv(aif,r);      % indicator-dilution theory
tac         = tac(1:numel(aif));     % trf should have same size as aif
keyboard
end


% % Patlak model
% function tac = GetTissueTac(aif,t,dt,mtt,cbf)
% cbf         = cbf/6000;              % ml/100ml/min = 1/(100*60s)
% tmin        = mtt; 
% r           = cbf *exp( -(t-mtt)); 
% r(t<tmin)   = cbf;                   %   ...residue function
% tac         = dt * conv(aif,r);      % indicator-dilution theory
% tac         = tac(1:numel(aif));     % trf should have same size as aif
% end

