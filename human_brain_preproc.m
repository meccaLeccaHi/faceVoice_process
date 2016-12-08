% human_brain_preproc.m
%
% begins parsing tdt data file and header into final data format
%
% apj
% last modified
% 12/8/16
%%%%%%%%%%%%%%%%
tic

%% read header file produced by ptb
pp.special = '/mnt/hbrl2/PetkovLab/Lazer_Morph/';
pp.printout = [pp.special '352L/results/trigs/'];
ptb_data = table2array(readtable([pp.special '352L/SPECIAL_mat/header_082416_1607.csv'],...
    'ReadVariableNames',false));

% split into header vars
headNms = ptb_data(1,:);
header = ptb_data(2:end,:);
header(cellfun(@isempty,header)) = {1}; % fill in empty 'SCALE' values

movNms = header(:,ismember(headNms,'MOVIE_NAME'));
varNms = {'VOICE' 'FACE' 'NOISE' 'LEVEL' 'IDENTITY' 'TRAJ'}; % define variable names

%% read physio file
tdt_data = load([pp.special '352L/SPECIAL_mat/352-007_SPECIALevents_DBT1.mat']);

% remove unneccessary fields
fields = fieldnames(tdt_data);
channels = fields(~cellfun(@isempty,regexp(fields,'LFPx*')));

% read analog photodiode signal, derive onset/offset times of stimuli
photodiode = tdt_data.Inpt_RZ2_chn002.dat;
flnkTm = 500;
minCross = find(photodiode > -0.06); % threshold crossings
IsoCross = find(diff(minCross)> 250); % only keep crossings longer than .25 secs
Starts = minCross(IsoCross(3:end)); % vector of indices (not times)
Stops = minCross(IsoCross(3:end)+1); % vector of indices (not times) 
maxTimeLen = max(Stops-Starts)+2*flnkTm; % for later pre-allocation
% figure; hist(Stops-Starts)
% mkdir('/mnt/hbrl2/PetkovLab/Lazer_Morph/352L/procData/')

% extract sampling rate (length = 1)
fs = tdt_data.(channels{1}).fs(1);      
        
% save dat arrays
for j = 1:length(channels)
    
    chan_str = ['LFPx_RZ2_chn' sprintf('%03d',j)];
    dat.fs = tdt_data.(chan_str).('fs')(1);

    dat.resps = nan(length(Starts),maxTimeLen+1); % preallocate space
    data = tdt_data.(chan_str).('dat');
    for i = 1:length(Starts) % trial loop
        trialSnip = data(Starts(i)-flnkTm:Stops(i)+flnkTm);
        dat.resps(i,1:length(trialSnip)) = trialSnip;
    end

    if j==1
        figure('Color','w','Visible','off')
        smp_inc = max(tdt_data.Inpt_RZ2_chn001.t)/(length(tdt_data.Inpt_RZ2_chn001.dat)-1);
        rsTime = 0:smp_inc:max(tdt_data.Inpt_RZ2_chn001.t);
        plot(rsTime,tdt_data.Inpt_RZ2_chn001.dat,'k'); hold on
        plot(rsTime,tdt_data.Inpt_RZ2_chn002.dat,'b');
        plot(rsTime(Starts),zeros(size(Starts))-0.05,'.g');
        plot(rsTime(Stops),zeros(size(Stops))-0.05,'.r');
        xlim([52 80])
        export_fig([pp.printout 'san_check_chan' sprintf('%03d',j) '.png'])
        close(gcf)
    end
    
    vars = zeros(length(movNms),length(varNms)); % preallocate variable matrix

    for i = 1:length(movNms)
        [~,name,~] = fileparts(movNms{i}); % get name substring
        
        % extract stim info
        if any(strfind(name,'aud'))
            vars(i,1) = 1;
        end
        if any(strfind(lower(name),'vi'))
            vars(i,2) = 1;
        end
        if any(strfind(name,'noisy'))
            vars(i,3) = 1;
        end
        
        identNum = regexp(name,'\d*','Match'); % identity

        %%%%% ADD RAD V TAN VARIABLE TO VARS BELOW
        
        % extract morph level, identity, & trajectory
        name = strsplit(name,'_');
        if regexp(name{1},'Average')
            vars(i,4) = 0;
            vars(i,5) = 0;
            vars(i,6) = 0;
        else
            vars(i,4) = str2double(name{2}); % morph level
            vars(i,5) = str2double(identNum{1}); % identity
            if regexp(name{1}(end-2:end),'rad'); % trajectory
                vars(i,6) = 1;
            else
                vars(i,6) = 2;
            end
        end
    end
    
    dat.headNms = [headNms varNms];
    dat.header	= [header num2cell(vars)];

    save([pp.special '352L/procData/352-007_chan' num2str(j) '.mat'],'-struct','dat')
end

keyboard

%% left off here
%% just need to wrap up this 3d array below

tic
eeglab_mat = [];

% set up eeglab 3D array - channels, timepoints, epochs
for s = 1:length(Starts) % epoch loop
    
    epoch_mat = nan(length(channels),maxTimeLen+1); % pre-allocate
    
    % create array for epoch N - rows are channels and columns are data points
    for c = 1:length(channels) % channel loop
        volt_data = tdt_data.(channels{c}).dat(Starts(s)-flnkTm:Stops(s)+flnkTm);
        epoch_mat(c,1:length(volt_data)) = volt_data';
    end
    
    % iteratively concatenate epoch arrays to create 3d array
    eeglab_mat = cat(3,eeglab_mat,epoch_mat);
    
end
toc

% import data into EEGlab array
EEG = pop_importdata('dataformat','array','data','eeglab_mat','setname','LazerMorph',...
    'srate',fs,'pnts',length(eeglab_mat),'nbchan',length(channels),'xmin',0);
% import epochs into EEGlab array
EEG = pop_importepoch(EEG, dat.header, dat.headNms);

% % figure out what to do with each of these later
toc