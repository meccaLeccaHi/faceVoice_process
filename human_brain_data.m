% human_brain_data.m
%
% begins plotting responses from tdt data and header files
%
% apj
% last modified
% 1/5/17
%%%%%%%%%%%%%%%%

%% set constants
% patient info
PATIENT = '352L';
BLOCK = '007';

% header file from ptb
HEADERFILE = 'header_082416_1607.csv';

% figure print resolution
PRS = 150;

% length of time before/after stimulus event to retrieve
FLNKTIM = 500;

%% set directories
PROJDIR = '/home/lab/Cloud2/movies/human/LazerMorph/';
DATADIR = '/mnt/hbrl2/PetkovLab/Lazer_Morph/';
PREPROCDIR = [DATADIR PATIENT '/procData/'];
RESULTSDIR = [DATADIR PATIENT '/results/'];

eeg_fname = dir([PREPROCDIR '*eeglab.set']);
EEG = pop_loadset([PREPROCDIR eeg_fname.name]);

chan_num = length(plot_resp(:,1));

freq_range = [3 100];%[1.2 100]
maxfreq = max(freq_range);
padratio = 2;     
alpha_val = 0.05; % prob. level for two-tailed bootstrap
% elec = 1;
% basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
% basemax = 0;

% average responses on each channels
plot_resp = nanmean(EEG.data,3)';

% plot mean resp on all channels
figure
plot(plot_resp)
title(['all ' num2str(chan_num) ' channels'])

% save figure
save_name = [RESULTSDIR 'erp/aveResp_' PATIENT '_allChan.png'];
save_plot(save_name)

for pl = 1:chan_num
    
%     tic

    % define limits for y-axis based on limits of the current channel
    EEGchan = pop_select(EEG,'channel',pl);
    yLimsChan = 1.1.*[min(min(EEGchan.data)) max(max(EEGchan.data))];

    %% build header
    hd.voice = cell2mat({EEG.epoch.VOICE});
    hd.face = cell2mat({EEG.epoch.FACE});
    hd.noise = cell2mat({EEG.epoch.NOISE});
    hd.levels = cell2mat({EEG.epoch.LEVEL});
    hd.idents = cell2mat({EEG.epoch.IDENTITY});
    hd.traj = cell2mat({EEG.epoch.TRAJ});
    
    %% subset data
    EEGchan.face = pop_select(EEG,'channel',pl,'trial',find(hd.face));
    EEGchan.voice = pop_select(EEG,'channel',pl,'trial',find(hd.voice));
    EEGchan.both = pop_select(EEG,'channel',pl,'trial',find(hd.face&hd.voice));
    EEGchan.noise = pop_select(EEG,'channel',pl,'trial',find(hd.noise));

    %% plot mean erp responses for each modality
    meanResp.face = nanmean(EEGchan.face.data,3);
    meanResp.voice = nanmean(EEGchan.voice.data,3);
    meanResp.both = nanmean(EEGchan.both.data,3);
    
    % remove singletons
    meanResp = structfun(@squeeze,meanResp, 'UniformOutput', false);
    
    keyboard
    
%     plotArray = {meanResp_audOnly; meanResp_visOnly; meanResp_audVis};
%     plotArray = {meanResp_face; meanResp_voice; meanResp_both};

    plotArray = cell2mat(struct2cell(meanResp));
    ttlArray = fieldnames(meanResp);

    yLims = 1.1.*[min(min(plotArray)) max(max(plotArray))];

    figure('Color','w'); % ,'Visible','off'
    for i = 1:3
        subplot(3,1,i)
        plot(plotArray(i,:))
        xlim([0 length(plotArray(i,:))])
        ylim(yLims)
        
        % indicate stimulus onset/offset
        line([FLNKTIM FLNKTIM*2; FLNKTIM FLNKTIM*2],[yLims' yLims'],'LineWidth',2,'Color','k')
        
        title(ttlArray{i})
    end
    
    % add notes to corner of figure
    uicontrol('Style','text','String',[PATIENT '_chan' sprintf('%03d',pl)],...
        'Units','normalized','Position',[0 0 .175 .03]); 
    
    save_name = [RESULTSDIR 'erp/aveMod_' PATIENT '_chan' sprintf('%03d',pl) '.png'];
%     save_plot(save_name)
    
    
    %% plot time-frequency plots
 

    % identity levels
    levels = unique(hd.levels);
    loop_lev = levels(2:end);
    level_len = length(loop_lev);
    REPORT.levels = levels; % add to report
    
    % identity numbers
    idents = unique(hd.idents);
    loop_id = idents(2:end);
    ident_len = length(loop_id);
    REPORT.identityNum = ident_len; % add to report
    
    % trajectory direction
    traj_inds = hd.traj==1;
    
    % set up figure
    rows = ident_len;
    cols = level_len*3 + 2;
    fHan_erp = figure('Units','normalized','Outerposition',[0 0 .7 .7],'Color','w'); % ,'Visible','off'
    fHan_spect = figure('Units','normalized','Outerposition',[0 0 .7 .7],'Color','w'); % ,'Visible','off'

%     mod_title = {'AUD','VIS','A/V'};
    ttlArray = fieldnames(meanResp);
        
    for i = 1:ident_len % identity loop
        
        for ii = 1:level_len % identity level loop
            
            ident_inds = hd.idents==loop_id(i); % create indices for identity i
            level_inds = hd.levels==loop_lev(ii); % create indices for identity level ii
            
            loop_inds = level_inds&ident_inds&traj_inds; % create logical array for current loop
            
            %% subset data - modality>noise>level
            EEGstim.face = pop_select(EEG,'trial',find(hd.face&~hd.noise&loop_inds),'channel',pl);
            EEGstim.voice = pop_select(EEG,'trial',find(hd.voice&~hd.noise&loop_inds),'channel',pl);
            EEGstim.noise = pop_select(EEG,'trial',find(hd.face&~hd.noise&loop_inds),'channel',pl);
            
            % convert structure to cell array
            resps = struct2cell(EEGstim);
            
%                 resp{1} = dat.resps(ind.voice&~ind.face&~ind.noise&loop_inds,:);
%                 resp{2} = dat.resps(~ind.voice&ind.face&~ind.noise&loop_inds,:);
%                 resp{3} = dat.resps(ind.voice&ind.face&~ind.noise&loop_inds,:);
            
            for iii = 1:length(ttlArray)
                
                %% plot erp
                set(0,'CurrentFigure',fHan_erp)
                subplot(rows,cols,ii + ((i-1)*cols) + ((iii-1)*(level_len+1)))
                plot_resp = nanmean(squeeze(resps{iii}.data),2); % average over stimuli
%                 plot_resp = plot_resp(isfinite(resps{1}.data)); % remove NaNs
                plot(plot_resp)
                %             axis off
                set(gca,'XLim',[1 length(plot_resp)],'YLim',yLimsChan,...
                    'Position',get(gca,'Position').*[1 1 1.1 1]) % set position
                
                 % indicate stimulus onset/offset
                 line([FLNKTIM FLNKTIM*2; FLNKTIM FLNKTIM*2],[yLimsChan' yLimsChan'],...
                     'LineWidth',2,'Color','k')
        
                %% set text labels
                if i~=ident_len||ii~=1
                    set(gca,'XLabel',[],'YLabel',[],'XTick',[],'YTick',[])
                else
                    set(gca,'XTick',[FLNKTIM numel(plot_resp)-FLNKTIM])
                    if iii==1
                        xlabel('ms')
                        uicontrol('Style','text','String',[PATIENT '_chan' sprintf('%03d',pl)],...
                            'Units','normalized','Position',[0 0 .08 .02]);
                    end
                end
                if i==1&&iii==1
                    title(loop_lev(ii))
                end
                if ii==1&&iii==1
                    ylabel(loop_id(i),'rot',0,'FontWeight','Bold')
                end
                if i==1&&ii==1
                    % plot second title
                    text(2300,2,ttlArray{iii},'FontSize',14,'FontWeight','Bold', ...
                        'HorizontalAlignment','Center');
                end
                
                %% plot spectrogram
%                 plot_resp = gpuArray(plot_resp);
                set(0,'CurrentFigure',fHan_spect)
                subplot(rows,cols,ii + ((i-1)*cols) + ((iii-1)*(level_len+1)))
                
                keyboard
                
                %%%%%%%%%%%%%%%%%% LEFT OFF HERE
                 %%%% REMEMBER TO REMOVE NANS BEFORE RUNNING NEWTIMEF
                
                 [ersp,itc,powbase,times,freqs] = ...
                     newtimef(plot_resp(isfinite(plot_resp)),sum(isfinite(plot_resp)),...
                     [1 numel(plot_resp)-FLNKTIM],EEG.srate,0,'plotitc','off', ...
                     'maxfreq',maxfreq,'freqs',freq_range,'padratio', padratio, ...
                     'plotphase', 'off',  'alpha', alpha_val);

                 
                 %%%% maybe start plotting the output of newtimef() below
                 %%%% to eliminate extraneuous labels...
        
%                 segmentLength = round(numel(plot_resp)/4.5);
%                 [y,f,t,p] = spectrogram(plot_resp,round(segmentLength/5), ...
%                     round(80/100*segmentLength/5),[],dat.fs,'yaxis','MinThreshold',-20);
%                 surf(t*1000,f,10*log10(abs(p)),'EdgeColor','none');
%                 axis xy; axis tight; colormap(jet); view(0,90); 
%                 ylim([0 200]) 
%                 whitebg(2,'k')
                
%% FIX>?
%                 % plot stim time
%                 line([FLNKTIM numel(plot_resp)-FLNKTIM],[max(ylim) max(ylim)]+6,'LineWidth',2.5,...
%                     'Color','k','clipping','off');
                
                % set text labels
                if i~=ident_len||ii~=1
                    set(gca,'XLabel',[],'YLabel',[],'XTick',[],'YTick',[])
                else
                    set(gca,'XTick',[FLNKTIM numel(plot_resp)-FLNKTIM])
                    xlabel('ms')
                    uicontrol('Style','text','String',[PATIENT '_chan' sprintf('%03d',pl)],...
                        'Units','normalized','Position',[0 0 .08 .02]);
                end
                if i==1&&iii==1
                    title(loop_lev(ii),'Color','k')
                end
                if ii==1&&iii==1
                    ylabel(loop_id(i),'rot',0,'FontWeight','Bold')
                end
                if i==1&&ii==1
                    % plot second title  320
                    text(2300,diff(ylim)*1.5,ttlArray{iii},'FontSize',14','FontWeight','Bold', ...
                        'HorizontalAlignment','Center','Color','k');
                end
                
                % set position
                set(gca,'Position',get(gca,'Position').*[1 1 1.05 1])                
            end
        end
    end
    
    save_name = [RESULTSDIR 'erp/erp_' PATIENT '_chan' sprintf('%03d',pl) '.png'];
    save_plot(save_name, fHan_erp)
    
    save_name = [RESULTSDIR 'spect/spect_' PATIENT '_chan' sprintf('%03d',pl) '.png'];
    save_plot(save_name, fHan_erp, PRS)
    
%     toc
end

% rsTime = 0:max(tdt_data.Inpt_RZ2_chn001.t)/(length(tdt_data.Inpt_RZ2_chn001.dat)-1):max(tdt_data.Inpt_RZ2_chn001.t);
