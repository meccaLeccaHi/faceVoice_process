% human_brain_data.m
%
% begins plotting responses from tdt data and header files
%
% apj
% last modified
% 12/13/16
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



% alpha_val = 0.05; %compute two-tailed bootstrap s, ignificance prob. level.
% freq_range = [3 100];%[1.2 100]
% maxfreq = max(freq_range);
% 
% padratio = 2;     
% 
% elec = 1;
% 
% %The figure only shows 500ms
% basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
% basemax = 0;
% 
% 
% 
% %  [ersp,itc,powbase,times,freqs,erspboot,itcboot,alltfX] = pop_newtimef(EEG, ...
% %             1, elec, [EEG.xmin EEG.xmax]*EEG.srate,[3 0.5], 'maxfreq',maxfreq, 'freqs',freq_range,'padratio', padratio, ...
% %             'plotphase', 'off', 'timesout', outtimes, 'alpha', alpha_val, 'naccu', 300, 'baseboot',1,'rmerp','off', ...
% %             'erspmax', maxersp, 'plotersp','on', 'plotitc',itcplot,'baseline',[basemin basemax],'marktimes',plot_marktimes);       
% 
% % [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
% %     newtimef(EEG.data, EEG.pnts, [EEG.xmin EEG.xmax], EEG.srate,0);
% 
% 
% plot_resp = nanmean(EEG.data,3);
% 
% figure
% plot(plot_resp)
% 
% % eeglab
% 
% figure;
% [ersp,itc,powbase,times,freqs] = ...
%     newtimef(plot_resp,length(plot_resp),[-FLNKTIM numel(plot_resp)-FLNKTIM], ...
%     EEG.srate,0,'plotitc','off');
                
                
                

preProcList = dir([PREPROCDIR PATIENT(1:3) '*chan*.mat']); % get list of pre-processed files
ppNameList = sort_nat({preProcList.name}'); % sort according to channel number

for pl = 1:length(ppNameList)
    
%     tic
    
    % load preProc file
    dat = load([PREPROCDIR ppNameList{pl}]);

    % grab header info
    ind.voice = cell2mat(dat.header(:,ismember(dat.headNms,'VOICE')));
    ind.face = cell2mat(dat.header(:,ismember(dat.headNms,'FACE')));
    ind.noise = cell2mat(dat.header(:,ismember(dat.headNms,'NOISE')));
    ind = structfun(@logical,ind, 'UniformOutput', false);

    hd.levels = cell2mat(dat.header(:,ismember(dat.headNms,'LEVEL')));
    hd.idents = cell2mat(dat.header(:,ismember(dat.headNms,'IDENTITY')));
    hd.traj = cell2mat(dat.header(:,ismember(dat.headNms,'TRAJ')));

    % grab responses
    subdat.voOnly = dat.resps(ind.voice&~ind.face&~ind.noise,:);
    subdat.faOnly = dat.resps(ind.face&~ind.voice&~ind.noise,:);
    subdat.voFa = dat.resps(ind.face&~ind.voice&~ind.noise,:);
%     subdat.voNoise = dat.resps(ind.voice&ind.noise,:);
%     subdat.faNoise = dat.resps(ind.face&ind.noise,:);
%     subdat.voFaNoise = dat.resps(ind.face&ind.voice&ind.noise,:);


    %% build header
    clear hd
    hd.voice = cell2mat({EEG.epoch.VOICE});
    hd.face = cell2mat({EEG.epoch.FACE});
    hd.noise = cell2mat({EEG.epoch.NOISE});

    % grab indices of each group of stimuli
    ind = structfun(@find,ind, 'UniformOutput', false);

    % add remaining variables to header
    hd.levels = cell2mat({EEG.epoch.LEVEL});
    hd.idents = cell2mat({EEG.epoch.IDENTITY});
    hd.traj = cell2mat({EEG.epoch.TRAJ});
    
    %% subset data
    EEGchan.face = pop_select(EEG,'channel',pl,'trial',ind.face);
    EEGchan.voice = pop_select(EEG,'channel',pl,'trial',ind.voice);
    EEGchan.both = pop_select(EEG,'channel',pl,'trial',hd.face&hd.voice);
    EEGchan.noise = pop_select(EEG,'channel',pl,'trial',ind.noise);

    %% plot mean erp responses for each modality
        meanResp_audOnly = nanmean(subdat.voOnly);
        meanResp_visOnly = nanmean(subdat.faOnly);
        meanResp_audVis = nanmean(subdat.voFa);

    meanResp_face = nanmean(EEGchan.face.data);
    meanResp_voice = nanmean(EEGchan.voice.data);
    meanResp_both = nanmean(EEGchan.both.data);
    
    keyboard
    
%     plotArray = {meanResp_audOnly; meanResp_visOnly; meanResp_audVis};
    plotArray = {meanResp_face; meanResp_voice; meanResp_both};

    ttlArray = {'face'; 'voice'; 'both'};
    
    yLims = 1.1*[min(min(dat.resps)) max(max(dat.resps))];

    figure('Color','w','Visible','off'); 
    for i = 1:3
        subplot(3,1,i)
        plot(plotArray{i})
        xlim([0 1000])
        ylim(yLims)
        line([FLNKTIM FLNKTIM],yLims,'LineWidth',2,'Color','k')
        line([FLNKTIM*3 FLNKTIM*3],yLims,'LineWidth',2,'Color','k')
        title(ttlArray{i})
    end
    uicontrol('Style','text','String',[PATIENT '_chan' sprintf('%03d',pl)],...
        'Units','normalized','Position',[0 0 .175 .03]); 
    
    save_name = [RESULTSDIR 'erp/aveMod_' PATIENT '_chan' sprintf('%03d',pl) '.png'];
    export_fig(save_name)
    fprintf('\nsaved: %s',save_name);
    close(gcf)
    
    %% plot time-frequency plots
    levels = unique(hd.levels);
    loop_lev = levels(2:end);
    level_len = length(loop_lev);
    idents = unique(hd.idents);
    loop_id = idents(2:end);
    ident_len = length(loop_id);
    traj_inds = hd.traj==1;
    
    rows = ident_len;
    cols = level_len*3 + 2;
    fHan_erp = figure('Units','normalized','Outerposition',[0 0 .7 .7],'Color','w','Visible','off'); % ,'Visible','off'
    fHan_spect = figure('Units','normalized','Outerposition',[0 0 .7 .7],'Color','w','Visible','off'); % ,'Visible','off'

    mod_title = {'AUD','VIS','A/V'};
    for i = 1:ident_len
        for ii = 1:level_len
            
            ident_inds = hd.idents==loop_id(i);
            level_inds = hd.levels==loop_lev(ii);

            resp{1} = dat.resps(ind.voice&~ind.face&~ind.noise&...
                level_inds&ident_inds&traj_inds,:);
            resp{2} = dat.resps(~ind.voice&ind.face&~ind.noise&...
                level_inds&ident_inds&traj_inds,:);
            resp{3} = dat.resps(ind.voice&ind.face&~ind.noise&...
                level_inds&ident_inds&traj_inds,:);
            
            for iii = 1:length(mod_title)
                %% plot erp
                set(0,'CurrentFigure',fHan_erp)
                subplot(rows,cols,ii + ((i-1)*cols) + ((iii-1)*(level_len+1)))
                plot_resp = resp{iii}(isfinite(resp{1})); % remove NaNs
                plot(plot_resp)
                %             axis off
                set(gca,'XLim',[1 length(plot_resp)],'YLim',yLims,...
                    'Position',get(gca,'Position').*[1 1 1.1 1]) % set position
                
                % plot stim time
                line([FLNKTIM FLNKTIM],yLims,'LineWidth',1,'Color','k')
                line([length(plot_resp) length(plot_resp)]-FLNKTIM,yLims,'LineWidth',1,'Color','k')
                
                % set text labels
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
                    text(2300,max(ylim)*2,mod_title{iii},'FontSize',14','FontWeight','Bold', ...
                        'HorizontalAlignment','Center');
                end
                
                %% plot spectrogram
%                 plot_resp = gpuArray(plot_resp);
                set(0,'CurrentFigure',fHan_spect)
                subplot(rows,cols,ii + ((i-1)*cols) + ((iii-1)*(level_len+1)))
                
                keyboard
                eeglab
                figure;
                [ersp,itc,powbase,times,freqs]=...
                    newtimef(plot_resp,length(plot_resp),[-FLNKTIM numel(plot_resp)-FLNKTIM],dat.fs,...
                    0,'plotitc','off');
                
                EEG = pop_importdata( 'dataformat', 'array', 'data', 'data', 'setname', 'Level', 'srate',srate, 'pnts',0, 'xmin',0, 'nbchan',0);


        
                segmentLength = round(numel(plot_resp)/4.5);
                [y,f,t,p] = spectrogram(plot_resp,round(segmentLength/5), ...
                    round(80/100*segmentLength/5),[],dat.fs,'yaxis','MinThreshold',-20);
                surf(t*1000,f,10*log10(abs(p)),'EdgeColor','none');
                axis xy; axis tight; colormap(jet); view(0,90); 
                ylim([0 200]) 
                whitebg(2,'k')
                
                % plot stim time
                line([FLNKTIM numel(plot_resp)-FLNKTIM],[200 200]+6,'LineWidth',2.5,...
                    'Color','k','clipping','off');
                
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
                    text(2300,diff(ylim)*1.5,mod_title{iii},'FontSize',14','FontWeight','Bold', ...
                        'HorizontalAlignment','Center','Color','k');
                end
                
                % set position
                set(gca,'Position',get(gca,'Position').*[1 1 1.05 1])                
            end
        end
    end
    
    save_name = [RESULTSDIR 'erp/erp_' PATIENT '_chan' sprintf('%03d',pl) '.png'];
    export_fig(fHan_erp,save_name,'-nocrop')
    fprintf('\nsaved: %s',save_name);
    close(fHan_erp)
    
    save_name = [RESULTSDIR 'spect/spect_' PATIENT '_chan' sprintf('%03d',pl) '.png'];
    export_fig(fHan_spect,save_name,'-nocrop',['-r' num2str(PRS)])
    fprintf('\nsaved: %s',save_name);
    close(fHan_spect)
    
%     toc
end

% rsTime = 0:max(tdt_data.Inpt_RZ2_chn001.t)/(length(tdt_data.Inpt_RZ2_chn001.dat)-1):max(tdt_data.Inpt_RZ2_chn001.t);
