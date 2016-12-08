% human_brain_data.m
%
% begins plotting responses from tdt data and header files
%
% apj
% last modified
% 12/7/16
%%%%%%%%%%%%%%%%

patient = '352L';
flnkTm = 500;
prs = 150;
pp.special = '/mnt/hbrl2/PetkovLab/Lazer_Morph/';
pp.preProc = [pp.special patient '/procData/'];
pp.cloud = '/home/lab/Cloud2/movies/human/LazerMorph/';
pp.aveMod = [pp.cloud patient '/results/aveMod/'];
pp.erp = [pp.cloud patient '/results/erp/'];
pp.spect = [pp.cloud patient '/results/spect/'];

preProcList = dir([pp.preProc patient(1:3) '*.mat']); % get list of pre-processed files
ppNameList = sort_nat({preProcList.name}'); % sort according to channel number

for pl = 1:length(ppNameList)
    
%     tic
    
    % load preProc file
    dat = load([pp.preProc ppNameList{pl}]);

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

 
    %% plot mean erp responses for each modality
    meanResp_audOnly = nanmean(subdat.voOnly);
    meanResp_visOnly = nanmean(subdat.faOnly);
    meanResp_audVis = nanmean(subdat.voFa);
    
    plotArray = {meanResp_audOnly; meanResp_visOnly; meanResp_audVis};
    ttlArray = {'audio'; 'visual'; 'audio/visual'};
    yLims = 1.1*[min(min(dat.resps)) max(max(dat.resps))];

    figure('Color','w','Visible','off'); 
    for i = 1:3
        subplot(3,1,i)
        plot(plotArray{i})
        xlim([0 1000])
        ylim(yLims)
        line([flnkTm flnkTm],yLims,'LineWidth',2,'Color','k')
        line([flnkTm*3 flnkTm*3],yLims,'LineWidth',2,'Color','k')
        title(ttlArray{i})
    end
    uicontrol('Style','text','String',[patient '_chan' sprintf('%03d',pl)],...
        'Units','normalized','Position',[0 0 .175 .03]); 
    
    save_name = [pp.aveMod 'aveMod_' patient '_chan' sprintf('%03d',pl) '.png'];
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
                line([flnkTm flnkTm],yLims,'LineWidth',1,'Color','k')
                line([length(plot_resp) length(plot_resp)]-flnkTm,yLims,'LineWidth',1,'Color','k')
                
                % set text labels
                if i~=ident_len||ii~=1
                    set(gca,'XLabel',[],'YLabel',[],'XTick',[],'YTick',[])
                else
                    set(gca,'XTick',[flnkTm numel(plot_resp)-flnkTm])
                    if iii==1
                        xlabel('ms')
                        uicontrol('Style','text','String',[patient '_chan' sprintf('%03d',pl)],...
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
                    newtimef(plot_resp,length(plot_resp),[-flnkTm numel(plot_resp)-flnkTm],dat.fs,...
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
                line([flnkTm numel(plot_resp)-flnkTm],[200 200]+6,'LineWidth',2.5,...
                    'Color','k','clipping','off');
                
                % set text labels
                if i~=ident_len||ii~=1
                    set(gca,'XLabel',[],'YLabel',[],'XTick',[],'YTick',[])
                else
                    set(gca,'XTick',[flnkTm numel(plot_resp)-flnkTm])
                    xlabel('ms')
                    uicontrol('Style','text','String',[patient '_chan' sprintf('%03d',pl)],...
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
    
    save_name = [pp.erp 'erp_' patient '_chan' sprintf('%03d',pl) '.png'];
    export_fig(fHan_erp,save_name,'-nocrop')
    fprintf('\nsaved: %s',save_name);
    close(fHan_erp)
    
    save_name = [pp.spect 'spect_' patient '_chan' sprintf('%03d',pl) '.png'];
    export_fig(fHan_spect,save_name,'-nocrop',['-r' num2str(prs)])
    fprintf('\nsaved: %s',save_name);
    close(fHan_spect)
    
%     toc
end

% rsTime = 0:max(tdt_data.Inpt_RZ2_chn001.t)/(length(tdt_data.Inpt_RZ2_chn001.dat)-1):max(tdt_data.Inpt_RZ2_chn001.t);
