% Extraction_LM_Ecog.m
% 12/7/2016 APJ


%This script extracts mat file in response to each sound and also generate
%the overall average ERP figure. It also extracts sound wave forms when
%audio_in=1. Adapted from 'Extraction_AGL_Ecog.m'


% clc
% close all
% clear all
% 
% linux_on = 1;

% set constants
DATADIR = '/home/lab/Cloud2/movies/human/LazerMorph/';
PATIENT = '352L';
BLOCK = '007';


% connection template
connectSumm = readtable(fullfile(DATADIR, PATIENT, ...
    [PATIENT '_conTemp.csv']));

% convert channels from strings of indices to literal vectors
Channel = cell(length(connectSumm{:,'TDTAcqChannel'}),1);
for i = 1:length(connectSumm{:,'TDTAcqChannel'})
    Channel{i} = eval(cell2mat(connectSumm{i,'TDTAcqChannel'}));
end
connectSumm = [connectSumm table(Channel)];  % concatenate variables


%%%%%%%%%%%%%%S
% 
% %L307_________________________________________
% HG = [184:191];
% FC =[193:224];
% ATL = [161:183];
% %STG = [65:153];
% STG = [100:153];
% aIns = [23:30];
% pIns = [33:40];
% Amg = [41:48];
% pHC = [49:60];
% pHCstrip = [7:16];
% %subTandFr [17:20;225:232]

%% INPUT
% % if linux_on
% %      datadir ='/mnt/hbrl2/PetkovLab/AGL/348R/';
% % %     datadir ='/home/windows_shares/lcn_fileserver/projects/Iowa/Newcastle/AGL_EcoG/data/2015-03-09/L307_AGLdata/';
% % else
% %     datadir ='V:\lcn\projects\Iowa\Newcastle\AGL_EcoG\data\2015-03-09\L307_AGLdata\';
% % end

% % %datafn = '348-014_SPECIALevents_DBT1';%data mat file
% % datafn = '348-014_SPECIALevents';%data mat file

% data file name
data_fname = [strjoin(regexp(PATIENT,['\d'],'match'),'') '-' BLOCK '_SPECIALevents'];

channels = [65:99];
audio_in = 0;


trial_on = 0; %Turn on figures for each trial

preevetime=1;
postevetime=5;

%stimID=[1:16];%Test
stimID=[1:8]; %Exposure

%%

cd ([datadir])
    
for ch = 1:length(channels)

    % fn=sprintf('LFPx_RZ2_chn%d', channels(ch));
    
    
    if audio_in
        fn=['Inpt_RZ2_chn00' num2str(channels(ch))];
    else
        if channels(ch)>=100
            fn=['LFPx_RZ2_chn' num2str(channels(ch))];
        elseif 10<=channels(ch) & channels(ch)<100
            fn=['LFPx_RZ2_chn0' num2str(channels(ch))];
        elseif channels(ch) < 10
            fn=['LFPx_RZ2_chn00' num2str(channels(ch))];
        end
    end
    
    display(['Processing '  fn]);
    %%
    %Loading events and data for each contact
    matObj=load(data_fname,'Epoc','Evnt', fn);
    
    data_dat = matObj.(fn).dat;
    data_t = matObj.(fn).t;%1 to 1010
    srate = matObj.(fn).fs(1);
    
    
    
    
    %%
    %     Evnt.time;%diff is 5.0738 s Inter-event-interval?
    %     Evnt.evnt; %1 to 16 (sequence ID)& 255 at the end; 10 repeats
    
    mkdir([datafn '/' fn]);
    cd ([datafn '/' fn]);
    
    figure;
    
    
    for stim_no=1:length(stimID)
        evnt=matObj.Evnt.evnt(matObj.Evnt.evnt==stim_no);
        evnt_time=matObj.Evnt.time(matObj.Evnt.evnt==stim_no);
        
        for ii=1:length(evnt_time) %ntrials
            %             if audio_in
            %                 t_index=int32(evnt_time(ii)*srate)-int32(preevetime*srate):int32(evnt_time(ii)*srate)+int32(postevetime*srate);
            %             else
            %t_index=int32(evnt_time(ii)*srate)-preevetime*srate:1000/srate:int32(evnt_time(ii)*srate)+postevetime*srate;
            t_index=int32(evnt_time(ii)*srate)-int32(preevetime*srate):int32(evnt_time(ii)*srate)+int32(postevetime*srate);
            %end
            t=0:1/srate:preevetime+postevetime; % in sec
            
            if ii==1
                X=NaN(1,length(t),ii);
            end
            
            
            %Check below
            X(1,1:length(t),ii)=data_dat(t_index);
            x=squeeze(X);
            
            if trial_on
                subplot(5,2,ii);
                plot(t-preevetime,x(:,ii));
                xlim([preevetime preevetime+postevetime]);
                title(num2str(ii));
            end
        end
        
        subplot(4,4,stim_no);
        
        plot(t-preevetime,nanmean(x,2));
        
        if stim_no==1
            title(['Grand average ERP: stim' num2str(stim_no)]);
            xlabel('Time from sequence onset(s)');
        else
            title(['stim' num2str(stim_no)]);
        end
        xlim([-preevetime postevetime]);
        
        if stim_no<10
            save(['snd0' num2str(stim_no) '.mat'],'X', 'postevetime','preevetime','stim_no','evnt','evnt_time','datafn','channels','fn');
        else
            save(['snd' num2str(stim_no) '.mat'],'X', 'postevetime','preevetime','stim_no','evnt','evnt_time','datafn','channels','fn');
        end
        clear X x
    end %end of stim_no
    saveas(gcf, 'ERP.png');
    close all
    clear data_dat data_t fn
end % end of channels

cd ([datadir])

clear all
