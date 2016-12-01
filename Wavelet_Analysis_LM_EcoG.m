%Wavelet_Analysis_EcoG.m
%2015-3-18 by  YK

%This script runs Time-Frequency decomposition and ITC on the data using EEGLAB
%wavelet command. For human EcoG data, lower alpha_val might be better to
%use (monkey, use alpha_val=0.01).


%From newtimef.m
% Outputs:
%            ersp   = (nfreqs,timesout) matrix of log spectral diffs from baseline
%                     (in dB log scale or absolute scale). Use the 'plot' output format
%                     above to output the ERSP as shown on the plot.
%            itc    = (nfreqs,timesout) matrix of complex inter-trial coherencies.
%                     itc is complex -- ITC magnitude is abs(itc); ITC phase in radians
%                     is angle(itc), or in deg phase(itc)*180/pi.
%          powbase  = baseline power spectrum. Note that even, when selecting the
%                     the 'trialbase' option, the average power spectrum is
%                     returned (not trial based). To obtain the baseline of
%                     each trial, recompute it manually using the tfdata
%                     output described below.
%            times  = vector of output times (spectral time window centers) (in ms).
%            freqs  = vector of frequency bin centers (in Hz).
%         erspboot  = (nfreqs,2) matrix of [lower upper] ERSP significance.
%          itcboot  = (nfreqs) matrix of [upper] abs(itc) threshold.
%           tfdata  = optional (nfreqs,timesout,trials) time/frequency decomposition
%                      of the single data trials. Values are complex.

clear all
clc;



Linuxformat = 1; %0: Windows; 1: Linux

num_lev = [1:8]; % [1:16] or 99 Stimulus numbers. Put 99 to analyse all sound response together.


alpha_val = 0.05; %compute two-tailed bootstrap s, ignificance prob. level.
freq_range = [3 100];%[1.2 100]

%The figure only shows 500ms
basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
basemax = 0;

epochmin = -1; %start trial relative to sound onset (s)
epochmax = 4; %Set 5 for humans, 4 for monkeys end trial relative to sound onset (s)

plot_marktimes = [0 413 563 976 1126 1539 1689 2102 2253 2666];



maxfreq = max(freq_range);

%% INPUT---------------------------------------------

% PathName=uigetdir('Choose a directory');
% cd (PathName);
dirnames = uipickfiles( 'Prompt','Pick directorie(s) that includes Level mat files');
ndir = size(dirnames,2);

for jj = 1:ndir
    %Prep for a directory for results
    loaddir = char(dirnames(jj));
    cd(loaddir);
    
    if Linuxformat == 1   
        savdir = [loaddir '/alpha' num2str(alpha_val) 'maxfreq' num2str(maxfreq) '_baseline' num2str(basemin) '-' num2str(basemax)];
     elseif Linuxformat == 0
        savdir = [loaddir '\alpha' num2str(alpha_val) 'maxfreq' num2str(maxfreq) '_baseline' num2str(basemin)  '-' num2str(basemax)];
    end
    mkdir(savdir);
    
    %Parameter settings (fixed for AGL EcoG data)    
    num_channels = 1;
    srate = 1000;
    padratio = 2;     
    Elecs = 1;
    maxersp = 6; %db in power
    outtimes = 1500; %200 for data points in time for cohrence and power analysis %Original 1500
    
    itcplot = 'on'; % on for ITC plots

    
    for level = 1:length(num_lev)
        if num_lev(level)<10
            load(['snd0' num2str(num_lev(level)) '.mat']);
        elseif num_lev(level)>=10
            load(['snd' num2str(num_lev(level)) '.mat']);
        end
        
        
        data =X(Elecs,1:(abs(epochmin)+epochmax)*1000,:);
        %         data=diff(data,1,2);
        
        elec=1;
        b=size(data,3); %b = number of trials
        in(:,1)=[1:b]';
        in(:,2)=abs(epochmin)*ones(b,1); % if I change here, it's stuck,. duration of the baseline before sound onset (sec)
        
        eeglab
        EEG = pop_importdata( 'dataformat', 'array', 'data', 'data', 'setname', 'Level', 'srate',srate, 'pnts',0, 'xmin',0, 'nbchan',0);
        EEG = eeg_checkset( EEG );
        EEG = pop_importepoch( EEG, 'in', { 'Epoch', 'stim'}, 'latencyfields',{ 'stim'}, 'timeunit',1, 'headerlines',0);
        EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  'stim'  }, [epochmin         epochmax], 'newname', 'Level epochs', 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        EEG = pop_rmbase( EEG, [basemin    0]);
        EEG = eeg_checkset( EEG );
        clear  X_C delay_data_comb delay_data_ind data in temp
        close all
        %set(gcf,'Position',[1 1 1600 1000]);
        
        
        %Figure settings-------------------------------------
        title([loaddir  'stim' num2str(num_lev(level))], 'Fontsize', 6);
        fprintf('\n\n');
        fprintf('Stim number  %d\n', num_lev(level))
        fprintf('\n\n');
        %fprintf('Electrode number (of 16) %d\n', elec)    
        
        
        [ersp,itc,powbase,times,freqs,erspboot,itcboot,alltfX] = pop_newtimef(EEG, ...
            1, elec, [EEG.xmin EEG.xmax]*srate,[3 0.5], 'maxfreq',maxfreq, 'freqs',freq_range,'padratio', padratio, ...
            'plotphase', 'off', 'timesout', outtimes, 'alpha', alpha_val, 'naccu', 300, 'baseboot',1,'rmerp','off', ...
            'erspmax', maxersp, 'plotersp','on', 'plotitc',itcplot,'baseline',[basemin basemax],'marktimes',plot_marktimes);
        
        save([savdir  filesep 'snd_' num2str(num_lev(level)) ],'ersp','itc','powbase','freqs','times','erspboot','itcboot','alltfX',...
            'X','preevetime','postevetime');

        clear alltfX times freqs erspboot itcboot ersp itc powbase X preevetime postevetime
        saveas(gcf, [ savdir filesep num2str(jj) '-snd' num2str(num_lev(level)) ], 'png')

    end
end

