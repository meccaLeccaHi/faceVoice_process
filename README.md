# faceVoice_process
#### Matlab scripts for analysis of human intra-cranial recordings
These scripts includes elements of EEGlab (https://sccn.ucsd.edu/eeglab/).

### human_brain_preproc.m
- De-noised digital signals are read in and parsed into trial epochs.
- Examples of trigger timing are produced for debugging.

### human_brain_data.m
- Parsed data is read iteratively to produce figures for each recording channel.
- Data from all channels are summarized in a report. 
