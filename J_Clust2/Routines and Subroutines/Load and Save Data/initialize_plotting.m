%Description: This .m file is a script which loads spike waveforms from the user dataset, calculates the spike features from the spike waveforms, 
%and plots the peak amplitudes of each spike on the two axes in the main GUI (channels axis and time axis)

%% Set sig samples based on start time

if ~handles.preload
    if handles.start_time == 0
        sig_sample_start = 1;
    else
        sig_sample_start = handles.Fs * handles.start_time;
    end
    sig_sample_end = handles.Fs * handles.end_time;
    if sig_sample_end > length(handles.filt_sig)
        sig_sample_end = length(handles.filt_sig);
        warning(['Last recorded sample is at ', num2str(length(handles.filt_sig) / handles.Fs), ' seconds']);
    end
else
    handles.first_spk = find(handles.ts > handles.start_time, 1);
    handles.last_spk = find(handles.ts > handles.end_time, 1);
    if isempty(handles.last_spk)
        handles.last_spk = length(handles.ts);
        warning(['Last recorded spike is at ', num2str(handles.ts(end)), ' seconds']);
    end   
end
%% Detect spikes and calculate features
if ~handles.preload
    disp('Loading Spike Waveforms...')
    [handles.waveforms, handles.ts, handles.num_spks, handles.new_spk_mrkr, handles.num_samples, handles.overlaps, handles.threshold] = ...
        spk_detection(handles.filt_sig(:,sig_sample_start:sig_sample_end), handles.Fs, handles.threshold);    
    disp('Calculating Spike Features...')
    handles.features = calc_features(handles.waveforms, handles.num_spks, handles.num_samples, handles.Fs);
    handles.ts = handles.ts + handles.start_time;
else
    disp('Calculating Spike Features...')
    handles.cur_ts = [handles.first_spk:handles.last_spk];
    handles.features = calc_features(handles.waveforms(:,:,handles.cur_ts), length(handles.cur_ts), handles.num_samples, handles.Fs);
    handles.cur_overlaps = handles.overlaps(handles.overlaps < length(handles.cur_ts));
end
    
disp('Done')
%% Initial plot of spike peak amps

handles.chan_disp1 = [1 2 3]; handles.feat_disp1 = 'Peak Amplitude';
handles.chan_disp2 = 1; handles.feat_disp2 = 'Peak Amplitude';
handles.feature_wires = handles.features{1};
handles.feature_time = handles.features{1};

plot_on_channel_scatter;
plot_on_time_scatter;

brush 'k' %sets default brushing color to black

