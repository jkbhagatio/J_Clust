function [Fs, uV_conversion, var_outs]  = load_sig(load_data_contents, load_data_val)

%Description: This .mfile loads tetrode data into the main J_Clust2 GUI. It allows the user to select from three different forms of data: 
%1) the raw signal from each tetrode channel, 2) the filtered signal from each tetrode channel, 3) the pre-detected spike waveforms and their
%respective timestamps. 
%
%Input: 'load_data_contents' & 'load_data_val' = user-defined (from GUI drop-down menu) form of data (either 1), 2), or 3) mentioned above)
%
%Output: 'Fs' = sampling rate, 'uV_conversion' = number for converting arbitraty units to uV value, 'var_outs' = depending on type of data file
%loaded - 1)filtered signal OR 1)spike waveforms, 2)spike timestamps, 3)number of spikes, 4)sample number of start of each spikes 5)number of 
%samples in spike waveform,  6)number of overlapped spike waveforms
%

[filename, filepath] = uigetfile;
if ~filepath
    Fs = [];
    uV_conversion = [];
    var_outs = cell(1,1);
    return;
end

addpath(genpath(filepath));
user_data = load(filename);
data_fieldnames = fieldnames(user_data);
 
%User system info input
user_sys_info_title = 'Recording System Info';
user_sys_info_prompt = {'Enter Sampling Rate (in Hz)', '(Optional) Enter Voltage Range (V)):','(Optional) Enter Number of Bits in ADC:', '(Optional) Enter gain (amplification):'};
user_sys_info = inputdlg(user_sys_info_prompt, user_sys_info_title, [1 50]);
Fs = str2num(user_sys_info{1}); 
if isempty(user_sys_info{2}) || isempty(user_sys_info{3}) || isempty(user_sys_info{4})
    uV_conversion = 1;
else
    adc_range = str2num(user_sys_info{2}); adc_bits = str2num(user_sys_info{3}); adc_gain = str2num(user_sys_info{4});
    uV_conversion = adc_range / (adc_bits * adc_gain) * 10^6; %convert adc values to uV
end

%load/calculate necessary variables from file

switch load_data_contents{load_data_val}
    case 'Raw Signal'
        raw_sig = user_data.(data_fieldnames{1});
        disp('Filtering Signal...')
        filt_sig = filter_tetrode(raw_sig, Fs);
        filt_sig = filt_sig * uV_conversion;
        var_outs{1} = filt_sig;  
        disp('Done')
    case 'Filtered Signal'
        filt_sig = user_data.(data_fieldnames(1));
        filt_sig = filt_sig * uV_conversion;
        var_outs{1} = filt_sig;
        disp('Filtered Signal loaded.')
    case 'Spike Waveforms'
        waveforms = user_data.(data_fieldnames{1});
        if size(waveforms,1) == 4
            ts = user_data.(data_fieldnames{2});
        else
            ts = user_data.(data_fieldnames{1});
            waveforms = user_data.(data_fieldnames{2});
        end
        if abs(min(waveforms(1,:,1))) > max(waveforms(1,:,1))
            waveforms = -waveforms;
        end
        num_spks = length(ts);
        new_spk_mrkr = ts * Fs;
        num_samples = size(waveforms,2);
        overlaps1 = find(diff(ts * 1000) < 1.5); %spks less than 1.5 ms apart,
        overlaps2 = unique([overlaps1; overlaps1 + 1]);  %merge to count both 1st and 2nd spks in an overlap, but remove duplicates
        var_outs{1} = waveforms; var_outs{2} = ts; var_outs{3} = num_spks; var_outs{4} = new_spk_mrkr; var_outs{5} = num_samples; var_outs{6} = overlaps2;
        disp('Spike Waveforms loaded')
end
