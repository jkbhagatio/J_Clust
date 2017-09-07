function [waveforms, ts, num_spks, new_spk_mrkr, num_samples, overlaps2, threshold] = spk_detection(filt_sig, Fs, uV_conversion, threshold)

%Description: This .mfile loads spike waveforms from a filtered tetrode signal. Using a novel approach to spike detection, it is able to capture
%multiple spike events within one spike window, which most online acquistion systems cannot do.
%
%Input: 'filt_sig' = filtered tetrode signal, 'Fs' = sampling rate, 'uV_conversion' = number for converting arbitraty units to uV value,
%'threshold' = calculated threshold above which spikes are detected(set as input argument because once threshold has been set, it should not change,
%even if user changes the "Set Time")
%
%Output: 'waveforms' = spike waveforms, 'ts' = spike timestamps,
%'num_spks' = number of spikes, 'new_spk_mrkr' = sample number where spike is detected, 'num_samples' = number of samples in spike waveform (1.5 ms),
%'overlaps2' = number of temporal spike waveform overlaps (within 1.5 ms), 'threshold' = calculated threshold above which spikes are detected
%

num_samples = ceil(Fs / 1000 * 1.5); %1.5 ms per spike waveform
filt_sig = filt_sig * uV_conversion;

if abs(min(filt_sig(:))) > max(filt_sig(:))
    filt_sig = -filt_sig;
end

if isempty(threshold)
    threshold = 5.5 * mean(median(abs(filt_sig / .6745),2));
end
[~, spk_datapoints] = find((filt_sig)>threshold);
spk_datapoints = unique(spk_datapoints); %all sample points in dataset with values > threshold
diff_spk_datapoints = diff(spk_datapoints);
new_spk_mrkr = spk_datapoints([1; find(diff_spk_datapoints > 1) + 1]); %1 to include first value in extra_datapoints, +1 for indexing purposes
num_spks = length(new_spk_mrkr);

pre_peak_samples = floor(1/3 * num_samples - 1);
post_peak_samples = ceil(2/3 * num_samples);
waveforms = zeros(4, num_samples, num_spks);

parfor i = 1:4
    for j = 2:num_spks-1 %throw out first and last spike detected in set in case full waveform can't be detected
        [peak, peakindx] = findpeaks(filt_sig(i, new_spk_mrkr(j):new_spk_mrkr(j)+pre_peak_samples),'NPeaks', 1); 
        if isempty(peakindx) || peak < filt_sig(i, new_spk_mrkr(j)) %this waveform starts off at peak
            peakindx = 1;
        end
        waveforms(i,:,j) = filt_sig(i, (new_spk_mrkr(j)+ peakindx - pre_peak_samples):new_spk_mrkr(j) + peakindx + post_peak_samples); %aligns to peak
    end
end

num_spks = num_spks - 2;
new_spk_mrkr(1) = []; new_spk_mrkr(end) = []; %throw out first and last spike detected in set in case full waveform can't be detected
waveforms(:,:,1) = []; waveforms(:,:,end) = [];

ts = new_spk_mrkr / Fs; %in seconds

overlaps1 = find(diff(ts * 1000) < 1.5); %spks less than 1.5 ms apart,    
overlaps2 = unique([overlaps1; overlaps1 + 1]);  %merge to count both 1st and 2nd spks in an overlap, but remove duplicates

end