function [waveforms, ts, num_spks, new_spk_mrkr, num_samples, overlaps2, threshold] = spk_detection(filt_sig, Fs, uV_conversion, threshold)

%Description: This .mfile loads spike waveforms from a filtered tetrode signal. Using a novel approach to spike detection, it is able to capture
%multiple spike events within one spike window, and both positive peak and negative peak spikes (which can both occur depending on location of 
%electrode relative to particular point of a neuron in extracellular space), which most online acquistion systems cannot do. 

%A negative peak value is more common and is clasically indicative of current flowing away from an electrode (presumably from extracellular space -> 
%intracellular space). A positive peak value is uncommon but may be indicative of an electrode placed in an area that is subbject to a dipole-like 
%phenomenon (e.g. an area near interesection of apical dendrites and soma at a moment when the soma is summating input currents), wherein a 
%return current essentially flips the polarity of the signal (spike).

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

%check to see if signal is oriented with most spikes having negative peaks - if not, flip polarity of signal

negs = sort(filt_sig(:), 'ascend');
poss = sort(filt_sig(:), 'descend');

if abs(mean(negs(1:100))) < mean(poss(1:100)) %then flip polarity of signal
    filt_sig = -filt_sig;
end

if isempty(threshold)
    threshold(1) = -5.5 * mean(median(abs(filt_sig / .6745),2)); %negative threshold
    threshold(2) = 5.5 * mean(median(abs(filt_sig / .6745), 2)); %positive threshold
end

[~, spk_datapoints1] = find((filt_sig)<threshold(1));
[~, spk_datapoints2] = find((filt_sig)>threshold(2));
spk_datapoints1 = unique(spk_datapoints1);
spk_datapoints2 = unique(spk_datapoints2);

parfor i=1:length(spk_datapoints2)
    close_updownspk_vec = sort(abs(spk_datapoints2(i) - spk_datapoints1));
    if close_updownspk_vec(1) < num_samples %if this positive peak is actually due to the repolarization/hyperpolarized period of a negative peak spike
        spk_datapoints2(i) = -1; %mark this spike for later removal
    end
end

spk_datapoints2(spk_datapoints2 == -1) = []; %remove all previously marked "bad" positive peak spikes

spk_datapoints = sort([spk_datapoints1; spk_datapoints2]); %all sample points in dataset with |values| > threshold
diff_spk_datapoints = diff(spk_datapoints);
new_spk_mrkr = spk_datapoints([1; find(diff_spk_datapoints > 1) + 1]); %1 to include first value in extra_datapoints, +1 for indexing purposes
num_spks = length(new_spk_mrkr);
[~, pos_pk_spks] = intersect(new_spk_mrkr, spk_datapoints2);

pre_peak_samples = floor(1/3 * num_samples - 1);
post_peak_samples = ceil(2/3 * num_samples);
waveforms = zeros(4, num_samples, num_spks);

parfor i = 1:4
    for j = 2:num_spks-1 %throw out first and last spike detected in set in case full waveform can't be detected
        if ~isempty(intersect(j, pos_pk_spks)) %we are dealing with a positive peak spike
            [peak, peakindx] = findpeaks(filt_sig(i, new_spk_mrkr(j):new_spk_mrkr(j)+pre_peak_samples),'NPeaks', 1);
            if isempty(peakindx) || peak < filt_sig(i, new_spk_mrkr(j)) %this waveform starts off at peak
                peakindx = 1;
            end
        else
            [peak, peakindx] = findpeaks(-filt_sig(i, new_spk_mrkr(j):new_spk_mrkr(j)+pre_peak_samples),'NPeaks', 1);
            if isempty(peakindx) || peak > filt_sig(i, new_spk_mrkr(j)) %this waveform starts off at peak
                peakindx = 1;
            end
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