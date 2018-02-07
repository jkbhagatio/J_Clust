function [waveforms, ts, num_spks, new_spk_mrkr, num_samples, overlaps2, threshold] = spk_detection(filt_sig, Fs, threshold)

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

%check to see if signal is oriented with most spikes having negative peaks - if not, flip polarity of signal

sort_sig = sort(filt_sig(:), 'ascend');
negs = sort_sig(1:100);
poss = sort_sig(end-99:end);

if abs(mean(negs)) < mean(poss) %then flip polarity of signal
    filt_sig = -filt_sig;
end

if isempty(threshold)
    threshold(1) = -5 * mean(median(abs(filt_sig / .6745),2)); %negative threshold
    threshold(2) = 5 * mean(median(abs(filt_sig / .6745), 2)); %positive threshold
end

[~, spk_datapoints1] = find((filt_sig)<threshold(1));
new_spk_mrkr1 = spk_datapoints1([1; find(diff(spk_datapoints1) > 1) + 1]);

new_spk_mrkr1_long = zeros(num_samples+1, length(new_spk_mrkr1));
for i =1:length(new_spk_mrkr1)
    new_spk_mrkr1_long(:,i) = [new_spk_mrkr1(i):new_spk_mrkr1(i)+num_samples];
end
new_spk_mrkr1_long = new_spk_mrkr1_long(:);

filt_sig2 = filt_sig;
filt_sig2(:,new_spk_mrkr1_long) = 0;
[~, spk_datapoints2] = find((filt_sig2) > threshold(2));
new_spk_mrkr2 = spk_datapoints2([1; find(diff(spk_datapoints2) > 1) + 1]);

new_spk_mrkr = sort([new_spk_mrkr1; new_spk_mrkr2]);
num_spks = length(new_spk_mrkr);

pre_peak_samples = floor(1/3 * num_samples - 1);
post_peak_samples = ceil(2/3 * num_samples);
waveforms = zeros(4, num_samples, num_spks);

parfor i = 1:4
    for j = 2:num_spks-1 %throw out first and last spike detected in set in case full waveform can't be detected
        deriv = sign(gradient(filt_sig(i, new_spk_mrkr(j):new_spk_mrkr(j)+pre_peak_samples)));
        peak_indx = find(deriv == -deriv(1), 1);
        if isempty(peak_indx) || abs(filt_sig(i, new_spk_mrkr(j)+peak_indx)) < abs(filt_sig(i, new_spk_mrkr(j))) %this waveform starts off at peak
            peak_indx = 1;
        end
        waveforms(i,:,j) = filt_sig(i, (new_spk_mrkr(j)+ peak_indx - pre_peak_samples):new_spk_mrkr(j) + peak_indx + post_peak_samples); %aligns to peak
    end
end

%throw out first and last spike detected in set in case full waveform can't be detected
num_spks = num_spks - 2;
new_spk_mrkr(1) = []; new_spk_mrkr(end) = []; 
waveforms(:,:,1) = []; waveforms(:,:,end) = [];

ts = new_spk_mrkr / Fs; %in seconds

overlaps1 = find(diff(ts * 1000) < 1.5); %spks less than 1.5 ms apart,    
overlaps2 = unique([overlaps1; overlaps1 + 1]);  %merge to count both 1st and 2nd spks in an overlap, but remove duplicates

end