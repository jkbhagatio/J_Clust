function [daub8_wavelets] = MyDWT(waveforms)

%where 'waveforms' = 2-d array (number of spikes x length of spike waveform),
%'n' = number of filter bank levels, %'wname' = mother wavelet
%(see 'wfilters' function) 'lf' = low-pass (aka rough/approximate) filter
%coefficients, 'hf' = high-pass(aka detail) filter coefficients

%Run In form of 
%1) MyDWT(...)
%or 2) MyDWT(... n) *wname will default to 'db4'
%or 3) MyDWT(... 'wname') *n will default to 5
%or 4) MyDWT(... n, 'wname')
%or 5) MyDWT(... n, lf, hf) -> allows for specification of filter
%coefficients without using a specific mother wavelet

dwtmode2('per');

num_spks = size(waveforms,1);
num_samples = size(waveforms, 2);

l_daub8 = wmaxlev(num_samples, 'db4');

num_coeffs = length(wavedec(waveforms(1,:), l_daub8, 'db4'));
daub8_wavelets = zeros(num_spks, num_coeffs);

for i = 1:num_spks
    daub8_wavelets(i,:) = wavedec(waveforms(i,:), l_daub8, 'db4');
end

[std_coeffs, std_coeffs_indx] = sort(std(daub8_wavelets), 'descend');

daub8_wavelets(:,:) = daub8_wavelets(:, std_coeffs_indx);
    
end




