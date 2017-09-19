function filt_sig = filter_tetrode(raw_sig, Fs)

%Description: This .mfile takes a raw signal from a tetrode and filters it using adaptive, FIR filter-building, based on the sampling rate of the data
%
%Input: 'raw_sig' = raw tetrode signal, 'Fs' = sampling rate
%
%Output: 'filt_sig' = filtered tetrode signal
%

Wn = [400 5500] * 2 / Fs; %cutoff frequencies where attenuation = -3dB
trans_bw = 350 / Fs; %width of transition band as fraction of Fs
n = ceil(4/trans_bw); %filter length using window method for Hamming window (2 * Fs / (f_co_1 - f_sb_1)) ***((M = 4 / transition bandwidth percentage of sampling rate))
B = fir1(n, Wn); %build filter using Hamming window (default)

filt_sig = zeros(4, length(raw_sig));

parfor i = 1:4
    filt_sig(i,:) = filtfilt(B, 1, double(raw_sig(i,:))); % runs filter
end

end

