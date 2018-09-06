function[B] = get_polar_pattern_from_FIR_wav(wavfile,n_fft,steering_vectors,fs)
% steering_vectors: [nFreq,nChans,nDirections]


% load in the wav file
[in,fs_in] = audioread(wavfile);
if fs_in~=fs
    error('Expected fs of input file to be %s Hz',num2str(fs));
end

if n_fft<size(in,1)
    error('fft is smaller than FIR length')
end

H = rfft(in,n_fft,1);                   %[nFreq,nChans]
B = sum(bsxfun(@times,H,steering_vectors),2);     %[nFreq,1,nDirections]

