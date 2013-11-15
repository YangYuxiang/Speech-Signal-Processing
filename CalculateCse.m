function cse = CalculateCse(x, fs)
%  This function computes cochlea-scaled spectral entropy of speech.
%  input:
%  x    vector of the speech 
%  fs   sampling frequency of the speech 
%  output:
%  cse     corresponding vector of CSE values
%
%  Functions required 
%   1. enframe.m
%   2. MakeERBFilters.m
%   3. ERBFilterBank.m
%   4. rfft.m
%   
%   Reference
%   Stilp, Christian E., and Keith R. Kluender. "Cochlea-scaled entropy,
%   not consonants, vowels, or time, best predicts speech intelligibility." 
%   Proceedings of the National Academy of Sciences 107.27 (2010): 12387-12392.
%
%    Written by Yuxiang Yang (dreamspark4ever@gmail.com) Nov. 2013
%

numChannels = 33;
vcDuration = 5;
nFFT = 2^nextpow2(0.016*fs);

x = x/norm(x);
frames = enframe(x, 0.016*fs);

numFrames = size(frames, 1);
mag = zeros(numChannels*((nFFT+2)/2), numFrames);
eucDist = zeros(1, numFrames-1);
numBoxcars = numFrames-1-(vcDuration-1);
cse = zeros(1, numBoxcars);

fcoefs = MakeERBFilters(fs, numChannels, 0);
for i = 1:numFrames
    y = ERBFilterBank(frames(i,:), fcoefs);
    yf = rfft(y, nFFT, 2);
    ym = abs(yf);
    temp = ym';
    mag(:,i) = temp(:);
end

for j = 1:numFrames-1
    eucDist(j) = norm(mag(:, j)-mag(:, j+1));
end

for k = 1:numBoxcars
    cse(k) = sum(eucDist(k:k+vcDuration-1));
end
