function psd=perogram(num_psd,window,overlap,seg_size,T,x,y)

% Classical periodogram auto/cross power spectral density estimate by Welch
% procedure, eqs. (5.39),(5.40),(5.43).                              Version 1.0
%
%    Auto:   psd=perogram(num_psd,window,overlap,seg_size,T,x)
%    Cross:  psd=perogram(num_psd,window,overlap,seg_size,T,x,y)
%
% num_psd  -- number of psd vector elements (must be power of 2);
%             frequency spacing between elements is F = 1/(num_psd*T)
%             Hz, with psd(1) corresponding to frequency  -1/(2*T)
%             Hz, psd(num_psd/2 + 1) corresponding to  0  Hz, and
%             psd(num_psd) corresponding to 1/(2*T) - F  Hz.
% window   -- window selection: 0 -- none, 1 -- Hamming, 2 -- Nuttall
% overlap  -- number of overlap samples between segments
% seg_size -- number of samples per segment (must be even)
% T        -- sample interval in seconds
% x        -- vector of data samples
% y        -- vector of data samples (cross PSD only);
%                 note that length(y) must equal length(x)
% psd      -- vector of num_psd auto/cross PSD values

if nargin == 7
    if length(x) ~= length(y)
        error('The input data vectors are not of equal lengths.')
    end
end
if seg_size > length(x)
    error(['Seg_size cannot exceed num_psd = ',int2str(num_psd)])
end
if (overlap < 0) | (overlap >= seg_size)
    error(['Out of range: 0 <= overlap <',int2str(seg_size)])
end
shift = seg_size - overlap;
num_segs = fix(( length(x) - seg_size )/shift) + 1;
psd = zeros(num_psd,1);
range = 1:seg_size;
s = seg_size - 1;
ph = 2*pi*(0:s)/s;

if     window == 0  wind = ones(seg_size,1);          % rectangular (aka boxcar)
elseif window == 1  wind = (.53836 - .46164*cos(ph))';                 % Hamming
elseif window == 2  wind=(.42323-.49755*cos(ph)+.07922*cos(2*ph))';    % Nuttall
else
    error('Window selection number is invalid.')
end
% Wind  could also be specified as one of the following from the MATLAB Signal
% Processing Toolbox:
%   wind = barlett(seg_size);
%   wind = blackman(seg_size);
%   wind = chebwin(seg_size,sidelobe level in db);
%   wind = hanning(seg_size);
%   wind = kaiser(seg_size,0.1102*(sidelobe level in db - 8.7));
%   wind = triang(seg_size);

for k=1:num_segs
    z = T*fft(wind.*x(range),num_psd);
    if nargin < 7
        psd = psd + real(z).^2 + imag(z).^2;
    else
        psd = psd + z.*(T*conj(fft(wind.*y(range),num_psd)));
    end
    range = range + shift;
end
pow_wind = sum(wind.^2)/seg_size;          % adjustment in gain for window power
psd = (1/(num_segs*pow_wind*seg_size*T))*fftshift(psd);
%
