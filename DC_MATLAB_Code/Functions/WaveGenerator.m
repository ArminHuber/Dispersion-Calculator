function [FrequencyRange,ExcitationMagnitude,ExcitationSignal,SampleRate,FourierTransformLength] = WaveGenerator(Frequency,FrequencyResolution,Cycles,SamplesPerCycle,Chirp,ChirpType,Window,WidthFactor)
N = 0:min(Frequency)/max(Frequency)/SamplesPerCycle:Cycles;
f = Frequency*1e3;
if  Chirp % generate the carrier wave
    if  ChirpType == 1 % linear
        t = N/mean(f);
        Y = cos(2*pi*(f(1)+diff(f)/t(end)/2*t).*t);
        % y = chirp(t,f(1),t(end),f(2),'linear');
    elseif ChirpType == 2 % quadratic
        t = 3*N/(2*f(1)+f(2));
        Y = cos(2*pi*(f(1)+diff(f)/t(end)^2/3*t.^2).*t);
        % y = chirp(t,f(1),t(end),f(2),'quadratic');
    elseif ChirpType == 3 % exponential
        t = N*(log(f(2))-log(f(1)))/(f(2)-f(1));
        k = f(2)/f(1);
        Y = cos(2*pi*f(1)*t(end)/log(k)*(k.^(t/t(end))-1));
        % y = chirp(t,f(1),t(end),f(2),'logarithmic');
    end
% figure,plot(t,Y),figure,plot(t,y)
else
    t = N/f;
    Y = cos(2*pi*N);
end
n = length(N)-1;
if  Window == 1 % Gaussian window function
    H = exp(-2*(WidthFactor*linspace(-.5,.5,n+1)).^2);
    % h = gausswin(n+1,WidthFactor);
elseif Window == 2 % Tukey
    x = 0:round(WidthFactor*(n+1)/2)-1;
    H = ones(1,n+1);
    H(1:length(x)) = (1-cos(2*pi*x/n/WidthFactor))/2;
    H(end-length(x)+1:end) = fliplr(H(1:length(x)));
    % h = tukeywin(n+1,WidthFactor);
elseif Window == 3 % Hann (equal to Tukey with WidthFactor = 1)
    H = sin(pi*(0:n)/n).^2;
    % h = hann(n+1);
elseif Window == 4 % Hamming
    H = .54-.46*cos(2*pi*(0:n)/n);
    % h = hamming(n+1);
elseif Window == 5 % Triangular
    H = linspace(0,1,n/2+1);
    H = [H fliplr(H(1:end-1))];
    % h = triang(n+1);
elseif Window == 6 % Rectangular
    H = ones(1,n+1);
    % h = rectwin(n+1);
end
% figure,plot(H),figure,plot(h)
ExcitationSignal = [t;Y.*H]; % time range of the excitation signal (s); apply the selected window function to the carrier wave to generate the wave packet
SampleRate = round(1/t(2)); % (Hz)
FourierTransformLength = round(max(Frequency)/FrequencyResolution*SamplesPerCycle); % determines the spectral resolution
if  mod(FourierTransformLength,2)
    FourierTransformLength = FourierTransformLength+1;
end
FrequencyRange = SampleRate/1e3*(0:FourierTransformLength/2)/FourierTransformLength; % the frequency range for the FFT (kHz)
Spectrum = fft(ExcitationSignal(2,:),FourierTransformLength)/FourierTransformLength; % compute the two-sided frequency spectrum
Spectrum(FourierTransformLength/2+2:end) = []; % compute the single-sided frequency spectrum
ExcitationMagnitude = [FrequencyRange;2*abs(Spectrum);angle(Spectrum)]; % frequencies, spectral amplitudes, and phases; multiply amplitudes by two to account for the skipped negative frequencies