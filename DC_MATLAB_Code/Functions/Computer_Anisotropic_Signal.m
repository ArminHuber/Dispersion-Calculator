% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2025 DLR
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% =========================================================================
function [ExcitationMagnitude,ExcitationSignal,ExcitationSpectrum,FourierTransformLength,Frequency,FrequencyLimitLow,FrequencyLimitHigh,FrequencyRange,PlotXAxis,Gate,XAxis,uSum,uSumALamb,uSumAShear,uSumBLamb,uSumBShear,uSumSLamb,uSumSShear] = Computer_Anisotropic_Signal(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,ALamb,AShear,BLamb,BShear,c,Material,Mode,SLamb,SShear,Layers,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,Delta,SamplesPerCycle,Cycles,Distance,FrequencyResolution,Frequency,MultiMode,SpectrumThreshold,TimeLimitFactor,ALambModes,AShearModes,BLambModes,BShearModes,SLambModes,SShearModes,Window,OutputWindowUI3)
SamplesX3 = 50; % through the total thickness

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
Memory = memory;
SamplesX3 = ceil(SamplesX3/Layers); % samples per layer
Height = Layers*(SamplesX3+1);
if  ~ToggleUpperFluid
    UpperFluid.Velocity = 1e-10;
    UpperFluid.Density = 1e-10;
end
if  ~ToggleLowerFluid
    LowerFluid.Velocity = 1e-10;
    LowerFluid.Density = 1e-10;
end
for m = 1:SuperLayerSize
    if  ~Decoupled
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
        A1=0;
    else
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
        a11=0;a12=0;a23=0;a33=0;a34=0;
    end
end
XAxis = [];
if  mod(Cycles,2)
    Sign = -1;
else
    Sign = 1;
end
CarrierWave = Sign*cos(2*pi*(0:1/SamplesPerCycle:Cycles)); % generate the carrier wave
if  Window == 1 % Gaussian window function
    H = exp(-20.*(-.5:1/(SamplesPerCycle*Cycles):.5).^2); % Michel's window function from "SineSignalSpectrum_Michel.m"
elseif Window == 2 % Hann
    H = sin(pi*(0:(SamplesPerCycle*Cycles))/(SamplesPerCycle*Cycles)).^2;
elseif Window == 3 % Hamming
    H = 25/46-21/46*cos(2*pi*(0:(SamplesPerCycle*Cycles))/(SamplesPerCycle*Cycles));
elseif Window == 4 % Triangular
    H = horzcat(linspace(0,1,SamplesPerCycle*Cycles/2+1));
    H = horzcat(H,fliplr(H(1:end-1)));
elseif Window == 5 % Rectangular
    H(1:SamplesPerCycle*Cycles+1) = 1;
end
ExcitationSignal = 1/Frequency/1e3*(0:Cycles/(SamplesPerCycle*Cycles):Cycles); % time range of the excitation signal (s)
ExcitationSignal(2,:) = CarrierWave.*H; % apply the selected window function to the carrier wave to generate the wave packet
SampleRate = round(SamplesPerCycle*Frequency*1e3); % (Hz)
FourierTransformLength = round(Frequency/FrequencyResolution*SamplesPerCycle); % determines the spectral resolution
FrequencyRange = SampleRate/1e3*(0:FourierTransformLength/2)/FourierTransformLength; % the frequency range for the FFT (kHz)
ExcitationMagnitude = abs(fft(ExcitationSignal(2,:),FourierTransformLength))/FourierTransformLength; % compute the two-sided frequency spectrum
ExcitationMagnitude(FourierTransformLength/2+2:end) = []; % compute the single-sided frequency spectrum
ExcitationMagnitude = 2*ExcitationMagnitude; % multiply the spectral amplitudes by two to account for the skipped negative frequencies
z1 = find(ExcitationMagnitude > SpectrumThreshold/100*(max(ExcitationMagnitude)),1); % find the index where the spectral amplitudes start to be higher than the SpectrumThreshold
z2 = find(ExcitationMagnitude > SpectrumThreshold/100*(max(ExcitationMagnitude)),1,'last'); % find the index where the spectral amplitudes end to be higher than the SpectrumThreshold
if  FrequencyRange(z1) == 0
    z1 = 2;
end
ExcitationSpectrum = FrequencyRange(z1:z2); % only this frequency range is considered for the superposition of the wave components
ExcitationSpectrum(2,:) = ExcitationMagnitude(z1:z2); % the corresponding spectral amplitudes
CoherenceTime = 1e-3/FrequencyResolution; % (s) the coherence time tc is the time after which we have a repetition of the temporal response, i.e., we will see a twin of the wave packet at each instant of time n*tc, n=1,2,..., away from the actual wave packet    
FrequencyLimitLow = 0;
FrequencyLimitHigh = 0;
PhaseVelocityALamb{1} = 0;
PhaseVelocitySLamb{1} = 0;
PhaseVelocityBLamb{1} = 0;
PhaseVelocityAShear{1} = 0;
PhaseVelocitySShear{1} = 0;
PhaseVelocityBShear{1} = 0;
uSum{1} = [];
uSumALamb{1} = [];
uSumSLamb{1} = [];
uSumBLamb{1} = [];
uSumAShear{1} = [];
uSumSShear{1} = [];
uSumBShear{1} = [];
if  ~MultiMode
    p = str2double(regexp(Mode,'\d*','match'))+1;
    if  ~contains(Mode,'SH') && Mode(1) == 'S'
        [~,q1] = min(abs(ExcitationSpectrum(1,1)-SLamb{p}(:,1))); % get the index of the starting frequency of the calculated frequency spectrum
        z = find(isapprox(diff(SLamb{p}(:,1)),SLamb{p}(end,1)-SLamb{p}(end-1,1),'loose'),1,'first'); % find the index of the first frequency from a phase velocity sweep where we have equidistant frequencies; this is to exclude the high phase velocity part of higher order modes obtained by frequency sweeps where we don't have equidistant frequencies; the energy velocity of this range is set to zero and can be identified by that way
        if  isempty(z)
            z = 1;
        end
        if  q1 < z % if the starting frequency is in the non-equidistant frequency range of higher order modes, we take the first frequency from the phase velocity sweep
            q1 = z;
        end
        FrequencyLimitLow = SLamb{p}(z,1); % used to indicate in the spectral amplitude plot below which frequency no data are available for higher order modes; it is the first frequency from a phase velocity sweep where we have equidistant frequencies; we indicate the range below this frequency with red shading, similarly as the range above the top frequency in the dispersion diagram 
        FrequencyLimitHigh = SLamb{p}(end,1); % used to indicate above which frequency no data are available; this can happen in attenuated cases
        [~,z1] = min(abs(ExcitationSpectrum(1,:)-SLamb{p}(q1,1))); % find the index in the calculated frequency spectrum corresponding the starting frequency (q1) determined above; we will shorten the frequency spectrum for the superposition of the wave components below
        [~,q2] = min(abs(ExcitationSpectrum(1,end)-SLamb{p}(:,1))); % get the index of the ending frequency of the calculated frequency spectrum
        [~,z2] = min(abs(ExcitationSpectrum(1,:)-SLamb{p}(q2,1))); % find the index corresponding the ending frequency (q2) determined above
        PhaseVelocity{1} = SLamb{p}(q1:q2,4)*1e3; % extract the phase velocities of the requested mode in the frequency range in question
        Attenuation{1} = SLamb{p}(q1:q2,7).*PhaseVelocity{1}./(SLamb{p}(q1:q2,1)*1e3); % get the attenuation in Np/wavelength
        EnergyVelocity = SLamb{p}(SLamb{p}(:,1) == Frequency,5)*1e3; % get the energy velocity of the wave packet
        ModeType = '';
    elseif ~contains(Mode,'SH') && Mode(1) == 'A'
        [~,q1] = min(abs(ExcitationSpectrum(1,1)-ALamb{p}(:,1)));
        z = find(isapprox(diff(ALamb{p}(:,1)),ALamb{p}(end,1)-ALamb{p}(end-1,1),'loose'),1,'first');
        if  isempty(z)
            z = 1;
        end
        if  q1 < z
            q1 = z;
        end
        FrequencyLimitLow = ALamb{p}(z,1);
        FrequencyLimitHigh = ALamb{p}(end,1);
        [~,z1] = min(abs(ExcitationSpectrum(1,:)-ALamb{p}(q1,1)));
        [~,q2] = min(abs(ExcitationSpectrum(1,end)-ALamb{p}(:,1)));
        [~,z2] = min(abs(ExcitationSpectrum(1,:)-ALamb{p}(q2,1)));
        PhaseVelocity{1} = ALamb{p}(q1:q2,4)*1e3;
        Attenuation{1} = ALamb{p}(q1:q2,7).*PhaseVelocity{1}./(ALamb{p}(q1:q2,1)*1e3);
        EnergyVelocity = ALamb{p}(ALamb{p}(:,1) == Frequency,5)*1e3;
        ModeType = '';
    elseif ~contains(Mode,'SH') && Mode(1) == 'B'
        [~,q1] = min(abs(ExcitationSpectrum(1,1)-BLamb{p}(:,1)));
        z = find(isapprox(diff(BLamb{p}(:,1)),BLamb{p}(end,1)-BLamb{p}(end-1,1),'loose'),1,'first');
        if  isempty(z)
            z = 1;
        end
        if  q1 < z
            q1 = z;
        end
        FrequencyLimitLow = BLamb{p}(z,1);
        FrequencyLimitHigh = BLamb{p}(end,1);
        [~,z1] = min(abs(ExcitationSpectrum(1,:)-BLamb{p}(q1,1)));
        [~,q2] = min(abs(ExcitationSpectrum(1,end)-BLamb{p}(:,1)));
        [~,z2] = min(abs(ExcitationSpectrum(1,:)-BLamb{p}(q2,1)));
        PhaseVelocity{1} = BLamb{p}(q1:q2,4)*1e3;
        Attenuation{1} = BLamb{p}(q1:q2,7).*PhaseVelocity{1}./(BLamb{p}(q1:q2,1)*1e3);
        EnergyVelocity = BLamb{p}(BLamb{p}(:,1) == Frequency,5)*1e3;
        ModeType = '';
    elseif  contains(Mode,'SH') && Mode(1) == 'S'
        [~,q1] = min(abs(ExcitationSpectrum(1,1)-SShear{p}(:,1)));
        z = find(isapprox(diff(SShear{p}(:,1)),SShear{p}(end,1)-SShear{p}(end-1,1),'loose'),1,'first');
        if  isempty(z)
            z = 1;
        end
        if  q1 < z
            q1 = z;
        end
        FrequencyLimitLow = SShear{p}(z,1);
        FrequencyLimitHigh = SShear{p}(end,1);
        [~,z1] = min(abs(ExcitationSpectrum(1,:)-SShear{p}(q1,1)));
        [~,q2] = min(abs(ExcitationSpectrum(1,end)-SShear{p}(:,1)));
        [~,z2] = min(abs(ExcitationSpectrum(1,:)-SShear{p}(q2,1)));
        PhaseVelocity{1} = SShear{p}(q1:q2,4)*1e3;
        Attenuation{1} = SShear{p}(q1:q2,7).*PhaseVelocity{1}./(SShear{p}(q1:q2,1)*1e3);
        EnergyVelocity = SShear{p}(SShear{p}(:,1) == Frequency,5)*1e3;
        ModeType = 'Shear';
    elseif contains(Mode,'SH') && Mode(1) == 'A'
        [~,q1] = min(abs(ExcitationSpectrum(1,1)-AShear{p-1}(:,1)));
        z = find(isapprox(diff(AShear{p-1}(:,1)),AShear{p-1}(end,1)-AShear{p-1}(end-1,1),'loose'),1,'first');
        if  isempty(z)
            z = 1;
        end
        if  q1 < z
            q1 = z;
        end
        FrequencyLimitLow = AShear{p-1}(z,1);
        FrequencyLimitHigh = AShear{p-1}(end,1);
        [~,z1] = min(abs(ExcitationSpectrum(1,:)-AShear{p-1}(q1,1)));
        [~,q2] = min(abs(ExcitationSpectrum(1,end)-AShear{p-1}(:,1)));
        [~,z2] = min(abs(ExcitationSpectrum(1,:)-AShear{p-1}(q2,1)));
        PhaseVelocity{1} = AShear{p-1}(q1:q2,4)*1e3;
        Attenuation{1} = AShear{p-1}(q1:q2,7).*PhaseVelocity{1}./(AShear{p-1}(q1:q2,1)*1e3);
        EnergyVelocity = AShear{p-1}(AShear{p-1}(:,1) == Frequency,5)*1e3;
        ModeType = 'Shear';
    elseif contains(Mode,'SH') && Mode(1) == 'B'
        [~,q1] = min(abs(ExcitationSpectrum(1,1)-BShear{p}(:,1)));
        z = find(isapprox(diff(BShear{p}(:,1)),BShear{p}(end,1)-BShear{p}(end-1,1),'loose'),1,'first');
        if  isempty(z)
            z = 1;
        end
        if  q1 < z
            q1 = z;
        end
        FrequencyLimitLow = BShear{p}(z,1);
        FrequencyLimitHigh = BShear{p}(end,1);
        [~,z1] = min(abs(ExcitationSpectrum(1,:)-BShear{p}(q1,1)));
        [~,q2] = min(abs(ExcitationSpectrum(1,end)-BShear{p}(:,1)));
        [~,z2] = min(abs(ExcitationSpectrum(1,:)-BShear{p}(q2,1)));
        PhaseVelocity{1} = BShear{p}(q1:q2,4)*1e3;
        Attenuation{1} = BShear{p}(q1:q2,7).*PhaseVelocity{1}./(BShear{p}(q1:q2,1)*1e3);
        EnergyVelocity = BShear{p}(BShear{p}(:,1) == Frequency,5)*1e3;
        ModeType = 'Shear';
    end
    Gate = (Distance/1e3/EnergyVelocity+[-.1*CoherenceTime .2*CoherenceTime])*1e6;
    TimeLimit = round(TimeLimitFactor*Distance/1e3/EnergyVelocity,6); % (s) the time limit for which the temporal response is calculated
    if  TimeLimit-ExcitationSignal(1,SamplesPerCycle*Cycles/2+1) < ExcitationSignal(1,SamplesPerCycle*Cycles/2+1)
        Time = -ExcitationSignal(1,SamplesPerCycle*Cycles/2+1):1/SampleRate:ExcitationSignal(1,SamplesPerCycle*Cycles/2+1);
        PlotXAxis = ExcitationSignal(1,:)*1e6; % (microsec)
    else
        Time = -ExcitationSignal(1,SamplesPerCycle*Cycles/2+1):1/SampleRate:TimeLimit-ExcitationSignal(1,SamplesPerCycle*Cycles/2+1);
        PlotXAxis = (0:1/SampleRate:TimeLimit)*1e6; % (microsec)
    end
    if  CoherenceTime < TimeLimit % if the coherence time is smaller than the calculated time range, we have to expect seeing unwanted twins of the wave packet
        for i = -1:100
            FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity/CoherenceTime),i); % (kHz) calculate the frequency resolution necessary to get a coherence time larger than the displayed time range
            CoherenceTime2 = 1e-3/FrequencyResolution2; % the new coherence time (s)
            if  FrequencyResolution2 > 0 && CoherenceTime2 > TimeLimit
                break
            end
        end
    else
        FrequencyResolution2 = FrequencyResolution;
    end
    [uSum,~,ExcitationSpectrumRange] = Computer_Anisotropic_Signal_Core(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,CoherenceTime,0,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,0,ModeType,0,PhaseVelocity,Attenuation,Time,TimeLimit,z1,z2,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
    if  isempty(uSum{1})
        return
    end
    if  CoherenceTime < TimeLimit
        String = ['TIME DOMAIN:',newline,...
        'Time limit:           ',num2str(TimeLimit*1e6),' micsec',newline,...
        'Coh. time (original): ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline,...
        'Coh. time (interp.):  ',num2str(CoherenceTime2*1e6,'%.0f'),' micsec',newline,...
        'Sample rate:          ',num2str(SampleRate/1e3),' kHz',newline,...
        'Samples:              ',num2str(length(Time)),newline,newline,...
        'FREQUENCY DOMAIN:',newline,...
        'Spectral range:  ',num2str(ExcitationSpectrumRange{1}(1,1)),' - ',num2str(ExcitationSpectrumRange{1}(1,end)),' kHz',newline,...
        'Res. (original): ',num2str(FrequencyResolution),' kHz',newline,...
        'Res. (interp.):  ',num2str(FrequencyResolution2),' kHz',newline,...
        'Frequencies:     ',num2str(length(ExcitationSpectrumRange{1}))];
    else
        String = ['TIME DOMAIN:',newline,...
        'Time limit:     ',num2str(TimeLimit*1e6),' micsec',newline,...
        'Coherence time: ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline,...
        'Sample rate:    ',num2str(SampleRate/1e3,'%.0f'),' kHz',newline,...
        'Samples:        ',num2str(length(Time)),newline,newline,...
        'FREQUENCY DOMAIN:',newline,...
        'Spectral range: ',num2str(ExcitationSpectrumRange{1}(1,1)),' - ',num2str(ExcitationSpectrumRange{1}(1,end)),' kHz',newline,...
        'Resolution:     ',num2str(FrequencyResolution),' kHz',newline,...
        'Frequencies:    ',num2str(length(ExcitationSpectrumRange{1}))];
    end
else
    EnergyVelocityALamb = NaN(1,10);
    EnergyVelocityAShear = NaN(1,10);
    EnergyVelocityBLamb = NaN(1,10);
    EnergyVelocityBShear = NaN(1,10);
    EnergyVelocitySLamb = NaN(1,10);
    EnergyVelocitySShear = NaN(1,10);
    ExcitationSpectrumRangeLimits = NaN(6,2);
    for p = 1:length(ALambModes)
        if  ALambModes(p)
            [~,q1] = min(abs(ExcitationSpectrum(1,1)-ALamb{p}(:,1)));
            z = find(isapprox(diff(ALamb{p}(:,1)),ALamb{p}(end,1)-ALamb{p}(end-1,1),'loose'),1,'first');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            [~,z1ALamb(p)] = min(abs(ExcitationSpectrum(1,:)-ALamb{p}(q1,1))); 
            [~,q2] = min(abs(ExcitationSpectrum(1,end)-ALamb{p}(:,1)));
            [~,z2ALamb(p)] = min(abs(ExcitationSpectrum(1,:)-ALamb{p}(q2,1)));
            PhaseVelocityALamb{p} = ALamb{p}(q1:q2,4)*1e3;
            AttenuationALamb{p} = ALamb{p}(q1:q2,7).*PhaseVelocityALamb{p}./(ALamb{p}(q1:q2,1)*1e3);
            EnergyVelocityALamb(p) = ALamb{p}(ALamb{p}(:,1) == Frequency,5)*1e3;
        end
    end
    for p = 1:length(AShearModes)
        if  AShearModes(p)
            [~,q1] = min(abs(ExcitationSpectrum(1,1)-AShear{p}(:,1)));
            z = find(isapprox(diff(AShear{p}(:,1)),AShear{p}(end,1)-AShear{p}(end-1,1),'loose'),1,'first');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            [~,z1AShear(p)] = min(abs(ExcitationSpectrum(1,:)-AShear{p}(q1,1)));
            [~,q2] = min(abs(ExcitationSpectrum(1,end)-AShear{p}(:,1)));
            [~,z2AShear(p)] = min(abs(ExcitationSpectrum(1,:)-AShear{p}(q2,1)));
            PhaseVelocityAShear{p} = AShear{p}(q1:q2,4)*1e3;
            AttenuationAShear{p} = AShear{p}(q1:q2,7).*PhaseVelocityAShear{p}./(AShear{p}(q1:q2,1)*1e3);
            EnergyVelocityAShear(p) = AShear{p}(AShear{p}(:,1) == Frequency,5)*1e3;
        end
    end
    for p = 1:length(BLambModes)
        if  BLambModes(p)
            [~,q1] = min(abs(ExcitationSpectrum(1,1)-BLamb{p}(:,1)));
            z = find(isapprox(diff(BLamb{p}(:,1)),BLamb{p}(end,1)-BLamb{p}(end-1,1),'loose'),1,'first');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            [~,z1BLamb(p)] = min(abs(ExcitationSpectrum(1,:)-BLamb{p}(q1,1))); 
            [~,q2] = min(abs(ExcitationSpectrum(1,end)-BLamb{p}(:,1)));
            [~,z2BLamb(p)] = min(abs(ExcitationSpectrum(1,:)-BLamb{p}(q2,1)));
            PhaseVelocityBLamb{p} = BLamb{p}(q1:q2,4)*1e3;
            AttenuationBLamb{p} = BLamb{p}(q1:q2,7).*PhaseVelocityBLamb{p}./(BLamb{p}(q1:q2,1)*1e3);
            EnergyVelocityBLamb(p) = BLamb{p}(BLamb{p}(:,1) == Frequency,5)*1e3;
        end
    end
    for p = 1:length(BShearModes)
        if  BShearModes(p)
            [~,q1] = min(abs(ExcitationSpectrum(1,1)-BShear{p}(:,1)));
            z = find(isapprox(diff(BShear{p}(:,1)),BShear{p}(end,1)-BShear{p}(end-1,1),'loose'),1,'first');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            [~,z1BShear(p)] = min(abs(ExcitationSpectrum(1,:)-BShear{p}(q1,1)));
            [~,q2] = min(abs(ExcitationSpectrum(1,end)-BShear{p}(:,1)));
            [~,z2BShear(p)] = min(abs(ExcitationSpectrum(1,:)-BShear{p}(q2,1)));
            PhaseVelocityBShear{p} = BShear{p}(q1:q2,4)*1e3;
            AttenuationBShear{p} = BShear{p}(q1:q2,7).*PhaseVelocityBShear{p}./(BShear{p}(q1:q2,1)*1e3);
            EnergyVelocityBShear(p) = BShear{p}(BShear{p}(:,1) == Frequency,5)*1e3;
        end
    end
    for p = 1:length(SLambModes)
        if  SLambModes(p)
            [~,q1] = min(abs(ExcitationSpectrum(1,1)-SLamb{p}(:,1)));
            z = find(isapprox(diff(SLamb{p}(:,1)),SLamb{p}(end,1)-SLamb{p}(end-1,1),'loose'),1,'first');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            [~,z1SLamb(p)] = min(abs(ExcitationSpectrum(1,:)-SLamb{p}(q1,1))); 
            [~,q2] = min(abs(ExcitationSpectrum(1,end)-SLamb{p}(:,1)));
            [~,z2SLamb(p)] = min(abs(ExcitationSpectrum(1,:)-SLamb{p}(q2,1)));
            PhaseVelocitySLamb{p} = SLamb{p}(q1:q2,4)*1e3;
            AttenuationSLamb{p} = SLamb{p}(q1:q2,7).*PhaseVelocitySLamb{p}./(SLamb{p}(q1:q2,1)*1e3);
            EnergyVelocitySLamb(p) = SLamb{p}(SLamb{p}(:,1) == Frequency,5)*1e3;
        end
    end
    for p = 1:length(SShearModes)
        if  SShearModes(p)
            [~,q1] = min(abs(ExcitationSpectrum(1,1)-SShear{p}(:,1)));
            z = find(isapprox(diff(SShear{p}(:,1)),SShear{p}(end,1)-SShear{p}(end-1,1),'loose'),1,'first');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            [~,z1SShear(p)] = min(abs(ExcitationSpectrum(1,:)-SShear{p}(q1,1)));
            [~,q2] = min(abs(ExcitationSpectrum(1,end)-SShear{p}(:,1)));
            [~,z2SShear(p)] = min(abs(ExcitationSpectrum(1,:)-SShear{p}(q2,1)));
            PhaseVelocitySShear{p} = SShear{p}(q1:q2,4)*1e3;
            AttenuationSShear{p} = SShear{p}(q1:q2,7).*PhaseVelocitySShear{p}./(SShear{p}(q1:q2,1)*1e3);
            EnergyVelocitySShear(p) = SShear{p}(SShear{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity = [min(EnergyVelocityALamb) max(EnergyVelocityALamb);min(EnergyVelocityAShear) max(EnergyVelocityAShear);min(EnergyVelocityBLamb) max(EnergyVelocityBLamb);min(EnergyVelocityBShear) max(EnergyVelocityBShear);min(EnergyVelocitySLamb) max(EnergyVelocitySLamb);min(EnergyVelocitySShear) max(EnergyVelocitySShear)];    
    EnergyVelocity = [min(EnergyVelocity(:,1)) max(EnergyVelocity(:,2))];
    Gate = [Distance/1e3/EnergyVelocity(2)-.1*CoherenceTime Distance/1e3/EnergyVelocity(1)+.2*CoherenceTime]*1e6;
    TimeLimit = round(TimeLimitFactor*Distance/1e3/EnergyVelocity(1),6);
    if  TimeLimit-ExcitationSignal(1,SamplesPerCycle*Cycles/2+1) < ExcitationSignal(1,SamplesPerCycle*Cycles/2+1)
        Time = -ExcitationSignal(1,SamplesPerCycle*Cycles/2+1):1/SampleRate:ExcitationSignal(1,SamplesPerCycle*Cycles/2+1);
        PlotXAxis = ExcitationSignal(1,:)*1e6;
    else
        Time = -ExcitationSignal(1,SamplesPerCycle*Cycles/2+1):1/SampleRate:TimeLimit-ExcitationSignal(1,SamplesPerCycle*Cycles/2+1);
        PlotXAxis = (0:1/SampleRate:TimeLimit)*1e6;
    end
    if  CoherenceTime < TimeLimit
        for i = -1:100
            FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity(1)/CoherenceTime),i);
            CoherenceTime2 = 1e-3/FrequencyResolution2;
            if  FrequencyResolution2 > 0 && CoherenceTime2 > TimeLimit
                break
            end
        end
    else
        FrequencyResolution2 = FrequencyResolution;
    end
    ModeTotal = length(EnergyVelocityALamb(~isnan(EnergyVelocityALamb)))+length(EnergyVelocityAShear(~isnan(EnergyVelocityAShear)))+length(EnergyVelocityBLamb(~isnan(EnergyVelocityBLamb)))+length(EnergyVelocityBShear(~isnan(EnergyVelocityBShear)))+length(EnergyVelocitySLamb(~isnan(EnergyVelocitySLamb)))+length(EnergyVelocitySShear(~isnan(EnergyVelocitySShear)));
    tic
    h1 = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating mode...'); % generate waitbar
    Counter = 0;
    if  any(cellfun(@any,PhaseVelocityALamb))
        [uSumALamb,Counter,ExcitationSpectrumRangeALamb] = Computer_Anisotropic_Signal_Core(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'',ModeTotal,PhaseVelocityALamb,AttenuationALamb,Time,TimeLimit,z1ALamb,z2ALamb,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            close(h1)
            return
        end
        for p = 1:length(ExcitationSpectrumRangeALamb)
            if  ~isempty(ExcitationSpectrumRangeALamb{p})
                ExcitationSpectrumRangeLimits(1,:) = [ExcitationSpectrumRangeALamb{p}(1,1) ExcitationSpectrumRangeALamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocityAShear))
        [uSumAShear,Counter,ExcitationSpectrumRangeAShear] = Computer_Anisotropic_Signal_Core(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Shear',ModeTotal,PhaseVelocityAShear,AttenuationAShear,Time,TimeLimit,z1AShear,z2AShear,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            close(h1)
            return
        end
        for p = 1:length(ExcitationSpectrumRangeAShear)
            if  ~isempty(ExcitationSpectrumRangeAShear{p})
                ExcitationSpectrumRangeLimits(2,:) = [ExcitationSpectrumRangeAShear{p}(1,1) ExcitationSpectrumRangeAShear{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocityBLamb))
        [uSumBLamb,Counter,ExcitationSpectrumRangeBLamb] = Computer_Anisotropic_Signal_Core(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'',ModeTotal,PhaseVelocityBLamb,AttenuationBLamb,Time,TimeLimit,z1BLamb,z2BLamb,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            close(h1)
            return
        end
        for p = 1:length(ExcitationSpectrumRangeBLamb)
            if  ~isempty(ExcitationSpectrumRangeBLamb{p})
                ExcitationSpectrumRangeLimits(3,:) = [ExcitationSpectrumRangeBLamb{p}(1,1) ExcitationSpectrumRangeBLamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocityBShear))
        [uSumBShear,Counter,ExcitationSpectrumRangeBShear] = Computer_Anisotropic_Signal_Core(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Shear',ModeTotal,PhaseVelocityBShear,AttenuationBShear,Time,TimeLimit,z1BShear,z2BShear,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            close(h1)
            return
        end
        for p = 1:length(ExcitationSpectrumRangeBShear)
            if  ~isempty(ExcitationSpectrumRangeBShear{p})
                ExcitationSpectrumRangeLimits(4,:) = [ExcitationSpectrumRangeBShear{p}(1,1) ExcitationSpectrumRangeBShear{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocitySLamb))
        [uSumSLamb,Counter,ExcitationSpectrumRangeSLamb] = Computer_Anisotropic_Signal_Core(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'',ModeTotal,PhaseVelocitySLamb,AttenuationSLamb,Time,TimeLimit,z1SLamb,z2SLamb,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            close(h1)
            return
        end
        for p = 1:length(ExcitationSpectrumRangeSLamb)
            if  ~isempty(ExcitationSpectrumRangeSLamb{p})
                ExcitationSpectrumRangeLimits(5,:) = [ExcitationSpectrumRangeSLamb{p}(1,1) ExcitationSpectrumRangeSLamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocitySShear))
        [uSumSShear,~,ExcitationSpectrumRangeSShear] = Computer_Anisotropic_Signal_Core(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Shear',ModeTotal,PhaseVelocitySShear,AttenuationSShear,Time,TimeLimit,z1SShear,z2SShear,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            close(h1)
            return
        end
        for p = 1:length(ExcitationSpectrumRangeSShear)
            if  ~isempty(ExcitationSpectrumRangeSShear{p})
                ExcitationSpectrumRangeLimits(6,:) = [ExcitationSpectrumRangeSShear{p}(1,1) ExcitationSpectrumRangeSShear{p}(1,end)];
                break
            end
        end
    end
    ExcitationSpectrumRangeLimits = [min(ExcitationSpectrumRangeLimits(:,1)) min(ExcitationSpectrumRangeLimits(:,2))];
    if  CoherenceTime < TimeLimit
        String = ['TIME DOMAIN:',newline,...
        'Time limit:           ',num2str(TimeLimit*1e6),' micsec',newline,...
        'Coh. time (original): ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline,...
        'Coh. time (interp.):  ',num2str(CoherenceTime2*1e6,'%.0f'),' micsec',newline,...
        'Sample rate:          ',num2str(SampleRate/1e3),' kHz',newline,...
        'Samples:              ',num2str(length(Time)),newline,newline,...
        'FREQUENCY DOMAIN:',newline,...
        'Spectral range:  ',num2str(ExcitationSpectrumRangeLimits(1)),' - ',num2str(ExcitationSpectrumRangeLimits(2)),' kHz',newline,...
        'Res. (original): ',num2str(FrequencyResolution),' kHz',newline,...
        'Res. (interp.):  ',num2str(FrequencyResolution2),' kHz',newline,...
        'Frequencies:     ',num2str(length(ExcitationSpectrumRangeLimits(1):FrequencyResolution2:ExcitationSpectrumRangeLimits(2)))];
    else
        String = ['TIME DOMAIN:',newline,...
        'Time limit:     ',num2str(TimeLimit*1e6),' micsec',newline,...
        'Coherence time: ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline,...
        'Sample rate:    ',num2str(SampleRate/1e3,'%.0f'),' kHz',newline,...
        'Samples:        ',num2str(length(Time)),newline,newline,...
        'FREQUENCY DOMAIN:',newline,...
        'Spectral range: ',num2str(ExcitationSpectrumRangeLimits(1)),' - ',num2str(ExcitationSpectrumRangeLimits(2)),' kHz',newline,...
        'Resolution:     ',num2str(FrequencyResolution),' kHz',newline,...
        'Frequencies:    ',num2str(length(ExcitationSpectrumRangeLimits(1):FrequencyResolution2:ExcitationSpectrumRangeLimits(2)))];
    end
    close(h1)
end
if  Gate(1) < 0
    Gate(1) = 0;
end
if  Gate(2) > PlotXAxis(end)
    Gate(2) = PlotXAxis(end);
end
XAxis = [0 TimeLimit*1e6];
OutputWindowUI3.String = String;
disp([String,newline,'----------------------------------'])