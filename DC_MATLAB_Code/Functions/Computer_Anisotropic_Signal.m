% =========================================================================
% Dispersion Calculator
% Copyright (C) 2018-2023 DLR
% Created by Armin Huber
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function [ExcitationMagnitude,ExcitationSignal,ExcitationSpectrum,FourierTransformLength,Frequency,FrequencyLimitLow,FrequencyLimitHigh,FrequencyRange,PhaseVelocityA,PhaseVelocityALamb,PhaseVelocityAShear,PhaseVelocityB,PhaseVelocityBLamb,PhaseVelocityBShear,PhaseVelocityS,PhaseVelocitySLamb,PhaseVelocitySShear,PlotXAxis,Gate,XAxis,uSum,uSumA,uSumALamb,uSumAShear,uSumB,uSumBLamb,uSumBShear,uSumS,uSumSLamb,uSumSShear] = Computer_Anisotropic_Signal(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,A,ALamb,AShear,B,BLamb,BShear,c,Material,Mode,S,SLamb,SShear,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Symmetric,Decoupled,SamplesPerCycle,Cycles,Distance,FrequencyResolution,Frequency,MultiMode,SpectrumThreshold,TimeLimitFactor,AModes,ALambModes,AShearModes,BModes,BLambModes,BShearModes,SModes,SLambModes,SShearModes,Window,OutputWindowUI3)
%#ok<*AGROW>
CarrierWave = imag(exp(1i*2*pi*(0:1/SamplesPerCycle:Cycles))); % generate the carrier wave
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
ExcitationSpectrum = FrequencyRange(z1:z2); % only this frequency range is considered for the superposition of the wave components
ExcitationSpectrum(2,:) = ExcitationMagnitude(z1:z2); % the corresponding spectral amplitudes
CoherenceTime = 1e-3/FrequencyResolution; % (s) the coherence time tc is the time after which we have a repetition of the temporal response, i.e., we will see a twin of the wave packet at each instant of time n*tc, n=1,2,..., away from the actual wave packet    
FrequencyLimitLow = 0;
FrequencyLimitHigh = 0;
PhaseVelocityA{1} = 0;
PhaseVelocityALamb{1} = 0;
PhaseVelocityAShear{1} = 0;
PhaseVelocityB{1} = 0;
PhaseVelocityBLamb{1} = 0;
PhaseVelocityBShear{1} = 0;
PhaseVelocityS{1} = 0;
PhaseVelocitySLamb{1} = 0;
PhaseVelocitySShear{1} = 0;
uSum{1} = [];
uSumA{1} = [];
uSumALamb{1} = [];
uSumAShear{1} = [];
uSumB{1} = [];
uSumBLamb{1} = [];
uSumBShear{1} = [];
uSumS{1} = [];
uSumSLamb{1} = [];
uSumSShear{1} = [];
if  ~MultiMode
    p = str2double(regexp(Mode,'\d*','match'))+1;
    if  Symmetric
        if  ~Decoupled
            if  Mode(1) == 'S'
                q1 = find(abs(ExcitationSpectrum(1,1)-S{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-S{p}(:,1)))); % get the index of the starting frequency of the calculated frequency spectrum
                z = find(ischange(S{p}(:,1),'linear'),1,'last'); % find the index of the first frequency from a phase velocity sweep where we have equidistant frequencies; this is to exclude the high phase velocity part of higher order modes obtained by frequency sweeps where we don't have equidistant frequencies; the energy velocity of this range is set to zero and can be identified by that way
                if  isempty(z)
                    z = 1;
                end
                if  q1 < z % if the starting frequency is in the non-equidistant frequency range of higher order modes, we take the first frequency from the phase velocity sweep
                    q1 = z;
                end
                FrequencyLimitLow = S{p}(z,1); % used to indicate in the spectral amplitude plot below which frequency no data are available for higher order modes; it is the first frequency from a phase velocity sweep where we have equidistant frequencies; we indicate the range below this frequency with red shading, similarly as the range above the top frequency in the dispersion diagram 
                FrequencyLimitHigh = S{p}(end,1); % used to indicate above which frequency no data are available; this can happen in attenuated cases
                z1 = find(abs(ExcitationSpectrum(1,:)-S{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-S{p}(q1,1)))); % find the index in the calculated frequency spectrum corresponding the starting frequency (q1) determined above; we will shorten the frequency spectrum for the superposition of the wave components below
                q2 = find(abs(ExcitationSpectrum(1,end)-S{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-S{p}(:,1)))); % get the index of the ending frequency of the calculated frequency spectrum
                z2 = find(abs(ExcitationSpectrum(1,:)-S{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-S{p}(q2,1)))); % find the index corresponding the ending frequency (q2) determined above
                PhaseVelocity{1} = S{p}(q1:q2,4)*1e3; % extract the phase velocities of the requested mode in the frequency range in question
                Attenuation{1} = S{p}(q1:q2,7).*PhaseVelocity{1}./(S{p}(q1:q2,1)*1e3); % get the attenuation in Np/wavelength
                EnergyVelocity = S{p}(S{p}(:,1) == Frequency,5)*1e3; % get the energy velocity of the wave packet
                ModeType = 'Coupled';
            elseif Mode(1) == 'A'
                q1 = find(abs(ExcitationSpectrum(1,1)-A{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-A{p}(:,1))));
                z = find(ischange(A{p}(:,1),'linear'),1,'last');
                if  isempty(z)
                    z = 1;
                end
                if  q1 < z
                    q1 = z;
                end
                FrequencyLimitLow = A{p}(z,1);
                FrequencyLimitHigh = A{p}(end,1);
                z1 = find(abs(ExcitationSpectrum(1,:)-A{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-A{p}(q1,1))));
                q2 = find(abs(ExcitationSpectrum(1,end)-A{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-A{p}(:,1))));
                z2 = find(abs(ExcitationSpectrum(1,:)-A{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-A{p}(q2,1))));
                PhaseVelocity{1} = A{p}(q1:q2,4)*1e3;
                Attenuation{1} = A{p}(q1:q2,7).*PhaseVelocity{1}./(A{p}(q1:q2,1)*1e3);
                EnergyVelocity = A{p}(A{p}(:,1) == Frequency,5)*1e3;
                ModeType = 'Coupled';
            end
        else
            if  ~contains(Mode,'SH')
                if  Mode(1) == 'S'
                    q1 = find(abs(ExcitationSpectrum(1,1)-SLamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-SLamb{p}(:,1))));
                    z = find(ischange(SLamb{p}(:,1),'linear'),1,'last');
                    if  isempty(z)
                        z = 1;
                    end
                    if  q1 < z
                        q1 = z;
                    end
                    FrequencyLimitLow = SLamb{p}(z,1);
                    FrequencyLimitHigh = SLamb{p}(end,1);
                    z1 = find(abs(ExcitationSpectrum(1,:)-SLamb{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-SLamb{p}(q1,1))));
                    q2 = find(abs(ExcitationSpectrum(1,end)-SLamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-SLamb{p}(:,1))));
                    z2 = find(abs(ExcitationSpectrum(1,:)-SLamb{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-SLamb{p}(q2,1))));
                    PhaseVelocity{1} = SLamb{p}(q1:q2,4)*1e3;
                    Attenuation{1} = SLamb{p}(q1:q2,6).*PhaseVelocity{1}./(SLamb{p}(q1:q2,1)*1e3);
                    EnergyVelocity = SLamb{p}(SLamb{p}(:,1) == Frequency,5)*1e3;
                    ModeType = 'Lamb';
                elseif Mode(1) == 'A'
                    q1 = find(abs(ExcitationSpectrum(1,1)-ALamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-ALamb{p}(:,1))));
                    z = find(ischange(ALamb{p}(:,1),'linear'),1,'last');
                    if  isempty(z)
                        z = 1;
                    end
                    if  q1 < z
                        q1 = z;
                    end
                    FrequencyLimitLow = ALamb{p}(z,1);
                    FrequencyLimitHigh = ALamb{p}(end,1);
                    z1 = find(abs(ExcitationSpectrum(1,:)-ALamb{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-ALamb{p}(q1,1))));
                    q2 = find(abs(ExcitationSpectrum(1,end)-ALamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-ALamb{p}(:,1))));
                    z2 = find(abs(ExcitationSpectrum(1,:)-ALamb{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-ALamb{p}(q2,1))));
                    PhaseVelocity{1} = ALamb{p}(q1:q2,4)*1e3;
                    Attenuation{1} = ALamb{p}(q1:q2,6).*PhaseVelocity{1}./(ALamb{p}(q1:q2,1)*1e3);
                    EnergyVelocity = ALamb{p}(ALamb{p}(:,1) == Frequency,5)*1e3;
                    ModeType = 'Lamb';
                end
            elseif contains(Mode,'SH')
                if  Mode(1) == 'S'
                    q1 = find(abs(ExcitationSpectrum(1,1)-SShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-SShear{p}(:,1))));
                    z = find(ischange(SShear{p}(:,1),'linear'),1,'last');
                    if  isempty(z)
                        z = 1;
                    end
                    if  q1 < z
                        q1 = z;
                    end
                    FrequencyLimitLow = SShear{p}(z,1);
                    FrequencyLimitHigh = SShear{p}(end,1);
                    z1 = find(abs(ExcitationSpectrum(1,:)-SShear{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-SShear{p}(q1,1))));
                    q2 = find(abs(ExcitationSpectrum(1,end)-SShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-SShear{p}(:,1))));
                    z2 = find(abs(ExcitationSpectrum(1,:)-SShear{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-SShear{p}(q2,1))));
                    PhaseVelocity{1} = SShear{p}(q1:q2,4)*1e3;
                    Attenuation{1} = SShear{p}(q1:q2,6).*PhaseVelocity{1}./(SShear{p}(q1:q2,1)*1e3);
                    EnergyVelocity = SShear{p}(SShear{p}(:,1) == Frequency,5)*1e3;
                    ModeType = 'Shear';
                elseif Mode(1) == 'A'
                    q1 = find(abs(ExcitationSpectrum(1,1)-AShear{p-1}(:,1)) == min(abs(ExcitationSpectrum(1,1)-AShear{p-1}(:,1))));
                    z = find(ischange(AShear{p-1}(:,1),'linear'),1,'last');
                    if  isempty(z)
                        z = 1;
                    end
                    if  q1 < z
                        q1 = z;
                    end
                    FrequencyLimitLow = AShear{p-1}(z,1);
                    FrequencyLimitHigh = AShear{p-1}(end,1);
                    z1 = find(abs(ExcitationSpectrum(1,:)-AShear{p-1}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-AShear{p-1}(q1,1))));
                    q2 = find(abs(ExcitationSpectrum(1,end)-AShear{p-1}(:,1)) == min(abs(ExcitationSpectrum(1,end)-AShear{p-1}(:,1))));
                    z2 = find(abs(ExcitationSpectrum(1,:)-AShear{p-1}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-AShear{p-1}(q2,1))));
                    PhaseVelocity{1} = AShear{p-1}(q1:q2,4)*1e3;
                    Attenuation{1} = AShear{p-1}(q1:q2,6).*PhaseVelocity{1}./(AShear{p-1}(q1:q2,1)*1e3);
                    EnergyVelocity = AShear{p-1}(AShear{p-1}(:,1) == Frequency,5)*1e3;
                    ModeType = 'Shear';
                end
            end
        end
    else
        if  ~Decoupled
            q1 = find(abs(ExcitationSpectrum(1,1)-B{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-B{p}(:,1))));
            z = find(ischange(B{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            FrequencyLimitLow = B{p}(z,1);
            FrequencyLimitHigh = B{p}(end,1);
            z1 = find(abs(ExcitationSpectrum(1,:)-B{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-B{p}(q1,1))));
            q2 = find(abs(ExcitationSpectrum(1,end)-B{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-B{p}(:,1))));
            z2 = find(abs(ExcitationSpectrum(1,:)-B{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-B{p}(q2,1))));
            PhaseVelocity{1} = B{p}(q1:q2,4)*1e3;
            Attenuation{1} = B{p}(q1:q2,7).*PhaseVelocity{1}./(B{p}(q1:q2,1)*1e3);
            EnergyVelocity = B{p}(B{p}(:,1) == Frequency,5)*1e3;
            ModeType = 'Coupled';            
        else
            if  ~contains(Mode,'SH')
                q1 = find(abs(ExcitationSpectrum(1,1)-BLamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-BLamb{p}(:,1))));
                z = find(ischange(BLamb{p}(:,1),'linear'),1,'last');
                if  isempty(z)
                    z = 1;
                end
                if  q1 < z
                    q1 = z;
                end
                FrequencyLimitLow = BLamb{p}(z,1);
                FrequencyLimitHigh = BLamb{p}(end,1);
                z1 = find(abs(ExcitationSpectrum(1,:)-BLamb{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-BLamb{p}(q1,1))));
                q2 = find(abs(ExcitationSpectrum(1,end)-BLamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-BLamb{p}(:,1))));
                z2 = find(abs(ExcitationSpectrum(1,:)-BLamb{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-BLamb{p}(q2,1))));
                PhaseVelocity{1} = BLamb{p}(q1:q2,4)*1e3;
                Attenuation{1} = BLamb{p}(q1:q2,6).*PhaseVelocity{1}./(BLamb{p}(q1:q2,1)*1e3);
                EnergyVelocity = BLamb{p}(BLamb{p}(:,1) == Frequency,5)*1e3;
                ModeType = 'Lamb';
            elseif contains(Mode,'SH')
                q1 = find(abs(ExcitationSpectrum(1,1)-BShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-BShear{p}(:,1))));
                z = find(ischange(BShear{p}(:,1),'linear'),1,'last');
                if  isempty(z)
                    z = 1;
                end
                if  q1 < z
                    q1 = z;
                end
                FrequencyLimitLow = BShear{p}(z,1);
                FrequencyLimitHigh = BShear{p}(end,1);
                z1 = find(abs(ExcitationSpectrum(1,:)-BShear{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-BShear{p}(q1,1))));
                q2 = find(abs(ExcitationSpectrum(1,end)-BShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-BShear{p}(:,1))));
                z2 = find(abs(ExcitationSpectrum(1,:)-BShear{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-BShear{p}(q2,1))));
                PhaseVelocity{1} = BShear{p}(q1:q2,4)*1e3;
                Attenuation{1} = BShear{p}(q1:q2,6).*PhaseVelocity{1}./(BShear{p}(q1:q2,1)*1e3);
                EnergyVelocity = BShear{p}(BShear{p}(:,1) == Frequency,5)*1e3;
                ModeType = 'Shear';
            end
        end
    end
    Gate = [Distance/1e3/EnergyVelocity-.1*CoherenceTime Distance/1e3/EnergyVelocity+.2*CoherenceTime]*1e6;
    TimeLimit = round(TimeLimitFactor*Distance/1e3/EnergyVelocity,6); % (s) the time limit for which the temporal response is calculated
    if  TimeLimit-ExcitationSignal(1,SamplesPerCycle*Cycles/2+1) < ExcitationSignal(1,SamplesPerCycle*Cycles/2+1)
        Time = -ExcitationSignal(1,SamplesPerCycle*Cycles/2+1):1/SampleRate:ExcitationSignal(1,SamplesPerCycle*Cycles/2+1);
        PlotXAxis = ExcitationSignal(1,:)*1e6; % (microsec)
    else
        Time = -ExcitationSignal(1,SamplesPerCycle*Cycles/2+1):1/SampleRate:TimeLimit-ExcitationSignal(1,SamplesPerCycle*Cycles/2+1);
        PlotXAxis = (0:1/SampleRate:TimeLimit)*1e6; % (microsec)
    end
    if  CoherenceTime < TimeLimit % if the coherence time is smaller than the calculated time range, we have to expect seeing unwanted twins of the wave packet
        if  FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity/CoherenceTime) >= 1 % calculate the frequency resolution necessary to get a coherence time larger than the displayed time range
            FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity/CoherenceTime)); % (kHz)
        elseif FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity/CoherenceTime) < 1 && FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity/CoherenceTime) >= .1
            FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity/CoherenceTime),1);
        else
            FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity/CoherenceTime),2);
        end
        CoherenceTime2 = 1e-3/FrequencyResolution2; % the new coherence time (s)
    else
        FrequencyResolution2 = FrequencyResolution;
    end
    [uSum,~,ExcitationSpectrumRange] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,0,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,0,ModeType,0,PhaseVelocity,Attenuation,Time,TimeLimit,z1,z2);
    if  isempty(uSum{1})
        return
    end
    if  CoherenceTime < TimeLimit
        String = ['TIME DOMAIN:',newline,...
        'Time limit:           ',num2str(TimeLimit*1e6),' micsec',newline,...
        'Coh. time (original): ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline,...
        'Coh. time (interp.):  ',num2str(CoherenceTime2*1e6,'%.0f'),' micsec',newline,...
        'Sample rate:          ',num2str(SampleRate/1e6),' MHz',newline,...
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
        'Sample rate:    ',num2str(SampleRate/1e6,'%.0f'),' MHz',newline,...
        'Samples:        ',num2str(length(Time)),newline,newline,...
        'FREQUENCY DOMAIN:',newline,...
        'Spectral range: ',num2str(ExcitationSpectrumRange{1}(1,1)),' - ',num2str(ExcitationSpectrumRange{1}(1,end)),' kHz',newline,...
        'Resolution:     ',num2str(FrequencyResolution),' kHz',newline,...
        'Frequencies:    ',num2str(length(ExcitationSpectrumRange{1}))];
    end
else
    EnergyVelocityA = NaN(1,10);
    EnergyVelocityALamb = NaN(1,10);
    EnergyVelocityAShear = NaN(1,10);
    EnergyVelocityB = NaN(1,10);
    EnergyVelocityBLamb = NaN(1,10);
    EnergyVelocityBShear = NaN(1,10);
    EnergyVelocityS = NaN(1,10);
    EnergyVelocitySLamb = NaN(1,10);
    EnergyVelocitySShear = NaN(1,10);
    ExcitationSpectrumRangeLimits(1:9,1:2) = NaN;
    for p = 1:length(AModes)
        if  AModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-A{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-A{p}(:,1))));
            z = find(ischange(A{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1A(p) = find(abs(ExcitationSpectrum(1,:)-A{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-A{p}(q1,1)))); 
            q2 = find(abs(ExcitationSpectrum(1,end)-A{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-A{p}(:,1))));
            z2A(p) = find(abs(ExcitationSpectrum(1,:)-A{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-A{p}(q2,1))));
            PhaseVelocityA{p} = A{p}(q1:q2,4)*1e3;
            AttenuationA{p} = A{p}(q1:q2,7).*PhaseVelocityA{p}./(A{p}(q1:q2,1)*1e3);
            EnergyVelocityA(p) = A{p}(A{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity = min(EnergyVelocityA);
    EnergyVelocity(2) = max(EnergyVelocityA);
    for p = 1:length(ALambModes)
        if  ALambModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-ALamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-ALamb{p}(:,1))));
            z = find(ischange(ALamb{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1ALamb(p) = find(abs(ExcitationSpectrum(1,:)-ALamb{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-ALamb{p}(q1,1)))); 
            q2 = find(abs(ExcitationSpectrum(1,end)-ALamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-ALamb{p}(:,1))));
            z2ALamb(p) = find(abs(ExcitationSpectrum(1,:)-ALamb{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-ALamb{p}(q2,1))));
            PhaseVelocityALamb{p} = ALamb{p}(q1:q2,4)*1e3;
            AttenuationALamb{p} = ALamb{p}(q1:q2,6).*PhaseVelocityALamb{p}./(ALamb{p}(q1:q2,1)*1e3);
            EnergyVelocityALamb(p) = ALamb{p}(ALamb{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity(2,1) = min(EnergyVelocityALamb);
    EnergyVelocity(2,2) = max(EnergyVelocityALamb);    
    for p = 1:length(AShearModes)
        if  AShearModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-AShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-AShear{p}(:,1))));
            z = find(ischange(AShear{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1AShear(p) = find(abs(ExcitationSpectrum(1,:)-AShear{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-AShear{p}(q1,1))));
            q2 = find(abs(ExcitationSpectrum(1,end)-AShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-AShear{p}(:,1))));
            z2AShear(p) = find(abs(ExcitationSpectrum(1,:)-AShear{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-AShear{p}(q2,1))));
            PhaseVelocityAShear{p} = AShear{p}(q1:q2,4)*1e3;
            AttenuationAShear{p} = AShear{p}(q1:q2,6).*PhaseVelocityAShear{p}./(AShear{p}(q1:q2,1)*1e3);
            EnergyVelocityAShear(p) = AShear{p}(AShear{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity(3,1) = min(EnergyVelocityAShear);
    EnergyVelocity(3,2) = max(EnergyVelocityAShear);    
    for p = 1:length(BModes)
        if  BModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-B{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-B{p}(:,1))));
            z = find(ischange(B{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1B(p) = find(abs(ExcitationSpectrum(1,:)-B{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-B{p}(q1,1)))); 
            q2 = find(abs(ExcitationSpectrum(1,end)-B{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-B{p}(:,1))));
            z2B(p) = find(abs(ExcitationSpectrum(1,:)-B{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-B{p}(q2,1))));
            PhaseVelocityB{p} = B{p}(q1:q2,4)*1e3;
            AttenuationB{p} = B{p}(q1:q2,7).*PhaseVelocityB{p}./(B{p}(q1:q2,1)*1e3);
            EnergyVelocityB(p) = B{p}(B{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity(4,1) = min(EnergyVelocityB);
    EnergyVelocity(4,2) = max(EnergyVelocityB);
    for p = 1:length(BLambModes)
        if  BLambModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-BLamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-BLamb{p}(:,1))));
            z = find(ischange(BLamb{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1BLamb(p) = find(abs(ExcitationSpectrum(1,:)-BLamb{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-BLamb{p}(q1,1)))); 
            q2 = find(abs(ExcitationSpectrum(1,end)-BLamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-BLamb{p}(:,1))));
            z2BLamb(p) = find(abs(ExcitationSpectrum(1,:)-BLamb{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-BLamb{p}(q2,1))));
            PhaseVelocityBLamb{p} = BLamb{p}(q1:q2,4)*1e3;
            AttenuationBLamb{p} = BLamb{p}(q1:q2,6).*PhaseVelocityBLamb{p}./(BLamb{p}(q1:q2,1)*1e3);
            EnergyVelocityBLamb(p) = BLamb{p}(BLamb{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity(5,1) = min(EnergyVelocityBLamb);
    EnergyVelocity(5,2) = max(EnergyVelocityBLamb);    
    for p = 1:length(BShearModes)
        if  BShearModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-BShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-BShear{p}(:,1))));
            z = find(ischange(BShear{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1BShear(p) = find(abs(ExcitationSpectrum(1,:)-BShear{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-BShear{p}(q1,1))));
            q2 = find(abs(ExcitationSpectrum(1,end)-BShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-BShear{p}(:,1))));
            z2BShear(p) = find(abs(ExcitationSpectrum(1,:)-BShear{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-BShear{p}(q2,1))));
            PhaseVelocityBShear{p} = BShear{p}(q1:q2,4)*1e3;
            AttenuationBShear{p} = BShear{p}(q1:q2,6).*PhaseVelocityBShear{p}./(BShear{p}(q1:q2,1)*1e3);
            EnergyVelocityBShear(p) = BShear{p}(BShear{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity(6,1) = min(EnergyVelocityBShear);
    EnergyVelocity(6,2) = max(EnergyVelocityBShear);    
    for p = 1:length(SModes)
        if  SModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-S{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-S{p}(:,1))));
            z = find(ischange(S{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1S(p) = find(abs(ExcitationSpectrum(1,:)-S{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-S{p}(q1,1)))); 
            q2 = find(abs(ExcitationSpectrum(1,end)-S{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-S{p}(:,1))));
            z2S(p) = find(abs(ExcitationSpectrum(1,:)-S{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-S{p}(q2,1))));
            PhaseVelocityS{p} = S{p}(q1:q2,4)*1e3;
            AttenuationS{p} = S{p}(q1:q2,7).*PhaseVelocityS{p}./(S{p}(q1:q2,1)*1e3);
            EnergyVelocityS(p) = S{p}(S{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity(7,1) = min(EnergyVelocityS);
    EnergyVelocity(7,2) = max(EnergyVelocityS);    
    for p = 1:length(SLambModes)
        if  SLambModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-SLamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-SLamb{p}(:,1))));
            z = find(ischange(SLamb{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1SLamb(p) = find(abs(ExcitationSpectrum(1,:)-SLamb{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-SLamb{p}(q1,1)))); 
            q2 = find(abs(ExcitationSpectrum(1,end)-SLamb{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-SLamb{p}(:,1))));
            z2SLamb(p) = find(abs(ExcitationSpectrum(1,:)-SLamb{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-SLamb{p}(q2,1))));
            PhaseVelocitySLamb{p} = SLamb{p}(q1:q2,4)*1e3;
            AttenuationSLamb{p} = SLamb{p}(q1:q2,6).*PhaseVelocitySLamb{p}./(SLamb{p}(q1:q2,1)*1e3);
            EnergyVelocitySLamb(p) = SLamb{p}(SLamb{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity(8,1) = min(EnergyVelocitySLamb);
    EnergyVelocity(8,2) = max(EnergyVelocitySLamb);
    for p = 1:length(SShearModes)
        if  SShearModes(p)
            q1 = find(abs(ExcitationSpectrum(1,1)-SShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,1)-SShear{p}(:,1))));
            z = find(ischange(SShear{p}(:,1),'linear'),1,'last');
            if  isempty(z)
                z = 1;
            end
            if  q1 < z
                q1 = z;
            end
            z1SShear(p) = find(abs(ExcitationSpectrum(1,:)-SShear{p}(q1,1)) == min(abs(ExcitationSpectrum(1,:)-SShear{p}(q1,1))));
            q2 = find(abs(ExcitationSpectrum(1,end)-SShear{p}(:,1)) == min(abs(ExcitationSpectrum(1,end)-SShear{p}(:,1))));
            z2SShear(p) = find(abs(ExcitationSpectrum(1,:)-SShear{p}(q2,1)) == min(abs(ExcitationSpectrum(1,:)-SShear{p}(q2,1))));
            PhaseVelocitySShear{p} = SShear{p}(q1:q2,4)*1e3;
            AttenuationSShear{p} = SShear{p}(q1:q2,6).*PhaseVelocitySShear{p}./(SShear{p}(q1:q2,1)*1e3);
            EnergyVelocitySShear(p) = SShear{p}(SShear{p}(:,1) == Frequency,5)*1e3;
        end
    end
    EnergyVelocity(9,1) = min(EnergyVelocitySShear);
    EnergyVelocity(9,2) = max(EnergyVelocitySShear);
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
        if  FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity(1)/CoherenceTime) >= 1
            FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity(1)/CoherenceTime));
        elseif FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity(1)/CoherenceTime) < 1 && FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity(1)/CoherenceTime) >= .1
            FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity(1)/CoherenceTime),1);
        else
            FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity(1)/CoherenceTime),2);
        end
        CoherenceTime2 = 1e-3/FrequencyResolution2;
    else
        FrequencyResolution2 = FrequencyResolution;
    end
    ModeTotal = length(EnergyVelocityA(~isnan(EnergyVelocityA)))+length(EnergyVelocityALamb(~isnan(EnergyVelocityALamb)))+length(EnergyVelocityAShear(~isnan(EnergyVelocityAShear)))+length(EnergyVelocityB(~isnan(EnergyVelocityB)))+length(EnergyVelocityBLamb(~isnan(EnergyVelocityBLamb)))+length(EnergyVelocityBShear(~isnan(EnergyVelocityBShear)))+length(EnergyVelocityS(~isnan(EnergyVelocityS)))+length(EnergyVelocitySLamb(~isnan(EnergyVelocitySLamb)))+length(EnergyVelocitySShear(~isnan(EnergyVelocitySShear)));
    tic
    h1 = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating mode...'); % generate waitbar
    Counter = 0;
    if  any(cellfun(@any,PhaseVelocityA))
        [uSumA,Counter,ExcitationSpectrumRangeA] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Coupled',ModeTotal,PhaseVelocityA,AttenuationA,Time,TimeLimit,z1A,z2A);
        for p = 1:length(ExcitationSpectrumRangeA)
            if  ~isempty(ExcitationSpectrumRangeA{p})
                ExcitationSpectrumRangeLimits(1,1:2) = [ExcitationSpectrumRangeA{p}(1,1) ExcitationSpectrumRangeA{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocityALamb))
        [uSumALamb,Counter,ExcitationSpectrumRangeALamb] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Lamb',ModeTotal,PhaseVelocityALamb,AttenuationALamb,Time,TimeLimit,z1ALamb,z2ALamb);
        for p = 1:length(ExcitationSpectrumRangeALamb)
            if  ~isempty(ExcitationSpectrumRangeALamb{p})
                ExcitationSpectrumRangeLimits(2,1:2) = [ExcitationSpectrumRangeALamb{p}(1,1) ExcitationSpectrumRangeALamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocityAShear))
        [uSumAShear,Counter,ExcitationSpectrumRangeAShear] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Shear',ModeTotal,PhaseVelocityAShear,AttenuationAShear,Time,TimeLimit,z1AShear,z2AShear);
        for p = 1:length(ExcitationSpectrumRangeAShear)
            if  ~isempty(ExcitationSpectrumRangeAShear{p})
                ExcitationSpectrumRangeLimits(3,1:2) = [ExcitationSpectrumRangeAShear{p}(1,1) ExcitationSpectrumRangeAShear{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocityB))
        [uSumB,Counter,ExcitationSpectrumRangeB] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Coupled',ModeTotal,PhaseVelocityB,AttenuationB,Time,TimeLimit,z1B,z2B);
        for p = 1:length(ExcitationSpectrumRangeB)
            if  ~isempty(ExcitationSpectrumRangeB{p})
                ExcitationSpectrumRangeLimits(4,1:2) = [ExcitationSpectrumRangeB{p}(1,1) ExcitationSpectrumRangeB{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocityBLamb))
        [uSumBLamb,Counter,ExcitationSpectrumRangeBLamb] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Lamb',ModeTotal,PhaseVelocityBLamb,AttenuationBLamb,Time,TimeLimit,z1BLamb,z2BLamb);
        for p = 1:length(ExcitationSpectrumRangeBLamb)
            if  ~isempty(ExcitationSpectrumRangeBLamb{p})
                ExcitationSpectrumRangeLimits(5,1:2) = [ExcitationSpectrumRangeBLamb{p}(1,1) ExcitationSpectrumRangeBLamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocityBShear))
        [uSumBShear,Counter,ExcitationSpectrumRangeBShear] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Shear',ModeTotal,PhaseVelocityBShear,AttenuationBShear,Time,TimeLimit,z1BShear,z2BShear);
        for p = 1:length(ExcitationSpectrumRangeBShear)
            if  ~isempty(ExcitationSpectrumRangeBShear{p})
                ExcitationSpectrumRangeLimits(6,1:2) = [ExcitationSpectrumRangeBShear{p}(1,1) ExcitationSpectrumRangeBShear{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocityS))
        [uSumS,Counter,ExcitationSpectrumRangeS] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Coupled',ModeTotal,PhaseVelocityS,AttenuationS,Time,TimeLimit,z1S,z2S);
        for p = 1:length(ExcitationSpectrumRangeS)
            if  ~isempty(ExcitationSpectrumRangeS{p})
                ExcitationSpectrumRangeLimits(7,1:2) = [ExcitationSpectrumRangeS{p}(1,1) ExcitationSpectrumRangeS{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocitySLamb))
        [uSumSLamb,Counter,ExcitationSpectrumRangeSLamb] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Lamb',ModeTotal,PhaseVelocitySLamb,AttenuationSLamb,Time,TimeLimit,z1SLamb,z2SLamb);
        for p = 1:length(ExcitationSpectrumRangeSLamb)
            if  ~isempty(ExcitationSpectrumRangeSLamb{p})
                ExcitationSpectrumRangeLimits(8,1:2) = [ExcitationSpectrumRangeSLamb{p}(1,1) ExcitationSpectrumRangeSLamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocitySShear))
        [uSumSShear,~,ExcitationSpectrumRangeSShear] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance/1e3,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,'Shear',ModeTotal,PhaseVelocitySShear,AttenuationSShear,Time,TimeLimit,z1SShear,z2SShear);
        for p = 1:length(ExcitationSpectrumRangeSShear)
            if  ~isempty(ExcitationSpectrumRangeSShear{p})
                ExcitationSpectrumRangeLimits(9,1:2) = [ExcitationSpectrumRangeSShear{p}(1,1) ExcitationSpectrumRangeSShear{p}(1,end)];
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
        'Sample rate:          ',num2str(SampleRate/1e6),' MHz',newline,...
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
        'Sample rate:    ',num2str(SampleRate/1e6,'%.0f'),' MHz',newline,...
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