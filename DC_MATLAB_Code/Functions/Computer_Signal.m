% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function [ExcitationMagnitude,ExcitationSignal,ExcitationSpectrum,FourierTransformLength,FrequencyLimitLow,FrequencyLimitHigh,FrequencyRange,PlotXAxis,Gate,XAxis,uSum,uSumALamb,uSumAShear,uSumBLamb,uSumBShear,uSumSLamb,uSumSShear] = Computer_Signal(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,ALamb,AShear,BLamb,BShear,c,Material,Mode,SLamb,SShear,Layers,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,Delta,SamplesPerCycle,Cycles,Distance,FrequencyResolution,Chirp,ChirpType,Frequency,MultiMode,SpectrumThreshold,TimeLimitFactor,ALambModes,AShearModes,BLambModes,BShearModes,SLambModes,SShearModes,Window,WidthFactor,OutputWindowUI3)
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
    if  Decoupled || (~Decoupled && DC(m))
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
        a11(m)=0;a12(m)=0;a23(m)=0;a33(m)=0;a34(m)=0;
    else
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
        A1(m)=0;
    end
end
XAxis = [];
if  ~MultiMode
    p = str2double(regexp(Mode,'\d*','match'))+1;
    if  ~contains(Mode,'SH') && Mode(1) == 'S'
        [X,EnergyVelocity,Direction] = EstimateEnergyVelocity(SLamb{p},Frequency);
    elseif ~contains(Mode,'SH') && Mode(1) == 'A'
        [X,EnergyVelocity,Direction] = EstimateEnergyVelocity(ALamb{p},Frequency);
    elseif ~contains(Mode,'SH') && Mode(1) == 'B'
        [X,EnergyVelocity,Direction] = EstimateEnergyVelocity(BLamb{p},Frequency);
    elseif contains(Mode,'SH') && Mode(1) == 'S'
        [X,EnergyVelocity,Direction] = EstimateEnergyVelocity(SShear{p},Frequency);
    elseif contains(Mode,'SH') && Mode(1) == 'A'
        [X,EnergyVelocity,Direction] = EstimateEnergyVelocity(AShear{p-1},Frequency);
    elseif contains(Mode,'SH') && Mode(1) == 'B'
        [X,EnergyVelocity,Direction] = EstimateEnergyVelocity(BShear{p},Frequency);
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
            [~,EnergyVelocityALamb(p),DirectionALamb(p)] = EstimateEnergyVelocity(ALamb{p},Frequency);
        end
    end
    for p = 1:length(AShearModes)
        if  AShearModes(p)
            [~,EnergyVelocityAShear(p),DirectionAShear(p)] = EstimateEnergyVelocity(AShear{p},Frequency);
        end
    end
    for p = 1:length(BLambModes)
        if  BLambModes(p)
            [~,EnergyVelocityBLamb(p),DirectionBLamb(p)] = EstimateEnergyVelocity(BLamb{p},Frequency);
        end
    end
    for p = 1:length(BShearModes)
        if  BShearModes(p)
            [~,EnergyVelocityBShear(p),DirectionBShear(p)] = EstimateEnergyVelocity(BShear{p},Frequency);
        end
    end
    for p = 1:length(SLambModes)
        if  SLambModes(p)
            [~,EnergyVelocitySLamb(p),DirectionSLamb(p)] = EstimateEnergyVelocity(SLamb{p},Frequency);
        end
    end
    for p = 1:length(SShearModes)
        if  SShearModes(p)
            [~,EnergyVelocitySShear(p),DirectionSShear(p)] = EstimateEnergyVelocity(SShear{p},Frequency);
        end
    end
    EnergyVelocity = [min(EnergyVelocityALamb) max(EnergyVelocityALamb);min(EnergyVelocityAShear) max(EnergyVelocityAShear);min(EnergyVelocityBLamb) max(EnergyVelocityBLamb);min(EnergyVelocityBShear) max(EnergyVelocityBShear);min(EnergyVelocitySLamb) max(EnergyVelocitySLamb);min(EnergyVelocitySShear) max(EnergyVelocitySShear)];    
    EnergyVelocity = [min(EnergyVelocity(:,1)) max(EnergyVelocity(:,2))];
end
CoherenceTime = 1e-3*min(Frequency)/max(Frequency)/FrequencyResolution; % (s) the coherence time tc is the time after which we have a repetition of the temporal response, i.e., we will see a twin of the wave packet at each instant of time n*tc, n=1,2,..., away from the actual wave packet
TimeLimit = round(TimeLimitFactor*Distance/1e3/EnergyVelocity(1),6); % (s) the time limit for which the temporal response is calculated
if  CoherenceTime < TimeLimit % if the coherence time is smaller than the calculated time range, we have to expect seeing unwanted twins of the wave packet
    for i = -1:100
        FrequencyResolution2 = round(FrequencyResolution/(TimeLimitFactor*Distance/1e3/EnergyVelocity(1)/CoherenceTime),i); % (kHz) calculate the frequency resolution necessary to get a coherence time larger than the displayed time range
        CoherenceTime2 = 1e-3/FrequencyResolution2; % the new coherence time (s)
        if  FrequencyResolution2 > 0 && CoherenceTime2 > TimeLimit
            break
        end
    end
else
    FrequencyResolution2 = FrequencyResolution;
end
[FrequencyRange,ExcitationMagnitude,ExcitationSignal,SampleRate,FourierTransformLength] = WaveGenerator(Frequency,FrequencyResolution2,Cycles,SamplesPerCycle,Chirp,ChirpType,Window,WidthFactor);
if  TimeLimit < ExcitationSignal(1,end)
    Time = ExcitationSignal(1,:);
else
    Time = 0:1/SampleRate:TimeLimit;
end
z1 = find(ExcitationMagnitude(2,:) > SpectrumThreshold/100*max(ExcitationMagnitude(2,:)),1); % find the index where the spectral amplitudes start to be higher than the SpectrumThreshold
z2 = find(ExcitationMagnitude(2,:) > SpectrumThreshold/100*max(ExcitationMagnitude(2,:)),1,'last'); % find the index where the spectral amplitudes end to be higher than the SpectrumThreshold
if  FrequencyRange(z1) == 0
    z1 = 2;
end
ExcitationSpectrum = ExcitationMagnitude(:,z1:z2); % only this frequency range is considered for the superposition of the wave components; the corresponding spectral amplitudes and phases
FrequencyLimitLow = 0;
FrequencyLimitHigh = 0;
PhaseVelocityALamb{1} = 0;
PhaseVelocitySLamb{1} = 0;
PhaseVelocityBLamb{1} = 0;
PhaseVelocityAShear{1} = 0;
PhaseVelocitySShear{1} = 0;
PhaseVelocityBShear{1} = 0;
uSum = {[]};
uSumALamb = {[]};
uSumSLamb = {[]};
uSumBLamb = {[]};
uSumAShear = {[]};
uSumSShear = {[]};
uSumBShear = {[]};
if  ~MultiMode
    [PhaseVelocity,Attenuation,SignFluids,ExcitationSpectrumRange,FrequencyLimitLow,FrequencyLimitHigh,ModeType] = ExtractData(X,ExcitationSpectrum,Mode);
    Gate = (Distance/1e3/EnergyVelocity+[-.1*CoherenceTime .2*CoherenceTime])*1e6;
    uSum = Computer(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,0,Distance/1e3,ExcitationSpectrumRange,0,ModeType,0,PhaseVelocity,Attenuation,SignFluids,Direction,Time,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
    if  isempty(uSum{1})
        return
    end
    if  CoherenceTime < TimeLimit
        String = ['TIME DOMAIN:',newline...
        'Time limit:           ',num2str(TimeLimit*1e6),' micsec',newline...
        'Coh. time (original): ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline...
        'Coh. time (interp.):  ',num2str(CoherenceTime2*1e6,'%.0f'),' micsec',newline...
        'Sample rate:          ',num2str(SampleRate/1e3),' kHz',newline...
        'Samples:              ',num2str(length(Time)),newline,newline...
        'FREQUENCY DOMAIN:',newline...
        'Range:  ',num2str(ExcitationSpectrumRange{1}(1)),' - ',num2str(ExcitationSpectrumRange{1}(1,end)),' kHz',newline...
        'Res. (original): ',num2str(FrequencyResolution),' kHz',newline...
        'Res. (interp.):  ',num2str(FrequencyResolution2),' kHz',newline...
        'Frequencies:     ',num2str(length(ExcitationSpectrumRange{1}))];
    else
        String = ['TIME DOMAIN:',newline...
        'Time limit:     ',num2str(TimeLimit*1e6),' micsec',newline...
        'Coherence time: ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline...
        'Sample rate:    ',num2str(SampleRate/1e3,'%.0f'),' kHz',newline...
        'Samples:        ',num2str(length(Time)),newline,newline...
        'FREQUENCY DOMAIN:',newline...
        'Range: ',num2str(ExcitationSpectrumRange{1}(1)),' - ',num2str(ExcitationSpectrumRange{1}(1,end)),' kHz',newline...
        'Resolution:     ',num2str(FrequencyResolution),' kHz',newline...
        'Frequencies:    ',num2str(length(ExcitationSpectrumRange{1}))];
    end
else
    for p = 1:length(ALambModes)
        if  ALambModes(p)
            [PhaseVelocityALamb(p),AttenuationALamb(p),SignFluidsALamb(p),ExcitationSpectrumRangeALamb(p),~,~,~] = ExtractData(ALamb{p},ExcitationSpectrum,'Lamb');
        end
    end
    for p = 1:length(AShearModes)
        if  AShearModes(p)
            [PhaseVelocityAShear(p),AttenuationAShear(p),~,ExcitationSpectrumRangeAShear(p),~,~,~] = ExtractData(AShear{p},ExcitationSpectrum,'SH');
        end
    end
    for p = 1:length(BLambModes)
        if  BLambModes(p)
            [PhaseVelocityBLamb(p),AttenuationBLamb(p),SignFluidsBLamb(p),ExcitationSpectrumRangeBLamb(p),~,~,~] = ExtractData(BLamb{p},ExcitationSpectrum,'Lamb');
        end
    end
    for p = 1:length(BShearModes)
        if  BShearModes(p)
            [PhaseVelocityBShear(p),AttenuationBShear(p),~,ExcitationSpectrumRangeBShear(p),~,~,~] = ExtractData(BShear{p},ExcitationSpectrum,'SH');
        end
    end
    for p = 1:length(SLambModes)
        if  SLambModes(p)
            [PhaseVelocitySLamb(p),AttenuationSLamb(p),SignFluidsSLamb(p),ExcitationSpectrumRangeSLamb(p),~,~,~] = ExtractData(SLamb{p},ExcitationSpectrum,'Lamb');
        end
    end
    for p = 1:length(SShearModes)
        if  SShearModes(p)
            [PhaseVelocitySShear(p),AttenuationSShear(p),~,ExcitationSpectrumRangeSShear(p),~,~,~] = ExtractData(SShear{p},ExcitationSpectrum,'SH');
        end
    end
    Gate = [Distance/1e3/EnergyVelocity(2)-.1*CoherenceTime Distance/1e3/EnergyVelocity(1)+.2*CoherenceTime]*1e6;
    ModeTotal = length(EnergyVelocityALamb(~isnan(EnergyVelocityALamb)))+length(EnergyVelocityAShear(~isnan(EnergyVelocityAShear)))+length(EnergyVelocityBLamb(~isnan(EnergyVelocityBLamb)))+length(EnergyVelocityBShear(~isnan(EnergyVelocityBShear)))+length(EnergyVelocitySLamb(~isnan(EnergyVelocitySLamb)))+length(EnergyVelocitySShear(~isnan(EnergyVelocitySShear)));
    tic
    h1 = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating mode...'); % generate waitbar
    Counter = 0;
    if  any(cellfun(@any,PhaseVelocityALamb))
        [uSumALamb,Counter] = Computer(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,Counter,Distance/1e3,ExcitationSpectrumRangeALamb,h1,'Lamb',ModeTotal,PhaseVelocityALamb,AttenuationALamb,SignFluidsALamb,DirectionALamb,Time,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            return
        end
        for p = 1:length(ExcitationSpectrumRangeALamb)
            if  ~isempty(ExcitationSpectrumRangeALamb{p})
                ExcitationSpectrumRangeLimits(1,:) = [ExcitationSpectrumRangeALamb{p}(1) ExcitationSpectrumRangeALamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocityAShear))
        [uSumAShear,Counter] = Computer(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,Counter,Distance/1e3,ExcitationSpectrumRangeAShear,h1,'Shear',ModeTotal,PhaseVelocityAShear,AttenuationAShear,0,DirectionAShear,Time,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            return
        end
        for p = 1:length(ExcitationSpectrumRangeAShear)
            if  ~isempty(ExcitationSpectrumRangeAShear{p})
                ExcitationSpectrumRangeLimits(2,:) = [ExcitationSpectrumRangeAShear{p}(1) ExcitationSpectrumRangeAShear{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocityBLamb))
        [uSumBLamb,Counter] = Computer(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,Counter,Distance/1e3,ExcitationSpectrumRangeBLamb,h1,'Lamb',ModeTotal,PhaseVelocityBLamb,AttenuationBLamb,SignFluidsBLamb,DirectionBLamb,Time,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            return
        end
        for p = 1:length(ExcitationSpectrumRangeBLamb)
            if  ~isempty(ExcitationSpectrumRangeBLamb{p})
                ExcitationSpectrumRangeLimits(3,:) = [ExcitationSpectrumRangeBLamb{p}(1) ExcitationSpectrumRangeBLamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocityBShear))
        [uSumBShear,Counter] = Computer(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,Counter,Distance/1e3,ExcitationSpectrumRangeBShear,h1,'Shear',ModeTotal,PhaseVelocityBShear,AttenuationBShear,0,DirectionBShear,Time,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            return
        end
        for p = 1:length(ExcitationSpectrumRangeBShear)
            if  ~isempty(ExcitationSpectrumRangeBShear{p})
                ExcitationSpectrumRangeLimits(4,:) = [ExcitationSpectrumRangeBShear{p}(1) ExcitationSpectrumRangeBShear{p}(1,end)];
                break
            end
        end
    end
    if  any(cellfun(@any,PhaseVelocitySLamb))
        [uSumSLamb,Counter] = Computer(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,Counter,Distance/1e3,ExcitationSpectrumRangeSLamb,h1,'Lamb',ModeTotal,PhaseVelocitySLamb,AttenuationSLamb,SignFluidsSLamb,DirectionSLamb,Time,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            return
        end
        for p = 1:length(ExcitationSpectrumRangeSLamb)
            if  ~isempty(ExcitationSpectrumRangeSLamb{p})
                ExcitationSpectrumRangeLimits(5,:) = [ExcitationSpectrumRangeSLamb{p}(1) ExcitationSpectrumRangeSLamb{p}(1,end)];
                break
            end
        end
    end    
    if  any(cellfun(@any,PhaseVelocitySShear))
        [uSumSShear] = Computer(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,Counter,Distance/1e3,ExcitationSpectrumRangeSShear,h1,'Shear',ModeTotal,PhaseVelocitySShear,AttenuationSShear,0,DirectionSShear,Time,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory);
        if  Stop
            return
        end
        for p = 1:length(ExcitationSpectrumRangeSShear)
            if  ~isempty(ExcitationSpectrumRangeSShear{p})
                ExcitationSpectrumRangeLimits(6,:) = [ExcitationSpectrumRangeSShear{p}(1) ExcitationSpectrumRangeSShear{p}(1,end)];
                break
            end
        end
    end
    ExcitationSpectrumRangeLimits = [min(ExcitationSpectrumRangeLimits(:,1)) min(ExcitationSpectrumRangeLimits(:,2))];
    if  CoherenceTime < TimeLimit
        String = ['TIME DOMAIN:',newline...
        'Time limit:           ',num2str(TimeLimit*1e6),' micsec',newline...
        'Coh. time (original): ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline...
        'Coh. time (interp.):  ',num2str(CoherenceTime2*1e6,'%.0f'),' micsec',newline...
        'Sample rate:          ',num2str(SampleRate/1e3),' kHz',newline...
        'Samples:              ',num2str(length(Time)),newline,newline...
        'FREQUENCY DOMAIN:',newline...
        'Range:  ',num2str(ExcitationSpectrumRangeLimits(1)),' - ',num2str(ExcitationSpectrumRangeLimits(2)),' kHz',newline...
        'Res. (original): ',num2str(FrequencyResolution),' kHz',newline...
        'Res. (interp.):  ',num2str(FrequencyResolution2),' kHz',newline...
        'Frequencies:     ',num2str(length(ExcitationSpectrumRangeLimits(1):FrequencyResolution2:ExcitationSpectrumRangeLimits(2)))];
    else
        String = ['TIME DOMAIN:',newline...
        'Time limit:     ',num2str(TimeLimit*1e6),' micsec',newline...
        'Coherence time: ',num2str(CoherenceTime*1e6,'%.0f'),' micsec',newline...
        'Sample rate:    ',num2str(SampleRate/1e3,'%.0f'),' kHz',newline...
        'Samples:        ',num2str(length(Time)),newline,newline...
        'FREQUENCY DOMAIN:',newline...
        'Range: ',num2str(ExcitationSpectrumRangeLimits(1)),' - ',num2str(ExcitationSpectrumRangeLimits(2)),' kHz',newline...
        'Resolution:     ',num2str(FrequencyResolution),' kHz',newline...
        'Frequencies:    ',num2str(length(ExcitationSpectrumRangeLimits(1):FrequencyResolution2:ExcitationSpectrumRangeLimits(2)))];
    end
    close(h1)
end
XAxis = [0 TimeLimit*1e6];
PlotXAxis = Time*1e6; % (micro sec)
if  Gate(1) < 0
    Gate(1) = 0;
end
if  Gate(2) > PlotXAxis(end)
    Gate(2) = PlotXAxis(end);
end
OutputWindowUI3.String = String;
disp([String,newline,'----------------------------------'])
end
function [X,EnergyVelocity,Direction] = EstimateEnergyVelocity(X,Frequency)
    [~,q] = min(abs(X(:,1)-mean(Frequency)));
    EnergyVelocity = X(q,5)*1e3; % estimate the energy velocity of the wave packet (m/s)
    if  EnergyVelocity > 0
        Direction = 1;
    else
        Direction = -1; % reverse propagation direction for backward propagating modes
    end
    if  EnergyVelocity < 200 % avoid excessively long time ranges in case of non-propagating modes
        EnergyVelocity = 200;
    end
end
function [PhaseVelocity,Attenuation,SignFluids,ExcitationSpectrumRange,FrequencyLimitLow,FrequencyLimitHigh,ModeType] = ExtractData(X,ExcitationSpectrum,Mode)
    [FrequencyLimitLow,FrequencyLimitHigh] = bounds(X(:,1)); % used to indicate frequencies where no data are available
    [~,q1] = min(abs(ExcitationSpectrum(1)-X(:,1))); % get the index of the starting frequency of the calculated frequency spectrum
    [~,q2] = min(abs(ExcitationSpectrum(1,end)-X(:,1))); % get the index of the ending frequency of the calculated frequency spectrum
    [~,z1] = min(abs(ExcitationSpectrum(1,:)-X(q1))); % find the index in the calculated frequency spectrum corresponding the starting frequency (q1) determined above; we will shorten the frequency spectrum for the superposition of the wave components below
    [~,z2] = min(abs(ExcitationSpectrum(1,:)-X(q2))); % find the index corresponding the ending frequency (q2) determined above
    Fit2 = fit(X(q1:q2,1),X(q1:q2,4)*1e3,'cubicspline'); % fit the phase velocity
    Fit3 = fit(X(q1:q2,1),X(q1:q2,7),'cubicspline'); % fit the attenuation
    ExcitationSpectrumRange = {ExcitationSpectrum(:,z1:z2)}; % extract from the frequency spectrum the range for which we have phase velocities
    PhaseVelocity = {Fit2(ExcitationSpectrumRange{1}(1,:)')}; % interpolate the phase velocity
    Attenuation = {Fit3(ExcitationSpectrumRange{1}(1,:)')}; % interpolate the attenuation
    if  contains(Mode,'SH')
        SignFluids = {0};
        ModeType = 'Shear';
    else
        Fit4 = fit(X(q1:q2,1),X(q1:q2,8),'cubicspline');
        SignFluids = {round(Fit4(ExcitationSpectrumRange{1}(1,:)'))};
        ModeType = 'Lamb';
    end
end
function [uSum,Counter] = Computer(FluidLoading,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,DC,I,I1,Counter,Distance,ExcitationSpectrumRange,h1,ModeType,ModeTotal,PhaseVelocity,Attenuation,SignFluids,Direction,Time,SamplesX3,Layers,Height,a1,a11,a12,a21,a22,a23,a31,a32,a33,a34,Memory)
    global Stop
    Stop = 0;
    for p = 1:length(PhaseVelocity)
        if  any(PhaseVelocity{p})
            Length = length(ExcitationSpectrumRange{p});
            if  Length*length(Time)*8 > .8*Memory.MaxPossibleArrayBytes % Bytes needed for u = real(exp(1i*(reshape(Wavenumber,[],1)*Distance-reshape(AngularFrequency,[],1).*Time)));
                errordlg('The data size exceeds your physical RAM! Decrease the parameter settings.','Error');
                return
            end
            PhaseVelocity{p} = reshape(PhaseVelocity{p},1,1,[]);
            Attenuation{p} = reshape(Attenuation{p},1,1,[]);
            if  strcmp(ModeType,'Lamb')
                SignFluids{p} = reshape(SignFluids{p},1,1,[]);
            end
            AngularFrequency = reshape(ExcitationSpectrumRange{p}(1,:),1,1,[])*pi*2e3;
            AngularFrequency2 = AngularFrequency.^2;
            k = Direction(p)*(AngularFrequency./PhaseVelocity{p}+1i*Attenuation{p});
            k2 = k.^2;
            A = cell(0);
            if  ~Decoupled
                k4 = k2.^2;
                k6 = k2.^3;
                v = zeros(Height,3,Length);
                sigma = zeros(Height,3,Length);
                for m = 1:SuperLayerSize
                    rw2 = Material{m}.Density*AngularFrequency2;
                    if  DC(m)
                        A2 = a21(m)*k2+a22(m)*rw2;
                        A3 = a31(m)*k4+a32(m)*rw2.*k2+rw2.^2;
                        d1 = sqrt(A2.^2-2*a1(m)*A3);
                        k32 = [d1-A2 -d1-A2]/a1(m);
                        k3 = sqrt(k32);
                        k3SH = sqrt((rw2-k2*c{m}(6,6))/c{m}(4,4));
                        W = (rw2-c{m}(1,1)*k2-c{m}(5,5)*k32)./((c{m}(1,3)+c{m}(5,5))*k.*k3);
                        e1 = k.*W+k3;
                        e3 = k3.*W;
                        D3 = 1i*(c{m}(1,3)*k+c{m}(3,3)*e3);
                        D4 = 1i*c{m}(4,4)*k3SH;
                        D5 = 1i*c{m}(5,5)*e1;
                        E = exp(1i*k3*LayerThicknesses(m));
                        ESH = exp(1i*k3SH*LayerThicknesses(m));
                        A{m,3} = [ones(1,2,Length) zeros(1,1,Length) E zeros(1,1,Length);zeros(1,2,Length) ones(1,1,Length) zeros(1,2,Length) ESH;W zeros(1,1,Length) -W.*E zeros(1,1,Length);E zeros(1,1,Length) ones(1,2,Length) zeros(1,1,Length);zeros(1,2,Length) ESH zeros(1,2,Length) ones(1,1,Length);W.*E zeros(1,1,Length) -W zeros(1,1,Length)]; % L2
                        A{m,1} = pagemrdivide([D3 zeros(1,1,Length) D3.*E zeros(1,1,Length);D5 zeros(1,1,Length) -D5.*E zeros(1,1,Length);zeros(1,2,Length) D4 zeros(1,2,Length) -D4.*ESH;D3.*E zeros(1,1,Length) D3 zeros(1,1,Length);D5.*E zeros(1,1,Length) -D5 zeros(1,1,Length);zeros(1,2,Length) D4.*ESH zeros(1,2,Length) -D4],A{m,3}); % L
                        A{m,5} = [k3 k3SH];
                        A{m,6} = [zeros(1,2,Length) ones(1,1,Length)];
                        A{m,7} = [W zeros(1,1,Length)];
                        A{m,8} = [1i*(c{m}(1,1)*k+c{m}(1,3)*e3) zeros(1,1,Length)]; % sigma11
                        A{m,9} = [D5 zeros(1,1,Length)]; % sigma13
                        A{m,10} = [zeros(1,2,Length) 1i*k*c{m}(6,6)]; % sigma12
                    else
                        r2w4 = rw2.^2;
                        A1 = a11(m)*k2+a12(m)*rw2;
                        A2 = a21(m)*k4+a22(m)*rw2.*k2+a23(m)*r2w4;
                        A3 = a31(m)*k6+a32(m)*rw2.*k4+a33(m)*r2w4.*k2+a34(m)*rw2.^3;
                        d1 = A1/3;
                        d2 = A2/3-d1.^2;
                        d3 = d1.^3-d1.*A2/2+A3/2;
                        d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
                        d5 = d2./d4;
                        d6 = (d5-d4)/2-d1;
                        d7 = (d5+d4)/2i*sqrt(3);
                        k32 = [d6+d7 d6-d7 d4-d5-d1];
                        k3 = sqrt(k32);
                        k3k = k3.*k;
                        m11 = c{m}(1,1)*k2+c{m}(5,5)*k32-rw2;
                        m22 = c{m}(6,6)*k2+c{m}(4,4)*k32-rw2;
                        m12 = c{m}(1,6)*k2+c{m}(4,5)*k32;
                        m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
                        m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
                        m1 = m13.*m22-m12.*m23;
                        V = (m11.*m23-m13.*m12)./m1;
                        W = (m11.*m22-m12.^2)./-m1;
                        e1 = k.*W+k3;
                        e2 = k3.*V;
                        e3 = k3.*W;
                        D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V).*k+c{m}(3,3)*e3);
                        D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
                        D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
                        E = exp(1i*k3*LayerThicknesses(m));
                        A{m,3} = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W]; % L2
                        A{m,1} = pagemrdivide([D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4],A{m,3}); % L
                        A{m,5} = k3;
                        A{m,6} = V;
                        A{m,7} = W;
                        A{m,8} = 1i*((c{m}(1,1)+c{m}(1,6)*V).*k+c{m}(1,3)*e3); % sigma11
                        A{m,9} = D5; % sigma13
                        A{m,10} = 1i*((c{m}(1,6)+c{m}(6,6)*V).*k+c{m}(3,6)*e3); % sigma12
                    end
                    A{m,4} = LayerThicknesses(m);
                    A{m,11} = DC(m);
                end
                A = repmat(A,Repetitions,1);
                M{1} = A{1};
                for m = 2:SuperLayerSize
                    M0 = A{m}(1:3,1:3,:)-M{1}(4:6,4:6,:);
                    M1 = pagemrdivide(M{1}(1:3,4:6,:),M0);
                    M2 = pagemrdivide(A{m}(4:6,1:3,:),M0);
                    M{1} = [M{1}(1:3,1:3,:)+pagemtimes(M1,M{1}(4:6,1:3,:)) -pagemtimes(M1,A{m}(1:3,4:6,:));pagemtimes(M2,M{1}(4:6,1:3,:)) A{m}(4:6,4:6,:)-pagemtimes(M2,A{m}(1:3,4:6,:))];
                end
                for m = 1:length(Pattern)
                    M0 = M{Pattern(m)}(1:3,1:3,:)-M{m}(4:6,4:6,:);
                    M1 = pagemrdivide(M{m}(1:3,4:6,:),M0);
                    M2 = pagemrdivide(M{Pattern(m)}(4:6,1:3,:),M0);
                    M{m+1} = [M{m}(1:3,1:3,:)+pagemtimes(M1,M{m}(4:6,1:3,:)) -pagemtimes(M1,M{Pattern(m)}(1:3,4:6,:));pagemtimes(M2,M{m}(4:6,1:3,:)) M{Pattern(m)}(4:6,4:6,:)-pagemtimes(M2,M{Pattern(m)}(1:3,4:6,:))];
                end
                if  SymmetricSystem
                    A = [A;flipud(A)];
                    M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
                    M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
                    M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
                    M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
                end
                if  FluidLoading
                    M{end} = pageinv(M{end});
                    k3Fu = SignFluids{p}.*sqrt(AngularFrequency2/UpperFluid.Velocity^2-k2);
                    DFu = 1i*UpperFluid.Density*AngularFrequency2./k; % sigma11, sigma22, sigma33 in the upper fluid
                    DFl = 1i*LowerFluid.Density*AngularFrequency2./k; % in the lower fluid
                    UFl = -(k3Fu./k+M{end}(3,1,:).*DFu)./(M{end}(3,4,:).*DFl);
                    uInterfaces = M{end}(1:3,1,:).*DFu+M{end}(1:3,4,:).*DFl.*UFl;
                    uInterfaces(:,Layers+1,:) = M{end}(4:6,1,:).*DFu+M{end}(4:6,4,:).*DFl.*UFl;
                else
                    Z1 = [pagemtimes(-M{end}(1:3,1:3,:),[ones(1,3,Length);A{1,6};-A{1,7}]) pagemtimes(-M{end}(1:3,4:6,:),[ones(1,3,Length);A{end,6};A{end,7}]);pagemtimes(M{end}(4:6,1:3,:),[ones(1,3,Length);A{1,6};-A{1,7}]) pagemtimes(M{end}(4:6,4:6,:),[ones(1,3,Length);A{end,6};A{end,7}])];
                    Z2 = [pagemtimes(M{end}(1:3,1:3,:),[ones(1,3,Length);A{1,6};A{1,7}]);pagemtimes(-M{end}(4:6,1:3,:),[ones(1,3,Length);A{1,6};A{1,7}])];
                    RT = pagemldivide(Z1,Z2(:,1,:));
                    uInterfaces = [ones(1,1,Length);A{1,6}(1,1,:);A{1,7}(1,1,:)]+pagemtimes([ones(1,3,Length);A{1,6};-A{1,7}],RT(1:3,1,:));
                    uInterfaces(:,Layers+1,:) = pagemtimes([ones(1,3,Length);A{end,6};A{end,7}],RT(4:6,1,:));
                end
                A{1,2} = A{1};
                for m = 2:Layers
                    M0 = A{m}(1:3,1:3,:)-A{m-1,2}(4:6,4:6,:);
                    M1 = pagemrdivide(A{m-1,2}(1:3,4:6,:),M0);
                    M2 = pagemrdivide(A{m}(4:6,1:3,:),M0);
                    A{m,2} = [A{m-1,2}(1:3,1:3,:)+pagemtimes(M1,A{m-1,2}(4:6,1:3,:)) -pagemtimes(M1,A{m}(1:3,4:6,:));pagemtimes(M2,A{m-1,2}(4:6,1:3,:)) A{m}(4:6,4:6,:)-pagemtimes(M2,A{m}(1:3,4:6,:))];
                end
                for m = Layers:-1:2
                    M0 = A{m}(1:3,1:3,:)-A{m-1,2}(4:6,4:6,:);
                    uInterfaces(:,m,:) = pagemtimes(pagemldivide(M0,A{m-1,2}(4:6,1:3,:)),uInterfaces(:,1,:))-pagemtimes(pagemldivide(M0,A{m}(1:3,4:6,:)),uInterfaces(:,m+1,:));
                end
                for m = 1:Layers
                    x3 = (0:A{m,4}/SamplesX3:A{m,4})';
                    U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                    if  m == 1
                        U0 = U;
                        x30 = x3;
                        x3Total = x3;
                    else
                        x3Total = [x3Total;x3Total(end)+x3];
                    end
                    r = (m-1)*length(x3)+1:m*length(x3);
                    E = exp(1i*[A{m,5}.*x3 A{m,5}.*(A{m,4}-x3)]);
                    if  A{m,11}
                        v(r,1,:) = -1i*AngularFrequency.*pagemtimes([E(:,1:2,:) zeros(length(x3),1,Length) E(:,4:5,:) zeros(length(x3),1,Length)],U); % v1
                    else
                        v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v1
                    end
                    v(r,2,:) = -1i*AngularFrequency.*pagemtimes([A{m,6} A{m,6}].*E,U); % v2
                    v(r,3,:) = -1i*AngularFrequency.*pagemtimes([A{m,7} -A{m,7}].*E,U); % v3
                    sigma(r,1,:) = pagemtimes([A{m,8} A{m,8}].*E,U); % sigma11
                    sigma(r,3,:) = pagemtimes([A{m,9} -A{m,9}].*E,U); % sigma13
                    sigma(r,2,:) = pagemtimes([A{m,10} A{m,10}].*E,U); % sigma12
                end
            elseif Decoupled && strcmp(ModeType,'Lamb')
                k4 = k2.^2;
                v = zeros(Height,2,Length);
                sigma = zeros(Height,2,Length);
                for m = 1:SuperLayerSize
                    rw2 = Material{m}.Density*AngularFrequency2;
                    A2 = a21(m)*k2+a22(m)*rw2;
                    A3 = a31(m)*k4+a32(m)*rw2.*k2+rw2.^2;
                    d1 = sqrt(A2.^2-2*a1(m)*A3);
                    k32 = [d1-A2 -d1-A2]/a1(m);
                    k3 = sqrt(k32);
                    W = (rw2-c{m}(1,1)*k2-c{m}(5,5)*k32)./((c{m}(1,3)+c{m}(5,5))*k.*k3);
                    e3 = k3.*W;
                    D3 = 1i*(c{m}(1,3)*k+c{m}(3,3)*e3);
                    D5 = 1i*c{m}(5,5)*(k.*W+k3);
                    E = exp(1i*k3*LayerThicknesses(m));
                    A{m,3} = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W]; % L2
                    A{m,1} = pagemrdivide([D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5],A{m,3}); % L
                    A{m,4} = LayerThicknesses(m);
                    A{m,5} = k3;
                    A{m,7} = W;
                    A{m,6} = 1i*(c{m}(1,1)*k+c{m}(1,3)*e3); % sigma11
                    A{m,8} = D5; % sigma13
                end
                A = repmat(A,Repetitions,1);
                M{1} = A{1};
                for m = 2:SuperLayerSize
                    M0 = A{m}(1:2,1:2,:)-M{1}(3:4,3:4,:);
                    M1 = pagemrdivide(M{1}(1:2,3:4,:),M0);
                    M2 = pagemrdivide(A{m}(3:4,1:2,:),M0);
                    M{1} = [M{1}(1:2,1:2,:)+pagemtimes(M1,M{1}(3:4,1:2,:)) -pagemtimes(M1,A{m}(1:2,3:4,:));pagemtimes(M2,M{1}(3:4,1:2,:)) A{m}(3:4,3:4,:)-pagemtimes(M2,A{m}(1:2,3:4,:))];
                end
                for m = 1:length(Pattern)
                    M0 = M{Pattern(m)}(1:2,1:2,:)-M{m}(3:4,3:4,:);
                    M1 = pagemrdivide(M{m}(1:2,3:4,:),M0);
                    M2 = pagemrdivide(M{Pattern(m)}(3:4,1:2,:),M0);
                    M{m+1} = [M{m}(1:2,1:2,:)+pagemtimes(M1,M{m}(3:4,1:2,:)) -pagemtimes(M1,M{Pattern(m)}(1:2,3:4,:));pagemtimes(M2,M{m}(3:4,1:2,:)) M{Pattern(m)}(3:4,3:4,:)-pagemtimes(M2,M{Pattern(m)}(1:2,3:4,:))];
                end
                if  SymmetricSystem
                    A = [A;flipud(A)];
                    M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
                    M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
                    M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
                    M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
                end
                if  FluidLoading
                    M{end} = pageinv(M{end});
                    k3Fu = SignFluids{p}.*sqrt(AngularFrequency2/UpperFluid.Velocity^2-k2);
                    DFu = 1i*UpperFluid.Density*AngularFrequency2./k; % sigma11, sigma22, sigma33 in the upper fluid
                    DFl = 1i*LowerFluid.Density*AngularFrequency2./k; % in the lower fluid
                    UFl = -(k3Fu./k+M{end}(2,1,:).*DFu)./(M{end}(2,3,:).*DFl);
                    uInterfaces = M{end}(1:2,1,:).*DFu+M{end}(1:2,3,:).*DFl.*UFl;
                    uInterfaces(:,Layers+1,:) = M{end}(3:4,1,:).*DFu+M{end}(3:4,3,:).*DFl.*UFl;
                else
                    Z1 = [pagemtimes(-M{end}(1:2,1:2,:),[ones(1,2,Length);-A{1,7}]) pagemtimes(-M{end}(1:2,3:4,:),[ones(1,2,Length);A{end,7}]);pagemtimes(M{end}(3:4,1:2,:),[ones(1,2,Length);-A{1,7}]) pagemtimes(M{end}(3:4,3:4,:),[ones(1,2,Length);A{end,7}])];
                    Z2 = [pagemtimes(M{end}(1:2,1:2,:),[ones(1,2,Length);A{1,7}]);pagemtimes(-M{end}(3:4,1:2,:),[ones(1,2,Length);A{1,7}])];
                    RT = pagemldivide(Z1,Z2(:,1,:));
                    uInterfaces = [ones(1,1,Length);A{1,7}(1,1,:)]+pagemtimes([ones(1,2,Length);-A{1,7}],RT(1:2,1,:));
                    uInterfaces(:,Layers+1,:) = pagemtimes([ones(1,2,Length);A{end,7}],RT(3:4,1,:));
                end
                A{1,2} = A{1};
                for m = 2:Layers
                    M0 = A{m}(1:2,1:2,:)-A{m-1,2}(3:4,3:4,:);
                    M1 = pagemrdivide(A{m-1,2}(1:2,3:4,:),M0);
                    M2 = pagemrdivide(A{m}(3:4,1:2,:),M0);
                    A{m,2} = [A{m-1,2}(1:2,1:2,:)+pagemtimes(M1,A{m-1,2}(3:4,1:2,:)) -pagemtimes(M1,A{m}(1:2,3:4,:));pagemtimes(M2,A{m-1,2}(3:4,1:2,:)) A{m}(3:4,3:4,:)-pagemtimes(M2,A{m}(1:2,3:4,:))];
                end
                for m = Layers:-1:2
                    M0 = A{m}(1:2,1:2,:)-A{m-1,2}(3:4,3:4,:);
                    uInterfaces(:,m,:) = pagemtimes(pagemldivide(M0,A{m-1,2}(3:4,1:2,:)),uInterfaces(:,1,:))-pagemtimes(pagemldivide(M0,A{m}(1:2,3:4,:)),uInterfaces(:,m+1,:));
                end
                for m = 1:Layers
                    x3 = (0:A{m,4}/SamplesX3:A{m,4})';
                    U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                    if  m == 1
                        U0 = U;
                        x30 = x3;
                        x3Total = x3;
                    else
                        x3Total = [x3Total;x3Total(end)+x3];
                    end
                    r = (m-1)*length(x3)+1:m*length(x3);
                    E = exp(1i*[A{m,5}.*x3 A{m,5}.*(A{m,4}-x3)]);
                    v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v1
                    v(r,2,:) = -1i*AngularFrequency.*pagemtimes([A{m,7} -A{m,7}].*E,U); % v3
                    sigma(r,1,:) = pagemtimes([A{m,6} A{m,6}].*E,U); % sigma11
                    sigma(r,2,:) = pagemtimes([A{m,8} -A{m,8}].*E,U); % sigma13
                end
            elseif Decoupled && strcmp(ModeType,'Shear') && DisplacementComponent == 2
                v = zeros(Height,1,Length);
                sigma = zeros(Height,1,Length);
                for m = 1:SuperLayerSize
                    k3 = sqrt((Material{m}.Density*AngularFrequency2-k2*c{m}(6,6))/c{m}(4,4));
                    if  k3(1) == 0
                        k3(1,1,:) = 1e-10;
                    end
                    E = exp(1i*k3*LayerThicknesses(m));
                    E2 = E.^2;
                    A{m,1} = 1i*k3*c{m}(4,4)./(E2-ones(1,1,Length)).*[-ones(1,1,Length)-E2 2*E;-2*E ones(1,1,Length)+E2]; % L
                    A{m,3} = [ones(1,1,Length) E;E ones(1,1,Length)]; % L2
                    A{m,4} = LayerThicknesses(m);
                    A{m,5} = k3;
                    A{m,6} = 1i*k*c{m}(6,6); % sigma12
                end
                A = repmat(A,Repetitions,1);
                M{1} = A{1};
                for m = 2:SuperLayerSize
                    M0 = A{m}(1,1,:)-M{1}(2,2,:);
                    M1 = M{1}(1,2,:)./M0;
                    M2 = A{m}(2,1,:)./M0;
                    M{1} = [M{1}(1,1,:)+M1.*M{1}(2,1,:) -M1.*A{m}(1,2,:);M2.*M{1}(2,1,:) A{m}(2,2,:)-M2.*A{m}(1,2,:)];
                end
                for m = 1:length(Pattern)
                    M0 = M{Pattern(m)}(1,1,:)-M{m}(2,2,:);
                    M1 = M{m}(1,2,:)./M0;
                    M2 = M{Pattern(m)}(2,1,:)./M0;
                    M{m+1} = [M{m}(1,1,:)+M1.*M{m}(2,1,:) -M1.*M{Pattern(m)}(1,2,:);M2.*M{m}(2,1,:) M{Pattern(m)}(2,2,:)-M2.*M{Pattern(m)}(1,2,:)];
                end
                if  SymmetricSystem
                    A = [A;flipud(A)];
                    M1 = -M{end}(1,2,:)./(2*M{end}(2,2,:));
                    M{end} = [M{end}(1,1,:)+M1.*M{end}(2,1,:) M1.*M{end}(2,1,:);-M1.*M{end}(2,1,:) -M{end}(1,1,:)-M1.*M{end}(2,1,:)];
                end
                Z1 = [M{end}(1,1,:) -M{end}(1,2,:);M{end}(2,:,:)]; % changed sign in Z1(1,1)!
                Z2 = [M{end}(1,1,:);-M{end}(2,1,:)];
                RT = pagemldivide(Z1,Z2);
                uInterfaces = ones(1,1,Length)+RT(1,1,:);
                uInterfaces(:,Layers+1,:) = RT(2,1,:);
                A{1,2} = A{1};
                for m = 2:Layers
                    M0 = A{m}(1,1,:)-A{m-1,2}(2,2,:);
                    M1 = A{m-1,2}(1,2,:)./M0;
                    M2 = A{m}(2,1,:)./M0;
                    A{m,2} = [A{m-1,2}(1,1,:)+M1.*A{m-1,2}(2,1,:) -M1.*A{m}(1,2,:);M2.*A{m-1,2}(2,1,:) A{m}(2,2,:)-M2.*A{m}(1,2,:)];
                end
                for m = Layers:-1:2
                    M0 = A{m}(1,1,:)-A{m-1,2}(2,2,:);
                    uInterfaces(:,m,:) = M0.\A{m-1,2}(2,1,:).*uInterfaces(:,1,:)-M0.\A{m}(1,2,:).*uInterfaces(:,m+1,:);
                end
                for m = 1:Layers
                    x3 = (0:A{m,4}/SamplesX3:A{m,4})';
                    U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]); 
                    if  m == 1
                        U0 = U;
                        x30 = x3;
                        x3Total = x3;
                    else
                        x3Total = [x3Total;x3Total(end)+x3];
                    end
                    r = (m-1)*length(x3)+1:m*length(x3);
                    E = exp(1i*A{m,5}.*[x3 (A{m,4}-x3)]);
                    v(r,1,:) = -1i*AngularFrequency.*pagemtimes(E,U); % v2
                    sigma(r,1,:) = pagemtimes([A{m,6} A{m,6}].*E,U); % sigma12
                end
            end
            if  strcmp(ModeType,'Shear') && DisplacementComponent == 1
                uSum{p} = zeros(1,length(Time));
            else
                PowerFlow = trapz(x3Total,-real(sum(sigma.*conj(v),2))/2);
                E = exp(1i*(k*Distance+[A{1,5}.*x30(1:2) A{1,5}.*(A{1,4}-x30(1:2))]));
                if  DisplacementComponent == 1 % u3
                    unorm = pagemtimes([A{1,7} -A{1,7}].*E,U0)./sqrt(PowerFlow);
                elseif DisplacementComponent == 2 % u1,u2
                    unorm = pagemtimes(E,U0)./sqrt(PowerFlow);
                end
                unorm = squeeze(abs(real(unorm(1,1,:).*exp(-1i*angle(unorm(2,1,:))))));
                uNorm = fillmissing(filloutliers(unorm,'spline','movmedian',5,'ThresholdFactor',1),'spline');
                u = real(exp(1i*(squeeze(k)*Distance-squeeze(AngularFrequency).*Time-ExcitationSpectrumRange{p}(3,:)')));
                u = u./max(u,[],2).*uNorm.*ExcitationSpectrumRange{p}(2,:)';
                uSum{p} = sum(u);
                if  all(isnan(uSum{p}))
                    uSum{p} = zeros(1,length(Time));
                end
% figure
% hold on
% plot(ExcitationSpectrumRange{p}(1,:),reshape(unorm,1,Length),'linewidth',4,'color','r')
% plot(ExcitationSpectrumRange{p}(1,:),reshape(uNorm,1,Length),'linewidth',1.5,'color','g')
% figure
% hold on
% for i = 1:length(PhaseVelocity{p})
% plot(u(i,:))
% end
            end
            Counter = Counter+1;
            if  ModeTotal > 0
                waitbar(Counter/ModeTotal,h1,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
            end
        end
        if  Stop
            close(h1)
            return
        end
    end
end