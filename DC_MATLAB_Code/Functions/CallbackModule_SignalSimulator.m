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
function a = CallbackModule_SignalSimulator(source,~,a,Tab3)
%#ok<*GVMIS> 
global Stop 
Stop = 0;
if  strcmp(source.Tag,'2') % Frequency
    a.Frequency3 = str2double(source.String);
    if  a.DataType3 == 1
        if  a.Frequency3 > a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1)
            a.Frequency3 = a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1);
            source.String = a.Frequency3;
            errordlg(['Selected frequency outside frequency range! Select between 0 and ',num2str(a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1)),' kHz.'],'Warning')    
        end
        [a.ALambModes1,a.AShearModes1,a.SLambModes1,a.SShearModes1,a.Frequency3] = ModeFinder_Isotropic(a.ALamb1,a.AShear1,a.SLamb1,a.SShear1,a.Frequency3);
    elseif a.DataType3 == 2
        if  a.Frequency3 > a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2)
            a.Frequency3 = a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2);
            source.String = a.Frequency3;
            errordlg(['Selected frequency outside frequency range! Select between 0 and ',num2str(a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2)),' kHz.'],'Warning')    
        end
        [a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Frequency3] = ModeFinder_Anisotropic(a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.S,a.SLamb2,a.SShear2,a.Frequency3);
    end
    source.String = a.Frequency3;
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3 == 1
        a.CalculateUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    Colorization(a,2)
    ShowModes_Signal(a)
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'3') % Cycles
    a.Cycles3 = str2double(source.String);
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3 == 1
        a.CalculateUI3.Enable = 'on';        
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'4') % Samples per cycle
    a.SamplesPerCycle3 = str2double(source.String);
    if  mod(a.SamplesPerCycle3,2) ~= 0
        a.SamplesPerCycle3 = ceil(a.SamplesPerCycle3)+1;
        source.String = a.SamplesPerCycle3;
        errordlg('Enter only even numbers of samples per cycle!','Warning')
    end
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3 == 1
        a.CalculateUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'5') % Window
    a.Window3 = source.Value;
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3 == 1
        a.CalculateUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'6') % Distance
    a.Distance3 = str2double(source.String);
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3 == 1
        a.CalculateUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'7') % Time limit
    a.TimeLimitFactor3 = str2double(source.String);
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3 == 1
        a.CalculateUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'8') % Spectrum threshold
    a.SpectrumThreshold3 = str2double(source.String);
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3 == 1
        a.CalculateUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'113') % Displacement component
    a.DisplacementComponent3 = source.Value;
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3 == 1
        a.CalculateUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'111') % Gate
    a.Gate3 = eval(source.String);
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'9') % Multi-mode
    a.MultiMode3 = source.Value;
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  source.Value == 1
        a.CalculateUI3.Enable = 'on';
        a.p2UI3.Title = 'Mode selection and magnification';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''off'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''off'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''off'';',i))
        end
    else
        a.CalculateUI3.Enable = 'off';
        a.p2UI3.Title = 'Mode selection';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Enable = ''on'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''on'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''on'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''on'';',i))
        end       
    end
    for i = 0:9
        eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
        eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
        eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
    end
    a.Plot_ALamb3(1:10) = 0;
    a.Plot_SLamb3(1:10) = 0;
    a.Plot_AShear3(1:10) = 0;
    a.Plot_SShear3(1:10) = 0;
    Colorization(a,2)
    ShowModes_Signal(a)
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
elseif strcmp(source.Tag,'92') % Calculate
    a.CalculateUI3.Enable = 'off';
    a.StopUI3.Enable = 'on';
    if  a.DataType3 == 1
        [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
    elseif a.DataType3 == 2
        [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
    end
    if  Stop == 0
        a.Mode3 = '';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 1;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 1;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 1;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 1;',i))
            eval(sprintf('a.ALamb%ubUI3.Enable = ''on'';',i))
            eval(sprintf('a.SLamb%ubUI3.Enable = ''on'';',i))
            eval(sprintf('a.AShear%ubUI3.Enable = ''on'';',i+1))
            eval(sprintf('a.SShear%ubUI3.Enable = ''on'';',i))
        end
        a.Plot_ALamb3(1:10) = 1;
        a.Plot_SLamb3(1:10) = 1;
        a.Plot_AShear3(1:10) = 1;
        a.Plot_SShear3(1:10) = 1;
        Colorization(a,2)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end        
        a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
        a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
        a.YAxisUI3.String = a.YAxis3;
        a.PlotUI3.Enable = 'on';
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
    a.StopUI3.Enable = 'off';
elseif strcmp(source.Tag,'109') % Stop
    Stop = 1;
    a.CalculateUI3.Enable = 'on';
    a.PlotUI3.Enable = 'off';
elseif strcmp(source.Tag,'10') % A0,B0 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 0
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A0';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A0';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B0';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
        end
    else
        a.Plot_ALamb3(1) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'11') % A0,B0 amplitude
    a.Amplitude_ALamb3(1) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'12') % A1,B1 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 1
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A1';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A1';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B1';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(2) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);        
    end
elseif strcmp(source.Tag,'13') % A1,B1 amplitude
    a.Amplitude_ALamb3(2) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'14') % A2,B2 plot       
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 2
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A2';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A2';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B2';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(3) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'15') % A2,B2 amplitude
    a.Amplitude_ALamb3(3) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'16') % A3,B3 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 3
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A3';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A3';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B3';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(4) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'17') % A3,B3 amplitude
    a.Amplitude_ALamb3(4) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'18') % A4,B4 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 4
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A4';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A4';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B4';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(5) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'19') % A4,B4 amplitude
    a.Amplitude_ALamb3(5) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'20') % A5,B5 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 5
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A5';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A5';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B5';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(6) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'21') % A5,B5 amplitude
    a.Amplitude_ALamb3(6) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'22') % A6,B6 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 6
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A6';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A6';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B6';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(7) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'23') % A6,B6 amplitude
    a.Amplitude_ALamb3(7) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'24') % A7,B7 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 7
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A7';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A7';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B7';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(8) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'25') % A7,B7 amplitude
    a.Amplitude_ALamb3(8) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'26') % A8,B8 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 8
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A8';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A8';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B8';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(9) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'27') % A8,B8 amplitude
    a.Amplitude_ALamb3(9) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'28') % A9,B9 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 9
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'A9';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'A9';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                    a.Mode3 = 'B9';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_ALamb3(10) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'29') % A9,B9 amplitude
    a.Amplitude_ALamb3(10) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'30') % S0,B'0 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 0
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S0';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S0';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH0';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(1) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'31') % S0,B'0 amplitude
    a.Amplitude_SLamb3(1) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'32') % S1,B'1 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 1
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S1';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S1';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH1';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(2) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'33') % S1,B'1 amplitude
    a.Amplitude_SLamb3(2) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);    
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'34') % S2,B'2 plot       
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 2
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S2';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S2';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH2';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(3) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'35') % S2,B'2 amplitude
    a.Amplitude_SLamb3(3) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'36') % S3,B'3 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 3
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S3';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S3';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH3';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(4) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'37') % S3,B'3 amplitude
    a.Amplitude_SLamb3(4) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'38') % S4,B'4 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 4
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S4';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S4';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH4';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(5) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'39') % S4,B'4 amplitude
    a.Amplitude_SLamb3(5) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'40') % S5,B'5 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 5
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S5';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S5';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH5';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(6) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'41') % S5,B'5 amplitude
    a.Amplitude_SLamb3(6) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'42') % S6,B'6 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 6
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S6';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S6';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH6';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(7) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'43') % S6,B'6 amplitude
    a.Amplitude_SLamb3(7) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'44') % S7,B'7 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 7
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S7';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S7';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH7';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(8) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'45') % S7,B'7 amplitude
    a.Amplitude_SLamb3(8) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'46') % S8,B'8 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 8
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S8';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S8';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH8';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(9) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'47') % S8,B'8 amplitude
    a.Amplitude_SLamb3(9) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'48') % S9,B'9 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 9
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            if  a.DataType3 == 1
                a.Mode3 = 'S9';
            elseif a.DataType3 == 2
                if  a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                    a.Mode3 = 'S9';
                elseif a.SuperLayerSize > 1 && a.SymmetricSystem == 0 && a.Decoupled == 1
                    a.Mode3 = 'BSH9';
                end
            end
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SLamb3(10) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'49') % S9,B'9 amplitude
    a.Amplitude_SLamb3(10) = str2double(source.String); 
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'50') % A'1 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 0
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH1';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(1) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'51') % A'1 amplitude
    a.Amplitude_AShear3(1) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'52') % A'2 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 1
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH2';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(2) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'53') % A'2 amplitude
    a.Amplitude_AShear3(2) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'54') % A'3 plot       
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 2
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH3';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(3) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'55') % A'3 amplitude
    a.Amplitude_AShear3(3) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'56') % A'4 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 3
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH4';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(4) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'57') % A'4 amplitude
    a.Amplitude_AShear3(4) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'58') % A'5 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 4
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH5';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(5) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'59') % A'5 amplitude
    a.Amplitude_AShear3(5) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'60') % A'6 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 5
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH6';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(6) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'61') % A'6 amplitude
    a.Amplitude_AShear3(6) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'62') % A'7 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 6
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH7';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(7) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'63') % A'7 amplitude
    a.Amplitude_AShear3(7) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'64') % A'8 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 7
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH8';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(8) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'65') % A'8 amplitude
    a.Amplitude_AShear3(8) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'66') % A'9 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 8
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH9';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(9) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'67') % A'9 amplitude
    a.Amplitude_AShear3(9) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);   
elseif strcmp(source.Tag,'68') % A'10 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 9
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value == 1
            a.Mode3 = 'ASH10';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_AShear3(10) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'69') % A'10 amplitude
    a.Amplitude_AShear3(10) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'70') % S'0 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 0
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH0';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SShear3(1) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'71') % S'0 amplitude
    a.Amplitude_SShear3(1) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'72') % S'1 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 1
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH1';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SShear3(2) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'73') % S'1 amplitude
    a.Amplitude_SShear3(2) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);    
elseif strcmp(source.Tag,'74') % S'2 plot       
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 2
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH2';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SShear3(3) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'75') % S'2 amplitude
    a.Amplitude_SShear3(3) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'76') % S'3 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 3
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH3';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SShear3(4) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'77') % S'3 amplitude
    a.Amplitude_SShear3(4) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'78') % S'4 plot
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 4
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH4';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SShear3(5) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'79') % S'4 amplitude
    a.Amplitude_SShear3(5) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'80') % S'5 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 5
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH5';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SShear3(6) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'81') % S'5 amplitude
    a.Amplitude_SShear3(6) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'82') % S'6 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 6
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH6';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
        end
    else
        a.Plot_SShear3(7) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'83') % S'6 amplitude
    a.Amplitude_SShear3(7) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'84') % S'7 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 7
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH7';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SShear3(8) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'85') % S'7 amplitude
    a.Amplitude_SShear3(8) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'86') % S'8 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 8
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH8';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1); 
        end
    else
        a.Plot_SShear3(9) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'87') % S'8 amplitude
    a.Amplitude_SShear3(9) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'88') % S'9 plot        
    if  a.MultiMode3 == 0
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 9
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value == 1
            a.Mode3 = 'SSH9';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Isotropic_Signal(a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.Thickness,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.c,a.Material2,a.Mode3,a.S,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
                a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
        end
    else
        a.Plot_SShear3(10) = source.Value;
        Colorization(a,1)
        if  a.DataType3 == 1
            a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
        elseif a.DataType3 == 2
            a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
        end
        a.YAxisUI3.String = a.YAxis3;
        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
    end
elseif strcmp(source.Tag,'89') % S'9 amplitude
    a.Amplitude_SShear3(10) = str2double(source.String);
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax_Isotropic(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumAShear3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax_Anisotropic(a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.Symmetric,a.Decoupled,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'90') % X-axis
    a.XAxis3 = eval(source.String);
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'112') % Y-axis
    a.YAxis3 = eval(source.String);
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'91') % Plot
    if  a.DataType3 == 1
        Signal_Isotropic_External(a.CropPlots3,a.FluidLoading1,a.Fluid1,a.DisplacementComponent3,a.Amplitude_ALamb3,a.Amplitude_AShear3,a.ALamb1,a.AShear1,a.Amplitude_SLamb3,a.Amplitude_SShear3,a.BoxLineWidth3,a.Directory3,a.Distance3,a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.ExportPlots3,a.FileName3,a.AxesTickFontSize3,a.AxesLabelFontSize3,a.TitleFontSize3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimit1,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.Title3,a.LineColors3,a.LineWidth3,a.Material1,a.Mode3,a.MultiMode3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.Plot_ALamb3,a.Plot_AShear3,a.Plot_SLamb3,a.Plot_SShear3,a.PDF3,a.PNG3,a.PlotXAxis3,a.PNGresolution3,a.Gate3,a.SLamb1,a.SShear1,a.Thickness,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3,a.XAxis3,a.YAxis3)
    elseif a.DataType3 == 2
        Signal_Anisotropic_External(a.Hybrid,a.CropPlots3,a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.S,a.SLamb2,a.SShear2,a.PhaseVelocityA3,a.PhaseVelocityALamb3,a.PhaseVelocityAShear3,a.PhaseVelocityB3,a.PhaseVelocityBLamb3,a.PhaseVelocityBShear3,a.PhaseVelocityS3,a.PhaseVelocitySLamb3,a.PhaseVelocitySShear3,a.PropagationAngle,a.LayupString1,a.SuperLayers,a.SymmetricSystem,a.Symmetric,a.SuperLayerSize,a.Decoupled,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.BoxLineWidth3,a.Directory3,a.Distance3,a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.ExportPlots3,a.FileName3,a.AxesTickFontSize3,a.AxesLabelFontSize3,a.TitleFontSize3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimit2,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.Title3,a.LineColors3,a.LineWidth3,a.Material2,a.Mode3,a.MultiMode3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.PDF3,a.PNG3,a.PlotXAxis3,a.PNGresolution3,a.Gate3,a.PlateThickness/1e3,a.uSum3,a.uSumA3,a.uSumALamb3,a.uSumAShear3,a.uSumB3,a.uSumBLamb3,a.uSumBShear3,a.uSumS3,a.uSumSLamb3,a.uSumSShear3,a.XAxis3,a.YAxis3)        
    end
elseif strcmp(source.Tag,'93') % Export plots
    a.ExportPlots3 = source.Value;
    if  source.Value == 1
        a.PlotUI3.String = 'Export';
    else
        a.PlotUI3.String = 'Plot';
    end
elseif strcmp(source.Tag,'110') % Crop plots
    a.CropPlots3 = source.Value;
elseif strcmp(source.Tag,'94') % PDF
    a.PDF3 = source.Value;
elseif strcmp(source.Tag,'95') % PNG
    a.PNG3 = source.Value;
elseif strcmp(source.Tag,'96') % PNG resolution
    a.PNGresolution3 = str2double(source.String);
elseif strcmp(source.Tag,'97') % Matlab
    Export_Signal(a,0,0,1)  
elseif strcmp(source.Tag,'98') % Excel
    Export_Signal(a,1,0,0) 
elseif strcmp(source.Tag,'99') % TXT
    Export_Signal(a,0,1,0) 
elseif strcmp(source.Tag,'100') % File name
    a.FileName3 = source.String;
elseif strcmp(source.Tag,'101') % Directory
    a.Directory3 = source.String;
elseif strcmp(source.Tag,'102') % Title
    if  a.DataType3 == 1
        a.Title3 = source.Value;
    elseif a.DataType3 == 2
        switch source.Value
        case 1
           a.Title3 = 2;
        case 2
           a.Title3 = 1;
        case 3
           a.Title3 = 0;
        end
    end
elseif strcmp(source.Tag,'103') % Box line width
    a.BoxLineWidth3 = str2double(source.String);
elseif strcmp(source.Tag,'104') % Curve line width
    a.LineWidth3 = str2double(source.String);
elseif strcmp(source.Tag,'105') % Title font size
    a.TitleFontSize3 = str2double(source.String);
elseif strcmp(source.Tag,'106') % Axes labels font size
    a.AxesLabelFontSize3 = str2double(source.String);
elseif strcmp(source.Tag,'107') % Axes ticks font size
    a.AxesTickFontSize3 = str2double(source.String);
elseif strcmp(source.Tag,'108') % Default
    if  a.DataType3 == 1
        a.Title3 = 1;
    elseif a.DataType3 == 2
        a.Title3 = 2;
    end
    a.TitleUI3.Value = 1;
    a.BoxLineWidth3 = .5;
    a.BoxLineWidthUI3.String = a.BoxLineWidth3;
    a.LineWidth3 = 1;
    a.LineWidthUI3.String = a.LineWidth3;
    a.TitleFontSize3 = 30;
    a.TitleFontSizeUI3.String = a.TitleFontSize3;
    a.AxesLabelFontSize3 = 30;
    a.AxesLabelFontSizeUI3.String = a.AxesLabelFontSize3;
    a.AxesTickFontSize3 = 24;
    a.AxesTickFontSizeUI3.String = a.AxesTickFontSize3;
end