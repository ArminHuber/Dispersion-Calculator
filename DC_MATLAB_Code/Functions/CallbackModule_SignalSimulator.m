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
function a = CallbackModule_SignalSimulator(source,~,a,Tab3)
%#ok<*GVMIS>
%#ok<*FXUP>
global Stop 
Stop = 0;
source.String = replace(source.String,',','.');
if  strcmp(source.Tag,'2') % Frequency
    a.Frequency3 = str2double(source.String);
    if  a.DataType3 == 1
        if  a.Frequency3 > a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1)
            a.Frequency3 = a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1);
            errordlg(['Selected frequency outside frequency range! Select between 0 and ',num2str(a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1)),' kHz.'],'Warning')    
        elseif a.Frequency3 < a.FrequencyResolution1
            a.Frequency3 = a.FrequencyResolution1;
        end
        [a.ALambModes1,a.SLambModes1,a.BLambModes1,a.AShearModes1,a.SShearModes1,a.BShearModes1,a.Frequency3] = ModeFinder(a.ALamb1,a.SLamb1,a.BLamb1,a.AShear1,a.SShear1,{0},a.Frequency3);
    elseif a.DataType3 == 2
        if  a.Frequency3 > a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2)
            a.Frequency3 = a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2);
            errordlg(['Selected frequency outside frequency range! Select between 0 and ',num2str(a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2)),' kHz.'],'Warning')    
        elseif a.Frequency3 < a.FrequencyResolution2
            a.Frequency3 = a.FrequencyResolution2;
        end
        [a.ALambModes2,a.SLambModes2,a.BLambModes2,a.AShearModes2,a.SShearModes2,a.BShearModes2,a.Frequency3] = ModeFinder(a.ALamb2,a.SLamb2,a.BLamb2,a.AShear2,a.SShear2,a.BShear2,a.Frequency3);
    end
    source.String = a.Frequency3;
    Colorization(a,2)
    ShowModes_Signal(a)
    HelperFunction1
elseif strcmp(source.Tag,'3') % Cycles
    a.Cycles3 = str2double(source.String);
    HelperFunction1
elseif strcmp(source.Tag,'4') % Samples per cycle
    a.SamplesPerCycle3 = str2double(source.String);
    if  mod(a.SamplesPerCycle3,2) ~= 0
        a.SamplesPerCycle3 = ceil(a.SamplesPerCycle3)+1;
        source.String = a.SamplesPerCycle3;
        errordlg('Enter only even numbers of samples per cycle!','Warning')
    end
    HelperFunction1
elseif strcmp(source.Tag,'5') % Window
    a.Window3 = source.Value;
    HelperFunction1
elseif strcmp(source.Tag,'6') % Distance
    a.Distance3 = str2double(source.String);
    HelperFunction1
elseif strcmp(source.Tag,'7') % Time limit
    a.TimeLimitFactor3 = str2double(source.String);
    HelperFunction1
elseif strcmp(source.Tag,'8') % Spectrum threshold
    a.SpectrumThreshold3 = str2double(source.String);
    HelperFunction1
elseif strcmp(source.Tag,'113') % Displacement component
    a.DisplacementComponent3 = source.Value;
    HelperFunction1
elseif strcmp(source.Tag,'111') % Gate
    a.Gate3 = eval(source.String);
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'9') % Multi-mode
    a.MultiMode3 = source.Value;
    if  source.Value
        a.p2UI3.Title = 'Mode selection and magnification';
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
    Colorization(a,2)
    ShowModes_Signal(a)
    HelperFunction1
elseif strcmp(source.Tag,'92') % Calculate/Stop
    if  source.Value && strcmp(source.String,'Calculate')
        try
            source.String = 'Stop';
            source.TooltipString = 'Stop calculating.';
            if  a.DataType3 == 1
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3,a.uSumBLamb3,a.uSumBShear3] = Computer_Isotropic_Signal(a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.UpperFluid1,a.LowerFluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.BLamb1,a.BLambModes1,a.Thickness1/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
            elseif a.DataType3 == 2
                [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumBLamb3,a.uSumBShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.c,a.Material2,a.Mode3,a.SLamb2,a.SShear2,str2double(a.LayerCountUI2.String),a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.Delta,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.ALambModes2,a.AShearModes2,a.BLambModes2,a.BShearModes2,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
            end
            source.Value = 0;
            source.String = 'Calculate';
            source.TooltipString = 'Calculate multiple modes at a time.';
            if  Stop
                return
            end
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
                a.YAxis3 = SignalMax(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            elseif a.DataType3 == 2
                a.YAxis3 = SignalMax(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
            end
            a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
            a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
            a.YAxisUI3.String = a.YAxis3;
            a.PlotUI3.Enable = 'on';
            a.CalculateUI3.Enable = 'off';
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
        catch
            source.Value = 0;
            source.String = 'Calculate';
            source.TooltipString = 'Calculate multiple modes at a time.';
        end
    else
        Stop = 1;
        a.PlotUI3.Enable = 'off';
    end
elseif strcmp(source.Tag,'10') % A0,B0 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 0
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A0';
            else
                a.Mode3 = 'B0';
            end
            HelperFunction2
        end
    else
        a.Plot_ALamb3(1) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'11') % A0,B0 amplitude
    a.Amplitude_ALamb3(1) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'12') % A1,B1 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 1
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A1';
            else
                a.Mode3 = 'B1';
            end
            HelperFunction2
        end
    else
        a.Plot_ALamb3(2) = source.Value;
        HelperFunction3        
    end
elseif strcmp(source.Tag,'13') % A1,B1 amplitude
    a.Amplitude_ALamb3(2) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'14') % A2,B2 plot       
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 2
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A2';
            else
                a.Mode3 = 'B2';
            end
            HelperFunction2 
        end
    else
        a.Plot_ALamb3(3) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'15') % A2,B2 amplitude
    a.Amplitude_ALamb3(3) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'16') % A3,B3 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 3
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A3';
            else
                a.Mode3 = 'B3';
            end
            HelperFunction2 
        end
    else
        a.Plot_ALamb3(4) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'17') % A3,B3 amplitude
    a.Amplitude_ALamb3(4) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'18') % A4,B4 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 4
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A4';
            else
                a.Mode3 = 'B4';
            end
            HelperFunction2 
        end
    else
        a.Plot_ALamb3(5) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'19') % A4,B4 amplitude
    a.Amplitude_ALamb3(5) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'20') % A5,B5 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 5
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A5';
            else
                a.Mode3 = 'B5';
            end
            HelperFunction2 
        end
    else
        a.Plot_ALamb3(6) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'21') % A5,B5 amplitude
    a.Amplitude_ALamb3(6) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'22') % A6,B6 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 6
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A6';
            else
                a.Mode3 = 'B6';
            end
            HelperFunction2 
        end
    else
        a.Plot_ALamb3(7) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'23') % A6,B6 amplitude
    a.Amplitude_ALamb3(7) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'24') % A7,B7 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 7
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A7';
            else
                a.Mode3 = 'B7';
            end
            HelperFunction2 
        end
    else
        a.Plot_ALamb3(8) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'25') % A7,B7 amplitude
    a.Amplitude_ALamb3(8) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'26') % A8,B8 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 8
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A8';
            else
                a.Mode3 = 'B8';
            end
            HelperFunction2 
        end
    else
        a.Plot_ALamb3(9) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'27') % A8,B8 amplitude
    a.Amplitude_ALamb3(9) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'28') % A9,B9 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            if  i ~= 9
                eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  (a.DataType3 == 1 && a.Symmetric1) || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'A9';
            else
                a.Mode3 = 'B9';
            end
            HelperFunction2 
        end
    else
        a.Plot_ALamb3(10) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'29') % A9,B9 amplitude
    a.Amplitude_ALamb3(10) = str2double(source.String);
    HelperFunction4   
elseif strcmp(source.Tag,'30') % S0,BSH0 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 0
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S0';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH0';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(1) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'31') % S0,BSH0 amplitude
    a.Amplitude_SLamb3(1) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'32') % S1,BSH1 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 1
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S1';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH1';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(2) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'33') % S1,BSH1 amplitude
    a.Amplitude_SLamb3(2) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'34') % S2,BSH2 plot       
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 2
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S2';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH2';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(3) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'35') % S2,BSH2 amplitude
    a.Amplitude_SLamb3(3) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'36') % S3,BSH3 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 3
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S3';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH3';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(4) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'37') % S3,BSH3 amplitude
    a.Amplitude_SLamb3(4) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'38') % S4,BSH4 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 4
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S4';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH4';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(5) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'39') % S4,BSH4 amplitude
    a.Amplitude_SLamb3(5) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'40') % S5,BSH5 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 5
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S5';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH5';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(6) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'41') % S5,BSH5 amplitude
    a.Amplitude_SLamb3(6) = str2double(source.String);
    HelperFunction4   
elseif strcmp(source.Tag,'42') % S6,BSH6 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 6
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S6';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH6';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(7) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'43') % S6,BSH6 amplitude
    a.Amplitude_SLamb3(7) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'44') % S7,BSH7 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 7
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S7';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH7';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(8) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'45') % S7,BSH7 amplitude
    a.Amplitude_SLamb3(8) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'46') % S8,BSH8 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 8
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S8';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH8';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(9) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'47') % S8,BSH8 amplitude
    a.Amplitude_SLamb3(9) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'48') % S9,BSH9 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            if  i ~= 9
                eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            end
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            if  a.DataType3 == 1 || (a.DataType3 == 2 && a.Symmetric2)
                a.Mode3 = 'S9';
            elseif a.DataType3 == 2 && ~a.Symmetric2 && a.Decoupled
                a.Mode3 = 'BSH9';
            end
            HelperFunction2 
        end
    else
        a.Plot_SLamb3(10) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'49') % S9,BSH9 amplitude
    a.Amplitude_SLamb3(10) = str2double(source.String); 
    HelperFunction4    
elseif strcmp(source.Tag,'50') % ASH1 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 0
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH1';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(1) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'51') % ASH1 amplitude
    a.Amplitude_AShear3(1) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'52') % ASH2 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 1
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH2';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(2) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'53') % ASH2 amplitude
    a.Amplitude_AShear3(2) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'54') % ASH3 plot       
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 2
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH3';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(3) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'55') % ASH3 amplitude
    a.Amplitude_AShear3(3) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'56') % ASH4 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 3
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH4';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(4) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'57') % ASH4 amplitude
    a.Amplitude_AShear3(4) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'58') % ASH5 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 4
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH5';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(5) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'59') % ASH5 amplitude
    a.Amplitude_AShear3(5) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'60') % ASH6 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 5
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH6';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(6) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'61') % ASH6 amplitude
    a.Amplitude_AShear3(6) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'62') % ASH7 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 6
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH7';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(7) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'63') % ASH7 amplitude
    a.Amplitude_AShear3(7) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'64') % ASH8 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 7
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH8';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(8) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'65') % ASH8 amplitude
    a.Amplitude_AShear3(8) = str2double(source.String);
    HelperFunction4   
elseif strcmp(source.Tag,'66') % ASH9 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 8
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH9';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(9) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'67') % ASH9 amplitude
    a.Amplitude_AShear3(9) = str2double(source.String);
    HelperFunction4   
elseif strcmp(source.Tag,'68') % ASH10 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            if  i ~= 9
                eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            end
            eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
        end
        if  source.Value
            a.Mode3 = 'ASH10';
            HelperFunction2 
        end
    else
        a.Plot_AShear3(10) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'69') % ASH10 amplitude
    a.Amplitude_AShear3(10) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'70') % SSH0 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 0
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH0';
            HelperFunction2 
        end
    else
        a.Plot_SShear3(1) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'71') % SSH0 amplitude
    a.Amplitude_SShear3(1) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'72') % SSH1 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 1
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH1';
            HelperFunction2 
        end
    else
        a.Plot_SShear3(2) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'73') % SSH1 amplitude
    a.Amplitude_SShear3(2) = str2double(source.String);
    HelperFunction4    
elseif strcmp(source.Tag,'74') % SSH2 plot       
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 2
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH2';
            HelperFunction2 
        end
    else
        a.Plot_SShear3(3) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'75') % SSH2 amplitude
    a.Amplitude_SShear3(3) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'76') % SSH3 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 3
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH3';
            HelperFunction2 
        end
    else
        a.Plot_SShear3(4) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'77') % SSH3 amplitude
    a.Amplitude_SShear3(4) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'78') % SSH4 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 4
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH4';
            HelperFunction2 
        end
    else
        a.Plot_SShear3(5) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'79') % SSH4 amplitude
    a.Amplitude_SShear3(5) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'80') % SSH5 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 5
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH5';
            HelperFunction2 
        end
    else
        a.Plot_SShear3(6) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'81') % SSH5 amplitude
    a.Amplitude_SShear3(6) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'82') % SSH6 plot
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 6
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH6';
            HelperFunction2
        end
    else
        a.Plot_SShear3(7) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'83') % SSH6 amplitude
    a.Amplitude_SShear3(7) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'84') % SSH7 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 7
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH7';
            HelperFunction2 
        end
    else
        a.Plot_SShear3(8) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'85') % SSH7 amplitude
    a.Amplitude_SShear3(8) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'86') % SSH8 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 8
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH8';
            HelperFunction2 
        end
    else
        a.Plot_SShear3(9) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'87') % SSH8 amplitude
    a.Amplitude_SShear3(9) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'88') % SSH9 plot        
    if  ~a.MultiMode3
        a.PlotUI3.Enable = 'on';
        for i = 0:9
            eval(sprintf('a.ALamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.SLamb%ubUI3.Value = 0;',i))
            eval(sprintf('a.AShear%ubUI3.Value = 0;',i+1))
            if  i ~= 9
                eval(sprintf('a.SShear%ubUI3.Value = 0;',i))
            end
        end
        if  source.Value
            a.Mode3 = 'SSH9';
            HelperFunction2
        end
    else
        a.Plot_SShear3(10) = source.Value;
        HelperFunction3
    end
elseif strcmp(source.Tag,'89') % SSH9 amplitude
    a.Amplitude_SShear3(10) = str2double(source.String);
    HelperFunction4
elseif strcmp(source.Tag,'90') % X-axis
    a.XAxis3 = eval(source.String);
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'112') % Y-axis
    a.YAxis3 = eval(source.String);
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
elseif strcmp(source.Tag,'91') % Plot
    if  a.DataType3 == 1
        Signal_External(0,a.FluidLoading1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.UpperFluid1,a.LowerFluid1,a.DisplacementComponent3,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,0,'',0,0,1,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.BoxLineWidth3,a.Directory3,a.Distance3,a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.ExportPlots3,a.FileName3,a.AxesTickFontSize3,a.AxesLabelFontSize3,a.TitleFontSize3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimit1,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.Title3,a.LineColors3,a.LineWidth3,a.Material1,a.Mode3,a.MultiMode3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.PDF3,a.PNG3,a.PlotXAxis3,a.PNGresolution3,a.Gate3,a.Thickness1,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumBLamb3,a.uSumBShear3,a.uSumSLamb3,a.uSumSShear3,a.XAxis3,a.YAxis3)
    elseif a.DataType3 == 2
        Signal_External(a.Hybrid,a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.PropagationAngle,a.LayupString1,a.SuperLayers,a.SymmetricSystem,a.Decoupled,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.BoxLineWidth3,a.Directory3,a.Distance3,a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.ExportPlots3,a.FileName3,a.AxesTickFontSize3,a.AxesLabelFontSize3,a.TitleFontSize3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimit2,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.Title3,a.LineColors3,a.LineWidth3,a.Material2,a.Mode3,a.MultiMode3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.PDF3,a.PNG3,a.PlotXAxis3,a.PNGresolution3,a.Gate3,a.PlateThickness/1e3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumBLamb3,a.uSumBShear3,a.uSumSLamb3,a.uSumSShear3,a.XAxis3,a.YAxis3)        
    end
elseif strcmp(source.Tag,'93') % Export plots
    a.ExportPlots3 = source.Value;
    if  source.Value
        a.PlotUI3.String = 'Export';
    else
        a.PlotUI3.String = 'Plot';
    end
elseif strcmp(source.Tag,'94') % PDF
    a.PDF3 = source.Value;
elseif strcmp(source.Tag,'95') % PNG
    a.PNG3 = source.Value;
elseif strcmp(source.Tag,'96') % PNG resolution
    a.PNGresolution3 = str2double(source.String);
elseif strcmp(source.Tag,'97') % Matlab
    if  a.DataType3 == 1
        Export_Signal(0,0,1,a.ALambModes1,a.SLambModes1,a.BLambModes1,a.AShearModes1,a.SShearModes1,a.BShearModes1,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.FrequencyRange3,a.ExcitationSignal3,a.ExcitationMagnitude3,a.Gate3,a.FourierTransformLength3,a.PlotXAxis3,a.MultiMode3,a.Mode3,a.Frequency3,a.Thickness1/1e3,a.Distance3,a.Directory3,a.FileName3)
    elseif a.DataType3 == 2
        Export_Signal(0,0,1,a.ALambModes2,a.SLambModes2,a.BLambModes2,a.AShearModes2,a.SShearModes2,a.BShearModes2,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.FrequencyRange3,a.ExcitationSignal3,a.ExcitationMagnitude3,a.Gate3,a.FourierTransformLength3,a.PlotXAxis3,a.MultiMode3,a.Mode3,a.Frequency3,a.PlateThickness/1e3,a.Distance3,a.Directory3,a.FileName3)
    end
elseif strcmp(source.Tag,'98') % Excel
    if  a.DataType3 == 1
        Export_Signal(1,0,0,a.ALambModes1,a.SLambModes1,a.BLambModes1,a.AShearModes1,a.SShearModes1,a.BShearModes1,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.FrequencyRange3,a.ExcitationSignal3,a.ExcitationMagnitude3,a.Gate3,a.FourierTransformLength3,a.PlotXAxis3,a.MultiMode3,a.Mode3,a.Frequency3,a.Thickness1/1e3,a.Distance3,a.Directory3,a.FileName3)
    elseif a.DataType3 == 2
        Export_Signal(1,0,0,a.ALambModes2,a.SLambModes2,a.BLambModes2,a.AShearModes2,a.SShearModes2,a.BShearModes2,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.FrequencyRange3,a.ExcitationSignal3,a.ExcitationMagnitude3,a.Gate3,a.FourierTransformLength3,a.PlotXAxis3,a.MultiMode3,a.Mode3,a.Frequency3,a.PlateThickness/1e3,a.Distance3,a.Directory3,a.FileName3)
    end
elseif strcmp(source.Tag,'99') % TXT
    if  a.DataType3 == 1
        Export_Signal(0,1,0,a.ALambModes1,a.SLambModes1,a.BLambModes1,a.AShearModes1,a.SShearModes1,a.BShearModes1,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.FrequencyRange3,a.ExcitationSignal3,a.ExcitationMagnitude3,a.Gate3,a.FourierTransformLength3,a.PlotXAxis3,a.MultiMode3,a.Mode3,a.Frequency3,a.Thickness1/1e3,a.Distance3,a.Directory3,a.FileName3)
    elseif a.DataType3 == 2
        Export_Signal(0,1,0,a.ALambModes2,a.SLambModes2,a.BLambModes2,a.AShearModes2,a.SShearModes2,a.BShearModes2,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.FrequencyRange3,a.ExcitationSignal3,a.ExcitationMagnitude3,a.Gate3,a.FourierTransformLength3,a.PlotXAxis3,a.MultiMode3,a.Mode3,a.Frequency3,a.PlateThickness/1e3,a.Distance3,a.Directory3,a.FileName3)
    end
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
function HelperFunction1
    a.Mode3 = '';
    a.PlotUI3.Enable = 'off';
    if  a.MultiMode3
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
end
function HelperFunction2
    if  a.DataType3 == 1
        [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumSLamb3,a.uSumSShear3,a.uSumBLamb3,a.uSumBShear3] = Computer_Isotropic_Signal(a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.UpperFluid1,a.LowerFluid1,a.DisplacementComponent3,a.ALamb1,a.ALambModes1,a.AShear1,a.AShearModes1,a.Cycles3,a.Distance3,a.Frequency3,a.FrequencyResolution1,a.Material1,a.Mode3,a.MultiMode3,a.SamplesPerCycle3,a.SLamb1,a.SLambModes1,a.SpectrumThreshold3,a.SShear1,a.SShearModes1,a.BLamb1,a.BLambModes1,a.Thickness1/1e3,a.TimeLimitFactor3,a.Window3,a.OutputWindowUI3);
        a.YAxis3 = SignalMax(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        [a.ExcitationMagnitude3,a.ExcitationSignal3,a.ExcitationSpectrum3,a.FourierTransformLength3,a.Frequency3,a.FrequencyLimitLow3,a.FrequencyLimitHigh3,a.FrequencyRange3,a.PlotXAxis3,a.Gate3,a.XAxis3,a.uSum3,a.uSumALamb3,a.uSumAShear3,a.uSumBLamb3,a.uSumBShear3,a.uSumSLamb3,a.uSumSShear3] = Computer_Anisotropic_Signal(a.FluidLoading2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.UpperFluid2,a.LowerFluid2,a.DisplacementComponent3,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.c,a.Material2,a.Mode3,a.SLamb2,a.SShear2,str2double(a.LayerCountUI2.String),a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.Delta,a.SamplesPerCycle3,a.Cycles3,a.Distance3,a.FrequencyResolution2,a.Frequency3,a.MultiMode3,a.SpectrumThreshold3,a.TimeLimitFactor3,a.ALambModes2,a.AShearModes2,a.BLambModes2,a.BShearModes2,a.SLambModes2,a.SShearModes2,a.Window3,a.OutputWindowUI3);
        a.YAxis3 = SignalMax(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    end
    a.GateUI3.String = ['[',num2str(a.Gate3(1),'%.1f'),' ',num2str(a.Gate3(2),'%.1f'),']'];
    a.XAxisUI3.String = ['[',num2str(a.XAxis3(1)),' ',num2str(a.XAxis3(2)),']'];
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
end
function HelperFunction3
    Colorization(a,1)
    HelperFunction4
end
function HelperFunction4
    if  a.DataType3 == 1
        a.YAxis3 = SignalMax(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    elseif a.DataType3 == 2
        a.YAxis3 = SignalMax(a.uSum3,a.uSumALamb3,a.uSumSLamb3,a.uSumBLamb3,a.uSumAShear3,a.uSumSShear3,a.uSumBShear3,a.Plot_ALamb3,a.Plot_SLamb3,a.Plot_AShear3,a.Plot_SShear3,a.Amplitude_ALamb3,a.Amplitude_SLamb3,a.Amplitude_AShear3,a.Amplitude_SShear3,a.MultiMode3);
    end
    a.YAxisUI3.String = a.YAxis3;
    [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
end
end