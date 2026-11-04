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
function a = CallbackModule_Anisotropic(source,a,Tab3)
%#ok<*AGROW>
%#ok<*GVMIS>
%#ok<*FXUP>
global Stop
Stop = 0;
source.String = replace(source.String,',','.');
if  strcmp(source.Tag,'63') % Fix
    a.Fix2 = source.Value;
    if  ~a.Fix2
        a.Search2 = 0;
        AdjustRange
    end
elseif strcmp(source.Tag,'1') % Propagation angle
    a.PropagationAngle = str2double(source.String);
    CheckDecoupling
    a.EffectiveLayupString2 = ['[',replace(num2str(a.LayerOrientations-a.PropagationAngle),whitespacePattern,'/'),']'];
    if  a.SuperLayers > 1
        a.EffectiveLayupString2 = append(a.EffectiveLayupString2,num2str(a.SuperLayers));
    end
    if  a.SymmetricSystem
        a.EffectiveLayupString2 = append(a.EffectiveLayupString2,'s');
    end
    a.EffectiveLayupUI2.String = a.EffectiveLayupString2;
    a.Search2 = 0;
elseif strcmp(source.Tag,'3') % Frequency limit
    a.FrequencyLimit2 = str2double(source.String);
    a.Frequency12 = a.FrequencyLimit2;
    a.Frequency1UI2.String = a.FrequencyLimit2;
    a.Frequency22 = a.FrequencyLimit2;
    a.Frequency2UI2.String = a.FrequencyLimit2;
    a.Search2 = 0;
    AdjustAxes
elseif strcmp(source.Tag,'4') % Frequency step
    a.FrequencyResolution2 = str2double(source.String);
elseif strcmp(source.Tag,'2') % Phase velocity limit
    a.PhaseVelocityLimit2 = str2double(source.String)*1e3;
    if  a.Quantity12 == 1
        a.YAxisUI2.String = ['[0 ',source.String,']'];
        a.YAxis2 = eval(a.YAxisUI2.String);
    end
    a.Search2 = 0;
elseif strcmp(source.Tag,'82') % Phase velocity accuracy
    a.Accuracy2 = str2double(source.String);
elseif strcmp(source.Tag,'84') % Matrix methods
    a.MatrixMethods2 = source.Value;
elseif strcmp(source.Tag,'76') % Force 2-D tracing
    a.Force2DTracing2 = source.Value;
    if  a.Force2DTracing2
        a.PhaseVelocityLimit2UI2.Enable = 'on';
        a.AttenuationLimitUI2.Enable = 'on';
        a.SweepsUI2.Enable = 'on';
        a.SweepsSectionsUI2.Enable = 'on';
        a.HigherOrderModesUI2.Enable = 'off';
    else
        a.PhaseVelocityLimit2UI2.Enable = 'off';
        a.AttenuationLimitUI2.Enable = 'off';
        a.SweepsUI2.Enable = 'off';
        a.SweepsSectionsUI2.Enable = 'off';
        a.HigherOrderModesUI2.Enable = 'on';
    end
elseif strcmp(source.Tag,'79') % Phase velocity limit 2
    a.PhaseVelocityLimit22 = str2double(source.String)*1e3;
elseif strcmp(source.Tag,'86') % Attenuation limit
    a.AttenuationLimit2 = str2double(source.String);
elseif strcmp(source.Tag,'87') % Sweeps
    a.Sweeps2 = str2double(source.String(source.Value));
elseif strcmp(source.Tag,'88') % Sweep sections
    a.SweepSections2 = str2double(source.String(source.Value));
elseif strcmp(source.Tag,'7') % Higher order modes
    a.HigherOrderModes2 = source.Value;
    a.Search2 = 0;
elseif strcmp(source.Tag,'9') % Symmetric modes/Torsional modes
    if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Circumferential')
        a.SymmetricModes2 = source.Value;
    else
        a.TorsionalModes2 = source.Value;
    end
elseif strcmp(source.Tag,'10') % Antisymmetric modes/Longitudinal modes
    if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Circumferential')
        a.AntisymmetricModes2 = source.Value;
    else
        a.LongitudinalModes2 = source.Value;
    end
elseif strcmp(source.Tag,'13') % Lamb modes/Flexural modes
    if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Circumferential')
        a.LambModes2 = source.Value;
    else
        a.FlexuralModes2 = source.Value;
        if  a.FlexuralModes2
            a.ShearHorizontalModes_FlexuralModeOrdersUI2.Enable = 'on';
        else
            a.ShearHorizontalModes_FlexuralModeOrdersUI2.Enable = 'off';
        end
    end
elseif strcmp(source.Tag,'14') % Shear horizontal modes/Flexural mode orders
    if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Circumferential')
        a.ShearHorizontalModes2 = source.Value;
    else
        if  rem(str2double(source.String),1) ~= 0 || str2double(source.String) < 1
            source.String = ceil(abs(str2double(source.String)));
            errordlg('The number of circumferential orders of flexural modes must be a positive integer.','Error');
        end
        a.FlexuralModeOrders2 = str2double(source.String);
        a.Search2 = 0;
    end
elseif strcmp(source.Tag,'15') % Step
    a.Step2 = str2double(source.String);
    a.Search2 = 0;
elseif strcmp(source.Tag,'16') % Search % frequency sweeps/phase velocity sweeps at high frequency (like in the "Trace modes" callback, but without mode selection restrictions: no if HigherOrderModes, LambModes, etc.)
    CheckDecoupling
    TransformStiffness
    FrequencyRange = 0:a.FrequencyResolution2:a.FrequencyLimit2;
    FrequencyRange(1) = a.FrequencyRangeStart2;
    if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution2
        FrequencyRange(1) = .1*a.FrequencyResolution2;
    end
    if  strcmp(a.Geometry2,'Plate')
        [a.MatrixMethodLimit2,a.EvanescenceLimit2,a.XS0,a.XA0,a.XSH0] = MatrixMethodLimitComputer_Plate(a.c,a.DeltaReal,a.Material2,a.Hybrid,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.Symmetric2,a.SymmetricSystem,a.Decoupled,a.DC,FrequencyRange);
    elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
        [a.MatrixMethodLimit2,a.EvanescenceLimit2,a.XL01,a.XF11,a.XT01] = MatrixMethodLimitComputer_Rod_Pipe(a.Geometry2,a.Material2,a.Layers2,a.Radius/1e3,FrequencyRange,a.LowerFluid2,a.ToggleLowerFluid2);
    elseif strcmp(a.Geometry2,'Circumferential')
        [a.MatrixMethodLimit2,a.EvanescenceLimit2,a.XC0,a.XCSH0] = MatrixMethodLimitComputer_Circumferential(a.Material2,a.Layers2,a.Radius/1e3,FrequencyRange);
    end
    FrequencySweeps
elseif strcmp(source.Tag,'17') % Trace modes/Stop trace
    if  source.Value
        CheckDecoupling
        TransformStiffness
        FrequencyRange = 0:a.FrequencyResolution2:a.FrequencyLimit2;
        FrequencyRangeF = FrequencyRange;
        FrequencyRange(1) = a.FrequencyRangeStart2;
        if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution2
            FrequencyRange(1) = .1*a.FrequencyResolution2;
        end
        if  strcmp(a.Geometry2,'Pipe') && a.ToggleLowerFluid2 && ~a.Sink2
            FrequencyRangeF(1) = .5*a.FrequencyResolution2;
        else
            FrequencyRangeF(1) = .1*a.FrequencyResolution2;
        end
        if  strcmp(a.Geometry2,'Plate')
            [a.MatrixMethodLimit2,a.EvanescenceLimit2,a.XS0,a.XA0,a.XSH0] = MatrixMethodLimitComputer_Plate(a.c,a.DeltaReal,a.Material2,a.Hybrid,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.Symmetric2,a.SymmetricSystem,a.Decoupled,a.DC,FrequencyRange);
        elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
            [a.MatrixMethodLimit2,a.EvanescenceLimit2,a.XL01,a.XF11,a.XT01] = MatrixMethodLimitComputer_Rod_Pipe(a.Geometry2,a.Material2,a.Layers2,a.Radius/1e3,FrequencyRange,a.LowerFluid2,a.ToggleLowerFluid2);
        elseif strcmp(a.Geometry2,'Circumferential')
            [a.MatrixMethodLimit2,a.EvanescenceLimit2,a.XC0,a.XCSH0] = MatrixMethodLimitComputer_Circumferential(a.Material2,a.Layers2,a.Radius/1e3,FrequencyRange);
        end
        if  ~a.Search2
            FrequencySweeps
        end
        if  a.Multithreading
            imshow('Multithreading22.png','Parent',a.a1UI2)
        end
        source.String = 'Stop trace';
        source.Tooltip = 'Stop tracing dispersion curves.';
        if  strcmp(a.Geometry2,'Plate')
            [a.SLamb2,a.ALamb2,a.BLamb2,a.SShear2,a.AShear2,a.BShear2] = Computer_Anisotropic_Plate(a.Multithreading,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Hybrid,a.FrequencyLimit2,a.Material2,a.PropagationAngle,a.PhaseVelocityStep2,a.FrequencyOffset2,a.LambPhaseVelocitySweepRange12,a.LambPhaseVelocitySweepRange22,a.Force2DTracing2,a.AttenuationLimit2,a.SearchWidth,a.SearchAreaSections,a.SearchAreaExtensions,a.Sweeps2,a.SweepSections2,a.LayerThicknesses/1e3,FrequencyRange,FrequencyRangeF,a.HSLamb2,a.HSShear2,a.HALamb2,a.HAShear2,a.HBLamb2,a.HBShear2,a.c,a.Delta,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.PlateThickness/1e3,a.SuperLayerSize,a.Decoupled,a.DC,a.Pattern,a.SymmetricSystem,a.Symmetric2,a.FrequencyResolution2,a.PhaseVelocityLimit2,a.PhaseVelocityLimit22,a.Accuracy2,a.PhaseVelocitySections2,a.FrequencySections2,a.HigherOrderModes2,a.SymmetricModes2,a.AntisymmetricModes2,a.LambModes2,a.ShearHorizontalModes2,a.MissingSamples2,a.BelowCutoffWidth2,a.CriticalSlope,a.GearSlopeRanges,a.GearFrequencyResolution,a.MatrixMethods2,a.MatrixMethodLimit2,a.XS0,a.XA0,a.XSH0);
            a.F2={[]};a.L2={[]};a.T2={[]};a.CLamb2={[]};a.CShear2={[]};
        elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
            [a.F2,a.L2,a.T2] = Computer_Isotropic_Rod_Pipe(a.Geometry2,a.Multithreading,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,a.FrequencyLimit2,a.Material2,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Symmetric2,a.PhaseVelocityStep2,a.FrequencyOffset2,a.LambPhaseVelocitySweepRange12,a.LambPhaseVelocitySweepRange22,a.Force2DTracing2,a.AttenuationLimit2,a.SearchWidth,a.SearchAreaSections,a.SearchAreaExtensions,a.Sweeps2,a.SweepSections2,FrequencyRange,FrequencyRangeF,a.HL2,a.HF2,a.HT2,a.FrequencyResolution2,a.PhaseVelocityLimit2,a.PhaseVelocityLimit22,a.Accuracy2,a.PhaseVelocitySections2,a.FrequencySections2,a.HigherOrderModes2,a.LongitudinalModes2,a.FlexuralModes2,a.TorsionalModes2,a.MissingSamples2,a.BelowCutoffWidth2,a.FlexuralModeOrders2,a.CriticalSlope,a.GearSlopeRanges,a.GearFrequencyResolution,a.LineColors2,a.IsotropicMultilayerCylinder,a.Layers2,a.Radius/1e3,a.XL01,a.XF11,a.XT01,a.MatrixMethods2,a.MatrixMethodLimit2);
            a.SLamb2={[]};a.ALamb2={[]};a.BLamb2={[]};a.SShear2={[]};a.AShear2={[]};a.BShear2={[]};a.CLamb2={[]};a.CShear2={[]};
        elseif strcmp(a.Geometry2,'Circumferential')
            [a.CLamb2,a.CShear2] = Computer_Isotropic_Circumferential(a.Multithreading,a.FrequencyLimit2,a.Material2,a.Radius(end)/1e3,a.InnerRadius/1e3,a.PhaseVelocityStep2,a.FrequencyOffset2,FrequencyRange,a.HCLamb2,a.HCShear2,a.FrequencyResolution2,a.PhaseVelocityLimit2,a.Accuracy2,a.PhaseVelocitySections2,a.FrequencySections2,a.HigherOrderModes2,a.LambModes2,a.ShearHorizontalModes2,a.MissingSamples2,a.BelowCutoffWidth2,a.IsotropicMultilayerCylinder,a.Layers2,a.Radius/1e3,a.XC0,a.XCSH0,a.MatrixMethods2,a.MatrixMethodLimit2);
            a.SLamb2={[]};a.ALamb2={[]};a.BLamb2={[]};a.SShear2={[]};a.AShear2={[]};a.BShear2={[]};a.F2={[]};a.L2={[]};a.T2={[]};
        end
        source.Value = 0;
        source.String = 'Trace modes';
        source.Tooltip = 'Start tracing the dispersion curves.';
        if  a.Multithreading
            imshow('Multithreading21.png','Parent',a.a1UI2)
        end
        if  Stop
            return
        end
        % if  (strcmp(a.Geometry2,'Plate') && ~a.Viscoelastic2 && ((a.LambModes2 && ~a.FluidLoading2) || ~a.LambModes2)) ||...
        %     (strcmp(a.Geometry2,'Rod') && ~a.Viscoelastic2 && ~a.FlexuralModes2 && ((a.LongitudinalModes2 && ~a.FluidLoading2) || ~a.LongitudinalModes2)) ||...
        %     (strcmp(a.Geometry2,'Pipe') && ~a.Viscoelastic2 && ~a.FlexuralModes2 && ~a.LongitudinalModes2)
        %     a.SamplesX3ceUI2.Enable = 'off';
        %     a.CalculateCeUI2.Enable = 'off';
        % else
        %     a.SamplesX3ceUI2.Enable = 'on';
        %     a.CalculateCeUI2.Enable = 'on';
        % end
        a.ModeNames2 = {''};
        if  strcmp(a.Geometry2,'Plate')
            if  ~isempty(a.ALamb2{1})
                for i = 0:length(a.ALamb2)-1
                    eval(sprintf('a.ModeNames2(i+1) = {''A%u''};',i));
                end
            end
            if  ~isempty(a.SLamb2{1})
                for i = 0:length(a.SLamb2)-1
                    eval(sprintf('a.ModeNames2(end+1) = {''S%u''};',i));
                end
            end
            if  ~isempty(a.BLamb2{1})
                for i = 0:length(a.BLamb2)-1
                    eval(sprintf('a.ModeNames2(i+1) = {''B%u''};',i));
                end
            end
            if  ~isempty(a.AShear2{1})
                for i = 1:length(a.AShear2)
                    eval(sprintf('a.ModeNames2(end+1) = {''ASH%u''};',i));
                end
            end
            if  ~isempty(a.SShear2{1})
                for i = 0:length(a.SShear2)-1
                    eval(sprintf('a.ModeNames2(end+1) = {''SSH%u''};',i));
                end
            end
            if  ~isempty(a.BShear2{1})
                for i = 0:length(a.BShear2)-1
                    eval(sprintf('a.ModeNames2(end+1) = {''BSH%u''};',i));
                end
            end
        elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
            if  ~isempty(a.T2{1})
                for i = 1:length(a.T2)
                    eval(sprintf('a.ModeNames2(i) = {''T(0,%u)''};',i));
                end
            end
            if  ~isempty(a.L2{1})
                for i = 1:length(a.L2)
                    eval(sprintf('a.ModeNames2(end+1) = {''L(0,%u)''};',i));
                end
            end
            if  ~isempty(a.F2{1})
                for n = 1:length(a.F2)
                    for i = 1:length(a.F2{n})
                        eval(sprintf('a.ModeNames2(end+1) = {''F(%u,%u)''};',n,i));
                    end
                end
            end
        elseif strcmp(a.Geometry2,'Circumferential')
            if  ~isempty(a.CLamb2{1})
                for i = 0:length(a.CLamb2)-1
                    eval(sprintf('a.ModeNames2(i+1) = {''C%u''};',i));
                end
            end
            if  ~isempty(a.CShear2{1})
                for i = 0:length(a.CShear2)-1
                    eval(sprintf('a.ModeNames2(end+1) = {''CSH%u''};',i));
                end
            end
        end
        a.ModeNames2(cellfun(@isempty,a.ModeNames2)) = [];
        a.Mode1UI2.Value = 1;
        a.Mode1UI2.String = a.ModeNames2;
        a.Mode12 = a.ModeNames2{1};
        a.Mode2UI2.Value = 1;
        a.Mode2UI2.String = a.ModeNames2;
        a.Mode22 = a.ModeNames2{1};
        a.ModeShapeSettingChanged2 = 1;

        if  strcmp(a.Geometry2,'Plate') % connection to Signal simulator
            a.Title3 = 2;
            a.TitleUI3.Style = 'popupmenu';
            a.TitleUI3.Value = 1;
            a.TitleUI3.String = {'with layup','without layup','no title'};
            a.TitleUI3.Tooltip = 'Choose to plot the title with the layup  included, without layup, or no title at all.';
            if  ~a.Hybrid
                a.DataUI3.String = a.Material2{1}.Name;
            elseif a.Hybrid
                a.DataUI3.String = 'Hybrid';
            end
            a.WaveGeneratorUI3.Enable = 'on';
            a.PlotUI3.Enable = 'off';
            a.DataType3 = 2;
            a.Mode3 = '';
            a.Frequency3 = a.FrequencyLimit2/10;
            if  a.Chirp3
                a.Frequency3(2) = 2*a.Frequency3(1);
            end
            if  a.Chirp3
                String = [num2str(a.Frequency3(1)),' - ',num2str(a.Frequency3(2)),' kHz, '];
                if  a.ChirpType3 == 1
                    String = append(String,'linear chirp');
                elseif a.ChirpType3 == 2
                    String = append(String,'quadratic chirp');
                elseif a.ChirpType3 == 3
                    String = append(String,'exponential chirp');
                end
                a.FrequencyDisplayUI3.String = String;
            else
                a.FrequencyDisplayUI3.String = [num2str(a.Frequency3),' kHz'];
            end
            [a.ALambModes2,a.SLambModes2,a.BLambModes2,a.AShearModes2,a.SShearModes2,a.BShearModes2] = ModeFinder(a.ALamb2,a.SLamb2,a.BLamb2,a.AShear2,a.SShear2,a.BShear2,mean(a.Frequency3));
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
                if  a.Symmetric2
                    eval(sprintf('a.ALamb%uaUI3.String = [''A'',num2str(i)];',i));
                    eval(sprintf('a.SLamb%uaUI3.String = [''S'',num2str(i)];',i));
                    eval(sprintf('a.SLamb%uaUI3.Position(3) = 16;',i));
                    eval(sprintf('a.SLamb%ubUI3.Position(1) = 88;',i));
                    eval(sprintf('a.SLamb%ucUI3.Position(1) = 103;',i));
                else
                    if  ~a.Decoupled
                        eval(sprintf('a.ALamb%uaUI3.String = [''B'',num2str(i)];',i));
                    else
                        eval(sprintf('a.ALamb%uaUI3.String = [''B'',num2str(i)];',i));
                        eval(sprintf('a.SLamb%uaUI3.String = [''BSH'',num2str(i)];',i));
                        eval(sprintf('a.SLamb%uaUI3.Position(3) = 30;',i));
                        eval(sprintf('a.SLamb%ubUI3.Position(1) = 102;',i));
                        eval(sprintf('a.SLamb%ucUI3.Position(1) = 117;',i));
                    end
                end
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
        else
            a.DataType3 = 2;
            a.DataUI3.String = '--';
            a.CalculateUI3.Enable = 'off';
            a.ALambModes2 = false;
            a.SLambModes2 = false;
            a.BLambModes2 = false;
            a.AShearModes2 = false;
            a.SShearModes2 = false;
            a.BShearModes2 = false;
            ShowModes_Signal(a)
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
        end
    
        a.CalculateCeUI2.Value = 1;
        a.CalculateCeUI2.String = 'Stop calculate';
        a.CalculateCeUI2.Tooltip = 'Stop calculating the energy velocity.';
        if  strcmp(a.Geometry2,'Plate')
            [a.SLamb2,a.ALamb2,a.BLamb2,a.SShear2,a.AShear2,a.BShear2] = Computer_Anisotropic_Plate_EnergyVelocity(a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.SLamb2,a.ALamb2,a.BLamb2,a.SShear2,a.AShear2,a.BShear2,a.c,a.Material2,a.Layers2,a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Delta,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SamplesX3ce2);
        elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
            [a.F2,a.L2,a.T2] = Computer_Isotropic_Rod_Pipe_EnergyVelocity(a.Geometry2,a.F2,a.L2,a.T2,a.Material2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,a.Layers2,a.Radius/1e3,a.SamplesX3ce2);
        elseif strcmp(a.Geometry2,'Circumferential')
            [a.CLamb2,a.CShear2] = Computer_Isotropic_Circumferential_EnergyVelocity_Multi(a.CLamb2,a.CShear2,a.Material2,a.Layers2,a.Radius/1e3,a.SamplesX3ce2);
        end
        a.CalculateCeUI2.Value = 0;
        a.CalculateCeUI2.String = 'Calculate ce';
        a.CalculateCeUI2.Tooltip = 'Start calculating the energy velocity.';
    else
        Stop = 1;
        if  a.Multithreading
            source.Value = 0;
            source.String = 'Trace modes';
            source.Tooltip = 'Start tracing the dispersion curves.';
            imshow('Multithreading21.png','Parent',a.a1UI2)
            cancelAll(a.Pool.FevalQueue)
        end
    end
elseif strcmp(source.Tag,'83') % Samples x3 (ce)
    a.SamplesX3ce2 = str2double(source.String);
elseif strcmp(source.Tag,'18') % Calculate ce/Stop calculate
    if  source.Value && strcmp(source.String,'Calculate ce')
        source.String = 'Stop calculate';
        source.Tooltip = 'Stop calculating the energy velocity.';
        if  strcmp(a.Geometry2,'Plate')
            [a.SLamb2,a.ALamb2,a.BLamb2,a.SShear2,a.AShear2,a.BShear2] = Computer_Anisotropic_Plate_EnergyVelocity(a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.SLamb2,a.ALamb2,a.BLamb2,a.SShear2,a.AShear2,a.BShear2,a.c,a.Material2,a.Layers2,a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Delta,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SamplesX3ce2);
        elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
            [a.F2,a.L2,a.T2] = Computer_Isotropic_Rod_Pipe_EnergyVelocity(a.Geometry2,a.F2,a.L2,a.T2,a.Material2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,a.Layers2,a.Radius/1e3,a.SamplesX3ce2);
        elseif strcmp(a.Geometry2,'Circumferential')
            [a.CLamb2,a.CShear2] = Computer_Isotropic_Circumferential_EnergyVelocity_Multi(a.CLamb2,a.CShear2,a.Material2,a.Layers2,a.Radius/1e3,a.SamplesX3ce2);
        end
        source.Value = 0;
        source.String = 'Calculate ce';
        source.Tooltip = 'Start calculating the energy velocity.';
    else
        Stop = 1;
    end
elseif strcmp(source.Tag,'19') % Quantity (Dispersion diagrams)
    a.Quantity12 = source.Value;
    if  a.Quantity12 < 3
        if  a.SuperLayerSize == 1
            a.Option1UI2.Enable = 'on';
        else
            a.Option1UI2.Enable = 'off';
        end
        a.Option1UI2.Style = 'checkbox';
        a.Option1UI2.String = ' ';
        a.Option1UI2.Value = a.BulkVelocities2;
        a.Option1UI2.Tooltip = 'Check this to show the bulk wave velocities.';
        a.Option1TextUI2.String = 'Bulk velocities';
    elseif a.Quantity12 == 3
        a.Option1UI2.Enable = 'on';
        a.Option1UI2.Style = 'edit';
        a.Option1UI2.String = a.Distance2;
        a.Option1UI2.Tooltip = 'Enter the propagation distance from the source to the sensor for which the propagation time shall be calculated.';
        a.Option1TextUI2.String = 'Distance (mm)';
    elseif a.Quantity12 == 4
        a.Option1UI2.Enable = 'on';
        a.Option1UI2.Style = 'popupmenu';
        a.Option1UI2.String = a.FieldNames_Fluid;
        a.Option1UI2.Value = find(strcmp(a.Couplant2.Name,a.Option1UI2.String));
        a.Option1UI2.Tooltip = 'Select for which surrounding medium the coincidence angle shall be displayed.';
        a.Option1TextUI2.String = 'Couplant';
    elseif a.Quantity12 > 4
        a.Option1UI2.Enable = 'off';
    end
    if  a.Quantity12 == 3
        a.XAxisModeTextUI2.String = 'Y-axis mode';
        a.XAxisModeUI2.Tooltip = 'Select the frequency''s dimension on the Y-axis.';
        a.XAxisTextUI2.String = ['X-axis (',char(181),'s)'];
        a.XAxisUI2.Tooltip = 'Enter which propagation time range shall be plotted.';
        if  a.XAxisMode2 == 1
            a.YAxisTextUI2.String = 'Y-axis (kHz)';
            a.YAxisUI2.Tooltip = 'Enter which frequency range shall be plotted.';
        elseif a.XAxisMode2 == 2
            a.YAxisTextUI2.String = 'Y-axis (MHz)';
            a.YAxisUI2.Tooltip = 'Enter which frequency range shall be plotted.';
        elseif a.XAxisMode2 == 3
            a.YAxisTextUI2.String = ['Y-axis (MHz',char(8901),'mm)'];
            a.YAxisUI2.Tooltip = 'Enter which frequency-thickness range shall be plotted.';
        end
    else
        a.XAxisModeTextUI2.String = 'X-axis mode';
        a.XAxisModeUI2.Tooltip = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode2 == 1
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisUI2.Tooltip = 'Enter which frequency range shall be plotted.';
        elseif a.XAxisMode2 == 2
            a.XAxisTextUI2.String = 'X-axis (MHz)';
            a.XAxisUI2.Tooltip = 'Enter which frequency range shall be plotted.';
        elseif a.XAxisMode2 == 3
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI2.Tooltip = 'Enter which frequency-thickness range shall be plotted.';
        end
    end
    if  a.Quantity12 == 1
        a.YAxisTextUI2.String = 'Y-axis (m/ms)';
        a.YAxisUI2.Tooltip = 'Enter which phase velocity range shall be plotted.';
    elseif a.Quantity12 == 2
        a.YAxisTextUI2.String = 'Y-axis (m/ms)';
        a.YAxisUI2.Tooltip = 'Enter which energy velocity range shall be plotted.';
    elseif a.Quantity12 == 4
        a.YAxisTextUI2.String = ['Y-axis (',char(176),')'];
        a.YAxisUI2.Tooltip = 'Enter which coincidence angle range shall be plotted.';
    elseif a.Quantity12 == 5
        if  a.XAxisMode2 == 3
            a.YAxisTextUI2.String = 'Y-axis (mm/mm)';
            a.YAxisUI2.Tooltip = 'Enter which wavelength per thickness range shall be plotted.';
        else
            a.YAxisTextUI2.String = 'Y-axis (mm)';
            a.YAxisUI2.Tooltip = 'Enter which wavelength range shall be plotted.';
        end
    elseif a.Quantity12 == 6
        if  a.XAxisMode2 == 3
            a.YAxisTextUI2.String = ['Y-axis (rad/mm',char(8901),'mm)'];
            a.YAxisUI2.Tooltip = 'Enter which wavenumber-thickness range shall be plotted.';
        else
            a.YAxisTextUI2.String = 'Y-axis (rad/mm)';
            a.YAxisUI2.Tooltip = 'Enter which wavenumber range shall be plotted.';
        end 
    elseif a.Quantity12 == 7
        if  a.XAxisMode2 == 3
            a.YAxisTextUI2.String = ['Y-axis (Np/m',char(8901),'mm)'];
            a.YAxisUI2.Tooltip = 'Enter which attenuation-thickness range shall be plotted.';
        else
            a.YAxisTextUI2.String = 'Y-axis (Np/m)';
            a.YAxisUI2.Tooltip = 'Enter which attenuation range shall be plotted.';
        end
    end
    AdjustAxes
elseif strcmp(source.Tag,'20') % option (Dispersion diagrams)
    if  a.Quantity12 < 3
        a.BulkVelocities2 = source.Value;
    elseif a.Quantity12 == 3
        a.Distance2 = str2double(source.String);
    elseif a.Quantity12 == 4
        a.Couplant2 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{source.Value});
    end
    if  a.Quantity12 == 3 || a.Quantity12 == 4
        AdjustAxes
    end
elseif strcmp(source.Tag,'21') % X-axis mode
    a.XAxisMode2 = source.Value;
    if  a.XAxisMode2 == 3
        if  a.Quantity12 == 3
            a.YAxisTextUI2.String = ['Y-axis (MHz',char(8901),'mm)'];
            a.YAxisUI2.Tooltip = 'Enter which frequency-thickness range shall be plotted.';
        else
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI2.Tooltip = 'Enter which frequency-thickness range shall be plotted.';
        end
        if  a.Quantity12 == 5
            a.YAxisTextUI2.String = 'Y-axis (mm/mm)';
            a.YAxisUI2.Tooltip = 'Enter which wavelength per thickness range shall be plotted.';
        elseif a.Quantity12 == 6
            a.YAxisTextUI2.String = ['Y-axis (rad/mm',char(8901),'mm)'];
            a.YAxisUI2.Tooltip = 'Enter which wavenumber-thickness range shall be plotted.';
        elseif a.Quantity12 == 7
            a.YAxisTextUI2.String = ['Y-axis (Np/m',char(8901),'mm)'];
            a.YAxisUI2.Tooltip = 'Enter which attenuation-thickness range shall be plotted.';
        end
    else
        if  a.Quantity12 == 3
            if  a.XAxisMode2 == 1
                a.YAxisTextUI2.String = 'Y-axis (kHz)';
            elseif a.XAxisMode2 == 2
                a.YAxisTextUI2.String = 'Y-axis (MHz)';
            end
            a.YAxisUI2.Tooltip = 'Enter which frequency range shall be plotted.';
        else
            if  a.XAxisMode2 == 1
                a.XAxisTextUI2.String = 'X-axis (kHz)';
            elseif a.XAxisMode2 == 2
                a.XAxisTextUI2.String = 'X-axis (MHz)';
            end
            a.XAxisUI2.Tooltip = 'Enter which frequency range shall be plotted.';
        end
        if  a.Quantity12 == 5
            a.YAxisTextUI2.String = 'Y-axis (mm)';
            a.YAxisUI2.Tooltip = 'Enter which wavelength range shall be plotted.';
        elseif a.Quantity12 == 6
            a.YAxisTextUI2.String = 'Y-axis (rad/mm)';
            a.YAxisUI2.Tooltip = 'Enter which wavenumber range shall be plotted.';
        elseif a.Quantity12 == 7
            a.YAxisTextUI2.String = 'Y-axis (Np/m)';
            a.YAxisUI2.Tooltip = 'Enter which attenuation range shall be plotted.';
        end
    end
    AdjustAxes
elseif strcmp(source.Tag,'22') % X-axis
    a.XAxis2 = eval(source.String);
    if  a.Quantity12 ~= 3
        if  a.XAxisMode2 == 1
            a.XAxis2 = eval(source.String);
        elseif a.XAxisMode2 == 2
            a.XAxis2 = eval(source.String)*1e3;
        elseif a.XAxisMode2 == 3
            if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Rod')
                a.XAxis2 = eval(source.String)*1e3/a.PlateThickness;
            elseif strcmp(a.Geometry2,'Pipe') || strcmp(a.Geometry2,'Circumferential')
                a.XAxis2 = eval(source.String)*1e3/a.LayupThickness;
            end
        end
    end
elseif strcmp(source.Tag,'23') % Y-axis
    a.YAxis2 = eval(source.String);
    if  a.Quantity12 == 3
        if  a.XAxisMode2 == 1
            a.YAxis2 = eval(source.String);
        elseif a.XAxisMode2 == 2
            a.YAxis2 = eval(source.String)*1e3;
        elseif a.XAxisMode2 == 3
            if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Rod')
                a.YAxis2 = eval(source.String)*1e3/a.PlateThickness;
            elseif strcmp(a.Geometry2,'Pipe') || strcmp(a.Geometry2,'Circumferential')
                a.YAxis2 = eval(source.String)*1e3/a.LayupThickness;
            end
        end
    end
elseif strcmp(source.Tag,'24') % Plot (Dispersion diagrams)
    try
        if  a.Quantity12 == 1
            DispersionDiagram_PhaseVelocity(a.Geometry2,a.Hybrid,a.LayerOrientations-a.PropagationAngle,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.BulkVelocities2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BoxLineWidth2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.FileName2,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.PNG2,a.PropagationAngle,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled) 
        elseif a.Quantity12 == 2
            if  ~a.Decoupled
                DispersionDiagram_EnergyVelocity2(2,a.Hybrid,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.BLamb2,a.BoxLineWidth2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.FileName2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.SLamb2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2)
                DispersionDiagram_EnergyVelocity2(1,a.Hybrid,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.BLamb2,a.BoxLineWidth2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.FileName2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.SLamb2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2)
                DispersionDiagram_EnergyVelocity(2,a.Geometry2,a.Hybrid,a.LayerOrientations-a.PropagationAngle,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.BulkVelocities2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BoxLineWidth2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.FileName2,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.PNG2,a.PropagationAngle,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
            end
            DispersionDiagram_EnergyVelocity(1,a.Geometry2,a.Hybrid,a.LayerOrientations-a.PropagationAngle,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.BulkVelocities2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BoxLineWidth2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.FileName2,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.PNG2,a.PropagationAngle,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
        elseif a.Quantity12 == 3
            DispersionDiagram_PropagationTime(a.Geometry2,a.Hybrid,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Distance2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BoxLineWidth2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.FileName2,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.PNG2,a.PropagationAngle,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
        elseif a.Quantity12 == 4
            DispersionDiagram_CoincidenceAngle(a.Geometry2,a.Hybrid,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Couplant2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BoxLineWidth2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.FileName2,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.PNG2,a.PropagationAngle,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
        elseif a.Quantity12 == 5
            DispersionDiagram_Wavelength(a.Geometry2,a.Hybrid,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BoxLineWidth2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.FileName2,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.PNG2,a.PropagationAngle,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
        elseif a.Quantity12 == 6
            DispersionDiagram_Wavenumber(a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BoxLineWidth2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.FileName2,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.PNG2,a.PropagationAngle,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
        elseif a.Quantity12 == 7
            DispersionDiagram_Attenuation(a.Geometry2,a.Hybrid,a.LayupStringTex2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BoxLineWidth2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.Title2,a.LambModes2,a.LineWidth2,a.Material2,a.PDF2,a.FileName2,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.PNG2,a.PropagationAngle,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric2,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
        end
    catch ME
        switch ME.identifier
        case 'MATLAB:nonExistentField'
            errordlg('You must trace dispersion curves before you can plot any data. Press ''Traces modes''.','Trace dispersion curves')
        otherwise
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Error')
        end
    end
elseif strcmp(source.Tag,'25') % Quantity (Through thickness profiles)
    a.Quantity22 = source.Value;
    if  a.Quantity22 == 1
        a.PhaseUI2.Enable = 'on';
        a.x11UI2.Enable = 'off';
        a.x12UI2.Enable = 'off';
        a.x13UI2.Enable = 'on';
        a.x22UI2.Enable = 'off';
        a.x23UI2.Enable = 'on';
        a.x33UI2.Enable = 'on';
        a.x13UI2.String = '1';
        a.x23UI2.String = '2';
        a.x33UI2.String = '3';
        a.x13UI2.Tooltip = 'Displacement u1';
        a.x23UI2.Tooltip = 'Displacement u2';
        a.x33UI2.Tooltip = 'Displacement u3';
    elseif a.Quantity22 == 2
        a.PhaseUI2.Enable = 'on';
        a.x11UI2.Enable = 'on';
        a.x12UI2.Enable = 'on';
        a.x13UI2.Enable = 'on';
        a.x22UI2.Enable = 'on';
        a.x23UI2.Enable = 'on';
        a.x33UI2.Enable = 'on';
        a.x13UI2.String = '13';
        a.x23UI2.String = '23';
        a.x33UI2.String = '33';
        a.x11UI2.Tooltip = ['Stress ',char(963),'11'];
        a.x12UI2.Tooltip = ['Stress ',char(963),'12'];
        a.x13UI2.Tooltip = ['Stress ',char(963),'13'];
        a.x22UI2.Tooltip = ['Stress ',char(963),'22'];
        a.x23UI2.Tooltip = ['Stress ',char(963),'23'];
        a.x33UI2.Tooltip = ['Stress ',char(963),'33'];
    elseif a.Quantity22 == 3
        a.PhaseUI2.Enable = 'on';
        a.x11UI2.Enable = 'on';
        a.x12UI2.Enable = 'on';
        a.x13UI2.Enable = 'on';
        a.x22UI2.Enable = 'on';
        a.x23UI2.Enable = 'on';
        a.x33UI2.Enable = 'on';
        a.x13UI2.String = '13';
        a.x23UI2.String = '23';
        a.x33UI2.String = '33';
        a.x11UI2.Tooltip = ['Strain ',char(949),'11'];
        a.x12UI2.Tooltip = ['Strain ',char(949),'12'];
        a.x13UI2.Tooltip = ['Strain ',char(949),'13'];
        a.x22UI2.Tooltip = ['Strain ',char(949),'22'];
        a.x23UI2.Tooltip = ['Strain ',char(949),'23'];
        a.x33UI2.Tooltip = ['Strain ',char(949),'33'];
    elseif a.Quantity22 == 4
        a.PhaseUI2.Enable = 'off';
        a.x11UI2.Enable = 'off';
        a.x12UI2.Enable = 'off';
        a.x13UI2.Enable = 'on';
        a.x22UI2.Enable = 'off';
        a.x23UI2.Enable = 'on';
        a.x33UI2.Enable = 'on';
        a.x13UI2.String = 'Estn';
        a.x23UI2.String = 'Ekin';
        a.x33UI2.String = 'Etot';
        a.x13UI2.Tooltip = 'Strain energy density';
        a.x23UI2.Tooltip = 'Kinetic energy density';
        a.x33UI2.Tooltip = 'Total energy density';
    elseif a.Quantity22 == 5
        a.PhaseUI2.Enable = 'off';
        a.x11UI2.Enable = 'off';
        a.x12UI2.Enable = 'off';
        a.x13UI2.Enable = 'on';
        a.x22UI2.Enable = 'off';
        a.x23UI2.Enable = 'on';
        a.x33UI2.Enable = 'on';
        a.x13UI2.String = '1';
        a.x23UI2.String = '2';
        a.x33UI2.String = '3';
        a.x13UI2.Tooltip = 'Power flow density p1';
        a.x23UI2.Tooltip = 'Power flow density p2';
        a.x33UI2.Tooltip = 'Power flow density p3';
    end
elseif strcmp(source.Tag,'26') % Mode (Through thickness profiles)
    a.Mode12 = a.ModeNames2{source.Value};
elseif strcmp(source.Tag,'27') % Frequency (Through thickness profiles)
    a.Frequency12 = str2double(source.String);
elseif strcmp(source.Tag,'28') % Samples per layer (Through thickness profiles)
    a.Samples12 = str2double(source.String);
elseif strcmp(source.Tag,'11') % Half-spaces (Through thickness profiles)
    a.HalfspacesNumber12 = str2double(source.String);
elseif strcmp(source.Tag,'12') % Half-spaces (Through thickness profiles)
    a.Halfspaces12 = source.Value;
    if  a.Halfspaces12
        a.HalfspacesNumber1UI2.Enable = 'on';
    else
        a.HalfspacesNumber1UI2.Enable = 'off';
    end
elseif strcmp(source.Tag,'81') % Phase
    a.Phase2 = source.Value;
elseif strcmp(source.Tag,'30') % 11
    a.Plot2(1) = source.Value;
elseif strcmp(source.Tag,'31') % 22
    a.Plot2(2) = source.Value;
elseif strcmp(source.Tag,'32') % 33
    a.Plot2(3) = source.Value;
elseif strcmp(source.Tag,'34') % 23
    a.Plot2(4) = source.Value;
elseif strcmp(source.Tag,'33') % 13
    a.Plot2(5) = source.Value;
elseif strcmp(source.Tag,'35') % 12
    a.Plot2(6) = source.Value; 
elseif strcmp(source.Tag,'29') % Plot (Through thickness profiles)
    try
        if  a.Quantity22 == 1
            u1Color = [1 0 0];
            u2Color = [.13 .55 .13];
            u3Color = [0 0 1];
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(1,1,0,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,u1Color,u2Color,u3Color,0,0,0,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2) 
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(1,1,0,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,u1Color,u2Color,u3Color,0,0,0,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(1,1,0,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,u1Color,u2Color,u3Color,0,0,0,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Phase2)
            end
        elseif a.Quantity22 == 2
            sigma33Color = [0 0 1];
            sigma11Color = [1 0 0];
            sigma22Color = [.13 .55 .13];
            sigma13Color = [0 0 0];
            sigma23Color = [1 0 1];
            sigma12Color = [0 1 1];
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(2,1,0,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,sigma11Color,sigma22Color,sigma33Color,sigma23Color,sigma13Color,sigma12Color,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(2,1,0,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,sigma11Color,sigma22Color,sigma33Color,sigma23Color,sigma13Color,sigma12Color,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(2,1,0,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,sigma11Color,sigma22Color,sigma33Color,sigma23Color,sigma13Color,sigma12Color,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Phase2)
            end
        elseif a.Quantity22 == 3
            epsilon33Color = [0 0 1];
            epsilon11Color = [1 0 0];
            epsilon22Color = [.13 .55 .13];
            epsilon13Color = [0 0 0];
            epsilon23Color = [1 0 1];
            epsilon12Color = [0 1 1];
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(3,1,0,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,epsilon11Color,epsilon22Color,epsilon33Color,epsilon23Color,epsilon13Color,epsilon12Color,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(3,1,0,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,epsilon11Color,epsilon22Color,epsilon33Color,epsilon23Color,epsilon13Color,epsilon12Color,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(3,1,0,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,epsilon11Color,epsilon22Color,epsilon33Color,epsilon23Color,epsilon13Color,epsilon12Color,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Phase2)
            end
        elseif a.Quantity22 == 4
            StrainEnergyDensityColor = [0 0 1];
            KineticEnergyDensityColor = [1 0 0];
            TotalEnergyDensityColor = [0 0 0];
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(4,1,0,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,StrainEnergyDensityColor,KineticEnergyDensityColor,TotalEnergyDensityColor,0,0,0,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,0)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(4,1,0,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,StrainEnergyDensityColor,KineticEnergyDensityColor,TotalEnergyDensityColor,0,0,0,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,0)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(4,1,0,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,StrainEnergyDensityColor,KineticEnergyDensityColor,TotalEnergyDensityColor,0,0,0,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,0)
            end
        elseif a.Quantity22 == 5
            P1Color = [1 0 0];
            P2Color = [.13 .55 .13];
            P3Color = [0 0 1];
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(5,1,0,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,P1Color,P2Color,P3Color,0,0,0,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,0)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(5,1,0,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,P1Color,P2Color,P3Color,0,0,0,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,0)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(5,1,0,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,P1Color,P2Color,P3Color,0,0,0,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,0)
            end
        end
    catch ME
        switch ME.identifier
        case 'MATLAB:nonExistentField'
            errordlg('You must trace dispersion curves before you can plot any data. Press ''Traces modes''.','Trace dispersion curves')
        otherwise
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Error')
        end
    end
elseif strcmp(source.Tag,'85') % Live-plot (Through thickness profiles)
    try
        ModeShapeLive(a.Geometry2,a.Hybrid,a.LayupStringTex2,a.PropagationAngle,a.FrequencyLimit2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.BLamb2,a.BShear2,a.CLamb2,a.CShear2,a.FlexuralModes2,a.FlexuralModeOrders2,a.LongitudinalModes2,a.TorsionalModes2,a.F2,a.L2,a.T2,a.LineColors2,a.LambModes2,a.Material2,a.Radius/1e3,a.PlateThickness/1e3,2*a.InnerRadius/1e3,a.ShearHorizontalModes2,a.SLamb2,a.SShear2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Decoupled,a.DC,a.c,a.Delta,a.Layers2,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.Pattern,a.LayerThicknesses/1e3)    
    catch ME
        switch ME.identifier
        case 'MATLAB:nonExistentField'
            errordlg('You must trace dispersion curves before you can plot any data. Press ''Traces modes''.','Trace dispersion curves')
        otherwise
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Error')
        end
    end
elseif strcmp(source.Tag,'36') % Mode (Mode shape)
    a.Mode22 = a.ModeNames2{source.Value};
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'45') % Plane
    if  source.Value == 1
        a.Plane2 = 'r-z';
    elseif source.Value == 2
        a.Plane2 = 'r-theta';
    end
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'37') % Frequency (Mode shape)
    a.Frequency22 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'38') % Wavelengths
    a.Length2 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'39') % Samples x1
    a.Samples22 = str2double(source.String);
    a.GridLine2 = 2*round(a.Samples22/80);
    if  ~a.GridLine2
        a.GridLine2 = 1;
    end
    a.GridLineUI2.String = a.GridLine2;
    a.ModeShapeSettingChanged2 = 1;
    AdjustSamples32
elseif strcmp(source.Tag,'40') % Samples per layer (Mode shape)
    a.Samples32 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'41') % Gain
    a.Gain2 = str2double(source.String);
elseif strcmp(source.Tag,'42') % Grid line
    a.GridLine2 = str2double(source.String);
elseif strcmp(source.Tag,'43') % Undistorted
    a.Undistorted2 = source.Value;
elseif strcmp(source.Tag,'6') || strcmp(source.Tag,'8') % Half-spaces (mode shape)
    if  strcmp(source.Tag,'6')
        a.HalfspacesNumber22 = str2double(source.String);
    elseif strcmp(source.Tag,'8')
        a.Halfspaces22 = source.Value;
        if  a.Halfspaces22
            a.HalfspacesNumber2UI2.Enable = 'on';
        else
            a.HalfspacesNumber2UI2.Enable = 'off';
        end
    end
    a.ModeShapeSettingChanged2 = 1;
    AdjustSamples32
elseif strcmp(source.Tag,'46') % Cycle duration
    a.CycleDuration2 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'47') % Frame rate
    a.FrameRate2 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'48') % Movie quality
    a.MovieQuality2 = str2double(source.String);
elseif strcmp(source.Tag,'49') % Animate
    a.Animate2 = source.Value;
elseif strcmp(source.Tag,'44') % Plot (Mode shape)
    try
        if  a.Animate2
            if  a.ModeShapeSettingChanged2
                a.ModeShapeSettingChanged2 = 0;
                if  strcmp(a.Geometry2,'Plate')
                    [a.u2_,a.x12_,a.x32_,a.p2_] = ModeShapeGridComputer_Plate(a.Material2,a.Layers2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.c,a.Delta,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.Frequency22,a.Length2,a.Mode22,a.PlateThickness/1e3,a.Samples22,a.Samples32,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces22,a.HalfspacesNumber22,a.CycleDuration2,a.FrameRate2);
                elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                    [a.u2_,a.ua2_,a.ub2_,a.x12_,a.r2_,a.Thetaa2_,a.Thetab2_,a.p2_] = ModeShapeGridComputer_Rod_Pipe(a.Geometry2,a.Plane2,a.Material2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,a.F2,a.L2,a.T2,a.Frequency22,a.Length2,a.Mode22,a.Layers2,a.Radius/1e3,a.Samples22,a.Samples32,a.Halfspaces22,a.HalfspacesNumber22,a.CycleDuration2,a.FrameRate2);
                elseif strcmp(a.Geometry2,'Circumferential')
                    [a.u2_,a.x12_,a.r2_,a.p2_] = ModeShapeGridComputer_Circumferential(a.Material2,a.CLamb2,a.CShear2,a.Frequency22,a.Length2,a.Mode22,a.Layers2,a.Radius/1e3,a.Samples22,a.Samples32,a.CycleDuration2,a.FrameRate2);
                end
            end
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeGrid_Plate(a.Hybrid,a.Layers2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Directory2,a.TitleFontSize2,a.AxesLabelFontSize2,a.Frequency22,a.GridLine2,a.Title2,a.LineWidth2,a.Material2,a.Mode22,a.FileName2,a.PlateThickness,a.PropagationAngle,a.Samples32,a.Gain2,a.u2_,a.Undistorted2,a.x12_,a.x32_,a.p2_,a.Halfspaces22,a.HalfspacesNumber22,0,0,0,0,a.FrameRate2,a.ExportPlots2,a.MovieQuality2)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeGrid_Rod_Pipe(a.Hybrid,a.LayupStringTex2,a.Geometry2,a.Plane2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Directory2,a.TitleFontSize2,a.AxesLabelFontSize2,a.Frequency22,a.GridLine2,a.Title2,a.LineWidth2,a.Material2,a.Mode22,a.FileName2,a.Radius(end),a.InnerRadius,a.Samples32,a.Gain2,a.u2_,a.ua2_,a.ub2_,a.Undistorted2,a.x12_,a.r2_,a.Thetaa2_,a.Thetab2_,a.p2_,a.Halfspaces22,a.HalfspacesNumber22,0,0,0,0,a.FrameRate2,a.ExportPlots2,a.MovieQuality2)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeGrid_Circumferential(a.Hybrid,a.LayupStringTex2,a.Directory2,a.TitleFontSize2,a.AxesLabelFontSize2,a.Frequency22,a.GridLine2,a.Title2,a.LineWidth2,a.Material2,a.Mode22,a.FileName2,a.Radius(end),a.InnerRadius,a.Gain2,a.u2_,a.Undistorted2,a.x12_,a.r2_,a.p2_,0,0,0,0,a.FrameRate2,a.ExportPlots2,a.MovieQuality2)
            end
        else
            if  strcmp(a.Geometry2,'Plate')
                [a.u2,a.x12,a.x32,a.p2] = ModeShapeGridComputer_Plate(a.Material2,a.Layers2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.c,a.Delta,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.Frequency22,a.Length2,a.Mode22,a.PlateThickness/1e3,a.Samples22,a.Samples32,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces22,a.HalfspacesNumber22,0,0);
                ModeShapeGrid_Plate(a.Hybrid,a.Layers2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Directory2,a.TitleFontSize2,a.AxesLabelFontSize2,a.Frequency22,a.GridLine2,a.Title2,a.LineWidth2,a.Material2,a.Mode22,a.FileName2,a.PlateThickness,a.PropagationAngle,a.Samples32,a.Gain2,a.u2,a.Undistorted2,a.x12,a.x32,a.p2,a.Halfspaces22,a.HalfspacesNumber22,a.PNGresolution2,a.ExportPlots2,a.PDF2,a.PNG2,0,0,0)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                [a.u2,a.ua2,a.ub2,a.x12,a.r2,a.Thetaa2,a.Thetab2,a.p2] = ModeShapeGridComputer_Rod_Pipe(a.Geometry2,a.Plane2,a.Material2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,a.F2,a.L2,a.T2,a.Frequency22,a.Length2,a.Mode22,a.Layers2,a.Radius/1e3,a.Samples22,a.Samples32,a.Halfspaces22,a.HalfspacesNumber22,0,0);
                ModeShapeGrid_Rod_Pipe(a.Hybrid,a.LayupStringTex2,a.Geometry2,a.Plane2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Directory2,a.TitleFontSize2,a.AxesLabelFontSize2,a.Frequency22,a.GridLine2,a.Title2,a.LineWidth2,a.Material2,a.Mode22,a.FileName2,a.Radius(end),a.InnerRadius,a.Samples32,a.Gain2,a.u2,a.ua2,a.ub2,a.Undistorted2,a.x12,a.r2,a.Thetaa2,a.Thetab2,a.p2,a.Halfspaces22,a.HalfspacesNumber22,a.PNGresolution2,a.ExportPlots2,a.PDF2,a.PNG2,0,0,0)
            elseif strcmp(a.Geometry2,'Circumferential')
                [a.u2,a.x12,a.r2,a.p2] = ModeShapeGridComputer_Circumferential(a.Material2,a.CLamb2,a.CShear2,a.Frequency22,a.Length2,a.Mode22,a.Layers2,a.Radius/1e3,a.Samples22,a.Samples32,0,0);
                ModeShapeGrid_Circumferential(a.Hybrid,a.LayupStringTex2,a.Directory2,a.TitleFontSize2,a.AxesLabelFontSize2,a.Frequency22,a.GridLine2,a.Title2,a.LineWidth2,a.Material2,a.Mode22,a.FileName2,a.Radius(end),a.InnerRadius,a.Gain2,a.u2,a.Undistorted2,a.x12,a.r2,a.p2,a.PNGresolution2,a.ExportPlots2,a.PDF2,a.PNG2,0,0,0)
            end
        end
    catch ME
        switch ME.identifier
        case 'MATLAB:nonExistentField'
            errordlg('You must trace dispersion curves before you can plot any data. Press ''Traces modes''.','Trace dispersion curves')
        otherwise
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Error')
        end
    end
elseif strcmp(source.Tag,'50') % Export plots
    a.ExportPlots2 = source.Value;
    if  a.ExportPlots2
        a.Plot1UI2.String = 'Export';
        a.Plot2UI2.String = 'Export';
        a.Plot3UI2.String = 'Export';
    else
        a.Plot1UI2.String = 'Plot';
        a.Plot2UI2.String = 'Plot';
        a.Plot3UI2.String = 'Plot';
    end
elseif strcmp(source.Tag,'51') % PDF
    a.PDF2 = source.Value;
elseif strcmp(source.Tag,'52') % PNG
    a.PNG2 = source.Value;
elseif strcmp(source.Tag,'53') % PNG resolution
    a.PNGresolution2 = str2double(source.String);
elseif strcmp(source.Tag,'54') % X-Axis mode (Export settings)
    a.XAxisMode22 = source.Value;
elseif strcmp(source.Tag,'55') % Arrange
    a.Arrange2 = source.Value;
elseif strcmp(source.Tag,'56') % Dispersion curves
    a.DispersionCurves2 = source.Value;
    if  a.DispersionCurves2
        a.XAxisMode2UI2.Enable = 'on';
        a.ArrangeUI2.Enable = 'on';
    else
        a.XAxisMode2UI2.Enable = 'off';
        a.ArrangeUI2.Enable = 'off';
    end
elseif strcmp(source.Tag,'57') % Through thickness
    a.ThroughThickness2 = source.Value;
elseif strcmp(source.Tag,'80') || strcmp(source.Tag,'58') || strcmp(source.Tag,'59') % Matlab Excel txt
    if  strcmp(source.Tag,'80') % mat
        DataFormat = 1;
    elseif strcmp(source.Tag,'58') % xlsx
        DataFormat = 2;
    elseif strcmp(source.Tag,'59') % txt
        DataFormat = 3;
    end
    try
        if  a.DispersionCurves2
            Export(a.Geometry2,DataFormat,a.SLamb2,a.ALamb2,a.BLamb2,a.SShear2,a.AShear2,a.BShear2,a.F2,a.L2,a.T2,a.CLamb2,a.CShear2,a.Arrange2,a.XAxisMode22,a.Distance2,a.Couplant2,a.PlateThickness,2*a.InnerRadius,a.Directory2,a.FileName2)
        end
        if  a.Quantity22 == 1
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(1,2,DataFormat,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2) 
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(1,2,DataFormat,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,0,0,0,0,0,0,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(1,2,DataFormat,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,0,0,0,0,0,0,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Phase2)
            end
        elseif a.Quantity22 == 2
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(2,2,DataFormat,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(2,2,DataFormat,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,0,0,0,0,0,0,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(2,2,DataFormat,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,0,0,0,0,0,0,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Phase2)
            end
        elseif a.Quantity22 == 3
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(3,2,DataFormat,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(3,2,DataFormat,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,0,0,0,0,0,0,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(3,2,DataFormat,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,0,0,0,0,0,0,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Phase2)
            end
        elseif a.Quantity22 == 4
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(4,2,DataFormat,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,0)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(4,2,DataFormat,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,0,0,0,0,0,0,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,0)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(4,2,DataFormat,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,0,0,0,0,0,0,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,0)
            end
        elseif a.Quantity22 == 5
            if  strcmp(a.Geometry2,'Plate')
                ModeShapeLines_Plate(5,2,DataFormat,a.Hybrid,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupStringTex2,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.ALamb2,a.AShear2,a.BLamb2,a.BShear2,a.SLamb2,a.SShear2,a.BoxLineWidth2,a.c,a.Delta,a.Material2,a.Layers2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.DC,a.Halfspaces12,a.HalfspacesNumber12,0)
            elseif strcmp(a.Geometry2,'Rod') || strcmp(a.Geometry2,'Pipe')
                ModeShapeLines_Rod_Pipe(5,2,DataFormat,a.Geometry2,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Sink2,0,0,0,0,0,0,a.F2,a.L2,a.T2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.Layers2,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,a.Halfspaces12,a.HalfspacesNumber12,0)
            elseif strcmp(a.Geometry2,'Circumferential')
                ModeShapeLines_Circumferential(5,2,DataFormat,a.Hybrid,a.LayupStringTex2,a.Plot2,a.PNGresolution2,a.Material2,0,0,0,0,0,0,a.CLamb2,a.CShear2,a.BoxLineWidth2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,1,a.Radius/1e3,a.Radius(end)/1e3,a.InnerRadius/1e3,a.Samples12,0)
            end
        end
    catch ME
        switch ME.identifier
        case 'MATLAB:nonExistentField'
            errordlg('You must trace dispersion curves before you can export any data. Press ''Traces modes''.','Trace dispersion curves')
        otherwise
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Error')
        end
    end
elseif strcmp(source.Tag,'60') % File name
    a.FileName2 = source.String;
elseif strcmp(source.Tag,'61') % Directory
    a.Directory2 = source.String;
elseif strcmp(source.Tag,'62') % Title
    if  source.Value == 1
        a.Title2 = 2;
    elseif source.Value == 2
        a.Title2 = 1;
    elseif source.Value == 3
        a.Title2 = 0;
    end
elseif strcmp(source.Tag,'64') % Legend location
    if  source.Value == 1
        a.LegendLocation2 = 'out';
    elseif source.Value == 2
        a.LegendLocation2 = 'in';
    end
elseif strcmp(source.Tag,'65') % Box line width
    a.BoxLineWidth2 = str2double(source.String);
elseif strcmp(source.Tag,'66') % Curve line width
    a.LineWidth2 = str2double(source.String);
elseif strcmp(source.Tag,'67') % S
    a.SColor2 = eval(source.String);
elseif strcmp(source.Tag,'68') % A
    a.AColor2 = eval(source.String);
elseif strcmp(source.Tag,'69') % B
    a.BColor2 = eval(source.String);
elseif strcmp(source.Tag,'73') % Title (Font size)
    a.TitleFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'74') % Axes labels
    a.AxesLabelFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'75') % Axes ticks
    a.AxesTickFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'77') % Legend
    a.LegendFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'78') % Default
    a.Title2 = 2;
    a.TitleUI2.Value = 1;
    a.LegendLocation2 = 'out';
    a.LegendLocationUI2.Value = 1;
    a.BoxLineWidth2 = .5;
    a.BoxLineWidthUI2.String = a.BoxLineWidth2;
    a.LineWidth2 = 1;
    a.LineWidthUI2.String = a.LineWidth2;
    a.SColor2 = [1 0 0];
    a.SColorUI2.String = '[1 0 0]';
    a.AColor2 = [0 0 1];
    a.AColorUI2.String = '[0 0 1]';
    a.BColor2 = [.5 0 1];
    a.BColorUI2.String = '[.5 0 1]';
    a.TitleFontSize2 = 30;
    a.TitleFontSizeUI2.String = a.TitleFontSize2;
    a.AxesLabelFontSize2 = 30;
    a.AxesLabelFontSizeUI2.String = a.AxesLabelFontSize2;
    a.AxesTickFontSize2 = 24;
    a.AxesTickFontSizeUI2.String = a.AxesTickFontSize2;
    a.LegendFontSize2 = 24;
    a.LegendFontSizeUI2.String = a.LegendFontSize2;
elseif strcmp(source.Tag,'1001') % OK in Specimen setup menu when Fix is on
    CheckDecoupling
    AdjustSamples32
elseif strcmp(source.Tag,'1002') % OK in Specimen setup menu when Fix is off
    AdjustRange
    CheckDecoupling
    AdjustSamples32
end
function AdjustRange
    if  a.IsotropicSinglelayer
        if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Rod') || (strcmp(a.Geometry2,'Pipe') && a.ToggleLowerFluid2 && ~a.Sink2)
            a.FrequencyLimit2 = round(a.Material2{1}.PlateVelocity*a.XRange2/a.PlateThickness,2,'significant');
        elseif strcmp(a.Geometry2,'Pipe') || strcmp(a.Geometry2,'Circumferential')
            a.FrequencyLimit2 = round(a.Material2{1}.PlateVelocity*a.XRange2/(a.PlateThickness-2*a.InnerRadius)*2,2,'significant');
        end
        a.PhaseVelocityLimit2 = round(a.Material2{1}.PlateVelocity*a.YRange2,2,'significant');
        a.Step2 = a.FrequencyLimit2/a.Steps1;
    else
        if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Rod') || (strcmp(a.Geometry2,'Pipe') && a.ToggleLowerFluid2 && ~a.Sink2)
            a.FrequencyLimit2 = round(1e3*a.XRange2/a.PlateThickness,2,'significant');
        elseif strcmp(a.Geometry2,'Pipe') || strcmp(a.Geometry2,'Circumferential')
            a.FrequencyLimit2 = round(1e3*a.XRange2/(a.PlateThickness-2*a.InnerRadius)*2,2,'significant');
        end
        a.PhaseVelocityLimit2 = 2e4;
        a.Step2 = a.FrequencyLimit2/a.Steps2;
    end
    a.FrequencyResolution2 = a.FrequencyLimit2/a.XSamples2;
    a.FrequencyLimitUI2.String = a.FrequencyLimit2;
    a.FrequencyResolutionUI2.String = a.FrequencyResolution2;
    a.PhaseVelocityLimitUI2.String = a.PhaseVelocityLimit2/1e3;
    a.StepUI2.String = a.Step2;
    a.Frequency12 = a.FrequencyLimit2;
    a.Frequency1UI2.String = a.FrequencyLimit2;
    a.Frequency22 = a.FrequencyLimit2;
    a.Frequency2UI2.String = a.FrequencyLimit2;
    AdjustAxes
end
function AdjustAxes
    if  strcmp(a.Geometry2,'Plate') || strcmp(a.Geometry2,'Rod')
        Thickness = a.PlateThickness;
    elseif strcmp(a.Geometry2,'Pipe') || strcmp(a.Geometry2,'Circumferential')
        Thickness = a.LayupThickness;
    end
    if  a.Quantity12 == 3
        if  a.XAxisMode2 == 1
            a.YAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        elseif a.XAxisMode2 == 2
            a.YAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.YAxis2 = eval(a.YAxisUI2.String)*1e3;
        elseif a.XAxisMode2 == 3
            a.YAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*Thickness),']'];
            a.YAxis2 = eval(a.YAxisUI2.String)*1e3/Thickness;
        end
    else
        if  a.XAxisMode2 == 1
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.XAxis2 = eval(a.XAxisUI2.String);
        elseif a.XAxisMode2 == 2
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
        elseif a.XAxisMode2 == 3
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*Thickness),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3/Thickness;
        end
    end
    if  a.Quantity12 == 1
        a.YAxisUI2.String = ['[0 ',num2str(a.PhaseVelocityLimit2/1e3),']'];
        a.YAxis2 = eval(a.YAxisUI2.String);
    elseif a.Quantity12 == 2
        if  a.IsotropicSinglelayer
            a.YAxisUI2.String = ['[0 ',num2str(ceil(a.Material2{1}.PlateVelocity/1e3)),']'];
        else
            a.YAxisUI2.String = '[0 11]';
        end
        a.YAxis2 = eval(a.YAxisUI2.String);
    elseif a.Quantity12 == 3
        if  a.IsotropicSinglelayer
            x = round(a.Distance2/a.Material2{1}.PlateVelocity*5e3,2,'significant');
        else
            x = round(2*a.Distance2,2,'significant');
        end
        a.XAxisUI2.String = ['[0 ',num2str(x),']'];
        a.XAxis2 = eval(a.XAxisUI2.String);
    elseif a.Quantity12 == 4
        if  a.IsotropicSinglelayer
            x = round(8e4/a.Material2{1}.PlateVelocity*a.Couplant2.Velocity/343);
        else
            x = round(30*a.Couplant2.Velocity/343);
        end
        if  x > 90
            x = 90;
        end
        a.YAxisUI2.String = ['[0 ',num2str(x),']'];
        a.YAxis2 = eval(a.YAxisUI2.String);
    elseif a.Quantity12 == 5
        if  a.IsotropicSinglelayer
            x = round(a.YRange2/1e3*a.Material2{1}.PlateVelocity*Thickness,2,'significant');
        end
        if  a.XAxisMode2 == 3
            if  a.IsotropicSinglelayer
                a.YAxisUI2.String = ['[0 ',num2str(x/Thickness),']'];
            else
                a.YAxisUI2.String = '[0 50]';
            end
            a.YAxis2 = eval(a.YAxisUI2.String);
        else
            if  a.IsotropicSinglelayer
                a.YAxisUI2.String = ['[0 ',num2str(x),']'];
            else
                a.YAxisUI2.String = ['[0 ',num2str(50*Thickness),']'];
            end
            a.YAxis2 = eval(a.YAxisUI2.String);
        end
    else
        if  a.Quantity12 == 6
            if  a.IsotropicSinglelayer
                x = round(2*pi*a.FrequencyLimit2/a.Material2{1}.RayleighVelocity,2,'significant');
            else
                x = round(2*pi*a.FrequencyLimit2/1.6e3,2,'significant');
            end
        elseif a.Quantity12 == 7
            if  a.Symmetric2
                if  strcmp(a.Material2{1}.Class,'Isotropic')
                    MaterialVelocity = a.Material2{1}.PlateVelocity;
                else
                    MaterialVelocity = a.Material2{1}.LongitudinalVelocity_1;
                end
                MaterialDensity = a.Material2{1}.Density;
                FluidVelocity = a.UpperFluid2.Velocity;
                FluidDensity = a.UpperFluid2.Density;
            else
                if  strcmp(a.Material2{1}.Class,'Isotropic')
                    UpperMaterialVelocity = a.Material2{1}.PlateVelocity;
                else
                    UpperMaterialVelocity = a.Material2{1}.LongitudinalVelocity_1;
                end
                if  strcmp(a.Material2{end}.Class,'Isotropic')
                    LowerMaterialVelocity = a.Material2{end}.PlateVelocity;
                else
                    LowerMaterialVelocity = a.Material2{end}.LongitudinalVelocity_1;
                end
                MaterialVelocity = (UpperMaterialVelocity+LowerMaterialVelocity)/2;
                MaterialDensity = (a.Material2{1}.Density+a.Material2{end}.Density)/2;
                if  a.ToggleUpperFluid2 && a.ToggleLowerFluid2
                    FluidVelocity = (a.UpperFluid2.Velocity+a.LowerFluid2.Velocity)/2;
                    FluidDensity = (a.UpperFluid2.Density+a.LowerFluid2.Density)/2;
                elseif a.ToggleUpperFluid2 && ~a.ToggleLowerFluid2
                    FluidVelocity = a.UpperFluid2.Velocity/2;
                    FluidDensity = a.UpperFluid2.Density/2;
                elseif ~a.ToggleUpperFluid2 && a.ToggleLowerFluid2
                    FluidVelocity = a.LowerFluid2.Velocity/2;
                    FluidDensity = a.LowerFluid2.Density/2;
                end
            end
            if  a.Viscoelastic2
                for i = 1:a.SuperLayerSize
                    if  ~isreal(a.Material2{i}.C)
                        xV = pi*(imag(a.Material2{i}.C(1,1))/real(a.Material2{i}.C(1,1))+imag(a.Material2{i}.C(6,6))/real(a.Material2{i}.C(6,6)));
                        break
                    end
                end
            end
            if  a.FluidLoading2 && a.Viscoelastic2
                x = round(1e4*(FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity)+xV)/Thickness,2,'significant');
            elseif a.FluidLoading2 && ~a.Viscoelastic2
                x = round(1e4*FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity)/Thickness,2,'significant');
            elseif ~a.FluidLoading2 && a.Viscoelastic2
                x = round(1e4*xV/Thickness,2,'significant');
            else
                x = 1;
            end
        end
        if  a.XAxisMode2 == 3
            a.YAxisUI2.String = ['[0 ',num2str(x*Thickness),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        else
            a.YAxisUI2.String = ['[0 ',num2str(x),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        end
    end
end
function CheckDecoupling
    a.DC = [];
    a.Phi = a.LayerOrientations-a.PropagationAngle;
    for i = 1:a.SuperLayerSize
        if  strcmp(a.MaterialClasses(i),'Isotropic')
            a.DC(i) = 1;
        elseif strcmp(a.MaterialClasses(i),'Cubic')
            if  mod(a.Phi(i),45)
                a.DC(i) = 0;
            else
                a.DC(i) = 1;
            end
        elseif strcmp(a.MaterialClasses(i),'Transversely isotropic') || strcmp(a.MaterialClasses(i),'Orthotropic')
            if  mod(a.Phi(i),90)
                a.DC(i) = 0;
            else
                a.DC(i) = 1;
            end
        end
    end
    if  all(a.DC)
        a.Decoupled = 1;
        a.ShearHorizontalModes_FlexuralModeOrdersUI2.Enable = 'on';
    else
        a.Decoupled = 0;
        a.ShearHorizontalModes_FlexuralModeOrdersUI2.Enable = 'off';
    end
end
function TransformStiffness
    for i = 1:a.SuperLayerSize % transformed layer stiffness matrices for orthotropic materials [2], p. 28; transverse isotropy is included
        C = conj(a.Material2{i}.C);
        s = sind(a.Phi(i));
        g = cosd(a.Phi(i));
        a.c{i}(1,1) = C(1,1)*g^4+C(2,2)*s^4+2*(C(1,2)+2*C(6,6))*s^2*g^2;
        a.c{i}(1,2) = (C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2+C(1,2);
        a.c{i}(1,3) = C(1,3)*g^2+C(2,3)*s^2;
        a.c{i}(1,6) = (C(1,2)+2*C(6,6)-C(1,1))*s*g^3+(C(2,2)-C(1,2)-2*C(6,6))*g*s^3;
        a.c{i}(2,2) = C(1,1)*s^4+C(2,2)*g^4+2*(C(1,2)+2*C(6,6))*s^2*g^2;
        a.c{i}(2,3) = C(2,3)*g^2+C(1,3)*s^2;
        a.c{i}(2,6) = (C(1,2)+2*C(6,6)-C(1,1))*g*s^3+(C(2,2)-C(1,2)-2*C(6,6))*s*g^3;
        a.c{i}(3,3) = C(3,3);
        a.c{i}(3,6) = (C(2,3)-C(1,3))*s*g;
        a.c{i}(4,4) = C(4,4)*g^2+C(5,5)*s^2;
        a.c{i}(4,5) = (C(4,4)-C(5,5))*s*g;
        a.c{i}(5,5) = C(5,5)*g^2+C(4,4)*s^2;
        a.c{i}(6,6) = C(6,6)+(C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2;
        a.Delta(i) = a.c{i}(3,3)*a.c{i}(4,4)*a.c{i}(5,5)-a.c{i}(3,3)*a.c{i}(4,5)^2;
        a.DeltaReal(i) = real(a.c{i}(3,3))*real(a.c{i}(4,4))*real(a.c{i}(5,5))-real(a.c{i}(3,3))*real(a.c{i}(4,5))^2;
    end
end
function AdjustSamples32
    if  strcmp(a.Geometry2,'Plate')
        if  a.FluidLoading2 && a.Halfspaces22
            if  a.ToggleUpperFluid2 && a.ToggleLowerFluid2
                a.Samples32 = round(a.Samples22/(2*a.HalfspacesNumber22+1)/a.Layers2/2);
            else
                a.Samples32 = round(a.Samples22/(a.HalfspacesNumber22+1)/a.Layers2/2);
            end
        else
            a.Samples32 = round(a.Samples22/a.Layers2/2);
        end
    elseif strcmp(a.Geometry2,'Rod')
        if  a.FluidLoading2 && a.Halfspaces22
            a.Samples32 = round(a.Samples22/(a.HalfspacesNumber22+1)/a.Layers2/2);
        else
            a.Samples32 = round(a.Samples22/a.Layers2/2);
        end
    elseif strcmp(a.Geometry2,'Pipe')
        if  a.ToggleUpperFluid2 && a.Halfspaces22 && a.ToggleLowerFluid2
            a.Samples32 = round(a.Samples22/(a.HalfspacesNumber22+2)/a.Layers2/2);
        elseif a.ToggleUpperFluid2 && a.Halfspaces22 && ~a.ToggleLowerFluid2
            a.Samples32 = round(a.Samples22/(a.HalfspacesNumber22+1)/a.Layers2/2);
        elseif ((a.ToggleUpperFluid2 && ~a.Halfspaces22) || ~a.ToggleUpperFluid2) && a.ToggleLowerFluid2
            a.Samples32 = round(a.Samples22/a.Layers2/4);
        else
            a.Samples32 = round(a.Samples22/a.Layers2/2);
        end
    elseif strcmp(a.Geometry2,'Circumferential')
        a.Samples32 = round(a.Samples22/a.Layers2/2);
    end
    if  ~a.Samples32
        a.Samples32 = 1;
    end
    a.Samples3UI2.String = a.Samples32;
end
function FrequencySweeps
    FrequencySweepRange = a.Step2:a.Step2:a.FrequencyLimit2;
    if  strcmp(a.Geometry2,'Plate')
        [a.HSLamb2,a.HSShear2,a.HALamb2,a.HAShear2,a.HBLamb2,a.HBShear2] = FrequencySweeper_Anisotropic_Plate(a.c,a.DeltaReal,a.Material2,FrequencySweepRange,a.PhaseVelocityLimit2,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.Pattern,a.SuperLayerSize,a.LayerThicknesses/1e3,a.Symmetric2,a.SymmetricSystem,a.Decoupled,a.DC,a.EvanescenceLimit2,a.OutputWindow1aUI2,a.OutputWindow1bUI2,a.OutputWindow2aUI2,a.OutputWindow2bUI2);
    elseif strcmp(a.Geometry2,'Rod')
        [a.HF2,a.HL2,a.HT2] = FrequencySweeper_Isotropic_Rod_Pipe_Multi(a.Geometry2,a.PhaseVelocityLimit2,FrequencySweepRange,a.FlexuralModeOrders2,a.Material2,a.Layers2,a.Radius/1e3,a.LowerFluid2,0,a.EvanescenceLimit2,a.OutputWindow1aUI2,a.OutputWindow1bUI2,a.OutputWindow2aUI2,a.OutputWindow2bUI2);
    elseif strcmp(a.Geometry2,'Pipe')
        if  a.Force2DTracing2 || a.Viscoelastic2 || a.ToggleUpperFluid2 || (a.ToggleLowerFluid2 && a.Sink2)
            [a.HF2,a.HL2,a.HT2] = FrequencySweeper_Isotropic_Rod_Pipe_Multi(a.Geometry2,a.PhaseVelocityLimit2,FrequencySweepRange,a.FlexuralModeOrders2,a.Material2,a.Layers2,a.Radius/1e3,a.LowerFluid2,0,a.EvanescenceLimit2,a.OutputWindow1aUI2,a.OutputWindow1bUI2,a.OutputWindow2aUI2,a.OutputWindow2bUI2);
            if  a.ToggleLowerFluid2 && ~a.Sink2
                [HF2,HL2] = FrequencySweeper_Isotropic_Rod_Pipe_Multi(a.Geometry2,a.PhaseVelocityLimit2,FrequencySweepRange,a.FlexuralModeOrders2,a.Material2,a.Layers2,a.Radius/1e3,a.LowerFluid2,a.ToggleLowerFluid2,a.EvanescenceLimit2,a.OutputWindow1aUI2,a.OutputWindow1bUI2,a.OutputWindow2aUI2,a.OutputWindow2bUI2);
                for n = 1:length(a.HF2)
                    for i = 1:length(a.HF2{n})
                        if  ~any(isapprox(HF2{n},a.HF2{n}(i),RelativeTolerance=.05))
                            HF2{n} = sort([HF2{n} a.HF2{n}(i)]);
                        end
                    end
                    a.HF2{n} = sort(HF2{n});
                end
                for i = 1:length(a.HL2)
                    if  ~any(isapprox(HL2,a.HL2(i),RelativeTolerance=.05))
                        HL2 = [HL2 a.HL2(i)];
                    end
                end
                a.HL2 = sort(HL2);
            end
        else
            [a.HF2,a.HL2,a.HT2] = FrequencySweeper_Isotropic_Rod_Pipe_Multi(a.Geometry2,a.PhaseVelocityLimit2,FrequencySweepRange,a.FlexuralModeOrders2,a.Material2,a.Layers2,a.Radius/1e3,a.LowerFluid2,a.ToggleLowerFluid2,a.EvanescenceLimit2,a.OutputWindow1aUI2,a.OutputWindow1bUI2,a.OutputWindow2aUI2,a.OutputWindow2bUI2);
        end
    elseif strcmp(a.Geometry2,'Circumferential')
        [a.HCLamb2,a.HCShear2] = FrequencySweeper_Isotropic_Circumferential_Multi(a.PhaseVelocityLimit2,FrequencySweepRange,a.Material2,a.Layers2,a.Radius/1e3,a.EvanescenceLimit2,a.OutputWindow1aUI2,a.OutputWindow1bUI2,a.OutputWindow2aUI2,a.OutputWindow2bUI2);
    end
    a.Search2 = 1;
end
end