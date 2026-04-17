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
function a = CallbackModule_Isotropic(source,~,a,Tab3)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
if  ~strcmp(source.Tag,'24') && ~strcmp(source.Tag,'34') % there are commas like in L(0,1)
    source.String = replace(source.String,',','.');
end
if  strcmp(source.Tag,'6') % Geometry
    a.Geometry1 = source.String(source.Value);
    switch source.Value
    case 1 % Plate
        a.UpperFluidTextUI1.String = 'Upper fluid';
        a.LowerFluidTextUI1.String = 'Lower fluid';
        a.ThicknessTextUI1.String = 'Thickness (mm)';
        a.ToggleUpperFluidUI1.Enable = 'on';
        a.ToggleLowerFluidUI1.Enable = 'on';
        if  a.ToggleUpperFluid1
            a.SelectUpperFluidUI1.Enable = 'on';
        end
        if  a.ToggleLowerFluid1
            a.SelectLowerFluidUI1.Enable = 'on';
        end
        a.SinkUI1.Enable = 'off';
        a.ThicknessInnerUI1.Enable = 'off';
        a.Symmetric_Torsional_ModesUI1.Enable = 'on';
        a.Antisymmetric_Longitudinal_ModesUI1.Enable = 'on';
        a.Symmetric_Torsional_ModesUI1.Tooltip = 'Check this in order to calculate the symmetric modes. These modes have a symmetric displacement pattern with respect to the middle plane of the plate.';
        a.Antisymmetric_Longitudinal_ModesUI1.Tooltip = 'Check this in order to calculate the antisymmetric modes. These modes have an antisymmetric displacement pattern with respect to the middle plane of the plate.';
        a.Lamb_Flexural_ModesUI1.Tooltip = 'Check this in order to calculate the Lamb wave modes. These modes show displacement only in the sagittal plane spanned by the propagation direction x1 and by the out-of-plane direction x3. These kind of waves are termed ''pure'' Lamb waves. Lamb waves are indicated by solid lines in the dispersion diagram.';
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.Tooltip = 'Check this in order to calculate the shear horizontal modes. These modes show displacement only perpendicular (x2) to the propagation direction x1, and are therefore termed ''pure'' modes. Shear horizontal waves are indicated by dashed lines in the dispersion diagram.';
        a.SamplesX3ceTextUI1.String = 'Samples x3';
        a.SamplesX3TextUI1.String = 'Samples x3';
        a.SamplesX1X3TextUI1.String = 'Samples x1,x3';
        a.SColorTextUI1.String = 'S';
        a.AColorTextUI1.String = 'A';
        a.BColorTextUI1.String = 'B';
    case 2 % Rod
        a.UpperFluidTextUI1.String = 'Outer fluid';
        a.ThicknessTextUI1.String = 'Outer diameter (mm)';
        a.ToggleUpperFluidUI1.Enable = 'on';
        a.ToggleLowerFluidUI1.Enable = 'off';
        if  a.ToggleUpperFluid1
            a.SelectUpperFluidUI1.Enable = 'on';
        end
        a.SelectLowerFluidUI1.Enable = 'off';
        a.SinkUI1.Enable = 'off';
        a.ThicknessInnerUI1.Enable = 'off';
        a.Symmetric_Torsional_ModesUI1.Enable = 'on';
        a.Antisymmetric_Longitudinal_ModesUI1.Enable = 'on';
        a.Symmetric_Torsional_ModesUI1.Tooltip = 'Check this in order to calculate the torsional modes. These modes show displacement only in the circumferential direction. Torsional waves are indicated by dashed lines in the dispersion diagram.';
        a.Antisymmetric_Longitudinal_ModesUI1.Tooltip = 'Check this in order to calculate the longitudinal modes. These modes have displacement components only in the axial and radial directions.';
        a.Lamb_Flexural_ModesUI1.Tooltip = 'Check this in order to calculate the flexural modes. These modes have displacement components in all three directions (axial, radial, and circumferential).';
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.Tooltip = 'For flexural modes, an infinite number n = 1,2,3,... of circumferential wavenumbers can fit into the circumference of a pipe or rod. For torsional and longitudinal modes, n = 0. For each family of n, an infinite number p = 1,2,3,... of higher order modes exists, similarly as we know it from the modes in a plate. The naming for modes in pipes and rods is T(0,p), L(0,p), F(n,p). Enter the maximum number of mode families n to be calculated.';
        a.SamplesX3ceTextUI1.String = 'Samples r';
        a.SamplesX3TextUI1.String = 'Samples r';
        a.SamplesX1X3TextUI1.String = 'Samples z,r';
        a.SColorTextUI1.String = 'T,L';
        a.AColorTextUI1.String = 'F';
    case 3 % Pipe
        a.UpperFluidTextUI1.String = 'Outer fluid';
        a.LowerFluidTextUI1.String = 'Inner fluid';
        a.ThicknessTextUI1.String = 'Outer diameter (mm)';
        a.ToggleUpperFluidUI1.Enable = 'on';
        a.ToggleLowerFluidUI1.Enable = 'on';
        if  a.ToggleUpperFluid1
            a.SelectUpperFluidUI1.Enable = 'on';
        end
        if  a.ToggleLowerFluid1
            a.SelectLowerFluidUI1.Enable = 'on';
            a.SinkUI1.Enable = 'on';
        else
            a.SinkUI1.Enable = 'off';
        end
        a.ThicknessInnerUI1.Enable = 'on';
        if  a.Thickness1 <= a.ThicknessInner1
            a.ThicknessInner1 = .8*a.Thickness1;
            a.ThicknessInnerUI1.String = a.ThicknessInner1;
        end
        a.Symmetric_Torsional_ModesUI1.Enable = 'on';
        a.Antisymmetric_Longitudinal_ModesUI1.Enable = 'on';
        a.Symmetric_Torsional_ModesUI1.Tooltip = 'Check this in order to calculate the torsional modes. These modes show displacement only in the circumferential direction. Torsional waves are indicated by dashed lines in the dispersion diagram.';
        a.Antisymmetric_Longitudinal_ModesUI1.Tooltip = 'Check this in order to calculate the longitudinal modes. These modes have displacement components only in the axial and radial directions.';
        a.Lamb_Flexural_ModesUI1.Tooltip = 'Check this in order to calculate the flexural modes. These modes have displacement components in all three directions (axial, radial, and circumferential).';
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.Tooltip = 'For flexural modes, an infinite number n = 1,2,3,... of circumferential wavenumbers can fit into the circumference of a pipe or rod. For torsional and longitudinal modes, n = 0. For each family of n, an infinite number p = 1,2,3,... of higher order modes exists, similarly as we know it from the modes in a plate. The naming for modes in pipes and rods is T(0,p), L(0,p), F(n,p). Enter the maximum number of mode families n to be calculated.';
        a.SamplesX3ceTextUI1.String = 'Samples r';
        a.SamplesX3TextUI1.String = 'Samples r';
        a.SamplesX1X3TextUI1.String = 'Samples z,r';
        a.SColorTextUI1.String = 'T,L';
        a.AColorTextUI1.String = 'F';
    case 4 % Circumferential
        a.ThicknessTextUI1.String = 'Outer diameter (mm)';
        a.ToggleUpperFluidUI1.Enable = 'off';
        a.ToggleLowerFluidUI1.Enable = 'off';
        a.SelectUpperFluidUI1.Enable = 'off';
        a.SelectLowerFluidUI1.Enable = 'off';
        a.SinkUI1.Enable = 'off';
        a.ThicknessInnerUI1.Enable = 'on';
        if  a.Thickness1 <= a.ThicknessInner1
            a.ThicknessInner1 = .8*a.Thickness1;
            a.ThicknessInnerUI1.String = a.ThicknessInner1;
        end
        a.Symmetric_Torsional_ModesUI1.Enable = 'off';
        a.Antisymmetric_Longitudinal_ModesUI1.Enable = 'off';
        a.Lamb_Flexural_ModesUI1.Tooltip = 'Check this in order to calculate the Lamb-type modes. These modes show displacement only in the sagittal plane spanned by the cirumferential direction and by the radial direction. Lamb-type waves are indicated by solid lines in the dispersion diagram.';
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.Tooltip = 'Check this in order to calculate the shear horizontal-type modes. These modes show displacement only in the axial direction. Shear horizontal-type waves are indicated by dashed lines in the dispersion diagram.';
        a.SamplesX3ceTextUI1.String = 'Samples r';
        a.SamplesX3TextUI1.String = 'Samples r';
        a.SamplesX1X3TextUI1.String = 'Samples z,r';
        a.BColorTextUI1.String = 'C';
    end
    if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Circumferential')
        a.Symmetric_Torsional_ModesTextUI1.String = 'Symmetric modes';
        a.Antisymmetric_Longitudinal_ModesTextUI1.String = 'Antisymmetric modes';
        a.Lamb_Flexural_ModesTextUI1.String = 'Lamb modes';
        a.ShearHorizontalModes_FlexuralModeOrdersTextUI1.String = 'Shear horizontal modes';
        a.Symmetric_Torsional_ModesUI1.Value = a.SymmetricModes1;
        a.Antisymmetric_Longitudinal_ModesUI1.Value = a.AntisymmetricModes1;
        a.Lamb_Flexural_ModesUI1.Value = a.LambModes1;
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.Enable = 'on';
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.Style = 'checkbox';
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.String = '';
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.Value = a.ShearHorizontalModes1;
        a.PlaneUI1.Enable = 'off';
    else
        a.Symmetric_Torsional_ModesTextUI1.String = 'Torsional modes';
        a.Antisymmetric_Longitudinal_ModesTextUI1.String = 'Longitudinal modes';
        a.Lamb_Flexural_ModesTextUI1.String = 'Flexural modes';
        a.ShearHorizontalModes_FlexuralModeOrdersTextUI1.String = 'Flexural mode orders';
        a.Symmetric_Torsional_ModesUI1.Value = a.TorsionalModes1;
        a.Antisymmetric_Longitudinal_ModesUI1.Value = a.LongitudinalModes1;
        a.Lamb_Flexural_ModesUI1.Value = a.FlexuralModes1;
        if  a.FlexuralModes1
            a.ShearHorizontalModes_FlexuralModeOrdersUI1.Enable = 'on';
        else
            a.ShearHorizontalModes_FlexuralModeOrdersUI1.Enable = 'off';  
        end
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.Style = 'edit';
        a.ShearHorizontalModes_FlexuralModeOrdersUI1.String = a.FlexuralModeOrders1;
        a.PlaneUI1.Enable = 'on';
    end
    if  (strcmp(a.Geometry1,'Rod') && a.ToggleUpperFluid1) || ((strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Pipe')) && (a.ToggleUpperFluid1 || a.ToggleLowerFluid1))
        a.FluidLoading1 = 1;
        if  strcmp(a.Geometry1,'Pipe') && ~a.ToggleUpperFluid1
            a.HalfspacesNumber1UI1.Enable = 'off';
            a.Halfspaces1UI1.Enable = 'off';
            a.HalfspacesNumber2UI1.Enable = 'off';
            a.Halfspaces2UI1.Enable = 'off';
        else
            if  a.Halfspaces11
                a.HalfspacesNumber1UI1.Enable = 'on';
            end
            a.Halfspaces1UI1.Enable = 'on';
            if  a.Halfspaces21
                a.HalfspacesNumber2UI1.Enable = 'on';
            end
            a.Halfspaces2UI1.Enable = 'on';
        end
        a.Couplant1 = a.UpperFluid1;
        if  a.Quantity11 == 4
            a.Option1UI1.Value = find(strcmp(a.Couplant1.Name,a.Option1UI1.String));
        end
    else
        a.FluidLoading1 = 0;
        a.HalfspacesNumber1UI1.Enable = 'off';
        a.Halfspaces1UI1.Enable = 'off';
        a.HalfspacesNumber2UI1.Enable = 'off';
        a.Halfspaces2UI1.Enable = 'off';
    end
    a.Search1 = 0;
    Enable2DTracingSettings
    if  ~a.Fix1
        AdjustBasicParameters
    end
    CheckSymmetry
    AdjustSamples31
elseif strcmp(source.Tag,'1') % Material
    a.Material1 = getfield(a.Materials.Isotropic,cell2mat(source.String(source.Value)));
    if  isreal(a.Material1.C)
        a.Viscoelastic1 = 0;
    else
        a.Viscoelastic1 = 1;
    end
    a.Search1 = 0;
    Enable2DTracingSettings
    if  ~a.Fix1
        a.PhaseVelocityLimit1 = round(a.YRange1*a.Material1.PlateVelocity,-3);
        a.PhaseVelocityLimitUI1.String = a.PhaseVelocityLimit1/1e3;
        AdjustBasicParameters
    end
elseif strcmp(source.Tag,'2') || strcmp(source.Tag,'90') % Thickness/Inner diameter
    if  strcmp(source.Tag,'2')
        a.Thickness1 = str2double(source.String);
    elseif strcmp(source.Tag,'90')
        a.ThicknessInner1 = str2double(source.String);
    end
    if  (strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')) && a.Thickness1 <= a.ThicknessInner1
        a.ThicknessInner1 = .8*a.Thickness1;
        a.ThicknessInnerUI1.String = a.ThicknessInner1;
    end
    a.Search1 = 0;
    if  ~a.Fix1
        AdjustFrequencyRange
    end
elseif strcmp(source.Tag,'7') || strcmp(source.Tag,'86') % Toggle upper fluid/Toggle lower fluid
    if  strcmp(source.Tag,'7')
        a.ToggleUpperFluid1 = source.Value;
        if  source.Value
            a.SelectUpperFluidUI1.Enable = 'on';
        else
            a.SelectUpperFluidUI1.Enable = 'off';
        end
    elseif strcmp(source.Tag,'86')
        a.ToggleLowerFluid1 = source.Value;
        if  source.Value
            a.SelectLowerFluidUI1.Enable = 'on';
            if  strcmp(a.Geometry1,'Pipe')
                a.SinkUI1.Enable = 'on';
            end
        else
            a.SelectLowerFluidUI1.Enable = 'off';
            if  strcmp(a.Geometry1,'Pipe')
                a.SinkUI1.Enable = 'off';
            end
        end
    end
    if  (strcmp(a.Geometry1,'Rod') && a.ToggleUpperFluid1) || ((strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Pipe')) && (a.ToggleUpperFluid1 || a.ToggleLowerFluid1))
        a.FluidLoading1 = 1;
        if  strcmp(a.Geometry1,'Pipe') && ~a.ToggleUpperFluid1
            a.HalfspacesNumber1UI1.Enable = 'off';
            a.Halfspaces1UI1.Enable = 'off';
            a.HalfspacesNumber2UI1.Enable = 'off';
            a.Halfspaces2UI1.Enable = 'off';
        else
            if  a.Halfspaces11
                a.HalfspacesNumber1UI1.Enable = 'on';
            end
            a.Halfspaces1UI1.Enable = 'on';
            if  a.Halfspaces21
                a.HalfspacesNumber2UI1.Enable = 'on';
            end
            a.Halfspaces2UI1.Enable = 'on';
        end
        a.Couplant1 = a.UpperFluid1;
        if  a.Quantity11 == 4
            a.Option1UI1.Value = find(strcmp(a.Couplant1.Name,a.Option1UI1.String));
        end
    else
        a.FluidLoading1 = 0;
        a.HalfspacesNumber1UI1.Enable = 'off';
        a.Halfspaces1UI1.Enable = 'off';
        a.HalfspacesNumber2UI1.Enable = 'off';
        a.Halfspaces2UI1.Enable = 'off';
    end
    a.Search1 = 0;
    Enable2DTracingSettings
    if  ~a.Fix1
        AdjustBasicParameters
    end
    CheckSymmetry
    AdjustSamples31
elseif strcmp(source.Tag,'87') || strcmp(source.Tag,'88') % Upper fluid/Lower fluid
    if  strcmp(source.Tag,'87')
        a.UpperFluid1 = getfield(a.Materials.Fluid,cell2mat(source.String(source.Value)));
        a.Couplant1 = a.UpperFluid1;
    elseif strcmp(source.Tag,'88')
        a.LowerFluid1 = getfield(a.Materials.Fluid,cell2mat(source.String(source.Value)));
    end
    if  a.Quantity11 == 4
        a.Option1UI1.Value = find(strcmp(a.Couplant1.Name,a.Option1UI1.String));
    end
    a.Search1 = 0;
    if  a.Quantity11 == 4 || a.Quantity11 == 7
        AdjustAxes
    end
    CheckSymmetry
elseif strcmp(source.Tag,'89') % Sink at center
    a.Sink1 = source.Value;
    a.Search1 = 0;
    Enable2DTracingSettings
    if  ~a.Fix1
        AdjustBasicParameters
    end
elseif strcmp(source.Tag,'61') % Fix
    a.Fix1 = source.Value;
    if  ~a.Fix1
        a.Search1 = 0;
        a.PhaseVelocityLimit1 = round(a.YRange1*a.Material1.PlateVelocity,-3);
        a.PhaseVelocityLimitUI1.String = a.PhaseVelocityLimit1/1e3;
        AdjustBasicParameters
    end
elseif strcmp(source.Tag,'4') % Frequency limit
    a.FrequencyLimit1 = str2double(source.String);
    a.Frequency11 = a.FrequencyLimit1;
    a.Frequency1UI1.String = a.FrequencyLimit1;
    a.Frequency21 = a.FrequencyLimit1;
    a.Frequency2UI1.String = a.FrequencyLimit1;
    a.Search1 = 0;
    AdjustAxes
elseif strcmp(source.Tag,'5') % Frequency step
    a.FrequencyResolution1 = str2double(source.String);
elseif strcmp(source.Tag,'3') % Phase velocity limit
    a.PhaseVelocityLimit1 = str2double(source.String)*1e3;
    if  a.Quantity11 == 1
        a.YAxisUI1.String = ['[0 ',source.String,']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    end
    a.Search1 = 0;
elseif strcmp(source.Tag,'85') % Phase velocity accuracy
    a.Accuracy1 = str2double(source.String);
elseif strcmp(source.Tag,'77') % Force 2-D tracing
    a.Force2DTracing1 = source.Value;
    if  source.Value
        a.PhaseVelocityLimit2UI1.Enable = 'on';
        a.AttenuationLimitUI1.Enable = 'on';
        a.SweepsUI1.Enable = 'on';
        a.SweepsSectionsUI1.Enable = 'on';
        a.HigherOrderModesUI1.Enable = 'off';
    else
        a.PhaseVelocityLimit2UI1.Enable = 'off';
        a.AttenuationLimitUI1.Enable = 'off';
        a.SweepsUI1.Enable = 'off';
        a.SweepsSectionsUI1.Enable = 'off';
        a.HigherOrderModesUI1.Enable = 'on';
    end
elseif strcmp(source.Tag,'79') % Phase velocity limit 2
    a.PhaseVelocityLimit21 = str2double(source.String)*1e3;
elseif strcmp(source.Tag,'93') % Attenuation limit
    a.AttenuationLimit1 = str2double(source.String);
elseif strcmp(source.Tag,'94') % Sweeps
    a.Sweeps1 = str2double(source.String(source.Value));
elseif strcmp(source.Tag,'95') % Sweep sections
    a.SweepSections1 = str2double(source.String(source.Value));
elseif strcmp(source.Tag,'8') % Higher order modes
    a.HigherOrderModes1 = source.Value;
    a.Search1 = 0;
elseif strcmp(source.Tag,'9') % Symmetric modes/Torsional modes
    if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Circumferential')
        a.SymmetricModes1 = source.Value;
    else
        a.TorsionalModes1 = source.Value;
    end
elseif strcmp(source.Tag,'10') % Antisymmetric modes/Longitudinal modes
    if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Circumferential')
        a.AntisymmetricModes1 = source.Value;
    else
        a.LongitudinalModes1 = source.Value;
    end
elseif strcmp(source.Tag,'11') % Lamb modes/Flexural modes
    if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Circumferential')
        a.LambModes1 = source.Value;
    else
        a.FlexuralModes1 = source.Value;
        if  source.Value
            a.ShearHorizontalModes_FlexuralModeOrdersUI1.Enable = 'on';
        else
            a.ShearHorizontalModes_FlexuralModeOrdersUI1.Enable = 'off';
        end
    end
elseif strcmp(source.Tag,'12') % Shear horizontal modes/Flexural mode orders
    if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Circumferential')
        a.ShearHorizontalModes1 = source.Value;
    else
        a.FlexuralModeOrders1 = str2double(source.String);
        if  rem(a.FlexuralModeOrders1,1) ~= 0 || a.FlexuralModeOrders1 < 1
            a.FlexuralModeOrders1 = ceil(abs(str2double(source.String)));
            source.String = a.FlexuralModeOrders1;
            errordlg('The number of circumferential orders of flexural modes must be a positive integer.','Error');
        end
        a.Search1 = 0;
    end
elseif strcmp(source.Tag,'13') % Step
    a.Step1 = str2double(source.String);
    a.Search1 = 0;
elseif strcmp(source.Tag,'14') % Search
    FrequencySweeps
elseif strcmp(source.Tag,'15') % Trace modes/Stop trace
    if  source.Value
        if  strcmp(a.Geometry1,'Circumferential')
            if  a.Viscoelastic1
                warndlg('The material is viscoelastic. Material damping is ignored for circumferential waves.','Warning');
            end
            if  a.ThicknessInner1/a.Thickness1 >= .98
                helpdlg(['The circumference closely resembles a flat plate (inner diameter/outer diameter >= 0.98). Consider calculating a ',num2str((a.Thickness1-a.ThicknessInner1)/2),' mm thick plate.'],'Note');
            end
        end
        FrequencyRange = 0:a.FrequencyResolution1:a.FrequencyLimit1;
        FrequencyRangeF = FrequencyRange;
        FrequencyRange(1) = a.FrequencyRangeStart1;
        if  a.FrequencyRangeStart1 > .1*a.FrequencyResolution1
            FrequencyRange(1) = .1*a.FrequencyResolution1;
        end
        if  strcmp(a.Geometry1,'Pipe') && a.ToggleLowerFluid1 && ~a.Sink1
            FrequencyRangeF(1) = .5*a.FrequencyResolution1;
        else
            FrequencyRangeF(1) = .1*a.FrequencyResolution1;
        end
        if  ~a.Search1
            FrequencySweeps
        end
        if  a.Multithreading
            imshow('Multithreading22.png','Parent',a.a1UI1)
        end
        source.String = 'Stop trace';
        source.Tooltip = 'Stop tracing dispersion curves.';
        if  strcmp(a.Geometry1,'Plate')
            [a.SLamb1,a.ALamb1,a.BLamb1,a.SShear1,a.AShear1] = Computer_Isotropic_Plate(a.Multithreading,a.Viscoelastic1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.FrequencyLimit1,a.Material1,a.Thickness1/2e3,a.Symmetric1,a.PhaseVelocityStep1,a.FrequencyOffset1,a.LambPhaseVelocitySweepRange11,a.LambPhaseVelocitySweepRange21,a.Force2DTracing1,a.AttenuationLimit1,a.SearchWidth,a.SearchAreaSections,a.SearchAreaExtensions,a.Sweeps1,a.SweepSections1,FrequencyRange,FrequencyRangeF,a.HSLamb1,a.HSShear1,a.HALamb1,a.HAShear1,a.HBLamb1,a.FrequencyResolution1,a.PhaseVelocityLimit1,a.PhaseVelocityLimit21,a.Accuracy1,a.PhaseVelocitySections1,a.FrequencySections1,a.HigherOrderModes1,a.SymmetricModes1,a.AntisymmetricModes1,a.LambModes1,a.ShearHorizontalModes1,a.MissingSamples1,a.BelowCutoffWidth1,a.CriticalSlope,a.GearSlopeRanges,a.GearFrequencyResolution);
            a.F1={[]};a.L1={[]};a.T1={[]};a.CLamb1={[]};a.CShear1={[]};
        elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
            [a.F1,a.L1,a.T1] = Computer_Isotropic_Rod_Pipe(a.Geometry1,a.Multithreading,a.Viscoelastic1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,a.FrequencyLimit1,a.Material1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Symmetric1,a.PhaseVelocityStep1,a.FrequencyOffset1,a.LambPhaseVelocitySweepRange11,a.LambPhaseVelocitySweepRange21,a.Force2DTracing1,a.AttenuationLimit1,a.SearchWidth,a.SearchAreaSections,a.SearchAreaExtensions,a.Sweeps1,a.SweepSections1,FrequencyRange,FrequencyRangeF,a.HL1,a.HF1,a.HT1,a.FrequencyResolution1,a.PhaseVelocityLimit1,a.PhaseVelocityLimit21,a.Accuracy1,a.PhaseVelocitySections1,a.FrequencySections1,a.HigherOrderModes1,a.LongitudinalModes1,a.FlexuralModes1,a.TorsionalModes1,a.MissingSamples1,a.BelowCutoffWidth1,a.FlexuralModeOrders1,a.CriticalSlope,a.GearSlopeRanges,a.GearFrequencyResolution,a.LineColors1);
            a.SLamb1={[]};a.ALamb1={[]};a.BLamb1={[]};a.SShear1={[]};a.AShear1={[]};a.CLamb1={[]};a.CShear1={[]};
        elseif strcmp(a.Geometry1,'Circumferential')
            [a.CLamb1,a.CShear1] = Computer_Isotropic_Circumferential(a.Multithreading,a.FrequencyLimit1,a.Material1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.PhaseVelocityStep1,a.FrequencyOffset1,FrequencyRange,a.HCLamb1,a.HCShear1,a.FrequencyResolution1,a.PhaseVelocityLimit1,a.Accuracy1,a.PhaseVelocitySections1,a.FrequencySections1,a.HigherOrderModes1,a.LambModes1,a.ShearHorizontalModes1,a.MissingSamples1,a.BelowCutoffWidth1);
            a.SLamb1={[]};a.ALamb1={[]};a.BLamb1={[]};a.SShear1={[]};a.AShear1={[]};a.F1={[]};a.L1={[]};a.T1={[]};
        end
        source.Value = 0;
        source.String = 'Trace modes';
        source.Tooltip = 'Start tracing the dispersion curves.';
        if  a.Multithreading
            imshow('Multithreading21.png','Parent',a.a1UI1)
        end
        if  Stop
            return
        end
        if  (strcmp(a.Geometry1,'Plate') && ~a.Viscoelastic1 && ((a.LambModes1 && ~a.FluidLoading1) || ~a.LambModes1)) ||...
            (strcmp(a.Geometry1,'Rod') && ~a.Viscoelastic1 && ~a.FlexuralModes1 && ((a.LongitudinalModes1 && ~a.FluidLoading1) || ~a.LongitudinalModes1)) ||...
            (strcmp(a.Geometry1,'Pipe') && ~a.Viscoelastic1 && ~a.FlexuralModes1 && ~a.LongitudinalModes1)
            a.SamplesX3ceUI1.Enable = 'off';
            a.CalculateCeUI1.Enable = 'off';
        else
            a.SamplesX3ceUI1.Enable = 'on';
            a.CalculateCeUI1.Enable = 'on';
        end
        a.ModeNames1 = {''};
        if  strcmp(a.Geometry1,'Plate')
            if  ~isempty(a.ALamb1{1})
                for i = 0:length(a.ALamb1)-1
                    eval(sprintf('a.ModeNames1(i+1) = {''A%u''};',i));
                end
            end
            if  ~isempty(a.SLamb1{1})
                for i = 0:length(a.SLamb1)-1
                    eval(sprintf('a.ModeNames1(end+1) = {''S%u''};',i));
                end
            end
            if  ~isempty(a.BLamb1{1})
                for i = 0:length(a.BLamb1)-1
                    eval(sprintf('a.ModeNames1(i+1) = {''B%u''};',i));
                end
            end
            if  ~isempty(a.AShear1{1})
                for i = 1:length(a.AShear1)
                    eval(sprintf('a.ModeNames1(end+1) = {''ASH%u''};',i));
                end
            end
            if  ~isempty(a.SShear1{1})
                for i = 0:length(a.SShear1)-1
                    eval(sprintf('a.ModeNames1(end+1) = {''SSH%u''};',i));
                end
            end                
        elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
            if  ~isempty(a.T1{1})
                for i = 1:length(a.T1)
                    eval(sprintf('a.ModeNames1(i) = {''T(0,%u)''};',i));
                end
            end
            if  ~isempty(a.L1{1})
                for i = 1:length(a.L1)
                    eval(sprintf('a.ModeNames1(end+1) = {''L(0,%u)''};',i));
                end
            end
            if  ~isempty(a.F1{1})
                for n = 1:length(a.F1)
                    for i = 1:length(a.F1{n})
                        eval(sprintf('a.ModeNames1(end+1) = {''F(%u,%u)''};',n,i));
                    end
                end
            end
        elseif strcmp(a.Geometry1,'Circumferential')
            if  ~isempty(a.CLamb1{1})
                for i = 0:length(a.CLamb1)-1
                    eval(sprintf('a.ModeNames1(i+1) = {''C%u''};',i));
                end
            end
            if  ~isempty(a.CShear1{1})
                for i = 0:length(a.CShear1)-1
                    eval(sprintf('a.ModeNames1(end+1) = {''CSH%u''};',i));
                end
            end
        end
        for i = 1:length(a.ModeNames1)
            if  isempty(a.ModeNames1{i})
                z(i) = 0;
            else
                z(i) = 1;
            end
        end
        a.ModeNames1(z == 0) = [];
        a.Mode1UI1.Value = 1;
        a.Mode1UI1.String = a.ModeNames1;
        a.Mode11 = a.ModeNames1{1};
        a.Mode2UI1.Value = 1;
        a.Mode2UI1.String = a.ModeNames1;
        a.Mode21 = a.ModeNames1{1};
        a.ModeShapeSettingChanged1 = 1;

        if  strcmp(a.Geometry1,'Plate') % connection to Signal simulator
            a.Title3 = 1;
            a.TitleUI3.Style = 'checkbox';
            a.TitleUI3.Value = 1;
            a.TitleUI3.String = '';
            a.TitleUI3.Tooltip = 'Check this in order to show the plot title.';
            a.DataUI3.String = a.Material1.Name;
            a.PlotUI3.Enable = 'off';
            a.DataType3 = 1;
            a.Mode3 = '';
            a.Frequency3 = a.FrequencyResolution1*1e2;
            a.FrequencyUI3.String = a.Frequency3;
            [a.ALambModes1,a.SLambModes1,a.BLambModes1,a.AShearModes1,a.SShearModes1,a.BShearModes1] = ModeFinder(a.ALamb1,a.SLamb1,a.BLamb1,a.AShear1,a.SShear1,{0},a.Frequency3);
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
                if  a.Symmetric1
                    eval(sprintf('a.ALamb%uaUI3.String = [''A'',num2str(i)];',i));
                else
                    eval(sprintf('a.ALamb%uaUI3.String = [''B'',num2str(i)];',i));
                end
                eval(sprintf('a.SLamb%uaUI3.String = [''S'',num2str(i)];',i));
                eval(sprintf('a.SLamb%uaUI3.Position(3) = 16;',i));
                eval(sprintf('a.SLamb%ubUI3.Position(1) = 88;',i));
                eval(sprintf('a.SLamb%ucUI3.Position(1) = 103;',i));
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
            a.DataType3 = 1;
            a.DataUI3.String = '--';
            a.CalculateUI3.Enable = 'off';
            a.ALambModes1 = false;
            a.SLambModes1 = false;
            a.BLambModes1 = false;
            a.AShearModes1 = false;
            a.SShearModes1 = false;
            ShowModes_Signal(a)
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
        end

        a.CalculateCeUI1.Value = 1;
        a.CalculateCeUI1.String = 'Stop calculate';
        a.CalculateCeUI1.Tooltip = 'Stop calculating the energy velocity.';
        if  strcmp(a.Geometry1,'Plate') % calculate ce
            [a.SLamb1,a.ALamb1,a.BLamb1,a.SShear1,a.AShear1] = Computer_Isotropic_Plate_EnergyVelocity(a.SLamb1,a.ALamb1,a.BLamb1,a.SShear1,a.AShear1,a.Material1,a.Force2DTracing1,a.Viscoelastic1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Thickness1/2e3,a.SamplesX3ce1);
        elseif strcmp(a.Geometry1,'Rod')
            [a.F1,a.L1,a.T1] = Computer_Isotropic_Rod_EnergyVelocity(a.F1,a.L1,a.T1,a.Material1,a.Force2DTracing1,a.Viscoelastic1,a.FluidLoading1,a.UpperFluid1,a.Thickness1/2e3,a.SamplesX3ce1);
        elseif strcmp(a.Geometry1,'Pipe')
            [a.F1,a.L1,a.T1] = Computer_Isotropic_Pipe_EnergyVelocity(a.F1,a.L1,a.T1,a.Material1,a.Force2DTracing1,a.Viscoelastic1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.SamplesX3ce1);
        elseif strcmp(a.Geometry1,'Circumferential')
            [a.CLamb1,a.CShear1] = Computer_Isotropic_Circumferential_EnergyVelocity(a.CLamb1,a.CShear1,a.Material1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.SamplesX3ce1);
        end
        a.CalculateCeUI1.Value = 0;
        a.CalculateCeUI1.String = 'Calculate ce';
        a.CalculateCeUI1.Tooltip = 'Start calculating the energy velocity.';
    else
        Stop = 1;
        if  a.Multithreading
            source.Value = 0;
            source.String = 'Trace modes';
            source.Tooltip = 'Start tracing the dispersion curves.';
            imshow('Multithreading21.png','Parent',a.a1UI1)
            cancelAll(a.Pool.FevalQueue)
        end
    end
elseif strcmp(source.Tag,'91') % Samples x3 (ce)
    a.SamplesX3ce1 = str2double(source.String);
elseif strcmp(source.Tag,'16') % Calculate ce/Stop calculate
    if  source.Value && strcmp(source.String,'Calculate ce')
        source.String = 'Stop calculate';
        source.Tooltip = 'Stop calculating the energy velocity.';
        if  strcmp(a.Geometry1,'Plate') % calculate ce
            [a.SLamb1,a.ALamb1,a.BLamb1,a.SShear1,a.AShear1] = Computer_Isotropic_Plate_EnergyVelocity(a.SLamb1,a.ALamb1,a.BLamb1,a.SShear1,a.AShear1,a.Material1,a.Force2DTracing1,a.Viscoelastic1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Thickness1/2e3,a.SamplesX3ce1);
        elseif strcmp(a.Geometry1,'Rod')
            [a.F1,a.L1,a.T1] = Computer_Isotropic_Rod_EnergyVelocity(a.F1,a.L1,a.T1,a.Material1,a.Force2DTracing1,a.Viscoelastic1,a.FluidLoading1,a.UpperFluid1,a.Thickness1/2e3,a.SamplesX3ce1);
        elseif strcmp(a.Geometry1,'Pipe')
            [a.F1,a.L1,a.T1] = Computer_Isotropic_Pipe_EnergyVelocity(a.F1,a.L1,a.T1,a.Material1,a.Force2DTracing1,a.Viscoelastic1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.SamplesX3ce1);
        elseif strcmp(a.Geometry1,'Circumferential')
            [a.CLamb1,a.CShear1] = Computer_Isotropic_Circumferential_EnergyVelocity(a.CLamb1,a.CShear1,a.Material1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.SamplesX3ce1);
        end
        source.Value = 0;
        source.String = 'Calculate ce';
        source.Tooltip = 'Start calculating the energy velocity.';
    else
        Stop = 1;
    end
elseif strcmp(source.Tag,'17') % Quantity (Dispersion diagrams)
    switch source.Value
    case 1
        a.Quantity11 = 1;
        a.Option1UI1.Enable = 'on';
        a.Option1UI1.Style = 'checkbox';
        a.Option1UI1.String = '';
        a.Option1UI1.Value = a.BulkVelocities1;
        a.Option1UI1.Tooltip = 'Check this to show the bulk wave velocities.';
        a.Option1TextUI1.String = 'Bulk velocities';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeUI1.Tooltip = 'Select the frequency''s dimension on the X-axis.';        
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';
        end
        a.YAxisTextUI1.String = 'Y-axis (m/ms)';
        a.YAxisUI1.Tooltip = 'Enter which phase velocity range shall be plotted.';
    case 2
        a.Quantity11 = 2;
        a.Option1UI1.Enable = 'on';
        a.Option1UI1.Style = 'checkbox';
        a.Option1UI1.String = '';
        a.Option1UI1.Value = a.BulkVelocities1;
        a.Option1UI1.Tooltip = 'Check this to show the bulk wave velocities.';
        a.Option1TextUI1.String = 'Bulk velocities';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeUI1.Tooltip = 'Select the frequency''s dimension on the X-axis.';        
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';          
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';           
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';          
        end
        a.YAxisTextUI1.String = 'Y-axis (m/ms)';
        a.YAxisUI1.Tooltip = 'Enter which energy velocity range shall be plotted.';
    case 3
        a.Quantity11 = 3;
        a.Option1UI1.Enable = 'on';
        a.Option1UI1.Style = 'edit';
        a.Option1UI1.String = a.Distance1;
        a.Option1UI1.Tooltip = 'Enter the propagation distance from the source to the sensor for which the propagation time shall be calculated.';
        a.Option1TextUI1.String = 'Distance (mm)';
        a.XAxisModeTextUI1.String = 'Y-axis mode';
        a.XAxisModeUI1.Tooltip = 'Select the frequency''s dimension on the Y-axis.';        
        a.XAxisTextUI1.String = ['X-axis (',char(181),'s)'];
        a.XAxisUI1.Tooltip = 'Enter which propagation time range shall be plotted.';
        if  a.XAxisMode1 == 1
            a.YAxisTextUI1.String = 'Y-axis (kHz)';
            a.YAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';           
        elseif a.XAxisMode1 == 2
            a.YAxisTextUI1.String = 'Y-axis (MHz)';
            a.YAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';         
        elseif a.XAxisMode1 == 3
            a.YAxisTextUI1.String = ['Y-axis (MHz',char(8901),'mm)'];
            a.YAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';         
        end
    case 4
        a.Quantity11 = 4;
        a.Option1UI1.Enable = 'on';
        a.Option1UI1.Style = 'popupmenu';
        a.Option1UI1.String = fieldnames(a.Materials.Fluid);
        a.Option1UI1.Value = find(strcmp(a.Couplant1.Name,a.Option1UI1.String));
        a.Option1UI1.Tooltip = 'Select for which surrounding medium the coincidence angle shall be displayed.';
        a.Option1TextUI1.String = 'Couplant';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeUI1.Tooltip = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';           
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';          
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';          
        end
        a.YAxisTextUI1.String = ['Y-axis (',char(176),')'];
        a.YAxisUI1.Tooltip = 'Enter which coincidence angle range shall be plotted.';
    case 5
        a.Quantity11 = 5;
        a.Option1UI1.Enable = 'off';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeUI1.Tooltip = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';          
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';          
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';         
        end
        if  a.XAxisMode1 == 3
            a.YAxisTextUI1.String = 'Y-axis (mm/mm)';
            a.YAxisUI1.Tooltip = 'Enter which wavelength per thickness range shall be plotted.';           
        else
            a.YAxisTextUI1.String = 'Y-axis (mm)';
            a.YAxisUI1.Tooltip = 'Enter which wavelength range shall be plotted.';           
        end
    case 6
        a.Quantity11 = 6;
        a.Option1UI1.Enable = 'off';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeUI1.Tooltip = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';           
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';          
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';          
        end        
        if  a.XAxisMode1 == 3
            a.YAxisTextUI1.String = ['Y-axis (rad/mm',char(8901),'mm)'];
            a.YAxisUI1.Tooltip = 'Enter which wavenumber-thickness range shall be plotted.'; 
        else
            a.YAxisTextUI1.String = 'Y-axis (rad/mm)';
            a.YAxisUI1.Tooltip = 'Enter which wavenumber range shall be plotted.';
        end
    case 7
        a.Quantity11 = 7;
        a.Option1UI1.Enable = 'off';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeUI1.Tooltip = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';          
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';          
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';           
        end
        if  a.XAxisMode1 == 3
            a.YAxisTextUI1.String = ['Y-axis (Np/m',char(8901),'mm)'];
            a.YAxisUI1.Tooltip = 'Enter which attenuation-thickness range shall be plotted.'; 
        else
            a.YAxisTextUI1.String = 'Y-axis (Np/m)';
            a.YAxisUI1.Tooltip = 'Enter which attenuation range shall be plotted.';
        end
    end
    AdjustAxes
elseif strcmp(source.Tag,'18') % option (Dispersion diagrams)
    if  a.Quantity11 < 3
        a.BulkVelocities1 = source.Value;
    elseif a.Quantity11 == 3
        a.Distance1 = str2double(source.String);
    elseif a.Quantity11 == 4
        E = fieldnames(a.Materials.Fluid);
        a.Couplant1 = getfield(a.Materials.Fluid,E{source.Value});
    end
    if  a.Quantity11 == 3 || a.Quantity11 == 4
        AdjustAxes
    end
elseif strcmp(source.Tag,'19') % X-Axis mode
    switch source.Value
    case 1
        if  a.Quantity11 ~= 3
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';
            val = eval(a.XAxisUI1.String);
            if  a.XAxisMode1 == 2
                a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
                a.XAxisUI1.String = ['[',num2str(val(1)*1e3),' ',num2str(val(2)*1e3),']'];
            elseif a.XAxisMode1 == 3
                if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                    a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness1;
                    a.XAxisUI1.String = ['[',num2str(val(1)*1e3/a.Thickness1),' ',num2str(val(2)*1e3/a.Thickness1),']'];
                elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                    a.XAxis1 = eval(a.XAxisUI1.String)*2e3/(a.Thickness1-a.ThicknessInner1);
                    a.XAxisUI1.String = ['[',num2str(val(1)*2e3/(a.Thickness1-a.ThicknessInner1)),' ',num2str(val(2)*2e3/(a.Thickness1-a.ThicknessInner1)),']'];
                end
            end
        else
            a.YAxisTextUI1.String = 'Y-axis (kHz)';
            a.YAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';
            val = eval(a.YAxisUI1.String);
            if  a.XAxisMode1 == 2
                a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
                a.YAxisUI1.String = ['[',num2str(val(1)*1e3),' ',num2str(val(2)*1e3),']'];
            elseif a.XAxisMode1 == 3
                if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                    a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness1;
                    a.YAxisUI1.String = ['[',num2str(val(1)*1e3/a.Thickness1),' ',num2str(val(2)*1e3/a.Thickness1),']'];
                elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                    a.YAxis1 = eval(a.YAxisUI1.String)*2e3/(a.Thickness1-a.ThicknessInner1);
                    a.YAxisUI1.String = ['[',num2str(val(1)*2e3/(a.Thickness1-a.ThicknessInner1)),' ',num2str(val(2)*2e3/(a.Thickness1-a.ThicknessInner1)),']'];
                end
            end
        end
        if  (a.Quantity11 == 5 || a.Quantity11 == 6 || a.Quantity11 == 7) && a.XAxisMode1 == 3
            if  a.Quantity11 == 5
                a.YAxisTextUI1.String = 'Y-axis (mm)';
                a.YAxisUI1.Tooltip = 'Enter which wavelength range shall be plotted.';
            elseif a.Quantity11 == 6
                a.YAxisTextUI1.String = 'Y-axis (rad/mm)';
                a.YAxisUI1.Tooltip = 'Enter which wavenumber range shall be plotted.';
            elseif a.Quantity11 == 7
                a.YAxisTextUI1.String = 'Y-axis (Np/m)';
                a.YAxisUI1.Tooltip = 'Enter which attenuation range shall be plotted.';
            end
            val = eval(a.YAxisUI1.String);
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                if  a.Quantity11 == 5
                    a.YAxis1 = eval(a.YAxisUI1.String)*a.Thickness1;
                    a.YAxisUI1.String = ['[',num2str(val(1)*a.Thickness1),' ',num2str(val(2)*a.Thickness1),']'];
                elseif a.Quantity11 == 6 || a.Quantity11 == 7
                    a.YAxis1 = eval(a.YAxisUI1.String)/a.Thickness1;
                    a.YAxisUI1.String = ['[',num2str(val(1)/a.Thickness1),' ',num2str(val(2)/a.Thickness1),']'];
                end
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                if  a.Quantity11 == 5
                    a.YAxis1 = eval(a.YAxisUI1.String)*(a.Thickness1-a.ThicknessInner1)/2;
                    a.YAxisUI1.String = ['[',num2str(val(1)*(a.Thickness1-a.ThicknessInner1)/2),' ',num2str(val(2)*(a.Thickness1-a.ThicknessInner1)/2),']'];
                elseif a.Quantity11 == 6 || a.Quantity11 == 7
                    a.YAxis1 = eval(a.YAxisUI1.String)/(a.Thickness1-a.ThicknessInner1)*2;
                    a.YAxisUI1.String = ['[',num2str(val(1)/(a.Thickness1-a.ThicknessInner1)*2),' ',num2str(val(2)/(a.Thickness1-a.ThicknessInner1)*2),']'];
                end
            end
        end
        a.XAxisMode1 = 1;
    case 2
        if  a.Quantity11 ~= 3
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';
            val = eval(a.XAxisUI1.String);
            if  a.XAxisMode1 == 1
                a.XAxis1 = eval(a.XAxisUI1.String);
                a.XAxisUI1.String = ['[',num2str(val(1)/1e3),' ',num2str(val(2)/1e3),']'];
            elseif a.XAxisMode1 == 3
                if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                    a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness1;
                    a.XAxisUI1.String = ['[',num2str(val(1)/a.Thickness1),' ',num2str(val(2)/a.Thickness1),']'];
                elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                    a.XAxis1 = eval(a.XAxisUI1.String)*2e3/(a.Thickness1-a.ThicknessInner1);
                    a.XAxisUI1.String = ['[',num2str(val(1)/(a.Thickness1-a.ThicknessInner1)*2),' ',num2str(val(2)/(a.Thickness1-a.ThicknessInner1)*2),']'];
                end
            end
        else
            a.YAxisTextUI1.String = 'Y-axis (MHz)';
            a.YAxisUI1.Tooltip = 'Enter which frequency range shall be plotted.';
            val = eval(a.YAxisUI1.String);
            if  a.XAxisMode1 == 1
                a.YAxis1 = eval(a.YAxisUI1.String);
                a.YAxisUI1.String = ['[',num2str(val(1)/1e3),' ',num2str(val(2)/1e3),']'];
            elseif a.XAxisMode1 == 3
                if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                    a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness1;
                    a.YAxisUI1.String = ['[',num2str(val(1)/a.Thickness1),' ',num2str(val(2)/a.Thickness1),']'];
                elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                    a.YAxis1 = eval(a.YAxisUI1.String)*2e3/(a.Thickness1-a.ThicknessInner1);
                    a.YAxisUI1.String = ['[',num2str(val(1)/(a.Thickness1-a.ThicknessInner1)*2),' ',num2str(val(2)/(a.Thickness1-a.ThicknessInner1)*2),']'];
                end
            end
        end
        if  (a.Quantity11 == 5 || a.Quantity11 == 6 || a.Quantity11 == 7) && a.XAxisMode1 == 3
            if  a.Quantity11 == 5
                a.YAxisTextUI1.String = 'Y-axis (mm)';
                a.YAxisUI1.Tooltip = 'Enter which wavelength range shall be plotted.';
            elseif a.Quantity11 == 6
                a.YAxisTextUI1.String = 'Y-axis (rad/mm)';
                a.YAxisUI1.Tooltip = 'Enter which wavenumber range shall be plotted.';
            elseif a.Quantity11 == 7
                a.YAxisTextUI1.String = 'Y-axis (Np/m)';
                a.YAxisUI1.Tooltip = 'Enter which attenuation range shall be plotted.';
            end
            val = eval(a.YAxisUI1.String);
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                if  a.Quantity11 == 5
                    a.YAxis1 = eval(a.YAxisUI1.String)*a.Thickness1;
                    a.YAxisUI1.String = ['[',num2str(val(1)*a.Thickness1),' ',num2str(val(2)*a.Thickness1),']'];
                elseif a.Quantity11 == 6 || a.Quantity11 == 7
                    a.YAxis1 = eval(a.YAxisUI1.String)/a.Thickness1;
                    a.YAxisUI1.String = ['[',num2str(val(1)/a.Thickness1),' ',num2str(val(2)/a.Thickness1),']'];
                end
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                if  a.Quantity11 == 5
                    a.YAxis1 = eval(a.YAxisUI1.String)*(a.Thickness1-a.ThicknessInner1)/2;
                    a.YAxisUI1.String = ['[',num2str(val(1)*(a.Thickness1-a.ThicknessInner1)/2),' ',num2str(val(2)*(a.Thickness1-a.ThicknessInner1)/2),']'];
                elseif a.Quantity11 == 6 || a.Quantity11 == 7
                    a.YAxis1 = eval(a.YAxisUI1.String)/(a.Thickness1-a.ThicknessInner1)*2;
                    a.YAxisUI1.String = ['[',num2str(val(1)/(a.Thickness1-a.ThicknessInner1)*2),' ',num2str(val(2)/(a.Thickness1-a.ThicknessInner1)*2),']'];
                end
            end
        end
        a.XAxisMode1 = 2;
    case 3
        if  a.Quantity11 ~= 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';
            val = eval(a.XAxisUI1.String);
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                if  a.XAxisMode1 == 1
                    a.XAxis1 = eval(a.XAxisUI1.String);
                    a.XAxisUI1.String = ['[',num2str(val(1)/1e3*a.Thickness1),' ',num2str(val(2)/1e3*a.Thickness1),']'];
                elseif a.XAxisMode1 == 2
                    a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
                    a.XAxisUI1.String = ['[',num2str(val(1)*a.Thickness1),' ',num2str(val(2)*a.Thickness1),']'];
                end
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                if  a.XAxisMode1 == 1
                    a.XAxis1 = eval(a.XAxisUI1.String);
                    a.XAxisUI1.String = ['[',num2str(val(1)/2e3*(a.Thickness1-a.ThicknessInner1)),' ',num2str(val(2)/2e3*(a.Thickness1-a.ThicknessInner1)),']'];
                elseif a.XAxisMode1 == 2
                    a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
                    a.XAxisUI1.String = ['[',num2str(val(1)*(a.Thickness1-a.ThicknessInner1)/2),' ',num2str(val(2)*(a.Thickness1-a.ThicknessInner1)/2),']'];
                end
            end
        else
            a.YAxisTextUI1.String = ['Y-axis (MHz',char(8901),'mm)'];
            a.YAxisUI1.Tooltip = 'Enter which frequency-thickness range shall be plotted.';
            val = eval(a.YAxisUI1.String);
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                if  a.XAxisMode1 == 1
                    a.YAxis1 = eval(a.YAxisUI1.String);
                    a.YAxisUI1.String = ['[',num2str(val(1)/1e3*a.Thickness1),' ',num2str(val(2)/1e3*a.Thickness1),']'];
                elseif a.XAxisMode1 == 2
                    a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
                    a.YAxisUI1.String = ['[',num2str(val(1)*a.Thickness1),' ',num2str(val(2)*a.Thickness1),']'];
                end
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                if  a.XAxisMode1 == 1
                    a.YAxis1 = eval(a.YAxisUI1.String);
                    a.YAxisUI1.String = ['[',num2str(val(1)/2e3*(a.Thickness1-a.ThicknessInner1)),' ',num2str(val(2)/2e3*(a.Thickness1-a.ThicknessInner1)),']'];
                elseif a.XAxisMode1 == 2
                    a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
                    a.YAxisUI1.String = ['[',num2str(val(1)*(a.Thickness1-a.ThicknessInner1)/2),' ',num2str(val(2)*(a.Thickness1-a.ThicknessInner1)/2),']'];
                end
            end
        end
        if  (a.Quantity11 == 5 || a.Quantity11 == 6 || a.Quantity11 == 7) && a.XAxisMode1 ~= 3 
            if  a.Quantity11 == 5
                a.YAxisTextUI1.String = 'Y-axis (mm/mm)';
                a.YAxisUI1.Tooltip = 'Enter which wavelength per thickness range shall be plotted.';
            elseif a.Quantity11 == 6
                a.YAxisTextUI1.String = ['Y-axis (rad/mm',char(8901),'mm)'];
                a.YAxisUI1.Tooltip = 'Enter which wavenumber-thickness range shall be plotted.';
            elseif a.Quantity11 == 7
                a.YAxisTextUI1.String = ['Y-axis (Np/m',char(8901),'mm)'];
                a.YAxisUI1.Tooltip = 'Enter which attenuation-thickness range shall be plotted.'; 
            end
            val = eval(a.YAxisUI1.String);
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                if  a.Quantity11 == 5
                    a.YAxis1 = eval(a.YAxisUI1.String)/a.Thickness1;
                    a.YAxisUI1.String = ['[',num2str(val(1)/a.Thickness1),' ',num2str(val(2)/a.Thickness1),']'];
                elseif a.Quantity11 == 6 || a.Quantity11 == 7
                    a.YAxis1 = eval(a.YAxisUI1.String)*a.Thickness1;
                    a.YAxisUI1.String = ['[',num2str(val(1)*a.Thickness1),' ',num2str(val(2)*a.Thickness1),']'];
                end
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                if  a.Quantity11 == 5
                    a.YAxis1 = eval(a.YAxisUI1.String)/(a.Thickness1-a.ThicknessInner1)*2;
                    a.YAxisUI1.String = ['[',num2str(val(1)/(a.Thickness1-a.ThicknessInner1)*2),' ',num2str(val(2)/(a.Thickness1-a.ThicknessInner1)*2),']'];
                elseif a.Quantity11 == 6 || a.Quantity11 == 7
                    a.YAxis1 = eval(a.YAxisUI1.String)*(a.Thickness1-a.ThicknessInner1)/2;
                    a.YAxisUI1.String = ['[',num2str(val(1)*(a.Thickness1-a.ThicknessInner1)/2),' ',num2str(val(2)*(a.Thickness1-a.ThicknessInner1)/2),']'];
                end
            end
        end
        a.XAxisMode1 = 3;
    end
elseif strcmp(source.Tag,'20') % X-Axis
    a.XAxis1 = eval(source.String);
    if  a.Quantity11 ~= 3
        if  a.XAxisMode1 == 1
            a.XAxis1 = eval(source.String);
        elseif a.XAxisMode1 == 2
            a.XAxis1 = eval(source.String)*1e3;
        elseif a.XAxisMode1 == 3
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                a.XAxis1 = eval(source.String)*1e3/a.Thickness1;
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                a.XAxis1 = eval(source.String)*2e3/(a.Thickness1-a.ThicknessInner1);
            end
        end
    end
elseif strcmp(source.Tag,'21') % Y-Axis
    a.YAxis1 = eval(source.String);
    if  a.Quantity11 == 3
        if  a.XAxisMode1 == 1
            a.YAxis1 = eval(source.String);
        elseif a.XAxisMode1 == 2
            a.YAxis1 = eval(source.String)*1e3;
        elseif a.XAxisMode1 == 3
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                a.YAxis1 = eval(source.String)*1e3/a.Thickness1;
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                a.YAxis1 = eval(source.String)*2e3/(a.Thickness1-a.ThicknessInner1);
            end
        end
    end
elseif strcmp(source.Tag,'22') % Plot (Dispersion diagrams)
    try
        if  a.Quantity11 == 1
            DispersionDiagram_PhaseVelocity(a.Geometry1,0,0,'',a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.BulkVelocities1,a.PNGresolution1,a.SColor1,a.AColor1,a.BColor1,a.ALamb1,a.AntisymmetricModes1,a.AShear1,a.BLamb1,a.BoxLineWidth1,{[]},a.CLamb1,a.CShear1,a.FlexuralModes1,a.FlexuralModeOrders1,a.LongitudinalModes1,a.TorsionalModes1,a.F1,a.L1,a.T1,a.LineColors1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.Title1,a.LambModes1,a.LineWidth1,a.Material1,a.PDF1,a.FileName1,a.Thickness1/1e3,a.ThicknessInner1/1e3,a.PNG1,0,a.ShearHorizontalModes1,a.SLamb1,a.SShear1,1,1,a.SymmetricModes1,0,a.Symmetric1,a.XAxis1,a.XAxisMode1,a.YAxis1,1)
        elseif a.Quantity11 == 2
            DispersionDiagram_EnergyVelocity(1,a.Geometry1,0,0,'',a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.BulkVelocities1,a.PNGresolution1,a.SColor1,a.AColor1,a.BColor1,a.ALamb1,a.AntisymmetricModes1,a.AShear1,a.BLamb1,a.BoxLineWidth1,{[]},a.CLamb1,a.CShear1,a.FlexuralModes1,a.FlexuralModeOrders1,a.LongitudinalModes1,a.TorsionalModes1,a.F1,a.L1,a.T1,a.LineColors1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.Title1,a.LambModes1,a.LineWidth1,a.Material1,a.PDF1,a.FileName1,a.Thickness1/1e3,a.ThicknessInner1/1e3,a.PNG1,0,a.ShearHorizontalModes1,a.SLamb1,a.SShear1,1,1,a.SymmetricModes1,0,a.Symmetric1,a.XAxis1,a.XAxisMode1,a.YAxis1,1)
        elseif a.Quantity11 == 3
            DispersionDiagram_PropagationTime(a.Geometry1,0,'',a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Distance1,a.PNGresolution1,a.SColor1,a.AColor1,a.BColor1,a.ALamb1,a.AntisymmetricModes1,a.AShear1,a.BLamb1,a.BoxLineWidth1,{[]},a.CLamb1,a.CShear1,a.FlexuralModes1,a.FlexuralModeOrders1,a.LongitudinalModes1,a.TorsionalModes1,a.F1,a.L1,a.T1,a.LineColors1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.Title1,a.LambModes1,a.LineWidth1,a.Material1,a.PDF1,a.FileName1,a.Thickness1/1e3,a.ThicknessInner1/1e3,a.PNG1,0,a.ShearHorizontalModes1,a.SLamb1,a.SShear1,1,1,a.SymmetricModes1,0,a.Symmetric1,a.XAxis1,a.XAxisMode1,a.YAxis1,1)
        elseif a.Quantity11 == 4
            DispersionDiagram_CoincidenceAngle(a.Geometry1,0,'',a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Couplant1,a.PNGresolution1,a.SColor1,a.AColor1,a.BColor1,a.ALamb1,a.AntisymmetricModes1,a.AShear1,a.BLamb1,a.BoxLineWidth1,{[]},a.CLamb1,a.CShear1,a.FlexuralModes1,a.FlexuralModeOrders1,a.LongitudinalModes1,a.TorsionalModes1,a.F1,a.L1,a.T1,a.LineColors1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.Title1,a.LambModes1,a.LineWidth1,a.Material1,a.PDF1,a.FileName1,a.Thickness1/1e3,a.ThicknessInner1/1e3,a.PNG1,0,a.ShearHorizontalModes1,a.SLamb1,a.SShear1,1,1,a.SymmetricModes1,0,a.Symmetric1,a.XAxis1,a.XAxisMode1,a.YAxis1,1)
        elseif a.Quantity11 == 5
            DispersionDiagram_Wavelength(a.Geometry1,0,'',a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.PNGresolution1,a.SColor1,a.AColor1,a.BColor1,a.ALamb1,a.AntisymmetricModes1,a.AShear1,a.BLamb1,a.BoxLineWidth1,{[]},a.CLamb1,a.CShear1,a.FlexuralModes1,a.FlexuralModeOrders1,a.LongitudinalModes1,a.TorsionalModes1,a.F1,a.L1,a.T1,a.LineColors1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.Title1,a.LambModes1,a.LineWidth1,a.Material1,a.PDF1,a.FileName1,a.Thickness1/1e3,a.ThicknessInner1/1e3,a.PNG1,0,a.ShearHorizontalModes1,a.SLamb1,a.SShear1,1,1,a.SymmetricModes1,0,a.Symmetric1,a.XAxis1,a.XAxisMode1,a.YAxis1,1)
        elseif a.Quantity11 == 6
            DispersionDiagram_Wavenumber(a.Geometry1,0,'',a.Viscoelastic1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.PNGresolution1,a.SColor1,a.AColor1,a.BColor1,a.ALamb1,a.AntisymmetricModes1,a.AShear1,a.BLamb1,a.BoxLineWidth1,{[]},a.CLamb1,a.CShear1,a.FlexuralModes1,a.FlexuralModeOrders1,a.LongitudinalModes1,a.TorsionalModes1,a.F1,a.L1,a.T1,a.LineColors1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.Title1,a.LambModes1,a.LineWidth1,a.Material1,a.PDF1,a.FileName1,a.Thickness1/1e3,a.ThicknessInner1/1e3,a.PNG1,0,a.ShearHorizontalModes1,a.SLamb1,a.SShear1,1,1,a.SymmetricModes1,0,a.Symmetric1,a.XAxis1,a.XAxisMode1,a.YAxis1,1)
        elseif a.Quantity11 == 7
            DispersionDiagram_Attenuation(a.Geometry1,0,'',a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.PNGresolution1,a.SColor1,a.AColor1,a.BColor1,a.ALamb1,a.AntisymmetricModes1,a.AShear1,a.BLamb1,a.BoxLineWidth1,{[]},a.CLamb1,a.CShear1,a.FlexuralModes1,a.FlexuralModeOrders1,a.LongitudinalModes1,a.TorsionalModes1,a.F1,a.L1,a.T1,a.LineColors1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.Title1,a.LambModes1,a.LineWidth1,a.Material1,a.PDF1,a.FileName1,a.Thickness1/1e3,a.ThicknessInner1/1e3,a.PNG1,0,a.ShearHorizontalModes1,a.SLamb1,a.SShear1,1,1,a.SymmetricModes1,0,a.Symmetric1,a.XAxis1,a.XAxisMode1,a.YAxis1,1)
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
elseif strcmp(source.Tag,'23') % Quantity (Through thickness profiles)
    switch source.Value
    case 1
        a.Quantity21 = 1;
        a.PhaseUI1.Enable = 'on';
        a.x11UI1.Enable = 'off';
        a.x12UI1.Enable = 'off';
        a.x13UI1.Enable = 'on';
        a.x22UI1.Enable = 'off';
        a.x23UI1.Enable = 'on';
        a.x33UI1.Enable = 'on';
        a.x13UI1.String = '1';
        a.x23UI1.String = '2';
        a.x33UI1.String = '3';
        a.x13UI1.Tooltip = 'Displacement u1';
        a.x23UI1.Tooltip = 'Displacement u2';
        a.x33UI1.Tooltip = 'Displacement u3';
    case 2
        a.Quantity21 = 2;
        a.PhaseUI1.Enable = 'on';
        a.x11UI1.Enable = 'on';
        a.x12UI1.Enable = 'on';
        a.x13UI1.Enable = 'on';
        a.x22UI1.Enable = 'on';
        a.x23UI1.Enable = 'on';
        a.x33UI1.Enable = 'on';       
        a.x13UI1.String = '13';
        a.x23UI1.String = '23';
        a.x33UI1.String = '33';
        a.x11UI1.Tooltip = ['Stress ',char(963),'11'];
        a.x12UI1.Tooltip = ['Stress ',char(963),'12'];
        a.x13UI1.Tooltip = ['Stress ',char(963),'13'];
        a.x22UI1.Tooltip = ['Stress ',char(963),'22'];
        a.x23UI1.Tooltip = ['Stress ',char(963),'23'];
        a.x33UI1.Tooltip = ['Stress ',char(963),'33'];
    case 3
        a.Quantity21 = 3;
        a.PhaseUI1.Enable = 'on';
        a.x11UI1.Enable = 'on';
        a.x12UI1.Enable = 'on';
        a.x13UI1.Enable = 'on';
        a.x22UI1.Enable = 'on';
        a.x23UI1.Enable = 'on';
        a.x33UI1.Enable = 'on';       
        a.x13UI1.String = '13';
        a.x23UI1.String = '23';
        a.x33UI1.String = '33';
        a.x11UI1.Tooltip = ['Strain ',char(949),'11'];
        a.x12UI1.Tooltip = ['Strain ',char(949),'12'];
        a.x13UI1.Tooltip = ['Strain ',char(949),'13'];
        a.x22UI1.Tooltip = ['Strain ',char(949),'22'];
        a.x23UI1.Tooltip = ['Strain ',char(949),'23'];
        a.x33UI1.Tooltip = ['Strain ',char(949),'33'];
    case 4
        a.Quantity21 = 4;
        a.PhaseUI1.Enable = 'off';
        a.x11UI1.Enable = 'off';
        a.x12UI1.Enable = 'off';
        a.x13UI1.Enable = 'on';
        a.x22UI1.Enable = 'off';
        a.x23UI1.Enable = 'on';
        a.x33UI1.Enable = 'on';
        a.x13UI1.String = 'Estn';
        a.x23UI1.String = 'Ekin';
        a.x33UI1.String = 'Etot';
        a.x13UI1.Tooltip = 'Strain energy density';
        a.x23UI1.Tooltip = 'Kinetic energy density';
        a.x33UI1.Tooltip = 'Total energy density';        
    case 5
        a.Quantity21 = 5;
        a.PhaseUI1.Enable = 'off';
        a.x11UI1.Enable = 'off';
        a.x12UI1.Enable = 'off';
        a.x13UI1.Enable = 'on';
        a.x22UI1.Enable = 'off';
        a.x23UI1.Enable = 'on';
        a.x33UI1.Enable = 'on';
        a.x13UI1.String = '1';
        a.x23UI1.String = '2';
        a.x33UI1.String = '3';
        a.x13UI1.Tooltip = 'Power flow density p1';
        a.x23UI1.Tooltip = 'Power flow density p2';
        a.x33UI1.Tooltip = 'Power flow density p3';                
    end
elseif strcmp(source.Tag,'24') % Mode (Through thickness profiles)
    a.Mode11 = a.ModeNames1{source.Value};    
elseif strcmp(source.Tag,'25') % Frequency (Through thickness profiles)
    a.Frequency11 = str2double(source.String);
elseif strcmp(source.Tag,'26') % Samples x3 (Through thickness profiles)
    a.Samples11 = str2double(source.String);
elseif strcmp(source.Tag,'80') % Half-spaces (Through thickness profiles)
    a.HalfspacesNumber11 = str2double(source.String);
elseif strcmp(source.Tag,'81') % Half-spaces (Through thickness profiles)
    a.Halfspaces11 = source.Value;
    if  source.Value
        a.HalfspacesNumber1UI1.Enable = 'on';
    else
        a.HalfspacesNumber1UI1.Enable = 'off';
    end
elseif strcmp(source.Tag,'84') % Phase
    a.Phase1 = source.Value;
elseif strcmp(source.Tag,'28') % 11
    a.Plot1(1) = source.Value;
elseif strcmp(source.Tag,'29') % 22
    a.Plot1(2) = source.Value;
elseif strcmp(source.Tag,'30') % 33
    a.Plot1(3) = source.Value;
elseif strcmp(source.Tag,'32') % 23
    a.Plot1(4) = source.Value;
elseif strcmp(source.Tag,'31') % 13
    a.Plot1(5) = source.Value;
elseif strcmp(source.Tag,'33') % 12
    a.Plot1(6) = source.Value;
elseif strcmp(source.Tag,'27') % Plot (Through thickness profiles)
    try
        if  a.Quantity21 == 1
            u1Color = [1 0 0];
            u2Color = [.13 .55 .13];
            u3Color = [0 0 1];
            if  strcmp(a.Geometry1,'Plate')
                ModeShapeLines(1,1,0,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,u1Color,u2Color,u3Color,0,0,0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
            elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                ModeShapeLines_Isotropic_Rod_Pipe(1,1,0,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,u1Color,u2Color,u3Color,0,0,0,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
            elseif strcmp(a.Geometry1,'Circumferential')
                ModeShapeLines_Isotropic_Circumferential(1,1,0,a.Plot1,a.PNGresolution1,a.Material1,u1Color,u2Color,u3Color,0,0,0,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Phase1)
            end
        elseif a.Quantity21 == 2
            sigma33Color = [0 0 1];
            sigma11Color = [1 0 0];
            sigma22Color = [.13 .55 .13];
            sigma13Color = [0 0 0];
            sigma23Color = [1 0 1];
            sigma12Color = [0 1 1];
            if  strcmp(a.Geometry1,'Plate')
                ModeShapeLines(2,1,0,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,sigma11Color,sigma22Color,sigma33Color,sigma23Color,sigma13Color,sigma12Color,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
            elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                ModeShapeLines_Isotropic_Rod_Pipe(2,1,0,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,sigma11Color,sigma22Color,sigma33Color,sigma23Color,sigma13Color,sigma12Color,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
            elseif strcmp(a.Geometry1,'Circumferential')
                ModeShapeLines_Isotropic_Circumferential(2,1,0,a.Plot1,a.PNGresolution1,a.Material1,sigma11Color,sigma22Color,sigma33Color,sigma23Color,sigma13Color,sigma12Color,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Phase1)
            end
        elseif a.Quantity21 == 3
            epsilon33Color = [0 0 1];
            epsilon11Color = [1 0 0];
            epsilon22Color = [.13 .55 .13];
            epsilon13Color = [0 0 0];
            epsilon23Color = [1 0 1];
            epsilon12Color = [0 1 1];
            if  strcmp(a.Geometry1,'Plate')
                ModeShapeLines(3,1,0,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,epsilon11Color,epsilon22Color,epsilon33Color,epsilon23Color,epsilon13Color,epsilon12Color,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
            elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                ModeShapeLines_Isotropic_Rod_Pipe(3,1,0,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,epsilon11Color,epsilon22Color,epsilon33Color,epsilon23Color,epsilon13Color,epsilon12Color,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
            elseif strcmp(a.Geometry1,'Circumferential')
                ModeShapeLines_Isotropic_Circumferential(3,1,0,a.Plot1,a.PNGresolution1,a.Material1,epsilon11Color,epsilon22Color,epsilon33Color,epsilon23Color,epsilon13Color,epsilon12Color,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Phase1)
            end
        elseif a.Quantity21 == 4
            StrainEnergyDensityColor = [0 0 1];
            KineticEnergyDensityColor = [1 0 0];
            TotalEnergyDensityColor = [0 0 0];
            if  strcmp(a.Geometry1,'Plate')
                ModeShapeLines(4,1,0,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,StrainEnergyDensityColor,KineticEnergyDensityColor,TotalEnergyDensityColor,0,0,0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,0)
            elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                ModeShapeLines_Isotropic_Rod_Pipe(4,1,0,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,StrainEnergyDensityColor,KineticEnergyDensityColor,TotalEnergyDensityColor,0,0,0,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
            elseif strcmp(a.Geometry1,'Circumferential')
                ModeShapeLines_Isotropic_Circumferential(4,1,0,a.Plot1,a.PNGresolution1,a.Material1,StrainEnergyDensityColor,KineticEnergyDensityColor,TotalEnergyDensityColor,0,0,0,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,0)
            end
        elseif a.Quantity21 == 5
            P1Color = [1 0 0];
            P2Color = [.13 .55 .13];
            P3Color = [0 0 1];
            if  strcmp(a.Geometry1,'Plate')
                ModeShapeLines(5,1,0,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,P1Color,P2Color,P3Color,0,0,0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,0)
            elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                ModeShapeLines_Isotropic_Rod_Pipe(5,1,0,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,P1Color,P2Color,P3Color,0,0,0,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
            elseif strcmp(a.Geometry1,'Circumferential')
                ModeShapeLines_Isotropic_Circumferential(5,1,0,a.Plot1,a.PNGresolution1,a.Material1,P1Color,P2Color,P3Color,0,0,0,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,0)
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
elseif strcmp(source.Tag,'92') % Live-plot (Through thickness profiles)
    try
        ModeShapeLive(a.Geometry1,0,'',0,a.FrequencyLimit1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,a.ALamb1,a.AntisymmetricModes1,a.AShear1,a.BLamb1,{[]},a.CLamb1,a.CShear1,a.FlexuralModes1,a.FlexuralModeOrders1,a.LongitudinalModes1,a.TorsionalModes1,a.F1,a.L1,a.T1,a.LineColors1,a.LambModes1,{a.Material1},a.Thickness1/1e3,a.ThicknessInner1/1e3,a.ShearHorizontalModes1,a.SLamb1,a.SShear1,1,1,a.SymmetricModes1,0,1,{conj(a.Material1.C)},0,1,0,[1 -1;-1 1],[],a.Thickness1/1e3)
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
elseif strcmp(source.Tag,'34') % Mode (Mode shape)
    a.Mode21 = a.ModeNames1{source.Value};
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'43') % Plane
    switch source.Value 
    case 1
        a.Plane1 = 'r-z';
    case 2
        a.Plane1 = 'r-theta';
    end
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'35') % Frequency (mode shape)
    a.Frequency21 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'36') % Wavelengths
    a.Length1 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'37') % Samples x1
    a.Samples21 = str2double(source.String);
    a.GridLine1 = 2*round(a.Samples21/80);
    if  a.GridLine1 == 0
        a.GridLine1 = 1;
    end
    a.GridLineUI1.String = a.GridLine1; 
    a.ModeShapeSettingChanged1 = 1;
    AdjustSamples31
elseif strcmp(source.Tag,'38') % Samples x3 (mode shape)
    a.Samples31 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'39') % Gain
    a.Gain1 = str2double(source.String);
elseif strcmp(source.Tag,'40') % Grid line
    a.GridLine1 = str2double(source.String);
elseif strcmp(source.Tag,'41') % Undistorted
    a.Undistorted1 = source.Value;
elseif strcmp(source.Tag,'82') || strcmp(source.Tag,'83') % Half-spaces (mode shape)
    if  strcmp(source.Tag,'82')
        a.HalfspacesNumber21 = str2double(source.String);
    elseif strcmp(source.Tag,'83')
        a.Halfspaces21 = source.Value;
        if  source.Value
            a.HalfspacesNumber2UI1.Enable = 'on';
        else
            a.HalfspacesNumber2UI1.Enable = 'off';
        end
    end
    a.ModeShapeSettingChanged1 = 1;
    AdjustSamples31
elseif strcmp(source.Tag,'44') % Cycle duration
    a.CycleDuration1 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'45') % Frame rate
    a.FrameRate1 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'46') % Movie quality
    a.MovieQuality1 = str2double(source.String);
elseif strcmp(source.Tag,'47') % Animate
    a.Animate1 = source.Value;
elseif strcmp(source.Tag,'42') % Plot (mode shape)
    try
        if  ~a.Animate1
            if  strcmp(a.Geometry1,'Plate')
                ModeShapeGrid(0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',{a.Material1},1,{conj(a.Material1.C)},0,a.PNGresolution1,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.Directory1,a.ExportPlots1,a.TitleFontSize1,a.AxesLabelFontSize1,a.Frequency21,a.GridLine1,a.Title1,a.Length1,a.LineWidth1,a.Mode21,a.FileName1,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples21,a.Samples31,a.Gain1,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,a.Undistorted1,1,a.Halfspaces21,a.HalfspacesNumber21)
            elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                ModeShapeGrid_Isotropic_Rod_Pipe(a.Geometry1,a.Plane1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,a.F1,a.L1,a.T1,a.Directory1,a.ExportPlots1,a.TitleFontSize1,a.AxesLabelFontSize1,a.Frequency21,a.GridLine1,a.Title1,a.Length1,a.LineWidth1,a.Mode21,a.FileName1,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples21,a.Samples31,a.Gain1,a.Undistorted1,a.Halfspaces21,a.HalfspacesNumber21)
            elseif strcmp(a.Geometry1,'Circumferential')
                ModeShapeGrid_Isotropic_Circumferential(a.PNGresolution1,a.Material1,a.CLamb1,a.CShear1,a.Directory1,a.ExportPlots1,a.TitleFontSize1,a.AxesLabelFontSize1,a.Frequency21,a.GridLine1,a.Title1,a.Length1,a.LineWidth1,a.Mode21,a.FileName1,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples21,a.Samples31,a.Gain1,a.Undistorted1)
            end
        else
            if  a.ModeShapeSettingChanged1
                a.ModeShapeSettingChanged1 = 0;
                if  strcmp(a.Geometry1,'Plate')
                    [a.Time1,a.u1,a.x11,a.x3Total,a.p1] = ModeShapeGridAnimationComputer({a.Material1},1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,{conj(a.Material1.C)},0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.CycleDuration1,a.FrameRate1,a.Frequency21,a.Length1,a.Mode21,a.Thickness1/1e3,a.Samples21,a.Samples31,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces21,a.HalfspacesNumber21);
                elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                    [a.Time1,a.u1,a.ua1,a.ub1,a.x11,a.r1,a.Thetaa1,a.Thetab1,a.p1] = ModeShapeGridAnimationComputer_Isotropic_Rod_Pipe(a.Geometry1,a.Plane1,a.Material1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,a.F1,a.L1,a.T1,a.CycleDuration1,a.FrameRate1,a.Frequency21,a.Length1,a.Mode21,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples21,a.Samples31,a.Halfspaces21,a.HalfspacesNumber21);
                elseif strcmp(a.Geometry1,'Circumferential')
                    [a.Time1,a.u1,a.x11,a.r1,a.p1] = ModeShapeGridAnimationComputer_Isotropic_Circumferential(a.Material1,a.CLamb1,a.CShear1,a.CycleDuration1,a.FrameRate1,a.Frequency21,a.Length1,a.Mode21,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples21,a.Samples31);
                end
            end
            if  strcmp(a.Geometry1,'Plate')
                ModeShapeGridAnimation(0,1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Directory1,a.TitleFontSize1,a.AxesLabelFontSize1,a.FrameRate1,a.Frequency21,a.GridLine1,a.Title1,a.LineWidth1,a.Material1,a.Mode21,a.ExportPlots1,a.FileName1,a.MovieQuality1,a.Thickness1,0,a.Samples31,a.Gain1,1,0,a.Time1,a.u1,a.Undistorted1,a.x11,a.x3Total,a.p1,a.Halfspaces21,a.HalfspacesNumber21)
            elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                ModeShapeGridAnimation_Rod_Pipe(a.Geometry1,a.Plane1,0,1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Directory1,a.TitleFontSize1,a.AxesLabelFontSize1,a.FrameRate1,a.Frequency21,a.GridLine1,a.Title1,a.LineWidth1,a.Material1,a.Mode21,a.ExportPlots1,a.FileName1,a.MovieQuality1,a.Thickness1/2,a.ThicknessInner1/2,0,a.Samples31,a.Gain1,1,0,a.Time1,a.u1,a.ua1,a.ub1,a.Undistorted1,a.x11,a.r1,a.Thetaa1,a.Thetab1,a.p1,a.Halfspaces21,a.HalfspacesNumber21)
            elseif strcmp(a.Geometry1,'Circumferential')
                ModeShapeGridAnimation_Circumferential(0,1,'',a.Directory1,a.TitleFontSize1,a.AxesLabelFontSize1,a.FrameRate1,a.Frequency21,a.GridLine1,a.Title1,a.LineWidth1,a.Material1,a.Mode21,a.ExportPlots1,a.FileName1,a.MovieQuality1,a.Thickness1/2,a.ThicknessInner1/2,0,a.Gain1,1,0,a.Time1,a.u1,a.Undistorted1,a.x11,a.r1,a.p1)
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
elseif strcmp(source.Tag,'48') % Export plots
    a.ExportPlots1 = source.Value;
    if  source.Value
        a.Plot1UI1.String = 'Export';
        a.Plot2UI1.String = 'Export';
        a.Plot3UI1.String = 'Export';
    else
        a.Plot1UI1.String = 'Plot';
        a.Plot2UI1.String = 'Plot';
        a.Plot3UI1.String = 'Plot';
    end
elseif strcmp(source.Tag,'49') % PDF
    a.PDF1 = source.Value;
elseif strcmp(source.Tag,'50') % PNG
    a.PNG1 = source.Value;
elseif strcmp(source.Tag,'51') % PNG resolution
    a.PNGresolution1 = str2double(source.String);
elseif strcmp(source.Tag,'52') % X-Axis mode (Export settings)
    a.XAxisMode21 = source.Value;
elseif strcmp(source.Tag,'53') % Arrange
    a.Arrange1 = source.Value;
elseif strcmp(source.Tag,'54') % Dispersion curves
    a.DispersionCurves1 = source.Value;
    switch source.Value
    case 1
        a.XAxisMode2UI1.Enable = 'on';
        a.ArrangeUI1.Enable = 'on';
    case 0
        a.XAxisMode2UI1.Enable = 'off';
        a.ArrangeUI1.Enable = 'off';
    end
elseif strcmp(source.Tag,'55') % Through thickness
    a.ThroughThickness1 = source.Value;
elseif strcmp(source.Tag,'78') || strcmp(source.Tag,'56') || strcmp(source.Tag,'57') % Matlab Excel txt
    if  strcmp(source.Tag,'78') % mat
        DataFormat = 1;
    elseif strcmp(source.Tag,'56') % xlsx
        DataFormat = 2;
    elseif strcmp(source.Tag,'57') % txt
        DataFormat = 3;
    end
    try
        if  a.DispersionCurves1
            Export(a.Geometry1,DataFormat,a.SLamb1,a.ALamb1,a.BLamb1,a.SShear1,a.AShear1,{[]},a.F1,a.L1,a.T1,a.CLamb1,a.CShear1,a.Arrange1,a.XAxisMode21,a.Distance1,a.Couplant1,a.Thickness1,a.ThicknessInner1,a.Directory1,a.FileName1)
        end
        if  a.ThroughThickness1
            if  a.Quantity21 == 1
                if  strcmp(a.Geometry1,'Plate')
                    ModeShapeLines(1,2,DataFormat,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
                elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                    ModeShapeLines_Isotropic_Rod_Pipe(1,2,DataFormat,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,0,0,0,0,0,0,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
                elseif strcmp(a.Geometry1,'Circumferential')
                    ModeShapeLines_Isotropic_Circumferential(1,2,DataFormat,a.Plot1,a.PNGresolution1,a.Material1,0,0,0,0,0,0,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Phase1)
                end
            elseif a.Quantity21 == 2
                if  strcmp(a.Geometry1,'Plate')
                    ModeShapeLines(2,2,DataFormat,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
                elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                    ModeShapeLines_Isotropic_Rod_Pipe(2,2,DataFormat,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,0,0,0,0,0,0,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
                elseif strcmp(a.Geometry1,'Circumferential')
                    ModeShapeLines_Isotropic_Circumferential(2,2,DataFormat,a.Plot1,a.PNGresolution1,a.Material1,0,0,0,0,0,0,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Phase1)
                end
            elseif a.Quantity21 == 3
                if  strcmp(a.Geometry1,'Plate')
                    ModeShapeLines(3,2,DataFormat,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
                elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                    ModeShapeLines_Isotropic_Rod_Pipe(3,2,DataFormat,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,0,0,0,0,0,0,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
                elseif strcmp(a.Geometry1,'Circumferential')
                    ModeShapeLines_Isotropic_Circumferential(3,2,DataFormat,a.Plot1,a.PNGresolution1,a.Material1,0,0,0,0,0,0,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Phase1)
                end
            elseif a.Quantity21 == 4
                if  strcmp(a.Geometry1,'Plate')
                    ModeShapeLines(4,2,DataFormat,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,0)
                elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                    ModeShapeLines_Isotropic_Rod_Pipe(4,2,DataFormat,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,0,0,0,0,0,0,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
                elseif strcmp(a.Geometry1,'Circumferential')
                    ModeShapeLines_Isotropic_Circumferential(4,2,DataFormat,a.Plot1,a.PNGresolution1,a.Material1,0,0,0,0,0,0,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,0)
                end
            elseif a.Quantity21 == 5
                if  strcmp(a.Geometry1,'Plate')
                    ModeShapeLines(5,2,DataFormat,0,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,'',a.Plot1,a.PNGresolution1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.BLamb1,{0},a.SLamb1,a.SShear1,a.BoxLineWidth1,{conj(a.Material1.C)},0,{a.Material1},1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/1e3,0,a.Samples11,0,[1 -1;-1 1],1,[],1,a.Thickness1/1e3,0,1,a.Halfspaces11,a.HalfspacesNumber11,0)
                elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe')
                    ModeShapeLines_Isotropic_Rod_Pipe(5,2,DataFormat,a.Geometry1,a.Plot1,a.PNGresolution1,a.Material1,a.FluidLoading1,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,0,0,0,0,0,0,a.F1,a.L1,a.T1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
                elseif strcmp(a.Geometry1,'Circumferential')
                    ModeShapeLines_Isotropic_Circumferential(5,2,DataFormat,a.Plot1,a.PNGresolution1,a.Material1,0,0,0,0,0,0,a.CLamb1,a.CShear1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.Samples11,0)
                end
            end
        end
    catch ME
        if  strcmp(ME.identifier,'MATLAB:nonExistentField')
            errordlg('You must trace dispersion curves before you can export any data. Press ''Traces modes''.','Trace dispersion curves')
        end
    end
elseif strcmp(source.Tag,'58') % File name
    a.FileName1 = source.String;
elseif strcmp(source.Tag,'59') % Directory
    a.Directory1 = source.String;
elseif strcmp(source.Tag,'60') % Title
    a.Title1 = source.Value;
elseif strcmp(source.Tag,'62') % Legend location
    switch source.Value
    case 1
        a.LegendLocation1 = 'out';
    case 2
        a.LegendLocation1 = 'in';
    end
elseif strcmp(source.Tag,'63') % Box line width
    a.BoxLineWidth1 = str2double(source.String);
elseif strcmp(source.Tag,'64') % Curve line width
    a.LineWidth1 = str2double(source.String);
elseif strcmp(source.Tag,'66') % S
    a.SColor1 = eval(source.String);
elseif strcmp(source.Tag,'67') % A
    a.AColor1 = eval(source.String);
elseif strcmp(source.Tag,'69') % B
    a.BColor1 = eval(source.String);
elseif strcmp(source.Tag,'71') % Title (Font size)
    a.TitleFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'72') % Axes labels
    a.AxesLabelFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'73') % Axes ticks
    a.AxesTickFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'75') % Legend
    a.LegendFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'76') % Default
    a.Title1 = 1;
    a.TitleUI1.Value = a.Title1;
    a.LegendLocation1 = 'out';
    a.LegendLocationUI1.Value = 1;
    a.BoxLineWidth1 = .5;
    a.BoxLineWidthUI1.String = a.BoxLineWidth1;
    a.LineWidth1 = 1;
    a.LineWidthUI1.String = a.LineWidth1;
    a.SColor1 = [1 0 0];
    a.SColorUI1.String = '[1 0 0]';
    a.AColor1 = [0 0 1];
    a.AColorUI1.String = '[0 0 1]';
    a.BColor1 = [.5 0 1];
    a.BColorUI1.String = '[.5 0 1]';
    a.TitleFontSize1 = 30;
    a.TitleFontSizeUI1.String = a.TitleFontSize1;
    a.AxesLabelFontSize1 = 30;
    a.AxesLabelFontSizeUI1.String = a.AxesLabelFontSize1;
    a.AxesTickFontSize1 = 24;
    a.AxesTickFontSizeUI1.String = a.AxesTickFontSize1;
    a.LegendFontSize1 = 24;
    a.LegendFontSizeUI1.String = a.LegendFontSize1;
end
function Enable2DTracingSettings
    if  strcmp(a.Geometry1,'Circumferential')
        a.Force2DTracingUI1.Enable = 'off';
        a.PhaseVelocityLimit2UI1.Enable = 'off';
        a.AttenuationLimitUI1.Enable = 'off';
        a.SweepsUI1.Enable = 'off';
        a.SweepsSectionsUI1.Enable = 'off';
        a.HigherOrderModesUI1.Enable = 'on';
    else
        if  a.Viscoelastic1 || ((strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')) && a.FluidLoading1) || (strcmp(a.Geometry1,'Pipe') && (a.ToggleUpperFluid1 || (a.ToggleLowerFluid1 && a.Sink1)))
            a.Force2DTracingUI1.Enable = 'off';
            a.PhaseVelocityLimit2UI1.Enable = 'on';
            a.AttenuationLimitUI1.Enable = 'on';
            a.SweepsUI1.Enable = 'on';
            a.SweepsSectionsUI1.Enable = 'on';
            a.HigherOrderModesUI1.Enable = 'off';
        else
            a.Force2DTracingUI1.Enable = 'on';
            if  a.Force2DTracing1
                a.PhaseVelocityLimit2UI1.Enable = 'on';
                a.AttenuationLimitUI1.Enable = 'on';
                a.SweepsUI1.Enable = 'on';
                a.SweepsSectionsUI1.Enable = 'on';
                a.HigherOrderModesUI1.Enable = 'off';
            else
                a.PhaseVelocityLimit2UI1.Enable = 'off';
                a.AttenuationLimitUI1.Enable = 'off';
                a.SweepsUI1.Enable = 'off';
                a.SweepsSectionsUI1.Enable = 'off';
                a.HigherOrderModesUI1.Enable = 'on';
            end
        end
    end
end
function AdjustBasicParameters
    if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Circumferential')
        a.XRange1 = 2; % 10000/10 determines how broad the dispersion diagram is by default
        a.XSamples1 = 1000; % determines how many samples a dispersion curve consists of by default
        a.BelowCutoffWidth1 = 50;  % (kHz/mm); if we exceed the allowed scanning width (in kHz/mm) below the cut-off frequency of a damped mode without finding it, the damped search stops, and only the nondamped tracing continues to have that data for the next mode; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one
    elseif strcmp(a.Geometry1,'Pipe')
        if  ~a.FluidLoading1 && ~a.Viscoelastic1
            a.XRange1 = 1; % 5000/5
            a.XSamples1 = 1000;
            a.BelowCutoffWidth1 = 50;
        else
            if  a.ToggleLowerFluid1 && ~a.Sink1
                a.XRange1 = 1; % 500/1.25
                a.BelowCutoffWidth1 = 12.5;
            else
                a.XRange1 = .4; % 2000/5
                a.BelowCutoffWidth1 = 50;
            end
            a.XSamples1 = 400;
        end
    end
    AdjustFrequencyRange
end
function AdjustFrequencyRange
    if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod') || (strcmp(a.Geometry1,'Pipe') && a.ToggleLowerFluid1 && ~a.Sink1)
        a.FrequencyLimit1 = round(a.Material1.PlateVelocity*a.XRange1/a.Thickness1,-2);
    elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
        a.FrequencyLimit1 = round(a.Material1.PlateVelocity*a.XRange1/(a.Thickness1-a.ThicknessInner1)*2,-2);
    end
    if  a.FrequencyLimit1 == 0
        for i = -1:100 %#ok<*FXUP>
            a.FrequencyLimit1 = round(a.Material1.PlateVelocity*a.XRange1/a.Thickness1,i);
            if  a.FrequencyLimit1 > 0
                break
            end
        end
    elseif a.FrequencyLimit1 >= 1e4
        a.FrequencyLimit1 = round(a.FrequencyLimit1,-3);
    end
    a.FrequencyResolution1 = a.FrequencyLimit1/a.XSamples1;
    a.FrequencyLimitUI1.String = a.FrequencyLimit1;
    a.FrequencyResolutionUI1.String = a.FrequencyResolution1;
    a.Step1 = a.FrequencyLimit1/a.Steps1;
    a.StepUI1.String = a.Step1;
    a.Frequency11 = a.FrequencyLimit1;
    a.Frequency1UI1.String = a.FrequencyLimit1;
    a.Frequency21 = a.FrequencyLimit1;
    a.Frequency2UI1.String = a.FrequencyLimit1;
    AdjustAxes
end
function AdjustAxes
    if  a.Quantity11 == 3
        if  a.XAxisMode1 == 1
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness1),']'];
                a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness1;
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/2e3*(a.Thickness1-a.ThicknessInner1)),']'];
                a.YAxis1 = eval(a.YAxisUI1.String)*2e3/(a.Thickness1-a.ThicknessInner1);
            end
        end
    else
        if  a.XAxisMode1 == 1
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness1),']'];
                a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness1;
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/2e3*(a.Thickness1-a.ThicknessInner1)),']'];
                a.XAxis1 = eval(a.XAxisUI1.String)*2e3/(a.Thickness1-a.ThicknessInner1);
            end
        end
    end
    if  a.Quantity11 == 1
        a.YAxisUI1.String = ['[0 ',num2str(a.PhaseVelocityLimit1/1e3),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    elseif a.Quantity11 == 2
        a.YAxisUI1.String = ['[0 ',num2str(ceil(a.Material1.PlateVelocity/1e3)),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    elseif a.Quantity11 == 3
        x = a.Distance1/a.Material1.PlateVelocity*5e3;
        if  x >= 1e4
            x = round(x,-3);
        elseif x >= 1e3 && x < 1e4
            x = round(x,-2);
        elseif x >= 1e2 && x < 1e3
            x = round(x,-1);
        elseif x < 1e2
            x = round(x);
        end
        if  x == 0
            x = 1;
        end
        a.XAxisUI1.String = ['[0 ',num2str(x),']'];
        a.XAxis1 = eval(a.XAxisUI1.String);
    elseif a.Quantity11 == 4
        x = round(8e4/a.Material1.PlateVelocity*a.Couplant1.Velocity/343);
        if  x > 90
            x = 90;
        end
        a.YAxisUI1.String = ['[0 ',num2str(x),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    elseif a.Quantity11 == 5
        if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
            x = a.YRange1*1e-3*a.Material1.PlateVelocity*a.Thickness1;
        elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
            x = a.YRange1*1e-3*a.Material1.PlateVelocity*(a.Thickness1-a.ThicknessInner1)/2;
        end
        if  x >= 1e3 && x < 1e4
            x = round(x,-2);
        elseif x >= 1e2 && x < 1e3
            x = round(x,-1);
        elseif x < 1e2
            x = round(x);
        end
        if  a.XAxisMode1 == 3
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                a.YAxisUI1.String = ['[0 ',num2str(x/a.Thickness1),']'];
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                a.YAxisUI1.String = ['[0 ',num2str(x/(a.Thickness1-a.ThicknessInner1)*2),']'];
            end
            a.YAxis1 = eval(a.YAxisUI1.String);
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        end
    elseif a.Quantity11 == 6 || a.Quantity11 == 7
        if  a.Quantity11 == 6
            x = 2*pi*a.FrequencyLimit1/a.Material1.RayleighVelocity;
        elseif a.Quantity11 == 7
            if  a.Symmetric1
                FluidVelocity = a.UpperFluid1.Velocity;
                FluidDensity = a.UpperFluid1.Density;
            else
                if  a.ToggleUpperFluid1 && a.ToggleLowerFluid1
                    FluidVelocity = .5*(a.UpperFluid1.Velocity+a.LowerFluid1.Velocity);
                    FluidDensity = .5*(a.UpperFluid1.Density+a.LowerFluid1.Density);
                elseif a.ToggleUpperFluid1 && ~a.ToggleLowerFluid1
                    FluidVelocity = .5*a.UpperFluid1.Velocity;
                    FluidDensity = .5*a.UpperFluid1.Density;
                elseif ~a.ToggleUpperFluid1 && a.ToggleLowerFluid1
                    FluidVelocity = .5*a.LowerFluid1.Velocity;
                    FluidDensity = .5*a.LowerFluid1.Density;
                end
            end
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                if  a.FluidLoading1 && a.Viscoelastic1
                    x = 1e4*(FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+a.Material1.Density*a.Material1.PlateVelocity)+a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness1;
                elseif a.FluidLoading1 && ~a.Viscoelastic1
                    x = 1e4*FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+a.Material1.Density*a.Material1.PlateVelocity)/a.Thickness1;
                elseif ~a.FluidLoading1 && a.Viscoelastic1
                    x = 1e4*(a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness1;
                else
                    x = 1;
                end
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                if  a.FluidLoading1 && a.Viscoelastic1
                    x = 1e4*(FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+a.Material1.Density*a.Material1.PlateVelocity)+a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/(a.Thickness1-a.ThicknessInner1)*2;
                elseif a.FluidLoading1 && ~a.Viscoelastic1
                    x = 1e4*FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+a.Material1.Density*a.Material1.PlateVelocity)/(a.Thickness1-a.ThicknessInner1)*2;
                elseif ~a.FluidLoading1 && a.Viscoelastic1
                    x = 1e4*(a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/(a.Thickness1-a.ThicknessInner1)*2;
                else
                    x = 1;
                end
            end
        end
        if  x >= 1e3 && x < 1e4
            x = 1e2*ceil(x/1e2);
        elseif x >= 1e2 && x < 1e3
            x = 1e1*ceil(x/1e1);
        elseif x >= 1e1 && x < 1e2
            x = ceil(x);
        elseif x < 1e1
            x = 1e-1*ceil(x/1e-1);
        end
        if  a.XAxisMode1 == 3
            if  strcmp(a.Geometry1,'Plate') || strcmp(a.Geometry1,'Rod')
                a.YAxisUI1.String = ['[0 ',num2str(x*a.Thickness1),']'];
            elseif strcmp(a.Geometry1,'Pipe') || strcmp(a.Geometry1,'Circumferential')
                a.YAxisUI1.String = ['[0 ',num2str(x*(a.Thickness1-a.ThicknessInner1)/2),']'];
            end
            a.YAxis1 = eval(a.YAxisUI1.String);
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        end
    end
end
function CheckSymmetry
    if  (~a.ToggleUpperFluid1 && ~a.ToggleLowerFluid1) || (a.ToggleUpperFluid1 && a.ToggleLowerFluid1 && strcmp(a.UpperFluid1.Name,a.LowerFluid1.Name))
        a.Symmetric1 = 1;
        if  strcmp(a.Geometry1,'Plate')
            a.Symmetric_Torsional_ModesUI1.Enable = 'on';
            a.Antisymmetric_Longitudinal_ModesUI1.Enable = 'on';
        end
    else
        a.Symmetric1 = 0;
        if  strcmp(a.Geometry1,'Plate')
            a.Symmetric_Torsional_ModesUI1.Enable = 'off';
            a.Antisymmetric_Longitudinal_ModesUI1.Enable = 'off';
        end
    end
end
function AdjustSamples31
    if  strcmp(a.Geometry1,'Plate')
        if  a.FluidLoading1 && a.Halfspaces21
            if  a.ToggleUpperFluid1 && a.ToggleLowerFluid1
                a.Samples31 = round(.5*a.Samples21/(2*a.HalfspacesNumber21+1));
            else
                a.Samples31 = round(.5*a.Samples21/(a.HalfspacesNumber21+1));
            end
        else
            a.Samples31 = round(.5*a.Samples21);
        end
    elseif strcmp(a.Geometry1,'Rod')
        if  a.FluidLoading1 && a.Halfspaces21
            a.Samples31 = round(.5*a.Samples21/(a.HalfspacesNumber21+1));
        else
            a.Samples31 = round(.5*a.Samples21);
        end
    elseif strcmp(a.Geometry1,'Pipe')
        if  a.ToggleUpperFluid1 && a.Halfspaces21 && a.ToggleLowerFluid1
            a.Samples31 = round(.5*a.Samples21/(a.HalfspacesNumber21+2));
        elseif a.ToggleUpperFluid1 && a.Halfspaces21 && ~a.ToggleLowerFluid1
            a.Samples31 = round(.5*a.Samples21/(a.HalfspacesNumber21+1));
        elseif ((a.ToggleUpperFluid1 && ~a.Halfspaces21) || ~a.ToggleUpperFluid1) && a.ToggleLowerFluid1
            a.Samples31 = round(.25*a.Samples21);
        else
            a.Samples31 = round(.5*a.Samples21);
        end
    elseif strcmp(a.Geometry1,'Circumferential')
        a.Samples31 = round(.5*a.Samples21);
    end
    if  a.Samples31 == 0
        a.Samples31 = 1;
    end
    a.Samples3UI1.String = a.Samples31;
end
function FrequencySweeps
    FrequencySweepRange = a.Step1:a.Step1:a.FrequencyLimit1;
    if  strcmp(a.Geometry1,'Plate')
        [a.HSLamb1,a.HALamb1,a.HBLamb1,a.HSShear1,a.HAShear1] = FrequencySweeper_Isotropic(FrequencySweepRange,a.Material1,a.PhaseVelocityLimit1,a.Thickness1/2e3,a.Symmetric1,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
    elseif strcmp(a.Geometry1,'Rod')
        [a.HF1,a.HL1,a.HT1] = FrequencySweeper_Isotropic_Rod_Pipe(a.Geometry1,a.Material1,a.Thickness1/2e3,0,a.UpperFluid1,a.LowerFluid1,0,0,0,FrequencySweepRange,a.FlexuralModeOrders1,a.PhaseVelocityLimit1,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
    elseif strcmp(a.Geometry1,'Pipe')
        if  a.Force2DTracing1 || a.Viscoelastic1 || a.ToggleUpperFluid1 || (a.ToggleLowerFluid1 && a.Sink1)
            [a.HF1,a.HL1,a.HT1] = FrequencySweeper_Isotropic_Rod_Pipe(a.Geometry1,a.Material1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.UpperFluid1,a.LowerFluid1,0,0,0,FrequencySweepRange,a.FlexuralModeOrders1,a.PhaseVelocityLimit1,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
            if  a.ToggleLowerFluid1 && ~a.Sink1
                [HF2,HL2,~] = FrequencySweeper_Isotropic_Rod_Pipe(a.Geometry1,a.Material1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,FrequencySweepRange,a.FlexuralModeOrders1,a.PhaseVelocityLimit1,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
                for n = 1:length(a.HF1)
                    for i = 1:length(a.HF1{n})
                        if  ~any(isapprox(HF2{n},a.HF1{n}(i),RelativeTolerance=.05))
                            HF2{n} = sort([HF2{n} a.HF1{n}(i)]);
                        end
                    end
                    a.HF1{n} = sort(HF2{n});
                end
                for i = 1:length(a.HL1)
                    if  ~any(isapprox(HL2,a.HL1(i),RelativeTolerance=.05))
                        HL2 = [HL2 a.HL1(i)];
                    end
                end
                a.HL1 = sort(HL2);
            end
        else
            if  ~a.ToggleUpperFluid1 && a.ToggleLowerFluid1 && ~a.Sink1 && ~a.Viscoelastic1
                [a.HF1,a.HL1,a.HT1] = FrequencySweeper_Isotropic_Rod_Pipe(a.Geometry1,a.Material1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.UpperFluid1,a.LowerFluid1,a.ToggleUpperFluid1,a.ToggleLowerFluid1,a.Sink1,FrequencySweepRange,a.FlexuralModeOrders1,a.PhaseVelocityLimit1,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
            else
                [a.HF1,a.HL1,a.HT1] = FrequencySweeper_Isotropic_Rod_Pipe(a.Geometry1,a.Material1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.UpperFluid1,a.LowerFluid1,0,0,0,FrequencySweepRange,a.FlexuralModeOrders1,a.PhaseVelocityLimit1,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
            end
        end
    elseif strcmp(a.Geometry1,'Circumferential')
        [a.HCLamb1,a.HCShear1] = FrequencySweeper_Isotropic_Circumferential(FrequencySweepRange,a.Material1,a.PhaseVelocityLimit1,a.Thickness1/2e3,a.ThicknessInner1/2e3,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
    end
    a.Search1 = 1;
end
end