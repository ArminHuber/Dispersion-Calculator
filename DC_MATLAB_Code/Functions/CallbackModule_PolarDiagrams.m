% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2024 DLR
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
function a = CallbackModule_PolarDiagrams(source,~,a)
%#ok<*AGROW>
%#ok<*GFLD>
%#ok<*GVMIS>
global Stop
Stop = 0;
source.String = replace(source.String,',','.');
if  strcmp(source.Tag,'1') % Frequency limit
    a.FrequencyLimit_Polar = str2double(source.String);
elseif strcmp(source.Tag,'2') % Frequency step
    a.FrequencyResolution_Polar = str2double(source.String);
elseif strcmp(source.Tag,'34') % Phase velocity accuracy
    a.PhaseVelocityResolution_Polar = str2double(source.String);
elseif strcmp(source.Tag,'33') % Propagation angle limit
    a.PropagationAngleMode_Polar = source.Value;
elseif strcmp(source.Tag,'3') % Propagation angle step
    a.PropagationAngleStep_Polar = str2double(source.String);
elseif strcmp(source.Tag,'4') % Phase velocity sections
    a.PhaseVelocitySections_Polar = str2double(source.String);
elseif strcmp(source.Tag,'5') % S1/B2
    a.S0_Polar = source.Value;
elseif strcmp(source.Tag,'6') % S0/B1
    a.SH0_Polar = source.Value;
elseif strcmp(source.Tag,'7') % A0/B0
    a.A0_Polar = source.Value;
elseif strcmp(source.Tag,'8') % Trace modes/Stop trace
    if  source.Value
        FrequencyRange = 0:a.FrequencyResolution_Polar:a.FrequencyLimit_Polar;
        FrequencyRange(1) = a.FrequencyRangeStart2;
        if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution_Polar
            FrequencyRange(1) = .1*a.FrequencyResolution_Polar;
        elseif a.FrequencyRangeStart2 < 1e-4*a.FrequencyResolution_Polar
            FrequencyRange(1) = 1e-4*a.FrequencyResolution_Polar;
        end
        if  a.PropagationAngleMode_Polar == 1
            PropagationAngle = 0:a.PropagationAngleStep_Polar:180;
        elseif a.PropagationAngleMode_Polar == 2
            PropagationAngle = 0:a.PropagationAngleStep_Polar:90;
        end
        for i = 1:length(PropagationAngle)
            Phi_Polar(i,:) = PropagationAngle(i)-a.LayerOrientations_Polar;
        end
        for i = 1:length(a.Material_Polar)
            if  strcmp(a.MaterialClasses_Polar(i),'Cubic')
                Phi_Polar(mod(Phi_Polar(:,i),45) == 0,i) = Phi_Polar(mod(Phi_Polar(:,i),45) == 0,i)+1e-1;
            elseif strcmp(a.MaterialClasses_Polar(i),'Transversely isotropic') || strcmp(a.MaterialClasses_Polar(i),'Orthotropic')
                Phi_Polar(mod(Phi_Polar(:,i),90) == 0,i) = Phi_Polar(mod(Phi_Polar(:,i),90) == 0,i)+1e-1;
            end
            if  isreal(a.Material_Polar{i}.C)
                Viscoelastic(i) = 0;
            else
                Viscoelastic(i) = 1;
                a.Material_Polar{i}.C = real(a.Material_Polar{i}.C);
            end
        end
        if  any(Viscoelastic)
            warndlg('A viscoelastic material is contained in the layup. The damping is ignored, i.e., only the real parts of the stiffness components are taken.','Warning');
        end
        if  a.Multithreading
            imshow('Multithreading22.png','Parent',a.a1UI4)
        end
        source.String = 'Stop trace';
        source.TooltipString = 'Stop tracing dispersion curves.';
        [a.A_Polar,a.c,a.a11,a.a12,a.a21,a.a22,a.a23,a.a31,a.a32,a.a33,a.a34] = Computer_Polar(a.Multithreading,a.Hybrid_Polar,a.MaterialClasses_Polar,a.A0_Polar,a.FrequencyLimit_Polar,FrequencyRange,a.LambPhaseVelocitySweepRange12,a.LambPhaseVelocitySweepRange22,a.LayerThicknesses_Polar/1e3,a.Material_Polar,a.PhaseVelocityResolution_Polar,a.PhaseVelocitySections_Polar,Phi_Polar,a.PlateThickness_Polar/1e3,PropagationAngle,a.SuperLayers_Polar,a.Pattern_Polar,a.SH0_Polar,a.S0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar,a.MissingSamples_Polar);
        source.Value = 0;
        source.String = 'Trace modes';
        source.TooltipString = 'Start tracing the dispersion curves.';
        if  Stop
            return
        end
        a.FrequencyUI4.Value = 1;
        a.FrequencyUI4.String = {fliplr(FrequencyRange)};
        a.Frequency_Polar = FrequencyRange(end);
        a.CalculateCeUI4.Value = 1;
        a.CalculateCeUI4.String = 'Stop calculate';
        a.CalculateCeUI4.TooltipString = 'Stop calculating the energy velocity.';
        a.A_Polar = Computer_Polar_EnergyVelocity(a.Multithreading,a.A_Polar,a.c,a.a11,a.a12,a.a21,a.a22,a.a23,a.a31,a.a32,a.a33,a.a34,a.Material_Polar,a.PlateThickness_Polar/1e3,a.SuperLayers_Polar,a.Pattern_Polar,a.SuperLayerSize_Polar,a.LayerThicknesses_Polar/1e3,a.SymmetricSystem_Polar,FrequencyRange,[1 1 -1;-1 -1 1;-1 -1 1],a.SamplesX3ce_Polar);
        a.CalculateCeUI4.Value = 0;
        a.CalculateCeUI4.String = 'Calculate ce';
        a.CalculateCeUI4.TooltipString = 'Start calculating the energy velocity.';
        if  a.Multithreading
            imshow('Multithreading21.png','Parent',a.a1UI4)
        end
    else
        Stop = 1;
        if  a.Multithreading
            source.Value = 0;
            source.String = 'Trace modes';
            source.TooltipString = 'Start tracing the dispersion curves.';
            imshow('Multithreading21.png','Parent',a.a1UI4)
            cancelAll(a.Pool.FevalQueue)
        end
    end
elseif strcmp(source.Tag,'35') % Samples x3 (ce)
    a.SamplesX3ce_Polar = str2double(source.String);
elseif strcmp(source.Tag,'9') % Calculate ce/Stop calculate
    if  source.Value && strcmp(source.String,'Calculate ce')
        if  a.Multithreading
            imshow('Multithreading22.png','Parent',a.a1UI4)
        end
        FrequencyRange = 0:a.FrequencyResolution_Polar:a.FrequencyLimit_Polar;
        FrequencyRange(1) = a.FrequencyRangeStart2;
        if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution_Polar
            FrequencyRange(1) = .1*a.FrequencyResolution_Polar;
        elseif a.FrequencyRangeStart2 < 1e-4*a.FrequencyResolution_Polar
            FrequencyRange(1) = 1e-4*a.FrequencyResolution_Polar;
        end
        source.String = 'Stop calculate';
        source.TooltipString = 'Stop calculating the energy velocity.';
        a.A_Polar = Computer_Polar_EnergyVelocity(a.Multithreading,a.A_Polar,a.c,a.a11,a.a12,a.a21,a.a22,a.a23,a.a31,a.a32,a.a33,a.a34,a.Material_Polar,a.PlateThickness_Polar/1e3,a.SuperLayers_Polar,a.Pattern_Polar,a.SuperLayerSize_Polar,a.LayerThicknesses_Polar/1e3,a.SymmetricSystem_Polar,FrequencyRange,[1 1 -1;-1 -1 1;-1 -1 1],a.SamplesX3ce_Polar);
        source.Value = 0;
        source.String = 'Calculate ce';
        source.TooltipString = 'Start calculating the energy velocity.';
        if  a.Multithreading
            imshow('Multithreading21.png','Parent',a.a1UI4)
        end
    else
        Stop = 1;
    end
elseif strcmp(source.Tag,'10') % Quantity
    switch source.Value
    case 1
        a.Quantity_Polar = 1;
        if  strcmp(a.LayerCountUI4.String,'1')
            a.Option1UI4.Enable = 'on';
            a.Option1UI4.Style = 'checkbox';
            a.Option1UI4.String = ' ';
            a.Option1UI4.Value = a.BulkVelocities_Polar;
            a.Option1UI4.TooltipString = 'Check this to show the bulk wave velocities.';
            a.Option1UI4.Position(3) = 50;
            a.Option1TextUI4.String = 'Bulk velocities';
            a.Option1TextUI4.Position(3) = 70;
        else
            a.Option1UI4.Enable = 'off';
        end
    case 2
        a.Quantity_Polar = 2;
        if  strcmp(a.LayerCountUI4.String,'1')
            a.Option1UI4.Enable = 'on';
            a.Option1UI4.Style = 'checkbox';
            a.Option1UI4.String = ' ';
            a.Option1UI4.Value = a.BulkVelocities_Polar;
            a.Option1UI4.TooltipString = 'Check this to show the bulk wave velocities.';
            a.Option1UI4.Position(3) = 50;
            a.Option1TextUI4.String = 'Bulk velocities';
            a.Option1TextUI4.Position(3) = 70;
        else
            a.Option1UI4.Enable = 'off';
        end
    case 3
        a.Quantity_Polar = 3;
        a.Option1UI4.Enable = 'on';
        a.Option1UI4.Style = 'edit';
        a.Option1UI4.String = a.Distance_Polar;
        a.Option1UI4.TooltipString = ['Enter the propagation distance from the source to the',newline,'sensor for which the propagation time shall be calculated.'];
        a.Option1UI4.Position(3) = 50;
        a.Option1TextUI4.String = 'Distance (mm)';
        a.Option1TextUI4.Position(3) = 71;
    case 4
        a.Quantity_Polar = 4;
        a.Option1UI4.Enable = 'on';
        a.Option1UI4.Style = 'popupmenu';
        a.Option1UI4.String = fieldnames(a.Materials.Fluid);
        a.Option1UI4.Value = find(strcmp(a.Couplant_Polar.Name,a.Option1UI4.String));
        a.Option1UI4.TooltipString = 'Select for which surrounding medium the coincidence angle shall be displayed.';
        a.Option1UI4.Position(3) = 130;
        a.Option1TextUI4.String = 'Couplant';
        a.Option1TextUI4.Position(3) = 44;
    case 5
        a.Quantity_Polar = 5;
        a.Option1UI4.Enable = 'off';
    case 6
        a.Quantity_Polar = 6;
        a.Option1UI4.Enable = 'off';
    end
elseif strcmp(source.Tag,'11') % option
    if  a.Quantity_Polar < 3
        a.BulkVelocities_Polar = source.Value;    
    elseif a.Quantity_Polar == 3
        a.Distance_Polar = str2double(source.String);
    elseif a.Quantity_Polar == 4
        E = fieldnames(a.Materials.Fluid);
        a.Couplant_Polar = getfield(a.Materials.Fluid,E{source.Value}); 
    end
elseif strcmp(source.Tag,'12') % Frequency
    FrequencyRange = a.FrequencyLimit_Polar:-a.FrequencyResolution_Polar:0;
    FrequencyRange(end) = a.FrequencyRangeStart2;
    if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution_Polar
        FrequencyRange(end) = .1*a.FrequencyResolution_Polar;
    elseif a.FrequencyRangeStart2 < 1e-4*a.FrequencyResolution_Polar
        FrequencyRange(end) = 1e-4*a.FrequencyResolution_Polar;
    end
    a.Frequency_Polar = FrequencyRange(source.Value);
elseif strcmp(source.Tag,'13') % Plot
    FrequencyRange = 0:a.FrequencyResolution_Polar:a.FrequencyLimit_Polar;
    FrequencyRange(1) = a.FrequencyRangeStart2;
    if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution_Polar
        FrequencyRange(1) = .1*a.FrequencyResolution_Polar;
    elseif a.FrequencyRangeStart2 < 1e-4*a.FrequencyResolution_Polar
        FrequencyRange(1) = 1e-4*a.FrequencyResolution_Polar;
    end
    if  a.PropagationAngleMode_Polar == 1
        PropagationAngle = 0:a.PropagationAngleStep_Polar:180;
    elseif a.PropagationAngleMode_Polar == 2
        PropagationAngle = 0:a.PropagationAngleStep_Polar:90;
    end
    if  a.Quantity_Polar == 1
        PhaseVelocity_Polar(a.Hybrid_Polar,a.LayupString_Polar,a.PropagationAngleStep_Polar,a.BulkVelocities_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
    elseif a.Quantity_Polar == 2
        EnergyVelocity_Polar(3,a.Hybrid_Polar,a.LayupString_Polar,a.PropagationAngleStep_Polar,a.BulkVelocities_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
        EnergyVelocity_Polar(2,a.Hybrid_Polar,a.LayupString_Polar,a.PropagationAngleStep_Polar,a.BulkVelocities_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)           
        EnergyVelocity_Polar(1,a.Hybrid_Polar,a.LayupString_Polar,a.PropagationAngleStep_Polar,a.BulkVelocities_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
    elseif a.Quantity_Polar == 3
        PropagationTime_Polar(2,a.Hybrid_Polar,a.LayupString_Polar,a.Distance_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
        PropagationTime_Polar(1,a.Hybrid_Polar,a.LayupString_Polar,a.Distance_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
    elseif a.Quantity_Polar == 4
        CoincidenceAngle_Polar(a.Hybrid_Polar,a.LayupString_Polar,a.Couplant_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
    elseif a.Quantity_Polar == 5
        Wavelength_Polar(a.Hybrid_Polar,a.LayupString_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
    elseif a.Quantity_Polar == 6
        Wavenumber_Polar(a.Hybrid_Polar,a.LayupString_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.SuperLayers_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
    end
elseif strcmp(source.Tag,'14') % Export
    a.ExportPlots_Polar = source.Value;
    if  source.Value
        a.PlotUI4.String = 'Export';
    else
        a.PlotUI4.String = 'Plot';
    end
elseif strcmp(source.Tag,'15') % PDF
    a.PDF_Polar = source.Value;
elseif strcmp(source.Tag,'16') % PNG
    a.PNG_Polar = source.Value;
elseif strcmp(source.Tag,'17') % PNG resolution
    a.PNGresolution_Polar = str2double(source.String);
elseif strcmp(source.Tag,'18') % Matlab
    Export_Polar(a,0,0,1)
elseif strcmp(source.Tag,'19') % Excel
    Export_Polar(a,1,0,0)
elseif strcmp(source.Tag,'20') % TXT
    Export_Polar(a,0,1,0)
elseif strcmp(source.Tag,'21') % File name
    a.FileName_Polar = source.String;
elseif strcmp(source.Tag,'22') % Directory
    a.Directory4 = source.String;
elseif strcmp(source.Tag,'23') % Title
    switch source.Value
    case 1
       a.Title_Polar = 2;
    case 2
       a.Title_Polar = 1;
    case 3
       a.Title_Polar = 0;
    end
elseif strcmp(source.Tag,'24') % Curve line width
    a.LineWidth_Polar = str2double(source.String);
elseif strcmp(source.Tag,'28') % Title (Font size)
    a.TitleFontSize_Polar = str2double(source.String);
elseif strcmp(source.Tag,'29') % Axes ticks
    a.AxesTickFontSize_Polar = str2double(source.String);
elseif strcmp(source.Tag,'30') % Mode labels
    a.ModeLabelFontSize_Polar = str2double(source.String);
elseif strcmp(source.Tag,'31') % Default
    a.Title_Polar = 2;
    a.TitleUI4.Value = 1;
    a.LineWidth_Polar = 1;
    a.LineWidthUI4.String = a.LineWidth_Polar;
    a.TitleFontSize_Polar = 30;
    a.TitleFontSizeUI4.String = a.TitleFontSize_Polar;
    a.AxesTickFontSize_Polar = 24;
    a.AxesTickFontSizeUI4.String = a.AxesTickFontSize_Polar;
    a.ModeLabelFontSize_Polar = 24;
    a.ModeLabelFontSizeUI4.String = a.ModeLabelFontSize_Polar;
end