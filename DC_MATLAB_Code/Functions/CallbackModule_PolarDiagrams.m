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
function a = CallbackModule_PolarDiagrams(source,a)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
source.String = replace(source.String,',','.');
if  strcmp(source.Tag,'1') % Frequency limit
    a.FrequencyLimit_Polar = str2double(source.String);
elseif strcmp(source.Tag,'2') % Frequency step
    a.FrequencyResolution_Polar = str2double(source.String);
elseif strcmp(source.Tag,'34') % Phase velocity accuracy
    a.Accuracy_Polar = str2double(source.String);
elseif strcmp(source.Tag,'33') % Propagation angle limit
    a.PropagationAngleMode_Polar = source.Value;
elseif strcmp(source.Tag,'3') % Propagation angle step
    a.PropagationAngleStep_Polar = str2double(source.String);
elseif strcmp(source.Tag,'4') % Phase velocity sections
    a.PhaseVelocitySections_Polar = str2double(source.String);
elseif strcmp(source.Tag,'84') % Matrix methods
    a.MatrixMethods_Polar = source.Value; 
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
        end
        if  a.PropagationAngleMode_Polar == 1
            PropagationAngle = 0:a.PropagationAngleStep_Polar:180;
        elseif a.PropagationAngleMode_Polar == 2
            PropagationAngle = 0:a.PropagationAngleStep_Polar:90;
        end
        for i = 1:length(PropagationAngle)
            Phi_Polar(i,:) = PropagationAngle(i)-a.LayerOrientations_Polar;
        end
        a.DC_Polar = zeros(size(Phi_Polar,1),length(a.Material_Polar));
        for i = 1:length(a.Material_Polar)
            for j = 1:size(Phi_Polar,1)
                if  strcmp(a.MaterialClasses_Polar(i),'Isotropic')
                    a.DC_Polar(j,i) = 1;
                elseif strcmp(a.MaterialClasses_Polar(i),'Cubic')
                    if  mod(Phi_Polar(j,i),45) % check whether the laminate consists only of 0deg, 45deg, and 90deg layers and the wave propagation is along 0deg, 45deg, or 90deg; then shear horizontal waves decouple from Lamb waves and they can be distinguished
                        a.DC_Polar(j,i) = 0;
                    else
                        a.DC_Polar(j,i) = 1;
                    end
                elseif strcmp(a.MaterialClasses_Polar(i),'Transversely isotropic') || strcmp(a.MaterialClasses_Polar(i),'Orthotropic')
                    if  mod(Phi_Polar(j,i),90) % check whether the laminate consists only of 0deg and 90deg layers
                        a.DC_Polar(j,i) = 0;
                    else
                        a.DC_Polar(j,i) = 1;
                    end
                end
            end
        end
        if  a.Multithreading
            imshow('Multithreading22.png','Parent',a.a1UI4)
        end
        source.String = 'Stop trace';
        source.Tooltip = 'Stop tracing dispersion curves.';
        [a.A_Polar,a.c,a.A1,a.a11,a.a12,a.a21,a.a22,a.a23,a.a31,a.a32,a.a33,a.a34] = Computer_Polar(a.Multithreading,a.A0_Polar,a.SH0_Polar,a.S0_Polar,a.Hybrid_Polar,a.FrequencyLimit_Polar,FrequencyRange,a.LambPhaseVelocitySweepRange12,a.LambPhaseVelocitySweepRange22,a.LayerThicknesses_Polar/1e3,a.Material_Polar,a.Accuracy_Polar,a.PhaseVelocitySections_Polar,Phi_Polar,a.PlateThickness_Polar/1e3,PropagationAngle,a.Pattern_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar,a.DC_Polar,a.MatrixMethods_Polar,a.MissingSamples_Polar);
        source.Value = 0;
        source.String = 'Trace modes';
        source.Tooltip = 'Start tracing the dispersion curves.';
        if  a.Multithreading
            imshow('Multithreading21.png','Parent',a.a1UI4)
        end
        if  Stop
            return
        end
        a.FrequencyUI4.Value = 1;
        a.FrequencyUI4.String = {fliplr(FrequencyRange)};
        a.Frequency_Polar = FrequencyRange(end);
        a.CalculateCeUI4.Value = 1;
        a.CalculateCeUI4.String = 'Stop calculate';
        a.CalculateCeUI4.Tooltip = 'Stop calculating the energy velocity.';
        a.A_Polar = Computer_Polar_EnergyVelocity(a.A_Polar,a.c,a.A1,a.a11,a.a12,a.a21,a.a22,a.a23,a.a31,a.a32,a.a33,a.a34,a.Material_Polar,a.Layers_Polar,a.SuperLayers_Polar,a.Pattern_Polar,a.SuperLayerSize_Polar,a.LayerThicknesses_Polar/1e3,a.SymmetricSystem_Polar,a.DC_Polar,FrequencyRange,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SamplesX3ce_Polar);
        a.CalculateCeUI4.Value = 0;
        a.CalculateCeUI4.String = 'Calculate ce';
        a.CalculateCeUI4.Tooltip = 'Start calculating the energy velocity.';
    else
        Stop = 1;
        if  a.Multithreading
            source.Value = 0;
            source.String = 'Trace modes';
            source.Tooltip = 'Start tracing the dispersion curves.';
            imshow('Multithreading21.png','Parent',a.a1UI4)
            cancelAll(a.Pool.FevalQueue)
        end
    end
elseif strcmp(source.Tag,'35') % Samples x3 (ce)
    a.SamplesX3ce_Polar = str2double(source.String);
elseif strcmp(source.Tag,'9') % Calculate ce/Stop calculate
    if  source.Value && strcmp(source.String,'Calculate ce')
        FrequencyRange = 0:a.FrequencyResolution_Polar:a.FrequencyLimit_Polar;
        FrequencyRange(1) = a.FrequencyRangeStart2;
        if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution_Polar
            FrequencyRange(1) = .1*a.FrequencyResolution_Polar;
        end
        source.String = 'Stop calculate';
        source.Tooltip = 'Stop calculating the energy velocity.';
        a.A_Polar = Computer_Polar_EnergyVelocity(a.A_Polar,a.c,a.A1,a.a11,a.a12,a.a21,a.a22,a.a23,a.a31,a.a32,a.a33,a.a34,a.Material_Polar,a.Layers_Polar,a.SuperLayers_Polar,a.Pattern_Polar,a.SuperLayerSize_Polar,a.LayerThicknesses_Polar/1e3,a.SymmetricSystem_Polar,a.DC_Polar,FrequencyRange,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SamplesX3ce_Polar);
        source.Value = 0;
        source.String = 'Calculate ce';
        source.Tooltip = 'Start calculating the energy velocity.';
    else
        Stop = 1;
    end
elseif strcmp(source.Tag,'10') % Quantity
    a.Quantity_Polar = source.Value;
    if  a.Quantity_Polar < 3
        if  a.SuperLayerSize_Polar == 1
            a.Option1UI4.Enable = 'on';
        else
            a.Option1UI4.Enable = 'off';
        end
        a.Option1UI4.Style = 'checkbox';
        a.Option1UI4.String = ' ';
        a.Option1UI4.Value = a.BulkVelocities_Polar;
        a.Option1UI4.Tooltip = 'Check this to show the bulk wave velocities.';
        a.Option1TextUI4.String = 'Bulk velocities';
    elseif a.Quantity_Polar == 3
        a.Option1UI4.Enable = 'on';
        a.Option1UI4.Style = 'edit';
        a.Option1UI4.String = a.Distance_Polar;
        a.Option1UI4.Tooltip = 'Enter the propagation distance from the source to the sensor for which the propagation time shall be calculated.';
        a.Option1TextUI4.String = 'Distance (mm)';
    elseif a.Quantity_Polar == 4
        a.Option1UI4.Enable = 'on';
        a.Option1UI4.Style = 'popupmenu';
        a.Option1UI4.String = a.FieldNames_Fluid;
        a.Option1UI4.Value = find(strcmp(a.Couplant_Polar.Name,a.Option1UI4.String));
        a.Option1UI4.Tooltip = 'Select for which surrounding medium the coincidence angle shall be displayed.';
        a.Option1TextUI4.String = 'Couplant';
    elseif a.Quantity_Polar > 4
        a.Option1UI4.Enable = 'off';
    end
elseif strcmp(source.Tag,'11') % option
    if  a.Quantity_Polar < 3
        a.BulkVelocities_Polar = source.Value;
    elseif a.Quantity_Polar == 3
        a.Distance_Polar = str2double(source.String);
    elseif a.Quantity_Polar == 4
        a.Couplant_Polar = getfield(a.Materials.Fluid,a.FieldNames_Fluid{source.Value});
    end
elseif strcmp(source.Tag,'12') % Frequency
    FrequencyRange = a.FrequencyLimit_Polar:-a.FrequencyResolution_Polar:0;
    FrequencyRange(end) = a.FrequencyRangeStart2;
    if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution_Polar
        FrequencyRange(end) = .1*a.FrequencyResolution_Polar;
    end
    a.Frequency_Polar = FrequencyRange(source.Value);
elseif strcmp(source.Tag,'13') % Plot
    try
        FrequencyRange = 0:a.FrequencyResolution_Polar:a.FrequencyLimit_Polar;
        FrequencyRange(1) = a.FrequencyRangeStart2;
        if  a.FrequencyRangeStart2 > .1*a.FrequencyResolution_Polar
            FrequencyRange(1) = .1*a.FrequencyResolution_Polar;
        end
        if  a.PropagationAngleMode_Polar == 1
            PropagationAngle = 0:a.PropagationAngleStep_Polar:180;
        elseif a.PropagationAngleMode_Polar == 2
            PropagationAngle = 0:a.PropagationAngleStep_Polar:90;
        end
        if  a.Quantity_Polar == 1
            DispersionDiagram_PhaseVelocity_Polar(a.Hybrid_Polar,a.LayupStringTex_Polar,a.PropagationAngleStep_Polar,a.BulkVelocities_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
        elseif a.Quantity_Polar == 2
            DispersionDiagram_EnergyVelocity_Polar(3,a.Hybrid_Polar,a.LayupStringTex_Polar,a.PropagationAngleStep_Polar,a.BulkVelocities_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
            DispersionDiagram_EnergyVelocity_Polar(2,a.Hybrid_Polar,a.LayupStringTex_Polar,a.PropagationAngleStep_Polar,a.BulkVelocities_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)           
            DispersionDiagram_EnergyVelocity_Polar(1,a.Hybrid_Polar,a.LayupStringTex_Polar,a.PropagationAngleStep_Polar,a.BulkVelocities_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
        elseif a.Quantity_Polar == 3
            DispersionDiagram_PropagationTime_Polar(2,a.Hybrid_Polar,a.LayupStringTex_Polar,a.Distance_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
            DispersionDiagram_PropagationTime_Polar(1,a.Hybrid_Polar,a.LayupStringTex_Polar,a.Distance_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
        elseif a.Quantity_Polar == 4
            DispersionDiagram_CoincidenceAngle_Polar(a.Hybrid_Polar,a.LayupStringTex_Polar,a.Couplant_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
        elseif a.Quantity_Polar == 5
            DispersionDiagram_Wavelength_Polar(a.Hybrid_Polar,a.LayupStringTex_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
        elseif a.Quantity_Polar == 6
            DispersionDiagram_Wavenumber_Polar(a.Hybrid_Polar,a.LayupStringTex_Polar,a.A_Polar,a.A0_Polar,a.Directory4,a.ExportPlots_Polar,a.FileName_Polar,a.AxesTickFontSize_Polar,a.TitleFontSize_Polar,a.ModeLabelFontSize_Polar,a.Frequency_Polar,FrequencyRange,a.Title_Polar,a.LineWidth_Polar,a.Material_Polar,a.PDF_Polar,a.PlateThickness_Polar/1e3,a.PNG_Polar,a.PNGresolution_Polar,PropagationAngle,a.PropagationAngleMode_Polar,a.S0_Polar,a.SH0_Polar,a.SuperLayerSize_Polar,a.SymmetricSystem_Polar)
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
elseif strcmp(source.Tag,'14') % Export
    a.ExportPlots_Polar = source.Value;
    if  a.ExportPlots_Polar
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
    try
        Export_Polar(a,1)
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
elseif strcmp(source.Tag,'19') % Excel
    try
        Export_Polar(a,2)
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
elseif strcmp(source.Tag,'20') % TXT
    try
        Export_Polar(a,3)
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
elseif strcmp(source.Tag,'21') % File name
    a.FileName_Polar = source.String;
elseif strcmp(source.Tag,'22') % Directory
    a.Directory4 = source.String;
elseif strcmp(source.Tag,'23') % Title
    if  source.Value == 1
        a.Title_Polar = 2;
    elseif source.Value == 2
        a.Title_Polar = 1;
    elseif source.Value == 3
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