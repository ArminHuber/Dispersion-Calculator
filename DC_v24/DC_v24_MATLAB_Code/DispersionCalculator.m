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
function DispersionCalculator
Version = '2.4';
Date = 'April 19, 2023';
Mode = 1; % 1: DC executed in MATLAB 2: DC executed as stand-alone; this switches the default directory for saving plots and data (see lines 39ff)

%#ok<*STRNU>
%#ok<*AGROW>
%#ok<*JAVFM>
%#ok<*GVMIS>

warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:illConditionedMatrix')
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
warning('off','MATLAB:ui:javaframe:PropertyToBeRemoved')

Dir = which('DispersionCalculator'); % get the full path and file name of "DispersionCalculator.m"
a.MaterialListDirectory = [Dir(1:length(Dir)-23),filesep,'Materials']; % get the path where the materials are stored
if  Mode == 1 % default directory if DC is executed in MATLAB
    a.Directory = Dir(1:length(Dir)-24-length(extractAfter(Dir(1:length(Dir)-23),asManyOfPattern(wildcardPattern+filesep))));
elseif Mode == 2 % default directory if DC is executed as stand-alone
    a.Directory = [extractBefore(Dir,filesep),filesep];
end

f1 = figure('NumberTitle','off','Name',['Dispersion Calculator v',Version],'Visible','off','MenuBar','none','Position',[0 0 1198 800],'CloseRequestFcn',@CloseRequest); % generate the main GUI
m1 = uimenu(f1,'Text','File'); % generate the menu bar
uimenu(m1,'Text','Open project','MenuSelectedFcn',@Open_Callback)
uimenu(m1,'Text','Save project','MenuSelectedFcn',@Save_Callback)
m2 = uimenu(f1,'Text','Materials');
uimenu(m2,'Text','Import','MenuSelectedFcn',@Import_Callback)
uimenu(m2,'Text','Export','MenuSelectedFcn',@Export_Callback)
m3 = uimenu(f1,'Text','Multicore');
uimenu(m3,'Text','Enable','MenuSelectedFcn',@Enable_Callback)
uimenu(m3,'Text','Disable','MenuSelectedFcn',@Disable_Callback)
uimenu(f1,'Text','Help','MenuSelectedFcn',@Help_Callback)
uimenu(f1,'Text','About','MenuSelectedFcn',@About_Callback)

jframe = get(gcf,'javaframe');
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png')))) % replace MATLAB logo by DC logo

tgroup = uitabgroup('Parent',f1); % generate the tab bar
Tab1 = uitab('Parent',tgroup,'Title','Isotropic');
Tab2 = uitab('Parent',tgroup,'Title','Anisotropic');
Tab3 = uitab('Parent',tgroup,'Title','Signal simulator');
Tab4 = uitab('Parent',tgroup,'Title','Polar diagrams');
Tab8 = uitab('Parent',tgroup,'Title','Bulk waves');
Tab6 = uitab('Parent',tgroup,'Title','Laminate stiffness');
Tab5 = uitab('Parent',tgroup,'Title','Material editor');
Tab7 = uitab('Parent',tgroup,'Title','Advanced','Backgroundcolor','white');

a.Materials = MaterialList; % load materials from "MaterialList_Fluid.txt", "MaterialList_Isotropic.txt", "MaterialList_TransverselyIsotropic.txt", and "MaterialList_Orthotropic.txt";

global Stop % global variable used to stop a calculation (if Stop == 1)
Stop = 0;

if  strcmp(Date(1:3),'Jan') % release string length
    if  length(Date) == 15
        ReleaseStringLength = 254;
    elseif length(Date) == 16
        ReleaseStringLength = 262;
    end    
elseif strcmp(Date(1:3),'Feb')
    if  length(Date) == 16
        ReleaseStringLength = 261;
    elseif length(Date) == 17
        ReleaseStringLength = 269;
    end
elseif strcmp(Date(1:3),'Mar')
    if  length(Date) == 13
        ReleaseStringLength = 243;
    elseif length(Date) == 14
        ReleaseStringLength = 251;
    end    
elseif strcmp(Date(1:3),'Apr')
    if  length(Date) == 13
        ReleaseStringLength = 231;
    elseif length(Date) == 14
        ReleaseStringLength = 239;
    end    
elseif strcmp(Date(1:3),'May')
    if  length(Date) == 11
        ReleaseStringLength = 229;
    elseif length(Date) == 12
        ReleaseStringLength = 237;
    end    
elseif strcmp(Date(1:3),'Jun')
    if  length(Date) == 12
        ReleaseStringLength = 234;
    elseif length(Date) == 13
        ReleaseStringLength = 242;
    end    
elseif strcmp(Date(1:3),'Jul')
    if  length(Date) == 12
        ReleaseStringLength = 228;
    elseif length(Date) == 13
        ReleaseStringLength = 236;
    end    
elseif strcmp(Date(1:3),'Aug')
    if  length(Date) == 14
        ReleaseStringLength = 248;
    elseif length(Date) == 15
        ReleaseStringLength = 256;
    end    
elseif strcmp(Date(1:3),'Sep')
    if  length(Date) == 17
        ReleaseStringLength = 275;
    elseif length(Date) == 18
        ReleaseStringLength = 283;
    end    
elseif strcmp(Date(1:3),'Oct')
    if  length(Date) == 15
        ReleaseStringLength = 256;
    elseif length(Date) == 16
        ReleaseStringLength = 264;
    end    
elseif strcmp(Date(1:3),'Nov')
    if  length(Date) == 16
        ReleaseStringLength = 270;
    elseif length(Date) == 17
        ReleaseStringLength = 278;
    end
elseif strcmp(Date(1:3),'Dec')
    if  length(Date) == 16
        ReleaseStringLength = 272;
    elseif length(Date) == 17
        ReleaseStringLength = 280;
    end    
end
if  length(Version) == 5
    ReleaseStringLength = ReleaseStringLength+12;
end

if  1 % advanced settings
    a.Multithreading = 0;
    
    % isotropic -----------------------------------------------------------
    a.YRange = 4; % determines how high the dispersion diagram is by default
    a.XRange1 = 2; % determines how broad the dispersion diagram is by default
    a.XSamples1 = 1e3; % determines how many samples a dispersion curve consists of by default
    a.Steps1 = 1e4; % determines how many samples the frequency sweep for the higher order modes contains

    a.FrequencyRangeStart1 = 1e-3; % start computation @ x kHz
    a.PhaseVelocityResolution1 = 1e-6; % (m/s)
    
    a.PhaseVelocitySections1 = 6;
    a.FrequencySections1 = 6;

    a.LambPhaseVelocitySweepRange11 = 2; % phase velocity search interval for Lamb waves for negative curvature; a higher value decreases the risk of missing a solution; a larger one might require a higher number of PhaseVelocitySections which increases processing time    
    a.LambPhaseVelocitySweepRange21 = 2; % for postive curvature

    a.PhaseVelocityStep1 = 100; % phase velocity steps for the frequency sweeps at high phase velocity to complete the modes (m/s) 
    a.FrequencyOffset1 = 20; % offset to low frequency for the frequency sweeps at high phase velocity for the modes (kHz/mm); its absolute must be large enough to find any solution but small enough to not reach a lower mode on the left; a larger SStep requires also a larger SOffset
    
    a.MissingSamples1 = 5; % if more than that number of sample points are missing, the tracing of a dispersion curve stops
    a.BelowCutoffWidth1 = 50;  % (kHz/mm); if we exceed the allowed scanning width (in kHz/mm) below the cut-off frequency of a damped mode without finding it, the damped search stops, and only the nondamped tracing continues to have that data for the next mode; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one
    
        % fluid loading and viscoelasticity settings
        a.FluidDensityThreshold1 = 200; % (g/cm3); above this value, zero attenuation solutions for A0 are discarded (in "Computer_Isotropic_ALamb_F.m") because these belong to the A_Scholte mode; below this value, A0 is very close to zero attenuation, and we don't discard the zero attenuation solutions anymore

        a.SearchWidthReal1 = [1.5 -1.5];
        a.SearchWidthImag1 = [2 .5];
        a.SearchAreaSections1 = 6;
        a.SearchAreaExtensions1 = 6;

    % anisotropic ----------------------------------------------
    a.XRange2 = 2; % determines how broad the dispersion diagram is by default (MHz*mm)
    a.XSamples2 = 200; % determines how many samples a dispersion curve consists of by default
    a.Steps2 = 2e3; % determines how many samples the frequency sweep for the higher order modes contains
    
    a.FrequencyRangeStart2 = 1e-3; % start computation @ x kHz 
    a.PhaseVelocityResolution2 = 1e-2; % (m/s)
    
    a.PhaseVelocitySections2 = 5;
    a.FrequencySections2 = 6;
    
    a.LambPhaseVelocitySweepRange12 = 10; % phase velocity search interval for Lamb waves for negative curvature; a higher value decreases the risk of missing a solution; a larger one might require a higher number of PhaseVelocitySections which increases processing time    
    a.LambPhaseVelocitySweepRange22 = 5; % for postive curvature
    a.ShearPhaseVelocitySweepRange2 = 2; % phase velocity search interval for shear horizontal waves; they have always positive curvature
    
    a.PhaseVelocityStep2 = 100; % phase velocity steps for the frequency sweeps at high phase velocity to complete the modes (m/s) 
    a.FrequencyOffset2 = 20; % offset to low frequency for the frequency sweeps at high phase velocity for the modes (kHz/mm); its absolute must be large enough to find any solution but small enough to not reach a lower mode on the left; a larger SStep requires also a larger SOffset
    
    a.MissingSamples2 = 5; % if more than that number of sample points are missing, the tracing of a dispersion curve stops
    a.BelowCutoffWidth2 = 50;  % (kHz/mm); if we exceed the allowed scanning width (in kHz/mm) below the cut-off frequency of a damped mode without finding it, the damped search stops, and only the nondamped tracing continues to have that data for the next mode; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one
    
        % fluid loading and viscoelasticity settings
        a.FluidDensityThreshold2 = 200; % (g/cm3); above this value, zero attenuation solutions for A0 are discarded (in "Computer_Isotropic_ALamb_F.m") because these belong to the A_Scholte mode; below this value, A0 is very close to zero attenuation, and we don't discard the zero attenuation solutions anymore

        a.SearchWidthReal2 = [1.5 -1.5];
        a.SearchWidthImag2 = [2 .5];
        a.SearchAreaSections2 = 5;
        a.SearchAreaExtensions2 = 3;

    % polar diagrams ----------------------------------------------
    a.XRange_Polar = .5; % determines how broad the dispersion diagram is by default (MHz*mm)
    a.XSamples_Polar = 50; % determines how many samples a dispersion curve consists of by default 
    a.MissingSamples_Polar = 5; % if more than that number of sample points are missing, the tracing of a dispersion curve stops
end
if  1 % default settings
    if  1 % Tab1_isotropic
        a.Directory1 = a.Directory;
        
        a.Detect1 = 0; % for isotropic tab: variable in order to check whether the higher order mode detection has been performed before their calculation       
        
        a.HSLamb1 = [];
        a.HSShear1 = [];
        a.HSScholte1 = [];
        a.HALamb1 = [];
        a.HAShear1 = [];
        a.HAScholte1 = [];

        a.FluidLoading1 = 0;
        E = fieldnames(a.Materials.Fluid);
        a.Fluid1 = getfield(a.Materials.Fluid,E{1});
        E = fieldnames(a.Materials.Isotropic);
        a.Material1 = getfield(a.Materials.Isotropic,E{1});
        if  isreal(a.Material1.C)
            a.Viscoelastic1 = 0;
        else
            a.Viscoelastic1 = 1;
        end

        a.Thickness = 1;

        a.PhaseVelocityLimit1 = round(a.YRange*a.Material1.PlateVelocity,-3);
        a.FrequencyLimit1 = round(a.Material1.PlateVelocity*a.XRange1/a.Thickness,-2);
        if  a.FrequencyLimit1 >= 1e4
            a.FrequencyLimit1 = round(a.FrequencyLimit1,-3);
        end
        a.FrequencyResolution1 = a.FrequencyLimit1/a.XSamples1;

        a.HigherOrderModes1 = 1;
        a.SymmetricModes1 = 1;
        a.AntisymmetricModes1 = 1;
        a.LambModes1 = 1;
        a.ShearHorizontalModes1 = 0;
        a.ScholteModes1 = 1;

        a.Step1 = a.FrequencyLimit1/a.Steps1;

        a.Quantity11 = 1;
        a.BulkVelocities1 = 0;
        a.Distance1 = 1e2;
        
        E = fieldnames(a.Materials.Fluid);
        a.Couplant1 = getfield(a.Materials.Fluid,E{1});

        a.XAxisMode1 = 1;
        a.XAxis1 = [0 a.FrequencyLimit1];
        a.YAxis1 = [0 a.PhaseVelocityLimit1/1e3];

        a.Quantity21 = 1;
        a.Frequency11 = a.FrequencyLimit1;
        a.Samples11 = 500;
        a.HalfspacesNumber11 = 1;
        a.Halfspaces11 = 0;
        a.Phase1 = 0;
        a.Plot1 = [1 1 1 1 1 1]; % displacement components - - 3 2 1 - / stress components 11 22 33 23 13 12

        a.ModeShapeSettingChanged1 = 1;        
        a.Frequency21 = a.FrequencyLimit1;
        a.Length1 = 2.5;
        a.Samples21 = 80;
        a.Samples31 = .5*a.Samples21;
        a.Scale1 = 1;
        a.GridLine1 = 2;
        a.Undistorted1 = 0;
        a.HalfspacesNumber21 = 1;
        a.Halfspaces21 = 0;

        a.Cycles1 = 1;
        a.CycleDuration1 = 1.5;
        a.FrameRate1 = 30;
        a.MovieQuality1 = 75;
        a.Animate1 = 0;

        a.ExportPlots1 = 0;
        a.CropPlots1 = 0;
        a.PDF1 = 1;
        a.PNG1 = 0;
        a.PNGresolution1 = 150;
        a.XAxisMode21 = 1;
        a.Arrange1 = 1;
        a.DispersionCurves1 = 1;
        a.ThroughThickness1 = 0;
        a.FileName1 = 'File name';

        a.Title1 = 1;
        a.ModeLabels1 = 1;
        a.LegendLocation1 = 'out';
        a.BoxLineWidth1 = .5;
        a.LineWidth1 = 1;
        a.SColor1 = [1 0 0];
        a.AColor1 = [0 0 1];
        a.S0Label = .05;
        a.SH0Label1 = .05;
        a.A0Label = .05;
        a.TitleFontSize1 = 30;
        a.AxesLabelFontSize1 = 30;
        a.AxesTickFontSize1 = 24;
        a.ModeLabelFontSize1 = 24;
        a.LegendFontSize1 = 24;

        a.MaterialName4A = a.Material1.Name;
    end
    if  1 % Tab2_anisotropic
        a.Directory2 = a.Directory;
        
        a.Detect2 = 0; % for anisotropic tab
        
        a.HS = [];
        a.HSLamb2 = [];
        a.HSShear2 = [];
        a.HSScholte2 = [];
        a.HA = [];
        a.HALamb2 = [];
        a.HAShear2 = [];
        a.HAScholte2 = [];
        a.HB = [];
        a.HBLamb = [];
        a.HBShear = [];
        a.HBScholte2 = [];

        a.FluidLoading2 = 0;
        a.ToggleUpperFluid2 = 0;
        a.ToggleLowerFluid2 = 0;
        a.SelectUpperFluidUI2Value = 1;
        a.SelectLowerFluidUI2Value = 1;
        a.SelectUpperFluidUI2Enable = 'off';
        a.SelectLowerFluidUI2Enable = 'off';
        E = fieldnames(a.Materials.Fluid);
        a.UpperFluid2 = getfield(a.Materials.Fluid,E{1});
        a.LowerFluid2 = getfield(a.Materials.Fluid,E{1});
        a.MaterialType2 = 1;
        a.MaterialTypeUI2Enable = 'on';
        E = fieldnames(a.Materials.Orthotropic);
        a.Material2{1} = getfield(a.Materials.Orthotropic,E{1});
        if  isreal(a.Material2{1}.C)
            a.Viscoelastic2 = 0;
        else
            a.Viscoelastic2 = 1;
        end
        a.MaterialNames{1} = a.Material2{1}.Name;
        a.MaterialClasses{1} = a.Material2{1}.Class;
        a.MaterialUI2Enable = 'on';
        
        a.MaterialUI2Value = 1;
        a.Hybrid = 0;
        a.UniformLayerThickness = 1;
        a.PlateThicknessUI2Enable = 'on';
        a.TablesUI2ColumnEditable = [true false false false false false true];
        a.PlateThickness = 1;
        a.SuperLayers = 1;
        a.SymmetricSystem = 0;
        a.Symmetric = 1;
        a.UnitCell = cell(400,6);
        a.UnitCell{1} = '0';
        a.UnitCell{1,2} = '1';
        a.UnitCell{1,3} = E{1};
        a.LayerOrientations = 0;
        a.LayerThicknesses = 1;
        
        a.LayupString1 = '0';
        a.EffectiveLayupString1 = '0';        

        a.PropagationAngle = 0;
        a.PhaseVelocityLimit2 = 20000;
        a.FrequencyLimit2 = 2000;
        a.FrequencyResolution2 = a.FrequencyLimit2/a.XSamples2;

        a.HigherOrderModes2 = 1;
        a.SymmetricModes2 = 1;
        a.AntisymmetricModes2 = 1;
        a.LambModes2 = 1;
        a.ShearHorizontalModes2 = 0;
        a.ScholteModes2 = 1;

        a.Step2 = a.FrequencyLimit2/a.Steps2;

        a.Quantity12 = 1;
        a.BulkVelocities2 = 0;
        a.Distance2 = 1e2;
        
        E = fieldnames(a.Materials.Fluid);
        a.Couplant2 = getfield(a.Materials.Fluid,E{1});
        
        a.XAxisMode2 = 1;
        a.XAxis2 = [0 a.FrequencyLimit2];
        a.YAxis2 = [0 a.PhaseVelocityLimit2/1e3];

        a.Quantity22 = 1;
        a.Frequency12 = a.FrequencyLimit2;
        a.Samples12 = a.Samples11;
        a.HalfspacesNumber12 = 1;
        a.Halfspaces12 = 0;
        a.Phase2 = 0;
        a.Plot2 = a.Plot1;

        a.ModeShapeSettingChanged2 = 1;
        a.Frequency22 = a.FrequencyLimit2;
        a.Length2 = 2.5;
        a.Samples22 = a.Samples21;
        a.Samples32 = .5*a.Samples22;
        a.Scale2 = 1;
        a.GridLine2 = 2;
        a.Undistorted2 = 0;
        a.HalfspacesNumber22 = 1;
        a.Halfspaces22 = 0;

        a.Cycles2 = 1;
        a.CycleDuration2 = 1.5;
        a.FrameRate2 = 30;
        a.MovieQuality2 = 75;
        a.Animate2 = 0;

        a.ExportPlots2 = 0;
        a.CropPlots2 = 0;
        a.PDF2 = 1;
        a.PNG2 = 0;
        a.PNGresolution2 = 150;
        a.XAxisMode22 = 1;
        a.Arrange2 = 1;
        a.DispersionCurves2 = 1;
        a.ThroughThickness2 = 0;
        a.FileName2 = 'File name';

        a.Title2 = 2;
        a.ModeLabels2 = 1;
        a.LegendLocation2 = 'out';
        a.BoxLineWidth2 = .5;
        a.LineWidth2 = 1;
        a.SColor2 = [1 0 0];
        a.AColor2 = [0 0 1];
        a.BColor = [.5 0 1];
        a.S0B1Label = .05;
        a.SH0Label2 = .05;
        a.A0B0Label = .05;
        a.TitleFontSize2 = 30;
        a.AxesLabelFontSize2 = 30;
        a.AxesTickFontSize2 = 24;
        a.ModeLabelFontSize2 = 24;
        a.LegendFontSize2 = 24;

        a.MaterialName4B = a.Material2{1}.Name;
    end
    if  1 % Tab3_signal simulator
        a.Directory3 = a.Directory;
        
        a.Cycles3 = 10;
        a.SamplesPerCycle3 = 20;
        a.Window3 = 1;
        a.Distance3 = 100;
        a.TimeLimitFactor3 = 2;
        a.SpectrumThreshold3 = .5;
        a.DisplacementComponent3 = 1;

        a.Gate3 = [0 0];
        a.MultiMode3 = 0;
        
        a.Plot_ALamb3(1:10) = 0;
        a.Plot_SLamb3(1:10) = 0;
        a.Plot_AShear3(1:10) = 0;
        a.Plot_SShear3(1:10) = 0;
        a.Amplitude_ALamb3(1:10) = 1;
        a.Amplitude_SLamb3(1:10) = 1;
        a.Amplitude_AShear3(1:10) = 1;
        a.Amplitude_SShear3(1:10) = 1;
        
        a.LineColors3 = [0 0 1;... blue (1)
            1 0 0;... red (2)
            .13 .55 .13;... green (3)
            0 0 0;... black (4)
            1 0 1;... magenta (5)
            0 1 1;... cyan (6)
            1 .7 0;... orange (7)
            .55 .27 .13;... brown (8)
            .5 0 1;... violet (9)
            .5 .5 .5]; % gray (10)

        a.ExportPlots3 = 0;
        a.CropPlots3 = 0;
        a.PDF3 = 1;
        a.PNG3 = 0;
        a.PNGresolution3 = 150;
        a.FileName3 = 'File name';

        a.Title3 = 1;
        a.BoxLineWidth3 = .5;
        a.LineWidth3 = 1;
        a.TitleFontSize3 = 30;
        a.AxesLabelFontSize3 = 30;
        a.AxesTickFontSize3 = 24;
    end
    if  1 % Tab4_polar diagrams
        a.Directory4 = a.Directory;
        
        a.MaterialType_Polar = 1;
        a.MaterialTypeUI4Enable = 'on';
        E = fieldnames(a.Materials.Orthotropic);
        a.Material_Polar{1} = getfield(a.Materials.Orthotropic,E{1}); 
        a.MaterialNames_Polar{1} = a.Material_Polar{1}.Name;
        a.MaterialClasses_Polar{1} = a.Material_Polar{1}.Class;
        a.MaterialUI4Enable = 'on';

        a.MaterialUI4Value = 1;
        a.Hybrid_Polar = 0;
        a.UniformLayerThickness_Polar = 1;
        a.PlateThicknessUI4Enable = 'on';
        a.TablesUI4ColumnEditable = [true false false false false false true];
        a.PlateThickness_Polar = 1;
        a.SuperLayers_Polar = 1;
        a.SymmetricSystem_Polar = 0;
        a.UnitCell_Polar = cell(400,6);
        a.UnitCell_Polar{1} = '0';
        a.UnitCell_Polar{1,2} = '1';
        a.UnitCell_Polar{1,3} = E{1};
        a.LayerOrientations_Polar = 0;
        a.LayerThicknesses_Polar = 1;
        
        a.LayupString_Polar = '0';

        a.FrequencyLimit_Polar = 500;
        a.FrequencyResolution_Polar = a.FrequencyLimit_Polar/a.XSamples_Polar;
        a.PropagationAngleMode_Polar = 1;
        a.PropagationAngleStep_Polar = 3;
        a.PhaseVelocitySections_Polar = 3;

        a.S0_Polar = 1;
        a.SH0_Polar = 1;
        a.A0_Polar = 1;

        a.Quantity_Polar = 1;
        a.BulkVelocities_Polar = 0;
        a.Distance_Polar = 100;
        
        E = fieldnames(a.Materials.Fluid);
        a.Couplant_Polar = getfield(a.Materials.Fluid,E{1});
        
        a.Frequency_Polar = 500;

        a.ExportPlots_Polar = 0;
        a.CropPlots_Polar = 0;
        a.PDF_Polar = 1;
        a.PNG_Polar = 0;
        a.PNGresolution_Polar = 150;
        a.FileName_Polar = 'File name';

        a.Title_Polar = 2;
        a.LineWidth_Polar = 1;
        a.SColor_Polar = [1 0 0];
        a.AColor_Polar = [0 0 1];
        a.BColor_Polar = [.5 0 1];
        a.TitleFontSize_Polar = 30;
        a.AxesTickFontSize_Polar = 24;
        a.ModeLabelFontSize_Polar = 24;
    end
    if  1 % Tab8_bulk waves
        a.Directory8 = a.Directory;
        
        a.MaterialType_Bulk = 1;
        E = fieldnames(a.Materials.Orthotropic);
        a.Material_Bulk = getfield(a.Materials.Orthotropic,E{1});
        a.Quantity_Bulk = 1;        
        a.ThetaStep_Bulk1 = 1;
        
        a.Plane = 13;
        a.Phi_Bulk11 = 0;
        
        a.PhiStep_Bulk = 1;
        a.Mode_Bulk = 'L';
        a.ViewPhi_Bulk1 = -37.5;
        a.ViewTheta_Bulk1 = 30;
        a.MarkerSize_Bulk = 20;
        a.ColorbarX_Bulk = .9;
        
        a.Phi_Bulk21 = 0;
        a.Theta_Bulk1 = 0;
        
        E = fieldnames(a.Materials.Fluid);
        a.Couplant_Bulk = getfield(a.Materials.Fluid,E{1});
        
        a.SolidType_Bulk = 1;
        E = fieldnames(a.Materials.Isotropic);
        a.Solid_Bulk = getfield(a.Materials.Isotropic,E{1});        
        a.Phi_Bulk2 = 0;
        a.Theta_Bulk2 = 2;
        a.ThetaStep_Bulk2 = .1;
        a.ViewPhi_Bulk2 = -37.5;
        a.ViewTheta_Bulk2 = 30;

        a.ExportPlots_Bulk = 0;
        a.CropPlots_Bulk = 0;
        a.PDF_Bulk = 1;
        a.PNG_Bulk = 0;
        a.PNGresolution_Bulk = 150;
        a.FileName_Bulk = 'File name';
        
        a.Title_Bulk = 1;
        a.BoxLineWidth_Bulk = .5;
        a.LineWidth_Bulk = 1;
        a.WaveVectorLineWidth_Bulk = 2;
        a.TitleFontSize_Bulk = 30;
        a.AxesLabelFontSize_Bulk = 30;
        a.AxesTickFontSize_Bulk = 24;
        a.ModeLabelFontSize_Bulk = 24;       
    end
    if  1 % Tab6_laminate stiffness
        a.Directory6 = a.Directory;
        
        a.MaterialType6 = 1;
        a.MaterialTypeUI6Enable = 'on';
        E = fieldnames(a.Materials.Orthotropic);
        a.Material3{1} = getfield(a.Materials.Orthotropic,E{1});
        a.MaterialNames3{1} = a.Material3{1}.Name;
        a.MaterialClasses3{1} = a.Material3{1}.Class;
        a.MaterialUI6Enable = 'on';

        a.MaterialUI6Value = 1;
        a.Hybrid6 = 0;
        a.UniformLayerThickness3 = 1;
        a.TablesUI6ColumnEditable = [true false false false false false true];
        a.UnitCell3 = cell(400,6);   
        a.UnitCell3{1} = '0';
        a.UnitCell3{1,2} = '1';
        a.UnitCell3{1,3} = E{1};
        a.LayerOrientations3 = 0;
        a.LayerThicknesses3 = 1;

        a.PropagationAngle3 = 0;

        a.Polar = [1 0 0 0 1 0 0 0 0 0 0 0 0];
    end
    if  1 % Tab5_material editor
        a.AttenuationUnit1 = 1;
        a.AtFrequency1 = 1000; % kHz
        a.MaterialTypeME = 1;
        a.MaterialName4C = a.Couplant1.Name;        
        a.Material1ME = a.Material1;
        a.Material2ME = a.Material2{1};
        a.Material3ME = a.Couplant1;
    end
end
if  1 % populate tabs
    if  1 % Tab1_isotropic
        uicontrol('Parent',Tab1,'Style','slider','Value',1,'Position',[10 10 25 150],'Callback',@SliderUI1_Callback);
        
        a.p1UI1 = uipanel('Parent',Tab1,'Title','Specimen','Units','pixels','Position',[10 610 225 145],'FontSize',10);
        uicontrol('Parent',a.p1UI1,'Style','text','String','Fluid-loading','Position',[10 105 62 13]);
        uicontrol('Parent',a.p1UI1,'Style','text','String','Fluid','Position',[10 75 24 13]);
        uicontrol('Parent',a.p1UI1,'Style','text','String','Material','Position',[10 45 39 13]);
        uicontrol('Parent',a.p1UI1,'Style','text','String','Thickness (mm)','Position',[10 15 78 13]);
        uicontrol('Parent',a.p1UI1,'Style','checkbox','Value',a.FluidLoading1,'TooltipString','Check this to enable fluid-loading.','Position',[160 100 50 23],'Callback',@CallbackUI1,'Tag','6');
        a.FluidUI1 = uicontrol('Parent',a.p1UI1,'Style','popupmenu','String',fieldnames(a.Materials.Fluid),'TooltipString','Select a fluid into which the plate is immersed.','Position',[60 70 150 23],'Enable','off','Callback',@CallbackUI1,'Tag','7');
        a.MaterialUI1 = uicontrol('Parent',a.p1UI1,'Style','popupmenu','String',fieldnames(a.Materials.Isotropic),'TooltipString','Select a material.','Position',[60 40 150 23],'Callback',@CallbackUI1,'Tag','1');
        uicontrol('Parent',a.p1UI1,'Style','edit','String',a.Thickness,'TooltipString','Enter the plate''s thickness.','Position',[160 10 50 23],'Callback',@CallbackUI1,'Tag','2');

        a.p2UI1 = uipanel('Parent',Tab1,'Title','Computational settings','Units','pixels','Position',[10 490 225 115],'FontSize',10);
        uicontrol('Parent',a.p2UI1,'Style','text','String','Phase velocity limit (m/ms)','Position',[10 75 128 13]);
        uicontrol('Parent',a.p2UI1,'Style','text','String','Frequency limit (kHz)','Position',[10 45 103 13]);
        uicontrol('Parent',a.p2UI1,'Style','text','String','Frequency step (kHz)','Position',[10 15 107 13]);
        a.PhaseVelocityLimitUI1 = uicontrol('Parent',a.p2UI1,'Style','edit','String',a.PhaseVelocityLimit1/1e3,'TooltipString',['Enter the phase velocity (Y-axis in the dispersion diagram)',newline,'up to which the dispersion curves shall be traced.'],'Position',[160 70 50 23],'Callback',@CallbackUI1,'Tag','3');
        a.FrequencyLimitUI1 = uicontrol('Parent',a.p2UI1,'Style','edit','String',a.FrequencyLimit1,'TooltipString',['Enter the frequency (X-axis in the dispersion diagram)',newline,'up to which the dispersion curves shall be traced.'],'Position',[160 40 50 23],'Callback',@CallbackUI1,'Tag','4');
        a.FrequencyResolutionUI1 = uicontrol('Parent',a.p2UI1,'Style','edit','String',a.FrequencyResolution1,'TooltipString','Enter the frequency step.','Position',[160 10 50 23],'Callback',@CallbackUI1,'Tag','5');

        a.p3UI1 = uipanel('Parent',Tab1,'Title','Mode selection','Units','pixels','Position',[10 280 225 205],'FontSize',10);
        uicontrol('Parent',a.p3UI1,'Style','text','String','Higher order modes','Position',[10 165 97 13]);
        uicontrol('Parent',a.p3UI1,'Style','text','String','Symmetric modes','Position',[10 135 87 13]);
        uicontrol('Parent',a.p3UI1,'Style','text','String','Antisymmetric modes','Position',[10 105 105 13]);
        uicontrol('Parent',a.p3UI1,'Style','text','String','Lamb modes','Position',[10 75 63 13]);
        uicontrol('Parent',a.p3UI1,'Style','text','String','Shear horizontal modes','Position',[10 45 117 13]);
        uicontrol('Parent',a.p3UI1,'Style','text','String','Scholte modes','Position',[10 15 73 13]);
        uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.HigherOrderModes1,'TooltipString','Check this in order to calculate the higher order modes in addition to the fundamental ones.','Position',[160 160 50 23],'Callback',@CallbackUI1,'Tag','8');
        uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.SymmetricModes1,'TooltipString',['Check this in order to calculate the symmetric modes. These modes have a',newline,'symmetric displacement pattern with respect to the middle plane of the plate.'],'Position',[160 130 50 23],'Callback',@CallbackUI1,'Tag','9');
        uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.AntisymmetricModes1,'TooltipString',['Check this in order to calculate the antisymmetric modes. These modes have an',newline,'antisymmetric displacement pattern with respect to the middle plane of the plate.'],'Position',[160 100 50 23],'Callback',@CallbackUI1,'Tag','10');
        uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.LambModes1,'TooltipString',['Check this in order to calculate the Lamb wave modes. In isotropic media, these modes',newline,'show displacement only in the sagittal plane spanned by the propagation direction x1',newline,'and by the out-of-plane direction x3. These kind of waves are termed ''pure'' Lamb',newline,'waves. Lamb waves are indicated by solid lines in the dispersion diagram.'],'Position',[160 70 50 23],'Callback',@CallbackUI1,'Tag','11');
        uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.ShearHorizontalModes1,'TooltipString',['Check this in order to calculate the shear horizontal modes. In isotropic media, these modes',newline,'show displacement only perpendicular (x2) to the propagation direction x1 and are therefore',newline,'termed ''pure'' modes. Shear horizontal waves are indicated by dashed lines in the dispersion',newline,'diagram.'],'Position',[160 40 50 23],'Callback',@CallbackUI1,'Tag','12');
        a.ScholteModesUI1 = uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.ScholteModes1,'TooltipString',['Check this in order to calculate the Scholte wave modes. Scholte waves are surface waves propagating',newline,'at an interface between a fluid and a solid. They are solved by the same characteristic function as the',newline,'Lamb waves, and similarly have displacement in the sagittal plane x1-x3. Scholte waves are indicated',newline,'by dashed-dotted lines in the dispersion diagram.'],'Position',[160 10 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','79');

        %------------------------------------------------------------------
        a.p4UI1 = uipanel('Parent',Tab1,'Title','Manually detect higher order modes','Units','pixels','Position',[245 690 223 65],'FontSize',10);
        uicontrol('Parent',a.p4UI1,'Style','text','String','Step (kHz)','Position',[10 20 53 13]);
        a.SteptUI1 = uicontrol('Parent',a.p4UI1,'Style','edit','String',a.Step1,'TooltipString',['Enter the step size for the frequency sweep. In general, a finer step increases the chance to',newline,'find mode cut-off frequencies, although exceptions from that rule might occur occasionally.',newline,'It is essential that all higher order modes are detected. Therefore, try smaller and larger step',newline,'sizes if you are not sure that you have found all modes.'],'Position',[75 15 50 23],'Callback',@CallbackUI1,'Tag','13');
        uicontrol('Parent',a.p4UI1,'Style','pushbutton','String','Detect','TooltipString',['Before you can trace the higher order modes, their cut-off frequencies',newline,'at the phase velocity limit must be detected. This is done automatically',newline,'upon pressing ''calculate'' if it has not already been done manually.',newline,'Sometimes the automatic search does not find all modes, then change',newline,'''Step'' and press ''Detect'' to find the missing modes.'],'Position',[140 10 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','14');

        a.OutputWindow1aUI1 = uicontrol('Parent',Tab1,'Style','text','String','','Position',[245 50 105 630],'BackgroundColor','white');
        a.OutputWindow1bUI1 = uicontrol('Parent',Tab1,'Style','text','String','','Position',[362 50 105 630],'BackgroundColor','white');
        a.OutputWindow2aUI1 = uicontrol('Parent',Tab1,'Style','text','String','','Position',[245 10 105 28],'BackgroundColor','white');
        a.OutputWindow2bUI1 = uicontrol('Parent',Tab1,'Style','text','String','','Position',[362 10 105 28],'BackgroundColor','white');

        %------------------------------------------------------------------
        a.c1UI1 = uicontrol('Parent',Tab1,'Style','pushbutton','String','Calculate','TooltipString','Start the dispersion curve tracing.','Position',[478 702 127 40],'FontSize',10,'Callback',@CallbackUI1,'Tag','15');
        a.c2UI1 = uicontrol('Parent',Tab1,'Style','pushbutton','String','Stop calculation','TooltipString','Stop the dispersion curve tracing.','Position',[615 702 127 40],'FontSize',10,'Callback',@CallbackUI1,'Tag','16');

        a.p5UI1 = uipanel('Parent',Tab1,'Title','Dispersion diagrams','Units','pixels','Position',[478 465 263 225],'FontSize',10);
        uicontrol('Parent',a.p5UI1,'Style','text','String','Quantity','Position',[10 185 42 13]);
        a.Option1TextUI1 = uicontrol('Parent',a.p5UI1,'Style','text','String','Bulk velocities','Position',[10 155 70 13]);
        a.XAxisModeTextUI1 = uicontrol('Parent',a.p5UI1,'Style','text','String','X-axis mode','Position',[10 125 62 13]);
        a.XAxisTextUI1 = uicontrol('Parent',a.p5UI1,'Style','text','String','X-axis (kHz)','Position',[10 95 62 13]);
        a.YAxisTextUI1 = uicontrol('Parent',a.p5UI1,'Style','text','String','Y-axis (m/ms)','Position',[10 65 70 13]);
        uicontrol('Parent',a.p5UI1,'Style','popupmenu','String',{'Phase velocity (m/ms)','Energy velocity (m/ms)',['Propagation time (',char(181),'s)'],['Coincidence angle (',char(176),')'],'Wavelength (mm)','Wavenumber (rad/mm)','Attenuation (Np/m)'},'TooltipString','Select which quantity to plot in the dispersion diagram.','Position',[117 180 130 23],'Callback',@CallbackUI1,'Tag','17');
        a.Option1UI1 = uicontrol('Parent',a.p5UI1,'Style','checkbox','Value',a.BulkVelocities1,'TooltipString','Check this to show the bulk wave velocities.','Position',[117 150 50 23],'Callback',@CallbackUI1,'Tag','18');
        a.XAxisModeUI1 = uicontrol('Parent',a.p5UI1,'Style','popupmenu','String',{'Frequency (kHz)','Frequency (MHz)',['f',char(8901),'d (MHz',char(8901),'mm)']},'TooltipString','Select the frequency''s dimension on the X-axis.','Position',[117 120 130 23],'Callback',@CallbackUI1,'Tag','19');
        a.XAxisUI1 = uicontrol('Parent',a.p5UI1,'Style','edit','String',['[0 ',num2str(a.FrequencyLimit1),']'],'TooltipString','Enter which frequency range shall be plotted.','Position',[117 90 75 23],'Callback',@CallbackUI1,'Tag','20');
        a.YAxisUI1 = uicontrol('Parent',a.p5UI1,'Style','edit','String',['[0 ',num2str(a.PhaseVelocityLimit1/1e3),']'],'TooltipString','Enter which phase velocity range shall be plotted.','Position',[117 60 75 23],'Callback',@CallbackUI1,'Tag','21');
        a.Plot1UI1 = uicontrol('Parent',a.p5UI1,'Style','pushbutton','String','Plot','TooltipString',['Plot the dispersion diagram. If you have checked ''Export plots''',newline,'in the export settings, the plot will be exported automatically.'],'Position',[117 15 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','22');

        a.p6UI1 = uipanel('Parent',Tab1,'Title','Through-thickness profiles','Units','pixels','Position',[751 465 224 290],'FontSize',10);
        uicontrol('Parent',a.p6UI1,'Style','text','String','Quantity','Position',[10 245 42 13]);
        uicontrol('Parent',a.p6UI1,'Style','text','String','Mode','Position',[10 215 28 13]);
        uicontrol('Parent',a.p6UI1,'Style','text','String','Frequency (kHz)','Position',[10 185 83 13]);
        a.Option2TextUI1 = uicontrol('Parent',a.p6UI1,'Style','text','String','Samples x3','Position',[10 155 58 13]);
        uicontrol('Parent',a.p6UI1,'Style','text','String','Half-spaces','Position',[10 125 61 13]);
        uicontrol('Parent',a.p6UI1,'Style','text','String','Phase','Position',[10 95 32 13]);
        uicontrol('Parent',a.p6UI1,'Style','popupmenu','String',{'Displacement','Stress','Strain','Energy density','Power flow density'},'TooltipString','Select which through-thickness quantity to plot.','Position',[119 240 90 23],'Callback',@CallbackUI1,'Tag','23');
        a.Mode1UI1 = uicontrol('Parent',a.p6UI1,'Style','popupmenu','String',{''},'TooltipString','Select the mode you want to analyze.','Position',[119 210 90 23],'Callback',@CallbackUI1,'Tag','24');
        a.Frequency1UI1 = uicontrol('Parent',a.p6UI1,'Style','edit','String',a.FrequencyLimit1,'TooltipString','Enter the frequency at which to analyze the selected mode.','Position',[119 180 50 23],'Callback',@CallbackUI1,'Tag','25');
        uicontrol('Parent',a.p6UI1,'Style','edit','String',a.Samples11,'TooltipString',['Enter the number of sample points over the plate''s thickness (x3)',newline,'at which the selected quantities are calculated.'],'Position',[119 150 50 23],'Callback',@CallbackUI1,'Tag','26');
        a.Halfspaces1UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','Value',a.Halfspaces11,'TooltipString','Check this to show the quantities in the upper and lower fluid.','Position',[89 120 20 23],'Enable','off','Callback',@CallbackUI1,'Tag','81');        
        a.HalfspacesNumber1UI1 = uicontrol('Parent',a.p6UI1,'Style','edit','String',a.HalfspacesNumber11,'TooltipString','Set the height of the half-spaces in plate thicknesses.','Position',[119 120 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','80');
        a.PhaseUI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','Value',a.Phase1,'TooltipString','Check this to plot also the phase of the field components.','Position',[119 90 20 23],'Callback',@CallbackUI1,'Tag','84');
        a.Plot2UI1 = uicontrol('Parent',a.p6UI1,'Style','pushbutton','String','Plot','TooltipString',['Plot the profile. If you have checked ''Export plots'' in the',newline,'export settings, the plot will be exported automatically.'],'Position',[119 15 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','27');
        a.x11UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','11','Value',a.Plot1(1),'TooltipString','','Position',[10 50 35 23],'Enable','off','Callback',@CallbackUI1,'Tag','28');
        a.x22UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','22','Value',a.Plot1(2),'TooltipString','','Position',[45 30 35 23],'Enable','off','Callback',@CallbackUI1,'Tag','29');
        a.x33UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','3','Value',a.Plot1(3),'TooltipString','Displacement u3','Position',[80 10 40 23],'Callback',@CallbackUI1,'Tag','30');
        a.x13UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','1','Value',a.Plot1(4),'TooltipString','Displacement u1','Position',[80 50 40 23],'Callback',@CallbackUI1,'Tag','31');
        a.x23UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','2','Value',a.Plot1(5),'TooltipString','Displacement u2','Position',[80 30 40 23],'Callback',@CallbackUI1,'Tag','32');
        a.x12UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','12','Value',a.Plot1(6),'TooltipString','','Position',[45 50 35 23],'Enable','off','Callback',@CallbackUI1,'Tag','33');
        
        a.p7UI1 = uipanel('Parent',Tab1,'Title','Mode shape','Units','pixels','Position',[478 195 497 265],'FontSize',10);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Mode','Position',[10 225 28 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Frequency (kHz)','Position',[10 195 83 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Wavelengths','Position',[10 165 65 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Samples x1','Position',[10 135 58 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Samples x3','Position',[10 105 58 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Scale','Position',[10 75 29 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Grid line','Position',[10 45 41 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Undistorted','Position',[10 15 57 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','String','Half-spaces','Position',[195 45 61 13]);
        a.Mode2UI1 = uicontrol('Parent',a.p7UI1,'Style','popupmenu','String',{''},'TooltipString','Select the mode you want to analyze.','Position',[117 220 90 23],'Callback',@CallbackUI1,'Tag','34');        
        a.Frequency2UI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.FrequencyLimit1,'TooltipString','Enter the frequency at which to analyze the selected mode.','Position',[117 190 50 23],'Callback',@CallbackUI1,'Tag','35');
        uicontrol('Parent',a.p7UI1,'Style','edit','String',a.Length1,'TooltipString','Specify how many wavelengths you want to display.','Position',[117 160 50 23],'Callback',@CallbackUI1,'Tag','36');
        uicontrol('Parent',a.p7UI1,'Style','edit','String',a.Samples21,'TooltipString',['Enter the number of sample points along the propagation',newline,'direction (x1) at which the displacement is calculated.'],'Position',[117 130 50 23],'Callback',@CallbackUI1,'Tag','37');
        a.Samples3UI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.Samples31,'TooltipString',['Enter the number of sample points over the plate''s',newline,'thickness (x3) at which the displacement is calculated.'],'Position',[117 100 50 23],'Callback',@CallbackUI1,'Tag','38');
        uicontrol('Parent',a.p7UI1,'Style','edit','String',a.Scale1,'TooltipString','Change the scaling of the displacement.','Position',[117 70 50 23],'Callback',@CallbackUI1,'Tag','39');
        a.GridLineUI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.GridLine1,'TooltipString','Draw a grid line at every ith sample point.','Position',[117 40 50 23],'Callback',@CallbackUI1,'Tag','40');
        uicontrol('Parent',a.p7UI1,'Style','checkbox','Value',a.Undistorted1,'TooltipString','Check this to draw the undistorted grid.','Position',[117 10 50 23],'Callback',@CallbackUI1,'Tag','41');
        a.HalfspacesNumber2UI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.HalfspacesNumber21,'TooltipString','Set the height of the half-spaces in plate thicknesses.','Position',[272 40 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','82');
        a.Halfspaces2UI1 = uicontrol('Parent',a.p7UI1,'Style','checkbox','Value',a.Halfspaces21,'TooltipString','Check this to show the quantities in the upper and lower fluid.','Position',[272 10 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','83');
        a.Plot3UI1 = uicontrol('Parent',a.p7UI1,'Style','pushbutton','String','Plot','TooltipString',['Plot the mode shape. If you have checked ''Animate'',',newline,'the mode shape will be animated. If you have',newline,'checked ''Export plots'' in the export settings, the',newline,'plot/movie will be exported automatically.'],'Position',[392 15 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','42');

        a.p8UI1 = uipanel('Parent',Tab1,'Title','Animation settings','Units','pixels','Position',[751 265 212 175],'FontSize',9);
        uicontrol('Parent',a.p8UI1,'Style','text','String','Cycles','Position',[10 135 35 13]);
        uicontrol('Parent',a.p8UI1,'Style','text','String','Cycle duration (s)','Position',[10 105 88 13]);
        uicontrol('Parent',a.p8UI1,'Style','text','String','Frame rate (Hz)','Position',[10 75 78 13]);
        uicontrol('Parent',a.p8UI1,'Style','text','String','Movie quality (0-100)','Position',[10 45 103 13]);
        uicontrol('Parent',a.p8UI1,'Style','text','String','Animate','Position',[10 15 41 13]);
        uicontrol('Parent',a.p8UI1,'Style','edit','String',a.Cycles1,'TooltipString',['Enter how many cycles shall be calculated. That cylce',newline,'will be repeated until you close the plot figure.'],'Position',[119 130 50 23],'Callback',@CallbackUI1,'Tag','43');
        uicontrol('Parent',a.p8UI1,'Style','edit','String',a.CycleDuration1,'TooltipString','Enter how long a cycle shall take.','Position',[119 100 50 23],'Callback',@CallbackUI1,'Tag','44');
        uicontrol('Parent',a.p8UI1,'Style','edit','String',a.FrameRate1,'TooltipString','Define the frame rate of the movie.','Position',[119 70 50 23],'Callback',@CallbackUI1,'Tag','45');
        uicontrol('Parent',a.p8UI1,'Style','edit','String',a.MovieQuality1,'TooltipString','Define the quality of the exported movie.','Position',[119 40 50 23],'Callback',@CallbackUI1,'Tag','46');
        uicontrol('Parent',a.p8UI1,'Style','checkbox','Value',a.Animate1,'TooltipString',['Check this in order to show the animated mode',newline,'shape upon pressing the plot button below.'],'Position',[119 10 70 23],'Callback',@CallbackUI1,'Tag','47');
        
        a.p9UI1 = uipanel('Parent',Tab1,'Title','Export settings','Units','pixels','Position',[478 10 497 180],'FontSize',10);
        uicontrol('Parent',a.p9UI1,'Style','text','String','Export plots','Position',[10 140 59 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','String','Crop plots','Position',[10 120 51 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','String','PDF','Position',[124 140 21 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','String','PNG','Position',[124 120 23 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','String','PNG resolution (dpi)','Position',[10 93 98 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','String','Dispersion curves','Position',[240 140 90 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','String','Through-thickness','Position',[240 120 92 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','String','Directory','Position',[10 20 46 13]);
        uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.ExportPlots1,'TooltipString',['Check this in order to export plots/movies automatically upon',newline,'pressing the respective plot button. You can also export plots',newline,'manually by using the ''File'' menu inside the plot figure.'],'Position',[80 135 20 23],'Callback',@CallbackUI1,'Tag','48');
        uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.CropPlots1,'TooltipString','Plots will be cropped tightly with minimal white space.','Position',[80 115 20 23],'Callback',@CallbackUI1,'Tag','77');
        uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.PDF1,'TooltipString','Check this to export a plot as pdf.','Position',[160 135 20 23],'Callback',@CallbackUI1,'Tag','49');
        uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.PNG1,'TooltipString','Check this to export a plot as png.','Position',[160 115 20 23],'Callback',@CallbackUI1,'Tag','50');
        uicontrol('Parent',a.p9UI1,'Style','edit','String',a.PNGresolution1,'TooltipString','Enter the resolution of the png image.','Position',[125 88 50 23],'Callback',@CallbackUI1,'Tag','51');
        a.XAxisMode2UI1 = uicontrol('Parent',a.p9UI1,'Style','popupmenu','String',{'Frequency (kHz)','Frequency (MHz)',['f',char(8901),'d (MHz',char(8901),'mm)']},'TooltipString','Select the frequency''s dimension to be exported.','Position',[370 130 110 23],'Callback',@CallbackUI1,'Tag','52');
        a.ArrangeUI1 = uicontrol('Parent',a.p9UI1,'Style','popupmenu','String',{'Horizontal arrangement','Vertical arrangement'},'TooltipString','Choose how to arrange the dispersion curve data in the Excel/txt.','Position',[370 95 110 23],'Callback',@CallbackUI1,'Tag','53');        
        uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.DispersionCurves1,'TooltipString','Check this to export the dispersion curves.','Position',[340 135 20 23],'Callback',@CallbackUI1,'Tag','54');
        uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.ThroughThickness1,'TooltipString','Check this to export the through-thickness profiles selected above.','Position',[340 115 20 23],'Callback',@CallbackUI1,'Tag','55');
        uicontrol('Parent',a.p9UI1,'Style','pushbutton','String','*.mat','TooltipString','Export the data as Matlab''s mat-file.','Position',[220 80 40 33],'Callback',@CallbackUI1,'Tag','78');
        uicontrol('Parent',a.p9UI1,'Style','pushbutton','String','*.xlsx','TooltipString','Export the data as Excel sheet.','Position',[270 80 40 33],'Callback',@CallbackUI1,'Tag','56');
        uicontrol('Parent',a.p9UI1,'Style','pushbutton','String','*.txt','TooltipString','Export the data as txt-file.','Position',[320 80 40 33],'Callback',@CallbackUI1,'Tag','57');
        uicontrol('Parent',a.p9UI1,'Style','edit','String',a.FileName1,'TooltipString',['Specify the name of the plots/movies to be exported.',newline,'The dispersion curve raw data are named automatically.'],'Position',[70 45 412 23],'Callback',@CallbackUI1,'Tag','58');
        uicontrol('Parent',a.p9UI1,'Style','edit','String',a.Directory,'TooltipString',['Specify the directory to which plots, movies, and',newline,'the dispersion curve raw data shall be exported.'],'Position',[70 15 412 23],'Callback',@CallbackUI1,'Tag','59');

        %------------------------------------------------------------------
        a.p10UI1 = uipanel('Parent',Tab1,'Title','Plot layout settings','Units','pixels','Position',[985 90 200 600],'FontSize',10);
        uicontrol('Parent',a.p10UI1,'Style','text','String','Title','Position',[10 560 21 13]);
        uicontrol('Parent',a.p10UI1,'Style','text','String','Mode labels','Position',[10 530 59 13]);
        uicontrol('Parent',a.p10UI1,'Style','text','String','Legend location','Position',[10 500 78 13]);
        uicontrol('Parent',a.p10UI1,'Style','text','String','Box line width','Position',[10 470 70 13]);
        uicontrol('Parent',a.p10UI1,'Style','text','String','Curve line width','Position',[10 440 80 13]);
        a.TitleUI1 = uicontrol('Parent',a.p10UI1,'Style','checkbox','Value',a.Title1,'TooltipString','Check this in order to show the plot title.','Position',[105 555 20 23],'Callback',@CallbackUI1,'Tag','60');
        a.ModeLabelsUI1 = uicontrol('Parent',a.p10UI1,'Style','checkbox','Value',a.ModeLabels1,'TooltipString',['Check this in order to show the mode',newline,'labels of the fundamental modes.'],'Position',[105 525 20 23],'Callback',@CallbackUI1,'Tag','61');
        a.LegendLocationUI1 = uicontrol('Parent',a.p10UI1,'Style','popupmenu','String',{'outside','inside'},'TooltipString','Determine the legend location in the through-thickness plots.','Position',[105 495 80 23],'Callback',@CallbackUI1,'Tag','62');
        a.BoxLineWidthUI1 = uicontrol('Parent',a.p10UI1,'Style','edit','String',a.BoxLineWidth1,'TooltipString','Enter the box line width.','Position',[105 465 50 23],'Callback',@CallbackUI1,'Tag','63');
        a.LineWidthUI1 = uicontrol('Parent',a.p10UI1,'Style','edit','String',a.LineWidth1,'TooltipString','Enter the curve line width.','Position',[105 435 50 23],'Callback',@CallbackUI1,'Tag','64');

        a.p11UI1 = uipanel('Parent',Tab1,'Title','Dispersion curve colors [R G B]','Units','pixels','Position',[995 400 180 85],'FontSize',9); 
        uicontrol('Parent',a.p11UI1,'Style','text','String','S','Position',[10 45 7 13]);
        uicontrol('Parent',a.p11UI1,'Style','text','String','A','Position',[10 15 8 13]);
        a.SColorUI1 = uicontrol('Parent',a.p11UI1,'Style','edit','String','[1 0 0]','TooltipString',['Specify the color of symmetric modes. You can compose any color',newline,'by entering the corresponding RGB values. Set numbers from 0 to 1.'],'Position',[95 40 50 23],'Callback',@CallbackUI1,'Tag','66');
        a.AColorUI1 = uicontrol('Parent',a.p11UI1,'Style','edit','String','[0 0 1]','TooltipString','Specify the color of antisymmetric modes.','Position',[95 10 50 23],'Callback',@CallbackUI1,'Tag','67');

        a.p12UI1 = uipanel('Parent',Tab1,'Title','Mode labels x-position (0-1)','Units','pixels','Position',[995 280 180 115],'FontSize',9); 
        uicontrol('Parent',a.p12UI1,'Style','text','String','S0','Position',[10 75 15 13]);
        uicontrol('Parent',a.p12UI1,'Style','text','String','SSH0','Position',[10 45 29 13]);
        uicontrol('Parent',a.p12UI1,'Style','text','String','A0','Position',[10 15 16 13]);
        a.S0LabelUI1 = uicontrol('Parent',a.p12UI1,'Style','edit','String',a.S0Label,'TooltipString',['Change the relative X-position of the fundamental symmetric Lamb wave',newline,'mode label. 0 corresponds to the left, 1 to the right X-axis limit.'],'Position',[95 70 50 23],'Callback',@CallbackUI1,'Tag','68');
        a.SH0LabelUI1 = uicontrol('Parent',a.p12UI1,'Style','edit','String',a.SH0Label1,'TooltipString',['Change the relative X-position of the fundamental shear horizontal wave',newline,'mode label.'],'Position',[95 40 50 23],'Callback',@CallbackUI1,'Tag','69');
        a.A0LabelUI1 = uicontrol('Parent',a.p12UI1,'Style','edit','String',a.A0Label,'TooltipString',['Change the relative X-position of the fundamental antisymmetric Lamb',newline,'wave mode label.'],'Position',[95 10 50 23],'Callback',@CallbackUI1,'Tag','70');

        a.p13UI1 = uipanel('Parent',Tab1,'Title','Font size','Units','pixels','Position',[995 100 180 175],'FontSize',9);
        uicontrol('Parent',a.p13UI1,'Style','text','String','Title','Position',[10 135 21 13]);
        uicontrol('Parent',a.p13UI1,'Style','text','String','Axes labels','Position',[10 105 59 13]);
        uicontrol('Parent',a.p13UI1,'Style','text','String','Axes ticks','Position',[10 75 53 13]);
        uicontrol('Parent',a.p13UI1,'Style','text','String','Mode labels','Position',[10 45 59 13]);
        uicontrol('Parent',a.p13UI1,'Style','text','String','Legend','Position',[10 15 38 13]);
        a.TitleFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.TitleFontSize1,'TooltipString','Set the title font size.','Position',[95 130 50 23],'Callback',@CallbackUI1,'Tag','71');
        a.AxesLabelFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.AxesLabelFontSize1,'TooltipString','Set the axes label font size.','Position',[95 100 50 23],'Callback',@CallbackUI1,'Tag','72');
        a.AxesTickFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.AxesTickFontSize1,'TooltipString','Set the axes tick label font size.','Position',[95 70 50 23],'Callback',@CallbackUI1,'Tag','73');
        a.ModeLabelFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.ModeLabelFontSize1,'TooltipString','Set the mode label font size.','Position',[95 40 50 23],'Callback',@CallbackUI1,'Tag','74');
        a.LegendFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.LegendFontSize1,'TooltipString','Set the legend font size.','Position',[95 10 50 23],'Callback',@CallbackUI1,'Tag','75');

        a.c3UI1 = uicontrol('Parent',Tab1,'Style','pushbutton','String','Default','TooltipString','Reset the plot layout settings to the default.','Position',[1052 45 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','76');    

        a.a1UI1 = axes('Parent',Tab1,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading20.png')
    end
    if  1 % Tab2_anisotropic
        a.p1UI2 = uipanel('Parent',Tab2,'Title','Specimen','Units','pixels','Position',[10 525 225 230],'FontSize',10);
        uicontrol('Parent',a.p1UI2,'Style','pushbutton','String','Edit','TooltipString','Press to define your specimen.','Position',[58 165 127 40],'FontSize',10,'Callback',@SpecimenSettingsUI2_Callback);
        uicontrol('Parent',a.p1UI2,'Style','text','String','Fluids','Position',[10 135 30 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','String','Material','Position',[10 105 39 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','String','Layup','Position',[10 75 32 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','String','Effective','Position',[10 45 45 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','String','Layers , d (mm)','Position',[10 15 78 13]);
        a.UpperFluidDisplayUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','vacuum','TooltipString','Fluid above the laminate.','Position',[58 130 70 23],'BackgroundColor','white');
        a.LowerFluidDisplayUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','vacuum','TooltipString','Fluid below the laminate.','Position',[140 130 70 23],'BackgroundColor','white');
        a.MaterialNameUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String',a.Material2{1}.Name,'Position',[58 100 152 23],'BackgroundColor','white');
        a.LayupUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','[0]','Position',[58 70 152 23],'BackgroundColor','white');
        a.EffectiveLayupUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','[0]','TooltipString','= Layup - Propagation angle','Position',[58 40 152 23],'BackgroundColor','white');
        a.LayerCountUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','1','TooltipString','Number of layers in the laminate.','Position',[100 10 50 23],'BackgroundColor','white');
        a.ThicknessCountUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String',a.PlateThickness,'TooltipString','Overall thickness of the laminate.','Position',[160 10 50 23],'BackgroundColor','white');

        a.p2UI2 = uipanel('Parent',Tab2,'Title','Computational settings','Units','pixels','Position',[10 375 225 145],'FontSize',10);
        uicontrol('Parent',a.p2UI2,'Style','text','String',['Propagation angle (',char(176),')'],'Position',[10 105 103 13]);
        uicontrol('Parent',a.p2UI2,'Style','text','String','Phase velocity limit (m/ms)','Position',[10 75 128 13]);
        uicontrol('Parent',a.p2UI2,'Style','text','String','Frequency limit (kHz)','Position',[10 45 103 13]);
        uicontrol('Parent',a.p2UI2,'Style','text','String','Frequency step (kHz)','Position',[10 15 107 13]);
        uicontrol('Parent',a.p2UI2,'Style','edit','String',a.PropagationAngle,'TooltipString',['Enter the angle with respect to the fiber orientations',newline,'along which the guided waves shall propagate.'],'Position',[160 100 50 23],'Callback',@CallbackUI2,'Tag','1');
        a.PhaseVelocityLimitUI2 = uicontrol('Parent',a.p2UI2,'Style','edit','String',a.PhaseVelocityLimit2/1e3,'TooltipString',['Enter the phase velocity (Y-axis in the dispersion diagram)',newline,'up to which the dispersion curves shall be traced.'],'Position',[160 70 50 23],'Callback',@CallbackUI2,'Tag','2');
        a.FrequencyLimitUI2 = uicontrol('Parent',a.p2UI2,'Style','edit','String',a.FrequencyLimit2,'TooltipString',['Enter the frequency (X-axis in the dispersion diagram)',newline,'up to which the dispersion curves shall be traced.'],'Position',[160 40 50 23],'Callback',@CallbackUI2,'Tag','3');
        a.FrequencyResolutionUI2 = uicontrol('Parent',a.p2UI2,'Style','edit','String',a.FrequencyResolution2,'TooltipString','Enter the frequency step.','Position',[160 10 50 23],'Callback',@CallbackUI2,'Tag','4');

        a.p3UI2 = uipanel('Parent',Tab2,'Title','Mode selection','Units','pixels','Position',[10 165 225 205],'FontSize',10);
        uicontrol('Parent',a.p3UI2,'Style','text','String','Higher order modes','Position',[10 165 97 13]);
        uicontrol('Parent',a.p3UI2,'Style','text','String','Symmetric modes','Position',[10 135 87 13]);
        uicontrol('Parent',a.p3UI2,'Style','text','String','Antisymmetric modes','Position',[10 105 105 13]);
        uicontrol('Parent',a.p3UI2,'Style','text','String','Lamb modes','Position',[10 75 63 13]);
        uicontrol('Parent',a.p3UI2,'Style','text','String','Shear horizontal modes','Position',[10 45 117 13]);
        uicontrol('Parent',a.p3UI2,'Style','text','String','Scholte modes','Position',[10 15 73 13]);
        uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.HigherOrderModes2,'TooltipString','Check this in order to calculate the higher order modes in addition to the fundamental ones.','Position',[160 160 50 23],'Callback',@CallbackUI2,'Tag','7');
        a.SymmetricModesUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.SymmetricModes2,'TooltipString',['If the layup is symmtric with respect to the middle plane, one can distinguish',newline,'symmetric and antisymmetric modes. Check this in order to calculate the',newline,'symmetric modes. These modes have a symmetric displacement pattern with',newline,'respect to the middle plane of the plate.'],'Position',[160 130 50 23],'Callback',@CallbackUI2,'Tag','9');
        a.AntisymmetricModesUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.AntisymmetricModes2,'TooltipString',['Check this in order to calculate the antisymmetric modes.',newline,'These modes have an antisymmetric displacement',newline,'pattern with respect to the middle plane of the plate.'],'Position',[160 100 50 23],'Callback',@CallbackUI2,'Tag','10');
        uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.LambModes2,'TooltipString',['Check this in order to calculate the pure Lamb wave modes. This option is available only',newline,'in the decoupled case. Lamb and shear horizontal waves are decoupled if the specimen',newline,'consists solely of 0',char(176),' and 90',char(176),' layers, and if the wave propagation is along either direction.',newline,'In the decoupled case, Lamb waves are termed ''pure'' modes. These Lamb waves show',newline,'displacement only in the sagittal plane spanned by the propagation direction x1 and by',newline,'the out-of-plane direction x3. In the coupled case, the modes are still called Lamb waves,',newline,'but they have a third displacement component perpendicular to the propagation direction',newline,'(shear horizontal, x2), and no distinction between Lamb and shear horizontal modes is',newline,'possible. Lamb waves are indicated by solid lines in the dispersion diagram.'],'Position',[160 70 50 23],'Callback',@CallbackUI2,'Tag','13');
        a.ShearHorizontalModesUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.ShearHorizontalModes2,'TooltipString',['The so-called ''pure'' shear horizontal modes can only exist in the decoupled case. Otherwise,',newline,'they cannot be distinguished from Lamb waves. The pure shear horizontal modes, however,',newline,'show displacement only perpendicular (x2) to the propagation direction. Shear horizontal',newline,'waves are indicated by dashed lines in the dispersion diagram.'],'Position',[160 40 50 23],'Callback',@CallbackUI2,'Tag','14');
        a.ScholteModesUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.ScholteModes2,'TooltipString',['Check this in order to calculate the Scholte wave modes. Scholte waves are surface waves propagating',newline,'at an interface between a fluid and a solid. They are solved by the same characteristic function as the',newline,'Lamb waves. Scholte waves are indicated by dashed-dotted lines in the dispersion diagram.'],'Position',[160 10 50 23],'Enable','off','Callback',@CallbackUI2,'Tag','5');        
        
        %------------------------------------------------------------------
        a.p7UI2 = uipanel('Parent',Tab2,'Title','Manually detect higher order modes','Units','pixels','Position',[245 690 223 65],'FontSize',10);
        uicontrol('Parent',a.p7UI2,'Style','text','String','Step (kHz)','Position',[10 20 53 13]);
        a.SteptUI2 = uicontrol('Parent',a.p7UI2,'Style','edit','String',a.Step2,'TooltipString',['Enter the step size for the frequency sweep. In general, a finer step increases the chance to',newline,'find mode cut-off frequencies, although exceptions from that rule might occur occasionally.',newline,'It is essential that all higher order modes are detected. Therefore, try smaller and larger step',newline,'sizes if you are not sure that you have found all modes.'],'Position',[75 15 50 23],'Callback',@CallbackUI2,'Tag','15');
        uicontrol('Parent',a.p7UI2,'Style','pushbutton','String','Detect','TooltipString',['Before you can trace the higher order modes, their cut-off frequencies',newline,'at the phase velocity limit must be detected. This is done automatically',newline,'upon pressing ''calculate'' if it has not already been done manually.',newline,'Sometimes the automatic search does not find all modes, then change',newline,'''Step'' and press ''Detect'' to find the missing modes.'],'Position',[140 10 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','16');

        a.OutputWindow1aUI2 = uicontrol('Parent',Tab2,'Style','text','String','','Position',[245 50 105 630],'BackgroundColor','white');
        a.OutputWindow1bUI2 = uicontrol('Parent',Tab2,'Style','text','String','','Position',[362 50 105 630],'BackgroundColor','white');
        a.OutputWindow2aUI2 = uicontrol('Parent',Tab2,'Style','text','String','','Position',[245 10 105 28],'BackgroundColor','white');
        a.OutputWindow2bUI2 = uicontrol('Parent',Tab2,'Style','text','String','','Position',[362 10 105 28],'BackgroundColor','white');

        %------------------------------------------------------------------
        a.c1UI2 = uicontrol('Parent',Tab2,'Style','pushbutton','String','Calculate','TooltipString','Start the dispersion curve tracing.','Position',[478 702 127 40],'FontSize',10,'Callback',@CallbackUI2,'Tag','17');
        a.c2UI2 = uicontrol('Parent',Tab2,'Style','pushbutton','String','Stop calculation','TooltipString','Stop the dispersion curve tracing.','Position',[615 702 127 40],'FontSize',10,'Callback',@CallbackUI2,'Tag','18');

        a.p8UI2 = uipanel('Parent',Tab2,'Title','Dispersion diagrams','Units','pixels','Position',[478 465 263 225],'FontSize',10);
        uicontrol('Parent',a.p8UI2,'Style','text','String','Quantity','Position',[10 185 42 13]);
        a.Option1TextUI2 = uicontrol('Parent',a.p8UI2,'Style','text','String','Bulk velocities','Position',[10 155 70 13]);
        a.XAxisModeTextUI2 = uicontrol('Parent',a.p8UI2,'Style','text','String','X-axis mode','Position',[10 125 62 13]);
        a.XAxisTextUI2 = uicontrol('Parent',a.p8UI2,'Style','text','String','X-axis (kHz)','Position',[10 95 62 13]);
        a.YAxisTextUI2 = uicontrol('Parent',a.p8UI2,'Style','text','String','Y-axis (m/ms)','Position',[10 65 70 13]);
        uicontrol('Parent',a.p8UI2,'Style','popupmenu','String',{'Phase velocity (m/ms)','Energy velocity (m/ms)',['Propagation time (',char(181),'s)'],['Coincidence angle (',char(176),')'],'Wavelength (mm)','Wavenumber (rad/mm)','Attenuation (Np/m)'},'TooltipString','Select which quantity to plot in the dispersion diagram.','Position',[117 180 130 23],'Callback',@CallbackUI2,'Tag','19');
        a.Option1UI2 = uicontrol('Parent',a.p8UI2,'Style','checkbox','Value',a.BulkVelocities2,'TooltipString','Check this to show the bulk wave velocities.','Position',[117 150 50 23],'Callback',@CallbackUI2,'Tag','20');
        a.XAxisModeUI2 = uicontrol('Parent',a.p8UI2,'Style','popupmenu','String',{'Frequency (kHz)','Frequency (MHz)',['f',char(8901),'d (MHz',char(8901),'mm)']},'TooltipString','Select the frequency''s dimension on the X-axis.','Position',[117 120 130 23],'Callback',@CallbackUI2,'Tag','21');
        a.XAxisUI2 = uicontrol('Parent',a.p8UI2,'Style','edit','String',['[0 ',num2str(a.FrequencyLimit2),']'],'TooltipString','Enter which frequency range shall be plotted.','Position',[117 90 75 23],'Callback',@CallbackUI2,'Tag','22');
        a.YAxisUI2 = uicontrol('Parent',a.p8UI2,'Style','edit','String',['[0 ',num2str(a.PhaseVelocityLimit2/1e3),']'],'TooltipString','Enter which phase velocity range shall be plotted.','Position',[117 60 75 23],'Callback',@CallbackUI2,'Tag','23');
        a.Plot1UI2 = uicontrol('Parent',a.p8UI2,'Style','pushbutton','String','Plot','TooltipString',['Plot the dispersion diagram. If you have checked ''Export plots''',newline,'in the export settings, the plot will be exported automatically.'],'Position',[117 15 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','24');

        a.p9UI2 = uipanel('Parent',Tab2,'Title','Through-thickness profiles','Units','pixels','Position',[751 465 224 290],'FontSize',10);
        uicontrol('Parent',a.p9UI2,'Style','text','String','Quantity','Position',[10 245 42 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','String','Mode','Position',[10 215 28 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','String','Frequency (kHz)','Position',[10 185 83 13]);
        a.Option2TextUI2 = uicontrol('Parent',a.p9UI2,'Style','text','String','Samples per layer','Position',[10 155 89 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','String','Half-spaces','Position',[10 125 61 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','String','Phase','Position',[10 95 32 13]);
        uicontrol('Parent',a.p9UI2,'Style','popupmenu','String',{'Displacement','Stress','Strain','Energy density','Power flow density'},'TooltipString','Select which through-thickness quantity to plot.','Position',[119 240 90 23],'Callback',@CallbackUI2,'Tag','25');
        a.Mode1UI2 = uicontrol('Parent',a.p9UI2,'Style','popupmenu','String',{''},'TooltipString','Select the mode you want to analyze.','Position',[119 210 90 23],'Callback',@CallbackUI2,'Tag','26');
        a.Frequency1UI2 = uicontrol('Parent',a.p9UI2,'Style','edit','String',a.FrequencyLimit2,'TooltipString','Enter the frequency at which to analyze the selected mode.','Position',[119 180 50 23],'Callback',@CallbackUI2,'Tag','27');
        a.Samples1UI2 = uicontrol('Parent',a.p9UI2,'Style','edit','String',a.Samples12,'TooltipString',['Enter the number of sample points per layer over the plate''s thickness',newline,'(x3) at which the selected quantities are calculated.'],'Position',[119 150 50 23],'Callback',@CallbackUI2,'Tag','28');
        a.Halfspaces1UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','Value',a.Halfspaces12,'TooltipString','Check this to show the quantities in the upper and lower fluid.','Position',[89 120 20 23],'Enable','off','Callback',@CallbackUI2,'Tag','12');         
        a.HalfspacesNumber1UI2 = uicontrol('Parent',a.p9UI2,'Style','edit','String',a.HalfspacesNumber12,'TooltipString','Set the height of the half-spaces in plate thicknesses.','Position',[119 120 50 23],'Enable','off','Callback',@CallbackUI2,'Tag','11');
        a.PhaseUI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','Value',a.Phase2,'TooltipString','Check this to plot also the phase of the field components.','Position',[119 90 20 23],'Callback',@CallbackUI2,'Tag','81');
        a.Plot2UI2 = uicontrol('Parent',a.p9UI2,'Style','pushbutton','String','Plot','TooltipString',['Plot the profile. If you have checked ''Export plots'' in the',newline,'export settings, the plot will be exported automatically.'],'Position',[119 15 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','29');
        a.x11UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','11','Value',a.Plot2(1),'TooltipString','','Position',[10 50 35 23],'Enable','off','Callback',@CallbackUI2,'Tag','30');
        a.x22UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','22','Value',a.Plot2(2),'TooltipString','','Position',[45 30 35 23],'Enable','off','Callback',@CallbackUI2,'Tag','31');
        a.x33UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','3','Value',a.Plot2(3),'TooltipString','Displacement u3','Position',[80 10 40 23],'Callback',@CallbackUI2,'Tag','32');
        a.x13UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','1','Value',a.Plot2(4),'TooltipString','Displacement u1','Position',[80 50 40 23],'Callback',@CallbackUI2,'Tag','33');
        a.x23UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','2','Value',a.Plot2(5),'TooltipString','Displacement u2','Position',[80 30 40 23],'Callback',@CallbackUI2,'Tag','34');
        a.x12UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','12','Value',a.Plot2(6),'TooltipString','','Position',[45 50 35 23],'Enable','off','Callback',@CallbackUI2,'Tag','35');        
        
        a.p10UI2 = uipanel('Parent',Tab2,'Title','Mode shape','Units','pixels','Position',[478 195 497 265],'FontSize',10);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Mode','Position',[10 225 28 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Frequency (kHz)','Position',[10 195 83 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Wavelengths','Position',[10 165 65 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Samples x1','Position',[10 135 58 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Samples per layer','Position',[10 105 89 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Scale','Position',[10 75 29 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Grid line','Position',[10 45 41 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Undistorted','Position',[10 15 57 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','String','Half-spaces','Position',[195 45 61 13]);
        a.Mode2UI2 = uicontrol('Parent',a.p10UI2,'Style','popupmenu','String',{''},'TooltipString','Select the mode you want to analyze.','Position',[117 220 90 23],'Callback',@CallbackUI2,'Tag','36');        
        a.Frequency2UI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Frequency22,'TooltipString','Enter the frequency at which to analyze the selected mode.','Position',[117 190 50 23],'Callback',@CallbackUI2,'Tag','37');
        uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Length2,'TooltipString','Specify how many wavelengths you want to display.','Position',[117 160 50 23],'Callback',@CallbackUI2,'Tag','38');
        uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Samples22,'TooltipString',['Enter the number of sample points along the propagation',newline,'direction (x1) at which the displacement is calculated.'],'Position',[117 130 50 23],'Callback',@CallbackUI2,'Tag','39');
        a.Samples3UI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Samples32,'TooltipString',['Enter the number of sample points per layer over the plate''s',newline,'thickness (x3) at which the displacement is calculated.'],'Position',[117 100 50 23],'Callback',@CallbackUI2,'Tag','40');
        uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Scale2,'TooltipString','Change the scaling of the displacement.','Position',[117 70 50 23],'Callback',@CallbackUI2,'Tag','41');
        a.GridLineUI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.GridLine2,'TooltipString','Draw a grid line at every ith sample point.','Position',[117 40 50 23],'Callback',@CallbackUI2,'Tag','42');
        uicontrol('Parent',a.p10UI2,'Style','checkbox','Value',a.Undistorted2,'TooltipString','Check this to draw the undistorted grid.','Position',[117 10 50 23],'Callback',@CallbackUI2,'Tag','43');
        a.HalfspacesNumber2UI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.HalfspacesNumber22,'TooltipString','Set the height of the half-spaces in plate thicknesses.','Position',[272 40 50 23],'Enable','off','Callback',@CallbackUI2,'Tag','6');
        a.Halfspaces2UI2 = uicontrol('Parent',a.p10UI2,'Style','checkbox','Value',a.Halfspaces22,'TooltipString','Check this to show the quantities in the upper and lower fluid.','Position',[272 10 50 23],'Enable','off','Callback',@CallbackUI2,'Tag','8');
        a.Plot3UI2 = uicontrol('Parent',a.p10UI2,'Style','pushbutton','String','Plot','TooltipString',['Plot the mode shape. If you have checked ''Animate'',',newline,'the mode shape will be animated. If you have',newline,'checked ''Export plots'' in the export settings, the',newline,'plot/movie will be exported automatically.'],'Position',[392 15 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','44');
       
        a.p11UI2 = uipanel('Parent',Tab2,'Title','Animation settings','Units','pixels','Position',[751 265 212 175],'FontSize',9);
        uicontrol('Parent',a.p11UI2,'Style','text','String','Cycles','Position',[10 135 35 13]);
        uicontrol('Parent',a.p11UI2,'Style','text','String','Cycle duration (s)','Position',[10 105 88 13]);
        uicontrol('Parent',a.p11UI2,'Style','text','String','Frame rate (Hz)','Position',[10 75 78 13]);
        uicontrol('Parent',a.p11UI2,'Style','text','String','Movie quality (0-100)','Position',[10 45 103 13]);
        uicontrol('Parent',a.p11UI2,'Style','text','String','Animate','Position',[10 15 41 13]);
        uicontrol('Parent',a.p11UI2,'Style','edit','String',a.Cycles2,'TooltipString',['Enter how many cycles shall be calculated. That cylce',newline,'will be repeated until you close the plot figure.'],'Position',[119 130 50 23],'Callback',@CallbackUI2,'Tag','45');
        uicontrol('Parent',a.p11UI2,'Style','edit','String',a.CycleDuration2,'TooltipString','Enter how long a cycle shall take.','Position',[119 100 50 23],'Callback',@CallbackUI2,'Tag','46');
        uicontrol('Parent',a.p11UI2,'Style','edit','String',a.FrameRate2,'TooltipString','Define the frame rate of the movie.','Position',[119 70 50 23],'Callback',@CallbackUI2,'Tag','47');
        uicontrol('Parent',a.p11UI2,'Style','edit','String',a.MovieQuality2,'TooltipString','Define the quality of the exported movie.','Position',[119 40 50 23],'Callback',@CallbackUI2,'Tag','48');
        uicontrol('Parent',a.p11UI2,'Style','checkbox','Value',a.Animate2,'TooltipString',['Check this in order to show the animated mode',newline,'shape upon pressing the plot button below.'],'Position',[119 10 70 23],'Callback',@CallbackUI2,'Tag','49');

        a.p12UI2 = uipanel('Parent',Tab2,'Title','Export settings','Units','pixels','Position',[478 10 497 180],'FontSize',10);
        uicontrol('Parent',a.p12UI2,'Style','text','String','Export plots','Position',[10 140 59 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','String','Crop plots','Position',[10 120 51 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','String','PDF','Position',[124 140 21 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','String','PNG','Position',[124 120 23 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','String','PNG resolution (dpi)','Position',[10 93 98 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','String','Dispersion curves','Position',[240 140 90 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','String','Through-thickness','Position',[240 120 92 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','String','Directory','Position',[10 20 46 13]);
        uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.ExportPlots2,'TooltipString',['Check this in order to export plots/movies automatically upon',newline,'pressing the respective plot button. You can also export plots',newline,'manually by using the ''File'' menu inside the plot figure.'],'Position',[80 135 20 23],'Callback',@CallbackUI2,'Tag','50');
        uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.CropPlots2,'TooltipString','Plots will be cropped tightly with minimal white space.','Position',[80 115 20 23],'Callback',@CallbackUI2,'Tag','79');
        uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.PDF2,'TooltipString','Check this to export a plot as pdf.','Position',[160 135 20 23],'Callback',@CallbackUI2,'Tag','51');
        uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.PNG2,'TooltipString','Check this to export a plot as png.','Position',[160 115 20 23],'Callback',@CallbackUI2,'Tag','52');
        uicontrol('Parent',a.p12UI2,'Style','edit','String',a.PNGresolution2,'TooltipString','Enter the resolution of the png image.','Position',[125 88 50 23],'Callback',@CallbackUI2,'Tag','53');
        
        a.XAxisMode2UI2 = uicontrol('Parent',a.p12UI2,'Style','popupmenu','String',{'Frequency (kHz)','Frequency (MHz)',['f',char(8901),'d (MHz',char(8901),'mm)']},'TooltipString','Select the frequency''s dimension to be exported.','Position',[370 130 110 23],'Callback',@CallbackUI2,'Tag','54');
        a.ArrangeUI2 = uicontrol('Parent',a.p12UI2,'Style','popupmenu','String',{'Horizontal arrangement','Vertical arrangement'},'TooltipString','Choose how to arrange the dispersion curve data in the Excel/txt.','Position',[370 95 110 23],'Callback',@CallbackUI2,'Tag','55');        
        uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.DispersionCurves2,'TooltipString','Check this to export the dispersion curves.','Position',[340 135 20 23],'Callback',@CallbackUI2,'Tag','56');
        uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.ThroughThickness2,'TooltipString','Check this to export the through-thickness profiles selected above.','Position',[340 115 20 23],'Callback',@CallbackUI2,'Tag','57');        
        uicontrol('Parent',a.p12UI2,'Style','pushbutton','String','*.mat','TooltipString','Export the data as Matlab''s mat-file.','Position',[220 80 40 33],'Callback',@CallbackUI2,'Tag','80');
        uicontrol('Parent',a.p12UI2,'Style','pushbutton','String','*.xlsx','TooltipString','Export the data as Excel sheet.','Position',[270 80 40 33],'Callback',@CallbackUI2,'Tag','58');
        uicontrol('Parent',a.p12UI2,'Style','pushbutton','String','*.txt','TooltipString','Export the data as txt-file.','Position',[320 80 40 33],'Callback',@CallbackUI2,'Tag','59');      
        uicontrol('Parent',a.p12UI2,'Style','edit','String',a.FileName2,'TooltipString',['Specify the name of the plots/movies to be exported.',newline,'The dispersion curve raw data are named automatically.'],'Position',[70 45 412 23],'Callback',@CallbackUI2,'Tag','60');
        uicontrol('Parent',a.p12UI2,'Style','edit','String',a.Directory,'TooltipString',['Specify the directory to which plots, movies, and',newline,'the dispersion curve raw data shall be exported.'],'Position',[70 15 412 23],'Callback',@CallbackUI2,'Tag','61');

        %------------------------------------------------------------------
        a.p13UI2 = uipanel('Parent',Tab2,'Title','Plot layout settings','Units','pixels','Position',[985 90 200 600],'FontSize',10);
        uicontrol('Parent',a.p13UI2,'Style','text','String','Title','Position',[10 560 21 13]);
        uicontrol('Parent',a.p13UI2,'Style','text','String','Mode labels','Position',[10 530 59 13]);
        uicontrol('Parent',a.p13UI2,'Style','text','String','Legend location','Position',[10 500 78 13]);
        uicontrol('Parent',a.p13UI2,'Style','text','String','Box line width','Position',[10 470 70 13]);
        uicontrol('Parent',a.p13UI2,'Style','text','String','Curve line width','Position',[10 440 80 13]);
        a.TitleUI2 = uicontrol('Parent',a.p13UI2,'Style','popupmenu','String',{'with layup','without layup','no title'},'TooltipString',['Choose to plot the title with the layup ',newline,'included, without layup, or no title at all.'],'Position',[105 555 80 23],'Callback',@CallbackUI2,'Tag','62');
        a.ModeLabelsUI2 = uicontrol('Parent',a.p13UI2,'Style','checkbox','Value',a.ModeLabels2,'TooltipString',['Check this in order to show the mode',newline,'labels of the fundamental modes.'],'Position',[105 525 20 23],'Callback',@CallbackUI2,'Tag','63');
        a.LegendLocationUI2 = uicontrol('Parent',a.p13UI2,'Style','popupmenu','String',{'outside','inside'},'TooltipString','Determine the legend location in the through-thickness plots.','Position',[105 495 80 23],'Callback',@CallbackUI2,'Tag','64');
        a.BoxLineWidthUI2 = uicontrol('Parent',a.p13UI2,'Style','edit','String',a.BoxLineWidth2,'TooltipString','Enter the box line width.','Position',[105 465 50 23],'Callback',@CallbackUI2,'Tag','65');
        a.LineWidthUI2 = uicontrol('Parent',a.p13UI2,'Style','edit','String',a.LineWidth2,'TooltipString','Enter the curve line width.','Position',[105 435 50 23],'Callback',@CallbackUI2,'Tag','66');

        a.p14UI2 = uipanel('Parent',Tab2,'Title','Dispersion curve colors [R G B]','Units','pixels','Position',[995 400 180 115],'FontSize',9);
        uicontrol('Parent',a.p14UI2,'Style','text','String','S','Position',[10 75 7 13]);
        uicontrol('Parent',a.p14UI2,'Style','text','String','A','Position',[10 45 8 13]);
        uicontrol('Parent',a.p14UI2,'Style','text','String','B','Position',[10 15 7 13]);
        a.SColorUI2 = uicontrol('Parent',a.p14UI2,'Style','edit','String','[1 0 0]','TooltipString',['Specify the color of symmetric modes. You can compose any color',newline,'by entering the corresponding RGB values. Set numbers from 0 to 1.'],'Position',[95 70 50 23],'Callback',@CallbackUI2,'Tag','67');
        a.AColorUI2 = uicontrol('Parent',a.p14UI2,'Style','edit','String','[0 0 1]','TooltipString','Specify the color of antisymmetric modes.','Position',[95 40 50 23],'Callback',@CallbackUI2,'Tag','68');
        a.BColorUI2 = uicontrol('Parent',a.p14UI2,'Style','edit','String','[.5 0 1]','TooltipString','Specify the color of nonsymmetric modes.','Position',[95 10 50 23],'Callback',@CallbackUI2,'Tag','69');

        a.p15UI2 = uipanel('Parent',Tab2,'Title','Mode labels x-position (0-1)','Units','pixels','Position',[995 280 180 115],'FontSize',9);
        uicontrol('Parent',a.p15UI2,'Style','text','String','S0/B1','Position',[10 75 31 13]);
        uicontrol('Parent',a.p15UI2,'Style','text','String','SSH0/BSH0','Position',[10 45 59 13]);
        uicontrol('Parent',a.p15UI2,'Style','text','String','A0/B0','Position',[10 15 32 13]);
        a.S0B1LabelUI2 = uicontrol('Parent',a.p15UI2,'Style','edit','String',a.S0B1Label,'TooltipString',['Change the relative X-position of the fundamental symmetric (S0)',newline,'or nonsymmetric (B1) Lamb wave mode label. 0 corresponds to the',newline,'left, 1 to the right X-axis limit.'],'Position',[95 70 50 23],'Callback',@CallbackUI2,'Tag','70');
        a.SH0LabelUI2 = uicontrol('Parent',a.p15UI2,'Style','edit','String',a.SH0Label2,'TooltipString',['Change the relative X-position of the fundamental shear horizontal wave',newline,'mode label.'],'Position',[95 40 50 23],'Callback',@CallbackUI2,'Tag','71');
        a.A0B0LabelUI2 = uicontrol('Parent',a.p15UI2,'Style','edit','String',a.A0B0Label,'TooltipString',['Change the relative X-position of the fundamental antisymmetric (A0)',newline,'or nonsymmetric (B0) Lamb wave mode label.'],'Position',[95 10 50 23],'Callback',@CallbackUI2,'Tag','72');

        a.p16UI2 = uipanel('Parent',Tab2,'Title','Font size','Units','pixels','Position',[995 100 180 175],'FontSize',9);
        uicontrol('Parent',a.p16UI2,'Style','text','String','Title','Position',[10 135 21 13]);
        uicontrol('Parent',a.p16UI2,'Style','text','String','Axes labels','Position',[10 105 59 13]);
        uicontrol('Parent',a.p16UI2,'Style','text','String','Axes ticks','Position',[10 75 53 13]);
        uicontrol('Parent',a.p16UI2,'Style','text','String','Mode labels','Position',[10 45 59 13]);
        uicontrol('Parent',a.p16UI2,'Style','text','String','Legend','Position',[10 15 38 13]);
        a.TitleFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.TitleFontSize2,'TooltipString','Set the title font size.','Position',[95 130 50 23],'Callback',@CallbackUI2,'Tag','73');
        a.AxesLabelFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.AxesLabelFontSize2,'TooltipString','Set the axes label font size.','Position',[95 100 50 23],'Callback',@CallbackUI2,'Tag','74');
        a.AxesTickFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.AxesTickFontSize2,'TooltipString','Set the axes tick label font size.','Position',[95 70 50 23],'Callback',@CallbackUI2,'Tag','75');
        a.ModeLabelFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.ModeLabelFontSize2,'TooltipString','Set the mode label font size.','Position',[95 40 50 23],'Callback',@CallbackUI2,'Tag','76');
        a.LegendFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.LegendFontSize2,'TooltipString','Set the legend font size.','Position',[95 10 50 23],'Callback',@CallbackUI2,'Tag','77');

        a.c3UI2 = uicontrol('Parent',Tab2,'Style','pushbutton','String','Default','TooltipString','Reset the plot layout settings to the default.','Position',[1052 45 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','78');    

        a.a1UI2 = axes('Parent',Tab2,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading20.png') 
    end
    if  1 % Tab3_signal simulator
        a.DataUI3 = uicontrol('Parent',Tab3,'Style','text','String','','TooltipString','The material to which the dispersion curves belonge to.','Position',[10 740 279 13],'BackgroundColor','white');
        
        a.p1UI3 = uipanel('Parent',Tab3,'Title','Computational settings','Units','pixels','Position',[10 420 279 315],'FontSize',10);         
        uicontrol('Parent',a.p1UI3,'Style','text','String','Frequency (kHz)','Position',[10 275 83 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','String','Cycles, Samples/cycle','Position',[10 245 111 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','String','Window','Position',[10 215 42 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','String','Distance (mm)','Position',[10 185 71 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','String',['n',char(8901),'Distance/ce'],'Position',[10 155 70 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','String',['Spectral threshold (',char(37),')'],'Position',[10 125 111 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','String','Displacement component','Position',[10 95 122 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','String',['Gate (',char(181),'s)'],'Position',[10 65 49 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','String','Multi-mode','Position',[10 25 53 13]);
        a.FrequencyUI3 = uicontrol('Parent',a.p1UI3,'Style','edit','String','','TooltipString','Enter the carrier wave frequency in the wave packet.','Position',[140 270 50 23],'Callback',@CallbackUI3,'Tag','2');
        uicontrol('Parent',a.p1UI3,'Style','edit','String',a.Cycles3,'TooltipString','Enter the number of cycles of the signal within the wave packet.','Position',[140 240 50 23],'Callback',@CallbackUI3,'Tag','3');
        uicontrol('Parent',a.p1UI3,'Style','edit','String',a.SamplesPerCycle3,'TooltipString',['Enter the samples per cycle of the excitation',newline,'signal. This determines the sample rate.'],'Position',[210 240 50 23],'Callback',@CallbackUI3,'Tag','4');
        uicontrol('Parent',a.p1UI3,'Style','popupmenu','String',{'Gauss','Hann','Hamming','Triangular'},'TooltipString',['Select the window function to be multiplied with the',newline,'carrier wave. This sets the shape of the wave packet.'],'Position',[140 210 120 23],'Callback',@CallbackUI3,'Tag','5');
        uicontrol('Parent',a.p1UI3,'Style','edit','String',a.Distance3,'TooltipString','Enter the propagation distance of the wave packet.','Position',[140 180 50 23],'Callback',@CallbackUI3,'Tag','6');
        uicontrol('Parent',a.p1UI3,'Style','edit','String',a.TimeLimitFactor3,'TooltipString',['Enter the temporal limit of the simulation in multiples ''n'' of the propagation',newline,'distance divided through the energy velocity of the mode in question. In',newline,'case of multiple modes, the lowest energy velocity among all contributing',newline,'modes is used. Due to dispersion, the wave packets spread with the',newline,'propagation time so that the value should always be greater than one.'],'Position',[140 150 50 23],'Callback',@CallbackUI3,'Tag','7');
        uicontrol('Parent',a.p1UI3,'Style','edit','String',a.SpectrumThreshold3,'TooltipString',['Enter the threshold of the spectral amplitudes, which shall be used for',newline,'the construction of the propagated wave packet. The value is in percent',newline,'of the maximal spectral amplitude of the excitation signal at the center',newline,'frequency defined above. A smaller value widens the frequency range',newline,'to be taken into acount.'],'Position',[140 120 50 23],'Callback',@CallbackUI3,'Tag','8');
        uicontrol('Parent',a.p1UI3,'Style','popupmenu','String',{'Out-of-plane (u3)','In-plane (u1/2)'},'TooltipString','Select which displacement component you want to calculate.','Position',[140 90 120 23],'Callback',@CallbackUI3,'Tag','113');
        a.GateUI3 = uicontrol('Parent',a.p1UI3,'Style','edit','String',['[',num2str(a.Gate3(1)),' ',num2str(a.Gate3(2)),']'],'TooltipString',['Enter the gate for the calculation of the propagated spectrum. FFT will',newline,'be performed on the temporal response between the gate''s limits. Set',newline,'the limits left and right of the wave packet you want to consider.'],'Position',[140 60 75 23],'Callback',@CallbackUI3,'Tag','111');
        uicontrol('Parent',a.p1UI3,'Style','checkbox','Value',a.MultiMode3,'TooltipString','Check this to simulate multiple modes at a time.','Position',[100 20 20 23],'Callback',@CallbackUI3,'Tag','9');
        a.CalculateUI3 = uicontrol('Parent',a.p1UI3,'Style','pushbutton','String','Calculate','TooltipString','Calculate multiple modes at a time.','Position',[140 15 65 33],'Enable','off','FontSize',10,'Callback',@CallbackUI3,'Tag','92');
        a.StopUI3 = uicontrol('Parent',a.p1UI3,'Style','pushbutton','String','Stop','TooltipString','Stop the calculation.','Position',[216 15 45 33],'Enable','off','FontSize',10,'Callback',@CallbackUI3,'Tag','109');
        
        a.p2UI3 = uipanel('Parent',Tab3,'Title','Mode selection','Units','pixels','Position',[10 185 279 230],'FontSize',10);
        a.ALamb0aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A0','Position',[10 190 16 13],'Visible','off');
        a.ALamb0bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 188 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','10');
        a.ALamb0cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 190 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','11');
        a.ALamb1aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A1','Position',[10 170 16 13],'Visible','off');
        a.ALamb1bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 168 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','12');
        a.ALamb1cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 170 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','13');
        a.ALamb2aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A2','Position',[10 150 16 13],'Visible','off');
        a.ALamb2bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 148 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','14');
        a.ALamb2cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 150 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','15');
        a.ALamb3aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A3','Position',[10 130 16 13],'Visible','off');
        a.ALamb3bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 128 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','16');
        a.ALamb3cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 130 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','17');
        a.ALamb4aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A4','Position',[10 110 16 13],'Visible','off');
        a.ALamb4bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 108 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','18');
        a.ALamb4cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 110 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','19');
        a.ALamb5aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A5','Position',[10 90 16 13],'Visible','off');
        a.ALamb5bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 88 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','20');
        a.ALamb5cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 90 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','21');
        a.ALamb6aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A6','Position',[10 70 16 13],'Visible','off');
        a.ALamb6bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 68 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','22');
        a.ALamb6cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 70 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','23');
        a.ALamb7aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A7','Position',[10 50 16 13],'Visible','off');
        a.ALamb7bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 48 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','24');
        a.ALamb7cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 50 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','25');
        a.ALamb8aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A8','Position',[10 30 16 13],'Visible','off');
        a.ALamb8bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 28 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','26');
        a.ALamb8cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 30 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','27');
        a.ALamb9aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','A9','Position',[10 10 16 13],'Visible','off');
        a.ALamb9bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 8 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','28');
        a.ALamb9cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[43 10 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','29');        

        a.SLamb0aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S0','Position',[70 190 16 13],'Visible','off');
        a.SLamb0bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 188 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','30');
        a.SLamb0cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 190 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','31');
        a.SLamb1aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S1','Position',[70 170 16 13],'Visible','off');
        a.SLamb1bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 168 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','32');
        a.SLamb1cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 170 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','33');
        a.SLamb2aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S2','Position',[70 150 16 13],'Visible','off');
        a.SLamb2bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 148 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','34');
        a.SLamb2cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 150 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','35');
        a.SLamb3aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S3','Position',[70 130 16 13],'Visible','off');
        a.SLamb3bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 128 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','36');
        a.SLamb3cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 130 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','37');
        a.SLamb4aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S4','Position',[70 110 16 13],'Visible','off');
        a.SLamb4bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 108 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','38');
        a.SLamb4cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 110 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','39');
        a.SLamb5aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S5','Position',[70 90 16 13],'Visible','off');
        a.SLamb5bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 88 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','40');
        a.SLamb5cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 90 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','41');
        a.SLamb6aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S6','Position',[70 70 16 13],'Visible','off');
        a.SLamb6bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 68 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','42');
        a.SLamb6cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 70 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','43');
        a.SLamb7aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S7','Position',[70 50 16 13],'Visible','off');
        a.SLamb7bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 48 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','44');
        a.SLamb7cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 50 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','45');
        a.SLamb8aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S8','Position',[70 30 16 13],'Visible','off');
        a.SLamb8bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 28 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','46');
        a.SLamb8cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 30 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','47');
        a.SLamb9aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','S9','Position',[70 10 16 13],'Visible','off');
        a.SLamb9bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 8 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','48');
        a.SLamb9cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[103 10 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','49');

        a.AShear1aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH1','Position',[129 190 30 13],'Visible','off');
        a.AShear1bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 188 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','50');
        a.AShear1cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 190 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','51'); 
        a.AShear2aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH2','Position',[129 170 30 13],'Visible','off');
        a.AShear2bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 168 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','52');
        a.AShear2cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 170 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','53');
        a.AShear3aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH3','Position',[129 150 30 13],'Visible','off');
        a.AShear3bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 148 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','54');
        a.AShear3cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 150 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','55');
        a.AShear4aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH4','Position',[129 130 30 13],'Visible','off');
        a.AShear4bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 128 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','56');
        a.AShear4cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 130 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','57');
        a.AShear5aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH5','Position',[129 110 30 13],'Visible','off');
        a.AShear5bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 108 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','58');
        a.AShear5cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 110 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','59');
        a.AShear6aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH6','Position',[129 90 30 13],'Visible','off');
        a.AShear6bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 88 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','60');
        a.AShear6cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 90 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','61');
        a.AShear7aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH7','Position',[129 70 30 13],'Visible','off');
        a.AShear7bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 68 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','62');
        a.AShear7cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 70 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','63');
        a.AShear8aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH8','Position',[129 50 30 13],'Visible','off');
        a.AShear8bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 48 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','64');
        a.AShear8cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 50 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','65');
        a.AShear9aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH9','Position',[129 30 30 13],'Visible','off');
        a.AShear9bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 28 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','66');
        a.AShear9cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 30 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','67');
        a.AShear10aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','ASH10','Position',[125 10 36 13],'Visible','off');
        a.AShear10bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 8 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','68');
        a.AShear10cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[175 10 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','69');        
        
        a.SShear0aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH0','Position',[201 190 30 13],'Visible','off');
        a.SShear0bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 188 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','70');
        a.SShear0cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 190 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','71');
        a.SShear1aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH1','Position',[201 170 30 13],'Visible','off');
        a.SShear1bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 168 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','72');
        a.SShear1cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 170 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','73');
        a.SShear2aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH2','Position',[201 150 30 13],'Visible','off');
        a.SShear2bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 148 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','74');
        a.SShear2cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 150 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','75');
        a.SShear3aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH3','Position',[201 130 30 13],'Visible','off');
        a.SShear3bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 128 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','76');
        a.SShear3cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 130 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','77');
        a.SShear4aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH4','Position',[201 110 30 13],'Visible','off');
        a.SShear4bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 108 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','78');
        a.SShear4cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 110 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','79');
        a.SShear5aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH5','Position',[201 90 30 13],'Visible','off');
        a.SShear5bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 88 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','80');
        a.SShear5cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 90 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','81');
        a.SShear6aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH6','Position',[201 70 30 13],'Visible','off');
        a.SShear6bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 68 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','82');
        a.SShear6cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 70 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','83');
        a.SShear7aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH7','Position',[201 50 30 13],'Visible','off');
        a.SShear7bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 48 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','84');
        a.SShear7cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 50 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','85');
        a.SShear8aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH8','Position',[201 30 30 13],'Visible','off');
        a.SShear8bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 28 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','86');
        a.SShear8cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 30 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','87');
        a.SShear9aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','String','SSH9','Position',[201 10 30 13],'Visible','off');
        a.SShear9bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 8 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','88');
        a.SShear9cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','TooltipString','Scale the amplitude.','Position',[247 10 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','89');
        
        a.OutputWindowUI3 = uicontrol('Parent',Tab3,'Style','text','String','','Position',[10 10 195 165],'BackgroundColor','white');
        a.c1UI3 = uicontrol('Parent',Tab3,'Style','text','String',['X-axis (',char(181),'s)'],'Position',[220 165 58 13]);
        a.c2UI3 = uicontrol('Parent',Tab3,'Style','text','String','Y-axis (nm)','Position',[220 115 59 13]);
        a.XAxisUI3 = uicontrol('Parent',Tab3,'Style','edit','String','','TooltipString','Enter which frequency range shall be plotted.','Position',[218 135 65 23],'Callback',@CallbackUI3,'Tag','90');
        a.YAxisUI3 = uicontrol('Parent',Tab3,'Style','edit','String','','TooltipString','Enter which displacement range shall be plotted.','Position',[218 85 65 23],'Callback',@CallbackUI3,'Tag','112'); 
        a.PlotUI3 = uicontrol('Parent',Tab3,'Style','pushbutton','String','Plot','TooltipString',['Plot the temporal and frequency responses. If you have checked ''Export',newline,'plots''in the export settings, the plot will be exported automatically.'],'Position',[218 40 65 33],'FontSize',10,'Callback',@CallbackUI3,'Tag','91');
        
        %------------------------------------------------------------------
        a.p3UI3 = uipanel('Parent',Tab3,'Title','Export settings','Units','pixels','Position',[296 10 446 180],'FontSize',10);
        uicontrol('Parent',a.p3UI3,'Style','text','String','Export plots','Position',[10 140 59 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','String','Crop plots','Position',[10 120 51 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','String','PDF','Position',[124 140 21 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','String','PNG','Position',[124 120 23 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','String','PNG resolution (dpi)','Position',[10 93 98 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','String','Directory','Position',[10 20 46 13]);
        uicontrol('Parent',a.p3UI3,'Style','checkbox','Value',a.ExportPlots3,'TooltipString',['Check this in order to export plots automatically upon',newline,'pressing the plot button. You can also export plots',newline,'manually by using the ''File'' menu inside the plot figure.'],'Position',[80 135 20 23],'Callback',@CallbackUI3,'Tag','93');
        uicontrol('Parent',a.p3UI3,'Style','checkbox','Value',a.CropPlots3,'TooltipString','Plots will be cropped tightly with minimal white space.','Position',[80 115 20 23],'Callback',@CallbackUI3,'Tag','110');
        uicontrol('Parent',a.p3UI3,'Style','checkbox','Value',a.PDF3,'TooltipString','Check this to export a plot as pdf.','Position',[160 135 20 23],'Callback',@CallbackUI3,'Tag','94');
        uicontrol('Parent',a.p3UI3,'Style','checkbox','Value',a.PNG3,'TooltipString','Check this to export a plot as png.','Position',[160 115 20 23],'Callback',@CallbackUI3,'Tag','95');
        uicontrol('Parent',a.p3UI3,'Style','edit','String',a.PNGresolution3,'TooltipString','Enter the resolution of the png image.','Position',[125 88 50 23],'Callback',@CallbackUI3,'Tag','96');
        uicontrol('Parent',a.p3UI3,'Style','pushbutton','String','*.mat','TooltipString','Export the data as Matlab''s mat-file.','Position',[220 80 40 33],'Callback',@CallbackUI3,'Tag','97');
        uicontrol('Parent',a.p3UI3,'Style','pushbutton','String','*.xlsx','TooltipString','Export the data as Excel sheet.','Position',[270 80 40 33],'Callback',@CallbackUI3,'Tag','98');
        uicontrol('Parent',a.p3UI3,'Style','pushbutton','String','*.txt','TooltipString','Export the data as txt-file.','Position',[320 80 40 33],'Callback',@CallbackUI3,'Tag','99');
        uicontrol('Parent',a.p3UI3,'Style','edit','String',a.FileName3,'TooltipString',['Specify the name of the plot to be exported. The',newline,'signal''s raw data are named automatically.'],'Position',[70 45 361 23],'Callback',@CallbackUI3,'Tag','100');
        uicontrol('Parent',a.p3UI3,'Style','edit','String',a.Directory,'TooltipString',['Specify the directory to which the plot and',newline,'the signal''s raw data shall be exported.'],'Position',[70 15 361 23],'Callback',@CallbackUI3,'Tag','101');

        %------------------------------------------------------------------
        a.p4UI3 = uipanel('Parent',Tab3,'Title','Plot layout settings','Units','pixels','Position',[750 10 360 180],'FontSize',10);
        uicontrol('Parent',a.p4UI3,'Style','text','String','Title','Position',[10 115 21 13]);
        uicontrol('Parent',a.p4UI3,'Style','text','String','Box line width','Position',[10 85 70 13]);
        uicontrol('Parent',a.p4UI3,'Style','text','String','Curve line width','Position',[10 55 80 13]);
        a.TitleUI3 = uicontrol('Parent',a.p4UI3,'Style','checkbox','Value',a.Title3,'TooltipString','Check this in order to show the plot title.','Position',[105 110 20 23],'Callback',@CallbackUI3,'Tag','102');
        a.BoxLineWidthUI3 = uicontrol('Parent',a.p4UI3,'Style','edit','String',a.BoxLineWidth3,'TooltipString','Enter the box line width.','Position',[105 80 50 23],'Callback',@CallbackUI3,'Tag','103');
        a.LineWidthUI3 = uicontrol('Parent',a.p4UI3,'Style','edit','String',a.LineWidth3,'TooltipString','Enter the curve line width.','Position',[105 50 50 23],'Callback',@CallbackUI3,'Tag','104');

        %------------------------------------------------------------------
        a.p5UI3 = uipanel('Parent',Tab3,'Title','Font size','Units','pixels','Position',[935 50 160 115],'FontSize',9);
        uicontrol('Parent',a.p5UI3,'Style','text','String','Title','Position',[10 75 21 13]);
        uicontrol('Parent',a.p5UI3,'Style','text','String','Axes labels','Position',[10 45 59 13]);
        uicontrol('Parent',a.p5UI3,'Style','text','String','Axes ticks','Position',[10 15 53 13]);
        a.TitleFontSizeUI3 = uicontrol('Parent',a.p5UI3,'Style','edit','String',a.TitleFontSize3,'TooltipString','Set the title font size.','Position',[95 70 50 23],'Callback',@CallbackUI3,'Tag','105');
        a.AxesLabelFontSizeUI3 = uicontrol('Parent',a.p5UI3,'Style','edit','String',a.AxesLabelFontSize3,'TooltipString','Set the axes label font size.','Position',[95 40 50 23],'Callback',@CallbackUI3,'Tag','106');
        a.AxesTickFontSizeUI3 = uicontrol('Parent',a.p5UI3,'Style','edit','String',a.AxesTickFontSize3,'TooltipString','Set the axes tick label font size.','Position',[95 10 50 23],'Callback',@CallbackUI3,'Tag','107');

        a.c3UI3 = uicontrol('Parent',Tab3,'Style','pushbutton','String','Default','TooltipString','Reset the plot layout settings to the default.','Position',[1120 149 65 33],'FontSize',10,'Callback',@CallbackUI3,'Tag','108');                

        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
    end
    if  1 % Tab4_polar diagrams
        a.p1UI4 = uipanel('Parent',Tab4,'Title','Specimen','Units','pixels','Position',[10 555 225 200],'FontSize',10);
        uicontrol('Parent',a.p1UI4,'Style','pushbutton','String','Edit','TooltipString','Press to define your specimen.','Position',[58 135 127 40],'FontSize',10,'Callback',@SpecimenSettingsUI4_Callback);
        uicontrol('Parent',a.p1UI4,'Style','text','String','Material','Position',[10 105 39 13]);
        uicontrol('Parent',a.p1UI4,'Style','text','String','Layup','Position',[10 75 32 13]);
        uicontrol('Parent',a.p1UI4,'Style','text','String','Layers','Position',[10 45 36 13]);
        uicontrol('Parent',a.p1UI4,'Style','text','String','Thickness (mm)','Position',[10 15 78 13]);
        a.MaterialNameUI4 = uicontrol('Parent',a.p1UI4,'Style','text','String',a.Material_Polar{1}.Name,'Position',[58 100 152 23],'BackgroundColor','white');
        a.LayupUI4 = uicontrol('Parent',a.p1UI4,'Style','text','String','[0]','Position',[58 70 152 23],'BackgroundColor','white');
        a.LayerCountUI4 = uicontrol('Parent',a.p1UI4,'Style','text','String','1','TooltipString','Number of layers in the laminate.','Position',[160 40 50 23],'BackgroundColor','white');
        a.ThicknessCountUI4 = uicontrol('Parent',a.p1UI4,'Style','text','String',a.PlateThickness_Polar,'TooltipString','Overall thickness of the laminate.','Position',[160 10 50 23],'BackgroundColor','white');

        a.p2UI4 = uipanel('Parent',Tab4,'Title','Computational settings','Units','pixels','Position',[10 375 225 175],'FontSize',10);
        uicontrol('Parent',a.p2UI4,'Style','text','String','Frequency limit (kHz)','Position',[10 135 103 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','String','Frequency step (kHz)','Position',[10 105 107 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','String',['Propagation angle limit (',char(176),')'],'Position',[10 75 123 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','String',['Propagation angle step (',char(176),')'],'Position',[10 45 127 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','String','Phase velocity sections (5^x)','Position',[10 15 144 13]);
        a.FrequencyLimitUI4 = uicontrol('Parent',a.p2UI4,'Style','edit','String',a.FrequencyLimit_Polar,'TooltipString','Enter the frequency up to which the dispersion curves shall be traced.','Position',[160 130 50 23],'Callback',@CallbackUI4,'Tag','1');
        a.FrequencyResolutionUI4 = uicontrol('Parent',a.p2UI4,'Style','edit','String',a.FrequencyResolution_Polar,'TooltipString','Enter the frequency step.','Position',[160 100 50 23],'Callback',@CallbackUI4,'Tag','2');
        uicontrol('Parent',a.p2UI4,'Style','popupmenu','String',{'180','90'},'TooltipString',['Select the propagation angle range. In case of 180 ',char(176),', it will be calculated',newline,'from 0 to 180 ',char(176),', and the result be copied to 180 to 360 ',char(176),'. In case of 90 ',char(176),',',newline,'it will be calculated from 0 to 90 ',char(176),', and the result be copied to the other',newline,'three quadrants. Notice that the latter approach is insufficient for unit',newline,'cells such as [0 45].'],'Position',[160 70 50 23],'Callback',@CallbackUI4,'Tag','33');
        uicontrol('Parent',a.p2UI4,'Style','edit','String',a.PropagationAngleStep_Polar,'TooltipString',['The polar dispersion curves are calculated for propagation',newline,'angles ranging from 0 to 90',char(176),'. Set the step.'],'Position',[160 40 50 23],'Callback',@CallbackUI4,'Tag','3');
        uicontrol('Parent',a.p2UI4,'Style','edit','String',a.PhaseVelocitySections_Polar,'TooltipString',['The fundamental dispersion curves are determined by sweeping the phase velocity at fixed',newline,'frequencies. The phase velocity sections as the power of five give the maximum number of',newline,'sections into which the phase velocity search interval is devided during the search for the',newline,'modal solution at a given frequency. A higher number increases the chance to find the',newline,'solution at the cost of processing time. Missing a solution is not necessarily critical since an',newline,'extrapolation routine replaces missing samples successfully, as long as not too many samples',newline,'are missing.'],'Position',[160 10 50 23],'Callback',@CallbackUI4,'Tag','4');

        a.p3UI4 = uipanel('Parent',Tab4,'Title','Mode selection','Units','pixels','Position',[10 255 225 115],'FontSize',10);
        uicontrol('Parent',a.p3UI4,'Style','text','String','S1/B2','Position',[10 75 32 13]);
        uicontrol('Parent',a.p3UI4,'Style','text','String','S0/B1','Position',[10 45 32 13]);
        uicontrol('Parent',a.p3UI4,'Style','text','String','A0/B0','Position',[10 15 32 13]);
        a.S0B1UI4 = uicontrol('Parent',a.p3UI4,'Style','checkbox','Value',a.S0_Polar,'TooltipString',['Check this to calculate the fundamental symmetric (S1)',newline,'or nonsymmetric (B2) mode.'],'Position',[160 70 50 23],'Callback',@CallbackUI4,'Tag','5');        
        a.SH0UI4 = uicontrol('Parent',a.p3UI4,'Style','checkbox','Value',a.SH0_Polar,'TooltipString',['Check this to calculate the fundamental symmetric (S0)',newline,'or nonsymmetric (B1) mode.'],'Position',[160 40 50 23],'Callback',@CallbackUI4,'Tag','6');
        a.A0B0UI4 = uicontrol('Parent',a.p3UI4,'Style','checkbox','Value',a.A0_Polar,'TooltipString',['Check this to calculate the fundamental antisymmetric (A0)',newline,'or nonsymmetric (B0) mode.'],'Position',[160 10 50 23],'Callback',@CallbackUI4,'Tag','7');

        %------------------------------------------------------------------
        a.c1UI4 = uicontrol('Parent',Tab4,'Style','pushbutton','String','Calculate','TooltipString','Start the dispersion curve tracing.','Position',[245 697 127 40],'FontSize',10,'Callback',@CallbackUI4,'Tag','8');
        a.c2UI4 = uicontrol('Parent',Tab4,'Style','pushbutton','String','Stop calculation','TooltipString','Stop the dispersion curve tracing.','Position',[382 697 127 40],'FontSize',10,'Callback',@CallbackUI4,'Tag','9');

        a.p4UI4 = uipanel('Parent',Tab4,'Title','Dispersion diagrams','Units','pixels','Position',[245 520 263 165],'FontSize',10);
        uicontrol('Parent',a.p4UI4,'Style','text','String','Quantity','Position',[10 125 42 13]);
        a.Option1TextUI4 = uicontrol('Parent',a.p4UI4,'Style','text','String','Bulk velocities','Position',[10 95 70 13]);
        uicontrol('Parent',a.p4UI4,'Style','text','String','Frequency (kHz)','Position',[10 65 83 13]);
        uicontrol('Parent',a.p4UI4,'Style','popupmenu','String',{'Phase velocity (m/ms)','Energy velocity (m/ms)',['Propagation time (',char(181),'s)'],['Coincidence angle (',char(176),')'],'Wavelength (mm)','Wavenumber (rad/mm)'},'TooltipString','Select which quantity to plot in the polar dispersion diagram.','Position',[117 120 130 23],'Callback',@CallbackUI4,'Tag','10');
        a.Option1UI4 = uicontrol('Parent',a.p4UI4,'Style','checkbox','Value',a.BulkVelocities_Polar,'TooltipString','Check this to show the bulk wave velocities.','Position',[117 90 50 23],'Callback',@CallbackUI4,'Tag','11');
        a.FrequencyUI4 = uicontrol('Parent',a.p4UI4,'Style','popupmenu','String',{''},'TooltipString','Select the frequency for which the polar dispersion curves shall be plotted.','Position',[117 60 65 23],'Callback',@CallbackUI4,'Tag','12');
        a.PlotUI4 = uicontrol('Parent',a.p4UI4,'Style','pushbutton','String','Plot','TooltipString',['Plot the polar dispersion diagram. If you have checked ''Export plots''',newline,'in the export settings, the plot will be exported automatically.'],'Position',[117 15 65 33],'FontSize',10,'Callback',@CallbackUI4,'Tag','13');

        a.p5UI4 = uipanel('Parent',Tab4,'Title','Export settings','Units','pixels','Position',[245 335 473 180],'FontSize',10);
        uicontrol('Parent',a.p5UI4,'Style','text','String','Export plots','Position',[10 140 59 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','String','Crop plots','Position',[10 120 51 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','String','PDF','Position',[124 140 21 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','String','PNG','Position',[124 120 23 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','String','PNG resolution (dpi)','Position',[10 93 98 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','String','Directory','Position',[10 20 46 13]);
        uicontrol('Parent',a.p5UI4,'Style','checkbox','Value',a.ExportPlots_Polar,'TooltipString',['Check this in order to export plots automatically upon',newline,'pressing the plot button. You can also export plots',newline,'manually by using the ''File'' menu inside the plot figure.'],'Position',[80 135 20 23],'Callback',@CallbackUI4,'Tag','14');
        uicontrol('Parent',a.p5UI4,'Style','checkbox','Value',a.CropPlots_Polar,'TooltipString','Plots will be cropped tightly with minimal white space.','Position',[80 115 20 23],'Callback',@CallbackUI4,'Tag','32');
        uicontrol('Parent',a.p5UI4,'Style','checkbox','Value',a.PDF_Polar,'TooltipString','Check this to export a plot as pdf.','Position',[160 135 20 23],'Callback',@CallbackUI4,'Tag','15');
        uicontrol('Parent',a.p5UI4,'Style','checkbox','Value',a.PNG_Polar,'TooltipString','Check this to export a plot as png.','Position',[160 115 20 23],'Callback',@CallbackUI4,'Tag','16');
        uicontrol('Parent',a.p5UI4,'Style','edit','String',a.PNGresolution_Polar,'TooltipString','Enter the resolution of the png image.','Position',[125 88 50 23],'Callback',@CallbackUI4,'Tag','17');
        uicontrol('Parent',a.p5UI4,'Style','pushbutton','String','*.mat','TooltipString','Export the data as Matlab''s mat-file.','Position',[220 80 40 33],'Callback',@CallbackUI4,'Tag','18');
        uicontrol('Parent',a.p5UI4,'Style','pushbutton','String','*.xlsx','TooltipString','Export the data as Excel sheet.','Position',[270 80 40 33],'Callback',@CallbackUI4,'Tag','19');
        uicontrol('Parent',a.p5UI4,'Style','pushbutton','String','*.txt','TooltipString','Export the data as txt-file.','Position',[320 80 40 33],'Callback',@CallbackUI4,'Tag','20');
        uicontrol('Parent',a.p5UI4,'Style','edit','String',a.FileName_Polar,'TooltipString',['Specify the name of the plot to be exported. The',newline,'dispersion curve raw data are named automatically.'],'Position',[70 45 385 23],'Callback',@CallbackUI4,'Tag','21');
        uicontrol('Parent',a.p5UI4,'Style','edit','String',a.Directory,'TooltipString',['Specify the directory to which the plot and the',newline,'dispersion curve raw data shall be exported.'],'Position',[70 15 385 23],'Callback',@CallbackUI4,'Tag','22');
        
        %------------------------------------------------------------------
        a.p6UI4 = uipanel('Parent',Tab4,'Title','Plot layout settings','Units','pixels','Position',[518 520 200 235],'FontSize',10);
        uicontrol('Parent',a.p6UI4,'Style','text','String','Title','Position',[10 175 21 13]);
        uicontrol('Parent',a.p6UI4,'Style','text','String','Curve line width','Position',[10 145 80 13]);
        a.TitleUI4 = uicontrol('Parent',a.p6UI4,'Style','popupmenu','String',{'with layup','without layup','no title'},'TooltipString',['Choose to plot the title with the layup ',newline,'included, without layup, or no title at all.'],'Position',[95 170 90 23],'Callback',@CallbackUI4,'Tag','23');
        a.LineWidthUI4 = uicontrol('Parent',a.p6UI4,'Style','edit','String',a.LineWidth_Polar,'TooltipString','Enter the curve line width.','Position',[95 140 50 23],'Callback',@CallbackUI4,'Tag','24');

        a.p7UI4 = uipanel('Parent',Tab4,'Title','Font size','Units','pixels','Position',[528 534 180 115],'FontSize',9);
        uicontrol('Parent',a.p7UI4,'Style','text','String','Title','Position',[10 75 21 13]);
        uicontrol('Parent',a.p7UI4,'Style','text','String','Axes ticks','Position',[10 45 53 13]);
        uicontrol('Parent',a.p7UI4,'Style','text','String','Mode labels','Position',[10 15 59 13]);
        a.TitleFontSizeUI4 = uicontrol('Parent',a.p7UI4,'Style','edit','String',a.TitleFontSize_Polar,'TooltipString','Set the title font size.','Position',[85 70 50 23],'Callback',@CallbackUI4,'Tag','28');
        a.AxesTickFontSizeUI4 = uicontrol('Parent',a.p7UI4,'Style','edit','String',a.AxesTickFontSize_Polar,'TooltipString','Set the axes tick label font size.','Position',[85 40 50 23],'Callback',@CallbackUI4,'Tag','29');
        a.ModeLabelFontSizeUI4 = uicontrol('Parent',a.p7UI4,'Style','edit','String',a.ModeLabelFontSize_Polar,'TooltipString','Set the mode label font size.','Position',[85 10 50 23],'Callback',@CallbackUI4,'Tag','30');

        a.c3UI4 = uicontrol('Parent',Tab4,'Style','pushbutton','String','Default','TooltipString','Reset the plot layout settings to the default.','Position',[730 693 65 33],'FontSize',10,'Callback',@CallbackUI4,'Tag','31');
    
        a.a1UI4 = axes('Parent',Tab4,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading20.png')    
    end
    if  1 % Tab8_bulk waves
        a.OutputWindowUI8 = uicontrol('Parent',Tab8,'Style','text','String','','Position',[10 10 220 737],'BackgroundColor','white');
        [a.BulkWaves,a.X,a.Y] = Computer_Isotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Theta_Bulk2,a.OutputWindowUI8);        
        a.OutputWindowUI8.String = ''; clc
        
        a.p1UI8 = uipanel('Parent',Tab8,'Title','Elastic waves in bulk material','Units','pixels','Position',[237 545 639 210],'FontSize',10);
        uicontrol('Parent',a.p1UI8,'Style','text','String','Class','Position',[10 175 29 13]);
        uicontrol('Parent',a.p1UI8,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic','Isotropic'},'Value',a.MaterialType_Bulk,'TooltipString','Select a material symmetry class.','Position',[10 145 130 23],'Callback',@CallbackUI8,'Tag','43');
        uicontrol('Parent',a.p1UI8,'Style','text','String','Material','Position',[10 125 39 13]);
        a.MaterialUI8 = uicontrol('Parent',a.p1UI8,'Style','popupmenu','String',fieldnames(a.Materials.Orthotropic),'TooltipString','Select a material.','Position',[10 95 130 23],'Callback',@CallbackUI8,'Tag','1');
        uicontrol('Parent',a.p1UI8,'Style','text','String','Quantity','Position',[10 75 42 13]);
        uicontrol('Parent',a.p1UI8,'Style','popupmenu','String',{'Phase velocity (m/ms)','Group velocity (m/ms)','Slowness (ms/m)',['Polarization skew (',char(176),')'],['Energy skew (',char(176),')']},'TooltipString','Select which quantity to plot.','Position',[10 45 130 23],'Callback',@CallbackUI8,'Tag','2');
        uicontrol('Parent',a.p1UI8,'Style','text','String',[char(916),char(920),' (',char(176),')'],'Position',[10 15 33 13]);
        uicontrol('Parent',a.p1UI8,'Style','edit','String',a.ThetaStep_Bulk1,'TooltipString','Enter the Theta step.','Position',[80 10 50 23],'Callback',@CallbackUI8,'Tag','3');
        
        a.p2UI8 = uipanel('Parent',Tab8,'Title','2-D profiles','Units','pixels','Position',[394 555 85 185],'FontSize',9);
        uicontrol('Parent',a.p2UI8,'Style','text','String','Plane','Position',[10 145 28 13]);
        uicontrol('Parent',a.p2UI8,'Style','text','String',[char(934),' (',char(176),')'],'Position',[10 90 25 13]);
        uicontrol('Parent',a.p2UI8,'Style','popupmenu','String',{'1-3','1-2'},'TooltipString',['Select to plot the 1-3-plane or the 1-2-plane. In the upper sketch to the',newline,'right, the x''3 and x''2-axes are swapped according to your selection.'],'Position',[10 115 50 23],'Callback',@CallbackUI8,'Tag','44');                        
        uicontrol('Parent',a.p2UI8,'Style','edit','String',a.Phi_Bulk11,'TooltipString','Enter Phi.','Position',[10 60 50 23],'Callback',@CallbackUI8,'Tag','4');
        a.Plot1UI8 = uicontrol('Parent',a.p2UI8,'Style','pushbutton','String','Plot','TooltipString',['Plot the profiles. If you have checked ''Export plots'' in the',newline,'export settings, the plot will be exported automatically.'],'Position',[10 15 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','5');

        a.p3UI8 = uipanel('Parent',Tab8,'Title','3-D surfaces and vectors','Units','pixels','Position',[489 555 377 185],'FontSize',9);
        uicontrol('Parent',a.p3UI8,'Style','text','String',[char(916),char(934),' (',char(176),')'],'Position',[10 145 33 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','String','Mode','Position',[10 65 28 13]);
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.PhiStep_Bulk,'TooltipString','Enter the Phi step.','Position',[50 140 50 23],'Callback',@CallbackUI8,'Tag','6');
        a.CalculateUI8 = uicontrol('Parent',a.p3UI8,'Style','pushbutton','String','Calculate','Position',[50 95 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','7');
        uicontrol('Parent',a.p3UI8,'Style','popupmenu','String',{'Longitudinal','Fast shear','Slow shear'},'TooltipString','Select which mode to plot.','Position',[50 60 80 23],'Callback',@CallbackUI8,'Tag','8');                        
        a.Plot2UI8 = uicontrol('Parent',a.p3UI8,'Style','pushbutton','String','Plot','TooltipString',['Plot the surface. If you have checked ''Export plots'' in the',newline,'export settings, the plot will be exported automatically.'],'Position',[50 15 65 33],'Enable','off','FontSize',10,'Callback',@CallbackUI8,'Tag','9');

        uicontrol('Parent',a.p3UI8,'Style','text','String',[char(934),' (',char(176),')'],'Position',[136 145 25 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','String',[char(920),' (',char(176),')'],'Position',[136 115 25 13]);
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.Phi_Bulk21,'TooltipString','Enter propagation direction Phi.','Position',[166 140 50 23],'Callback',@CallbackUI8,'Tag','10');
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.Theta_Bulk1,'TooltipString','Enter propagation direction Theta.','Position',[166 110 50 23],'Callback',@CallbackUI8,'Tag','11');
        a.Plot3UI8 = uicontrol('Parent',a.p3UI8,'Style','pushbutton','String','Plot','TooltipString',['Plot the vectors. If you have checked ''Export plots'' in the',newline,'export settings, the plot will be exported automatically.'],'Position',[166 15 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','12');
                        
        uicontrol('Parent',a.p3UI8,'Style','text','String',['View ',char(934),' (',char(176),')'],'Position',[246 145 54 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','String',['View ',char(920),' (',char(176),')'],'Position',[246 115 54 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','String','Marker size','Position',[246 85 58 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','String','Colorbar','Position',[246 55 43 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','String','x-pos. (0-1)','Position',[246 40 61 13]);
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.ViewPhi_Bulk1,'TooltipString','Enter the view angle Phi.','Position',[310 140 50 23],'Callback',@CallbackUI8,'Tag','13');
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.ViewTheta_Bulk1,'TooltipString','Enter the view angle Theta.','Position',[310 110 50 23],'Callback',@CallbackUI8,'Tag','14');
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.MarkerSize_Bulk,'TooltipString','Enter the size of the markers.','Position',[310 80 50 23],'Callback',@CallbackUI8,'Tag','15');
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.ColorbarX_Bulk,'TooltipString','Shift the colorbar left-right.','Position',[310 50 50 23],'Callback',@CallbackUI8,'Tag','16');

        %------------------------------------------------------------------
        a.p4UI8 = uipanel('Parent',Tab8,'Title','Bulk waves on interfaces','Units','pixels','Position',[237 275 245 265],'FontSize',10);
        uicontrol('Parent',a.p4UI8,'Style','text','String','Fluid','Position',[10 225 24 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','String','Class','Position',[10 195 29 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','String','Solid','Position',[10 165 25 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','String',[char(934),' (',char(176),')'],'Position',[10 135 25 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','String',[char(920),'i (',char(176),')'],'Position',[10 105 27 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','String',[char(916),char(920),' (',char(176),')'],'Position',[10 75 33 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','String',['View ',char(934),' (',char(176),')'],'Position',[10 45 54 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','String',['View ',char(920),' (',char(176),')'],'Position',[10 15 54 13]);
        
        a.CouplantUI8 = uicontrol('Parent',a.p4UI8,'Style','popupmenu','String',fieldnames(a.Materials.Fluid),'TooltipString',['Select a default fluid from which a',newline,' plane wave impinges on the solid.'],'Position',[80 220 150 23],'Callback',@CallbackUI8,'Tag','17');
        uicontrol('Parent',a.p4UI8,'Style','popupmenu','String',{'Isotropic','Cubic','Transversely isotropic','Orthotropic'},'TooltipString','Select the material type of the solid.','Position',[80 190 150 23],'Callback',@CallbackUI8,'Tag','19');
        a.SolidUI8 = uicontrol('Parent',a.p4UI8,'Style','popupmenu','String',fieldnames(a.Materials.Isotropic),'TooltipString','Select a material.','Position',[80 160 150 23],'Callback',@CallbackUI8,'Tag','20');
        a.Phi2UI8 = uicontrol('Parent',a.p4UI8,'Style','edit','String',a.Phi_Bulk2,'TooltipString','Enter Phi.','Position',[80 130 50 23],'Enable','off','Callback',@CallbackUI8,'Tag','21');
        uicontrol('Parent',a.p4UI8,'Style','edit','String',a.Theta_Bulk2,'TooltipString','Enter the incidence angle Theta.','Position',[80 100 50 23],'Callback',@CallbackUI8,'Tag','22');
        a.ThetaStep2UI8 = uicontrol('Parent',a.p4UI8,'Style','edit','String',a.ThetaStep_Bulk2,'TooltipString','Enter the Theta step for plotting the reflection and transmission coefficients.','Position',[80 70 50 23],'Callback',@CallbackUI8,'Tag','23');
        a.ViewPhi2UI8 = uicontrol('Parent',a.p4UI8,'Style','edit','String',a.ViewPhi_Bulk2,'TooltipString','Enter the view angle Phi.','Position',[80 40 50 23],'Enable','off','Callback',@CallbackUI8,'Tag','24');
        a.ViewTheta2UI8 = uicontrol('Parent',a.p4UI8,'Style','edit','String',a.ViewTheta_Bulk2,'TooltipString','Enter the view angle Theta.','Position',[80 10 50 23],'Enable','off','Callback',@CallbackUI8,'Tag','25');
        a.Plot2DUI8 = uicontrol('Parent',a.p4UI8,'Style','pushbutton','String','Plot 2-D','TooltipString',['Plot the bulk waves. If you have checked ''Export plots'' in ',newline,' the export settings, the plot will be exported automatically.'],'Position',[165 90 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','26');        
        a.Plot3DUI8 = uicontrol('Parent',a.p4UI8,'Style','pushbutton','String','Plot 3-D','TooltipString',['Plot the bulk waves. If you have checked ''Export plots'' in',newline,' the export settings, the plot will be exported automatically.'],'Position',[165 50 65 33],'FontSize',10,'Enable','off','Callback',@CallbackUI8,'Tag','27');
        a.PlotRTUI8 = uicontrol('Parent',a.p4UI8,'Style','pushbutton','String','Plot R,T','TooltipString',['Plot reflected and transmitted energy',newline,'coefficients versus incidence angle.'],'Position',[165 10 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','46');

        %------------------------------------------------------------------
        a.p5UI8 = uipanel('Parent',Tab8,'Title','Plot layout settings','Units','pixels','Position',[492 215 170 325],'FontSize',10);
        uicontrol('Parent',a.p5UI8,'Style','text','String','Title','Position',[10 285 21 13]);
        a.TitleUI8 = uicontrol('Parent',a.p5UI8,'Style','checkbox','Value',a.Title_Bulk,'TooltipString','Check this in order to show the plot title.','Position',[95 280 20 23],'Callback',@CallbackUI8,'Tag','28');      

        a.p6UI8 = uipanel('Parent',Tab8,'Title','Line width','Units','pixels','Position',[502 375 150 115],'FontSize',9);
        uicontrol('Parent',a.p6UI8,'Style','text','String','Box','Position',[10 75 21 13]);
        uicontrol('Parent',a.p6UI8,'Style','text','String','Profile','Position',[10 45 32 13]);
        uicontrol('Parent',a.p6UI8,'Style','text','String','Bulk wave','Position',[10 15 53 13]);
        a.BoxLineWidthUI8 = uicontrol('Parent',a.p6UI8,'Style','edit','String',a.BoxLineWidth_Bulk,'TooltipString','Enter the box line width.','Position',[85 70 50 23],'Callback',@CallbackUI8,'Tag','29');
        a.LineWidthUI8 = uicontrol('Parent',a.p6UI8,'Style','edit','String',a.LineWidth_Bulk,'TooltipString','Enter the profile line width.','Position',[85 40 50 23],'Callback',@CallbackUI8,'Tag','30');
        a.WaveVectorLineWidthUI8 = uicontrol('Parent',a.p6UI8,'Style','edit','String',a.WaveVectorLineWidth_Bulk,'TooltipString','Enter the bulk wave line width.','Position',[85 10 50 23],'Callback',@CallbackUI8,'Tag','31');
        
        a.p7UI8 = uipanel('Parent',Tab8,'Title','Font size','Units','pixels','Position',[502 225 150 145],'FontSize',9);
        uicontrol('Parent',a.p7UI8,'Style','text','String','Title','Position',[10 105 21 13]);
        uicontrol('Parent',a.p7UI8,'Style','text','String','Axes labels','Position',[10 75 59 13]);
        uicontrol('Parent',a.p7UI8,'Style','text','String','Axes ticks','Position',[10 45 53 13]);
        uicontrol('Parent',a.p7UI8,'Style','text','String','Mode labels','Position',[10 15 59 13]);
        a.TitleFontSizeUI8 = uicontrol('Parent',a.p7UI8,'Style','edit','String',a.TitleFontSize_Bulk,'TooltipString','Set the title font size.','Position',[85 100 50 23],'Callback',@CallbackUI8,'Tag','32');
        a.AxesLabelFontSizeUI8 = uicontrol('Parent',a.p7UI8,'Style','edit','String',a.AxesLabelFontSize_Bulk,'TooltipString','Set the axes label font size.','Position',[85 70 50 23],'Callback',@CallbackUI8,'Tag','33');
        a.AxesTickFontSizeUI8 = uicontrol('Parent',a.p7UI8,'Style','edit','String',a.AxesTickFontSize_Bulk,'TooltipString','Set the axes tick label font size.','Position',[85 40 50 23],'Callback',@CallbackUI8,'Tag','34');
        a.ModeLabelFontSizeUI8 = uicontrol('Parent',a.p7UI8,'Style','edit','String',a.ModeLabelFontSize_Bulk,'TooltipString','Set the mode label font size.','Position',[85 10 50 23],'Callback',@CallbackUI8,'Tag','35');        

        %------------------------------------------------------------------
        a.p8UI8 = uipanel('Parent',Tab8,'Title','Export settings','Units','pixels','Position',[237 70 425 140],'FontSize',10);
        uicontrol('Parent',a.p8UI8,'Style','text','String','Export plots','Position',[10 100 59 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','String','Crop plots','Position',[10 80 51 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','String','PDF','Position',[124 100 21 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','String','PNG','Position',[124 80 23 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','String','PNG resolution (dpi)','Position',[235 80 98 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','String','Directory','Position',[10 20 46 13]);
        uicontrol('Parent',a.p8UI8,'Style','checkbox','Value',a.ExportPlots_Bulk,'TooltipString',['Check this in order to export plots automatically upon',newline,'pressing the respective plot button. You can also export plots',newline,'manually by using the ''File'' menu inside the plot figure.'],'Position',[80 95 20 23],'Callback',@CallbackUI8,'Tag','36');
        uicontrol('Parent',a.p8UI8,'Style','checkbox','Value',a.CropPlots_Bulk,'TooltipString','Plots will be cropped tightly with minimal white space.','Position',[80 75 20 23],'Callback',@CallbackUI8,'Tag','45');
        uicontrol('Parent',a.p8UI8,'Style','checkbox','Value',a.PDF_Bulk,'TooltipString','Check this to export a plot as pdf.','Position',[160 95 20 23],'Callback',@CallbackUI8,'Tag','37');
        uicontrol('Parent',a.p8UI8,'Style','checkbox','Value',a.PNG_Bulk,'TooltipString','Check this to export a plot as png.','Position',[160 75 20 23],'Callback',@CallbackUI8,'Tag','38');
        uicontrol('Parent',a.p8UI8,'Style','edit','String',a.PNGresolution_Bulk,'TooltipString','Enter the resolution of the png image.','Position',[350 75 50 23],'Callback',@CallbackUI8,'Tag','39');
        uicontrol('Parent',a.p8UI8,'Style','edit','String',a.FileName_Bulk,'TooltipString','Specify the name of the plots to be exported.','Position',[70 45 330 23],'Callback',@CallbackUI8,'Tag','40');
        uicontrol('Parent',a.p8UI8,'Style','edit','String',a.Directory,'TooltipString','Specify the directory to which plots shall be exported.','Position',[70 15 330 23],'Callback',@CallbackUI8,'Tag','41');        

        a.c1UI8 = uicontrol('Parent',Tab8,'Style','pushbutton','String','Default','TooltipString','Reset the plot layout settings to the default.','Position',[675 500 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','42');    

        a.a1UI8 = axes('Parent',Tab8,'Units','pixels','Position',[880 390 312 367]);
        imshow('Bulk1.png')
        a.a2UI8 = axes('Parent',Tab8,'Units','pixels','Position',[806 5 386 381]);
        imshow('Bulk2.png')        
    end
    if  1 % Tab6_laminate stiffness
        a.p1UI6 = uipanel('Parent',Tab6,'Title','','Units','pixels','Position',[0 0 1198 1000],'FontSize',10);
        
        uicontrol('Parent',a.p1UI6,'Style','text','String','Specimen','Position',[20 565+170 58 15],'FontSize',9);
        uicontrol('Parent',a.p1UI6,'Style','pushbutton','String','Edit','TooltipString','Press to define your specimen.','Position',[81.5 565+125 127 40],'FontSize',10,'Callback',@SpecimenSettingsUI6_Callback);
        uicontrol('Parent',a.p1UI6,'Style','text','String','Material','Position',[20 565+95 39 13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','Unit cell','Position',[20 565+65 39 13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String',['Azimuthal angle (',char(176),')'],'Position',[20 565+35 93 13]);
        a.MaterialNameUI6 = uicontrol('Parent',a.p1UI6,'Style','text','String',a.Material3{1}.Name,'Position',[81.5 565+90 190 23],'BackgroundColor','white');
        a.LayupUI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','[0]','Position',[81.5 565+60 190 23],'BackgroundColor','white');
        uicontrol('Parent',a.p1UI6,'Style','edit','String',a.PropagationAngle,'TooltipString','Enter the azimuthal angle with respect to the fiber orientations.','Position',[140 565+30 50 23],'Callback',@CallbackUI6,'Tag','2');

        uicontrol('Parent',a.p1UI6,'Style','text','String','Laminate stiffness components (homogenized stiffness tensor) (GPa)','Position',[280+40 495+240 385 15],'FontSize',9);
        a.C11UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+40 495+210 60 23],'BackgroundColor','white');
        a.C12UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+105 495+210 60 23],'BackgroundColor','white');
        a.C13UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+170 495+210 60 23],'BackgroundColor','white');
        uicontrol('Parent',a.p1UI6,'Style','text','String','0','Position',[280+235 495+210 60 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','0','Position',[280+300 495+210 60 23],'BackgroundColor',[.88 .88 .88]);
        a.C16UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+365 495+210 60 23],'BackgroundColor','white');
        a.C22UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+105 495+180 60 23],'BackgroundColor','white');
        a.C23UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+170 495+180 60 23],'BackgroundColor','white');
        uicontrol('Parent',a.p1UI6,'Style','text','String','0','Position',[280+235 495+180 60 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','0','Position',[280+300 495+180 60 23],'BackgroundColor',[.88 .88 .88]);
        a.C26UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+365 495+180 60 23],'BackgroundColor','white');          
        a.C33UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+170 495+150 60 23],'BackgroundColor','white');
        uicontrol('Parent',a.p1UI6,'Style','text','String','0','Position',[280+235 495+150 60 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','0','Position',[280+300 495+150 60 23],'BackgroundColor',[.88 .88 .88]);
        a.C36UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+365 495+150 60 23],'BackgroundColor','white');        
        a.C44UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+235 495+120 60 23],'BackgroundColor','white');
        a.C45UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+300 495+120 60 23],'BackgroundColor','white');
        uicontrol('Parent',a.p1UI6,'Style','text','String','0','Position',[280+365 495+120 60 23],'BackgroundColor',[.88 .88 .88]);        
        a.C55UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+300 495+90 60 23],'BackgroundColor','white');
        uicontrol('Parent',a.p1UI6,'Style','text','String','0','Position',[280+365 495+90 61 23],'BackgroundColor',[.88 .88 .88]);         
        a.C66UI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','','Position',[280+365 495+60 61 23],'BackgroundColor','white');      

        uicontrol('Parent',a.p1UI6,'Style','pushbutton','String','*.mat','TooltipString','Export the data as Matlab''s mat-file.','Position',[750 590 40 33],'Callback',@CallbackUI6,'Tag','3');
        uicontrol('Parent',a.p1UI6,'Style','pushbutton','String','*.xlsx','TooltipString','Export the data as Excel sheet.','Position',[800 590 40 33],'Callback',@CallbackUI6,'Tag','4');
        uicontrol('Parent',a.p1UI6,'Style','pushbutton','String','*.txt','TooltipString','Export the data as txt-file.','Position',[850 590 40 33],'Callback',@CallbackUI6,'Tag','5');
        uicontrol('Parent',a.p1UI6,'Style','edit','String',a.Directory,'TooltipString',['Specify the directory to which the laminate',newline,'stiffness matrix shall be exported.'],'Position',[750 495+60 410 23],'Callback',@CallbackUI6,'Tag','6');        

        uicontrol('Parent',a.p1UI6,'Style','text','String','C11','Position',[20,130+365,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C12','Position',[20,130+335,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C13','Position',[20,130+305,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C16','Position',[20,130+275,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C22','Position',[20,130+245,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C23','Position',[20,130+215,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C26','Position',[20,130+185,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C33','Position',[20,130+155,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C36','Position',[20,130+125,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C44','Position',[20,130+95,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C45','Position',[20,130+65,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C55','Position',[20,130+35,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','text','String','C66','Position',[20,130+5,21,13]);
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(1),'Position',[60 130+360 15 23],'Callback',@CallbackUI6,'Tag','7','backgroundcolor','r');
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(2),'Position',[60 130+330 15 23],'Callback',@CallbackUI6,'Tag','8','backgroundcolor',[.13 .55 .13]);
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(3),'Position',[60 130+300 15 23],'Callback',@CallbackUI6,'Tag','9','backgroundcolor','b');
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(4),'Position',[60 130+270 15 23],'Callback',@CallbackUI6,'Tag','10','backgroundcolor','k');
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(5),'Position',[60 130+240 15 23],'Callback',@CallbackUI6,'Tag','11','backgroundcolor','m');
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(6),'Position',[60 130+210 15 23],'Callback',@CallbackUI6,'Tag','12','backgroundcolor','c');
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(7),'Position',[60 130+180 15 23],'Callback',@CallbackUI6,'Tag','13','backgroundcolor',[1 .7 0]);
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(8),'Position',[60 130+150 15 23],'Callback',@CallbackUI6,'Tag','14','backgroundcolor',[.55 .27 .13]);
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(9),'Position',[60 130+120 15 23],'Callback',@CallbackUI6,'Tag','15','backgroundcolor',[.5 0 1]);
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(10),'Position',[60 130+90 15 23],'Callback',@CallbackUI6,'Tag','16','backgroundcolor',[.5 .5 .5]);
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(11),'Position',[60 130+60 15 23],'Callback',@CallbackUI6,'Tag','17','backgroundcolor','r');
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(12),'Position',[60 130+30 15 23],'Callback',@CallbackUI6,'Tag','18','backgroundcolor',[.13 .55 .13]);
        uicontrol('Parent',a.p1UI6,'Style','checkbox','Value',a.Polar(13),'Position',[60 130+0 15 23],'Callback',@CallbackUI6,'Tag','19','backgroundcolor','b');
        
        a.CLaminate = Computer_LaminateStiffness(a.Material3,a.LayerOrientations3,a.LayerThicknesses3,a.PropagationAngle3);
        a.CLP = Computer_LaminateStiffnessPolar(a.Material3,a.LayerOrientations3,a.LayerThicknesses3,1);
        
        a.C11UI6.String = a.CLaminate(1,1)/1e9;
        a.C12UI6.String = a.CLaminate(1,2)/1e9;
        a.C13UI6.String = a.CLaminate(1,3)/1e9;
        a.C16UI6.String = a.CLaminate(1,6)/1e9;
        a.C22UI6.String = a.CLaminate(2,2)/1e9;
        a.C23UI6.String = a.CLaminate(2,3)/1e9;
        a.C26UI6.String = a.CLaminate(2,6)/1e9;
        a.C33UI6.String = a.CLaminate(3,3)/1e9;
        a.C36UI6.String = a.CLaminate(3,6)/1e9;
        a.C44UI6.String = a.CLaminate(4,4)/1e9;
        a.C45UI6.String = a.CLaminate(4,5)/1e9;
        a.C55UI6.String = a.CLaminate(5,5)/1e9;
        a.C66UI6.String = a.CLaminate(6,6)/1e9;
        
        a.h1 = LaminateStiffness_Internal(a); 
    end
    if  1 % Tab5_material editor
        a.p1UI5 = uipanel('Parent',Tab5,'Title','Isotropic materials','Units','pixels','Position',[10 220 440 535],'FontSize',10);
        uicontrol('Parent',a.p1UI5,'Style','pushbutton','String','?','Position',[395 485 23 23],'FontSize',12,'FontWeight','bold','foregroundcolor',[1 1 1],'backgroundcolor',[.4 .4 .4],'Callback',@CallbackUI5,'Tag','68');
        uicontrol('Parent',a.p1UI5,'Style','text','String','Material','Position',[10 180+310 39 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Mass density (kg/m3)','Position',[10 180+280 105 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Engineering constants','Position',[10 160+255 126 15],'FontSize',9);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Real part','Position',[60 160+235 45 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Imaginary part','Position',[135 160+235 70 15]);        
        uicontrol('Parent',a.p1UI5,'Style','text','String','E (GPa)','Position',[10 180+190 39 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','v','Position',[10 180+160 6 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Stiffness components (GPa)','Position',[250 160+255 157 15],'FontSize',9);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Real part','Position',[280 160+235 45 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Imaginary part','Position',[355 160+235 70 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','C11','Position',[250 180+190 21 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','C66','Position',[250 180+160 21 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Bulk waves','Position',[10 180+115 64 15],'FontSize',9);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Longitudinal velocity (m/s)','Position',[10 180+90 127 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Shear velocity (m/s)','Position',[10 180+60 99 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','Attenuation unit','Position',[10 180+30 77 13]);
        a.LongitudinalAttenuationTextUI5 = uicontrol('Parent',a.p1UI5,'Style','text','String',['Longitudinal attenuation (Np/',char(955),')'],'Position',[10 180+0 148 13]);
        a.TransverseAttenuationTextUI5 = uicontrol('Parent',a.p1UI5,'Style','text','String',['Shear attenuation (Np/',char(955),')'],'Position',[10 180-30 120 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','String','At frequency (kHz)','Position',[10 180-60 95 13]);
        a.Material1UI5 = uicontrol('Parent',a.p1UI5,'Style','popupmenu','String',fieldnames(a.Materials.Isotropic),'TooltipString','Select a material to display its parameters.','Position',[135 180+305 220 23],'Callback',@CallbackUI5,'Tag','1');
        a.Density1UI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.Density,'TooltipString','Enter the mass density.','Position',[135 180+275 65 23],'Callback',@CallbackUI5,'Tag','2');
        a.YoungsModulusUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',real(a.Material1.YoungsModulus)/1e9,'TooltipString','Enter the real part of Young''s modulus.','Position',[60 160+205 65 23],'Callback',@CallbackUI5,'Tag','3');
        a.YoungsModulusImagUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',imag(a.Material1.YoungsModulus)/1e9,'TooltipString','Enter the imaginary part of Young''s modulus.','Position',[135 160+205 65 23],'Callback',@CallbackUI5,'Tag','4');        
        a.PoissonsNumberUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',real(a.Material1.PoissonsNumber),'TooltipString','Enter the real part of Poisson''s ratio.','Position',[60 160+175 65 23],'Callback',@CallbackUI5,'Tag','5');
        a.PoissonsNumberImagUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',imag(a.Material1.PoissonsNumber),'TooltipString','Enter the imaginary part of Poisson''s ratio.','Position',[135 160+175 65 23],'Callback',@CallbackUI5,'Tag','6');
        a.C11IsoUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',real(a.Material1.C(1,1))/1e9,'TooltipString','Enter the real stiffness component.','Position',[280 160+205 65 23],'Callback',@CallbackUI5,'Tag','64');
        a.C11ImagIsoUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',imag(a.Material1.C(1,1))/1e9,'TooltipString','Enter the imaginary stiffness component.','Position',[355 160+205 65 23],'Callback',@CallbackUI5,'Tag','65');
        a.C66IsoUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',real(a.Material1.C(6,6))/1e9,'TooltipString','Enter the real stiffness component.','Position',[280 160+175 65 23],'Callback',@CallbackUI5,'Tag','66');        
        a.C66ImagIsoUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',imag(a.Material1.C(6,6))/1e9,'TooltipString','Enter the imaginary stiffness component.','Position',[355 160+175 65 23],'Callback',@CallbackUI5,'Tag','67');
        a.LongitudinalVelocityUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.LongitudinalVelocity,'TooltipString','Enter the longitudinal bulk velocity.','Position',[165 180+85 65 23],'Callback',@CallbackUI5,'Tag','7');
        a.TransverseVelocityUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.TransverseVelocity,'TooltipString','Enter the shear bulk velocity.','Position',[165 180+55 65 23],'Callback',@CallbackUI5,'Tag','8');
        uicontrol('Parent',a.p1UI5,'Style','popupmenu','String',{['Np/',char(955)],['dB/',char(955)],'Np/m','dB/m'},'TooltipString',['Select in which units you want to enter the attenuation of the bulk waves. ',newline,'If you select ''Np/m'' or ''dB/m'', you need to enter at which frequency the',newline,'attenuation was measured. The DC assumes a linear increase of attenuation',newline,'with frequency, i.e., a damping loss which is constant per wavelength',newline,'(hysteretic damping).'],'Position',[165 180+25 65 23],'Callback',@CallbackUI5,'Tag','9');
        a.LongitudinalAttenuationUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.LongitudinalAttenuation,'TooltipString','The attenuation of longitudinal bulk waves.','Position',[165 180-5 65 23],'Callback',@CallbackUI5,'Tag','10');
        a.TransverseAttenuationUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.LongitudinalAttenuation,'TooltipString','The attenuation of shear bulk waves.','Position',[165 180-35 65 23],'Callback',@CallbackUI5,'Tag','11');
        a.AtFrequency1UI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.AtFrequency1,'TooltipString','Enter at which frequency the above attenuations were measured.','Position',[165 180-65 65 23],'Enable','off','Callback',@CallbackUI5,'Tag','12');

        uicontrol('Parent',a.p1UI5,'Style','text','String','New material''s name','Position',[10 110-45 102 13]);
        a.Name1UI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.MaterialName4A,'TooltipString',['Enter a new material''s name and press ''Save'' or enter',newline,'an available material''s name and press ''Delete''.'],'Position',[135 110-50 220 23],'Callback',@CallbackUI5,'Tag','13');
        uicontrol('Parent',a.p1UI5,'Style','pushbutton','String','Save material','TooltipString','Save the material to the isotropic materials list.','Position',[135 110-95 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','14');
        uicontrol('Parent',a.p1UI5,'Style','pushbutton','String','Delete material','TooltipString','Remove the material from the isotropic materials list.','Position',[255 110-95 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','15');

        %------------------------------------------------------------------        
        a.p2UI5 = uipanel('Parent',Tab5,'Title','Anisotropic materials','Units','pixels','Position',[460 10 725 745],'FontSize',10);
        uicontrol('Parent',a.p2UI5,'Style','pushbutton','String','?','Position',[680 695 23 23],'FontSize',12,'FontWeight','bold','foregroundcolor',[1 1 1],'backgroundcolor',[.4 .4 .4],'Callback',@CallbackUI5,'Tag','69');
        uicontrol('Parent',a.p2UI5,'Style','text','String','Class','Position',[10 175+525 29 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','Material','Position',[10 175+495 39 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','Mass density (kg/m3)','Position',[10 175+465 105 13]);
        
        uicontrol('Parent',a.p2UI5,'Style','text','String','Engineering constants (GPa)','Position',[10 155+440 126 15],'FontSize',9);
        uicontrol('Parent',a.p2UI5,'Style','text','String','Real part','Position',[75 155+420 45 15]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','Imaginary part','Position',[150 155+420 70 15]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','E1 (GPa)','Position',[10 155+395 49 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','E2 (GPa)','Position',[10 155+365 49 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','E3 (GPa)','Position',[10 155+335 49 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','G12 (GPa)','Position',[10 155+305 57 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','G13 (GPa)','Position',[10 155+275 57 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','G23 (GPa)','Position',[10 155+245 57 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','v12','Position',[10 155+215 20 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','v13','Position',[10 155+185 20 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','v23','Position',[10 155+155 20 13]);
        uicontrol('Parent',a.p2UI5,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic'},'Value',a.MaterialTypeME,'TooltipString','Select a material symmetry class.','Position',[135 155+540 220 23],'Callback',@CallbackUI5,'Tag','16');
        a.Material2UI5 = uicontrol('Parent',a.p2UI5,'Style','popupmenu','String',fieldnames(a.Materials.Orthotropic),'TooltipString','Select a material to display its parameters.','Position',[135 155+510 220 23],'Callback',@CallbackUI5,'Tag','17');
        a.Density2UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',a.Material2{1}.Density,'TooltipString','Enter the mass density.','Position',[135 155+480 65 23],'Callback',@CallbackUI5,'Tag','18');        
        a.E1UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.E1)/1e9,'TooltipString','Enter the real part of Young''s modulus in the 1-direction.','Position',[75 155+390 65 23],'Callback',@CallbackUI5,'Tag','19');
        a.E2UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.E2)/1e9,'TooltipString','Enter the real part of Young''s modulus in the 2-direction.','Position',[75 155+360 65 23],'Callback',@CallbackUI5,'Tag','20');
        a.E3UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.E3)/1e9,'TooltipString','Enter the real part of Young''s modulus in the 3-direction.','Position',[75 155+330 65 23],'Callback',@CallbackUI5,'Tag','21');
        a.G12UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.G12)/1e9,'TooltipString','Enter the real part of the shear modulus in the 2-direction on the 2-3-plane.','Position',[75 155+300 65 23],'Callback',@CallbackUI5,'Tag','22');
        a.G13UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.G13)/1e9,'TooltipString','Enter the real part of the shear modulus in the 3-direction on the 2-3-plane.','Position',[75 155+270 65 23],'Callback',@CallbackUI5,'Tag','23');
        a.G23UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.G23)/1e9,'TooltipString','Enter the real part of the shear modulus in the 3-direction on the 1-3-plane.','Position',[75 155+240 65 23],'Callback',@CallbackUI5,'Tag','24');
        a.v12UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.v12),'TooltipString',['Enter the real part of Poisson''s ratio corresponding to a contraction',newline,'in the 2-direction when an extension is applied in the 1-direction.'],'Position',[75 155+210 65 23],'Callback',@CallbackUI5,'Tag','25');
        a.v13UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.v13),'TooltipString',['Enter the real part of Poisson''s ratio corresponding to a contraction',newline,'in the 3-direction when an extension is applied in the 1-direction.'],'Position',[75 155+180 65 23],'Callback',@CallbackUI5,'Tag','26');
        a.v23UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.v23),'TooltipString',['Enter the real part of Poisson''s ratio corresponding to a contraction',newline,'in the 3-direction when an extension is applied in the 2-direction.'],'Position',[75 155+150 65 23],'Callback',@CallbackUI5,'Tag','27');
        a.E1ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.E1)/1e9,'TooltipString','Enter the imaginary part of Young''s modulus in the 1-direction.','Position',[150 155+390 65 23],'Callback',@CallbackUI5,'Tag','28');
        a.E2ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.E2)/1e9,'TooltipString','Enter the imaginary part of Young''s modulus in the 2-direction.','Position',[150 155+360 65 23],'Callback',@CallbackUI5,'Tag','29');
        a.E3ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.E3)/1e9,'TooltipString','Enter the imaginary part of Young''s modulus in the 3-direction.','Position',[150 155+330 65 23],'Callback',@CallbackUI5,'Tag','30');
        a.G12ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.G12)/1e9,'TooltipString','Enter the imaginary part of the shear modulus in the 2-direction on the 2-3-plane.','Position',[150 155+300 65 23],'Callback',@CallbackUI5,'Tag','31');
        a.G13ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.G13)/1e9,'TooltipString','Enter the imaginary part of the shear modulus in the 3-direction on the 2-3-plane.','Position',[150 155+270 65 23],'Callback',@CallbackUI5,'Tag','32');
        a.G23ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.G23)/1e9,'TooltipString','Enter the imaginary part of the shear modulus in the 3-direction on the 1-3-plane.','Position',[150 155+240 65 23],'Callback',@CallbackUI5,'Tag','33');
        a.v12ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.v12),'TooltipString',['Enter the imaginary part of Poisson''s ratio corresponding to a contraction',newline,'in the 2-direction when an extension is applied in the 1-direction.'],'Position',[150 155+210 65 23],'Callback',@CallbackUI5,'Tag','34');
        a.v13ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.v13),'TooltipString',['Enter the imaginary part of Poisson''s ratio corresponding to a contraction',newline,'in the 3-direction when an extension is applied in the 1-direction.'],'Position',[150 155+180 65 23],'Callback',@CallbackUI5,'Tag','35');
        a.v23ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.v23),'TooltipString',['Enter the imaginary part of Poisson''s ratio corresponding to a contraction',newline,'in the 3-direction when an extension is applied in the 2-direction.'],'Position',[150 155+150 65 23],'Callback',@CallbackUI5,'Tag','36');

        uicontrol('Parent',a.p2UI5,'Style','text','String','Stiffness components (GPa)','Position',[230+40 155+440 157 15],'FontSize',9);
        uicontrol('Parent',a.p2UI5,'Style','text','String','Real part','Position',[230+40 155+420 45 15]);
        a.C11UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(1,1))/1e9,'TooltipString','Enter the real stiffness component in the 1-direction.','Position',[230+40 155+390 65 23],'Callback',@CallbackUI5,'Tag','37');
        a.C12UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(1,2))/1e9,'TooltipString','Enter the real stiffness component.','Position',[230+115 155+390 65 23],'Callback',@CallbackUI5,'Tag','38');
        a.C13UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(1,3))/1e9,'Position',[230+190 155+390 65 23],'Callback',@CallbackUI5,'Tag','39');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 155+390 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 155+390 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+390 65 23],'BackgroundColor',[.88 .88 .88]);
        a.C22UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(2,2))/1e9,'TooltipString','Enter the real stiffness component in the 2-direction.','Position',[230+115 155+360 65 23],'Callback',@CallbackUI5,'Tag','40');
        a.C23UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(2,3))/1e9,'TooltipString','Enter the real stiffness component.','Position',[230+190 155+360 65 23],'Callback',@CallbackUI5,'Tag','41');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 155+360 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 155+360 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+360 65 23],'BackgroundColor',[.88 .88 .88]);
        a.C33UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(3,3))/1e9,'TooltipString','Enter the real stiffness component in the 3-direction.','Position',[230+190 155+330 65 23],'Callback',@CallbackUI5,'Tag','42');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 155+330 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 155+330 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+330 65 23],'BackgroundColor',[.88 .88 .88]);        
        a.C44UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(4,4))/1e9,'TooltipString','Enter the real stiffness component.','Position',[230+265 155+300 65 23],'Callback',@CallbackUI5,'Tag','43');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 155+300 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+300 65 23],'BackgroundColor',[.88 .88 .88]);        
        a.C55UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(5,5))/1e9,'TooltipString','Enter the real stiffness component.','Position',[230+340 155+270 65 23],'Callback',@CallbackUI5,'Tag','44');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+270 65 23],'BackgroundColor',[.88 .88 .88]);         
        a.C66UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(6,6))/1e9,'TooltipString','Enter the real stiffness component.','Position',[230+415 155+240 65 23],'Callback',@CallbackUI5,'Tag','45');

        uicontrol('Parent',a.p2UI5,'Style','text','String','Imaginary part','Position',[230+40 200+180 70 15]);
        a.C11ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(1,1))/1e9,'TooltipString','Enter the imaginary stiffness component in the 1-direction.','Position',[230+40 200+150 65 23],'Callback',@CallbackUI5,'Tag','46');
        a.C12ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(1,2))/1e9,'TooltipString','Enter the imaginary stiffness component.','Position',[230+115 200+150 65 23],'Callback',@CallbackUI5,'Tag','47');
        a.C13ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(1,3))/1e9,'Position',[230+190 200+150 65 23],'Callback',@CallbackUI5,'Tag','48');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 200+150 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 200+150 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+150 65 23],'BackgroundColor',[.88 .88 .88]);
        a.C22ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(2,2))/1e9,'TooltipString','Enter the imaginary stiffness component in the 2-direction.','Position',[230+115 200+120 65 23],'Callback',@CallbackUI5,'Tag','49');
        a.C23ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(2,3))/1e9,'TooltipString','Enter the imaginary stiffness component.','Position',[230+190 200+120 65 23],'Callback',@CallbackUI5,'Tag','50');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 200+120 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 200+120 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+120 65 23],'BackgroundColor',[.88 .88 .88]);
        a.C33ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(3,3))/1e9,'TooltipString','Enter the imaginary stiffness component in the 3-direction.','Position',[230+190 200+90 65 23],'Callback',@CallbackUI5,'Tag','51');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 200+90 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 200+90 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+90 65 23],'BackgroundColor',[.88 .88 .88]);        
        a.C44ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(4,4))/1e9,'TooltipString','Enter the imaginary stiffness component.','Position',[230+265 200+60 65 23],'Callback',@CallbackUI5,'Tag','52');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 200+60 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+60 65 23],'BackgroundColor',[.88 .88 .88]);        
        a.C55ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(5,5))/1e9,'TooltipString','Enter the imaginary stiffness component.','Position',[230+340 200+30 65 23],'Callback',@CallbackUI5,'Tag','53');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+30 65 23],'BackgroundColor',[.88 .88 .88]);         
        a.C66ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(6,6))/1e9,'TooltipString','Enter the imaginary stiffness component.','Position',[230+415 200 65 23],'Callback',@CallbackUI5,'Tag','54');

        uicontrol('Parent',a.p2UI5,'Style','text','String','Bulk wave velocities (m/s)','Position',[10 120+95 143 15],'FontSize',9);
        uicontrol('Parent',a.p2UI5,'Style','text','String','1','Position',[154 120+90 4 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','2','Position',[228 120+90 6 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','3','Position',[302 120+90 6 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','Longitudinal velocity','Position',[10 120+65 99 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','Fast shear velocity','Position',[10 120+35 94 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','Slow shear velocity','Position',[10 120+5 98 13]);
        a.LongitudinalVelocity_1UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.LongitudinalVelocity_1,'TooltipString',['The phase velocity of longitudinal bulk',newline,'waves propagating in the 1-direction.'],'Position',[125 120+60 60 23],'BackgroundColor','white');
        a.FastShearVelocity_1UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.FastShearVelocity_1,'TooltipString',['The phase velocity of fast shear bulk',newline,'waves propagating in the 1-direction.'],'Position',[125 120+30 60 23],'BackgroundColor','white');
        a.SlowShearVelocity_1UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.SlowShearVelocity_1,'TooltipString',['The phase velocity of slow shear bulk',newline,'waves propagating in the 1-direction.'],'Position',[125 120+0 60 23],'BackgroundColor','white');
        a.LongitudinalVelocity_2UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.LongitudinalVelocity_2,'TooltipString',['The phase velocity of longitudinal bulk',newline,'waves propagating in the 2-direction.'],'Position',[200 120+60 60 23],'BackgroundColor','white');
        a.FastShearVelocity_2UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.FastShearVelocity_2,'TooltipString',['The phase velocity of fast shear bulk',newline,'waves propagating in the 2-direction.'],'Position',[200 120+30 60 23],'BackgroundColor','white');
        a.SlowShearVelocity_2UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.SlowShearVelocity_2,'TooltipString',['The phase velocity of slow shear bulk',newline,'waves propagating in the 2-direction.'],'Position',[200 120+0 60 23],'BackgroundColor','white');
        a.LongitudinalVelocity_3UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.LongitudinalVelocity_3,'TooltipString',['The phase velocity of longitudinal bulk',newline,'waves propagating in the 3-direction.'],'Position',[275 120+60 60 23],'BackgroundColor','white');
        a.FastShearVelocity_3UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.FastShearVelocity_3,'TooltipString',['The phase velocity of fast shear bulk',newline,'waves propagating in the 3-direction.'],'Position',[275 120+30 60 23],'BackgroundColor','white');
        a.SlowShearVelocity_3UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.SlowShearVelocity_3,'TooltipString',['The phase velocity of slow shear bulk',newline,'waves propagating in the 3-direction.'],'Position',[275 120+0 60 23],'BackgroundColor','white');

        uicontrol('Parent',a.p2UI5,'Style','text','String','New material''s name','Position',[10 65 102 13]);
        a.Name2UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',a.MaterialName4B,'TooltipString',['Enter a new material''s name and press ''Save'' or enter',newline,'an available material''s name and press ''Delete''.'],'Position',[125 60 220 23],'Callback',@CallbackUI5,'Tag','55');
        uicontrol('Parent',a.p2UI5,'Style','pushbutton','String','Save material','TooltipString','Save the material to the corresponding material list.','Position',[125 15 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','56');
        uicontrol('Parent',a.p2UI5,'Style','pushbutton','String','Delete material','TooltipString','Remove the material from the corresponding material list.','Position',[245 15 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','57');
        
        %------------------------------------------------------------------
        a.p3UI5 = uipanel('Parent',Tab5,'Title','Fluids','Units','pixels','Position',[10 10 440 205],'FontSize',10);
        uicontrol('Parent',a.p3UI5,'Style','text','String','Fluid','Position',[10 160 24 13]);
        uicontrol('Parent',a.p3UI5,'Style','text','String','Mass density (kg/m3)','Position',[10 130 105 13]);
        uicontrol('Parent',a.p3UI5,'Style','text','String','Velocity (m/s)','Position',[10 100 69 13]);
        a.Material3UI5 = uicontrol('Parent',a.p3UI5,'Style','popupmenu','String',fieldnames(a.Materials.Fluid),'TooltipString','Select a fluid to display its parameters.','Position',[135 155 220 23],'Callback',@CallbackUI5,'Tag','58');
        a.Density3UI5 = uicontrol('Parent',a.p3UI5,'Style','edit','String',a.Couplant1.Density,'TooltipString','Enter the mass density.','Position',[135 125 65 23],'Callback',@CallbackUI5,'Tag','59');
        a.VelocityUI5 = uicontrol('Parent',a.p3UI5,'Style','edit','String',a.Couplant1.Velocity,'TooltipString','Enter the phase velocity of ultrasonic waves.','Position',[135 95 65 23],'Callback',@CallbackUI5,'Tag','60');

        uicontrol('Parent',a.p3UI5,'Style','text','String','New fluid''s name','Position',[10 65 85 13]);
        a.Name3UI5 = uicontrol('Parent',a.p3UI5,'Style','edit','String',a.MaterialName4C,'TooltipString',['Enter a new fluid''s name and press ''Save'' or enter',newline,'an available fluid''s name and press ''Delete''.'],'Position',[135 60 220 23],'Callback',@CallbackUI5,'Tag','61');
        uicontrol('Parent',a.p3UI5,'Style','pushbutton','String','Save fluid','TooltipString','Save the fluid to the fluids list.','Position',[135 15 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','62');
        uicontrol('Parent',a.p3UI5,'Style','pushbutton','String','Delete fluid','TooltipString','Remove the fluid from the fluids list.','Position',[255 15 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','63');
    end
    if  1 % Tab7_advanced
        a.c1UI7 = uicontrol('Parent',Tab7,'Style','text','String','isotropic','Position',[280 740 51 15],'FontSize',10,'Backgroundcolor','white');
        a.c2UI7 = uicontrol('Parent',Tab7,'Style','text','String','anisotropic','Position',[360 740 65 15],'FontSize',10,'Backgroundcolor','white');
        
        a.p1UI7 = uipanel('Parent',Tab7,'Title','Phase velocity sweeps','Units','pixels','Position',[10 555 430 180],'FontSize',10,'Backgroundcolor','white');    
        uicontrol('Parent',a.p1UI7,'Style','text','String','Phase velocity resolution (m/s)','Position',[10 135 150 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p1UI7,'Style','text','String','Phase velocity sections (5^x)','Position',[10 105 144 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p1UI7,'Style','text','String','Lamb wave search width for negative curvature','Position',[10 75 237 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p1UI7,'Style','text','String','Lamb wave search width for positive curvature','Position',[10 45 233 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p1UI7,'Style','text','String','Shear horizontal wave search width','Position',[10 15 179 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.PhaseVelocityResolution1,'TooltipString','Enter the resolution obtained during the phase velocity sweeps.','Position',[270 130 50 23],'Callback',@CallbackUI7,'Tag','1');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.PhaseVelocityResolution2,'TooltipString','Enter the resolution obtained during the phase velocity sweeps.','Position',[350 130 50 23],'Callback',@CallbackUI7,'Tag','2');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.PhaseVelocitySections1,'TooltipString',['The dispersion curves are determined by sweeping the phase velocity at fixed frequencies for',newline,'the most part. The phase velocity sections as the power of five give the maximum number of',newline,'sections into which the phase velocity search interval is devided during the search for the',newline,'modal solution at a given frequency. A higher number increases the chance to find the',newline,'solution at the cost of processing time. Missing a solution is not necessarily critical since an',newline,'extrapolation routine replaces missing samples successfully, as long as not too many samples',newline,'are missing.'],'Position',[270 100 50 23],'Callback',@CallbackUI7,'Tag','18');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.PhaseVelocitySections2,'TooltipString',['The dispersion curves are determined by sweeping the phase velocity at fixed frequencies for',newline,'the most part. The phase velocity sections as the power of five give the maximum number of',newline,'sections into which the phase velocity search interval is devided during the search for the',newline,'modal solution at a given frequency. A higher number increases the chance to find the',newline,'solution at the cost of processing time. Missing a solution is not necessarily critical since an',newline,'extrapolation routine replaces missing samples successfully, as long as not too many samples',newline,'are missing.'],'Position',[350 100 50 23],'Callback',@CallbackUI7,'Tag','20');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.LambPhaseVelocitySweepRange11,'TooltipString',['This value determines the phase velocity search interval width for',newline,'Lamb waves. Please read the manual for more information.'],'Position',[270 70 50 23],'Callback',@CallbackUI7,'Tag','3');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.LambPhaseVelocitySweepRange21,'TooltipString',['This value determines the phase velocity search interval width for',newline,'Lamb waves. Please read the manual for more information.'],'Position',[270 40 50 23],'Callback',@CallbackUI7,'Tag','4');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.LambPhaseVelocitySweepRange12,'TooltipString',['This value determines the phase velocity search interval width for',newline,'Lamb waves. Please read the manual for more information.'],'Position',[350 70 50 23],'Callback',@CallbackUI7,'Tag','5');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.LambPhaseVelocitySweepRange22,'TooltipString',['This value determines the phase velocity search interval width for',newline,'Lamb waves. Please read the manual for more information.'],'Position',[350 40 50 23],'Callback',@CallbackUI7,'Tag','6');
        uicontrol('Parent',a.p1UI7,'Style','edit','String',a.ShearPhaseVelocitySweepRange2,'TooltipString',['This value determines the phase velocity search interval width for shear',newline,'horizontal waves. Please read the manual for more information.'],'Position',[350 10 50 23],'Callback',@CallbackUI7,'Tag','7');
        
        a.p2UI7 = uipanel('Parent',Tab7,'Title','Frequency sweeps to complete dispersion curves at high phase velocity','Units','pixels','Position',[10 425 430 120],'FontSize',10,'Backgroundcolor','white');
        uicontrol('Parent',a.p2UI7,'Style','text','String','Frequency sections (5^x)','Position',[10 75 126 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p2UI7,'Style','text','String','Phase velocity step (m/s)','Position',[10 45 124 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p2UI7,'Style','text','String','Search interval (kHz/mm)','Position',[10 15 123 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p2UI7,'Style','edit','String',a.FrequencySections1,'TooltipString',['Frequency sweeps at certain phase velocities are performed to complete the higher order',newline,'dispersion curves at high phase velocities, if necessary. The frequency sections work similarly',newline,'as the phase velocity sections, but with respect to the frequency search interval.'],'Position',[270 70 50 23],'Callback',@CallbackUI7,'Tag','19');
        uicontrol('Parent',a.p2UI7,'Style','edit','String',a.FrequencySections2,'TooltipString',['Frequency sweeps at certain phase velocities are performed to complete the higher order',newline,'dispersion curves at high phase velocities, if necessary. The frequency sections work similarly',newline,'as the phase velocity sections, but with respect to the frequency search interval.'],'Position',[350 70 50 23],'Callback',@CallbackUI7,'Tag','21');
        uicontrol('Parent',a.p2UI7,'Style','edit','String',a.PhaseVelocityStep1,'TooltipString','Enter the phase velocity step.','Position',[270 40 50 23],'Callback',@CallbackUI7,'Tag','8');
        uicontrol('Parent',a.p2UI7,'Style','edit','String',a.PhaseVelocityStep2,'TooltipString','Enter the phase velocity step.','Position',[350 40 50 23],'Callback',@CallbackUI7,'Tag','12');
        uicontrol('Parent',a.p2UI7,'Style','edit','String',a.FrequencyOffset1,'TooltipString','Enter the frequency search interval width.','Position',[270 10 50 23],'Callback',@CallbackUI7,'Tag','10');
        uicontrol('Parent',a.p2UI7,'Style','edit','String',a.FrequencyOffset2,'TooltipString','Enter the frequency search interval width.','Position',[350 10 50 23],'Callback',@CallbackUI7,'Tag','15');
  
        a.c3UI7 = uicontrol('Parent',Tab7,'Style','text','String','isotropic','Position',[628 740 51 15],'FontSize',10,'Backgroundcolor','white');
        a.c4UI7 = uicontrol('Parent',Tab7,'Style','text','String','anisotropic','Position',[720 740 65 15],'FontSize',10,'Backgroundcolor','white');
        
        a.p3UI7 = uipanel('Parent',Tab7,'Title','Fluid-loading and viscoelasticity settings','Units','pixels','Position',[450 585 370 150],'FontSize',10,'Backgroundcolor','white');    
        uicontrol('Parent',a.p3UI7,'Style','text','String','Real part search width','Position',[10 105 112 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p3UI7,'Style','text','String','Imaginary part search width','Position',[10 75 137 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p3UI7,'Style','text','String','Search area sections','Position',[10 45 106 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p3UI7,'Style','text','String','Search area extensions','Position',[10 15 118 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p3UI7,'Style','edit','String',['[',num2str(a.SearchWidthReal1(1)),' ',num2str(a.SearchWidthReal1(2)),']'],'TooltipString',['This determines the search range in the real wavenumber',newline,'dimension. For more information, read the manual.'],'Position',[170 100 70 23],'Callback',@CallbackUI7,'Tag','22');
        uicontrol('Parent',a.p3UI7,'Style','edit','String',['[',num2str(a.SearchWidthReal2(1)),' ',num2str(a.SearchWidthReal2(2)),']'],'TooltipString',['This determines the search range in the real wavenumber',newline,'dimension. For more information, read the manual.'],'Position',[270 100 70 23],'Callback',@CallbackUI7,'Tag','26');
        uicontrol('Parent',a.p3UI7,'Style','edit','String',['[',num2str(a.SearchWidthImag1(1)),' ',num2str(a.SearchWidthImag1(2)),']'],'TooltipString',['This determines the search range in the imaginary wavenumber',newline,'dimension. For more information, read the manual.'],'Position',[170 70 70 23],'Callback',@CallbackUI7,'Tag','23');
        uicontrol('Parent',a.p3UI7,'Style','edit','String',['[',num2str(a.SearchWidthImag2(1)),' ',num2str(a.SearchWidthImag2(2)),']'],'TooltipString',['This determines the search range in the imaginary wavenumber',newline,'dimension. For more information, read the manual.'],'Position',[270 70 70 23],'Callback',@CallbackUI7,'Tag','27');
        uicontrol('Parent',a.p3UI7,'Style','edit','String',a.SearchAreaSections1,'TooltipString',['This value determines at how many grid points of the search area, spanned by the real and imaginary',newline,'wavenumber parts, the characteristic function is evaluated for finding a modal solution.'],'Position',[170 40 70 23],'Callback',@CallbackUI7,'Tag','24');
        uicontrol('Parent',a.p3UI7,'Style','edit','String',a.SearchAreaSections2,'TooltipString',['This value determines at how many grid points of the search area, spanned by the real and imaginary',newline,'wavenumber parts, the characteristic function is evaluated for finding a modal solution.'],'Position',[270 40 70 23],'Callback',@CallbackUI7,'Tag','28');
        uicontrol('Parent',a.p3UI7,'Style','edit','String',a.SearchAreaExtensions1,'TooltipString',['This value determines how many times the search area, spanned by the real and imaginary',newline,'wavenumber parts, is increased in size in case a modal solution was not found previously.'],'Position',[170 10 70 23],'Callback',@CallbackUI7,'Tag','25');
        uicontrol('Parent',a.p3UI7,'Style','edit','String',a.SearchAreaExtensions2,'TooltipString',['This value determines how many times the search area, spanned by the real and imaginary',newline,'wavenumber parts, is increased in size in case a modal solution was not found previously.'],'Position',[270 10 70 23],'Callback',@CallbackUI7,'Tag','29');
        
        uicontrol('Parent',Tab7,'Style','text','String','DLR','Position',[220-53 50 40 23],'FontSize',13.25,'FontWeight','bold','Foregroundcolor',[.4 .4 .4],'Backgroundcolor','white');
        axes('Parent',Tab7,'Units','pixels','Position',[220-95 40 83 100]);
        imshow('DLR_Logo_gray.jpg')
        axes('Parent',Tab7,'Units','pixels','Position',[310 -80 1000 1000*.692]);
        imshow('Earth.jpg')
    end
end

f1.Units = 'normalized';
movegui(f1,'center') % center the GUI
f1.Visible = 'on'; % make the GUI visible

function Open_Callback(~,~)
    [Path,File] = uigetfile('*.mat');
    if  Path ~= 0
        h = msgbox('Opening project in a new DC instance...');
        delete(f1)
        load(fullfile(File,Path)); %#ok<LOAD>
        jframe = get(gcf,'javaframe');
        jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))))
        close(h)
        if  a.Multithreading == 1
            try
                p = gcp('nocreate');
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to start parallel pool')
                return
            end
            if  isempty(p)
                try
                    h = msgbox('Starting parallel pool...');
                    parpool(2,'IdleTimeout',180)
                    close(h)
                catch ME
                    close(h)
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to start parallel pool')
                    return
                end
            else
                msgbox('A parallel pool is already active.')
            end
        end
    end
end
function Save_Callback(~,~)
    uisave('a',fullfile(a.Directory,'Project'))
end
function Import_Callback(~,~)
    [FileName,Path] = uigetfile('title','Select one or multiple material lists','*.txt','MultiSelect','on');
    if  Path ~= 0
        if  ischar(FileName)
            File = [Path,FileName];
            fileID = fopen(File);
            n = numel(strfind(fgets(fileID),' '))+1;
            fclose(fileID);
            if  n == 11
                copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_Orthotropic.txt'])
                String = ['Orthotropic materials imported.',newline];
            elseif n == 7
                copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_TransverselyIsotropic.txt'])
                String = ['Transversely isotropic materials imported.',newline];
            elseif n == 5
                copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_Cubic.txt'])
                String = ['Cubic materials imported.',newline];
            elseif n == 6
                copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_Isotropic.txt'])
                String = ['Isotropic materials imported.',newline];
            elseif n == 3
                copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_Fluid.txt'])
                String = ['Fluids imported.',newline];
            end
        else
            String = '';
            for i = 1:length(FileName)
                File = [Path,FileName{i}];
                fileID = fopen(File);
                n = numel(strfind(fgets(fileID),' '))+1;
                fclose(fileID);
                if  n == 11
                    copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_Orthotropic.txt'])
                    String = append(String,'Orthotropic materials imported.',newline);
                elseif n == 7
                    copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_TransverselyIsotropic.txt'])
                    String = append(String,'Transversely isotropic materials imported.',newline);
                elseif n == 5
                    copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_Cubic.txt'])
                    String = append(String,'Cubic materials imported.',newline);
                elseif n == 6
                    copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_Isotropic.txt'])
                    String = append(String,'Isotropic materials imported.',newline);
                elseif n == 3
                    copyfile(File,[a.MaterialListDirectory,filesep,'MaterialList_Fluid.txt'])
                    String = append(String,'Fluids imported.',newline);
                end
            end
        end
        msgbox([String,newline,'Please restart the DC after importing all material lists to have them available.'],'Info')
    end
end
function Export_Callback(~,~)
    try
        Path = uigetdir('title','Select folder where to save the material lists');
        if  Path ~= 0
            copyfile(a.MaterialListDirectory,Path)
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export material lists')
        return
    end
end
function Enable_Callback(~,~)
    try
        p = gcp('nocreate');
    catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to start parallel pool')
        return
    end
    if  isempty(p)
        try
            h = msgbox('Starting parallel pool...');
            parpool(2,'IdleTimeout',180)
            a.Multithreading = 1;
            delete(a.a1UI1)
            delete(a.a1UI2)
            delete(a.a1UI4)
            a.a1UI1 = axes('Parent',Tab1,'Units','pixels','Position',[1130 705 34 34]);
            a.a1UI2 = axes('Parent',Tab2,'Units','pixels','Position',[1130 705 34 34]);
            a.a1UI4 = axes('Parent',Tab4,'Units','pixels','Position',[1130 705 34 34]);
            imshow('Multithreading21.png','Parent',a.a1UI1)
            imshow('Multithreading21.png','Parent',a.a1UI2)
            imshow('Multithreading21.png','Parent',a.a1UI4)
            close(h)
        catch ME
            close(h)
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to start parallel pool')
            return
        end
    else
        msgbox('A parallel pool is already active.')
        a.Multithreading = 1;
        delete(a.a1UI1)
        delete(a.a1UI2)
        delete(a.a1UI4)
        a.a1UI1 = axes('Parent',Tab1,'Units','pixels','Position',[1130 705 34 34]);
        a.a1UI2 = axes('Parent',Tab2,'Units','pixels','Position',[1130 705 34 34]);
        a.a1UI4 = axes('Parent',Tab4,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading21.png','Parent',a.a1UI1)
        imshow('Multithreading21.png','Parent',a.a1UI2)
        imshow('Multithreading21.png','Parent',a.a1UI4)
    end
end
function Disable_Callback(~,~)
    try
        p = gcp('nocreate');
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to shut down parallel pool')
        return
    end
    if  isempty(p)
        msgbox('There is no active parallel pool to shut down.')
        a.Multithreading = 0;
        delete(a.a1UI1)
        delete(a.a1UI2)
        delete(a.a1UI4)
        a.a1UI1 = axes('Parent',Tab1,'Units','pixels','Position',[1130 705 34 34]);
        a.a1UI2 = axes('Parent',Tab2,'Units','pixels','Position',[1130 705 34 34]);
        a.a1UI4 = axes('Parent',Tab4,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading20.png','Parent',a.a1UI1)
        imshow('Multithreading20.png','Parent',a.a1UI2)
        imshow('Multithreading20.png','Parent',a.a1UI4)
    else
        h = msgbox('Shutting down parallel pool...');
        delete(p)
        a.Multithreading = 0;
        delete(a.a1UI1)
        delete(a.a1UI2)
        delete(a.a1UI4)
        a.a1UI1 = axes('Parent',Tab1,'Units','pixels','Position',[1130 705 34 34]);
        a.a1UI2 = axes('Parent',Tab2,'Units','pixels','Position',[1130 705 34 34]);
        a.a1UI4 = axes('Parent',Tab4,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading20.png','Parent',a.a1UI1)
        imshow('Multithreading20.png','Parent',a.a1UI2)
        imshow('Multithreading20.png','Parent',a.a1UI4)
        close(h)
    end
end
function Help_Callback(~,~)
    winopen 'DispersionCalculator_Manual.pdf'   
end
function About_Callback(~,~)
    f3 = figure('NumberTitle','off','Name','About','Visible','off','MenuBar','none','Position',[0 0 520 680],'color','w');

    jframe = get(gcf,'javaframe');
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
  
    f3.Units = 'normalized';
    movegui(f3,'center')
    f3.Visible = 'on';

    uicontrol('Parent',f3,'Style','text','String','Dispersion Calculator','Position',[130 340+272 182 23],'FontSize',14,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','2','Position',[310 340+272 10 23],'FontSize',10,'FontWeight','bold','Foregroundcolor',[1 .788 .055],'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String',['Version ',Version,' - compiled ',Date],'Position',[130 340+238 ReleaseStringLength 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String',['Copyright ',char(169),' 2018-20',Date(end-1:end),' DLR'],'Position',[130 340+220 186 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','Dispersion Calculator was created by Armin Huber','Position',[130 340+186 330 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','Center for Lightweight Production Technology (ZLP)','Position',[130 340+168 337 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','Institute of Structures and Design','Position',[130 340+150 218 18],'FontSize',11,'Backgroundcolor','white');

    uicontrol('Parent',f3,'Style','text','String','Deutsches Zentrum','Position',[130 340+101 173 23],'FontSize',13.25,'FontWeight','bold','Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String',['f',char(252),'r Luft- und Raumfahrt'],'Position',[130 340+79 205 23],'FontSize',13.25,'FontWeight','bold','Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','DLR','Position',[130-53 340+79 40 23],'FontSize',13.25,'FontWeight','bold','Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','German Aerospace Center','Position',[130 340+56.5 215 23],'FontSize',13.25,'Backgroundcolor','white');

    uicontrol('Parent',f3,'Style','text','String','For more information contact Armin Huber at','Position',[130 340+18 291 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','armin.huber@dlr.de','Position',[130 340 131 18],'FontSize',11,'Backgroundcolor','white');

    axes('Parent',f3,'Units','pixels','Position',[130-80 340+205 60 60]);
    imshow('DC_Logo60.png')
    axes('Parent',f3,'Units','pixels','Position',[130-95 340+69 83 100]);
    imshow('DLR_Logo_black.jpg')
    axes('Parent',f3,'Units','pixels','Position',[176 150 166 166]);
    imshow('DC_Logo_Earth.png') 
    
    uicontrol('Parent',f3,'Style','text','String','Thanks to','Position',[130 110 66 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','Michael Lowe, Michel Castaings,','Position',[130 92 213 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','Stanislav Rokhlin, Victor Giurgiutiu,','Position',[130 74 227 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','Marc Deschamps, Eric Ducasse,','Position',[130 56 219 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','Markus Sause.','Position',[130 38 99 18],'FontSize',11,'Backgroundcolor','white');
end
function CloseRequest(~,~)
    Selection = questdlg('Do you really want to close the Dispersion Calculator?','Close request','Yes','No','Yes'); 
    switch Selection 
        case 'Yes'
            delete(gcf)
    end
end
function SliderUI1_Callback(source,~,~)
    a.p1UI1.Position(2) = 610+(source.Value-1)*115;
    a.p2UI1.Position(2) = 490+(source.Value-1)*115;
    a.p3UI1.Position(2) = 280+(source.Value-1)*115;
    a.p4UI1.Position(2) = 690+(source.Value-1)*115;
    a.OutputWindow1aUI1.Position(2) = 50+(source.Value-1)*115;
    a.OutputWindow1bUI1.Position(2) = 50+(source.Value-1)*115;
    a.OutputWindow2aUI1.Position(2) = 10+(source.Value-1)*115;
    a.OutputWindow2bUI1.Position(2) = 10+(source.Value-1)*115;
    a.c1UI1.Position(2) = 702+(source.Value-1)*115;
    a.c2UI1.Position(2) = 702+(source.Value-1)*115;
    a.p5UI1.Position(2) = 465+(source.Value-1)*115;
    a.p6UI1.Position(2) = 465+(source.Value-1)*115;
    a.p7UI1.Position(2) = 195+(source.Value-1)*115;
    a.p8UI1.Position(2) = 265+(source.Value-1)*115;
    a.p9UI1.Position(2) = 10+(source.Value-1)*115;
    a.p10UI1.Position(2) = 90+(source.Value-1)*115;
    a.p11UI1.Position(2) = 400+(source.Value-1)*115;
    a.p12UI1.Position(2) = 280+(source.Value-1)*115;
    a.p13UI1.Position(2) = 100+(source.Value-1)*115;
    a.c3UI1.Position(2) = 45+(source.Value-1)*115;
    a.a1UI1.Position(2) = 705+(source.Value-1)*115;
    
    a.p1UI2.Position(2) = 525+(source.Value-1)*115;
    a.p2UI2.Position(2) = 375+(source.Value-1)*115;
    a.p3UI2.Position(2) = 165+(source.Value-1)*115;
    a.p7UI2.Position(2) = 690+(source.Value-1)*115;
    a.OutputWindow1aUI2.Position(2) = 50+(source.Value-1)*115;
    a.OutputWindow1bUI2.Position(2) = 50+(source.Value-1)*115;
    a.OutputWindow2aUI2.Position(2) = 10+(source.Value-1)*115;
    a.OutputWindow2bUI2.Position(2) = 10+(source.Value-1)*115;
    a.c1UI2.Position(2) = 702+(source.Value-1)*115;
    a.c2UI2.Position(2) = 702+(source.Value-1)*115;
    a.p8UI2.Position(2) = 465+(source.Value-1)*115;
    a.p9UI2.Position(2) = 465+(source.Value-1)*115;
    a.p10UI2.Position(2) = 195+(source.Value-1)*115;
    a.p11UI2.Position(2) = 265+(source.Value-1)*115;
    a.p12UI2.Position(2) = 10+(source.Value-1)*115;
    a.p13UI2.Position(2) = 90+(source.Value-1)*115;
    a.p14UI2.Position(2) = 400+(source.Value-1)*115;
    a.p15UI2.Position(2) = 280+(source.Value-1)*115;
    a.p16UI2.Position(2) = 100+(source.Value-1)*115;
    a.c3UI2.Position(2) = 45+(source.Value-1)*115;
    a.a1UI2.Position(2) = 705+(source.Value-1)*115;

    a.DataUI3.Position(2) = 740+(source.Value-1)*115;
    a.p1UI3.Position(2) = 420+(source.Value-1)*115;
    a.p2UI3.Position(2) = 185+(source.Value-1)*115;
    a.p3UI3.Position(2) = 10+(source.Value-1)*115;
    a.OutputWindowUI3.Position(2) = 10+(source.Value-1)*115;
    a.c1UI3.Position(2) = 165+(source.Value-1)*115;
    a.c2UI3.Position(2) = 115+(source.Value-1)*115;
    a.XAxisUI3.Position(2) = 135+(source.Value-1)*115;
    a.YAxisUI3.Position(2) = 85+(source.Value-1)*115;
    a.PlotUI3.Position(2) = 40+(source.Value-1)*115;
    a.p4UI3.Position(2) = 10+(source.Value-1)*115;
    a.p5UI3.Position(2) = 50+(source.Value-1)*115;
    a.c3UI3.Position(2) = 149+(source.Value-1)*115;
    if  a.MultiMode3 == 0
        a.h2.Position(2) = 515+(source.Value-1)*115;
    elseif a.MultiMode3 == 1
        a.h2.Position(2) = 575+(source.Value-1)*115;
    end
    a.h3.Position(2) = 230+(source.Value-1)*115;
    a.h4.Position(2) = 230+(source.Value-1)*115;
    a.h5.Position(2) = 455+(source.Value-1)*115;
    
    a.p1UI4.Position(2) = 555+(source.Value-1)*115;
    a.p2UI4.Position(2) = 375+(source.Value-1)*115;
    a.p3UI4.Position(2) = 255+(source.Value-1)*115;
    a.c1UI4.Position(2) = 702+(source.Value-1)*115;
    a.c2UI4.Position(2) = 702+(source.Value-1)*115;
    a.p4UI4.Position(2) = 520+(source.Value-1)*115;
    a.p5UI4.Position(2) = 335+(source.Value-1)*115;
    a.p6UI4.Position(2) = 520+(source.Value-1)*115;
    a.p7UI4.Position(2) = 534+(source.Value-1)*115;
    a.c3UI4.Position(2) = 693+(source.Value-1)*115;
    a.a1UI4.Position(2) = 705+(source.Value-1)*115;
    
    a.OutputWindowUI8.Position(2) = 10+(source.Value-1)*115;
    a.p1UI8.Position(2) = 545+(source.Value-1)*115;
    a.p2UI8.Position(2) = 555+(source.Value-1)*115;
    a.p3UI8.Position(2) = 555+(source.Value-1)*115;
    a.p4UI8.Position(2) = 275+(source.Value-1)*115;
    a.p5UI8.Position(2) = 215+(source.Value-1)*115;
    a.p6UI8.Position(2) = 375+(source.Value-1)*115;
    a.p7UI8.Position(2) = 225+(source.Value-1)*115;
    a.p8UI8.Position(2) = 70+(source.Value-1)*115;
    a.a1UI8.Position(2) = 390+(source.Value-1)*115;
    a.a2UI8.Position(2) = 5+(source.Value-1)*115;
    a.c1UI8.Position(2) = 500+(source.Value-1)*115;
    
    a.p1UI6.Position(2) = (source.Value-1)*115;
    
    a.p1UI5.Position(2) = 215+(source.Value-1)*115;
    a.p2UI5.Position(2) = 5+(source.Value-1)*115;
    a.p3UI5.Position(2) = 5+(source.Value-1)*115;
    
    a.c1UI7.Position(2) = 740+(source.Value-1)*115;
    a.c2UI7.Position(2) = 740+(source.Value-1)*115;
    a.p1UI7.Position(2) = 555+(source.Value-1)*115;
    a.p2UI7.Position(2) = 425+(source.Value-1)*115;
    a.c3UI7.Position(2) = 740+(source.Value-1)*115;
    a.c4UI7.Position(2) = 740+(source.Value-1)*115;
    a.p3UI7.Position(2) = 585+(source.Value-1)*115;
end
function SpecimenSettingsUI2_Callback(~,~) % Tab2_anisotropic
    f2 = figure('NumberTitle','off','Name','Specimen_Anisotropic','Visible','off','MenuBar','none','Position',[0 0 910 400]);
    jframe = get(gcf,'javaframe');
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
    f2.Units = 'normalized';
    movegui(f2,'center')
    f2.Visible = 'on';
        
    uicontrol('Parent',f2,'Style','pushbutton','String','Open','TooltipString','Open an existing specimen definition.','Position',[20 350 95 33],'FontSize',10,'Callback',@OpenUI2_Callback);
    uicontrol('Parent',f2,'Style','pushbutton','String','Save','TooltipString','Save the current specimen definition.','Position',[20+130 350 95 33],'FontSize',10,'Callback',@SaveUI2_Callback);
    uicontrol('Parent',f2,'Style','pushbutton','String','Reset','TooltipString','Reset specimen definition.','Position',[20+260 350 95 33],'FontSize',10,'Callback',@ResetUI2_Callback);    
    
    uicontrol('Parent',f2,'Style','text','String','Upper fluid','Position',[20 70+245 54 13])
    uicontrol('Parent',f2,'Style','text','String','Lower fluid','Position',[20 70+215 57 13])
    uicontrol('Parent',f2,'Style','text','String','Hybrid','Position',[20 70+185 33 13])
    uicontrol('Parent',f2,'Style','text','String','Class','Position',[20 70+155 29 13])
    uicontrol('Parent',f2,'Style','text','String','Material','Position',[20 70+125 39 13]);
    uicontrol('Parent',f2,'Style','text','String','Uniform layer thickness','Position',[20 70+95 115 13]);
    uicontrol('Parent',f2,'Style','text','String','Total thickness (mm)','Position',[20 70+65 101 13]);
    uicontrol('Parent',f2,'Style','text','String','Unit cell repetitions','Position',[20 70+35 92 13]);
    uicontrol('Parent',f2,'Style','text','String','Symmetric system','Position',[20 70+5 90 13]);
       
    a.ToggleUpperFluidUI2 = uicontrol('Parent',f2,'Style','checkbox','Value',a.ToggleUpperFluid2,'TooltipString','Toggle upper fluid. If unchecked, the upper half-space will be vacuum.','Position',[150-60  70+240 50 23],'Callback',@ToggleUpperFluidUI2_Callback);
    a.ToggleLowerFluidUI2 = uicontrol('Parent',f2,'Style','checkbox','Value',a.ToggleLowerFluid2,'TooltipString','Toggle lower fluid. If unchecked, the lower half-space will be vacuum.','Position',[150-60 70+210 50 23],'Callback',@ToggleLowerFluidUI2_Callback);
    a.SelectUpperFluidUI2 = uicontrol('Parent',f2,'Style','popupmenu','Value',a.SelectUpperFluidUI2Value,'String',fieldnames(a.Materials.Fluid),'TooltipString','Select the fluid filling the half-space above the laminate.','Position',[150-30 70+240 140 23],'Enable',a.SelectUpperFluidUI2Enable,'Callback',@SelectUpperFluidUI2_Callback);
    a.SelectLowerFluidUI2 = uicontrol('Parent',f2,'Style','popupmenu','Value',a.SelectLowerFluidUI2Value,'String',fieldnames(a.Materials.Fluid),'TooltipString','Select the fluid filling the half-space below the laminate.','Position',[150-30 70+210 140 23],'Enable',a.SelectLowerFluidUI2Enable,'Callback',@SelectLowerFluidUI2_Callback);
    a.HybridUI2 = uicontrol('Parent',f2,'Style','checkbox','Value',a.Hybrid,'TooltipString','The layup may contain different materials.','Position',[150-60 70+180 50 23],'Callback',@HybridUI2_Callback);
    a.MaterialTypeUI2 = uicontrol('Parent',f2,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic','Isotropic'},'Value',a.MaterialType2,'TooltipString','Select a material symmetry class.','Enable',a.MaterialTypeUI2Enable,'Position',[150-60 70+150 170 23],'Callback',@MaterialTypeUI2_Callback);
    if  a.MaterialType2 == 1
        a.MaterialUI2 = uicontrol('Parent',f2,'Style','popupmenu','Value',a.MaterialUI2Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Orthotropic),'Enable',a.MaterialUI2Enable,'Position',[150-60 70+120 170 23],'Callback',@MaterialUI2_Callback);
    elseif a.MaterialType2 == 2
        a.MaterialUI2 = uicontrol('Parent',f2,'Style','popupmenu','Value',a.MaterialUI2Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.TransverselyIsotropic),'Enable',a.MaterialUI2Enable,'Position',[150-60 70+120 170 23],'Callback',@MaterialUI2_Callback);
    elseif a.MaterialType2 == 3
        a.MaterialUI2 = uicontrol('Parent',f2,'Style','popupmenu','Value',a.MaterialUI2Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Cubic),'Enable',a.MaterialUI2Enable,'Position',[150-60 70+120 170 23],'Callback',@MaterialUI2_Callback);
    elseif a.MaterialType2 == 4
        a.MaterialUI2 = uicontrol('Parent',f2,'Style','popupmenu','Value',a.MaterialUI2Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Isotropic),'Enable',a.MaterialUI2Enable,'Position',[150-60 70+120 170 23],'Callback',@MaterialUI2_Callback);
    end
    a.UniformLayerThicknessUI2 = uicontrol('Parent',f2,'Style','checkbox','Value',a.UniformLayerThickness,'TooltipString',['Check this if the layup has uniform layer thicknesses. Then, the Dispersion Calculator',newline,'deduces the layer thicknesses from the overall plate thickness as defined below, and',newline,'you do not need to enter the individual layer thicknesses into the table to the right.'],'Position',[150 70+90 50 23],'Callback',@UniformLayerThicknessUI2_Callback);
    a.PlateThicknessUI2 = uicontrol('Parent',f2,'Style','edit','String',a.PlateThickness,'TooltipString',['Enter the overall plate thickness. This is available',newline,'only if you have checked ''Uniform layer thickness''.'],'Enable',a.PlateThicknessUI2Enable,'Position',[150 70+60 50 23],'Callback',@PlateThicknessUI2_Callback);
    a.SuperLayersUI2 = uicontrol('Parent',f2,'Style','edit','String',a.SuperLayers,'TooltipString',['Enter the number of repetitions of the unit cell',newline,'contained in the laminate. Use only integers.'],'Position',[150 70+30 50 23],'Callback',@SuperLayersUI2_Callback);
    a.SymmetricSystemUI2 = uicontrol('Parent',f2,'Style','checkbox','Value',a.SymmetricSystem,'TooltipString',['Check this if the laminate is symmetric about its middle plane.',newline,'The Dispersion Calculator extends the laminate accordingly.'],'Position',[150 70+0 50 23],'Callback',@SymmetricSystemUI2_Callback);    

    uicontrol('Parent',f2,'Style','pushbutton','String','OK','TooltipString','Accept the current specimen definition and continue.','Position',[20 20 95 33],'FontSize',10,'Callback',@OKUI2_Callback);
    uicontrol('Parent',f2,'Style','pushbutton','String','Cancel','TooltipString','Discard the current specimen definition.','Position',[20+130 20 95 33],'FontSize',10,'Callback',@CancelUI2_Callback);
    
    uicontrol('Parent',f2,'Style','text','String','Unit cell','Position',[280 20+297 46 13],'FontSize',9);
    a.TablesUI2 = uitable(f2,'ColumnName',{['Phi (',char(176),')'],'d (mm)','Orthotropic','Trans. iso.','Cubic','Isotropic','Delete'},'ColumnWidth',{50 50 100 100 100 100 50},'ColumnFormat',({[] [] fieldnames(a.Materials.Orthotropic)' fieldnames(a.Materials.TransverselyIsotropic)' fieldnames(a.Materials.Cubic)' fieldnames(a.Materials.Isotropic)' 'logical'}),'TooltipString',['Enter the fiber orientation of every layer into the left column, and if you don''t have checked',newline,'''Uniform layer thickness'', enter also the corresponding layer thicknesses into the second column.'],'Data',a.UnitCell,'ColumnEditable',a.TablesUI2ColumnEditable,'RowStriping','off','Position',[280 20 613 291],'CellEditCallback',@TablesUI2_Callback);

    function OpenUI2_Callback(~,~)
        [Path,File] = uigetfile('*.mat');
        if  Path ~= 0
            D4 = load(fullfile(File,Path));
            a.ToggleUpperFluid2 = D4.D5.ToggleUpperFluid;
            a.ToggleLowerFluid2 = D4.D5.ToggleLowerFluid;
            a.UpperFluid2 = D4.D5.UpperFluid;
            a.LowerFluid2 = D4.D5.LowerFluid;
            a.Hybrid = D4.D5.Hybrid;
            a.MaterialType2 = D4.D5.MaterialType;
            a.Material2 = D4.D5.Material;
            a.MaterialNames = D4.D5.MaterialNames;
            a.MaterialClasses = D4.D5.MaterialClasses;
            a.UniformLayerThickness = D4.D5.UniformLayerThickness;
            a.PlateThickness = D4.D5.PlateThickness;
            a.SuperLayers = D4.D5.SuperLayers;
            a.SymmetricSystem = D4.D5.SymmetricSystem;
            a.UnitCell = D4.D5.UnitCell;
            a.LayerOrientations = D4.D5.LayerOrientations;
            a.LayerThicknesses = D4.D5.LayerThicknesses;
            a.ToggleUpperFluidUI2.Value = a.ToggleUpperFluid2;
            a.ToggleLowerFluidUI2.Value = a.ToggleLowerFluid2;
            a.SelectUpperFluidUI2.Value = find(strcmp(fieldnames(a.Materials.Fluid),a.UpperFluid2.Name));
            a.SelectLowerFluidUI2.Value = find(strcmp(fieldnames(a.Materials.Fluid),a.LowerFluid2.Name));
            a.SelectUpperFluidUI2Value = a.SelectUpperFluidUI2.Value;
            a.SelectLowerFluidUI2Value = a.SelectLowerFluidUI2.Value;            
            if  a.ToggleUpperFluid2 == 1
                a.SelectUpperFluidUI2.Enable = 'on';
                a.SelectUpperFluidUI2Enable = 'on';
            else
                a.SelectUpperFluidUI2.Enable = 'off';
                a.SelectUpperFluidUI2Enable = 'off';
            end
            if  a.ToggleLowerFluid2 == 1
                a.SelectLowerFluidUI2.Enable = 'on';
                a.SelectLowerFluidUI2Enable = 'on';
            else
                a.SelectLowerFluidUI2.Enable = 'off';
                a.SelectLowerFluidUI2Enable = 'off';
            end
            a.HybridUI2.Value = a.Hybrid;
            a.MaterialTypeUI2.Value = a.MaterialType2;
            if  a.Hybrid == 0
                if  a.MaterialType2 == 1
                    a.MaterialUI2.String = fieldnames(a.Materials.Orthotropic);
                    a.MaterialUI2.Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material2{1}.Name));
                elseif a.MaterialType2 == 2
                    a.MaterialUI2.String = fieldnames(a.Materials.TransverselyIsotropic);
                    a.MaterialUI2.Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material2{1}.Name));
                elseif a.MaterialType2 == 3
                    a.MaterialUI2.String = fieldnames(a.Materials.Cubic);
                    a.MaterialUI2.Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material2{1}.Name));
                elseif a.MaterialType2 == 4
                    a.MaterialUI2.String = fieldnames(a.Materials.Isotropic);
                    a.MaterialUI2.Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material2{1}.Name));
                end
            end
            a.MaterialUI2Value = a.MaterialUI2.Value;
            if  a.Hybrid == 1
                a.MaterialTypeUI2.Enable = 'off';
                a.MaterialTypeUI2Enable = 'off';
                a.MaterialUI2.Enable = 'off';
                a.MaterialUI2Enable = 'off';
            else
                a.MaterialTypeUI2.Enable = 'on';
                a.MaterialTypeUI2Enable = 'on';
                a.MaterialUI2.Enable = 'on';
                a.MaterialUI2Enable = 'on';
            end
            a.UniformLayerThicknessUI2.Value = a.UniformLayerThickness;
            if  a.UniformLayerThickness == 1
                a.PlateThicknessUI2.Enable = 'on';
                a.PlateThicknessUI2Enable = 'on';
                if  a.Hybrid == 1
                    a.TablesUI2.ColumnEditable = [true false true true true true true];
                    a.TablesUI2ColumnEditable = [true false true true true true true];
                else
                    a.TablesUI2.ColumnEditable = [true false false false false false true];
                    a.TablesUI2ColumnEditable = [true false false false false false true];
                end
            else
                a.PlateThicknessUI2.Enable = 'off';
                a.PlateThicknessUI2Enable = 'off';
                if  a.Hybrid == 1
                    a.TablesUI2.ColumnEditable = [true true true true true true true];
                    a.TablesUI2ColumnEditable = [true true true true true true true];
                else
                    a.TablesUI2.ColumnEditable = [true true false false false false true];
                    a.TablesUI2ColumnEditable = [true true false false false false true];
                end
            end
            a.PlateThicknessUI2.String = a.PlateThickness;
            a.SuperLayersUI2.String = a.SuperLayers;
            a.SymmetricSystemUI2.Value = a.SymmetricSystem;
            a.TablesUI2.Data = a.UnitCell;
        end
    end
    function SaveUI2_Callback(~,~)
        D5.ToggleUpperFluid = a.ToggleUpperFluid2;
        D5.ToggleLowerFluid = a.ToggleLowerFluid2;
        D5.UpperFluid = a.UpperFluid2;
        D5.LowerFluid = a.LowerFluid2;
        D5.Hybrid = a.Hybrid;
        D5.MaterialType = a.MaterialType2;
        D5.Material = a.Material2;
        D5.MaterialNames = a.MaterialNames;
        D5.MaterialClasses = a.MaterialClasses;
        D5.UniformLayerThickness = a.UniformLayerThickness;
        D5.PlateThickness = a.PlateThickness;
        D5.SuperLayers = a.SuperLayers;
        D5.SymmetricSystem = a.SymmetricSystem;
        D5.UnitCell = a.TablesUI2.Data;
        D5.LayerOrientations = a.LayerOrientations;
        D5.LayerThicknesses = a.LayerThicknesses;
        uisave('D5',fullfile(a.Directory,'Specimen'))
    end
    function ResetUI2_Callback(~,~)
        a.ToggleUpperFluid2 = 0;
        a.ToggleLowerFluid2 = 0;
        a.SelectUpperFluidUI2Value = 1;
        a.SelectLowerFluidUI2Value = 1;
        a.SelectUpperFluidUI2Enable = 'off';
        a.SelectLowerFluidUI2Enable = 'off';
        E = fieldnames(a.Materials.Fluid);
        a.UpperFluid2 = getfield(a.Materials.Fluid,E{1});
        a.LowerFluid2 = getfield(a.Materials.Fluid,E{1});
        a.Hybrid = 0;
        a.MaterialType2 = 1;
        a.MaterialTypeUI2Enable = 'on';
        E = fieldnames(a.Materials.Orthotropic);
        a.Material2{1} = getfield(a.Materials.Orthotropic,E{1});
        a.MaterialNames{1} = a.Material2{1}.Name;
        a.MaterialClasses{1} = a.Material2{1}.Class;
        if  length(a.Material2) > 1
            a.Material2(2:end) = [];
            a.MaterialNames(2:end) = [];
            a.MaterialClasses(2:end) = [];
        end
        a.MaterialUI2Value = 1;
        a.MaterialUI2Enable = 'on';
        a.UniformLayerThickness = 1;
        a.PlateThicknessUI2Enable = 'on';
        a.TablesUI2.ColumnEditable = [true false false false false false true];
        a.TablesUI2ColumnEditable = [true false false false false false true];
        a.PlateThickness = 1;
        a.SuperLayers = 1;
        a.SymmetricSystem = 0;
        a.UnitCell = cell(400,6);
        a.UnitCell{1} = '0';
        a.UnitCell{1,2} = '1';
        a.UnitCell{1,3} = E{1};
        a.LayerOrientations = 0;
        a.LayerThicknesses = 1;
        
        a.ToggleUpperFluidUI2.Value = a.ToggleUpperFluid2;
        a.ToggleLowerFluidUI2.Value = a.ToggleLowerFluid2;
        a.SelectUpperFluidUI2.Value = a.SelectUpperFluidUI2Value;
        a.SelectLowerFluidUI2.Value = a.SelectLowerFluidUI2Value;
        a.SelectUpperFluidUI2.Enable = a.SelectUpperFluidUI2Enable;
        a.SelectLowerFluidUI2.Enable = a.SelectLowerFluidUI2Enable;
        a.HybridUI2.Value = a.Hybrid;
        a.MaterialTypeUI2.Value = a.MaterialType2;
        a.MaterialTypeUI2.Enable = a.MaterialTypeUI2Enable;
        a.MaterialUI2.Value = a.MaterialUI2Value;
        a.MaterialUI2.String = E;
        a.MaterialUI2.Enable = a.MaterialUI2Enable;
        a.UniformLayerThicknessUI2.Value = a.UniformLayerThickness;
        a.PlateThicknessUI2.String = a.PlateThickness;
        a.SuperLayersUI2.String = a.SuperLayers;
        a.SymmetricSystemUI2.Value = a.SymmetricSystem;
        a.TablesUI2.Data = a.UnitCell;
    end
    function ToggleUpperFluidUI2_Callback(source,~,~)
        a.ToggleUpperFluid2 = source.Value;
        if  source.Value == 1
            a.SelectUpperFluidUI2.Enable = 'on';
            a.SelectUpperFluidUI2Enable = 'on';
        else
            a.SelectUpperFluidUI2.Enable = 'off';
            a.SelectUpperFluidUI2Enable = 'off';
        end
    end
    function ToggleLowerFluidUI2_Callback(source,~,~)
        a.ToggleLowerFluid2 = source.Value;
        if  source.Value == 1
            a.SelectLowerFluidUI2.Enable = 'on';
            a.SelectLowerFluidUI2Enable = 'on';
        else
            a.SelectLowerFluidUI2.Enable = 'off';
            a.SelectLowerFluidUI2Enable = 'off';
        end
    end
    function SelectUpperFluidUI2_Callback(source,~,~)
        a.UpperFluid2 = getfield(a.Materials.Fluid,cell2mat(source.String(source.Value)));
        a.SelectUpperFluidUI2Value = find(strcmp(fieldnames(a.Materials.Fluid),a.UpperFluid2.Name));
    end
    function SelectLowerFluidUI2_Callback(source,~,~)
        a.LowerFluid2 = getfield(a.Materials.Fluid,cell2mat(source.String(source.Value)));
        a.SelectLowerFluidUI2Value = find(strcmp(fieldnames(a.Materials.Fluid),a.LowerFluid2.Name));
    end
    function HybridUI2_Callback(source,~,~)
        a.Hybrid = source.Value;
        if  source.Value == 1
            a.MaterialTypeUI2.Enable = 'off';
            a.MaterialTypeUI2Enable = 'off';
            a.MaterialUI2.Enable = 'off';
            a.MaterialUI2Enable = 'off';
            if  a.UniformLayerThickness == 1
                a.TablesUI2.ColumnEditable = [true false true true true true true];
                a.TablesUI2ColumnEditable = [true false true true true true true];
            else
                a.TablesUI2.ColumnEditable = [true true true true true true true];
                a.TablesUI2ColumnEditable = [true true true true true true true];
            end
        else
            a.MaterialTypeUI2.Enable = 'on';
            a.MaterialTypeUI2Enable = 'on';
            a.MaterialUI2.Enable = 'on';
            a.MaterialUI2Enable = 'on';
            if  a.UniformLayerThickness == 1
                a.TablesUI2.ColumnEditable = [true false false false false false true];
                a.TablesUI2ColumnEditable = [true false false false false false true];
            else
                a.TablesUI2.ColumnEditable = [true true false false false false true];
                a.TablesUI2ColumnEditable = [true true false false false false true];
            end
            
            if  a.MaterialType2 == 1
                [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(a.MaterialUI2.String(a.MaterialUI2.Value))));
                a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material2{1}.Name));
                a.TablesUI2.Data(1:length(a.LayerOrientations),4:6) = {''};
            elseif a.MaterialType2 == 2
                [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(a.MaterialUI2.String(a.MaterialUI2.Value))));
                a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material2{1}.Name));
                a.TablesUI2.Data(1:length(a.LayerOrientations),3) = {''};
                a.TablesUI2.Data(1:length(a.LayerOrientations),5:6) = {''};
            elseif a.MaterialType2 == 3
                [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Cubic,cell2mat(a.MaterialUI2.String(a.MaterialUI2.Value))));
                a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material2{1}.Name));
                a.TablesUI2.Data(1:length(a.LayerOrientations),3:4) = {''};
                a.TablesUI2.Data(1:length(a.LayerOrientations),6) = {''};
            elseif a.MaterialType2 == 4
                [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Isotropic,cell2mat(a.MaterialUI2.String(a.MaterialUI2.Value))));
                a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material2{1}.Name));
                a.TablesUI2.Data(1:length(a.LayerOrientations),3:5) = {''};
            end
            a.TablesUI2.Data(1:length(a.LayerOrientations),a.MaterialType2+2) = {a.Material2{1}.Name};
        end
    end
    function MaterialTypeUI2_Callback(source,~,~)
        a.MaterialType2 = source.Value;
        a.MaterialUI2.Value = 1;
        switch source.Value
        case 1
            a.MaterialUI2.String = fieldnames(a.Materials.Orthotropic);
            E = fieldnames(a.Materials.Orthotropic);
            [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Orthotropic,E{1}));
            a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material2{1}.Name));
            a.TablesUI2.Data(1:length(a.LayerOrientations),4:6) = {''};
        case 2
            a.MaterialUI2.String = fieldnames(a.Materials.TransverselyIsotropic);
            E = fieldnames(a.Materials.TransverselyIsotropic);
            [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.TransverselyIsotropic,E{1}));
            a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material2{1}.Name));
            a.TablesUI2.Data(1:length(a.LayerOrientations),3) = {''};
            a.TablesUI2.Data(1:length(a.LayerOrientations),5:6) = {''};
        case 3
            a.MaterialUI2.String = fieldnames(a.Materials.Cubic);
            E = fieldnames(a.Materials.Cubic);
            [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Cubic,E{1}));
            a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material2{1}.Name));
            a.TablesUI2.Data(1:length(a.LayerOrientations),3:4) = {''};
            a.TablesUI2.Data(1:length(a.LayerOrientations),6) = {''};
        case 4
            a.MaterialUI2.String = fieldnames(a.Materials.Isotropic);
            E = fieldnames(a.Materials.Isotropic);
            [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Isotropic,E{1}));
            a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material2{1}.Name));
            a.TablesUI2.Data(1:length(a.LayerOrientations),3:5) = {''};
        end
        a.TablesUI2.Data(1:length(a.LayerOrientations),a.MaterialType2+2) = {a.Material2{1}.Name};
    end
    function MaterialUI2_Callback(source,~,~)
        if  a.MaterialType2 == 1
            [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material2{1}.Name));
            a.TablesUI2.Data(1:length(a.LayerOrientations),4:6) = {''};
        elseif a.MaterialType2 == 2
            [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material2{1}.Name));
            a.TablesUI2.Data(1:length(a.LayerOrientations),3) = {''};
            a.TablesUI2.Data(1:length(a.LayerOrientations),5:6) = {''};
        elseif a.MaterialType2 == 3
            [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Cubic,cell2mat(source.String(source.Value))));
            a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material2{1}.Name));
            a.TablesUI2.Data(1:length(a.LayerOrientations),3:4) = {''};
            a.TablesUI2.Data(1:length(a.LayerOrientations),6) = {''};
        elseif a.MaterialType2 == 4
            [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Isotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material2{1}.Name));
            a.TablesUI2.Data(1:length(a.LayerOrientations),3:5) = {''};
        end
        a.TablesUI2.Data(1:length(a.LayerOrientations),a.MaterialType2+2) = {a.Material2{1}.Name};        
    end
    function UniformLayerThicknessUI2_Callback(source,~,~)
        a.UniformLayerThickness = source.Value;
        if  source.Value == 1
            a.PlateThicknessUI2.Enable = 'on';
            a.PlateThicknessUI2Enable = 'on';
            if  a.Hybrid == 1
                a.TablesUI2.ColumnEditable = [true false true true true true true];
                a.TablesUI2ColumnEditable = [true false true true true true true];
            else
                a.TablesUI2.ColumnEditable = [true false false false false false true];
                a.TablesUI2ColumnEditable = [true false false false false false true];
            end
            if  a.SymmetricSystem == 1
                a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(2*length(a.LayerOrientations)*a.SuperLayers);
            else
                a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(length(a.LayerOrientations)*a.SuperLayers);
            end
            a.TablesUI2.Data(1:length(a.LayerOrientations),2) = {a.LayerThicknesses(1)};
        else
            a.PlateThicknessUI2.Enable = 'off';
            a.PlateThicknessUI2Enable = 'off';
            if  a.Hybrid == 1
                a.TablesUI2.ColumnEditable = [true true true true true true true];
                a.TablesUI2ColumnEditable = [true true true true true true true];
            else
                a.TablesUI2.ColumnEditable = [true true false false false false true];
                a.TablesUI2ColumnEditable = [true true false false false false true];
            end
        end
    end
    function PlateThicknessUI2_Callback(source,~,~)
        a.PlateThickness = str2double(source.String);
        if  a.SymmetricSystem == 1
            a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(2*length(a.LayerOrientations)*a.SuperLayers);
        else
            a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(length(a.LayerOrientations)*a.SuperLayers);
        end
        a.TablesUI2.Data(1:length(a.LayerOrientations),2) = {a.LayerThicknesses(1)};
    end
    function SuperLayersUI2_Callback(source,~,~)
        a.SuperLayers = str2double(source.String);
        if  a.UniformLayerThickness == 1
            if  a.SymmetricSystem == 1
                a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(2*length(a.LayerOrientations)*a.SuperLayers);
            else
                a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(length(a.LayerOrientations)*a.SuperLayers);
            end
            a.TablesUI2.Data(1:length(a.LayerOrientations),2) = {a.LayerThicknesses(1)};
        else
            if  a.SymmetricSystem == 1
                a.PlateThickness = 2*a.SuperLayers*sum(a.LayerThicknesses);
            else
                a.PlateThickness = a.SuperLayers*sum(a.LayerThicknesses);
            end
            a.PlateThicknessUI2.String = a.PlateThickness;
        end
    end
    function SymmetricSystemUI2_Callback(source,~,~)
        a.SymmetricSystem = source.Value;
        if  a.UniformLayerThickness == 1
            if  a.SymmetricSystem == 1
                a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(2*length(a.LayerOrientations)*a.SuperLayers);
            else
                a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(length(a.LayerOrientations)*a.SuperLayers);
            end
            a.TablesUI2.Data(1:length(a.LayerOrientations),2) = {a.LayerThicknesses(1)};
        else
            if  a.SymmetricSystem == 1
                a.PlateThickness = 2*a.SuperLayers*sum(a.LayerThicknesses);
            else
                a.PlateThickness = a.SuperLayers*sum(a.LayerThicknesses);
            end
            a.PlateThicknessUI2.String = a.PlateThickness;
        end
    end
    function TablesUI2_Callback(hObject,callbackdata)
        if  callbackdata.Indices(2) == 1
            if  ~isempty(callbackdata.EditData)
                a.LayerOrientations(callbackdata.Indices(1)) = str2double(callbackdata.EditData);
            else
                errordlg('Invalid input!','Error');
                hObject.Data(callbackdata.Indices(1),1) = {callbackdata.PreviousData};
                return
            end
            if  a.UniformLayerThickness == 1
                if  a.SymmetricSystem == 1
                    a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(2*length(a.LayerOrientations)*a.SuperLayers);
                else
                    a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(length(a.LayerOrientations)*a.SuperLayers);
                end
                a.TablesUI2.Data(1:length(a.LayerOrientations),2) = {a.LayerThicknesses(1)};
            else
                if  a.SymmetricSystem == 1
                    a.PlateThickness = 2*a.SuperLayers*sum(a.LayerThicknesses);
                else
                    a.PlateThickness = a.SuperLayers*sum(a.LayerThicknesses);
                end
                a.PlateThicknessUI2.String = a.PlateThickness;
            end
            if  a.Hybrid == 0
                [a.Material2{1:length(a.LayerOrientations)}] = deal(a.Material2{1});
                hObject.Data(1:length(a.LayerOrientations),a.MaterialType2+2) = {a.Material2{1}.Name};
                if  a.MaterialType2 == 1
                    hObject.Data(1:length(a.LayerOrientations),4:6) = {''};
                elseif a.MaterialType2 == 2
                    hObject.Data(1:length(a.LayerOrientations),3) = {''};
                    hObject.Data(1:length(a.LayerOrientations),5:6) = {''};
                elseif a.MaterialType2 == 3
                    hObject.Data(1:length(a.LayerOrientations),3:4) = {''};
                    hObject.Data(1:length(a.LayerOrientations),6) = {''}; 
                elseif a.MaterialType2 == 4
                    hObject.Data(1:length(a.LayerOrientations),3:5) = {''};
                end
            end
        elseif callbackdata.Indices(2) == 2
            if  ~isempty(callbackdata.EditData)
                a.LayerThicknesses(callbackdata.Indices(1)) = str2double(callbackdata.EditData);
            else
                errordlg('Invalid input!','Error');
                hObject.Data(callbackdata.Indices(1),2) = {callbackdata.PreviousData};
                return
            end
            if  a.SymmetricSystem == 1
                a.PlateThickness = 2*a.SuperLayers*sum(a.LayerThicknesses);
            else
                a.PlateThickness = a.SuperLayers*sum(a.LayerThicknesses);
            end
            a.PlateThicknessUI2.String = a.PlateThickness;
            if  a.Hybrid == 0
                [a.Material2{1:length(a.LayerThicknesses)}] = deal(a.Material2{1});
                hObject.Data(1:length(a.LayerThicknesses),a.MaterialType2+2) = {a.Material2{1}.Name};
                if  a.MaterialType2 == 1
                    hObject.Data(1:length(a.LayerThicknesses),4:6) = {''};
                elseif a.MaterialType2 == 2
                    hObject.Data(1:length(a.LayerThicknesses),3) = {''};
                    hObject.Data(1:length(a.LayerThicknesses),5:6) = {''};
                elseif a.MaterialType2 == 3
                    hObject.Data(1:length(a.LayerThicknesses),3:4) = {''};
                    hObject.Data(1:length(a.LayerThicknesses),6) = {''}; 
                elseif a.MaterialType2 == 4
                    hObject.Data(1:length(a.LayerThicknesses),3:5) = {''};
                end
            end
        elseif callbackdata.Indices(2) == 3
            a.Material2{callbackdata.Indices(1)} = getfield(a.Materials.Orthotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),4:6) = {''};
        elseif callbackdata.Indices(2) == 4
            a.Material2{callbackdata.Indices(1)} = getfield(a.Materials.TransverselyIsotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3) = {''};
            hObject.Data(callbackdata.Indices(1),5:6) = {''};
        elseif callbackdata.Indices(2) == 5
            a.Material2{callbackdata.Indices(1)} = getfield(a.Materials.Cubic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3:4) = {''};
            hObject.Data(callbackdata.Indices(1),6) = {''};
        elseif callbackdata.Indices(2) == 6
            a.Material2{callbackdata.Indices(1)} = getfield(a.Materials.Isotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3:5) = {''};
        elseif callbackdata.Indices(2) == 7 
            if  length(a.LayerOrientations) > 1 || length(a.LayerThicknesses) > 1 || length(a.Material2) > 1
                if  length(a.LayerOrientations) >= callbackdata.Indices(1)
                    a.LayerOrientations(callbackdata.Indices(1)) = [];
                end
                if  length(a.LayerThicknesses) >= callbackdata.Indices(1)
                    a.LayerThicknesses(callbackdata.Indices(1)) = [];
                end
                if  length(a.Material2) >= callbackdata.Indices(1)
                    a.Material2(callbackdata.Indices(1)) = [];
                    a.Material2(cellfun(@isempty,a.Material2)) = [];
                end
                if  length(a.MaterialNames) >= callbackdata.Indices(1)
                    a.MaterialNames(callbackdata.Indices(1)) = [];
                end
                if  length(a.MaterialClasses) >= callbackdata.Indices(1)
                    a.MaterialClasses(callbackdata.Indices(1)) = [];
                end
                hObject.Data(callbackdata.Indices(1),:) = [];
                hObject.Data = vertcat(hObject.Data,cell(1,6));
            else
                errordlg('Do not delete the only layer!','Error');
                return
            end
            if  a.UniformLayerThickness == 1
                if  a.SymmetricSystem == 1
                    a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(2*length(a.LayerOrientations)*a.SuperLayers);
                else
                    a.LayerThicknesses(1:length(a.LayerOrientations)) = a.PlateThickness/(length(a.LayerOrientations)*a.SuperLayers);
                end
                a.TablesUI2.Data(1:length(a.LayerOrientations),2) = {a.LayerThicknesses(1)};
            else
                if  a.SymmetricSystem == 1
                    a.PlateThickness = 2*a.SuperLayers*sum(a.LayerThicknesses);
                else
                    a.PlateThickness = a.SuperLayers*sum(a.LayerThicknesses);
                end
                a.PlateThicknessUI2.String = a.PlateThickness;
            end
            if  a.Hybrid == 0
                [a.Material2{1:length(a.LayerThicknesses)}] = deal(a.Material2{1});
                hObject.Data(1:length(a.LayerThicknesses),a.MaterialType2+2) = {a.Material2{1}.Name};
                if  a.MaterialType2 == 1
                    hObject.Data(1:length(a.LayerThicknesses),4:6) = {''};
                elseif a.MaterialType2 == 2
                    hObject.Data(1:length(a.LayerThicknesses),3) = {''};
                    hObject.Data(1:length(a.LayerThicknesses),5:6) = {''};
                elseif a.MaterialType2 == 3
                    hObject.Data(1:length(a.LayerThicknesses),3:4) = {''};
                    hObject.Data(1:length(a.LayerThicknesses),6) = {''}; 
                elseif a.MaterialType2 == 4
                    hObject.Data(1:length(a.LayerThicknesses),3:5) = {''};
                end
            end
        end
        try
            for i = 1:length(a.Material2)
                a.MaterialNames{i} = a.Material2{i}.Name;
                a.MaterialClasses{i} = a.Material2{i}.Class;
            end
        catch
            
        end
    end
    function OKUI2_Callback(~,~)
        try
            for i = 1:length(a.Material2)
                a.MaterialNames{i} = a.Material2{i}.Name;
                a.MaterialClasses{i} = a.Material2{i}.Class;
            end
        catch
            errordlg('A layer in between misses a material!','Error');
            return  
        end
        if  length(a.LayerOrientations) ~= length(a.LayerThicknesses) && a.UniformLayerThickness == 0
            errordlg('If the laminate does not have uniform layer thicknesses, you must enter every layer''s thickness in the ''d (mm)'' column! All columns must have the same number of entries.','Error');
            return
        elseif length(a.LayerOrientations) ~= length(a.Material2) && a.Hybrid == 1
            errordlg('If you set up a hybrid layup, you must assign a material to every layer defined in the ''Phi'' column and vice versa!','Error');
            return
        elseif length(a.LayerOrientations) == 1 && a.SuperLayers > 1
            errordlg('Do not define multiple repetitions if the laminate consists of only one layer! Increase the layer thickness instead.','Error');
            a.SuperLayers = 1;
            a.SuperLayersUI2.String = 1;
            a.LayerThicknesses = a.PlateThickness;
            a.TablesUI2.Data(1,2) = {a.LayerThicknesses(1)};
            return
        elseif length(a.LayerOrientations) == 1 && a.SymmetricSystem == 1
            errordlg('Do not define a symmetric layup if the laminate consists of only one layer! Set the double layer thickness instead.','Error');
            a.SymmetricSystem = 0;
            a.SymmetricSystemUI2.Value = 0;
            a.LayerThicknesses = a.PlateThickness;
            a.TablesUI2.Data(1,2) = {a.LayerThicknesses(1)};
            return
        elseif length(a.LayerOrientations) > 1 && all(strcmp(a.MaterialClasses,'Isotropic')) && all(strcmp(a.MaterialNames(1),a.MaterialNames))
            errordlg('Do not define a unit cell containing solely layers of the same isotropic material! Replace by one layer of greater thickness. This will be the same because the DC assumes rigid bonding between the layers. It is also more efficient since the computational expense scales with the number of layers.','Error');
            return
        elseif length(a.LayerOrientations) > 1 && all(a.LayerOrientations(1) == a.LayerOrientations) && all(strcmp(a.MaterialNames(1),a.MaterialNames))
            errordlg('Do not define a unit cell containing solely layers of the same material and orientation! Replace by one layer of greater thickness. This will be the same because the DC assumes rigid bonding between the layers. It is also more efficient since the computational expense scales with the number of layers.','Error');
            return
        elseif a.SymmetricSystem == 0 && length(a.LayerOrientations) > 1 &&...
            (mod(length(a.LayerOrientations),2) == 0 && all(fliplr(a.LayerOrientations(1:length(a.LayerOrientations)/2)) == a.LayerOrientations(length(a.LayerOrientations)/2+1:end)) && all(fliplr(a.LayerThicknesses(1:length(a.LayerThicknesses)/2)) == a.LayerThicknesses(length(a.LayerThicknesses)/2+1:end)) && all(strcmp(fliplr(a.MaterialNames(1:length(a.MaterialNames)/2)),a.MaterialNames(length(a.MaterialNames)/2+1:end))) ||...
            mod(length(a.LayerOrientations),2) ~= 0 && all(fliplr(a.LayerOrientations(1:length(a.LayerOrientations)/2-.5)) == a.LayerOrientations(length(a.LayerOrientations)/2+1.5:end)) && all(fliplr(a.LayerThicknesses(1:length(a.LayerThicknesses)/2-.5)) == a.LayerThicknesses(length(a.LayerThicknesses)/2+1.5:end)) && all(strcmp(fliplr(a.MaterialNames(1:length(a.MaterialNames)/2-.5)),a.MaterialNames(length(a.MaterialNames)/2+1.5:end))))
            errordlg('Do not define a symmetric unit cell! Cut the unit cell in half and check ''Symmetric system'' instead. This will be more robust and more efficient.','Error');
            return
        elseif all(strcmp(a.MaterialNames(1),a.MaterialNames)) && a.Hybrid == 1
            if  length(a.LayerOrientations) > 1
                errordlg('It is no hybrid if all layers are the same material!','Error');
            else
                errordlg('It is no hybrid if you have only one layer!','Error');
            end
            a.Hybrid = 0;
            a.HybridUI2.Value = 0;
            
            a.MaterialTypeUI2.Enable = 'on';
            a.MaterialTypeUI2Enable = 'on';
            a.MaterialUI2.Enable = 'on';
            a.MaterialUI2Enable = 'on';
            if  a.UniformLayerThickness == 1
                a.TablesUI2.ColumnEditable = [true false false false false false true];
                a.TablesUI2ColumnEditable = [true false false false false false true];
            else
                a.TablesUI2.ColumnEditable = [true true false false false false true];
                a.TablesUI2ColumnEditable = [true true false false false false true];
            end
            
            if  a.MaterialType2 == 1
                [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(a.MaterialUI2.String(a.MaterialUI2.Value))));
                a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material2{1}.Name));
                a.TablesUI2.Data(1:length(a.LayerOrientations),4:6) = {''};
            elseif a.MaterialType2 == 2
                [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(a.MaterialUI2.String(a.MaterialUI2.Value))));
                a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material2{1}.Name));
                a.TablesUI2.Data(1:length(a.LayerOrientations),3) = {''};
                a.TablesUI2.Data(1:length(a.LayerOrientations),5:6) = {''};
            elseif a.MaterialType2 == 3
                [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Cubic,cell2mat(a.MaterialUI2.String(a.MaterialUI2.Value))));
                a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material2{1}.Name));
                a.TablesUI2.Data(1:length(a.LayerOrientations),3:4) = {''};
                a.TablesUI2.Data(1:length(a.LayerOrientations),6) = {''};
            elseif a.MaterialType2 == 4
                [a.Material2{1:length(a.LayerOrientations)}] = deal(getfield(a.Materials.Isotropic,cell2mat(a.MaterialUI2.String(a.MaterialUI2.Value))));
                a.MaterialUI2Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material2{1}.Name));
                a.TablesUI2.Data(1:length(a.LayerOrientations),3:5) = {''};
            end
            a.TablesUI2.Data(1:length(a.LayerOrientations),a.MaterialType2+2) = {a.Material2{1}.Name};
            return
        end
        a.UnitCell = a.TablesUI2.Data;
        a.ThicknessCountUI2.String = a.PlateThickness;
        if  a.ToggleUpperFluid2 == 1 || a.ToggleLowerFluid2 == 1
            a.FluidLoading2 = 1;
            a.ScholteModesUI2.Enable = 'on';
            a.HalfspacesNumber1UI2.Enable = 'on';
            a.Halfspaces1UI2.Enable = 'on';
            a.HalfspacesNumber2UI2.Enable = 'on';
            a.Halfspaces2UI2.Enable = 'on';
            a.Couplant2 = a.UpperFluid2;
            if  a.Quantity12 == 4
                a.Option1UI2.Value = find(strcmp(a.Couplant2.Name,a.Option1UI2.String));
                x = round(30*a.Couplant2.Velocity/343);
                if  x > 90
                    x = 90;
                end
                a.YAxisUI2.String = ['[0 ',num2str(x),']'];
                a.YAxis2 = eval(a.YAxisUI2.String);
            end
        else
            a.FluidLoading2 = 0;
            a.ScholteModesUI2.Enable = 'off';
            a.HalfspacesNumber1UI2.Enable = 'off';
            a.Halfspaces1UI2.Enable = 'off';
            a.HalfspacesNumber2UI2.Enable = 'off';
            a.Halfspaces2UI2.Enable = 'off';
        end
        if  a.ToggleUpperFluid2 == 1
            a.UpperFluidDisplayUI2.String = a.UpperFluid2.Name;
        else
            a.UpperFluidDisplayUI2.String = 'vacuum';
        end
        if  a.ToggleLowerFluid2 == 1
            a.LowerFluidDisplayUI2.String = a.LowerFluid2.Name;
        else
            a.LowerFluidDisplayUI2.String = 'vacuum';
        end
        if  a.SymmetricSystem == 1
            a.LayerCountUI2.String = 2*length(a.LayerOrientations)*a.SuperLayers;
        else
            a.LayerCountUI2.String = length(a.LayerOrientations)*a.SuperLayers;
        end
        if  (a.ToggleUpperFluid2 == 0 && a.ToggleLowerFluid2 == 0 && (length(a.LayerOrientations) == 1 || a.SymmetricSystem == 1)) || (a.ToggleUpperFluid2 == 1 && a.ToggleLowerFluid2 == 1 && strcmp(a.UpperFluid2.Name,a.LowerFluid2.Name) && (length(a.LayerOrientations) == 1 || a.SymmetricSystem == 1))
            a.Symmetric = 1;
        else
            a.Symmetric = 0;
        end
        if  a.Symmetric == 1
            a.SymmetricModesUI2.Enable = 'on';
            a.AntisymmetricModesUI2.Enable = 'on';
        else
            a.SymmetricModesUI2.Enable = 'off';
            a.AntisymmetricModesUI2.Enable = 'off';
        end
        if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
            a.PhaseVelocityLimit2 = round(a.YRange*a.Material2{1}.PlateVelocity,-3);
            a.FrequencyLimit2 = round(a.Material2{1}.PlateVelocity*a.XRange1/a.PlateThickness,-2);
            if  a.FrequencyLimit2 == 0
                a.FrequencyLimit2 = 50;
            elseif a.FrequencyLimit2 >= 1e4
                a.FrequencyLimit2 = round(a.FrequencyLimit2,-3);
            end
            a.FrequencyResolution2 = a.FrequencyLimit2/a.XSamples1;
            a.Step2 = a.FrequencyLimit2/a.Steps1;
        else
            a.PhaseVelocityLimit2 = 20e3;
            a.FrequencyLimit2 = 1e3*round(a.XRange2/a.PlateThickness,1);
            if  a.FrequencyLimit2 == 0
                a.FrequencyLimit2 = 50;
            end
            a.FrequencyResolution2 = a.FrequencyLimit2/a.XSamples2;
            a.Step2 = a.FrequencyLimit2/a.Steps2;
        end
        a.PhaseVelocityLimitUI2.String = a.PhaseVelocityLimit2/1e3;
        a.FrequencyLimitUI2.String = a.FrequencyLimit2;
        a.FrequencyResolutionUI2.String = a.FrequencyResolution2;
        a.SteptUI2.String = a.Step2;
        a.Frequency12 = a.FrequencyLimit2;
        a.Frequency1UI2.String = a.FrequencyLimit2;
        a.Frequency22 = a.FrequencyLimit2;
        a.Frequency2UI2.String = a.FrequencyLimit2;
        Phi = a.LayerOrientations-a.PropagationAngle;
        for i = 1:length(a.MaterialClasses)
            if  strcmp(a.MaterialClasses(i),'Isotropic')
                DC(i) = 1;
            elseif strcmp(a.MaterialClasses(i),'Cubic')
                if  mod(Phi(i),45) == 0
                    DC(i) = 1;
                else
                    DC(i) = 0;
                end
            elseif strcmp(a.MaterialClasses(i),'Transversely isotropic') || strcmp(a.MaterialClasses(i),'Orthotropic')
                if  mod(Phi(i),90) == 0
                    DC(i) = 1;
                else
                    DC(i) = 0;
                end
            end
            if  isreal(a.Material2{i}.C)
                a.Viscoelastic2(i) = 0;
            else
                a.Viscoelastic2(i) = 1;
            end
        end
        if  all(DC)
            a.ShearHorizontalModesUI2.Enable = 'on';
        else
            a.ShearHorizontalModesUI2.Enable = 'off';
        end
        if  any(a.Viscoelastic2)
            a.Viscoelastic2 = 1;
        else
            a.Viscoelastic2 = 0;
        end        
        if  a.Quantity12 == 1
            if  strcmp(a.LayerCountUI2.String,'1')
                a.Option1UI2.Enable = 'on';
                a.Option1UI2.Style = 'checkbox';
                a.Option1UI2.String = ' ';
                a.Option1UI2.Value = a.BulkVelocities2;
                a.Option1UI2.TooltipString = 'Check this to show the bulk wave velocities.';
                a.Option1UI2.Position(3) = 50;
                a.Option1TextUI2.String = 'Bulk velocities';
                a.Option1TextUI2.Position(3) = 70;
            else
                a.Option1UI2.Enable = 'off';
            end
            a.XAxisModeTextUI2.String = 'X-axis mode';
            a.XAxisModeTextUI2.Position(3) = 62;
            a.XAxisModeUI2.TooltipString = 'Select the frequency''s dimension on the X-axis.';
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            a.YAxisTextUI2.String = 'Y-axis (m/ms)';
            a.YAxisTextUI2.Position(3) = 70;
            a.YAxisUI2.TooltipString = 'Enter which phase velocity range shall be plotted.';
            if  a.XAxisMode2 == 1
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
                a.XAxis2 = eval(a.XAxisUI2.String);
            elseif a.XAxisMode2 == 2
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
            elseif a.XAxisMode2 == 3
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;                    
            end
            a.YAxisUI2.String = ['[0 ',num2str(a.PhaseVelocityLimit2/1e3),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        elseif a.Quantity12 == 2
            if  strcmp(a.LayerCountUI2.String,'1')
                a.Option1UI2.Enable = 'on';
                a.Option1UI2.Style = 'checkbox';
                a.Option1UI2.String = ' ';
                a.Option1UI2.Value = a.BulkVelocities2;
                a.Option1UI2.TooltipString = 'Check this to show the bulk wave velocities.';
                a.Option1UI2.Position(3) = 50;
                a.Option1TextUI2.String = 'Bulk velocities';
                a.Option1TextUI2.Position(3) = 70;
            else
                a.Option1UI2.Enable = 'off';
            end
            a.XAxisModeTextUI2.String = 'X-axis mode';
            a.XAxisModeTextUI2.Position(3) = 62;
            a.XAxisModeUI2.TooltipString = 'Select the frequency''s dimension on the X-axis.';
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            a.YAxisTextUI2.String = 'Y-axis (m/ms)';
            a.YAxisTextUI2.Position(3) = 70;
            a.YAxisUI2.TooltipString = 'Enter which phase velocity range shall be plotted.';
            if  a.XAxisMode2 == 1
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
                a.XAxis2 = eval(a.XAxisUI2.String);
            elseif a.XAxisMode2 == 2
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
            elseif a.XAxisMode2 == 3
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;                    
            end
            if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                a.YAxisUI2.String = ['[0 ',num2str(ceil(a.Material2{1}.PlateVelocity/1e3)),']'];
            else
                a.YAxisUI2.String = '[0 11]';
            end
            a.YAxis2 = eval(a.YAxisUI2.String);
        elseif a.Quantity12 == 3
            if  a.XAxisMode2 == 1
                a.YAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
                a.YAxis2 = eval(a.YAxisUI2.String);
            elseif a.XAxisMode2 == 2
                a.YAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
                a.YAxis2 = eval(a.YAxisUI2.String)*1e3;
            elseif a.XAxisMode2 == 3
                a.YAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
                a.YAxis2 = eval(a.YAxisUI2.String)*1e3/a.PlateThickness;                    
            end
            if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                x = a.Distance2/a.Material2{1}.PlateVelocity*5e3;
            else
                x = 2*a.Distance2;
            end
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
            a.XAxisUI2.String = ['[0 ',num2str(x),']'];
            a.XAxis2 = eval(a.XAxisUI2.String);
        elseif a.Quantity12 == 4
            if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                x = round(8e4/a.Material2{1}.PlateVelocity*a.Couplant2.Velocity/343);
            else
                x = round(30*a.Couplant2.Velocity/343);
            end
            if  x > 90
                x = 90;
            end
            a.YAxisUI2.String = ['[0 ',num2str(x),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
            if  a.XAxisMode2 == 1
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
                a.XAxis2 = eval(a.XAxisUI2.String);
            elseif a.XAxisMode2 == 2
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
            elseif a.XAxisMode2 == 3
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;                    
            end
        elseif a.Quantity12 == 5
            if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                x = a.YRange*1e-3*a.Material2{1}.PlateVelocity*a.PlateThickness;
                if  x >= 1e3 && x < 1e4
                    x = round(x,-2);
                elseif x >= 1e2 && x < 1e3
                    x = round(x,-1);
                elseif x < 1e2
                    x = round(x);
                end
            end
            if  a.XAxisMode2 == 3
                if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                    a.YAxisUI2.String = ['[0 ',num2str(x/a.PlateThickness),']'];
                else
                    a.YAxisUI2.String = '[0 50]';
                end
                a.YAxis2 = eval(a.YAxisUI2.String);                
            else
                if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                    a.YAxisUI2.String = ['[0 ',num2str(x),']'];
                else
                    a.YAxisUI2.String = ['[0 ',num2str(50*a.PlateThickness),']'];
                end
                a.YAxis2 = eval(a.YAxisUI2.String);
            end
            if  a.XAxisMode2 == 1
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
                a.XAxis2 = eval(a.XAxisUI2.String);
            elseif a.XAxisMode2 == 2
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
            elseif a.XAxisMode2 == 3
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;                    
            end
        elseif a.Quantity12 == 6
            if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                x = 2*pi*a.FrequencyLimit2/a.Material2{1}.RayleighVelocity;
            else
                x = 2*pi*a.FrequencyLimit2/1.6e3;
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
            if  a.XAxisMode2 == 3
                a.YAxisUI2.String = ['[0 ',num2str(x*a.PlateThickness),']'];
                a.YAxis2 = eval(a.YAxisUI2.String);
            else
                a.YAxisUI2.String = ['[0 ',num2str(x),']'];
                a.YAxis2 = eval(a.YAxisUI2.String);
            end            
            if  a.XAxisMode2 == 1
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
                a.XAxis2 = eval(a.XAxisUI2.String);            
            elseif a.XAxisMode2 == 2
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3;            
            elseif a.XAxisMode2 == 3
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;             
            end 
        elseif a.Quantity12 == 7
            if  a.Symmetric == 1
                Fluid = a.UpperFluid2;
                if  strcmp(a.Material2{1}.Class,'Isotropic')
                    MaterialVelocity = a.Material2{1}.PlateVelocity;
                else
                    MaterialVelocity = a.Material2{1}.LongitudinalVelocity_1;
                end
                MaterialDensity = a.Material2{1}.Density;
                FluidVelocity = Fluid.Velocity;
                FluidDensity = Fluid.Density;
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
                MaterialVelocity = .5*(UpperMaterialVelocity+LowerMaterialVelocity);
                MaterialDensity = .5*(a.Material2{1}.Density+a.Material2{end}.Density);
                if  a.ToggleUpperFluid2 == 1 && a.ToggleLowerFluid2 == 1
                    FluidVelocity = .5*(a.UpperFluid2.Velocity+a.LowerFluid2.Velocity);
                    FluidDensity = .5*(a.UpperFluid2.Density+a.LowerFluid2.Density);
                elseif a.ToggleUpperFluid2 == 1 && a.ToggleLowerFluid2 == 0
                    FluidVelocity = .5*a.UpperFluid2.Velocity;
                    FluidDensity = .5*a.UpperFluid2.Density;
                elseif a.ToggleUpperFluid2 == 0 && a.ToggleLowerFluid2 == 1
                    FluidVelocity = .5*a.LowerFluid2.Velocity;
                    FluidDensity = .5*a.LowerFluid2.Density;
                end
            end
            if  a.Viscoelastic2 == 1
                for i = 1:length(a.Material2)
                    if  ~isreal(a.Material2{i}.C)
                        xV = pi*(imag(a.Material2{i}.C(1,1))/real(a.Material2{i}.C(1,1))+imag(a.Material2{i}.C(6,6))/real(a.Material2{i}.C(6,6)));
                        break
                    end
                end
            end
            if  a.FluidLoading2 == 1 && a.Viscoelastic2 == 1
                x = 1e4*(FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity)+xV)/a.PlateThickness;
            elseif a.FluidLoading2 == 1 && a.Viscoelastic2 == 0
                x = 1e4*FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity)/a.PlateThickness;
            elseif a.FluidLoading2 == 0 && a.Viscoelastic2 == 1
                x = 1e4*xV/a.PlateThickness;
            else
                x = 1;
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
            if  a.XAxisMode2 == 3
                a.YAxisUI2.String = ['[0 ',num2str(x*a.PlateThickness),']'];
                a.YAxis2 = eval(a.YAxisUI2.String);
            else
                a.YAxisUI2.String = ['[0 ',num2str(x),']'];
                a.YAxis2 = eval(a.YAxisUI2.String);
            end            
            if  a.XAxisMode2 == 1
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
                a.XAxis2 = eval(a.XAxisUI2.String);            
            elseif a.XAxisMode2 == 2
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3;            
            elseif a.XAxisMode2 == 3
                a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;             
            end
        end
        a.Samples12 = round(a.Samples11/str2double(a.LayerCountUI2.String));
        a.Samples1UI2.String = a.Samples12;
        if  a.FluidLoading2 == 1 && a.Halfspaces22 == 1
            if  a.ToggleUpperFluid2 == 1 && a.ToggleLowerFluid2 == 1
                a.Samples32 = round(.5*a.Samples22/((2*a.HalfspacesNumber22+1)*str2double(a.LayerCountUI2.String)));
            else
                a.Samples32 = round(.5*a.Samples22/((a.HalfspacesNumber22+1)*str2double(a.LayerCountUI2.String)));
            end
        else
            a.Samples32 = round(.5*a.Samples22/str2double(a.LayerCountUI2.String));
        end
        if  a.Samples32 == 0
            a.Samples32 = 1;
        end
        a.Samples3UI2.String = a.Samples32;
        if  a.Hybrid == 0
            a.MaterialNameUI2.String = a.Material2{1}.Name;
        elseif a.Hybrid == 1
            a.MaterialNameUI2.String = 'Hybrid';
        end
        a.LayupString1 = char(join(split(num2str(a.LayerOrientations)),'/'));
        a.EffectiveLayupString1 = char(join(split(num2str(a.LayerOrientations-a.PropagationAngle)),'/'));
        if  a.SymmetricSystem == 0
            if  a.SuperLayers == 1
                a.LayupUI2.String = ['[',a.LayupString1,']'];
                a.EffectiveLayupUI2.String = ['[',a.EffectiveLayupString1,']'];
            elseif a.SuperLayers > 1
                a.LayupUI2.String = ['[',a.LayupString1,']',num2str(a.SuperLayers)];
                a.EffectiveLayupUI2.String = ['[',a.EffectiveLayupString1,']',num2str(a.SuperLayers)];
            end
        else
            if  a.SuperLayers == 1
                a.LayupUI2.String = ['[',a.LayupString1,']s'];
                a.EffectiveLayupUI2.String = ['[',a.EffectiveLayupString1,']s'];
            elseif a.SuperLayers > 1
                a.LayupUI2.String = ['[',a.LayupString1,']',num2str(a.SuperLayers),'s'];
                a.EffectiveLayupUI2.String = ['[',a.EffectiveLayupString1,']',num2str(a.SuperLayers),'s'];
            end
        end
        close(f2)
        String = 'Unit cell:';
        for i = 1:length(a.MaterialNames)
            String = append(String,newline,num2str(i),': ',a.MaterialNames{i},' (',a.MaterialClasses{i},')');
        end
        disp([String,newline,'-----------------------------------'])
    end
    function CancelUI2_Callback(~,~)
        close(f2)
    end
    a.Detect2 = 0;
end
function SpecimenSettingsUI4_Callback(~,~) % Tab4_polar diagrams
    f3 = figure('NumberTitle','off','Name','Specimen_Polar diagrams','Visible','off','MenuBar','none','Position',[0 0 910 340]);
    jframe = get(gcf,'javaframe');
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
    f3.Units = 'normalized';
    movegui(f3,'center')
    f3.Visible = 'on';
        
    uicontrol('Parent',f3,'Style','pushbutton','String','Open','TooltipString','Open an existing specimen definition.','Position',[20 290 95 33],'FontSize',10,'Callback',@OpenUI4_Callback);
    uicontrol('Parent',f3,'Style','pushbutton','String','Save','TooltipString','Save the current specimen definition.','Position',[20+130 290 95 33],'FontSize',10,'Callback',@SaveUI4_Callback);
    uicontrol('Parent',f3,'Style','pushbutton','String','Reset','TooltipString','Reset specimen definition.','Position',[20+260 290 95 33],'FontSize',10,'Callback',@ResetUI4_Callback);    
    
    uicontrol('Parent',f3,'Style','text','String','Hybrid','Position',[20 70+185 33 13])
    uicontrol('Parent',f3,'Style','text','String','Class','Position',[20 70+155 29 13])
    uicontrol('Parent',f3,'Style','text','String','Material','Position',[20 70+125 39 13]);
    uicontrol('Parent',f3,'Style','text','String','Uniform layer thickness','Position',[20 70+95 115 13]);
    uicontrol('Parent',f3,'Style','text','String','Total thickness (mm)','Position',[20 70+65 101 13]);
    uicontrol('Parent',f3,'Style','text','String','Unit cell repetitions','Position',[20 70+35 92 13]);
    uicontrol('Parent',f3,'Style','text','String','Symmetric system','Position',[20 70+5 90 13]);
    a.HybridUI4 = uicontrol('Parent',f3,'Style','checkbox','Value',a.Hybrid_Polar,'TooltipString','The layup may contain different materials','Position',[150 70+180 50 23],'Callback',@HybridUI4_Callback);
    a.MaterialTypeUI4 = uicontrol('Parent',f3,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic','Isotropic'},'Value',a.MaterialType_Polar,'TooltipString','Select a material symmetry class.','Enable',a.MaterialTypeUI4Enable,'Position',[150-75 70+150 170 23],'Callback',@MaterialTypeUI4_Callback);
    if  a.MaterialType_Polar == 1
        a.MaterialUI4 = uicontrol('Parent',f3,'Style','popupmenu','Value',a.MaterialUI4Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Orthotropic),'Enable',a.MaterialUI4Enable,'Position',[150-75 70+120 170 23],'Callback',@MaterialUI4_Callback);
    elseif a.MaterialType_Polar == 2
        a.MaterialUI4 = uicontrol('Parent',f3,'Style','popupmenu','Value',a.MaterialUI4Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.TransverselyIsotropic),'Enable',a.MaterialUI4Enable,'Position',[150-75 70+120 170 23],'Callback',@MaterialUI4_Callback);
    elseif a.MaterialType_Polar == 3
        a.MaterialUI4 = uicontrol('Parent',f3,'Style','popupmenu','Value',a.MaterialUI4Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Cubic),'Enable',a.MaterialUI4Enable,'Position',[150-75 70+120 170 23],'Callback',@MaterialUI4_Callback);
    elseif a.MaterialType_Polar == 4
        a.MaterialUI4 = uicontrol('Parent',f3,'Style','popupmenu','Value',a.MaterialUI4Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Isotropic),'Enable',a.MaterialUI4Enable,'Position',[150-75 70+120 170 23],'Callback',@MaterialUI4_Callback);
    end
    a.UniformLayerThicknessUI4 = uicontrol('Parent',f3,'Style','checkbox','Value',a.UniformLayerThickness_Polar,'TooltipString',['Check this if the layup has uniform layer thicknesses. Then, the Dispersion Calculator',newline,'deduces the layer thicknesses from the overall plate thickness as defined below, and',newline,'you do not need to enter the individual layer thicknesses into the table to the right.'],'Position',[150 70+90 50 23],'Callback',@UniformLayerThicknessUI4_Callback);
    a.PlateThicknessUI4 = uicontrol('Parent',f3,'Style','edit','String',a.PlateThickness_Polar,'TooltipString',['Enter the overall plate thickness. This is available',newline,'only if you have checked ''Uniform layer thickness''.'],'Enable',a.PlateThicknessUI4Enable,'Position',[150 70+60 50 23],'Callback',@PlateThicknessUI4_Callback);
    a.SuperLayersUI4 = uicontrol('Parent',f3,'Style','edit','String',a.SuperLayers_Polar,'TooltipString',['Enter the number of repetitions of the unit cell',newline,'contained in the laminate. Use only integers.'],'Position',[150 70+30 50 23],'Callback',@SuperLayersUI4_Callback);
    a.SymmetricSystemUI4 = uicontrol('Parent',f3,'Style','checkbox','Value',a.SymmetricSystem_Polar,'TooltipString',['Check this if the laminate is symmetric about its middle plane.',newline,'The Dispersion Calculator extends the laminate accordingly.'],'Position',[150 70+0 50 23],'Callback',@SymmetricSystemUI4_Callback);    

    uicontrol('Parent',f3,'Style','pushbutton','String','OK','TooltipString','Accept the current specimen definition and continue.','Position',[20 20 95 33],'FontSize',10,'Callback',@OKUI4_Callback);
    uicontrol('Parent',f3,'Style','pushbutton','String','Cancel','TooltipString','Discard the current specimen definition.','Position',[20+130 20 95 33],'FontSize',10,'Callback',@CancelUI4_Callback);
    
    uicontrol('Parent',f3,'Style','text','String','Unit cell','Position',[280 20+246 46 13],'FontSize',9);
    a.TablesUI4 = uitable(f3,'ColumnName',{['Phi (',char(176),')'],'d (mm)','Orthotropic','Trans. iso.','Cubic','Isotropic','Delete'},'ColumnWidth',{50 50 100 100 100 100 50},'ColumnFormat',({[] [] fieldnames(a.Materials.Orthotropic)' fieldnames(a.Materials.TransverselyIsotropic)' fieldnames(a.Materials.Cubic)' fieldnames(a.Materials.Isotropic)' 'logical'}),'TooltipString',['Enter the fiber orientation of every layer into the left column, and if you don''t have checked',newline,'''Uniform layer thickness'', enter also the corresponding layer thicknesses into the second column.'],'Data',a.UnitCell_Polar,'ColumnEditable',a.TablesUI4ColumnEditable,'RowStriping','off','Position',[280 20 613 237],'CellEditCallback',@TablesUI4_Callback);

    function OpenUI4_Callback(~,~)
        [Path,File] = uigetfile('*.mat');
        if  Path ~= 0
            D4 = load(fullfile(File,Path));
            a.Hybrid_Polar = D4.D5.Hybrid;
            a.MaterialType_Polar = D4.D5.MaterialType;
            a.Material_Polar = D4.D5.Material;
            a.MaterialNames_Polar = D4.D5.MaterialNames;
            a.MaterialClasses_Polar = D4.D5.MaterialClasses;
            a.UniformLayerThickness_Polar = D4.D5.UniformLayerThickness;
            a.PlateThickness_Polar = D4.D5.PlateThickness;
            a.SuperLayers_Polar = D4.D5.SuperLayers;
            a.SymmetricSystem_Polar = D4.D5.SymmetricSystem;
            a.UnitCell_Polar = D4.D5.UnitCell;
            a.LayerOrientations_Polar = D4.D5.LayerOrientations;
            a.LayerThicknesses_Polar = D4.D5.LayerThicknesses;

            a.HybridUI4.Value = a.Hybrid_Polar;
            a.MaterialTypeUI4.Value = a.MaterialType_Polar;
            if  a.Hybrid_Polar == 0
                if  a.MaterialType_Polar == 1
                    a.MaterialUI4.String = fieldnames(a.Materials.Orthotropic);
                    a.MaterialUI4.Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material_Polar{1}.Name));
                elseif a.MaterialType_Polar == 2
                    a.MaterialUI4.String = fieldnames(a.Materials.TransverselyIsotropic);
                    a.MaterialUI4.Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material_Polar{1}.Name));
                elseif a.MaterialType_Polar == 3
                    a.MaterialUI4.String = fieldnames(a.Materials.Cubic);
                    a.MaterialUI4.Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material_Polar{1}.Name));
                elseif a.MaterialType_Polar == 4
                    a.MaterialUI4.String = fieldnames(a.Materials.Isotropic);
                    a.MaterialUI4.Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material_Polar{1}.Name));
                end
            end
            a.MaterialUI4Value = a.MaterialUI4.Value;
            if  a.Hybrid_Polar == 1
                a.MaterialTypeUI4.Enable = 'off';
                a.MaterialTypeUI4Enable = 'off';
                a.MaterialUI4.Enable = 'off';
                a.MaterialUI4Enable = 'off';
            else
                a.MaterialTypeUI4.Enable = 'on';
                a.MaterialTypeUI4Enable = 'on';
                a.MaterialUI4.Enable = 'on';
                a.MaterialUI4Enable = 'on';
            end
            a.UniformLayerThicknessUI4.Value = a.UniformLayerThickness_Polar;
            if  a.UniformLayerThickness_Polar == 1
                a.PlateThicknessUI4.Enable = 'on';
                a.PlateThicknessUI4Enable = 'on';
                if  a.Hybrid_Polar == 1
                    a.TablesUI4.ColumnEditable = [true false true true true true true];
                    a.TablesUI4ColumnEditable = [true false true true true true true];
                else
                    a.TablesUI4.ColumnEditable = [true false false false false false true];
                    a.TablesUI4ColumnEditable = [true false false false false false true];
                end
            else
                a.PlateThicknessUI4.Enable = 'off';
                a.PlateThicknessUI4Enable = 'off';
                if  a.Hybrid_Polar == 1
                    a.TablesUI4.ColumnEditable = [true true true true true true true];
                    a.TablesUI4ColumnEditable = [true true true true true true true];
                else
                    a.TablesUI4.ColumnEditable = [true true false false false false true];
                    a.TablesUI4ColumnEditable = [true true false false false false true];
                end
            end
            a.PlateThicknessUI4.String = a.PlateThickness_Polar;
            a.SuperLayersUI4.String = a.SuperLayers_Polar;
            a.SymmetricSystemUI4.Value = a.SymmetricSystem_Polar;
            a.TablesUI4.Data = a.UnitCell_Polar;
        end
    end
    function SaveUI4_Callback(~,~)
        D5.Hybrid = a.Hybrid_Polar;
        D5.MaterialType = a.MaterialType_Polar;
        D5.Material = a.Material_Polar;
        D5.MaterialNames = a.MaterialNames_Polar;
        D5.MaterialClasses = a.MaterialClasses_Polar;
        D5.UniformLayerThickness = a.UniformLayerThickness_Polar;
        D5.PlateThickness = a.PlateThickness_Polar;
        D5.SuperLayers = a.SuperLayers_Polar;
        D5.SymmetricSystem = a.SymmetricSystem_Polar;
        D5.UnitCell = a.TablesUI4.Data;
        D5.LayerOrientations = a.LayerOrientations_Polar;
        D5.LayerThicknesses = a.LayerThicknesses_Polar;
        uisave('D5',fullfile(a.Directory,'Specimen'))
    end
    function ResetUI4_Callback(~,~)
        a.Hybrid_Polar = 0;
        a.MaterialType_Polar = 1;
        a.MaterialTypeUI4Enable = 'on';
        E = fieldnames(a.Materials.Orthotropic);
        a.Material_Polar{1} = getfield(a.Materials.Orthotropic,E{1});
        a.MaterialNames_Polar{1} = a.Material_Polar{1}.Name;
        a.MaterialClasses_Polar{1} = a.Material_Polar{1}.Class;
        if  length(a.Material_Polar) > 1
            a.Material_Polar(2:end) = [];
            a.MaterialNames_Polar(2:end) = [];
            a.MaterialClasses_Polar(2:end) = [];
        end
        a.MaterialUI4Value = 1;
        a.MaterialUI4Enable = 'on';
        a.UniformLayerThickness_Polar = 1;
        a.PlateThicknessUI4Enable = 'on';
        a.TablesUI4.ColumnEditable = [true false false false false false true];
        a.TablesUI4ColumnEditable = [true false false false false false true];
        a.PlateThickness_Polar = 1;
        a.SuperLayers_Polar = 1;
        a.SymmetricSystem_Polar = 0;
        a.UnitCell_Polar = cell(400,6);
        a.UnitCell_Polar{1} = '0';
        a.UnitCell_Polar{1,2} = '1';
        a.UnitCell_Polar{1,3} = E{1};
        a.LayerOrientations_Polar = 0;
        a.LayerThicknesses_Polar = 1;
        
        a.HybridUI4.Value = a.Hybrid_Polar;
        a.MaterialTypeUI4.Value = a.MaterialType_Polar;
        a.MaterialTypeUI4.Enable = a.MaterialTypeUI4Enable;
        a.MaterialUI4.Value = a.MaterialUI4Value;
        a.MaterialUI4.String = E;
        a.MaterialUI4.Enable = a.MaterialUI4Enable;
        a.UniformLayerThicknessUI4.Value = a.UniformLayerThickness_Polar;
        a.PlateThicknessUI4.String = a.PlateThickness_Polar;
        a.SuperLayersUI4.String = a.SuperLayers_Polar;
        a.SymmetricSystemUI4.Value = a.SymmetricSystem_Polar;
        a.TablesUI4.Data = a.UnitCell_Polar;
    end
    function HybridUI4_Callback(source,~,~)
        a.Hybrid_Polar = source.Value;
        if  source.Value == 1
            a.MaterialTypeUI4.Enable = 'off';
            a.MaterialTypeUI4Enable = 'off';
            a.MaterialUI4.Enable = 'off';
            a.MaterialUI4Enable = 'off';
            if  a.UniformLayerThickness_Polar == 1
                a.TablesUI4.ColumnEditable = [true false true true true true true];
                a.TablesUI4ColumnEditable = [true false true true true true true];
            else
                a.TablesUI4.ColumnEditable = [true true true true true true true];
                a.TablesUI4ColumnEditable = [true true true true true true true];
            end
        else
            a.MaterialTypeUI4.Enable = 'on';
            a.MaterialTypeUI4Enable = 'on';
            a.MaterialUI4.Enable = 'on';
            a.MaterialUI4Enable = 'on';
            if  a.UniformLayerThickness_Polar == 1
                a.TablesUI4.ColumnEditable = [true false false false false false true];
                a.TablesUI4ColumnEditable = [true false false false false false true];
            else
                a.TablesUI4.ColumnEditable = [true true false false false false true];
                a.TablesUI4ColumnEditable = [true true false false false false true];
            end
            
            if  a.MaterialType_Polar == 1
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(a.MaterialUI4.String(a.MaterialUI4.Value))));
                a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material_Polar{1}.Name));
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),4:6) = {''};
            elseif a.MaterialType_Polar == 2
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(a.MaterialUI4.String(a.MaterialUI4.Value))));
                a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material_Polar{1}.Name));
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3) = {''};
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),5:6) = {''};
            elseif a.MaterialType_Polar == 3
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Cubic,cell2mat(a.MaterialUI4.String(a.MaterialUI4.Value))));
                a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material_Polar{1}.Name));
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3:4) = {''};
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),6) = {''};
            elseif a.MaterialType_Polar == 4
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Isotropic,cell2mat(a.MaterialUI4.String(a.MaterialUI4.Value))));
                a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material_Polar{1}.Name));
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3:5) = {''};
            end
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),a.MaterialType_Polar+2) = {a.Material_Polar{1}.Name};
        end
    end
    function MaterialTypeUI4_Callback(source,~,~)
        a.MaterialType_Polar = source.Value;
        a.MaterialUI4.Value = 1;
        switch source.Value
        case 1
            a.MaterialUI4.String = fieldnames(a.Materials.Orthotropic);
            E = fieldnames(a.Materials.Orthotropic);
            [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Orthotropic,E{1}));
            a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material_Polar{1}.Name));
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),4:6) = {''};
        case 2
            a.MaterialUI4.String = fieldnames(a.Materials.TransverselyIsotropic);
            E = fieldnames(a.Materials.TransverselyIsotropic);
            [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.TransverselyIsotropic,E{1}));
            a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material_Polar{1}.Name));
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3) = {''};
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),5:6) = {''};
        case 3
            a.MaterialUI4.String = fieldnames(a.Materials.Cubic);
            E = fieldnames(a.Materials.Cubic);
            [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Cubic,E{1}));
            a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material_Polar{1}.Name));
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3:4) = {''};
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),6) = {''};
        case 4
            a.MaterialUI4.String = fieldnames(a.Materials.Isotropic);
            E = fieldnames(a.Materials.Isotropic);
            [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Isotropic,E{1}));
            a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material_Polar{1}.Name));
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3:5) = {''};
        end
        a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),a.MaterialType_Polar+2) = {a.Material_Polar{1}.Name};
    end
    function MaterialUI4_Callback(source,~,~)
        if  a.MaterialType_Polar == 1
            [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material_Polar{1}.Name));
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),4:6) = {''};
        elseif a.MaterialType_Polar == 2
            [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material_Polar{1}.Name));
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3) = {''};
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),5:6) = {''};
        elseif a.MaterialType_Polar == 3
            [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Cubic,cell2mat(source.String(source.Value))));
            a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material_Polar{1}.Name));
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3:4) = {''};
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),6) = {''};
        elseif a.MaterialType_Polar == 4
            [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Isotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material_Polar{1}.Name));
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3:5) = {''};
        end
        a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),a.MaterialType_Polar+2) = {a.Material_Polar{1}.Name};        
    end
    function UniformLayerThicknessUI4_Callback(source,~,~)
        a.UniformLayerThickness_Polar = source.Value;
        if  source.Value == 1
            a.PlateThicknessUI4.Enable = 'on';
            a.PlateThicknessUI4Enable = 'on';
            if  a.Hybrid_Polar == 1
                a.TablesUI4.ColumnEditable = [true false true true true true true];
                a.TablesUI4ColumnEditable = [true false true true true true true];
            else
                a.TablesUI4.ColumnEditable = [true false false false false false true];
                a.TablesUI4ColumnEditable = [true false false false false false true];
            end
            if  a.SymmetricSystem_Polar == 1
                a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(2*length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
            else
                a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
            end
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),2) = {a.LayerThicknesses_Polar(1)};
        else
            a.PlateThicknessUI4.Enable = 'off';
            a.PlateThicknessUI4Enable = 'off';
            if  a.Hybrid_Polar == 1
                a.TablesUI4.ColumnEditable = [true true true true true true true];
                a.TablesUI4ColumnEditable = [true true true true true true true];
            else
                a.TablesUI4.ColumnEditable = [true true false false false false true];
                a.TablesUI4ColumnEditable = [true true false false false false true];
            end
        end
    end
    function PlateThicknessUI4_Callback(source,~,~)
        a.PlateThickness_Polar = str2double(source.String);
        if  a.SymmetricSystem_Polar == 1
            a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(2*length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
        else
            a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
        end
        a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),2) = {a.LayerThicknesses_Polar(1)};
    end
    function SuperLayersUI4_Callback(source,~,~)
        a.SuperLayers_Polar = str2double(source.String);
        if  a.UniformLayerThickness_Polar == 1
            if  a.SymmetricSystem_Polar == 1
                a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(2*length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
            else
                a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
            end
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),2) = {a.LayerThicknesses_Polar(1)};
        else
            if  a.SymmetricSystem_Polar == 1
                a.PlateThickness_Polar = 2*a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
            else
                a.PlateThickness_Polar = a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
            end
            a.PlateThicknessUI4.String = a.PlateThickness_Polar;
        end
    end
    function SymmetricSystemUI4_Callback(source,~,~)
        a.SymmetricSystem_Polar = source.Value;
        if  a.UniformLayerThickness_Polar == 1
            if  a.SymmetricSystem_Polar == 1
                a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(2*length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
            else
                a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
            end
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),2) = {a.LayerThicknesses_Polar(1)};
        else
            if  a.SymmetricSystem_Polar == 1
                a.PlateThickness_Polar = 2*a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
            else
                a.PlateThickness_Polar = a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
            end
            a.PlateThicknessUI4.String = a.PlateThickness_Polar;
        end
    end
    function TablesUI4_Callback(hObject,callbackdata)
        if  callbackdata.Indices(2) == 1
            if  ~isempty(callbackdata.EditData)
                a.LayerOrientations_Polar(callbackdata.Indices(1)) = str2double(callbackdata.EditData);
            else
                errordlg('Invalid input!','Error');
                hObject.Data(callbackdata.Indices(1),1) = {callbackdata.PreviousData};
                return
            end
            if  a.UniformLayerThickness_Polar == 1
                if  a.SymmetricSystem_Polar == 1
                    a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(2*length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
                else
                    a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
                end
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),2) = {a.LayerThicknesses_Polar(1)};
            else
                if  a.SymmetricSystem_Polar == 1
                    a.PlateThickness_Polar = 2*a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
                else
                    a.PlateThickness_Polar = a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
                end
                a.PlateThicknessUI4.String = a.PlateThickness_Polar;
            end
            if  a.Hybrid_Polar == 0
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(a.Material_Polar{1});
                hObject.Data(1:length(a.LayerOrientations_Polar),a.MaterialType_Polar+2) = {a.Material_Polar{1}.Name};
                if  a.MaterialType_Polar == 1
                    hObject.Data(1:length(a.LayerOrientations_Polar),4:6) = {''};
                elseif a.MaterialType_Polar == 2
                    hObject.Data(1:length(a.LayerOrientations_Polar),3) = {''};
                    hObject.Data(1:length(a.LayerOrientations_Polar),5:6) = {''};
                elseif a.MaterialType_Polar == 3
                    hObject.Data(1:length(a.LayerOrientations_Polar),3:4) = {''};
                    hObject.Data(1:length(a.LayerOrientations_Polar),6) = {''}; 
                elseif a.MaterialType_Polar == 4
                    hObject.Data(1:length(a.LayerOrientations_Polar),3:5) = {''};
                end
            end
        elseif callbackdata.Indices(2) == 2
            if  ~isempty(callbackdata.EditData)
                a.LayerThicknesses_Polar(callbackdata.Indices(1)) = str2double(callbackdata.EditData);
            else
                errordlg('Invalid input!','Error');
                hObject.Data(callbackdata.Indices(1),1) = {callbackdata.PreviousData};
                return
            end
            if  a.SymmetricSystem_Polar == 1
                a.PlateThickness_Polar = 2*a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
            else
                a.PlateThickness_Polar = a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
            end
            a.PlateThicknessUI4.String = a.PlateThickness_Polar;
            if  a.Hybrid_Polar == 0
                [a.Material_Polar{1:length(a.LayerThicknesses_Polar)}] = deal(a.Material_Polar{1});
                hObject.Data(1:length(a.LayerThicknesses_Polar),a.MaterialType_Polar+2) = {a.Material_Polar{1}.Name};
                if  a.MaterialType_Polar == 1
                    hObject.Data(1:length(a.LayerThicknesses_Polar),4:6) = {''};
                elseif a.MaterialType_Polar == 2
                    hObject.Data(1:length(a.LayerThicknesses_Polar),3) = {''};
                    hObject.Data(1:length(a.LayerThicknesses_Polar),5:6) = {''};
                elseif a.MaterialType_Polar == 3
                    hObject.Data(1:length(a.LayerThicknesses_Polar),3:4) = {''};
                    hObject.Data(1:length(a.LayerThicknesses_Polar),6) = {''}; 
                elseif a.MaterialType_Polar == 4
                    hObject.Data(1:length(a.LayerThicknesses_Polar),3:5) = {''};
                end
            end
        elseif callbackdata.Indices(2) == 3
            a.Material_Polar{callbackdata.Indices(1)} = getfield(a.Materials.Orthotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),4:6) = {''};
        elseif callbackdata.Indices(2) == 4
            a.Material_Polar{callbackdata.Indices(1)} = getfield(a.Materials.TransverselyIsotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3) = {''};
            hObject.Data(callbackdata.Indices(1),5:6) = {''};
        elseif callbackdata.Indices(2) == 5
            a.Material_Polar{callbackdata.Indices(1)} = getfield(a.Materials.Cubic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3:4) = {''};
            hObject.Data(callbackdata.Indices(1),6) = {''};
        elseif callbackdata.Indices(2) == 6
            a.Material_Polar{callbackdata.Indices(1)} = getfield(a.Materials.Isotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3:5) = {''};
        elseif callbackdata.Indices(2) == 7
            if  length(a.LayerOrientations_Polar) > 1 || length(a.LayerThicknesses_Polar) > 1 || length(a.Material_Polar) > 1
                if  length(a.LayerOrientations_Polar) >= callbackdata.Indices(1)
                    a.LayerOrientations_Polar(callbackdata.Indices(1)) = [];
                end
                if  length(a.LayerThicknesses_Polar) >= callbackdata.Indices(1)
                    a.LayerThicknesses_Polar(callbackdata.Indices(1)) = [];
                end
                if  length(a.Material_Polar) >= callbackdata.Indices(1)
                    a.Material_Polar(callbackdata.Indices(1)) = [];
                    a.Material_Polar(cellfun(@isempty,a.Material_Polar)) = [];
                end
                if  length(a.MaterialNames_Polar) >= callbackdata.Indices(1)
                    a.MaterialNames_Polar(callbackdata.Indices(1)) = [];
                end
                if  length(a.MaterialClasses_Polar) >= callbackdata.Indices(1)
                    a.MaterialClasses_Polar(callbackdata.Indices(1)) = [];
                end
                hObject.Data(callbackdata.Indices(1),:) = [];
                hObject.Data = vertcat(hObject.Data,cell(1,6));
            else
                errordlg('Do not delete the only layer!','Error');
                return
            end
            if  a.UniformLayerThickness_Polar == 1
                if  a.SymmetricSystem_Polar == 1
                    a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(2*length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
                else
                    a.LayerThicknesses_Polar(1:length(a.LayerOrientations_Polar)) = a.PlateThickness_Polar/(length(a.LayerOrientations_Polar)*a.SuperLayers_Polar);
                end
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),2) = {a.LayerThicknesses_Polar(1)};
            else
                if  a.SymmetricSystem_Polar == 1
                    a.PlateThickness_Polar = 2*a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
                else
                    a.PlateThickness_Polar = a.SuperLayers_Polar*sum(a.LayerThicknesses_Polar);
                end
                a.PlateThicknessUI4.String = a.PlateThickness_Polar;
            end
            if  a.Hybrid_Polar == 0
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(a.Material_Polar{1});
                hObject.Data(1:length(a.LayerOrientations_Polar),a.MaterialType_Polar+2) = {a.Material_Polar{1}.Name};
                if  a.MaterialType_Polar == 1
                    hObject.Data(1:length(a.LayerOrientations_Polar),4:6) = {''};
                elseif a.MaterialType_Polar == 2
                    hObject.Data(1:length(a.LayerOrientations_Polar),3) = {''};
                    hObject.Data(1:length(a.LayerOrientations_Polar),5:6) = {''};
                elseif a.MaterialType_Polar == 3
                    hObject.Data(1:length(a.LayerOrientations_Polar),3:4) = {''};
                    hObject.Data(1:length(a.LayerOrientations_Polar),6) = {''}; 
                elseif a.MaterialType_Polar == 4
                    hObject.Data(1:length(a.LayerOrientations_Polar),3:5) = {''};
                end
            end
        end
        try
            for i = 1:length(a.Material_Polar)
                a.MaterialNames_Polar{i} = a.Material_Polar{i}.Name;
                a.MaterialClasses_Polar{i} = a.Material_Polar{i}.Class;
            end
        catch
            
        end
    end
    function OKUI4_Callback(~,~)
        try
            for i = 1:length(a.Material_Polar)
                a.MaterialNames_Polar{i} = a.Material_Polar{i}.Name;
                a.MaterialClasses_Polar{i} = a.Material_Polar{i}.Class;
            end
        catch
            errordlg('A layer in between misses a material!','Error');
            return  
        end    
        if  length(a.LayerOrientations_Polar) ~= length(a.LayerThicknesses_Polar) && a.UniformLayerThickness_Polar == 0
            errordlg('If the laminate does not have uniform layer thicknesses, you must enter every layer''s thickness in the ''d (mm)'' column! All columns must have the same number of entries.','Error');
            return
        elseif length(a.LayerOrientations_Polar) ~= length(a.Material_Polar) && a.Hybrid_Polar == 1
            errordlg('If you set up a hybrid layup, you must assign a material to every layer defined in the ''Phi'' column and vice versa!','Error');
            return
        elseif length(a.LayerOrientations_Polar) == 1 && a.SuperLayers_Polar > 1
            errordlg('Do not define multiple repetitions if the laminate consists of only one layer! Increase the layer thickness instead.','Error');
            a.SuperLayers_Polar = 1;
            a.SuperLayersUI4.String = 1;
            a.LayerThicknesses_Polar = a.PlateThickness_Polar;
            a.TablesUI4.Data(1,2) = {a.LayerThicknesses_Polar(1)};
            return
        elseif length(a.LayerOrientations_Polar) == 1 && a.SymmetricSystem_Polar == 1
            errordlg('Do not define a symmetric layup if the laminate consists of only one layer! Set the double layer thickness instead.','Error');
            a.SymmetricSystem_Polar = 0;
            a.SymmetricSystemUI4.Value = 0;
            a.LayerThicknesses_Polar = a.PlateThickness_Polar;
            a.TablesUI4.Data(1,2) = {a.LayerThicknesses_Polar(1)};
            return  
        elseif length(a.LayerOrientations_Polar) > 1 && all(strcmp(a.MaterialClasses_Polar,'Isotropic')) && all(strcmp(a.MaterialNames_Polar(1),a.MaterialNames_Polar))
            errordlg('Do not define a unit cell containing solely layers of the same isotropic material! Replace by one layer of greater thickness. This will be the same because the DC assumes rigid bonding between the layers. It is also more efficient since the computational expense scales with the number of layers.','Error');
            return
        elseif length(a.LayerOrientations_Polar) > 1 && all(a.LayerOrientations_Polar(1) == a.LayerOrientations_Polar) && all(strcmp(a.MaterialNames_Polar(1),a.MaterialNames_Polar))
            errordlg('Do not define a unit cell containing solely layers of the same material and orientation! Replace by one layer of greater thickness. This will be the same because the DC assumes rigid bonding between the layers. It is also more efficient since the computational expense scales with the number of layers.','Error');
            return
        elseif a.SymmetricSystem_Polar == 0 && length(a.LayerOrientations_Polar) > 1 &&...
            (mod(length(a.LayerOrientations_Polar),2) == 0 && all(fliplr(a.LayerOrientations_Polar(1:length(a.LayerOrientations_Polar)/2)) == a.LayerOrientations_Polar(length(a.LayerOrientations_Polar)/2+1:end)) && all(fliplr(a.LayerThicknesses_Polar(1:length(a.LayerThicknesses_Polar)/2)) == a.LayerThicknesses_Polar(length(a.LayerThicknesses_Polar)/2+1:end)) && all(strcmp(fliplr(a.MaterialNames_Polar(1:length(a.MaterialNames_Polar)/2)),a.MaterialNames_Polar(length(a.MaterialNames_Polar)/2+1:end))) ||...
            mod(length(a.LayerOrientations_Polar),2) ~= 0 && all(fliplr(a.LayerOrientations_Polar(1:length(a.LayerOrientations_Polar)/2-.5)) == a.LayerOrientations_Polar(length(a.LayerOrientations_Polar)/2+1.5:end)) && all(fliplr(a.LayerThicknesses_Polar(1:length(a.LayerThicknesses_Polar)/2-.5)) == a.LayerThicknesses_Polar(length(a.LayerThicknesses_Polar)/2+1.5:end)) && all(strcmp(fliplr(a.MaterialNames_Polar(1:length(a.MaterialNames_Polar)/2-.5)),a.MaterialNames_Polar(length(a.MaterialNames_Polar)/2+1.5:end))))
            errordlg('Do not define a symmetric unit cell! Cut the unit cell in half and check ''Symmetric system'' instead. This will be more robust and more efficient.','Error');
            return
        elseif all(strcmp(a.MaterialNames_Polar(1),a.MaterialNames_Polar)) && a.Hybrid_Polar == 1
            if  length(a.LayerOrientations_Polar) > 1
                errordlg('It is no hybrid if all layers are the same material!','Error');
            else
                errordlg('It is no hybrid if you have only one layer!','Error');
            end
            a.Hybrid_Polar = 0;
            a.HybridUI4.Value = 0;
            
            a.MaterialTypeUI4.Enable = 'on';
            a.MaterialTypeUI4Enable = 'on';
            a.MaterialUI4.Enable = 'on';
            a.MaterialUI4Enable = 'on';
            if  a.UniformLayerThickness_Polar == 1
                a.TablesUI4.ColumnEditable = [true false false false false false true];
                a.TablesUI4ColumnEditable = [true false false false false false true];
            else
                a.TablesUI4.ColumnEditable = [true true false false false false true];
                a.TablesUI4ColumnEditable = [true true false false false false true];
            end
            
            if  a.MaterialType_Polar == 1
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(a.MaterialUI4.String(a.MaterialUI4.Value))));
                a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material_Polar{1}.Name));
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),4:6) = {''};
            elseif a.MaterialType_Polar == 2
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(a.MaterialUI4.String(a.MaterialUI4.Value))));
                a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material_Polar{1}.Name));
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3) = {''};
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),5:6) = {''};
            elseif a.MaterialType_Polar == 3
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Cubic,cell2mat(a.MaterialUI4.String(a.MaterialUI4.Value))));
                a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material_Polar{1}.Name));
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3:4) = {''};
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),6) = {''};
            elseif a.MaterialType_Polar == 4
                [a.Material_Polar{1:length(a.LayerOrientations_Polar)}] = deal(getfield(a.Materials.Isotropic,cell2mat(a.MaterialUI4.String(a.MaterialUI4.Value))));
                a.MaterialUI4Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material_Polar{1}.Name));
                a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),3:5) = {''};
            end
            a.TablesUI4.Data(1:length(a.LayerOrientations_Polar),a.MaterialType_Polar+2) = {a.Material_Polar{1}.Name};
            
            return
        elseif ~any(a.LayerOrientations_Polar == 0)
            errordlg('The layup must contain at least one 0 deg layer!','Error');
            return
        end
        a.UnitCell_Polar = a.TablesUI4.Data;
        a.ThicknessCountUI4.String = a.PlateThickness_Polar;
        if  a.SymmetricSystem_Polar == 1
            a.LayerCountUI4.String = 2*length(a.LayerOrientations_Polar)*a.SuperLayers_Polar;
        else
            a.LayerCountUI4.String = length(a.LayerOrientations_Polar)*a.SuperLayers_Polar;
        end
        a.FrequencyLimit_Polar = 1e3*round(a.XRange_Polar/a.PlateThickness_Polar,1);
        if  a.FrequencyLimit_Polar == 0
            a.FrequencyLimit_Polar = 50;
        end
        a.FrequencyLimitUI4.String = a.FrequencyLimit_Polar;
        a.FrequencyResolution_Polar = a.FrequencyLimit_Polar/a.XSamples_Polar;
        a.FrequencyResolutionUI4.String = a.FrequencyResolution_Polar;
        a.Frequency_Polar = a.FrequencyLimit_Polar;
        a.FrequencyUI4.String = a.FrequencyLimit_Polar;
        if  a.Quantity_Polar == 1
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
        elseif a.Quantity_Polar == 2
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
        end
        if  a.Hybrid_Polar == 0
            a.MaterialNameUI4.String = a.Material_Polar{1}.Name;
        elseif a.Hybrid_Polar == 1
            a.MaterialNameUI4.String = 'Hybrid';
        end
        a.LayupString_Polar = char(join(split(num2str(a.LayerOrientations_Polar)),'/'));
        a.EffectiveLayupString1 = char(join(split(num2str(a.LayerOrientations_Polar-a.PropagationAngle)),'/'));
        if  a.SymmetricSystem_Polar == 0
            if  a.SuperLayers_Polar == 1
                a.LayupUI4.String = ['[',a.LayupString_Polar,']'];
                a.EffectiveLayupUI4.String = ['[',a.EffectiveLayupString1,']'];
            elseif a.SuperLayers_Polar > 1
                a.LayupUI4.String = ['[',a.LayupString_Polar,']',num2str(a.SuperLayers_Polar)];
                a.EffectiveLayupUI4.String = ['[',a.EffectiveLayupString1,']',num2str(a.SuperLayers_Polar)];
            end
        else
            if  a.SuperLayers_Polar == 1
                a.LayupUI4.String = ['[',a.LayupString_Polar,']s'];
                a.EffectiveLayupUI4.String = ['[',a.EffectiveLayupString1,']s'];
            elseif a.SuperLayers_Polar > 1
                a.LayupUI4.String = ['[',a.LayupString_Polar,']',num2str(a.SuperLayers_Polar),'s'];
                a.EffectiveLayupUI4.String = ['[',a.EffectiveLayupString1,']',num2str(a.SuperLayers_Polar),'s'];
            end
        end
        close(f3)
        String = 'Unit cell:';
        for i = 1:length(a.MaterialNames_Polar)
            String = append(String,newline,num2str(i),': ',a.MaterialNames_Polar{i},' (',a.MaterialClasses_Polar{i},')');
        end
        disp([String,newline,'-----------------------------------'])
    end
    function CancelUI4_Callback(~,~)
        close(f3)
    end
end
function SpecimenSettingsUI6_Callback(~,~) % Tab6_laminate stiffness
    f4 = figure('NumberTitle','off','Name','Specimen_Laminate stiffness','Visible','off','MenuBar','none','Position',[0 0 910 340]);
    jframe = get(gcf,'javaframe');
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
    f4.Units = 'normalized';
    movegui(f4,'center')
    f4.Visible = 'on';
        
    uicontrol('Parent',f4,'Style','pushbutton','String','Open','TooltipString','Open an existing specimen definition.','Position',[20 290 95 33],'FontSize',10,'Callback',@OpenUI6_Callback);
    uicontrol('Parent',f4,'Style','pushbutton','String','Reset','TooltipString','Reset specimen definition.','Position',[20+260 290 95 33],'FontSize',10,'Callback',@ResetUI6_Callback);    
    
    uicontrol('Parent',f4,'Style','text','String','Hybrid','Position',[20 70+185 33 13])
    uicontrol('Parent',f4,'Style','text','String','Class','Position',[20 70+155 29 13])
    uicontrol('Parent',f4,'Style','text','String','Material','Position',[20 70+125 39 13]);
    uicontrol('Parent',f4,'Style','text','String','Uniform layer thickness','Position',[20 70+95 115 13]);
    a.HybridUI6 = uicontrol('Parent',f4,'Style','checkbox','Value',a.Hybrid6,'TooltipString','The layup may contain different materials','Position',[150 70+180 50 23],'Callback',@HybridUI6_Callback);
    a.MaterialTypeUI6 = uicontrol('Parent',f4,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic','Isotropic'},'Value',a.MaterialType6,'TooltipString','Select a material symmetry class.','Enable',a.MaterialTypeUI6Enable,'Position',[150-75 70+150 170 23],'Callback',@MaterialTypeUI6_Callback);
    if  a.MaterialType6 == 1
        a.MaterialUI6 = uicontrol('Parent',f4,'Style','popupmenu','Value',a.MaterialUI6Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Orthotropic),'Enable',a.MaterialUI6Enable,'Position',[150-75 70+120 170 23],'Callback',@MaterialUI6_Callback);
    elseif a.MaterialType6 == 2
        a.MaterialUI6 = uicontrol('Parent',f4,'Style','popupmenu','Value',a.MaterialUI6Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.TransverselyIsotropic),'Enable',a.MaterialUI6Enable,'Position',[150-75 70+120 170 23],'Callback',@MaterialUI6_Callback);
    elseif a.MaterialType6 == 3
        a.MaterialUI6 = uicontrol('Parent',f4,'Style','popupmenu','Value',a.MaterialUI6Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Cubic),'Enable',a.MaterialUI6Enable,'Position',[150-75 70+120 170 23],'Callback',@MaterialUI6_Callback);
    elseif a.MaterialType6 == 4
        a.MaterialUI6 = uicontrol('Parent',f4,'Style','popupmenu','Value',a.MaterialUI6Value,'TooltipString','Select a material.','String',fieldnames(a.Materials.Isotropic),'Enable',a.MaterialUI6Enable,'Position',[150-75 70+120 170 23],'Callback',@MaterialUI6_Callback);
    end
    a.UniformLayerThicknessUI6 = uicontrol('Parent',f4,'Style','checkbox','Value',a.UniformLayerThickness3,'TooltipString',['Check this if the layup has uniform layer thicknesses. Then, the Dispersion Calculator',newline,'deduces the layer thicknesses from the overall plate thickness as defined below, and',newline,'you do not need to enter the individual layer thicknesses into the table to the right.'],'Position',[150 70+90 50 23],'Callback',@UniformLayerThicknessUI6_Callback);

    uicontrol('Parent',f4,'Style','pushbutton','String','OK','TooltipString','Accept the current specimen definition and continue.','Position',[20 20 95 33],'FontSize',10,'Callback',@OKUI6_Callback);
    uicontrol('Parent',f4,'Style','pushbutton','String','Cancel','TooltipString','Discard the current specimen definition.','Position',[20+130 20 95 33],'FontSize',10,'Callback',@CancelUI6_Callback);
    
    uicontrol('Parent',f4,'Style','text','String','Unit cell','Position',[280 20+246 46 13],'FontSize',9);
    a.TablesUI6 = uitable(f4,'ColumnName',{['Phi (',char(176),')'],'d (mm)','Orthotropic','Trans. iso.','Cubic','Isotropic','Delete'},'ColumnWidth',{50 50 100 100 100 100 50},'ColumnFormat',({[] [] fieldnames(a.Materials.Orthotropic)' fieldnames(a.Materials.TransverselyIsotropic)' fieldnames(a.Materials.Cubic)' fieldnames(a.Materials.Isotropic)' 'logical'}),'TooltipString',['Enter the fiber orientation of every layer into the left column, and if you don''t have checked',newline,'''Uniform layer thickness'', enter also the corresponding layer thicknesses into the second column.'],'Data',a.UnitCell3,'ColumnEditable',a.TablesUI6ColumnEditable,'RowStriping','off','Position',[280 20 613 237],'CellEditCallback',@TablesUI6_Callback);

    function OpenUI6_Callback(~,~)
        [Path,File] = uigetfile('*.mat');
        if  Path ~= 0
            D4 = load(fullfile(File,Path));
            a.Hybrid6 = D4.D5.Hybrid;
            a.MaterialType6 = D4.D5.MaterialType;
            a.Material3 = D4.D5.Material;
            a.MaterialNames3 = D4.D5.MaterialNames;
            a.MaterialClasses3 = D4.D5.MaterialClasses;
            a.UniformLayerThickness3 = D4.D5.UniformLayerThickness;
            a.UnitCell3 = D4.D5.UnitCell;
            a.LayerOrientations3 = D4.D5.LayerOrientations;
            a.LayerThicknesses3 = D4.D5.LayerThicknesses;

            a.HybridUI6.Value = a.Hybrid6;
            a.MaterialTypeUI6.Value = a.MaterialType6;
            if  a.Hybrid6 == 0
                if  a.MaterialType6 == 1
                    a.MaterialUI6.String = fieldnames(a.Materials.Orthotropic);
                    a.MaterialUI6.Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material3{1}.Name));
                elseif a.MaterialType6 == 2
                    a.MaterialUI6.String = fieldnames(a.Materials.TransverselyIsotropic);
                    a.MaterialUI6.Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material3{1}.Name));
                elseif a.MaterialType6 == 3
                    a.MaterialUI6.String = fieldnames(a.Materials.Cubic);
                    a.MaterialUI6.Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material3{1}.Name));
                elseif a.MaterialType6 == 4
                    a.MaterialUI6.String = fieldnames(a.Materials.Isotropic);
                    a.MaterialUI6.Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material3{1}.Name));
                end
            end
            a.MaterialUI6Value = a.MaterialUI6.Value;
            if  a.Hybrid6 == 1
                a.MaterialTypeUI6.Enable = 'off';
                a.MaterialTypeUI6Enable = 'off';
                a.MaterialUI6.Enable = 'off';
                a.MaterialUI6Enable = 'off';
            else
                a.MaterialTypeUI6.Enable = 'on';
                a.MaterialTypeUI6Enable = 'on';
                a.MaterialUI6.Enable = 'on';
                a.MaterialUI6Enable = 'on';
            end
            a.UniformLayerThicknessUI6.Value = a.UniformLayerThickness3;
            if  a.UniformLayerThickness3 == 1
                if  a.Hybrid6 == 1
                    a.TablesUI6.ColumnEditable = [true false true true true true true];
                    a.TablesUI6ColumnEditable = [true false true true true true true];
                else
                    a.TablesUI6.ColumnEditable = [true false false false false false true];
                    a.TablesUI6ColumnEditable = [true false false false false false true];
                end
            else
                if  a.Hybrid6 == 1
                    a.TablesUI6.ColumnEditable = [true true true true true true true];
                    a.TablesUI6ColumnEditable = [true true true true true true true];
                else
                    a.TablesUI6.ColumnEditable = [true true false false false false true];
                    a.TablesUI6ColumnEditable = [true true false false false false true];
                end
            end
            a.TablesUI6.Data = a.UnitCell3;
        end
    end
    function ResetUI6_Callback(~,~)
        a.Hybrid6 = 0;
        a.MaterialType6 = 1;
        a.MaterialTypeUI6Enable = 'on';
        E = fieldnames(a.Materials.Orthotropic);
        a.Material3{1} = getfield(a.Materials.Orthotropic,E{1});
        a.MaterialNames3{1} = a.Material3{1}.Name;
        a.MaterialClasses3{1} = a.Material3{1}.Class;
        if  length(a.Material3) > 1
            a.Material3(2:end) = [];
            a.MaterialNames3(2:end) = [];
            a.MaterialClasses3(2:end) = [];
        end
        a.MaterialUI6Value = 1;
        a.MaterialUI6Enable = 'on';
        a.UniformLayerThickness3 = 1;
        a.TablesUI6.ColumnEditable = [true false false false false false true];
        a.TablesUI6ColumnEditable = [true false false false false false true];
        a.UnitCell3 = cell(400,6);
        a.UnitCell3{1} = '0';
        a.UnitCell3{1,2} = '1';
        a.UnitCell3{1,3} = E{1};
        a.LayerOrientations3 = 0;
        a.LayerThicknesses3 = 1;
        
        a.HybridUI6.Value = a.Hybrid6;
        a.MaterialTypeUI6.Value = a.MaterialType6;
        a.MaterialTypeUI6.Enable = a.MaterialTypeUI6Enable;
        a.MaterialUI6.Value = a.MaterialUI6Value;
        a.MaterialUI6.String = E;
        a.MaterialUI6.Enable = a.MaterialUI6Enable;
        a.UniformLayerThicknessUI6.Value = a.UniformLayerThickness3;
        a.TablesUI6.Data = a.UnitCell3;
    end
    function HybridUI6_Callback(source,~,~)
        a.Hybrid6 = source.Value;
        if  source.Value == 1
            a.MaterialTypeUI6.Enable = 'off';
            a.MaterialTypeUI6Enable = 'off';
            a.MaterialUI6.Enable = 'off';
            a.MaterialUI6Enable = 'off';
            if  a.UniformLayerThickness3 == 1
                a.TablesUI6.ColumnEditable = [true false true true true true true];
                a.TablesUI6ColumnEditable = [true false true true true true true];
            else
                a.TablesUI6.ColumnEditable = [true true true true true true true];
                a.TablesUI6ColumnEditable = [true true true true true true true];
            end
        else
            a.MaterialTypeUI6.Enable = 'on';
            a.MaterialTypeUI6Enable = 'on';
            a.MaterialUI6.Enable = 'on';
            a.MaterialUI6Enable = 'on';
            if  a.UniformLayerThickness3 == 1
                a.TablesUI6.ColumnEditable = [true false false false false false true];
                a.TablesUI6ColumnEditable = [true false false false false false true];
            else
                a.TablesUI6.ColumnEditable = [true true false false false false true];
                a.TablesUI6ColumnEditable = [true true false false false false true];
            end
            if  a.MaterialType6 == 1
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(a.MaterialUI6.String(a.MaterialUI6.Value))));
                a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material3{1}.Name));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),4:6) = {''};
            elseif a.MaterialType6 == 2
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(a.MaterialUI6.String(a.MaterialUI6.Value))));
                a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material3{1}.Name));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),3) = {''};
                a.TablesUI6.Data(1:length(a.LayerOrientations3),5:6) = {''};
            elseif a.MaterialType6 == 3
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Cubic,cell2mat(a.MaterialUI6.String(a.MaterialUI6.Value))));
                a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material3{1}.Name));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),3:4) = {''};
                a.TablesUI6.Data(1:length(a.LayerOrientations3),6) = {''};
            elseif a.MaterialType6 == 4
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Isotropic,cell2mat(a.MaterialUI6.String(a.MaterialUI6.Value))));
                a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material3{1}.Name));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),3:5) = {''};
            end
            a.TablesUI6.Data(1:length(a.LayerOrientations3),a.MaterialType6+2) = {a.Material3{1}.Name};
        end
    end
    function MaterialTypeUI6_Callback(source,~,~)
        a.MaterialType6 = source.Value;
        a.MaterialUI6.Value = 1;
        switch source.Value
        case 1
            a.MaterialUI6.String = fieldnames(a.Materials.Orthotropic);
            E = fieldnames(a.Materials.Orthotropic);
            [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Orthotropic,E{1}));
            a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material3{1}.Name));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),4:6) = {''};
        case 2
            a.MaterialUI6.String = fieldnames(a.Materials.TransverselyIsotropic);
            E = fieldnames(a.Materials.TransverselyIsotropic);
            [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.TransverselyIsotropic,E{1}));
            a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material3{1}.Name));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),3) = {''};
            a.TablesUI6.Data(1:length(a.LayerOrientations3),5:6) = {''};
        case 3
            a.MaterialUI6.String = fieldnames(a.Materials.Cubic);
            E = fieldnames(a.Materials.Cubic);
            [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Cubic,E{1}));
            a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material3{1}.Name));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),3:4) = {''};
            a.TablesUI6.Data(1:length(a.LayerOrientations3),6) = {''};
        case 4
            a.MaterialUI6.String = fieldnames(a.Materials.Isotropic);
            E = fieldnames(a.Materials.Isotropic);
            [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Isotropic,E{1}));
            a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material3{1}.Name));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),3:5) = {''};
        end
        a.TablesUI6.Data(1:length(a.LayerOrientations3),a.MaterialType6+2) = {a.Material3{1}.Name};
    end
    function MaterialUI6_Callback(source,~,~)
        if  a.MaterialType6 == 1
            [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material3{1}.Name));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),4:6) = {''};
        elseif a.MaterialType6 == 2
            [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material3{1}.Name));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),3) = {''};
            a.TablesUI6.Data(1:length(a.LayerOrientations3),5:6) = {''};
        elseif a.MaterialType6 == 3
            [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Cubic,cell2mat(source.String(source.Value))));
            a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material3{1}.Name));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),3:4) = {''};
            a.TablesUI6.Data(1:length(a.LayerOrientations3),6) = {''};
        elseif a.MaterialType6 == 4
            [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Isotropic,cell2mat(source.String(source.Value))));
            a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material3{1}.Name));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),3:5) = {''};
        end
        a.TablesUI6.Data(1:length(a.LayerOrientations3),a.MaterialType6+2) = {a.Material3{1}.Name};        
    end
    function UniformLayerThicknessUI6_Callback(source,~,~)
        a.UniformLayerThickness3 = source.Value;
        if  source.Value == 1
            if  a.Hybrid6 == 1
                a.TablesUI6.ColumnEditable = [true false true true true true true];
                a.TablesUI6ColumnEditable = [true false true true true true true];
            else
                a.TablesUI6.ColumnEditable = [true false false false false false true];
                a.TablesUI6ColumnEditable = [true false false false false false true];
            end
            a.LayerThicknesses3(1:length(a.LayerOrientations3)) = sum(a.LayerThicknesses3)/(length(a.LayerOrientations3));
            a.TablesUI6.Data(1:length(a.LayerOrientations3),2) = {a.LayerThicknesses3(1)};
        else
            if  a.Hybrid6 == 1
                a.TablesUI6.ColumnEditable = [true true true true true true true];
                a.TablesUI6ColumnEditable = [true true true true true true true];
            else
                a.TablesUI6.ColumnEditable = [true true false false false false true];
                a.TablesUI6ColumnEditable = [true true false false false false true];
            end
        end
    end
    function TablesUI6_Callback(hObject,callbackdata)
        if  callbackdata.Indices(2) == 1
            if  ~isempty(callbackdata.EditData)
                a.LayerOrientations3(callbackdata.Indices(1)) = str2double(callbackdata.EditData);
            else
                errordlg('Invalid input!','Error');
                hObject.Data(callbackdata.Indices(1),1) = {callbackdata.PreviousData};
                return
            end
            if  a.UniformLayerThickness3 == 1
                a.LayerThicknesses3(1:length(a.LayerOrientations3)) = sum(a.LayerThicknesses3)/(length(a.LayerOrientations3));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),2) = {a.LayerThicknesses3(1)};
            end
            if  a.Hybrid6 == 0
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(a.Material3{1});
                hObject.Data(1:length(a.LayerOrientations3),a.MaterialType6+2) = {a.Material3{1}.Name};
                if  a.MaterialType6 == 1
                    hObject.Data(1:length(a.LayerOrientations3),4:6) = {''};
                elseif a.MaterialType6 == 2
                    hObject.Data(1:length(a.LayerOrientations3),3) = {''};
                    hObject.Data(1:length(a.LayerOrientations3),5:6) = {''};
                elseif a.MaterialType6 == 3
                    hObject.Data(1:length(a.LayerOrientations3),3:4) = {''};
                    hObject.Data(1:length(a.LayerOrientations3),6) = {''}; 
                elseif a.MaterialType6 == 4
                    hObject.Data(1:length(a.LayerOrientations3),3:5) = {''};
                end
            end
        elseif callbackdata.Indices(2) == 2
            if  ~isempty(callbackdata.EditData)
                a.LayerThicknesses3(callbackdata.Indices(1)) = str2double(callbackdata.EditData);
            else
                errordlg('Invalid input!','Error');
                hObject.Data(callbackdata.Indices(1),1) = {callbackdata.PreviousData};
                return
            end
            if  a.Hybrid6 == 0
                [a.Material3{1:length(a.LayerThicknesses3)}] = deal(a.Material3{1});
                hObject.Data(1:length(a.LayerThicknesses3),a.MaterialType6+2) = {a.Material3{1}.Name};
                if  a.MaterialType6 == 1
                    hObject.Data(1:length(a.LayerThicknesses3),4:6) = {''};
                elseif a.MaterialType6 == 2
                    hObject.Data(1:length(a.LayerThicknesses3),3) = {''};
                    hObject.Data(1:length(a.LayerThicknesses3),5:6) = {''};
                elseif a.MaterialType6 == 3
                    hObject.Data(1:length(a.LayerThicknesses3),3:4) = {''};
                    hObject.Data(1:length(a.LayerThicknesses3),6) = {''}; 
                elseif a.MaterialType6 == 4
                    hObject.Data(1:length(a.LayerThicknesses3),3:5) = {''};
                end
            end
        elseif callbackdata.Indices(2) == 3
            a.Material3{callbackdata.Indices(1)} = getfield(a.Materials.Orthotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),4:6) = {''};
        elseif callbackdata.Indices(2) == 4
            a.Material3{callbackdata.Indices(1)} = getfield(a.Materials.TransverselyIsotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3) = {''};
            hObject.Data(callbackdata.Indices(1),5:6) = {''};
        elseif callbackdata.Indices(2) == 5
            a.Material3{callbackdata.Indices(1)} = getfield(a.Materials.Cubic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3:4) = {''};
            hObject.Data(callbackdata.Indices(1),6) = {''};
        elseif callbackdata.Indices(2) == 6
            a.Material3{callbackdata.Indices(1)} = getfield(a.Materials.Isotropic,callbackdata.EditData);
            hObject.Data(callbackdata.Indices(1),3:5) = {''};
        elseif callbackdata.Indices(2) == 7
            if  length(a.LayerOrientations3) > 1 || length(a.LayerThicknesses3) > 1 || length(a.Material3) > 1
                if  length(a.LayerOrientations3) >= callbackdata.Indices(1)
                    a.LayerOrientations3(callbackdata.Indices(1)) = [];
                end
                if  length(a.LayerThicknesses3) >= callbackdata.Indices(1)
                    a.LayerThicknesses3(callbackdata.Indices(1)) = [];
                end
                if  length(a.Material3) >= callbackdata.Indices(1)
                    a.Material3(callbackdata.Indices(1)) = [];
                    a.Material3(cellfun(@isempty,a.Material3)) = [];
                end
                if  length(a.MaterialNames3) >= callbackdata.Indices(1)
                    a.MaterialNames3(callbackdata.Indices(1)) = [];
                end
                if  length(a.MaterialClasses3) >= callbackdata.Indices(1)
                    a.MaterialClasses3(callbackdata.Indices(1)) = [];
                end
                hObject.Data(callbackdata.Indices(1),:) = [];
                hObject.Data = vertcat(hObject.Data,cell(1,6));
            else
                errordlg('Do not delete the only layer!','Error');
                return
            end
            if  a.UniformLayerThickness3 == 1
                a.LayerThicknesses3(1:length(a.LayerOrientations3)) = sum(a.LayerThicknesses3)/(length(a.LayerOrientations3));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),2) = {a.LayerThicknesses3(1)};
            end
            if  a.Hybrid6 == 0
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(a.Material3{1});
                hObject.Data(1:length(a.LayerOrientations3),a.MaterialType6+2) = {a.Material3{1}.Name};
                if  a.MaterialType6 == 1
                    hObject.Data(1:length(a.LayerOrientations3),4:6) = {''};
                elseif a.MaterialType6 == 2
                    hObject.Data(1:length(a.LayerOrientations3),3) = {''};
                    hObject.Data(1:length(a.LayerOrientations3),5:6) = {''};
                elseif a.MaterialType6 == 3
                    hObject.Data(1:length(a.LayerOrientations3),3:4) = {''};
                    hObject.Data(1:length(a.LayerOrientations3),6) = {''}; 
                elseif a.MaterialType6 == 4
                    hObject.Data(1:length(a.LayerOrientations3),3:5) = {''};
                end
            end
        end
        try
            for i = 1:length(a.Material3)
                a.MaterialNames3{i} = a.Material3{i}.Name;
                a.MaterialClasses3{i} = a.Material3{i}.Class;
            end
        catch
            
        end
    end
    function OKUI6_Callback(~,~)
        try
            for i = 1:length(a.Material3)
                a.MaterialNames3{i} = a.Material3{i}.Name;
                a.MaterialClasses3{i} = a.Material3{i}.Class;
            end
        catch
            errordlg('A layer in between misses a material!','Error');
            return
        end
        if  length(a.LayerOrientations3) ~= length(a.LayerThicknesses3) && a.UniformLayerThickness3 == 0
            errordlg('If the laminate does not have uniform layer thicknesses, you must enter every layer''s thickness in the ''d (mm)'' column! All columns must have the same number of entries.','Error');
            return
        elseif length(a.LayerOrientations3) ~= length(a.Material3) && a.Hybrid6 == 1
            errordlg('If you set up a hybrid layup, you must assign a material to every layer defined in the ''Phi'' column and vice versa!','Error');
            return
        elseif all(strcmp(a.MaterialNames3(1),a.MaterialNames3)) && a.Hybrid6 == 1
            if  length(a.LayerOrientations3) > 1
                errordlg('It is no hybrid if all layers are the same material!','Error');
            else
                errordlg('It is no hybrid if you have only one layer!','Error');
            end
            a.Hybrid6 = 0;
            a.HybridUI6.Value = 0;
            
            a.MaterialTypeUI6.Enable = 'on';
            a.MaterialTypeUI6Enable = 'on';
            a.MaterialUI6.Enable = 'on';
            a.MaterialUI6Enable = 'on';
            if  a.UniformLayerThickness3 == 1
                a.TablesUI6.ColumnEditable = [true false false false false false true];
                a.TablesUI6ColumnEditable = [true false false false false false true];
            else
                a.TablesUI6.ColumnEditable = [true true false false false false true];
                a.TablesUI6ColumnEditable = [true true false false false false true];
            end
            
            if  a.MaterialType6 == 1
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Orthotropic,cell2mat(a.MaterialUI6.String(a.MaterialUI6.Value))));
                a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Orthotropic),a.Material3{1}.Name));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),4:6) = {''};
            elseif a.MaterialType6 == 2
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.TransverselyIsotropic,cell2mat(a.MaterialUI6.String(a.MaterialUI6.Value))));
                a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.TransverselyIsotropic),a.Material3{1}.Name));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),3) = {''};
                a.TablesUI6.Data(1:length(a.LayerOrientations3),5:6) = {''};
            elseif a.MaterialType6 == 3
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Cubic,cell2mat(a.MaterialUI6.String(a.MaterialUI6.Value))));
                a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Cubic),a.Material3{1}.Name));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),3:4) = {''};
                a.TablesUI6.Data(1:length(a.LayerOrientations3),6) = {''};
            elseif a.MaterialType6 == 4
                [a.Material3{1:length(a.LayerOrientations3)}] = deal(getfield(a.Materials.Isotropic,cell2mat(a.MaterialUI6.String(a.MaterialUI6.Value))));
                a.MaterialUI6Value = find(strcmp(fieldnames(a.Materials.Isotropic),a.Material3{1}.Name));
                a.TablesUI6.Data(1:length(a.LayerOrientations3),3:5) = {''};
            end
            a.TablesUI6.Data(1:length(a.LayerOrientations3),a.MaterialType6+2) = {a.Material3{1}.Name};
            
            return
        end        
        a.UnitCell3 = a.TablesUI6.Data;
        if  a.Hybrid6 == 0
            a.MaterialNameUI6.String = a.Material3{1}.Name;
        elseif a.Hybrid6 == 1
            a.MaterialNameUI6.String = 'Hybrid';
        end
        a.LayupUI6.String = ['[',char(join(split(num2str(a.LayerOrientations3)),'/')),']'];
        close(f4)
        String = 'Unit cell:';
        for i = 1:length(a.MaterialNames3)
            String = append(String,newline,num2str(i),': ',a.MaterialNames3{i},' (',a.MaterialClasses3{i},')');
            if  ~isreal(a.Material3{i}.C)
                a.Material3{i}.C = real(a.Material3{i}.C);
            end
        end
        disp([String,newline,'-----------------------------------'])
        
        a.CLaminate = Computer_LaminateStiffness(a.Material3,a.LayerOrientations3,a.LayerThicknesses3,a.PropagationAngle3);
        a.CLP = Computer_LaminateStiffnessPolar(a.Material3,a.LayerOrientations3,a.LayerThicknesses3,1);
        
        a.C11UI6.String = a.CLaminate(1,1)/1e9;
        a.C12UI6.String = a.CLaminate(1,2)/1e9;
        a.C13UI6.String = a.CLaminate(1,3)/1e9;
        a.C16UI6.String = a.CLaminate(1,6)/1e9;
        a.C22UI6.String = a.CLaminate(2,2)/1e9;
        a.C23UI6.String = a.CLaminate(2,3)/1e9;
        a.C26UI6.String = a.CLaminate(2,6)/1e9;
        a.C33UI6.String = a.CLaminate(3,3)/1e9;
        a.C36UI6.String = a.CLaminate(3,6)/1e9;
        a.C44UI6.String = a.CLaminate(4,4)/1e9;
        a.C45UI6.String = a.CLaminate(4,5)/1e9;
        a.C55UI6.String = a.CLaminate(5,5)/1e9;
        a.C66UI6.String = a.CLaminate(6,6)/1e9;
        
        a.h1 = LaminateStiffness_Internal(a);
    end
    function CancelUI6_Callback(~,~)
        close(f4)
    end
end
function CallbackUI1(source,eventdata) % Tab1_isotropic
    a = CallbackModule_Isotropic(source,eventdata,a,Tab1,Tab3); 
end
function CallbackUI2(source,eventdata) % Tab2_anisotropic
    a = CallbackModule_Anisotropic(source,eventdata,a,Tab2,Tab3);
end
function CallbackUI3(source,eventdata) % Tab3_signal simulator
    a = CallbackModule_SignalSimulator(source,eventdata,a,Tab3);
end
function CallbackUI4(source,eventdata) % Tab4_polar diagrams
    a = CallbackModule_PolarDiagrams(source,eventdata,a,Tab4);
end
function CallbackUI8(source,eventdata) % Tab8_bulk waves
    a = CallbackModule_BulkWaves(source,eventdata,a);
end
function CallbackUI6(source,eventdata) % Tab6_laminate stiffness
    a = CallbackModule_LaminateStiffness(source,eventdata,a);
end
function CallbackUI5(source,eventdata) % Tab5_material editor
    a = CallbackModule_MaterialEditor(source,eventdata,a);
end
function CallbackUI7(source,eventdata) % Tab7_advanced
    a = CallbackModule_Advanced(source,eventdata,a);
end
end