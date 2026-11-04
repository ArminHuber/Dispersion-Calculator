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
function DispersionCalculator
Version = '3.3';
Date = 'November 4, 2026';
Mode = 1; % 1: DC executed in MATLAB 2: DC executed as stand-alone; this switches the default directory for saving plots and data (see lines 55ff)

%#ok<*STRNU>
%#ok<*AGROW>
%#ok<*GVMIS>

warning('off','MATLAB:nearlySingularMatrix') % turn warnings off
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:illConditionedMatrix')

try % make "DC@GitHub" blue to indicate that a new version is available
    X = webread('https://github.com/ArminHuber/Dispersion-Calculator');
    z = strfind(X,'Dispersion Calculator');
    F = extractBetween(X(z(1):z(1)+30),'Calculator v',' ');
    if  strcmp(Version,F{1})
        UpdateColor = [.1294 .1294 .1294];
    else
        UpdateColor = 'b';
    end
catch
    UpdateColor = [.1294 .1294 .1294];
end

Dir = which('DispersionCalculator'); % get the full path and file name of "DispersionCalculator.m"
a.MaterialListDirectory = [Dir(1:length(Dir)-23),filesep,'Materials']; % get the path where the materials are stored
if  Mode == 1 % default directory if DC is executed in MATLAB
    a.Directory = Dir(1:length(Dir)-24-length(extractAfter(Dir(1:length(Dir)-23),asManyOfPattern(wildcardPattern+filesep))));
elseif Mode == 2 % default directory if DC is executed as stand-alone
    a.Directory = [extractBefore(Dir,filesep),filesep];
end

f1 = figure('Icon',which('DC_Logo16.png'),'NumberTitle','off','Name',['Dispersion Calculator v',Version],'Visible','off','MenuBar','none','Position',[0 0 1198 800],'CloseRequestFcn',@CloseRequest); % generate the main GUI
m1 = uimenu(f1,'Text','File'); % generate the menu bar
uimenu(m1,'Text','New project','MenuSelectedFcn',@New_Callback)
uimenu(m1,'Text','Open project','MenuSelectedFcn',@Open_Callback)
uimenu(m1,'Text','Save project','MenuSelectedFcn',@Save_Callback)
uimenu(m1,'Text','Quit','Separator','on','MenuSelectedFcn',@CloseRequest)
m2 = uimenu(f1,'Text','Materials');
uimenu(m2,'Text','Import','MenuSelectedFcn',@Import_Callback)
uimenu(m2,'Text','Export','MenuSelectedFcn',@Export_Callback)
uimenu(m2,'Text','Open directory','MenuSelectedFcn',@OpenDirectory_Callback)
m3 = uimenu(f1,'Text','Multicore');
uimenu(m3,'Text','Enable','MenuSelectedFcn',@Enable_Callback)
uimenu(m3,'Text','Disable','MenuSelectedFcn',@Disable_Callback)
uimenu(f1,'Text','Help','MenuSelectedFcn',@Help_Callback)
uimenu(f1,'Text','DC@GitHub','MenuSelectedFcn',@Homepage_Callback,'ForegroundColor',UpdateColor)
uimenu(f1,'Text','About','MenuSelectedFcn',@About_Callback)

tgroup = uitabgroup('Parent',f1); % generate the tab bar
Tab1 = uitab('Parent',tgroup,'Title','Isotropic');
Tab2 = uitab('Parent',tgroup,'Title','Anisotropic');
Tab3 = uitab('Parent',tgroup,'Title','Signal simulator');
Tab4 = uitab('Parent',tgroup,'Title','Polar diagrams');
Tab8 = uitab('Parent',tgroup,'Title','Bulk waves');
Tab6 = uitab('Parent',tgroup,'Title','Laminate stiffness');
Tab5 = uitab('Parent',tgroup,'Title','Material editor');
Tab7 = uitab('Parent',tgroup,'Title','Advanced','Backgroundcolor','white');

try
    a.Materials = MaterialList; % load materials from "MaterialList_Fluid.txt", "MaterialList_Isotropic.txt", "MaterialList_TransverselyIsotropic.txt", and "MaterialList_Orthotropic.txt";
catch ME
    switch ME.identifier
    case 'MATLAB:UndefinedFunction'
        errordlg('The DC files are not on the MATLAB path. Right-click on the DC folder -> Add to Path -> Selected Folders and Subfolders.','Files are not on the MATLAB path')        
    case 'MATLAB:setfield'
        try
            errordlg(['Invalid material name ',extractAfter(ME.message,':'),newline,newline,['Open the corresponding material list and fix the invalid material name. Valid names must start with a letter and may contain only letters, numbers, and the underscore character. The length may not exceed ',num2str(namelengthmax),' characters.'],newline,newline,'Notice: The recommended way to add, edit, and delete materials is to use the ''Material editor'' inside DC. Do not edit the material list txt-files directly.'],'Invalid material name')
            pause(2)
            winopen(a.MaterialListDirectory)
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Error')
        end
    otherwise
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Error')
    end
    return
end

global Stop % global variable used to stop a calculation
Stop = 0;

a.FieldNames_Orthotropic = fieldnames(a.Materials.Orthotropic);
a.FieldNames_TransverselyIsotropic = fieldnames(a.Materials.TransverselyIsotropic);
a.FieldNames_Cubic = fieldnames(a.Materials.Cubic);
a.FieldNames_Isotropic = fieldnames(a.Materials.Isotropic);
a.FieldNames_Fluid = fieldnames(a.Materials.Fluid);

if  1 % advanced settings
    a.Multithreading = 0;

    % 2-D tracing settings
    a.SearchWidth = [1.5 -1.5];
    a.SearchAreaSections = 5;
    a.SearchAreaExtensions = 3;
    
    a.CriticalSlope = 1000; % (m/s) if the slope of the phase velocity curve exceeds this value, a smaller frequency step is used
    a.GearSlopeRanges = [0 1 2 4 Inf]*a.CriticalSlope;
    a.GearFrequencyResolution = [1 2 4];
    
    % isotropic -----------------------------------------------------------
    a.YRange1 = 4; % determines how high the dispersion diagram is by default
    a.XRange1 = 2; % determines how broad the dispersion diagram is by default
    a.XSamples1 = 1000; % determines how many samples a dispersion curve consists of by default
    a.Steps1 = 10000; % determines how many samples the frequency sweep for the higher order modes contains

    a.FrequencyRangeStart1 = 1e-3; % start computation @ x kHz
    
    a.PhaseVelocitySections1 = 6;
    a.FrequencySections1 = 6;

    a.LambPhaseVelocitySweepRange11 = 2; % phase velocity search interval for Lamb waves for negative curvature; a higher value decreases the risk of missing a solution; a larger one might require a higher number of PhaseVelocitySections which increases processing time    
    a.LambPhaseVelocitySweepRange21 = 2; % for postive curvature

    a.PhaseVelocityStep1 = 100; % phase velocity steps for the frequency sweeps at high phase velocity to complete the modes (m/s) 
    a.FrequencyOffset1 = 20; % offset to low frequency for the frequency sweeps at high phase velocity for the modes (kHz/mm); its absolute must be large enough to find any solution but small enough to not reach a lower mode on the left; a larger SStep requires also a larger SOffset
    
    a.MissingSamples1 = 3; % if more than that number of sample points are missing, the tracing of a dispersion curve stops
    a.BelowCutoffWidth1 = 50;  % (kHz/mm); if we exceed the allowed scanning width (in kHz/mm) below the cut-off frequency of a damped mode without finding it, the damped search stops, and only the nondamped tracing continues to have that data for the next mode; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one

    % anisotropic ----------------------------------------------
    a.YRange2 = 4; % determines how high the dispersion diagram is by default (for isotropic single layers)
    a.XRange2 = 4; % determines how broad the dispersion diagram is by default (MHz*mm)
    a.XSamples2 = 400; % determines how many samples a dispersion curve consists of by default
    a.Steps2 = 4000; % determines how many samples the frequency sweep for the higher order modes contains
    
    a.FrequencyRangeStart2 = 1e-3; % start computation @ x kHz 
    
    a.PhaseVelocitySections2 = 5;
    a.FrequencySections2 = 6;
    
    a.LambPhaseVelocitySweepRange12 = 10; % phase velocity search interval for Lamb waves for negative curvature; a higher value decreases the risk of missing a solution; a larger one might require a higher number of PhaseVelocitySections which increases processing time    
    a.LambPhaseVelocitySweepRange22 = 5; % for postive curvature
    
    a.PhaseVelocityStep2 = 100; % phase velocity steps for the frequency sweeps at high phase velocity to complete the modes (m/s) 
    a.FrequencyOffset2 = 20; % offset to low frequency for the frequency sweeps at high phase velocity for the modes (kHz/mm); its absolute must be large enough to find any solution but small enough to not reach a lower mode on the left; a larger SStep requires also a larger SOffset
    
    a.MissingSamples2 = 3; % if more than that number of sample points are missing, the tracing of a dispersion curve stops
    a.BelowCutoffWidth2 = 50;  % (kHz/mm); if we exceed the allowed scanning width (in kHz/mm) below the cut-off frequency of a damped mode without finding it, the damped search stops, and only the nondamped tracing continues to have that data for the next mode; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one

    % polar diagrams ----------------------------------------------
    a.XRange_Polar = .5; % determines how broad the dispersion diagram is by default (MHz*mm)
    a.XSamples_Polar = 50; % determines how many samples a dispersion curve consists of by default 
    a.MissingSamples_Polar = 3; % if more than that number of sample points are missing, the tracing of a dispersion curve stops
end
if  1 % default settings
    if  1 % Tab1_isotropic
        a.Directory1 = a.Directory;
        
        a.Search1 = 0; % for isotropic tab: variable in order to check whether the higher order mode detection has been performed before their calculation       

        a.HSLamb1 = [];
        a.HALamb1 = [];
        a.HBLamb1 = [];
        a.HSShear1 = [];
        a.HAShear1 = [];
        a.HL1 = [];
        a.HF1 = [];
        a.HT1 = [];
        a.HCLamb1 = [];
        a.HCShear1 = [];

        a.Geometry1 = 'Plate';
        a.FluidLoading1 = 0;
        a.ToggleUpperFluid1 = 0;
        a.ToggleLowerFluid1 = 0;
        a.UpperFluid1 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        a.LowerFluid1 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        a.Sink1 = 0;
        a.Material1 = getfield(a.Materials.Isotropic,a.FieldNames_Isotropic{1});
        if  isreal(a.Material1.C)
            a.Viscoelastic1 = 0;
        else
            a.Viscoelastic1 = 1;
        end
        a.Thickness1 = 1;
        a.ThicknessInner1 = .8*a.Thickness1;
        a.Symmetric1 = 1;

        a.Fix1 = 0;
        a.PhaseVelocityLimit1 = round(a.Material1.PlateVelocity*a.YRange1,2,'significant');
        a.Accuracy1 = 1e-6;
        a.FrequencyLimit1 = round(a.Material1.PlateVelocity*a.XRange1/a.Thickness1,2,'significant');
        a.FrequencyResolution1 = a.FrequencyLimit1/a.XSamples1;

        a.Force2DTracing1 = 0;
        a.PhaseVelocityLimit21 = 1000*1e3;
        a.AttenuationLimit1 = 1;
        a.Sweeps1 = 1;
        a.SweepSections1 = 256;

        a.HigherOrderModes1 = 1;
            a.SymmetricModes1 = 1;
            a.AntisymmetricModes1 = 1;
            a.LambModes1 = 1;
            a.ShearHorizontalModes1 = 1;

            a.TorsionalModes1 = 1;
            a.LongitudinalModes1 = 1;
            a.FlexuralModes1 = 1;
            a.FlexuralModeOrders1 = 1;

        a.LineColors1 = [0 0 1 % blue (1) for flexural modes
            .13 .55 .13 % green (2)
            0 0 0 % black (3)
            1 0 1 % magenta (4)
            0 1 1 % cyan (5)
            1 .7 0 % orange (6)
            .55 .27 .13 % brown (7)
            .5 0 1 % violet (8)
            .5 .5 .5 % gray (9)
            1 0 0]; % red (10)

        a.Step1 = a.FrequencyLimit1/a.Steps1;

        a.SamplesX3ce1 = 50;

        a.Quantity11 = 1;
        a.BulkVelocities1 = 0;
        a.Distance1 = 1e2;
        
        a.Couplant1 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});

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
        a.Plane1 = 'r-z';
        a.Frequency21 = a.FrequencyLimit1;
        a.Length1 = 2.5;
        a.Samples21 = 80;
        a.Samples31 = .5*a.Samples21;
        a.Gain1 = 1;
        a.GridLine1 = 2;
        a.Undistorted1 = 0;
        a.HalfspacesNumber21 = 1;
        a.Halfspaces21 = 0;

        a.CycleDuration1 = 1.5;
        a.FrameRate1 = 30;
        a.MovieQuality1 = 75;
        a.Animate1 = 0;

        a.ExportPlots1 = 0;
        a.PDF1 = 1;
        a.PNG1 = 0;
        a.PNGresolution1 = 150;
        a.XAxisMode21 = 1;
        a.Arrange1 = 1;
        a.DispersionCurves1 = 1;
        a.ThroughThickness1 = 0;
        a.FileName1 = 'File name';

        a.Title1 = 1;
        a.LegendLocation1 = 'out';
        a.BoxLineWidth1 = .5;
        a.LineWidth1 = 1;
        a.SColor1 = [1 0 0];
        a.AColor1 = [0 0 1];
        a.BColor1 = [.5 0 1];
        a.TitleFontSize1 = 30;
        a.AxesLabelFontSize1 = 30;
        a.AxesTickFontSize1 = 24;
        a.LegendFontSize1 = 24;

        a.MaterialName4A = a.Material1.Name;
    end
    if  1 % Tab2_anisotropic
        a.Directory2 = a.Directory;
        
        a.Search2 = 0; % for anisotropic tab
        
        a.HSLamb2 = [];
        a.HALamb2 = [];
        a.HBLamb2 = [];
        a.HSShear2 = [];
        a.HAShear2 = [];
        a.HBShear2 = [];
        a.HL2 = [];
        a.HF2 = [];
        a.HT2 = [];
        a.HCLamb2 = [];
        a.HCShear2 = [];

        a.Geometry2 = 'Plate';
        a.Radius = 0;
        a.InnerRadius = 1;
        a.Sink2 = 0;
        a.FluidLoading2 = 0;
        a.ToggleUpperFluid2 = 0;
        a.ToggleLowerFluid2 = 0;
        a.UpperFluid2 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        a.LowerFluid2 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        a.MaterialType2 = 1;
        a.Material2 = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})};
        if  isreal(a.Material2{1}.C)
            a.Viscoelastic2 = 0;
        else
            a.Viscoelastic2 = 1;
        end
        a.IsotropicSinglelayer = 0;
        a.IsotropicMultilayerCylinder = 0;
        a.MaterialNames = {a.Material2{1}.Name};
        a.MaterialClasses = {a.Material2{1}.Class};
        a.Hybrid = 0;
        a.UniformLayerThickness = 1;
        a.LayupThickness = 1;
        a.PlateThickness = a.LayupThickness;
        a.SuperLayers = 1;
        a.Pattern = [];
        a.SuperLayerSize = 1;
        a.SymmetricSystem = 0;
        a.Symmetric2 = 1;
        a.Layers2 = 1;
        a.UnitCell = cell(99,7);
        a.UnitCell{1} = 0;
        a.UnitCell{1,2} = 1;
        a.UnitCell{1,3} = a.FieldNames_Orthotropic{1};
        a.LayerOrientations = 0;
        a.LayerThicknesses = 1;
        a.LayupString2 = '[0]';
        a.LayupStringTex2 = '[0]';
        a.EffectiveLayupString2 = '[0]';

        a.Fix2 = 0;
        a.PropagationAngle = 0;
        a.PhaseVelocityLimit2 = 20e3;
        a.Accuracy2 = 1e-6;
        a.FrequencyLimit2 = round(1e3*a.XRange2/a.LayupThickness,2,'significant');
        a.FrequencyResolution2 = a.FrequencyLimit2/a.XSamples2;
        a.MatrixMethods2 = 1;

        a.Force2DTracing2 = 0;
        a.PhaseVelocityLimit22 = 1000*1e3;
        a.AttenuationLimit2 = 1;
        a.Sweeps2 = 1;
        a.SweepSections2 = 256;

        a.HigherOrderModes2 = 1;
            a.SymmetricModes2 = 1;
            a.AntisymmetricModes2 = 1;
            a.LambModes2 = 1;
            a.ShearHorizontalModes2 = 1;
    
            a.TorsionalModes2 = 1;
            a.LongitudinalModes2 = 1;
            a.FlexuralModes2 = 1;
            a.FlexuralModeOrders2 = 1;

        a.LineColors2 = [0 0 1 % blue (1) for flexural modes
            .13 .55 .13 % green (2)
            0 0 0 % black (3)
            1 0 1 % magenta (4)
            0 1 1 % cyan (5)
            1 .7 0 % orange (6)
            .55 .27 .13 % brown (7)
            .5 0 1 % violet (8)
            .5 .5 .5 % gray (9)
            1 0 0]; % red (10)        

        a.Step2 = a.FrequencyLimit2/a.Steps2;

        a.SamplesX3ce2 = 50;

        a.Quantity12 = 1;
        a.BulkVelocities2 = 0;
        a.Distance2 = 1e2;
        
        a.Couplant2 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        
        a.XAxisMode2 = 1;
        a.XAxis2 = [0 a.FrequencyLimit2];
        a.YAxis2 = [0 a.PhaseVelocityLimit2/1e3];

        a.Quantity22 = 1;
        a.Frequency12 = a.FrequencyLimit2;
        a.Samples12_0 = 500;
        a.Samples12 = a.Samples12_0;
        a.HalfspacesNumber12 = 1;
        a.Halfspaces12 = 0;
        a.Phase2 = 0;
        a.Plot2 = [1 1 1 1 1 1]; % displacement components - - 3 2 1 - / stress components 11 22 33 23 13 12

        a.ModeShapeSettingChanged2 = 1;
        a.Plane2 = 'r-z';
        a.Frequency22 = a.FrequencyLimit2;
        a.Length2 = 2.5;
        a.Samples22 = 80;
        a.Samples32 = .5*a.Samples22;
        a.Gain2 = 1;
        a.GridLine2 = 2;
        a.Undistorted2 = 0;
        a.HalfspacesNumber22 = 1;
        a.Halfspaces22 = 0;

        a.CycleDuration2 = 1.5;
        a.FrameRate2 = 30;
        a.MovieQuality2 = 75;
        a.Animate2 = 0;

        a.ExportPlots2 = 0;
        a.PDF2 = 1;
        a.PNG2 = 0;
        a.PNGresolution2 = 150;
        a.XAxisMode22 = 1;
        a.Arrange2 = 1;
        a.DispersionCurves2 = 1;
        a.ThroughThickness2 = 0;
        a.FileName2 = 'File name';

        a.Title2 = 2;
        a.LegendLocation2 = 'out';
        a.BoxLineWidth2 = .5;
        a.LineWidth2 = 1;
        a.SColor2 = [1 0 0];
        a.AColor2 = [0 0 1];
        a.BColor2 = [.5 0 1];
        a.TitleFontSize2 = 30;
        a.AxesLabelFontSize2 = 30;
        a.AxesTickFontSize2 = 24;
        a.LegendFontSize2 = 24;

        a.MaterialName4B = a.Material2{1}.Name;
    end
    if  1 % Tab3_signal simulator
        a.Directory3 = a.Directory;
        
        a.Chirp3 = 0;
        a.ChirpType3 = 1;
        a.Cycles3 = 10;
        a.SamplesPerCycle3 = 20;
        a.Window3 = 1;
        a.WidthFactorGauss3 = 3.16;
        a.WidthFactorTukey3 = .5;
        a.WidthFactor3 = a.WidthFactorGauss3;
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
        
        a.LineColors3 = [0 0 1; % blue (1)
            1 0 0; % red (2)
            .13 .55 .13; % green (3)
            0 0 0; % black (4)
            1 0 1; % magenta (5)
            0 1 1; % cyan (6)
            1 .7 0; % orange (7)
            .55 .27 .13; % brown (8)
            .5 0 1; % violet (9)
            .5 .5 .5]; % gray (10)

        a.ExportPlots3 = 0;
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
        a.Material_Polar = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})}; 
        a.MaterialNames_Polar = {a.Material_Polar{1}.Name};
        a.MaterialClasses_Polar = {a.Material_Polar{1}.Class};
        a.Hybrid_Polar = 0;
        a.UniformLayerThickness_Polar = 1;
        a.LayupThickness_Polar = 1;
        a.PlateThickness_Polar = a.LayupThickness_Polar;
        a.SuperLayers_Polar = 1;
        a.Pattern_Polar = [];
        a.SuperLayerSize_Polar = 1;
        a.SymmetricSystem_Polar = 0;
        a.Layers_Polar = 1;
        a.UnitCell_Polar = cell(99,7);
        a.UnitCell_Polar{1} = 0;
        a.UnitCell_Polar{1,2} = 1;
        a.UnitCell_Polar{1,3} = a.FieldNames_Orthotropic{1};
        a.LayerOrientations_Polar = 0;
        a.LayerThicknesses_Polar = 1;
        a.LayupString_Polar = '[0]';
        a.LayupStringTex_Polar = '[0]';

        a.FrequencyLimit_Polar = 500;
        a.FrequencyResolution_Polar = a.FrequencyLimit_Polar/a.XSamples_Polar;
        a.Accuracy_Polar = 1e-4;
        a.PropagationAngleMode_Polar = 1;
        a.PropagationAngleStep_Polar = 3;
        a.PhaseVelocitySections_Polar = 5;
        a.MatrixMethods_Polar = 1;

        a.S0_Polar = 1;
        a.SH0_Polar = 1;
        a.A0_Polar = 1;

        a.SamplesX3ce_Polar = 50;

        a.Quantity_Polar = 1;
        a.BulkVelocities_Polar = 0;
        a.Distance_Polar = 100;
        
        a.Couplant_Polar = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        
        a.Frequency_Polar = 500;

        a.ExportPlots_Polar = 0;
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
        a.Material_Bulk = getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1});
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
        
        a.Couplant_Bulk = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        
        a.SolidType_Bulk = 1;
        a.Solid_Bulk = getfield(a.Materials.Isotropic,a.FieldNames_Isotropic{1});        
        a.Phi_Bulk2 = 0;
        a.Theta_Bulk2 = 2;
        a.ThetaStep_Bulk2 = .1;
        a.ViewPhi_Bulk2 = -37.5;
        a.ViewTheta_Bulk2 = 30;

        a.ExportPlots_Bulk = 0;
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
        a.Material6 = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})};
        a.MaterialNames6 = {a.Material6{1}.Name};
        a.MaterialClasses6 = {a.Material6{1}.Class};
        a.Hybrid6 = 0;
        a.UniformLayerThickness6 = 1;
        a.LayupThickness6 = 1;
        a.SuperLayerSize6 = 1;
        a.UnitCell6 = cell(99,7);
        a.UnitCell6{1} = 0;
        a.UnitCell6{1,2} = 1;
        a.UnitCell6{1,3} = a.FieldNames_Orthotropic{1};
        a.LayerOrientations6 = 0;
        a.LayerThicknesses6 = 1;
        a.LayupString6 = '[0]';

        a.PropagationAngle6 = 0;

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
        a.p1UI1 = uipanel('Parent',Tab1,'Title','Specimen','Units','pixels','Position',[10 520 225 235],'FontSize',10);
        uicontrol('Parent',a.p1UI1,'Style','text','HorizontalAlignment','left','String','Geometry','Position',[10 195 49 13])
        uicontrol('Parent',a.p1UI1,'Style','text','HorizontalAlignment','left','String','Material','Position',[10 165 39 13])
        a.ThicknessTextUI1 = uicontrol('Parent',a.p1UI1,'Style','text','HorizontalAlignment','left','String','Thickness (mm)','Position',[10 135 100 13]); % thickness 78, Outer diameter
        uicontrol('Parent',a.p1UI1,'Style','text','HorizontalAlignment','left','String','Inner diameter (mm)','Position',[10 105 97 13])
        a.UpperFluidTextUI1 = uicontrol('Parent',a.p1UI1,'Style','text','HorizontalAlignment','left','String','Upper fluid','Position',[10 75 54 13]);
        a.LowerFluidTextUI1 = uicontrol('Parent',a.p1UI1,'Style','text','HorizontalAlignment','left','String','Lower fluid','Position',[10 45 57 13]);
        uicontrol('Parent',a.p1UI1,'Style','text','HorizontalAlignment','left','String','Sink at center','Position',[10 15 68 13]);
        a.GeometryUI1 = uicontrol('Parent',a.p1UI1,'Style','popupmenu','String',{'Plate','Rod','Pipe','Circumferential'},'Tooltip','Select a geometry.','Position',[80 190 130 23],'Callback',@CallbackUI1,'Tag','6');
        a.MaterialUI1 = uicontrol('Parent',a.p1UI1,'Style','popupmenu','String',a.FieldNames_Isotropic,'Tooltip','Select a material.','Position',[80 160 130 23],'Callback',@CallbackUI1,'Tag','1');
        a.ThicknessUI1 = uicontrol('Parent',a.p1UI1,'Style','edit','String',a.Thickness1,'Tooltip','Enter the plate''s thickness/outer pipe/rod diameter.','Position',[160 130 50 23],'Callback',@CallbackUI1,'Tag','2');
        a.ThicknessInnerUI1 = uicontrol('Parent',a.p1UI1,'Style','edit','String',a.ThicknessInner1,'Tooltip','Enter the inner pipe diameter.','Position',[160 100 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','90');
        a.ToggleUpperFluidUI1 = uicontrol('Parent',a.p1UI1,'Style','checkbox','Value',a.ToggleUpperFluid1,'Tooltip','Toggle upper fluid. If unchecked, the upper half-space will be vacuum.','Position',[80 70 50 23],'Callback',@CallbackUI1,'Tag','7');
        a.ToggleLowerFluidUI1 = uicontrol('Parent',a.p1UI1,'Style','checkbox','Value',a.ToggleLowerFluid1,'Tooltip','Toggle lower fluid. If unchecked, the lower half-space will be vacuum.','Position',[80 40 50 23],'Callback',@CallbackUI1,'Tag','86');
        a.SelectUpperFluidUI1 = uicontrol('Parent',a.p1UI1,'Style','popupmenu','String',a.FieldNames_Fluid,'Tooltip','Select the fluid in the upper half-space.','Position',[110 70 100 23],'Enable','off','Callback',@CallbackUI1,'Tag','87');
        a.SelectLowerFluidUI1 = uicontrol('Parent',a.p1UI1,'Style','popupmenu','String',a.FieldNames_Fluid,'Tooltip','Select the fluid in the lower half-space.','Position',[110 40 100 23],'Enable','off','Callback',@CallbackUI1,'Tag','88');
        a.SinkUI1 = uicontrol('Parent',a.p1UI1,'Style','checkbox','Value',a.Sink1,'Tooltip','Put a sink at the center of the pipe. This absorbs the bulk waves travelling from one side of the inner pipe wall through the fluid- filling to the other side, and therefore suppresses the fluid modes.','Position',[80 10 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','89');
        
        a.p2UI1 = uipanel('Parent',Tab1,'Title','Computational settings','Units','pixels','Position',[10 370 225 145],'FontSize',10);
        uicontrol('Parent',Tab1,'Style','togglebutton','String','Fix','Tooltip','Fix the computational settings such that they are not adjusted automatically.','Position',[212 495 23 23],'Callback',@CallbackUI1,'Tag','61');
        uicontrol('Parent',a.p2UI1,'Style','text','HorizontalAlignment','left','String','Frequency limit (kHz)','Position',[10 105 103 13]);
        uicontrol('Parent',a.p2UI1,'Style','text','HorizontalAlignment','left','String','Frequency step (kHz)','Position',[10 75 107 13]);
        uicontrol('Parent',a.p2UI1,'Style','text','HorizontalAlignment','left','String','Phase velocity limit (m/ms)','Position',[10 45 128 13]);
        uicontrol('Parent',a.p2UI1,'Style','text','HorizontalAlignment','left','String','Phase velocity accuracy (m/s)','Position',[10 15 149 13]);
        a.FrequencyLimitUI1 = uicontrol('Parent',a.p2UI1,'Style','edit','String',a.Frequency21,'Tooltip','Enter the frequency (X-axis in the dispersion diagram) up to which the dispersion curves shall be traced.','Position',[160 100 50 23],'Callback',@CallbackUI1,'Tag','4');
        a.FrequencyResolutionUI1 = uicontrol('Parent',a.p2UI1,'Style','edit','String',a.FrequencyResolution1,'Tooltip','Enter the frequency step.','Position',[160 70 50 23],'Callback',@CallbackUI1,'Tag','5');        
        a.PhaseVelocityLimitUI1 = uicontrol('Parent',a.p2UI1,'Style','edit','String',a.PhaseVelocityLimit1/1e3,'Tooltip','Enter the phase velocity (Y-axis in the dispersion diagram) up to which the dispersion curves shall be traced in case there is no attenuation and 2-D tracing is not forced.','Position',[160 40 50 23],'Callback',@CallbackUI1,'Tag','3');
        a.AccuracyUI1 = uicontrol('Parent',a.p2UI1,'Style','edit','String',a.Accuracy1,'Tooltip','Bisections are performed until the modal solutions have reached this accuracy.','Position',[160 10 50 23],'Callback',@CallbackUI1,'Tag','85');

        a.p14UI1 = uipanel('Parent',Tab1,'Title','2-D tracing settings','Units','pixels','Position',[10 190 225 175],'FontSize',10);
        uicontrol('Parent',a.p14UI1,'Style','text','HorizontalAlignment','left','String','Force 2-D tracing','Position',[10 135 145 13]);
        uicontrol('Parent',a.p14UI1,'Style','text','HorizontalAlignment','left','String','Phase velocity limit 2 (m/ms)','Position',[10 105 145 13]);
        uicontrol('Parent',a.p14UI1,'Style','text','HorizontalAlignment','left','String','Attenuation limit (k_im/k_re)','Position',[10 75 145 13]);
        uicontrol('Parent',a.p14UI1,'Style','text','HorizontalAlignment','left','String','Sweeps','Position',[10 45 145 13]);
        uicontrol('Parent',a.p14UI1,'Style','text','HorizontalAlignment','left','String','Sweep sections','Position',[10 15 145 13]);
        a.Force2DTracingUI1 = uicontrol('Parent',a.p14UI1,'Style','checkbox','Value',a.Force2DTracing1,'Tooltip','','Position',[160 130 50 23],'Callback',@CallbackUI1,'Tag','77');
        a.PhaseVelocityLimit2UI1 = uicontrol('Parent',a.p14UI1,'Style','edit','String',a.PhaseVelocityLimit21/1e3,'Tooltip','Enter the phase velocity (Y-axis in the dispersion diagram) up to which the dispersion curves shall be traced in case there is attenuation or 2-D tracing is forced.','Position',[160 100 50 23],'Callback',@CallbackUI1,'Tag','79');
        a.AttenuationLimitUI1 = uicontrol('Parent',a.p14UI1,'Style','edit','String',a.AttenuationLimit1,'Tooltip','Define the maximum attenuation covered by the initial sweeps for finding modal solutions as the ratio of the imaginary part of the wavenumber to the real part of the wavenumber. Hence, increasing this ratio increases the attenuation search range. However, note that once DC has found a solution and starts tracing a dispersion curve, DC tries to continue tracing until either ''Frequency limit'' or ''Phase velocity limit 2'' are reached. DC does not stop tracing when a dispersion curve exceeds ''Attenuation limit''.','Position',[160 70 50 23],'Callback',@CallbackUI1,'Tag','93');
        a.SweepsUI1 = uicontrol('Parent',a.p14UI1,'Style','popupmenu','String',{'1','2','4','8'},'Tooltip','Define the number of extra initial sweeps for finding modal solutions (DC first performs sweeps at the cut-off frequencies found in the ''Search higher order modes'' panel). Dispersion curve tracing starts from the solutions found by those sweeps. Performing more extra sweeps increases the likelihood to find all modes. By default, DC performs two extra sweeps close to zero frequency, one sweep close to the ''Frequency limit'', and two sweeps evenly distributed along the frequency range. The number of those evenly distributed sweeps is two multiplied by ''Sweeps''. Hence, beside the three fixed sweeps, you can distribute a maximum number of 2*8 = 16 sweeps.','Position',[150 40 60 23],'Callback',@CallbackUI1,'Tag','94');
        a.SweepsSectionsUI1 = uicontrol('Parent',a.p14UI1,'Style','popupmenu','String',{'256','512','1024','2048'},'Tooltip','Define the number of samples in the real and imaginary wavenumber parts at which the dispersion equation is evaluated. For example, setting ''Sweep sections'' to 256 means that the dispersion equation is calculated for a quadratic grid of 256 x 256 wavenumbers (ignoring those where k_imag/k_real is larger than ''Attenuation limit'') on both the positive and negative attenuation sides, resulting in a total grid of 256 x 512. Increasing ''Sweep sections'' improves the chance to find weakly pronounced minima (modal solutions) in the dispersion equation amplitude.','Position',[150 10 60 23],'Callback',@CallbackUI1,'Tag','95');
        if  a.Viscoelastic1
            a.Force2DTracingUI1.Enable = 'off';
            a.PhaseVelocityLimit2UI1.Enable = 'on';
            a.AttenuationLimitUI1.Enable = 'on';
            a.SweepsUI1.Enable = 'on';
            a.SweepsSectionsUI1.Enable = 'on';
        else
            a.Force2DTracingUI1.Enable = 'on';
            a.PhaseVelocityLimit2UI1.Enable = 'off';
            a.AttenuationLimitUI1.Enable = 'off';
            a.SweepsUI1.Enable = 'off';
            a.SweepsSectionsUI1.Enable = 'off';
        end

        a.p3UI1 = uipanel('Parent',Tab1,'Title','Mode selection','Units','pixels','Position',[10 10 225 175],'FontSize',10);
        uicontrol('Parent',a.p3UI1,'Style','text','HorizontalAlignment','left','String','Higher order modes','Position',[10 135 97 13]);
        a.Symmetric_Torsional_ModesTextUI1 = uicontrol('Parent',a.p3UI1,'Style','text','HorizontalAlignment','left','String','Symmetric modes','Position',[10 105 87 13]);
        a.Antisymmetric_Longitudinal_ModesTextUI1 = uicontrol('Parent',a.p3UI1,'Style','text','HorizontalAlignment','left','String','Antisymmetric modes','Position',[10 75 105 13]);
        a.Lamb_Flexural_ModesTextUI1 = uicontrol('Parent',a.p3UI1,'Style','text','HorizontalAlignment','left','String','Lamb modes','Position',[10 45 75 13]);
        a.ShearHorizontalModes_FlexuralModeOrdersTextUI1 = uicontrol('Parent',a.p3UI1,'Style','text','HorizontalAlignment','left','String','Shear horizontal modes','Position',[10 15 117 13]);
        a.HigherOrderModesUI1 = uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.HigherOrderModes1,'Tooltip','Check this in order to calculate the higher order modes in addition to the fundamental ones.','Position',[160 130 50 23],'Callback',@CallbackUI1,'Tag','8');
        a.Symmetric_Torsional_ModesUI1 = uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.SymmetricModes1,'Tooltip','Check this in order to calculate the symmetric modes. These modes have a symmetric displacement pattern with respect to the middle plane of the plate.','Position',[160 100 50 23],'Callback',@CallbackUI1,'Tag','9');
        a.Antisymmetric_Longitudinal_ModesUI1 = uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.AntisymmetricModes1,'Tooltip','Check this in order to calculate the antisymmetric modes. These modes have an antisymmetric displacement pattern with respect to the middle plane of the plate.','Position',[160 70 50 23],'Callback',@CallbackUI1,'Tag','10');
        a.Lamb_Flexural_ModesUI1 = uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.LambModes1,'Tooltip','Check this in order to calculate the Lamb wave modes. These modes show displacement only in the sagittal plane spanned by the propagation direction x1 and by the out-of-plane direction x3. These kind of waves are termed ''pure'' Lamb waves. Lamb waves are indicated by solid lines in the dispersion diagram.','Position',[160 40 50 23],'Callback',@CallbackUI1,'Tag','11');
        a.ShearHorizontalModes_FlexuralModeOrdersUI1 = uicontrol('Parent',a.p3UI1,'Style','checkbox','Value',a.ShearHorizontalModes1,'Tooltip','Check this in order to calculate the shear horizontal modes. These modes show displacement only perpendicular (x2) to the propagation direction x1, and are therefore termed ''pure'' modes. Shear horizontal waves are indicated by dashed lines in the dispersion diagram.','Position',[160 10 50 23],'Callback',@CallbackUI1,'Tag','12');
        if  a.Viscoelastic1
            a.HigherOrderModesUI1.Enable = 'off';
        else
            a.HigherOrderModesUI1.Enable = 'on';
        end

        %------------------------------------------------------------------
        a.p4UI1 = uipanel('Parent',Tab1,'Title','Search higher order modes','Units','pixels','Position',[245 690 223 65],'FontSize',10);
        uicontrol('Parent',a.p4UI1,'Style','text','HorizontalAlignment','left','String','Step (kHz)','Position',[10 20 53 13]);
        a.StepUI1 = uicontrol('Parent',a.p4UI1,'Style','edit','String',a.Step1,'Tooltip','Enter the step size for the frequency sweep. In general, a finer step increases the chance to find mode cut-off frequencies, although exceptions from that rule might occur occasionally. It is essential that all higher order modes are detected. Therefore, try smaller and larger step sizes if you are not sure that you have found all modes.','Position',[75 15 50 23],'Callback',@CallbackUI1,'Tag','13');
        uicontrol('Parent',a.p4UI1,'Style','pushbutton','String','Search','Tooltip','Before you can trace the higher order modes, their cut-off frequencies at the phase velocity limit must be detected. This is done automatically upon pressing ''calculate'' if it has not already been done manually. Sometimes the automatic search does not find all modes, then change ''Step'' and press ''Search'' to find the missing modes.','Position',[140 10 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','14');

        a.OutputWindow1aUI1 = uicontrol('Parent',Tab1,'Style','text','String','','Position',[245 50 105 630],'BackgroundColor','white');
        a.OutputWindow1bUI1 = uicontrol('Parent',Tab1,'Style','text','String','','Position',[362 50 105 630],'BackgroundColor','white');
        a.OutputWindow2aUI1 = uicontrol('Parent',Tab1,'Style','text','String','','Position',[245 10 105 28],'BackgroundColor','white');
        a.OutputWindow2bUI1 = uicontrol('Parent',Tab1,'Style','text','String','','Position',[362 10 105 28],'BackgroundColor','white');

        %------------------------------------------------------------------
        a.p12UI1 = uipanel('Parent',Tab1,'Title','Trace modes and calculate energy velocity','Units','pixels','Position',[478 635 263 120],'FontSize',10);
        a.SamplesX3ceTextUI1 = uicontrol('Parent',a.p12UI1,'Style','text','HorizontalAlignment','left','String','Samples x3','Position',[30 40 58 13]);
        a.SamplesX3ceUI1 = uicontrol('Parent',a.p12UI1,'Style','edit','String',a.SamplesX3ce1,'Tooltip','Enter the number of through-thickness samples for the energy velocity calculation. The more samples the more accurate the result. The accuracy decreases with an increasing frequency-thickness product because the mode shapes become more complex. Therefore, choose a sufficiently high number of samples.','Position',[34 15 50 23],'Callback',@CallbackUI1,'Tag','91');        
        a.TraceModesUI1 = uicontrol('Parent',a.p12UI1,'Style','togglebutton','String','Trace modes','Tooltip','Start tracing the dispersion curves.','Position',[117 55 130 40],'FontSize',10,'Callback',@CallbackUI1,'BusyAction','cancel','Tag','15');
        a.CalculateCeUI1 = uicontrol('Parent',a.p12UI1,'Style','togglebutton','String','Calculate ce','Tooltip','Start calculating the energy velocity.','Position',[117 10 130 40],'FontSize',10,'Callback',@CallbackUI1,'BusyAction','cancel','Tag','16');

        a.p5UI1 = uipanel('Parent',Tab1,'Title','Dispersion diagrams','Units','pixels','Position',[478 405 263 225],'FontSize',10);
        uicontrol('Parent',a.p5UI1,'Style','text','HorizontalAlignment','left','String','Quantity','Position',[10 185 42 13]);
        a.Option1TextUI1 = uicontrol('Parent',a.p5UI1,'Style','text','HorizontalAlignment','left','String','Bulk velocities','Position',[10 155 71 13]);
        a.XAxisModeTextUI1 = uicontrol('Parent',a.p5UI1,'Style','text','HorizontalAlignment','left','String','X-axis mode','Position',[10 125 63 13]);
        a.XAxisTextUI1 = uicontrol('Parent',a.p5UI1,'Style','text','HorizontalAlignment','left','String','X-axis (kHz)','Position',[10 95 86 13]);
        a.YAxisTextUI1 = uicontrol('Parent',a.p5UI1,'Style','text','HorizontalAlignment','left','String','Y-axis (m/ms)','Position',[10 65 101 13]);
        a.Quantity1UI1 = uicontrol('Parent',a.p5UI1,'Style','popupmenu','String',{'Phase velocity (m/ms)','Energy velocity (m/ms)',['Propagation time (',char(181),'s)'],['Coincidence angle (',char(176),')'],'Wavelength (mm)','Wavenumber (rad/mm)','Attenuation (Np/m)'},'Tooltip','Select which quantity to plot in the dispersion diagram.','Position',[117 180 130 23],'Callback',@CallbackUI1,'Tag','17');
        a.Option1UI1 = uicontrol('Parent',a.p5UI1,'Style','checkbox','Value',a.BulkVelocities1,'Tooltip','Check this to show the bulk wave velocities.','Position',[117 150 130 23],'Callback',@CallbackUI1,'Tag','18');
        a.XAxisModeUI1 = uicontrol('Parent',a.p5UI1,'Style','popupmenu','String',{'Frequency (kHz)','Frequency (MHz)',['f',char(8901),'d (MHz',char(8901),'mm)']},'Tooltip','Select the frequency''s dimension on the X-axis.','Position',[117 120 130 23],'Callback',@CallbackUI1,'Tag','19');
        a.XAxisUI1 = uicontrol('Parent',a.p5UI1,'Style','edit','String',['[0 ',num2str(a.FrequencyLimit1),']'],'Tooltip','Enter which frequency range shall be plotted.','Position',[117 90 75 23],'Callback',@CallbackUI1,'Tag','20');
        a.YAxisUI1 = uicontrol('Parent',a.p5UI1,'Style','edit','String',['[0 ',num2str(a.PhaseVelocityLimit1/1e3),']'],'Tooltip','Enter which phase velocity range shall be plotted.','Position',[117 60 75 23],'Callback',@CallbackUI1,'Tag','21');
        a.Plot1UI1 = uicontrol('Parent',a.p5UI1,'Style','pushbutton','String','Plot','Tooltip','Plot the dispersion diagram. If you have checked ''Export plots'' in the export settings, the plot will be exported automatically.','Position',[117 15 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','22');

        a.p6UI1 = uipanel('Parent',Tab1,'Title','Through-thickness profiles','Units','pixels','Position',[751 405 224 290],'FontSize',10);
        uicontrol('Parent',a.p6UI1,'Style','text','HorizontalAlignment','left','String','Quantity','Position',[10 245 42 13]);
        uicontrol('Parent',a.p6UI1,'Style','text','HorizontalAlignment','left','String','Mode','Position',[10 215 28 13]);
        uicontrol('Parent',a.p6UI1,'Style','text','HorizontalAlignment','left','String','Frequency (kHz)','Position',[10 185 83 13]);
        a.SamplesX3TextUI1 = uicontrol('Parent',a.p6UI1,'Style','text','HorizontalAlignment','left','String','Samples x3','Position',[10 155 58 13]);
        uicontrol('Parent',a.p6UI1,'Style','text','HorizontalAlignment','left','String','Half-spaces','Position',[10 125 61 13]);
        uicontrol('Parent',a.p6UI1,'Style','text','HorizontalAlignment','left','String','Phase','Position',[10 95 32 13]);
        a.Quantity2UI1 = uicontrol('Parent',a.p6UI1,'Style','popupmenu','String',{'Displacement','Stress','Strain','Energy density','Power flow density'},'Tooltip','Select which through-thickness quantity to plot.','Position',[119 240 90 23],'Callback',@CallbackUI1,'Tag','23');
        a.Mode1UI1 = uicontrol('Parent',a.p6UI1,'Style','popupmenu','String',{''},'Tooltip','Select the mode you want to analyze.','Position',[119 210 90 23],'Callback',@CallbackUI1,'Tag','24');
        a.Frequency1UI1 = uicontrol('Parent',a.p6UI1,'Style','edit','String',a.FrequencyLimit1,'Tooltip','Enter the frequency at which to analyze the selected mode.','Position',[119 180 50 23],'Callback',@CallbackUI1,'Tag','25');
        a.Samples1UI1 = uicontrol('Parent',a.p6UI1,'Style','edit','String',a.Samples11,'Tooltip','Enter the number of sample points over the plate''s thickness (x3) at which the selected quantities are calculated.','Position',[119 150 50 23],'Callback',@CallbackUI1,'Tag','26');
        a.Halfspaces1UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','Value',a.Halfspaces11,'Tooltip','Check this to show the quantities in the upper and lower fluid.','Position',[89 120 20 23],'Enable','off','Callback',@CallbackUI1,'Tag','81');        
        a.HalfspacesNumber1UI1 = uicontrol('Parent',a.p6UI1,'Style','edit','String',a.HalfspacesNumber11,'Tooltip','Set the height of the half-spaces in plate thicknesses.','Position',[119 120 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','80');
        a.PhaseUI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','Value',a.Phase1,'Tooltip','Check this to plot also the phase of the field components.','Position',[119 90 20 23],'Callback',@CallbackUI1,'Tag','84');
        a.Plot2UI1 = uicontrol('Parent',a.p6UI1,'Style','pushbutton','String','Plot','Tooltip','Plot the profile. If you have checked ''Export plots'' in the export settings, the plot will be exported automatically.','Position',[119 50 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','27');
        a.Plot2UIa1 = uicontrol('Parent',a.p6UI1,'Style','pushbutton','String','Live-plot','Tooltip','The profiles are updating as you are sliding along the dispersion curves.','Position',[119 15 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','92');
        a.x11UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','11','Value',a.Plot1(1),'Tooltip','','Position',[10 60 35 23],'Enable','off','Callback',@CallbackUI1,'Tag','28');
        a.x22UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','22','Value',a.Plot1(2),'Tooltip','','Position',[45 40 35 23],'Enable','off','Callback',@CallbackUI1,'Tag','29');
        a.x33UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','3','Value',a.Plot1(3),'Tooltip','Displacement u3','Position',[80 20 40 23],'Callback',@CallbackUI1,'Tag','30');
        a.x13UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','1','Value',a.Plot1(4),'Tooltip','Displacement u1','Position',[80 60 40 23],'Callback',@CallbackUI1,'Tag','31');
        a.x23UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','2','Value',a.Plot1(5),'Tooltip','Displacement u2','Position',[80 40 40 23],'Callback',@CallbackUI1,'Tag','32');
        a.x12UI1 = uicontrol('Parent',a.p6UI1,'Style','checkbox','String','12','Value',a.Plot1(6),'Tooltip','','Position',[45 60 35 23],'Enable','off','Callback',@CallbackUI1,'Tag','33');
        
        a.p7UI1 = uipanel('Parent',Tab1,'Title','Mode shape','Units','pixels','Position',[478 195 497 205],'FontSize',10);
        uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Mode, plane','Position',[10 165 60 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Frequency (kHz)','Position',[10 135 83 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Wavelengths','Position',[10 105 65 13]);
        a.SamplesX1X3TextUI1 = uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Samples x1,x3','Position',[10 75 73 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Half-spaces','Position',[10 45 61 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Grid line every','Position',[10 15 72 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Gain','Position',[185 15 24 13]);        
        uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Animate','Position',[282 45 41 13]);
        uicontrol('Parent',a.p7UI1,'Style','text','HorizontalAlignment','left','String','Undistorted','Position',[282 15 57 13]);
        a.Mode2UI1 = uicontrol('Parent',a.p7UI1,'Style','popupmenu','String',{''},'Tooltip','Select the mode you want to analyze.','Position',[117 160 90 23],'Callback',@CallbackUI1,'Tag','34');
        a.PlaneUI1 = uicontrol('Parent',a.p7UI1,'Style','popupmenu','String',{'r-z',['r-',char(952)]},'Tooltip','Select in which plane you want to plot.','Position',[215 160 48 23],'Enable','off','Callback',@CallbackUI1,'Tag','43');
        a.Frequency2UI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.FrequencyLimit1,'Tooltip','Enter the frequency at which to analyze the selected mode.','Position',[117 130 50 23],'Callback',@CallbackUI1,'Tag','35');
        a.LengthUI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.Length1,'Tooltip','Specify how many wavelengths you want to display.','Position',[117 100 50 23],'Callback',@CallbackUI1,'Tag','36');
        a.Samples2UI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.Samples21,'Tooltip','Enter the number of sample points along the propagation direction (x1) at which the displacement is calculated.','Position',[117 70 50 23],'Callback',@CallbackUI1,'Tag','37');
        a.Samples3UI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.Samples31,'Tooltip','Enter the number of sample points over the plate''s thickness (x3) at which the displacement is calculated.','Position',[177 70 50 23],'Callback',@CallbackUI1,'Tag','38');
        a.Halfspaces2UI1 = uicontrol('Parent',a.p7UI1,'Style','checkbox','Value',a.Halfspaces21,'Tooltip','Check this to show the quantities in the upper and lower fluid.','Position',[87 40 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','83');
        a.HalfspacesNumber2UI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.HalfspacesNumber21,'Tooltip','Set the height of the half-spaces in plate thicknesses.','Position',[117 40 50 23],'Enable','off','Callback',@CallbackUI1,'Tag','82');        
        a.GridLineUI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.GridLine1,'Tooltip','Draw a grid line at every ith sample point.','Position',[117 10 50 23],'Callback',@CallbackUI1,'Tag','40');
        a.GainUI1 = uicontrol('Parent',a.p7UI1,'Style','edit','String',a.Gain1,'Tooltip','Magnify the displacement.','Position',[220 10 50 23],'Callback',@CallbackUI1,'Tag','39');
        a.AnimateUI1 = uicontrol('Parent',a.p7UI1,'Style','checkbox','Value',a.Animate1,'Tooltip','Check this in order to show the animated mode shape upon pressing the plot button below.','Position',[352 40 70 23],'Callback',@CallbackUI1,'Tag','47');
        a.UndistortedUI1 = uicontrol('Parent',a.p7UI1,'Style','checkbox','Value',a.Undistorted1,'Tooltip','Check this to draw the undistorted grid.','Position',[352 10 50 23],'Callback',@CallbackUI1,'Tag','41');
        a.Plot3UI1 = uicontrol('Parent',a.p7UI1,'Style','pushbutton','String','Plot','Tooltip','Plot the mode shape. If you have checked ''Animate'', the mode shape will be animated. If you have checked ''Export plots'' in the export settings, the plot/movie will be exported automatically.','Position',[392 15 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','42');

        a.p8UI1 = uipanel('Parent',Tab1,'Title','Animation settings','Units','pixels','Position',[751 265 212 115],'FontSize',9);
        uicontrol('Parent',a.p8UI1,'Style','text','HorizontalAlignment','left','String','Cycle duration (s)','Position',[10 75 88 13]);
        uicontrol('Parent',a.p8UI1,'Style','text','HorizontalAlignment','left','String','Frame rate (Hz)','Position',[10 45 78 13]);
        uicontrol('Parent',a.p8UI1,'Style','text','HorizontalAlignment','left','String','Movie quality (0-100)','Position',[10 15 103 13]);
        a.CycleDurationUI1 = uicontrol('Parent',a.p8UI1,'Style','edit','String',a.CycleDuration1,'Tooltip','Enter how long a cycle shall take.','Position',[119 70 50 23],'Callback',@CallbackUI1,'Tag','44');
        a.FrameRateUI1 = uicontrol('Parent',a.p8UI1,'Style','edit','String',a.FrameRate1,'Tooltip','Define the frame rate of the movie.','Position',[119 40 50 23],'Callback',@CallbackUI1,'Tag','45');
        a.MovieQualityUI1 = uicontrol('Parent',a.p8UI1,'Style','edit','String',a.MovieQuality1,'Tooltip','Define the quality of the exported movie.','Position',[119 10 50 23],'Callback',@CallbackUI1,'Tag','46');
        
        a.p9UI1 = uipanel('Parent',Tab1,'Title','Export settings','Units','pixels','Position',[478 10 497 180],'FontSize',10);
        uicontrol('Parent',a.p9UI1,'Style','text','HorizontalAlignment','left','String','Export plots','Position',[10 140 59 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','HorizontalAlignment','left','String','PDF','Position',[124 140 26 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','HorizontalAlignment','left','String','PNG','Position',[124 120 28 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','HorizontalAlignment','left','String','PNG resolution (dpi)','Position',[10 93 98 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','HorizontalAlignment','left','String','Dispersion curves','Position',[240 140 90 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','HorizontalAlignment','left','String','Through-thickness','Position',[240 120 92 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','HorizontalAlignment','left','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p9UI1,'Style','text','HorizontalAlignment','left','String','Directory','Position',[10 20 46 13]);
        a.ExportPlotsUI1 = uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.ExportPlots1,'Tooltip','Check this in order to export plots/movies automatically upon pressing the respective plot button. You can also export plots manually by using the ''File'' menu inside the plot figure.','Position',[80 135 20 23],'Callback',@CallbackUI1,'Tag','48');
        a.PDFUI1 = uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.PDF1,'Tooltip','Check this to export a plot as pdf.','Position',[160 135 20 23],'Callback',@CallbackUI1,'Tag','49');
        a.PNGUI1 = uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.PNG1,'Tooltip','Check this to export a plot as png.','Position',[160 115 20 23],'Callback',@CallbackUI1,'Tag','50');
        a.PNGresolutionUI1 = uicontrol('Parent',a.p9UI1,'Style','edit','String',a.PNGresolution1,'Tooltip','Enter the resolution of the png image.','Position',[125 88 50 23],'Callback',@CallbackUI1,'Tag','51');
        a.XAxisMode2UI1 = uicontrol('Parent',a.p9UI1,'Style','popupmenu','String',{'Frequency (kHz)','Frequency (MHz)',['f',char(8901),'d (MHz',char(8901),'mm)']},'Tooltip','Select the frequency''s dimension to be exported.','Position',[370 130 110 23],'Callback',@CallbackUI1,'Tag','52');
        a.ArrangeUI1 = uicontrol('Parent',a.p9UI1,'Style','popupmenu','String',{'Horizontal arrangement','Vertical arrangement'},'Tooltip','Choose how to arrange the dispersion curve data in the Excel/txt.','Position',[370 95 110 23],'Callback',@CallbackUI1,'Tag','53');        
        a.DispersionCurvesUI1 = uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.DispersionCurves1,'Tooltip','Check this to export the dispersion curves.','Position',[340 135 20 23],'Callback',@CallbackUI1,'Tag','54');
        a.ThroughThicknessUI1 = uicontrol('Parent',a.p9UI1,'Style','checkbox','Value',a.ThroughThickness1,'Tooltip','Check this to export the through-thickness profiles selected above.','Position',[340 115 20 23],'Callback',@CallbackUI1,'Tag','55');
        uicontrol('Parent',a.p9UI1,'Style','pushbutton','String','*.mat','Tooltip','Export the data as Matlab''s mat-file.','Position',[220 80 40 33],'Callback',@CallbackUI1,'Tag','78');
        uicontrol('Parent',a.p9UI1,'Style','pushbutton','String','*.xlsx','Tooltip','Export the data as Excel sheet.','Position',[270 80 40 33],'Callback',@CallbackUI1,'Tag','56');
        uicontrol('Parent',a.p9UI1,'Style','pushbutton','String','*.txt','Tooltip','Export the data as txt-file.','Position',[320 80 40 33],'Callback',@CallbackUI1,'Tag','57');
        a.FileNameUI1 = uicontrol('Parent',a.p9UI1,'Style','edit','String',a.FileName1,'Tooltip','Specify the name of the plots/movies to be exported. The dispersion curve raw data are named automatically.','Position',[70 45 412 23],'Callback',@CallbackUI1,'Tag','58');
        a.DirectoryUI1 = uicontrol('Parent',a.p9UI1,'Style','edit','String',a.Directory,'Tooltip','Specify the directory to which plots, movies, and the dispersion curve raw data shall be exported.','Position',[70 15 412 23],'Callback',@CallbackUI1,'Tag','59');

        %------------------------------------------------------------------
        a.p10UI1 = uipanel('Parent',Tab1,'Title','Plot layout settings','Units','pixels','Position',[985 275 200 420],'FontSize',10);
        uicontrol('Parent',a.p10UI1,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 380 21 13]);
        uicontrol('Parent',a.p10UI1,'Style','text','HorizontalAlignment','left','String','Legend location','Position',[10 350 78 13]);
        uicontrol('Parent',a.p10UI1,'Style','text','HorizontalAlignment','left','String','Box line width','Position',[10 320 70 13]);
        uicontrol('Parent',a.p10UI1,'Style','text','HorizontalAlignment','left','String','Curve line width','Position',[10 290 80 13]);
        a.TitleUI1 = uicontrol('Parent',a.p10UI1,'Style','checkbox','Value',a.Title1,'Tooltip','Check this in order to show the plot title.','Position',[105 375 20 23],'Callback',@CallbackUI1,'Tag','60');
        a.LegendLocationUI1 = uicontrol('Parent',a.p10UI1,'Style','popupmenu','String',{'outside','inside'},'Tooltip','Determine the legend location in the through-thickness plots.','Position',[105 345 80 23],'Callback',@CallbackUI1,'Tag','62');
        a.BoxLineWidthUI1 = uicontrol('Parent',a.p10UI1,'Style','edit','String',a.BoxLineWidth1,'Tooltip','Enter the box line width.','Position',[105 315 50 23],'Callback',@CallbackUI1,'Tag','63');
        a.LineWidthUI1 = uicontrol('Parent',a.p10UI1,'Style','edit','String',a.LineWidth1,'Tooltip','Enter the curve line width.','Position',[105 285 50 23],'Callback',@CallbackUI1,'Tag','64');

        a.p11UI1 = uipanel('Parent',Tab1,'Title','Dispersion curve colors [R G B]','Units','pixels','Position',[995 435 180 115],'FontSize',9); 
        a.SColorTextUI1 = uicontrol('Parent',a.p11UI1,'Style','text','HorizontalAlignment','left','String','S','Position',[10 75 17 13]);
        a.AColorTextUI1 = uicontrol('Parent',a.p11UI1,'Style','text','HorizontalAlignment','left','String','A','Position',[10 45 8 13]);
        a.BColorTextUI1 = uicontrol('Parent',a.p11UI1,'Style','text','HorizontalAlignment','left','String','B','Position',[10 15 7 13]);
        a.SColorUI1 = uicontrol('Parent',a.p11UI1,'Style','edit','String','[1 0 0]','Tooltip','Specify the color of symmetric modes. You can compose any color by entering the corresponding RGB values. Set numbers from 0 to 1.','Position',[95 70 50 23],'Callback',@CallbackUI1,'Tag','66');
        a.AColorUI1 = uicontrol('Parent',a.p11UI1,'Style','edit','String','[0 0 1]','Tooltip','Specify the color of antisymmetric modes.','Position',[95 40 50 23],'Callback',@CallbackUI1,'Tag','67');
        a.BColorUI1 = uicontrol('Parent',a.p11UI1,'Style','edit','String','[.5 0 1]','Tooltip','Specify the color of nonsymmetric modes.','Position',[95 10 50 23],'Callback',@CallbackUI1,'Tag','69');

        a.p13UI1 = uipanel('Parent',Tab1,'Title','Font size','Units','pixels','Position',[995 285 180 145],'FontSize',9);
        uicontrol('Parent',a.p13UI1,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 105 21 13]);
        uicontrol('Parent',a.p13UI1,'Style','text','HorizontalAlignment','left','String','Axes labels','Position',[10 75 59 13]);
        uicontrol('Parent',a.p13UI1,'Style','text','HorizontalAlignment','left','String','Axes ticks','Position',[10 45 53 13]);
        uicontrol('Parent',a.p13UI1,'Style','text','HorizontalAlignment','left','String','Legend','Position',[10 15 38 13]);
        a.TitleFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.TitleFontSize1,'Tooltip','Set the title font size.','Position',[95 100 50 23],'Callback',@CallbackUI1,'Tag','71');
        a.AxesLabelFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.AxesLabelFontSize1,'Tooltip','Set the axes label font size.','Position',[95 70 50 23],'Callback',@CallbackUI1,'Tag','72');
        a.AxesTickFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.AxesTickFontSize1,'Tooltip','Set the axes tick label font size.','Position',[95 40 50 23],'Callback',@CallbackUI1,'Tag','73');
        a.LegendFontSizeUI1 = uicontrol('Parent',a.p13UI1,'Style','edit','String',a.LegendFontSize1,'Tooltip','Set the legend font size.','Position',[95 10 50 23],'Callback',@CallbackUI1,'Tag','75');

        a.c3UI1 = uicontrol('Parent',Tab1,'Style','pushbutton','String','Default','Tooltip','Reset the plot layout settings to the default.','Position',[1052 230 65 33],'FontSize',10,'Callback',@CallbackUI1,'Tag','76');    

        uicontrol('Parent',Tab1,'Style','text','HorizontalAlignment','left','String','CPU cores','Position',[1035 715 53 13]);
        a.c4UI1 = uicontrol('Parent',Tab1,'Style','text','String','0','Position',[1095 715 30 13],'BackgroundColor','white');
        a.a1UI1 = axes('Parent',Tab1,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading20.png')
    end
    if  1 % Tab2_anisotropic
        a.p1UI2 = uipanel('Parent',Tab2,'Title','Specimen','Units','pixels','Position',[10 580 225 175],'FontSize',10);
        uicontrol('Parent',a.p1UI2,'Style','pushbutton','String','Edit','Tooltip','Press to define your specimen.','Position',[105 110 105 40],'FontSize',10,'Callback',@SpecimenSettingsUI2_Callback);
        uicontrol('Parent',a.p1UI2,'Style','text','HorizontalAlignment','left','String','Geometry','Position',[10 130 60 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','HorizontalAlignment','left','String','Fluids','Position',[10 90 30 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','HorizontalAlignment','left','String','Material','Position',[10 70 39 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','HorizontalAlignment','left','String','Layup','Position',[10 50 32 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','HorizontalAlignment','left','String','Effective','Position',[10 30 45 13]);
        uicontrol('Parent',a.p1UI2,'Style','text','HorizontalAlignment','left','String','Layers , d (mm)','Position',[10 10 78 13]);
        a.GeometryDisplayUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String',a.Geometry2,'Tooltip','Geometry.','Position',[10 110 80 13],'BackgroundColor','white');
        a.Fluid1DisplayUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','vacuum','Tooltip','Fluid above the layup.','Position',[58 90 70 13],'BackgroundColor','white');
        a.Fluid2DisplayUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','vacuum','Tooltip','Fluid below the layup.','Position',[140 90 70 13],'BackgroundColor','white');
        a.MaterialNameUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String',a.Material2{1}.Name,'Position',[58 70 152 13],'BackgroundColor','white');
        a.LayupUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','[0]','Position',[58 50 152 13],'BackgroundColor','white');
        a.EffectiveLayupUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','[0]','Tooltip','= Layup - Propagation angle','Position',[58 30 152 13],'BackgroundColor','white');
        a.LayerCountUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String','1','Tooltip','Number of layers in the layup.','Position',[100 10 50 13],'BackgroundColor','white');
        a.ThicknessCountUI2 = uicontrol('Parent',a.p1UI2,'Style','text','String',a.LayupThickness,'Tooltip','Overall thickness of the layup.','Position',[160 10 50 13],'BackgroundColor','white');

        a.p2UI2 = uipanel('Parent',Tab2,'Title','Computational settings','Units','pixels','Position',[10 370 225 205],'FontSize',10);
        uicontrol('Parent',Tab2,'Style','togglebutton','String','Fix','Tooltip','Fix the computational settings such that they are not adjusted automatically.','Position',[212 555 23 23],'Callback',@CallbackUI2,'Tag','63');
        uicontrol('Parent',a.p2UI2,'Style','text','HorizontalAlignment','left','String',['Propagation angle (',char(176),')'],'Position',[10 165 103 13]);
        uicontrol('Parent',a.p2UI2,'Style','text','HorizontalAlignment','left','String','Frequency limit (kHz)','Position',[10 135 103 13]);
        uicontrol('Parent',a.p2UI2,'Style','text','HorizontalAlignment','left','String','Frequency step (kHz)','Position',[10 105 107 13]);
        uicontrol('Parent',a.p2UI2,'Style','text','HorizontalAlignment','left','String','Phase velocity limit (m/ms)','Position',[10 75 128 13]);
        uicontrol('Parent',a.p2UI2,'Style','text','HorizontalAlignment','left','String','Phase velocity accuracy (m/s)','Position',[10 45 149 13]);
        uicontrol('Parent',a.p2UI2,'Style','text','HorizontalAlignment','left','String','Matrix methods','Position',[10 15 75 13]);
        a.PropagationAngleUI2 = uicontrol('Parent',a.p2UI2,'Style','edit','String',a.PropagationAngle,'Tooltip','Enter the angle with respect to the fiber orientations along which the guided waves shall propagate.','Position',[160 160 50 23],'Callback',@CallbackUI2,'Tag','1');
        a.FrequencyLimitUI2 = uicontrol('Parent',a.p2UI2,'Style','edit','String',a.FrequencyLimit2,'Tooltip','Enter the frequency (X-axis in the dispersion diagram) up to which the dispersion curves shall be traced.','Position',[160 130 50 23],'Callback',@CallbackUI2,'Tag','3');
        a.FrequencyResolutionUI2 = uicontrol('Parent',a.p2UI2,'Style','edit','String',a.FrequencyResolution2,'Tooltip','Enter the frequency step.','Position',[160 100 50 23],'Callback',@CallbackUI2,'Tag','4');        
        a.PhaseVelocityLimitUI2 = uicontrol('Parent',a.p2UI2,'Style','edit','String',a.PhaseVelocityLimit2/1e3,'Tooltip','Enter the phase velocity (Y-axis in the dispersion diagram) up to which the dispersion curves shall be traced.','Position',[160 70 50 23],'Callback',@CallbackUI2,'Tag','2');
        a.AccuracyUI2 = uicontrol('Parent',a.p2UI2,'Style','edit','String',a.Accuracy2,'Tooltip','Bisections are performed until the modal solutions have reached this accuracy.','Position',[160 40 50 23],'Callback',@CallbackUI2,'Tag','82');
        a.MatrixMethodsUI2 = uicontrol('Parent',a.p2UI2,'Style','popupmenu','String',{'SMM+TMM','SMM'},'Tooltip','Choose if the transfer matrix method (TMM) shall be used instead of the stiffness matrix method (SMM) whenever there are stable conditions for TMM. Stable conditions for TMM are above a diagonal line starting in the origin of the frequency- phase velocity dispersion diagram. SMM is used below the diagonal. DC determines the slope of the diagonal before starting dispersion curve tracing. TMM is faster and more sensitive to quasi-shear horizontal modes than SMM. Another advantage is that its dispersion equation amplitude does not have maxima, unlike the case with SMM. This makes finding the minima (modal solutions) easier. However, sometimes SMM alone brings better results. Try both options in difficult cases.','Position',[120 10 90 23],'Callback',@CallbackUI2,'Tag','84');

        a.p17UI2 = uipanel('Parent',Tab2,'Title','2-D tracing settings','Units','pixels','Position',[10 190 225 175],'FontSize',10);
        uicontrol('Parent',a.p17UI2,'Style','text','HorizontalAlignment','left','String','Force 2-D tracing','Position',[10 135 145 13]);
        uicontrol('Parent',a.p17UI2,'Style','text','HorizontalAlignment','left','String','Phase velocity limit 2 (m/ms)','Position',[10 105 145 13]);
        uicontrol('Parent',a.p17UI2,'Style','text','HorizontalAlignment','left','String','Attenuation limit (k_im/k_re)','Position',[10 75 145 13]);
        uicontrol('Parent',a.p17UI2,'Style','text','HorizontalAlignment','left','String','Sweeps','Position',[10 45 145 13]);
        uicontrol('Parent',a.p17UI2,'Style','text','HorizontalAlignment','left','String','Sweep sections','Position',[10 15 145 13]);
        a.Force2DTracingUI2 = uicontrol('Parent',a.p17UI2,'Style','checkbox','Value',a.Force2DTracing2,'Tooltip','','Position',[160 130 50 23],'Callback',@CallbackUI2,'Tag','76');
        a.PhaseVelocityLimit2UI2 = uicontrol('Parent',a.p17UI2,'Style','edit','String',a.PhaseVelocityLimit22/1e3,'Tooltip','Enter the phase velocity (Y-axis in the dispersion diagram) up to which the dispersion curves shall be traced in case there is attenuation or 2-D tracing is forced.','Position',[160 100 50 23],'Callback',@CallbackUI2,'Tag','79');
        a.AttenuationLimitUI2 = uicontrol('Parent',a.p17UI2,'Style','edit','String',a.AttenuationLimit2,'Tooltip','Define the maximum attenuation covered by the initial sweeps for finding modal solutions as the ratio of the imaginary part of the wavenumber to the real part of the wavenumber. Hence, increasing this ratio increases the attenuation search range. However, note that once DC has found a solution and starts tracing a dispersion curve, DC tries to continue tracing until either ''Frequency limit'' or ''Phase velocity limit 2'' are reached. DC does not stop tracing when a dispersion curve exceeds ''Attenuation limit''.','Position',[160 70 50 23],'Callback',@CallbackUI2,'Tag','86');
        a.SweepsUI2 = uicontrol('Parent',a.p17UI2,'Style','popupmenu','String',{'1','2','4','8'},'Tooltip','Define the number of extra initial sweeps for finding modal solutions (DC first performs sweeps at the cut-off frequencies found in the ''Search higher order modes'' panel). Dispersion curve tracing starts from the solutions found by those sweeps. Performing more extra sweeps increases the likelihood to find all modes. By default, DC performs two extra sweeps close to zero frequency, one sweep close to the ''Frequency limit'', and two sweeps evenly distributed along the frequency range. The number of those evenly distributed sweeps is two multiplied by ''Sweeps''. Hence, beside the three fixed sweeps, you can distribute a maximum number of 2*8 = 16 sweeps.','Position',[150 40 60 23],'Callback',@CallbackUI2,'Tag','87');
        a.SweepsSectionsUI2 = uicontrol('Parent',a.p17UI2,'Style','popupmenu','String',{'256','512','1024','2048'},'Tooltip','Define the number of samples in the real and imaginary wavenumber parts at which the dispersion equation is evaluated. For example, setting ''Sweep sections'' to 256 means that the dispersion equation is calculated for a quadratic grid of 256 x 256 wavenumbers (ignoring those where k_imag/k_real is larger than ''Attenuation limit'') on both the positive and negative attenuation sides, resulting in a total grid of 256 x 512. Increasing ''Sweep sections'' improves the chance to find weakly pronounced minima (modal solutions) in the dispersion equation amplitude.','Position',[150 10 60 23],'Callback',@CallbackUI2,'Tag','88');
        if  a.Viscoelastic2
            a.Force2DTracingUI2.Enable = 'off';
            a.PhaseVelocityLimit2UI2.Enable = 'on';
            a.AttenuationLimitUI2.Enable = 'on';
            a.SweepsUI2.Enable = 'on';
            a.SweepsSectionsUI2.Enable = 'on';
        else
            a.Force2DTracingUI2.Enable = 'on';
            a.PhaseVelocityLimit2UI2.Enable = 'off';
            a.AttenuationLimitUI2.Enable = 'off';
            a.SweepsUI2.Enable = 'off';
            a.SweepsSectionsUI2.Enable = 'off';
        end

        a.p3UI2 = uipanel('Parent',Tab2,'Title','Mode selection','Units','pixels','Position',[10 10 225 175],'FontSize',10);
        uicontrol('Parent',a.p3UI2,'Style','text','HorizontalAlignment','left','String','Higher order modes','Position',[10 135 97 13]);
        a.Symmetric_Torsional_ModesTextUI2 = uicontrol('Parent',a.p3UI2,'Style','text','HorizontalAlignment','left','String','Symmetric modes','Position',[10 105 87 13]);
        a.Antisymmetric_Longitudinal_ModesTextUI2 = uicontrol('Parent',a.p3UI2,'Style','text','HorizontalAlignment','left','String','Antisymmetric modes','Position',[10 75 105 13]);
        a.Lamb_Flexural_ModesTextUI2 = uicontrol('Parent',a.p3UI2,'Style','text','HorizontalAlignment','left','String','Lamb modes','Position',[10 45 75 13]);
        a.ShearHorizontalModes_FlexuralModeOrdersTextUI2 = uicontrol('Parent',a.p3UI2,'Style','text','HorizontalAlignment','left','String','Shear horizontal modes','Position',[10 15 117 13]);
        a.HigherOrderModesUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.HigherOrderModes2,'Tooltip','Check this in order to calculate the higher order modes in addition to the fundamental ones.','Position',[160 130 50 23],'Callback',@CallbackUI2,'Tag','7');
        a.Symmetric_Torsional_ModesUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.SymmetricModes2,'Tooltip','If the layup is symmtric with respect to the middle plane, one can distinguish symmetric and antisymmetric modes. Check this in order to calculate the symmetric modes. These modes have a symmetric displacement pattern with respect to the middle plane of the plate.','Position',[160 100 50 23],'Callback',@CallbackUI2,'Tag','9');
        a.Antisymmetric_Longitudinal_ModesUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.AntisymmetricModes2,'Tooltip','Check this in order to calculate the antisymmetric modes. These modes have an antisymmetric displacement pattern with respect to the middle plane of the plate.','Position',[160 70 50 23],'Callback',@CallbackUI2,'Tag','10');
        a.Lamb_Flexural_ModesUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.LambModes2,'Tooltip',['Check this in order to calculate the pure Lamb wave modes. This option is available only in the decoupled case. Lamb and shear horizontal waves are decoupled if the specimen consists solely of 0',char(176),' and 90',char(176),' layers, and if the wave propagation is along either direction. In the decoupled case, Lamb waves are termed ''pure'' modes. These Lamb waves show displacement only in the sagittal plane spanned by the propagation direction x1 and by the out-of-plane direction x3. In the coupled case, the modes are still called Lamb waves, but they have a third displacement component perpendicular to the propagation direction (shear horizontal, x2), and no distinction between Lamb and shear horizontal modes is possible. Lamb waves are indicated by solid lines in the dispersion diagram.'],'Position',[160 40 50 23],'Callback',@CallbackUI2,'Tag','13');
        a.ShearHorizontalModes_FlexuralModeOrdersUI2 = uicontrol('Parent',a.p3UI2,'Style','checkbox','Value',a.ShearHorizontalModes2,'Tooltip','The so-called ''pure'' shear horizontal modes can only exist in the decoupled case. Otherwise, they cannot be distinguished from Lamb waves. The pure shear horizontal modes, however, show displacement only perpendicular (x2) to the propagation direction. Shear horizontal waves are indicated by dashed lines in the dispersion diagram.','Position',[160 10 50 23],'Callback',@CallbackUI2,'Tag','14');
        if  a.Viscoelastic2
            a.HigherOrderModesUI2.Enable = 'off';
        else
            a.HigherOrderModesUI2.Enable = 'on';
        end

        %------------------------------------------------------------------
        a.p7UI2 = uipanel('Parent',Tab2,'Title','Search higher order modes','Units','pixels','Position',[245 690 223 65],'FontSize',10);
        uicontrol('Parent',a.p7UI2,'Style','text','HorizontalAlignment','left','String','Step (kHz)','Position',[10 20 53 13]);
        a.StepUI2 = uicontrol('Parent',a.p7UI2,'Style','edit','String',a.Step2,'Tooltip','Enter the step size for the frequency sweep. In general, a finer step increases the chance to find mode cut-off frequencies, although exceptions from that rule might occur occasionally. It is essential that all higher order modes are detected. Therefore, try smaller and larger step sizes if you are not sure that you have found all modes.','Position',[75 15 50 23],'Callback',@CallbackUI2,'Tag','15');
        uicontrol('Parent',a.p7UI2,'Style','pushbutton','String','Search','Tooltip','Before you can trace the higher order modes, their cut-off frequencies at the phase velocity limit must be detected. This is done automatically upon pressing ''calculate'' if it has not already been done manually. Sometimes the automatic search does not find all modes, then change ''Step'' and press ''Search'' to find the missing modes.','Position',[140 10 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','16');

        a.OutputWindow1aUI2 = uicontrol('Parent',Tab2,'Style','text','String','','Position',[245 50 105 630],'BackgroundColor','white');
        a.OutputWindow1bUI2 = uicontrol('Parent',Tab2,'Style','text','String','','Position',[362 50 105 630],'BackgroundColor','white');
        a.OutputWindow2aUI2 = uicontrol('Parent',Tab2,'Style','text','String','','Position',[245 10 105 28],'BackgroundColor','white');
        a.OutputWindow2bUI2 = uicontrol('Parent',Tab2,'Style','text','String','','Position',[362 10 105 28],'BackgroundColor','white');

        %------------------------------------------------------------------
        a.p15UI2 = uipanel('Parent',Tab2,'Title','Trace modes and calculate energy velocity','Units','pixels','Position',[478 635 263 120],'FontSize',10);
        a.SamplesX3ceTextUI2 = uicontrol('Parent',a.p15UI2,'Style','text','HorizontalAlignment','left','String','Samples x3','Position',[30 40 58 13]);
        a.SamplesX3ceUI2 = uicontrol('Parent',a.p15UI2,'Style','edit','String',a.SamplesX3ce2,'Tooltip','Enter the number of through-thickness samples for the energy velocity calculation. The more samples the more accurate the result. The accuracy decreases with an increasing frequency-thickness product because the mode shapes become more complex. Therefore, choose a sufficiently high number of samples.','Position',[34 15 50 23],'Callback',@CallbackUI2,'Tag','83');        
        a.TraceModesUI2 = uicontrol('Parent',a.p15UI2,'Style','togglebutton','String','Trace modes','Tooltip','Start tracing the dispersion curves.','Position',[117 55 130 40],'FontSize',10,'Callback',@CallbackUI2,'BusyAction','cancel','Tag','17');
        a.CalculateCeUI2 = uicontrol('Parent',a.p15UI2,'Style','togglebutton','String','Calculate ce','Tooltip','Start calculating the energy velocity.','Position',[117 10 130 40],'FontSize',10,'Callback',@CallbackUI2,'BusyAction','cancel','Tag','18');

        a.p8UI2 = uipanel('Parent',Tab2,'Title','Dispersion diagrams','Units','pixels','Position',[478 405 263 225],'FontSize',10);
        uicontrol('Parent',a.p8UI2,'Style','text','HorizontalAlignment','left','String','Quantity','Position',[10 185 42 13]);
        a.Option1TextUI2 = uicontrol('Parent',a.p8UI2,'Style','text','HorizontalAlignment','left','String','Bulk velocities','Position',[10 155 71 13]);
        a.XAxisModeTextUI2 = uicontrol('Parent',a.p8UI2,'Style','text','HorizontalAlignment','left','String','X-axis mode','Position',[10 125 63 13]);
        a.XAxisTextUI2 = uicontrol('Parent',a.p8UI2,'Style','text','HorizontalAlignment','left','String','X-axis (kHz)','Position',[10 95 86 13]);
        a.YAxisTextUI2 = uicontrol('Parent',a.p8UI2,'Style','text','HorizontalAlignment','left','String','Y-axis (m/ms)','Position',[10 65 101 13]);
        a.Quantity1UI2 = uicontrol('Parent',a.p8UI2,'Style','popupmenu','String',{'Phase velocity (m/ms)','Energy velocity (m/ms)',['Propagation time (',char(181),'s)'],['Coincidence angle (',char(176),')'],'Wavelength (mm)','Wavenumber (rad/mm)','Attenuation (Np/m)'},'Tooltip','Select which quantity to plot in the dispersion diagram.','Position',[117 180 130 23],'Callback',@CallbackUI2,'Tag','19');
        a.Option1UI2 = uicontrol('Parent',a.p8UI2,'Style','checkbox','Value',a.BulkVelocities2,'Tooltip','Check this to show the bulk wave velocities.','Position',[117 150 130 23],'Callback',@CallbackUI2,'Tag','20');
        a.XAxisModeUI2 = uicontrol('Parent',a.p8UI2,'Style','popupmenu','String',{'Frequency (kHz)','Frequency (MHz)',['f',char(8901),'d (MHz',char(8901),'mm)']},'Tooltip','Select the frequency''s dimension on the X-axis.','Position',[117 120 130 23],'Callback',@CallbackUI2,'Tag','21');
        a.XAxisUI2 = uicontrol('Parent',a.p8UI2,'Style','edit','String',['[0 ',num2str(a.FrequencyLimit2),']'],'Tooltip','Enter which frequency range shall be plotted.','Position',[117 90 75 23],'Callback',@CallbackUI2,'Tag','22');
        a.YAxisUI2 = uicontrol('Parent',a.p8UI2,'Style','edit','String',['[0 ',num2str(a.PhaseVelocityLimit2/1e3),']'],'Tooltip','Enter which phase velocity range shall be plotted.','Position',[117 60 75 23],'Callback',@CallbackUI2,'Tag','23');
        a.Plot1UI2 = uicontrol('Parent',a.p8UI2,'Style','pushbutton','String','Plot','Tooltip','Plot the dispersion diagram. If you have checked ''Export plots'' in the export settings, the plot will be exported automatically.','Position',[117 15 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','24');

        a.p9UI2 = uipanel('Parent',Tab2,'Title','Through-thickness profiles','Units','pixels','Position',[751 405 224 290],'FontSize',10);
        uicontrol('Parent',a.p9UI2,'Style','text','HorizontalAlignment','left','String','Quantity','Position',[10 245 42 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','HorizontalAlignment','left','String','Mode','Position',[10 215 28 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','HorizontalAlignment','left','String','Frequency (kHz)','Position',[10 185 83 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','HorizontalAlignment','left','String','Samples per layer','Position',[10 155 89 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','HorizontalAlignment','left','String','Half-spaces','Position',[10 125 61 13]);
        uicontrol('Parent',a.p9UI2,'Style','text','HorizontalAlignment','left','String','Phase','Position',[10 95 32 13]);
        a.Quantity2UI2 = uicontrol('Parent',a.p9UI2,'Style','popupmenu','String',{'Displacement','Stress','Strain','Energy density','Power flow density'},'Tooltip','Select which through-thickness quantity to plot.','Position',[119 240 90 23],'Callback',@CallbackUI2,'Tag','25');
        a.Mode1UI2 = uicontrol('Parent',a.p9UI2,'Style','popupmenu','String',{''},'Tooltip','Select the mode you want to analyze.','Position',[119 210 90 23],'Callback',@CallbackUI2,'Tag','26');
        a.Frequency1UI2 = uicontrol('Parent',a.p9UI2,'Style','edit','String',a.FrequencyLimit2,'Tooltip','Enter the frequency at which to analyze the selected mode.','Position',[119 180 50 23],'Callback',@CallbackUI2,'Tag','27');
        a.Samples1UI2 = uicontrol('Parent',a.p9UI2,'Style','edit','String',a.Samples12,'Tooltip','Enter the number of sample points per layer over the plate''s thickness (x3) at which the selected quantities are calculated.','Position',[119 150 50 23],'Callback',@CallbackUI2,'Tag','28');
        a.Halfspaces1UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','Value',a.Halfspaces12,'Tooltip','Check this to show the quantities in the upper and lower fluid.','Position',[89 120 20 23],'Enable','off','Callback',@CallbackUI2,'Tag','12');         
        a.HalfspacesNumber1UI2 = uicontrol('Parent',a.p9UI2,'Style','edit','String',a.HalfspacesNumber12,'Tooltip','Set the height of the half-spaces in plate thicknesses.','Position',[119 120 50 23],'Enable','off','Callback',@CallbackUI2,'Tag','11');
        a.PhaseUI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','Value',a.Phase2,'Tooltip','Check this to plot also the phase of the field components.','Position',[119 90 20 23],'Callback',@CallbackUI2,'Tag','81');
        a.Plot2UI2 = uicontrol('Parent',a.p9UI2,'Style','pushbutton','String','Plot','Tooltip','Plot the profile. If you have checked ''Export plots'' in the export settings, the plot will be exported automatically.','Position',[119 50 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','29');
        a.Plot2UIa2 = uicontrol('Parent',a.p9UI2,'Style','pushbutton','String','Live-plot','Tooltip','The profiles are updating as you are sliding along the dispersion curves.','Position',[119 15 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','85');
        a.x11UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','11','Value',a.Plot2(1),'Tooltip','','Position',[10 60 35 23],'Enable','off','Callback',@CallbackUI2,'Tag','30');
        a.x22UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','22','Value',a.Plot2(2),'Tooltip','','Position',[45 40 35 23],'Enable','off','Callback',@CallbackUI2,'Tag','31');
        a.x33UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','3','Value',a.Plot2(3),'Tooltip','Displacement u3','Position',[80 20 40 23],'Callback',@CallbackUI2,'Tag','32');
        a.x13UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','1','Value',a.Plot2(4),'Tooltip','Displacement u1','Position',[80 60 40 23],'Callback',@CallbackUI2,'Tag','33');
        a.x23UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','2','Value',a.Plot2(5),'Tooltip','Displacement u2','Position',[80 40 40 23],'Callback',@CallbackUI2,'Tag','34');
        a.x12UI2 = uicontrol('Parent',a.p9UI2,'Style','checkbox','String','12','Value',a.Plot2(6),'Tooltip','','Position',[45 60 35 23],'Enable','off','Callback',@CallbackUI2,'Tag','35');        
        
        a.p10UI2 = uipanel('Parent',Tab2,'Title','Mode shape','Units','pixels','Position',[478 195 497 205],'FontSize',10);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Mode, plane','Position',[10 165 60 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Frequency (kHz)','Position',[10 135 83 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Wavelengths','Position',[10 105 65 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Samples per layer','Position',[10 75 89 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Half-spaces','Position',[10 45 61 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Grid line every','Position',[10 15 72 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Gain','Position',[185 15 24 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Animate','Position',[282 45 41 13]);
        uicontrol('Parent',a.p10UI2,'Style','text','HorizontalAlignment','left','String','Undistorted','Position',[282 15 57 13]);
        a.Mode2UI2 = uicontrol('Parent',a.p10UI2,'Style','popupmenu','String',{''},'Tooltip','Select the mode you want to analyze.','Position',[117 160 90 23],'Callback',@CallbackUI2,'Tag','36');
        a.PlaneUI2 = uicontrol('Parent',a.p10UI2,'Style','popupmenu','String',{'r-z',['r-',char(952)]},'Tooltip','Select in which plane you want to plot.','Position',[215 160 48 23],'Enable','off','Callback',@CallbackUI2,'Tag','45');
        a.Frequency2UI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Frequency22,'Tooltip','Enter the frequency at which to analyze the selected mode.','Position',[117 130 50 23],'Callback',@CallbackUI2,'Tag','37');
        a.LengthUI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Length2,'Tooltip','Specify how many wavelengths you want to display.','Position',[117 100 50 23],'Callback',@CallbackUI2,'Tag','38');
        a.Samples2UI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Samples22,'Tooltip','Enter the number of sample points along the propagation direction (x1) at which the displacement is calculated.','Position',[117 70 50 23],'Callback',@CallbackUI2,'Tag','39');
        a.Samples3UI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Samples32,'Tooltip','Enter the number of sample points over the plate''s thickness (x3) at which the displacement is calculated.','Position',[177 70 50 23],'Callback',@CallbackUI2,'Tag','40');
        a.Halfspaces2UI2 = uicontrol('Parent',a.p10UI2,'Style','checkbox','Value',a.Halfspaces22,'Tooltip','Check this to show the quantities in the upper and lower fluid.','Position',[87 40 50 23],'Enable','off','Callback',@CallbackUI2,'Tag','8');
        a.HalfspacesNumber2UI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.HalfspacesNumber22,'Tooltip','Set the height of the half-spaces in plate thicknesses.','Position',[117 40 50 23],'Enable','off','Callback',@CallbackUI2,'Tag','6');        
        a.GridLineUI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.GridLine2,'Tooltip','Draw a grid line at every ith sample point.','Position',[117 10 50 23],'Callback',@CallbackUI2,'Tag','42');
        a.GainUI2 = uicontrol('Parent',a.p10UI2,'Style','edit','String',a.Gain2,'Tooltip','Magnify the displacement.','Position',[220 10 50 23],'Callback',@CallbackUI2,'Tag','41');
        a.AnimateUI2 = uicontrol('Parent',a.p10UI2,'Style','checkbox','Value',a.Animate2,'Tooltip','Check this in order to show the animated mode shape upon pressing the plot button below.','Position',[352 40 70 23],'Callback',@CallbackUI2,'Tag','49');
        a.UndistortedUI2 = uicontrol('Parent',a.p10UI2,'Style','checkbox','Value',a.Undistorted2,'Tooltip','Check this to draw the undistorted grid.','Position',[352 10 50 23],'Callback',@CallbackUI2,'Tag','43');
        a.Plot3UI2 = uicontrol('Parent',a.p10UI2,'Style','pushbutton','String','Plot','Tooltip','Plot the mode shape. If you have checked ''Animate'', the mode shape will be animated. If you have checked ''Export plots'' in the export settings, the plot/movie will be exported automatically.','Position',[392 15 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','44');

        a.p11UI2 = uipanel('Parent',Tab2,'Title','Animation settings','Units','pixels','Position',[751 265 212 115],'FontSize',9);
        uicontrol('Parent',a.p11UI2,'Style','text','HorizontalAlignment','left','String','Cycle duration (s)','Position',[10 75 88 13]);
        uicontrol('Parent',a.p11UI2,'Style','text','HorizontalAlignment','left','String','Frame rate (Hz)','Position',[10 45 78 13]);
        uicontrol('Parent',a.p11UI2,'Style','text','HorizontalAlignment','left','String','Movie quality (0-100)','Position',[10 15 103 13]);
        a.CycleDurationUI2 = uicontrol('Parent',a.p11UI2,'Style','edit','String',a.CycleDuration2,'Tooltip','Enter how long a cycle shall take.','Position',[119 70 50 23],'Callback',@CallbackUI2,'Tag','46');
        a.FrameRateUI2 = uicontrol('Parent',a.p11UI2,'Style','edit','String',a.FrameRate2,'Tooltip','Define the frame rate of the movie.','Position',[119 40 50 23],'Callback',@CallbackUI2,'Tag','47');
        a.MovieQualityUI2 = uicontrol('Parent',a.p11UI2,'Style','edit','String',a.MovieQuality2,'Tooltip','Define the quality of the exported movie.','Position',[119 10 50 23],'Callback',@CallbackUI2,'Tag','48');

        a.p12UI2 = uipanel('Parent',Tab2,'Title','Export settings','Units','pixels','Position',[478 10 497 180],'FontSize',10);
        uicontrol('Parent',a.p12UI2,'Style','text','HorizontalAlignment','left','String','Export plots','Position',[10 140 59 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','HorizontalAlignment','left','String','PDF','Position',[124 140 26 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','HorizontalAlignment','left','String','PNG','Position',[124 120 28 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','HorizontalAlignment','left','String','PNG resolution (dpi)','Position',[10 93 98 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','HorizontalAlignment','left','String','Dispersion curves','Position',[240 140 90 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','HorizontalAlignment','left','String','Through-thickness','Position',[240 120 92 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','HorizontalAlignment','left','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p12UI2,'Style','text','HorizontalAlignment','left','String','Directory','Position',[10 20 46 13]);
        a.ExportPlotsUI2 = uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.ExportPlots2,'Tooltip','Check this in order to export plots/movies automatically upon pressing the respective plot button. You can also export plots manually by using the ''File'' menu inside the plot figure.','Position',[80 135 20 23],'Callback',@CallbackUI2,'Tag','50');
        a.PDFUI2 = uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.PDF2,'Tooltip','Check this to export a plot as pdf.','Position',[160 135 20 23],'Callback',@CallbackUI2,'Tag','51');
        a.PNGUI2 = uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.PNG2,'Tooltip','Check this to export a plot as png.','Position',[160 115 20 23],'Callback',@CallbackUI2,'Tag','52');
        a.PNGresolutionUI2 = uicontrol('Parent',a.p12UI2,'Style','edit','String',a.PNGresolution2,'Tooltip','Enter the resolution of the png image.','Position',[125 88 50 23],'Callback',@CallbackUI2,'Tag','53');
        a.XAxisMode2UI2 = uicontrol('Parent',a.p12UI2,'Style','popupmenu','String',{'Frequency (kHz)','Frequency (MHz)',['f',char(8901),'d (MHz',char(8901),'mm)']},'Tooltip','Select the frequency''s dimension to be exported.','Position',[370 130 110 23],'Callback',@CallbackUI2,'Tag','54');
        a.ArrangeUI2 = uicontrol('Parent',a.p12UI2,'Style','popupmenu','String',{'Horizontal arrangement','Vertical arrangement'},'Tooltip','Choose how to arrange the dispersion curve data in the Excel/txt.','Position',[370 95 110 23],'Callback',@CallbackUI2,'Tag','55');        
        a.DispersionCurvesUI2 = uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.DispersionCurves2,'Tooltip','Check this to export the dispersion curves.','Position',[340 135 20 23],'Callback',@CallbackUI2,'Tag','56');
        a.ThroughThicknessUI2 = uicontrol('Parent',a.p12UI2,'Style','checkbox','Value',a.ThroughThickness2,'Tooltip','Check this to export the through-thickness profiles selected above.','Position',[340 115 20 23],'Callback',@CallbackUI2,'Tag','57');        
        uicontrol('Parent',a.p12UI2,'Style','pushbutton','String','*.mat','Tooltip','Export the data as Matlab''s mat-file.','Position',[220 80 40 33],'Callback',@CallbackUI2,'Tag','80');
        uicontrol('Parent',a.p12UI2,'Style','pushbutton','String','*.xlsx','Tooltip','Export the data as Excel sheet.','Position',[270 80 40 33],'Callback',@CallbackUI2,'Tag','58');
        uicontrol('Parent',a.p12UI2,'Style','pushbutton','String','*.txt','Tooltip','Export the data as txt-file.','Position',[320 80 40 33],'Callback',@CallbackUI2,'Tag','59');      
        a.FileNameUI2 = uicontrol('Parent',a.p12UI2,'Style','edit','String',a.FileName2,'Tooltip','Specify the name of the plots/movies to be exported. The dispersion curve raw data are named automatically.','Position',[70 45 412 23],'Callback',@CallbackUI2,'Tag','60');
        a.DirectoryUI2 = uicontrol('Parent',a.p12UI2,'Style','edit','String',a.Directory,'Tooltip','Specify the directory to which plots, movies, and the dispersion curve raw data shall be exported.','Position',[70 15 412 23],'Callback',@CallbackUI2,'Tag','61');

        %------------------------------------------------------------------
        a.p13UI2 = uipanel('Parent',Tab2,'Title','Plot layout settings','Units','pixels','Position',[985 275 200 420],'FontSize',10);
        uicontrol('Parent',a.p13UI2,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 380 21 13]);
        uicontrol('Parent',a.p13UI2,'Style','text','HorizontalAlignment','left','String','Legend location','Position',[10 350 78 13]);
        uicontrol('Parent',a.p13UI2,'Style','text','HorizontalAlignment','left','String','Box line width','Position',[10 320 70 13]);
        uicontrol('Parent',a.p13UI2,'Style','text','HorizontalAlignment','left','String','Curve line width','Position',[10 290 80 13]);
        a.TitleUI2 = uicontrol('Parent',a.p13UI2,'Style','popupmenu','String',{'with layup','without layup','no title'},'Tooltip','Choose to plot the title with the layup  included, without layup, or no title at all.','Position',[105 375 80 23],'Callback',@CallbackUI2,'Tag','62');
        a.LegendLocationUI2 = uicontrol('Parent',a.p13UI2,'Style','popupmenu','String',{'outside','inside'},'Tooltip','Determine the legend location in the through-thickness plots.','Position',[105 345 80 23],'Callback',@CallbackUI2,'Tag','64');
        a.BoxLineWidthUI2 = uicontrol('Parent',a.p13UI2,'Style','edit','String',a.BoxLineWidth2,'Tooltip','Enter the box line width.','Position',[105 315 50 23],'Callback',@CallbackUI2,'Tag','65');
        a.LineWidthUI2 = uicontrol('Parent',a.p13UI2,'Style','edit','String',a.LineWidth2,'Tooltip','Enter the curve line width.','Position',[105 285 50 23],'Callback',@CallbackUI2,'Tag','66');

        a.p14UI2 = uipanel('Parent',Tab2,'Title','Dispersion curve colors [R G B]','Units','pixels','Position',[995 435 180 115],'FontSize',9);
        a.SColorTextUI2 = uicontrol('Parent',a.p14UI2,'Style','text','HorizontalAlignment','left','String','S','Position',[10 75 17 13]);
        a.AColorTextUI2 = uicontrol('Parent',a.p14UI2,'Style','text','HorizontalAlignment','left','String','A','Position',[10 45 8 13]);
        a.BColorTextUI2 = uicontrol('Parent',a.p14UI2,'Style','text','HorizontalAlignment','left','String','B','Position',[10 15 7 13]);
        a.SColorUI2 = uicontrol('Parent',a.p14UI2,'Style','edit','String','[1 0 0]','Tooltip','Specify the color of symmetric modes. You can compose any color by entering the corresponding RGB values. Set numbers from 0 to 1.','Position',[95 70 50 23],'Callback',@CallbackUI2,'Tag','67');
        a.AColorUI2 = uicontrol('Parent',a.p14UI2,'Style','edit','String','[0 0 1]','Tooltip','Specify the color of antisymmetric modes.','Position',[95 40 50 23],'Callback',@CallbackUI2,'Tag','68');
        a.BColorUI2 = uicontrol('Parent',a.p14UI2,'Style','edit','String','[.5 0 1]','Tooltip','Specify the color of nonsymmetric modes.','Position',[95 10 50 23],'Callback',@CallbackUI2,'Tag','69');

        a.p16UI2 = uipanel('Parent',Tab2,'Title','Font size','Units','pixels','Position',[995 285 180 145],'FontSize',9);
        uicontrol('Parent',a.p16UI2,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 105 21 13]);
        uicontrol('Parent',a.p16UI2,'Style','text','HorizontalAlignment','left','String','Axes labels','Position',[10 75 59 13]);
        uicontrol('Parent',a.p16UI2,'Style','text','HorizontalAlignment','left','String','Axes ticks','Position',[10 45 53 13]);
        uicontrol('Parent',a.p16UI2,'Style','text','HorizontalAlignment','left','String','Legend','Position',[10 15 38 13]);
        a.TitleFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.TitleFontSize2,'Tooltip','Set the title font size.','Position',[95 100 50 23],'Callback',@CallbackUI2,'Tag','73');
        a.AxesLabelFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.AxesLabelFontSize2,'Tooltip','Set the axes label font size.','Position',[95 70 50 23],'Callback',@CallbackUI2,'Tag','74');
        a.AxesTickFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.AxesTickFontSize2,'Tooltip','Set the axes tick label font size.','Position',[95 40 50 23],'Callback',@CallbackUI2,'Tag','75');
        a.LegendFontSizeUI2 = uicontrol('Parent',a.p16UI2,'Style','edit','String',a.LegendFontSize2,'Tooltip','Set the legend font size.','Position',[95 10 50 23],'Callback',@CallbackUI2,'Tag','77');

        a.c3UI2 = uicontrol('Parent',Tab2,'Style','pushbutton','String','Default','Tooltip','Reset the plot layout settings to the default.','Position',[1052 230 65 33],'FontSize',10,'Callback',@CallbackUI2,'Tag','78');    

        uicontrol('Parent',Tab2,'Style','text','HorizontalAlignment','left','String','CPU cores','Position',[1035 715 53 13]);
        a.c4UI2 = uicontrol('Parent',Tab2,'Style','text','String','0','Position',[1095 715 30 13],'BackgroundColor','white');
        a.a1UI2 = axes('Parent',Tab2,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading20.png') 
    end
    if  1 % Tab3_signal simulator
        a.DataUI3 = uicontrol('Parent',Tab3,'Style','text','String','','Tooltip','The material to which the dispersion curves belonge to.','Position',[10 740 279 13],'BackgroundColor','white');
        
        a.p1UI3 = uipanel('Parent',Tab3,'Title','Computational settings','Units','pixels','Position',[10 420 279 315],'FontSize',10);
        a.WaveGeneratorUI3 = uicontrol('Parent',a.p1UI3,'Style','pushbutton','String','Generate wave','Tooltip','Press to generate a wave form.','Position',[140 226 127 40],'Enable','off','FontSize',10,'Callback',@WaveGeneratorUI3_Callback);
        
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String','Frequency','Position',[10 273 60 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String','Cycles','Position',[10 253 40 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String','Samples/cycle','Position',[10 233 80 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String','Window','Position',[10 213 50 13]);
        a.FrequencyDisplayUI3 = uicontrol('Parent',a.p1UI3,'Style','text','String','','Tooltip','The carrier wave frequency in the wave packet and the chirp type (if selected).','Position',[70 273 197 13],'BackgroundColor','white');
        a.CyclesDisplayUI3 = uicontrol('Parent',a.p1UI3,'Style','text','String',a.Cycles3,'Tooltip','The number of cycles of the signal within the wave packet.','Position',[90 253 40 13],'BackgroundColor','white');
        a.SamplesPerCycleDisplayUI3 = uicontrol('Parent',a.p1UI3,'Style','text','String',a.SamplesPerCycle3,'Tooltip','The samples per cycle in the excitation signal.','Position',[90 233 40 13],'BackgroundColor','white');
        a.WindowDisplayUI3 = uicontrol('Parent',a.p1UI3,'Style','text','String','Gauss','Tooltip','The window function to be multiplied with the carrier wave.','Position',[60 213 70 13],'BackgroundColor','white');
        
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String','Distance (mm)','Position',[10 185 71 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String',['n',char(8901),'Distance/ce'],'Position',[10 155 70 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String',['Spectral threshold (',char(37),')'],'Position',[10 125 111 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String','Displacement component','Position',[10 95 122 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String',['Gate (',char(181),'s)'],'Position',[10 65 49 13]);
        uicontrol('Parent',a.p1UI3,'Style','text','HorizontalAlignment','left','String','Multi-mode','Position',[10 25 53 13]);
        a.DistanceUI3 = uicontrol('Parent',a.p1UI3,'Style','edit','String',a.Distance3,'Tooltip','Enter the propagation distance of the wave packet.','Position',[140 180 50 23],'Callback',@CallbackUI3,'Tag','6');
        a.TimeLimitFactorUI3 = uicontrol('Parent',a.p1UI3,'Style','edit','String',a.TimeLimitFactor3,'Tooltip','Enter the temporal limit of the simulation in multiples ''n'' of the propagation distance divided through the energy velocity of the mode in question. In case of multiple modes, the lowest energy velocity among all contributing modes is used. Due to dispersion, the wave packets spread with the propagation time so that the value should always be greater than one.','Position',[140 150 50 23],'Callback',@CallbackUI3,'Tag','7');
        a.SpectrumThresholdUI3 = uicontrol('Parent',a.p1UI3,'Style','edit','String',a.SpectrumThreshold3,'Tooltip','Enter the threshold of the spectral amplitudes, which shall be used for the construction of the propagated wave packet. The value is in percent of the maximal spectral amplitude of the excitation signal at the center frequency defined above. A smaller value widens the frequency range to be taken into acount.','Position',[140 120 50 23],'Callback',@CallbackUI3,'Tag','8');
        a.DisplacementComponentUI3 = uicontrol('Parent',a.p1UI3,'Style','popupmenu','String',{'Out-of-plane (u3)','In-plane (u1/2)'},'Tooltip','Select which displacement component you want to calculate.','Position',[140 90 120 23],'Callback',@CallbackUI3,'Tag','113');
        a.GateUI3 = uicontrol('Parent',a.p1UI3,'Style','edit','String',['[',num2str(a.Gate3(1)),' ',num2str(a.Gate3(2)),']'],'Tooltip','Enter the gate for the calculation of the propagated spectrum. FFT will be performed on the temporal response between the gate''s limits. Set the limits left and right of the wave packet you want to consider.','Position',[140 60 75 23],'Callback',@CallbackUI3,'Tag','111');
        a.MultiModeUI3 = uicontrol('Parent',a.p1UI3,'Style','checkbox','Value',a.MultiMode3,'Tooltip','Check this to simulate multiple modes at a time.','Position',[100 20 20 23],'Callback',@CallbackUI3,'Tag','9');
        a.CalculateUI3 = uicontrol('Parent',a.p1UI3,'Style','togglebutton','String','Calculate','Tooltip','Calculate multiple modes at a time.','Position',[140 15 75 33],'Enable','off','FontSize',10,'Callback',@CallbackUI3,'BusyAction','cancel','Tag','92');

        a.p2UI3 = uipanel('Parent',Tab3,'Title','Mode selection','Units','pixels','Position',[10 185 279 230],'FontSize',10);
        a.ALamb0aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A0','Position',[10 190 16 13],'Visible','off');
        a.ALamb0bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 188 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','10');
        a.ALamb0cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 190 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','11');
        a.ALamb1aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A1','Position',[10 170 16 13],'Visible','off');
        a.ALamb1bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 168 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','12');
        a.ALamb1cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 170 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','13');
        a.ALamb2aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A2','Position',[10 150 16 13],'Visible','off');
        a.ALamb2bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 148 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','14');
        a.ALamb2cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 150 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','15');
        a.ALamb3aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A3','Position',[10 130 16 13],'Visible','off');
        a.ALamb3bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 128 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','16');
        a.ALamb3cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 130 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','17');
        a.ALamb4aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A4','Position',[10 110 16 13],'Visible','off');
        a.ALamb4bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 108 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','18');
        a.ALamb4cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 110 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','19');
        a.ALamb5aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A5','Position',[10 90 16 13],'Visible','off');
        a.ALamb5bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 88 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','20');
        a.ALamb5cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 90 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','21');
        a.ALamb6aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A6','Position',[10 70 16 13],'Visible','off');
        a.ALamb6bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 68 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','22');
        a.ALamb6cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 70 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','23');
        a.ALamb7aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A7','Position',[10 50 16 13],'Visible','off');
        a.ALamb7bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 48 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','24');
        a.ALamb7cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 50 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','25');
        a.ALamb8aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A8','Position',[10 30 16 13],'Visible','off');
        a.ALamb8bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 28 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','26');
        a.ALamb8cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 30 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','27');
        a.ALamb9aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','A9','Position',[10 10 16 13],'Visible','off');
        a.ALamb9bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[28 8 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','28');
        a.ALamb9cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[43 10 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','29');        

        a.SLamb0aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S0','Position',[70 190 16 13],'Visible','off');
        a.SLamb0bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 188 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','30');
        a.SLamb0cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 190 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','31');
        a.SLamb1aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S1','Position',[70 170 16 13],'Visible','off');
        a.SLamb1bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 168 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','32');
        a.SLamb1cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 170 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','33');
        a.SLamb2aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S2','Position',[70 150 16 13],'Visible','off');
        a.SLamb2bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 148 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','34');
        a.SLamb2cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 150 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','35');
        a.SLamb3aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S3','Position',[70 130 16 13],'Visible','off');
        a.SLamb3bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 128 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','36');
        a.SLamb3cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 130 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','37');
        a.SLamb4aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S4','Position',[70 110 16 13],'Visible','off');
        a.SLamb4bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 108 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','38');
        a.SLamb4cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 110 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','39');
        a.SLamb5aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S5','Position',[70 90 16 13],'Visible','off');
        a.SLamb5bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 88 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','40');
        a.SLamb5cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 90 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','41');
        a.SLamb6aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S6','Position',[70 70 16 13],'Visible','off');
        a.SLamb6bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 68 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','42');
        a.SLamb6cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 70 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','43');
        a.SLamb7aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S7','Position',[70 50 16 13],'Visible','off');
        a.SLamb7bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 48 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','44');
        a.SLamb7cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 50 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','45');
        a.SLamb8aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S8','Position',[70 30 16 13],'Visible','off');
        a.SLamb8bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 28 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','46');
        a.SLamb8cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 30 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','47');
        a.SLamb9aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','S9','Position',[70 10 16 13],'Visible','off');
        a.SLamb9bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[88 8 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','48');
        a.SLamb9cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[103 10 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','49');

        a.AShear1aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH1','Position',[129 190 30 13],'Visible','off');
        a.AShear1bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 188 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','50');
        a.AShear1cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 190 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','51'); 
        a.AShear2aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH2','Position',[129 170 30 13],'Visible','off');
        a.AShear2bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 168 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','52');
        a.AShear2cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 170 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','53');
        a.AShear3aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH3','Position',[129 150 30 13],'Visible','off');
        a.AShear3bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 148 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','54');
        a.AShear3cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 150 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','55');
        a.AShear4aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH4','Position',[129 130 30 13],'Visible','off');
        a.AShear4bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 128 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','56');
        a.AShear4cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 130 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','57');
        a.AShear5aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH5','Position',[129 110 30 13],'Visible','off');
        a.AShear5bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 108 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','58');
        a.AShear5cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 110 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','59');
        a.AShear6aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH6','Position',[129 90 30 13],'Visible','off');
        a.AShear6bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 88 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','60');
        a.AShear6cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 90 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','61');
        a.AShear7aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH7','Position',[129 70 30 13],'Visible','off');
        a.AShear7bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 68 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','62');
        a.AShear7cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 70 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','63');
        a.AShear8aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH8','Position',[129 50 30 13],'Visible','off');
        a.AShear8bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 48 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','64');
        a.AShear8cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 50 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','65');
        a.AShear9aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH9','Position',[129 30 30 13],'Visible','off');
        a.AShear9bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 28 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','66');
        a.AShear9cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 30 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','67');
        a.AShear10aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','ASH10','Position',[125 10 36 13],'Visible','off');
        a.AShear10bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[160 8 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','68');
        a.AShear10cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[175 10 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','69');        
        
        a.SShear0aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH0','Position',[201 190 30 13],'Visible','off');
        a.SShear0bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 188 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','70');
        a.SShear0cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 190 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','71');
        a.SShear1aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH1','Position',[201 170 30 13],'Visible','off');
        a.SShear1bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 168 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','72');
        a.SShear1cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 170 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','73');
        a.SShear2aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH2','Position',[201 150 30 13],'Visible','off');
        a.SShear2bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 148 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','74');
        a.SShear2cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 150 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','75');
        a.SShear3aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH3','Position',[201 130 30 13],'Visible','off');
        a.SShear3bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 128 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','76');
        a.SShear3cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 130 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','77');
        a.SShear4aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH4','Position',[201 110 30 13],'Visible','off');
        a.SShear4bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 108 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','78');
        a.SShear4cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 110 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','79');
        a.SShear5aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH5','Position',[201 90 30 13],'Visible','off');
        a.SShear5bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 88 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','80');
        a.SShear5cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 90 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','81');
        a.SShear6aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH6','Position',[201 70 30 13],'Visible','off');
        a.SShear6bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 68 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','82');
        a.SShear6cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 70 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','83');
        a.SShear7aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH7','Position',[201 50 30 13],'Visible','off');
        a.SShear7bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 48 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','84');
        a.SShear7cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 50 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','85');
        a.SShear8aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH8','Position',[201 30 30 13],'Visible','off');
        a.SShear8bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 28 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','86');
        a.SShear8cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 30 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','87');
        a.SShear9aUI3 = uicontrol('Parent',a.p2UI3,'Style','text','HorizontalAlignment','left','String','SSH9','Position',[201 10 30 13],'Visible','off');
        a.SShear9bUI3 = uicontrol('Parent',a.p2UI3,'Style','checkbox','Value',0,'Position',[232 8 15 17],'Visible','off','Callback',@CallbackUI3,'Tag','88');
        a.SShear9cUI3 = uicontrol('Parent',a.p2UI3,'Style','edit','String','1','Tooltip','Scale the amplitude.','Position',[247 10 20 13],'Visible','off','Callback',@CallbackUI3,'Tag','89');
        
        a.OutputWindowUI3 = uicontrol('Parent',Tab3,'Style','text','String','','Position',[10 10 195 165],'BackgroundColor','white');
        a.c1UI3 = uicontrol('Parent',Tab3,'Style','text','HorizontalAlignment','left','String',['X-axis (',char(181),'s)'],'Position',[220 165 58 13]);
        a.c2UI3 = uicontrol('Parent',Tab3,'Style','text','HorizontalAlignment','left','String','Y-axis (nm)','Position',[220 115 59 13]);
        a.XAxisUI3 = uicontrol('Parent',Tab3,'Style','edit','String','','Tooltip','Enter which frequency range shall be plotted.','Position',[218 135 65 23],'Callback',@CallbackUI3,'Tag','90');
        a.YAxisUI3 = uicontrol('Parent',Tab3,'Style','edit','String','','Tooltip','Enter which displacement range shall be plotted.','Position',[218 85 65 23],'Callback',@CallbackUI3,'Tag','112'); 
        a.PlotUI3 = uicontrol('Parent',Tab3,'Style','pushbutton','String','Plot','Tooltip','Plot the temporal and frequency responses. If you have checked ''Export plots''in the export settings, the plot will be exported automatically.','Position',[218 40 65 33],'FontSize',10,'Callback',@CallbackUI3,'Tag','91');
        
        %------------------------------------------------------------------
        a.p3UI3 = uipanel('Parent',Tab3,'Title','Export settings','Units','pixels','Position',[296 10 446 180],'FontSize',10);
        uicontrol('Parent',a.p3UI3,'Style','text','HorizontalAlignment','left','String','Export plots','Position',[10 140 59 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','HorizontalAlignment','left','String','PDF','Position',[124 140 26 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','HorizontalAlignment','left','String','PNG','Position',[124 120 28 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','HorizontalAlignment','left','String','PNG resolution (dpi)','Position',[10 93 98 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','HorizontalAlignment','left','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p3UI3,'Style','text','HorizontalAlignment','left','String','Directory','Position',[10 20 46 13]);
        a.ExportPlotsUI3 = uicontrol('Parent',a.p3UI3,'Style','checkbox','Value',a.ExportPlots3,'Tooltip','Check this in order to export plots automatically upon pressing the plot button. You can also export plots manually by using the ''File'' menu inside the plot figure.','Position',[80 135 20 23],'Callback',@CallbackUI3,'Tag','93');
        a.PDFUI3 = uicontrol('Parent',a.p3UI3,'Style','checkbox','Value',a.PDF3,'Tooltip','Check this to export a plot as pdf.','Position',[160 135 20 23],'Callback',@CallbackUI3,'Tag','94');
        a.PNGUI3 = uicontrol('Parent',a.p3UI3,'Style','checkbox','Value',a.PNG3,'Tooltip','Check this to export a plot as png.','Position',[160 115 20 23],'Callback',@CallbackUI3,'Tag','95');
        a.PNGresolutionUI3 = uicontrol('Parent',a.p3UI3,'Style','edit','String',a.PNGresolution3,'Tooltip','Enter the resolution of the png image.','Position',[125 88 50 23],'Callback',@CallbackUI3,'Tag','96');
        uicontrol('Parent',a.p3UI3,'Style','pushbutton','String','*.mat','Tooltip','Export the data as Matlab''s mat-file.','Position',[220 80 40 33],'Callback',@CallbackUI3,'Tag','97');
        uicontrol('Parent',a.p3UI3,'Style','pushbutton','String','*.xlsx','Tooltip','Export the data as Excel sheet.','Position',[270 80 40 33],'Callback',@CallbackUI3,'Tag','98');
        uicontrol('Parent',a.p3UI3,'Style','pushbutton','String','*.txt','Tooltip','Export the data as txt-file.','Position',[320 80 40 33],'Callback',@CallbackUI3,'Tag','99');
        a.FileNameUI3 = uicontrol('Parent',a.p3UI3,'Style','edit','String',a.FileName3,'Tooltip','Specify the name of the plot to be exported. The signal''s raw data are named automatically.','Position',[70 45 361 23],'Callback',@CallbackUI3,'Tag','100');
        a.DirectoryUI3 = uicontrol('Parent',a.p3UI3,'Style','edit','String',a.Directory,'Tooltip','Specify the directory to which the plot and the signal''s raw data shall be exported.','Position',[70 15 361 23],'Callback',@CallbackUI3,'Tag','101');

        %------------------------------------------------------------------
        a.p4UI3 = uipanel('Parent',Tab3,'Title','Plot layout settings','Units','pixels','Position',[750 10 360 180],'FontSize',10);
        uicontrol('Parent',a.p4UI3,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 115 21 13]);
        uicontrol('Parent',a.p4UI3,'Style','text','HorizontalAlignment','left','String','Box line width','Position',[10 85 70 13]);
        uicontrol('Parent',a.p4UI3,'Style','text','HorizontalAlignment','left','String','Curve line width','Position',[10 55 80 13]);
        a.TitleUI3 = uicontrol('Parent',a.p4UI3,'Style','checkbox','Value',a.Title3,'Tooltip','Check this in order to show the plot title.','Position',[105 110 65 23],'Callback',@CallbackUI3,'Tag','102');
        a.BoxLineWidthUI3 = uicontrol('Parent',a.p4UI3,'Style','edit','String',a.BoxLineWidth3,'Tooltip','Enter the box line width.','Position',[105 80 50 23],'Callback',@CallbackUI3,'Tag','103');
        a.LineWidthUI3 = uicontrol('Parent',a.p4UI3,'Style','edit','String',a.LineWidth3,'Tooltip','Enter the curve line width.','Position',[105 50 50 23],'Callback',@CallbackUI3,'Tag','104');

        %------------------------------------------------------------------
        a.p5UI3 = uipanel('Parent',Tab3,'Title','Font size','Units','pixels','Position',[935 50 160 115],'FontSize',9);
        uicontrol('Parent',a.p5UI3,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 75 21 13]);
        uicontrol('Parent',a.p5UI3,'Style','text','HorizontalAlignment','left','String','Axes labels','Position',[10 45 59 13]);
        uicontrol('Parent',a.p5UI3,'Style','text','HorizontalAlignment','left','String','Axes ticks','Position',[10 15 53 13]);
        a.TitleFontSizeUI3 = uicontrol('Parent',a.p5UI3,'Style','edit','String',a.TitleFontSize3,'Tooltip','Set the title font size.','Position',[95 70 50 23],'Callback',@CallbackUI3,'Tag','105');
        a.AxesLabelFontSizeUI3 = uicontrol('Parent',a.p5UI3,'Style','edit','String',a.AxesLabelFontSize3,'Tooltip','Set the axes label font size.','Position',[95 40 50 23],'Callback',@CallbackUI3,'Tag','106');
        a.AxesTickFontSizeUI3 = uicontrol('Parent',a.p5UI3,'Style','edit','String',a.AxesTickFontSize3,'Tooltip','Set the axes tick label font size.','Position',[95 10 50 23],'Callback',@CallbackUI3,'Tag','107');

        a.c3UI3 = uicontrol('Parent',Tab3,'Style','pushbutton','String','Default','Tooltip','Reset the plot layout settings to the default.','Position',[1120 149 65 33],'FontSize',10,'Callback',@CallbackUI3,'Tag','108');                

        [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
    end
    if  1 % Tab4_polar diagrams
        a.p1UI4 = uipanel('Parent',Tab4,'Title','Specimen','Units','pixels','Position',[10 555 225 200],'FontSize',10);
        uicontrol('Parent',a.p1UI4,'Style','pushbutton','String','Edit','Tooltip','Press to define your specimen.','Position',[58 135 127 40],'FontSize',10,'Callback',@SpecimenSettingsUI4_Callback);
        uicontrol('Parent',a.p1UI4,'Style','text','HorizontalAlignment','left','String','Material','Position',[10 105 39 13]);
        uicontrol('Parent',a.p1UI4,'Style','text','HorizontalAlignment','left','String','Layup','Position',[10 75 32 13]);
        uicontrol('Parent',a.p1UI4,'Style','text','HorizontalAlignment','left','String','Layers','Position',[10 45 36 13]);
        uicontrol('Parent',a.p1UI4,'Style','text','HorizontalAlignment','left','String','Thickness (mm)','Position',[10 15 78 13]);
        a.MaterialNameUI4 = uicontrol('Parent',a.p1UI4,'Style','text','String',a.Material_Polar{1}.Name,'Position',[58 100 152 23],'BackgroundColor','white');
        a.LayupUI4 = uicontrol('Parent',a.p1UI4,'Style','text','String','[0]','Position',[58 70 152 23],'BackgroundColor','white');
        a.LayerCountUI4 = uicontrol('Parent',a.p1UI4,'Style','text','String','1','Tooltip','Number of layers in the layup.','Position',[160 40 50 23],'BackgroundColor','white');
        a.ThicknessCountUI4 = uicontrol('Parent',a.p1UI4,'Style','text','String',a.PlateThickness_Polar,'Tooltip','Overall thickness of the layup.','Position',[160 10 50 23],'BackgroundColor','white');

        a.p2UI4 = uipanel('Parent',Tab4,'Title','Computational settings','Units','pixels','Position',[10 315 225 235],'FontSize',10);
        uicontrol('Parent',a.p2UI4,'Style','text','HorizontalAlignment','left','String','Frequency limit (kHz)','Position',[10 195 103 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','HorizontalAlignment','left','String','Frequency step (kHz)','Position',[10 165 107 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','HorizontalAlignment','left','String','Phase velocity accuracy (m/s)','Position',[10 135 149 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','HorizontalAlignment','left','String',['Propagation angle limit (',char(176),')'],'Position',[10 105 123 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','HorizontalAlignment','left','String',['Propagation angle step (',char(176),')'],'Position',[10 75 127 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','HorizontalAlignment','left','String','Phase velocity sections (5^x)','Position',[10 45 144 13]);
        uicontrol('Parent',a.p2UI4,'Style','text','HorizontalAlignment','left','String','Matrix methods','Position',[10 15 75 13]);
        a.FrequencyLimitUI4 = uicontrol('Parent',a.p2UI4,'Style','edit','String',a.FrequencyLimit_Polar,'Tooltip','Enter the frequency up to which the dispersion curves shall be traced.','Position',[160 190 50 23],'Callback',@CallbackUI4,'Tag','1');
        a.FrequencyResolutionUI4 = uicontrol('Parent',a.p2UI4,'Style','edit','String',a.FrequencyResolution_Polar,'Tooltip','Enter the frequency step.','Position',[160 160 50 23],'Callback',@CallbackUI4,'Tag','2');
        a.AccuracyUI4 = uicontrol('Parent',a.p2UI4,'Style','edit','String',a.Accuracy_Polar,'Tooltip','Bisections are performed until the modal solutions have reached this accuracy.','Position',[160 130 50 23],'Callback',@CallbackUI4,'Tag','34');
        a.PropagationAngleModeUI4 = uicontrol('Parent',a.p2UI4,'Style','popupmenu','String',{'180','90'},'Tooltip',['Select the propagation angle range. In case of 180 ',char(176),', it will be calculated from 0 to 180 ',char(176),', and the result be copied to 180 to 360 ',char(176),'. In case of 90 ',char(176),', it will be calculated from 0 to 90 ',char(176),', and the result be copied to the other three quadrants. Notice that the latter approach is insufficient for unit cells such as [0 45].'],'Position',[160 100 50 23],'Callback',@CallbackUI4,'Tag','33');
        a.PropagationAngleStepUI4 = uicontrol('Parent',a.p2UI4,'Style','edit','String',a.PropagationAngleStep_Polar,'Tooltip',['The polar dispersion curves are calculated for propagation angles ranging from 0 to 90',char(176),'. Set the step.'],'Position',[160 70 50 23],'Callback',@CallbackUI4,'Tag','3');
        a.PhaseVelocitySectionsUI4 = uicontrol('Parent',a.p2UI4,'Style','edit','String',a.PhaseVelocitySections_Polar,'Tooltip','The fundamental dispersion curves are determined by sweeping the phase velocity at fixed frequencies. The phase velocity sections as the power of five give the maximum number of sections into which the phase velocity search interval is devided during the search for the modal solution at a given frequency. A higher number increases the chance to find the solution at the cost of processing time. Missing a solution is not necessarily critical since an extrapolation routine replaces missing samples successfully, as long as not too many samples are missing.','Position',[160 40 50 23],'Callback',@CallbackUI4,'Tag','4');
        a.MatrixMethodsUI4 = uicontrol('Parent',a.p2UI4,'Style','popupmenu','String',{'SMM+TMM','SMM'},'Tooltip','Choose if the transfer matrix method (TMM) shall be used instead of the stiffness matrix method (SMM) whenever there are stable conditions for TMM. Stable conditions for TMM are above a diagonal line starting in the origin of the frequency- phase velocity dispersion diagram. SMM is used below the diagonal. DC determines the slope of the diagonal before starting dispersion curve tracing. TMM is faster and more sensitive to quasi-shear horizontal modes than SMM. Another advantage is that its dispersion equation amplitude does not have maxima, unlike the case with SMM. This makes finding the minima (modal solutions) easier. However, sometimes SMM alone brings better results. Try both options in difficult cases.','Position',[120 10 90 23],'Callback',@CallbackUI4,'Tag','84');

        a.p3UI4 = uipanel('Parent',Tab4,'Title','Mode selection','Units','pixels','Position',[10 195 225 115],'FontSize',10);
        uicontrol('Parent',a.p3UI4,'Style','text','HorizontalAlignment','left','String','S1/B2','Position',[10 75 32 13]);
        uicontrol('Parent',a.p3UI4,'Style','text','HorizontalAlignment','left','String','S0/B1','Position',[10 45 32 13]);
        uicontrol('Parent',a.p3UI4,'Style','text','HorizontalAlignment','left','String','A0/B0','Position',[10 15 32 13]);
        a.S0B1UI4 = uicontrol('Parent',a.p3UI4,'Style','checkbox','Value',a.S0_Polar,'Tooltip','Check this to calculate the fundamental symmetric (S1) or nonsymmetric (B2) mode.','Position',[160 70 50 23],'Callback',@CallbackUI4,'Tag','5');        
        a.SH0UI4 = uicontrol('Parent',a.p3UI4,'Style','checkbox','Value',a.SH0_Polar,'Tooltip','Check this to calculate the fundamental symmetric (S0) or nonsymmetric (B1) mode.','Position',[160 40 50 23],'Callback',@CallbackUI4,'Tag','6');
        a.A0B0UI4 = uicontrol('Parent',a.p3UI4,'Style','checkbox','Value',a.A0_Polar,'Tooltip','Check this to calculate the fundamental antisymmetric (A0) or nonsymmetric (B0) mode.','Position',[160 10 50 23],'Callback',@CallbackUI4,'Tag','7');

        %------------------------------------------------------------------
        a.p8UI4 = uipanel('Parent',Tab4,'Title','Trace modes and calculate energy velocity','Units','pixels','Position',[245 635 263 120],'FontSize',10);
        a.SamplesX3ceTextUI4 = uicontrol('Parent',a.p8UI4,'Style','text','HorizontalAlignment','left','String','Samples x3','Position',[30 40 58 13]);
        a.SamplesX3ceUI4 = uicontrol('Parent',a.p8UI4,'Style','edit','String',a.SamplesX3ce_Polar,'Tooltip','Enter the number of through-thickness samples for the energy velocity calculation. The more samples the more accurate the result. The accuracy decreases with an increasing frequency-thickness product because the mode shapes become more complex. Therefore, choose a sufficiently high number of samples.','Position',[34 15 50 23],'Callback',@CallbackUI4,'Tag','35');        
        a.TraceModesUI4 = uicontrol('Parent',a.p8UI4,'Style','togglebutton','String','Trace modes','Tooltip','Start tracing the dispersion curves.','Position',[117 55 130 40],'FontSize',10,'Callback',@CallbackUI4,'BusyAction','cancel','Tag','8');
        a.CalculateCeUI4 = uicontrol('Parent',a.p8UI4,'Style','togglebutton','String','Calculate ce','Tooltip','Start calculating the energy velocity.','Position',[117 10 130 40],'FontSize',10,'Callback',@CallbackUI4,'BusyAction','cancel','Tag','9');

        a.p4UI4 = uipanel('Parent',Tab4,'Title','Dispersion diagrams','Units','pixels','Position',[245 465 263 165],'FontSize',10);
        uicontrol('Parent',a.p4UI4,'Style','text','HorizontalAlignment','left','String','Quantity','Position',[10 125 42 13]);
        a.Option1TextUI4 = uicontrol('Parent',a.p4UI4,'Style','text','HorizontalAlignment','left','String','Bulk velocities','Position',[10 95 71 13]);
        uicontrol('Parent',a.p4UI4,'Style','text','HorizontalAlignment','left','String','Frequency (kHz)','Position',[10 65 83 13]);
        a.QuantityUI4 = uicontrol('Parent',a.p4UI4,'Style','popupmenu','String',{'Phase velocity (m/ms)','Energy velocity (m/ms)',['Propagation time (',char(181),'s)'],['Coincidence angle (',char(176),')'],'Wavelength (mm)','Wavenumber (rad/mm)'},'Tooltip','Select which quantity to plot in the polar dispersion diagram.','Position',[117 120 130 23],'Callback',@CallbackUI4,'Tag','10');
        a.Option1UI4 = uicontrol('Parent',a.p4UI4,'Style','checkbox','Value',a.BulkVelocities_Polar,'Tooltip','Check this to show the bulk wave velocities.','Position',[117 90 130 23],'Callback',@CallbackUI4,'Tag','11');
        a.FrequencyUI4 = uicontrol('Parent',a.p4UI4,'Style','popupmenu','String',{''},'Tooltip','Select the frequency for which the polar dispersion curves shall be plotted.','Position',[117 60 65 23],'Callback',@CallbackUI4,'Tag','12');
        a.PlotUI4 = uicontrol('Parent',a.p4UI4,'Style','pushbutton','String','Plot','Tooltip','Plot the polar dispersion diagram. If you have checked ''Export plots'' in the export settings, the plot will be exported automatically.','Position',[117 15 65 33],'FontSize',10,'Callback',@CallbackUI4,'Tag','13');

        a.p5UI4 = uipanel('Parent',Tab4,'Title','Export settings','Units','pixels','Position',[245 280 473 180],'FontSize',10);
        uicontrol('Parent',a.p5UI4,'Style','text','HorizontalAlignment','left','String','Export plots','Position',[10 140 59 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','HorizontalAlignment','left','String','PDF','Position',[124 140 26 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','HorizontalAlignment','left','String','PNG','Position',[124 120 28 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','HorizontalAlignment','left','String','PNG resolution (dpi)','Position',[10 93 98 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','HorizontalAlignment','left','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p5UI4,'Style','text','HorizontalAlignment','left','String','Directory','Position',[10 20 46 13]);
        a.ExportPlots_PolarUI4 = uicontrol('Parent',a.p5UI4,'Style','checkbox','Value',a.ExportPlots_Polar,'Tooltip','Check this in order to export plots automatically upon pressing the plot button. You can also export plots manually by using the ''File'' menu inside the plot figure.','Position',[80 135 20 23],'Callback',@CallbackUI4,'Tag','14');
        a.PDF_PolarUI4 = uicontrol('Parent',a.p5UI4,'Style','checkbox','Value',a.PDF_Polar,'Tooltip','Check this to export a plot as pdf.','Position',[160 135 20 23],'Callback',@CallbackUI4,'Tag','15');
        a.PNG_PolarUI4 = uicontrol('Parent',a.p5UI4,'Style','checkbox','Value',a.PNG_Polar,'Tooltip','Check this to export a plot as png.','Position',[160 115 20 23],'Callback',@CallbackUI4,'Tag','16');
        a.PNGresolution_PolarUI4 = uicontrol('Parent',a.p5UI4,'Style','edit','String',a.PNGresolution_Polar,'Tooltip','Enter the resolution of the png image.','Position',[125 88 50 23],'Callback',@CallbackUI4,'Tag','17');
        uicontrol('Parent',a.p5UI4,'Style','pushbutton','String','*.mat','Tooltip','Export the data as Matlab''s mat-file.','Position',[220 80 40 33],'Callback',@CallbackUI4,'Tag','18');
        uicontrol('Parent',a.p5UI4,'Style','pushbutton','String','*.xlsx','Tooltip','Export the data as Excel sheet.','Position',[270 80 40 33],'Callback',@CallbackUI4,'Tag','19');
        uicontrol('Parent',a.p5UI4,'Style','pushbutton','String','*.txt','Tooltip','Export the data as txt-file.','Position',[320 80 40 33],'Callback',@CallbackUI4,'Tag','20');
        a.FileName_PolarUI4 = uicontrol('Parent',a.p5UI4,'Style','edit','String',a.FileName_Polar,'Tooltip','Specify the name of the plot to be exported. The dispersion curve raw data are named automatically.','Position',[70 45 385 23],'Callback',@CallbackUI4,'Tag','21');
        a.DirectoryUI4 = uicontrol('Parent',a.p5UI4,'Style','edit','String',a.Directory,'Tooltip','Specify the directory to which the plot and the dispersion curve raw data shall be exported.','Position',[70 15 385 23],'Callback',@CallbackUI4,'Tag','22');
        
        %------------------------------------------------------------------
        a.p6UI4 = uipanel('Parent',Tab4,'Title','Plot layout settings','Units','pixels','Position',[518 520 200 235],'FontSize',10);
        uicontrol('Parent',a.p6UI4,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 175 21 13]);
        uicontrol('Parent',a.p6UI4,'Style','text','HorizontalAlignment','left','String','Curve line width','Position',[10 145 80 13]);
        a.TitleUI4 = uicontrol('Parent',a.p6UI4,'Style','popupmenu','String',{'with layup','without layup','no title'},'Tooltip','Choose to plot the title with the layup  included, without layup, or no title at all.','Position',[95 170 90 23],'Callback',@CallbackUI4,'Tag','23');
        a.LineWidthUI4 = uicontrol('Parent',a.p6UI4,'Style','edit','String',a.LineWidth_Polar,'Tooltip','Enter the curve line width.','Position',[95 140 50 23],'Callback',@CallbackUI4,'Tag','24');

        a.p7UI4 = uipanel('Parent',Tab4,'Title','Font size','Units','pixels','Position',[528 534 180 115],'FontSize',9);
        uicontrol('Parent',a.p7UI4,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 75 21 13]);
        uicontrol('Parent',a.p7UI4,'Style','text','HorizontalAlignment','left','String','Axes ticks','Position',[10 45 53 13]);
        uicontrol('Parent',a.p7UI4,'Style','text','HorizontalAlignment','left','String','Mode labels','Position',[10 15 59 13]);
        a.TitleFontSizeUI4 = uicontrol('Parent',a.p7UI4,'Style','edit','String',a.TitleFontSize_Polar,'Tooltip','Set the title font size.','Position',[85 70 50 23],'Callback',@CallbackUI4,'Tag','28');
        a.AxesTickFontSizeUI4 = uicontrol('Parent',a.p7UI4,'Style','edit','String',a.AxesTickFontSize_Polar,'Tooltip','Set the axes tick label font size.','Position',[85 40 50 23],'Callback',@CallbackUI4,'Tag','29');
        a.ModeLabelFontSizeUI4 = uicontrol('Parent',a.p7UI4,'Style','edit','String',a.ModeLabelFontSize_Polar,'Tooltip','Set the mode label font size.','Position',[85 10 50 23],'Callback',@CallbackUI4,'Tag','30');

        a.c3UI4 = uicontrol('Parent',Tab4,'Style','pushbutton','String','Default','Tooltip','Reset the plot layout settings to the default.','Position',[730 693 65 33],'FontSize',10,'Callback',@CallbackUI4,'Tag','31');
    
        uicontrol('Parent',Tab4,'Style','text','HorizontalAlignment','left','String','CPU cores','Position',[1035 715 53 13]);
        a.c4UI4 = uicontrol('Parent',Tab4,'Style','text','String','0','Position',[1095 715 30 13],'BackgroundColor','white');
        a.a1UI4 = axes('Parent',Tab4,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading20.png')    
    end
    if  1 % Tab8_bulk waves
        a.OutputWindowUI8 = uicontrol('Parent',Tab8,'Style','text','String','','Position',[10 10 220 737],'BackgroundColor','white');
        [a.BulkWaves,a.X,a.Y] = Computer_SnellsLaw_Isotropic(a.Couplant_Bulk,a.Solid_Bulk,a.Theta_Bulk2,a.OutputWindowUI8);        
        a.OutputWindowUI8.String = ''; clc
        
        a.p1UI8 = uipanel('Parent',Tab8,'Title','Elastic waves in bulk material','Units','pixels','Position',[237 545 639 210],'FontSize',10);
        uicontrol('Parent',a.p1UI8,'Style','text','HorizontalAlignment','left','String','Class','Position',[10 175 29 13]);
        uicontrol('Parent',a.p1UI8,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic','Isotropic'},'Value',a.MaterialType_Bulk,'Tooltip','Select a material symmetry class.','Position',[10 145 130 23],'Callback',@CallbackUI8,'Tag','43');
        uicontrol('Parent',a.p1UI8,'Style','text','HorizontalAlignment','left','String','Material','Position',[10 125 39 13]);
        a.MaterialUI8 = uicontrol('Parent',a.p1UI8,'Style','popupmenu','String',a.FieldNames_Orthotropic,'Tooltip','Select a material.','Position',[10 95 130 23],'Callback',@CallbackUI8,'Tag','1');
        uicontrol('Parent',a.p1UI8,'Style','text','HorizontalAlignment','left','String','Quantity','Position',[10 75 42 13]);
        uicontrol('Parent',a.p1UI8,'Style','popupmenu','String',{'Phase velocity (m/ms)','Group velocity (m/ms)','Slowness (ms/m)',['Polarization skew (',char(176),')'],['Energy skew (',char(176),')']},'Tooltip','Select which quantity to plot.','Position',[10 45 130 23],'Callback',@CallbackUI8,'Tag','2');
        uicontrol('Parent',a.p1UI8,'Style','text','HorizontalAlignment','left','String',[char(916),char(920),' (',char(176),')'],'Position',[10 15 33 13]);
        uicontrol('Parent',a.p1UI8,'Style','edit','String',a.ThetaStep_Bulk1,'Tooltip','Enter the Theta step.','Position',[80 10 50 23],'Callback',@CallbackUI8,'Tag','3');
        
        a.p2UI8 = uipanel('Parent',Tab8,'Title','2-D profiles','Units','pixels','Position',[394 550 85 185],'FontSize',9);
        uicontrol('Parent',a.p2UI8,'Style','text','HorizontalAlignment','left','String','Plane','Position',[10 145 28 13]);
        uicontrol('Parent',a.p2UI8,'Style','text','HorizontalAlignment','left','String',[char(934),' (',char(176),')'],'Position',[10 90 25 13]);
        uicontrol('Parent',a.p2UI8,'Style','popupmenu','String',{'1-3','1-2'},'Tooltip','Select to plot the 1-3-plane or the 1-2-plane. In the upper sketch to the right, the x''3 and x''2-axes are swapped according to your selection.','Position',[10 115 50 23],'Callback',@CallbackUI8,'Tag','44');                        
        uicontrol('Parent',a.p2UI8,'Style','edit','String',a.Phi_Bulk11,'Tooltip','Enter Phi.','Position',[10 60 50 23],'Callback',@CallbackUI8,'Tag','4');
        a.Plot1UI8 = uicontrol('Parent',a.p2UI8,'Style','pushbutton','String','Plot','Tooltip','Plot the profiles. If you have checked ''Export plots'' in the export settings, the plot will be exported automatically.','Position',[10 15 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','5');

        a.p3UI8 = uipanel('Parent',Tab8,'Title','3-D surfaces and vectors','Units','pixels','Position',[489 550 377 185],'FontSize',9);
        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String',[char(916),char(934),' (',char(176),')'],'Position',[10 145 33 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String','Mode','Position',[10 65 28 13]);
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.PhiStep_Bulk,'Tooltip','Enter the Phi step.','Position',[50 140 50 23],'Callback',@CallbackUI8,'Tag','6');
        a.CalculateUI8 = uicontrol('Parent',a.p3UI8,'Style','pushbutton','String','Calculate','Position',[50 95 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','7');
        uicontrol('Parent',a.p3UI8,'Style','popupmenu','String',{'Longitudinal','Fast shear','Slow shear'},'Tooltip','Select which mode to plot.','Position',[50 60 80 23],'Callback',@CallbackUI8,'Tag','8');                        
        a.Plot2UI8 = uicontrol('Parent',a.p3UI8,'Style','pushbutton','String','Plot','Tooltip','Plot the surface. If you have checked ''Export plots'' in the export settings, the plot will be exported automatically.','Position',[50 15 65 33],'Enable','off','FontSize',10,'Callback',@CallbackUI8,'Tag','9');

        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String',[char(934),' (',char(176),')'],'Position',[136 145 25 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String',[char(920),' (',char(176),')'],'Position',[136 115 25 13]);
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.Phi_Bulk21,'Tooltip','Enter propagation direction Phi.','Position',[166 140 50 23],'Callback',@CallbackUI8,'Tag','10');
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.Theta_Bulk1,'Tooltip','Enter propagation direction Theta.','Position',[166 110 50 23],'Callback',@CallbackUI8,'Tag','11');
        a.Plot3UI8 = uicontrol('Parent',a.p3UI8,'Style','pushbutton','String','Plot','Tooltip','Plot the vectors. If you have checked ''Export plots'' in the export settings, the plot will be exported automatically.','Position',[166 15 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','12');
                        
        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String',['View ',char(934),' (',char(176),')'],'Position',[246 145 54 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String',['View ',char(920),' (',char(176),')'],'Position',[246 115 54 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String','Marker size','Position',[246 85 58 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String','Colorbar','Position',[246 55 43 13]);
        uicontrol('Parent',a.p3UI8,'Style','text','HorizontalAlignment','left','String','x-pos. (0-1)','Position',[246 40 61 13]);
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.ViewPhi_Bulk1,'Tooltip','Enter the view angle Phi.','Position',[310 140 50 23],'Callback',@CallbackUI8,'Tag','13');
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.ViewTheta_Bulk1,'Tooltip','Enter the view angle Theta.','Position',[310 110 50 23],'Callback',@CallbackUI8,'Tag','14');
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.MarkerSize_Bulk,'Tooltip','Enter the size of the markers.','Position',[310 80 50 23],'Callback',@CallbackUI8,'Tag','15');
        uicontrol('Parent',a.p3UI8,'Style','edit','String',a.ColorbarX_Bulk,'Tooltip','Shift the colorbar left-right.','Position',[310 50 50 23],'Callback',@CallbackUI8,'Tag','16');

        %------------------------------------------------------------------
        a.p4UI8 = uipanel('Parent',Tab8,'Title','Bulk waves on interfaces','Units','pixels','Position',[237 275 245 265],'FontSize',10);
        uicontrol('Parent',a.p4UI8,'Style','text','HorizontalAlignment','left','String','Fluid','Position',[10 225 24 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','HorizontalAlignment','left','String','Class','Position',[10 195 29 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','HorizontalAlignment','left','String','Solid','Position',[10 165 25 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','HorizontalAlignment','left','String',[char(934),' (',char(176),')'],'Position',[10 135 25 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','HorizontalAlignment','left','String',[char(920),'i (',char(176),')'],'Position',[10 105 27 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','HorizontalAlignment','left','String',[char(916),char(920),' (',char(176),')'],'Position',[10 75 33 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','HorizontalAlignment','left','String',['View ',char(934),' (',char(176),')'],'Position',[10 45 54 13]);
        uicontrol('Parent',a.p4UI8,'Style','text','HorizontalAlignment','left','String',['View ',char(920),' (',char(176),')'],'Position',[10 15 54 13]);
        
        a.CouplantUI8 = uicontrol('Parent',a.p4UI8,'Style','popupmenu','String',a.FieldNames_Fluid,'Tooltip','Select a default fluid from which a  plane wave impinges on the solid.','Position',[80 220 150 23],'Callback',@CallbackUI8,'Tag','17');
        uicontrol('Parent',a.p4UI8,'Style','popupmenu','String',{'Isotropic','Cubic','Transversely isotropic','Orthotropic'},'Tooltip','Select the material type of the solid.','Position',[80 190 150 23],'Callback',@CallbackUI8,'Tag','19');
        a.SolidUI8 = uicontrol('Parent',a.p4UI8,'Style','popupmenu','String',a.FieldNames_Isotropic,'Tooltip','Select a material.','Position',[80 160 150 23],'Callback',@CallbackUI8,'Tag','20');
        a.Phi2UI8 = uicontrol('Parent',a.p4UI8,'Style','edit','String',a.Phi_Bulk2,'Tooltip','Enter Phi.','Position',[80 130 50 23],'Enable','off','Callback',@CallbackUI8,'Tag','21');
        uicontrol('Parent',a.p4UI8,'Style','edit','String',a.Theta_Bulk2,'Tooltip','Enter the incidence angle Theta.','Position',[80 100 50 23],'Callback',@CallbackUI8,'Tag','22');
        a.ThetaStep2UI8 = uicontrol('Parent',a.p4UI8,'Style','edit','String',a.ThetaStep_Bulk2,'Tooltip','Enter the Theta step for plotting the reflection and transmission coefficients.','Position',[80 70 50 23],'Callback',@CallbackUI8,'Tag','23');
        a.ViewPhi2UI8 = uicontrol('Parent',a.p4UI8,'Style','edit','String',a.ViewPhi_Bulk2,'Tooltip','Enter the view angle Phi.','Position',[80 40 50 23],'Enable','off','Callback',@CallbackUI8,'Tag','24');
        a.ViewTheta2UI8 = uicontrol('Parent',a.p4UI8,'Style','edit','String',a.ViewTheta_Bulk2,'Tooltip','Enter the view angle Theta.','Position',[80 10 50 23],'Enable','off','Callback',@CallbackUI8,'Tag','25');
        a.Plot2DUI8 = uicontrol('Parent',a.p4UI8,'Style','pushbutton','String','Plot 2-D','Tooltip','Plot the bulk waves. If you have checked ''Export plots'' in   the export settings, the plot will be exported automatically.','Position',[165 90 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','26');        
        a.Plot3DUI8 = uicontrol('Parent',a.p4UI8,'Style','pushbutton','String','Plot 3-D','Tooltip','Plot the bulk waves. If you have checked ''Export plots'' in  the export settings, the plot will be exported automatically.','Position',[165 50 65 33],'FontSize',10,'Enable','off','Callback',@CallbackUI8,'Tag','27');
        a.PlotRTUI8 = uicontrol('Parent',a.p4UI8,'Style','pushbutton','String','Plot R,T','Tooltip','Plot reflected and transmitted energy coefficients versus incidence angle.','Position',[165 10 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','46');

        %------------------------------------------------------------------
        a.p5UI8 = uipanel('Parent',Tab8,'Title','Plot layout settings','Units','pixels','Position',[492 215 170 325],'FontSize',10);
        uicontrol('Parent',a.p5UI8,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 285 21 13]);
        a.TitleUI8 = uicontrol('Parent',a.p5UI8,'Style','checkbox','Value',a.Title_Bulk,'Tooltip','Check this in order to show the plot title.','Position',[95 280 20 23],'Callback',@CallbackUI8,'Tag','28');      

        a.p6UI8 = uipanel('Parent',Tab8,'Title','Line width','Units','pixels','Position',[502 375 150 115],'FontSize',9);
        uicontrol('Parent',a.p6UI8,'Style','text','HorizontalAlignment','left','String','Box','Position',[10 75 21 13]);
        uicontrol('Parent',a.p6UI8,'Style','text','HorizontalAlignment','left','String','Profile','Position',[10 45 32 13]);
        uicontrol('Parent',a.p6UI8,'Style','text','HorizontalAlignment','left','String','Bulk wave','Position',[10 15 53 13]);
        a.BoxLineWidthUI8 = uicontrol('Parent',a.p6UI8,'Style','edit','String',a.BoxLineWidth_Bulk,'Tooltip','Enter the box line width.','Position',[85 70 50 23],'Callback',@CallbackUI8,'Tag','29');
        a.LineWidthUI8 = uicontrol('Parent',a.p6UI8,'Style','edit','String',a.LineWidth_Bulk,'Tooltip','Enter the profile line width.','Position',[85 40 50 23],'Callback',@CallbackUI8,'Tag','30');
        a.WaveVectorLineWidthUI8 = uicontrol('Parent',a.p6UI8,'Style','edit','String',a.WaveVectorLineWidth_Bulk,'Tooltip','Enter the bulk wave line width.','Position',[85 10 50 23],'Callback',@CallbackUI8,'Tag','31');
        
        a.p7UI8 = uipanel('Parent',Tab8,'Title','Font size','Units','pixels','Position',[502 225 150 145],'FontSize',9);
        uicontrol('Parent',a.p7UI8,'Style','text','HorizontalAlignment','left','String','Title','Position',[10 105 21 13]);
        uicontrol('Parent',a.p7UI8,'Style','text','HorizontalAlignment','left','String','Axes labels','Position',[10 75 59 13]);
        uicontrol('Parent',a.p7UI8,'Style','text','HorizontalAlignment','left','String','Axes ticks','Position',[10 45 53 13]);
        uicontrol('Parent',a.p7UI8,'Style','text','HorizontalAlignment','left','String','Mode labels','Position',[10 15 59 13]);
        a.TitleFontSizeUI8 = uicontrol('Parent',a.p7UI8,'Style','edit','String',a.TitleFontSize_Bulk,'Tooltip','Set the title font size.','Position',[85 100 50 23],'Callback',@CallbackUI8,'Tag','32');
        a.AxesLabelFontSizeUI8 = uicontrol('Parent',a.p7UI8,'Style','edit','String',a.AxesLabelFontSize_Bulk,'Tooltip','Set the axes label font size.','Position',[85 70 50 23],'Callback',@CallbackUI8,'Tag','33');
        a.AxesTickFontSizeUI8 = uicontrol('Parent',a.p7UI8,'Style','edit','String',a.AxesTickFontSize_Bulk,'Tooltip','Set the axes tick label font size.','Position',[85 40 50 23],'Callback',@CallbackUI8,'Tag','34');
        a.ModeLabelFontSizeUI8 = uicontrol('Parent',a.p7UI8,'Style','edit','String',a.ModeLabelFontSize_Bulk,'Tooltip','Set the mode label font size.','Position',[85 10 50 23],'Callback',@CallbackUI8,'Tag','35');        

        %------------------------------------------------------------------
        a.p8UI8 = uipanel('Parent',Tab8,'Title','Export settings','Units','pixels','Position',[237 70 425 140],'FontSize',10);
        uicontrol('Parent',a.p8UI8,'Style','text','HorizontalAlignment','left','String','Export plots','Position',[10 100 59 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','HorizontalAlignment','left','String','PDF','Position',[124 100 26 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','HorizontalAlignment','left','String','PNG','Position',[124 80 28 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','HorizontalAlignment','left','String','PNG resolution (dpi)','Position',[235 80 98 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','HorizontalAlignment','left','String','File name','Position',[10 50 47 13]);
        uicontrol('Parent',a.p8UI8,'Style','text','HorizontalAlignment','left','String','Directory','Position',[10 20 46 13]);
        uicontrol('Parent',a.p8UI8,'Style','checkbox','Value',a.ExportPlots_Bulk,'Tooltip','Check this in order to export plots automatically upon pressing the respective plot button. You can also export plots manually by using the ''File'' menu inside the plot figure.','Position',[80 95 20 23],'Callback',@CallbackUI8,'Tag','36');
        uicontrol('Parent',a.p8UI8,'Style','checkbox','Value',a.PDF_Bulk,'Tooltip','Check this to export a plot as pdf.','Position',[160 95 20 23],'Callback',@CallbackUI8,'Tag','37');
        uicontrol('Parent',a.p8UI8,'Style','checkbox','Value',a.PNG_Bulk,'Tooltip','Check this to export a plot as png.','Position',[160 75 20 23],'Callback',@CallbackUI8,'Tag','38');
        uicontrol('Parent',a.p8UI8,'Style','edit','String',a.PNGresolution_Bulk,'Tooltip','Enter the resolution of the png image.','Position',[350 75 50 23],'Callback',@CallbackUI8,'Tag','39');
        uicontrol('Parent',a.p8UI8,'Style','edit','String',a.FileName_Bulk,'Tooltip','Specify the name of the plots to be exported.','Position',[70 45 330 23],'Callback',@CallbackUI8,'Tag','40');
        uicontrol('Parent',a.p8UI8,'Style','edit','String',a.Directory,'Tooltip','Specify the directory to which plots shall be exported.','Position',[70 15 330 23],'Callback',@CallbackUI8,'Tag','41');        

        a.c1UI8 = uicontrol('Parent',Tab8,'Style','pushbutton','String','Default','Tooltip','Reset the plot layout settings to the default.','Position',[675 500 65 33],'FontSize',10,'Callback',@CallbackUI8,'Tag','42');    

        a.a1UI8 = axes('Parent',Tab8,'Units','pixels','Position',[880 390 312 367]);
        imshow('Bulk1.png')
        a.a2UI8 = axes('Parent',Tab8,'Units','pixels','Position',[806 5 386 381]);
        imshow('Bulk2.png')        
    end
    if  1 % Tab6_laminate stiffness
        a.p1UI6 = uipanel('Parent',Tab6,'Title','','Units','pixels','Position',[0 0 1198 1000],'FontSize',10);
        
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','Specimen','Position',[20 565+170 58 15],'FontSize',9);
        uicontrol('Parent',a.p1UI6,'Style','pushbutton','String','Edit','Tooltip','Press to define your specimen.','Position',[81.5 565+125 127 40],'FontSize',10,'Callback',@SpecimenSettingsUI6_Callback);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','Material','Position',[20 565+95 39 13]);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','Unit cell','Position',[20 565+65 39 13]);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String',['Azimuthal angle (',char(176),')'],'Position',[20 565+35 93 13]);
        a.MaterialNameUI6 = uicontrol('Parent',a.p1UI6,'Style','text','String',a.Material6{1}.Name,'Position',[81.5 565+90 190 23],'BackgroundColor','white');
        a.LayupUI6 = uicontrol('Parent',a.p1UI6,'Style','text','String','[0]','Position',[81.5 565+60 190 23],'BackgroundColor','white');
        uicontrol('Parent',a.p1UI6,'Style','edit','String',a.PropagationAngle,'Tooltip','Enter the azimuthal angle with respect to the fiber orientations.','Position',[140 565+30 50 23],'Callback',@CallbackUI6,'Tag','2');

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

        uicontrol('Parent',a.p1UI6,'Style','pushbutton','String','*.mat','Tooltip','Export the data as Matlab''s mat-file.','Position',[750 590 40 33],'Callback',@CallbackUI6,'Tag','3');
        uicontrol('Parent',a.p1UI6,'Style','pushbutton','String','*.xlsx','Tooltip','Export the data as Excel sheet.','Position',[800 590 40 33],'Callback',@CallbackUI6,'Tag','4');
        uicontrol('Parent',a.p1UI6,'Style','pushbutton','String','*.txt','Tooltip','Export the data as txt-file.','Position',[850 590 40 33],'Callback',@CallbackUI6,'Tag','5');
        uicontrol('Parent',a.p1UI6,'Style','edit','String',a.Directory,'Tooltip','Specify the directory to which the layup stiffness matrix shall be exported.','Position',[750 495+60 410 23],'Callback',@CallbackUI6,'Tag','6');        

        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C11','Position',[20,130+365,21,13],'foregroundcolor','r');
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C12','Position',[20,130+335,21,13],'foregroundcolor',[.13 .55 .13]);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C13','Position',[20,130+305,21,13],'foregroundcolor','b');
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C16','Position',[20,130+275,21,13],'foregroundcolor','k');
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C22','Position',[20,130+245,21,13],'foregroundcolor','m');
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C23','Position',[20,130+215,21,13],'foregroundcolor','c');
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C26','Position',[20,130+185,21,13],'foregroundcolor',[1 .7 0]);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C33','Position',[20,130+155,21,13],'foregroundcolor',[.55 .27 .13]);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C36','Position',[20,130+125,21,13],'foregroundcolor',[.5 0 1]);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C44','Position',[20,130+95,21,13],'foregroundcolor',[.5 .5 .5]);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C45','Position',[20,130+65,21,13],'foregroundcolor','r');
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C55','Position',[20,130+35,21,13],'foregroundcolor',[.13 .55 .13]);
        uicontrol('Parent',a.p1UI6,'Style','text','HorizontalAlignment','left','String','C66','Position',[20,130+5,21,13],'foregroundcolor','b');
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
        
        a.CLaminate = Computer_LaminateStiffness(1,a.Material6,a.LayerOrientations6,a.LayerThicknesses6,a.PropagationAngle6);
        a.CLP = Computer_LaminateStiffness(2,a.Material6,a.LayerOrientations6,a.LayerThicknesses6,1);
        
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
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Material','Position',[10 180+310 39 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Mass density (kg/m3)','Position',[10 180+280 105 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Engineering constants','Position',[10 160+255 126 15],'FontSize',9);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Real part','Position',[60 160+235 45 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Imaginary part','Position',[135 160+235 70 15]);        
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','E (GPa)','Position',[10 180+190 39 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','v','Position',[10 180+160 6 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Stiffness components (GPa)','Position',[250 160+255 157 15],'FontSize',9);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Real part','Position',[280 160+235 45 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Imaginary part','Position',[355 160+235 70 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','C11','Position',[250 180+190 21 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','C66','Position',[250 180+160 21 15]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Bulk waves','Position',[10 180+115 64 15],'FontSize',9);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Longitudinal velocity (m/s)','Position',[10 180+90 127 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Shear velocity (m/s)','Position',[10 180+60 99 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','Attenuation unit','Position',[10 180+30 77 13]);
        a.LongitudinalAttenuationTextUI5 = uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String',['Longitudinal attenuation (Np/',char(955),')'],'Position',[10 180+0 150 13]);
        a.TransverseAttenuationTextUI5 = uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String',['Shear attenuation (Np/',char(955),')'],'Position',[10 180-30 122 13]);
        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','At frequency (kHz)','Position',[10 180-60 95 13]);
        a.Material1UI5 = uicontrol('Parent',a.p1UI5,'Style','popupmenu','String',a.FieldNames_Isotropic,'Tooltip','Select a material to display its parameters.','Position',[135 180+305 220 23],'Callback',@CallbackUI5,'Tag','1');
        a.Density1UI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.Density,'Tooltip','Enter the mass density.','Position',[135 180+275 65 23],'Callback',@CallbackUI5,'Tag','2');
        a.YoungsModulusUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',real(a.Material1.YoungsModulus)/1e9,'Tooltip','Enter the real part of Young''s modulus.','Position',[60 160+205 65 23],'Callback',@CallbackUI5,'Tag','3');
        a.YoungsModulusImagUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',imag(a.Material1.YoungsModulus)/1e9,'Tooltip','Enter the imaginary part of Young''s modulus.','Position',[135 160+205 65 23],'Callback',@CallbackUI5,'Tag','4');        
        a.PoissonsNumberUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',real(a.Material1.PoissonsNumber),'Tooltip','Enter the real part of Poisson''s ratio.','Position',[60 160+175 65 23],'Callback',@CallbackUI5,'Tag','5');
        a.PoissonsNumberImagUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',imag(a.Material1.PoissonsNumber),'Tooltip','Enter the imaginary part of Poisson''s ratio.','Position',[135 160+175 65 23],'Callback',@CallbackUI5,'Tag','6');
        a.C11IsoUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',real(a.Material1.C(1,1))/1e9,'Tooltip','Enter the real stiffness component.','Position',[280 160+205 65 23],'Callback',@CallbackUI5,'Tag','64');
        a.C11ImagIsoUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',imag(a.Material1.C(1,1))/1e9,'Tooltip','Enter the imaginary stiffness component.','Position',[355 160+205 65 23],'Callback',@CallbackUI5,'Tag','65');
        a.C66IsoUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',real(a.Material1.C(6,6))/1e9,'Tooltip','Enter the real stiffness component.','Position',[280 160+175 65 23],'Callback',@CallbackUI5,'Tag','66');        
        a.C66ImagIsoUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',imag(a.Material1.C(6,6))/1e9,'Tooltip','Enter the imaginary stiffness component.','Position',[355 160+175 65 23],'Callback',@CallbackUI5,'Tag','67');
        a.LongitudinalVelocityUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.LongitudinalVelocity,'Tooltip','Enter the longitudinal bulk velocity.','Position',[165 180+85 65 23],'Callback',@CallbackUI5,'Tag','7');
        a.TransverseVelocityUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.TransverseVelocity,'Tooltip','Enter the shear bulk velocity.','Position',[165 180+55 65 23],'Callback',@CallbackUI5,'Tag','8');
        uicontrol('Parent',a.p1UI5,'Style','popupmenu','String',{['Np/',char(955)],['dB/',char(955)],'Np/m','dB/m'},'Tooltip','Select in which units you want to enter the attenuation of the bulk waves.  If you select ''Np/m'' or ''dB/m'', you need to enter at which frequency the attenuation was measured. DC assumes a linear increase of attenuation with frequency, i.e., a damping loss which is constant per wavelength (hysteretic damping).','Position',[165 180+25 65 23],'Callback',@CallbackUI5,'Tag','9');
        a.LongitudinalAttenuationUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.LongitudinalAttenuation,'Tooltip','The attenuation of longitudinal bulk waves.','Position',[165 180-5 65 23],'Callback',@CallbackUI5,'Tag','10');
        a.TransverseAttenuationUI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.Material1.LongitudinalAttenuation,'Tooltip','The attenuation of shear bulk waves.','Position',[165 180-35 65 23],'Callback',@CallbackUI5,'Tag','11');
        a.AtFrequency1UI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.AtFrequency1,'Tooltip','Enter at which frequency the above attenuations were measured.','Position',[165 180-65 65 23],'Enable','off','Callback',@CallbackUI5,'Tag','12');

        uicontrol('Parent',a.p1UI5,'Style','text','HorizontalAlignment','left','String','New material''s name','Position',[10 110-45 102 13]);
        a.Name1UI5 = uicontrol('Parent',a.p1UI5,'Style','edit','String',a.MaterialName4A,'Tooltip','Enter a new material''s name and press ''Save'' or enter an available material''s name and press ''Delete''.','Position',[135 110-50 220 23],'Callback',@CallbackUI5,'Tag','13');
        uicontrol('Parent',a.p1UI5,'Style','pushbutton','String','Save material','Tooltip','Save the material to the isotropic materials list.','Position',[135 110-95 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','14');
        uicontrol('Parent',a.p1UI5,'Style','pushbutton','String','Delete material','Tooltip','Remove the material from the isotropic materials list.','Position',[255 110-95 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','15');

        %------------------------------------------------------------------        
        a.p2UI5 = uipanel('Parent',Tab5,'Title','Anisotropic materials','Units','pixels','Position',[460 10 725 745],'FontSize',10);
        uicontrol('Parent',a.p2UI5,'Style','pushbutton','String','?','Position',[680 695 23 23],'FontSize',12,'FontWeight','bold','foregroundcolor',[1 1 1],'backgroundcolor',[.4 .4 .4],'Callback',@CallbackUI5,'Tag','69');
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Class','Position',[10 175+525 29 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Material','Position',[10 175+495 39 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Mass density (kg/m3)','Position',[10 175+465 105 13]);
        
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Engineering constants (GPa)','Position',[10 155+440 126 15],'FontSize',9);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Real part','Position',[75 155+420 45 15]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Imaginary part','Position',[150 155+420 70 15]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','E1 (GPa)','Position',[10 155+395 49 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','E2 (GPa)','Position',[10 155+365 49 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','E3 (GPa)','Position',[10 155+335 49 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','G12 (GPa)','Position',[10 155+305 57 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','G13 (GPa)','Position',[10 155+275 57 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','G23 (GPa)','Position',[10 155+245 57 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','v12','Position',[10 155+215 20 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','v13','Position',[10 155+185 20 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','v23','Position',[10 155+155 20 13]);
        uicontrol('Parent',a.p2UI5,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic'},'Value',a.MaterialTypeME,'Tooltip','Select a material symmetry class.','Position',[135 155+540 220 23],'Callback',@CallbackUI5,'Tag','16');
        a.Material2UI5 = uicontrol('Parent',a.p2UI5,'Style','popupmenu','String',a.FieldNames_Orthotropic,'Tooltip','Select a material to display its parameters.','Position',[135 155+510 220 23],'Callback',@CallbackUI5,'Tag','17');
        a.Density2UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',a.Material2{1}.Density,'Tooltip','Enter the mass density.','Position',[135 155+480 65 23],'Callback',@CallbackUI5,'Tag','18');        
        a.E1UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.E1)/1e9,'Tooltip','Enter the real part of Young''s modulus in the 1-direction.','Position',[75 155+390 65 23],'Callback',@CallbackUI5,'Tag','19');
        a.E2UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.E2)/1e9,'Tooltip','Enter the real part of Young''s modulus in the 2-direction.','Position',[75 155+360 65 23],'Callback',@CallbackUI5,'Tag','20');
        a.E3UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.E3)/1e9,'Tooltip','Enter the real part of Young''s modulus in the 3-direction.','Position',[75 155+330 65 23],'Callback',@CallbackUI5,'Tag','21');
        a.G12UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.G12)/1e9,'Tooltip','Enter the real part of the shear modulus in the 2-direction on the 2-3-plane.','Position',[75 155+300 65 23],'Callback',@CallbackUI5,'Tag','22');
        a.G13UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.G13)/1e9,'Tooltip','Enter the real part of the shear modulus in the 3-direction on the 2-3-plane.','Position',[75 155+270 65 23],'Callback',@CallbackUI5,'Tag','23');
        a.G23UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.G23)/1e9,'Tooltip','Enter the real part of the shear modulus in the 3-direction on the 1-3-plane.','Position',[75 155+240 65 23],'Callback',@CallbackUI5,'Tag','24');
        a.v12UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.v12),'Tooltip','Enter the real part of Poisson''s ratio corresponding to a contraction in the 2-direction when an extension is applied in the 1-direction.','Position',[75 155+210 65 23],'Callback',@CallbackUI5,'Tag','25');
        a.v13UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.v13),'Tooltip','Enter the real part of Poisson''s ratio corresponding to a contraction in the 3-direction when an extension is applied in the 1-direction.','Position',[75 155+180 65 23],'Callback',@CallbackUI5,'Tag','26');
        a.v23UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.v23),'Tooltip','Enter the real part of Poisson''s ratio corresponding to a contraction in the 3-direction when an extension is applied in the 2-direction.','Position',[75 155+150 65 23],'Callback',@CallbackUI5,'Tag','27');
        a.E1ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.E1)/1e9,'Tooltip','Enter the imaginary part of Young''s modulus in the 1-direction.','Position',[150 155+390 65 23],'Callback',@CallbackUI5,'Tag','28');
        a.E2ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.E2)/1e9,'Tooltip','Enter the imaginary part of Young''s modulus in the 2-direction.','Position',[150 155+360 65 23],'Callback',@CallbackUI5,'Tag','29');
        a.E3ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.E3)/1e9,'Tooltip','Enter the imaginary part of Young''s modulus in the 3-direction.','Position',[150 155+330 65 23],'Callback',@CallbackUI5,'Tag','30');
        a.G12ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.G12)/1e9,'Tooltip','Enter the imaginary part of the shear modulus in the 2-direction on the 2-3-plane.','Position',[150 155+300 65 23],'Callback',@CallbackUI5,'Tag','31');
        a.G13ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.G13)/1e9,'Tooltip','Enter the imaginary part of the shear modulus in the 3-direction on the 2-3-plane.','Position',[150 155+270 65 23],'Callback',@CallbackUI5,'Tag','32');
        a.G23ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.G23)/1e9,'Tooltip','Enter the imaginary part of the shear modulus in the 3-direction on the 1-3-plane.','Position',[150 155+240 65 23],'Callback',@CallbackUI5,'Tag','33');
        a.v12ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.v12),'Tooltip','Enter the imaginary part of Poisson''s ratio corresponding to a contraction in the 2-direction when an extension is applied in the 1-direction.','Position',[150 155+210 65 23],'Callback',@CallbackUI5,'Tag','34');
        a.v13ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.v13),'Tooltip','Enter the imaginary part of Poisson''s ratio corresponding to a contraction in the 3-direction when an extension is applied in the 1-direction.','Position',[150 155+180 65 23],'Callback',@CallbackUI5,'Tag','35');
        a.v23ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.v23),'Tooltip','Enter the imaginary part of Poisson''s ratio corresponding to a contraction in the 3-direction when an extension is applied in the 2-direction.','Position',[150 155+150 65 23],'Callback',@CallbackUI5,'Tag','36');

        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Stiffness components (GPa)','Position',[230+40 155+440 157 15],'FontSize',9);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Real part','Position',[230+40 155+420 45 15]);
        a.C11UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(1,1))/1e9,'Tooltip','Enter the real stiffness component in the 1-direction.','Position',[230+40 155+390 65 23],'Callback',@CallbackUI5,'Tag','37');
        a.C12UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(1,2))/1e9,'Tooltip','Enter the real stiffness component.','Position',[230+115 155+390 65 23],'Callback',@CallbackUI5,'Tag','38');
        a.C13UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(1,3))/1e9,'Position',[230+190 155+390 65 23],'Callback',@CallbackUI5,'Tag','39');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 155+390 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 155+390 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+390 65 23],'BackgroundColor',[.88 .88 .88]);
        a.C22UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(2,2))/1e9,'Tooltip','Enter the real stiffness component in the 2-direction.','Position',[230+115 155+360 65 23],'Callback',@CallbackUI5,'Tag','40');
        a.C23UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(2,3))/1e9,'Tooltip','Enter the real stiffness component.','Position',[230+190 155+360 65 23],'Callback',@CallbackUI5,'Tag','41');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 155+360 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 155+360 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+360 65 23],'BackgroundColor',[.88 .88 .88]);
        a.C33UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(3,3))/1e9,'Tooltip','Enter the real stiffness component in the 3-direction.','Position',[230+190 155+330 65 23],'Callback',@CallbackUI5,'Tag','42');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 155+330 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 155+330 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+330 65 23],'BackgroundColor',[.88 .88 .88]);        
        a.C44UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(4,4))/1e9,'Tooltip','Enter the real stiffness component.','Position',[230+265 155+300 65 23],'Callback',@CallbackUI5,'Tag','43');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 155+300 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+300 65 23],'BackgroundColor',[.88 .88 .88]);        
        a.C55UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(5,5))/1e9,'Tooltip','Enter the real stiffness component.','Position',[230+340 155+270 65 23],'Callback',@CallbackUI5,'Tag','44');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 155+270 65 23],'BackgroundColor',[.88 .88 .88]);         
        a.C66UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',real(a.Material2{1}.C(6,6))/1e9,'Tooltip','Enter the real stiffness component.','Position',[230+415 155+240 65 23],'Callback',@CallbackUI5,'Tag','45');

        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Imaginary part','Position',[230+40 200+180 70 15]);
        a.C11ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(1,1))/1e9,'Tooltip','Enter the imaginary stiffness component in the 1-direction.','Position',[230+40 200+150 65 23],'Callback',@CallbackUI5,'Tag','46');
        a.C12ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(1,2))/1e9,'Tooltip','Enter the imaginary stiffness component.','Position',[230+115 200+150 65 23],'Callback',@CallbackUI5,'Tag','47');
        a.C13ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(1,3))/1e9,'Position',[230+190 200+150 65 23],'Callback',@CallbackUI5,'Tag','48');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 200+150 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 200+150 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+150 65 23],'BackgroundColor',[.88 .88 .88]);
        a.C22ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(2,2))/1e9,'Tooltip','Enter the imaginary stiffness component in the 2-direction.','Position',[230+115 200+120 65 23],'Callback',@CallbackUI5,'Tag','49');
        a.C23ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(2,3))/1e9,'Tooltip','Enter the imaginary stiffness component.','Position',[230+190 200+120 65 23],'Callback',@CallbackUI5,'Tag','50');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 200+120 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 200+120 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+120 65 23],'BackgroundColor',[.88 .88 .88]);
        a.C33ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(3,3))/1e9,'Tooltip','Enter the imaginary stiffness component in the 3-direction.','Position',[230+190 200+90 65 23],'Callback',@CallbackUI5,'Tag','51');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+265 200+90 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 200+90 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+90 65 23],'BackgroundColor',[.88 .88 .88]);        
        a.C44ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(4,4))/1e9,'Tooltip','Enter the imaginary stiffness component.','Position',[230+265 200+60 65 23],'Callback',@CallbackUI5,'Tag','52');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+340 200+60 65 23],'BackgroundColor',[.88 .88 .88]);
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+60 65 23],'BackgroundColor',[.88 .88 .88]);        
        a.C55ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(5,5))/1e9,'Tooltip','Enter the imaginary stiffness component.','Position',[230+340 200+30 65 23],'Callback',@CallbackUI5,'Tag','53');
        uicontrol('Parent',a.p2UI5,'Style','text','String','0','Position',[230+415 200+30 65 23],'BackgroundColor',[.88 .88 .88]);         
        a.C66ImagUI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',imag(a.Material2{1}.C(6,6))/1e9,'Tooltip','Enter the imaginary stiffness component.','Position',[230+415 200 65 23],'Callback',@CallbackUI5,'Tag','54');

        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Bulk wave velocities (m/s)','Position',[10 120+95 143 15],'FontSize',9);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','1','Position',[154 120+90 4 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','2','Position',[228 120+90 6 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','3','Position',[302 120+90 6 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Longitudinal velocity','Position',[10 120+65 99 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Fast shear velocity','Position',[10 120+35 94 13]);
        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','Slow shear velocity','Position',[10 120+5 98 13]);
        a.LongitudinalVelocity_1UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.LongitudinalVelocity_1,'Tooltip','The phase velocity of longitudinal bulk waves propagating in the 1-direction.','Position',[125 120+60 60 23],'BackgroundColor','white');
        a.FastShearVelocity_1UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.FastShearVelocity_1,'Tooltip','The phase velocity of fast shear bulk waves propagating in the 1-direction.','Position',[125 120+30 60 23],'BackgroundColor','white');
        a.SlowShearVelocity_1UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.SlowShearVelocity_1,'Tooltip','The phase velocity of slow shear bulk waves propagating in the 1-direction.','Position',[125 120+0 60 23],'BackgroundColor','white');
        a.LongitudinalVelocity_2UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.LongitudinalVelocity_2,'Tooltip','The phase velocity of longitudinal bulk waves propagating in the 2-direction.','Position',[200 120+60 60 23],'BackgroundColor','white');
        a.FastShearVelocity_2UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.FastShearVelocity_2,'Tooltip','The phase velocity of fast shear bulk waves propagating in the 2-direction.','Position',[200 120+30 60 23],'BackgroundColor','white');
        a.SlowShearVelocity_2UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.SlowShearVelocity_2,'Tooltip','The phase velocity of slow shear bulk waves propagating in the 2-direction.','Position',[200 120+0 60 23],'BackgroundColor','white');
        a.LongitudinalVelocity_3UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.LongitudinalVelocity_3,'Tooltip','The phase velocity of longitudinal bulk waves propagating in the 3-direction.','Position',[275 120+60 60 23],'BackgroundColor','white');
        a.FastShearVelocity_3UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.FastShearVelocity_3,'Tooltip','The phase velocity of fast shear bulk waves propagating in the 3-direction.','Position',[275 120+30 60 23],'BackgroundColor','white');
        a.SlowShearVelocity_3UI5 = uicontrol('Parent',a.p2UI5,'Style','text','String',a.Material2{1}.SlowShearVelocity_3,'Tooltip','The phase velocity of slow shear bulk waves propagating in the 3-direction.','Position',[275 120+0 60 23],'BackgroundColor','white');

        uicontrol('Parent',a.p2UI5,'Style','text','HorizontalAlignment','left','String','New material''s name','Position',[10 65 102 13]);
        a.Name2UI5 = uicontrol('Parent',a.p2UI5,'Style','edit','String',a.MaterialName4B,'Tooltip','Enter a new material''s name and press ''Save'' or enter an available material''s name and press ''Delete''.','Position',[125 60 220 23],'Callback',@CallbackUI5,'Tag','55');
        uicontrol('Parent',a.p2UI5,'Style','pushbutton','String','Save material','Tooltip','Save the material to the corresponding material list.','Position',[125 15 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','56');
        uicontrol('Parent',a.p2UI5,'Style','pushbutton','String','Delete material','Tooltip','Remove the material from the corresponding material list.','Position',[245 15 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','57');
        
        %------------------------------------------------------------------
        a.p3UI5 = uipanel('Parent',Tab5,'Title','Fluids','Units','pixels','Position',[10 10 440 205],'FontSize',10);
        uicontrol('Parent',a.p3UI5,'Style','text','HorizontalAlignment','left','String','Fluid','Position',[10 160 24 13]);
        uicontrol('Parent',a.p3UI5,'Style','text','HorizontalAlignment','left','String','Mass density (kg/m3)','Position',[10 130 105 13]);
        uicontrol('Parent',a.p3UI5,'Style','text','HorizontalAlignment','left','String','Velocity (m/s)','Position',[10 100 69 13]);
        a.Material3UI5 = uicontrol('Parent',a.p3UI5,'Style','popupmenu','String',a.FieldNames_Fluid,'Tooltip','Select a fluid to display its parameters.','Position',[135 155 220 23],'Callback',@CallbackUI5,'Tag','58');
        a.Density3UI5 = uicontrol('Parent',a.p3UI5,'Style','edit','String',a.Couplant1.Density,'Tooltip','Enter the mass density.','Position',[135 125 65 23],'Callback',@CallbackUI5,'Tag','59');
        a.VelocityUI5 = uicontrol('Parent',a.p3UI5,'Style','edit','String',a.Couplant1.Velocity,'Tooltip','Enter the phase velocity of ultrasonic waves.','Position',[135 95 65 23],'Callback',@CallbackUI5,'Tag','60');

        uicontrol('Parent',a.p3UI5,'Style','text','HorizontalAlignment','left','String','New fluid''s name','Position',[10 65 85 13]);
        a.Name3UI5 = uicontrol('Parent',a.p3UI5,'Style','edit','String',a.MaterialName4C,'Tooltip','Enter a new fluid''s name and press ''Save'' or enter an available fluid''s name and press ''Delete''.','Position',[135 60 220 23],'Callback',@CallbackUI5,'Tag','61');
        uicontrol('Parent',a.p3UI5,'Style','pushbutton','String','Save fluid','Tooltip','Save the fluid to the fluids list.','Position',[135 15 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','62');
        uicontrol('Parent',a.p3UI5,'Style','pushbutton','String','Delete fluid','Tooltip','Remove the fluid from the fluids list.','Position',[255 15 100 33],'FontSize',10,'Callback',@CallbackUI5,'Tag','63');
    end
    if  1 % Tab7_advanced
        a.c1UI7 = uicontrol('Parent',Tab7,'Style','text','HorizontalAlignment','left','String','Isotropic','Position',[280 740 54 15],'FontSize',10,'Backgroundcolor','white');
        a.c2UI7 = uicontrol('Parent',Tab7,'Style','text','HorizontalAlignment','left','String','Anisotropic','Position',[360 740 68 15],'FontSize',10,'Backgroundcolor','white');
        
        a.p1UI7 = uipanel('Parent',Tab7,'Title','Phase velocity sweeps','Units','pixels','Position',[10 615 430 120],'FontSize',10,'Backgroundcolor','white');    
        uicontrol('Parent',a.p1UI7,'Style','text','HorizontalAlignment','left','String','Phase velocity sections (5^x)','Position',[10 75 144 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p1UI7,'Style','text','HorizontalAlignment','left','String','Search width for negative curvature','Position',[10 45 237 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p1UI7,'Style','text','HorizontalAlignment','left','String','Search width for positive curvature','Position',[10 15 233 13],'Backgroundcolor','white');
        a.PhaseVelocitySections1UI7 = uicontrol('Parent',a.p1UI7,'Style','edit','String',a.PhaseVelocitySections1,'Tooltip','The dispersion curves are determined by sweeping the phase velocity at fixed frequencies for the most part. The phase velocity sections as the power of five give the maximum number of sections into which the phase velocity search interval is devided during the search for the modal solution at a given frequency. A higher number increases the chance to find the solution at the cost of processing time. Missing a solution is not necessarily critical since an extrapolation routine replaces missing samples successfully, as long as not too many samples are missing.','Position',[270 70 50 23],'Callback',@CallbackUI7,'Tag','18');
        a.PhaseVelocitySections2UI7 = uicontrol('Parent',a.p1UI7,'Style','edit','String',a.PhaseVelocitySections2,'Tooltip','The dispersion curves are determined by sweeping the phase velocity at fixed frequencies for the most part. The phase velocity sections as the power of five give the maximum number of sections into which the phase velocity search interval is devided during the search for the modal solution at a given frequency. A higher number increases the chance to find the solution at the cost of processing time. Missing a solution is not necessarily critical since an extrapolation routine replaces missing samples successfully, as long as not too many samples are missing.','Position',[350 70 50 23],'Callback',@CallbackUI7,'Tag','20');
        a.LambPhaseVelocitySweepRange11UI7 = uicontrol('Parent',a.p1UI7,'Style','edit','String',a.LambPhaseVelocitySweepRange11,'Tooltip','This value determines the phase velocity search interval width for Lamb waves. Please read the manual for more information.','Position',[270 40 50 23],'Callback',@CallbackUI7,'Tag','3');
        a.LambPhaseVelocitySweepRange21UI7 = uicontrol('Parent',a.p1UI7,'Style','edit','String',a.LambPhaseVelocitySweepRange21,'Tooltip','This value determines the phase velocity search interval width for Lamb waves. Please read the manual for more information.','Position',[270 10 50 23],'Callback',@CallbackUI7,'Tag','4');
        a.LambPhaseVelocitySweepRange12UI7 = uicontrol('Parent',a.p1UI7,'Style','edit','String',a.LambPhaseVelocitySweepRange12,'Tooltip','This value determines the phase velocity search interval width for Lamb waves. Please read the manual for more information.','Position',[350 40 50 23],'Callback',@CallbackUI7,'Tag','5');
        a.LambPhaseVelocitySweepRange22UI7 = uicontrol('Parent',a.p1UI7,'Style','edit','String',a.LambPhaseVelocitySweepRange22,'Tooltip','This value determines the phase velocity search interval width for Lamb waves. Please read the manual for more information.','Position',[350 10 50 23],'Callback',@CallbackUI7,'Tag','6');
        
        a.p2UI7 = uipanel('Parent',Tab7,'Title','Frequency sweeps to complete dispersion curves at high phase velocity','Units','pixels','Position',[10 485 430 120],'FontSize',10,'Backgroundcolor','white');
        uicontrol('Parent',a.p2UI7,'Style','text','HorizontalAlignment','left','String','Frequency sections (5^x)','Position',[10 75 126 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p2UI7,'Style','text','HorizontalAlignment','left','String','Phase velocity step (m/s)','Position',[10 45 124 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p2UI7,'Style','text','HorizontalAlignment','left','String','Search interval (kHz/mm)','Position',[10 15 123 13],'Backgroundcolor','white');
        a.FrequencySections1UI7 = uicontrol('Parent',a.p2UI7,'Style','edit','String',a.FrequencySections1,'Tooltip','Frequency sweeps at certain phase velocities are performed to complete the higher order dispersion curves at high phase velocities, if necessary. The frequency sections work similarly as the phase velocity sections, but with respect to the frequency search interval.','Position',[270 70 50 23],'Callback',@CallbackUI7,'Tag','19');
        a.FrequencySections2UI7 = uicontrol('Parent',a.p2UI7,'Style','edit','String',a.FrequencySections2,'Tooltip','Frequency sweeps at certain phase velocities are performed to complete the higher order dispersion curves at high phase velocities, if necessary. The frequency sections work similarly as the phase velocity sections, but with respect to the frequency search interval.','Position',[350 70 50 23],'Callback',@CallbackUI7,'Tag','21');
        a.PhaseVelocityStep1UI7 = uicontrol('Parent',a.p2UI7,'Style','edit','String',a.PhaseVelocityStep1,'Tooltip','Enter the phase velocity step.','Position',[270 40 50 23],'Callback',@CallbackUI7,'Tag','8');
        a.PhaseVelocityStep2UI7 = uicontrol('Parent',a.p2UI7,'Style','edit','String',a.PhaseVelocityStep2,'Tooltip','Enter the phase velocity step.','Position',[350 40 50 23],'Callback',@CallbackUI7,'Tag','12');
        a.FrequencyOffset1UI7 = uicontrol('Parent',a.p2UI7,'Style','edit','String',a.FrequencyOffset1,'Tooltip','Enter the frequency search interval width.','Position',[270 10 50 23],'Callback',@CallbackUI7,'Tag','10');
        a.FrequencyOffset2UI7 = uicontrol('Parent',a.p2UI7,'Style','edit','String',a.FrequencyOffset2,'Tooltip','Enter the frequency search interval width.','Position',[350 10 50 23],'Callback',@CallbackUI7,'Tag','15');

        a.p3UI7 = uipanel('Parent',Tab7,'Title','2-D tracing settings','Units','pixels','Position',[450 615 270 120],'FontSize',10,'Backgroundcolor','white');    
        uicontrol('Parent',a.p3UI7,'Style','text','HorizontalAlignment','left','String','Search width','Position',[10 75 112 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p3UI7,'Style','text','HorizontalAlignment','left','String','Search area sections','Position',[10 45 106 13],'Backgroundcolor','white');
        uicontrol('Parent',a.p3UI7,'Style','text','HorizontalAlignment','left','String','Search area extensions','Position',[10 15 118 13],'Backgroundcolor','white');
        a.SearchWidthUI7 = uicontrol('Parent',a.p3UI7,'Style','edit','String',['[',num2str(a.SearchWidth(1)),' ',num2str(a.SearchWidth(2)),']'],'Tooltip','This determines the search range. For more information, read the manual.','Position',[170 70 70 23],'Callback',@CallbackUI7,'Tag','22');
        a.SearchAreaSectionsUI7 = uicontrol('Parent',a.p3UI7,'Style','edit','String',a.SearchAreaSections,'Tooltip','This value determines at how many grid points of the search area, spanned by the real and imaginary wavenumber parts, the characteristic function is evaluated for finding a modal solution.','Position',[170 40 70 23],'Callback',@CallbackUI7,'Tag','23');
        a.SearchAreaExtensionsUI7 = uicontrol('Parent',a.p3UI7,'Style','edit','String',a.SearchAreaExtensions,'Tooltip','This value determines how many times the search area, spanned by the real and imaginary wavenumber parts, is increased in size in case a modal solution was not found previously.','Position',[170 10 70 23],'Callback',@CallbackUI7,'Tag','24');
        
        uicontrol('Parent',Tab7,'Style','text','HorizontalAlignment','left','String','DLR','Position',[170-53 50 40 23],'FontSize',13.25,'FontWeight','bold','Foregroundcolor',[.4 .4 .4],'Backgroundcolor','white');
        axes('Parent',Tab7,'Units','pixels','Position',[170-95 40 83 100]);
        imshow('DLR_Logo_gray.jpg')
        axes('Parent',Tab7,'Units','pixels','Position',[210 -50 1000 1000*.692]);
        imshow('Earth.jpg')
    end
end

f1.Units = 'normalized';
movegui(f1,'center') % center the GUI
f1.Visible = 'on'; % make the GUI visible

function New_Callback(~,~)
    Selection = questdlg('Do you really want to create a new project? The current project will be closed.','Confirm','Ok','Cancel','Ok'); 
    switch Selection
    case 'Ok'
        delete(f1)
        h = msgbox('Opening new DC instance. Please wait...');
        DispersionCalculator
        close(h)
    end
end
function Open_Callback(~,~)
    [Path,File] = uigetfile('*.mat');
    if  Path ~= 0
        h = msgbox('Loading data. Please wait...');
        load(fullfile(File,Path),'M');
        Fields = fieldnames(M);
        for i = 1:length(Fields)
            if  isstruct(M.(Fields{i})) && isfield(M.(Fields{i}),'Style')
                a.(Fields{i}).Style = M.(Fields{i}).Style;
                a.(Fields{i}).Value = M.(Fields{i}).Value;
                a.(Fields{i}).String = M.(Fields{i}).String;
                a.(Fields{i}).Enable = M.(Fields{i}).Enable;
                a.(Fields{i}).Position = M.(Fields{i}).Position;
            else
                a.(Fields{i}) = M.(Fields{i});
            end
        end
        clear M
        close(h)
        try
            Colorization(a,1)
            ShowModes_Signal(a)
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,1);
        catch
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
        end
        if  a.Multithreading
            try
                a.Pool = gcp('nocreate');
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to start parallel pool')
                return
            end
            if  isempty(a.Pool)
                h = msgbox('Starting parallel pool...');
                a.Pool = parpool('IdleTimeout',180);
                parfevalOnAll(@warning,0,'off','all')
                close(h)
            end
            imshow('Multithreading21.png','Parent',a.a1UI1)
            imshow('Multithreading21.png','Parent',a.a1UI2)
            imshow('Multithreading21.png','Parent',a.a1UI4)
        else
            imshow('Multithreading20.png','Parent',a.a1UI1)
            imshow('Multithreading20.png','Parent',a.a1UI2)
            imshow('Multithreading20.png','Parent',a.a1UI4)
        end
    end
end
function Save_Callback(~,~)
    h = msgbox('Gathering data. Please wait...');
    Fields = fieldnames(a);
    for i = 1:length(Fields)
        if  ~isobject(a.(Fields{i}))
            M.(Fields{i}) = a.(Fields{i});
        elseif isobject(a.(Fields{i})) && isvalid(a.(Fields{i})) && ~strcmp(Fields{i},'Pool') && strcmp(a.(Fields{i}).Type,'uicontrol')
            M.(Fields{i}).Style = a.(Fields{i}).Style;
            M.(Fields{i}).Value = a.(Fields{i}).Value;
            M.(Fields{i}).String = a.(Fields{i}).String;
            M.(Fields{i}).Enable = a.(Fields{i}).Enable;
            M.(Fields{i}).Position = a.(Fields{i}).Position;
        end
    end
    close(h)
    uisave('M',fullfile(a.Directory,'Project'))
    clear M
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
        msgbox([String,newline,'Please restart DC (File -> New project) after importing all material lists to have them available.'],'Info')
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
function OpenDirectory_Callback(~,~)
    try
        winopen(a.MaterialListDirectory)
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export material lists')
        return
    end
end
function Enable_Callback(~,~)
    try
        a.Pool = gcp('nocreate');
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to start parallel pool')
        return
    end
    if  isempty(a.Pool)
        h = msgbox('Starting parallel pool...');
        a.Pool = parpool('IdleTimeout',180);
        parfevalOnAll(@warning,0,'off','all')
        close(h)
    end
    a.Multithreading = 1;
    a.c4UI1.String = a.Pool.NumWorkers;
    a.c4UI2.String = a.Pool.NumWorkers;
    a.c4UI4.String = a.Pool.NumWorkers;
    imshow('Multithreading21.png','Parent',a.a1UI1)
    imshow('Multithreading21.png','Parent',a.a1UI2)
    imshow('Multithreading21.png','Parent',a.a1UI4)
end
function Disable_Callback(~,~)
    a.Multithreading = 0;
    a.c4UI1.String = '0';
    a.c4UI2.String = '0';
    a.c4UI4.String = '0';
    imshow('Multithreading20.png','Parent',a.a1UI1)
    imshow('Multithreading20.png','Parent',a.a1UI2)
    imshow('Multithreading20.png','Parent',a.a1UI4)
end
function Help_Callback(~,~)
    winopen 'DispersionCalculator_Manual.pdf'
end
function Homepage_Callback(~,~)
    try
        web('https://github.com/ArminHuber/Dispersion-Calculator')
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to open url')
        return
    end
end
function About_Callback(~,~)
    f3 = figure('Icon',which('DC_Logo16.png'),'NumberTitle','off','Name','About','Visible','off','MenuBar','none','Position',[0 0 520 780 ],'color','w');
    f3.Units = 'normalized';
    movegui(f3,'center')
    f3.Visible = 'on';

    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Dispersion Calculator','Position',[130 440+272 182 23],'FontSize',14,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','3','Position',[310 440+272 10 23],'FontSize',10,'FontWeight','bold','Foregroundcolor',[1 .5 0],'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String',['Version ',Version,' - compiled ',Date],'Position',[130 440+238 300 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String',['Copyright ',char(169),' 2018-20',Date(end-1:end),' DLR'],'Position',[130 440+220 200 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Dispersion Calculator was created by Armin Huber','Position',[130 440+186 350 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Center for Lightweight Production Technology (ZLP)','Position',[130 440+168 350 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Institute of Structures and Design','Position',[130 440+150 250 18],'FontSize',11,'Backgroundcolor','white');

    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Deutsches Zentrum','Position',[130 440+101 200 23],'FontSize',13.25,'FontWeight','bold','Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String',['f',char(252),'r Luft- und Raumfahrt'],'Position',[130 440+79 250 23],'FontSize',13.25,'FontWeight','bold','Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','DLR','Position',[130-53 440+79 50 23],'FontSize',13.25,'FontWeight','bold','Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','German Aerospace Center','Position',[130 440+56.5 250 23],'FontSize',13.25,'Backgroundcolor','white');

    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','For more information contact Armin Huber at','Position',[130 458 300 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','armin.huber@dlr.de','Position',[130 440 150 18],'FontSize',11,'Backgroundcolor','white');

    axes('Parent',f3,'Units','pixels','Position',[130-80 645 60 60]);
    imshow('DC_Logo60.png')
    axes('Parent',f3,'Units','pixels','Position',[130-95 509 83 100]);
    imshow('DLR_Logo_black.jpg')
    axes('Parent',f3,'Units','pixels','Position',[176 250 166 166]);
    imshow('DC_Logo_Earth.png') 
    
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Thanks to','Position',[130 210 100 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Michael Lowe, Michel Castaings,','Position',[130 192 250 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Stanislav Rokhlin, Victor Giurgiutiu,','Position',[130 174 250 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Marc Deschamps, Eric Ducasse,','Position',[130 156 250 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Markus Sause.','Position',[130 138 150 18],'FontSize',11,'Backgroundcolor','white');

    uicontrol('Parent',f3,'Style','text','String','Dispersion Calculator is dedicated to','Position',[130 80 260 18],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f3,'Style','text','String','Sir Patrick Stewart','Position',[130 55 260 18],'FontSize',11,'Backgroundcolor','white');
end
function CloseRequest(~,~)
    Selection = questdlg('Do you really want to close Dispersion Calculator?','Close request','Yes','No','Yes'); 
    switch Selection
    case 'Yes'
        delete(f1)
    end
end
function SpecimenSettingsUI2_Callback(~,~) % Tab2_anisotropic
    Geometry = a.Geometry2;
    Radius = a.Radius;
    InnerRadius = a.InnerRadius;
    Sink = a.Sink2;
    FluidLoading = a.FluidLoading2;
    if  strcmp(Geometry,'Plate')
        ToggleFluid1 = a.ToggleUpperFluid2;
        ToggleFluid2 = a.ToggleLowerFluid2;
        Fluid1 = a.UpperFluid2;
        Fluid2 = a.LowerFluid2;
    else
        ToggleFluid1 = a.ToggleLowerFluid2;
        ToggleFluid2 = a.ToggleUpperFluid2;
        Fluid1 = a.LowerFluid2;
        Fluid2 = a.UpperFluid2;
    end
    Couplant = a.Couplant2;
    MaterialType = a.MaterialType2;
    Material = a.Material2;
    MaterialNames = a.MaterialNames;
    MaterialClasses = a.MaterialClasses;
    Viscoelastic = a.Viscoelastic2;
    IsotropicSinglelayer = a.IsotropicSinglelayer;
    IsotropicMultilayerCylinder = a.IsotropicMultilayerCylinder;
    Hybrid = a.Hybrid;
    UniformLayerThickness = a.UniformLayerThickness;
    LayupThickness = a.LayupThickness;
    PlateThickness = a.PlateThickness;
    SuperLayers = a.SuperLayers;
    Pattern = a.Pattern;
    SuperLayerSize = a.SuperLayerSize;
    SymmetricSystem = a.SymmetricSystem;
    Symmetric = a.Symmetric2;
    Layers = a.Layers2;
    LayerOrientations = a.LayerOrientations;
    LayerThicknesses = a.LayerThicknesses;
    LayupString = a.LayupString2;
    LayupStringTex = a.LayupStringTex2;
    EffectiveLayupString = a.EffectiveLayupString2;
    UnitCell = a.UnitCell;

    Material = Material(1:SuperLayerSize);
    if  ~strcmp(Geometry,'Plate') && SymmetricSystem && SuperLayers == 1
        LayerThicknesses = [LayerThicknesses(1:SuperLayerSize-1) LayerThicknesses(SuperLayerSize)/2];
    else
        LayerThicknesses = LayerThicknesses(1:SuperLayerSize);
    end

    f2 = figure('Icon',which('DC_Logo16.png'),'NumberTitle','off','Name','Specimen','Visible','off','MenuBar','none','Position',[0 0 1010 460]);
    f2.Units = 'normalized';
    movegui(f2,'center')
    f2.Visible = 'on';

    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Geometry','Position',[20 435 49 13])
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Inner radius (mm)','Position',[20 405 100 13])
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Sink at center','Position',[20 375 68 13]);
    Fluid1TextUI2 = uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','Position',[20 345 54 13]);
    Fluid2TextUI2 = uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','Position',[20 285 57 13]);
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Layup','Position',[20 315 54 13])
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Hybrid','Position',[20 255 33 13])
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Class','Position',[20 225 29 13])
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Material','Position',[20 195 39 13]);
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Uniform layer thickness','Position',[20 165 115 13]);
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Layup thickness (mm)','Position',[20 135 101 13]);
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Unit cell repetitions','Position',[20 105 92 13]);
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Symmetric unit cell','Position',[20 75 90 13]);

    GeometryUI = uicontrol('Parent',f2,'Style','popupmenu','Value',find(strcmp(Geometry,{'Plate','Rod','Pipe','Circumferential'})),'String',{'Plate','Rod','Pipe','Circumferential'},'Tooltip','Select a geometry.','Position',[120 430 140 23],'Callback',@Callback,'Tag','1');
    InnerRadiusUI = uicontrol('Parent',f2,'Style','edit','String',InnerRadius,'Tooltip','Enter the inner pipe radius.','Position',[120 400 50 23],'Enable','off','Callback',@Callback,'Tag','2');
    SinkUI = uicontrol('Parent',f2,'Style','checkbox','Value',Sink,'Tooltip','Put a sink at the center of the pipe. This absorbs the bulk waves travelling from one side of the inner pipe wall through the fluid-filling to the other side, and therefore suppresses the fluid modes.','Position',[90 370 50 23],'Enable','off','Callback',@Callback,'Tag','3');
    ToggleFluid1UI = uicontrol('Parent',f2,'Style','checkbox','Value',ToggleFluid1,'Position',[90 340 50 23],'Callback',@Callback,'Tag','4');
    ToggleFluid2UI = uicontrol('Parent',f2,'Style','checkbox','Value',ToggleFluid2,'Position',[90 280 50 23],'Callback',@Callback,'Tag','5');
    SelectFluid1UI = uicontrol('Parent',f2,'Style','popupmenu','Value',find(strcmp(Fluid1.Name,a.FieldNames_Fluid)),'String',a.FieldNames_Fluid,'Position',[120 340 140 23],'Callback',@Callback,'Tag','6');
    SelectFluid2UI = uicontrol('Parent',f2,'Style','popupmenu','Value',find(strcmp(Fluid2.Name,a.FieldNames_Fluid)),'String',a.FieldNames_Fluid,'Position',[120 280 140 23],'Callback',@Callback,'Tag','7');
    LayupUI = uicontrol('Parent',f2,'Style','text','String',LayupString,'Position',[120 316 140 13],'BackgroundColor','white');
    HybridUI = uicontrol('Parent',f2,'Style','checkbox','Value',Hybrid,'Tooltip','The layup may contain different materials.','Position',[90 250 50 23],'Callback',@Callback,'Tag','8');
    MaterialClassUI = uicontrol('Parent',f2,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic','Isotropic'},'Value',MaterialType,'Tooltip','Select a material symmetry class.','Position',[90 220 170 23],'Callback',@Callback,'Tag','9');
    MaterialUI = uicontrol('Parent',f2,'Style','popupmenu','String',{''},'Tooltip','Select a material.','Position',[90 190 170 23],'Callback',@Callback,'Tag','10');
    UniformLayerThicknessUI = uicontrol('Parent',f2,'Style','checkbox','Value',UniformLayerThickness,'Tooltip','Check this if the layup has uniform layer thicknesses. Dispersion Calculator deduces the layer thicknesses from the total layup thickness as defined below, and you do not need to enter the individual layer thicknesses into the table to the right.','Position',[150 160 50 23],'Callback',@Callback,'Tag','11');
    LayupThicknessUI = uicontrol('Parent',f2,'Style','edit','String',LayupThickness,'Tooltip','Enter the total layup thickness. This is available only if you have checked ''Uniform layer thickness''.','Position',[150 130 50 23],'Callback',@Callback,'Tag','12');
    SuperLayersUI = uicontrol('Parent',f2,'Style','edit','String',SuperLayers,'Tooltip','Enter the number of repetitions of the unit cell contained in the layup. Use only integers.','Position',[150 100 50 23],'Callback',@Callback,'Tag','13');
    SymmetricSystemUI = uicontrol('Parent',f2,'Style','checkbox','Value',SymmetricSystem,'Tooltip','Check this if the layup is symmetric about its middle plane. Dispersion Calculator extends the layup accordingly.','Position',[150 70 50 23],'Callback',@Callback,'Tag','14');    

    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Unit cell','Position',[280 370 46 13],'FontSize',9);
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','Note: When setting up a plate, the first layer is the top layer.','Position',[380 385 500 13])
    uicontrol('Parent',f2,'Style','text','HorizontalAlignment','left','String','When setting up a rod, pipe, or circumference, the first layers is the innermost layer.','Position',[408 370 500 13])
    TableUI = uitable(f2,'ColumnName',{['Phi (',char(176),')'],'d (mm)','Orthotropic','Trans. iso.','Cubic','Isotropic','Delete'},'ColumnWidth',{60 60 120 120 120 120 55},'ColumnFormat',({[] [] a.FieldNames_Orthotropic' a.FieldNames_TransverselyIsotropic' a.FieldNames_Cubic' a.FieldNames_Isotropic' 'logical'}),'Tooltip','Enter the fiber orientation of every layer into the left column, and if you don''t have checked ''Uniform layer thickness'', enter also the corresponding layer thicknesses into the second column.','Data',UnitCell,'RowStriping','off','Position',[280 20 715 340],'CellEditCallback',@TableCallback);

    uicontrol('Parent',f2,'Style','pushbutton','String','Reset','Tooltip','Reset specimen definition.','Position',[280 420 95 33],'FontSize',10,'Callback',@Callback,'Tag','15'); 
    uicontrol('Parent',f2,'Style','pushbutton','String','OK','Tooltip','Accept the current specimen definition and continue.','Position',[20 20 95 33],'FontSize',10,'Callback',@Callback,'Tag','16');
    uicontrol('Parent',f2,'Style','pushbutton','String','Cancel','Tooltip','Discard the current specimen definition.','Position',[150 20 95 33],'FontSize',10,'Callback',@Callback,'Tag','17');

    GeometryAdjustments
    TableColumnEditable

    if  Hybrid
        MaterialClassUI.Enable = 'off';
        MaterialUI.Enable = 'off';
        MaterialUI.Value = 1;
    else
        MaterialClassUI.Enable = 'on';
        MaterialUI.Enable = 'on';
        if  MaterialType == 1
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Orthotropic));
        elseif MaterialType == 2
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_TransverselyIsotropic));
        elseif MaterialType == 3
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Cubic));
        elseif MaterialType == 4
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Isotropic));
        end
    end
    if  MaterialType == 1
        MaterialUI.String = a.FieldNames_Orthotropic;
    elseif MaterialType == 2
        MaterialUI.String = a.FieldNames_TransverselyIsotropic;
    elseif MaterialType == 3
        MaterialUI.String = a.FieldNames_Cubic;
    elseif MaterialType == 4
        MaterialUI.String = a.FieldNames_Isotropic;
    end
    if  UniformLayerThickness
        LayupThicknessUI.Enable = 'on';
    else
        LayupThicknessUI.Enable = 'off';
    end

    function GeometryAdjustments
        if  strcmp(Geometry,'Plate')
            Fluid1TextUI2.String = 'Upper fluid';
            Fluid2TextUI2.String = 'Lower fluid';
            InnerRadiusUI.Enable = 'off';
            SinkUI.Enable = 'off';
            ToggleFluid1UI.Enable = 'on';
            ToggleFluid2UI.Enable = 'on';
            if  ToggleFluid1
                SelectFluid1UI.Enable = 'on';
            else
                SelectFluid1UI.Enable = 'off';
            end
            if  ToggleFluid2
                SelectFluid2UI.Enable = 'on';
            else
                SelectFluid2UI.Enable = 'off';
            end
            ToggleFluid1UI.Tooltip = 'Toggle the upper fluid. If unchecked, the upper half-space will be vacuum.';
            ToggleFluid2UI.Tooltip = 'Toggle the lower fluid. If unchecked, the lower half-space will be vacuum.';
            SelectFluid1UI.Tooltip = 'Select the fluid in the upper half-space.';
            SelectFluid2UI.Tooltip = 'Select the fluid in the lower half-space.';
        elseif strcmp(Geometry,'Rod')
            Fluid2TextUI2.String = 'Outer fluid';
            InnerRadiusUI.Enable = 'off';
            SinkUI.Enable = 'off';
            ToggleFluid1UI.Enable = 'off';
            ToggleFluid2UI.Enable = 'on';
            SelectFluid1UI.Enable = 'off';
            if  ToggleFluid2
                SelectFluid2UI.Enable = 'on';
            else
                SelectFluid2UI.Enable = 'off';
            end
            ToggleFluid2UI.Tooltip = 'Toggle the outer fluid. If unchecked, the rod will be surrounded by vacuum.';
            SelectFluid2UI.Tooltip = 'Select the fluid surrounding the rod.';
        elseif strcmp(Geometry,'Pipe')
            Fluid1TextUI2.String = 'Inner fluid';
            Fluid2TextUI2.String = 'Outer fluid';
            InnerRadiusUI.Enable = 'on';
            ToggleFluid1UI.Enable = 'on';
            ToggleFluid2UI.Enable = 'on';
            if  ToggleFluid1
                SelectFluid1UI.Enable = 'on';
                SinkUI.Enable = 'on';
            else
                SelectFluid1UI.Enable = 'off';
                SinkUI.Enable = 'off';
            end
            if  ToggleFluid2
                SelectFluid2UI.Enable = 'on';
            else
                SelectFluid2UI.Enable = 'off';
            end
            ToggleFluid1UI.Tooltip = 'Toggle the inner fluid. If unchecked, the pipe will contain vacuum inside.';
            ToggleFluid2UI.Tooltip = 'Toggle the outer fluid. If unchecked, the pipe will be surrounded by vacuum.';
            SelectFluid1UI.Tooltip = 'Select the fluid inside the pipe.';
            SelectFluid2UI.Tooltip = 'Select the fluid surrounding the pipe.';
        elseif strcmp(Geometry,'Circumferential')
            Fluid1TextUI2.String = 'Inner fluid';
            Fluid2TextUI2.String = 'Outer fluid';
            InnerRadiusUI.Enable = 'on';
            SinkUI.Enable = 'off';
            ToggleFluid1UI.Enable = 'off';
            ToggleFluid2UI.Enable = 'off';
            SelectFluid1UI.Enable = 'off';
            SelectFluid2UI.Enable = 'off';
            ToggleFluid1UI.Tooltip = 'Toggle the inner fluid. If unchecked, the pipe will contain vacuum inside.';
            ToggleFluid2UI.Tooltip = 'Toggle the outer fluid. If unchecked, the pipe will be surrounded by vacuum.';
            SelectFluid1UI.Tooltip = 'Select the fluid inside the pipe.';
            SelectFluid2UI.Tooltip = 'Select the fluid surrounding the pipe.';
        end
    end
    function TableColumnEditable
        if  Hybrid && UniformLayerThickness
            TableUI.ColumnEditable = [true false true true true true true];
        elseif Hybrid && ~UniformLayerThickness
            TableUI.ColumnEditable = [true true true true true true true];
        elseif ~Hybrid && UniformLayerThickness
            TableUI.ColumnEditable = [true false false false false false true];
        elseif ~Hybrid && ~UniformLayerThickness
            TableUI.ColumnEditable = [true true false false false false true];
        end
if  Hybrid && ~strcmp(Geometry,'Plate')
    TableUI.ColumnEditable(3:5) = false;
end
    end
    function DistributeMaterials
        if  MaterialType == 1
            Material(1:SuperLayerSize) = {getfield(a.Materials.Orthotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 2
            Material(1:SuperLayerSize) = {getfield(a.Materials.TransverselyIsotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 3
            Material(1:SuperLayerSize) = {getfield(a.Materials.Cubic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 4
            Material(1:SuperLayerSize) = {getfield(a.Materials.Isotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        end
    end
    function ClearTableColumns
        if  MaterialType == 1
            TableUI.Data(1:SuperLayerSize,4:6) = {''};
        elseif MaterialType == 2
            TableUI.Data(1:SuperLayerSize,[3 5 6]) = {''};
        elseif MaterialType == 3
            TableUI.Data(1:SuperLayerSize,[3 4 6]) = {''};
        elseif MaterialType == 4
            TableUI.Data(1:SuperLayerSize,3:5) = {''};
        end
        TableUI.Data(1:SuperLayerSize,MaterialType+2) = {Material{1}.Name};
    end
    function SetThicknesses
        if  UniformLayerThickness
            if  SymmetricSystem
                LayerThicknesses(1:SuperLayerSize) = LayupThickness/SuperLayerSize/SuperLayers/2;
            else
                LayerThicknesses(1:SuperLayerSize) = LayupThickness/SuperLayerSize/SuperLayers;
            end
            TableUI.Data(1:SuperLayerSize,2) = {LayerThicknesses(1)};
        else
            LayupThickness = SuperLayers*sum(LayerThicknesses);
            if  SymmetricSystem
                LayupThickness = LayupThickness*2;
            end
            LayupThicknessUI.String = LayupThickness;
        end
    end
    function MakeLayupString
        LayupString = ['[',replace(num2str(LayerOrientations),whitespacePattern,'/'),']'];
        LayupStringTex = [' ',LayupString];
        EffectiveLayupString = ['[',replace(num2str(LayerOrientations-a.PropagationAngle),whitespacePattern,'/'),']'];
        if  SuperLayers > 1
            LayupString = [LayupString,num2str(SuperLayers)];
            LayupStringTex = [LayupStringTex,'$_{',num2str(SuperLayers),'}$'];
            EffectiveLayupString = [EffectiveLayupString,num2str(SuperLayers)];
        end
        if  SymmetricSystem
            LayupString = [LayupString,'s'];
            LayupStringTex = [LayupStringTex,'$_{\mathrm s}$'];
            LayupStringTex = erase(LayupStringTex,'}$$_{');
            EffectiveLayupString = [EffectiveLayupString,'s'];
        end
    end
    function Reset
        Geometry = 'Plate';
        Radius = 0;
        InnerRadius = 1;
        Sink = 0;
        FluidLoading = 0;
        ToggleFluid1 = 0;
        ToggleFluid2 = 0;
        Fluid1 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        Fluid2 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
        Couplant = Fluid1;
        MaterialType = 1;
        Material = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})};
        MaterialNames = {Material{1}.Name};
        MaterialClasses = {Material{1}.Class};
        if  isreal(Material{1}.C)
            Viscoelastic = 0;
        else
            Viscoelastic = 1;
        end
        IsotropicSinglelayer = 0;
        IsotropicMultilayerCylinder = 0;
        Hybrid = 0;
        UniformLayerThickness = 1;
        LayupThickness = 1;
        PlateThickness = LayupThickness;
        SuperLayers = 1;
        Pattern = [];
        SuperLayerSize = 1;
        SymmetricSystem = 0;
        Symmetric = 1;
        Layers = 1;
        LayerOrientations = 0;
        LayerThicknesses = 1;
        LayupString = '[0]';
        LayupStringTex = '[0]';
        EffectiveLayupString = '[0]';
        UnitCell = cell(99,7);
        UnitCell{1} = 0;
        UnitCell{1,2} = 1;
        UnitCell{1,3} = a.FieldNames_Orthotropic{1};
        Fluid1TextUI2.String = 'Upper fluid';
        Fluid2TextUI2.String = 'Lower fluid';
        GeometryUI.Value = 1;
        InnerRadiusUI.String = InnerRadius;
        InnerRadiusUI.Enable = 'off';
        SinkUI.Value = Sink;
        SinkUI.Enable = 'off';
        ToggleFluid1UI.Value = ToggleFluid1;
        ToggleFluid2UI.Value = ToggleFluid2;
        ToggleFluid1UI.Enable = 'on';
        ToggleFluid2UI.Enable = 'on';
        ToggleFluid1UI.Tooltip = 'Toggle the upper fluid. If unchecked, the upper half-space will be vacuum.';
        ToggleFluid2UI.Tooltip = 'Toggle the lower fluid. If unchecked, the lower half-space will be vacuum.';
        SelectFluid1UI.Value = 1;
        SelectFluid2UI.Value = 1;
        SelectFluid1UI.Enable = 'off';
        SelectFluid2UI.Enable = 'off';
        SelectFluid1UI.Tooltip = 'Select the fluid in the upper half-space.';
        SelectFluid2UI.Tooltip = 'Select the fluid in the lower half-space.';
        LayupUI.String = LayupString;
        HybridUI.Value = Hybrid;
        MaterialClassUI.Value = MaterialType;
        MaterialClassUI.Enable = 'on';
        MaterialUI.String = a.FieldNames_Orthotropic;
        MaterialUI.Value = 1;
        MaterialUI.Enable = 'on';
        UniformLayerThicknessUI.Value = UniformLayerThickness;
        LayupThicknessUI.String = LayupThickness;
        LayupThicknessUI.Enable = 'on';
        SuperLayersUI.String = SuperLayers;
        SymmetricSystemUI.Value = SymmetricSystem;
        TableUI.ColumnEditable = [true false false false false false true];
        TableUI.Data = UnitCell;
    end
function Reset2
    Radius = 1;
    InnerRadius = 1;
    Sink = 0;
    FluidLoading = 0;
    ToggleFluid1 = 0;
    ToggleFluid2 = 0;
    Fluid1 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
    Fluid2 = getfield(a.Materials.Fluid,a.FieldNames_Fluid{1});
    Couplant = Fluid1;
    MaterialType = 4;
    Material = {getfield(a.Materials.Isotropic,a.FieldNames_Isotropic{1})};
    MaterialNames = {Material{1}.Name};
    MaterialClasses = {Material{1}.Class};
    if  isreal(Material{1}.C)
        Viscoelastic = 0;
    else
        Viscoelastic = 1;
    end
    IsotropicSinglelayer = 1;
    IsotropicMultilayerCylinder = 1;
    Hybrid = 0;
    UniformLayerThickness = 1;
    LayupThickness = 1;
    PlateThickness = 2*Radius;
    SuperLayers = 1;
    Pattern = [];
    SuperLayerSize = 1;
    SymmetricSystem = 0;
    Symmetric = 1;
    Layers = 1;
    LayerOrientations = 0;
    LayerThicknesses = 1;
    LayupString = '[0]';
    LayupStringTex = '[0]';
    EffectiveLayupString = '[0]';
    UnitCell = cell(99,7);
    UnitCell{1} = 0;
    UnitCell{1,2} = 1;
    UnitCell{1,6} = a.FieldNames_Isotropic{1};
    InnerRadiusUI.String = InnerRadius;
    SinkUI.Value = Sink;
    ToggleFluid1UI.Value = ToggleFluid1;
    ToggleFluid2UI.Value = ToggleFluid2;
    SelectFluid1UI.Value = 1;
    SelectFluid2UI.Value = 1;
    LayupUI.String = LayupString;
    HybridUI.Value = Hybrid;
    MaterialClassUI.Value = MaterialType;
    MaterialClassUI.Enable = 'on';
    MaterialUI.String = a.FieldNames_Isotropic;
    MaterialUI.Value = 1;
    MaterialUI.Enable = 'on';
    UniformLayerThicknessUI.Value = UniformLayerThickness;
    LayupThicknessUI.String = LayupThickness;
    LayupThicknessUI.Enable = 'on';
    SuperLayersUI.String = SuperLayers;
    SymmetricSystemUI.Value = SymmetricSystem;
    TableUI.ColumnEditable = [true false false false false false true];
    TableUI.Data = UnitCell;
end
    function DisplayProperties
        String = ['Geometry:            ',Geometry,newline...
                  'LayerOrientations:   ',regexprep(num2str(LayerOrientations),'\s+',' '),newline...
                  'LayerThicknesses:    ',regexprep(num2str(LayerThicknesses),'\s+',' '),newline...
                  'Layers:              ',num2str(SuperLayerSize),newline...
                  'Total layers:        ',num2str(Layers),newline...
                  'LayupThickness:      ',num2str(LayupThickness),newline...
                  'Waveguide thickness: ',num2str(PlateThickness),newline];
        if  ~strcmp(Geometry,'Plate')
            String = append(String,'Radius:              ',regexprep(num2str(Radius'),'\s+',' '),newline);
        end
        for i = 1:length(Material)
            String = append(String,num2str(i),' ',Material{i}.Name,' ',Material{i}.Class,newline);
        end
        disp([String,'---------------------------'])
    end
    function Callback(source,~)
        source.String = replace(source.String,',','.');
        if  strcmp(source.Tag,'1') % Geometry
if  strcmp(Geometry,'Plate') && source.Value > 1
    Reset2
elseif ~strcmp(Geometry,'Plate') && source.Value == 1
    Reset
end
            Geometry = source.String{source.Value};
            GeometryAdjustments
        elseif strcmp(source.Tag,'2') % Inner radius
            if  str2double(source.String) <= 0
                source.String = 1;
                errordlg('The inner pipe radius must be greater than zero. Choose ''Rod'' as a geometry to model a zero inner radius.','Error');
            end
            InnerRadius = str2double(source.String);
        elseif strcmp(source.Tag,'3') % Sink at center
            Sink = source.Value;
        elseif strcmp(source.Tag,'4') % Toggle fluid 1
            ToggleFluid1 = source.Value;
            if  ToggleFluid1
                SelectFluid1UI.Enable = 'on';
                if  strcmp(Geometry,'Pipe')
                    SinkUI.Enable = 'on';
                end
            else
                SelectFluid1UI.Enable = 'off';
                if  strcmp(Geometry,'Pipe')
                    SinkUI.Enable = 'off';
                end
            end
        elseif strcmp(source.Tag,'5') % Toggle fluid 2
            ToggleFluid2 = source.Value;
            if  ToggleFluid2
                SelectFluid2UI.Enable = 'on';
            else
                SelectFluid2UI.Enable = 'off';
            end
        elseif strcmp(source.Tag,'6') % Select fluid 1
            Fluid1 = getfield(a.Materials.Fluid,cell2mat(source.String(source.Value)));
            if  strcmp(Geometry,'Plate')
                Couplant = Fluid1;
            end
        elseif strcmp(source.Tag,'7') % Select fluid 2
            Fluid2 = getfield(a.Materials.Fluid,cell2mat(source.String(source.Value)));
            if  strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
                Couplant = Fluid2;
            end
        elseif strcmp(source.Tag,'8') % Hybrid
            Hybrid = source.Value;
            if  Hybrid
                MaterialClassUI.Enable = 'off';
                MaterialUI.Enable = 'off';
            else
                MaterialClassUI.Enable = 'on';
                MaterialUI.Enable = 'on';
                DistributeMaterials
                ClearTableColumns
            end
            TableColumnEditable
        elseif strcmp(source.Tag,'9') % Class
if  ~strcmp(Geometry,'Plate') && source.Value < 4
    source.Value = 4;
    errordlg('Rod, Pipe, and Circumferential can consist of isotropic layers only.','Error');
    return
end
            MaterialType = source.Value;
            MaterialUI.Value = 1;
            if  MaterialType == 1
                MaterialUI.String = a.FieldNames_Orthotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})};
            elseif MaterialType == 2
                MaterialUI.String = a.FieldNames_TransverselyIsotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.TransverselyIsotropic,a.FieldNames_TransverselyIsotropic{1})};
            elseif MaterialType == 3
                MaterialUI.String = a.FieldNames_Cubic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Cubic,a.FieldNames_Cubic{1})};
            elseif MaterialType == 4
                MaterialUI.String = a.FieldNames_Isotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Isotropic,a.FieldNames_Isotropic{1})};
            end
            ClearTableColumns
        elseif strcmp(source.Tag,'10') % Material
            DistributeMaterials
            ClearTableColumns
        elseif strcmp(source.Tag,'11') % Uniform layer thickness
            UniformLayerThickness = source.Value;
            if  UniformLayerThickness
                LayupThicknessUI.Enable = 'on';
                SetThicknesses
            else
                LayupThicknessUI.Enable = 'off';
            end
            TableColumnEditable
        elseif strcmp(source.Tag,'12') % Layup thickness
            if  str2double(source.String) <= 0
                source.String = 1;
                errordlg('The layup thickness must be greater than zero.','Error');
            end
            LayupThickness = str2double(source.String);
            SetThicknesses
        elseif strcmp(source.Tag,'13') % Unit cell repetitions
            if  rem(str2double(source.String),1) ~= 0 || str2double(source.String) < 1
                source.String = ceil(abs(str2double(source.String)));
                errordlg('Unit cell repetitions must be a positive integer.','Error');
            end
            SuperLayers = str2double(source.String);
            SetThicknesses
            MakeLayupString
            LayupUI.String = LayupString;
        elseif strcmp(source.Tag,'14') % Symmetric unit cell
            SymmetricSystem = source.Value;
            SetThicknesses
            MakeLayupString
            LayupUI.String = LayupString;
        elseif strcmp(source.Tag,'15') % Reset
            Reset
        elseif strcmp(source.Tag,'16') % OK
            MaterialNames = {};
            MaterialClasses = {};
            for i = 1:SuperLayerSize
                MaterialNames{i} = Material{i}.Name;
                MaterialClasses{i} = Material{i}.Class;
            end
            if  SuperLayerSize == 1 && SuperLayers > 1
                errordlg('Do not define multiple repetitions if the layup consists of only one layer! Increase the layer thickness instead.','Error');
                SuperLayers = 1;
                SuperLayersUI.String = 1;
                LayerThicknesses = LayupThickness;
                TableUI.Data{1,2} = LayerThicknesses(1);
                MakeLayupString
                LayupUI.String = LayupString;
                return
            elseif SuperLayerSize == 1 && SymmetricSystem
                errordlg('Do not define a symmetric layup if the layup consists of only one layer! Set double layer thickness instead.','Error');
                SymmetricSystem = 0;
                SymmetricSystemUI.Value = 0;
                LayerThicknesses = LayupThickness;
                TableUI.Data{1,2} = LayerThicknesses(1);
                MakeLayupString
                LayupUI.String = LayupString;
                return
            elseif SuperLayerSize > 1 && all(strcmp(MaterialClasses,'Isotropic')) && all(strcmp(MaterialNames(1),MaterialNames))
                errordlg('Do not define a unit cell containing solely layers of the same isotropic material! Replace by one layer of greater thickness. This will be the same because DC assumes rigid bonding between the layers. It is also more efficient since the computational expense scales with the number of layers.','Error');
                return
            elseif SuperLayerSize > 1 && all(LayerOrientations(1) == LayerOrientations) && all(strcmp(MaterialNames(1),MaterialNames))
                errordlg('Do not define a unit cell containing solely layers of the same material and orientation! Replace by one layer of greater thickness. This will be the same because DC assumes rigid bonding between the layers. It is also more efficient since the computational expense scales with the number of layers.','Error');
                return
            elseif strcmp(Geometry,'Plate') && ~SymmetricSystem && SuperLayerSize > 1 &&...
                (~mod(SuperLayerSize,2) && all(fliplr(LayerOrientations(1:SuperLayerSize/2)) == LayerOrientations(SuperLayerSize/2+1:end)) && all(fliplr(LayerThicknesses(1:SuperLayerSize/2)) == LayerThicknesses(SuperLayerSize/2+1:end)) && all(strcmp(fliplr(MaterialNames(1:SuperLayerSize/2)),MaterialNames(SuperLayerSize/2+1:end))) ||...
                mod(SuperLayerSize,2) && all(fliplr(LayerOrientations(1:SuperLayerSize/2-.5)) == LayerOrientations(SuperLayerSize/2+1.5:end)) && all(fliplr(LayerThicknesses(1:SuperLayerSize/2-.5)) == LayerThicknesses(SuperLayerSize/2+1.5:end)) && all(strcmp(fliplr(MaterialNames(1:SuperLayerSize/2-.5)),MaterialNames(SuperLayerSize/2+1.5:end))))
                errordlg('Unlike in rods, pipes, and circumferentials, do not define a symmetric unit cell in a plate! Cut the unit cell in half and check ''Symmetric unit cell'' instead. This will be more robust and more efficient.','Error');
                return
            elseif all(strcmp(MaterialNames(1),MaterialNames)) && Hybrid
                if  SuperLayerSize > 1
                    errordlg('It is no hybrid if all layers are the same material!','Error');
                else
                    errordlg('It is no hybrid if you have only one layer!','Error');
                end
                Hybrid = 0;
                HybridUI.Value = 0;
                MaterialClassUI.Enable = 'on';
                MaterialUI.Enable = 'on';
                DistributeMaterials
                ClearTableColumns
                TableColumnEditable
                return
            end
            UnitCell = TableUI.Data;
            for i = 1:SuperLayerSize
                if  isreal(Material{i}.C)
                    Viscoelastic(i) = 0;
                else
                    Viscoelastic(i) = 1;
                end
            end
            if  any(Viscoelastic)
                Viscoelastic = 1;
            else
                Viscoelastic = 0;
            end
            if  SuperLayerSize == 1 && strcmp(MaterialClasses,'Isotropic')
                IsotropicSinglelayer = 1;
            else
                IsotropicSinglelayer = 0;
            end
            if  all(strcmp(MaterialClasses,'Isotropic')) && ~strcmp(Geometry,'Plate')
                IsotropicMultilayerCylinder = 1;
            else
                IsotropicMultilayerCylinder = 0;
            end
            if  strcmp(Geometry,'Plate')
                Layers = SuperLayerSize*SuperLayers;
                if  SymmetricSystem
                    Layers = Layers*2;
                end
                PlateThickness = LayupThickness;
                n = 2.^(0:log2(SuperLayers));
                Pattern = [0:length(n)-1;n];
                while Pattern(2,end) < SuperLayers
                    for i = length(n):-1:1
                        if  Pattern(2,end)+n(i) <= SuperLayers
                            Pattern(1:2,end+1) = [i;Pattern(2,end)+n(i)];
                            break
                        end
                    end
                end
                Pattern = Pattern(1,2:end);
            else
                Material = repmat(Material(1:SuperLayerSize),1,SuperLayers);
                LayerThicknesses = repmat(LayerThicknesses(1:SuperLayerSize),1,SuperLayers);
                if  SymmetricSystem
                    Material = [Material fliplr(Material(1:end-1))];
                    LayerThicknesses = [LayerThicknesses(1:end-1) 2*LayerThicknesses(end) fliplr(LayerThicknesses(1:end-1))];
                end
                Layers = length(LayerThicknesses);
                if  strcmp(Geometry,'Rod')
                    Radius = [0 cumsum(LayerThicknesses)]';
                else
                    Radius = [InnerRadius InnerRadius+cumsum(LayerThicknesses)]';
                end
                PlateThickness = 2*Radius(end);
            end
            if  (strcmp(Geometry,'Rod') && ToggleFluid2) || ((strcmp(Geometry,'Plate') || strcmp(Geometry,'Pipe')) && (ToggleFluid1 || ToggleFluid2))
                FluidLoading = 1;
                if  strcmp(Geometry,'Pipe') && ~ToggleFluid2
                    a.HalfspacesNumber1UI2.Enable = 'off';
                    a.Halfspaces1UI2.Enable = 'off';
                    a.HalfspacesNumber2UI2.Enable = 'off';
                    a.Halfspaces2UI2.Enable = 'off';
                else
                    if  a.Halfspaces12
                        a.HalfspacesNumber1UI2.Enable = 'on';
                    end
                    a.Halfspaces1UI2.Enable = 'on';
                    if  a.Halfspaces22
                        a.HalfspacesNumber2UI2.Enable = 'on';
                    end
                    a.Halfspaces2UI2.Enable = 'on';
                end
                if  a.Quantity12 == 4
                    a.Option1UI2.Value = find(strcmp(Couplant.Name,a.Option1UI2.String));
                    if  IsotropicSinglelayer
                        x = round(8e4/Material{1}.PlateVelocity*Couplant.Velocity/343);
                    else
                        x = round(30*Couplant.Velocity/343);
                    end
                    if  x > 90
                        x = 90;
                    end
                    a.YAxisUI2.String = ['[0 ',num2str(x),']'];
                    a.YAxis2 = eval(a.YAxisUI2.String);
                end
            else
                FluidLoading = 0;
                a.HalfspacesNumber1UI2.Enable = 'off';
                a.Halfspaces1UI2.Enable = 'off';
                a.HalfspacesNumber2UI2.Enable = 'off';
                a.Halfspaces2UI2.Enable = 'off';
            end
            if  (SuperLayerSize == 1 || SymmetricSystem) && ((~ToggleFluid1 && ~ToggleFluid2) || (ToggleFluid1 && ToggleFluid2 && strcmp(Fluid1.Name,Fluid2.Name)))
                Symmetric = 1;
            else
                Symmetric = 0;
            end
            if  IsotropicSinglelayer
                if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod') || strcmp(Geometry,'Circumferential')
                    a.XRange2 = 2;
                    a.XSamples2 = 1000;
                    a.BelowCutoffWidth2 = 50;
                elseif strcmp(Geometry,'Pipe')
                    if  ~FluidLoading && ~Viscoelastic
                        a.XRange2 = 1;
                        a.XSamples2 = 1000;
                        a.BelowCutoffWidth2 = 50;
                    else
                        if  ToggleFluid1 && ~Sink
                            a.XRange2 = 1;
                            a.BelowCutoffWidth2 = 12.5;
                        else
                            a.XRange2 = .4;
                            a.BelowCutoffWidth2 = 50;
                        end
                        a.XSamples2 = 400;
                    end
                end
            else
                if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod') || strcmp(Geometry,'Circumferential')
                    a.XRange2 = 4;
                    a.XSamples2 = 400;
                    a.BelowCutoffWidth2 = 50;
                elseif strcmp(Geometry,'Pipe')
                    if  ~FluidLoading && ~Viscoelastic
                        a.XRange2 = 2;
                        a.XSamples2 = 400;
                        a.BelowCutoffWidth2 = 50;
                    else
                        if  ToggleFluid1 && ~Sink
                            a.XRange2 = 2;
                            a.BelowCutoffWidth2 = 12.5;
                        else
                            a.XRange2 = .8;
                            a.BelowCutoffWidth2 = 50;
                        end
                        a.XSamples2 = 200;
                    end
                end
            end
            if  ToggleFluid1 && ~strcmp(Geometry,'Rod') && ~strcmp(Geometry,'Circumferential')
                a.Fluid1DisplayUI2.String = Fluid1.Name;
            else
                if  strcmp(Geometry,'Rod')
                    a.Fluid1DisplayUI2.String = '--';
                else
                    a.Fluid1DisplayUI2.String = 'vacuum';
                end
            end
            if  ToggleFluid2 && ~strcmp(Geometry,'Circumferential')
                a.Fluid2DisplayUI2.String = Fluid2.Name;
            else
                a.Fluid2DisplayUI2.String = 'vacuum';
            end
            if  Hybrid
                a.MaterialNameUI2.String = 'Hybrid';
            else
                a.MaterialNameUI2.String = Material{1}.Name;
            end
            MakeLayupString
            a.GeometryDisplayUI2.String = Geometry;
            a.LayupUI2.String = LayupString;
            if  strcmp(Geometry,'Plate')
                a.PropagationAngleUI2.Enable = 'on';
                a.EffectiveLayupUI2.String = EffectiveLayupString;
            else
                a.PropagationAngleUI2.Enable = 'off';
                a.EffectiveLayupUI2.String = LayupString;
            end
            if  ~strcmp(Geometry,'Plate') && SymmetricSystem
                a.LayerCountUI2.String = Layers+1;
            else
                a.LayerCountUI2.String = Layers;
            end
            if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod')
                a.ThicknessCountUI2.String = PlateThickness;
            else
                a.ThicknessCountUI2.String = [num2str(PlateThickness),char(215),num2str(LayupThickness)];
            end
            if  a.Quantity12 < 3
                if  Layers == 1
                    a.Option1UI2.Enable = 'on';
                else
                    a.Option1UI2.Enable = 'off';
                end
            end
            a.Samples12 = round(a.Samples12_0/Layers);
            a.Samples1UI2.String = a.Samples12;
            if  strcmp(Geometry,'Plate')
                if  Symmetric
                    a.Symmetric_Torsional_ModesUI2.Enable = 'on';
                    a.Antisymmetric_Longitudinal_ModesUI2.Enable = 'on';
                else
                    a.Symmetric_Torsional_ModesUI2.Enable = 'off';
                    a.Antisymmetric_Longitudinal_ModesUI2.Enable = 'off';
                end
                a.Symmetric_Torsional_ModesUI2.Tooltip = 'If the layup is symmtric with respect to the middle plane, one can distinguish symmetric and antisymmetric modes. Check this in order to calculate the symmetric modes. These modes have a symmetric displacement pattern with respect to the middle plane of the plate.';
                a.Antisymmetric_Longitudinal_ModesUI2.Tooltip = 'Check this in order to calculate the antisymmetric modes. These modes have an antisymmetric displacement pattern with respect to the middle plane of the plate.';
                a.Lamb_Flexural_ModesUI2.Tooltip = ['Check this in order to calculate the pure Lamb wave modes. This option is available only in the decoupled case. Lamb and shear horizontal waves are decoupled if the specimen consists solely of 0',char(176),' and 90',char(176),' layers, and if the wave propagation is along either direction. In the decoupled case, Lamb waves are termed ''pure'' modes. These Lamb waves show displacement only in the sagittal plane spanned by the propagation direction x1 and by the out-of-plane direction x3. In the coupled case, the modes are still called Lamb waves, but they have a third displacement component perpendicular to the propagation direction (shear horizontal, x2), and no distinction between Lamb and shear horizontal modes is possible. Lamb waves are indicated by solid lines in the dispersion diagram.'];
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.Tooltip = 'The so-called ''pure'' shear horizontal modes can only exist in the decoupled case. Otherwise, they cannot be distinguished from Lamb waves. The pure shear horizontal modes, however, show displacement only perpendicular (x2) to the propagation direction. Shear horizontal waves are indicated by dashed lines in the dispersion diagram.';
                a.SamplesX3ceTextUI2.String = 'Samples x3';
                a.SColorTextUI2.String = 'S';
                a.AColorTextUI2.String = 'A';
                a.BColorTextUI2.String = 'B';
            elseif strcmp(Geometry,'Rod')
                a.Symmetric_Torsional_ModesUI2.Enable = 'on';
                a.Antisymmetric_Longitudinal_ModesUI2.Enable = 'on';
                a.Symmetric_Torsional_ModesUI2.Tooltip = 'Check this in order to calculate the torsional modes. These modes show displacement only in the circumferential direction. Torsional waves are indicated by dashed lines in the dispersion diagram.';
                a.Antisymmetric_Longitudinal_ModesUI2.Tooltip = 'Check this in order to calculate the longitudinal modes. These modes have displacement components only in the axial and radial directions.';
                a.Lamb_Flexural_ModesUI2.Tooltip = 'Check this in order to calculate the flexural modes. These modes have displacement components in all three directions (axial, radial, and circumferential).';
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.Tooltip = 'For flexural modes, an infinite number n = 1,2,3,... of circumferential wavenumbers can fit into the circumference of a pipe or rod. For torsional and longitudinal modes, n = 0. For each family of n, an infinite number p = 1,2,3,... of higher order modes exists, similarly as we know it from the modes in a plate. The naming for modes in pipes and rods is T(0,p), L(0,p), F(n,p). Enter the maximum number of mode families n to be calculated.';
                a.SamplesX3ceTextUI2.String = 'Samples r';
                a.SColorTextUI2.String = 'T,L';
                a.AColorTextUI2.String = 'F';
            elseif strcmp(Geometry,'Pipe')
                a.Symmetric_Torsional_ModesUI2.Enable = 'on';
                a.Antisymmetric_Longitudinal_ModesUI2.Enable = 'on';
                a.Symmetric_Torsional_ModesUI2.Tooltip = 'Check this in order to calculate the torsional modes. These modes show displacement only in the circumferential direction. Torsional waves are indicated by dashed lines in the dispersion diagram.';
                a.Antisymmetric_Longitudinal_ModesUI2.Tooltip = 'Check this in order to calculate the longitudinal modes. These modes have displacement components only in the axial and radial directions.';
                a.Lamb_Flexural_ModesUI2.Tooltip = 'Check this in order to calculate the flexural modes. These modes have displacement components in all three directions (axial, radial, and circumferential).';
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.Tooltip = 'For flexural modes, an infinite number n = 1,2,3,... of circumferential wavenumbers can fit into the circumference of a pipe or rod. For torsional and longitudinal modes, n = 0. For each family of n, an infinite number p = 1,2,3,... of higher order modes exists, similarly as we know it from the modes in a plate. The naming for modes in pipes and rods is T(0,p), L(0,p), F(n,p). Enter the maximum number of mode families n to be calculated.';
                a.SamplesX3ceTextUI2.String = 'Samples r';
                a.SColorTextUI2.String = 'T,L';
                a.AColorTextUI2.String = 'F';
            elseif strcmp(Geometry,'Circumferential')
                a.Symmetric_Torsional_ModesUI2.Enable = 'off';
                a.Antisymmetric_Longitudinal_ModesUI2.Enable = 'off';
                a.Lamb_Flexural_ModesUI2.Tooltip = 'Check this in order to calculate the Lamb-type modes. These modes show displacement only in the sagittal plane spanned by the cirumferential direction and by the radial direction. Lamb-type waves are indicated by solid lines in the dispersion diagram.';
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.Tooltip = 'Check this in order to calculate the shear horizontal-type modes. These modes show displacement only in the axial direction. Shear horizontal-type waves are indicated by dashed lines in the dispersion diagram.';
                a.SamplesX3ceTextUI2.String = 'Samples r';
                a.BColorTextUI2.String = 'C';
            end
            if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Circumferential')
                a.Symmetric_Torsional_ModesTextUI2.String = 'Symmetric modes';
                a.Antisymmetric_Longitudinal_ModesTextUI2.String = 'Antisymmetric modes';
                a.Lamb_Flexural_ModesTextUI2.String = 'Lamb modes';
                a.ShearHorizontalModes_FlexuralModeOrdersTextUI2.String = 'Shear horizontal modes';
                a.Symmetric_Torsional_ModesUI2.Value = a.SymmetricModes2;
                a.Antisymmetric_Longitudinal_ModesUI2.Value = a.AntisymmetricModes2;
                a.Lamb_Flexural_ModesUI2.Value = a.LambModes2;
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.Enable = 'on';
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.Style = 'checkbox';
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.String = '';
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.Value = a.ShearHorizontalModes2;
                a.PlaneUI2.Enable = 'off';
            else
                a.Symmetric_Torsional_ModesTextUI2.String = 'Torsional modes';
                a.Antisymmetric_Longitudinal_ModesTextUI2.String = 'Longitudinal modes';
                a.Lamb_Flexural_ModesTextUI2.String = 'Flexural modes';
                a.ShearHorizontalModes_FlexuralModeOrdersTextUI2.String = 'Flexural mode orders';
                a.Symmetric_Torsional_ModesUI2.Value = a.TorsionalModes2;
                a.Antisymmetric_Longitudinal_ModesUI2.Value = a.LongitudinalModes2;
                a.Lamb_Flexural_ModesUI2.Value = a.FlexuralModes2;
                if  a.FlexuralModes2
                    a.ShearHorizontalModes_FlexuralModeOrdersUI2.Enable = 'on';
                else
                    a.ShearHorizontalModes_FlexuralModeOrdersUI2.Enable = 'off';  
                end
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.Style = 'edit';
                a.ShearHorizontalModes_FlexuralModeOrdersUI2.String = a.FlexuralModeOrders2;
                a.PlaneUI2.Enable = 'on';
            end
            if  strcmp(Geometry,'Circumferential')
                a.Force2DTracingUI2.Enable = 'off';
                a.PhaseVelocityLimit2UI2.Enable = 'off';
                a.AttenuationLimitUI2.Enable = 'off';
                a.SweepsUI2.Enable = 'off';
                a.SweepsSectionsUI2.Enable = 'off';
                a.HigherOrderModesUI2.Enable = 'on';
            else
                if  Viscoelastic || ((strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod')) && FluidLoading) || (strcmp(Geometry,'Pipe') && (ToggleFluid2 || (ToggleFluid1 && Sink)))
                    a.Force2DTracingUI2.Enable = 'off';
                    a.PhaseVelocityLimit2UI2.Enable = 'on';
                    a.AttenuationLimitUI2.Enable = 'on';
                    a.SweepsUI2.Enable = 'on';
                    a.SweepsSectionsUI2.Enable = 'on';
                    a.HigherOrderModesUI2.Enable = 'off';
                else
                    a.Force2DTracingUI2.Enable = 'on';
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
                end
            end
            a.Geometry2 = Geometry;
            a.Radius = Radius;
            a.InnerRadius = InnerRadius;
            a.Sink2 = Sink;
            a.FluidLoading2 = FluidLoading;
            if  strcmp(Geometry,'Plate')
                a.ToggleUpperFluid2 = ToggleFluid1;
                a.ToggleLowerFluid2 = ToggleFluid2;
                a.UpperFluid2 = Fluid1;
                a.LowerFluid2 = Fluid2;
            else
                a.ToggleLowerFluid2 = ToggleFluid1;
                a.ToggleUpperFluid2 = ToggleFluid2;
                a.LowerFluid2 = Fluid1;
                a.UpperFluid2 = Fluid2;
            end
            a.Couplant2 = Couplant;
            a.MaterialType2 = MaterialType;
            a.Material2 = Material;
            a.MaterialNames = MaterialNames;
            a.MaterialClasses = MaterialClasses;
            a.Viscoelastic2 = Viscoelastic;
            a.IsotropicSinglelayer = IsotropicSinglelayer;
            a.IsotropicMultilayerCylinder = IsotropicMultilayerCylinder;
            a.Hybrid = Hybrid;
            a.UniformLayerThickness = UniformLayerThickness;
            a.LayupThickness = LayupThickness;
            a.PlateThickness = PlateThickness;
            a.SuperLayers = SuperLayers;
            a.Pattern = Pattern;
            a.SuperLayerSize = SuperLayerSize;
            a.SymmetricSystem = SymmetricSystem;
            a.Symmetric2 = Symmetric;
            a.Layers2 = Layers;
            a.LayerOrientations = LayerOrientations;
            a.LayerThicknesses = LayerThicknesses;
            a.LayupString2 = LayupString;
            a.LayupStringTex2 = LayupStringTex;
            a.EffectiveLayupString2 = EffectiveLayupString;
            a.UnitCell = UnitCell;
            a.Search2 = 0;
            if  a.Fix2
                source.Tag = '1001';
            else
                source.Tag = '1002';
            end
            a = CallbackModule_Anisotropic(source,a);
            if  strcmp(Geometry,'Circumferential')
                if  Viscoelastic
                    warndlg('The material is viscoelastic. Material damping is ignored for circumferential waves.','Warning');
                end
                if  InnerRadius/Radius(end) >= .98
                    helpdlg(['The circumference closely resembles a flat plate (inner diameter/outer diameter >= 0.98). Consider calculating a ',num2str(LayupThickness),' mm thick plate.'],'Note');
                end
            end
            DisplayProperties
            close(f2)
        elseif strcmp(source.Tag,'17') % Cancel
            close(f2)
        end
    end
    function TableCallback(hObject,callbackdata)
        if  callbackdata.Indices(2) < 7
            if  callbackdata.Indices(2) < 3
                if  isempty(callbackdata.EditData)
                    hObject.Data{callbackdata.Indices(1),callbackdata.Indices(2)} = callbackdata.PreviousData;
                    errordlg('Invalid input!','Error');
                    return
                end
                hObject.Data{callbackdata.Indices(1),callbackdata.Indices(2)} = str2double(replace(callbackdata.EditData,',','.'));
                if  callbackdata.Indices(2) == 1
                    LayerOrientations(callbackdata.Indices(1)) = hObject.Data{callbackdata.Indices(1),1};
                elseif callbackdata.Indices(2) == 2
                    if  hObject.Data{callbackdata.Indices(1),2} <= 0
                        hObject.Data{callbackdata.Indices(1),2} = 1;
                        errordlg('The layer thickness must be greater than zero.','Error');
                    end
                    LayerThicknesses(callbackdata.Indices(1)) = hObject.Data{callbackdata.Indices(1),2};
                end
            end
            z = cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),1));
            hObject.Data(z,1) = {0};
            LayerOrientations(z) = 0;
            if  ~UniformLayerThickness
                z = cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),2));
                hObject.Data(z,2) = {1};
                LayerThicknesses(z) = 1;
            end
            SuperLayerSize = length(LayerOrientations);
            if  Hybrid
                if  callbackdata.Indices(2) == 3
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Orthotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),4:6) = {''};
                elseif callbackdata.Indices(2) == 4
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.TransverselyIsotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),[3 5 6]) = {''};
                elseif callbackdata.Indices(2) == 5
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Cubic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),[3 4 6]) = {''};
                elseif callbackdata.Indices(2) == 6
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Isotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),3:5) = {''};
                end
                z = all(cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),3:6)),2);
                for i = 1:callbackdata.Indices(1)
                    if  z(i)
if  ~strcmp(Geometry,'Plate')
    Material{i} = getfield(a.Materials.Isotropic,a.FieldNames_Isotropic{1});
    hObject.Data{i,6} = Material{i}.Name;
else
                        Material{i} = getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1});
                        hObject.Data{i,3} = Material{i}.Name;
end
                    end
                end
            else
                Material(1:SuperLayerSize) = Material(1);
                hObject.Data(1:SuperLayerSize,MaterialType+2) = {Material{1}.Name};
                if  MaterialType == 1
                    hObject.Data(1:SuperLayerSize,4:6) = {''};
                elseif MaterialType == 2
                    hObject.Data(1:SuperLayerSize,[3 5 6]) = {''};
                elseif MaterialType == 3
                    hObject.Data(1:SuperLayerSize,[3 4 6]) = {''};
                elseif MaterialType == 4
                    hObject.Data(1:SuperLayerSize,3:5) = {''};
                end
            end
        else
            if  SuperLayerSize > 1
                if  SuperLayerSize >= callbackdata.Indices(1)
                    LayerOrientations(callbackdata.Indices(1)) = [];
                    LayerThicknesses(callbackdata.Indices(1)) = [];
                    Material(callbackdata.Indices(1)) = [];
                end
                SuperLayerSize = length(LayerOrientations);
                hObject.Data(callbackdata.Indices(1),:) = [];
                hObject.Data = [hObject.Data;cell(1,7)];
            else
                errordlg('Do not delete the only layer!','Error');
                return
            end
        end
        SetThicknesses
        MakeLayupString
        LayupUI.String = LayupString;
    end
end
function WaveGeneratorUI3_Callback(~,~) % Tab3_signal simulator
    Frequency = a.Frequency3;
    Cycles = a.Cycles3;
    SamplesPerCycle = a.SamplesPerCycle3;
    Chirp = a.Chirp3;
    ChirpType = a.ChirpType3;
    Window = a.Window3;
    WidthFactorGauss = a.WidthFactorGauss3;
    WidthFactorTukey = a.WidthFactorTukey3;
    WidthFactor = a.WidthFactor3;
    if  a.DataType3 == 1
        ALambModes = a.ALambModes1;
        SLambModes = a.SLambModes1;
        BLambModes = a.BLambModes1;
        AShearModes = a.AShearModes1;
        SShearModes = a.SShearModes1;
        BShearModes = a.BShearModes1;
    elseif a.DataType3 == 2
        ALambModes = a.ALambModes2;
        SLambModes = a.SLambModes2;
        BLambModes = a.BLambModes2;
        AShearModes = a.AShearModes2;
        SShearModes = a.SShearModes2;
        BShearModes = a.BShearModes2;
    end

    f5 = figure('Icon',which('DC_Logo16.png'),'NumberTitle','off','Name','Wave generator','Visible','off','MenuBar','none','Position',[0 0 930 260]);
    f5.Units = 'normalized';
    movegui(f5,'center')
    f5.Visible = 'on';

    uicontrol('Parent',f5,'Style','text','HorizontalAlignment','left','String','Chirp/Type','Position',[20 225 83 13]);
    uicontrol('Parent',f5,'Style','text','HorizontalAlignment','left','String','Frequency (kHz)','Position',[20 195 83 13]);
    uicontrol('Parent',f5,'Style','text','HorizontalAlignment','left','String','Cycles','Position',[20 165 111 13]);
    uicontrol('Parent',f5,'Style','text','HorizontalAlignment','left','String','Samples per cycle','Position',[20 135 111 13]);
    uicontrol('Parent',f5,'Style','text','HorizontalAlignment','left','String','Window function','Position',[20 105 100 13]);
    uicontrol('Parent',f5,'Style','text','HorizontalAlignment','left','String','Width factor','Position',[20 75 100 13]);
    uicontrol('Parent',f5,'Style','checkbox','Value',Chirp,'Tooltip','Enable a chirp signal. Due to the dispersive nature of guided waves, the frequency components in a wave packet travel at different speeds, leading to a spreading of the wave packet. In a chirp signal, the frequency increases (up-chirp) or decreases (down-chirp) with time. With a properly designed chirp signal, undesirable dispersive effects in guided waves can be compensated for such that the wave packet maintains its shape while propagating large distances.','Position',[110 220 20 23],'Callback',@Callback,'Tag','1');
    ChirpTypeUI = uicontrol('Parent',f5,'Style','popupmenu','Value',ChirpType,'String',{'Linear','Quadratic','Exponential'},'Tooltip','Select how the frequency changes with time in the chirp signal.','Position',[140 220 120 23],'Callback',@Callback,'Tag','2');
    FrequencyUI = uicontrol('Parent',f5,'Style','edit','String','','Tooltip','Enter the carrier wave frequency in the wave packet. When you have enabled ''Chirp'', you must define the instantaneous frequencies at the start and at the end of the signal. The frequency can increase (up-chirp) or decrease (down-chirp) with time.','Position',[140 190 120 23],'Callback',@Callback,'Tag','3');
    uicontrol('Parent',f5,'Style','edit','String',Cycles,'Tooltip','Enter the number of cycles of the signal within the wave packet.','Position',[140 160 50 23],'Callback',@Callback,'Tag','4');
    uicontrol('Parent',f5,'Style','edit','String',SamplesPerCycle,'Tooltip','Enter the samples per cycle in the excitation signal. This determines the sample rate. In a chirp signal, the number corresponds to the shortest cycle.','Position',[140 130 50 23],'Callback',@Callback,'Tag','5');
    uicontrol('Parent',f5,'Style','popupmenu','Value',Window,'String',{'Gauss','Tukey','Hann','Hamming','Triangular'},'Tooltip','Select the window function to be multiplied with the carrier wave. This sets the shape of the wave packet.','Position',[140 100 120 23],'Callback',@Callback,'Tag','6');
    WidthFactorUI = uicontrol('Parent',f5,'Style','edit','String',WidthFactor,'Tooltip','In case of a Gaussian window, the width is inversely proportional to Width factor. The default value is 3.16. In case of a Tukey window, this value sets the plateau width. A value of 1 returns a Hann window, and 0 gives a rectangular window. The default value is 0.5.','Position',[140 70 50 23],'Callback',@Callback,'Tag','9');
    uicontrol('Parent',f5,'Style','pushbutton','String','OK','Tooltip','Accept the current settings and continue.','Position',[20 20 95 33],'FontSize',10,'Callback',@Callback,'Tag','7');
    uicontrol('Parent',f5,'Style','pushbutton','String','Cancel','Tooltip','Discard the current settings.','Position',[140 20 95 33],'FontSize',10,'Callback',@Callback,'Tag','8');

    if  Chirp
        ChirpTypeUI.Enable = 'on';
        FrequencyUI.String = ['[',num2str(Frequency(1)),' ',num2str(Frequency(2)),']'];
    else
        ChirpTypeUI.Enable = 'off';
        FrequencyUI.String = Frequency;
    end
    if  Window < 3
        WidthFactorUI.Enable = 'on';
    else
        WidthFactorUI.Enable = 'off';
    end
    
    ax1 = axes('Parent',f5,'Units','pixels','Position',[330 40 250 200]); 
    ax1.Box = 'on';
    ax1.FontSize = 9;
    ax1.XLabel.FontSize = 10;
    ax1.YLabel.FontSize = 10;
    ax1.XLabel.String = ['Time (',char(181),'s)'];
    ax1.YAxis.Label.String  = 'Amplitude (a.u.)';
    
    ax2 = axes('Parent',f5,'Units','pixels','Position',[650 40 250 200]);
    ax2.Box = 'on';
    ax2.FontSize = 9;
    ax2.Title.FontSize = 10;
    ax2.XLabel.FontSize = 10;
    ax2.YLabel.FontSize = 10;
    ax2.XLabel.String = 'Frequency (kHz)';
    ax2.YLabel.String = 'Spectral amplitude';
    
    l1 = line(ax1,0,0,'LineWidth',.5,'Color','r');
    l2 = line(ax2,0,0,'LineWidth',.5,'Color','r');

    PlotWave

    function PlotWave
        if  a.DataType3 == 1
            [FrequencyRange,ExcitationMagnitude,ExcitationSignal] = WaveGenerator(Frequency,a.FrequencyResolution1,Cycles,SamplesPerCycle,Chirp,ChirpType,Window,WidthFactor);
        elseif a.DataType3 == 2
            [FrequencyRange,ExcitationMagnitude,ExcitationSignal] = WaveGenerator(Frequency,a.FrequencyResolution2,Cycles,SamplesPerCycle,Chirp,ChirpType,Window,WidthFactor);
        end
        ax1.XLim = [0 ExcitationSignal(1,end)*1e6];
        ax2.XLim = [0 2*mean(Frequency)];
        l1.XData = ExcitationSignal(1,:)*1e6;
        l1.YData = ExcitationSignal(2,:);
        l2.XData = FrequencyRange;
        l2.YData = ExcitationMagnitude(2,:);
    end
    function Callback(source,~)
        source.String = replace(source.String,',','.');
        if  strcmp(source.Tag,'1') % Chirp
            Chirp = source.Value;
            if  Chirp
                Frequency(2) = 2*Frequency(1);
                FrequencyUI.String = ['[',num2str(Frequency(1)),' ',num2str(Frequency(2)),']'];               
                ChirpTypeUI.Enable = 'on';
            else
                Frequency(2) = [];
                FrequencyUI.String = Frequency;
                ChirpTypeUI.Enable = 'off';
            end
        elseif strcmp(source.Tag,'2') % Chirp type
            ChirpType = source.Value;
        elseif strcmp(source.Tag,'3') % Frequency
            Frequency = eval(source.String);
            if  Chirp
                source.String = ['[',num2str(Frequency(1)),' ',num2str(Frequency(2)),']'];
            else
                source.String = Frequency;
            end
        elseif strcmp(source.Tag,'4') % Cycles
            Cycles = str2double(source.String);
        elseif strcmp(source.Tag,'5') % Samples per cycle
            SamplesPerCycle = str2double(source.String);
            if  mod(SamplesPerCycle,2)
                SamplesPerCycle = ceil(SamplesPerCycle)+1;
                source.String = SamplesPerCycle;
                errordlg('Enter only even numbers of samples per cycle!','Warning')
            end
        elseif strcmp(source.Tag,'6') % Window function
            Window = source.Value;
            if  Window == 1
                WidthFactor = WidthFactorGauss;
            elseif Window == 2
                WidthFactor = WidthFactorTukey;
            end
            WidthFactorUI.String = WidthFactor;
            if  Window < 3
                WidthFactorUI.Enable = 'on';
            else
                WidthFactorUI.Enable = 'off';
            end
        elseif strcmp(source.Tag,'9') % Width factor
            WidthFactor = str2double(source.String);
            if  Window == 1
                WidthFactorGauss = WidthFactor;
            elseif Window == 2
                WidthFactorTukey = WidthFactor;
            end
        elseif strcmp(source.Tag,'7') % OK
            if  a.DataType3 == 1
                if  ~isbetween(mean(Frequency),a.FrequencyResolution1,a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1))
                    if  isscalar(Frequency)
                        errordlg(['Selected frequency is out of range! Select between ',num2str(a.FrequencyResolution1),' and ',num2str(a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1)),' kHz.'],'Error');
                    else
                        errordlg(['Selected frequencies are out of range! Mean value must be between ',num2str(a.FrequencyResolution1),' and ',num2str(a.FrequencyLimit1-mod(a.FrequencyLimit1,a.FrequencyResolution1)),' kHz.'],'Error');
                    end
                    return
                end
                [ALambModes,SLambModes,BLambModes,AShearModes,SShearModes,BShearModes] = ModeFinder(a.ALamb1,a.SLamb1,a.BLamb1,a.AShear1,a.SShear1,{0},mean(Frequency));
            elseif a.DataType3 == 2
                if  ~isbetween(mean(Frequency),a.FrequencyResolution2,a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2))
                    if  isscalar(Frequency)
                        errordlg(['Selected frequency is out of range! Select between ',num2str(a.FrequencyResolution2),' and ',num2str(a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2)),' kHz.'],'Error');
                    else
                        errordlg(['Selected frequencies are out of range! Mean value must be between ',num2str(a.FrequencyResolution2),' and ',num2str(a.FrequencyLimit2-mod(a.FrequencyLimit2,a.FrequencyResolution2)),' kHz.'],'Error');
                    end
                    return
                end
                [ALambModes,SLambModes,BLambModes,AShearModes,SShearModes,BShearModes] = ModeFinder(a.ALamb2,a.SLamb2,a.BLamb2,a.AShear2,a.SShear2,a.BShear2,mean(Frequency));
            end
            if  Chirp
                String = [num2str(Frequency(1)),' - ',num2str(Frequency(2)),' kHz, '];
                if  ChirpType == 1
                    String = append(String,'linear chirp');
                elseif ChirpType == 2
                    String = append(String,'quadratic chirp');
                elseif ChirpType == 3
                    String = append(String,'exponential chirp');
                end
                a.FrequencyDisplayUI3.String = String;
            else
                a.FrequencyDisplayUI3.String = [num2str(Frequency),' kHz'];
            end
            a.CyclesDisplayUI3.String = Cycles;
            a.SamplesPerCycleDisplayUI3.String = SamplesPerCycle;
            if  Window == 1
                a.WindowDisplayUI3.String = 'Gauss';
            elseif Window == 2
                a.WindowDisplayUI3.String = 'Tukey';
            elseif Window == 3
                a.WindowDisplayUI3.String = 'Hann';
            elseif Window == 4
                a.WindowDisplayUI3.String = 'Hamming';
            elseif Window == 5
                a.WindowDisplayUI3.String = 'Triangular';
            end
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
            a.Frequency3 = Frequency;
            a.Cycles3 = Cycles;
            a.SamplesPerCycle3 = SamplesPerCycle;
            a.Chirp3 = Chirp;
            a.ChirpType3 = ChirpType;
            a.Window3 = Window;
            a.WidthFactorGauss3 = WidthFactorGauss;
            a.WidthFactorTukey3 = WidthFactorTukey;
            a.WidthFactor3 = WidthFactor;
            if  a.DataType3 == 1
                a.ALambModes1 = ALambModes;
                a.SLambModes1 = SLambModes;
                a.BLambModes1 = BLambModes;
                a.AShearModes1 = AShearModes;
                a.SShearModes1 = SShearModes;
                a.BShearModes1 = BShearModes;
            elseif a.DataType3 == 2
                a.ALambModes2 = ALambModes;
                a.SLambModes2 = SLambModes;
                a.BLambModes2 = BLambModes;
                a.AShearModes2 = AShearModes;
                a.SShearModes2 = SShearModes;
                a.BShearModes2 = BShearModes;
            end
            [a.h2,a.h3,a.h4,a.h5] = Signal_Internal(a,Tab3,0);
            Colorization(a,2)
            ShowModes_Signal(a)
            close(f5)
        elseif strcmp(source.Tag,'8') % Cancel
            close(f5)
        end
        if  isvalid(source) && str2double(source.Tag) ~= 8
            PlotWave
        end
    end
end
function SpecimenSettingsUI4_Callback(~,~) % Tab4_polar diagrams
    MaterialType = a.MaterialType_Polar;
    Material = a.Material_Polar;
    MaterialNames = a.MaterialNames_Polar;
    MaterialClasses = a.MaterialClasses_Polar;
    Hybrid = a.Hybrid_Polar;
    UniformLayerThickness = a.UniformLayerThickness_Polar;
    LayupThickness = a.LayupThickness_Polar;
    PlateThickness = a.PlateThickness_Polar;
    SuperLayers = a.SuperLayers_Polar;
    Pattern = a.Pattern_Polar;
    SuperLayerSize = a.SuperLayerSize_Polar;
    SymmetricSystem = a.SymmetricSystem_Polar;
    Layers = a.Layers_Polar;
    LayerOrientations = a.LayerOrientations_Polar;
    LayerThicknesses = a.LayerThicknesses_Polar;
    LayupString = a.LayupString_Polar;
    LayupStringTex = a.LayupStringTex_Polar;
    UnitCell = a.UnitCell_Polar;

    f3 = figure('Icon',which('DC_Logo16.png'),'NumberTitle','off','Name','Specimen','Visible','off','MenuBar','none','Position',[0 0 1010 370]);
    f3.Units = 'normalized';
    movegui(f3,'center')
    f3.Visible = 'on';

    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Layup','Position',[20 285 54 13])
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Hybrid','Position',[20 255 33 13])
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Class','Position',[20 225 29 13])
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Material','Position',[20 195 39 13]);
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Uniform layer thickness','Position',[20 165 115 13]);
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Layup thickness (mm)','Position',[20 135 101 13]);
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Unit cell repetitions','Position',[20 105 92 13]);
    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Symmetric unit cell','Position',[20 75 90 13]);

    LayupUI = uicontrol('Parent',f3,'Style','text','String',LayupString,'Position',[90 286 170 13],'BackgroundColor','white');
    HybridUI = uicontrol('Parent',f3,'Style','checkbox','Value',Hybrid,'Tooltip','The layup may contain different materials.','Position',[90 250 50 23],'Callback',@Callback,'Tag','8');
    MaterialClassUI = uicontrol('Parent',f3,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic','Isotropic'},'Value',MaterialType,'Tooltip','Select a material symmetry class.','Position',[90 220 170 23],'Callback',@Callback,'Tag','9');
    MaterialUI = uicontrol('Parent',f3,'Style','popupmenu','String',{''},'Tooltip','Select a material.','Position',[90 190 170 23],'Callback',@Callback,'Tag','10');
    UniformLayerThicknessUI = uicontrol('Parent',f3,'Style','checkbox','Value',UniformLayerThickness,'Tooltip','Check this if the layup has uniform layer thicknesses. Dispersion Calculator deduces the layer thicknesses from the total layup thickness as defined below, and you do not need to enter the individual layer thicknesses into the table to the right.','Position',[150 160 50 23],'Callback',@Callback,'Tag','11');
    LayupThicknessUI = uicontrol('Parent',f3,'Style','edit','String',LayupThickness,'Tooltip','Enter the total layup thickness. This is available only if you have checked ''Uniform layer thickness''.','Position',[150 130 50 23],'Callback',@Callback,'Tag','12');
    SuperLayersUI = uicontrol('Parent',f3,'Style','edit','String',SuperLayers,'Tooltip','Enter the number of repetitions of the unit cell contained in the layup. Use only integers.','Position',[150 100 50 23],'Callback',@Callback,'Tag','13');
    SymmetricSystemUI = uicontrol('Parent',f3,'Style','checkbox','Value',SymmetricSystem,'Tooltip','Check this if the layup is symmetric about its middle plane. Dispersion Calculator extends the layup accordingly.','Position',[150 70 50 23],'Callback',@Callback,'Tag','14');    

    uicontrol('Parent',f3,'Style','text','HorizontalAlignment','left','String','Unit cell','Position',[280 310 46 13],'FontSize',9);
    TableUI = uitable(f3,'ColumnName',{['Phi (',char(176),')'],'d (mm)','Orthotropic','Trans. iso.','Cubic','Isotropic','Delete'},'ColumnWidth',{60 60 120 120 120 120 55},'ColumnFormat',({[] [] a.FieldNames_Orthotropic' a.FieldNames_TransverselyIsotropic' a.FieldNames_Cubic' a.FieldNames_Isotropic' 'logical'}),'Tooltip','Enter the fiber orientation of every layer into the left column, and if you don''t have checked ''Uniform layer thickness'', enter also the corresponding layer thicknesses into the second column.','Data',UnitCell,'RowStriping','off','Position',[280 20 715 280],'CellEditCallback',@TableCallback);

    uicontrol('Parent',f3,'Style','pushbutton','String','Reset','Tooltip','Reset specimen definition.','Position',[280 330 95 33],'FontSize',10,'Callback',@Callback,'Tag','15'); 
    uicontrol('Parent',f3,'Style','pushbutton','String','OK','Tooltip','Accept the current specimen definition and continue.','Position',[20 20 95 33],'FontSize',10,'Callback',@Callback,'Tag','16');
    uicontrol('Parent',f3,'Style','pushbutton','String','Cancel','Tooltip','Discard the current specimen definition.','Position',[150 20 95 33],'FontSize',10,'Callback',@Callback,'Tag','17');

    TableColumnEditable

    if  Hybrid
        MaterialClassUI.Enable = 'off';
        MaterialUI.Enable = 'off';
        MaterialUI.Value = 1;
    else
        MaterialClassUI.Enable = 'on';
        MaterialUI.Enable = 'on';
        if  MaterialType == 1
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Orthotropic));
        elseif MaterialType == 2
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_TransverselyIsotropic));
        elseif MaterialType == 3
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Cubic));
        elseif MaterialType == 4
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Isotropic));
        end
    end
    if  MaterialType == 1
        MaterialUI.String = a.FieldNames_Orthotropic;
    elseif MaterialType == 2
        MaterialUI.String = a.FieldNames_TransverselyIsotropic;
    elseif MaterialType == 3
        MaterialUI.String = a.FieldNames_Cubic;
    elseif MaterialType == 4
        MaterialUI.String = a.FieldNames_Isotropic;
    end
    if  UniformLayerThickness
        LayupThicknessUI.Enable = 'on';
    else
        LayupThicknessUI.Enable = 'off';
    end

    function TableColumnEditable
        if  Hybrid && UniformLayerThickness
            TableUI.ColumnEditable = [true false true true true true true];
        elseif Hybrid && ~UniformLayerThickness
            TableUI.ColumnEditable = [true true true true true true true];
        elseif ~Hybrid && UniformLayerThickness
            TableUI.ColumnEditable = [true false false false false false true];
        elseif ~Hybrid && ~UniformLayerThickness
            TableUI.ColumnEditable = [true true false false false false true];
        end
    end
    function DistributeMaterials
        if  MaterialType == 1
            Material(1:SuperLayerSize) = {getfield(a.Materials.Orthotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 2
            Material(1:SuperLayerSize) = {getfield(a.Materials.TransverselyIsotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 3
            Material(1:SuperLayerSize) = {getfield(a.Materials.Cubic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 4
            Material(1:SuperLayerSize) = {getfield(a.Materials.Isotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        end
    end
    function ClearTableColumns
        if  MaterialType == 1
            TableUI.Data(1:SuperLayerSize,4:6) = {''};
        elseif MaterialType == 2
            TableUI.Data(1:SuperLayerSize,[3 5 6]) = {''};
        elseif MaterialType == 3
            TableUI.Data(1:SuperLayerSize,[3 4 6]) = {''};
        elseif MaterialType == 4
            TableUI.Data(1:SuperLayerSize,3:5) = {''};
        end
        TableUI.Data(1:SuperLayerSize,MaterialType+2) = {Material{1}.Name};
    end
    function SetThicknesses
        if  UniformLayerThickness
            if  SymmetricSystem
                LayerThicknesses(1:SuperLayerSize) = LayupThickness/SuperLayerSize/SuperLayers/2;
            else
                LayerThicknesses(1:SuperLayerSize) = LayupThickness/SuperLayerSize/SuperLayers;
            end
            TableUI.Data(1:SuperLayerSize,2) = {LayerThicknesses(1)};
        else
            LayupThickness = SuperLayers*sum(LayerThicknesses);
            if  SymmetricSystem
                LayupThickness = LayupThickness*2;
            end
            LayupThicknessUI.String = LayupThickness;
        end
    end
    function MakeLayupString
        LayupString = ['[',replace(num2str(LayerOrientations),whitespacePattern,'/'),']'];
        LayupStringTex = [' ',LayupString];
        if  SuperLayers > 1
            LayupString = [LayupString,num2str(SuperLayers)];
            LayupStringTex = [LayupStringTex,'$_{',num2str(SuperLayers),'}$'];
        end
        if  SymmetricSystem
            LayupString = [LayupString,'s'];
            LayupStringTex = [LayupStringTex,'$_{\mathrm s}$'];
            LayupStringTex = erase(LayupStringTex,'}$$_{');
        end
    end
    function DisplayProperties
        String = ['LayerOrientations:   ',regexprep(num2str(LayerOrientations),'\s+',' '),newline...
                  'LayerThicknesses:    ',regexprep(num2str(LayerThicknesses),'\s+',' '),newline...
                  'Layers:              ',num2str(SuperLayerSize),newline...
                  'Total layers:        ',num2str(Layers),newline...
                  'LayupThickness:      ',num2str(LayupThickness),newline...
                  'Waveguide thickness: ',num2str(PlateThickness),newline];
        for i = 1:length(Material)
            String = append(String,num2str(i),' ',Material{i}.Name,' ',Material{i}.Class,newline);
        end
        disp([String,'---------------------------'])
    end
    function Callback(source,~)
        source.String = replace(source.String,',','.');
        if  strcmp(source.Tag,'8') % Hybrid
            Hybrid = source.Value;
            if  Hybrid
                MaterialClassUI.Enable = 'off';
                MaterialUI.Enable = 'off';
            else
                MaterialClassUI.Enable = 'on';
                MaterialUI.Enable = 'on';
                DistributeMaterials
                ClearTableColumns
            end
            TableColumnEditable
        elseif strcmp(source.Tag,'9') % Class
            MaterialType = source.Value;
            MaterialUI.Value = 1;
            if  MaterialType == 1
                MaterialUI.String = a.FieldNames_Orthotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})};
            elseif MaterialType == 2
                MaterialUI.String = a.FieldNames_TransverselyIsotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.TransverselyIsotropic,a.FieldNames_TransverselyIsotropic{1})};
            elseif MaterialType == 3
                MaterialUI.String = a.FieldNames_Cubic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Cubic,a.FieldNames_Cubic{1})};
            elseif MaterialType == 4
                MaterialUI.String = a.FieldNames_Isotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Isotropic,a.FieldNames_Isotropic{1})};
            end
            ClearTableColumns
        elseif strcmp(source.Tag,'10') % Material
            DistributeMaterials
            ClearTableColumns
        elseif strcmp(source.Tag,'11') % Uniform layer thickness
            UniformLayerThickness = source.Value;
            if  UniformLayerThickness
                LayupThicknessUI.Enable = 'on';
                SetThicknesses
            else
                LayupThicknessUI.Enable = 'off';
            end
            TableColumnEditable
        elseif strcmp(source.Tag,'12') % Layup thickness
            if  str2double(source.String) <= 0
                source.String = 1;
                errordlg('The layup thickness must be greater than zero.','Error');
            end
            LayupThickness = str2double(source.String);
            SetThicknesses
        elseif strcmp(source.Tag,'13') % Unit cell repetitions
            if  rem(str2double(source.String),1) ~= 0 || str2double(source.String) < 1
                source.String = ceil(abs(str2double(source.String)));
                errordlg('Unit cell repetitions must be a positive integer.','Error');
            end
            SuperLayers = str2double(source.String);
            SetThicknesses
            MakeLayupString
            LayupUI.String = LayupString;
        elseif strcmp(source.Tag,'14') % Symmetric unit cell
            SymmetricSystem = source.Value;
            SetThicknesses
            MakeLayupString
            LayupUI.String = LayupString;
        elseif strcmp(source.Tag,'15') % Reset
            MaterialType = 1;
            Material = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})};
            MaterialNames = {Material{1}.Name};
            MaterialClasses = {Material{1}.Class};
            Hybrid = 0;
            UniformLayerThickness = 1;
            LayupThickness = 1;
            PlateThickness = LayupThickness;
            SuperLayers = 1;
            Pattern = [];
            SuperLayerSize = 1;
            SymmetricSystem = 0;
            Layers = 1;
            LayerOrientations = 0;
            LayerThicknesses = 1;
            LayupString = '[0]';
            LayupStringTex = '[0]';
            UnitCell = cell(99,7);
            UnitCell{1} = 0;
            UnitCell{1,2} = 1;
            UnitCell{1,3} = a.FieldNames_Orthotropic{1};
            LayupUI.String = LayupString;
            HybridUI.Value = Hybrid;
            MaterialClassUI.Value = MaterialType;
            MaterialClassUI.Enable = 'on';
            MaterialUI.String = a.FieldNames_Orthotropic;
            MaterialUI.Value = 1;
            MaterialUI.Enable = 'on';
            UniformLayerThicknessUI.Value = UniformLayerThickness;
            LayupThicknessUI.String = LayupThickness;
            LayupThicknessUI.Enable = 'on';
            SuperLayersUI.String = SuperLayers;
            SymmetricSystemUI.Value = SymmetricSystem;
            TableUI.ColumnEditable = [true false false false false false true];
            TableUI.Data = UnitCell;
        elseif strcmp(source.Tag,'16') % OK
            MaterialNames = {};
            MaterialClasses = {};
            for i = 1:SuperLayerSize
                MaterialNames{i} = Material{i}.Name;
                MaterialClasses{i} = Material{i}.Class;
            end
            if  SuperLayerSize == 1 && SuperLayers > 1
                errordlg('Do not define multiple repetitions if the layup consists of only one layer! Increase the layer thickness instead.','Error');
                SuperLayers = 1;
                SuperLayersUI.String = 1;
                LayerThicknesses = LayupThickness;
                TableUI.Data{1,2} = LayerThicknesses(1);
                MakeLayupString
                LayupUI.String = LayupString;
                return
            elseif SuperLayerSize == 1 && SymmetricSystem
                errordlg('Do not define a symmetric layup if the layup consists of only one layer! Set double layer thickness instead.','Error');
                SymmetricSystem = 0;
                SymmetricSystemUI.Value = 0;
                LayerThicknesses = LayupThickness;
                TableUI.Data{1,2} = LayerThicknesses(1);
                MakeLayupString
                LayupUI.String = LayupString;
                return
            elseif SuperLayerSize > 1 && all(strcmp(MaterialClasses,'Isotropic')) && all(strcmp(MaterialNames(1),MaterialNames))
                errordlg('Do not define a unit cell containing solely layers of the same isotropic material! Replace by one layer of greater thickness. This will be the same because DC assumes rigid bonding between the layers. It is also more efficient since the computational expense scales with the number of layers.','Error');
                return
            elseif SuperLayerSize > 1 && all(LayerOrientations(1) == LayerOrientations) && all(strcmp(MaterialNames(1),MaterialNames))
                errordlg('Do not define a unit cell containing solely layers of the same material and orientation! Replace by one layer of greater thickness. This will be the same because DC assumes rigid bonding between the layers. It is also more efficient since the computational expense scales with the number of layers.','Error');
                return
            elseif ~SymmetricSystem && SuperLayerSize > 1 &&...
                (~mod(SuperLayerSize,2) && all(fliplr(LayerOrientations(1:SuperLayerSize/2)) == LayerOrientations(SuperLayerSize/2+1:end)) && all(fliplr(LayerThicknesses(1:SuperLayerSize/2)) == LayerThicknesses(SuperLayerSize/2+1:end)) && all(strcmp(fliplr(MaterialNames(1:SuperLayerSize/2)),MaterialNames(SuperLayerSize/2+1:end))) ||...
                mod(SuperLayerSize,2) && all(fliplr(LayerOrientations(1:SuperLayerSize/2-.5)) == LayerOrientations(SuperLayerSize/2+1.5:end)) && all(fliplr(LayerThicknesses(1:SuperLayerSize/2-.5)) == LayerThicknesses(SuperLayerSize/2+1.5:end)) && all(strcmp(fliplr(MaterialNames(1:SuperLayerSize/2-.5)),MaterialNames(SuperLayerSize/2+1.5:end))))
                errordlg('Do not define a symmetric unit cell! Cut the unit cell in half and check ''Symmetric unit cell'' instead. This will be more robust and more efficient.','Error');
                return
            elseif all(strcmp(MaterialNames(1),MaterialNames)) && Hybrid
                if  SuperLayerSize > 1
                    errordlg('It is no hybrid if all layers are the same material!','Error');
                else
                    errordlg('It is no hybrid if you have only one layer!','Error');
                end
                Hybrid = 0;
                HybridUI.Value = 0;
                MaterialClassUI.Enable = 'on';
                MaterialUI.Enable = 'on';
                DistributeMaterials
                ClearTableColumns
                TableColumnEditable
                return
            elseif all(LayerOrientations)
                errordlg('The layup must contain at least one 0 deg layer!','Error');
                return
            end
            UnitCell = TableUI.Data;
            for i = 1:SuperLayerSize
                if  ~isreal(Material{i}.C)
                    warndlg('A viscoelastic material is contained in the layup. The damping is ignored, i.e., only the real parts of the stiffness components are taken.','Warning');
                    break
                end
            end
            Layers = SuperLayerSize*SuperLayers;
            if  SymmetricSystem
                Layers = Layers*2;
            end
            PlateThickness = LayupThickness;
            n = 2.^(0:log2(SuperLayers));
            Pattern = [0:length(n)-1;n];
            while Pattern(2,end) < SuperLayers
                for i = length(n):-1:1
                    if  Pattern(2,end)+n(i) <= SuperLayers
                        Pattern(1:2,end+1) = [i;Pattern(2,end)+n(i)];
                        break
                    end
                end
            end
            Pattern = Pattern(1,2:end);
            if  Hybrid
                a.MaterialNameUI4.String = 'Hybrid';
            else
                a.MaterialNameUI4.String = Material{1}.Name;
            end
            MakeLayupString
            a.LayupUI4.String = LayupString;
            a.LayerCountUI4.String = Layers;
            a.ThicknessCountUI4.String = PlateThickness;
            if  a.Quantity_Polar < 3
                if  Layers == 1
                    a.Option1UI4.Enable = 'on';
                else
                    a.Option1UI4.Enable = 'off';
                end
            end
            a.FrequencyLimit_Polar = round(1e3*a.XRange_Polar/PlateThickness,2,'significant');
            a.FrequencyLimitUI4.String = a.FrequencyLimit_Polar;
            a.FrequencyResolution_Polar = a.FrequencyLimit_Polar/a.XSamples_Polar;
            a.FrequencyResolutionUI4.String = a.FrequencyResolution_Polar;
            a.Frequency_Polar = a.FrequencyLimit_Polar;
            a.FrequencyUI4.String = a.FrequencyLimit_Polar;
            a.MaterialType_Polar = MaterialType;
            a.Material_Polar = Material;
            a.MaterialNames_Polar = MaterialNames;
            a.MaterialClasses_Polar = MaterialClasses;
            a.Hybrid_Polar = Hybrid;
            a.UniformLayerThickness_Polar = UniformLayerThickness;
            a.LayupThickness_Polar = LayupThickness;
            a.PlateThickness_Polar = PlateThickness;
            a.SuperLayers_Polar = SuperLayers;
            a.Pattern_Polar = Pattern;
            a.SuperLayerSize_Polar = SuperLayerSize;
            a.SymmetricSystem_Polar = SymmetricSystem;
            a.Layers_Polar = Layers;
            a.LayerOrientations_Polar = LayerOrientations;
            a.LayerThicknesses_Polar = LayerThicknesses;
            a.LayupString_Polar = LayupString;
            a.LayupStringTex_Polar = LayupStringTex;
            a.UnitCell_Polar = UnitCell;
            DisplayProperties
            close(f3)
        elseif strcmp(source.Tag,'17') % Cancel
            close(f3)
        end
    end
    function TableCallback(hObject,callbackdata)
        if  callbackdata.Indices(2) < 7
            if  callbackdata.Indices(2) < 3
                if  isempty(callbackdata.EditData)
                    hObject.Data{callbackdata.Indices(1),callbackdata.Indices(2)} = callbackdata.PreviousData;
                    errordlg('Invalid input!','Error');
                    return
                end
                hObject.Data{callbackdata.Indices(1),callbackdata.Indices(2)} = str2double(replace(callbackdata.EditData,',','.'));
                if  callbackdata.Indices(2) == 1
                    LayerOrientations(callbackdata.Indices(1)) = hObject.Data{callbackdata.Indices(1),1};
                elseif callbackdata.Indices(2) == 2
                    if  hObject.Data{callbackdata.Indices(1),2} <= 0
                        hObject.Data{callbackdata.Indices(1),2} = 1;
                        errordlg('The layer thickness must be greater than zero.','Error');
                    end
                    LayerThicknesses(callbackdata.Indices(1)) = hObject.Data{callbackdata.Indices(1),2};
                end
            end
            z = cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),1));
            hObject.Data(z,1) = {0};
            LayerOrientations(z) = 0;
            if  ~UniformLayerThickness
                z = cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),2));
                hObject.Data(z,2) = {1};
                LayerThicknesses(z) = 1;
            end
            SuperLayerSize = length(LayerOrientations);
            if  Hybrid
                if  callbackdata.Indices(2) == 3
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Orthotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),4:6) = {''};
                elseif callbackdata.Indices(2) == 4
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.TransverselyIsotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),[3 5 6]) = {''};
                elseif callbackdata.Indices(2) == 5
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Cubic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),[3 4 6]) = {''};
                elseif callbackdata.Indices(2) == 6
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Isotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),3:5) = {''};
                end
                z = all(cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),3:6)),2);
                for i = 1:callbackdata.Indices(1)
                    if  z(i)
                        Material{i} = getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1});
                        hObject.Data{i,3} = Material{i}.Name;
                    end
                end
            else
                Material(1:SuperLayerSize) = Material(1);
                hObject.Data(1:SuperLayerSize,MaterialType+2) = {Material{1}.Name};
                if  MaterialType == 1
                    hObject.Data(1:SuperLayerSize,4:6) = {''};
                elseif MaterialType == 2
                    hObject.Data(1:SuperLayerSize,[3 5 6]) = {''};
                elseif MaterialType == 3
                    hObject.Data(1:SuperLayerSize,[3 4 6]) = {''};
                elseif MaterialType == 4
                    hObject.Data(1:SuperLayerSize,3:5) = {''};
                end
            end
        else
            if  SuperLayerSize > 1
                if  SuperLayerSize >= callbackdata.Indices(1)
                    LayerOrientations(callbackdata.Indices(1)) = [];
                    LayerThicknesses(callbackdata.Indices(1)) = [];
                    Material(callbackdata.Indices(1)) = [];
                end
                SuperLayerSize = length(LayerOrientations);
                hObject.Data(callbackdata.Indices(1),:) = [];
                hObject.Data = [hObject.Data;cell(1,7)];
            else
                errordlg('Do not delete the only layer!','Error');
                return
            end
        end
        SetThicknesses
        MakeLayupString
        LayupUI.String = LayupString;
    end
end
function SpecimenSettingsUI6_Callback(~,~) % Tab6_laminate stiffness
    MaterialType = a.MaterialType6;
    Material = a.Material6;
    MaterialNames = a.MaterialNames6;
    MaterialClasses = a.MaterialClasses6;
    Hybrid = a.Hybrid6;
    UniformLayerThickness = a.UniformLayerThickness6;
    LayupThickness = a.LayupThickness6;
    SuperLayerSize = a.SuperLayerSize6;
    LayerOrientations = a.LayerOrientations6;
    LayerThicknesses = a.LayerThicknesses6;
    LayupString = a.LayupString6;
    UnitCell = a.UnitCell6;

    f4 = figure('Icon',which('DC_Logo16.png'),'NumberTitle','off','Name','Specimen','Visible','off','MenuBar','none','Position',[0 0 1010 370]);
    f4.Units = 'normalized';
    movegui(f4,'center')
    f4.Visible = 'on';

    uicontrol('Parent',f4,'Style','text','HorizontalAlignment','left','String','Layup','Position',[20 285 54 13])
    uicontrol('Parent',f4,'Style','text','HorizontalAlignment','left','String','Hybrid','Position',[20 255 33 13])
    uicontrol('Parent',f4,'Style','text','HorizontalAlignment','left','String','Class','Position',[20 225 29 13])
    uicontrol('Parent',f4,'Style','text','HorizontalAlignment','left','String','Material','Position',[20 195 39 13]);
    uicontrol('Parent',f4,'Style','text','HorizontalAlignment','left','String','Uniform layer thickness','Position',[20 165 115 13]);

    LayupUI = uicontrol('Parent',f4,'Style','text','String',LayupString,'Position',[90 286 170 13],'BackgroundColor','white');
    HybridUI = uicontrol('Parent',f4,'Style','checkbox','Value',Hybrid,'Tooltip','The layup may contain different materials.','Position',[90 250 50 23],'Callback',@Callback,'Tag','8');
    MaterialClassUI = uicontrol('Parent',f4,'Style','popupmenu','String',{'Orthotropic','Transversely isotropic','Cubic','Isotropic'},'Value',MaterialType,'Tooltip','Select a material symmetry class.','Position',[90 220 170 23],'Callback',@Callback,'Tag','9');
    MaterialUI = uicontrol('Parent',f4,'Style','popupmenu','String',{''},'Tooltip','Select a material.','Position',[90 190 170 23],'Callback',@Callback,'Tag','10');
    UniformLayerThicknessUI = uicontrol('Parent',f4,'Style','checkbox','Value',UniformLayerThickness,'Tooltip','Check this if the layup has uniform layer thicknesses. Dispersion Calculator deduces the layer thicknesses from the total layup thickness as defined below, and you do not need to enter the individual layer thicknesses into the table to the right.','Position',[150 160 50 23],'Callback',@Callback,'Tag','11');

    uicontrol('Parent',f4,'Style','text','HorizontalAlignment','left','String','Unit cell','Position',[280 310 46 13],'FontSize',9);
    TableUI = uitable(f4,'ColumnName',{['Phi (',char(176),')'],'d (mm)','Orthotropic','Trans. iso.','Cubic','Isotropic','Delete'},'ColumnWidth',{60 60 120 120 120 120 55},'ColumnFormat',({[] [] a.FieldNames_Orthotropic' a.FieldNames_TransverselyIsotropic' a.FieldNames_Cubic' a.FieldNames_Isotropic' 'logical'}),'Tooltip','Enter the fiber orientation of every layer into the left column, and if you don''t have checked ''Uniform layer thickness'', enter also the corresponding layer thicknesses into the second column.','Data',UnitCell,'RowStriping','off','Position',[280 20 715 280],'CellEditCallback',@TableCallback);

    uicontrol('Parent',f4,'Style','pushbutton','String','Reset','Tooltip','Reset specimen definition.','Position',[280 330 95 33],'FontSize',10,'Callback',@Callback,'Tag','15'); 
    uicontrol('Parent',f4,'Style','pushbutton','String','OK','Tooltip','Accept the current specimen definition and continue.','Position',[20 20 95 33],'FontSize',10,'Callback',@Callback,'Tag','16');
    uicontrol('Parent',f4,'Style','pushbutton','String','Cancel','Tooltip','Discard the current specimen definition.','Position',[150 20 95 33],'FontSize',10,'Callback',@Callback,'Tag','17');

    TableColumnEditable

    if  Hybrid
        MaterialClassUI.Enable = 'off';
        MaterialUI.Enable = 'off';
        MaterialUI.Value = 1;
    else
        MaterialClassUI.Enable = 'on';
        MaterialUI.Enable = 'on';
        if  MaterialType == 1
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Orthotropic));
        elseif MaterialType == 2
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_TransverselyIsotropic));
        elseif MaterialType == 3
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Cubic));
        elseif MaterialType == 4
            MaterialUI.Value = find(strcmp(Material{1}.Name,a.FieldNames_Isotropic));
        end
    end
    if  MaterialType == 1
        MaterialUI.String = a.FieldNames_Orthotropic;
    elseif MaterialType == 2
        MaterialUI.String = a.FieldNames_TransverselyIsotropic;
    elseif MaterialType == 3
        MaterialUI.String = a.FieldNames_Cubic;
    elseif MaterialType == 4
        MaterialUI.String = a.FieldNames_Isotropic;
    end

    function TableColumnEditable
        if  Hybrid && UniformLayerThickness
            TableUI.ColumnEditable = [true false true true true true true];
        elseif Hybrid && ~UniformLayerThickness
            TableUI.ColumnEditable = [true true true true true true true];
        elseif ~Hybrid && UniformLayerThickness
            TableUI.ColumnEditable = [true false false false false false true];
        elseif ~Hybrid && ~UniformLayerThickness
            TableUI.ColumnEditable = [true true false false false false true];
        end
    end
    function DistributeMaterials
        if  MaterialType == 1
            Material(1:SuperLayerSize) = {getfield(a.Materials.Orthotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 2
            Material(1:SuperLayerSize) = {getfield(a.Materials.TransverselyIsotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 3
            Material(1:SuperLayerSize) = {getfield(a.Materials.Cubic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        elseif MaterialType == 4
            Material(1:SuperLayerSize) = {getfield(a.Materials.Isotropic,cell2mat(MaterialUI.String(MaterialUI.Value)))};
        end
    end
    function ClearTableColumns
        if  MaterialType == 1
            TableUI.Data(1:SuperLayerSize,4:6) = {''};
        elseif MaterialType == 2
            TableUI.Data(1:SuperLayerSize,[3 5 6]) = {''};
        elseif MaterialType == 3
            TableUI.Data(1:SuperLayerSize,[3 4 6]) = {''};
        elseif MaterialType == 4
            TableUI.Data(1:SuperLayerSize,3:5) = {''};
        end
        TableUI.Data(1:SuperLayerSize,MaterialType+2) = {Material{1}.Name};
    end
    function SetThicknesses
        if  UniformLayerThickness
            LayerThicknesses(1:SuperLayerSize) = LayupThickness/SuperLayerSize;
            TableUI.Data(1:SuperLayerSize,2) = {LayerThicknesses(1)};
        else
            LayupThickness = sum(LayerThicknesses);
        end
    end
    function MakeLayupString
        LayupString = ['[',replace(num2str(LayerOrientations),whitespacePattern,'/'),']'];
    end
    function DisplayProperties
        String = ['LayerOrientations:   ',regexprep(num2str(LayerOrientations),'\s+',' '),newline...
                  'LayerThicknesses:    ',regexprep(num2str(LayerThicknesses),'\s+',' '),newline...
                  'Layers:              ',num2str(SuperLayerSize),newline...
                  'LayupThickness:      ',num2str(LayupThickness),newline];
        for i = 1:length(Material)
            String = append(String,num2str(i),' ',Material{i}.Name,' ',Material{i}.Class,newline);
        end
        disp([String,'---------------------------'])
    end
    function Callback(source,~)
        source.String = replace(source.String,',','.');
        if  strcmp(source.Tag,'8') % Hybrid
            Hybrid = source.Value;
            if  Hybrid
                MaterialClassUI.Enable = 'off';
                MaterialUI.Enable = 'off';
            else
                MaterialClassUI.Enable = 'on';
                MaterialUI.Enable = 'on';
                DistributeMaterials
                ClearTableColumns
            end
            TableColumnEditable
        elseif strcmp(source.Tag,'9') % Class
            MaterialType = source.Value;
            MaterialUI.Value = 1;
            if  MaterialType == 1
                MaterialUI.String = a.FieldNames_Orthotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})};
            elseif MaterialType == 2
                MaterialUI.String = a.FieldNames_TransverselyIsotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.TransverselyIsotropic,a.FieldNames_TransverselyIsotropic{1})};
            elseif MaterialType == 3
                MaterialUI.String = a.FieldNames_Cubic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Cubic,a.FieldNames_Cubic{1})};
            elseif MaterialType == 4
                MaterialUI.String = a.FieldNames_Isotropic;
                Material(1:SuperLayerSize) = {getfield(a.Materials.Isotropic,a.FieldNames_Isotropic{1})};
            end
            ClearTableColumns
        elseif strcmp(source.Tag,'10') % Material
            DistributeMaterials
            ClearTableColumns
        elseif strcmp(source.Tag,'11') % Uniform layer thickness
            UniformLayerThickness = source.Value;
            if  UniformLayerThickness
                SetThicknesses
            end
            TableColumnEditable
        elseif strcmp(source.Tag,'15') % Reset
            MaterialType = 1;
            Material = {getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1})};
            MaterialNames = {Material{1}.Name};
            MaterialClasses = {Material{1}.Class};
            Hybrid = 0;
            UniformLayerThickness = 1;
            LayupThickness = 1;
            SuperLayerSize = 1;
            LayerOrientations = 0;
            LayerThicknesses = 1;
            LayupString = '[0]';
            UnitCell = cell(99,7);
            UnitCell{1} = 0;
            UnitCell{1,2} = 1;
            UnitCell{1,3} = a.FieldNames_Orthotropic{1};
            LayupUI.String = LayupString;
            HybridUI.Value = Hybrid;
            MaterialClassUI.Value = MaterialType;
            MaterialClassUI.Enable = 'on';
            MaterialUI.String = a.FieldNames_Orthotropic;
            MaterialUI.Value = 1;
            MaterialUI.Enable = 'on';
            UniformLayerThicknessUI.Value = UniformLayerThickness;
            TableUI.ColumnEditable = [true false false false false false true];
            TableUI.Data = UnitCell;
        elseif strcmp(source.Tag,'16') % OK
            MaterialNames = {};
            MaterialClasses = {};
            for i = 1:SuperLayerSize
                MaterialNames{i} = Material{i}.Name;
                MaterialClasses{i} = Material{i}.Class;
            end
            if  all(strcmp(MaterialNames(1),MaterialNames)) && Hybrid
                if  SuperLayerSize > 1
                    errordlg('It is no hybrid if all layers are the same material!','Error');
                else
                    errordlg('It is no hybrid if you have only one layer!','Error');
                end
                Hybrid = 0;
                HybridUI.Value = 0;
                MaterialClassUI.Enable = 'on';
                MaterialUI.Enable = 'on';
                DistributeMaterials
                ClearTableColumns
                TableColumnEditable
                return
            end
            UnitCell = TableUI.Data;
            if  Hybrid
                a.MaterialNameUI6.String = 'Hybrid';
            else
                a.MaterialNameUI6.String = Material{1}.Name;
            end
            MakeLayupString
            a.LayupUI6.String = LayupString;
            a.MaterialType6 = MaterialType;
            a.Material6 = Material;
            a.MaterialNames6 = MaterialNames;
            a.MaterialClasses6 = MaterialClasses;
            a.Hybrid6 = Hybrid;
            a.UniformLayerThickness6 = UniformLayerThickness;
            a.LayupThickness6 = LayupThickness;
            a.SuperLayerSize6 = SuperLayerSize;
            a.LayerOrientations6 = LayerOrientations;
            a.LayerThicknesses6 = LayerThicknesses;
            a.LayupString6 = LayupString;
            a.UnitCell6 = UnitCell;
            DisplayProperties
            close(f4)
            a.CLaminate = Computer_LaminateStiffness(1,a.Material6,a.LayerOrientations6,a.LayerThicknesses6,a.PropagationAngle6);
            a.CLP = Computer_LaminateStiffness(2,a.Material6,a.LayerOrientations6,a.LayerThicknesses6,1);
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
        elseif strcmp(source.Tag,'17') % Cancel
            close(f4)
        end
    end
    function TableCallback(hObject,callbackdata)
        if  callbackdata.Indices(2) < 7
            if  callbackdata.Indices(2) < 3
                if  isempty(callbackdata.EditData)
                    hObject.Data{callbackdata.Indices(1),callbackdata.Indices(2)} = callbackdata.PreviousData;
                    errordlg('Invalid input!','Error');
                    return
                end
                hObject.Data{callbackdata.Indices(1),callbackdata.Indices(2)} = str2double(replace(callbackdata.EditData,',','.'));
                if  callbackdata.Indices(2) == 1
                    LayerOrientations(callbackdata.Indices(1)) = hObject.Data{callbackdata.Indices(1),1};
                elseif callbackdata.Indices(2) == 2
                    if  hObject.Data{callbackdata.Indices(1),2} <= 0
                        hObject.Data{callbackdata.Indices(1),2} = 1;
                        errordlg('The layer thickness must be greater than zero.','Error');
                    end
                    LayerThicknesses(callbackdata.Indices(1)) = hObject.Data{callbackdata.Indices(1),2};
                end
            end
            z = cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),1));
            hObject.Data(z,1) = {0};
            LayerOrientations(z) = 0;
            if  ~UniformLayerThickness
                z = cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),2));
                hObject.Data(z,2) = {1};
                LayerThicknesses(z) = 1;
            end
            SuperLayerSize = length(LayerOrientations);
            if  Hybrid
                if  callbackdata.Indices(2) == 3
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Orthotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),4:6) = {''};
                elseif callbackdata.Indices(2) == 4
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.TransverselyIsotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),[3 5 6]) = {''};
                elseif callbackdata.Indices(2) == 5
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Cubic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),[3 4 6]) = {''};
                elseif callbackdata.Indices(2) == 6
                    Material{callbackdata.Indices(1)} = getfield(a.Materials.Isotropic,callbackdata.EditData);
                    hObject.Data(callbackdata.Indices(1),3:5) = {''};
                end
                z = all(cellfun(@isempty,hObject.Data(1:callbackdata.Indices(1),3:6)),2);
                for i = 1:callbackdata.Indices(1)
                    if  z(i)
                        Material{i} = getfield(a.Materials.Orthotropic,a.FieldNames_Orthotropic{1});
                        hObject.Data{i,3} = Material{i}.Name;
                    end
                end
            else
                Material(1:SuperLayerSize) = Material(1);
                hObject.Data(1:SuperLayerSize,MaterialType+2) = {Material{1}.Name};
                if  MaterialType == 1
                    hObject.Data(1:SuperLayerSize,4:6) = {''};
                elseif MaterialType == 2
                    hObject.Data(1:SuperLayerSize,[3 5 6]) = {''};
                elseif MaterialType == 3
                    hObject.Data(1:SuperLayerSize,[3 4 6]) = {''};
                elseif MaterialType == 4
                    hObject.Data(1:SuperLayerSize,3:5) = {''};
                end
            end
        else
            if  SuperLayerSize > 1
                if  SuperLayerSize >= callbackdata.Indices(1)
                    LayerOrientations(callbackdata.Indices(1)) = [];
                    LayerThicknesses(callbackdata.Indices(1)) = [];
                    Material(callbackdata.Indices(1)) = [];
                end
                SuperLayerSize = length(LayerOrientations);
                hObject.Data(callbackdata.Indices(1),:) = [];
                hObject.Data = [hObject.Data;cell(1,7)];
            else
                errordlg('Do not delete the only layer!','Error');
                return
            end
        end
        SetThicknesses
        MakeLayupString
        LayupUI.String = LayupString;
    end
end
function CallbackUI1(source,~) % Tab1_isotropic
    a = CallbackModule_Isotropic(source,a,Tab3); 
end
function CallbackUI2(source,~) % Tab2_anisotropic
    a = CallbackModule_Anisotropic(source,a,Tab3);
end
function CallbackUI3(source,~) % Tab3_signal simulator
    a = CallbackModule_SignalSimulator(source,a,Tab3);
end
function CallbackUI4(source,~) % Tab4_polar diagrams
    a = CallbackModule_PolarDiagrams(source,a);
end
function CallbackUI8(source,~) % Tab8_bulk waves
    a = CallbackModule_BulkWaves(source,a);
end
function CallbackUI6(source,~) % Tab6_laminate stiffness
    a = CallbackModule_LaminateStiffness(source,a);
end
function CallbackUI5(source,~) % Tab5_material editor
    a = CallbackModule_MaterialEditor(source,a);
end
function CallbackUI7(source,~) % Tab7_advanced
    a = CallbackModule_Advanced(source,a);
end
end