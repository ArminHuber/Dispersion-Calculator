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
function a = CallbackModule_BulkWaves(source,~,a)
%#ok<*GFLD>
%#ok<*NASGU>
source.String = replace(source.String,',','.');
if  strcmp(source.Tag,'43') % Class
    a.MaterialType_Bulk = source.Value;
    a.MaterialUI8.Value = 1;
    switch source.Value
    case 1
        a.MaterialUI8.String = fieldnames(a.Materials.Orthotropic);
        E = fieldnames(a.Materials.Orthotropic);
        a.Material_Bulk = getfield(a.Materials.Orthotropic,E{1});
    case 2
        a.MaterialUI8.String = fieldnames(a.Materials.TransverselyIsotropic);
        E = fieldnames(a.Materials.TransverselyIsotropic);
        a.Material_Bulk = getfield(a.Materials.TransverselyIsotropic,E{1});
    case 3
        a.MaterialUI8.String = fieldnames(a.Materials.Cubic);
        E = fieldnames(a.Materials.Cubic);
        a.Material_Bulk = getfield(a.Materials.Cubic,E{1});
    case 4
        a.MaterialUI8.String = fieldnames(a.Materials.Isotropic);
        E = fieldnames(a.Materials.Isotropic);
        a.Material_Bulk = getfield(a.Materials.Isotropic,E{1});
    end
    a.CalculateUI8.Enable = 'on';
    a.Plot2UI8.Enable = 'off';
elseif strcmp(source.Tag,'1') % Material
    if  a.MaterialType_Bulk == 1
        E = fieldnames(a.Materials.Orthotropic); 
        a.Material_Bulk = getfield(a.Materials.Orthotropic,cell2mat(source.String(source.Value)));        
    elseif a.MaterialType_Bulk == 2
        E = fieldnames(a.Materials.TransverselyIsotropic); 
        a.Material_Bulk = getfield(a.Materials.TransverselyIsotropic,cell2mat(source.String(source.Value))); 
    elseif a.MaterialType_Bulk == 3
        E = fieldnames(a.Materials.Cubic); 
        a.Material_Bulk = getfield(a.Materials.Cubic,cell2mat(source.String(source.Value)));
    elseif a.MaterialType_Bulk == 4
        E = fieldnames(a.Materials.Isotropic); 
        a.Material_Bulk = getfield(a.Materials.Isotropic,cell2mat(source.String(source.Value)));
    end
    a.CalculateUI8.Enable = 'on';
    a.Plot2UI8.Enable = 'off';
elseif strcmp(source.Tag,'2') % Quantity
    a.Quantity_Bulk = source.Value;
elseif strcmp(source.Tag,'3') % Delta Theta (Elastic waves in bulk material)
    a.ThetaStep_Bulk1 = str2double(source.String);
    a.CalculateUI8.Enable = 'on';
    a.Plot2UI8.Enable = 'off';
elseif strcmp(source.Tag,'44') % Plane
    switch source.Value
    case 1
        a.Plane = 13;
    case 2
        a.Plane = 12;
    end
elseif strcmp(source.Tag,'4') % Phi (2-D profiles)
    a.Phi_Bulk11 = str2double(source.String);
elseif strcmp(source.Tag,'5') % Plot (2-D profiles)
    if  a.Quantity_Bulk == 1
        PhaseVelocity_BulkWave2D(a.Material_Bulk,a.Plane,a.Phi_Bulk11,a.ThetaStep_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.LineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.AxesTickFontSize_Bulk,a.ModeLabelFontSize_Bulk,a.PNGresolution_Bulk)
    elseif a.Quantity_Bulk == 2
        GroupVelocity_BulkWave2D(a.Material_Bulk,a.Plane,a.Phi_Bulk11,a.ThetaStep_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.LineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.AxesTickFontSize_Bulk,a.ModeLabelFontSize_Bulk,a.PNGresolution_Bulk)
    elseif a.Quantity_Bulk == 3
        Slowness_BulkWave2D(a.Material_Bulk,a.Plane,a.Phi_Bulk11,a.ThetaStep_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.LineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.AxesTickFontSize_Bulk,a.ModeLabelFontSize_Bulk,a.PNGresolution_Bulk)        
    elseif a.Quantity_Bulk == 4
        PolarizationSkew_BulkWave2D(a.Material_Bulk,a.Plane,a.Phi_Bulk11,a.ThetaStep_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.LineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.AxesTickFontSize_Bulk,a.ModeLabelFontSize_Bulk,a.PNGresolution_Bulk)
    elseif a.Quantity_Bulk == 5
        EnergySkew_BulkWave2D(a.Material_Bulk,a.Plane,a.Phi_Bulk11,a.ThetaStep_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.LineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.AxesTickFontSize_Bulk,a.ModeLabelFontSize_Bulk,a.PNGresolution_Bulk)
    end
elseif strcmp(source.Tag,'6') % Delta Phi
    a.PhiStep_Bulk = str2double(source.String);
    a.CalculateUI8.Enable = 'on';
    a.Plot2UI8.Enable = 'off';
elseif strcmp(source.Tag,'7') % Calculate
    [a.PhaseVelocityCartesian,a.GroupVelocityCartesian,a.SlownessCartesian,a.PolarizationSkewCartesian,a.EnergySkewCartesian,a.X3D] = Computer_BulkWaves(a.Material_Bulk,a.PhiStep_Bulk,a.ThetaStep_Bulk1,a.OutputWindowUI8);
    a.CalculateUI8.Enable = 'off';
    a.Plot2UI8.Enable = 'on';
elseif strcmp(source.Tag,'8') % Mode
    switch source.Value
    case 1
       a.Mode_Bulk = 'L';
    case 2
       a.Mode_Bulk = 'Sf';
    case 3
       a.Mode_Bulk = 'Ss';
    end
elseif strcmp(source.Tag,'9') % Plot 1 (3-D surfaces and vectors)
    if  a.Quantity_Bulk == 1
        PhaseVelocity_BulkWave3D(a.Material_Bulk,a.Title_Bulk,a.TitleFontSize_Bulk,a.BoxLineWidth_Bulk,a.Mode_Bulk,a.MarkerSize_Bulk,a.ColorbarX_Bulk,a.ViewPhi_Bulk1,a.ViewTheta_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.AxesTickFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.PNGresolution_Bulk,a.PhaseVelocityCartesian,a.X3D)        
    elseif a.Quantity_Bulk == 2
        GroupVelocity_BulkWave3D(a.Material_Bulk,a.Title_Bulk,a.TitleFontSize_Bulk,a.BoxLineWidth_Bulk,a.Mode_Bulk,a.MarkerSize_Bulk,a.ColorbarX_Bulk,a.ViewPhi_Bulk1,a.ViewTheta_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.AxesTickFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.PNGresolution_Bulk,a.GroupVelocityCartesian,a.X3D)
    elseif a.Quantity_Bulk == 3
        Slowness_BulkWave3D(a.Material_Bulk,a.Title_Bulk,a.TitleFontSize_Bulk,a.BoxLineWidth_Bulk,a.Mode_Bulk,a.MarkerSize_Bulk,a.ColorbarX_Bulk,a.ViewPhi_Bulk1,a.ViewTheta_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.AxesTickFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.PNGresolution_Bulk,a.SlownessCartesian,a.X3D)        
    elseif a.Quantity_Bulk == 4
        PolarizationSkew_BulkWave3D(a.Material_Bulk,a.Title_Bulk,a.TitleFontSize_Bulk,a.BoxLineWidth_Bulk,a.Mode_Bulk,a.MarkerSize_Bulk,a.ColorbarX_Bulk,a.ViewPhi_Bulk1,a.ViewTheta_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.AxesTickFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.PNGresolution_Bulk,a.PolarizationSkewCartesian,a.X3D)
    elseif a.Quantity_Bulk == 5
        EnergySkew_BulkWave3D(a.Material_Bulk,a.Title_Bulk,a.TitleFontSize_Bulk,a.BoxLineWidth_Bulk,a.Mode_Bulk,a.MarkerSize_Bulk,a.ColorbarX_Bulk,a.ViewPhi_Bulk1,a.ViewTheta_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.AxesTickFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.PNGresolution_Bulk,a.EnergySkewCartesian,a.X3D)
    end
elseif strcmp(source.Tag,'13') % View Phi (3-D surfaces and vectors)
    a.ViewPhi_Bulk1 = str2double(source.String);
elseif strcmp(source.Tag,'14') % View Theta (3-D surfaces and vectors)
    a.ViewTheta_Bulk1 = str2double(source.String);
elseif strcmp(source.Tag,'15') % Marker size
    a.MarkerSize_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'16') % Colorbar x-pos.
    a.ColorbarX_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'10') % Phi (3-D surfaces and vectors)
    a.Phi_Bulk21 = str2double(source.String);
elseif strcmp(source.Tag,'11') % Theta (3-D surfaces and vectors)
    a.Theta_Bulk1 = str2double(source.String);
elseif strcmp(source.Tag,'12') % Plot 2 (3-D surfaces and vectors)
    BulkWaves3D(a.ModeLabelFontSize_Bulk,a.Material_Bulk,a.AxesLabelFontSize_Bulk,a.BoxLineWidth_Bulk,a.ViewTheta_Bulk1,a.ViewPhi_Bulk1,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.TitleFontSize_Bulk,a.AxesTickFontSize_Bulk,a.WaveVectorLineWidth_Bulk,a.PNGresolution_Bulk,a.Phi_Bulk21,a.Theta_Bulk1)
elseif strcmp(source.Tag,'17') % Fluid
    E = fieldnames(a.Materials.Fluid);
    a.Couplant_Bulk = getfield(a.Materials.Fluid,E{source.Value});
    a.FluidVelocityUI8.String = a.Couplant_Bulk.Velocity/1e3;
    if  a.SolidType_Bulk == 1
        [a.BulkWaves,a.X,a.Y] = Computer_Isotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Theta_Bulk2,a.OutputWindowUI8);
    elseif a.SolidType_Bulk > 1
        [a.BulkWaves,a.X,a.Y] = Computer_Anisotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Phi_Bulk2,a.Theta_Bulk2,a.OutputWindowUI8);
    end
elseif strcmp(source.Tag,'19') % Solid type
    a.SolidType_Bulk = source.Value;
    a.SolidUI8.Value = 1;
    switch source.Value
    case 1
        a.SolidUI8.String = fieldnames(a.Materials.Isotropic);
        E = fieldnames(a.Materials.Isotropic);
        a.Solid_Bulk = getfield(a.Materials.Isotropic,E{1});
        a.Phi2UI8.Enable = 'off';
        a.ViewPhi2UI8.Enable = 'off';
        a.ViewTheta2UI8.Enable = 'off';
        a.Plot3DUI8.Enable = 'off';
    case 2
        a.SolidUI8.String = fieldnames(a.Materials.Cubic);
        E = fieldnames(a.Materials.Cubic);
        a.Solid_Bulk = getfield(a.Materials.Cubic,E{1});
        a.Phi2UI8.Enable = 'on';
        a.ViewPhi2UI8.Enable = 'on';
        a.ViewTheta2UI8.Enable = 'on';
        a.Plot3DUI8.Enable = 'on';
    case 3
        a.SolidUI8.String = fieldnames(a.Materials.TransverselyIsotropic);
        E = fieldnames(a.Materials.TransverselyIsotropic);
        a.Solid_Bulk = getfield(a.Materials.TransverselyIsotropic,E{1});
        a.Phi2UI8.Enable = 'on';
        a.ViewPhi2UI8.Enable = 'on';
        a.ViewTheta2UI8.Enable = 'on';
        a.Plot3DUI8.Enable = 'on';
    case 4
        a.SolidUI8.String = fieldnames(a.Materials.Orthotropic);
        E = fieldnames(a.Materials.Orthotropic);
        a.Solid_Bulk = getfield(a.Materials.Orthotropic,E{1});
        a.Phi2UI8.Enable = 'on';
        a.ViewPhi2UI8.Enable = 'on';
        a.ViewTheta2UI8.Enable = 'on';
        a.Plot3DUI8.Enable = 'on';
    end
    if  a.SolidType_Bulk == 1
        [a.BulkWaves,a.X,a.Y] = Computer_Isotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Theta_Bulk2,a.OutputWindowUI8);
    elseif a.SolidType_Bulk > 1
        [a.BulkWaves,a.X,a.Y] = Computer_Anisotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Phi_Bulk2,a.Theta_Bulk2,a.OutputWindowUI8);
    end
elseif strcmp(source.Tag,'20') % Solid
    if  a.SolidType_Bulk == 1 
    	a.Solid_Bulk = getfield(a.Materials.Isotropic,cell2mat(source.String(source.Value)));
    elseif a.SolidType_Bulk == 2
        a.Solid_Bulk = getfield(a.Materials.Cubic,cell2mat(source.String(source.Value)));
    elseif a.SolidType_Bulk == 3
        a.Solid_Bulk = getfield(a.Materials.TransverselyIsotropic,cell2mat(source.String(source.Value)));
    elseif a.SolidType_Bulk == 4
        a.Solid_Bulk = getfield(a.Materials.Orthotropic,cell2mat(source.String(source.Value)));
    end
    if  a.SolidType_Bulk == 1
        [a.BulkWaves,a.X,a.Y] = Computer_Isotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Theta_Bulk2,a.OutputWindowUI8);
    elseif a.SolidType_Bulk > 1
        [a.BulkWaves,a.X,a.Y] = Computer_Anisotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Phi_Bulk2,a.Theta_Bulk2,a.OutputWindowUI8);
    end
elseif strcmp(source.Tag,'21') % Phi (Bulk waves on interfaces)
    a.Phi_Bulk2 = str2double(source.String);
    [a.BulkWaves,a.X,a.Y] = Computer_Anisotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Phi_Bulk2,a.Theta_Bulk2,a.OutputWindowUI8);
elseif strcmp(source.Tag,'22') % Theta (Bulk waves on interfaces)
    a.Theta_Bulk2 = str2double(source.String);
    if  a.SolidType_Bulk == 1
        [a.BulkWaves,a.X,a.Y] = Computer_Isotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Theta_Bulk2,a.OutputWindowUI8);
    elseif a.SolidType_Bulk > 1
        [a.BulkWaves,a.X,a.Y] = Computer_Anisotropic_SnellsLaw(a.Couplant_Bulk,a.Solid_Bulk,a.Phi_Bulk2,a.Theta_Bulk2,a.OutputWindowUI8);
    end
elseif strcmp(source.Tag,'23') % Delta Theta (Bulk waves on interfaces)
    a.ThetaStep_Bulk2 = str2double(source.String);
elseif strcmp(source.Tag,'24') % View Phi (Bulk waves on interfaces)
    a.ViewPhi_Bulk2 = str2double(source.String);
elseif strcmp(source.Tag,'25') % View Theta (Bulk waves on interfaces)
    a.ViewTheta_Bulk2 = str2double(source.String);
elseif strcmp(source.Tag,'26') % Plot 2-D
    if  a.SolidType_Bulk == 1
        SnellsLaw2D_Isotropic(a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,1,a.LineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesTickFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.BulkWaves,a.ModeLabelFontSize_Bulk,a.WaveVectorLineWidth_Bulk,a.PNGresolution_Bulk,a.X,a.Y,a.Couplant_Bulk,a.Solid_Bulk,a.Theta_Bulk2)
    elseif a.SolidType_Bulk > 1
        SnellsLaw2D_Anisotropic(a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,1,a.LineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesTickFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.BulkWaves,a.ModeLabelFontSize_Bulk,a.WaveVectorLineWidth_Bulk,a.PNGresolution_Bulk,a.X,a.Y,a.Couplant_Bulk,a.Solid_Bulk,a.Phi_Bulk2,a.Theta_Bulk2)
    end
elseif strcmp(source.Tag,'27') % Plot 3-D
    SnellsLaw3D_Anisotropic(a.ModeLabelFontSize_Bulk,a.Solid_Bulk,a.AxesLabelFontSize_Bulk,a.BoxLineWidth_Bulk,a.ViewTheta_Bulk2,a.ViewPhi_Bulk2,a.ExportPlots_Bulk,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.TitleFontSize_Bulk,a.AxesTickFontSize_Bulk,a.BulkWaves,a.WaveVectorLineWidth_Bulk,a.PNGresolution_Bulk,a.Phi_Bulk2,a.Theta_Bulk2)
elseif strcmp(source.Tag,'46') % Plot R,T
    if  a.SolidType_Bulk == 1
        EnergyScattering_Isotropic(a.Couplant_Bulk,a.Solid_Bulk,a.ExportPlots_Bulk,a.Theta_Bulk2,a.ThetaStep_Bulk2,2,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.LineWidth_Bulk,a.BoxLineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.AxesTickFontSize_Bulk,a.PNGresolution_Bulk)
    elseif a.SolidType_Bulk > 1
        EnergyScattering_Anisotropic(a.Couplant_Bulk,a.Solid_Bulk,a.ExportPlots_Bulk,a.Phi_Bulk2,a.Theta_Bulk2,a.ThetaStep_Bulk2,2,a.PDF_Bulk,a.PNG_Bulk,a.FileName_Bulk,a.Directory8,a.Title_Bulk,a.LineWidth_Bulk,a.BoxLineWidth_Bulk,a.TitleFontSize_Bulk,a.AxesLabelFontSize_Bulk,a.AxesTickFontSize_Bulk,a.PNGresolution_Bulk)    
    end
elseif strcmp(source.Tag,'36') % Export plots
    a.ExportPlots_Bulk = source.Value;
    if  source.Value
        a.Plot1UI8.String = 'Export';
        a.Plot2UI8.String = 'Export';
        a.Plot3UI8.String = 'Export';
        a.Plot2DUI8.String = 'Export 2-D';
        a.Plot3DUI8.String = 'Export 3-D';
        a.PlotRTUI8.String = 'Export R,T';
    else
        a.Plot1UI8.String = 'Plot';
        a.Plot2UI8.String = 'Plot';
        a.Plot3UI8.String = 'Plot';
        a.Plot2DUI8.String = 'Plot 2-D';
        a.Plot3DUI8.String = 'Plot 3-D';
        a.PlotRTUI8.String = 'Plot R,T';
    end
elseif strcmp(source.Tag,'37') % PDF
    a.PDF_Bulk = source.Value;
elseif strcmp(source.Tag,'38') % PNG
    a.PNG_Bulk = source.Value;
elseif strcmp(source.Tag,'39') % PNG resolution
    a.PNGresolution_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'40') % File name
    a.FileName_Bulk = source.String;
elseif strcmp(source.Tag,'41') % Directory
    a.Directory8 = source.String;
elseif strcmp(source.Tag,'28') % Title
    a.Title_Bulk = source.Value;
elseif strcmp(source.Tag,'29') % Box
    a.BoxLineWidth_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'30') % Profile
    a.LineWidth_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'31') % Bulk wave
    a.WaveVectorLineWidth_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'32') % Title (Font size)
    a.TitleFontSize_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'33') % Axes labels
    a.AxesLabelFontSize_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'34') % Axes ticks
    a.AxesTickFontSize_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'35') % Mode labels
    a.ModeLabelFontSize_Bulk = str2double(source.String);
elseif strcmp(source.Tag,'42') % Default
    a.Title_Bulk = 1;
    a.TitleUI8.Value = a.Title_Bulk;
    a.BoxLineWidth_Bulk = .5;
    a.BoxLineWidthUI8.String = a.BoxLineWidth_Bulk;
    a.LineWidth_Bulk = 1;
    a.LineWidthUI8.String = a.LineWidth_Bulk;
    a.WaveVectorLineWidth_Bulk = 2;
    a.WaveVectorLineWidthUI8.String = a.WaveVectorLineWidth_Bulk;
    a.TitleFontSize_Bulk = 30;
    a.TitleFontSizeUI8.String = a.TitleFontSize_Bulk;
    a.AxesLabelFontSize_Bulk = 30;
    a.AxesLabelFontSizeUI8.String = a.AxesLabelFontSize_Bulk;
    a.AxesTickFontSize_Bulk = 24;
    a.AxesTickFontSizeUI8.String = a.AxesTickFontSize_Bulk;
    a.ModeLabelFontSize_Bulk = 24;
    a.ModeLabelFontSizeUI8.String = a.ModeLabelFontSize_Bulk;
end