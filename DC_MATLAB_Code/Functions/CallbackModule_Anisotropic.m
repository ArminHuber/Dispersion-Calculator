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
function a = CallbackModule_Anisotropic(source,~,a,Tab2,Tab3)
%#ok<*AGROW>
%#ok<*GFLD>
%#ok<*GVMIS>
global Stop  
Stop = 0;
if  strcmp(source.Tag,'1') % Propagation angle
    a.PropagationAngle = str2double(source.String);
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
    end
    if  all(DC)
        a.ShearHorizontalModesUI2.Enable = 'on';
    else
        a.ShearHorizontalModesUI2.Enable = 'off';
    end
    a.Detect2 = 0;
    a.EffectiveLayupString1 = char(join(split(num2str(a.LayerOrientations-a.PropagationAngle)),'/'));
    if  a.SymmetricSystem == 0
        if  a.SuperLayers == 1
            a.EffectiveLayupUI2.String = ['[',a.EffectiveLayupString1,']'];
        elseif a.SuperLayers > 1
            a.EffectiveLayupUI2.String = ['[',a.EffectiveLayupString1,']',num2str(a.SuperLayers)];
        end
    else
        if  a.SuperLayers == 1
            a.EffectiveLayupUI2.String = ['[',a.EffectiveLayupString1,']s'];
        elseif a.SuperLayers > 1
            a.EffectiveLayupUI2.String = ['[',a.EffectiveLayupString1,']',num2str(a.SuperLayers),'s'];
        end
    end
elseif strcmp(source.Tag,'2') % Phase velocity limit
    a.PhaseVelocityLimit2 = str2double(source.String)*1e3;

    if  a.Quantity12 == 1
        a.YAxisUI2.String = ['[0 ',source.String,']'];
        a.YAxis2 = eval(a.YAxisUI2.String);
    end

    a.Detect2 = 0;
elseif strcmp(source.Tag,'3') % Frequency limit
    a.FrequencyLimit2 = str2double(source.String);
    
    if  a.Quantity12 == 3
        if  a.XAxisMode2 == 1
            a.YAxisUI2.String = ['[0 ',source.String,']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        elseif a.XAxisMode2 == 2
            a.YAxisUI2.String = ['[0 ',num2str(eval(source.String)/1e3),']'];
            a.YAxis2 = eval(a.YAxisUI2.String)*1e3;
        elseif a.XAxisMode2 == 3
            a.YAxisUI2.String = ['[0 ',num2str(eval(source.String)/1e3*a.PlateThickness),']'];
            a.YAxis2 = eval(a.YAxisUI2.String)*1e3/a.PlateThickness;                    
        end 
    else
        if  a.XAxisMode2 == 1
            a.XAxisUI2.String = ['[0 ',source.String,']'];
            a.XAxis2 = eval(a.XAxisUI2.String);
        elseif a.XAxisMode2 == 2
            a.XAxisUI2.String = ['[0 ',num2str(eval(source.String)/1e3),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
        elseif a.XAxisMode2 == 3
            a.XAxisUI2.String = ['[0 ',num2str(eval(source.String)/1e3*a.PlateThickness),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;                    
        end       
    end
    
    if  a.Quantity12 == 6
        x = 2*pi*a.FrequencyLimit2/1.6e3;
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
    end       
    
    a.Frequency12 = a.FrequencyLimit2;
    a.Frequency1UI2.String = a.FrequencyLimit2;

    a.Frequency22 = a.FrequencyLimit2;
    a.Frequency2UI2.String = a.FrequencyLimit2;
    
    a.Detect2 = 0;
elseif strcmp(source.Tag,'4') % Frequency step
    a.FrequencyResolution2 = str2double(source.String);
elseif strcmp(source.Tag,'7') % Higher order modes
    a.HigherOrderModes2 = source.Value;
elseif strcmp(source.Tag,'9') % Symmetric modes
    a.SymmetricModes2 = source.Value;
elseif strcmp(source.Tag,'10') % Antiymmetric modes
    a.AntisymmetricModes2 = source.Value;
elseif strcmp(source.Tag,'13') % Lamb modes
    a.LambModes2 = source.Value;
elseif strcmp(source.Tag,'14') % Shear horizontal modes
    a.ShearHorizontalModes2 = source.Value;
elseif strcmp(source.Tag,'5') % Scholte modes
    a.ScholteModes2 = source.Value;   
elseif strcmp(source.Tag,'15') % Step
    a.Step2 = str2double(source.String);
    a.Detect2 = 0;
elseif strcmp(source.Tag,'16') % Detect
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
    end
    if  all(DC)
        a.Decoupled = 1;
    else
        a.Decoupled = 0;
        for i = 1:length(a.MaterialClasses)
            if  DC(i)
                Phi(i) = Phi(i)+1e-3;
            end
        end
    end    
    a.SuperLayerSize = length(Phi);
    for i = 1:a.SuperLayerSize % transformed layer stiffness matrices for orthotropic materials [2], p. 28; transverse isotropy is included
        C = conj(a.Material2{i}.C);
        s = sind(Phi(i));
        g = cosd(Phi(i));
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
        DeltaReal(i) = real(a.c{i}(3,3))*real(a.c{i}(4,4))*real(a.c{i}(5,5))-real(a.c{i}(3,3))*real(a.c{i}(4,5))^2;
        if  strcmp(a.MaterialClasses(i),'Isotropic') && a.Decoupled == 0
            a.c{i}(1,6) = 1e3;
            a.c{i}(2,6) = 1e3;
            a.c{i}(3,6) = 1e3;
            a.c{i}(4,5) = 1e3;
        end
    end
    [a.HS,a.HSLamb2,a.HSShear2,a.HA,a.HALamb2,a.HAShear2,a.HB,a.HBLamb,a.HBShear] = FrequencySweeper_Anisotropic(a.c,DeltaReal,a.Material2,100/a.PlateThickness:a.Step2:a.FrequencyLimit2,a.PhaseVelocityLimit2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.Symmetric,a.SymmetricSystem,a.Decoupled,a.OutputWindow1aUI2,a.OutputWindow1bUI2,a.OutputWindow2aUI2,a.OutputWindow2bUI2);
    if  a.FluidLoading2 == 1 && a.Symmetric == 1
        [a.HSScholte2,a.HAScholte2,a.HBScholte2] = FrequencySweeper_Anisotropic_Scholte(a.c,DeltaReal,a.Material2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,100/a.PlateThickness:a.Step2:a.FrequencyLimit2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.Symmetric,a.SymmetricSystem,a.Decoupled);
    else
        a.HSScholte2 = [];
        a.HAScholte2 = [];
        a.HBScholte2 = [];
    end
    a.Detect2 = 1;
elseif strcmp(source.Tag,'17') % Calculate
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
    end
    if  all(DC)
        a.Decoupled = 1;
    else
        a.Decoupled = 0;
        for i = 1:length(a.MaterialClasses)
            if  DC(i)
                Phi(i) = Phi(i)+1e-3;
            end
        end
    end
    a.SuperLayerSize = length(Phi);
    for i = 1:a.SuperLayerSize % transformed layer stiffness matrices for orthotropic materials [2], p. 28; transverse isotropy is included
        C = conj(a.Material2{i}.C);
        s = sind(Phi(i));
        g = cosd(Phi(i));
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
        Delta(i) = a.c{i}(3,3)*a.c{i}(4,4)*a.c{i}(5,5)-a.c{i}(3,3)*a.c{i}(4,5)^2;
        DeltaReal(i) = real(a.c{i}(3,3))*real(a.c{i}(4,4))*real(a.c{i}(5,5))-real(a.c{i}(3,3))*real(a.c{i}(4,5))^2;
        if  strcmp(a.MaterialClasses(i),'Isotropic') && a.Decoupled == 0
            a.c{i}(1,6) = 1e3;
            a.c{i}(2,6) = 1e3;
            a.c{i}(3,6) = 1e3;
            a.c{i}(4,5) = 1e3;
        end
    end
    if  a.Detect2 == 0 && a.HigherOrderModes2 == 1 && (a.LambModes2 == 1 || a.ShearHorizontalModes2 == 1)
        [a.HS,a.HSLamb2,a.HSShear2,a.HA,a.HALamb2,a.HAShear2,a.HB,a.HBLamb,a.HBShear] = FrequencySweeper_Anisotropic(a.c,DeltaReal,a.Material2,100/a.PlateThickness:a.Step2:a.FrequencyLimit2,a.PhaseVelocityLimit2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.Symmetric,a.SymmetricSystem,a.Decoupled,a.OutputWindow1aUI2,a.OutputWindow1bUI2,a.OutputWindow2aUI2,a.OutputWindow2bUI2);
        a.Detect2 = 1;
    end
    if  a.FluidLoading2 == 1 && a.HigherOrderModes2 == 1 && a.ScholteModes2 == 1 && a.Symmetric == 1
        [a.HSScholte2,a.HAScholte2,a.HBScholte2] = FrequencySweeper_Anisotropic_Scholte(a.c,DeltaReal,a.Material2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,100/a.PlateThickness:a.Step2:a.FrequencyLimit2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.Symmetric,a.SymmetricSystem,a.Decoupled);
    else
        a.HSScholte2 = [];
        a.HAScholte2 = [];
        a.HBScholte2 = [];
    end
    FrequencyRange = 0:a.FrequencyResolution2:a.FrequencyLimit2;
    FrequencyRangeF = FrequencyRange;
    FrequencyRange(1) = a.FrequencyRangeStart2;
    FrequencyRangeF(1) = .1*a.FrequencyResolution2;
    if  a.FluidLoading2 == 1 || a.Viscoelastic2 == 1
        [FLambF,FScholte] = PhaseVelocitySweeper_Anisotropic_F(a.c,Delta,a.Material2,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.FluidDensityThreshold2,FrequencyRangeF(1),[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.Symmetric,a.SymmetricSystem,a.Decoupled);
    else
        FLambF = [];
        FScholte = [];
    end    
    if  a.Symmetric == 1 && ((a.Decoupled == 0 && a.SymmetricModes2 == 1 && a.AntisymmetricModes2 == 1 && a.Multithreading == 1)  || (a.Decoupled == 1 && a.LambModes2 == 1 && a.SymmetricModes2 == 1 && a.AntisymmetricModes2 == 1 && a.Multithreading == 1))
        delete(a.a1UI2)
        a.a1UI2 = axes('Parent',Tab2,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading22.png','Parent',a.a1UI2)
    end
    [a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2] = Computer_Anisotropic(a.Multithreading,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.FluidDensityThreshold2,a.Hybrid,a.FrequencyLimit2,a.Material2,a.PropagationAngle,a.PhaseVelocityStep2,a.FrequencyOffset2,a.ShearPhaseVelocitySweepRange2,a.LambPhaseVelocitySweepRange12,a.LambPhaseVelocitySweepRange22,a.SearchWidthReal2,a.SearchWidthImag2,a.SearchAreaSections2,a.SearchAreaExtensions2,a.LayerThicknesses/1e3,FrequencyRange,FrequencyRangeF,FLambF,FScholte,a.HS,a.HSLamb2,a.HSShear2,a.HSScholte2,a.HA,a.HALamb2,a.HAShear2,a.HAScholte2,a.HB,a.HBLamb,a.HBShear,a.HBScholte2,a.c,Delta,[1 1 -1;-1 -1 1;-1 -1 1],[1 -1;-1 1],a.PlateThickness/1e3,a.SuperLayerSize,a.Decoupled,a.SuperLayers,a.SymmetricSystem,a.Symmetric,a.FrequencyResolution2,a.PhaseVelocityLimit2,a.PhaseVelocityResolution2,a.PhaseVelocitySections2,a.FrequencySections2,a.HigherOrderModes2,a.SymmetricModes2,a.AntisymmetricModes2,a.LambModes2,a.ShearHorizontalModes2,a.ScholteModes2,a.MissingSamples2,a.BelowCutoffWidth2);
    if  Stop == 1
        return
    end
    [a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2] = Computer_Anisotropic_EnergyVelocity(a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.c,a.Material2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Symmetric,a.Decoupled,Delta);
    clear Computer_Anisotropic_SLamb_F Computer_Anisotropic_ALamb_F Computer_Anisotropic_BLamb_F Computer_Anisotropic_S_F Computer_Anisotropic_A_F Computer_Anisotropic_B_F Computer_Anisotropic_AScholte_V Computer_Anisotropic_SScholte_V Computer_Anisotropic_AScholte_Coupled_V Computer_Anisotropic_SScholte_Coupled_V Computer_Anisotropic_AShear_V Computer_Anisotropic_SShear_V Computer_Anisotropic_BShear_V
    if  a.Symmetric == 1 && ((a.Decoupled == 0 && a.SymmetricModes2 == 1 && a.AntisymmetricModes2 == 1 && a.Multithreading == 1)  || (a.Decoupled == 1 && a.LambModes2 == 1 && a.SymmetricModes2 == 1 && a.AntisymmetricModes2 == 1 && a.Multithreading == 1))
        delete(a.a1UI2)
        a.a1UI2 = axes('Parent',Tab2,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading21.png','Parent',a.a1UI2)
    end
    a.ModeNames2 = {''};
    if  a.Symmetric == 1
        if  a.Decoupled == 0
            if  ~isempty(a.A{1})
                for i = 0:length(a.A)-1
                    eval(sprintf('a.ModeNames2(i+1) = {''A%u''};',i));
                end
            end
            if  ~isempty(a.S{1})
                for i = 0:length(a.S)-1
                    eval(sprintf('a.ModeNames2(length(a.A)+i+1) = {''S%u''};',i));
                end
            end
            if  ~isempty(a.AScholte2{1})
                for i = 0:length(a.AScholte2)-1
                    eval(sprintf('a.ModeNames2(length(a.A)+length(a.S)+i+1) = {''AScholte%u''};',i));
                end
            end
            if  ~isempty(a.SScholte2{1})
                for i = 0:length(a.SScholte2)-1
                    eval(sprintf('a.ModeNames2(length(a.A)+length(a.S)+length(a.AScholte2)+i+1) = {''SScholte%u''};',i));
                end
            end
        elseif a.Decoupled == 1 
            if  ~isempty(a.ALamb2{1})
                for i = 0:length(a.ALamb2)-1
                    eval(sprintf('a.ModeNames2(i+1) = {''A%u''};',i));
                end
            end
            if  ~isempty(a.SLamb2{1})
                for i = 0:length(a.SLamb2)-1
                    eval(sprintf('a.ModeNames2(length(a.ALamb2)+i+1) = {''S%u''};',i));
                end
            end
            if  ~isempty(a.AShear2{1})
                for i = 1:length(a.AShear2)
                    eval(sprintf('a.ModeNames2(length(a.ALamb2)+length(a.SLamb2)+i) = {''ASH%u''};',i));
                end
            end
            if  ~isempty(a.SShear2{1})
                for i = 0:length(a.SShear2)-1
                    eval(sprintf('a.ModeNames2(length(a.ALamb2)+length(a.SLamb2)+length(a.AShear2)+i+1) = {''SSH%u''};',i));
                end
            end
            if  ~isempty(a.AScholte2{1})
                for i = 0:length(a.AScholte2)-1
                    eval(sprintf('a.ModeNames2(length(a.ALamb2)+length(a.SLamb2)+length(a.AShear2)+length(a.SShear2)+i+1) = {''AScholte%u''};',i));
                end
            end
            if  ~isempty(a.SScholte2{1})
                for i = 0:length(a.SScholte2)-1
                    eval(sprintf('a.ModeNames2(length(a.ALamb2)+length(a.SLamb2)+length(a.AShear2)+length(a.SShear2)+length(a.AScholte2)+i+1) = {''SScholte%u''};',i));
                end
            end
        end
    else
        if  a.Decoupled == 0 
            if  ~isempty(a.B{1})
                for i = 0:length(a.B)-1
                    eval(sprintf('a.ModeNames2(i+1) = {''B%u''};',i));
                end
            end
            if  ~isempty(a.BScholte2{1})
                for i = 0:length(a.BScholte2)-1
                    eval(sprintf('a.ModeNames2(length(a.B)+i+1) = {''BScholte%u''};',i));
                end
            end
        elseif a.Decoupled == 1
            if  ~isempty(a.BLamb{1})
                for i = 0:length(a.BLamb)-1
                    eval(sprintf('a.ModeNames2(i+1) = {''B%u''};',i));
                end
            end
            if  a.SuperLayerSize > 1 && a.SymmetricSystem == 0
                if  ~isempty(a.BShear{1})
                    for i = 0:length(a.BShear)-1
                        eval(sprintf('a.ModeNames2(length(a.BLamb)+i+1) = {''BSH%u''};',i));
                    end
                end
                if  ~isempty(a.BScholte2{1})
                    for i = 0:length(a.BScholte2)-1
                        eval(sprintf('a.ModeNames2(length(a.BLamb)+length(a.BShear)+i+1) = {''BScholte%u''};',i));
                    end
                end
            elseif a.SuperLayerSize == 1 || a.SymmetricSystem == 1
                if  ~isempty(a.AShear2{1})
                    for i = 1:length(a.AShear2)
                        eval(sprintf('a.ModeNames2(length(a.BLamb)+i) = {''ASH%u''};',i));
                    end
                end
                if  ~isempty(a.SShear2{1})
                    for i = 0:length(a.SShear2)-1
                        eval(sprintf('a.ModeNames2(length(a.BLamb)+length(a.AShear2)+i+1) = {''SSH%u''};',i));
                    end
                end
                if  ~isempty(a.BScholte2{1})
                    for i = 0:length(a.BScholte2)-1
                        eval(sprintf('a.ModeNames2(length(a.BLamb)+length(a.AShear2)+length(a.SShear2)+i+1) = {''BScholte%u''};',i));
                    end
                end
            end
        end
    end
    for i = 1:length(a.ModeNames2)
        if  isempty(a.ModeNames2{i})
            z(i) = 0;
        else
            z(i) = 1;
        end
    end
    a.ModeNames2(z == 0) = [];
    a.Mode1UI2.Value = 1;
    a.Mode1UI2.String = a.ModeNames2;
    a.Mode12 = a.ModeNames2{1};
    a.Mode2UI2.Value = 1;
    a.Mode2UI2.String = a.ModeNames2;
    a.Mode22 = a.ModeNames2{1};
    a.ModeShapeSettingChanged2 = 1;
    
    a.Title3 = 2;
    a.TitleUI3.Style = 'popupmenu';
    a.TitleUI3.Value = 1;
    a.TitleUI3.String = {'with layup','without layup','no title'};
    a.TitleUI3.TooltipString = ['Choose to plot the title with the layup ',newline,'included, without layup, or no title at all.'];
    a.TitleUI3.Position(3) = 65;
    if  a.Hybrid == 0
        a.DataUI3.String = a.Material2{1}.Name;
    elseif a.Hybrid == 1
        a.DataUI3.String = 'Hybrid';
    end
    a.PlotUI3.Enable = 'off';
    a.DataType3 = 2;
    a.Mode3 = '';
    a.Frequency3 = a.FrequencyResolution2*50;
    a.FrequencyUI3.String = a.Frequency3;
    [a.AModes,a.ALambModes2,a.AShearModes2,a.BModes,a.BLambModes,a.BShearModes,a.SModes,a.SLambModes2,a.SShearModes2,a.Frequency3] = ModeFinder_Anisotropic(a.A,a.ALamb2,a.AShear2,a.B,a.BLamb,a.BShear,a.S,a.SLamb2,a.SShear2,a.Frequency3);
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
        if  a.Symmetric == 1
            eval(sprintf('a.ALamb%uaUI3.String = [''A'',num2str(i)];',i));
            eval(sprintf('a.SLamb%uaUI3.String = [''S'',num2str(i)];',i));
            eval(sprintf('a.SLamb%uaUI3.Position(3) = 16;',i));
            eval(sprintf('a.SLamb%ubUI3.Position(1) = 88;',i));
            eval(sprintf('a.SLamb%ucUI3.Position(1) = 103;',i));
        elseif a.Symmetric == 0
            if  a.Decoupled == 0
                eval(sprintf('a.ALamb%uaUI3.String = [''B'',num2str(i)];',i));
            elseif a.Decoupled == 1
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
elseif strcmp(source.Tag,'18') % Stop
    Stop = 1;
    clear Computer_Anisotropic_SLamb_F Computer_Anisotropic_ALamb_F Computer_Anisotropic_BLamb_F Computer_Anisotropic_S_F Computer_Anisotropic_A_F Computer_Anisotropic_B_F Computer_Anisotropic_AScholte_V Computer_Anisotropic_SScholte_V Computer_Anisotropic_AScholte_Coupled_V Computer_Anisotropic_SScholte_Coupled_V Computer_Anisotropic_AShear_V Computer_Anisotropic_SShear_V Computer_Anisotropic_BShear_V
elseif strcmp(source.Tag,'19') % Quantity (Dispersion diagrams)
    switch source.Value
    case 1
        a.Quantity12 = 1;
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
        if  a.XAxisMode2 == 1
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.XAxis2 = eval(a.XAxisUI2.String);
        elseif a.XAxisMode2 == 2
            a.XAxisTextUI2.String = 'X-axis (MHz)';
            a.XAxisTextUI2.Position(3) = 65;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
        elseif a.XAxisMode2 == 3
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI2.Position(3) = 86;
            a.XAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;
        end
        a.YAxisTextUI2.String = 'Y-axis (m/ms)';
        a.YAxisTextUI2.Position(3) = 70;
        a.YAxisUI2.TooltipString = 'Enter which phase velocity range shall be plotted.';
            
        a.YAxisUI2.String = ['[0 ',num2str(a.PhaseVelocityLimit2/1e3),']'];
        a.YAxis2 = eval(a.YAxisUI2.String);
    case 2
        a.Quantity12 = 2;
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
        if  a.XAxisMode2 == 1
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.XAxis2 = eval(a.XAxisUI2.String);            
        elseif a.XAxisMode2 == 2
            a.XAxisTextUI2.String = 'X-axis (MHz)';
            a.XAxisTextUI2.Position(3) = 65;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3;            
        elseif a.XAxisMode2 == 3
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI2.Position(3) = 86;
            a.XAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;            
        end
        a.YAxisTextUI2.String = 'Y-axis (m/ms)';
        a.YAxisTextUI2.Position(3) = 70;
        a.YAxisUI2.TooltipString = 'Enter which phase velocity range shall be plotted.';

        if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
            a.YAxisUI2.String = ['[0 ',num2str(ceil(a.Material2{1}.PlateVelocity/1e3)),']'];
        else
            a.YAxisUI2.String = '[0 11]';
        end
        a.YAxis2 = eval(a.YAxisUI2.String);
    case 3        
        a.Quantity12 = 3;
        a.Option1UI2.Enable = 'on';
        a.Option1UI2.Style = 'edit';
        a.Option1UI2.String = a.Distance2;
        a.Option1UI2.TooltipString = ['Enter the propagation distance from the source to the',newline,'sensor for which the propagation time shall be calculated.'];
        a.Option1UI2.Position(3) = 75;
        a.Option1TextUI2.String = 'Distance (mm)';
        a.Option1TextUI2.Position(3) = 71;
        a.XAxisModeTextUI2.String = 'Y-axis mode';
        a.XAxisModeTextUI2.Position(3) = 63;
        a.XAxisModeUI2.TooltipString = 'Select the frequency''s dimension on the Y-axis.';
        a.XAxisTextUI2.String = ['X-axis (',char(181),'s)'];
        a.XAxisTextUI2.Position(3) = 56;
        a.XAxisUI2.TooltipString = 'Enter which propagation time range shall be plotted.';
        if  a.XAxisMode2 == 1
            a.YAxisTextUI2.String = 'Y-axis (kHz)';
            a.YAxisTextUI2.Position(3) = 63;
            a.YAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.YAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);             
        elseif a.XAxisMode2 == 2
            a.YAxisTextUI2.String = 'Y-axis (MHz)';
            a.YAxisTextUI2.Position(3) = 66;
            a.YAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.YAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.YAxis2 = eval(a.YAxisUI2.String)*1e3;            
        elseif a.XAxisMode2 == 3
            a.YAxisTextUI2.String = ['Y-axis (MHz',char(8901),'mm)'];
            a.YAxisTextUI2.Position(3) = 87;
            a.YAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
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
    case 4
        a.Quantity12 = 4;
        a.Option1UI2.Enable = 'on';
        a.Option1UI2.Style = 'popupmenu';
        a.Option1UI2.String = fieldnames(a.Materials.Fluid);
        a.Option1UI2.Value = find(strcmp(a.Couplant2.Name,a.Option1UI2.String));
        a.Option1UI2.TooltipString = 'Select for which surrounding medium the coincidence angle shall be displayed.';
        a.Option1UI2.Position(3) = 130;
        a.Option1TextUI2.String = 'Couplant';
        a.Option1TextUI2.Position(3) = 44;
        a.XAxisModeTextUI2.String = 'X-axis mode';
        a.XAxisModeTextUI2.Position(3) = 62;
        a.XAxisModeUI2.TooltipString = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode2 == 1
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.XAxis2 = eval(a.XAxisUI2.String);            
        elseif a.XAxisMode2 == 2
            a.XAxisTextUI2.String = 'X-axis (MHz)';
            a.XAxisTextUI2.Position(3) = 65;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3;            
        elseif a.XAxisMode2 == 3
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI2.Position(3) = 86;
            a.XAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;            
        end
        a.YAxisTextUI2.String = ['Y-axis (',char(176),')'];
        a.YAxisTextUI2.Position(3) = 49;
        a.YAxisUI2.TooltipString = 'Enter which coincidence angle range shall be plotted.';

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
    case 5
        a.Quantity12 = 5;
        a.Option1UI2.Enable = 'off';
        a.XAxisModeTextUI2.String = 'X-axis mode';
        a.XAxisModeTextUI2.Position(3) = 62;
        a.XAxisModeUI2.TooltipString = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode2 == 1
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.XAxis2 = eval(a.XAxisUI2.String);            
        elseif a.XAxisMode2 == 2
            a.XAxisTextUI2.String = 'X-axis (MHz)';
            a.XAxisTextUI2.Position(3) = 65;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3;            
        elseif a.XAxisMode2 == 3
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI2.Position(3) = 86;
            a.XAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;             
        end
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
            a.YAxisTextUI2.String = 'Y-axis (mm/mm)';
            a.YAxisTextUI2.Position(3) = 80;
            a.YAxisUI2.TooltipString = 'Enter which wavelength per thickness range shall be plotted.';

            if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                a.YAxisUI2.String = ['[0 ',num2str(x/a.PlateThickness),']'];
            else
                a.YAxisUI2.String = '[0 50]';
            end
            a.YAxis2 = eval(a.YAxisUI2.String);            
        else
            a.YAxisTextUI2.String = 'Y-axis (mm)';
            a.YAxisTextUI2.Position(3) = 61;
            a.YAxisUI2.TooltipString = 'Enter which wavelength range shall be plotted.';

            if  length(a.LayerOrientations) == 1 && strcmp(a.MaterialClasses,'Isotropic')
                a.YAxisUI2.String = ['[0 ',num2str(x),']'];
            else
                a.YAxisUI2.String = ['[0 ',num2str(50*a.PlateThickness),']'];
            end
            a.YAxis2 = eval(a.YAxisUI2.String);
        end
    case 6
        a.Quantity12 = 6;
        a.Option1UI2.Enable = 'off';
        a.XAxisModeTextUI2.String = 'X-axis mode';
        a.XAxisModeTextUI2.Position(3) = 62;
        a.XAxisModeUI2.TooltipString = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode2 == 1
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.XAxis2 = eval(a.XAxisUI2.String);            
        elseif a.XAxisMode2 == 2
            a.XAxisTextUI2.String = 'X-axis (MHz)';
            a.XAxisTextUI2.Position(3) = 65;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3;            
        elseif a.XAxisMode2 == 3
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI2.Position(3) = 86;
            a.XAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;             
        end
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
            a.YAxisTextUI2.String = ['Y-axis (rad/mm',char(8901),'mm)'];
            a.YAxisTextUI2.Position(3) = 101;
            a.YAxisUI2.TooltipString = 'Enter which wavenumber-thickness range shall be plotted.'; 
            
            a.YAxisUI2.String = ['[0 ',num2str(x*a.PlateThickness),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        else
            a.YAxisTextUI2.String = 'Y-axis (rad/mm)';
            a.YAxisTextUI2.Position(3) = 80;
            a.YAxisUI2.TooltipString = 'Enter which wavenumber range shall be plotted.';
            
            a.YAxisUI2.String = ['[0 ',num2str(x),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        end 
    case 7
        a.Quantity12 = 7;
        a.Option1UI2.Enable = 'off';
        a.XAxisModeTextUI2.String = 'X-axis mode';
        a.XAxisModeTextUI2.Position(3) = 62;
        a.XAxisModeUI2.TooltipString = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode2 == 1
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2),']'];
            a.XAxis2 = eval(a.XAxisUI2.String);            
        elseif a.XAxisMode2 == 2
            a.XAxisTextUI2.String = 'X-axis (MHz)';
            a.XAxisTextUI2.Position(3) = 65;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3;            
        elseif a.XAxisMode2 == 3
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI2.Position(3) = 86;
            a.XAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI2.String = ['[0 ',num2str(a.FrequencyLimit2/1e3*a.PlateThickness),']'];
            a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;             
        end
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
            a.YAxisTextUI2.String = ['Y-axis (Np/m',char(8901),'mm)'];
            a.YAxisTextUI2.Position(3) = 90;
            a.YAxisUI2.TooltipString = 'Enter which attenuation-thickness range shall be plotted.'; 
            
            a.YAxisUI2.String = ['[0 ',num2str(x*a.PlateThickness),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        else
            a.YAxisTextUI2.String = 'Y-axis (Np/m)';
            a.YAxisTextUI2.Position(3) = 69;
            a.YAxisUI2.TooltipString = 'Enter which attenuation range shall be plotted.';
            
            a.YAxisUI2.String = ['[0 ',num2str(x),']'];
            a.YAxis2 = eval(a.YAxisUI2.String);
        end
    end
elseif strcmp(source.Tag,'20') % option (Dispersion diagrams)
    if  a.Quantity12 < 3
        a.BulkVelocities2 = source.Value;
    elseif a.Quantity12 == 3
        a.Distance2 = str2double(source.String);
        
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
        E = fieldnames(a.Materials.Fluid);
        switch source.Value
        case source.Value
            a.Couplant2 = getfield(a.Materials.Fluid,E{source.Value}); 
        end

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
    end
elseif strcmp(source.Tag,'21') % X-axis mode
    switch source.Value
    case 1
        if  a.Quantity12 ~= 3
            a.XAxisTextUI2.String = 'X-axis (kHz)';
            a.XAxisTextUI2.Position(3) = 62;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';

            val = eval(a.XAxisUI2.String);
            if  a.XAxisMode2 == 2
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
                a.XAxisUI2.String = ['[',num2str(val(1)*1e3),' ',num2str(val(2)*1e3),']'];
            elseif a.XAxisMode2 == 3
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;
                a.XAxisUI2.String = ['[',num2str(val(1)*1e3/a.PlateThickness),' ',num2str(val(2)*1e3/a.PlateThickness),']'];
            end            
        elseif a.Quantity12 == 3
            a.YAxisTextUI2.String = 'Y-axis (kHz)';
            a.YAxisTextUI2.Position(3) = 63;
            a.YAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            val = eval(a.YAxisUI2.String);
            if  a.XAxisMode2 == 2
                a.YAxis2 = eval(a.YAxisUI2.String)*1e3;
                a.YAxisUI2.String = ['[',num2str(val(1)*1e3),' ',num2str(val(2)*1e3),']'];
            elseif a.XAxisMode2 == 3
                a.YAxis2 = eval(a.YAxisUI2.String)*1e3/a.PlateThickness;
                a.YAxisUI2.String = ['[',num2str(val(1)*1e3/a.PlateThickness),' ',num2str(val(2)*1e3/a.PlateThickness),']'];
            end
        end
        if  (a.Quantity12 == 5 || a.Quantity12 == 6 || a.Quantity12 == 7) && a.XAxisMode2 == 3
            if  a.Quantity12 == 5
                a.YAxisTextUI2.String = 'Y-axis (mm)';
                a.YAxisTextUI2.Position(3) = 61;
                a.YAxisUI2.TooltipString = 'Enter which wavelength range shall be plotted.';
            elseif a.Quantity12 == 6
                a.YAxisTextUI2.String = 'Y-axis (rad/mm)';
                a.YAxisTextUI2.Position(3) = 80;                
                a.YAxisUI2.TooltipString = 'Enter which wavenumber range shall be plotted.';
            elseif a.Quantity12 == 7
                a.YAxisTextUI2.String = 'Y-axis (Np/m)';
                a.YAxisTextUI2.Position(3) = 69;
                a.YAxisUI2.TooltipString = 'Enter which attenuation range shall be plotted.';
            end

            val = eval(a.YAxisUI2.String);
            if  a.Quantity12 == 5
                a.YAxis2 = eval(a.YAxisUI2.String)*a.PlateThickness;
                a.YAxisUI2.String = ['[',num2str(val(1)*a.PlateThickness),' ',num2str(val(2)*a.PlateThickness),']'];
            elseif a.Quantity12 == 6 || a.Quantity12 == 7
                a.YAxis2 = eval(a.YAxisUI2.String)/a.PlateThickness;
                a.YAxisUI2.String = ['[',num2str(val(1)/a.PlateThickness),' ',num2str(val(2)/a.PlateThickness),']'];                
            end
        end        
        
        a.XAxisMode2 = 1;
    case 2
        if  a.Quantity12 ~= 3
            a.XAxisTextUI2.String = 'X-axis (MHz)';
            a.XAxisTextUI2.Position(3) = 65;
            a.XAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            val = eval(a.XAxisUI2.String);
            if  a.XAxisMode2 == 1
                a.XAxis2 = eval(a.XAxisUI2.String);
                a.XAxisUI2.String = ['[',num2str(val(1)/1e3),' ',num2str(val(2)/1e3),']'];
            elseif a.XAxisMode2 == 3
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3/a.PlateThickness;
                a.XAxisUI2.String = ['[',num2str(val(1)/a.PlateThickness),' ',num2str(val(2)/a.PlateThickness),']'];
            end
        elseif a.Quantity12 == 3
            a.YAxisTextUI2.String = 'Y-axis (MHz)';
            a.YAxisTextUI2.Position(3) = 66;
            a.YAxisUI2.TooltipString = 'Enter which frequency range shall be plotted.';
            
            val = eval(a.YAxisUI2.String);
            if  a.XAxisMode2 == 1
                a.YAxis2 = eval(a.YAxisUI2.String);
                a.YAxisUI2.String = ['[',num2str(val(1)/1e3),' ',num2str(val(2)/1e3),']'];
            elseif a.XAxisMode2 == 3
                a.YAxis2 = eval(a.YAxisUI2.String)*1e3/a.PlateThickness;
                a.YAxisUI2.String = ['[',num2str(val(1)/a.PlateThickness),' ',num2str(val(2)/a.PlateThickness),']'];
            end
        end
        if  (a.Quantity12 == 5 || a.Quantity12 == 6 || a.Quantity12 == 7) && a.XAxisMode2 == 3
            if  a.Quantity12 == 5
                a.YAxisTextUI2.String = 'Y-axis (mm)';
                a.YAxisTextUI2.Position(3) = 61;
                a.YAxisUI2.TooltipString = 'Enter which wavelength range shall be plotted.';
            elseif a.Quantity12 == 6
                a.YAxisTextUI2.String = 'Y-axis (rad/mm)';
                a.YAxisTextUI2.Position(3) = 80;                
                a.YAxisUI2.TooltipString = 'Enter which wavenumber range shall be plotted.';
            elseif a.Quantity12 == 7
                a.YAxisTextUI2.String = 'Y-axis (Np/m)';
                a.YAxisTextUI2.Position(3) = 69;
                a.YAxisUI2.TooltipString = 'Enter which attenuation range shall be plotted.';
            end

            val = eval(a.YAxisUI2.String);
            if  a.Quantity12 == 5
                a.YAxis2 = eval(a.YAxisUI2.String)*a.PlateThickness;
                a.YAxisUI2.String = ['[',num2str(val(1)*a.PlateThickness),' ',num2str(val(2)*a.PlateThickness),']'];
            elseif a.Quantity12 == 6 || a.Quantity12 == 7
                a.YAxis2 = eval(a.YAxisUI2.String)/a.PlateThickness;
                a.YAxisUI2.String = ['[',num2str(val(1)/a.PlateThickness),' ',num2str(val(2)/a.PlateThickness),']'];                
            end
        end
        
        a.XAxisMode2 = 2;
    case 3
        if  a.Quantity12 ~= 3
            a.XAxisTextUI2.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI2.Position(3) = 86;
            a.XAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            val = eval(a.XAxisUI2.String);
            if  a.XAxisMode2 == 1
                a.XAxis2 = eval(a.XAxisUI2.String);
                a.XAxisUI2.String = ['[',num2str(val(1)/1e3*a.PlateThickness),' ',num2str(val(2)/1e3*a.PlateThickness),']'];
            elseif a.XAxisMode2 == 2
                a.XAxis2 = eval(a.XAxisUI2.String)*1e3;
                a.XAxisUI2.String = ['[',num2str(val(1)*a.PlateThickness),' ',num2str(val(2)*a.PlateThickness),']'];
            end            
        elseif a.Quantity12 == 3
            a.YAxisTextUI2.String = ['Y-axis (MHz',char(8901),'mm)'];
            a.YAxisTextUI2.Position(3) = 87;
            a.YAxisUI2.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            val = eval(a.YAxisUI2.String);
            if  a.XAxisMode2 == 1
                a.YAxis2 = eval(a.YAxisUI2.String);
                a.YAxisUI2.String = ['[',num2str(val(1)/1e3*a.PlateThickness),' ',num2str(val(2)/1e3*a.PlateThickness),']'];
            elseif a.XAxisMode2 == 2
                a.YAxis2 = eval(a.YAxisUI2.String)*1e3;
                a.YAxisUI2.String = ['[',num2str(val(1)*a.PlateThickness),' ',num2str(val(2)*a.PlateThickness),']'];
            end
        end
        if  (a.Quantity12 == 5 || a.Quantity12 == 6 || a.Quantity12 == 7) && a.XAxisMode2 ~= 3 
            if  a.Quantity12 == 5
                a.YAxisTextUI2.String = 'Y-axis (mm/mm)';
                a.YAxisTextUI2.Position(3) = 80;                
                a.YAxisUI2.TooltipString = 'Enter which wavelength per thickness range shall be plotted.';
            elseif a.Quantity12 == 6
                a.YAxisTextUI2.String = ['Y-axis (rad/mm',char(8901),'mm)'];
                a.YAxisTextUI2.Position(3) = 101; 
                a.YAxisUI2.TooltipString = 'Enter which wavenumber-thickness range shall be plotted.';
            elseif a.Quantity12 == 7
                a.YAxisTextUI2.String = ['Y-axis (Np/m',char(8901),'mm)'];
                a.YAxisTextUI2.Position(3) = 90;
                a.YAxisUI2.TooltipString = 'Enter which attenuation-thickness range shall be plotted.';
            end

            val = eval(a.YAxisUI2.String);
            if  a.Quantity12 == 5
                a.YAxis2 = eval(a.YAxisUI2.String)/a.PlateThickness;
                a.YAxisUI2.String = ['[',num2str(val(1)/a.PlateThickness),' ',num2str(val(2)/a.PlateThickness),']'];
            elseif a.Quantity12 == 6 || a.Quantity12 == 7
                a.YAxis2 = eval(a.YAxisUI2.String)*a.PlateThickness;
                a.YAxisUI2.String = ['[',num2str(val(1)*a.PlateThickness),' ',num2str(val(2)*a.PlateThickness),']'];                
            end
        end        

        a.XAxisMode2 = 3;
    end
elseif strcmp(source.Tag,'22') % X-axis
    a.XAxis2 = eval(source.String);
    
    if  a.Quantity12 ~= 3
        if  a.XAxisMode2 == 1
            a.XAxis2 = eval(source.String);
        elseif a.XAxisMode2 == 2
            a.XAxis2 = eval(source.String)*1e3;
        elseif a.XAxisMode2 == 3
            a.XAxis2 = eval(source.String)*1e3/a.PlateThickness;
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
            a.YAxis2 = eval(source.String)*1e3/a.PlateThickness;
        end          
    end
elseif strcmp(source.Tag,'24') % Plot (Dispersion diagrams)
    if  a.Quantity12 == 1
        PhaseVelocity_Anisotropic(a.Hybrid,a.CropPlots2,a.LayerOrientations-a.PropagationAngle,a.LayupString1,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.BulkVelocities2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BoxLineWidth2,a.BShear,a.BScholte2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.Title2,a.HigherOrderModes2,a.LambModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.FileName2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.ShearHorizontalModes2,a.ScholteModes2,a.SLamb2,a.SShear2,a.SScholte2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled) 
    elseif a.Quantity12 == 2
        if  a.Decoupled == 0
            EnergyVelocitySkewAngle_Anisotropic(a.Hybrid,a.CropPlots2,a.LayupString1,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.AScholte2,a.AntisymmetricModes2,a.B,a.BScholte2,a.BoxLineWidth2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.FileName2,a.Title2,a.HigherOrderModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.SScholte2,a.ScholteModes2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled) 
            EnergyVelocityAbsolut_Anisotropic(a.Hybrid,a.CropPlots2,a.LayupString1,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.AScholte2,a.AntisymmetricModes2,a.B,a.BScholte2,a.BoxLineWidth2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.FileName2,a.Title2,a.HigherOrderModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.SScholte2,a.ScholteModes2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)            
            EnergyVelocity2_Anisotropic(a.Hybrid,a.CropPlots2,a.LayerOrientations-a.PropagationAngle,a.LayupString1,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.BulkVelocities2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.AScholte2,a.AntisymmetricModes2,a.B,a.BScholte2,a.BoxLineWidth2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.FileName2,a.Title2,a.HigherOrderModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.SScholte2,a.ScholteModes2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
        end
        EnergyVelocity1_Anisotropic(a.Hybrid,a.CropPlots2,a.LayerOrientations-a.PropagationAngle,a.LayupString1,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.BulkVelocities2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BoxLineWidth2,a.BShear,a.BScholte2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.FileName2,a.Title2,a.HigherOrderModes2,a.LambModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.ShearHorizontalModes2,a.ScholteModes2,a.SLamb2,a.SShear2,a.SScholte2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
    elseif a.Quantity12 == 3
        PropagationTime_Anisotropic(a.Hybrid,a.CropPlots2,a.LayupString1,a.FluidLoading2,a.LowerFluid2,a.UpperFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.Distance2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BoxLineWidth2,a.BShear,a.BScholte2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.FileName2,a.Title2,a.HigherOrderModes2,a.LambModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.ShearHorizontalModes2,a.ScholteModes2,a.SLamb2,a.SShear2,a.SScholte2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
    elseif a.Quantity12 == 4
        CoincidenceAngle_Anisotropic(a.Hybrid,a.CropPlots2,a.LayupString1,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BoxLineWidth2,a.BShear,a.BScholte2,a.Couplant2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.Title2,a.HigherOrderModes2,a.LambModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.FileName2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.ShearHorizontalModes2,a.ScholteModes2,a.SLamb2,a.SShear2,a.SScholte2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)         
    elseif a.Quantity12 == 5
        Wavelength_Anisotropic(a.Hybrid,a.CropPlots2,a.LayupString1,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BoxLineWidth2,a.BShear,a.BScholte2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.Title2,a.HigherOrderModes2,a.LambModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.FileName2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.ShearHorizontalModes2,a.ScholteModes2,a.SLamb2,a.SShear2,a.SScholte2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)         
    elseif a.Quantity12 == 6
        Wavenumber_Anisotropic(a.Hybrid,a.CropPlots2,a.LayupString1,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BoxLineWidth2,a.BShear,a.BScholte2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.Title2,a.HigherOrderModes2,a.LambModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.FileName2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.ShearHorizontalModes2,a.ScholteModes2,a.SLamb2,a.SShear2,a.SScholte2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
    elseif a.Quantity12 == 7
        Attenuation_Anisotropic(a.Hybrid,a.CropPlots2,a.LayupString1,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.PNGresolution2,a.SColor2,a.AColor2,a.BColor,a.A,a.ALamb2,a.AntisymmetricModes2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BoxLineWidth2,a.BShear,a.BScholte2,a.Directory2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.ModeLabelFontSize2,a.Title2,a.HigherOrderModes2,a.LambModes2,a.LineWidth2,a.Material2,a.ModeLabels2,a.S0B1Label,a.SH0Label2,a.A0B0Label,a.PDF2,a.FileName2,a.PlateThickness/1e3,a.PNG2,a.PropagationAngle,a.S,a.ShearHorizontalModes2,a.ScholteModes2,a.SLamb2,a.SShear2,a.SScholte2,a.SuperLayers,a.SuperLayerSize,a.SymmetricModes2,a.SymmetricSystem,a.Symmetric,a.XAxis2,a.XAxisMode2,a.YAxis2,a.Decoupled)
    end
elseif strcmp(source.Tag,'25') % Quantity (Through thickness profiles)
    switch source.Value
    case 1
        a.Quantity22 = 1;
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
        a.x13UI2.TooltipString = 'Displacement u1';
        a.x23UI2.TooltipString = 'Displacement u2';
        a.x33UI2.TooltipString = 'Displacement u3';        
    case 2
        a.Quantity22 = 2;
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
        a.x11UI2.TooltipString = ['Stress ',char(963),'11'];
        a.x12UI2.TooltipString = ['Stress ',char(963),'12'];
        a.x13UI2.TooltipString = ['Stress ',char(963),'13'];
        a.x22UI2.TooltipString = ['Stress ',char(963),'22'];
        a.x23UI2.TooltipString = ['Stress ',char(963),'23'];
        a.x33UI2.TooltipString = ['Stress ',char(963),'33'];        
    case 3
        a.Quantity22 = 3;
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
        a.x11UI2.TooltipString = ['Strain ',char(949),'11'];
        a.x12UI2.TooltipString = ['Strain ',char(949),'12'];
        a.x13UI2.TooltipString = ['Strain ',char(949),'13'];
        a.x22UI2.TooltipString = ['Strain ',char(949),'22'];
        a.x23UI2.TooltipString = ['Strain ',char(949),'23'];
        a.x33UI2.TooltipString = ['Strain ',char(949),'33'];
    case 4
        a.Quantity22 = 4;
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
        a.x13UI2.TooltipString = 'Strain energy density';
        a.x23UI2.TooltipString = 'Kinetic energy density';
        a.x33UI2.TooltipString = 'Total energy density';        
    case 5
        a.Quantity22 = 5;
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
        a.x13UI2.TooltipString = 'Power flow density p1';
        a.x23UI2.TooltipString = 'Power flow density p2';
        a.x33UI2.TooltipString = 'Power flow density p3';        
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
    if  a.Quantity22 == 1
        u1Color = [1 0 0];
        u2Color =  [.13 .55 .13];
        u3Color = [0 0 1];
        ModeShapeLines_Anisotropic(1,1,0,0,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,u1Color,u2Color,u3Color,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2) 
    elseif a.Quantity22 == 2 
        sigma33Color = [0 0 1];
        sigma11Color = [1 0 0];
        sigma22Color = [.13 .55 .13];
        sigma13Color = [0 0 0];
        sigma23Color = [1 0 1];
        sigma12Color = [0 1 1];              
        ModeShapeLines_Anisotropic(2,1,0,0,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,sigma11Color,sigma22Color,sigma33Color,sigma23Color,sigma13Color,sigma12Color,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
    elseif a.Quantity22 == 3
        epsilon33Color = [0 0 1];
        epsilon11Color = [1 0 0];
        epsilon22Color = [.13 .55 .13];
        epsilon13Color = [0 0 0];
        epsilon23Color = [1 0 1];
        epsilon12Color = [0 1 1];
        ModeShapeLines_Anisotropic(3,1,0,0,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,epsilon11Color,epsilon22Color,epsilon33Color,epsilon23Color,epsilon13Color,epsilon12Color,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
    elseif a.Quantity22 == 4
        StrainEnergyDensityColor = [0 0 1];
        KineticEnergyDensityColor = [1 0 0];
        TotalEnergyDensityColor = [0 0 0];        
        ModeShapeLines_Anisotropic(4,1,0,0,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,StrainEnergyDensityColor,KineticEnergyDensityColor,TotalEnergyDensityColor,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,0)
    elseif a.Quantity22 == 5
        P1Color = [1 0 0];
        P2Color = [.13 .55 .13];
        P3Color = [0 0 1];
        ModeShapeLines_Anisotropic(5,1,0,0,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,P1Color,P2Color,P3Color,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,0)
    end
elseif strcmp(source.Tag,'36') % Mode (Mode shape)
    a.Mode22 = a.ModeNames2{source.Value};
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'37') % Frequency (Mode shape) 
    a.Frequency22 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'38') % Wavelengths
    a.Length2 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'39') % Samples x1
    a.Samples22 = str2double(source.String);
     
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
    a.GridLine2 = 2*round(a.Samples22/80);
    a.GridLineUI2.String = a.GridLine2;

    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'40') % Samples per layer (Mode shape) 
    a.Samples32 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'41') % Scale
    a.Scale2 = str2double(source.String);
elseif strcmp(source.Tag,'42') % Grid line
    a.GridLine2 = str2double(source.String);
elseif strcmp(source.Tag,'43') % Undistorted
    a.Undistorted2 = source.Value;
elseif strcmp(source.Tag,'6') % Half-spaces (mode shape)
    a.HalfspacesNumber22 = str2double(source.String);
    
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
    
    a.ModeShapeSettingChanged2 = 1;
elseif strcmp(source.Tag,'8') % Half-spaces (mode shape)
    a.Halfspaces22 = source.Value;
    
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
    
    a.ModeShapeSettingChanged2 = 1;   
elseif strcmp(source.Tag,'45') % Cycles  
    a.Cycles2 = str2double(source.String);
    a.ModeShapeSettingChanged2 = 1;
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
    if  a.Animate2 == 0
        ModeShapeImage_Anisotropic(a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,str2double(a.LayerCountUI2.String),a.Material2,a.c,a.PNGresolution2,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.Directory2,a.ExportPlots2,a.TitleFontSize2,a.AxesLabelFontSize2,a.Frequency22,a.GridLine2,a.Title2,a.Length2,a.LineWidth2,a.Mode22,a.FileName2,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples22,a.Samples32,a.Scale2,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Undistorted2,a.Decoupled,a.Halfspaces22,a.HalfspacesNumber22)
    else
        if  a.ModeShapeSettingChanged2 == 1
            [a.Time2,a.u2,a.x12,a.x32,a.x3Total,a.p2] = ModeShapeAnimationComputer_Anisotropic(a.Material2,str2double(a.LayerCountUI2.String),a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.c,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.CycleDuration2,a.Cycles2,a.FrameRate2,a.Frequency22,a.Length2,a.Mode22,a.PlateThickness/1e3,a.Samples22,a.Samples32,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces22,a.HalfspacesNumber22);
            a.ModeShapeSettingChanged2 = 0;
        end
        ModeShapeAnimation_Anisotropic(a.Hybrid,str2double(a.LayerCountUI2.String),a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.LayupString1,a.Directory2,a.TitleFontSize2,a.AxesLabelFontSize2,a.FrameRate2,a.Frequency22,a.GridLine2,a.Title2,a.LineWidth2,a.Material2,a.Mode22,a.ExportPlots2,a.FileName2,a.MovieQuality2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples22,a.Samples32,a.Scale2,a.SuperLayers,a.SuperLayerSize,a.SymmetricSystem,a.Time2,a.u2,a.Undistorted2,a.x12,a.x32,a.x3Total,a.p2,a.Halfspaces22,a.HalfspacesNumber22)
    end
elseif strcmp(source.Tag,'50') % Export plots
    a.ExportPlots2 = source.Value;

    if  source.Value == 1
        a.Plot1UI2.String = 'Export';
        a.Plot2UI2.String = 'Export';
        a.Plot3UI2.String = 'Export';
    else
        a.Plot1UI2.String = 'Plot';
        a.Plot2UI2.String = 'Plot';
        a.Plot3UI2.String = 'Plot';
    end
elseif strcmp(source.Tag,'79') % Crop plots
    a.CropPlots2 = source.Value;
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
    switch source.Value
    case 1
        a.XAxisMode2UI2.Enable = 'on';
        a.ArrangeUI2.Enable = 'on';
    case 0
        a.XAxisMode2UI2.Enable = 'off';
        a.ArrangeUI2.Enable = 'off';
    end
elseif strcmp(source.Tag,'57') % Through thickness
    a.ThroughThickness2 = source.Value;
elseif strcmp(source.Tag,'80') % Matlab
    if  a.DispersionCurves2 == 1
        Export_Anisotropic(a,0,0,1)
    end
    if  a.ThroughThickness2 == 1
        if  a.Quantity22 == 1
            ModeShapeLines_Anisotropic(1,0,1,0,0,1,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 2
            ModeShapeLines_Anisotropic(2,0,1,0,0,1,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 3
            ModeShapeLines_Anisotropic(3,0,1,0,0,1,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 4
            ModeShapeLines_Anisotropic(4,0,1,0,0,1,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,0)
        elseif a.Quantity22 == 5
            ModeShapeLines_Anisotropic(5,0,1,0,0,1,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,0)
        end
    end
elseif strcmp(source.Tag,'58') % Excel
    if  a.DispersionCurves2 == 1
        Export_Anisotropic(a,1,0,0)
    end
    if  a.ThroughThickness2 == 1
        if  a.Quantity22 == 1
            ModeShapeLines_Anisotropic(1,0,1,1,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 2
            ModeShapeLines_Anisotropic(2,0,1,1,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 3
            ModeShapeLines_Anisotropic(3,0,1,1,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 4
            ModeShapeLines_Anisotropic(4,0,1,1,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,0)
        elseif a.Quantity22 == 5
            ModeShapeLines_Anisotropic(5,0,1,1,0,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,0)
        end
    end
elseif strcmp(source.Tag,'59') % TXT
    if  a.DispersionCurves2 == 1
        Export_Anisotropic(a,0,1,0)
    end
    if  a.ThroughThickness2 == 1
        if  a.Quantity22 == 1
            ModeShapeLines_Anisotropic(1,0,1,0,1,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 2
            ModeShapeLines_Anisotropic(2,0,1,0,1,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 3
            ModeShapeLines_Anisotropic(3,0,1,0,1,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,a.Phase2)
        elseif a.Quantity22 == 4
            ModeShapeLines_Anisotropic(4,0,1,0,1,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,0)
        elseif a.Quantity22 == 5
            ModeShapeLines_Anisotropic(5,0,1,0,1,0,a.Hybrid,a.Viscoelastic2,a.FluidLoading2,a.UpperFluid2,a.LowerFluid2,a.ToggleUpperFluid2,a.ToggleLowerFluid2,a.CropPlots2,a.LayupString1,a.Plot2,a.PNGresolution2,0,0,0,0,0,0,a.A,a.ALamb2,a.AShear2,a.AScholte2,a.B,a.BLamb,a.BShear,a.BScholte2,a.S,a.SLamb2,a.SShear2,a.SScholte2,a.BoxLineWidth2,a.c,a.Material2,a.Directory2,a.FileName2,a.ExportPlots2,a.AxesTickFontSize2,a.AxesLabelFontSize2,a.TitleFontSize2,a.LegendFontSize2,a.Frequency12,a.Title2,a.LegendLocation2,a.LineWidth2,a.Mode12,a.PDF2,a.PNG2,a.PlateThickness/1e3,a.PropagationAngle,a.Samples12,a.SuperLayers,a.SuperLayerSize,a.LayerThicknesses/1e3,a.SymmetricSystem,a.Decoupled,a.Halfspaces12,a.HalfspacesNumber12,0)
        end
    end
elseif strcmp(source.Tag,'60') % File name
    a.FileName2 = source.String;
elseif strcmp(source.Tag,'61') % Directory
    a.Directory2 = source.String;
elseif strcmp(source.Tag,'62') % Title
    switch source.Value
    case 1
       a.Title2 = 2;
    case 2
       a.Title2 = 1;
    case 3
       a.Title2 = 0;
    end
elseif strcmp(source.Tag,'63') % Mode labels
    a.ModeLabels2 = source.Value;
elseif strcmp(source.Tag,'64') % Legend location
    switch source.Value
    case 1
       a.LegendLocation2 = 'out';
    case 2
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
    a.BColor = eval(source.String);
elseif strcmp(source.Tag,'70') % S0/B1
    a.S0B1Label = str2double(source.String);
elseif strcmp(source.Tag,'71') % S'0/B'0
    a.SH0Label2 = str2double(source.String);
elseif strcmp(source.Tag,'72') % A0/B0
    a.A0B0Label = str2double(source.String);
elseif strcmp(source.Tag,'73') % Title (Font size)
    a.TitleFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'74') % Axes labels
    a.AxesLabelFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'75') % Axes ticks
    a.AxesTickFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'76') % Mode labels
    a.ModeLabelFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'77') % Legend
    a.LegendFontSize2 = str2double(source.String);
elseif strcmp(source.Tag,'78') % Default
    a.Title2 = 2;
    a.TitleUI2.Value = 1;
    a.ModeLabels2 = 1;
    a.ModeLabelsUI2.Value = a.ModeLabels2;
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
    a.BColor = [.5 0 1];
    a.BColorUI2.String = '[.5 0 1]';    
    a.S0B1Label = .05;
    a.S0B1LabelUI2.String = a.S0B1Label;
    a.SH0Label2 = .05;
    a.SH0LabelUI2.String = a.SH0Label2;
    a.A0B0Label = .05;
    a.A0B0LabelUI2.String = a.A0B0Label;
    a.TitleFontSize2 = 30;
    a.TitleFontSizeUI2.String = a.TitleFontSize2;
    a.AxesLabelFontSize2 = 30;
    a.AxesLabelFontSizeUI2.String = a.AxesLabelFontSize2;
    a.AxesTickFontSize2 = 24;
    a.AxesTickFontSizeUI2.String = a.AxesTickFontSize2;
    a.ModeLabelFontSize2 = 24;
    a.ModeLabelFontSizeUI2.String = a.ModeLabelFontSize2;
    a.LegendFontSize2 = 24;
    a.LegendFontSizeUI2.String = a.LegendFontSize2;
end