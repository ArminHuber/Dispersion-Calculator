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
function a = CallbackModule_Isotropic(source,~,a,Tab1,Tab3)
%#ok<*AGROW>
%#ok<*GFLD>
%#ok<*GVMIS> 
global Stop 
Stop = 0;
if  strcmp(source.Tag,'6') % Fluid loading
    a.FluidLoading1 = source.Value;
    if  source.Value == 1
        a.FluidUI1.Enable = 'on';
        a.ScholteModesUI1.Enable = 'on';
        a.HalfspacesNumber1UI1.Enable = 'on';
        a.Halfspaces1UI1.Enable = 'on';
        a.HalfspacesNumber2UI1.Enable = 'on';
        a.Halfspaces2UI1.Enable = 'on';
        
        a.Couplant1 = a.Fluid1;
        if  a.Quantity11 == 4
            a.Option1UI1.Value = find(strcmp(a.Couplant1.Name,a.Option1UI1.String));

            x = round(8e4/a.Material1.PlateVelocity*a.Couplant1.Velocity/343);
            if  x > 90
                x = 90;
            end
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        end
    else
        a.FluidUI1.Enable = 'off';
        a.ScholteModesUI1.Enable = 'off';
        a.HalfspacesNumber1UI1.Enable = 'off';
        a.Halfspaces1UI1.Enable = 'off';
        a.HalfspacesNumber2UI1.Enable = 'off';
        a.Halfspaces2UI1.Enable = 'off';
    end
    if  a.Quantity11 == 7
        if  a.FluidLoading1 == 1 && a.Viscoelastic1 == 1
            x = 1e4*(a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)+a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
        elseif a.FluidLoading1 == 1 && a.Viscoelastic1 == 0
            x = 1e4*a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)/a.Thickness;
        elseif a.FluidLoading1 == 0 && a.Viscoelastic1 == 1
            x = 1e4*(a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
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
        if  a.XAxisMode1 == 3
            a.YAxisUI1.String = ['[0 ',num2str(x*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        end
    end
elseif strcmp(source.Tag,'7') % Fluid
    a.Fluid1 = getfield(a.Materials.Fluid,cell2mat(source.String(source.Value)));
    
    a.Couplant1 = a.Fluid1;
    if  a.Quantity11 == 4
        a.Option1UI1.Value = find(strcmp(a.Couplant1.Name,a.Option1UI1.String));
        
        x = round(8e4/a.Material1.PlateVelocity*a.Couplant1.Velocity/343);
        if  x > 90
            x = 90;
        end
        a.YAxisUI1.String = ['[0 ',num2str(x),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    elseif a.Quantity11 == 7
        if  a.Viscoelastic1 == 1
            x = 1e4*(a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)+a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
        elseif a.Viscoelastic1 == 0
            x = 1e4*a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)/a.Thickness;
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
            a.YAxisUI1.String = ['[0 ',num2str(x*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        end
    end    
elseif strcmp(source.Tag,'1') % Material
    a.Material1 = getfield(a.Materials.Isotropic,cell2mat(source.String(source.Value)));

    if  isreal(a.Material1.C)
        a.Viscoelastic1 = 0;
    else
        a.Viscoelastic1 = 1;
    end

    a.PhaseVelocityLimit1 = round(a.YRange*a.Material1.PlateVelocity,-3);
    a.PhaseVelocityLimitUI1.String = a.PhaseVelocityLimit1/1e3;
    
    a.FrequencyLimit1 = round(a.Material1.PlateVelocity*a.XRange1/a.Thickness,-2);
    if  a.FrequencyLimit1 == 0
        a.FrequencyLimit1 = 50;
    elseif a.FrequencyLimit1 >= 1e4
        a.FrequencyLimit1 = round(a.FrequencyLimit1,-3);
    end
    a.FrequencyLimitUI1.String = a.FrequencyLimit1;

    a.FrequencyResolution1 = a.FrequencyLimit1/a.XSamples1;
    a.FrequencyResolutionUI1.String = a.FrequencyResolution1;
    
    a.Step1 = a.FrequencyLimit1/a.Steps1;
    a.SteptUI1.String = a.Step1;
    
    a.Frequency11 = a.FrequencyLimit1;
    a.Frequency1UI1.String = a.FrequencyLimit1;
    
    a.Frequency21 = a.FrequencyLimit1;
    a.Frequency2UI1.String = a.FrequencyLimit1;

    if  a.Quantity11 == 1
        a.YAxisUI1.String = ['[0 ',num2str(a.PhaseVelocityLimit1/1e3),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);

        if  a.XAxisMode1 == 1
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;                    
        end
    elseif a.Quantity11 == 2
        a.YAxisUI1.String = ['[0 ',num2str(ceil(a.Material1.PlateVelocity/1e3)),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);

        if  a.XAxisMode1 == 1
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;                    
        end
    elseif a.Quantity11 == 3
        if  a.XAxisMode1 == 1
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness;                    
        end
        
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

        if  a.XAxisMode1 == 1
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;                    
        end
    elseif a.Quantity11 == 5
        x = a.YRange*1e-3*a.Material1.PlateVelocity*a.Thickness;
        if  x >= 1e3 && x < 1e4
            x = round(x,-2);
        elseif x >= 1e2 && x < 1e3
            x = round(x,-1);
        elseif x < 1e2
            x = round(x);
        end
        if  a.XAxisMode1 == 3
            a.YAxisUI1.String = ['[0 ',num2str(x/a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        end
        
        if  a.XAxisMode1 == 1
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;                    
        end
    elseif a.Quantity11 == 6 || a.Quantity11 == 7
        if  a.Quantity11 == 6 
            x = 2*pi*a.FrequencyLimit1/a.Material1.RayleighVelocity;
        elseif a.Quantity11 == 7
            if  a.FluidLoading1 == 1 && a.Viscoelastic1 == 1
                x = 1e4*(a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)+a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
            elseif a.FluidLoading1 == 1 && a.Viscoelastic1 == 0
                x = 1e4*a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)/a.Thickness;
            elseif a.FluidLoading1 == 0 && a.Viscoelastic1 == 1
                x = 1e4*(a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
            else
                x = 1;
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
            a.YAxisUI1.String = ['[0 ',num2str(x*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        end
        
        if  a.XAxisMode1 == 1
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;                    
        end
    end
    
    a.Detect1 = 0;
elseif strcmp(source.Tag,'2') % Thickness
    a.Thickness = str2double(source.String);
    
    a.FrequencyLimit1 = round(a.Material1.PlateVelocity*a.XRange1/a.Thickness,-2);
    if  a.FrequencyLimit1 == 0
        a.FrequencyLimit1 = 50;
    elseif a.FrequencyLimit1 >= 1e4
        a.FrequencyLimit1 = round(a.FrequencyLimit1,-3);
    end
    a.FrequencyLimitUI1.String = a.FrequencyLimit1;
        
    a.Step1 = a.FrequencyLimit1/a.Steps1;
    a.SteptUI1.String = a.Step1;    

    if  a.Quantity11 == 3
        if  a.XAxisMode1 == 1
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness;                    
        end 
    else
        if  a.XAxisMode1 == 1
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;                    
        end       
    end
    
    if  a.Quantity11 == 5
        x = a.YRange*1e-3*a.Material1.PlateVelocity*a.Thickness;
        if  x >= 1e3 && x < 1e4
            x = round(x,-2);
        elseif x >= 1e2 && x < 1e3
            x = round(x,-1);
        elseif x < 1e2
            x = round(x);
        end        
        if  a.XAxisMode1 == 3
            a.YAxisUI1.String = ['[0 ',num2str(x/a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        end
    elseif a.Quantity11 == 6 || a.Quantity11 == 7
        if  a.Quantity11 == 6 
            x = 2*pi*a.FrequencyLimit1/a.Material1.RayleighVelocity;           
        elseif a.Quantity11 == 7
            if  a.FluidLoading1 == 1 && a.Viscoelastic1 == 1
                x = 1e4*(a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)+a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
            elseif a.FluidLoading1 == 1 && a.Viscoelastic1 == 0
                x = 1e4*a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)/a.Thickness;
            elseif a.FluidLoading1 == 0 && a.Viscoelastic1 == 1
                x = 1e4*(a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
            else
                x = 1;
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
            a.YAxisUI1.String = ['[0 ',num2str(x*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        end        
    end
    
    a.FrequencyResolution1 = a.FrequencyLimit1/a.XSamples1;
    a.FrequencyResolutionUI1.String = a.FrequencyResolution1;
    
    a.Frequency11 = a.FrequencyLimit1;
    a.Frequency1UI1.String = a.FrequencyLimit1;
    
    a.Frequency21 = a.FrequencyLimit1;
    a.Frequency2UI1.String = a.FrequencyLimit1;
    
    a.Detect1 = 0;
elseif strcmp(source.Tag,'3') % Phase velocity limit
    a.PhaseVelocityLimit1 = str2double(source.String)*1e3;
    
    if  a.Quantity11 == 1
        a.YAxisUI1.String = ['[0 ',source.String,']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    end
    
    a.Detect1 = 0;
elseif strcmp(source.Tag,'4') % Frequency limit
    a.FrequencyLimit1 = str2double(source.String);

    if  a.Quantity11 == 3
        if  a.XAxisMode1 == 1
            a.YAxisUI1.String = ['[0 ',source.String,']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.YAxisUI1.String = ['[0 ',num2str(eval(source.String)/1e3),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.YAxisUI1.String = ['[0 ',num2str(eval(source.String)/1e3*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness;                    
        end 
    else
        if  a.XAxisMode1 == 1
            a.XAxisUI1.String = ['[0 ',source.String,']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisUI1.String = ['[0 ',num2str(eval(source.String)/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.XAxisUI1.String = ['[0 ',num2str(eval(source.String)/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;                    
        end       
    end
    if  a.Quantity11 == 6
        x = 2*pi*a.FrequencyLimit1/a.Material1.RayleighVelocity;
        if x >= 1e3 && x < 1e4
            x = 1e2*ceil(x/1e2);
        elseif x >= 1e2 && x < 1e3
            x = 1e1*ceil(x/1e1);
        elseif x >= 1e1 && x < 1e2
            x = ceil(x);
        elseif x < 1e1
            x = 1e-1*ceil(x/1e-1);
        end
        if  a.XAxisMode1 == 3
            a.YAxisUI1.String = ['[0 ',num2str(x*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        else
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        end
    end
    
    a.Frequency11 = a.FrequencyLimit1;
    a.Frequency1UI1.String = a.FrequencyLimit1;
    
    a.Frequency21 = a.FrequencyLimit1;
    a.Frequency2UI1.String = a.FrequencyLimit1;
    
    a.Detect1 = 0;
elseif strcmp(source.Tag,'5') % Frequency step
    a.FrequencyResolution1 = str2double(source.String);
elseif strcmp(source.Tag,'8') % Higher order modes
    a.HigherOrderModes1 = source.Value;
elseif strcmp(source.Tag,'9') % Symmetric modes
    a.SymmetricModes1 = source.Value;
elseif strcmp(source.Tag,'10') % Antisymmetric modes
    a.AntisymmetricModes1 = source.Value;
elseif strcmp(source.Tag,'11') % Lamb modes
    a.LambModes1 = source.Value;
elseif strcmp(source.Tag,'12') % Shear horizontal modes
    a.ShearHorizontalModes1 = source.Value;
elseif strcmp(source.Tag,'79') % Scholte modes
    a.ScholteModes1 = source.Value;
elseif strcmp(source.Tag,'13') % Step
    a.Step1 = str2double(source.String);
    a.Detect1 = 0;
elseif strcmp(source.Tag,'14') % Detect
    [a.HSLamb1,a.HALamb1,a.HSShear1,a.HAShear1] = FrequencySweeper_Isotropic(1:a.Step1:a.FrequencyLimit1,a.Material1,a.PhaseVelocityLimit1,a.Thickness/2e3,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
    if  a.FluidLoading1 == 1 && a.Material1.TransverseVelocity < a.Fluid1.Velocity
        [a.HSScholte1,a.HAScholte1] = FrequencySweeper_Isotropic_Scholte(a.Material1,a.Thickness/2e3,a.Fluid1,1:a.Step1:a.FrequencyLimit1);
    else
        a.HSScholte1 = [];
        a.HAScholte1 = [];
    end
    a.Detect1 = 1;
elseif strcmp(source.Tag,'15') % Calculate
    if  a.Detect1 == 0 && a.HigherOrderModes1 == 1 && (a.LambModes1 == 1 || a.ShearHorizontalModes1 == 1)
        [a.HSLamb1,a.HALamb1,a.HSShear1,a.HAShear1] = FrequencySweeper_Isotropic(1:a.Step1:a.FrequencyLimit1,a.Material1,a.PhaseVelocityLimit1,a.Thickness/2e3,a.OutputWindow1aUI1,a.OutputWindow1bUI1,a.OutputWindow2aUI1,a.OutputWindow2bUI1);
        a.Detect1 = 1;
    end
    if  a.FluidLoading1 == 1 && a.HigherOrderModes1 == 1 && a.ScholteModes1 == 1 && a.Material1.TransverseVelocity < a.Fluid1.Velocity
        [a.HSScholte1,a.HAScholte1] = FrequencySweeper_Isotropic_Scholte(a.Material1,a.Thickness/2e3,a.Fluid1,1:a.Step1:a.FrequencyLimit1);
    else
        a.HSScholte1 = [];
        a.HAScholte1 = [];
    end
    FrequencyRange = 0:a.FrequencyResolution1:a.FrequencyLimit1;
    FrequencyRangeF = FrequencyRange;
    FrequencyRange(1) = a.FrequencyRangeStart1;
    FrequencyRangeF(1) = .1*a.FrequencyResolution1;
    if  a.FluidLoading1 == 1 || a.Viscoelastic1 == 1
        [FALambF,FAScholte] = PhaseVelocitySweeper_Isotropic_F(a.Material1,a.Thickness/2e3,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,a.FluidDensityThreshold1,FrequencyRangeF(1));
    else
        FALambF = [];
        FAScholte = [];
    end
    if  a.Multithreading == 1 && a.LambModes1 == 1 && a.SymmetricModes1 == 1 && a.AntisymmetricModes1 == 1
        delete(a.a1UI1)
        a.a1UI1 = axes('Parent',Tab1,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading22.png','Parent',a.a1UI1)
    end
    [a.ALamb1,a.AShear1,a.AScholte1,a.SLamb1,a.SShear1,a.SScholte1] = Computer_Isotropic(a.Multithreading,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,a.FrequencyLimit1,a.Material1,a.Thickness/2e3,a.PhaseVelocityStep1,a.FrequencyOffset1,a.LambPhaseVelocitySweepRange11,a.LambPhaseVelocitySweepRange21,a.SearchWidthReal1,a.SearchWidthImag1,a.SearchAreaSections1,a.SearchAreaExtensions1,FrequencyRange,FrequencyRangeF,FALambF,FAScholte,a.HSLamb1,a.HSShear1,a.HSScholte1,a.HALamb1,a.HAShear1,a.HAScholte1,a.FrequencyResolution1,a.PhaseVelocityLimit1,a.PhaseVelocityResolution1,a.PhaseVelocitySections1,a.FrequencySections1,a.HigherOrderModes1,a.SymmetricModes1,a.AntisymmetricModes1,a.LambModes1,a.ShearHorizontalModes1,a.ScholteModes1,a.MissingSamples1,a.BelowCutoffWidth1);
    if  Stop == 1
        return
    end
    [a.ALamb1,a.AShear1,a.AScholte1,a.SLamb1,a.SShear1,a.SScholte1] = Computer_Isotropic_EnergyVelocity(a.ALamb1,a.AShear1,a.AScholte1,a.SLamb1,a.SShear1,a.SScholte1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,a.Thickness/1e3);
    clear Computer_Isotropic_SLamb_F Computer_Isotropic_ALamb_F Computer_Isotropic_SScholte_V Computer_Isotropic_AScholte_V Computer_Isotropic_SShear_V Computer_Isotropic_AShear_V        
    if  a.Multithreading == 1 && a.LambModes1 == 1 && a.SymmetricModes1 == 1 && a.AntisymmetricModes1 == 1
        delete(a.a1UI1)
        a.a1UI1 = axes('Parent',Tab1,'Units','pixels','Position',[1130 705 34 34]);
        imshow('Multithreading21.png','Parent',a.a1UI1)
    end
    a.ModeNames1 = {''};
    if  ~isempty(a.ALamb1{1})
        for i = 0:length(a.ALamb1)-1
            eval(sprintf('a.ModeNames1(i+1) = {''A%u''};',i));
        end
    end
    if  ~isempty(a.SLamb1{1})
        for i = 0:length(a.SLamb1)-1
            eval(sprintf('a.ModeNames1(length(a.ALamb1)+i+1) = {''S%u''};',i));
        end
    end
    if  ~isempty(a.AShear1{1})
        for i = 1:length(a.AShear1)
            eval(sprintf('a.ModeNames1(length(a.ALamb1)+length(a.SLamb1)+i) = {''ASH%u''};',i));
        end
    end
    if  ~isempty(a.SShear1{1})
        for i = 0:length(a.SShear1)-1
            eval(sprintf('a.ModeNames1(length(a.ALamb1)+length(a.SLamb1)+length(a.AShear1)+i+1) = {''SSH%u''};',i));
        end
    end
    if  ~isempty(a.AScholte1{1})
        for i = 0:length(a.AScholte1)-1
            eval(sprintf('a.ModeNames1(length(a.ALamb1)+length(a.SLamb1)+length(a.AShear1)+length(a.SShear1)+i+1) = {''AScholte%u''};',i));
        end
    end
    if  ~isempty(a.SScholte1{1})
        for i = 0:length(a.SScholte1)-1
            eval(sprintf('a.ModeNames1(length(a.ALamb1)+length(a.SLamb1)+length(a.AShear1)+length(a.SShear1)+length(a.AScholte1)+i+1) = {''SScholte%u''};',i));
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

    a.Title3 = 1;
    a.TitleUI3.Style = 'checkbox';
    a.TitleUI3.Value = 1;
    a.TitleUI3.String = '';
    a.TitleUI3.TooltipString = 'Check this in order to show the plot title.';
    a.TitleUI3.Position(3) = 20;
    a.DataUI3.String = a.Material1.Name;
    a.PlotUI3.Enable = 'off';
    a.DataType3 = 1;
    a.Mode3 = '';
    a.Frequency3 = a.FrequencyResolution1*1e2;
    a.FrequencyUI3.String = a.Frequency3;
    [a.ALambModes1,a.AShearModes1,a.SLambModes1,a.SShearModes1,a.Frequency3] = ModeFinder_Isotropic(a.ALamb1,a.AShear1,a.SLamb1,a.SShear1,a.Frequency3);    
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
        eval(sprintf('a.ALamb%uaUI3.String = [''A'',num2str(i)];',i));
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
elseif strcmp(source.Tag,'16') % Stop
    Stop = 1;
    clear Computer_Isotropic_SLamb_F Computer_Isotropic_ALamb_F Computer_Isotropic_SScholte_V Computer_Isotropic_AScholte_V Computer_Isotropic_SShear_V Computer_Isotropic_AShear_V
elseif strcmp(source.Tag,'17') % Quantity (Dispersion diagrams)
    switch source.Value
    case 1
        a.Quantity11 = 1;
        a.Option1UI1.Enable = 'on';
        a.Option1UI1.Style = 'checkbox';
        a.Option1UI1.String = ' ';
        a.Option1UI1.Value = a.BulkVelocities1;
        a.Option1UI1.TooltipString = 'Check this to show the bulk wave velocities.';
        a.Option1UI1.Position(3) = 50;
        a.Option1TextUI1.String = 'Bulk velocities';
        a.Option1TextUI1.Position(3) = 70;
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeTextUI1.Position(3) = 62;
        a.XAxisModeUI1.TooltipString = 'Select the frequency''s dimension on the X-axis.';        
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisTextUI1.Position(3) = 62;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisTextUI1.Position(3) = 65;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI1.Position(3) = 86;
            a.XAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;
        end
        a.YAxisTextUI1.String = 'Y-axis (m/ms)';
        a.YAxisTextUI1.Position(3) = 70;
        a.YAxisUI1.TooltipString = 'Enter which phase velocity range shall be plotted.';
        
        a.YAxisUI1.String = ['[0 ',num2str(a.PhaseVelocityLimit1/1e3),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    case 2
        a.Quantity11 = 2;
        a.Option1UI1.Enable = 'on';
        a.Option1UI1.Style = 'checkbox';
        a.Option1UI1.String = ' ';
        a.Option1UI1.Value = a.BulkVelocities1;
        a.Option1UI1.TooltipString = 'Check this to show the bulk wave velocities.';
        a.Option1UI1.Position(3) = 50;
        a.Option1TextUI1.String = 'Bulk velocities';
        a.Option1TextUI1.Position(3) = 70;
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeTextUI1.Position(3) = 62;
        a.XAxisModeUI1.TooltipString = 'Select the frequency''s dimension on the X-axis.';        
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisTextUI1.Position(3) = 62;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);            
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisTextUI1.Position(3) = 65;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;            
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI1.Position(3) = 86;
            a.XAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;            
        end
        a.YAxisTextUI1.String = 'Y-axis (m/ms)';
        a.YAxisTextUI1.Position(3) = 70;
        a.YAxisUI1.TooltipString = 'Enter which energy velocity range shall be plotted.';
        
        a.YAxisUI1.String = ['[0 ',num2str(ceil(a.Material1.PlateVelocity/1e3)),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    case 3
        a.Quantity11 = 3;
        a.Option1UI1.Enable = 'on';
        a.Option1UI1.Style = 'edit';
        a.Option1UI1.String = a.Distance1;
        a.Option1UI1.TooltipString = ['Enter the propagation distance from the source to the',newline,'sensor for which the propagation time shall be calculated.'];
        a.Option1UI1.Position(3) = 75;
        a.Option1TextUI1.String = 'Distance (mm)';
        a.Option1TextUI1.Position(3) = 71;
        a.XAxisModeTextUI1.String = 'Y-axis mode';
        a.XAxisModeTextUI1.Position(3) = 63;
        a.XAxisModeUI1.TooltipString = 'Select the frequency''s dimension on the Y-axis.';        
        a.XAxisTextUI1.String = ['X-axis (',char(181),'s)'];
        a.XAxisTextUI1.Position(3) = 56;
        a.XAxisUI1.TooltipString = 'Enter which propagation time range shall be plotted.';
        if  a.XAxisMode1 == 1
            a.YAxisTextUI1.String = 'Y-axis (kHz)';
            a.YAxisTextUI1.Position(3) = 63;
            a.YAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);             
        elseif a.XAxisMode1 == 2
            a.YAxisTextUI1.String = 'Y-axis (MHz)';
            a.YAxisTextUI1.Position(3) = 66;
            a.YAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3;            
        elseif a.XAxisMode1 == 3
            a.YAxisTextUI1.String = ['Y-axis (MHz',char(8901),'mm)'];
            a.YAxisTextUI1.Position(3) = 87;
            a.YAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.YAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness;            
        end
        
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
    case 4
        a.Quantity11 = 4;
        a.Option1UI1.Enable = 'on';
        a.Option1UI1.Style = 'popupmenu';
        a.Option1UI1.String = fieldnames(a.Materials.Fluid);
        a.Option1UI1.Value = find(strcmp(a.Couplant1.Name,a.Option1UI1.String));
        a.Option1UI1.TooltipString = 'Select for which surrounding medium the coincidence angle shall be displayed.';
        a.Option1UI1.Position(3) = 130;
        a.Option1TextUI1.String = 'Couplant';
        a.Option1TextUI1.Position(3) = 44;
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeTextUI1.Position(3) = 62;
        a.XAxisModeUI1.TooltipString = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisTextUI1.Position(3) = 62;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);            
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisTextUI1.Position(3) = 65;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;            
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI1.Position(3) = 86;
            a.XAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;            
        end
        a.YAxisTextUI1.String = ['Y-axis (',char(176),')'];
        a.YAxisTextUI1.Position(3) = 49;
        a.YAxisUI1.TooltipString = 'Enter which coincidence angle range shall be plotted.';
        
        x = round(8e4/a.Material1.PlateVelocity*a.Couplant1.Velocity/343);
        if  x > 90
            x = 90;
        end
        a.YAxisUI1.String = ['[0 ',num2str(x),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    case 5
        a.Quantity11 = 5;
        a.Option1UI1.Enable = 'off';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeTextUI1.Position(3) = 62;
        a.XAxisModeUI1.TooltipString = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisTextUI1.Position(3) = 62;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);            
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisTextUI1.Position(3) = 65;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;            
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI1.Position(3) = 86;
            a.XAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;             
        end
        x = a.YRange*1e-3*a.Material1.PlateVelocity*a.Thickness;
        if  x >= 1e3 && x < 1e4
            x = round(x,-2);
        elseif x >= 1e2 && x < 1e3
            x = round(x,-1);
        elseif x < 1e2
            x = round(x);
        end
        if  a.XAxisMode1 == 3
            a.YAxisTextUI1.String = 'Y-axis (mm/mm)';
            a.YAxisTextUI1.Position(3) = 80;
            a.YAxisUI1.TooltipString = 'Enter which wavelength per thickness range shall be plotted.'; 
            
            a.YAxisUI1.String = ['[0 ',num2str(x/a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        else
            a.YAxisTextUI1.String = 'Y-axis (mm)';
            a.YAxisTextUI1.Position(3) = 61;
            a.YAxisUI1.TooltipString = 'Enter which wavelength range shall be plotted.';
            
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);            
        end
    case 6
        a.Quantity11 = 6;
        a.Option1UI1.Enable = 'off';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeTextUI1.Position(3) = 62;
        a.XAxisModeUI1.TooltipString = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisTextUI1.Position(3) = 62;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);            
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisTextUI1.Position(3) = 65;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;            
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI1.Position(3) = 86;
            a.XAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;             
        end        
        x = 2*pi*a.FrequencyLimit1/a.Material1.RayleighVelocity;
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
            a.YAxisTextUI1.String = ['Y-axis (rad/mm',char(8901),'mm)'];
            a.YAxisTextUI1.Position(3) = 101;
            a.YAxisUI1.TooltipString = 'Enter which wavenumber-thickness range shall be plotted.'; 
            
            a.YAxisUI1.String = ['[0 ',num2str(x*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        else
            a.YAxisTextUI1.String = 'Y-axis (rad/mm)';
            a.YAxisTextUI1.Position(3) = 80;
            a.YAxisUI1.TooltipString = 'Enter which wavenumber range shall be plotted.';
            
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        end
    case 7
        a.Quantity11 = 7;
        a.Option1UI1.Enable = 'off';
        a.XAxisModeTextUI1.String = 'X-axis mode';
        a.XAxisModeTextUI1.Position(3) = 62;
        a.XAxisModeUI1.TooltipString = 'Select the frequency''s dimension on the X-axis.';
        if  a.XAxisMode1 == 1
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisTextUI1.Position(3) = 62;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1),']'];
            a.XAxis1 = eval(a.XAxisUI1.String);            
        elseif a.XAxisMode1 == 2
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisTextUI1.Position(3) = 65;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3;            
        elseif a.XAxisMode1 == 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI1.Position(3) = 86;
            a.XAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            a.XAxisUI1.String = ['[0 ',num2str(a.FrequencyLimit1/1e3*a.Thickness),']'];
            a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;             
        end
        if  a.FluidLoading1 == 1 && a.Viscoelastic1 == 1
            x = 1e4*(a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)+a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
        elseif a.FluidLoading1 == 1 && a.Viscoelastic1 == 0
            x = 1e4*a.Fluid1.Density*a.Fluid1.Velocity/(a.Fluid1.Density*a.Fluid1.Velocity+a.Material1.Density*a.Material1.PlateVelocity)/a.Thickness;
        elseif a.FluidLoading1 == 0 && a.Viscoelastic1 == 1
            x = 1e4*(a.Material1.LongitudinalAttenuation+a.Material1.TransverseAttenuation)/a.Thickness;
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
        if  a.XAxisMode1 == 3
            a.YAxisTextUI1.String = ['Y-axis (Np/m',char(8901),'mm)'];
            a.YAxisTextUI1.Position(3) = 90;
            a.YAxisUI1.TooltipString = 'Enter which attenuation-thickness range shall be plotted.'; 
            
            a.YAxisUI1.String = ['[0 ',num2str(x*a.Thickness),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        else
            a.YAxisTextUI1.String = 'Y-axis (Np/m)';
            a.YAxisTextUI1.Position(3) = 69;
            a.YAxisUI1.TooltipString = 'Enter which attenuation range shall be plotted.';
            
            a.YAxisUI1.String = ['[0 ',num2str(x),']'];
            a.YAxis1 = eval(a.YAxisUI1.String);
        end
    end
elseif strcmp(source.Tag,'18') % option (Dispersion diagrams)
    if  a.Quantity11 < 3
        a.BulkVelocities1 = source.Value;
    elseif a.Quantity11 == 3
        a.Distance1 = str2double(source.String);
        
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
        E = fieldnames(a.Materials.Fluid);
        switch source.Value
        case source.Value
            a.Couplant1 = getfield(a.Materials.Fluid,E{source.Value});
        end
        
        x = round(8e4/a.Material1.PlateVelocity*a.Couplant1.Velocity/343);
        if  x > 90
            x = 90;
        end
        a.YAxisUI1.String = ['[0 ',num2str(x),']'];
        a.YAxis1 = eval(a.YAxisUI1.String);
    end
elseif strcmp(source.Tag,'19') % X-Axis mode
    switch source.Value
    case 1
        if  a.Quantity11 ~= 3
            a.XAxisTextUI1.String = 'X-axis (kHz)';
            a.XAxisTextUI1.Position(3) = 62;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';

            val = eval(a.XAxisUI1.String);
            if  a.XAxisMode1 == 2
                a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
                a.XAxisUI1.String = ['[',num2str(val(1)*1e3),' ',num2str(val(2)*1e3),']'];
            elseif a.XAxisMode1 == 3
                a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;
                a.XAxisUI1.String = ['[',num2str(val(1)*1e3/a.Thickness),' ',num2str(val(2)*1e3/a.Thickness),']'];
            end            
        elseif a.Quantity11 == 3
            a.YAxisTextUI1.String = 'Y-axis (kHz)';
            a.YAxisTextUI1.Position(3) = 63;
            a.YAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            val = eval(a.YAxisUI1.String);
            if  a.XAxisMode1 == 2
                a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
                a.YAxisUI1.String = ['[',num2str(val(1)*1e3),' ',num2str(val(2)*1e3),']'];
            elseif a.XAxisMode1 == 3
                a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness;
                a.YAxisUI1.String = ['[',num2str(val(1)*1e3/a.Thickness),' ',num2str(val(2)*1e3/a.Thickness),']'];
            end
        end
        if  (a.Quantity11 == 5 || a.Quantity11 == 6 || a.Quantity11 == 7) && a.XAxisMode1 == 3
            if  a.Quantity11 == 5
                a.YAxisTextUI1.String = 'Y-axis (mm)';
                a.YAxisTextUI1.Position(3) = 61;                
                a.YAxisUI1.TooltipString = 'Enter which wavelength range shall be plotted.';
            elseif a.Quantity11 == 6
                a.YAxisTextUI1.String = 'Y-axis (rad/mm)';
                a.YAxisTextUI1.Position(3) = 80;                
                a.YAxisUI1.TooltipString = 'Enter which wavenumber range shall be plotted.';
            elseif a.Quantity11 == 7
                a.YAxisTextUI1.String = 'Y-axis (Np/m)';
                a.YAxisTextUI1.Position(3) = 69;
                a.YAxisUI1.TooltipString = 'Enter which attenuation range shall be plotted.';
            end

            val = eval(a.YAxisUI1.String);
            if  a.Quantity11 == 5
                a.YAxis1 = eval(a.YAxisUI1.String)*a.Thickness;
                a.YAxisUI1.String = ['[',num2str(val(1)*a.Thickness),' ',num2str(val(2)*a.Thickness),']'];
            elseif a.Quantity11 == 6 || a.Quantity11 == 7
                a.YAxis1 = eval(a.YAxisUI1.String)/a.Thickness;
                a.YAxisUI1.String = ['[',num2str(val(1)/a.Thickness),' ',num2str(val(2)/a.Thickness),']'];
            end
        end
        
        a.XAxisMode1 = 1;
    case 2
        if  a.Quantity11 ~= 3
            a.XAxisTextUI1.String = 'X-axis (MHz)';
            a.XAxisTextUI1.Position(3) = 65;
            a.XAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            val = eval(a.XAxisUI1.String);
            if  a.XAxisMode1 == 1
                a.XAxis1 = eval(a.XAxisUI1.String);
                a.XAxisUI1.String = ['[',num2str(val(1)/1e3),' ',num2str(val(2)/1e3),']'];
            elseif a.XAxisMode1 == 3
                a.XAxis1 = eval(a.XAxisUI1.String)*1e3/a.Thickness;
                a.XAxisUI1.String = ['[',num2str(val(1)/a.Thickness),' ',num2str(val(2)/a.Thickness),']'];
            end
        elseif a.Quantity11 == 3
            a.YAxisTextUI1.String = 'Y-axis (MHz)';
            a.YAxisTextUI1.Position(3) = 66;
            a.YAxisUI1.TooltipString = 'Enter which frequency range shall be plotted.';
            
            val = eval(a.YAxisUI1.String);
            if  a.XAxisMode1 == 1
                a.YAxis1 = eval(a.YAxisUI1.String);
                a.YAxisUI1.String = ['[',num2str(val(1)/1e3),' ',num2str(val(2)/1e3),']'];
            elseif a.XAxisMode1 == 3
                a.YAxis1 = eval(a.YAxisUI1.String)*1e3/a.Thickness;
                a.YAxisUI1.String = ['[',num2str(val(1)/a.Thickness),' ',num2str(val(2)/a.Thickness),']'];
            end
        end
        if  (a.Quantity11 == 5 || a.Quantity11 == 6 || a.Quantity11 == 7) && a.XAxisMode1 == 3
            if  a.Quantity11 == 5
                a.YAxisTextUI1.String = 'Y-axis (mm)';
                a.YAxisTextUI1.Position(3) = 61;                
                a.YAxisUI1.TooltipString = 'Enter which wavelength range shall be plotted.';
            elseif a.Quantity11 == 6
                a.YAxisTextUI1.String = 'Y-axis (rad/mm)';
                a.YAxisTextUI1.Position(3) = 80;                
                a.YAxisUI1.TooltipString = 'Enter which wavenumber range shall be plotted.';
            elseif a.Quantity11 == 7
                a.YAxisTextUI1.String = 'Y-axis (Np/m)';
                a.YAxisTextUI1.Position(3) = 69;
                a.YAxisUI1.TooltipString = 'Enter which attenuation range shall be plotted.';
            end

            val = eval(a.YAxisUI1.String);
            if  a.Quantity11 == 5
                a.YAxis1 = eval(a.YAxisUI1.String)*a.Thickness;
                a.YAxisUI1.String = ['[',num2str(val(1)*a.Thickness),' ',num2str(val(2)*a.Thickness),']'];
            elseif a.Quantity11 == 6 || a.Quantity11 == 7
                a.YAxis1 = eval(a.YAxisUI1.String)/a.Thickness;
                a.YAxisUI1.String = ['[',num2str(val(1)/a.Thickness),' ',num2str(val(2)/a.Thickness),']'];                
            end           
        end
        
        a.XAxisMode1 = 2;
    case 3
        if  a.Quantity11 ~= 3
            a.XAxisTextUI1.String = ['X-axis (MHz',char(8901),'mm)'];
            a.XAxisTextUI1.Position(3) = 86;
            a.XAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            val = eval(a.XAxisUI1.String);
            if  a.XAxisMode1 == 1
                a.XAxis1 = eval(a.XAxisUI1.String);
                a.XAxisUI1.String = ['[',num2str(val(1)/1e3*a.Thickness),' ',num2str(val(2)/1e3*a.Thickness),']'];
            elseif a.XAxisMode1 == 2
                a.XAxis1 = eval(a.XAxisUI1.String)*1e3;
                a.XAxisUI1.String = ['[',num2str(val(1)*a.Thickness),' ',num2str(val(2)*a.Thickness),']'];
            end            
        elseif a.Quantity11 == 3
            a.YAxisTextUI1.String = ['Y-axis (MHz',char(8901),'mm)'];
            a.YAxisTextUI1.Position(3) = 87;
            a.YAxisUI1.TooltipString = 'Enter which frequency-thickness range shall be plotted.';
            
            val = eval(a.YAxisUI1.String);
            if  a.XAxisMode1 == 1
                a.YAxis1 = eval(a.YAxisUI1.String);
                a.YAxisUI1.String = ['[',num2str(val(1)/1e3*a.Thickness),' ',num2str(val(2)/1e3*a.Thickness),']'];
            elseif a.XAxisMode1 == 2
                a.YAxis1 = eval(a.YAxisUI1.String)*1e3;
                a.YAxisUI1.String = ['[',num2str(val(1)*a.Thickness),' ',num2str(val(2)*a.Thickness),']'];
            end
        end
        if  (a.Quantity11 == 5 || a.Quantity11 == 6 || a.Quantity11 == 7) && a.XAxisMode1 ~= 3 
            if  a.Quantity11 == 5
                a.YAxisTextUI1.String = 'Y-axis (mm/mm)';
                a.YAxisTextUI1.Position(3) = 80;
                a.YAxisUI1.TooltipString = 'Enter which wavelength per thickness range shall be plotted.';
            elseif a.Quantity11 == 6
                a.YAxisTextUI1.String = ['Y-axis (rad/mm',char(8901),'mm)'];
                a.YAxisTextUI1.Position(3) = 101;
                a.YAxisUI1.TooltipString = 'Enter which wavenumber-thickness range shall be plotted.';
            elseif a.Quantity11 == 7
                a.YAxisTextUI1.String = ['Y-axis (Np/m',char(8901),'mm)'];
                a.YAxisTextUI1.Position(3) = 90;
                a.YAxisUI1.TooltipString = 'Enter which attenuation-thickness range shall be plotted.'; 
            end

            val = eval(a.YAxisUI1.String);
            if  a.Quantity11 == 5
                a.YAxis1 = eval(a.YAxisUI1.String)/a.Thickness;
                a.YAxisUI1.String = ['[',num2str(val(1)/a.Thickness),' ',num2str(val(2)/a.Thickness),']'];
            elseif a.Quantity11 == 6 || a.Quantity11 == 7
                a.YAxis1 = eval(a.YAxisUI1.String)*a.Thickness;
                a.YAxisUI1.String = ['[',num2str(val(1)*a.Thickness),' ',num2str(val(2)*a.Thickness),']'];                
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
            a.XAxis1 = eval(source.String)*1e3/a.Thickness;
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
            a.YAxis1 = eval(source.String)*1e3/a.Thickness;
        end          
    end    
elseif strcmp(source.Tag,'22') % Plot (Dispersion diagrams) 
    if  a.Quantity11 == 1
        PhaseVelocity_Isotropic(a.CropPlots1,a.FluidLoading1,a.Fluid1,a.BulkVelocities1,a.HigherOrderModes1,a.PNGresolution1,a.S0Label,a.SH0Label1,a.A0Label,a.LambModes1,a.ShearHorizontalModes1,a.ScholteModes1,a.SColor1,a.AColor1,a.ALamb1,a.AShear1,a.AScholte1,a.SymmetricModes1,a.AntisymmetricModes1,a.BoxLineWidth1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.ModeLabelFontSize1,a.Title1,a.LineWidth1,a.Material1,a.ModeLabels1,a.PDF1,a.FileName1,a.Thickness,a.PNG1,a.SLamb1,a.SShear1,a.SScholte1,a.XAxis1,a.XAxisMode1,a.YAxis1)
    elseif a.Quantity11 == 2
        EnergyVelocity_Isotropic(a.CropPlots1,a.FluidLoading1,a.Fluid1,a.BulkVelocities1,a.HigherOrderModes1,a.PNGresolution1,a.S0Label,a.SH0Label1,a.A0Label,a.LambModes1,a.ShearHorizontalModes1,a.ScholteModes1,a.SColor1,a.AColor1,a.ALamb1,a.AShear1,a.AScholte1,a.SymmetricModes1,a.AntisymmetricModes1,a.BoxLineWidth1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.ModeLabelFontSize1,a.Title1,a.LineWidth1,a.Material1,a.ModeLabels1,a.PDF1,a.FileName1,a.Thickness,a.PNG1,a.SLamb1,a.SShear1,a.SScholte1,a.XAxis1,a.XAxisMode1,a.YAxis1)
    elseif a.Quantity11 == 3
        PropagationTime_Isotropic(a.CropPlots1,a.FluidLoading1,a.Fluid1,a.Distance1,a.HigherOrderModes1,a.PNGresolution1,a.S0Label,a.SH0Label1,a.A0Label,a.LambModes1,a.ShearHorizontalModes1,a.ScholteModes1,a.SColor1,a.AColor1,a.ALamb1,a.AShear1,a.AScholte1,a.SymmetricModes1,a.AntisymmetricModes1,a.BoxLineWidth1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.ModeLabelFontSize1,a.Title1,a.LineWidth1,a.Material1,a.ModeLabels1,a.PDF1,a.FileName1,a.Thickness,a.PNG1,a.SLamb1,a.SShear1,a.SScholte1,a.XAxis1,a.XAxisMode1,a.YAxis1)
    elseif a.Quantity11 == 4
        CoincidenceAngle_Isotropic(a.CropPlots1,a.FluidLoading1,a.Fluid1,a.HigherOrderModes1,a.PNGresolution1,a.Couplant1,a.S0Label,a.SH0Label1,a.A0Label,a.LambModes1,a.ShearHorizontalModes1,a.ScholteModes1,a.SColor1,a.AColor1,a.ALamb1,a.AShear1,a.AScholte1,a.SymmetricModes1,a.AntisymmetricModes1,a.BoxLineWidth1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.ModeLabelFontSize1,a.Title1,a.LineWidth1,a.Material1,a.ModeLabels1,a.PDF1,a.FileName1,a.Thickness,a.PNG1,a.SLamb1,a.SShear1,a.SScholte1,a.XAxis1,a.XAxisMode1,a.YAxis1)
    elseif a.Quantity11 == 5
        Wavelength_Isotropic(a.CropPlots1,a.FluidLoading1,a.Fluid1,a.HigherOrderModes1,a.PNGresolution1,a.S0Label,a.SH0Label1,a.A0Label,a.LambModes1,a.ShearHorizontalModes1,a.ScholteModes1,a.SColor1,a.AColor1,a.ALamb1,a.AShear1,a.AScholte1,a.SymmetricModes1,a.AntisymmetricModes1,a.BoxLineWidth1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.ModeLabelFontSize1,a.Title1,a.LineWidth1,a.Material1,a.ModeLabels1,a.PDF1,a.FileName1,a.Thickness,a.PNG1,a.SLamb1,a.SShear1,a.SScholte1,a.XAxis1,a.XAxisMode1,a.YAxis1)
    elseif a.Quantity11 == 6
        Wavenumber_Isotropic(a.CropPlots1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,a.HigherOrderModes1,a.PNGresolution1,a.S0Label,a.SH0Label1,a.A0Label,a.LambModes1,a.ShearHorizontalModes1,a.ScholteModes1,a.SColor1,a.AColor1,a.ALamb1,a.AShear1,a.AScholte1,a.SymmetricModes1,a.AntisymmetricModes1,a.BoxLineWidth1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.ModeLabelFontSize1,a.Title1,a.LineWidth1,a.Material1,a.ModeLabels1,a.PDF1,a.FileName1,a.Thickness,a.PNG1,a.SLamb1,a.SShear1,a.SScholte1,a.XAxis1,a.XAxisMode1,a.YAxis1)
    elseif a.Quantity11 == 7
        Attenuation_Isotropic(a.CropPlots1,a.FluidLoading1,a.Fluid1,a.HigherOrderModes1,a.PNGresolution1,a.S0Label,a.A0Label,a.LambModes1,a.ShearHorizontalModes1,a.ScholteModes1,a.SColor1,a.AColor1,a.ALamb1,a.AShear1,a.AScholte1,a.SymmetricModes1,a.AntisymmetricModes1,a.BoxLineWidth1,a.Directory1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.ModeLabelFontSize1,a.Title1,a.LineWidth1,a.Material1,a.ModeLabels1,a.PDF1,a.FileName1,a.Thickness,a.PNG1,a.SLamb1,a.SShear1,a.SScholte1,a.XAxis1,a.XAxisMode1,a.YAxis1)
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
        a.x13UI1.TooltipString = 'Displacement u1';
        a.x23UI1.TooltipString = 'Displacement u2';
        a.x33UI1.TooltipString = 'Displacement u3';
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
        a.x11UI1.TooltipString = ['Stress ',char(963),'11'];
        a.x12UI1.TooltipString = ['Stress ',char(963),'12'];
        a.x13UI1.TooltipString = ['Stress ',char(963),'13'];
        a.x22UI1.TooltipString = ['Stress ',char(963),'22'];
        a.x23UI1.TooltipString = ['Stress ',char(963),'23'];
        a.x33UI1.TooltipString = ['Stress ',char(963),'33'];
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
        a.x11UI1.TooltipString = ['Strain ',char(949),'11'];
        a.x12UI1.TooltipString = ['Strain ',char(949),'12'];
        a.x13UI1.TooltipString = ['Strain ',char(949),'13'];
        a.x22UI1.TooltipString = ['Strain ',char(949),'22'];
        a.x23UI1.TooltipString = ['Strain ',char(949),'23'];
        a.x33UI1.TooltipString = ['Strain ',char(949),'33'];
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
        a.x13UI1.TooltipString = 'Strain energy density';
        a.x23UI1.TooltipString = 'Kinetic energy density';
        a.x33UI1.TooltipString = 'Total energy density';        
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
        a.x13UI1.TooltipString = 'Power flow density p1';
        a.x23UI1.TooltipString = 'Power flow density p2';
        a.x33UI1.TooltipString = 'Power flow density p3';                
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
    if  a.Quantity21 == 1
        u1Color = [1 0 0];
        u2Color = [.13 .55 .13];
        u3Color = [0 0 1];
        ModeShapeLines_Isotropic(1,1,0,0,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,u1Color,u2Color,u3Color,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
    elseif a.Quantity21 == 2
        sigma33Color = [0 0 1];
        sigma11Color = [1 0 0];
        sigma22Color = [.13 .55 .13];
        sigma13Color = [0 0 0];
        sigma23Color = [1 0 1];
        sigma12Color = [0 1 1]; 
        ModeShapeLines_Isotropic(2,1,0,0,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,sigma11Color,sigma22Color,sigma33Color,sigma23Color,sigma13Color,sigma12Color,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
    elseif a.Quantity21 == 3
        epsilon33Color = [0 0 1];
        epsilon11Color = [1 0 0];
        epsilon22Color = [.13 .55 .13];
        epsilon13Color = [0 0 0];
        epsilon23Color = [1 0 1];
        epsilon12Color = [0 1 1];
        ModeShapeLines_Isotropic(3,1,0,0,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,epsilon11Color,epsilon22Color,epsilon33Color,epsilon23Color,epsilon13Color,epsilon12Color,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
    elseif a.Quantity21 == 4
        StrainEnergyDensityColor = [0 0 1];
        KineticEnergyDensityColor = [1 0 0];
        TotalEnergyDensityColor = [0 0 0];
        ModeShapeLines_Isotropic(4,1,0,0,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,StrainEnergyDensityColor,KineticEnergyDensityColor,TotalEnergyDensityColor,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
    elseif a.Quantity21 == 5
        P1Color = [1 0 0];
        P2Color = [.13 .55 .13];
        P3Color = [0 0 1];
        ModeShapeLines_Isotropic(5,1,0,0,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,P1Color,P2Color,P3Color,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
    end
elseif strcmp(source.Tag,'34') % Mode (mode shape)
    a.Mode21 = a.ModeNames1{source.Value};
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'35') % Frequency (mode shape)
    a.Frequency21 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'36') % Wavelengths
    a.Length1 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'37') % Samples x1
    a.Samples21 = str2double(source.String);
    
    if  a.FluidLoading1 == 1 && a.Halfspaces21 == 1
        a.Samples31 = round(.5*a.Samples21/(2*a.HalfspacesNumber21+1));
    else
        a.Samples31 = round(.5*a.Samples21);
    end
    if  a.Samples31 == 0
        a.Samples31 = 1;
    end
    a.Samples3UI1.String = a.Samples31;
    a.GridLine1 = 2*round(a.Samples21/80);
    a.GridLineUI1.String = a.GridLine1; 
    
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'38') % Samples x3 (mode shape)
    a.Samples31 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'39') % Scale
    a.Scale1 = str2double(source.String);
elseif strcmp(source.Tag,'40') % Grid line
    a.GridLine1 = str2double(source.String);
elseif strcmp(source.Tag,'41') % Undistorted
    a.Undistorted1 = source.Value;
elseif strcmp(source.Tag,'82') % Half-spaces (mode shape)
    a.HalfspacesNumber21 = str2double(source.String);
    
    if  a.FluidLoading1 == 1 && a.Halfspaces21 == 1
        a.Samples31 = round(.5*a.Samples21/(2*a.HalfspacesNumber21+1));
    else
        a.Samples31 = round(.5*a.Samples21);
    end
    if  a.Samples31 == 0
        a.Samples31 = 1;
    end
    a.Samples3UI1.String = a.Samples31;
    
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'83') % Half-spaces (mode shape)
    a.Halfspaces21 = source.Value;
    
    if  a.FluidLoading1 == 1 && a.Halfspaces21 == 1
        a.Samples31 = round(.5*a.Samples21/(2*a.HalfspacesNumber21+1));
    else
        a.Samples31 = round(.5*a.Samples21);
    end
    if  a.Samples31 == 0
        a.Samples31 = 1;
    end
    a.Samples3UI1.String = a.Samples31;
    
    a.ModeShapeSettingChanged1 = 1;
elseif strcmp(source.Tag,'43') % Cycles
    a.Cycles1 = str2double(source.String);
    a.ModeShapeSettingChanged1 = 1;
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
    if  a.Animate1 == 0
        ModeShapeImage_Isotropic(a.CropPlots1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,a.ALamb1,a.AShear1,a.AScholte1,a.Directory1,a.ExportPlots1,a.TitleFontSize1,a.AxesLabelFontSize1,a.Frequency21,a.GridLine1,a.Title1,a.Length1,a.LineWidth1,a.Mode21,a.FileName1,a.PDF1,a.PNG1,a.Thickness/1e3,a.Samples21,a.Samples31,a.Scale1,a.SLamb1,a.SShear1,a.SScholte1,a.Undistorted1,a.Halfspaces21,a.HalfspacesNumber21)
    else
        if  a.ModeShapeSettingChanged1 == 1
            [a.Time1,a.u1,a.x11,a.x31,a.p1] = ModeShapeAnimationComputer_Isotropic(a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,a.ALamb1,a.AShear1,a.AScholte1,a.CycleDuration1,a.Cycles1,a.FrameRate1,a.Frequency21,a.Length1,a.Mode21,a.Thickness/1e3,a.Samples21,a.Samples31,a.SLamb1,a.SShear1,a.SScholte1,a.Halfspaces21,a.HalfspacesNumber21);
            a.ModeShapeSettingChanged1 = 0;
        end
        ModeShapeAnimation_Isotropic(a.Directory1,a.FluidLoading1,a.Fluid1,a.TitleFontSize1,a.AxesLabelFontSize1,a.FrameRate1,a.Frequency21,a.GridLine1,a.Title1,a.LineWidth1,a.Material1,a.Mode21,a.ExportPlots1,a.FileName1,a.MovieQuality1,a.Thickness,a.Samples21,a.Samples31,a.Scale1,a.Time1,a.u1,a.Undistorted1,a.x11,a.x31,a.p1,a.Halfspaces21,a.HalfspacesNumber21)
    end
elseif strcmp(source.Tag,'48') % Export plots
    a.ExportPlots1 = source.Value;

    if  source.Value == 1
        a.Plot1UI1.String = 'Export';
        a.Plot2UI1.String = 'Export';
        a.Plot3UI1.String = 'Export';
    else
        a.Plot1UI1.String = 'Plot';
        a.Plot2UI1.String = 'Plot';
        a.Plot3UI1.String = 'Plot';
    end
elseif strcmp(source.Tag,'77') % Crop plots
    a.CropPlots1 = source.Value;
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
elseif strcmp(source.Tag,'78') % Matlab
    if  a.DispersionCurves1 == 1
        Export_Isotropic(a,0,0,1)
    end
    if  a.ThroughThickness1 == 1
        if  a.Quantity21 == 1
            ModeShapeLines_Isotropic(1,0,1,0,0,1,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 2
            ModeShapeLines_Isotropic(2,0,1,0,0,1,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 3
            ModeShapeLines_Isotropic(3,0,1,0,0,1,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 4
            ModeShapeLines_Isotropic(4,0,1,0,0,1,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
        elseif a.Quantity21 == 5
            ModeShapeLines_Isotropic(5,0,1,0,0,1,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
        end
    end
elseif strcmp(source.Tag,'56') % Excel
    if  a.DispersionCurves1 == 1
        Export_Isotropic(a,1,0,0)
    end
    if  a.ThroughThickness1 == 1
        if  a.Quantity21 == 1
            ModeShapeLines_Isotropic(1,0,1,1,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 2
            ModeShapeLines_Isotropic(2,0,1,1,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 3
            ModeShapeLines_Isotropic(3,0,1,1,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 4
            ModeShapeLines_Isotropic(4,0,1,1,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
        elseif a.Quantity21 == 5
            ModeShapeLines_Isotropic(5,0,1,1,0,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
        end
    end
elseif strcmp(source.Tag,'57') % TXT
    if  a.DispersionCurves1 == 1
        Export_Isotropic(a,0,1,0)
    end
    if  a.ThroughThickness1 == 1
        if  a.Quantity21 == 1
            ModeShapeLines_Isotropic(1,0,1,0,1,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 2
            ModeShapeLines_Isotropic(2,0,1,0,1,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 3
            ModeShapeLines_Isotropic(3,0,1,0,1,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,a.Phase1)
        elseif a.Quantity21 == 4
            ModeShapeLines_Isotropic(4,0,1,0,1,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
        elseif a.Quantity21 == 5
            ModeShapeLines_Isotropic(5,0,1,0,1,0,a.CropPlots1,a.Plot1,a.PNGresolution1,a.Material1,a.Viscoelastic1,a.FluidLoading1,a.Fluid1,0,0,0,0,0,0,a.ALamb1,a.AShear1,a.AScholte1,a.BoxLineWidth1,a.Directory1,a.FileName1,a.ExportPlots1,a.AxesTickFontSize1,a.AxesLabelFontSize1,a.TitleFontSize1,a.LegendFontSize1,a.Frequency11,a.Title1,a.LegendLocation1,a.LineWidth1,a.Mode11,a.PDF1,a.PNG1,a.Thickness/1e3,a.SLamb1,a.SShear1,a.SScholte1,a.Samples11,a.Halfspaces11,a.HalfspacesNumber11,0)
        end
    end
elseif strcmp(source.Tag,'58') % File name
    a.FileName1 = source.String;
elseif strcmp(source.Tag,'59') % Directory
    a.Directory1 = source.String;
elseif strcmp(source.Tag,'60') % Title
    a.Title1 = source.Value;
elseif strcmp(source.Tag,'61') % Mode labels 
    a.ModeLabels1 = source.Value;
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
elseif strcmp(source.Tag,'68') % S0
    a.S0Label = str2double(source.String);
elseif strcmp(source.Tag,'69') % S'0
    a.SH0Label1 = str2double(source.String);
elseif strcmp(source.Tag,'70') % A0
    a.A0Label = str2double(source.String);
elseif strcmp(source.Tag,'71') % Title (Font size)
    a.TitleFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'72') % Axes labels
    a.AxesLabelFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'73') % Axes ticks
    a.AxesTickFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'74') % Mode labels
    a.ModeLabelFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'75') % Legend
    a.LegendFontSize1 = str2double(source.String);
elseif strcmp(source.Tag,'76') % Default
    a.Title1 = 1;
    a.TitleUI1.Value = a.Title1;
    a.ModeLabels1 = 1;
    a.ModeLabelsUI1.Value = a.ModeLabels1;
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
    a.S0Label = .05;
    a.S0LabelUI1.String = a.S0Label;
    a.SH0Label1 = .05;
    a.SH0LabelUI1.String = a.SH0Label1;
    a.A0Label = .05;
    a.A0LabelUI1.String = a.A0Label;
    a.TitleFontSize1 = 30;
    a.TitleFontSizeUI1.String = a.TitleFontSize1;
    a.AxesLabelFontSize1 = 30;
    a.AxesLabelFontSizeUI1.String = a.AxesLabelFontSize1;
    a.AxesTickFontSize1 = 24;
    a.AxesTickFontSizeUI1.String = a.AxesTickFontSize1;
    a.ModeLabelFontSize1 = 24;
    a.ModeLabelFontSizeUI1.String = a.ModeLabelFontSize1;
    a.LegendFontSize1 = 24;
    a.LegendFontSizeUI1.String = a.LegendFontSize1;
end