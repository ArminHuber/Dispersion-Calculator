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
function a = CallbackModule_MaterialEditor(source,~,a)
%#ok<*GFLD>
%#ok<*AGROW>
%#ok<*FXUP>
source.String = replace(source.String,',','.');
if  strcmp(source.Tag,'1') % Material (isotropic)
    a.Material1ME = getfield(a.Materials.Isotropic,cell2mat(source.String(source.Value)));
    a.Density1UI5.String = a.Material1ME.Density;
    a.YoungsModulusUI5.String = real(a.Material1ME.YoungsModulus_complex)/1e9;
    a.YoungsModulusImagUI5.String = imag(a.Material1ME.YoungsModulus_complex)/1e9;
    a.PoissonsNumberUI5.String = real(a.Material1ME.PoissonsNumber_complex);
    a.PoissonsNumberImagUI5.String = imag(a.Material1ME.PoissonsNumber_complex);
    a.C11IsoUI5.String = real(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/1e9;
    a.C11ImagIsoUI5.String = imag(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/1e9;
    a.C66IsoUI5.String = real(a.Material1ME.Mu_complex)/1e9;
    a.C66ImagIsoUI5.String = imag(a.Material1ME.Mu_complex)/1e9;
    a.LongitudinalVelocityUI5.String = a.Material1ME.LongitudinalVelocity;
    a.TransverseVelocityUI5.String = a.Material1ME.TransverseVelocity;
    if  a.AttenuationUnit1 == 1 % Np/Lambda
        a.LongitudinalAttenuationUI5.String = a.Material1ME.LongitudinalAttenuation;
        a.TransverseAttenuationUI5.String = a.Material1ME.TransverseAttenuation;
    elseif a.AttenuationUnit1 == 2 % dB/Lambda
        a.LongitudinalAttenuationUI5.String = 20*log10(exp(a.Material1ME.LongitudinalAttenuation));
        a.TransverseAttenuationUI5.String = 20*log10(exp(a.Material1ME.TransverseAttenuation));
    elseif a.AttenuationUnit1 == 3 % Np/m
        a.LongitudinalAttenuationUI5.String = a.Material1ME.LongitudinalAttenuation*a.AtFrequency1*1e3/a.Material1ME.LongitudinalVelocity;
        a.TransverseAttenuationUI5.String = a.Material1ME.TransverseAttenuation*a.AtFrequency1*1e3/a.Material1ME.TransverseVelocity;
    elseif a.AttenuationUnit1 == 4 % dB/m
        a.LongitudinalAttenuationUI5.String = 20*log10(exp(a.Material1ME.LongitudinalAttenuation))*a.AtFrequency1*1e3/a.Material1ME.LongitudinalVelocity;
        a.TransverseAttenuationUI5.String = 20*log10(exp(a.Material1ME.TransverseAttenuation))*a.AtFrequency1*1e3/a.Material1ME.TransverseVelocity;
    end
    a.MaterialName4A = a.Material1ME.Name;
    a.Name1UI5.String = a.MaterialName4A;
elseif strcmp(source.Tag,'2') % Density (isotropic)
    a.Material1ME.Density = str2double(source.String);
    
    a.Material1ME.LongitudinalVelocity = sqrt((a.Material1ME.Lambda+2*a.Material1ME.Mu)/a.Material1ME.Density);
    a.Material1ME.TransverseVelocity = sqrt(a.Material1ME.Mu/a.Material1ME.Density);
    a.Material1ME.PlateVelocity = sqrt(a.Material1ME.YoungsModulus/(a.Material1ME.Density*(1-a.Material1ME.PoissonsNumber^2)));
    a.Material1ME.CylinderVelocity = sqrt(a.Material1ME.YoungsModulus/a.Material1ME.Density);
    a.Material1ME.RayleighVelocity = a.Material1ME.TransverseVelocity*((.87+1.12*a.Material1ME.PoissonsNumber)./(1+a.Material1ME.PoissonsNumber));
    a.Material1ME.LongitudinalVelocity_complex = a.Material1ME.LongitudinalVelocity/(1+1i*a.Material1ME.LongitudinalAttenuation/(2*pi));
    a.Material1ME.TransverseVelocity_complex = a.Material1ME.TransverseVelocity/(1+1i*a.Material1ME.TransverseAttenuation/(2*pi));

    a.LongitudinalVelocityUI5.String = a.Material1ME.LongitudinalVelocity; 
    a.TransverseVelocityUI5.String = a.Material1ME.TransverseVelocity;
elseif strcmp(source.Tag,'3') || strcmp(source.Tag,'4') || strcmp(source.Tag,'5') || strcmp(source.Tag,'6') % Young's modulus and Poisson's ratio
    if  strcmp(source.Tag,'3')
        a.Material1ME.YoungsModulus_complex = str2double(source.String)*1e9+1i*imag(a.Material1ME.YoungsModulus_complex);
    elseif strcmp(source.Tag,'4')
        a.Material1ME.YoungsModulus_complex = real(a.Material1ME.YoungsModulus_complex)+1i*str2double(source.String)*1e9;
    elseif strcmp(source.Tag,'5')
        a.Material1ME.PoissonsNumber_complex = str2double(source.String)+1i*imag(a.Material1ME.PoissonsNumber_complex);
    elseif strcmp(source.Tag,'6')
        a.Material1ME.PoissonsNumber_complex = real(a.Material1ME.PoissonsNumber_complex)+1i*str2double(source.String);
    end

    a.Material1ME.Lambda_complex = a.Material1ME.YoungsModulus_complex*a.Material1ME.PoissonsNumber_complex/((1+a.Material1ME.PoissonsNumber_complex)*(1-2*a.Material1ME.PoissonsNumber_complex));
    a.Material1ME.Mu_complex = a.Material1ME.YoungsModulus_complex/(2*(1+a.Material1ME.PoissonsNumber_complex));
    a.Material1ME.LongitudinalVelocity_complex = conj(sqrt((a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/a.Material1ME.Density));
    a.Material1ME.TransverseVelocity_complex = conj(sqrt((a.Material1ME.Mu_complex)/a.Material1ME.Density));
    a.Material1ME.LongitudinalAttenuation = -2*pi*imag(a.Material1ME.LongitudinalVelocity_complex)/real(a.Material1ME.LongitudinalVelocity_complex);
    a.Material1ME.TransverseAttenuation = -2*pi*imag(a.Material1ME.TransverseVelocity_complex)/real(a.Material1ME.TransverseVelocity_complex);
    a.Material1ME.LongitudinalVelocity = real(a.Material1ME.LongitudinalVelocity_complex)-a.Material1ME.LongitudinalAttenuation*imag(a.Material1ME.LongitudinalVelocity_complex)/(2*pi);
    a.Material1ME.TransverseVelocity = real(a.Material1ME.TransverseVelocity_complex)-a.Material1ME.TransverseAttenuation*imag(a.Material1ME.TransverseVelocity_complex)/(2*pi);
    a.Material1ME.Mu = a.Material1ME.Density*a.Material1ME.TransverseVelocity^2;
    a.Material1ME.Lambda = a.Material1ME.Density*a.Material1ME.LongitudinalVelocity^2-2*a.Material1ME.Mu;
    a.Material1ME.YoungsModulus = a.Material1ME.Mu*((3*a.Material1ME.Lambda+2*a.Material1ME.Mu)/(a.Material1ME.Lambda+a.Material1ME.Mu));
    a.Material1ME.PoissonsNumber = a.Material1ME.Lambda/(2*(a.Material1ME.Lambda+a.Material1ME.Mu));
    a.Material1ME.PlateVelocity = sqrt(a.Material1ME.YoungsModulus/(a.Material1ME.Density*(1-a.Material1ME.PoissonsNumber^2)));
    a.Material1ME.CylinderVelocity = sqrt(a.Material1ME.YoungsModulus/a.Material1ME.Density);
    a.Material1ME.RayleighVelocity = a.Material1ME.TransverseVelocity*((.87+1.12*a.Material1ME.PoissonsNumber)./(1+a.Material1ME.PoissonsNumber));    
    
    a.C11IsoUI5.String = real(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/1e9;
    a.C11ImagIsoUI5.String = imag(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/1e9;
    a.C66IsoUI5.String = real(a.Material1ME.Mu_complex)/1e9;
    a.C66ImagIsoUI5.String = imag(a.Material1ME.Mu_complex)/1e9;
    if  a.AttenuationUnit1 == 1 % Np/Lambda
        a.LongitudinalAttenuationUI5.String = a.Material1ME.LongitudinalAttenuation;
        a.TransverseAttenuationUI5.String = a.Material1ME.TransverseAttenuation;
    elseif a.AttenuationUnit1 == 2 % dB/Lambda
        a.LongitudinalAttenuationUI5.String = 20*log10(exp(a.Material1ME.LongitudinalAttenuation));
        a.TransverseAttenuationUI5.String = 20*log10(exp(a.Material1ME.TransverseAttenuation));
    elseif a.AttenuationUnit1 == 3 % Np/m
        a.LongitudinalAttenuationUI5.String = a.Material1ME.LongitudinalAttenuation*a.AtFrequency1*1e3/a.Material1ME.LongitudinalVelocity;
        a.TransverseAttenuationUI5.String = a.Material1ME.TransverseAttenuation*a.AtFrequency1*1e3/a.Material1ME.TransverseVelocity;
    elseif a.AttenuationUnit1 == 4 % dB/m
        a.LongitudinalAttenuationUI5.String = 20*log10(exp(a.Material1ME.LongitudinalAttenuation))*a.AtFrequency1*1e3/a.Material1ME.LongitudinalVelocity;
        a.TransverseAttenuationUI5.String = 20*log10(exp(a.Material1ME.TransverseAttenuation))*a.AtFrequency1*1e3/a.Material1ME.TransverseVelocity;
    end
    a.LongitudinalVelocityUI5.String = a.Material1ME.LongitudinalVelocity;
    a.TransverseVelocityUI5.String = a.Material1ME.TransverseVelocity;
elseif strcmp(source.Tag,'64') || strcmp(source.Tag,'65') || strcmp(source.Tag,'66') || strcmp(source.Tag,'67') % C11 and C66 (isotropic)
    if  strcmp(source.Tag,'64')
        C11 = str2double(source.String)*1e9+1i*imag(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex);
    elseif strcmp(source.Tag,'65')
        C11 = real(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)+1i*str2double(source.String)*1e9;
    elseif strcmp(source.Tag,'66')
        C11 = a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex;
        a.Material1ME.Mu_complex = str2double(source.String)*1e9+1i*imag(a.Material1ME.Mu_complex); 
    elseif strcmp(source.Tag,'67')
        C11 = a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex;        
        a.Material1ME.Mu_complex = real(a.Material1ME.Mu_complex)+1i*str2double(source.String)*1e9;
    end

    a.Material1ME.Lambda_complex = C11-2*a.Material1ME.Mu_complex;
    a.Material1ME.YoungsModulus_complex = a.Material1ME.Mu_complex*((3*a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/(a.Material1ME.Lambda_complex+a.Material1ME.Mu_complex));
    a.Material1ME.PoissonsNumber_complex = a.Material1ME.Lambda_complex/(2*(a.Material1ME.Lambda_complex+a.Material1ME.Mu_complex));
    a.Material1ME.LongitudinalVelocity_complex = conj(sqrt((a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/a.Material1ME.Density));
    a.Material1ME.TransverseVelocity_complex = conj(sqrt((a.Material1ME.Mu_complex)/a.Material1ME.Density));
    a.Material1ME.LongitudinalAttenuation = -2*pi*imag(a.Material1ME.LongitudinalVelocity_complex)/real(a.Material1ME.LongitudinalVelocity_complex);
    a.Material1ME.TransverseAttenuation = -2*pi*imag(a.Material1ME.TransverseVelocity_complex)/real(a.Material1ME.TransverseVelocity_complex);
    a.Material1ME.LongitudinalVelocity = real(a.Material1ME.LongitudinalVelocity_complex)-a.Material1ME.LongitudinalAttenuation*imag(a.Material1ME.LongitudinalVelocity_complex)/(2*pi);
    a.Material1ME.TransverseVelocity = real(a.Material1ME.TransverseVelocity_complex)-a.Material1ME.TransverseAttenuation*imag(a.Material1ME.TransverseVelocity_complex)/(2*pi);
    a.Material1ME.Mu = a.Material1ME.Density*a.Material1ME.TransverseVelocity^2;
    a.Material1ME.Lambda = a.Material1ME.Density*a.Material1ME.LongitudinalVelocity^2-2*a.Material1ME.Mu;
    a.Material1ME.YoungsModulus = a.Material1ME.Mu*((3*a.Material1ME.Lambda+2*a.Material1ME.Mu)/(a.Material1ME.Lambda+a.Material1ME.Mu));
    a.Material1ME.PoissonsNumber = a.Material1ME.Lambda/(2*(a.Material1ME.Lambda+a.Material1ME.Mu));
    a.Material1ME.PlateVelocity = sqrt(a.Material1ME.YoungsModulus/(a.Material1ME.Density*(1-a.Material1ME.PoissonsNumber^2)));
    a.Material1ME.CylinderVelocity = sqrt(a.Material1ME.YoungsModulus/a.Material1ME.Density);
    a.Material1ME.RayleighVelocity = a.Material1ME.TransverseVelocity*((.87+1.12*a.Material1ME.PoissonsNumber)./(1+a.Material1ME.PoissonsNumber));    
    
    a.YoungsModulusUI5.String = real(a.Material1ME.YoungsModulus_complex)/1e9;
    a.YoungsModulusImagUI5.String = imag(a.Material1ME.YoungsModulus_complex)/1e9;
    a.PoissonsNumberUI5.String = real(a.Material1ME.PoissonsNumber_complex);
    a.PoissonsNumberImagUI5.String = imag(a.Material1ME.PoissonsNumber_complex);
    if  a.AttenuationUnit1 == 1 % Np/Lambda
        a.LongitudinalAttenuationUI5.String = a.Material1ME.LongitudinalAttenuation;
        a.TransverseAttenuationUI5.String = a.Material1ME.TransverseAttenuation;
    elseif a.AttenuationUnit1 == 2 % dB/Lambda
        a.LongitudinalAttenuationUI5.String = 20*log10(exp(a.Material1ME.LongitudinalAttenuation));
        a.TransverseAttenuationUI5.String = 20*log10(exp(a.Material1ME.TransverseAttenuation));
    elseif a.AttenuationUnit1 == 3 % Np/m
        a.LongitudinalAttenuationUI5.String = a.Material1ME.LongitudinalAttenuation*a.AtFrequency1*1e3/a.Material1ME.LongitudinalVelocity;
        a.TransverseAttenuationUI5.String = a.Material1ME.TransverseAttenuation*a.AtFrequency1*1e3/a.Material1ME.TransverseVelocity;
    elseif a.AttenuationUnit1 == 4 % dB/m
        a.LongitudinalAttenuationUI5.String = 20*log10(exp(a.Material1ME.LongitudinalAttenuation))*a.AtFrequency1*1e3/a.Material1ME.LongitudinalVelocity;
        a.TransverseAttenuationUI5.String = 20*log10(exp(a.Material1ME.TransverseAttenuation))*a.AtFrequency1*1e3/a.Material1ME.TransverseVelocity;
    end
    a.LongitudinalVelocityUI5.String = a.Material1ME.LongitudinalVelocity;
    a.TransverseVelocityUI5.String = a.Material1ME.TransverseVelocity;
elseif strcmp(source.Tag,'7') || strcmp(source.Tag,'8') % Longitudinal velocity and shear velocity
    if  strcmp(source.Tag,'7')
        a.Material1ME.LongitudinalVelocity = str2double(source.String);
    elseif strcmp(source.Tag,'8')
        a.Material1ME.TransverseVelocity = str2double(source.String);
    end
     
    a.Material1ME.Lambda = a.Material1ME.Density*(a.Material1ME.LongitudinalVelocity^2-2*a.Material1ME.TransverseVelocity^2);
    a.Material1ME.Mu = a.Material1ME.Density*a.Material1ME.TransverseVelocity^2;
    a.Material1ME.YoungsModulus = a.Material1ME.Mu*(3*a.Material1ME.Lambda+2*a.Material1ME.Mu)/(a.Material1ME.Lambda+a.Material1ME.Mu);
    a.Material1ME.PoissonsNumber = a.Material1ME.Lambda/(2*(a.Material1ME.Lambda+a.Material1ME.Mu));
    a.Material1ME.LongitudinalVelocity_complex = a.Material1ME.LongitudinalVelocity/(1+1i*a.Material1ME.LongitudinalAttenuation/(2*pi));
    a.Material1ME.TransverseVelocity_complex = a.Material1ME.TransverseVelocity/(1+1i*a.Material1ME.TransverseAttenuation/(2*pi));
    a.Material1ME.PlateVelocity = sqrt(a.Material1ME.YoungsModulus/(a.Material1ME.Density*(1-a.Material1ME.PoissonsNumber^2)));
    a.Material1ME.CylinderVelocity = sqrt(a.Material1ME.YoungsModulus/a.Material1ME.Density); 
    a.Material1ME.RayleighVelocity = a.Material1ME.TransverseVelocity*((.87+1.12*a.Material1ME.PoissonsNumber)/(1+a.Material1ME.PoissonsNumber));
    Mu_real = a.Material1ME.Density*a.Material1ME.TransverseVelocity^2*4*pi^2*(4*pi^2-a.Material1ME.TransverseAttenuation^2)/(a.Material1ME.TransverseAttenuation^2+4*pi^2)^2;
    Mu_imag = a.Material1ME.Density*a.Material1ME.TransverseVelocity^2*4*pi^2*(4*pi*a.Material1ME.TransverseAttenuation)/(a.Material1ME.TransverseAttenuation^2+4*pi^2)^2;
    Lambda_real = a.Material1ME.Density*a.Material1ME.LongitudinalVelocity^2*4*pi^2*(4*pi^2-a.Material1ME.LongitudinalAttenuation^2)/(a.Material1ME.LongitudinalAttenuation^2+4*pi^2)^2-2*Mu_real;
    Lambda_imag = a.Material1ME.Density*a.Material1ME.LongitudinalVelocity^2*4*pi^2*(4*pi*a.Material1ME.LongitudinalAttenuation)/(a.Material1ME.LongitudinalAttenuation^2+4*pi^2)^2-2*Mu_imag;
    a.Material1ME.Lambda_complex = Lambda_real+1i*Lambda_imag;
    a.Material1ME.Mu_complex = Mu_real+1i*Mu_imag;
    a.Material1ME.YoungsModulus_complex = a.Material1ME.Mu_complex*((3*a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/(a.Material1ME.Lambda_complex+a.Material1ME.Mu_complex));
    a.Material1ME.PoissonsNumber_complex = a.Material1ME.Lambda_complex/(2*(a.Material1ME.Lambda_complex+a.Material1ME.Mu_complex));

    a.YoungsModulusUI5.String = real(a.Material1ME.YoungsModulus_complex)/1e9;
    a.YoungsModulusImagUI5.String = imag(a.Material1ME.YoungsModulus_complex)/1e9;
    a.PoissonsNumberUI5.String = real(a.Material1ME.PoissonsNumber_complex);
    a.PoissonsNumberImagUI5.String = imag(a.Material1ME.PoissonsNumber_complex);
    a.C11IsoUI5.String = real(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/1e9;
    a.C11ImagIsoUI5.String = imag(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/1e9;
    a.C66IsoUI5.String = real(a.Material1ME.Mu_complex)/1e9;
    a.C66ImagIsoUI5.String = imag(a.Material1ME.Mu_complex)/1e9;
elseif strcmp(source.Tag,'9') % Attenuation unit
    a.AttenuationUnit1 = source.Value;
    switch source.Value
    case 1 % Np/Lambda
        a.LongitudinalAttenuationUI5.String = a.Material1ME.LongitudinalAttenuation;
        a.TransverseAttenuationUI5.String = a.Material1ME.TransverseAttenuation;

        a.LongitudinalAttenuationTextUI5.String = ['Longitudinal attenuation (Np/',char(955),')'];
        a.TransverseAttenuationTextUI5.String = ['Shear attenuation (Np/',char(955),')'];
        a.LongitudinalAttenuationTextUI5.Position(3) = 148;
        a.TransverseAttenuationTextUI5.Position(3) = 120;
        a.AtFrequency1UI5.Enable = 'off';
    case 2 % dB/Lambda
        a.LongitudinalAttenuationUI5.String = 20*log10(exp(a.Material1ME.LongitudinalAttenuation));
        a.TransverseAttenuationUI5.String = 20*log10(exp(a.Material1ME.TransverseAttenuation));

        a.LongitudinalAttenuationTextUI5.String = ['Longitudinal attenuation (dB/',char(955),')'];
        a.TransverseAttenuationTextUI5.String = ['Shear attenuation (dB/',char(955),')'];
        a.LongitudinalAttenuationTextUI5.Position(3) = 148;
        a.TransverseAttenuationTextUI5.Position(3) = 120;
        a.AtFrequency1UI5.Enable = 'off';
    case 3 % Np/m
        a.LongitudinalAttenuationUI5.String = a.Material1ME.LongitudinalAttenuation*a.AtFrequency1*1e3/a.Material1ME.LongitudinalVelocity;
        a.TransverseAttenuationUI5.String = a.Material1ME.TransverseAttenuation*a.AtFrequency1*1e3/a.Material1ME.TransverseVelocity;

        a.LongitudinalAttenuationTextUI5.String = 'Longitudinal attenuation (Np/m)';
        a.TransverseAttenuationTextUI5.String = 'Shear attenuation (Np/m)';
        a.LongitudinalAttenuationTextUI5.Position(3) = 150;
        a.TransverseAttenuationTextUI5.Position(3) = 122;
        a.AtFrequency1UI5.Enable = 'on';
    case 4 % dB/m
        a.LongitudinalAttenuationUI5.String = 20*log10(exp(a.Material1ME.LongitudinalAttenuation))*a.AtFrequency1*1e3/a.Material1ME.LongitudinalVelocity;
        a.TransverseAttenuationUI5.String = 20*log10(exp(a.Material1ME.TransverseAttenuation))*a.AtFrequency1*1e3/a.Material1ME.TransverseVelocity;

        a.LongitudinalAttenuationTextUI5.String = 'Longitudinal attenuation (dB/m)';
        a.TransverseAttenuationTextUI5.String = 'Shear attenuation (dB/m)';
        a.LongitudinalAttenuationTextUI5.Position(3) = 150;
        a.TransverseAttenuationTextUI5.Position(3) = 122;
        a.AtFrequency1UI5.Enable = 'on';
    end
elseif strcmp(source.Tag,'10') || strcmp(source.Tag,'11') % Longitudinal attenuation and shear attenuation
    if  strcmp(source.Tag,'10')
        if  a.AttenuationUnit1 == 1 % Np/Lambda
            a.Material1ME.LongitudinalAttenuation = str2double(source.String);
        elseif a.AttenuationUnit1 == 2 % dB/Lambda
            a.Material1ME.LongitudinalAttenuation = log(10)/20*str2double(source.String);
        elseif a.AttenuationUnit1 == 3 % Np/m
            a.Material1ME.LongitudinalAttenuation = str2double(source.String)*a.Material1ME.LongitudinalVelocity/(a.AtFrequency1*1e3);
        elseif a.AttenuationUnit1 == 4 % dB/m
            a.Material1ME.LongitudinalAttenuation = log(10)/20*str2double(source.String)*a.Material1ME.LongitudinalVelocity/(a.AtFrequency1*1e3);
        end
    elseif strcmp(source.Tag,'11')
        if  a.AttenuationUnit1 == 1 % Np/Lambda
            a.Material1ME.TransverseAttenuation = str2double(source.String);
        elseif a.AttenuationUnit1 == 2 % dB/Lambda
            a.Material1ME.TransverseAttenuation = log(10)/20*str2double(source.String);
        elseif a.AttenuationUnit1 == 3 % Np/m
            a.Material1ME.TransverseAttenuation = str2double(source.String)*a.Material1ME.TransverseVelocity/(a.AtFrequency1*1e3);
        elseif a.AttenuationUnit1 == 4 % dB/m
            a.Material1ME.TransverseAttenuation = log(10)/20*str2double(source.String)*a.Material1ME.TransverseVelocity/(a.AtFrequency1*1e3);
        end
    end

    a.Material1ME.LongitudinalVelocity_complex = a.Material1ME.LongitudinalVelocity/(1+1i*a.Material1ME.LongitudinalAttenuation/(2*pi));
    a.Material1ME.TransverseVelocity_complex = a.Material1ME.TransverseVelocity/(1+1i*a.Material1ME.TransverseAttenuation/(2*pi));
    Mu_real = a.Material1ME.Density*a.Material1ME.TransverseVelocity^2*4*pi^2*(4*pi^2-a.Material1ME.TransverseAttenuation^2)/(a.Material1ME.TransverseAttenuation^2+4*pi^2)^2;
    Mu_imag = a.Material1ME.Density*a.Material1ME.TransverseVelocity^2*4*pi^2*(4*pi*a.Material1ME.TransverseAttenuation)/(a.Material1ME.TransverseAttenuation^2+4*pi^2)^2;
    Lambda_real = a.Material1ME.Density*a.Material1ME.LongitudinalVelocity^2*4*pi^2*(4*pi^2-a.Material1ME.LongitudinalAttenuation^2)/(a.Material1ME.LongitudinalAttenuation^2+4*pi^2)^2-2*Mu_real;
    Lambda_imag = a.Material1ME.Density*a.Material1ME.LongitudinalVelocity^2*4*pi^2*(4*pi*a.Material1ME.LongitudinalAttenuation)/(a.Material1ME.LongitudinalAttenuation^2+4*pi^2)^2-2*Mu_imag;
    a.Material1ME.Lambda_complex = Lambda_real+1i*Lambda_imag;
    a.Material1ME.Mu_complex = Mu_real+1i*Mu_imag;
    a.Material1ME.YoungsModulus_complex = a.Material1ME.Mu_complex*((3*a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/(a.Material1ME.Lambda_complex+a.Material1ME.Mu_complex));
    a.Material1ME.PoissonsNumber_complex = a.Material1ME.Lambda_complex/(2*(a.Material1ME.Lambda_complex+a.Material1ME.Mu_complex));

    a.C11IsoUI5.String = real(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/1e9;
    a.C11ImagIsoUI5.String = imag(a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex)/1e9;
    a.C66IsoUI5.String = real(a.Material1ME.Mu_complex)/1e9;
    a.C66ImagIsoUI5.String = imag(a.Material1ME.Mu_complex)/1e9;
    a.YoungsModulusUI5.String = real(a.Material1ME.YoungsModulus_complex)/1e9;
    a.YoungsModulusImagUI5.String = imag(a.Material1ME.YoungsModulus_complex)/1e9;
    a.PoissonsNumberUI5.String = real(a.Material1ME.PoissonsNumber_complex);
    a.PoissonsNumberImagUI5.String = imag(a.Material1ME.PoissonsNumber_complex);
elseif strcmp(source.Tag,'12') % At frequency (isotropic)
    a.Material1ME.LongitudinalAttenuation = a.Material1ME.LongitudinalAttenuation*a.AtFrequency1/str2double(source.String);
    a.Material1ME.TransverseAttenuation = a.Material1ME.TransverseAttenuation*a.AtFrequency1/str2double(source.String);

    a.AtFrequency1 = str2double(source.String);
elseif strcmp(source.Tag,'13') % Material name (isotropic)
    a.MaterialName4A = source.String;
elseif strcmp(source.Tag,'14') % Save (isotropic)
    if  any(strcmp(a.MaterialName4A,fieldnames(a.Materials.Isotropic)))
        Selection = questdlg(['Overwrite ',a.MaterialName4A,'?'],'Confirm','Ok','Cancel','Ok'); 
        switch Selection 
        case 'Ok'
            Save1
        end
    else
        Save1
    end
elseif strcmp(source.Tag,'15') % Delete (isotropic)
    Selection = questdlg(['Delete ',a.MaterialName4A,'?'],'Confirm','Ok','Cancel','Ok');
    switch Selection 
    case 'Ok'
        a.Materials.Isotropic = rmfield(a.Materials.Isotropic,a.MaterialName4A);
        a.Material1UI5.String = fieldnames(a.Materials.Isotropic);
        a.MaterialUI1.String = fieldnames(a.Materials.Isotropic);
        a.Material1UI5.Value = 1;
        a.MaterialUI1.Value = 1;
        if  a.MaterialType_Bulk == 4
            a.MaterialUI8.String = fieldnames(a.Materials.Isotropic);
            a.MaterialUI8.Value = 1;
        end
        if  a.SolidType_Bulk == 1
            a.SolidUI8.String = fieldnames(a.Materials.Isotropic);
            a.SolidUI8.Value = 1;
        end
        E = fieldnames(a.Materials.Isotropic);
        for i = 1:length(fieldnames(a.Materials.Isotropic))
            F = struct2cell(getfield(a.Materials.Isotropic,E{i}))';
            G(i,1:6) = [F(1) F(3:4) F(6) F(8:9)];
        end
        writetable(cell2table(G(:,1:6)),fullfile(a.MaterialListDirectory,'MaterialList_Isotropic.txt'),'delimiter',' ','WriteVariableNames',false)
    end
elseif strcmp(source.Tag,'16') % Class
    a.MaterialTypeME = source.Value;
    a.Material2UI5.Value = 1;
    switch source.Value
    case 1
        a.Material2UI5.String = fieldnames(a.Materials.Orthotropic);
        E = fieldnames(a.Materials.Orthotropic);
        a.Material2ME = getfield(a.Materials.Orthotropic,E{1});        
        
        a.E2UI5.Style = 'edit';
        a.E3UI5.Style = 'edit';
        a.G13UI5.Style = 'edit';
        a.G23UI5.Style = 'edit';
        a.v13UI5.Style = 'edit';
        a.v23UI5.Style = 'edit';
        a.E2UI5.String = real(a.Material2ME.E2)/1e9;
        a.E3UI5.String = real(a.Material2ME.E3)/1e9;
        a.G13UI5.String = real(a.Material2ME.G13)/1e9;
        a.G23UI5.String = real(a.Material2ME.G23)/1e9;
        a.v13UI5.String = real(a.Material2ME.v13);
        a.v23UI5.String = real(a.Material2ME.v23);
        a.E2UI5.BackgroundColor = [1 1 1];
        a.E3UI5.BackgroundColor = [1 1 1];
        a.G13UI5.BackgroundColor = [1 1 1];
        a.G23UI5.BackgroundColor = [1 1 1];
        a.v13UI5.BackgroundColor = [1 1 1];
        a.v23UI5.BackgroundColor = [1 1 1];

        a.E2ImagUI5.Style = 'edit';
        a.E3ImagUI5.Style = 'edit';
        a.G13ImagUI5.Style = 'edit';
        a.G23ImagUI5.Style = 'edit';
        a.v13ImagUI5.Style = 'edit';
        a.v23ImagUI5.Style = 'edit';
        a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
        a.E3ImagUI5.String = imag(a.Material2ME.E3)/1e9;
        a.G13ImagUI5.String = imag(a.Material2ME.G13)/1e9;
        a.G23ImagUI5.String = imag(a.Material2ME.G23)/1e9;
        a.v13ImagUI5.String = imag(a.Material2ME.v13);
        a.v23ImagUI5.String = imag(a.Material2ME.v23);
        a.E2ImagUI5.BackgroundColor = [1 1 1];
        a.E3ImagUI5.BackgroundColor = [1 1 1];
        a.G13ImagUI5.BackgroundColor = [1 1 1];
        a.G23ImagUI5.BackgroundColor = [1 1 1];
        a.v13ImagUI5.BackgroundColor = [1 1 1];
        a.v23ImagUI5.BackgroundColor = [1 1 1];
        
        a.C13UI5.Style = 'edit';
        a.C22UI5.Style = 'edit';
        a.C23UI5.Style = 'edit';
        a.C33UI5.Style = 'edit';
        a.C44UI5.Style = 'edit';
        a.C55UI5.Style = 'edit';
        a.C13UI5.String = real(a.Material2ME.C(1,3))/1e9;
        a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
        a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;
        a.C33UI5.String = real(a.Material2ME.C(3,3))/1e9;
        a.C44UI5.String = real(a.Material2ME.C(4,4))/1e9;
        a.C55UI5.String = real(a.Material2ME.C(5,5))/1e9;
        a.C13UI5.BackgroundColor = [1 1 1];
        a.C22UI5.BackgroundColor = [1 1 1];
        a.C23UI5.BackgroundColor = [1 1 1];
        a.C33UI5.BackgroundColor = [1 1 1];
        a.C44UI5.BackgroundColor = [1 1 1];
        a.C55UI5.BackgroundColor = [1 1 1];

        a.C13ImagUI5.Style = 'edit';
        a.C22ImagUI5.Style = 'edit';
        a.C23ImagUI5.Style = 'edit';
        a.C33ImagUI5.Style = 'edit';
        a.C44ImagUI5.Style = 'edit';
        a.C55ImagUI5.Style = 'edit';
        a.C13ImagUI5.String = imag(a.Material2ME.C(1,3))/1e9;
        a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
        a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
        a.C33ImagUI5.String = imag(a.Material2ME.C(3,3))/1e9;
        a.C44ImagUI5.String = imag(a.Material2ME.C(4,4))/1e9;
        a.C55ImagUI5.String = imag(a.Material2ME.C(5,5))/1e9;
        a.C13ImagUI5.BackgroundColor = [1 1 1];
        a.C22ImagUI5.BackgroundColor = [1 1 1];
        a.C23ImagUI5.BackgroundColor = [1 1 1];
        a.C33ImagUI5.BackgroundColor = [1 1 1];
        a.C44ImagUI5.BackgroundColor = [1 1 1];
        a.C55ImagUI5.BackgroundColor = [1 1 1];
    case 2
        a.Material2UI5.String = fieldnames(a.Materials.TransverselyIsotropic);
        E = fieldnames(a.Materials.TransverselyIsotropic);
        a.Material2ME = getfield(a.Materials.TransverselyIsotropic,E{1});
        
        a.E2UI5.Style = 'edit';
        a.E3UI5.Style = 'text';
        a.G13UI5.Style = 'text';
        a.G23UI5.Style = 'text';
        a.v13UI5.Style = 'text';
        a.v23UI5.Style = 'edit';
        a.E2UI5.String = real(a.Material2ME.E2)/1e9;
        a.E3UI5.String = 'E2';
        a.G13UI5.String = 'G12';
        a.G23UI5.String = '.5E2/(1+v23)';
        a.v13UI5.String = 'v12';
        a.v23UI5.String = real(a.Material2ME.v23);
        a.E2UI5.BackgroundColor = [1 1 1];
        a.E3UI5.BackgroundColor = [.88 .88 .88];
        a.G13UI5.BackgroundColor = [.88 .88 .88];
        a.G23UI5.BackgroundColor = [.88 .88 .88];
        a.v13UI5.BackgroundColor = [.88 .88 .88];
        a.v23UI5.BackgroundColor = [1 1 1];

        a.E2ImagUI5.Style = 'edit';
        a.E3ImagUI5.Style = 'text';
        a.G13ImagUI5.Style = 'text';
        a.G23ImagUI5.Style = 'text';
        a.v13ImagUI5.Style = 'text';
        a.v23ImagUI5.Style = 'edit';
        a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
        a.E3ImagUI5.String = 'E2';
        a.G13ImagUI5.String = 'G12';
        a.G23ImagUI5.String = '.5E2/(1+v23)';
        a.v13ImagUI5.String = 'v12';
        a.v23ImagUI5.String = imag(a.Material2ME.v23);
        a.E2ImagUI5.BackgroundColor = [1 1 1];
        a.E3ImagUI5.BackgroundColor = [.88 .88 .88];
        a.G13ImagUI5.BackgroundColor = [.88 .88 .88];
        a.G23ImagUI5.BackgroundColor = [.88 .88 .88];
        a.v13ImagUI5.BackgroundColor = [.88 .88 .88];
        a.v23ImagUI5.BackgroundColor = [1 1 1];
        
        a.C13UI5.Style = 'text';
        a.C22UI5.Style = 'edit';
        a.C23UI5.Style = 'edit';
        a.C33UI5.Style = 'text';
        a.C44UI5.Style = 'text';
        a.C55UI5.Style = 'text';
        a.C13UI5.String = 'C12';
        a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
        a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;
        a.C33UI5.String = 'C22';
        a.C44UI5.String = '.5(C22-C23)';
        a.C55UI5.String = 'C66';
        a.C13UI5.BackgroundColor = [.88 .88 .88];
        a.C22UI5.BackgroundColor = [1 1 1];
        a.C23UI5.BackgroundColor = [1 1 1];
        a.C33UI5.BackgroundColor = [.88 .88 .88];
        a.C44UI5.BackgroundColor = [.88 .88 .88];
        a.C55UI5.BackgroundColor = [.88 .88 .88];

        a.C13ImagUI5.Style = 'text';
        a.C22ImagUI5.Style = 'edit';
        a.C23ImagUI5.Style = 'edit';
        a.C33ImagUI5.Style = 'text';
        a.C44ImagUI5.Style = 'text';
        a.C55ImagUI5.Style = 'text';
        a.C13ImagUI5.String = 'C12';
        a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
        a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
        a.C33ImagUI5.String = 'C22';
        a.C44ImagUI5.String = '.5(C22-C23)';
        a.C55ImagUI5.String = 'C66';
        a.C13ImagUI5.BackgroundColor = [.88 .88 .88];
        a.C22ImagUI5.BackgroundColor = [1 1 1];
        a.C23ImagUI5.BackgroundColor = [1 1 1];
        a.C33ImagUI5.BackgroundColor = [.88 .88 .88];
        a.C44ImagUI5.BackgroundColor = [.88 .88 .88];
        a.C55ImagUI5.BackgroundColor = [.88 .88 .88];
    case 3
        a.Material2UI5.String = fieldnames(a.Materials.Cubic);
        E = fieldnames(a.Materials.Cubic);
        a.Material2ME = getfield(a.Materials.Cubic,E{1});

        a.E2UI5.Style = 'text';
        a.E3UI5.Style = 'text';
        a.G13UI5.Style = 'text';
        a.G23UI5.Style = 'text';
        a.v13UI5.Style = 'text';
        a.v23UI5.Style = 'text';
        a.E2UI5.String = 'E1';
        a.E3UI5.String = 'E1';
        a.G13UI5.String = 'G12';
        a.G23UI5.String = 'G12';
        a.v13UI5.String = 'v12';
        a.v23UI5.String = 'v12';
        a.E2UI5.BackgroundColor = [.88 .88 .88];
        a.E3UI5.BackgroundColor = [.88 .88 .88];
        a.G13UI5.BackgroundColor = [.88 .88 .88];
        a.G23UI5.BackgroundColor = [.88 .88 .88];
        a.v13UI5.BackgroundColor = [.88 .88 .88];
        a.v23UI5.BackgroundColor = [.88 .88 .88];

        a.E2ImagUI5.Style = 'text';
        a.E3ImagUI5.Style = 'text';
        a.G13ImagUI5.Style = 'text';
        a.G23ImagUI5.Style = 'text';
        a.v13ImagUI5.Style = 'text';
        a.v23ImagUI5.Style = 'text';
        a.E2ImagUI5.String = 'E1';
        a.E3ImagUI5.String = 'E1';
        a.G13ImagUI5.String = 'G12';
        a.G23ImagUI5.String = 'G12';
        a.v13ImagUI5.String = 'v12';
        a.v23ImagUI5.String = 'v12';
        a.E2ImagUI5.BackgroundColor = [.88 .88 .88];
        a.E3ImagUI5.BackgroundColor = [.88 .88 .88];
        a.G13ImagUI5.BackgroundColor = [.88 .88 .88];
        a.G23ImagUI5.BackgroundColor = [.88 .88 .88];
        a.v13ImagUI5.BackgroundColor = [.88 .88 .88];
        a.v23ImagUI5.BackgroundColor = [.88 .88 .88];        
        
        a.C13UI5.Style = 'text';
        a.C22UI5.Style = 'text';
        a.C23UI5.Style = 'text';
        a.C33UI5.Style = 'text';
        a.C44UI5.Style = 'text';
        a.C55UI5.Style = 'text';
        a.C13UI5.String = 'C12';
        a.C22UI5.String = 'C11';
        a.C23UI5.String = 'C12';
        a.C33UI5.String = 'C11';
        a.C44UI5.String = 'C66';
        a.C55UI5.String = 'C66';
        a.C13UI5.BackgroundColor = [.88 .88 .88];
        a.C22UI5.BackgroundColor = [.88 .88 .88];
        a.C23UI5.BackgroundColor = [.88 .88 .88];
        a.C33UI5.BackgroundColor = [.88 .88 .88];
        a.C44UI5.BackgroundColor = [.88 .88 .88];
        a.C55UI5.BackgroundColor = [.88 .88 .88];

        a.C13ImagUI5.Style = 'text';
        a.C22ImagUI5.Style = 'text';
        a.C23ImagUI5.Style = 'text';
        a.C33ImagUI5.Style = 'text';
        a.C44ImagUI5.Style = 'text';
        a.C55ImagUI5.Style = 'text';
        a.C13ImagUI5.String = 'C12';
        a.C22ImagUI5.String = 'C11';
        a.C23ImagUI5.String = 'C12';
        a.C33ImagUI5.String = 'C11';
        a.C44ImagUI5.String = 'C66';
        a.C55ImagUI5.String = 'C66';
        a.C13ImagUI5.BackgroundColor = [.88 .88 .88];
        a.C22ImagUI5.BackgroundColor = [.88 .88 .88];
        a.C23ImagUI5.BackgroundColor = [.88 .88 .88];
        a.C33ImagUI5.BackgroundColor = [.88 .88 .88];
        a.C44ImagUI5.BackgroundColor = [.88 .88 .88];
        a.C55ImagUI5.BackgroundColor = [.88 .88 .88];
    end
    a.Density2UI5.String = a.Material2ME.Density;

    a.E1UI5.String = real(a.Material2ME.E1)/1e9;
    a.G12UI5.String = real(a.Material2ME.G12)/1e9;
    a.v12UI5.String = real(a.Material2ME.v12);
    a.C11UI5.String = real(a.Material2ME.C(1,1))/1e9;
    a.C12UI5.String = real(a.Material2ME.C(1,2))/1e9;
    a.C66UI5.String = real(a.Material2ME.C(6,6))/1e9;

    a.E1ImagUI5.String = imag(a.Material2ME.E1)/1e9;
    a.G12ImagUI5.String = imag(a.Material2ME.G12)/1e9;
    a.v12ImagUI5.String = imag(a.Material2ME.v12);
    a.C11ImagUI5.String = imag(a.Material2ME.C(1,1))/1e9;
    a.C12ImagUI5.String = imag(a.Material2ME.C(1,2))/1e9;
    a.C66ImagUI5.String = imag(a.Material2ME.C(6,6))/1e9;

    a.MaterialName4B = a.Material2ME.Name;
    a.Name2UI5.String = a.MaterialName4B;
    
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.FastShearVelocity_1UI5.String = a.Material2ME.FastShearVelocity_1;
    a.SlowShearVelocity_1UI5.String = a.Material2ME.SlowShearVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.FastShearVelocity_2UI5.String = a.Material2ME.FastShearVelocity_2;
    a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
    a.FastShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_3;
    a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;
elseif strcmp(source.Tag,'17') % Material (anisotropic)
    if  a.MaterialTypeME == 1
        a.Material2ME = getfield(a.Materials.Orthotropic,cell2mat(source.String(source.Value)));
    elseif a.MaterialTypeME == 2
        a.Material2ME = getfield(a.Materials.TransverselyIsotropic,cell2mat(source.String(source.Value)));
    elseif a.MaterialTypeME == 3
        a.Material2ME = getfield(a.Materials.Cubic,cell2mat(source.String(source.Value)));
    end
    a.Density2UI5.String = a.Material2ME.Density;

    a.E1UI5.String = real(a.Material2ME.E1)/1e9;
    a.G12UI5.String = real(a.Material2ME.G12)/1e9;
    a.v12UI5.String = real(a.Material2ME.v12);
    a.C11UI5.String = real(a.Material2ME.C(1,1))/1e9;
    a.C12UI5.String = real(a.Material2ME.C(1,2))/1e9;
    a.C66UI5.String = real(a.Material2ME.C(6,6))/1e9;

    a.E1ImagUI5.String = imag(a.Material2ME.E1)/1e9;
    a.G12ImagUI5.String = imag(a.Material2ME.G12)/1e9;
    a.v12ImagUI5.String = imag(a.Material2ME.v12);
    a.C11ImagUI5.String = imag(a.Material2ME.C(1,1))/1e9;
    a.C12ImagUI5.String = imag(a.Material2ME.C(1,2))/1e9;
    a.C66ImagUI5.String = imag(a.Material2ME.C(6,6))/1e9;
    
    a.MaterialName4B = a.Material2ME.Name;
    a.Name2UI5.String = a.MaterialName4B;
    if  a.MaterialTypeME == 1
        a.E2UI5.String = real(a.Material2ME.E2)/1e9;
        a.E3UI5.String = real(a.Material2ME.E3)/1e9;
        a.G13UI5.String = real(a.Material2ME.G13)/1e9;
        a.G23UI5.String = real(a.Material2ME.G23)/1e9;
        a.v13UI5.String = real(a.Material2ME.v13);
        a.v23UI5.String = real(a.Material2ME.v23);
        a.C13UI5.String = real(a.Material2ME.C(1,3))/1e9;
        a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;
        a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
        a.C33UI5.String = real(a.Material2ME.C(3,3))/1e9;
        a.C44UI5.String = real(a.Material2ME.C(4,4))/1e9;
        a.C55UI5.String = real(a.Material2ME.C(5,5))/1e9;

        a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
        a.E3ImagUI5.String = imag(a.Material2ME.E3)/1e9;
        a.G13ImagUI5.String = imag(a.Material2ME.G13)/1e9;
        a.G23ImagUI5.String = imag(a.Material2ME.G23)/1e9;
        a.v13ImagUI5.String = imag(a.Material2ME.v13);
        a.v23ImagUI5.String = imag(a.Material2ME.v23);
        a.C13ImagUI5.String = imag(a.Material2ME.C(1,3))/1e9;
        a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
        a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
        a.C33ImagUI5.String = imag(a.Material2ME.C(3,3))/1e9;
        a.C44ImagUI5.String = imag(a.Material2ME.C(4,4))/1e9;
        a.C55ImagUI5.String = imag(a.Material2ME.C(5,5))/1e9;
    elseif  a.MaterialTypeME == 2
        a.E2UI5.String = real(a.Material2ME.E2)/1e9;
        a.v23UI5.String = real(a.Material2ME.v23);
        a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;
        a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;

        a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
        a.v23ImagUI5.String = imag(a.Material2ME.v23);
        a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
        a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
    end
    
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.FastShearVelocity_1UI5.String = a.Material2ME.FastShearVelocity_1;
    a.SlowShearVelocity_1UI5.String = a.Material2ME.SlowShearVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.FastShearVelocity_2UI5.String = a.Material2ME.FastShearVelocity_2;
    a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
    a.FastShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_3;
    a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;
elseif strcmp(source.Tag,'18') % Density (anisotropic)
    a.Material2ME.Density = str2double(source.String);
    
    a.Material2ME.LongitudinalVelocity_1 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
    if  a.MaterialTypeME == 1
        a.Material2ME.FastShearVelocity_1 = sqrt(real(a.Material2ME.C(5,5))/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_1 = sqrt(real(a.Material2ME.C(6,6))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.FastShearVelocity_2 = sqrt(real(a.Material2ME.C(6,6))/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_2 = sqrt(real(a.Material2ME.C(4,4))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = sqrt(real(a.Material2ME.C(3,3))/a.Material2ME.Density);
        a.Material2ME.FastShearVelocity_3 = sqrt(real(a.Material2ME.C(5,5))/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_3 = sqrt(real(a.Material2ME.C(4,4))/a.Material2ME.Density);
    elseif a.MaterialTypeME == 2
        a.Material2ME.FastShearVelocity_1 = sqrt(real(a.Material2ME.C(5,5))/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_1 = sqrt(real(a.Material2ME.C(6,6))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.FastShearVelocity_2 = sqrt(real(a.Material2ME.C(6,6))/a.Material2ME.Density); 
        a.Material2ME.SlowShearVelocity_2 = sqrt(.5*(real(a.Material2ME.C(2,2))-real(a.Material2ME.C(2,3)))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_2;
        a.Material2ME.FastShearVelocity_3 = a.Material2ME.FastShearVelocity_2;
        a.Material2ME.SlowShearVelocity_3 = a.Material2ME.SlowShearVelocity_2;
    elseif a.MaterialTypeME == 3
        a.Material2ME.FastShearVelocity_1 = sqrt(real(a.Material2ME.G12)/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_1 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.LongitudinalVelocity_2 = a.Material2ME.LongitudinalVelocity_1;
        a.Material2ME.FastShearVelocity_2 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.SlowShearVelocity_2 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_1;
        a.Material2ME.FastShearVelocity_3 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.SlowShearVelocity_3 = a.Material2ME.FastShearVelocity_1;
    end
    
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.FastShearVelocity_1UI5.String = a.Material2ME.FastShearVelocity_1;
    a.SlowShearVelocity_1UI5.String = a.Material2ME.SlowShearVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.FastShearVelocity_2UI5.String = a.Material2ME.FastShearVelocity_2;
    a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
    a.FastShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_3;
    a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;
elseif strcmp(source.Tag,'19') || strcmp(source.Tag,'28') % E1
    if  strcmp(source.Tag,'19')
        a.Material2ME.E1 = str2double(source.String)*1e9+1i*imag(a.Material2ME.E1);
    elseif strcmp(source.Tag,'28')
        a.Material2ME.E1 = real(a.Material2ME.E1)+1i*str2double(source.String)*1e9;    
    end

    a.Material2ME.v21 = a.Material2ME.E2*a.Material2ME.v12/a.Material2ME.E1;
    if  a.MaterialTypeME == 1
        a.Material2ME.v31 = a.Material2ME.E3*a.Material2ME.v13/a.Material2ME.E1;
        a.Material2ME.v32 = a.Material2ME.E3*a.Material2ME.v23/a.Material2ME.E2;
        D = (1-a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v13*a.Material2ME.v31-a.Material2ME.v23*a.Material2ME.v32-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v31);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v23*a.Material2ME.v32)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v13*a.Material2ME.v32+a.Material2ME.v12)*a.Material2ME.E2/D; 
        a.Material2ME.C(1,3) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v13)*a.Material2ME.E3/D;
        a.Material2ME.C(2,2) = (1-a.Material2ME.v13*a.Material2ME.v31)*a.Material2ME.E2/D;
        a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v13+a.Material2ME.v23)*a.Material2ME.E3/D;
        a.Material2ME.C(3,3) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E3/D;
        
        a.C13UI5.String = real(a.Material2ME.C(1,3))/1e9;
        a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
        a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;
        a.C33UI5.String = real(a.Material2ME.C(3,3))/1e9;

        a.C13ImagUI5.String = imag(a.Material2ME.C(1,3))/1e9;
        a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
        a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
        a.C33ImagUI5.String = imag(a.Material2ME.C(3,3))/1e9;

        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = sqrt(real(a.Material2ME.C(3,3))/a.Material2ME.Density);
    elseif a.MaterialTypeME == 2
        D = (1-2*a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v23^2-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v21);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v23^2)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v12)*a.Material2ME.E2/D;
        a.Material2ME.C(2,2) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E2/D;
        a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v12+a.Material2ME.v23)*a.Material2ME.E2/D;
        
        a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
        a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;

        a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
        a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;

        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_2;
    elseif a.MaterialTypeME == 3
        D = (1-3*a.Material2ME.v12^2-2*a.Material2ME.v12^3);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v12^2)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v12^2+a.Material2ME.v12)*a.Material2ME.E1/D;
        
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_2;
    end
    a.C11UI5.String = real(a.Material2ME.C(1,1))/1e9;
    a.C12UI5.String = real(a.Material2ME.C(1,2))/1e9;

    a.C11ImagUI5.String = imag(a.Material2ME.C(1,1))/1e9;
    a.C12ImagUI5.String = imag(a.Material2ME.C(1,2))/1e9;

    a.Material2ME.LongitudinalVelocity_1 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
elseif strcmp(source.Tag,'20') || strcmp(source.Tag,'29') % E2
    if  strcmp(source.Tag,'20')
        a.Material2ME.E2 = str2double(source.String)*1e9+1i*imag(a.Material2ME.E2);
    elseif strcmp(source.Tag,'29')
        a.Material2ME.E2 = real(a.Material2ME.E2)+1i*str2double(source.String)*1e9;    
    end

    a.Material2ME.v21 = a.Material2ME.E2*a.Material2ME.v12/a.Material2ME.E1;
    if  a.MaterialTypeME == 1
        a.Material2ME.v31 = a.Material2ME.E3*a.Material2ME.v13/a.Material2ME.E1;
        a.Material2ME.v32 = a.Material2ME.E3*a.Material2ME.v23/a.Material2ME.E2;
        D = (1-a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v13*a.Material2ME.v31-a.Material2ME.v23*a.Material2ME.v32-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v31);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v23*a.Material2ME.v32)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v13*a.Material2ME.v32+a.Material2ME.v12)*a.Material2ME.E2/D; 
        a.Material2ME.C(1,3) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v13)*a.Material2ME.E3/D;
        a.Material2ME.C(2,2) = (1-a.Material2ME.v13*a.Material2ME.v31)*a.Material2ME.E2/D;
        a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v13+a.Material2ME.v23)*a.Material2ME.E3/D;
        a.Material2ME.C(3,3) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E3/D;
        
        a.C13UI5.String = real(a.Material2ME.C(1,3))/1e9;
        a.C33UI5.String = real(a.Material2ME.C(3,3))/1e9;

        a.C13ImagUI5.String = imag(a.Material2ME.C(1,3))/1e9;
        a.C33ImagUI5.String = imag(a.Material2ME.C(3,3))/1e9;
        
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = sqrt(real(a.Material2ME.C(3,3))/a.Material2ME.Density);
    elseif a.MaterialTypeME == 2
        a.Material2ME.E3 = a.Material2ME.E2;
        a.Material2ME.G23 = .5*a.Material2ME.E2/(1+a.Material2ME.v23);
        
        D = (1-2*a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v23^2-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v21);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v23^2)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v12)*a.Material2ME.E2/D;
        a.Material2ME.C(2,2) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E2/D;
        a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v12+a.Material2ME.v23)*a.Material2ME.E2/D;
        
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_2;
        a.Material2ME.SlowShearVelocity_2 = sqrt(.5*(real(a.Material2ME.C(2,2))-real(a.Material2ME.C(2,3)))/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_3 = a.Material2ME.SlowShearVelocity_2;
        a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
        a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;
    end
    a.C11UI5.String = real(a.Material2ME.C(1,1))/1e9;
    a.C12UI5.String = real(a.Material2ME.C(1,2))/1e9;
    a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
    a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;

    a.C11ImagUI5.String = imag(a.Material2ME.C(1,1))/1e9;
    a.C12ImagUI5.String = imag(a.Material2ME.C(1,2))/1e9;
    a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
    a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
    
    a.Material2ME.LongitudinalVelocity_1 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
elseif strcmp(source.Tag,'21') || strcmp(source.Tag,'30') % E3
    if  strcmp(source.Tag,'21')
        a.Material2ME.E3 = str2double(source.String)*1e9+1i*imag(a.Material2ME.E3);
    elseif strcmp(source.Tag,'30')
        a.Material2ME.E3 = real(a.Material2ME.E3)+1i*str2double(source.String)*1e9;    
    end

    a.Material2ME.v21 = a.Material2ME.E2*a.Material2ME.v12/a.Material2ME.E1;
    a.Material2ME.v31 = a.Material2ME.E3*a.Material2ME.v13/a.Material2ME.E1;
    a.Material2ME.v32 = a.Material2ME.E3*a.Material2ME.v23/a.Material2ME.E2;
    D = (1-a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v13*a.Material2ME.v31-a.Material2ME.v23*a.Material2ME.v32-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v31);
    a.Material2ME.C(1,1) = (1-a.Material2ME.v23*a.Material2ME.v32)*a.Material2ME.E1/D;
    a.Material2ME.C(1,2) = (a.Material2ME.v13*a.Material2ME.v32+a.Material2ME.v12)*a.Material2ME.E2/D; 
    a.Material2ME.C(1,3) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v13)*a.Material2ME.E3/D;
    a.Material2ME.C(2,2) = (1-a.Material2ME.v13*a.Material2ME.v31)*a.Material2ME.E2/D;
    a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v13+a.Material2ME.v23)*a.Material2ME.E3/D;
    a.Material2ME.C(3,3) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E3/D;

    a.C11UI5.String = real(a.Material2ME.C(1,1))/1e9;
    a.C12UI5.String = real(a.Material2ME.C(1,2))/1e9;
    a.C13UI5.String = real(a.Material2ME.C(1,3))/1e9;
    a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
    a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;
    a.C33UI5.String = real(a.Material2ME.C(3,3))/1e9;

    a.C11ImagUI5.String = imag(a.Material2ME.C(1,1))/1e9;
    a.C12ImagUI5.String = imag(a.Material2ME.C(1,2))/1e9;
    a.C13ImagUI5.String = imag(a.Material2ME.C(1,3))/1e9;
    a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
    a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
    a.C33ImagUI5.String = imag(a.Material2ME.C(3,3))/1e9;

    a.Material2ME.LongitudinalVelocity_1 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
    a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
    a.Material2ME.LongitudinalVelocity_3 = sqrt(real(a.Material2ME.C(3,3))/a.Material2ME.Density);
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
elseif strcmp(source.Tag,'22') || strcmp(source.Tag,'31') % G12
    if  strcmp(source.Tag,'22')
        a.Material2ME.G12 = str2double(source.String)*1e9+1i*imag(a.Material2ME.G12);
    elseif strcmp(source.Tag,'31')
        a.Material2ME.G12 = real(a.Material2ME.G12)+1i*str2double(source.String)*1e9;    
    end
    
    a.Material2ME.C(6,6) = a.Material2ME.G12;
    a.C66UI5.String = real(a.Material2ME.C(6,6))/1e9;
    a.C66ImagUI5.String = imag(a.Material2ME.C(6,6))/1e9;
    if  a.MaterialTypeME == 1
        a.Material2ME.SlowShearVelocity_1 = sqrt(real(a.Material2ME.G12)/a.Material2ME.Density);
        a.Material2ME.FastShearVelocity_2 = a.Material2ME.SlowShearVelocity_1;
    elseif a.MaterialTypeME == 2
        a.Material2ME.G13 = a.Material2ME.G12;
        a.Material2ME.C(5,5) = a.Material2ME.G12;
        
        a.Material2ME.FastShearVelocity_1 = sqrt(real(a.Material2ME.G12)/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_1 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.FastShearVelocity_2 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.FastShearVelocity_3 = a.Material2ME.FastShearVelocity_1;
        a.FastShearVelocity_1UI5.String = a.Material2ME.FastShearVelocity_1;
        a.FastShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_3;
    elseif a.MaterialTypeME == 3
        a.Material2ME.FastShearVelocity_1 = sqrt(real(a.Material2ME.G12)/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_1 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.FastShearVelocity_2 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.SlowShearVelocity_2 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.FastShearVelocity_3 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.SlowShearVelocity_3 = a.Material2ME.FastShearVelocity_1;
        a.FastShearVelocity_1UI5.String = a.Material2ME.FastShearVelocity_1;
        a.SlowShearVelocity_2UI5.String = a.Material2ME.FastShearVelocity_1;
        a.FastShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_1;
        a.SlowShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_1;
    end
    a.SlowShearVelocity_1UI5.String = a.Material2ME.SlowShearVelocity_1;
    a.FastShearVelocity_2UI5.String = a.Material2ME.FastShearVelocity_2;   
elseif strcmp(source.Tag,'23') || strcmp(source.Tag,'32') % G13
    if  strcmp(source.Tag,'23')
        a.Material2ME.G13 = str2double(source.String)*1e9+1i*imag(a.Material2ME.G13);
    elseif strcmp(source.Tag,'32')
        a.Material2ME.G13 = real(a.Material2ME.G13)+1i*str2double(source.String)*1e9;    
    end
    
    a.Material2ME.C(5,5) = a.Material2ME.G13;
    a.C55UI5.String = real(a.Material2ME.C(5,5))/1e9;
    a.C55ImagUI5.String = imag(a.Material2ME.C(5,5))/1e9;
    
    a.Material2ME.FastShearVelocity_1 = sqrt(real(a.Material2ME.G13)/a.Material2ME.Density);
    a.Material2ME.FastShearVelocity_3 = a.Material2ME.FastShearVelocity_1;
    a.FastShearVelocity_1UI5.String = a.Material2ME.FastShearVelocity_1;
    a.FastShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_3;      
elseif strcmp(source.Tag,'24') || strcmp(source.Tag,'33') % G23
    if  strcmp(source.Tag,'24')
        a.Material2ME.G23 = str2double(source.String)*1e9+1i*imag(a.Material2ME.G23);
    elseif strcmp(source.Tag,'33')
        a.Material2ME.G23 = real(a.Material2ME.G23)+1i*str2double(source.String)*1e9;    
    end
    
    a.Material2ME.C(4,4) = a.Material2ME.G23;
    a.C44UI5.String = real(a.Material2ME.C(4,4))/1e9;
    a.C44ImagUI5.String = imag(a.Material2ME.C(4,4))/1e9;
    
    a.Material2ME.SlowShearVelocity_2 = sqrt(real(a.Material2ME.G23)/a.Material2ME.Density);
    a.Material2ME.SlowShearVelocity_3 = a.Material2ME.SlowShearVelocity_2;
    a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
    a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;  
elseif strcmp(source.Tag,'25') || strcmp(source.Tag,'34') % v12
    if  strcmp(source.Tag,'25')
        a.Material2ME.v12 = str2double(source.String)+1i*imag(a.Material2ME.v12);
    elseif strcmp(source.Tag,'34')
        a.Material2ME.v12 = real(a.Material2ME.v12)+1i*str2double(source.String);    
    end
    
    a.Material2ME.v21 = a.Material2ME.E2*a.Material2ME.v12/a.Material2ME.E1;
    if  a.MaterialTypeME == 1
        a.Material2ME.v31 = a.Material2ME.E3*a.Material2ME.v13/a.Material2ME.E1;
        a.Material2ME.v32 = a.Material2ME.E3*a.Material2ME.v23/a.Material2ME.E2;
        D = (1-a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v13*a.Material2ME.v31-a.Material2ME.v23*a.Material2ME.v32-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v31);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v23*a.Material2ME.v32)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v13*a.Material2ME.v32+a.Material2ME.v12)*a.Material2ME.E2/D; 
        a.Material2ME.C(1,3) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v13)*a.Material2ME.E3/D;
        a.Material2ME.C(2,2) = (1-a.Material2ME.v13*a.Material2ME.v31)*a.Material2ME.E2/D;
        a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v13+a.Material2ME.v23)*a.Material2ME.E3/D;
        a.Material2ME.C(3,3) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E3/D;
        
        a.C13UI5.String = real(a.Material2ME.C(1,3))/1e9;
        a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
        a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;
        a.C33UI5.String = real(a.Material2ME.C(3,3))/1e9;

        a.C13ImagUI5.String = imag(a.Material2ME.C(1,3))/1e9;
        a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
        a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
        a.C33ImagUI5.String = imag(a.Material2ME.C(3,3))/1e9;
        
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = sqrt(real(a.Material2ME.C(3,3))/a.Material2ME.Density);
    elseif a.MaterialTypeME == 2
        a.Material2ME.v13 = a.Material2ME.v12;
        
        D = (1-2*a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v23^2-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v21);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v23^2)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v12)*a.Material2ME.E2/D;
        a.Material2ME.C(2,2) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E2/D;
        a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v12+a.Material2ME.v23)*a.Material2ME.E2/D;

        a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
        a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;

        a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
        a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;

        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_2;
    elseif a.MaterialTypeME == 3
        D = (1-3*a.Material2ME.v12^2-2*a.Material2ME.v12^3);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v12^2)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v12^2+a.Material2ME.v12)*a.Material2ME.E1/D;
        
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_2;
    end
    a.C11UI5.String = real(a.Material2ME.C(1,1))/1e9;
    a.C12UI5.String = real(a.Material2ME.C(1,2))/1e9;

    a.C11ImagUI5.String = imag(a.Material2ME.C(1,1))/1e9;
    a.C12ImagUI5.String = imag(a.Material2ME.C(1,2))/1e9;
    
    a.Material2ME.LongitudinalVelocity_1 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;   
elseif strcmp(source.Tag,'26') || strcmp(source.Tag,'35') % v13
    if  strcmp(source.Tag,'26')
        a.Material2ME.v13 = str2double(source.String)+1i*imag(a.Material2ME.v13);
    elseif strcmp(source.Tag,'35')
        a.Material2ME.v13 = real(a.Material2ME.v13)+1i*str2double(source.String);    
    end
    
    a.Material2ME.v21 = a.Material2ME.E2*a.Material2ME.v12/a.Material2ME.E1;
    a.Material2ME.v31 = a.Material2ME.E3*a.Material2ME.v13/a.Material2ME.E1;
    a.Material2ME.v32 = a.Material2ME.E3*a.Material2ME.v23/a.Material2ME.E2;
    D = (1-a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v13*a.Material2ME.v31-a.Material2ME.v23*a.Material2ME.v32-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v31);
    a.Material2ME.C(1,1) = (1-a.Material2ME.v23*a.Material2ME.v32)*a.Material2ME.E1/D;
    a.Material2ME.C(1,2) = (a.Material2ME.v13*a.Material2ME.v32+a.Material2ME.v12)*a.Material2ME.E2/D;
    a.Material2ME.C(1,3) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v13)*a.Material2ME.E3/D;
    a.Material2ME.C(2,2) = (1-a.Material2ME.v13*a.Material2ME.v31)*a.Material2ME.E2/D;
    a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v13+a.Material2ME.v23)*a.Material2ME.E3/D;
    a.Material2ME.C(3,3) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E3/D;
    
    a.C11UI5.String = real(a.Material2ME.C(1,1))/1e9;
    a.C12UI5.String = real(a.Material2ME.C(1,2))/1e9;
    a.C13UI5.String = real(a.Material2ME.C(1,3))/1e9;
    a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
    a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;
    a.C33UI5.String = real(a.Material2ME.C(3,3))/1e9;

    a.C11ImagUI5.String = imag(a.Material2ME.C(1,1))/1e9;
    a.C12ImagUI5.String = imag(a.Material2ME.C(1,2))/1e9;
    a.C13ImagUI5.String = imag(a.Material2ME.C(1,3))/1e9;
    a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
    a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
    a.C33ImagUI5.String = imag(a.Material2ME.C(3,3))/1e9;
    
    a.Material2ME.LongitudinalVelocity_1 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
    a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
    a.Material2ME.LongitudinalVelocity_3 = sqrt(real(a.Material2ME.C(3,3))/a.Material2ME.Density);
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
elseif strcmp(source.Tag,'27') || strcmp(source.Tag,'36') % v23
    if  strcmp(source.Tag,'27')
        a.Material2ME.v23 = str2double(source.String)+1i*imag(a.Material2ME.v23);
    elseif strcmp(source.Tag,'36')
        a.Material2ME.v23 = real(a.Material2ME.v23)+1i*str2double(source.String);    
    end
    
    a.Material2ME.v21 = a.Material2ME.E2*a.Material2ME.v12/a.Material2ME.E1;
    if  a.MaterialTypeME == 1
        a.Material2ME.v31 = a.Material2ME.E3*a.Material2ME.v13/a.Material2ME.E1;
        a.Material2ME.v32 = a.Material2ME.E3*a.Material2ME.v23/a.Material2ME.E2;
        D = (1-a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v13*a.Material2ME.v31-a.Material2ME.v23*a.Material2ME.v32-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v31);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v23*a.Material2ME.v32)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v13*a.Material2ME.v32+a.Material2ME.v12)*a.Material2ME.E2/D; 
        a.Material2ME.C(1,3) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v13)*a.Material2ME.E3/D;
        a.Material2ME.C(2,2) = (1-a.Material2ME.v13*a.Material2ME.v31)*a.Material2ME.E2/D;
        a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v13+a.Material2ME.v23)*a.Material2ME.E3/D;
        a.Material2ME.C(3,3) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E3/D;
        
        a.C13UI5.String = real(a.Material2ME.C(1,3))/1e9;
        a.C33UI5.String = real(a.Material2ME.C(3,3))/1e9;
        
        a.C13ImagUI5.String = imag(a.Material2ME.C(1,3))/1e9;
        a.C33ImagUI5.String = imag(a.Material2ME.C(3,3))/1e9;
        
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = sqrt(real(a.Material2ME.C(3,3))/a.Material2ME.Density);
    elseif a.MaterialTypeME == 2
        a.Material2ME.G23 = .5*a.Material2ME.E2/(1+a.Material2ME.v23);
        
        D = (1-2*a.Material2ME.v12*a.Material2ME.v21-a.Material2ME.v23^2-2*a.Material2ME.v12*a.Material2ME.v23*a.Material2ME.v21);
        a.Material2ME.C(1,1) = (1-a.Material2ME.v23^2)*a.Material2ME.E1/D;
        a.Material2ME.C(1,2) = (a.Material2ME.v12*a.Material2ME.v23+a.Material2ME.v12)*a.Material2ME.E2/D;
        a.Material2ME.C(2,2) = (1-a.Material2ME.v12*a.Material2ME.v21)*a.Material2ME.E2/D;
        a.Material2ME.C(2,3) = (a.Material2ME.v21*a.Material2ME.v12+a.Material2ME.v23)*a.Material2ME.E2/D;
        
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_2;
        a.Material2ME.SlowShearVelocity_2 = sqrt(.5*(real(a.Material2ME.C(2,2))-real(a.Material2ME.C(2,3)))/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_3 = a.Material2ME.SlowShearVelocity_2;
        a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
        a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;
    end
    a.C11UI5.String = real(a.Material2ME.C(1,1))/1e9;
    a.C12UI5.String = real(a.Material2ME.C(1,2))/1e9;
    a.C22UI5.String = real(a.Material2ME.C(2,2))/1e9;
    a.C23UI5.String = real(a.Material2ME.C(2,3))/1e9;

    a.C11ImagUI5.String = imag(a.Material2ME.C(1,1))/1e9;
    a.C12ImagUI5.String = imag(a.Material2ME.C(1,2))/1e9;
    a.C22ImagUI5.String = imag(a.Material2ME.C(2,2))/1e9;
    a.C23ImagUI5.String = imag(a.Material2ME.C(2,3))/1e9;
    
    a.Material2ME.LongitudinalVelocity_1 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
elseif strcmp(source.Tag,'37') || strcmp(source.Tag,'46') % C11
    if  strcmp(source.Tag,'37')
        a.Material2ME.C(1,1) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(1,1));
    elseif strcmp(source.Tag,'46')
        a.Material2ME.C(1,1) = real(a.Material2ME.C(1,1))+1i*str2double(source.String)*1e9;
    end
    
    if  a.MaterialTypeME == 1
        D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,3)^2+a.Material2ME.C(3,3)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)*a.Material2ME.C(3,3)-2*a.Material2ME.C(1,2)*a.Material2ME.C(1,3)*a.Material2ME.C(2,3);
        a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.E2 = D/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
        a.Material2ME.E3 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        a.Material2ME.v23 = (a.Material2ME.C(1,2)*a.Material2ME.C(1,3)-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
    
        a.E2UI5.String = real(a.Material2ME.E2)/1e9;
        a.E3UI5.String = real(a.Material2ME.E3)/1e9;
        a.v23UI5.String = real(a.Material2ME.v23);

        a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
        a.E3ImagUI5.String = imag(a.Material2ME.E3)/1e9;
        a.v23ImagUI5.String = imag(a.Material2ME.v23);
    elseif a.MaterialTypeME == 2
        D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,2)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)^2-2*a.Material2ME.C(1,2)^2*a.Material2ME.C(2,3);
        a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)^2);
        a.Material2ME.E2 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        a.Material2ME.E3 = a.Material2ME.E2;
        a.Material2ME.v23 = (a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        
        a.E2UI5.String = real(a.Material2ME.E2)/1e9;
        a.v23UI5.String = real(a.Material2ME.v23);

        a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
        a.v23ImagUI5.String = imag(a.Material2ME.v23);
    elseif a.MaterialTypeME == 3
        D = 3*a.Material2ME.C(1,1)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)^3-2*a.Material2ME.C(1,2)^3;
        a.Material2ME.E1 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)^2);
        a.Material2ME.v12 = (a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(1,2))/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)^2);

        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
        a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;
        a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_2;
        
        a.v12UI5.String = real(a.Material2ME.v12);
        a.v12ImagUI5.String = imag(a.Material2ME.v12);
    end
    a.E1UI5.String = real(a.Material2ME.E1)/1e9;
    a.E1ImagUI5.String = imag(a.Material2ME.E1)/1e9;
    
    a.Material2ME.LongitudinalVelocity_1 = sqrt(real(a.Material2ME.C(1,1))/a.Material2ME.Density);
    a.LongitudinalVelocity_1UI5.String = a.Material2ME.LongitudinalVelocity_1;
elseif strcmp(source.Tag,'38') || strcmp(source.Tag,'47') % C12
    if  strcmp(source.Tag,'38')
        a.Material2ME.C(1,2) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(1,2));
    elseif strcmp(source.Tag,'47')
        a.Material2ME.C(1,2) = real(a.Material2ME.C(1,2))+1i*str2double(source.String)*1e9;
    end
    
    if  a.MaterialTypeME == 1
        D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,3)^2+a.Material2ME.C(3,3)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)*a.Material2ME.C(3,3)-2*a.Material2ME.C(1,2)*a.Material2ME.C(1,3)*a.Material2ME.C(2,3);
        a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.E2 = D/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
        a.Material2ME.E3 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        a.Material2ME.v12 = (a.Material2ME.C(1,3)*a.Material2ME.C(2,3)-a.Material2ME.C(1,2)*a.Material2ME.C(3,3))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.v13 = (a.Material2ME.C(1,2)*a.Material2ME.C(2,3)-a.Material2ME.C(1,3)*a.Material2ME.C(2,2))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.v23 = (a.Material2ME.C(1,2)*a.Material2ME.C(1,3)-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
    
        a.E2UI5.String = real(a.Material2ME.E2)/1e9;
        a.E3UI5.String = real(a.Material2ME.E3)/1e9;
        a.v13UI5.String = real(a.Material2ME.v13);
        a.v23UI5.String = real(a.Material2ME.v23);

        a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
        a.E3ImagUI5.String = imag(a.Material2ME.E3)/1e9;
        a.v13ImagUI5.String = imag(a.Material2ME.v13);
        a.v23ImagUI5.String = imag(a.Material2ME.v23);
    elseif a.MaterialTypeME == 2
        D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,2)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)^2-2*a.Material2ME.C(1,2)^2*a.Material2ME.C(2,3);    
        a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)^2);
        a.Material2ME.E2 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        a.Material2ME.E3 = a.Material2ME.E2;
        a.Material2ME.v12 = (a.Material2ME.C(1,2)*a.Material2ME.C(2,3)-a.Material2ME.C(1,2)*a.Material2ME.C(2,2))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)^2);
        a.Material2ME.v13 = a.Material2ME.v12;
        a.Material2ME.v23 = (a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        
        a.E2UI5.String = real(a.Material2ME.E2)/1e9;
        a.v23UI5.String = real(a.Material2ME.v23);

        a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
        a.v23ImagUI5.String = imag(a.Material2ME.v23);
    elseif a.MaterialTypeME == 3
        D = 3*a.Material2ME.C(1,1)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)^3-2*a.Material2ME.C(1,2)^3;    
        a.Material2ME.E1 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)^2);
        a.Material2ME.v12 = (a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(1,2))/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)^2);
    end
    a.E1UI5.String = real(a.Material2ME.E1)/1e9;
    a.v12UI5.String = real(a.Material2ME.v12);

    a.E1ImagUI5.String = imag(a.Material2ME.E1)/1e9;
    a.v12ImagUI5.String = imag(a.Material2ME.v12);
elseif strcmp(source.Tag,'39') || strcmp(source.Tag,'48') % C13
    if  strcmp(source.Tag,'39')
        a.Material2ME.C(1,3) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(1,3));
    elseif strcmp(source.Tag,'48')
        a.Material2ME.C(1,3) = real(a.Material2ME.C(1,3))+1i*str2double(source.String)*1e9;
    end

    D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,3)^2+a.Material2ME.C(3,3)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)*a.Material2ME.C(3,3)-2*a.Material2ME.C(1,2)*a.Material2ME.C(1,3)*a.Material2ME.C(2,3);
    a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
    a.Material2ME.E2 = D/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
    a.Material2ME.E3 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
    a.Material2ME.v12 = (a.Material2ME.C(1,3)*a.Material2ME.C(2,3)-a.Material2ME.C(1,2)*a.Material2ME.C(3,3))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
    a.Material2ME.v13 = (a.Material2ME.C(1,2)*a.Material2ME.C(2,3)-a.Material2ME.C(1,3)*a.Material2ME.C(2,2))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
    a.Material2ME.v23 = (a.Material2ME.C(1,2)*a.Material2ME.C(1,3)-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));

    a.E1UI5.String = real(a.Material2ME.E1)/1e9;
    a.E2UI5.String = real(a.Material2ME.E2)/1e9;
    a.E3UI5.String = real(a.Material2ME.E3)/1e9;
    a.v12UI5.String = real(a.Material2ME.v12);
    a.v13UI5.String = real(a.Material2ME.v13);
    a.v23UI5.String = real(a.Material2ME.v23);

    a.E1ImagUI5.String = imag(a.Material2ME.E1)/1e9;
    a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
    a.E3ImagUI5.String = imag(a.Material2ME.E3)/1e9;
    a.v12ImagUI5.String = imag(a.Material2ME.v12);
    a.v13ImagUI5.String = imag(a.Material2ME.v13);
    a.v23ImagUI5.String = imag(a.Material2ME.v23);
elseif strcmp(source.Tag,'40') || strcmp(source.Tag,'49') % C22
    if  strcmp(source.Tag,'40')
        a.Material2ME.C(2,2) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(2,2));
    elseif strcmp(source.Tag,'49')
        a.Material2ME.C(2,2) = real(a.Material2ME.C(2,2))+1i*str2double(source.String)*1e9;
    end
    
    if  a.MaterialTypeME == 1
        D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,3)^2+a.Material2ME.C(3,3)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)*a.Material2ME.C(3,3)-2*a.Material2ME.C(1,2)*a.Material2ME.C(1,3)*a.Material2ME.C(2,3);
        a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.E2 = D/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
        a.Material2ME.E3 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        a.Material2ME.v12 = (a.Material2ME.C(1,3)*a.Material2ME.C(2,3)-a.Material2ME.C(1,2)*a.Material2ME.C(3,3))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.v13 = (a.Material2ME.C(1,2)*a.Material2ME.C(2,3)-a.Material2ME.C(1,3)*a.Material2ME.C(2,2))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
    
        a.E3UI5.String = real(a.Material2ME.E3)/1e9;
        a.v13UI5.String = real(a.Material2ME.v13);

        a.E3ImagUI5.String = imag(a.Material2ME.E3)/1e9;
        a.v13ImagUI5.String = imag(a.Material2ME.v13);
        
        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
    elseif a.MaterialTypeME == 2
        a.Material2ME.G23 = .5*(a.Material2ME.C(2,2)-a.Material2ME.C(2,3));
        
        D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,2)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)^2-2*a.Material2ME.C(1,2)^2*a.Material2ME.C(2,3);    
        a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)^2);
        a.Material2ME.E2 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        a.Material2ME.E3 = a.Material2ME.E2;
        a.Material2ME.v12 = (a.Material2ME.C(1,2)*a.Material2ME.C(2,3)-a.Material2ME.C(1,2)*a.Material2ME.C(2,2))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)^2);
        a.Material2ME.v13 = a.Material2ME.v12;
        a.Material2ME.v23 = (a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));

        a.Material2ME.LongitudinalVelocity_2 = sqrt(real(a.Material2ME.C(2,2))/a.Material2ME.Density);
        a.Material2ME.LongitudinalVelocity_3 = a.Material2ME.LongitudinalVelocity_2;
        a.Material2ME.SlowShearVelocity_2 = sqrt(.5*(real(a.Material2ME.C(2,2))-real(a.Material2ME.C(2,3)))/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_3 = a.Material2ME.SlowShearVelocity_2;
        a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
        a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
        a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;
    end
    a.E1UI5.String = real(a.Material2ME.E1)/1e9;
    a.E2UI5.String = real(a.Material2ME.E2)/1e9;
    a.v12UI5.String = real(a.Material2ME.v12);
    a.v23UI5.String = real(a.Material2ME.v23);

    a.E1ImagUI5.String = imag(a.Material2ME.E1)/1e9;
    a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
    a.v12ImagUI5.String = imag(a.Material2ME.v12);
    a.v23ImagUI5.String = imag(a.Material2ME.v23);

    a.LongitudinalVelocity_2UI5.String = a.Material2ME.LongitudinalVelocity_2;   
elseif strcmp(source.Tag,'41') || strcmp(source.Tag,'50') % C23
    if  strcmp(source.Tag,'41')
        a.Material2ME.C(2,3) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(2,3));
    elseif strcmp(source.Tag,'50')
        a.Material2ME.C(2,3) = real(a.Material2ME.C(2,3))+1i*str2double(source.String)*1e9;
    end
    
    if  a.MaterialTypeME == 1
        D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,3)^2+a.Material2ME.C(3,3)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)*a.Material2ME.C(3,3)-2*a.Material2ME.C(1,2)*a.Material2ME.C(1,3)*a.Material2ME.C(2,3);
        a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.E2 = D/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
        a.Material2ME.E3 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        a.Material2ME.v12 = (a.Material2ME.C(1,3)*a.Material2ME.C(2,3)-a.Material2ME.C(1,2)*a.Material2ME.C(3,3))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.v13 = (a.Material2ME.C(1,2)*a.Material2ME.C(2,3)-a.Material2ME.C(1,3)*a.Material2ME.C(2,2))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
        a.Material2ME.v23 = (a.Material2ME.C(1,2)*a.Material2ME.C(1,3)-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
    
        a.E3UI5.String = real(a.Material2ME.E3)/1e9;
        a.v13UI5.String = real(a.Material2ME.v13);

        a.E3ImagUI5.String = imag(a.Material2ME.E3)/1e9;
        a.v13ImagUI5.String = imag(a.Material2ME.v13);
    elseif a.MaterialTypeME == 2
        a.Material2ME.G23 = .5*(a.Material2ME.C(2,2)-a.Material2ME.C(2,3));
        
        D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,2)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)^2-2*a.Material2ME.C(1,2)^2*a.Material2ME.C(2,3);    
        a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)^2);
        a.Material2ME.E2 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        a.Material2ME.E3 = a.Material2ME.E2;
        a.Material2ME.v12 = (a.Material2ME.C(1,2)*a.Material2ME.C(2,3)-a.Material2ME.C(1,2)*a.Material2ME.C(2,2))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)^2);
        a.Material2ME.v13 = a.Material2ME.v12;
        a.Material2ME.v23 = (a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
        
        a.Material2ME.SlowShearVelocity_2 = sqrt(.5*(real(a.Material2ME.C(2,2))-real(a.Material2ME.C(2,3)))/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_3 = a.Material2ME.SlowShearVelocity_2;
        a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
        a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;
    end
    a.E1UI5.String = real(a.Material2ME.E1)/1e9;
    a.E2UI5.String = real(a.Material2ME.E2)/1e9;
    a.v12UI5.String = real(a.Material2ME.v12);
    a.v23UI5.String = real(a.Material2ME.v23);    

    a.E1ImagUI5.String = imag(a.Material2ME.E1)/1e9;
    a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
    a.v12ImagUI5.String = imag(a.Material2ME.v12);
    a.v23ImagUI5.String = imag(a.Material2ME.v23);
elseif strcmp(source.Tag,'42') || strcmp(source.Tag,'51') % C33
    if  strcmp(source.Tag,'42')
        a.Material2ME.C(3,3) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(3,3));
    elseif strcmp(source.Tag,'51')
        a.Material2ME.C(3,3) = real(a.Material2ME.C(3,3))+1i*str2double(source.String)*1e9;
    end
    
    D = a.Material2ME.C(1,1)*a.Material2ME.C(2,3)^2+a.Material2ME.C(2,2)*a.Material2ME.C(1,3)^2+a.Material2ME.C(3,3)*a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2)*a.Material2ME.C(3,3)-2*a.Material2ME.C(1,2)*a.Material2ME.C(1,3)*a.Material2ME.C(2,3);
    a.Material2ME.E1 = D/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
    a.Material2ME.E2 = D/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));
    a.Material2ME.E3 = D/(a.Material2ME.C(1,2)^2-a.Material2ME.C(1,1)*a.Material2ME.C(2,2));
    a.Material2ME.v12 = (a.Material2ME.C(1,3)*a.Material2ME.C(2,3)-a.Material2ME.C(1,2)*a.Material2ME.C(3,3))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
    a.Material2ME.v13 = (a.Material2ME.C(1,2)*a.Material2ME.C(2,3)-a.Material2ME.C(1,3)*a.Material2ME.C(2,2))/(a.Material2ME.C(2,3)^2-a.Material2ME.C(2,2)*a.Material2ME.C(3,3));
    a.Material2ME.v23 = (a.Material2ME.C(1,2)*a.Material2ME.C(1,3)-a.Material2ME.C(1,1)*a.Material2ME.C(2,3))/(a.Material2ME.C(1,3)^2-a.Material2ME.C(1,1)*a.Material2ME.C(3,3));

    a.E1UI5.String = real(a.Material2ME.E1)/1e9;
    a.E2UI5.String = real(a.Material2ME.E2)/1e9;
    a.E3UI5.String = real(a.Material2ME.E3)/1e9;
    a.v12UI5.String = real(a.Material2ME.v12);
    a.v13UI5.String = real(a.Material2ME.v13);
    a.v23UI5.String = real(a.Material2ME.v23);

    a.E1ImagUI5.String = imag(a.Material2ME.E1)/1e9;
    a.E2ImagUI5.String = imag(a.Material2ME.E2)/1e9;
    a.E3ImagUI5.String = imag(a.Material2ME.E3)/1e9;
    a.v12ImagUI5.String = imag(a.Material2ME.v12);
    a.v13ImagUI5.String = imag(a.Material2ME.v13);
    a.v23ImagUI5.String = imag(a.Material2ME.v23);
    
    a.Material2ME.LongitudinalVelocity_3 = sqrt(real(a.Material2ME.C(3,3))/a.Material2ME.Density);
    a.LongitudinalVelocity_3UI5.String = a.Material2ME.LongitudinalVelocity_3;
elseif strcmp(source.Tag,'43') || strcmp(source.Tag,'52') % C44
    if  strcmp(source.Tag,'43')
        a.Material2ME.C(4,4) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(4,4));
    elseif strcmp(source.Tag,'52')
        a.Material2ME.C(4,4) = real(a.Material2ME.C(4,4))+1i*str2double(source.String)*1e9;
    end
    
    a.Material2ME.G23 = a.Material2ME.C(4,4);
    a.G23UI5.String = real(a.Material2ME.G23)/1e9;
    a.G23ImagUI5.String = imag(a.Material2ME.G23)/1e9;
    
    a.Material2ME.SlowShearVelocity_2 = sqrt(real(a.Material2ME.C(4,4))/a.Material2ME.Density);
    a.Material2ME.SlowShearVelocity_3 = a.Material2ME.SlowShearVelocity_2;
    a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_2;
    a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_3;
elseif strcmp(source.Tag,'44') || strcmp(source.Tag,'53') % C55
    if  strcmp(source.Tag,'44')
        a.Material2ME.C(5,5) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(5,5));
    elseif strcmp(source.Tag,'53')
        a.Material2ME.C(5,5) = real(a.Material2ME.C(5,5))+1i*str2double(source.String)*1e9;
    end
    
    a.Material2ME.G13 = a.Material2ME.C(5,5);
    a.G13UI5.String = real(a.Material2ME.G13)/1e9;
    a.G13ImagUI5.String = imag(a.Material2ME.G13)/1e9;
    
    a.Material2ME.FastShearVelocity_1 = sqrt(real(a.Material2ME.C(5,5))/a.Material2ME.Density);
    a.Material2ME.FastShearVelocity_3 = a.Material2ME.FastShearVelocity_1;
    a.FastShearVelocity_1UI5.String = a.Material2ME.FastShearVelocity_1;
    a.FastShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_3;
elseif strcmp(source.Tag,'45') || strcmp(source.Tag,'54') % C66
    if  strcmp(source.Tag,'45')
        a.Material2ME.C(6,6) = str2double(source.String)*1e9+1i*imag(a.Material2ME.C(6,6));
    elseif strcmp(source.Tag,'54')
        a.Material2ME.C(6,6) = real(a.Material2ME.C(6,6))+1i*str2double(source.String)*1e9;
    end
    
    a.Material2ME.G12 = a.Material2ME.C(6,6);
    a.G12UI5.String = real(a.Material2ME.G12)/1e9;
    a.G12ImagUI5.String = imag(a.Material2ME.G12)/1e9;
    if  a.MaterialTypeME == 1
        a.Material2ME.SlowShearVelocity_1 = sqrt(real(a.Material2ME.C(6,6))/a.Material2ME.Density);
        a.Material2ME.FastShearVelocity_2 = a.Material2ME.SlowShearVelocity_1;
    elseif a.MaterialTypeME == 2
        a.Material2ME.G13 = a.Material2ME.C(6,6);
        a.Material2ME.C(5,5) = a.Material2ME.C(6,6);
        
        a.Material2ME.FastShearVelocity_1 = sqrt(real(a.Material2ME.G12)/a.Material2ME.Density);
        a.Material2ME.SlowShearVelocity_1 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.FastShearVelocity_2 = a.Material2ME.FastShearVelocity_1;
        a.Material2ME.FastShearVelocity_3 = a.Material2ME.FastShearVelocity_1;
        a.FastShearVelocity_1UI5.String = a.Material2ME.FastShearVelocity_1;
        a.FastShearVelocity_3UI5.String = a.Material2ME.FastShearVelocity_3;
    elseif a.MaterialTypeME == 3
        a.Material2ME.SlowShearVelocity_1 = sqrt(real(a.Material2ME.G12)/a.Material2ME.Density);
        a.Material2ME.FastShearVelocity_2 = a.Material2ME.SlowShearVelocity_1;
        a.FastShearVelocity_1UI5.String = a.Material2ME.SlowShearVelocity_1;
        a.SlowShearVelocity_2UI5.String = a.Material2ME.SlowShearVelocity_1;
        a.FastShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_1;
        a.SlowShearVelocity_3UI5.String = a.Material2ME.SlowShearVelocity_1;
    end
    a.SlowShearVelocity_1UI5.String = a.Material2ME.SlowShearVelocity_1;
    a.FastShearVelocity_2UI5.String = a.Material2ME.FastShearVelocity_2;
elseif strcmp(source.Tag,'55') % Material name (anisotropic)
    a.MaterialName4B = source.String;
elseif strcmp(source.Tag,'56') % Save (anisotropic)
    if  any(strcmp(a.MaterialName4B,fieldnames(a.Materials.Cubic))) || any(strcmp(a.MaterialName4B,fieldnames(a.Materials.TransverselyIsotropic))) || any(strcmp(a.MaterialName4B,fieldnames(a.Materials.Orthotropic)))
        Selection = questdlg(['Overwrite ',a.MaterialName4B,'?'],'Confirm','Ok','Cancel','Ok'); 
        switch Selection
        case 'Ok'
            Save2
        end
    else
        Save2
    end
elseif strcmp(source.Tag,'57') % Delete (anisotropic)
    Selection = questdlg(['Delete ',a.MaterialName4B,'?'],'Confirm','Ok','Cancel','Ok');
    a.Material2UI5.Value = 1;
    switch Selection
    case 'Ok'
        if  a.MaterialTypeME == 1
            a.Materials.Orthotropic = rmfield(a.Materials.Orthotropic,a.MaterialName4B);
            a.Material2UI5.String = fieldnames(a.Materials.Orthotropic);
            if  a.MaterialType2 == 1
                a.MaterialUI2Value = 1;
            end
            if  a.MaterialType_Polar == 1
                a.MaterialUI4Value = 1;
            end
            if  a.MaterialType6 == 1
                a.MaterialUI6Value = 1;
            end
            if  a.MaterialType_Bulk == 1
                a.MaterialUI8.String = fieldnames(a.Materials.Orthotropic);
                a.MaterialUI8.Value = 1;
            end
            if  a.SolidType_Bulk == 4
                a.SolidUI8.String = fieldnames(a.Materials.Orthotropic);
                a.SolidUI8.Value = 1;
            end
            E = fieldnames(a.Materials.Orthotropic);
            for i = 1:length(fieldnames(a.Materials.Orthotropic))
                F = struct2cell(getfield(a.Materials.Orthotropic,E{i}))';
                G(i,1:11) = [F(1) F(3:12)];
            end
            writetable(cell2table(G),fullfile(a.MaterialListDirectory,'MaterialList_Orthotropic.txt'),'delimiter',' ','WriteVariableNames',false)
        elseif a.MaterialTypeME == 2 
            a.Materials.TransverselyIsotropic = rmfield(a.Materials.TransverselyIsotropic,a.MaterialName4B);
            a.Material2UI5.String = fieldnames(a.Materials.TransverselyIsotropic);
            if  a.MaterialType2 == 2
                a.MaterialUI2Value = 1;
            end
            if  a.MaterialType_Polar == 2
                a.MaterialUI4Value = 1;
            end
            if  a.MaterialType6 == 2
                a.MaterialUI6Value = 1;
            end
            if  a.MaterialType_Bulk == 2
                a.MaterialUI8.String = fieldnames(a.Materials.TransverselyIsotropic);
                a.MaterialUI8.Value = 1;
            end   
            if  a.SolidType_Bulk == 3
                a.SolidUI8.String = fieldnames(a.Materials.TransverselyIsotropic);
                a.SolidUI8.Value = 1;
            end
            E = fieldnames(a.Materials.TransverselyIsotropic);
            for i = 1:length(fieldnames(a.Materials.TransverselyIsotropic))
                F = struct2cell(getfield(a.Materials.TransverselyIsotropic,E{i}))';
                G(i,1:7) = [F(1) F(3:5) F(7) F(10) F(12)];
            end
            writetable(cell2table(G),fullfile(a.MaterialListDirectory,'MaterialList_TransverselyIsotropic.txt'),'delimiter',' ','WriteVariableNames',false)
        elseif a.MaterialTypeME == 3 
            a.Materials.Cubic = rmfield(a.Materials.Cubic,a.MaterialName4B);
            a.Material2UI5.String = fieldnames(a.Materials.Cubic);
            if  a.MaterialType2 == 3
                a.MaterialUI2Value = 1;
            end
            if  a.MaterialType_Polar == 3
                a.MaterialUI4Value = 1;
            end
            if  a.MaterialType6 == 3
                a.MaterialUI6Value = 1;
            end
            if  a.MaterialType_Bulk == 3
                a.MaterialUI8.String = fieldnames(a.Materials.Cubic);
                a.MaterialUI8.Value = 1;
            end
            if  a.SolidType_Bulk == 2
                a.SolidUI8.String = fieldnames(a.Materials.Cubic);
                a.SolidUI8.Value = 1;
            end
            E = fieldnames(a.Materials.Cubic);
            for i = 1:length(fieldnames(a.Materials.Cubic))
                F = struct2cell(getfield(a.Materials.Cubic,E{i}))';
                G(i,1:5) = [F(1) F(3:4) F(7) F(10)];
            end
            writetable(cell2table(G),fullfile(a.MaterialListDirectory,'MaterialList_Cubic.txt'),'delimiter',' ','WriteVariableNames',false)
        end
    end
elseif strcmp(source.Tag,'58') % Fluid
    a.Material3ME = getfield(a.Materials.Fluid,cell2mat(source.String(source.Value)));
    a.Density3UI5.String = a.Material3ME.Density;
    a.VelocityUI5.String = a.Material3ME.Velocity;
    a.MaterialName4C = a.Material3ME.Name;
    a.Name3UI5.String = a.MaterialName4C;
elseif strcmp(source.Tag,'59') % Density (fluid)
    a.Material3ME.Density = str2double(source.String);
elseif strcmp(source.Tag,'60') % Velocity
    a.Material3ME.Velocity = str2double(source.String);
elseif strcmp(source.Tag,'61') % Medium name
    a.MaterialName4C = source.String;
elseif strcmp(source.Tag,'62') % Save (fluid)
    if  any(strcmp(a.MaterialName4C,fieldnames(a.Materials.Fluid)))
        Selection = questdlg(['Overwrite ',a.MaterialName4C,'?'],'Confirm','Ok','Cancel','Ok'); 
        switch Selection 
        case 'Ok'
            Save3
        end
    else
        Save3
    end
elseif strcmp(source.Tag,'63') % Delete (fluid)
    Selection = questdlg(['Delete ',a.MaterialName4C,'?'],'Confirm','Ok','Cancel','Ok'); 
    switch Selection 
    case 'Ok'
        a.Materials.Fluid = rmfield(a.Materials.Fluid,a.MaterialName4C);
        a.FluidUI1.String = fieldnames(a.Materials.Fluid);
        if  a.Quantity11 == 4
            a.Option1UI1.String = fieldnames(a.Materials.Fluid);
        end
        a.Material3UI5.String = fieldnames(a.Materials.Fluid);
        a.CouplantUI8.String = fieldnames(a.Materials.Fluid);
        a.FluidUI1.Value = 1;
        a.Option1UI1.Value = 1;
        a.Material3UI5.Value = 1;
        a.CouplantUI8.Value = 1;
        E = fieldnames(a.Materials.Fluid);
        for i = 1:length(fieldnames(a.Materials.Fluid))
            F = struct2cell(getfield(a.Materials.Fluid,E{i}))';
            G(i,1:3) = [F(1) F(3:4)];
        end
        writetable(cell2table(G(:,1:3)),fullfile(a.MaterialListDirectory,'MaterialList_Fluid.txt'),'delimiter',' ','WriteVariableNames',false)
    end
elseif strcmp(source.Tag,'68') % Explanation screen (isotropic)
    f6 = figure('NumberTitle','off','Name','Help: Isotropic materials','Visible','off','MenuBar','none','Position',[0 0 610 410],'color','w');

    jframe = get(gcf,'javaframe'); %#ok<*JAVFM> 
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
  
    f6.Units = 'normalized';
    movegui(f6,'center')
    f6.Visible = 'on';

    uicontrol('Parent',f6,'Style','text','String','In the Isotropic materials panel, you can enter isotropic material parameters. After','Position',[10 370 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','entering the mass density of your material, choose one of three categories to enter','Position',[10 370-1*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','the remaining parameters. Either enter the engineering constants or the stiffness','Position',[10 370-2*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','components or the bulk wave characteristics. When editing the parameters in one','Position',[10 370-3*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','category, those in the other two are updated automatically.','Position',[10 370-4*18 590 23],'FontSize',11,'Backgroundcolor','white');

    uicontrol('Parent',f6,'Style','text','String','The engineering constants and stiffness components have real and imaginary parts,','Position',[10 370-6*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','where the real parts describe the elastic behavior in the usual way, and the','Position',[10 370-7*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','imaginary parts describe the damping. If all imaginary parts are zero, the material is','Position',[10 370-8*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','perfectly elastic, i.e., without material damping. Otherwise it is viscoelastic so that','Position',[10 370-9*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','elastic waves are damped. Notice that the real parts may never be zero.','Position',[10 370-10*18 590 23],'FontSize',11,'Backgroundcolor','white');

    uicontrol('Parent',f6,'Style','text','String','If you decide to enter the bulk wave characteristics, you need to enter the wave','Position',[10 370-12*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','speeds and, in case you want to consider material damping, you also need to enter','Position',[10 370-13*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','the corresponding damping coefficients. You can choose to enter them in Nepers or','Position',[10 370-14*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','Decibels per wavelength or per meter.','Position',[10 370-15*18 590 23],'FontSize',11,'Backgroundcolor','white');

    uicontrol('Parent',f6,'Style','text','String','Notice the Materials menu in the menu bar. It is useful for example when you want to','Position',[10 370-17*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','transfer your custom materials to a newer DC version. Use Materials => Export to ','Position',[10 370-18*18 590 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f6,'Style','text','String','export the material lists. After installing the new DC version, use Materials => Import.','Position',[10 370-18*19 590 23],'FontSize',11,'Backgroundcolor','white');
elseif strcmp(source.Tag,'69') % Explanation screen (anisotropic)
    f5 = figure('NumberTitle','off','Name','Help: Anisotropic materials','Visible','off','MenuBar','none','Position',[0 0 800 500],'color','w');

    jframe = get(gcf,'javaframe'); %#ok<*JAVFM> 
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
  
    f5.Units = 'normalized';
    movegui(f5,'center')
    f5.Visible = 'on';

    uicontrol('Parent',f5,'Style','text','String','In the Anisotropic materials panel, you can enter orthotropic, transversely isotropic, and cubic materials. Make','Position',[10 460 780 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','your symmetry class selection in the Class drop-down menu. Compared to orthotropic materials, the number of','Position',[10 460-1*18 780 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','independent entries reduces from 9 to 5 for transversely isotropic materials to 3 for cubic ones.','Position',[10 460-2*18 780 23],'FontSize',11,'Backgroundcolor','white');

    uicontrol('Parent',f5,'Style','text','String','After entering the mass density of your material, choose one of two ways to enter the remaining parameters.','Position',[10 460-4*18 780 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','Either enter the engineering constants or the stiffness components. When editing the engineering constants, the','Position',[10 460-5*18 780 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','stiffness components are updated automatically and vice versa. The bulk wave velocities are updated as well,','Position',[10 460-6*18 780 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','but they cannot be edited directly.','Position',[10 460-7*18 780 23],'FontSize',11,'Backgroundcolor','white');

    uicontrol('Parent',f5,'Style','text','String','Be aware of the axis convention used by the DC.','Position',[395 460-9*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','Consider the unidirectional composite ply sketched','Position',[395 460-10*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','on the left. The primed coordinate system ','Position',[395 460-11*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','x''i = (x''1,x''2,x''3) represents the crystallographic axes','Position',[395 460-12*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','of the composite. The carbon fibers are aligned in','Position',[395 460-13*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','the x''1-direction. The x''2-direction is in-plane and','Position',[395 460-14*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','perpendicular to the fibers and x''3 is out-of-plane.','Position',[395 460-15*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','The engineering constants and stiffness components','Position',[395 460-16*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','are defined in x''i. For instance, C11 has the highest','Position',[395 460-17*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','value among all stiffness components, as it gives the','Position',[395 460-18*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','stiffness in the fiber direction. The global coordinate','Position',[395 460-19*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','system xi = (x1,x2,x3) is the reference frame of wave','Position',[395 460-20*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','propagation. It is rotated about the x''3-axis through','Position',[395 460-21*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','the propagation angle Phi, which you can adjust in','Position',[395 460-22*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','the Anisotropic tab. Wave propagation takes place in','Position',[395 460-23*18 390 23],'FontSize',11,'Backgroundcolor','white');
    uicontrol('Parent',f5,'Style','text','String','the x1-direction.','Position',[395 460-24*18 390 23],'FontSize',11,'Backgroundcolor','white');

    axes('Parent',f5,'Units','pixels','Position',[20 20 371 300]);
    imshow('CoordinateSystems.png')
end
function Save1
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'Name',a.MaterialName4A);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'Class','Isotropic');
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'Density',a.Material1ME.Density);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'YoungsModulus',a.Material1ME.YoungsModulus);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'YoungsModulus_complex',a.Material1ME.YoungsModulus_complex);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'PoissonsNumber',a.Material1ME.PoissonsNumber);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'PoissonsNumber_complex',a.Material1ME.PoissonsNumber_complex);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'LongitudinalAttenuation',a.Material1ME.LongitudinalAttenuation);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'TransverseAttenuation',a.Material1ME.TransverseAttenuation);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'Lambda',a.Material1ME.Lambda);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'Lambda_complex',a.Material1ME.Lambda_complex);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'Mu',a.Material1ME.Mu);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'Mu_complex',a.Material1ME.Mu_complex);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'LongitudinalVelocity',a.Material1ME.LongitudinalVelocity);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'LongitudinalVelocity_complex',a.Material1ME.LongitudinalVelocity_complex);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'TransverseVelocity',a.Material1ME.TransverseVelocity);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'TransverseVelocity_complex',a.Material1ME.TransverseVelocity_complex);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'PlateVelocity',a.Material1ME.PlateVelocity);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'CylinderVelocity',a.Material1ME.CylinderVelocity);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'RayleighVelocity',a.Material1ME.RayleighVelocity);
    a.Materials.Isotropic = setfield(a.Materials.Isotropic,{1,1},a.MaterialName4A,{1,1},'C',[a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex a.Material1ME.Lambda_complex a.Material1ME.Lambda_complex 0 0 0;0 a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex a.Material1ME.Lambda_complex 0 0 0;0 0 a.Material1ME.Lambda_complex+2*a.Material1ME.Mu_complex 0 0 0;0 0 0 a.Material1ME.Mu_complex 0 0;0 0 0 0 a.Material1ME.Mu_complex 0;0 0 0 0 0 a.Material1ME.Mu_complex]);
    a.Materials.Isotropic = orderfields(a.Materials.Isotropic);
    a.Material1UI5.String = fieldnames(a.Materials.Isotropic);
    a.MaterialUI1.String = fieldnames(a.Materials.Isotropic);
    if  a.MaterialType_Bulk == 4
        a.MaterialUI8.String = fieldnames(a.Materials.Isotropic);
    end
    if  a.SolidType_Bulk == 1
        a.SolidUI8.String = fieldnames(a.Materials.Isotropic);
    end
    E = fieldnames(a.Materials.Isotropic);
    for i = 1:length(fieldnames(a.Materials.Isotropic))
        F = struct2cell(getfield(a.Materials.Isotropic,E{i}))';
        G(i,1:6) = [F(1) F(3:4) F(6) F(8:9)];
    end
    writetable(cell2table(G(:,1:6)),fullfile(a.MaterialListDirectory,'MaterialList_Isotropic.txt'),'delimiter',' ','WriteVariableNames',false)
end
function Save2
    if  a.MaterialTypeME == 1
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'Name',a.MaterialName4B);
        if  a.Material2ME.C(3,3) == a.Material2ME.C(1,1) && a.Material2ME.C(2,2) == a.Material2ME.C(1,1) && a.Material2ME.C(1,3) == a.Material2ME.C(1,2) && a.Material2ME.C(2,3) == a.Material2ME.C(1,2) && a.Material2ME.C(4,4) == a.Material2ME.C(6,6) && a.Material2ME.C(5,5) == a.Material2ME.C(6,6)
            if  abs(a.Material2ME.C(1,1)-a.Material2ME.C(1,2)-2*a.Material2ME.C(4,4))/a.Material2ME.C(1,1) < 5e-3
                a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'Class','Isotropic');
            else
                a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'Class','Cubic');
            end
        elseif a.Material2ME.C(1,3) == a.Material2ME.C(1,2) && a.Material2ME.C(3,3) == a.Material2ME.C(2,2) && a.Material2ME.C(5,5) == a.Material2ME.C(6,6)
            a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'Class','Transversely isotropic');
        else
            a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'Class','Orthotropic');
        end
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'Density',a.Material2ME.Density);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'E1',a.Material2ME.E1);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'E2',a.Material2ME.E2);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'E3',a.Material2ME.E3);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'G12',a.Material2ME.G12);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'G13',a.Material2ME.G13);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'G23',a.Material2ME.G23);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'v12',a.Material2ME.v12);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'v13',a.Material2ME.v13);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'v23',a.Material2ME.v23);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_1',a.Material2ME.LongitudinalVelocity_1);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_2',a.Material2ME.LongitudinalVelocity_2);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_3',a.Material2ME.LongitudinalVelocity_3);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_1',a.Material2ME.FastShearVelocity_1);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_2',a.Material2ME.FastShearVelocity_2);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_3',a.Material2ME.FastShearVelocity_3);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_1',a.Material2ME.SlowShearVelocity_1);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_2',a.Material2ME.SlowShearVelocity_2);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_3',a.Material2ME.SlowShearVelocity_3);
        a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'C',[a.Material2ME.C(1,1) a.Material2ME.C(1,2) a.Material2ME.C(1,3) 0 0 0;0 a.Material2ME.C(2,2) a.Material2ME.C(2,3) 0 0 0;0 0 a.Material2ME.C(3,3) 0 0 0;0 0 0 a.Material2ME.G23 0 0;0 0 0 0 a.Material2ME.G13 0;0 0 0 0 0 a.Material2ME.G12]);
        if  a.Material2ME.C(3,3) == a.Material2ME.C(1,1) && a.Material2ME.C(2,2) == a.Material2ME.C(1,1) && a.Material2ME.C(1,3) == a.Material2ME.C(1,2) && a.Material2ME.C(2,3) == a.Material2ME.C(1,2) && a.Material2ME.C(4,4) == a.Material2ME.C(6,6) && a.Material2ME.C(5,5) == a.Material2ME.C(6,6) && abs(a.Material2ME.C(1,1)-a.Material2ME.C(1,2)-2*a.Material2ME.C(4,4))/a.Material2ME.C(1,1) < 5e-3
            a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'PlateVelocity',sqrt(a.Material2ME.E1/(a.Material2ME.Density*(1-a.Material2ME.v12^2))));
            a.Materials.Orthotropic = setfield(a.Materials.Orthotropic,{1,1},a.MaterialName4B,{1,1},'CylinderVelocity',sqrt(a.Material2ME.E1/a.Material2ME.Density));
        end
        a.Materials.Orthotropic = orderfields(a.Materials.Orthotropic);
        a.Material2UI5.String = fieldnames(a.Materials.Orthotropic);
        if  a.MaterialType_Bulk == 1
            a.MaterialUI8.String = fieldnames(a.Materials.Orthotropic);
        end
        if  a.SolidType_Bulk == 4
            a.SolidUI8.String = fieldnames(a.Materials.Orthotropic);
        end
        E = fieldnames(a.Materials.Orthotropic);
        for i = 1:length(fieldnames(a.Materials.Orthotropic))
            F = struct2cell(getfield(a.Materials.Orthotropic,E{i}))';
            G(i,1:11) = [F(1) F(3:12)];
        end
        writetable(cell2table(G),fullfile(a.MaterialListDirectory,'MaterialList_Orthotropic.txt'),'delimiter',' ','WriteVariableNames',false)
    elseif a.MaterialTypeME == 2
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'Name',a.MaterialName4B);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'Class','Transversely isotropic');
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'Density',a.Material2ME.Density);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'E1',a.Material2ME.E1);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'E2',a.Material2ME.E2);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'E3',a.Material2ME.E3);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'G12',a.Material2ME.G12);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'G13',a.Material2ME.G13);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'G23',a.Material2ME.G23);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'v12',a.Material2ME.v12);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'v13',a.Material2ME.v13);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'v23',a.Material2ME.v23);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_1',a.Material2ME.LongitudinalVelocity_1);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_2',a.Material2ME.LongitudinalVelocity_2);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_3',a.Material2ME.LongitudinalVelocity_3);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_1',a.Material2ME.FastShearVelocity_1);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_2',a.Material2ME.FastShearVelocity_2);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_3',a.Material2ME.FastShearVelocity_3);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_1',a.Material2ME.SlowShearVelocity_1);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_2',a.Material2ME.SlowShearVelocity_2);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_3',a.Material2ME.SlowShearVelocity_3);
        a.Materials.TransverselyIsotropic = setfield(a.Materials.TransverselyIsotropic,{1,1},a.MaterialName4B,{1,1},'C',[a.Material2ME.C(1,1) a.Material2ME.C(1,2) a.Material2ME.C(1,2) 0 0 0;0 a.Material2ME.C(2,2) a.Material2ME.C(2,3) 0 0 0;0 0 a.Material2ME.C(2,2) 0 0 0;0 0 0 .5*(a.Material2ME.C(2,2)-a.Material2ME.C(2,3)) 0 0;0 0 0 0 a.Material2ME.G12 0;0 0 0 0 0 a.Material2ME.G12]);
        a.Materials.TransverselyIsotropic = orderfields(a.Materials.TransverselyIsotropic);
        a.Material2UI5.String = fieldnames(a.Materials.TransverselyIsotropic);
        a.MaterialUI8.String = fieldnames(a.Materials.TransverselyIsotropic);
        if  a.MaterialType_Bulk == 2
            a.MaterialUI8.String = fieldnames(a.Materials.TransverselyIsotropic);
        end   
        if  a.SolidType_Bulk == 3
            a.SolidUI8.String = fieldnames(a.Materials.TransverselyIsotropic);
        end
        E = fieldnames(a.Materials.TransverselyIsotropic);
        for i = 1:length(fieldnames(a.Materials.TransverselyIsotropic))
            F = struct2cell(getfield(a.Materials.TransverselyIsotropic,E{i}))';
            G(i,1:7) = [F(1) F(3:5) F(7) F(10) F(12)];
        end
        writetable(cell2table(G),fullfile(a.MaterialListDirectory,'MaterialList_TransverselyIsotropic.txt'),'delimiter',' ','WriteVariableNames',false)
    elseif a.MaterialTypeME == 3
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'Name',a.MaterialName4B);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'Class','Cubic');
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'Density',a.Material2ME.Density);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'E1',a.Material2ME.E1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'E2',a.Material2ME.E1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'E3',a.Material2ME.E1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'G12',a.Material2ME.G12);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'G13',a.Material2ME.G12);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'G23',a.Material2ME.G12);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'v12',a.Material2ME.v12);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'v13',a.Material2ME.v12);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'v23',a.Material2ME.v12);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_1',a.Material2ME.LongitudinalVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_2',a.Material2ME.LongitudinalVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'LongitudinalVelocity_3',a.Material2ME.LongitudinalVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_1',a.Material2ME.FastShearVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_2',a.Material2ME.FastShearVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'FastShearVelocity_3',a.Material2ME.FastShearVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_1',a.Material2ME.FastShearVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_2',a.Material2ME.FastShearVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'SlowShearVelocity_3',a.Material2ME.FastShearVelocity_1);
        a.Materials.Cubic = setfield(a.Materials.Cubic,{1,1},a.MaterialName4B,{1,1},'C',[a.Material2ME.C(1,1) a.Material2ME.C(1,2) a.Material2ME.C(1,2) 0 0 0;0 a.Material2ME.C(1,1) a.Material2ME.C(1,2) 0 0 0;0 0 a.Material2ME.C(1,1) 0 0 0;0 0 0 a.Material2ME.G12 0 0;0 0 0 0 a.Material2ME.G12 0;0 0 0 0 0 a.Material2ME.G12]);
        a.Materials.Cubic = orderfields(a.Materials.Cubic);
        a.Material2UI5.String = fieldnames(a.Materials.Cubic);
        a.MaterialUI8.String = fieldnames(a.Materials.Cubic);
        if  a.MaterialType_Bulk == 3
            a.MaterialUI8.String = fieldnames(a.Materials.Cubic);
        end
        if  a.SolidType_Bulk == 2
            a.SolidUI8.String = fieldnames(a.Materials.Cubic);
        end
        E = fieldnames(a.Materials.Cubic);
        for i = 1:length(fieldnames(a.Materials.Cubic))
            F = struct2cell(getfield(a.Materials.Cubic,E{i}))';
            G(i,1:5) = [F(1) F(3:4) F(7) F(10)];
        end
        writetable(cell2table(G),fullfile(a.MaterialListDirectory,'MaterialList_Cubic.txt'),'delimiter',' ','WriteVariableNames',false)
    end
end
function Save3
    a.Materials.Fluid = setfield(a.Materials.Fluid,{1,1},a.MaterialName4C,{1,1},'Name',a.MaterialName4C);
    a.Materials.Fluid = setfield(a.Materials.Fluid,{1,1},a.MaterialName4C,{1,1},'Class','Fluid');
    a.Materials.Fluid = setfield(a.Materials.Fluid,{1,1},a.MaterialName4C,{1,1},'Density',a.Material3ME.Density);
    a.Materials.Fluid = setfield(a.Materials.Fluid,{1,1},a.MaterialName4C,{1,1},'Velocity',a.Material3ME.Velocity);
    a.Materials.Fluid = orderfields(a.Materials.Fluid);
    a.FluidUI1.String = fieldnames(a.Materials.Fluid);
    if  a.Quantity11 == 4
        a.Option1UI1.String = fieldnames(a.Materials.Fluid);
    end
    a.Material3UI5.String = fieldnames(a.Materials.Fluid);
    a.CouplantUI8.String = fieldnames(a.Materials.Fluid);
    E = fieldnames(a.Materials.Fluid);
    for i = 1:length(fieldnames(a.Materials.Fluid))
        F = struct2cell(getfield(a.Materials.Fluid,E{i}))';
        G(i,1:3) = [F(1) F(3:4)];
    end
    writetable(cell2table(G(:,1:3)),fullfile(a.MaterialListDirectory,'MaterialList_Fluid.txt'),'delimiter',' ','WriteVariableNames',false)
end
end