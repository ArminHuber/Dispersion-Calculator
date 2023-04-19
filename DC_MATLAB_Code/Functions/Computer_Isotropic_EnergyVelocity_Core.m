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
function [X,Counter] = Computer_Isotropic_EnergyVelocity_Core(X,ModeType,ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,Fluid,Thickness)
SamplesX3 = 50;

Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
if  ~FluidLoading
    Fluid.Density = 1e-10;
    Fluid.Velocity = 1e-10;
end
x3 = 0:Thickness/SamplesX3:Thickness;
for p = 1:length(X)
    if  strcmp(ModeType,'Lamb') || strcmp(ModeType,'Scholte')
        for q = 1:height(X{p})
            PhaseVelocity = X{p}(q,4)*1e3;
            Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(q,1)*1e3;
            k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
            k3(1) = sqrt(AngularFrequency^2/Material.LongitudinalVelocity_complex^2-k^2);
            k3(2) = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
            k3Fluid = sqrt(AngularFrequency^2/Fluid.Velocity^2-k^2);
            if  Viscoelastic && strcmp(ModeType,'Scholte')
                k3Fluid = -k3Fluid;
            end
            W = (Material.Density*AngularFrequency^2-(Lambda+2*Mu)*k^2-Mu*k3.^2)./((Lambda+Mu)*k*k3); 
            WFluid = k3Fluid/k;
            D1 = 1i*(k*(Lambda+2*Mu)+Lambda*k3.*W); % sigma11
            D3 = 1i*(k*Lambda+(Lambda+2*Mu)*k3.*W); % sigma33
            D5 = 1i*Mu*(k3+k*W); % sigma13
            DFluid = 1i*Fluid.Density*AngularFrequency^2/k; % sigma11, sigma22, sigma33 in the fluid
            E1 = 1i*k; % epsilon11
            E3 = 1i*k3.*W; % epsilon33
            E5 = 1i*(k3+k*W); % epsilon13
            E = exp(1i*k3*Thickness);
            Z1 = [-W(2) W.*E -WFluid 0;-D3(2) -D3.*E DFluid 0;-D5(2) D5.*E 0 0;-W(2)*E(2) W 0 WFluid;-D3(2)*E(2) -D3 0 DFluid];
            Z2 = [W(1);D3(1);D5(1);W(1)*E(1);D3(1)*E(1)];
            U = Z1\Z2;
            E = [exp(1i*k3.*x3') exp(1i*k3.*(Thickness-x3)')];
            v(:,1) = -1i*AngularFrequency*(E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4));
            v(:,3) = -1i*AngularFrequency*(W(1)*E(:,1)+W(2)*U(1)*E(:,2)-W(1)*U(2)*E(:,3)-W(2)*U(3)*E(:,4));
            sigma(:,1) = D1(1)*E(:,1)+D1(2)*U(1)*E(:,2)+D1(1)*U(2)*E(:,3)+D1(2)*U(3)*E(:,4); % sigma11
            sigma(:,3) = D3(1)*E(:,1)+D3(2)*U(1)*E(:,2)+D3(1)*U(2)*E(:,3)+D3(2)*U(3)*E(:,4); % sigma33
            sigma(:,5) = D5(1)*E(:,1)+D5(2)*U(1)*E(:,2)-D5(1)*U(2)*E(:,3)-D5(2)*U(3)*E(:,4); % sigma13
            epsilon(:,1) = E1*(E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4)); % epsilon11
            epsilon(:,3) = E3(1)*E(:,1)+E3(2)*U(1)*E(:,2)+E3(1)*U(2)*E(:,3)+E3(2)*U(3)*E(:,4); % epsilon33
            epsilon(:,5) = E5(1)*E(:,1)+E5(2)*U(1)*E(:,2)-E5(1)*U(2)*E(:,3)-E5(2)*U(3)*E(:,4); % epsilon13
            StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+real(epsilon(:,1)).*real(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,3)).*real(sigma(:,3))+imag(epsilon(:,5)).*imag(sigma(:,5))+real(epsilon(:,5)).*real(sigma(:,5)));
            KineticEnergyDensity = .5*Material.Density*(abs(v(:,1)).^2+abs(v(:,3)).^2);
            PowerFlowDensity = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
            PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
            TotalEnergy = x3(2)*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
            X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % ce1
        end
    elseif contains(ModeType,'Shear')
        for q = 1:height(X{p})
            PhaseVelocity = X{p}(q,4)*1e3;
            Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(q,1)*1e3;
            k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
            k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
            D4 = 1i*k3*Mu; % sigma23
            D6 = 1i*k*Mu; % sigma12
            E4 = 1i*k3; % epsilon23
            E6 = 1i*k; % epsilon12
            E = [exp(1i*k3*x3') exp(1i*k3*(Thickness-x3)')];
            if  strcmp(ModeType,'SShear')
                v(:,2) = -1i*AngularFrequency*(E(:,1)+E(:,2));
                sigma(:,4) = D4*E(:,1)-D4*E(:,2); % sigma23
                sigma(:,6) = D6*E(:,1)+D6*E(:,2); % sigma12
                epsilon(:,4) = E4*E(:,1)-E4*E(:,2); % epsilon23
                epsilon(:,6) = E6*E(:,1)+E6*E(:,2); % epsilon12
            elseif strcmp(ModeType,'AShear')
                v(:,2) = -1i*AngularFrequency*(E(:,1)-E(:,2));
                sigma(:,4) = D4*E(:,1)+D4*E(:,2); % sigma23
                sigma(:,6) = D6*E(:,1)-D6*E(:,2); % sigma12
                epsilon(:,4) = E4*E(:,1)+E4*E(:,2); % epsilon23
                epsilon(:,6) = E6*E(:,1)-E6*E(:,2); % epsilon12
            end
            StrainEnergyDensity = .5*(imag(epsilon(:,4)).*imag(sigma(:,4))+real(epsilon(:,4)).*real(sigma(:,4))+imag(epsilon(:,6)).*imag(sigma(:,6))+real(epsilon(:,6)).*real(sigma(:,6)));
            KineticEnergyDensity = .5*Material.Density*abs(v(:,2)).^2;
            PowerFlowDensity = -.5*real(sigma(:,6).*conj(v(:,2)));
            PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
            TotalEnergy = x3(2)*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
            X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % ce1
        end
    end
    X{p}(:,5) = filloutliers(X{p}(:,5),'spline','movmedian',5,'ThresholdFactor',1);
    Counter = Counter+1;
    if  ModeTotal > 0
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
    end
end