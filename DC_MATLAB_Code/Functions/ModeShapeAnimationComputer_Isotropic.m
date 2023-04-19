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
function [Time,u,x1,x3,p] = ModeShapeAnimationComputer_Isotropic(Material,Viscoelastic,FluidLoading,Fluid,ALamb,AShear,AScholte,CycleDuration,Cycles,FrameRate,Frequency,Length,Mode,Thickness,SamplesX1,SamplesX3,SLamb,SShear,SScholte,ShowHalfSpaces,HalfSpaces)
%#ok<*AGROW>
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
x3 = 0:Thickness/SamplesX3:Thickness;
Time = 0:Cycles/(Frequency*1e3*Cycles*CycleDuration*FrameRate):Cycles/(Frequency*1e3);
disp(['Movie duration: ',num2str(Cycles*CycleDuration),' s',newline,'Frames @ ',num2str(FrameRate),'/s: ',num2str(Cycles*CycleDuration*FrameRate)]);
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'SH')
    if  ~contains(Mode,'Scholte') && Mode(1) == 'S'
        q = find(SLamb{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SLamb{p}(:,1))) || Frequency < min(SLamb{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SLamb{p}(:,1))),' and ',num2str(ceil(max(SLamb{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(SLamb{p}(:,1)-Frequency) == min(abs(SLamb{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = SLamb{p}(q,1);
            end
        end
        PhaseVelocity = SLamb{p}(q,4)*1e3;
        Attenuation = SLamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif ~contains(Mode,'Scholte') && Mode(1) == 'A'
        q = find(ALamb{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(ALamb{p}(:,1))) || Frequency < min(ALamb{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(ALamb{p}(:,1))),' and ',num2str(ceil(max(ALamb{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(ALamb{p}(:,1)-Frequency) == min(abs(ALamb{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = ALamb{p}(q,1);
            end
        end
        PhaseVelocity = ALamb{p}(q,4)*1e3;
        Attenuation = ALamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'S' 
        q = find(SScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SScholte{p}(:,1))) || Frequency < min(SScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SScholte{p}(:,1))),' and ',num2str(ceil(max(SScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(SScholte{p}(:,1)-Frequency) == min(abs(SScholte{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = SScholte{p}(q,1);
            end
        end
        PhaseVelocity = SScholte{p}(q,4)*1e3;
        Attenuation = SScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'A'
        q = find(AScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(AScholte{p}(:,1))) || Frequency < min(AScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AScholte{p}(:,1))),' and ',num2str(ceil(max(AScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(AScholte{p}(:,1)-Frequency) == min(abs(AScholte{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = AScholte{p}(q,1);
            end
        end
        PhaseVelocity = AScholte{p}(q,4)*1e3;
        Attenuation = AScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    end
    if  ~FluidLoading
        Fluid.Density = 1e-10;
        Fluid.Velocity = 1e-10;
    end
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
    k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
    k3(1) = sqrt(AngularFrequency^2/Material.LongitudinalVelocity_complex^2-k^2);
    k3(2) = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
    k3Fluid = sqrt(AngularFrequency^2/Fluid.Velocity^2-k^2);
    if  Viscoelastic && contains(Mode,'Scholte')
        k3Fluid = -k3Fluid;
    end
    W = (Material.Density*AngularFrequency^2-(Lambda+2*Mu)*k^2-Mu*k3.^2)./((Lambda+Mu)*k*k3); 
    WFluid = k3Fluid/k;
    D3 = 1i*(k*Lambda+(Lambda+2*Mu)*k3.*W); % sigma33
    D5 = 1i*Mu*(k3+k*W); % sigma13
    DFluid = 1i*Fluid.Density*AngularFrequency^2/k; % sigma11, sigma22, sigma33 in the fluid
    E = exp(1i*k3*Thickness);
    Z1 = [-W(2) W.*E -WFluid 0;-D3(2) -D3.*E DFluid 0;-D5(2) D5.*E 0 0;-W(2)*E(2) W 0 WFluid;-D3(2)*E(2) -D3 0 DFluid];
    Z2 = [W(1);D3(1);D5(1);W(1)*E(1);D3(1)*E(1)];
    U = Z1\Z2;
    for i = 1:length(x1)
        for j = 1:length(x3)
            E = [exp(1i*(k*x1(i)+k3.*x3(j)-AngularFrequency*Time')) exp(1i*(k*x1(i)+k3.*(Thickness-x3(j))-AngularFrequency*Time'))];
            u{j,i}(:,1) = E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4);
            u{j,i}(:,2) = W(1)*E(:,1)+W(2)*U(1)*E(:,2)-W(1)*U(2)*E(:,3)-W(2)*U(3)*E(:,4);
        end
    end
    if  FluidLoading && ShowHalfSpaces
        x3Fluid = 0:Thickness/SamplesX3:HalfSpaces*Thickness;
        for i = 1:length(x1)
            for j = 1:length(x3Fluid)
                EFluid = exp(1i*(k*x1(i)+k3Fluid*x3Fluid(j)-AngularFrequency*Time));
                uFluid0{j,i}(:,1) = U(4)*EFluid;
                uFluid0{j,i}(:,2) = -WFluid*U(4)*EFluid;
                uFluid1{j,i}(:,1) = U(5)*EFluid;
                uFluid1{j,i}(:,2) = WFluid*U(5)*EFluid;
            end
        end
    end
elseif contains(Mode,'SH')
    if  Mode(1) == 'S'
        q = find(SShear{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SShear{p}(:,1))) || Frequency < min(SShear{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SShear{p}(:,1))),' and ',num2str(ceil(max(SShear{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(SShear{p}(:,1)-Frequency) == min(abs(SShear{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = SShear{p}(q,1);
            end
        end
        PhaseVelocity = SShear{p}(q,4)*1e3;
        Attenuation = SShear{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
        for i = 1:length(x1)
            for j = 1:length(x3)
                E = [exp(1i*(k*x1(i)+k3*x3(j)-AngularFrequency*Time')) exp(1i*(k*x1(i)+k3*(Thickness-x3(j))-AngularFrequency*Time'))];
                u{j,i}(:,1) = E(:,1)+E(:,2);
                u{j,i}(:,2) = 0;
            end
        end
    elseif Mode(1) == 'A'
        q = find(AShear{p-1}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(AShear{p-1}(:,1))) || Frequency < min(AShear{p-1}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AShear{p-1}(:,1))),' and ',num2str(ceil(max(AShear{p-1}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(AShear{p-1}(:,1)-Frequency) == min(abs(AShear{p-1}(:,1)-Frequency)));
                q = q(1);
                Frequency = AShear{p-1}(q,1);
            end
        end
        PhaseVelocity = AShear{p-1}(q,4)*1e3;
        Attenuation = AShear{p-1}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
        for i = 1:length(x1)
            for j = 1:length(x3)
                E = [exp(1i*(k*x1(i)+k3*x3(j)-AngularFrequency*Time')) exp(1i*(k*x1(i)+k3*(Thickness-x3(j))-AngularFrequency*Time'))];
                u{j,i}(:,1) = E(:,1)-E(:,2);
                u{j,i}(:,2) = 0;
            end
        end
    end
    if  FluidLoading && ShowHalfSpaces
        x3Fluid = 0:Thickness/SamplesX3:HalfSpaces*Thickness;
        for i = 1:length(x1)
            for j = 1:length(x3Fluid)
                uFluid0{j,i}(length(Time),2) = 0;
                uFluid1{j,i}(length(Time),2) = 0;
            end
        end
    end
end
if  FluidLoading && ShowHalfSpaces
    x3 = horzcat(-fliplr(x3Fluid),x3,x3(end)+x3Fluid);
    u = vertcat(flipud(uFluid0),u,uFluid1);
end