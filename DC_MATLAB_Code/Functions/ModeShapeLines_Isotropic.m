% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2025 DLR
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
function ModeShapeLines_Isotropic(FunctionMode,MakeFigure,ExportData,XSLX,TXT,MAT,Plot,PNGresolution,Material,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Color1,Color2,Color3,Color4,Color5,Color6,ALamb,AShear,AScholte,BLamb,BScholte,BoxLineWidth,Directory,FileName,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeLegend,Frequency,HeadLine,LegendLocation,LineWidth,Mode,PDF,PNG,Thickness,SLamb,SShear,SScholte,SamplesX3,ShowHalfSpaces,HalfSpaces,Phase)
%#ok<*FNDSB>
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
x3 = (0:Thickness/SamplesX3:Thickness)'; % (m)
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'SH')
    if  ~contains(Mode,'Scholte') && Mode(1) == 'S'
        q = find(SLamb{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SLamb{p}(:,1))) || Frequency < min(SLamb{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SLamb{p}(:,1))),' and ',num2str(ceil(max(SLamb{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(SLamb{p}(:,1)-Frequency));
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
                [~,q] = min(abs(ALamb{p}(:,1)-Frequency));
                Frequency = ALamb{p}(q,1);
            end
        end
        PhaseVelocity = ALamb{p}(q,4)*1e3;
        Attenuation = ALamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif ~contains(Mode,'Scholte') && Mode(1) == 'B'
        q = find(BLamb{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(BLamb{p}(:,1))) || Frequency < min(BLamb{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BLamb{p}(:,1))),' and ',num2str(ceil(max(BLamb{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(BLamb{p}(:,1)-Frequency));
                Frequency = BLamb{p}(q,1);
            end
        end
        PhaseVelocity = BLamb{p}(q,4)*1e3;
        Attenuation = BLamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'S' 
        q = find(SScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SScholte{p}(:,1))) || Frequency < min(SScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SScholte{p}(:,1))),' and ',num2str(ceil(max(SScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(SScholte{p}(:,1)-Frequency));
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
                [~,q] = min(abs(AScholte{p}(:,1)-Frequency));
                Frequency = AScholte{p}(q,1);
            end
        end
        PhaseVelocity = AScholte{p}(q,4)*1e3;
        Attenuation = AScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'B'
        q = find(BScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(BScholte{p}(:,1))) || Frequency < min(BScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BScholte{p}(:,1))),' and ',num2str(ceil(max(BScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(BScholte{p}(:,1)-Frequency));
                Frequency = BScholte{p}(q,1);
            end
        end
        PhaseVelocity = BScholte{p}(q,4)*1e3;
        Attenuation = BScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    end
    if  ~ToggleUpperFluid
        UpperFluid.Velocity = 1e-10;
        UpperFluid.Density = 1e-10;
    end
    if  ~ToggleLowerFluid
        LowerFluid.Velocity = 1e-10;
        LowerFluid.Density = 1e-10;
    end
    AngularFrequency = 2*pi*Frequency*1e3;
    AngularFrequency2 = AngularFrequency^2;
    k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
    k2 = k^2;
    k3(1) = sqrt(AngularFrequency^2/Material.LongitudinalVelocity_complex^2-k2);
    k3(2) = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k2);
    k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-k2);
    k3LowerFluid = sqrt(AngularFrequency2/LowerFluid.Velocity^2-k2);
    if  contains(Mode,'Scholte') && Attenuation ~= 0
        if  PhaseVelocity < UpperFluid.Velocity
            k3UpperFluid = -k3UpperFluid;
        end
        if  PhaseVelocity < LowerFluid.Velocity
            k3LowerFluid = -k3LowerFluid;
        end
    end
    W = (Material.Density*AngularFrequency^2-(Lambda+2*Mu)*k2-Mu*k3.^2)./((Lambda+Mu)*k*k3); 
    WUpperFluid = k3UpperFluid/k;
    WLowerFluid = k3LowerFluid/k;
    D1 = 1i*(k*(Lambda+2*Mu)+Lambda*k3.*W); % sigma11
    D2 = 1i*Lambda*(k+k3.*W); % sigma22
    D3 = 1i*(k*Lambda+(Lambda+2*Mu)*k3.*W); % sigma33
    D5 = 1i*Mu*(k3+k*W); % sigma13
    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/k; % sigma11, sigma22, sigma33 in the upper fluid
    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/k; % in the lower fluid
    E1 = 1i*k; % epsilon11
    E3 = 1i*k3.*W; % epsilon33
    E5 = 1i*(k3+k*W); % epsilon13
    E = exp(1i*k3*Thickness);
    Z1 = [-W(2) W.*E -WUpperFluid 0;-D3(2) -D3.*E DUpperFluid 0;-D5(2) D5.*E 0 0;-W(2)*E(2) W 0 WLowerFluid;-D3(2)*E(2) -D3 0 DLowerFluid];
    Z2 = [W(1);D3(1);D5(1);W(1)*E(1);D3(1)*E(1)];
    U = Z1\Z2;
    E = [exp(1i*k3.*x3) exp(1i*k3.*(Thickness-x3))];
    u(:,1) = E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4);
    u(:,3) = W(1)*E(:,1)+W(2)*U(1)*E(:,2)-W(1)*U(2)*E(:,3)-W(2)*U(3)*E(:,4);
    v = -1i*AngularFrequency*u;
    sigma(:,1) = D1(1)*E(:,1)+D1(2)*U(1)*E(:,2)+D1(1)*U(2)*E(:,3)+D1(2)*U(3)*E(:,4); % sigma11
    sigma(:,2) = D2(1)*E(:,1)+D2(2)*U(1)*E(:,2)+D2(1)*U(2)*E(:,3)+D2(2)*U(3)*E(:,4); % sigma22
    sigma(:,3) = D3(1)*E(:,1)+D3(2)*U(1)*E(:,2)+D3(1)*U(2)*E(:,3)+D3(2)*U(3)*E(:,4); % sigma33
    sigma(:,5) = D5(1)*E(:,1)+D5(2)*U(1)*E(:,2)-D5(1)*U(2)*E(:,3)-D5(2)*U(3)*E(:,4); % sigma13
    sigma(:,6) = 0;
    epsilon(:,1) = E1*(E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4)); % epsilon11
    epsilon(:,3) = E3(1)*E(:,1)+E3(2)*U(1)*E(:,2)+E3(1)*U(2)*E(:,3)+E3(2)*U(3)*E(:,4); % epsilon33
    epsilon(:,5) = E5(1)*E(:,1)+E5(2)*U(1)*E(:,2)-E5(1)*U(2)*E(:,3)-E5(2)*U(3)*E(:,4); % epsilon13
    epsilon(:,6) = 0;
    StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,5).*conj(sigma(:,5)));
    KineticEnergyDensity = .5*Material.Density*(abs(v(:,1)).^2+abs(v(:,3)).^2);
    PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
    PowerFlowDensity(:,3) = -.5*real(sigma(:,5).*conj(v(:,1))+sigma(:,3).*conj(v(:,3)));
    if  Phase
        sigmaPhase5 = rad2deg(angle(sigma(:,5)*exp(-1i*angle(u(1,1)))));
        epsilonPhase5 = rad2deg(angle(epsilon(:,5)*exp(-1i*angle(u(1,1)))));
    end
    sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5))); % zero in the fluid, so we make phase shift now    
    epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
elseif contains(Mode,'SH')
    if  Mode(1) == 'S'
        q = find(SShear{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SShear{p}(:,1))) || Frequency < min(SShear{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SShear{p}(:,1))),' and ',num2str(ceil(max(SShear{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(SShear{p}(:,1)-Frequency));
                Frequency = SShear{p}(q,1);
            end
        end
        PhaseVelocity = SShear{p}(q,4)*1e3;
        Attenuation = SShear{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
        D4 = 1i*k3*Mu; % sigma23
        D6 = 1i*k*Mu; % sigma12
        E4 = 1i*k3; % epsilon23
        E6 = 1i*k; % epsilon12
        E = [exp(1i*k3*x3) exp(1i*k3*(Thickness-x3))];
        u(:,2) = E(:,1)+E(:,2);
        v = -1i*AngularFrequency*u;
        sigma(:,4) = D4*E(:,1)-D4*E(:,2); % sigma23
        sigma(:,6) = D6*E(:,1)+D6*E(:,2); % sigma12
        epsilon(:,4) = E4*E(:,1)-E4*E(:,2); % epsilon23
        epsilon(:,6) = E6*E(:,1)+E6*E(:,2); % epsilon12
    elseif Mode(1) == 'A'
        q = find(AShear{p-1}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(AShear{p-1}(:,1))) || Frequency < min(AShear{p-1}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AShear{p-1}(:,1))),' and ',num2str(ceil(max(AShear{p-1}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(AShear{p-1}(:,1)-Frequency));
                Frequency = AShear{p-1}(q,1);
            end
        end
        PhaseVelocity = AShear{p-1}(q,4)*1e3;
        Attenuation = AShear{p-1}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
        D4 = 1i*k3*Mu; % sigma23
        D6 = 1i*k*Mu; % sigma12
        E4 = 1i*k3; % epsilon23
        E6 = 1i*k; % epsilon12
        E = [exp(1i*k3*x3) exp(1i*k3*(Thickness-x3))];
        u(:,2) = E(:,1)-E(:,2);
        v = -1i*AngularFrequency*u;
        sigma(:,4) = D4*E(:,1)+D4*E(:,2); % sigma23
        sigma(:,6) = D6*E(:,1)-D6*E(:,2); % sigma12
        epsilon(:,4) = E4*E(:,1)+E4*E(:,2); % epsilon23
        epsilon(:,6) = E6*E(:,1)-E6*E(:,2); % epsilon12
    end
    u(:,3) = 0;
    StrainEnergyDensity = .5*real(epsilon(:,4).*conj(sigma(:,4))+epsilon(:,6).*conj(sigma(:,6)));
    KineticEnergyDensity = .5*Material.Density*abs(v(:,2)).^2;
    PowerFlowDensity(:,1) = -.5*real(sigma(:,6).*conj(v(:,2)));
    PowerFlowDensity(:,3) = -.5*real(sigma(:,4).*conj(v(:,2)));
    if  Phase
        uPhase = rad2deg(angle(u));
        sigmaPhase = rad2deg(angle(sigma));
        epsilonPhase = rad2deg(angle(epsilon));
    end
    u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
    sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
    sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
    epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
    epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
end
PowerFlow = trapz(x3,PowerFlowDensity(:,1));
if  FluidLoading && ShowHalfSpaces
    x3Fluid = (0:Thickness/SamplesX3:HalfSpaces*Thickness)';
    E1Fluid = 1i*k; % epsilon11
    if  ToggleUpperFluid
        if  ~contains(Mode,'SH')
            EUpperFluid = exp(1i*k3UpperFluid*x3Fluid);
            E3UpperFluid = 1i*k3UpperFluid*WUpperFluid; % epsilon33
            uUpperFluid(:,1) = U(4)*EUpperFluid;
            uUpperFluid(:,3) = -WUpperFluid*U(4)*EUpperFluid;
            vUpperFluid = -1i*AngularFrequency*uUpperFluid;
            sigmaUpperFluid(:,1) = DUpperFluid*U(4)*EUpperFluid; % sigma11
            sigmaUpperFluid(:,2) = sigmaUpperFluid(:,1); % sigma22
            sigmaUpperFluid(:,3) = sigmaUpperFluid(:,1); % sigma33
            sigmaUpperFluid(:,6) = 0;
            epsilonUpperFluid(:,1) = E1Fluid*U(4)*EUpperFluid; % epsilon11
            epsilonUpperFluid(:,3) = E3UpperFluid*U(4)*EUpperFluid; % epsilon33
            epsilonUpperFluid(:,6) = 0;
            StrainEnergyDensityUpperFluid = .5*real(epsilonUpperFluid(:,1).*conj(sigmaUpperFluid(:,1))+epsilonUpperFluid(:,3).*conj(sigmaUpperFluid(:,3)));
            KineticEnergyDensityUpperFluid = .5*(UpperFluid.Density*(abs(vUpperFluid(:,1)).^2+abs(vUpperFluid(:,3)).^2));
            PowerFlowDensityUpperFluid(:,1) = -.5*real(sigmaUpperFluid(:,1).*conj(vUpperFluid(:,1)));
            PowerFlowDensityUpperFluid(:,3) = -.5*real(sigmaUpperFluid(:,3).*conj(vUpperFluid(:,3)));
        elseif contains(Mode,'SH')
            uUpperFluid(length(x3Fluid),3) = 0;
            sigmaUpperFluid(length(x3Fluid),6) = 0;
            epsilonUpperFluid(length(x3Fluid),6) = 0;
            StrainEnergyDensityUpperFluid(length(x3Fluid),1) = 0;
            KineticEnergyDensityUpperFluid(length(x3Fluid),1) = 0;
            PowerFlowDensityUpperFluid(length(x3Fluid),3) = 0;
            if  Phase
                uPhase = vertcat(zeros(length(x3Fluid),3),uPhase);
                sigmaPhase = vertcat(zeros(length(x3Fluid),6),sigmaPhase);
                epsilonPhase = vertcat(zeros(length(x3Fluid),6),epsilonPhase);
            end
        end
        x3 = vertcat(-flipud(x3Fluid),x3);
        u = vertcat(flipud(uUpperFluid),u);
        sigma = vertcat(flipud(sigmaUpperFluid),sigma);
        epsilon = vertcat(flipud(epsilonUpperFluid),epsilon);
        StrainEnergyDensity = vertcat(flipud(StrainEnergyDensityUpperFluid),StrainEnergyDensity);
        KineticEnergyDensity = vertcat(flipud(KineticEnergyDensityUpperFluid),KineticEnergyDensity);
        PowerFlowDensity = vertcat(flipud(PowerFlowDensityUpperFluid),PowerFlowDensity);
    end
    if  ToggleLowerFluid
        if  ~contains(Mode,'SH')
            ELowerFluid = exp(1i*k3LowerFluid*x3Fluid);
            E3LowerFluid = 1i*k3LowerFluid*WLowerFluid; % epsilon33
            uLowerFluid(:,1) = U(5)*ELowerFluid;
            uLowerFluid(:,3) = WLowerFluid*U(5)*ELowerFluid;
            vLowerFluid = -1i*AngularFrequency*uLowerFluid;
            sigmaLowerFluid(:,1) = DLowerFluid*U(5)*ELowerFluid; % sigma11
            sigmaLowerFluid(:,2) = sigmaLowerFluid(:,1); % sigma22
            sigmaLowerFluid(:,3) = sigmaLowerFluid(:,1); % sigma33
            sigmaLowerFluid(:,6) = 0;
            epsilonLowerFluid(:,1) = E1Fluid*U(5)*ELowerFluid; % epsilon11
            epsilonLowerFluid(:,3) = E3LowerFluid*U(5)*ELowerFluid; % epsilon33
            epsilonLowerFluid(:,6) = 0;
            StrainEnergyDensityLowerFluid = .5*real(epsilonLowerFluid(:,1).*conj(sigmaLowerFluid(:,1))+epsilonLowerFluid(:,3).*conj(sigmaLowerFluid(:,3)));
            KineticEnergyDensityLowerFluid = .5*(LowerFluid.Density*(abs(vLowerFluid(:,1)).^2+abs(vLowerFluid(:,3)).^2));
            PowerFlowDensityLowerFluid(:,1) = -.5*real(sigmaLowerFluid(:,1).*conj(vLowerFluid(:,1)));
            PowerFlowDensityLowerFluid(:,3) = -.5*real(sigmaLowerFluid(:,3).*conj(vLowerFluid(:,3)));
        elseif contains(Mode,'SH')
            uLowerFluid(length(x3Fluid),3) = 0;
            sigmaLowerFluid(length(x3Fluid),6) = 0;
            epsilonLowerFluid(length(x3Fluid),6) = 0;
            StrainEnergyDensityLowerFluid(length(x3Fluid),1) = 0;
            KineticEnergyDensityLowerFluid(length(x3Fluid),1) = 0;
            PowerFlowDensityLowerFluid(length(x3Fluid),3) = 0;
            if  Phase
                uPhase = vertcat(uPhase,zeros(length(x3Fluid),3));
                sigmaPhase = vertcat(sigmaPhase,zeros(length(x3Fluid),6));
                epsilonPhase = vertcat(epsilonPhase,zeros(length(x3Fluid),6));
            end
        end
        x3 = vertcat(x3,x3(end)+x3Fluid);
        u = vertcat(u,uLowerFluid);
        sigma = vertcat(sigma,sigmaLowerFluid);
        epsilon = vertcat(epsilon,epsilonLowerFluid);
        StrainEnergyDensity = vertcat(StrainEnergyDensity,StrainEnergyDensityLowerFluid);
        KineticEnergyDensity = vertcat(KineticEnergyDensity,KineticEnergyDensityLowerFluid);
        PowerFlowDensity = vertcat(PowerFlowDensity,PowerFlowDensityLowerFluid);
    end
end
if  ~contains(Mode,'SH')
    if  Phase
        if  FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid
                uPhase = rad2deg(angle(u*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                sigmaPhase(length(x3Fluid)+1:length(x3Fluid)+length(sigmaPhase5),5) = sigmaPhase5;
                epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(length(x3Fluid)+1,1)))));
                epsilonPhase(length(x3Fluid)+1:length(x3Fluid)+length(epsilonPhase5),5) = epsilonPhase5;
            else
                uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
                sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
                sigmaPhase(1:length(sigmaPhase5),5) = sigmaPhase5;
                epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
                epsilonPhase(1:length(epsilonPhase5),5) = epsilonPhase5;
            end
        else
            uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
            sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
            sigmaPhase(:,5) = sigmaPhase5;
            epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
            epsilonPhase(:,5) = epsilonPhase5;
        end
    end
    u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
    u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
    sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
    sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
    sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
    epsilon(:,1) = epsilon(:,1)*exp(-1i*angle(epsilon(2,1)));
    epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
end
u = u/sqrt(PowerFlow);
sigma = sigma/sqrt(PowerFlow);
epsilon = epsilon/sqrt(PowerFlow);
if  real(u(2,1)) < 0
    u = -u;
end
if  real(sigma(2,1)) < 0
    sigma = -sigma;
end
if  real(epsilon(2,1)) < 0
    epsilon = -epsilon;
end
if  Phase
    uPhase(find(round(uPhase) == -180)) = 180;
    sigmaPhase(find(round(sigmaPhase) == -180)) = 180;
    epsilonPhase(find(round(epsilonPhase) == -180)) = 180;
end
StrainEnergyDensity = .5*StrainEnergyDensity/PowerFlow;
KineticEnergyDensity = .5*KineticEnergyDensity/PowerFlow;
TotalEnergyDensity = StrainEnergyDensity+KineticEnergyDensity;
PowerFlowDensity = PowerFlowDensity/PowerFlow;
if  ExportData
    if  FunctionMode == 1
        if  ~Phase
            Table = table('Size',[length(x3) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','u1 (nm)','u2 (nm)','u3 (nm)'});
            Table(:,1) = num2cell(-1e3*x3);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
        else
            Table = table('Size',[length(x3) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','u1 (nm)','u2 (nm)','u3 (nm)','Phase u1 (deg)','Phase u2 (deg)','Phase u3 (deg)'});
            Table(:,1) = num2cell(-1e3*x3);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
            Table(:,5) = num2cell(uPhase(:,1));
            Table(:,6) = num2cell(uPhase(:,2));
            Table(:,7) = num2cell(uPhase(:,3));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.txt']));
            end
            if  MAT
                Displacement = Table;
                save(fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.mat']),'Displacement')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end  
    elseif FunctionMode == 2
        if  ~Phase
            Table = table('Size',[length(x3) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','sigma11 (kPa)','sigma22 (kPa)','sigma33 (kPa)','sigma23 (kPa)','sigma13 (kPa)','sigma12 (kPa)'});
            Table(:,1) = num2cell(-1e3*x3);
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
        else
            Table = table('Size',[length(x3) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','sigma11 (kPa)','sigma22 (kPa)','sigma33 (kPa)','sigma23 (kPa)','sigma13 (kPa)','sigma12 (kPa)','Phase sigma11 (deg)','Phase sigma22 (deg)','Phase sigma33 (deg)','Phase sigma23 (deg)','Phase sigma13 (deg)','Phase sigma12 (deg)'});
            Table(:,1) = num2cell(-1e3*x3);
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
            Table(:,8) = num2cell(sigmaPhase(:,1));
            Table(:,9) = num2cell(sigmaPhase(:,2));
            Table(:,10) = num2cell(sigmaPhase(:,3));
            Table(:,11) = num2cell(sigmaPhase(:,4));
            Table(:,12) = num2cell(sigmaPhase(:,5));
            Table(:,13) = num2cell(sigmaPhase(:,6));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.txt']));
            end
            if  MAT
                Stress = Table;
                save(fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.mat']),'Stress')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end   
    elseif FunctionMode == 3
        if  ~Phase
            Table = table('Size',[length(x3) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','epsilon11','epsilon22','epsilon33','epsilon23','epsilon13','epsilon12'});
            Table(:,1) = num2cell(-1e3*x3);
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
        else
            Table = table('Size',[length(x3) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','epsilon11','epsilon22','epsilon33','epsilon23','epsilon13','epsilon12','Phase epsilon11 (deg)','Phase epsilon22 (deg)','Phase epsilon33 (deg)','Phase epsilon23 (deg)','Phase epsilon13 (deg)','Phase epsilon12 (deg)'});
            Table(:,1) = num2cell(-1e3*x3);
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
            Table(:,8) = num2cell(epsilonPhase(:,1));
            Table(:,9) = num2cell(epsilonPhase(:,2));
            Table(:,10) = num2cell(epsilonPhase(:,3));
            Table(:,11) = num2cell(epsilonPhase(:,4));
            Table(:,12) = num2cell(epsilonPhase(:,5));
            Table(:,13) = num2cell(epsilonPhase(:,6));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.txt']));
            end
            if  MAT
                Strain = Table;
                save(fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.mat']),'Strain')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end   
    elseif FunctionMode == 4
        Table = table('Size',[length(x3) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','Estrain (J/m2)','Ekin (J/m2)','Etotal (J/m2)'});
        Table(:,1) = num2cell(-1e3*x3);
        Table(:,2) = num2cell(StrainEnergyDensity);
        Table(:,3) = num2cell(KineticEnergyDensity);
        Table(:,4) = num2cell(TotalEnergyDensity);
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.txt']));
            end
            if  MAT
                EnergyDensity = Table;
                save(fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.mat']),'EnergyDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end        
    elseif FunctionMode == 5
        Table = table('Size',[length(x3) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','P1 (W/m)','P2 (W/m)','P3 (W/m)'});
        Table(:,1) = num2cell(-1e3*x3);
        Table(:,2) = num2cell(PowerFlowDensity(:,1));
        Table(:,3) = num2cell(PowerFlowDensity(:,2));
        Table(:,4) = num2cell(PowerFlowDensity(:,3));
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.txt']));
            end
            if  MAT
                PowerFlowDensity = Table;
                save(fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*Thickness),'MHzmm.mat']),'PowerFlowDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    end
end
if  MakeFigure
    if  Phase
        if  FunctionMode == 1 
            f = figure('Name','Displacement phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        elseif FunctionMode == 2
            f = figure('Name','Stress phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        elseif FunctionMode == 3
            f = figure('Name','Strain phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        end   
        datacursormode on
        x = xline(0,'Color',[.6 .6 .6]); 
        hasbehavior(x,'legend',false);
        hold on
        if  FunctionMode == 1
            if  Plot(3)
                plot(uPhase(:,3),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),-x3*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),-x3*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),-x3*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),-x3*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),-x3*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),-x3*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end
        end    
        ax = gca;
        ax.Box = 'on';
        ax.LineWidth = BoxLineWidth;
        ax.FontSize = FontSizeAxes;
        ax.Title.Interpreter = 'latex';
        ax.Title.FontSize = FontSizeHeadLine;
        if  HeadLine
            if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
                ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
            elseif contains(Mode,'SH')
                ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
            elseif contains(Mode,'Scholte')
                ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
            end
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Thickness*1e3),'\,mm ',replace(Material.Name,'_','\_'),' plate'];
            if  FluidLoading
                if  ToggleUpperFluid && ToggleLowerFluid
                    String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
                end
            end
            ax.Title.String = String;
        end
        ax.XLabel.Interpreter = 'latex';
        ax.XLabel.FontSize = FontSizeAxesLabels;
        ax.XLim = max(abs(ax.XLim))*[-1 1];
        ax.YLabel.Interpreter = 'latex';
        ax.YLabel.FontSize = FontSizeAxesLabels;
        ax.YLabel.String = '$x_3$ (mm)';
        ax.TickLabelInterpreter = 'latex';
        ax.YLim = -1e3*[x3(end) x3(1)];
        if  FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid && ToggleLowerFluid
                if  strcmp(UpperFluid.Name,LowerFluid.Name)
                    UpperFaceAlpha = .2;
                    LowerFaceAlpha = .2;
                else
                    if  UpperFluid.Density*UpperFluid.Velocity > LowerFluid.Density*LowerFluid.Velocity
                        UpperFaceAlpha = .2;
                        LowerFaceAlpha = .1;
                    else
                        UpperFaceAlpha = .1;
                        LowerFaceAlpha = .2;
                    end
                end
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',UpperFaceAlpha,'EdgeColor','none')
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -Thickness*1e3 -Thickness*1e3],'b','FaceAlpha',LowerFaceAlpha,'EdgeColor','none')
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none')
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -Thickness*1e3 -Thickness*1e3],'b','FaceAlpha',.2,'EdgeColor','none')
            end
        end
        if  FunctionMode == 1
            if  Plot(3)
                plot(uPhase(:,3),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),-x3*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),-x3*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),-x3*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),-x3*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),-x3*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),-x3*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end
        end
        if  FunctionMode == 1
            ax.XLabel.String = 'Displacement phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$u_3$','$u_1$','$u_2$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
            else
                LegendNames = {'Out-of-plane ($u_3$)','In-plane ($u_1$)','Shear horizontal ($u_2$)'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')
            end
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        elseif FunctionMode == 2
            ax.XLabel.String = 'Stress phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$\sigma_{33}$','$\sigma_{11}$','$\sigma_{22}$','$\sigma_{13}$','$\sigma_{23}$','$\sigma_{12}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Out-of-plane ($\sigma_{33}$)','In-plane ($\sigma_{11}$)','In-plane ($\sigma_{22}$)','Shear ($\sigma_{13}$)','Shear ($\sigma_{23}$)','Shear ($\sigma_{12}$)'}; 
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
            end
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        elseif FunctionMode == 3
            ax.XLabel.String = 'Strain phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$\varepsilon_{33}$','$\varepsilon_{11}$','$\varepsilon_{22}$','$\varepsilon_{13}$','$\varepsilon_{23}$','$\varepsilon_{12}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Out-of-plane ($\varepsilon_{33}$)','In-plane ($\varepsilon_{11}$)','In-plane ($\varepsilon_{22}$)','Shear ($\varepsilon_{13}$)','Shear ($\varepsilon_{23}$)','Shear ($\varepsilon_{12}$)'}; 
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
            end 
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        end
        tb = axtoolbar('default');
        tb.Visible = 'on';
        d = datacursormode(f);
        d.Interpreter = 'latex';
        d.UpdateFcn = @CursorPhase;
    end
    if  FunctionMode == 1 
        f = figure('Name','Displacement','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 2
        f = figure('Name','Stress','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 3
        f = figure('Name','Strain','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 4
        f = figure('Name','Energy density','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 5 
        f = figure('Name','Power flow density','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    end   
    datacursormode on
    x = xline(0,'Color',[.6 .6 .6]); 
    hasbehavior(x,'legend',false);
    hold on
    if  FunctionMode == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),-x3*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),-x3*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),-x3*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
    end    
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.Interpreter = 'latex';
    ax.Title.FontSize = FontSizeHeadLine;
    if  HeadLine
        if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
            ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
        elseif contains(Mode,'SH')
            ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
        elseif contains(Mode,'Scholte')
            ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
        end
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Thickness*1e3),'\,mm ',replace(Material.Name,'_','\_'),' plate'];
        if  FluidLoading
            if  ToggleUpperFluid && ToggleLowerFluid
                String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
            end
        end
        ax.Title.String = String;
    end
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.XLim = max(abs(ax.XLim))*[-1 1];
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.String = '$x_3$ (mm)';    
    ax.TickLabelInterpreter = 'latex';    
    ax.YLim = -1e3*[x3(end) x3(1)];
    if  FluidLoading && ShowHalfSpaces
        if  ToggleUpperFluid && ToggleLowerFluid
            if  strcmp(UpperFluid.Name,LowerFluid.Name)
                UpperFaceAlpha = .2;
                LowerFaceAlpha = .2;
            else
                if  UpperFluid.Density*UpperFluid.Velocity > LowerFluid.Density*LowerFluid.Velocity
                    UpperFaceAlpha = .2;
                    LowerFaceAlpha = .1;
                else
                    UpperFaceAlpha = .1;
                    LowerFaceAlpha = .2;
                end
            end
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',UpperFaceAlpha,'EdgeColor','none')
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -Thickness*1e3 -Thickness*1e3],'b','FaceAlpha',LowerFaceAlpha,'EdgeColor','none')
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none')
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -Thickness*1e3 -Thickness*1e3],'b','FaceAlpha',.2,'EdgeColor','none')
        end
    end
    if  FunctionMode == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,-x3*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),-x3*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),-x3*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),-x3*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),-x3*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),-x3*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),-x3*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
    end
    if  FunctionMode == 1
        ax.XLabel.String = 'Displacement (nm)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$u_3$','$u_1$','$u_2$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
        else
            LegendNames = {'Out-of-plane ($u_3$)','In-plane ($u_1$)','Shear horizontal ($u_2$)'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')
        end
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 2
        ax.XLabel.String = 'Stress (kPa)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$\sigma_{33}$','$\sigma_{11}$','$\sigma_{22}$','$\sigma_{13}$','$\sigma_{23}$','$\sigma_{12}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($\sigma_{33}$)','In-plane ($\sigma_{11}$)','In-plane ($\sigma_{22}$)','Shear ($\sigma_{13}$)','Shear ($\sigma_{23}$)','Shear ($\sigma_{12}$)'}; 
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
        end
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 3
        ax.XLabel.String = 'Strain';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$\varepsilon_{33}$','$\varepsilon_{11}$','$\varepsilon_{22}$','$\varepsilon_{13}$','$\varepsilon_{23}$','$\varepsilon_{12}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($\varepsilon_{33}$)','In-plane ($\varepsilon_{11}$)','In-plane ($\varepsilon_{22}$)','Shear ($\varepsilon_{13}$)','Shear ($\varepsilon_{23}$)','Shear ($\varepsilon_{12}$)'}; 
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
        end 
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 4
        ax.XLabel.String = 'Energy density (J/m$^2$)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$E_\mathrm{strain}$','$E_\mathrm{kin}$','$E_\mathrm{total}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Strain energy density','Kinetic energy density','Total energy density'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')        
        end       
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end 
    elseif FunctionMode == 5
        ax.XLabel.String = 'Power flow density (W/m)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$p_3$','$p_1$','$p_2$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($p_3$)','In-plane ($p_1$)','Shear horizontal ($p_2$)'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')        
        end
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    end
    tb = axtoolbar('default');
    tb.Visible = 'on';
    d = datacursormode(f);
    d.Interpreter = 'latex';
    d.UpdateFcn = @Cursor;
end
function output_txt = Cursor(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = replace(UpperFluid.Name,'_','\_');
        else
            output_txt = replace(LowerFluid.Name,'_','\_');
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$u_1$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$u_2$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$u_3$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\sigma_{33}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\sigma_{13}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\sigma_{23}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\sigma_{11}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\sigma_{22}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\sigma_{12}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varepsilon_{33}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varepsilon_{13}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varepsilon_{23}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varepsilon_{11}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varepsilon_{22}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varepsilon_{12}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 4
            if  event_obj.Target.Color == Color1
                output_txt = {['$E_\mathrm{strain}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$E_\mathrm{kin}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$E_\mathrm{total}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 5 
            if  event_obj.Target.Color == Color1
                output_txt = {['$p_1$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$p_2$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$p_3$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        end
    end
end
function output_txt = CursorPhase(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = replace(UpperFluid.Name,'_','\_');
        else
            output_txt = replace(LowerFluid.Name,'_','\_');
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$\varphi(u_1)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(u_2)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$\varphi(u_3)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\sigma_{33})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\sigma_{13})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\sigma_{23})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\sigma_{11})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\sigma_{22})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\sigma_{12})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\varepsilon_{33})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\varepsilon_{13})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\varepsilon_{23})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\varepsilon_{11})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\varepsilon_{22})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\varepsilon_{12})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        end
    end
end
end