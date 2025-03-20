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
function [Time,u,x1,x3,p] = ModeShapeGridAnimationComputer_Isotropic(Material,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,ALamb,AShear,AScholte,BLamb,BScholte,CycleDuration,FrameRate,Frequency,Length,Mode,Thickness,SamplesX1,SamplesX3,SLamb,SShear,SScholte,ShowHalfSpaces,HalfSpaces)
%#ok<*AGROW>
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
x3 = 0:Thickness/SamplesX3:Thickness; % (m)
Time = (0:1/(Frequency*1e3*CycleDuration*FrameRate):1/(Frequency*1e3))';
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
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
    AngularFrequency2 = AngularFrequency^2;
    k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
    k2 = k^2;
    k3(1) = sqrt(AngularFrequency2/Material.LongitudinalVelocity_complex^2-k2);
    k3(2) = sqrt(AngularFrequency2/Material.TransverseVelocity_complex^2-k2);
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
    W = (Material.Density*AngularFrequency2-(Lambda+2*Mu)*k2-Mu*k3.^2)./((Lambda+Mu)*k*k3); 
    WUpperFluid = k3UpperFluid/k;
    WLowerFluid = k3LowerFluid/k;
    D3 = 1i*(k*Lambda+(Lambda+2*Mu)*k3.*W); % sigma33
    D5 = 1i*Mu*(k3+k*W); % sigma13
    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/k; % sigma11, sigma22, sigma33 in the upper fluid
    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/k; % in the lower fluid
    E = exp(1i*k3*Thickness);
    Z1 = [-W(2) W.*E -WUpperFluid 0;-D3(2) -D3.*E DUpperFluid 0;-D5(2) D5.*E 0 0;-W(2)*E(2) W 0 WLowerFluid;-D3(2)*E(2) -D3 0 DLowerFluid];
    Z2 = [W(1);D3(1);D5(1);W(1)*E(1);D3(1)*E(1)];
    U = Z1\Z2;
    for i = 1:length(x1)
        for j = 1:length(x3)
            E = [exp(1i*(k*x1(i)+k3.*x3(j)-AngularFrequency*Time)) exp(1i*(k*x1(i)+k3.*(Thickness-x3(j))-AngularFrequency*Time))];
            u{j,i}(:,1) = E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4);
            u{j,i}(:,2) = W(1)*E(:,1)+W(2)*U(1)*E(:,2)-W(1)*U(2)*E(:,3)-W(2)*U(3)*E(:,4);
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
                [~,q] = min(abs(SShear{p}(:,1)-Frequency));
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
                E = [exp(1i*(k*x1(i)+k3*x3(j)-AngularFrequency*Time)) exp(1i*(k*x1(i)+k3*(Thickness-x3(j))-AngularFrequency*Time))];
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
                [~,q] = min(abs(AShear{p-1}(:,1)-Frequency));
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
                E = [exp(1i*(k*x1(i)+k3*x3(j)-AngularFrequency*Time)) exp(1i*(k*x1(i)+k3*(Thickness-x3(j))-AngularFrequency*Time))];
                u{j,i}(:,1) = E(:,1)-E(:,2);
                u{j,i}(:,2) = 0;
            end
        end
    end
end
if  FluidLoading && ShowHalfSpaces
    x3Fluid = 0:Thickness/SamplesX3:HalfSpaces*Thickness;
    if  ToggleUpperFluid
        if  ~contains(Mode,'SH')
            for i = 1:length(x1)
                for j = 1:length(x3Fluid)
                    EUpperFluid = exp(1i*(k*x1(i)+k3UpperFluid*x3Fluid(j)-AngularFrequency*Time));
                    uUpperFluid{j,i}(:,1) = U(4)*EUpperFluid;
                    uUpperFluid{j,i}(:,2) = -WUpperFluid*U(4)*EUpperFluid;
                end
            end
        elseif contains(Mode,'SH')
            for i = 1:length(x1)
                for j = 1:length(x3Fluid)
                    uUpperFluid{j,i}(length(Time),2) = 0;
                end
            end
        end
        x3 = horzcat(-fliplr(x3Fluid),x3);
        u = vertcat(flipud(uUpperFluid),u);
    end
    if  ToggleLowerFluid
        if  ~contains(Mode,'SH')
            for i = 1:length(x1)
                for j = 1:length(x3Fluid)
                    ELowerFluid = exp(1i*(k*x1(i)+k3LowerFluid*x3Fluid(j)-AngularFrequency*Time));
                    uLowerFluid{j,i}(:,1) = U(5)*ELowerFluid;
                    uLowerFluid{j,i}(:,2) = WLowerFluid*U(5)*ELowerFluid;
                end
            end
        elseif contains(Mode,'SH')
            for i = 1:length(x1)
                for j = 1:length(x3Fluid)
                    uLowerFluid{j,i}(length(Time),2) = 0;
                end
            end
        end
        x3 = horzcat(x3,x3(end)+x3Fluid);
        u = vertcat(u,uLowerFluid);
    end  
end