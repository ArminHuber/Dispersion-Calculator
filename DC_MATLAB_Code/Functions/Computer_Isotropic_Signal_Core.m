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
function [uSum,Counter,ExcitationSpectrumRange] = Computer_Isotropic_Signal_Core(ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,CoherenceTime,Counter,Distance,ExcitationSpectrum,SamplesPerCycle,Cycles,FrequencyResolution,FrequencyResolution2,h1,Material,ModeType,ModeTotal,PhaseVelocity,Attenuation,Thickness,Time,TimeLimit,z1,z2)
SamplesX3 = 50; % samples per layer

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
if  ~ToggleUpperFluid
    UpperFluid.Velocity = 1e-10;
    UpperFluid.Density = 1e-10;
end
if  ~ToggleLowerFluid
    LowerFluid.Velocity = 1e-10;
    LowerFluid.Density = 1e-10;
end
x3 = (0:Thickness/SamplesX3:Thickness)'; % (m)
Half = Thickness/2;
T = (-Half:2*Half/SamplesX3:Half)';
uSum{length(PhaseVelocity)} = [];
ExcitationSpectrumRange{length(PhaseVelocity)} = [];
Time0 = Time(1:SamplesPerCycle*Cycles);
% figure
% hold on
for p = 1:length(PhaseVelocity)
    if  Stop
        return
    end
    if  any(PhaseVelocity{p})
        ExcitationSpectrumRange{p} = ExcitationSpectrum(:,z1(p):z2(p)); % extract from the frequency spectrum the range for which we have phase velocities
        if  CoherenceTime < TimeLimit % if the coherence time is smaller than the calculated time range, we have to expect seeing unwanted twins of the wave packet
            Fit1 = fit(ExcitationSpectrumRange{p}(1,:)',ExcitationSpectrumRange{p}(2,:)','cubicspline'); % fit the spectral amplitudes
            Fit2 = fit(ExcitationSpectrumRange{p}(1,:)',PhaseVelocity{p},'cubicspline'); % fit the phase velocity
            Fit3 = fit(ExcitationSpectrumRange{p}(1,:)',Attenuation{p},'cubicspline'); % fit the attenuation
            ExcitationSpectrumRange{p} = ExcitationSpectrumRange{p}(1,1):FrequencyResolution2:ExcitationSpectrumRange{p}(1,end); % generate the new frequency range with smaller steps
            ExcitationSpectrumRange{p}(2,:) = Fit1(ExcitationSpectrumRange{p}(1,:))/FrequencyResolution*FrequencyResolution2; % interpolate the spectral amplitudes
            PhaseVelocity{p} = Fit2(ExcitationSpectrumRange{p}(1,:)); % interpolate the phase velocity
            Attenuation{p} = Fit3(ExcitationSpectrumRange{p}(1,:)); % interpolate the attenuation
        end
        try
            u = zeros(length(ExcitationSpectrumRange{p}),length(Time));
            u0 = zeros(length(ExcitationSpectrumRange{p}),length(Time0));
        catch ME
            msgbox(['IDENTIFIER: ',ME.identifier,newline,'MESSAGE: ',ME.message,newline,'FILE: ',ME.stack.file,newline])
            errordlg('The data size exceeds your system''s RAM! Decrease the parameter settings.','Error');
            if  exist('h1','var') && ModeTotal > 0
                close(h1)
            end
            return
        end
        if  contains(ModeType,'Lamb') % this cannot be done with the potential approach because the phase information is lost; by contrast, phase information is present in the bulk wave approach used for Lamb waves below; for SH waves, the potential approach can be used
            for i = 1:length(PhaseVelocity{p})
                AngularFrequency = 2*pi*ExcitationSpectrumRange{p}(1,i)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                k = AngularFrequency/PhaseVelocity{p}(i)*(1+1i*Attenuation{p}(i)/(2*pi));
                k2 = k^2;
                k3(1) = sqrt(AngularFrequency2/Material.LongitudinalVelocity_complex^2-k2);
                k3(2) = sqrt(AngularFrequency2/Material.TransverseVelocity_complex^2-k2);
                k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-k2);
                k3LowerFluid = sqrt(AngularFrequency2/LowerFluid.Velocity^2-k2);
                W = (Material.Density*AngularFrequency2-(Lambda+2*Mu)*k2-Mu*k3.^2)./((Lambda+Mu)*k*k3);
                WUpperFluid = k3UpperFluid/k;
                WLowerFluid = k3LowerFluid/k;
                D1 = 1i*(k*(Lambda+2*Mu)+Lambda*k3.*W); % sigma11
                D3 = 1i*(k*Lambda+(Lambda+2*Mu)*k3.*W); % sigma33
                D5 = 1i*Mu*(k3+k*W); % sigma13
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/k; % sigma11, sigma22, sigma33 in the upper fluid
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/k; % in the lower fluid
                E = exp(1i*k3*Thickness);
                Z1 = [-W(2) W.*E -WUpperFluid 0;-D3(2) -D3.*E DUpperFluid 0;-D5(2) D5.*E 0 0;-W(2)*E(2) W 0 WLowerFluid;-D3(2)*E(2) -D3 0 DLowerFluid];
                Z2 = [W(1);D3(1);D5(1);W(1)*E(1);D3(1)*E(1)];
                U = Z1\Z2;
                E = [exp(1i*k3.*x3) exp(1i*k3.*(Thickness-x3))];
                v(:,1) = -1i*AngularFrequency*(E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4));
                v(:,3) = -1i*AngularFrequency*(W(1)*E(:,1)+W(2)*U(1)*E(:,2)-W(1)*U(2)*E(:,3)-W(2)*U(3)*E(:,4));
                sigma(:,1) = D1(1)*E(:,1)+D1(2)*U(1)*E(:,2)+D1(1)*U(2)*E(:,3)+D1(2)*U(3)*E(:,4); % sigma11
                sigma(:,5) = D5(1)*E(:,1)+D5(2)*U(1)*E(:,2)-D5(1)*U(2)*E(:,3)-D5(2)*U(3)*E(:,4); % sigma13
                PowerFlowDensity = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
                PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                E = [exp(1i*(k*Distance+k3.*x3(1:2))) exp(1i*(k*Distance+k3.*(Thickness-x3(1:2))))];
                if  strcmp(ModeType,'ALamb') && PhaseVelocity{p}(1) > Material.LongitudinalVelocity && PhaseVelocity{p}(end) < Material.LongitudinalVelocity
                    Shift = -1i*angle(k3(1));
                else
                    Shift = 0;
                end
                if  DisplacementComponent == 1 % u3
                    unorm = W(1)*E(:,1)+W(2)*U(1)*E(:,2)-W(1)*U(2)*E(:,3)-W(2)*U(3)*E(:,4);
                    u0(i,:) = real(...
                         W(1)*     exp(Shift-1i*AngularFrequency*Time0)...
                        +W(2)*U(1)*exp(Shift-1i*AngularFrequency*Time0)...
                        -W(1)*U(2)*exp(Shift+1i*(k3(1)*Thickness-AngularFrequency*Time0))...
                        -W(2)*U(3)*exp(Shift+1i*(k3(2)*Thickness-AngularFrequency*Time0)))/sqrt(PowerFlow); 
                    u(i,:) = real(...
                         W(1)*     exp(Shift+1i*(k*Distance-AngularFrequency*Time))...
                        +W(2)*U(1)*exp(Shift+1i*(k*Distance-AngularFrequency*Time))...
                        -W(1)*U(2)*exp(Shift+1i*(k*Distance+k3(1)*Thickness-AngularFrequency*Time))...
                        -W(2)*U(3)*exp(Shift+1i*(k*Distance+k3(2)*Thickness-AngularFrequency*Time)))/sqrt(PowerFlow); 
                elseif DisplacementComponent == 2 % u1
                    unorm = E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4);
                    u0(i,:) = real(...
                              exp(Shift-1i*AngularFrequency*Time0)...
                        +U(1)*exp(Shift-1i*AngularFrequency*Time0)...
                        +U(2)*exp(Shift+1i*(k3(1)*Thickness-AngularFrequency*Time0))...
                        +U(3)*exp(Shift+1i*(k3(2)*Thickness-AngularFrequency*Time0)))/sqrt(PowerFlow);
                    u(i,:) = real(...
                              exp(Shift+1i*(k*Distance-AngularFrequency*Time))...
                        +U(1)*exp(Shift+1i*(k*Distance-AngularFrequency*Time))...
                        +U(2)*exp(Shift+1i*(k*Distance+k3(1)*Thickness-AngularFrequency*Time))...
                        +U(3)*exp(Shift+1i*(k*Distance+k3(2)*Thickness-AngularFrequency*Time)))/sqrt(PowerFlow); 
                end
                uNorm(i) = abs(real(unorm(1)*exp(-1i*angle(unorm(2))))/sqrt(PowerFlow));
                if  i > 1 && max(u0(i-1,:)+u0(i,:)) < max(u0(i-1,:)-u0(i,:))
                    u(i,:) = -u(i,:);
                    u0(i,:) = -u0(i,:);
                end
            end
            uNormSmooth = filloutliers(uNorm,'spline','movmedian',5,'ThresholdFactor',1);
% figure
% hold on
% plot(uNorm,'linewidth',4,'color','r')
% plot(uNormSmooth,'linewidth',1.5,'color','g')
            for i = 1:length(PhaseVelocity{p})
                u(i,:) = u(i,:)*uNormSmooth(i)/max(u(i,:))*ExcitationSpectrumRange{p}(2,i);
% plot(u0(i,:))
            end
        elseif strcmp(ModeType,'SShear') && DisplacementComponent == 2
            for i = 1:length(PhaseVelocity{p})
                AngularFrequency = 2*pi*ExcitationSpectrumRange{p}(1,i)*1e3; % potential approach
                k = AngularFrequency/PhaseVelocity{p}(i)*(1+1i*Attenuation{p}(i)/(2*pi));
                u2 = cos(sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2)*T);
                v2 = -1i*AngularFrequency*u2;
                sigma6 = Mu*1i*k*u2;
                PowerFlowDensity = -.5*real(sigma6.*conj(v2));
                PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                u(i,:) = real(u2(1)*exp(1i*(k*Distance-AngularFrequency*Time)))/sqrt(PowerFlow)*ExcitationSpectrumRange{p}(2,i);

                % AngularFrequency = 2*pi*ExcitationSpectrumRange{p}(1,i)*1e3; % bulk wave approach
                % k = AngularFrequency/PhaseVelocity{p}(i)*(1+1i*Attenuation{p}(i)/(2*pi));
                % k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
                % D6 = 1i*k*Mu; % sigma12
                % E = [exp(1i*k3*x3) exp(1i*k3*(Thickness-x3))];
                % v = -1i*AngularFrequency*(E(:,1)+E(:,2));
                % sigma = D6*E(:,1)+D6*E(:,2); % sigma12
                % PowerFlowDensity = -.5*real(sigma.*conj(v));
                % PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                % u(i,:) = real(exp(1i*(k*Distance-AngularFrequency*Time))+exp(1i*(k*Distance+k3*Thickness-AngularFrequency*Time)))/sqrt(PowerFlow)*ExcitationSpectrumRange{p}(2,i);
            end
        elseif strcmp(ModeType,'AShear') && DisplacementComponent == 2
            for i = 1:length(PhaseVelocity{p})
                AngularFrequency = 2*pi*ExcitationSpectrumRange{p}(1,i)*1e3; % potential approach
                k = AngularFrequency/PhaseVelocity{p}(i)*(1+1i*Attenuation{p}(i)/(2*pi));
                u2 = sin(sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2)*T);
                v2 = -1i*AngularFrequency*u2;
                sigma6 = Mu*1i*k*u2;
                PowerFlowDensity = -.5*real(sigma6.*conj(v2));
                PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                u(i,:) = real(u2(1)*exp(1i*(k*Distance-AngularFrequency*Time)))/sqrt(PowerFlow)*ExcitationSpectrumRange{p}(2,i);

                % AngularFrequency = 2*pi*ExcitationSpectrumRange{p}(1,i)*1e3; % bulk wave approach
                % k = AngularFrequency/PhaseVelocity{p}(i)*(1+1i*Attenuation{p}(i)/(2*pi));
                % k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
                % D6 = 1i*k*Mu; % sigma12
                % E = [exp(1i*k3*x3) exp(1i*k3*(Thickness-x3))];
                % v = -1i*AngularFrequency*(E(:,1)-E(:,2));
                % sigma = D6*E(:,1)-D6*E(:,2); % sigma12
                % PowerFlowDensity = -.5*real(sigma.*conj(v));
                % PowerFlow = x3(2)*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                % u(i,:) = real(exp(1i*(k*Distance-AngularFrequency*Time))-exp(1i*(k*Distance+k3*Thickness-AngularFrequency*Time)))/sqrt(PowerFlow)*ExcitationSpectrumRange{p}(2,i);
            end
        end
        if  ~(contains(ModeType,'Shear') && DisplacementComponent == 1)
            uSum{p} = sum(u);
        else
            uSum{p} = u(1,:);
        end
        Counter = Counter+1;
        if  ModeTotal > 0
            waitbar(Counter/ModeTotal,h1,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
        end
    end
end