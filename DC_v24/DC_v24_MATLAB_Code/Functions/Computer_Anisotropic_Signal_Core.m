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

% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function [uSum,Counter,ExcitationSpectrumRange] = Computer_Anisotropic_Signal_Core(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,CoherenceTime,Counter,Distance,ExcitationSpectrum,FrequencyResolution,FrequencyResolution2,h1,ModeType,ModeTotal,PhaseVelocity,Attenuation,Time,TimeLimit,z1,z2)
SamplesX3 = 50; % samples per layer

%#ok<*AGROW>
%#ok<*GVMIS>
%#ok<*MINV> 
global Stop 
Stop = 0;
I = [1 1 -1;-1 -1 1;-1 -1 1];
I1 = [1 -1;-1 1];
if  ~ToggleUpperFluid
    UpperFluid.Velocity = 1e-10;
    UpperFluid.Density = 1e-10;
end
if  ~ToggleLowerFluid
    LowerFluid.Velocity = 1e-10;
    LowerFluid.Density = 1e-10;
end
uSum{length(PhaseVelocity)} = [];
ExcitationSpectrumRange{length(PhaseVelocity)} = [];
for m = 1:SuperLayerSize
    if  strcmp(ModeType,'Coupled')
        Delta(m) = c{m}(3,3)*c{m}(4,4)*c{m}(5,5)-c{m}(3,3)*c{m}(4,5)^2;
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
    elseif strcmp(ModeType,'Lamb')
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
    end
end
% figure
% hold on
for p = 1:length(PhaseVelocity)
    if  Stop == 1
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
        catch ME
            msgbox(['IDENTIFIER: ',ME.identifier,newline,'MESSAGE: ',ME.message,newline,'FILE: ',ME.stack.file,newline])
            errordlg('The data size exceeds your system''s RAM! Decrease the parameter settings.','Error');
            if  exist('h1','var') && ModeTotal > 0
                close(h1)
            end
            return
        end
        if  strcmp(ModeType,'Coupled')
            for i = 1:length(PhaseVelocity{p})
                Layup = cell(0);
                AngularFrequency = 2*pi*ExcitationSpectrumRange{p}(1,i)*1e3;
                WaveNumber = AngularFrequency/PhaseVelocity{p}(i)*(1+1i*Attenuation{p}(i)/(2*pi));
                WaveNumber2 = WaveNumber^2;
                WaveNumber4 = WaveNumber^4;
                WaveNumber6 = WaveNumber^6;
                for m = 1:SuperLayerSize
                    rw2 = Material{m}.Density*AngularFrequency^2;
                    r2w4 = Material{m}.Density^2*AngularFrequency^4;
                    A1 = a11(m)*WaveNumber2+a12(m)*rw2;
                    A2 = a21(m)*WaveNumber4+a22(m)*rw2*WaveNumber2+a23(m)*r2w4;
                    A3 = a31(m)*WaveNumber6+a32(m)*rw2*WaveNumber4+a33(m)*r2w4*WaveNumber2+a34(m)*Material{m}.Density^3*AngularFrequency^6;
                    k3a = A2/3-A1^2/9;
                    k3b = A1^3/27-A1*A2/6+A3/2;
                    k3c = (sqrt(k3b^2+k3a^3)-k3b)^(1/3);
                    k3d = k3a/(2*k3c)-k3c/2;
                    k3e = k3a/k3c;
                    k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
                    k3(1) = sqrt(k3d-k3f-A1/3);
                    k3(2) = sqrt(k3d+k3f-A1/3);
                    k3(3) = -sqrt(k3c-k3e-A1/3);
                    k32 = k3.^2;
                    m11 = c{m}(1,1)*WaveNumber2-rw2+c{m}(5,5)*k32;
                    m12 = c{m}(1,6)*WaveNumber2+c{m}(4,5)*k32;
                    m13 = (c{m}(1,3)+c{m}(5,5))*WaveNumber*k3;
                    m22 = c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k32;
                    m23 = (c{m}(3,6)+c{m}(4,5))*WaveNumber*k3;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D3 = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber*V+c{m}(3,3)*k3.*W); % sigma33
                    D4 = 1i*(c{m}(4,5)*(k3+WaveNumber*W)+c{m}(4,4)*k3.*V); % sigma23
                    D5 = 1i*(c{m}(5,5)*(k3+WaveNumber*W)+c{m}(4,5)*k3.*V); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                    Layup{m,7} = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W]; % L2
                    Layup{m,1} = L1/Layup{m,7}; % L
                    Layup{m,3} = k3;
                    Layup{m,4} = V;
                    Layup{m,5} = W;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,8} = 1i*(c{m}(1,1)*WaveNumber+c{m}(1,6)*WaveNumber*V+c{m}(1,3)*k3.*W); % sigma11
                    Layup{m,9} = D5; % sigma13
                    Layup{m,10} = 1i*(c{m}(1,6)*WaveNumber+c{m}(6,6)*WaveNumber*V+c{m}(3,6)*k3.*W); % sigma12
                end
                Layup = repmat(Layup,Repetitions,1);
                if  SymmetricSystem
                    Layup = vertcat(Layup,flipud(Layup));
                end
                M = Layup{1};
                for m = 2:SuperLayerSize
                    N = inv(Layup{m,1}(1:3,1:3)-M(4:6,4:6));
                    M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*Layup{m,1}(1:3,4:6);Layup{m,1}(4:6,1:3)*N*M(4:6,1:3) Layup{m,1}(4:6,4:6)-Layup{m,1}(4:6,1:3)*N*Layup{m,1}(1:3,4:6)];
                end
                MM{1} = M;
                for m = 2:Repetitions
                    N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                    MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                end
                if  FluidLoading
                    G = inv(MM{end});
                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2);
                    WUpperFluid = k3UpperFluid/WaveNumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber; % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber; % in the lower fluid
                    ULowerFluid = -(WUpperFluid+G(3,1)*DUpperFluid)/(G(3,4)*DLowerFluid);
                    uInterfaces = [G(1,1)*DUpperFluid+G(1,4)*DLowerFluid*ULowerFluid;G(2,1)*DUpperFluid+G(2,4)*DLowerFluid*ULowerFluid;G(3,1)*DUpperFluid+G(3,4)*DLowerFluid*ULowerFluid];
                    if  SymmetricSystem
                        uInterfaces(:,2*Repetitions*SuperLayerSize+1) = [G(4,1)*DUpperFluid+G(4,4)*DLowerFluid*ULowerFluid;G(5,1)*DUpperFluid+G(5,4)*DLowerFluid*ULowerFluid;G(6,1)*DUpperFluid+G(6,4)*DLowerFluid*ULowerFluid];
                    else
                        uInterfaces(:,Repetitions*SuperLayerSize+1) = [G(4,1)*DUpperFluid+G(4,4)*DLowerFluid*ULowerFluid;G(5,1)*DUpperFluid+G(5,4)*DLowerFluid*ULowerFluid;G(6,1)*DUpperFluid+G(6,4)*DLowerFluid*ULowerFluid];
                    end
                else
                    Z1 = [-MM{end}(1:3,1:3)*[1 1 1;Layup{1,4};-Layup{1,5}] -MM{end}(1:3,4:6)*[1 1 1;Layup{end,4};Layup{end,5}];MM{end}(4:6,1:3)*[1 1 1;Layup{1,4};-Layup{1,5}] MM{end}(4:6,4:6)*[1 1 1;Layup{end,4};Layup{end,5}]];
                    Z2 = [MM{end}(1:3,1:3)*[1 1 1;Layup{1,4};Layup{1,5}];-MM{end}(4:6,1:3)*[1 1 1;Layup{1,4};Layup{1,5}]];
                    RT = Z1\Z2(:,1);
                    uInterfaces = [1;Layup{1,4}(1);Layup{1,5}(1)]+[1 1 1;Layup{1,4};-Layup{1,5}]*RT(1:3);
                    if  SymmetricSystem
                        uInterfaces(:,2*Repetitions*SuperLayerSize+1) = [1 1 1;Layup{end,4};Layup{end,5}]*RT(4:6);
                    else
                        uInterfaces(:,Repetitions*SuperLayerSize+1) = [1 1 1;Layup{end,4};Layup{end,5}]*RT(4:6);
                    end
                end
                Layup{1,2} = Layup{1};
                for m = 2:size(Layup,1)
                    N = inv(Layup{m,1}(1:3,1:3)-Layup{m-1,2}(4:6,4:6));
                    Layup{m,2} = [Layup{m-1,2}(1:3,1:3)+Layup{m-1,2}(1:3,4:6)*N*Layup{m-1,2}(4:6,1:3) -Layup{m-1,2}(1:3,4:6)*N*Layup{m,1}(1:3,4:6);Layup{m,1}(4:6,1:3)*N*Layup{m-1,2}(4:6,1:3) Layup{m,1}(4:6,4:6)-Layup{m,1}(4:6,1:3)*N*Layup{m,1}(1:3,4:6)];
                end
                for m = size(Layup,1):-1:2
                    N = inv(Layup{m,1}(1:3,1:3)-Layup{m-1,2}(4:6,4:6));
                    uInterfaces(:,m) = N*Layup{m-1,2}(4:6,1:3)*uInterfaces(:,1)-N*Layup{m,1}(1:3,4:6)*uInterfaces(:,m+1);
                end
                for m = 1:size(Layup,1)
                    U = Layup{m,7}\[uInterfaces(:,m);uInterfaces(:,m+1)];
                    x3 = 0:Layup{m,6}/SamplesX3:Layup{m,6};
                    if  m == 1
                        U0 = U;
                        x30 = x3;
                        x3Total = x3;
                    else
                        x3Total = horzcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}.*x3') exp(1i*Layup{m,3}.*(Layup{m,6}-x3)')];
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = -1i*AngularFrequency*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6));
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency*(Layup{m,4}(1)*U(1)*E(:,1)+Layup{m,4}(2)*U(2)*E(:,2)+Layup{m,4}(3)*U(3)*E(:,3)+Layup{m,4}(1)*U(4)*E(:,4)+Layup{m,4}(2)*U(5)*E(:,5)+Layup{m,4}(3)*U(6)*E(:,6));
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = -1i*AngularFrequency*(Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)+Layup{m,5}(3)*U(3)*E(:,3)-Layup{m,5}(1)*U(4)*E(:,4)-Layup{m,5}(2)*U(5)*E(:,5)-Layup{m,5}(3)*U(6)*E(:,6));
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,8}(1)*U(1)*E(:,1)+Layup{m,8}(2)*U(2)*E(:,2)+Layup{m,8}(3)*U(3)*E(:,3)+Layup{m,8}(1)*U(4)*E(:,4)+Layup{m,8}(2)*U(5)*E(:,5)+Layup{m,8}(3)*U(6)*E(:,6); % sigma11
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)+Layup{m,9}(3)*U(3)*E(:,3)-Layup{m,9}(1)*U(4)*E(:,4)-Layup{m,9}(2)*U(5)*E(:,5)-Layup{m,9}(3)*U(6)*E(:,6); % sigma13
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,10}(1)*U(1)*E(:,1)+Layup{m,10}(2)*U(2)*E(:,2)+Layup{m,10}(3)*U(3)*E(:,3)+Layup{m,10}(1)*U(4)*E(:,4)+Layup{m,10}(2)*U(5)*E(:,5)+Layup{m,10}(3)*U(6)*E(:,6); % sigma12
                end
                PowerFlowDensity = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                E = [exp(1i*(WaveNumber*Distance+Layup{1,3}.*x30(1:2)')) exp(1i*(WaveNumber*Distance+Layup{1,3}.*(Layup{1,6}-x30(1:2))'))];
                if  DisplacementComponent == 1 % u3
                    unorm = Layup{1,5}(1)*U0(1)*E(:,1)+Layup{1,5}(2)*U0(2)*E(:,2)+Layup{1,5}(3)*U0(3)*E(:,3)-Layup{1,5}(1)*U0(4)*E(:,4)-Layup{1,5}(2)*U0(5)*E(:,5)-Layup{1,5}(3)*U0(6)*E(:,6);
                    if  SuperLayerSize > 1 || (SuperLayerSize == 1 && ~strcmp(Material{1}.Class,'Cubic'))
                        u(i,:) = real(...
                             abs(Layup{1,5}(1)*U0(1))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(Layup{1,5}(2)*U0(2))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(Layup{1,5}(3)*U0(3))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(Layup{1,5}(1)*U0(4))*exp(1i*(WaveNumber*Distance+Layup{1,3}(1)*LayerThicknesses(1)-AngularFrequency*Time'))...
                            +abs(Layup{1,5}(2)*U0(5))*exp(1i*(WaveNumber*Distance+Layup{1,3}(2)*LayerThicknesses(1)-AngularFrequency*Time'))...
                            +abs(Layup{1,5}(3)*U0(6))*exp(1i*(WaveNumber*Distance+Layup{1,3}(3)*LayerThicknesses(1)-AngularFrequency*Time')))';
                    end
                elseif DisplacementComponent == 2 % u1
                    unorm = U0(1)*E(:,1)+U0(2)*E(:,2)+U0(3)*E(:,3)+U0(4)*E(:,4)+U0(5)*E(:,5)+U0(6)*E(:,6);
                    if  SuperLayerSize > 1 || (SuperLayerSize == 1 && ~strcmp(Material{1}.Class,'Cubic'))
                        u(i,:) = real(...
                             abs(U0(1))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(U0(2))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(U0(3))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(U0(4))*exp(1i*(WaveNumber*Distance+Layup{1,3}(1)*LayerThicknesses(1)-AngularFrequency*Time'))...
                            +abs(U0(5))*exp(1i*(WaveNumber*Distance+Layup{1,3}(2)*LayerThicknesses(1)-AngularFrequency*Time'))...
                            +abs(U0(6))*exp(1i*(WaveNumber*Distance+Layup{1,3}(3)*LayerThicknesses(1)-AngularFrequency*Time')))';
                    end
                end
                uNorm(i) = abs(real(unorm(1)*exp(-1i*angle(unorm(2))))/sqrt(PowerFlow));
                if  SuperLayerSize == 1 && strcmp(Material{1}.Class,'Cubic')
                    ki(i) = WaveNumber;
                    AngularFrequencyi(i) = AngularFrequency;
                    UNorm(:,i) = abs(U0(2:6)/U0(1));
                    if  DisplacementComponent == 1
                        WNorm(:,i) = abs(Layup{1,5}(2:3)/Layup{1,5}(1));
                    end
                end
            end
            uNorm(isnan(uNorm)) = 0;
            uNormSmooth = filloutliers(uNorm,'spline','movmedian',5,'ThresholdFactor',1);
            if  SuperLayerSize == 1 && strcmp(Material{1}.Class,'Cubic')
                UNorm(:,isnan(UNorm(1,:))) = 0;
                UNormSmooth(1,1:length(UNorm)) = filloutliers(UNorm(1,:),'spline','movmedian',10,'ThresholdFactor',1);
                UNormSmooth(2,1:length(UNorm)) = filloutliers(UNorm(2,:),'spline','movmedian',10,'ThresholdFactor',1);
                UNormSmooth(3,1:length(UNorm)) = filloutliers(UNorm(3,:),'spline','movmedian',10,'ThresholdFactor',1);
                UNormSmooth(4,1:length(UNorm)) = filloutliers(UNorm(4,:),'spline','movmedian',10,'ThresholdFactor',1);
                UNormSmooth(5,1:length(UNorm)) = filloutliers(UNorm(5,:),'spline','movmedian',10,'ThresholdFactor',1);
                if  DisplacementComponent == 1
                    WNormSmooth(1,:) = filloutliers(WNorm(1,:),'spline','movmedian',10,'ThresholdFactor',1);
                    WNormSmooth(2,:) = filloutliers(WNorm(2,:),'spline','movmedian',10,'ThresholdFactor',1);
                end
            end
% figure('name','uNorm')
% hold on
% plot(uNorm,'linewidth',4,'color','r')
% plot(uNormSmooth,'linewidth',1.5,'color','g')
            for i = 1:length(PhaseVelocity{p})
                if  SuperLayerSize == 1 && strcmp(Material{1}.Class,'Cubic')
                    if  DisplacementComponent == 1 % u3
                        u(i,:) = real(...
                                                               exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +WNormSmooth(1,i)*UNormSmooth(1,i)*exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +WNormSmooth(2,i)*UNormSmooth(2,i)*exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +                 UNormSmooth(3,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(1)*LayerThicknesses(1)-AngularFrequencyi(i)*Time'))...
                            +WNormSmooth(1,i)*UNormSmooth(4,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(2)*LayerThicknesses(1)-AngularFrequencyi(i)*Time'))...
                            +WNormSmooth(2,i)*UNormSmooth(5,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(3)*LayerThicknesses(1)-AngularFrequencyi(i)*Time')))';
                    elseif DisplacementComponent == 2 % u1
                        u(i,:) = real(...
                                              exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +UNormSmooth(1,i)*exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +UNormSmooth(2,i)*exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +UNormSmooth(3,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(1)*LayerThicknesses(1)-AngularFrequencyi(i)*Time'))...
                            +UNormSmooth(4,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(2)*LayerThicknesses(1)-AngularFrequencyi(i)*Time'))...
                            +UNormSmooth(5,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(3)*LayerThicknesses(1)-AngularFrequencyi(i)*Time')))';
                    end
                end
                u(i,:) = u(i,:)*uNormSmooth(i)/max(u(i,:))*ExcitationSpectrumRange{p}(2,i);
% plot(u(i,1:120))
            end       
        elseif strcmp(ModeType,'Lamb')
            for i = 1:length(PhaseVelocity{p})
                Layup = cell(0);
                AngularFrequency = 2*pi*ExcitationSpectrumRange{p}(1,i)*1e3;
                WaveNumber = AngularFrequency/PhaseVelocity{p}(i)*(1+1i*Attenuation{p}(i)/(2*pi));
                WaveNumber2 = WaveNumber^2;
                WaveNumber4 = WaveNumber^4;
                for m = 1:SuperLayerSize
                    rw2 = Material{m}.Density*AngularFrequency^2;
                    A2 = a21(m)*WaveNumber2+a22(m)*rw2;
                    A3 = a31(m)*WaveNumber4+a32(m)*rw2*WaveNumber2+Material{m}.Density^2*AngularFrequency^4;
                    k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    W = (rw2-c{m}(1,1)*WaveNumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*WaveNumber*k3);
                    D3 = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,3)*k3.*W); % sigma33
                    D5 = 1i*(c{m}(5,5)*(k3+WaveNumber*W)); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
                    L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                    Layup{m,7} = [1 1 E;W -W.*E;E 1 1;W.*E -W]; % L2
                    Layup{m,1} = L1/Layup{m,7}; % L
                    Layup{m,3} = k3;
                    Layup{m,5} = W;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,8} = 1i*(c{m}(1,1)*WaveNumber+c{m}(1,3)*k3.*W); % sigma11
                    Layup{m,9} = D5; % sigma13
                end
                Layup = repmat(Layup,Repetitions,1);
                if  SymmetricSystem
                    Layup = vertcat(Layup,flipud(Layup));
                end
                M = Layup{1};
                for m = 2:SuperLayerSize
                    N = inv(Layup{m,1}(1:2,1:2)-M(3:4,3:4));
                    M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*Layup{m,1}(1:2,3:4);Layup{m,1}(3:4,1:2)*N*M(3:4,1:2) Layup{m,1}(3:4,3:4)-Layup{m,1}(3:4,1:2)*N*Layup{m,1}(1:2,3:4)];
                end
                MM{1} = M;
                for m = 2:Repetitions
                    N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                    MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                end
                if  FluidLoading
                    G = inv(MM{end});
                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2);
                    WUpperFluid = k3UpperFluid/WaveNumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber; % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber; % in the lower fluid
                    ULowerFluid = -(WUpperFluid+G(2,1)*DUpperFluid)/(G(2,3)*DLowerFluid);
                    uInterfaces = [G(1,1)*DUpperFluid+G(1,3)*DLowerFluid*ULowerFluid;G(2,1)*DUpperFluid+G(2,3)*DLowerFluid*ULowerFluid];
                    if  SymmetricSystem
                        uInterfaces(:,2*Repetitions*SuperLayerSize+1) = [G(3,1)*DUpperFluid+G(3,3)*DLowerFluid*ULowerFluid;G(4,1)*DUpperFluid+G(4,3)*DLowerFluid*ULowerFluid];
                    else
                        uInterfaces(:,Repetitions*SuperLayerSize+1) = [G(3,1)*DUpperFluid+G(3,3)*DLowerFluid*ULowerFluid;G(4,1)*DUpperFluid+G(4,3)*DLowerFluid*ULowerFluid];
                    end
                else
                    Z1 = [-MM{end}(1:2,1:2)*[1 1;-Layup{1,5}] -MM{end}(1:2,3:4)*[1 1;Layup{end,5}];MM{end}(3:4,1:2)*[1 1;-Layup{1,5}] MM{end}(3:4,3:4)*[1 1;Layup{end,5}]];
                    Z2 = [MM{end}(1:2,1:2)*[1 1;Layup{1,5}];-MM{end}(3:4,1:2)*[1 1;Layup{1,5}]];
                    RT = Z1\Z2(:,1);
                    uInterfaces = [1;Layup{1,5}(1)]+[1 1;-Layup{1,5}]*RT(1:2);
                    if  SymmetricSystem
                        uInterfaces(:,2*Repetitions*SuperLayerSize+1) = [1 1;Layup{end,5}]*RT(3:4);
                    else
                        uInterfaces(:,Repetitions*SuperLayerSize+1) = [1 1;Layup{end,5}]*RT(3:4);
                    end
                end
                Layup{1,2} = Layup{1};
                for m = 2:size(Layup,1)
                    N = inv(Layup{m,1}(1:2,1:2)-Layup{m-1,2}(3:4,3:4));
                    Layup{m,2} = [Layup{m-1,2}(1:2,1:2)+Layup{m-1,2}(1:2,3:4)*N*Layup{m-1,2}(3:4,1:2) -Layup{m-1,2}(1:2,3:4)*N*Layup{m,1}(1:2,3:4);Layup{m,1}(3:4,1:2)*N*Layup{m-1,2}(3:4,1:2) Layup{m,1}(3:4,3:4)-Layup{m,1}(3:4,1:2)*N*Layup{m,1}(1:2,3:4)];
                end
                for m = size(Layup,1):-1:2
                    N = inv(Layup{m,1}(1:2,1:2)-Layup{m-1,2}(3:4,3:4));
                    uInterfaces(:,m) = N*Layup{m-1,2}(3:4,1:2)*uInterfaces(:,1)-N*Layup{m,1}(1:2,3:4)*uInterfaces(:,m+1);
                end
                for m = 1:size(Layup,1)
                    U = Layup{m,7}\[uInterfaces(:,m);uInterfaces(:,m+1)];
                    x3 = 0:Layup{m,6}/SamplesX3:Layup{m,6};
                    if  m == 1
                        U0 = U;
                        x30 = x3;
                        x3Total = x3;
                    else
                        x3Total = horzcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}.*x3') exp(1i*Layup{m,3}.*(Layup{m,6}-x3)')];
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = -1i*AngularFrequency*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4));
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = -1i*AngularFrequency*(Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)-Layup{m,5}(1)*U(3)*E(:,3)-Layup{m,5}(2)*U(4)*E(:,4));
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,8}(1)*U(1)*E(:,1)+Layup{m,8}(2)*U(2)*E(:,2)+Layup{m,8}(1)*U(3)*E(:,3)+Layup{m,8}(2)*U(4)*E(:,4); % sigma11
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)-Layup{m,9}(1)*U(3)*E(:,3)-Layup{m,9}(2)*U(4)*E(:,4); % sigma13
                end
                PowerFlowDensity = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                E = [exp(1i*(WaveNumber*Distance+Layup{1,3}.*x30(1:2)')) exp(1i*(WaveNumber*Distance+Layup{1,3}.*(Layup{1,6}-x30(1:2))'))];
                if  DisplacementComponent == 1 % u3
                    unorm = Layup{1,5}(1)*U0(1)*E(:,1)+Layup{1,5}(2)*U0(2)*E(:,2)-Layup{1,5}(1)*U0(3)*E(:,3)-Layup{1,5}(2)*U0(4)*E(:,4);
                    if  SuperLayerSize > 1
                        u(i,:) = real(...
                             abs(Layup{1,5}(1)*U0(1))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(Layup{1,5}(2)*U0(2))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(Layup{1,5}(1)*U0(3))*exp(1i*(WaveNumber*Distance+Layup{1,3}(1)*LayerThicknesses(1)-AngularFrequency*Time'))...
                            +abs(Layup{1,5}(2)*U0(4))*exp(1i*(WaveNumber*Distance+Layup{1,3}(2)*LayerThicknesses(1)-AngularFrequency*Time')))';
                    end
                elseif DisplacementComponent == 2 % u1
                    unorm = U0(1)*E(:,1)+U0(2)*E(:,2)+U0(3)*E(:,3)+U0(4)*E(:,4);
                    if  SuperLayerSize > 1
                        u(i,:) = real(...
                             abs(U0(1))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(U0(2))*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                            +abs(U0(3))*exp(1i*(WaveNumber*Distance+Layup{1,3}(1)*LayerThicknesses(1)-AngularFrequency*Time'))...
                            +abs(U0(4))*exp(1i*(WaveNumber*Distance+Layup{1,3}(2)*LayerThicknesses(1)-AngularFrequency*Time')))';
                    end
                end
                uNorm(i) = abs(real(unorm(1)*exp(-1i*angle(unorm(2))))/sqrt(PowerFlow));
                if  SuperLayerSize == 1
                    ki(i) = WaveNumber;
                    AngularFrequencyi(i) = AngularFrequency;
                    UNorm(:,i) = abs(U0(2:4)/U0(1));
                    if  DisplacementComponent == 1
                        WNorm(i) = abs(Layup{1,5}(2)/Layup{1,5}(1));
                    end
                end
            end
            uNorm(isnan(uNorm)) = 0;
            uNormSmooth = filloutliers(uNorm,'spline','movmedian',5,'ThresholdFactor',1);
            if  SuperLayerSize == 1
                UNorm(:,isnan(UNorm(1,:))) = 0;
                UNormSmooth(1,1:length(UNorm)) = filloutliers(UNorm(1,:),'spline','movmedian',10,'ThresholdFactor',1);
                UNormSmooth(2,1:length(UNorm)) = filloutliers(UNorm(2,:),'spline','movmedian',10,'ThresholdFactor',1);
                UNormSmooth(3,1:length(UNorm)) = filloutliers(UNorm(3,:),'spline','movmedian',10,'ThresholdFactor',1);
                if  DisplacementComponent == 1
                    WNormSmooth = filloutliers(WNorm,'spline','movmedian',10,'ThresholdFactor',1);
                end
            end
% figure('name','uNorm')
% hold on
% plot(uNorm,'linewidth',4,'color','r')
% plot(uNormSmooth,'linewidth',1.5,'color','g')
            for i = 1:length(PhaseVelocity{p})
                if  SuperLayerSize == 1
                    if  DisplacementComponent == 1 % u3
                        u(i,:) = real(...
                                                             exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +WNormSmooth(i)*UNormSmooth(1,i)*exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +               UNormSmooth(2,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(1)*LayerThicknesses(1)-AngularFrequencyi(i)*Time'))...
                            +WNormSmooth(i)*UNormSmooth(3,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(2)*LayerThicknesses(1)-AngularFrequencyi(i)*Time')))';
                    elseif DisplacementComponent == 2 % u1
                        u(i,:) = real(...
                                              exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +UNormSmooth(1,i)*exp(1i*(ki(i)*Distance-AngularFrequencyi(i)*Time'))...
                            +UNormSmooth(2,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(1)*LayerThicknesses(1)-AngularFrequencyi(i)*Time'))...
                            +UNormSmooth(3,i)*exp(1i*(ki(i)*Distance+Layup{1,3}(2)*LayerThicknesses(1)-AngularFrequencyi(i)*Time')))';
                    end
                end
                u(i,:) = u(i,:)*uNormSmooth(i)/max(u(i,:))*ExcitationSpectrumRange{p}(2,i);
% plot(u(i,:))
            end
        elseif strcmp(ModeType,'Shear') && DisplacementComponent == 2
            for i = 1:length(PhaseVelocity{p})
                Layup = cell(0);
                AngularFrequency = 2*pi*ExcitationSpectrumRange{p}(1,i)*1e3;
                WaveNumber = AngularFrequency/PhaseVelocity{p}(i)*(1+1i*Attenuation{p}(i)/(2*pi));
                for m = 1:SuperLayerSize
                    k3 = sqrt((Material{m}.Density*AngularFrequency^2-WaveNumber^2*c{m}(6,6))/c{m}(4,4));
                    if  k3 == 0
                        k3 = 1e-10;
                    end
                    E = exp(1i*k3*LayerThicknesses(m));
                    Layup{m,1} = 1i*k3*c{m}(4,4)/(E^2-1)*[-1-E^2 2*E;-2*E 1+E^2];
                    Layup{m,3} = k3;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,7} = [1 E;E 1];
                    Layup{m,8} = 1i*WaveNumber*c{m}(6,6); % sigma12
                end
                Layup = repmat(Layup,Repetitions,1);
                if  SymmetricSystem
                    Layup = vertcat(Layup,flipud(Layup));
                end
                M = Layup{1};
                for m = 2:SuperLayerSize
                    N = inv(Layup{m,1}(1,1)-M(2,2));
                    M = [M(1,1)+M(1,2)*N*M(2,1) -M(1,2)*N*Layup{m,1}(1,2);Layup{m,1}(2,1)*N*M(2,1) Layup{m,1}(2,2)-Layup{m,1}(2,1)*N*Layup{m,1}(1,2)];
                end
                MM{1} = M;
                for m = 2:Repetitions
                    N = inv(MM{1}(1,1)-MM{m-1}(2,2));
                    MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{1}(1,2);MM{1}(2,1)*N*MM{m-1}(2,1) MM{1}(2,2)-MM{1}(2,1)*N*MM{1}(1,2)];
                end
                if  SymmetricSystem
                    N = inv(-2*MM{end}(2,2));
                    MM{end} = [MM{end}(1,1)+MM{end}(1,2)*N*MM{end}(2,1) MM{end}(1,2)*N*MM{end}(2,1);-MM{end}(1,2)*N*MM{end}(2,1) -MM{end}(1,1)-MM{end}(1,2)*N*MM{end}(2,1)];
                end
                Z1 = [MM{end}(1,1) -MM{end}(1,2);MM{end}(2,1) MM{end}(2,2)]; % changed sign in Z1(1,1)!
                Z2 = [MM{end}(1,1);-MM{end}(2,1)];
                RT = Z1\Z2;
                uInterfaces = 1+RT(1);
                if  SymmetricSystem
                    uInterfaces(:,2*Repetitions*SuperLayerSize+1) = RT(2);
                else
                    uInterfaces(:,Repetitions*SuperLayerSize+1) = RT(2);
                end
                Layup{1,2} = Layup{1};
                for m = 2:size(Layup,1)
                    N = inv(Layup{m,1}(1,1)-Layup{m-1,2}(2,2));
                    Layup{m,2} = [Layup{m-1,2}(1,1)+Layup{m-1,2}(1,2)*N*Layup{m-1,2}(2,1) -Layup{m-1,2}(1,2)*N*Layup{m,1}(1,2);Layup{m,1}(2,1)*N*Layup{m-1,2}(2,1) Layup{m,1}(2,2)-Layup{m,1}(2,1)*N*Layup{m,1}(1,2)];
                end
                for m = size(Layup,1):-1:2
                    N = inv(Layup{m,1}(1,1)-Layup{m-1,2}(2,2));
                    uInterfaces(:,m) = N*Layup{m-1,2}(2,1)*uInterfaces(:,1)-N*Layup{m,1}(1,2)*uInterfaces(:,m+1);
                end
                for m = 1:size(Layup,1)
                    U = Layup{m,7}\[uInterfaces(:,m);uInterfaces(:,m+1)];
                    x3 = 0:Layup{m,6}/SamplesX3:Layup{m,6};
                    if  m == 1
                        U0 = U;
                        x30 = x3;
                        x3Total = x3;
                    else
                        x3Total = horzcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}*x3') exp(1i*Layup{m,3}*(Layup{m,6}-x3)')];
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3)) = -1i*AngularFrequency*(U(1)*E(:,1)+U(2)*E(:,2));
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3)) = Layup{m,8}(1)*U(1)*E(:,1)+Layup{m,8}(1)*U(2)*E(:,2); % sigma12
                end
                PowerFlowDensity = -.5*real(sigma.*conj(v));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                E = [exp(1i*(WaveNumber*Distance+Layup{1,3}.*x30(1:2)')) exp(1i*(WaveNumber*Distance+Layup{1,3}.*(Layup{1,6}-x30(1:2))'))];
                unorm = U0(1)*E(:,1)+U0(2)*E(:,2);
                uNorm(i) = abs(real(unorm(1)*exp(-1i*angle(unorm(2))))/sqrt(PowerFlow));
                u(i,:) = real(...
                     U0(1)*exp(1i*(WaveNumber*Distance-AngularFrequency*Time'))...
                    +U0(2)*exp(1i*(WaveNumber*Distance+k3*LayerThicknesses(1)-AngularFrequency*Time')))';
            end
            uNormSmooth = filloutliers(uNorm,'spline','movmedian',5,'ThresholdFactor',1);
% figure('name','uNorm')
% hold on
% plot(uNorm,'linewidth',4,'color','r')
% plot(uNormSmooth,'linewidth',1.5,'color','g')
            for i = 1:length(PhaseVelocity{p})
                u(i,:) = u(i,:)*uNormSmooth(i)/max(u(i,:))*ExcitationSpectrumRange{p}(2,i);
% plot(u(i,:))
            end
        end
        if  ~(strcmp(ModeType,'Shear') && DisplacementComponent == 1)
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