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
function [X,Counter] = Computer_Anisotropic_EnergyVelocity_Core(X,ModeType,ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta)
Layers = Repetitions*length(LayerThicknesses);
if  SymmetricSystem
    Layers = 2*Layers;
end
SamplesX3 = ceil(50/Layers); % samples per layer

%#ok<*AGROW>
%#ok<*MINV>
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
for m = 1:SuperLayerSize
    if  contains(ModeType,'Coupled')
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
    elseif strcmp(ModeType,'Lamb') || strcmp(ModeType,'Scholte')
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
    end
end
for p = 1:length(X)
    if  contains(ModeType,'Coupled')
        for q = 1:height(X{p})
            Layup = cell(0);
            AngularFrequency = 2*pi*X{p}(q,1)*1e3;
            Attenuation = X{p}(q,7)*X{p}(q,4)/X{p}(q,1); % Np/wavelength
            WaveNumber = AngularFrequency/(X{p}(q,4)*1e3)*(1+1i*Attenuation/(2*pi));
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
                Layup{m,8} = Material{m}.Density;
                Layup{m,9} = 1i*(c{m}(1,1)*WaveNumber+c{m}(1,6)*WaveNumber*V+c{m}(1,3)*k3.*W); % sigma11
                Layup{m,10} = 1i*(c{m}(1,2)*WaveNumber+c{m}(2,6)*WaveNumber*V+c{m}(2,3)*k3.*W); % sigma22
                Layup{m,11} = D3; % sigma33
                Layup{m,12} = D4; % sigma23
                Layup{m,13} = D5; % sigma13
                Layup{m,14} = 1i*(c{m}(1,6)*WaveNumber+c{m}(6,6)*WaveNumber*V+c{m}(3,6)*k3.*W); % sigma12
                Layup{m,15} = 1i*WaveNumber; % epsilon11
                Layup{m,16} = 1i*k3.*W; % epsilon33
                Layup{m,17} = 1i*k3.*V; % epsilon23
                Layup{m,18} = 1i*(k3+WaveNumber*W); % epsilon13
                Layup{m,19} = 1i*WaveNumber*V; % epsilon12
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
                if  Viscoelastic && strcmp(ModeType,'ScholteCoupled')
                    k3UpperFluid = -k3UpperFluid;
                end
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
                    x3Total = x3;
                else
                    x3Total = horzcat(x3Total,x3+x3Total(end));
                end
                E = [exp(1i*Layup{m,3}.*x3') exp(1i*Layup{m,3}.*(Layup{m,6}-x3)')];
                v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = -1i*AngularFrequency*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6));
                v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency*(Layup{m,4}(1)*U(1)*E(:,1)+Layup{m,4}(2)*U(2)*E(:,2)+Layup{m,4}(3)*U(3)*E(:,3)+Layup{m,4}(1)*U(4)*E(:,4)+Layup{m,4}(2)*U(5)*E(:,5)+Layup{m,4}(3)*U(6)*E(:,6));
                v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = -1i*AngularFrequency*(Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)+Layup{m,5}(3)*U(3)*E(:,3)-Layup{m,5}(1)*U(4)*E(:,4)-Layup{m,5}(2)*U(5)*E(:,5)-Layup{m,5}(3)*U(6)*E(:,6));
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)+Layup{m,9}(3)*U(3)*E(:,3)+Layup{m,9}(1)*U(4)*E(:,4)+Layup{m,9}(2)*U(5)*E(:,5)+Layup{m,9}(3)*U(6)*E(:,6); % sigma11
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = Layup{m,10}(1)*U(1)*E(:,1)+Layup{m,10}(2)*U(2)*E(:,2)+Layup{m,10}(3)*U(3)*E(:,3)+Layup{m,10}(1)*U(4)*E(:,4)+Layup{m,10}(2)*U(5)*E(:,5)+Layup{m,10}(3)*U(6)*E(:,6); % sigma22
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(2)*U(2)*E(:,2)+Layup{m,11}(3)*U(3)*E(:,3)+Layup{m,11}(1)*U(4)*E(:,4)+Layup{m,11}(2)*U(5)*E(:,5)+Layup{m,11}(3)*U(6)*E(:,6); % sigma33
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,12}(1)*U(1)*E(:,1)+Layup{m,12}(2)*U(2)*E(:,2)+Layup{m,12}(3)*U(3)*E(:,3)-Layup{m,12}(1)*U(4)*E(:,4)-Layup{m,12}(2)*U(5)*E(:,5)-Layup{m,12}(3)*U(6)*E(:,6); % sigma23
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,13}(1)*U(1)*E(:,1)+Layup{m,13}(2)*U(2)*E(:,2)+Layup{m,13}(3)*U(3)*E(:,3)-Layup{m,13}(1)*U(4)*E(:,4)-Layup{m,13}(2)*U(5)*E(:,5)-Layup{m,13}(3)*U(6)*E(:,6); % sigma13
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,14}(1)*U(1)*E(:,1)+Layup{m,14}(2)*U(2)*E(:,2)+Layup{m,14}(3)*U(3)*E(:,3)+Layup{m,14}(1)*U(4)*E(:,4)+Layup{m,14}(2)*U(5)*E(:,5)+Layup{m,14}(3)*U(6)*E(:,6); % sigma12
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,15}*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6)); % epsilon11
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,16}(1)*U(1)*E(:,1)+Layup{m,16}(2)*U(2)*E(:,2)+Layup{m,16}(3)*U(3)*E(:,3)+Layup{m,16}(1)*U(4)*E(:,4)+Layup{m,16}(2)*U(5)*E(:,5)+Layup{m,16}(3)*U(6)*E(:,6); % epsilon33
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,17}(1)*U(1)*E(:,1)+Layup{m,17}(2)*U(2)*E(:,2)+Layup{m,17}(3)*U(3)*E(:,3)-Layup{m,17}(1)*U(4)*E(:,4)-Layup{m,17}(2)*U(5)*E(:,5)-Layup{m,17}(3)*U(6)*E(:,6); % epsilon23
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,18}(1)*U(1)*E(:,1)+Layup{m,18}(2)*U(2)*E(:,2)+Layup{m,18}(3)*U(3)*E(:,3)-Layup{m,18}(1)*U(4)*E(:,4)-Layup{m,18}(2)*U(5)*E(:,5)-Layup{m,18}(3)*U(6)*E(:,6); % epsilon13
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,19}(1)*U(1)*E(:,1)+Layup{m,19}(2)*U(2)*E(:,2)+Layup{m,19}(3)*U(3)*E(:,3)+Layup{m,19}(1)*U(4)*E(:,4)+Layup{m,19}(2)*U(5)*E(:,5)+Layup{m,19}(3)*U(6)*E(:,6); % epsilon12  
                KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,8}*(abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1)).^2+abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2)).^2+abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3)).^2);
            end
            StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+real(epsilon(:,1)).*real(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,3)).*real(sigma(:,3))+imag(epsilon(:,4)).*imag(sigma(:,4))+real(epsilon(:,4)).*real(sigma(:,4))+imag(epsilon(:,5)).*imag(sigma(:,5))+real(epsilon(:,5)).*real(sigma(:,5))+imag(epsilon(:,6)).*imag(sigma(:,6))+real(epsilon(:,6)).*real(sigma(:,6)));
            PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)));
            PowerFlowDensity(:,2) = -.5*real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)));
            PowerFlow = [trapz(x3Total,PowerFlowDensity(:,1)) trapz(x3Total,PowerFlowDensity(:,2))];
            TotalEnergy = .5*(trapz(x3Total,StrainEnergyDensity)+trapz(x3Total,KineticEnergyDensity));
            X{p}(q,5:6) = PowerFlow/TotalEnergy/1e3; % ce1,ce2
        end
    elseif strcmp(ModeType,'Lamb') || strcmp(ModeType,'Scholte')
        for q = 1:height(X{p})
            Layup = cell(0);
            AngularFrequency = 2*pi*X{p}(q,1)*1e3;
            Attenuation = X{p}(q,6)*X{p}(q,4)/X{p}(q,1); % Np/wavelength
            WaveNumber = AngularFrequency/(X{p}(q,4)*1e3)*(1+1i*Attenuation/(2*pi));
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
                Layup{m,8} = Material{m}.Density;
                Layup{m,9} = 1i*(c{m}(1,1)*WaveNumber+c{m}(1,3)*k3.*W); % sigma11
                Layup{m,11} = D3; % sigma33
                Layup{m,12} = D5; % sigma13
                Layup{m,13} = 1i*WaveNumber; % epsilon11
                Layup{m,14} = 1i*k3.*W; % epsilon33
                Layup{m,15} = 1i*(k3+WaveNumber*W); % epsilon13
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
                if  Viscoelastic && strcmp(ModeType,'Scholte')
                    k3UpperFluid = -k3UpperFluid;
                end
                WUpperFluid = k3UpperFluid/WaveNumber;
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber;
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber;
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
                    x3Total = x3;
                else
                    x3Total = horzcat(x3Total,x3+x3Total(end));
                end
                E = [exp(1i*Layup{m,3}.*x3') exp(1i*Layup{m,3}.*(Layup{m,6}-x3)')];
                v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = -1i*AngularFrequency*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4));
                v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = -1i*AngularFrequency*(Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)-Layup{m,5}(1)*U(3)*E(:,3)-Layup{m,5}(2)*U(4)*E(:,4));
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)+Layup{m,9}(1)*U(3)*E(:,3)+Layup{m,9}(2)*U(4)*E(:,4); % sigma11
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(2)*U(2)*E(:,2)+Layup{m,11}(1)*U(3)*E(:,3)+Layup{m,11}(2)*U(4)*E(:,4); % sigma33
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,12}(1)*U(1)*E(:,1)+Layup{m,12}(2)*U(2)*E(:,2)-Layup{m,12}(1)*U(3)*E(:,3)-Layup{m,12}(2)*U(4)*E(:,4); % sigma13
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,13}*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)); % epsilon11
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,14}(1)*U(1)*E(:,1)+Layup{m,14}(2)*U(2)*E(:,2)+Layup{m,14}(1)*U(3)*E(:,3)+Layup{m,14}(2)*U(4)*E(:,4); % epsilon33
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,15}(1)*U(1)*E(:,1)+Layup{m,15}(2)*U(2)*E(:,2)-Layup{m,15}(1)*U(3)*E(:,3)-Layup{m,15}(2)*U(4)*E(:,4); % epsilon13
                KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,8}*(abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1)).^2+abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3)).^2);
            end
            StrainEnergyDensity = .5*(imag(epsilon(:,1)).*imag(sigma(:,1))+real(epsilon(:,1)).*real(sigma(:,1))+imag(epsilon(:,3)).*imag(sigma(:,3))+real(epsilon(:,3)).*real(sigma(:,3))+imag(epsilon(:,5)).*imag(sigma(:,5))+real(epsilon(:,5)).*real(sigma(:,5)));
            PowerFlowDensity = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
            PowerFlow = trapz(x3Total,PowerFlowDensity);
            TotalEnergy = .5*(trapz(x3Total,StrainEnergyDensity)+trapz(x3Total,KineticEnergyDensity));
            X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % ce1
        end
    elseif strcmp(ModeType,'Shear')
        for q = 1:height(X{p})
            Layup = cell(0);
            PhaseVelocity = X{p}(q,4)*1e3;
            Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(q,1)*1e3;
            WaveNumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
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
                Layup{m,8} = 1i*k3*c{m}(4,4); % sigma23
                Layup{m,9} = 1i*WaveNumber*c{m}(6,6); % sigma12
                Layup{m,10} = 1i*k3; % epsilon23
                Layup{m,11} = 1i*WaveNumber; % epsilon12
                Layup{m,12} = Material{m}.Density;
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
                    x3Total = x3;
                else
                    x3Total = horzcat(x3Total,x3+x3Total(end));
                end
                E = [exp(1i*Layup{m,3}*x3') exp(1i*Layup{m,3}*(Layup{m,6}-x3)')];
                v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency*(U(1)*E(:,1)+U(2)*E(:,2));
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,8}(1)*U(1)*E(:,1)-Layup{m,8}(1)*U(2)*E(:,2); % sigma23
                sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(1)*U(2)*E(:,2); % sigma12
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,10}(1)*U(1)*E(:,1)-Layup{m,10}(1)*U(2)*E(:,2); % epsilon23
                epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(1)*U(2)*E(:,2); % epsilon12
                KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,12}*abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2)).^2;
            end
            StrainEnergyDensity = .5*(imag(epsilon(:,4)).*imag(sigma(:,4))+real(epsilon(:,4)).*real(sigma(:,4))+imag(epsilon(:,6)).*imag(sigma(:,6))+real(epsilon(:,6)).*real(sigma(:,6)));
            PowerFlowDensity = -.5*real(sigma(:,6).*conj(v(:,2)));
            PowerFlow = trapz(x3Total,PowerFlowDensity);
            TotalEnergy = .5*(trapz(x3Total,StrainEnergyDensity)+trapz(x3Total,KineticEnergyDensity));
            X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % ce1
        end
    end
    X{p}(:,5) = filloutliers(X{p}(:,5),'spline','movmedian',5,'ThresholdFactor',1);
    if  contains(ModeType,'Coupled')
        X{p}(:,6) = filloutliers(X{p}(:,6),'spline','movmedian',5,'ThresholdFactor',1);
    end
    Counter = Counter+1;
    if  ModeTotal > 0
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
    end
end