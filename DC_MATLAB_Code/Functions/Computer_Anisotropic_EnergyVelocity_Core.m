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
function [X,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,X,ModeType,ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3)
Layers = Repetitions*length(LayerThicknesses);
if  SymmetricSystem
    Layers = 2*Layers;
end
SamplesX3 = ceil(SamplesX3/Layers); % samples per layer

%#ok<*AGROW>
%#ok<*MINV>
%#ok<*PFTUSW>
%#ok<*PFBNS>
%#ok<*GVMIS>
global Stop 
Stop = 0;
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
        if  Multithreading
            EnergyVelocity = [0 0];
            PhaseVelocity = X{p}(:,4)*1e3;
            Attenuation = X{p}(:,7).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                v = 0;
                sigma = 0;
                epsilon = 0;
                KineticEnergyDensity = 0;
                PowerFlowDensity = 0;
                x3Total = 0;
                Layup = cell(0);
                AngularFrequency2 = AngularFrequency(q)^2;
                Wavenumber = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                Wavenumber2 = Wavenumber^2;
                Wavenumber4 = Wavenumber2^2;
                Wavenumber6 = Wavenumber2^3;
                for m = 1:SuperLayerSize
                    rw2 = Material{m}.Density*AngularFrequency2;
                    r2w4 = rw2^2;
                    A1 = a11(m)*Wavenumber2+a12(m)*rw2;
                    A2 = a21(m)*Wavenumber4+a22(m)*rw2*Wavenumber2+a23(m)*r2w4;
                    A3 = a31(m)*Wavenumber6+a32(m)*rw2*Wavenumber4+a33(m)*r2w4*Wavenumber2+a34(m)*rw2^3;
                    d1 = A2/3-A1^2/9;
                    d2 = A1^3/27-A1*A2/6+A3/2;
                    d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                    d4 = d1/(2*d3)-d3/2;
                    d5 = d1/d3;
                    d6 = (sqrt(3)*(d3+d5)*1i)/2;
                    k3 = [sqrt(d4-d6-A1/3) sqrt(d4+d6-A1/3) -sqrt(d3-d5-A1/3)];
                    k32 = k3.^2;
                    m11 = c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k32;
                    m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
                    m13 = (c{m}(1,3)+c{m}(5,5))*Wavenumber*k3;
                    m22 = c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k32;
                    m23 = (c{m}(3,6)+c{m}(4,5))*Wavenumber*k3;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber*V+c{m}(3,3)*k3.*W); % sigma33
                    D4 = 1i*(c{m}(4,5)*(k3+Wavenumber*W)+c{m}(4,4)*k3.*V); % sigma23
                    D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)+c{m}(4,5)*k3.*V); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                    Layup{m,7} = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W]; % L2
                    Layup{m,1} = L1/Layup{m,7}; % L
                    Layup{m,3} = k3;
                    Layup{m,4} = V;
                    Layup{m,5} = W;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,8} = Material{m}.Density;
                    Layup{m,9} = 1i*(c{m}(1,1)*Wavenumber+c{m}(1,6)*Wavenumber*V+c{m}(1,3)*k3.*W); % sigma11
                    Layup{m,10} = 1i*(c{m}(1,2)*Wavenumber+c{m}(2,6)*Wavenumber*V+c{m}(2,3)*k3.*W); % sigma22
                    Layup{m,11} = D3; % sigma33
                    Layup{m,12} = D4; % sigma23
                    Layup{m,13} = D5; % sigma13
                    Layup{m,14} = 1i*(c{m}(1,6)*Wavenumber+c{m}(6,6)*Wavenumber*V+c{m}(3,6)*k3.*W); % sigma12
                    Layup{m,15} = 1i*Wavenumber; % epsilon11
                    Layup{m,16} = 1i*k3.*W; % epsilon33
                    Layup{m,17} = 1i*k3.*V; % epsilon23
                    Layup{m,18} = 1i*(k3+Wavenumber*W); % epsilon13
                    Layup{m,19} = 1i*Wavenumber*V; % epsilon12
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
                MM = {M};
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                    MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                end
                if  FluidLoading
                    G = inv(MM{end});
                    k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
                    if  strcmp(ModeType,'Scholte') && Attenuation(q) ~= 0 && PhaseVelocity(q) < UpperFluid.Velocity
                        k3UpperFluid = -k3UpperFluid;
                    end
                    WUpperFluid = k3UpperFluid/Wavenumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/Wavenumber; % in the lower fluid
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
                    x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
                    if  m == 1
                        x3Total = x3;
                    else
                        x3Total = vertcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}.*x3) exp(1i*Layup{m,3}.*(Layup{m,6}-x3))];
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = -1i*AngularFrequency(q)*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)+U(5)*E(:,5)+U(6)*E(:,6));
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency(q)*(Layup{m,4}(1)*U(1)*E(:,1)+Layup{m,4}(2)*U(2)*E(:,2)+Layup{m,4}(3)*U(3)*E(:,3)+Layup{m,4}(1)*U(4)*E(:,4)+Layup{m,4}(2)*U(5)*E(:,5)+Layup{m,4}(3)*U(6)*E(:,6));
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = -1i*AngularFrequency(q)*(Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)+Layup{m,5}(3)*U(3)*E(:,3)-Layup{m,5}(1)*U(4)*E(:,4)-Layup{m,5}(2)*U(5)*E(:,5)-Layup{m,5}(3)*U(6)*E(:,6));
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
                StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,4).*conj(sigma(:,4))+epsilon(:,5).*conj(sigma(:,5))+epsilon(:,6).*conj(sigma(:,6)));
                PowerFlowDensity(1:height(v),1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)));
                PowerFlowDensity(1:height(v),2) = -.5*real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                TotalEnergy = trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)/2;
                EnergyVelocity(q,:) = PowerFlow/TotalEnergy/1e3; % ce1,ce2 (m/ms)
            end
            X{p}(:,5:6) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                Layup = cell(0);
                PhaseVelocity = X{p}(q,4)*1e3;
                Attenuation = X{p}(q,7)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                Wavenumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
                Wavenumber2 = Wavenumber^2;
                Wavenumber4 = Wavenumber2^2;
                Wavenumber6 = Wavenumber2^3;
                for m = 1:SuperLayerSize
                    rw2 = Material{m}.Density*AngularFrequency2;
                    r2w4 = rw2^2;
                    A1 = a11(m)*Wavenumber2+a12(m)*rw2;
                    A2 = a21(m)*Wavenumber4+a22(m)*rw2*Wavenumber2+a23(m)*r2w4;
                    A3 = a31(m)*Wavenumber6+a32(m)*rw2*Wavenumber4+a33(m)*r2w4*Wavenumber2+a34(m)*rw2^3;
                    d1 = A2/3-A1^2/9;
                    d2 = A1^3/27-A1*A2/6+A3/2;
                    d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                    d4 = d1/(2*d3)-d3/2;
                    d5 = d1/d3;
                    d6 = (sqrt(3)*(d3+d5)*1i)/2;
                    k3(1) = sqrt(d4-d6-A1/3);
                    k3(2) = sqrt(d4+d6-A1/3);
                    k3(3) = -sqrt(d3-d5-A1/3);
                    k32 = k3.^2;
                    m11 = c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k32;
                    m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
                    m13 = (c{m}(1,3)+c{m}(5,5))*Wavenumber*k3;
                    m22 = c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k32;
                    m23 = (c{m}(3,6)+c{m}(4,5))*Wavenumber*k3;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber*V+c{m}(3,3)*k3.*W); % sigma33
                    D4 = 1i*(c{m}(4,5)*(k3+Wavenumber*W)+c{m}(4,4)*k3.*V); % sigma23
                    D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)+c{m}(4,5)*k3.*V); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                    Layup{m,7} = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W]; % L2
                    Layup{m,1} = L1/Layup{m,7}; % L
                    Layup{m,3} = k3;
                    Layup{m,4} = V;
                    Layup{m,5} = W;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,8} = Material{m}.Density;
                    Layup{m,9} = 1i*(c{m}(1,1)*Wavenumber+c{m}(1,6)*Wavenumber*V+c{m}(1,3)*k3.*W); % sigma11
                    Layup{m,10} = 1i*(c{m}(1,2)*Wavenumber+c{m}(2,6)*Wavenumber*V+c{m}(2,3)*k3.*W); % sigma22
                    Layup{m,11} = D3; % sigma33
                    Layup{m,12} = D4; % sigma23
                    Layup{m,13} = D5; % sigma13
                    Layup{m,14} = 1i*(c{m}(1,6)*Wavenumber+c{m}(6,6)*Wavenumber*V+c{m}(3,6)*k3.*W); % sigma12
                    Layup{m,15} = 1i*Wavenumber; % epsilon11
                    Layup{m,16} = 1i*k3.*W; % epsilon33
                    Layup{m,17} = 1i*k3.*V; % epsilon23
                    Layup{m,18} = 1i*(k3+Wavenumber*W); % epsilon13
                    Layup{m,19} = 1i*Wavenumber*V; % epsilon12
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
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                    MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                end
                if  FluidLoading
                    G = inv(MM{end});
                    k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
                    if  strcmp(ModeType,'Scholte') && Attenuation ~= 0 && PhaseVelocity < UpperFluid.Velocity
                        k3UpperFluid = -k3UpperFluid;
                    end
                    WUpperFluid = k3UpperFluid/Wavenumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/Wavenumber; % in the lower fluid
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
                    x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
                    if  m == 1
                        x3Total = x3;
                    else
                        x3Total = vertcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}.*x3) exp(1i*Layup{m,3}.*(Layup{m,6}-x3))];
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
                StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,4).*conj(sigma(:,4))+epsilon(:,5).*conj(sigma(:,5))+epsilon(:,6).*conj(sigma(:,6)));
                PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)));
                PowerFlowDensity(:,2) = -.5*real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                TotalEnergy = trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)/2;
                X{p}(q,5:6) = PowerFlow/TotalEnergy/1e3; % ce1,ce2 (m/ms)
            end
        end
    elseif strcmp(ModeType,'Lamb') || strcmp(ModeType,'Scholte')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            Attenuation = X{p}(:,6).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                v = 0;
                sigma = 0;
                epsilon = 0;
                KineticEnergyDensity = 0;
                x3Total = 0;
                Layup = cell(0);
                AngularFrequency2 = AngularFrequency(q)^2;
                Wavenumber = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                Wavenumber2 = Wavenumber^2;
                Wavenumber4 = Wavenumber2^2;
                for m = 1:SuperLayerSize
                    rw2 = Material{m}.Density*AngularFrequency2;
                    A2 = a21(m)*Wavenumber2+a22(m)*rw2;
                    A3 = a31(m)*Wavenumber4+a32(m)*rw2*Wavenumber2+rw2^2;
                    k3 = [sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m)) sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m))];
                    W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber*k3);
                    D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W); % sigma33
                    D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
                    L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                    Layup{m,7} = [1 1 E;W -W.*E;E 1 1;W.*E -W]; % L2
                    Layup{m,1} = L1/Layup{m,7}; % L
                    Layup{m,3} = k3;
                    Layup{m,5} = W;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,8} = Material{m}.Density;
                    Layup{m,9} = 1i*(c{m}(1,1)*Wavenumber+c{m}(1,3)*k3.*W); % sigma11
                    Layup{m,11} = D3; % sigma33
                    Layup{m,12} = D5; % sigma13
                    Layup{m,13} = 1i*Wavenumber; % epsilon11
                    Layup{m,14} = 1i*k3.*W; % epsilon33
                    Layup{m,15} = 1i*(k3+Wavenumber*W); % epsilon13
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
                MM = {M};
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                    MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                end
                if  FluidLoading
                    G = inv(MM{end});
                    k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
                    if  strcmp(ModeType,'Scholte') && Attenuation(q) ~= 0 && PhaseVelocity(q) < UpperFluid.Velocity
                        k3UpperFluid = -k3UpperFluid;
                    end
                    WUpperFluid = k3UpperFluid/Wavenumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/Wavenumber;
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/Wavenumber;
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
                    x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
                    if  m == 1
                        x3Total = x3;
                    else
                        x3Total = vertcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}.*x3) exp(1i*Layup{m,3}.*(Layup{m,6}-x3))];
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = -1i*AngularFrequency(q)*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4));
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = -1i*AngularFrequency(q)*(Layup{m,5}(1)*U(1)*E(:,1)+Layup{m,5}(2)*U(2)*E(:,2)-Layup{m,5}(1)*U(3)*E(:,3)-Layup{m,5}(2)*U(4)*E(:,4));
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(2)*U(2)*E(:,2)+Layup{m,9}(1)*U(3)*E(:,3)+Layup{m,9}(2)*U(4)*E(:,4); % sigma11
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(2)*U(2)*E(:,2)+Layup{m,11}(1)*U(3)*E(:,3)+Layup{m,11}(2)*U(4)*E(:,4); % sigma33
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,12}(1)*U(1)*E(:,1)+Layup{m,12}(2)*U(2)*E(:,2)-Layup{m,12}(1)*U(3)*E(:,3)-Layup{m,12}(2)*U(4)*E(:,4); % sigma13
                    epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = Layup{m,13}*(U(1)*E(:,1)+U(2)*E(:,2)+U(3)*E(:,3)+U(4)*E(:,4)); % epsilon11
                    epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3) = Layup{m,14}(1)*U(1)*E(:,1)+Layup{m,14}(2)*U(2)*E(:,2)+Layup{m,14}(1)*U(3)*E(:,3)+Layup{m,14}(2)*U(4)*E(:,4); % epsilon33
                    epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),5) = Layup{m,15}(1)*U(1)*E(:,1)+Layup{m,15}(2)*U(2)*E(:,2)-Layup{m,15}(1)*U(3)*E(:,3)-Layup{m,15}(2)*U(4)*E(:,4); % epsilon13
                    KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,8}*(abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1)).^2+abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),3)).^2);
                end
                StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,5).*conj(sigma(:,5)));
                PowerFlowDensity = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                TotalEnergy = trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)/2;
                EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % ce1 (m/ms)
            end
            X{p}(:,5) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                Layup = cell(0);
                PhaseVelocity = X{p}(q,4)*1e3;
                Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                Wavenumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
                Wavenumber2 = Wavenumber^2;
                Wavenumber4 = Wavenumber2^2;
                for m = 1:SuperLayerSize
                    rw2 = Material{m}.Density*AngularFrequency2;
                    A2 = a21(m)*Wavenumber2+a22(m)*rw2;
                    A3 = a31(m)*Wavenumber4+a32(m)*rw2*Wavenumber2+rw2^2;
                    k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber*k3);
                    D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W); % sigma33
                    D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
                    L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                    Layup{m,7} = [1 1 E;W -W.*E;E 1 1;W.*E -W]; % L2
                    Layup{m,1} = L1/Layup{m,7}; % L
                    Layup{m,3} = k3;
                    Layup{m,5} = W;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,8} = Material{m}.Density;
                    Layup{m,9} = 1i*(c{m}(1,1)*Wavenumber+c{m}(1,3)*k3.*W); % sigma11
                    Layup{m,11} = D3; % sigma33
                    Layup{m,12} = D5; % sigma13
                    Layup{m,13} = 1i*Wavenumber; % epsilon11
                    Layup{m,14} = 1i*k3.*W; % epsilon33
                    Layup{m,15} = 1i*(k3+Wavenumber*W); % epsilon13
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
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                    MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                end
                if  FluidLoading
                    G = inv(MM{end});
                    k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-Wavenumber2);
                    if  strcmp(ModeType,'Scholte') && Attenuation ~= 0 && PhaseVelocity < UpperFluid.Velocity
                        k3UpperFluid = -k3UpperFluid;
                    end
                    WUpperFluid = k3UpperFluid/Wavenumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/Wavenumber;
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/Wavenumber;
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
                    x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
                    if  m == 1
                        x3Total = x3;
                    else
                        x3Total = vertcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}.*x3) exp(1i*Layup{m,3}.*(Layup{m,6}-x3))];
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
                StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,5).*conj(sigma(:,5)));
                PowerFlowDensity = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                TotalEnergy = trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)/2;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % ce1 (m/ms)
            end
        end
    elseif strcmp(ModeType,'Shear')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            Attenuation = X{p}(:,6).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                v = 0;
                sigma = 0;
                epsilon = 0;
                KineticEnergyDensity = 0;
                x3Total = 0;
                Layup = cell(0);
                AngularFrequency2 = AngularFrequency(q)^2;
                Wavenumber = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                Wavenumber2 = Wavenumber^2;
                for m = 1:SuperLayerSize
                    k3 = sqrt((Material{m}.Density*AngularFrequency2-Wavenumber2*c{m}(6,6))/c{m}(4,4));
                    if  k3 == 0
                        k3 = 1e-10;
                    end
                    E = exp(1i*k3*LayerThicknesses(m));
                    Layup{m,1} = 1i*k3*c{m}(4,4)/(E^2-1)*[-1-E^2 2*E;-2*E 1+E^2];
                    Layup{m,3} = k3;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,7} = [1 E;E 1];
                    Layup{m,8} = 1i*k3*c{m}(4,4); % sigma23
                    Layup{m,9} = 1i*Wavenumber*c{m}(6,6); % sigma12
                    Layup{m,10} = 1i*k3; % epsilon23
                    Layup{m,11} = 1i*Wavenumber; % epsilon12
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
                MM = {M};
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1,1)-MM{m-1}(2,2));
                    MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{m-1}(1,2);MM{m-1}(2,1)*N*MM{m-1}(2,1) MM{m-1}(2,2)-MM{m-1}(2,1)*N*MM{m-1}(1,2)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1,1)-MM{m-1}(2,2));
                    MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{Pattern(m)}(1,2);MM{Pattern(m)}(2,1)*N*MM{m-1}(2,1) MM{Pattern(m)}(2,2)-MM{Pattern(m)}(2,1)*N*MM{Pattern(m)}(1,2)];
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
                    x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
                    if  m == 1
                        x3Total = x3;
                    else
                        x3Total = vertcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}*x3) exp(1i*Layup{m,3}*(Layup{m,6}-x3))];
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency(q)*(U(1)*E(:,1)+U(2)*E(:,2));
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,8}(1)*U(1)*E(:,1)-Layup{m,8}(1)*U(2)*E(:,2); % sigma23
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(1)*U(2)*E(:,2); % sigma12
                    epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,10}(1)*U(1)*E(:,1)-Layup{m,10}(1)*U(2)*E(:,2); % epsilon23
                    epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(1)*U(2)*E(:,2); % epsilon12
                    KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,12}*abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2)).^2;
                end
                StrainEnergyDensity = .5*real(epsilon(:,4).*conj(sigma(:,4))+epsilon(:,6).*conj(sigma(:,6)));
                PowerFlowDensity = -.5*real(sigma(:,6).*conj(v(:,2)));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                TotalEnergy = trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)/2;
                EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % ce1 (m/ms)
            end
            X{p}(:,5) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                Layup = cell(0);
                PhaseVelocity = X{p}(q,4)*1e3;
                Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                Wavenumber = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
                Wavenumber2 = Wavenumber^2;
                for m = 1:SuperLayerSize
                    k3 = sqrt((Material{m}.Density*AngularFrequency2-Wavenumber2*c{m}(6,6))/c{m}(4,4));
                    if  k3 == 0
                        k3 = 1e-10;
                    end
                    E = exp(1i*k3*LayerThicknesses(m));
                    Layup{m,1} = 1i*k3*c{m}(4,4)/(E^2-1)*[-1-E^2 2*E;-2*E 1+E^2];
                    Layup{m,3} = k3;
                    Layup{m,6} = LayerThicknesses(m);
                    Layup{m,7} = [1 E;E 1];
                    Layup{m,8} = 1i*k3*c{m}(4,4); % sigma23
                    Layup{m,9} = 1i*Wavenumber*c{m}(6,6); % sigma12
                    Layup{m,10} = 1i*k3; % epsilon23
                    Layup{m,11} = 1i*Wavenumber; % epsilon12
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
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1,1)-MM{m-1}(2,2));
                    MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{m-1}(1,2);MM{m-1}(2,1)*N*MM{m-1}(2,1) MM{m-1}(2,2)-MM{m-1}(2,1)*N*MM{m-1}(1,2)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1,1)-MM{m-1}(2,2));
                    MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{Pattern(m)}(1,2);MM{Pattern(m)}(2,1)*N*MM{m-1}(2,1) MM{Pattern(m)}(2,2)-MM{Pattern(m)}(2,1)*N*MM{Pattern(m)}(1,2)];
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
                    x3 = (0:Layup{m,6}/SamplesX3:Layup{m,6})';
                    if  m == 1
                        x3Total = x3;
                    else
                        x3Total = vertcat(x3Total,x3+x3Total(end));
                    end
                    E = [exp(1i*Layup{m,3}*x3) exp(1i*Layup{m,3}*(Layup{m,6}-x3))];
                    v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2) = -1i*AngularFrequency*(U(1)*E(:,1)+U(2)*E(:,2));
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,8}(1)*U(1)*E(:,1)-Layup{m,8}(1)*U(2)*E(:,2); % sigma23
                    sigma((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,9}(1)*U(1)*E(:,1)+Layup{m,9}(1)*U(2)*E(:,2); % sigma12
                    epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),4) = Layup{m,10}(1)*U(1)*E(:,1)-Layup{m,10}(1)*U(2)*E(:,2); % epsilon23
                    epsilon((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),6) = Layup{m,11}(1)*U(1)*E(:,1)+Layup{m,11}(1)*U(2)*E(:,2); % epsilon12
                    KineticEnergyDensity((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),1) = .5*Layup{m,12}*abs(v((m-1)*length(x3)+1:(m-1)*length(x3)+length(x3),2)).^2;
                end
                StrainEnergyDensity = .5*real(epsilon(:,4).*conj(sigma(:,4))+epsilon(:,6).*conj(sigma(:,6)));
                PowerFlowDensity = -.5*real(sigma(:,6).*conj(v(:,2)));
                PowerFlow = trapz(x3Total,PowerFlowDensity);
                TotalEnergy = trapz(x3Total,StrainEnergyDensity+KineticEnergyDensity)/2;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % ce1 (m/ms)
            end
        end
    end
    X{p}(:,5) = fillmissing(filloutliers(X{p}(:,5),'spline','movmedian',5,'ThresholdFactor',1),'spline');
    if  contains(ModeType,'Coupled')
        X{p}(:,6) = fillmissing(filloutliers(X{p}(:,6),'spline','movmedian',5,'ThresholdFactor',1),'spline');
    end
    Counter = Counter+1;
    if  ModeTotal > 0
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
    end
    if  Stop
        return
    end
end