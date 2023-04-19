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
function [HSScholte,HAScholte,HBScholte] = FrequencySweeper_Anisotropic_Scholte(c,Delta,Material,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,SweepRange,Repetitions,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled)
Resolution = 1e-5; % (kHz)

%#ok<*AGROW>
%#ok<*MINV>
HSScholte = [];
HAScholte = [];
HBScholte = [];
if  length(SweepRange) < 2
    return
end
I = [1 1 -1;-1 -1 1;-1 -1 1];
I1 = [1 -1;-1 1];
for m = 1:SuperLayerSize
    c{m} = real(c{m});
    if  ~Decoupled
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
    else
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
    end
end
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
else
    Bisections = 1;
end
AngularFrequency = 2*pi*SweepRange*1e3;
XRough = [];
XFine = [];
if  Symmetric
    Fluid = UpperFluid;
    PhaseVelocity = (1-1e-4)*Fluid.Velocity;
    WaveNumber = AngularFrequency/PhaseVelocity;
    WaveNumber2 = WaveNumber.^2;
    WaveNumber4 = WaveNumber.^4;
    if  ~Decoupled
        WaveNumber6 = WaveNumber.^6;
        for i = 1:length(WaveNumber)
            for m = 1:SuperLayerSize
                rw2 = Material{m}.Density*AngularFrequency(i)^2;
                r2w4 = Material{m}.Density^2*AngularFrequency(i)^4;
                A1 = a11(m)*WaveNumber2(i)+a12(m)*rw2;
                A2 = a21(m)*WaveNumber4(i)+a22(m)*rw2*WaveNumber2(i)+a23(m)*r2w4;
                A3 = a31(m)*WaveNumber6(i)+a32(m)*rw2*WaveNumber4(i)+a33(m)*r2w4*WaveNumber2(i)+a34(m)*Material{m}.Density^3*AngularFrequency(i)^6;
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
                m11 = c{m}(1,1)*WaveNumber2(i)-rw2+c{m}(5,5)*k32;
                m12 = c{m}(1,6)*WaveNumber2(i)+c{m}(4,5)*k32;
                m13 = (c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3;
                m22 = c{m}(6,6)*WaveNumber2(i)-rw2+c{m}(4,4)*k32;
                m23 = (c{m}(3,6)+c{m}(4,5))*WaveNumber(i)*k3;
                V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                D3 = 1i*(c{m}(1,3)*WaveNumber(i)+c{m}(3,6)*WaveNumber(i)*V+c{m}(3,3)*k3.*W); % sigma33
                D4 = 1i*(c{m}(4,5)*(k3+WaveNumber(i)*W)+c{m}(4,4)*k3.*V); % sigma23
                D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(i)*W)+c{m}(4,5)*k3.*V); % sigma13
                E = exp(1i*k3*LayerThicknesses(m));
%                 L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
%                 L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(3);L1(1,4)=D3(1)*E(1);L1(1,5)=D3(2)*E(2);L1(1,6)=D3(3)*E(3);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=D5(3);L1(2,4)=-D5(1)*E(1);L1(2,5)=-D5(2)*E(2);L1(2,6)=-D5(3)*E(3);L1(3,1)=D4(1);L1(3,2)=D4(2);L1(3,3)=D4(3);L1(3,4)=-D4(1)*E(1);L1(3,5)=-D4(2)*E(2);L1(3,6)=-D4(3)*E(3);L1(4,1)=D3(1)*E(1);L1(4,2)=D3(2)*E(2);L1(4,3)=D3(3)*E(3);L1(4,4)=D3(1);L1(4,5)=D3(2);L1(4,6)=D3(3);L1(5,1)=D5(1)*E(1);L1(5,2)=D5(2)*E(2);L1(5,3)=D5(3)*E(3);L1(5,4)=-D5(1);L1(5,5)=-D5(2);L1(5,6)=-D5(3);L1(6,1)=D4(1)*E(1);L1(6,2)=D4(2)*E(2);L1(6,3)=D4(3)*E(3);L1(6,4)=-D4(1);L1(6,5)=-D4(2);L1(6,6)=-D4(3);
                L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(1);L2(2,2)=V(2);L2(2,3)=V(3);L2(2,4)=V(1)*E(1);L2(2,5)=V(2)*E(2);L2(2,6)=V(3)*E(3);L2(3,1)=W(1);L2(3,2)=W(2);L2(3,3)=W(3);L2(3,4)=-W(1)*E(1);L2(3,5)=-W(2)*E(2);L2(3,6)=-W(3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(1)*E(1);L2(5,2)=V(2)*E(2);L2(5,3)=V(3)*E(3);L2(5,4)=V(1);L2(5,5)=V(2);L2(5,6)=V(3);L2(6,1)=W(1)*E(1);L2(6,2)=W(2)*E(2);L2(6,3)=W(3)*E(3);L2(6,4)=-W(1);L2(6,5)=-W(2);L2(6,6)=-W(3);
                L{m} = L1/L2;
            end
            M = L{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
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
            G = inv(MM{end});
            k3Fluid = sqrt(AngularFrequency(i)^2/Fluid.Velocity^2-WaveNumber2(i));
            WFluid = k3Fluid/WaveNumber(i);
            DFluid = 1i*Fluid.Density*Fluid.Velocity^2*(WaveNumber(i)+k3Fluid*WFluid); % sigma11, sigma22, sigma33 in the fluid
            MFluid = WFluid/DFluid;
            Y(i) = (G(3,1)+MFluid)*(G(6,4)-MFluid)-G(3,4)*G(6,1);
            if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i))) && sign(real(Y(i-2))) ~= sign(real(Y(i)))
                XRough(end+1) = SweepRange(i-1);
            end
        end
% figure,plot(SweepRange,real(Y),SweepRange,abs(real(Y))),yline(0)
        if  isempty(XRough)
            disp('No higher order Scholte modes found!')
            return
        end
        for p = 1:length(XRough)
            Frequency = [XRough(p)-(SweepRange(2)-SweepRange(1)) XRough(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                AngularFrequency = 2*pi*Frequency*1e3;
                WaveNumber = AngularFrequency/PhaseVelocity;
                WaveNumber2 = WaveNumber.^2;
                WaveNumber4 = WaveNumber.^4;
                WaveNumber6 = WaveNumber.^6;
                for i = 1:length(WaveNumber)
                    for m = 1:SuperLayerSize
                        rw2 = Material{m}.Density*AngularFrequency(i)^2;
                        r2w4 = Material{m}.Density^2*AngularFrequency(i)^4;
                        A1 = a11(m)*WaveNumber2(i)+a12(m)*rw2;
                        A2 = a21(m)*WaveNumber4(i)+a22(m)*rw2*WaveNumber2(i)+a23(m)*r2w4;
                        A3 = a31(m)*WaveNumber6(i)+a32(m)*rw2*WaveNumber4(i)+a33(m)*r2w4*WaveNumber2(i)+a34(m)*Material{m}.Density^3*AngularFrequency(i)^6;
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
                        m11 = c{m}(1,1)*WaveNumber2(i)-rw2+c{m}(5,5)*k32;
                        m12 = c{m}(1,6)*WaveNumber2(i)+c{m}(4,5)*k32;
                        m13 = (c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3;
                        m22 = c{m}(6,6)*WaveNumber2(i)-rw2+c{m}(4,4)*k32;
                        m23 = (c{m}(3,6)+c{m}(4,5))*WaveNumber(i)*k3;
                        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                        D3 = 1i*(c{m}(1,3)*WaveNumber(i)+c{m}(3,6)*WaveNumber(i)*V+c{m}(3,3)*k3.*W); % sigma33
                        D4 = 1i*(c{m}(4,5)*(k3+WaveNumber(i)*W)+c{m}(4,4)*k3.*V); % sigma23
                        D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(i)*W)+c{m}(4,5)*k3.*V); % sigma13
                        E = exp(1i*k3*LayerThicknesses(m));
                        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                        L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                        L{m} = L1/L2;
                    end
                    M = L{1};
                    for m = 2:SuperLayerSize
                        N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                        M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
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
                    G = inv(MM{end});
                    k3Fluid = sqrt(AngularFrequency(i)^2/Fluid.Velocity^2-WaveNumber2(i));
                    WFluid = k3Fluid/WaveNumber(i);
                    DFluid = 1i*Fluid.Density*Fluid.Velocity^2*(WaveNumber(i)+k3Fluid*WFluid); % sigma11, sigma22, sigma33 in the fluid
                    MFluid = WFluid/DFluid;                    
                    Y(i) = abs((G(3,1)+MFluid)*(G(6,4)-MFluid)-G(3,4)*G(6,1));
                    if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                        if  o == Bisections
                            XFine(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        for p = 1:length(XFine)
            AngularFrequency = 2*pi*XFine(p)*1e3;
            WaveNumber = AngularFrequency/PhaseVelocity;
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
                L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                L{m} = L1/L2;
            end
            M = L{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
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
            G = inv(MM{end});
            k3Fluid = sqrt(AngularFrequency^2/Fluid.Velocity^2-WaveNumber2);
            WFluid = k3Fluid/WaveNumber;
            DFluid = 1i*Fluid.Density*Fluid.Velocity^2*(WaveNumber+k3Fluid*WFluid); % sigma11, sigma22, sigma33 in the fluid            
            UFluid = -(WFluid+G(3,1)*DFluid)/(G(3,4)*DFluid);
            if  sign(real(UFluid)) == 1 % this works in the unattenuated case only
                HSScholte(end+1) = XFine(p);
            else
                HAScholte(end+1) = XFine(p);
            end
        end        
    else
        for i = 1:length(WaveNumber)
            for m = 1:SuperLayerSize
                rw2 = Material{m}.Density*AngularFrequency(i)^2;
                A2 = a21(m)*WaveNumber2(i)+a22(m)*rw2;
                A3 = a31(m)*WaveNumber4(i)+a32(m)*rw2*WaveNumber2(i)+Material{m}.Density^2*AngularFrequency(i)^4;
                k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                W = (rw2-c{m}(1,1)*WaveNumber2(i)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3);
                D3 = 1i*(c{m}(1,3)*WaveNumber(i)+c{m}(3,3)*k3.*W); % sigma33
                D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(i)*W)); % sigma13
                E = exp(1i*k3*LayerThicknesses(m));
%                 L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
%                 L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(1)*E(1);L1(1,4)=D3(2)*E(2);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=-D5(1)*E(1);L1(2,4)=-D5(2)*E(2);L1(3,1)=D3(1)*E(1);L1(3,2)=D3(2)*E(2);L1(3,3)=D3(1);L1(3,4)=D3(2);L1(4,1)=D5(1)*E(1);L1(4,2)=D5(2)*E(2);L1(4,3)=-D5(1);L1(4,4)=-D5(2);
                L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(1);L2(2,2)=W(2);L2(2,3)=-W(1)*E(1);L2(2,4)=-W(2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(1)*E(1);L2(4,2)=W(2)*E(2);L2(4,3)=-W(1);L2(4,4)=-W(2);
                L{m} = L1/L2;
            end
            M = L{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
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
            G = inv(MM{end});
            k3Fluid = sqrt(AngularFrequency(i)^2/Fluid.Velocity^2-WaveNumber2(i));
            WFluid = k3Fluid/WaveNumber(i);
            DFluid = 1i*Fluid.Density*Fluid.Velocity^2*(WaveNumber(i)+k3Fluid*WFluid); % sigma11, sigma22, sigma33 in the fluid
            MFluid = WFluid/DFluid;            
            Y(i) = (G(2,1)+MFluid)*(G(4,3)-MFluid)-G(2,3)*G(4,1);
            if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i))) && sign(real(Y(i-2))) ~= sign(real(Y(i)))
                XRough(end+1) = SweepRange(i-1);
            end
        end
% figure,plot(SweepRange,real(Y),SweepRange,abs(real(Y))),yline(0)
        if  isempty(XRough)
            disp('No higher order Scholte modes found!')
            return
        end
        for p = 1:length(XRough)
            Frequency = [XRough(p)-(SweepRange(2)-SweepRange(1)) XRough(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                AngularFrequency = 2*pi*Frequency*1e3;
                WaveNumber = AngularFrequency/PhaseVelocity;
                WaveNumber2 = WaveNumber.^2;
                WaveNumber4 = WaveNumber.^4;
                for i = 1:length(WaveNumber)
                    for m = 1:SuperLayerSize
                        rw2 = Material{m}.Density*AngularFrequency(i)^2;
                        A2 = a21(m)*WaveNumber2(i)+a22(m)*rw2;
                        A3 = a31(m)*WaveNumber4(i)+a32(m)*rw2*WaveNumber2(i)+Material{m}.Density^2*AngularFrequency(i)^4;
                        k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                        k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                        W = (rw2-c{m}(1,1)*WaveNumber2(i)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3);
                        D3 = 1i*(c{m}(1,3)*WaveNumber(i)+c{m}(3,3)*k3.*W); % sigma33
                        D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(i)*W)); % sigma13
                        E = exp(1i*k3*LayerThicknesses(m));
                        L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                        L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                        L{m} = L1/L2;
                    end
                    M = L{1};
                    for m = 2:SuperLayerSize
                        N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                        M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
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
                    G = inv(MM{end});
                    k3Fluid = sqrt(AngularFrequency(i)^2/Fluid.Velocity^2-WaveNumber2(i));
                    WFluid = k3Fluid/WaveNumber(i);
                    DFluid = 1i*Fluid.Density*Fluid.Velocity^2*(WaveNumber(i)+k3Fluid*WFluid); % sigma11, sigma22, sigma33 in the fluid
                    MFluid = WFluid/DFluid;                    
                    Y(i) = abs((G(2,1)+MFluid)*(G(4,3)-MFluid)-G(2,3)*G(4,1));
                    if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i) 
                        if  o == Bisections
                            XFine(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        for p = 1:length(XFine)
            AngularFrequency = 2*pi*XFine(p)*1e3;
            WaveNumber = AngularFrequency/PhaseVelocity;
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
                L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                L{m} = L1/L2;
            end
            M = L{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
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
            G = inv(MM{end});
            k3Fluid = sqrt(AngularFrequency^2/Fluid.Velocity^2-WaveNumber2);
            WFluid = k3Fluid/WaveNumber;
            DFluid = 1i*Fluid.Density*Fluid.Velocity^2*(WaveNumber+k3Fluid*WFluid); % sigma11, sigma22, sigma33 in the fluid            
            UFluid = -(WFluid+G(2,1)*DFluid)/(G(2,3)*DFluid);
            if  sign(real(UFluid)) == 1 % this works in the unattenuated case only
                HSScholte(end+1) = XFine(p);
            else
                HAScholte(end+1) = XFine(p);
            end
        end
    end
    if  isempty(HAScholte)
        HSScholte = [];
    else
        HSScholte(HSScholte < HAScholte(1)) = [];
    end
    if  isempty(HSScholte) && isempty(HAScholte)
        disp('No higher order Scholte modes found!')
        return
    end
    String = ['Frq. @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode         Frq.(kHz)'];
    if  any(HAScholte)
        for i = 1:length(HAScholte)
            if  i < 10
                if  HAScholte(i) < 1e2
                    String = append(String,newline,'AScholte',num2str(i),'      ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e2 && HAScholte(i) < 1e3
                    String = append(String,newline,'AScholte',num2str(i),'     ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e3 && HAScholte(i) < 1e4
                    String = append(String,newline,'AScholte',num2str(i),'    ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e4
                    String = append(String,newline,'AScholte',num2str(i),'   ',num2str(HAScholte(i),'%.3f'));
                end
            else
                if  HAScholte(i) < 1e2
                    String = append(String,newline,'AScholte',num2str(i),'     ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e2 && HAScholte(i) < 1e3
                    String = append(String,newline,'AScholte',num2str(i),'    ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e3 && HAScholte(i) < 1e4
                    String = append(String,newline,'AScholte',num2str(i),'   ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e4
                    String = append(String,newline,'AScholte',num2str(i),'  ',num2str(HAScholte(i),'%.3f'));
                end
            end
        end
    end
    String = append(String,newline);
    if  any(HSScholte)
        for i = 1:length(HSScholte)
            if  i < 10
                if  HSScholte(i) < 1e2
                    String = append(String,newline,'SScholte',num2str(i),'      ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e2 && HSScholte(i) < 1e3
                    String = append(String,newline,'SScholte',num2str(i),'     ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e3 && HSScholte(i) < 1e4
                    String = append(String,newline,'SScholte',num2str(i),'    ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e4
                    String = append(String,newline,'SScholte',num2str(i),'   ',num2str(HSScholte(i),'%.3f'));
                end
            else
                if  HSScholte(i) < 1e2
                    String = append(String,newline,'SScholte',num2str(i),'     ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e2 && HSScholte(i) < 1e3
                    String = append(String,newline,'SScholte',num2str(i),'    ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e3 && HSScholte(i) < 1e4
                    String = append(String,newline,'SScholte',num2str(i),'   ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e4
                    String = append(String,newline,'SScholte',num2str(i),'  ',num2str(HSScholte(i),'%.3f'));
                end
            end
        end
    end
    String = append(String,newline,newline);
    if  any(HAScholte)
        String = append(String,'AScholte: ',num2str(length(HAScholte)));
    end
    String = append(String,newline);
    if  any(HSScholte)
        String = append(String,'SScholte: ',num2str(length(HSScholte)));
    end
else
%     if  ToggleUpperFluid && ToggleLowerFluid
%         if  strcmp(UpperFluid.Name,LowerFluid.Name)
%             FluidVelocity = UpperFluid.Velocity;
%         else
%             if  UpperFluid.Velocity > LowerFluid.Velocity
%                 FluidVelocity = UpperFluid.Velocity;
%             else
%                 FluidVelocity = LowerFluid.Velocity;
%             end
%         end
%     elseif ToggleUpperFluid && ~ToggleLowerFluid
%         FluidVelocity = UpperFluid.Velocity;
%     elseif ~ToggleUpperFluid && ToggleLowerFluid
%         FluidVelocity = LowerFluid.Velocity;
%     end
%     PhaseVelocity = (1-1e-4)*FluidVelocity;
%     WaveNumber = AngularFrequency./PhaseVelocity;
%     if  ~Decoupled
%         for i = 1:length(WaveNumber)
%             for m = 1:SuperLayerSize
%                 A1 = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)*WaveNumber(i)^2-c{m}(1,3)^2*c{m}(4,4)*WaveNumber(i)^2+2*c{m}(1,3)*c{m}(3,6)*c{m}(4,5)*WaveNumber(i)^2-2*c{m}(1,3)*c{m}(4,4)*c{m}(5,5)*WaveNumber(i)^2+2*c{m}(1,3)*c{m}(4,5)^2*WaveNumber(i)^2-2*c{m}(1,6)*c{m}(3,3)*c{m}(4,5)*WaveNumber(i)^2-c{m}(3,3)*c{m}(4,4)*Material{m}.Density*AngularFrequency(i)^2+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^2-c{m}(3,3)*c{m}(5,5)*Material{m}.Density*AngularFrequency(i)^2-c{m}(3,6)^2*c{m}(5,5)*WaveNumber(i)^2-c{m}(4,4)*c{m}(5,5)*Material{m}.Density*AngularFrequency(i)^2+c{m}(4,5)^2*Material{m}.Density*AngularFrequency(i)^2)/Delta(m);
%                 A2 = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)*WaveNumber(i)^4-c{m}(1,1)*c{m}(3,3)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(1,1)*c{m}(3,6)^2*WaveNumber(i)^4-2*c{m}(1,1)*c{m}(3,6)*c{m}(4,5)*WaveNumber(i)^4+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)*WaveNumber(i)^4-c{m}(1,1)*c{m}(4,4)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(1,1)*c{m}(4,5)^2*WaveNumber(i)^4-c{m}(1,3)^2*c{m}(6,6)*WaveNumber(i)^4+c{m}(1,3)^2*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+2*c{m}(1,3)*c{m}(1,6)*c{m}(3,6)*WaveNumber(i)^4+2*c{m}(1,3)*c{m}(1,6)*c{m}(4,5)*WaveNumber(i)^4-2*c{m}(1,3)*c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^4+2*c{m}(1,3)*c{m}(5,5)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(1,6)^2*c{m}(3,3)*WaveNumber(i)^4+2*c{m}(1,6)*c{m}(3,6)*c{m}(5,5)*WaveNumber(i)^4+2*c{m}(1,6)*c{m}(4,5)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(3,3)*c{m}(6,6)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+c{m}(3,3)*Material{m}.Density^2*AngularFrequency(i)^4+c{m}(3,6)^2*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+2*c{m}(3,6)*c{m}(4,5)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(4,4)*c{m}(5,5)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+c{m}(4,4)*Material{m}.Density^2*AngularFrequency(i)^4+c{m}(4,5)^2*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+c{m}(5,5)*Material{m}.Density^2*AngularFrequency(i)^4)/Delta(m);
%                 A3 = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^6-c{m}(1,1)*c{m}(5,5)*WaveNumber(i)^4*Material{m}.Density*AngularFrequency(i)^2-c{m}(1,1)*c{m}(6,6)*WaveNumber(i)^4*Material{m}.Density*AngularFrequency(i)^2+c{m}(1,1)*WaveNumber(i)^2*Material{m}.Density^2*AngularFrequency(i)^4-c{m}(1,6)^2*c{m}(5,5)*WaveNumber(i)^6+c{m}(1,6)^2*WaveNumber(i)^4*Material{m}.Density*AngularFrequency(i)^2-c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^4*Material{m}.Density*AngularFrequency(i)^2+c{m}(5,5)*WaveNumber(i)^2*Material{m}.Density^2*AngularFrequency(i)^4+c{m}(6,6)*WaveNumber(i)^2*Material{m}.Density^2*AngularFrequency(i)^4-Material{m}.Density^3*AngularFrequency(i)^6)/Delta(m);
%                 k3(1) = ((-A1^2/9+A2/3)/(2*(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))-(((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)/2-(3^(1/2)*((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)+(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))*1i)/2-A1/3)^(1/2);
%                 k3(2) = ((-A1^2/9+A2/3)/(2*(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))-(((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)/2+(3^(1/2)*((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)+(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))*1i)/2-A1/3)^(1/2);
%                 k3(3) = -((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)-A1/3-(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))^(1/2);
%                 V = ((c{m}(1,1)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(5,5)*k3.^2).*((c{m}(3,6)+c{m}(4,5))*WaveNumber(i)*k3)-((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3).*(c{m}(1,6)*WaveNumber(i)^2+c{m}(4,5)*k3.^2))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3).*(c{m}(6,6)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(4,4)*k3.^2)-(c{m}(1,6)*WaveNumber(i)^2+c{m}(4,5)*k3.^2).*((c{m}(3,6)+c{m}(4,5))*WaveNumber(i)*k3));
%                 W = ((c{m}(1,1)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(5,5)*k3.^2).*(c{m}(6,6)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(4,4)*k3.^2)-(c{m}(1,6)*WaveNumber(i)^2+c{m}(4,5)*k3.^2).^2)./((c{m}(1,6)*WaveNumber(i)^2+c{m}(4,5)*k3.^2).*((c{m}(3,6)+c{m}(4,5))*WaveNumber(i)*k3)-((c{m}(6,6)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(4,4)*k3.^2).*((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3)));
%                 D3 = 1i*(c{m}(1,3)*WaveNumber(i)+c{m}(3,6)*WaveNumber(i)*V+c{m}(3,3)*k3.*W); % sigma33
%                 D4 = 1i*(c{m}(4,5)*(k3+WaveNumber(i)*W)+c{m}(4,4)*k3.*V); % sigma23
%                 D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(i)*W)+c{m}(4,5)*k3.*V); % sigma13
%                 E = exp(1i*k3*LayerThicknesses(m));
%                 L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
%                 L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
%                 L{m} = L1/L2;
%             end
%             M = L{1};
%             for m = 2:SuperLayerSize
%                 N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
%                 M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
%             end
%             MM{1} = M;
%             for m = 2:Repetitions
%                 N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
%                 MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
%             end
%             if  SymmetricSystem
%                 N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
%                 MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
%             end
%             G = inv(MM{end});
%             if  ToggleUpperFluid && ToggleLowerFluid
%                 k3UpperFluid = sqrt(AngularFrequency(i)^2/UpperFluid.Velocity^2-WaveNumber(i)^2);
%                 k3LowerFluid = sqrt(AngularFrequency(i)^2/LowerFluid.Velocity^2-WaveNumber(i)^2);
%                 WUpperFluid = k3UpperFluid/WaveNumber(i);
%                 WLowerFluid = k3LowerFluid/WaveNumber(i);
%                 DUpperFluid = 1i*UpperFluid.Density*UpperFluid.Velocity^2*(WaveNumber(i)+k3UpperFluid*WUpperFluid); % sigma11, sigma22, sigma33 in the upper fluid
%                 DLowerFluid = 1i*LowerFluid.Density*LowerFluid.Velocity^2*(WaveNumber(i)+k3LowerFluid*WLowerFluid); % in the lower fluid
%                 Y(i) = (WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid;
%             elseif ToggleUpperFluid && ~ToggleLowerFluid
%                 k3UpperFluid = sqrt(AngularFrequency(i)^2/UpperFluid.Velocity^2-WaveNumber(i)^2);
%                 WUpperFluid = k3UpperFluid/WaveNumber(i);
%                 DUpperFluid = 1i*UpperFluid.Density*UpperFluid.Velocity^2*(WaveNumber(i)+k3UpperFluid*WUpperFluid);
%                 Y(i) = WUpperFluid+G(3,1)*DUpperFluid;
%             elseif ~ToggleUpperFluid && ToggleLowerFluid
%                 k3LowerFluid = sqrt(AngularFrequency(i)^2/LowerFluid.Velocity^2-WaveNumber(i)^2);
%                 WLowerFluid = k3LowerFluid/WaveNumber(i);
%                 DLowerFluid = 1i*LowerFluid.Density*LowerFluid.Velocity^2*(WaveNumber(i)+k3LowerFluid*WLowerFluid);
%                 Y(i) = WLowerFluid-G(6,4)*DLowerFluid;
%             end
%             if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
%                 XRough(end+1) = SweepRange(i-1);
%             end
%         end
% % figure,plot(SweepRange,real(Y),SweepRange,abs(real(Y))),yline(0)
%         if  isempty(XRough)
%             disp('No higher order Scholte modes found!')
%             return
%         end
%         for p = 1:length(XRough)
%             Frequency = [XRough(p)-(SweepRange(2)-SweepRange(1)) XRough(p)+(SweepRange(2)-SweepRange(1))];
%             for o = 1:Bisections
%                 Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
%                 AngularFrequency = 2*pi*Frequency*1e3;
%                 WaveNumber = AngularFrequency/PhaseVelocity;
%                 for i = 1:length(WaveNumber)
%                     for m = 1:SuperLayerSize
%                         A1 = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)*WaveNumber(i)^2-c{m}(1,3)^2*c{m}(4,4)*WaveNumber(i)^2+2*c{m}(1,3)*c{m}(3,6)*c{m}(4,5)*WaveNumber(i)^2-2*c{m}(1,3)*c{m}(4,4)*c{m}(5,5)*WaveNumber(i)^2+2*c{m}(1,3)*c{m}(4,5)^2*WaveNumber(i)^2-2*c{m}(1,6)*c{m}(3,3)*c{m}(4,5)*WaveNumber(i)^2-c{m}(3,3)*c{m}(4,4)*Material{m}.Density*AngularFrequency(i)^2+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^2-c{m}(3,3)*c{m}(5,5)*Material{m}.Density*AngularFrequency(i)^2-c{m}(3,6)^2*c{m}(5,5)*WaveNumber(i)^2-c{m}(4,4)*c{m}(5,5)*Material{m}.Density*AngularFrequency(i)^2+c{m}(4,5)^2*Material{m}.Density*AngularFrequency(i)^2)/Delta(m);
%                         A2 = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)*WaveNumber(i)^4-c{m}(1,1)*c{m}(3,3)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(1,1)*c{m}(3,6)^2*WaveNumber(i)^4-2*c{m}(1,1)*c{m}(3,6)*c{m}(4,5)*WaveNumber(i)^4+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)*WaveNumber(i)^4-c{m}(1,1)*c{m}(4,4)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(1,1)*c{m}(4,5)^2*WaveNumber(i)^4-c{m}(1,3)^2*c{m}(6,6)*WaveNumber(i)^4+c{m}(1,3)^2*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+2*c{m}(1,3)*c{m}(1,6)*c{m}(3,6)*WaveNumber(i)^4+2*c{m}(1,3)*c{m}(1,6)*c{m}(4,5)*WaveNumber(i)^4-2*c{m}(1,3)*c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^4+2*c{m}(1,3)*c{m}(5,5)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(1,6)^2*c{m}(3,3)*WaveNumber(i)^4+2*c{m}(1,6)*c{m}(3,6)*c{m}(5,5)*WaveNumber(i)^4+2*c{m}(1,6)*c{m}(4,5)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(3,3)*c{m}(6,6)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+c{m}(3,3)*Material{m}.Density^2*AngularFrequency(i)^4+c{m}(3,6)^2*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+2*c{m}(3,6)*c{m}(4,5)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(4,4)*c{m}(5,5)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+c{m}(4,4)*Material{m}.Density^2*AngularFrequency(i)^4+c{m}(4,5)^2*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2-c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^2*Material{m}.Density*AngularFrequency(i)^2+c{m}(5,5)*Material{m}.Density^2*AngularFrequency(i)^4)/Delta(m);
%                         A3 = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^6-c{m}(1,1)*c{m}(5,5)*WaveNumber(i)^4*Material{m}.Density*AngularFrequency(i)^2-c{m}(1,1)*c{m}(6,6)*WaveNumber(i)^4*Material{m}.Density*AngularFrequency(i)^2+c{m}(1,1)*WaveNumber(i)^2*Material{m}.Density^2*AngularFrequency(i)^4-c{m}(1,6)^2*c{m}(5,5)*WaveNumber(i)^6+c{m}(1,6)^2*WaveNumber(i)^4*Material{m}.Density*AngularFrequency(i)^2-c{m}(5,5)*c{m}(6,6)*WaveNumber(i)^4*Material{m}.Density*AngularFrequency(i)^2+c{m}(5,5)*WaveNumber(i)^2*Material{m}.Density^2*AngularFrequency(i)^4+c{m}(6,6)*WaveNumber(i)^2*Material{m}.Density^2*AngularFrequency(i)^4-Material{m}.Density^3*AngularFrequency(i)^6)/Delta(m);
%                         k3(1) = ((-A1^2/9+A2/3)/(2*(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))-(((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)/2-(3^(1/2)*((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)+(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))*1i)/2-A1/3)^(1/2);
%                         k3(2) = ((-A1^2/9+A2/3)/(2*(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))-(((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)/2+(3^(1/2)*((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)+(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))*1i)/2-A1/3)^(1/2);
%                         k3(3) = -((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)-A1/3-(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))^(1/2);
%                         V = ((c{m}(1,1)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(5,5)*k3.^2).*((c{m}(3,6)+c{m}(4,5))*WaveNumber(i)*k3)-((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3).*(c{m}(1,6)*WaveNumber(i)^2+c{m}(4,5)*k3.^2))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3).*(c{m}(6,6)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(4,4)*k3.^2)-(c{m}(1,6)*WaveNumber(i)^2+c{m}(4,5)*k3.^2).*((c{m}(3,6)+c{m}(4,5))*WaveNumber(i)*k3));
%                         W = ((c{m}(1,1)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(5,5)*k3.^2).*(c{m}(6,6)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(4,4)*k3.^2)-(c{m}(1,6)*WaveNumber(i)^2+c{m}(4,5)*k3.^2).^2)./((c{m}(1,6)*WaveNumber(i)^2+c{m}(4,5)*k3.^2).*((c{m}(3,6)+c{m}(4,5))*WaveNumber(i)*k3)-((c{m}(6,6)*WaveNumber(i)^2-Material{m}.Density*AngularFrequency(i)^2+c{m}(4,4)*k3.^2).*((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3)));
%                         D3 = 1i*(c{m}(1,3)*WaveNumber(i)+c{m}(3,6)*WaveNumber(i)*V+c{m}(3,3)*k3.*W); % sigma33
%                         D4 = 1i*(c{m}(4,5)*(k3+WaveNumber(i)*W)+c{m}(4,4)*k3.*V); % sigma23
%                         D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(i)*W)+c{m}(4,5)*k3.*V); % sigma13
%                         E = exp(1i*k3*LayerThicknesses(m));
%                         L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
%                         L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
%                         L{m} = L1/L2;
%                     end
%                     M = L{1};
%                     for m = 2:SuperLayerSize
%                         N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
%                         M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
%                     end
%                     MM{1} = M;
%                     for m = 2:Repetitions
%                         N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
%                         MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
%                     end
%                     if  SymmetricSystem
%                         N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
%                         MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
%                     end
%                     G = inv(MM{end});
%                     if  ToggleUpperFluid && ToggleLowerFluid
%                         k3UpperFluid = sqrt(AngularFrequency(i)^2/UpperFluid.Velocity^2-WaveNumber(i)^2);
%                         k3LowerFluid = sqrt(AngularFrequency(i)^2/LowerFluid.Velocity^2-WaveNumber(i)^2);
%                         WUpperFluid = k3UpperFluid/WaveNumber(i);
%                         WLowerFluid = k3LowerFluid/WaveNumber(i);
%                         DUpperFluid = 1i*UpperFluid.Density*UpperFluid.Velocity^2*(WaveNumber(i)+k3UpperFluid*WUpperFluid); % sigma11, sigma22, sigma33 in the upper fluid
%                         DLowerFluid = 1i*LowerFluid.Density*LowerFluid.Velocity^2*(WaveNumber(i)+k3LowerFluid*WLowerFluid); % in the lower fluid
%                         Y(i) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
%                     elseif ToggleUpperFluid && ~ToggleLowerFluid
%                         k3UpperFluid = sqrt(AngularFrequency(i)^2/UpperFluid.Velocity^2-WaveNumber(i)^2);
%                         WUpperFluid = k3UpperFluid/WaveNumber(i);
%                         DUpperFluid = 1i*UpperFluid.Density*UpperFluid.Velocity^2*(WaveNumber(i)+k3UpperFluid*WUpperFluid);
%                         Y(i) = abs(WUpperFluid+G(3,1)*DUpperFluid);
%                     elseif ~ToggleUpperFluid && ToggleLowerFluid
%                         k3LowerFluid = sqrt(AngularFrequency(i)^2/LowerFluid.Velocity^2-WaveNumber(i)^2);
%                         WLowerFluid = k3LowerFluid/WaveNumber(i);
%                         DLowerFluid = 1i*LowerFluid.Density*LowerFluid.Velocity^2*(WaveNumber(i)+k3LowerFluid*WLowerFluid);
%                         Y(i) = abs(WLowerFluid-G(6,4)*DLowerFluid);
%                     end
%                     if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)  
%                         if  o == Bisections
%                             HBScholte(end+1) = Frequency(i-1);
%                         end
%                         Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
%                         break
%                     end
%                 end
%             end
%         end
%     else
%         for i = 1:length(WaveNumber)
%             for m = 1:SuperLayerSize
%                 A2 = (c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2)*WaveNumber(i)^2-(c{m}(3,3)+c{m}(5,5))*Material{m}.Density*AngularFrequency(i)^2;
%                 A3 = c{m}(1,1)*c{m}(5,5)*WaveNumber(i)^4-(c{m}(1,1)+c{m}(5,5))*Material{m}.Density*AngularFrequency(i)^2*WaveNumber(i)^2+Material{m}.Density^2*AngularFrequency(i)^4;
%                 k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
%                 k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
%                 W = (Material{m}.Density*AngularFrequency(i)^2-c{m}(1,1)*WaveNumber(i)^2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3);
%                 D3 = 1i*(c{m}(1,3)*WaveNumber(i)+c{m}(3,3)*k3.*W); % sigma33
%                 D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(i)*W)); % sigma13
%                 E = exp(1i*k3*LayerThicknesses(m));
%                 L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
%                 L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
%                 L{m} = L1/L2;
%             end
%             M = L{1};
%             for m = 2:SuperLayerSize
%                 N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
%                 M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
%             end
%             MM{1} = M;
%             for m = 2:Repetitions
%                 N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
%                 MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
%             end
%             if  SymmetricSystem
%                 N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
%                 MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
%             end
%             G = inv(MM{end});
%             if  ToggleUpperFluid && ToggleLowerFluid
%                 k3UpperFluid = sqrt(AngularFrequency(i)^2/UpperFluid.Velocity^2-WaveNumber(i)^2);
%                 k3LowerFluid = sqrt(AngularFrequency(i)^2/LowerFluid.Velocity^2-WaveNumber(i)^2);
%                 WUpperFluid = k3UpperFluid/WaveNumber(i);
%                 WLowerFluid = k3LowerFluid/WaveNumber(i);
%                 DUpperFluid = 1i*UpperFluid.Density*UpperFluid.Velocity^2*(WaveNumber(i)+k3UpperFluid*WUpperFluid); % sigma11, sigma22, sigma33 in the upper fluid
%                 DLowerFluid = 1i*LowerFluid.Density*LowerFluid.Velocity^2*(WaveNumber(i)+k3LowerFluid*WLowerFluid); % in the lower fluid
%                 Y(i) = (WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid;
%             elseif ToggleUpperFluid && ~ToggleLowerFluid
%                 k3UpperFluid = sqrt(AngularFrequency(i)^2/UpperFluid.Velocity^2-WaveNumber(i)^2);
%                 WUpperFluid = k3UpperFluid/WaveNumber(i);
%                 DUpperFluid = 1i*UpperFluid.Density*UpperFluid.Velocity^2*(WaveNumber(i)+k3UpperFluid*WUpperFluid);
%                 Y(i) = WUpperFluid+G(2,1)*DUpperFluid;
%             elseif ~ToggleUpperFluid && ToggleLowerFluid
%                 k3LowerFluid = sqrt(AngularFrequency(i)^2/LowerFluid.Velocity^2-WaveNumber(i)^2);
%                 WLowerFluid = k3LowerFluid/WaveNumber(i);
%                 DLowerFluid = 1i*LowerFluid.Density*LowerFluid.Velocity^2*(WaveNumber(i)+k3LowerFluid*WLowerFluid);
%                 Y(i) = WLowerFluid-G(4,3)*DLowerFluid;
%             end
%             if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
%                 XRough(end+1) = SweepRange(i-1);
%             end
%         end
% % figure,plot(SweepRange,real(Y),SweepRange,abs(real(Y))),yline(0)
%         if  isempty(XRough)
%             disp('No higher order Scholte modes found!')
%             return
%         end
%         for p = 1:length(XRough)
%             Frequency = [XRough(p)-(SweepRange(2)-SweepRange(1)) XRough(p)+(SweepRange(2)-SweepRange(1))];
%             for o = 1:Bisections
%                 Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
%                 AngularFrequency = 2*pi*Frequency*1e3;
%                 WaveNumber = AngularFrequency/PhaseVelocity;
%                 for i = 1:length(WaveNumber)
%                     for m = 1:SuperLayerSize
%                         A2 = (c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2)*WaveNumber(i)^2-(c{m}(3,3)+c{m}(5,5))*Material{m}.Density*AngularFrequency(i)^2;
%                         A3 = c{m}(1,1)*c{m}(5,5)*WaveNumber(i)^4-(c{m}(1,1)+c{m}(5,5))*Material{m}.Density*AngularFrequency(i)^2*WaveNumber(i)^2+Material{m}.Density^2*AngularFrequency(i)^4;
%                         k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
%                         k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
%                         W = (Material{m}.Density*AngularFrequency(i)^2-c{m}(1,1)*WaveNumber(i)^2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*WaveNumber(i)*k3);
%                         D3 = 1i*(c{m}(1,3)*WaveNumber(i)+c{m}(3,3)*k3.*W); % sigma33
%                         D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(i)*W)); % sigma13
%                         E = exp(1i*k3*LayerThicknesses(m));
%                         L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
%                         L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
%                         L{m} = L1/L2;
%                     end
%                     M = L{1};
%                     for m = 2:SuperLayerSize
%                         N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
%                         M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
%                     end
%                     MM{1} = M;
%                     for m = 2:Repetitions
%                         N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
%                         MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
%                     end
%                     if  SymmetricSystem
%                         N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
%                         MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
%                     end
%                     G = inv(MM{end});
%                     if  ToggleUpperFluid && ToggleLowerFluid
%                         k3UpperFluid = sqrt(AngularFrequency(i)^2/UpperFluid.Velocity^2-WaveNumber(i)^2);
%                         k3LowerFluid = sqrt(AngularFrequency(i)^2/LowerFluid.Velocity^2-WaveNumber(i)^2);
%                         WUpperFluid = k3UpperFluid/WaveNumber(i);
%                         WLowerFluid = k3LowerFluid/WaveNumber(i);
%                         DUpperFluid = 1i*UpperFluid.Density*UpperFluid.Velocity^2*(WaveNumber(i)+k3UpperFluid*WUpperFluid); % sigma11, sigma22, sigma33 in the upper fluid
%                         DLowerFluid = 1i*LowerFluid.Density*LowerFluid.Velocity^2*(WaveNumber(i)+k3LowerFluid*WLowerFluid); % in the lower fluid
%                         Y(i) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
%                     elseif ToggleUpperFluid && ~ToggleLowerFluid
%                         k3UpperFluid = sqrt(AngularFrequency(i)^2/UpperFluid.Velocity^2-WaveNumber(i)^2);
%                         WUpperFluid = k3UpperFluid/WaveNumber(i);
%                         DUpperFluid = 1i*UpperFluid.Density*UpperFluid.Velocity^2*(WaveNumber(i)+k3UpperFluid*WUpperFluid);
%                         Y(i) = abs(WUpperFluid+G(2,1)*DUpperFluid);
%                     elseif ~ToggleUpperFluid && ToggleLowerFluid
%                         k3LowerFluid = sqrt(AngularFrequency(i)^2/LowerFluid.Velocity^2-WaveNumber(i)^2);
%                         WLowerFluid = k3LowerFluid/WaveNumber(i);
%                         DLowerFluid = 1i*LowerFluid.Density*LowerFluid.Velocity^2*(WaveNumber(i)+k3LowerFluid*WLowerFluid);
%                         Y(i) = abs(WLowerFluid-G(4,3)*DLowerFluid);
%                     end
%                     if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
%                         if  o == Bisections
%                             HBScholte(end+1) = Frequency(i-1);
%                         end
%                         Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
%                         break
%                     end
%                 end
%             end
%         end
%     end
%     if  isempty(HBScholte)
%         disp('No higher order Scholte modes found!')
%         return
%     end
%     String = ['Frq. @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode         Frq.(kHz)'];
%     if  any(HBScholte)
%         for i = 1:length(HBScholte)
%             if  i < 9
%                 if  HBScholte(i) < 1e2
%                     String = append(String,newline,'BScholte',num2str(i+1),'      ',num2str(HBScholte(i),'%.3f'));
%                 elseif HBScholte(i) >= 1e2 && HBScholte(i) < 1e3
%                     String = append(String,newline,'BScholte',num2str(i+1),'     ',num2str(HBScholte(i),'%.3f'));
%                 elseif HBScholte(i) >= 1e3 && HBScholte(i) < 1e4
%                     String = append(String,newline,'BScholte',num2str(i+1),'    ',num2str(HBScholte(i),'%.3f'));
%                 elseif HBScholte(i) >= 1e4
%                     String = append(String,newline,'BScholte',num2str(i+1),'   ',num2str(HBScholte(i),'%.3f'));
%                 end
%             else
%                 if  HBScholte(i) < 1e2
%                     String = append(String,newline,'BScholte',num2str(i+1),'     ',num2str(HBScholte(i),'%.3f'));
%                 elseif HBScholte(i) >= 1e2 && HBScholte(i) < 1e3
%                     String = append(String,newline,'BScholte',num2str(i+1),'    ',num2str(HBScholte(i),'%.3f'));
%                 elseif HBScholte(i) >= 1e3 && HBScholte(i) < 1e4
%                     String = append(String,newline,'BScholte',num2str(i+1),'   ',num2str(HBScholte(i),'%.3f'));
%                 elseif HBScholte(i) >= 1e4
%                     String = append(String,newline,'BScholte',num2str(i+1),'  ',num2str(HBScholte(i),'%.3f'));
%                 end
%             end
%         end
%     end
%     String = append(String,newline,newline);
%     if  any(HBScholte)
%         String = append(String,'BScholte: ',num2str(length(HBScholte)));
%     end
end
disp([String,newline,'----------------'])