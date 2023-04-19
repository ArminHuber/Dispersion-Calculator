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
function [HS,HSLamb,HSShear,HA,HALamb,HAShear,HB,HBLamb,HBShear] = FrequencySweeper_Anisotropic(c,Delta,Material,SweepRange,PhaseVelocityLimit,Repetitions,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled,OutputWindow1aUI2,OutputWindow1bUI2,OutputWindow2aUI2,OutputWindow2bUI2)
Resolution = 1e-5; % (kHz)

%#ok<*AGROW>
%#ok<*MINV>
HS = [];
HSLamb = [];
HSShear = [];
HA = [];
HALamb = [];
HAShear = [];
HB = [];
HBLamb = [];
HBShear = [];
if  length(SweepRange) < 2
    return
end
I = [1 1 -1;-1 -1 1;-1 -1 1];
I1 = [1 -1;-1 1];
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
else
    Bisections = 1;
end
for m = 1:SuperLayerSize
    c{m} = real(c{m});
    if  ~Decoupled
        A1 = (2*(c{m}(1,3)*c{m}(4,5)^2+c{m}(1,3)*c{m}(3,6)*c{m}(4,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)-c{m}(1,3)*c{m}(4,4)*c{m}(5,5))+c{m}(1,1)*c{m}(3,3)*c{m}(4,4)-c{m}(1,3)^2*c{m}(4,4)-c{m}(3,6)^2*c{m}(5,5)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-(c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)+c{m}(4,4)*c{m}(5,5)-c{m}(4,5)^2)*Material{m}.Density*PhaseVelocityLimit^2)/Delta(m);
        A2 = (2*(c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)+c{m}(1,6)*c{m}(3,6)*c{m}(5,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5))+c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,6)^2*c{m}(3,3)-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)+(2*(c{m}(1,6)*c{m}(4,5)+c{m}(1,3)*c{m}(5,5)+c{m}(3,6)*c{m}(4,5))+c{m}(1,3)^2+c{m}(3,6)^2+c{m}(4,5)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)-c{m}(5,5)*c{m}(6,6))*Material{m}.Density*PhaseVelocityLimit^2+(c{m}(3,3)+c{m}(4,4)+c{m}(5,5))*Material{m}.Density^2*PhaseVelocityLimit^4)/Delta(m);
        A3 = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5)-(c{m}(1,1)*c{m}(5,5)+c{m}(1,1)*c{m}(6,6)+c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2)*Material{m}.Density*PhaseVelocityLimit^2+(c{m}(1,1)+c{m}(5,5)+c{m}(6,6))*Material{m}.Density^2*PhaseVelocityLimit^4-Material{m}.Density^3*PhaseVelocityLimit^6)/Delta(m);
        Alphaa = A2/3-A1^2/9;
        Alphab = A1^3/27-A1*A2/6+A3/2;
        Alphac = (sqrt(Alphab^2+Alphaa^3)-Alphab)^(1/3);
        Alphad = Alphaa/(2*Alphac)-Alphac/2;
        Alphae = Alphaa/Alphac;
        Alphaf = (sqrt(3)*(Alphac+Alphae)*1i)/2;
        Alpha(m,1) = sqrt(Alphad-Alphaf-A1/3);
        Alpha(m,2) = sqrt(Alphad+Alphaf-A1/3);
        Alpha(m,3) = -sqrt(Alphac-Alphae-A1/3);
        m11 = c{m}(1,1)-Material{m}.Density*PhaseVelocityLimit^2+c{m}(5,5)*Alpha(m,:).^2;
        m12 = c{m}(1,6)+c{m}(4,5)*Alpha(m,:).^2;
        m13 = (c{m}(1,3)+c{m}(5,5))*Alpha(m,:);
        m22 = c{m}(6,6)-Material{m}.Density*PhaseVelocityLimit^2+c{m}(4,4)*Alpha(m,:).^2;
        m23 = (c{m}(3,6)+c{m}(4,5))*Alpha(m,:);
        V(m,:) = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
        W(m,:) = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
        D3(m,:) = c{m}(1,3)+c{m}(3,6)*V(m,:)+c{m}(3,3)*Alpha(m,:).*W(m,:);
        D4(m,:) = c{m}(4,5)*(Alpha(m,:)+W(m,:))+c{m}(4,4)*Alpha(m,:).*V(m,:);
        D5(m,:) = c{m}(5,5)*(Alpha(m,:)+W(m,:))+c{m}(4,5)*Alpha(m,:).*V(m,:);
    else
        A1 = 2*c{m}(3,3)*c{m}(5,5);
        A2 = c{m}(5,5)*(c{m}(5,5)-Material{m}.Density*PhaseVelocityLimit^2)+c{m}(3,3)*(c{m}(1,1)-Material{m}.Density*PhaseVelocityLimit^2)-(c{m}(1,3)+c{m}(5,5))^2;
        A3 = (c{m}(5,5)-Material{m}.Density*PhaseVelocityLimit^2)*(c{m}(1,1)-Material{m}.Density*PhaseVelocityLimit^2);
        Alpha(m,1) = sqrt((-A2+sqrt(A2^2-2*A1*A3))/A1);
        Alpha(m,2) = sqrt((-A2-sqrt(A2^2-2*A1*A3))/A1);
        W(m,:) = (Material{m}.Density*PhaseVelocityLimit^2-c{m}(1,1)-c{m}(5,5)*Alpha(m,:).^2)./((c{m}(1,3)+c{m}(5,5))*Alpha(m,:));
        D3(m,:) = c{m}(1,3)+c{m}(3,3)*Alpha(m,:).*W(m,:);
        D5(m,:) = c{m}(5,5)*(Alpha(m,:)+W(m,:));
        AlphaSH(m) = sqrt((Material{m}.Density*PhaseVelocityLimit^2-c{m}(6,6))/c{m}(4,4));
        D(m) = AlphaSH(m)*c{m}(4,4);
    end
end
if  SuperLayerSize == 1 || SymmetricSystem
    if  ~Decoupled
        XRough = [];
        XFine = [];
        for i = 1:length(SweepRange)
            WaveNumber = 2*pi*SweepRange(i)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
%                 L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
%                 L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,3);L1(1,4)=D3(m,1)*E(1);L1(1,5)=D3(m,2)*E(2);L1(1,6)=D3(m,3)*E(3);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=D5(m,3);L1(2,4)=-D5(m,1)*E(1);L1(2,5)=-D5(m,2)*E(2);L1(2,6)=-D5(m,3)*E(3);L1(3,1)=D4(m,1);L1(3,2)=D4(m,2);L1(3,3)=D4(m,3);L1(3,4)=-D4(m,1)*E(1);L1(3,5)=-D4(m,2)*E(2);L1(3,6)=-D4(m,3)*E(3);L1(4,1)=D3(m,1)*E(1);L1(4,2)=D3(m,2)*E(2);L1(4,3)=D3(m,3)*E(3);L1(4,4)=D3(m,1);L1(4,5)=D3(m,2);L1(4,6)=D3(m,3);L1(5,1)=D5(m,1)*E(1);L1(5,2)=D5(m,2)*E(2);L1(5,3)=D5(m,3)*E(3);L1(5,4)=-D5(m,1);L1(5,5)=-D5(m,2);L1(5,6)=-D5(m,3);L1(6,1)=D4(m,1)*E(1);L1(6,2)=D4(m,2)*E(2);L1(6,3)=D4(m,3)*E(3);L1(6,4)=-D4(m,1);L1(6,5)=-D4(m,2);L1(6,6)=-D4(m,3);
                L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(m,1);L2(2,2)=V(m,2);L2(2,3)=V(m,3);L2(2,4)=V(m,1)*E(1);L2(2,5)=V(m,2)*E(2);L2(2,6)=V(m,3)*E(3);L2(3,1)=W(m,1);L2(3,2)=W(m,2);L2(3,3)=W(m,3);L2(3,4)=-W(m,1)*E(1);L2(3,5)=-W(m,2)*E(2);L2(3,6)=-W(m,3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(m,1)*E(1);L2(5,2)=V(m,2)*E(2);L2(5,3)=V(m,3)*E(3);L2(5,4)=V(m,1);L2(5,5)=V(m,2);L2(5,6)=V(m,3);L2(6,1)=W(m,1)*E(1);L2(6,2)=W(m,2)*E(2);L2(6,3)=W(m,3)*E(3);L2(6,4)=-W(m,1);L2(6,5)=-W(m,2);L2(6,6)=-W(m,3);
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
            Y(i) = det(MM{end});
            if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i)))
                XRough(end+1) = SweepRange(i-1);
            end
        end
% figure('name','Lamb'),plot(SweepRange,real(Y),SweepRange,abs(real(Y))),yline(0)
        if  isempty(XRough)
            String = 'No higher order modes found!';
            OutputWindow1aUI2.String = String;
            OutputWindow1bUI2.String = '';
            OutputWindow2aUI2.String = '';
            OutputWindow2bUI2.String = '';
            disp(String)
            return
        end
        for p = 1:length(XRough)
            Frequency = [XRough(p)-(SweepRange(2)-SweepRange(1)) XRough(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                for i = 1:length(Frequency)
                    WaveNumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
%                         L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
%                         L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                        L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,3);L1(1,4)=D3(m,1)*E(1);L1(1,5)=D3(m,2)*E(2);L1(1,6)=D3(m,3)*E(3);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=D5(m,3);L1(2,4)=-D5(m,1)*E(1);L1(2,5)=-D5(m,2)*E(2);L1(2,6)=-D5(m,3)*E(3);L1(3,1)=D4(m,1);L1(3,2)=D4(m,2);L1(3,3)=D4(m,3);L1(3,4)=-D4(m,1)*E(1);L1(3,5)=-D4(m,2)*E(2);L1(3,6)=-D4(m,3)*E(3);L1(4,1)=D3(m,1)*E(1);L1(4,2)=D3(m,2)*E(2);L1(4,3)=D3(m,3)*E(3);L1(4,4)=D3(m,1);L1(4,5)=D3(m,2);L1(4,6)=D3(m,3);L1(5,1)=D5(m,1)*E(1);L1(5,2)=D5(m,2)*E(2);L1(5,3)=D5(m,3)*E(3);L1(5,4)=-D5(m,1);L1(5,5)=-D5(m,2);L1(5,6)=-D5(m,3);L1(6,1)=D4(m,1)*E(1);L1(6,2)=D4(m,2)*E(2);L1(6,3)=D4(m,3)*E(3);L1(6,4)=-D4(m,1);L1(6,5)=-D4(m,2);L1(6,6)=-D4(m,3);
                        L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(m,1);L2(2,2)=V(m,2);L2(2,3)=V(m,3);L2(2,4)=V(m,1)*E(1);L2(2,5)=V(m,2)*E(2);L2(2,6)=V(m,3)*E(3);L2(3,1)=W(m,1);L2(3,2)=W(m,2);L2(3,3)=W(m,3);L2(3,4)=-W(m,1)*E(1);L2(3,5)=-W(m,2)*E(2);L2(3,6)=-W(m,3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(m,1)*E(1);L2(5,2)=V(m,2)*E(2);L2(5,3)=V(m,3)*E(3);L2(5,4)=V(m,1);L2(5,5)=V(m,2);L2(5,6)=V(m,3);L2(6,1)=W(m,1)*E(1);L2(6,2)=W(m,2)*E(2);L2(6,3)=W(m,3)*E(3);L2(6,4)=-W(m,1);L2(6,5)=-W(m,2);L2(6,6)=-W(m,3);
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
                    Y(i) = det(MM{end});
                    if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i)))
                        if  o == Bisections && sign(real(Y(i-2))) ~= sign(real(Y(i)))
                            XFine(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        for p = 1:length(XFine)
            WaveNumber = 2*pi*XFine(p)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
                L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                if  SuperLayerSize > 1
                    L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                    L{m} = L1/L2;
                end
            end
            if  SuperLayerSize == 1
                U = L1(2:6,2:6)\-L1(2:6,1);
                u1(1) = 1+U(1)+U(2)+U(3)*exp(1i*WaveNumber*Alpha(1)*LayerThicknesses)+U(4)*exp(1i*WaveNumber*Alpha(2)*LayerThicknesses)+U(5)*exp(1i*WaveNumber*Alpha(3)*LayerThicknesses);
                u1(2) = exp(1i*WaveNumber*Alpha(1)*LayerThicknesses)+U(1)*exp(1i*WaveNumber*Alpha(2)*LayerThicknesses)+U(2)*exp(1i*WaveNumber*Alpha(3)*LayerThicknesses)+U(3)+U(4)+U(5);
            else
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
                Z1 = [-MM{end}(1:3,1:3)*[1 1 1;V(1,:);-W(1,:)] -MM{end}(1:3,4:6)*[1 1 1;V(1,:);W(1,:)];MM{end}(4:6,1:3)*[1 1 1;V(1,:);-W(1,:)] MM{end}(4:6,4:6)*[1 1 1;V(1,:);W(1,:)]];
                Z2 = [MM{end}(1:3,1:3)*[1 1 1;V(1,:);W(1,:)];-MM{end}(4:6,1:3)*[1 1 1;V(1,:);W(1,:)]];
                RT = Z1\Z2(:,1);
                u1 = [1+[1 1 1]*RT(1:3) [1 1 1]*RT(4:6)];
            end
            u1 = u1*exp(-1i*angle(u1(1)));
            if  sign(real(u1(1))) == sign(real(u1(2)))
                HS(end+1) = XFine(p);
            else
                HA(end+1) = XFine(p);
            end
        end
        String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
        if  Symmetric
            if  any(HA)
                for i = 1:length(HA)
                    if  i < 10
                        if  HA(i) < 1e2
                            String = append(String,newline,'A',num2str(i),'       ',num2str(HA(i),'%.3f'));
                        elseif HA(i) >= 1e2 && HA(i) < 1e3
                            String = append(String,newline,'A',num2str(i),'      ',num2str(HA(i),'%.3f'));
                        elseif HA(i) >= 1e3 && HA(i) < 1e4
                            String = append(String,newline,'A',num2str(i),'     ',num2str(HA(i),'%.3f'));
                        elseif HA(i) >= 1e4
                            String = append(String,newline,'A',num2str(i),'    ',num2str(HA(i),'%.3f'));
                        end
                    else
                        if  HA(i) < 1e2
                            String = append(String,newline,'A',num2str(i),'      ',num2str(HA(i),'%.3f'));
                        elseif HA(i) >= 1e2 && HA(i) < 1e3
                            String = append(String,newline,'A',num2str(i),'     ',num2str(HA(i),'%.3f'));
                        elseif HA(i) >= 1e3 && HA(i) < 1e4
                            String = append(String,newline,'A',num2str(i),'    ',num2str(HA(i),'%.3f'));
                        elseif HA(i) >= 1e4
                            String = append(String,newline,'A',num2str(i),'   ',num2str(HA(i),'%.3f'));
                        end
                    end
                end
            end
        else
            HB = sort(horzcat(HA,HS));
            if  any(HB)
                for i = 1:length(HB)
                    if  i < 8
                        if  HB(i) < 1e2
                            String = append(String,newline,'B',num2str(i+2),'       ',num2str(HB(i),'%.3f'));
                        elseif HB(i) >= 1e2 && HB(i) < 1e3
                            String = append(String,newline,'B',num2str(i+2),'      ',num2str(HB(i),'%.3f'));
                        elseif HB(i) >= 1e3 && HB(i) < 1e4
                            String = append(String,newline,'B',num2str(i+2),'     ',num2str(HB(i),'%.3f'));
                        elseif HB(i) >= 1e4
                            String = append(String,newline,'B',num2str(i+2),'    ',num2str(HB(i),'%.3f'));
                        end
                    else
                        if  HB(i) < 1e2
                            String = append(String,newline,'B',num2str(i+2),'      ',num2str(HB(i),'%.3f'));
                        elseif HB(i) >= 1e2 && HB(i) < 1e3
                            String = append(String,newline,'B',num2str(i+2),'     ',num2str(HB(i),'%.3f'));
                        elseif HB(i) >= 1e3 && HB(i) < 1e4
                            String = append(String,newline,'B',num2str(i+2),'    ',num2str(HB(i),'%.3f'));
                        elseif HB(i) >= 1e4
                            String = append(String,newline,'B',num2str(i+2),'   ',num2str(HB(i),'%.3f'));
                        end
                    end
                end
            end
        end
        OutputWindow1aUI2.String = String;
        disp(String)
        String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
        if  any(HS) && Symmetric
            for i = 1:length(HS)
                if  i < 9
                    if  HS(i) < 1e2
                        String = append(String,newline,'S',num2str(i+1),'       ',num2str(HS(i),'%.3f'));
                    elseif HS(i) >= 1e2 && HS(i) < 1e3
                        String = append(String,newline,'S',num2str(i+1),'      ',num2str(HS(i),'%.3f'));
                    elseif HS(i) >= 1e3 && HS(i) < 1e4
                        String = append(String,newline,'S',num2str(i+1),'     ',num2str(HS(i),'%.3f'));
                    elseif HS(i) >= 1e4
                        String = append(String,newline,'S',num2str(i+1),'    ',num2str(HS(i),'%.3f'));
                    end
                else
                    if  HS(i) < 1e2
                        String = append(String,newline,'S',num2str(i+1),'      ',num2str(HS(i),'%.3f'));
                    elseif HS(i) >= 1e2 && HS(i) < 1e3
                        String = append(String,newline,'S',num2str(i+1),'     ',num2str(HS(i),'%.3f'));
                    elseif HS(i) >= 1e3 && HS(i) < 1e4
                        String = append(String,newline,'S',num2str(i+1),'    ',num2str(HS(i),'%.3f'));
                    elseif HS(i) >= 1e4
                        String = append(String,newline,'S',num2str(i+1),'   ',num2str(HS(i),'%.3f'));
                    end
                end
            end
        end
        OutputWindow1bUI2.String = String;
        disp(extractAfter(String,')'))
        disp(' ')
        if  Symmetric
            if  any(HA)
                disp(['A: ',num2str(length(HA))])
                OutputWindow2aUI2.String = ['A: ',num2str(length(HA))];
            else
                OutputWindow2aUI2.String = '';
            end
        else
            if  any(HB)
                disp(['B: ',num2str(length(HB))])
                OutputWindow2aUI2.String = ['B: ',num2str(length(HB))];
            else
                OutputWindow2aUI2.String = '';
            end
        end
        if  any(HS) && Symmetric
            disp(['S: ',num2str(length(HS))])
            OutputWindow2bUI2.String = ['S: ',num2str(length(HS))];
        else
            OutputWindow2bUI2.String = '';
        end
    else
        XRough = [];
        XRoughSSH = [];
        XRoughASH = [];
        XFine = [];
        for i = 1:length(SweepRange)
            WaveNumber = 2*pi*SweepRange(i)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
%                 L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
%                 L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,1)*E(1);L1(1,4)=D3(m,2)*E(2);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=-D5(m,1)*E(1);L1(2,4)=-D5(m,2)*E(2);L1(3,1)=D3(m,1)*E(1);L1(3,2)=D3(m,2)*E(2);L1(3,3)=D3(m,1);L1(3,4)=D3(m,2);L1(4,1)=D5(m,1)*E(1);L1(4,2)=D5(m,2)*E(2);L1(4,3)=-D5(m,1);L1(4,4)=-D5(m,2);
                L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(m,1);L2(2,2)=W(m,2);L2(2,3)=-W(m,1)*E(1);L2(2,4)=-W(m,2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(m,1)*E(1);L2(4,2)=W(m,2)*E(2);L2(4,3)=-W(m,1);L2(4,4)=-W(m,2);
                L{m} = L1/L2;
                Gamma = AlphaSH(m)*WaveNumber*LayerThicknesses(m);
                LSH{m} = [cos(Gamma) 1i*sin(Gamma)/D(m);1i*D(m)*sin(Gamma) cos(Gamma)];
            end
            M = L{1};
            MSH = LSH{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                MSH = MSH*LSH{m};
            end
            if  Repetitions > 1
                MSH = MSH^Repetitions;
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
            Y(i) = det(MM{end});
            if  SuperLayerSize == 1
                YSSH(i) = sin(Gamma/2);
                YASH(i) = cos(Gamma/2);
            else
                YSSH(i) = MSH(2,1);
                YASH(i) = MSH(2,2);
            end
            if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i)))
                XRough(end+1) = SweepRange(i-1);
            end
            if  i > 2 && abs(YSSH(i-1)) < abs(YSSH(i-2)) && abs(YSSH(i-1)) < abs(YSSH(i))
                XRoughSSH(end+1) = SweepRange(i-1);
            end
            if  i > 2 && abs(YASH(i-1)) < abs(YASH(i-2)) && abs(YASH(i-1)) < abs(YASH(i))
                XRoughASH(end+1) = SweepRange(i-1);
            end
        end
% figure('name','Lamb'),plot(SweepRange,real(Y),SweepRange,abs(real(Y))),yline(0)
% figure('name','SSH'),plot(SweepRange,abs(YSSH))
% figure('name','ASH'),plot(SweepRange,abs(YASH))
        if  isempty(XRough) && isempty(XRoughSSH) && isempty(XRoughASH)
            String = 'No higher order modes found!';
            OutputWindow1aUI2.String = String;
            OutputWindow1bUI2.String = '';
            OutputWindow2aUI2.String = '';
            OutputWindow2bUI2.String = '';
            disp(String)
            return
        end
        for p = 1:length(XRough)
            Frequency = [XRough(p)-(SweepRange(2)-SweepRange(1)) XRough(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                for i = 1:length(Frequency)
                    WaveNumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
%                         L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
%                         L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                        L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,1)*E(1);L1(1,4)=D3(m,2)*E(2);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=-D5(m,1)*E(1);L1(2,4)=-D5(m,2)*E(2);L1(3,1)=D3(m,1)*E(1);L1(3,2)=D3(m,2)*E(2);L1(3,3)=D3(m,1);L1(3,4)=D3(m,2);L1(4,1)=D5(m,1)*E(1);L1(4,2)=D5(m,2)*E(2);L1(4,3)=-D5(m,1);L1(4,4)=-D5(m,2);
                        L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(m,1);L2(2,2)=W(m,2);L2(2,3)=-W(m,1)*E(1);L2(2,4)=-W(m,2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(m,1)*E(1);L2(4,2)=W(m,2)*E(2);L2(4,3)=-W(m,1);L2(4,4)=-W(m,2);
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
                    Y(i) = det(MM{end});
                    if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i)))
                        if  o == Bisections && sign(real(Y(i-2))) ~= sign(real(Y(i)))
                            XFine(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        for p = 1:length(XFine)
            WaveNumber = 2*pi*XFine(p)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
                L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                if  SuperLayerSize > 1
                    L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                    L{m} = L1/L2;
                end
            end
            if  SuperLayerSize == 1
                U = L1(2:4,2:4)\-L1(2:4,1);
                u1(1) = 1+U(1)+U(2)*exp(1i*WaveNumber*Alpha(1)*LayerThicknesses)+U(3)*exp(1i*WaveNumber*Alpha(2)*LayerThicknesses);
                u1(2) = exp(1i*WaveNumber*Alpha(1)*LayerThicknesses)+U(1)*exp(1i*WaveNumber*Alpha(2)*LayerThicknesses)+U(2)+U(3);
            else
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
                Z1 = [-MM{end}(1:2,1:2)*[1 1;-W(1,:)] -MM{end}(1:2,3:4)*[1 1;W(1,:)];MM{end}(3:4,1:2)*[1 1;-W(1,:)] MM{end}(3:4,3:4)*[1 1;W(1,:)]];
                Z2 = [MM{end}(1:2,1:2)*[1 1;W(1,:)];-MM{end}(3:4,1:2)*[1 1;W(1,:)]];
                RT = Z1\Z2(:,1);
                u1 = [1+[1 1]*RT(1:2) [1 1]*RT(3:4)];
            end
            u1 = u1*exp(-1i*angle(u1(1)));
            if  sign(real(u1(1))) == sign(real(u1(2)))
                HSLamb(end+1) = XFine(p);
            else
                HALamb(end+1) = XFine(p);
            end
        end
        for p = 1:length(XRoughSSH)
            Frequency = [XRoughSSH(p)-(SweepRange(2)-SweepRange(1)) XRoughSSH(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                for i = 1:length(Frequency)
                    WaveNumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        Gamma = AlphaSH(m)*WaveNumber*LayerThicknesses(m);
                        L{m} = [cos(Gamma) 1i*sin(Gamma)/D(m);1i*D(m)*sin(Gamma) cos(Gamma)];
                    end
                    M = L{1};
                    for m = 2:SuperLayerSize
                        M = M*L{m};
                    end
                    if  Repetitions > 1
                        M = M^Repetitions;
                    end
                    Y(i) = M(2,1);
                    if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
                        if  o == Bisections
                            HSShear(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        for p = 1:length(XRoughASH)
            Frequency = [XRoughASH(p)-(SweepRange(2)-SweepRange(1)) XRoughASH(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                for i = 1:length(Frequency)
                    WaveNumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        Gamma = AlphaSH(m)*WaveNumber*LayerThicknesses(m);
                        L{m} = [cos(Gamma) 1i*sin(Gamma)/D(m);1i*D(m)*sin(Gamma) cos(Gamma)];
                    end
                    M = L{1};
                    for m = 2:SuperLayerSize
                        M = M*L{m};
                    end
                    if  Repetitions > 1
                        M = M^Repetitions;
                    end
                    if  SuperLayerSize == 1
                        Y(i) = cos(Gamma/2);
                    else
                        Y(i) = M(2,2);
                    end
                    if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i)) 
                        if  o == Bisections
                            HAShear(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
        if  Symmetric
            if  any(HALamb)
                for i = 1:length(HALamb)
                    if  i < 10
                        if  HALamb(i) < 1e2
                            String = append(String,newline,'A',num2str(i),'       ',num2str(HALamb(i),'%.3f'));
                        elseif HALamb(i) >= 1e2 && HALamb(i) < 1e3
                            String = append(String,newline,'A',num2str(i),'      ',num2str(HALamb(i),'%.3f'));
                        elseif HALamb(i) >= 1e3 && HALamb(i) < 1e4
                            String = append(String,newline,'A',num2str(i),'     ',num2str(HALamb(i),'%.3f'));
                        elseif HALamb(i) >= 1e4
                            String = append(String,newline,'A',num2str(i),'    ',num2str(HALamb(i),'%.3f'));
                        end
                    else
                        if  HALamb(i) < 1e2
                            String = append(String,newline,'A',num2str(i),'      ',num2str(HALamb(i),'%.3f'));
                        elseif HALamb(i) >= 1e2 && HALamb(i) < 1e3
                            String = append(String,newline,'A',num2str(i),'     ',num2str(HALamb(i),'%.3f'));
                        elseif HALamb(i) >= 1e3 && HALamb(i) < 1e4
                            String = append(String,newline,'A',num2str(i),'    ',num2str(HALamb(i),'%.3f'));
                        elseif HALamb(i) >= 1e4
                            String = append(String,newline,'A',num2str(i),'   ',num2str(HALamb(i),'%.3f'));
                        end
                    end
                end
            end
        else
            HBLamb = sort(horzcat(HALamb,HSLamb));
            if  any(HBLamb)
                for i = 1:length(HBLamb)
                    if  i < 9
                        if  HBLamb(i) < 1e2
                            String = append(String,newline,'B',num2str(i+1),'       ',num2str(HBLamb(i),'%.3f'));
                        elseif HBLamb(i) >= 1e2 && HBLamb(i) < 1e3
                            String = append(String,newline,'B',num2str(i+1),'      ',num2str(HBLamb(i),'%.3f'));
                        elseif HBLamb(i) >= 1e3 && HBLamb(i) < 1e4
                            String = append(String,newline,'B',num2str(i+1),'     ',num2str(HBLamb(i),'%.3f'));
                        elseif HBLamb(i) >= 1e4
                            String = append(String,newline,'B',num2str(i+1),'    ',num2str(HBLamb(i),'%.3f'));
                        end
                    else
                        if  HBLamb(i) < 1e2
                            String = append(String,newline,'B',num2str(i+1),'      ',num2str(HBLamb(i),'%.3f'));
                        elseif HBLamb(i) >= 1e2 && HBLamb(i) < 1e3
                            String = append(String,newline,'B',num2str(i+1),'     ',num2str(HBLamb(i),'%.3f'));
                        elseif HBLamb(i) >= 1e3 && HBLamb(i) < 1e4
                            String = append(String,newline,'B',num2str(i+1),'    ',num2str(HBLamb(i),'%.3f'));
                        elseif HBLamb(i) >= 1e4
                            String = append(String,newline,'B',num2str(i+1),'   ',num2str(HBLamb(i),'%.3f'));
                        end
                    end
                end
            end
        end
        String = append(String,newline,' ');
        if  any(HAShear)
            for i = 1:length(HAShear)
                if  i < 10
                    if  HAShear(i) < 1e2
                        String = append(String,newline,'ASH',num2str(i),'     ',num2str(HAShear(i),'%.3f'));
                    elseif HAShear(i) >= 1e2 && HAShear(i) < 1e3
                        String = append(String,newline,'ASH',num2str(i),'    ',num2str(HAShear(i),'%.3f'));
                    elseif HAShear(i) >= 1e3 && HAShear(i) < 1e4
                        String = append(String,newline,'ASH',num2str(i),'   ',num2str(HAShear(i),'%.3f'));
                    elseif HAShear(i) >= 1e4
                        String = append(String,newline,'ASH',num2str(i),'  ',num2str(HAShear(i),'%.3f'));
                    end
                else
                    if  HAShear(i) < 1e2
                        String = append(String,newline,'ASH',num2str(i),'    ',num2str(HAShear(i),'%.3f'));
                    elseif HAShear(i) >= 1e2 && HAShear(i) < 1e3
                        String = append(String,newline,'ASH',num2str(i),'   ',num2str(HAShear(i),'%.3f'));
                    elseif HAShear(i) >= 1e3 && HAShear(i) < 1e4
                        String = append(String,newline,'ASH',num2str(i),'  ',num2str(HAShear(i),'%.3f'));
                    elseif HAShear(i) >= 1e4
                        String = append(String,newline,'ASH',num2str(i),' ',num2str(HAShear(i),'%.3f'));
                    end
                end
            end
        end
        OutputWindow1aUI2.String = String;
        disp(String)
        String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
        if  any(HSLamb) && Symmetric
            for i = 1:length(HSLamb)
                if  i < 10
                    if  HSLamb(i) < 1e2
                        String = append(String,newline,'S',num2str(i),'       ',num2str(HSLamb(i),'%.3f'));
                    elseif HSLamb(i) >= 1e2 && HSLamb(i) < 1e3
                        String = append(String,newline,'S',num2str(i),'      ',num2str(HSLamb(i),'%.3f'));
                    elseif HSLamb(i) >= 1e3 && HSLamb(i) < 1e4
                        String = append(String,newline,'S',num2str(i),'     ',num2str(HSLamb(i),'%.3f'));
                    elseif HSLamb(i) >= 1e4
                        String = append(String,newline,'S',num2str(i),'    ',num2str(HSLamb(i),'%.3f'));
                    end
                else
                    if  HSLamb(i) < 1e2
                        String = append(String,newline,'S',num2str(i),'      ',num2str(HSLamb(i),'%.3f'));
                    elseif HSLamb(i) >= 1e2 && HSLamb(i) < 1e3
                        String = append(String,newline,'S',num2str(i),'     ',num2str(HSLamb(i),'%.3f'));
                    elseif HSLamb(i) >= 1e3 && HSLamb(i) < 1e4
                        String = append(String,newline,'S',num2str(i),'    ',num2str(HSLamb(i),'%.3f'));
                    elseif HSLamb(i) >= 1e4
                        String = append(String,newline,'S',num2str(i),'   ',num2str(HSLamb(i),'%.3f'));
                    end
                end
            end
        end
        String = append(String,newline,' ');
        if  any(HSShear)
            for i = 1:length(HSShear)
                if  i < 10
                    if  HSShear(i) < 1e2
                        String = append(String,newline,'SSH',num2str(i),'     ',num2str(HSShear(i),'%.3f'));
                    elseif HSShear(i) >= 1e2 && HSShear(i) < 1e3
                        String = append(String,newline,'SSH',num2str(i),'    ',num2str(HSShear(i),'%.3f'));
                    elseif HSShear(i) >= 1e3 && HSShear(i) < 1e4
                        String = append(String,newline,'SSH',num2str(i),'   ',num2str(HSShear(i),'%.3f'));
                    elseif HSShear(i) >= 1e4
                        String = append(String,newline,'SSH',num2str(i),'  ',num2str(HSShear(i),'%.3f'));
                    end
                else
                    if  HSShear(i) < 1e2
                        String = append(String,newline,'SSH',num2str(i),'    ',num2str(HSShear(i),'%.3f'));
                    elseif HSShear(i) >= 1e2 && HSShear(i) < 1e3
                        String = append(String,newline,'SSH',num2str(i),'   ',num2str(HSShear(i),'%.3f'));
                    elseif HSShear(i) >= 1e3 && HSShear(i) < 1e4
                        String = append(String,newline,'SSH',num2str(i),'  ',num2str(HSShear(i),'%.3f'));
                    elseif HSShear(i) >= 1e4
                        String = append(String,newline,'SSH',num2str(i),' ',num2str(HSShear(i),'%.3f'));
                    end
                end
            end
        end
        OutputWindow1bUI2.String = String;
        disp(extractAfter(String,')'))
        disp(' ')
        String = '';
        if  Symmetric
            if  any(HALamb)
                String = append(String,'A:   ',num2str(length(HALamb)));
            end
        else
            if  any(HBLamb)
                String = append(String,'B:   ',num2str(length(HBLamb)));
            end            
        end
        if  any(HAShear)
            String = append(String,newline,'ASH: ',num2str(length(HAShear)));
        end
        OutputWindow2aUI2.String = String;
        disp(String)
        String = '';
        if  any(HSLamb) && Symmetric
            String = append(String,'S:   ',num2str(length(HSLamb)));
        end
        if  any(HSShear)
            String = append(String,newline,'SSH: ',num2str(length(HSShear)));
        end
        OutputWindow2bUI2.String = String;
        disp(String)
    end
elseif SuperLayerSize > 1 && ~SymmetricSystem
    if  ~Decoupled
        XRough = [];
        for i = 1:length(SweepRange)
            WaveNumber = 2*pi*SweepRange(i)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
%                 L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
%                 L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,3);L1(1,4)=D3(m,1)*E(1);L1(1,5)=D3(m,2)*E(2);L1(1,6)=D3(m,3)*E(3);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=D5(m,3);L1(2,4)=-D5(m,1)*E(1);L1(2,5)=-D5(m,2)*E(2);L1(2,6)=-D5(m,3)*E(3);L1(3,1)=D4(m,1);L1(3,2)=D4(m,2);L1(3,3)=D4(m,3);L1(3,4)=-D4(m,1)*E(1);L1(3,5)=-D4(m,2)*E(2);L1(3,6)=-D4(m,3)*E(3);L1(4,1)=D3(m,1)*E(1);L1(4,2)=D3(m,2)*E(2);L1(4,3)=D3(m,3)*E(3);L1(4,4)=D3(m,1);L1(4,5)=D3(m,2);L1(4,6)=D3(m,3);L1(5,1)=D5(m,1)*E(1);L1(5,2)=D5(m,2)*E(2);L1(5,3)=D5(m,3)*E(3);L1(5,4)=-D5(m,1);L1(5,5)=-D5(m,2);L1(5,6)=-D5(m,3);L1(6,1)=D4(m,1)*E(1);L1(6,2)=D4(m,2)*E(2);L1(6,3)=D4(m,3)*E(3);L1(6,4)=-D4(m,1);L1(6,5)=-D4(m,2);L1(6,6)=-D4(m,3);
                L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(m,1);L2(2,2)=V(m,2);L2(2,3)=V(m,3);L2(2,4)=V(m,1)*E(1);L2(2,5)=V(m,2)*E(2);L2(2,6)=V(m,3)*E(3);L2(3,1)=W(m,1);L2(3,2)=W(m,2);L2(3,3)=W(m,3);L2(3,4)=-W(m,1)*E(1);L2(3,5)=-W(m,2)*E(2);L2(3,6)=-W(m,3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(m,1)*E(1);L2(5,2)=V(m,2)*E(2);L2(5,3)=V(m,3)*E(3);L2(5,4)=V(m,1);L2(5,5)=V(m,2);L2(5,6)=V(m,3);L2(6,1)=W(m,1)*E(1);L2(6,2)=W(m,2)*E(2);L2(6,3)=W(m,3)*E(3);L2(6,4)=-W(m,1);L2(6,5)=-W(m,2);L2(6,6)=-W(m,3);
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
            Y(i) = det(MM{end});
            if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i)))
                XRough(end+1) = SweepRange(i-1);
            end
        end
% figure('name','Lamb'),plot(SweepRange,real(Y),SweepRange,abs(real(Y))),yline(0)
        if  isempty(XRough)
            String = 'No higher order modes found!';
            OutputWindow1aUI2.String = String;
            OutputWindow1bUI2.String = '';
            OutputWindow2aUI2.String = '';
            OutputWindow2bUI2.String = '';
            disp(String)
            return
        end
        for p = 1:length(XRough)
            Frequency = [XRough(p)-(SweepRange(2)-SweepRange(1)) XRough(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                for i = 1:length(Frequency)
                    WaveNumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
%                         L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
%                         L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                        L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,3);L1(1,4)=D3(m,1)*E(1);L1(1,5)=D3(m,2)*E(2);L1(1,6)=D3(m,3)*E(3);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=D5(m,3);L1(2,4)=-D5(m,1)*E(1);L1(2,5)=-D5(m,2)*E(2);L1(2,6)=-D5(m,3)*E(3);L1(3,1)=D4(m,1);L1(3,2)=D4(m,2);L1(3,3)=D4(m,3);L1(3,4)=-D4(m,1)*E(1);L1(3,5)=-D4(m,2)*E(2);L1(3,6)=-D4(m,3)*E(3);L1(4,1)=D3(m,1)*E(1);L1(4,2)=D3(m,2)*E(2);L1(4,3)=D3(m,3)*E(3);L1(4,4)=D3(m,1);L1(4,5)=D3(m,2);L1(4,6)=D3(m,3);L1(5,1)=D5(m,1)*E(1);L1(5,2)=D5(m,2)*E(2);L1(5,3)=D5(m,3)*E(3);L1(5,4)=-D5(m,1);L1(5,5)=-D5(m,2);L1(5,6)=-D5(m,3);L1(6,1)=D4(m,1)*E(1);L1(6,2)=D4(m,2)*E(2);L1(6,3)=D4(m,3)*E(3);L1(6,4)=-D4(m,1);L1(6,5)=-D4(m,2);L1(6,6)=-D4(m,3);
                        L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(m,1);L2(2,2)=V(m,2);L2(2,3)=V(m,3);L2(2,4)=V(m,1)*E(1);L2(2,5)=V(m,2)*E(2);L2(2,6)=V(m,3)*E(3);L2(3,1)=W(m,1);L2(3,2)=W(m,2);L2(3,3)=W(m,3);L2(3,4)=-W(m,1)*E(1);L2(3,5)=-W(m,2)*E(2);L2(3,6)=-W(m,3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(m,1)*E(1);L2(5,2)=V(m,2)*E(2);L2(5,3)=V(m,3)*E(3);L2(5,4)=V(m,1);L2(5,5)=V(m,2);L2(5,6)=V(m,3);L2(6,1)=W(m,1)*E(1);L2(6,2)=W(m,2)*E(2);L2(6,3)=W(m,3)*E(3);L2(6,4)=-W(m,1);L2(6,5)=-W(m,2);L2(6,6)=-W(m,3);
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
                    Y(i) = det(MM{end});
                    if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i)))
                        if  o == Bisections && sign(real(Y(i-2))) ~= sign(real(Y(i)))
                            HB(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
        if  any(HB)
            for i = 1:length(HB)
                if  i < 8
                    if  HB(i) < 1e2
                        String = append(String,newline,'B',num2str(i+2),'       ',num2str(HB(i),'%.3f'));
                    elseif HB(i) >= 1e2 && HB(i) < 1e3
                        String = append(String,newline,'B',num2str(i+2),'      ',num2str(HB(i),'%.3f'));
                    elseif HB(i) >= 1e3 && HB(i) < 1e4
                        String = append(String,newline,'B',num2str(i+2),'     ',num2str(HB(i),'%.3f'));
                    elseif HB(i) >= 1e4
                        String = append(String,newline,'B',num2str(i+2),'    ',num2str(HB(i),'%.3f'));
                    end
                else
                    if  HB(i) < 1e2
                        String = append(String,newline,'B',num2str(i+2),'      ',num2str(HB(i),'%.3f'));
                    elseif HB(i) >= 1e2 && HB(i) < 1e3
                        String = append(String,newline,'B',num2str(i+2),'     ',num2str(HB(i),'%.3f'));
                    elseif HB(i) >= 1e3 && HB(i) < 1e4
                        String = append(String,newline,'B',num2str(i+2),'    ',num2str(HB(i),'%.3f'));
                    elseif HB(i) >= 1e4
                        String = append(String,newline,'B',num2str(i+2),'   ',num2str(HB(i),'%.3f'));
                    end
                end
            end
        end
        OutputWindow1aUI2.String = String;
        disp(String)
        disp(' ')
        if  any(HB)
            disp(['B: ',num2str(length(HB))])
            OutputWindow2aUI2.String = ['B: ',num2str(length(HB))];
        end
        OutputWindow1bUI2.String = '';
        OutputWindow2bUI2.String = '';
    else
        XRough = [];
        XRoughSH = [];
        for i = 1:length(SweepRange)
            WaveNumber = 2*pi*SweepRange(i)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
%                 L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
%                 L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,1)*E(1);L1(1,4)=D3(m,2)*E(2);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=-D5(m,1)*E(1);L1(2,4)=-D5(m,2)*E(2);L1(3,1)=D3(m,1)*E(1);L1(3,2)=D3(m,2)*E(2);L1(3,3)=D3(m,1);L1(3,4)=D3(m,2);L1(4,1)=D5(m,1)*E(1);L1(4,2)=D5(m,2)*E(2);L1(4,3)=-D5(m,1);L1(4,4)=-D5(m,2);
                L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(m,1);L2(2,2)=W(m,2);L2(2,3)=-W(m,1)*E(1);L2(2,4)=-W(m,2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(m,1)*E(1);L2(4,2)=W(m,2)*E(2);L2(4,3)=-W(m,1);L2(4,4)=-W(m,2);
                L{m} = L1/L2;
                Gamma = AlphaSH(m)*WaveNumber*LayerThicknesses(m);
                LSH{m} = [cos(Gamma) 1i*sin(Gamma)/D(m);1i*D(m)*sin(Gamma) cos(Gamma)];
            end
            M = L{1};
            MSH = LSH{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                MSH = MSH*LSH{m};
            end
            if  Repetitions > 1
                MSH = MSH^Repetitions;
            end
            MM{1} = M;
            for m = 2:Repetitions
                N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
            end
            Y(i) = det(MM{end});
            YSH(i) = MSH(2,1);
            if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i)))    
                XRough(end+1) = SweepRange(i-1);
            end
            if  i > 2 && abs(YSH(i-1)) < abs(YSH(i-2)) && abs(YSH(i-1)) < abs(YSH(i))
                XRoughSH(end+1) = SweepRange(i-1);
            end
        end
% figure('name','Lamb'),plot(SweepRange,real(Y),SweepRange,abs(real(Y))),yline(0)
% figure('name','SH'),plot(SweepRange,abs(YSH))
        if  isempty(XRough) && isempty(XRoughSH)
            String = 'No higher order modes found!';
            OutputWindow1aUI2.String = String;
            OutputWindow1bUI2.String = '';
            OutputWindow2aUI2.String = '';
            OutputWindow2bUI2.String = '';
            disp(String)
            return
        end
        Y = [];
        for p = 1:length(XRough)
            Frequency = [XRough(p)-(SweepRange(2)-SweepRange(1)) XRough(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                for i = 1:length(Frequency)
                    WaveNumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
%                         L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
%                         L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                        L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,1)*E(1);L1(1,4)=D3(m,2)*E(2);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=-D5(m,1)*E(1);L1(2,4)=-D5(m,2)*E(2);L1(3,1)=D3(m,1)*E(1);L1(3,2)=D3(m,2)*E(2);L1(3,3)=D3(m,1);L1(3,4)=D3(m,2);L1(4,1)=D5(m,1)*E(1);L1(4,2)=D5(m,2)*E(2);L1(4,3)=-D5(m,1);L1(4,4)=-D5(m,2);
                        L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(m,1);L2(2,2)=W(m,2);L2(2,3)=-W(m,1)*E(1);L2(2,4)=-W(m,2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(m,1)*E(1);L2(4,2)=W(m,2)*E(2);L2(4,3)=-W(m,1);L2(4,4)=-W(m,2);
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
                    Y(i) = det(MM{end});
                    if  i > 2 && abs(real(Y(i-1))) < abs(real(Y(i-2))) && abs(real(Y(i-1))) < abs(real(Y(i)))
                        if  o == Bisections && sign(real(Y(i-2))) ~= sign(real(Y(i)))
                            HBLamb(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        for p = 1:length(XRoughSH)
            Frequency = [XRoughSH(p)-(SweepRange(2)-SweepRange(1)) XRoughSH(p)+(SweepRange(2)-SweepRange(1))];
            for o = 1:Bisections
                Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
                for i = 1:length(Frequency)
                    WaveNumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        Gamma = AlphaSH(m)*WaveNumber*LayerThicknesses(m);
                        L{m} = [cos(Gamma) 1i*sin(Gamma)/D(m);1i*D(m)*sin(Gamma) cos(Gamma)];
                    end
                    M = L{1};
                    for m = 2:SuperLayerSize
                        M = M*L{m};
                    end
                    if  Repetitions > 1
                        M = M^Repetitions;
                    end
                    Y(i) = M(2,1);
                    if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
                        if  o == Bisections
                            HBShear(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
        if  any(HBLamb)
            for i = 1:length(HBLamb)
                if  i < 9
                    if  HBLamb(i) < 1e2
                        String = append(String,newline,'B',num2str(i+1),'       ',num2str(HBLamb(i),'%.3f'));
                    elseif HBLamb(i) >= 1e2 && HBLamb(i) < 1e3
                        String = append(String,newline,'B',num2str(i+1),'      ',num2str(HBLamb(i),'%.3f'));
                    elseif HBLamb(i) >= 1e3 && HBLamb(i) < 1e4
                        String = append(String,newline,'B',num2str(i+1),'     ',num2str(HBLamb(i),'%.3f'));
                    elseif HBLamb(i) >= 1e4
                        String = append(String,newline,'B',num2str(i+1),'    ',num2str(HBLamb(i),'%.3f'));
                    end
                else
                    if  HBLamb(i) < 1e2
                        String = append(String,newline,'B',num2str(i+1),'      ',num2str(HBLamb(i),'%.3f'));
                    elseif HBLamb(i) >= 1e2 && HBLamb(i) < 1e3
                        String = append(String,newline,'B',num2str(i+1),'     ',num2str(HBLamb(i),'%.3f'));
                    elseif HBLamb(i) >= 1e3 && HBLamb(i) < 1e4
                        String = append(String,newline,'B',num2str(i+1),'    ',num2str(HBLamb(i),'%.3f'));
                    elseif HBLamb(i) >= 1e4
                        String = append(String,newline,'B',num2str(i+1),'   ',num2str(HBLamb(i),'%.3f'));
                    end
                end
            end
        end
        OutputWindow1aUI2.String = String;
        disp(String)
        String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
        if  any(HBShear)
            for i = 1:length(HBShear)
                if  i < 10
                    if  HBShear(i) < 1e2
                        String = append(String,newline,'BSH',num2str(i),'     ',num2str(HBShear(i),'%.3f'));
                    elseif HBShear(i) >= 1e2 && HBShear(i) < 1e3
                        String = append(String,newline,'BSH',num2str(i),'    ',num2str(HBShear(i),'%.3f'));
                    elseif HBShear(i) >= 1e3 && HBShear(i) < 1e4
                        String = append(String,newline,'BSH',num2str(i),'   ',num2str(HBShear(i),'%.3f'));
                    elseif HBShear(i) >= 1e4
                        String = append(String,newline,'BSH',num2str(i),'  ',num2str(HBShear(i),'%.3f'));
                    end
                else
                    if  HBShear(i) < 1e2
                        String = append(String,newline,'BSH',num2str(i),'    ',num2str(HBShear(i),'%.3f'));
                    elseif HBShear(i) >= 1e2 && HBShear(i) < 1e3
                        String = append(String,newline,'BSH',num2str(i),'   ',num2str(HBShear(i),'%.3f'));
                    elseif HBShear(i) >= 1e3 && HBShear(i) < 1e4
                        String = append(String,newline,'BSH',num2str(i),'  ',num2str(HBShear(i),'%.3f'));
                    elseif HBShear(i) >= 1e4
                        String = append(String,newline,'BSH',num2str(i),' ',num2str(HBShear(i),'%.3f'));
                    end
                end
            end
        end
        OutputWindow1bUI2.String = String;
        disp(extractAfter(String,')'))
        disp(' ')
        if  any(HBLamb)
            disp(['B:   ',num2str(length(HBLamb))])
            OutputWindow2aUI2.String = ['B:   ',num2str(length(HBLamb))];
        else
            OutputWindow2aUI2.String = '';
        end
        if  any(HBShear)
            disp(['BSH: ',num2str(length(HBShear))])
            OutputWindow2bUI2.String = ['BSH: ',num2str(length(HBShear))];
        else
            OutputWindow2bUI2.String = '';
        end
    end
end
disp('----------------')