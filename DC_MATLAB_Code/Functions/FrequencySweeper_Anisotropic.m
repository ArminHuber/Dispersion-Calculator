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
function [HS,HSLamb,HSShear,HA,HALamb,HAShear,HB,HBLamb,HBShear] = FrequencySweeper_Anisotropic(c,Delta,Material,SweepRange,PhaseVelocityLimit,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled,OutputWindow1aUI2,OutputWindow1bUI2,OutputWindow2aUI2,OutputWindow2bUI2)
TMM = true; % use TMM instead of SMM when a minimum number of bulk waves is propagating (not evanescent)
TMMcondition_Coupled = 1*3*SuperLayerSize; % minimum number of pairs of bulk waves which must be propagating in a superlayer
TMMcondition_Decoupled = 1*2*SuperLayerSize;
Resolution = 1e-5; % (kHz)

%#ok<*AGROW>
%#ok<*MINV>
%#ok<*UNRCH>
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
        d1 = A2/3-A1^2/9;
        d2 = A1^3/27-A1*A2/6+A3/2;
        d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
        d4 = d1/(2*d3)-d3/2;
        d5 = d1/d3;
        d6 = (sqrt(3)*(d3+d5)*1i)/2;
        Alpha(m,1) = sqrt(d4-d6-A1/3);
        Alpha(m,2) = sqrt(d4+d6-A1/3);
        Alpha(m,3) = -sqrt(d3-d5-A1/3);
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
        Propagating(3*m-2:3*m) = abs(imag(Alpha(m,:))) < 1e-10;
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
        DSH(m) = AlphaSH(m)*c{m}(4,4);
        Propagating(2*m-1:2*m) = abs(imag(Alpha(m,:))) < 1e-10;
    end
end
if  TMM && ((~Decoupled && numel(find(Propagating)) >= TMMcondition_Coupled) || (Decoupled && numel(find(Propagating)) >= TMMcondition_Decoupled))
    UseTMM = true;
else
    UseTMM = false;
end
if  SuperLayerSize == 1 || SymmetricSystem
    if  ~Decoupled
        XRough = [];
        XFine = [];
        for i = 1:length(SweepRange)
            Wavenumber = 2*pi*SweepRange(i)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                if  UseTMM
                    D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                    % R = [1 1 1 1 1 1;V(m,:) V(m,:);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:);D4(m,:) -D4(m,:)];
                    R(1,1)=1;R(1,2)=1;R(1,3)=1;R(1,4)=1;R(1,5)=1;R(1,6)=1;R(2,1)=V(m,1);R(2,2)=V(m,2);R(2,3)=V(m,3);R(2,4)=V(m,1);R(2,5)=V(m,2);R(2,6)=V(m,3);R(3,1)=W(m,1);R(3,2)=W(m,2);R(3,3)=W(m,3);R(3,4)=-W(m,1);R(3,5)=-W(m,2);R(3,6)=-W(m,3);R(4,1)=D3(m,1);R(4,2)=D3(m,2);R(4,3)=D3(m,3);R(4,4)=D3(m,1);R(4,5)=D3(m,2);R(4,6)=D3(m,3);R(5,1)=D5(m,1);R(5,2)=D5(m,2);R(5,3)=D5(m,3);R(5,4)=-D5(m,1);R(5,5)=-D5(m,2);R(5,6)=-D5(m,3);R(6,1)=D4(m,1);R(6,2)=D4(m,2);R(6,3)=D4(m,3);R(6,4)=-D4(m,1);R(6,5)=-D4(m,2);R(6,6)=-D4(m,3);
                    L{m} = R*D/R;
                else
                    E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                    % L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                    % L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                    L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,3);L1(1,4)=D3(m,1)*E(1);L1(1,5)=D3(m,2)*E(2);L1(1,6)=D3(m,3)*E(3);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=D5(m,3);L1(2,4)=-D5(m,1)*E(1);L1(2,5)=-D5(m,2)*E(2);L1(2,6)=-D5(m,3)*E(3);L1(3,1)=D4(m,1);L1(3,2)=D4(m,2);L1(3,3)=D4(m,3);L1(3,4)=-D4(m,1)*E(1);L1(3,5)=-D4(m,2)*E(2);L1(3,6)=-D4(m,3)*E(3);L1(4,1)=D3(m,1)*E(1);L1(4,2)=D3(m,2)*E(2);L1(4,3)=D3(m,3)*E(3);L1(4,4)=D3(m,1);L1(4,5)=D3(m,2);L1(4,6)=D3(m,3);L1(5,1)=D5(m,1)*E(1);L1(5,2)=D5(m,2)*E(2);L1(5,3)=D5(m,3)*E(3);L1(5,4)=-D5(m,1);L1(5,5)=-D5(m,2);L1(5,6)=-D5(m,3);L1(6,1)=D4(m,1)*E(1);L1(6,2)=D4(m,2)*E(2);L1(6,3)=D4(m,3)*E(3);L1(6,4)=-D4(m,1);L1(6,5)=-D4(m,2);L1(6,6)=-D4(m,3);
                    L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(m,1);L2(2,2)=V(m,2);L2(2,3)=V(m,3);L2(2,4)=V(m,1)*E(1);L2(2,5)=V(m,2)*E(2);L2(2,6)=V(m,3)*E(3);L2(3,1)=W(m,1);L2(3,2)=W(m,2);L2(3,3)=W(m,3);L2(3,4)=-W(m,1)*E(1);L2(3,5)=-W(m,2)*E(2);L2(3,6)=-W(m,3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(m,1)*E(1);L2(5,2)=V(m,2)*E(2);L2(5,3)=V(m,3)*E(3);L2(5,4)=V(m,1);L2(5,5)=V(m,2);L2(5,6)=V(m,3);L2(6,1)=W(m,1)*E(1);L2(6,2)=W(m,2)*E(2);L2(6,3)=W(m,3)*E(3);L2(6,4)=-W(m,1);L2(6,5)=-W(m,2);L2(6,6)=-W(m,3);
                    L{m} = L1/L2;
                end
            end
            M = L{1};
            if  UseTMM
                for m = 2:SuperLayerSize
                    M = M*L{m};
                end
                if  Repetitions > 1
                    M = M^Repetitions;
                end
                if  SymmetricSystem
                    M2 = L{end};
                    for m = SuperLayerSize-1:-1:1
                        M2 = M2*L{m};
                    end
                    M = M*M2^Repetitions;
                end
                Y(i) = imag(det(M(4:6,1:3)));
            else
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                    M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
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
                Y(i) = real(det(MM{end}));
            end
            if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
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
                    Wavenumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        if  UseTMM
                            D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                            % R = [1 1 1 1 1 1;V(m,:) V(m,:);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:);D4(m,:) -D4(m,:)];
                            R(1,1)=1;R(1,2)=1;R(1,3)=1;R(1,4)=1;R(1,5)=1;R(1,6)=1;R(2,1)=V(m,1);R(2,2)=V(m,2);R(2,3)=V(m,3);R(2,4)=V(m,1);R(2,5)=V(m,2);R(2,6)=V(m,3);R(3,1)=W(m,1);R(3,2)=W(m,2);R(3,3)=W(m,3);R(3,4)=-W(m,1);R(3,5)=-W(m,2);R(3,6)=-W(m,3);R(4,1)=D3(m,1);R(4,2)=D3(m,2);R(4,3)=D3(m,3);R(4,4)=D3(m,1);R(4,5)=D3(m,2);R(4,6)=D3(m,3);R(5,1)=D5(m,1);R(5,2)=D5(m,2);R(5,3)=D5(m,3);R(5,4)=-D5(m,1);R(5,5)=-D5(m,2);R(5,6)=-D5(m,3);R(6,1)=D4(m,1);R(6,2)=D4(m,2);R(6,3)=D4(m,3);R(6,4)=-D4(m,1);R(6,5)=-D4(m,2);R(6,6)=-D4(m,3);
                            L{m} = R*D/R;
                        else
                            E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                            % L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                            % L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                            L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,3);L1(1,4)=D3(m,1)*E(1);L1(1,5)=D3(m,2)*E(2);L1(1,6)=D3(m,3)*E(3);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=D5(m,3);L1(2,4)=-D5(m,1)*E(1);L1(2,5)=-D5(m,2)*E(2);L1(2,6)=-D5(m,3)*E(3);L1(3,1)=D4(m,1);L1(3,2)=D4(m,2);L1(3,3)=D4(m,3);L1(3,4)=-D4(m,1)*E(1);L1(3,5)=-D4(m,2)*E(2);L1(3,6)=-D4(m,3)*E(3);L1(4,1)=D3(m,1)*E(1);L1(4,2)=D3(m,2)*E(2);L1(4,3)=D3(m,3)*E(3);L1(4,4)=D3(m,1);L1(4,5)=D3(m,2);L1(4,6)=D3(m,3);L1(5,1)=D5(m,1)*E(1);L1(5,2)=D5(m,2)*E(2);L1(5,3)=D5(m,3)*E(3);L1(5,4)=-D5(m,1);L1(5,5)=-D5(m,2);L1(5,6)=-D5(m,3);L1(6,1)=D4(m,1)*E(1);L1(6,2)=D4(m,2)*E(2);L1(6,3)=D4(m,3)*E(3);L1(6,4)=-D4(m,1);L1(6,5)=-D4(m,2);L1(6,6)=-D4(m,3);
                            L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(m,1);L2(2,2)=V(m,2);L2(2,3)=V(m,3);L2(2,4)=V(m,1)*E(1);L2(2,5)=V(m,2)*E(2);L2(2,6)=V(m,3)*E(3);L2(3,1)=W(m,1);L2(3,2)=W(m,2);L2(3,3)=W(m,3);L2(3,4)=-W(m,1)*E(1);L2(3,5)=-W(m,2)*E(2);L2(3,6)=-W(m,3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(m,1)*E(1);L2(5,2)=V(m,2)*E(2);L2(5,3)=V(m,3)*E(3);L2(5,4)=V(m,1);L2(5,5)=V(m,2);L2(5,6)=V(m,3);L2(6,1)=W(m,1)*E(1);L2(6,2)=W(m,2)*E(2);L2(6,3)=W(m,3)*E(3);L2(6,4)=-W(m,1);L2(6,5)=-W(m,2);L2(6,6)=-W(m,3);
                            L{m} = L1/L2;
                        end
                    end
                    M = L{1};
                    if  UseTMM
                        for m = 2:SuperLayerSize
                            M = M*L{m};
                        end
                        if  Repetitions > 1
                            M = M^Repetitions;
                        end
                        if  SymmetricSystem
                            M2 = L{end};
                            for m = SuperLayerSize-1:-1:1
                                M2 = M2*L{m};
                            end
                            M = M*M2^Repetitions;
                        end
                        Y(i) = imag(det(M(4:6,1:3)));
                    else
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                            M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
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
                        Y(i) = real(det(MM{end}));
                    end
                    if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
                        if  o == Bisections && sign(Y(i-2)) ~= sign(Y(i))
                            XFine(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        for p = 1:length(XFine)
            Wavenumber = 2*pi*XFine(p)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                if  SuperLayerSize > 1
                    L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                    L{m} = L1/L2;
                end
            end
            if  SuperLayerSize == 1
                U = L1(2:6,2:6)\-L1(2:6,1);
                u1(1) = 1+U(1)+U(2)+U(3)*exp(1i*Wavenumber*Alpha(1)*LayerThicknesses)+U(4)*exp(1i*Wavenumber*Alpha(2)*LayerThicknesses)+U(5)*exp(1i*Wavenumber*Alpha(3)*LayerThicknesses);
                u1(2) = exp(1i*Wavenumber*Alpha(1)*LayerThicknesses)+U(1)*exp(1i*Wavenumber*Alpha(2)*LayerThicknesses)+U(2)*exp(1i*Wavenumber*Alpha(3)*LayerThicknesses)+U(3)+U(4)+U(5);
            else
                M = L{1};
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                    M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
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
            Wavenumber = 2*pi*SweepRange(i)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                if  UseTMM
                    D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                    % R = [1 1 1 1;W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:)];
                    R(1,1)=1;R(1,2)=1;R(1,3)=1;R(1,4)=1;R(2,1)=W(m,1);R(2,2)=W(m,2);R(2,3)=-W(m,1);R(2,4)=-W(m,2);R(3,1)=D3(m,1);R(3,2)=D3(m,2);R(3,3)=D3(m,1);R(3,4)=D3(m,2);R(4,1)=D5(m,1);R(4,2)=D5(m,2);R(4,3)=-D5(m,1);R(4,4)=-D5(m,2);
                    L{m} = R*D/R;
                else
                    E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                    % L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                    % L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                    L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,1)*E(1);L1(1,4)=D3(m,2)*E(2);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=-D5(m,1)*E(1);L1(2,4)=-D5(m,2)*E(2);L1(3,1)=D3(m,1)*E(1);L1(3,2)=D3(m,2)*E(2);L1(3,3)=D3(m,1);L1(3,4)=D3(m,2);L1(4,1)=D5(m,1)*E(1);L1(4,2)=D5(m,2)*E(2);L1(4,3)=-D5(m,1);L1(4,4)=-D5(m,2);
                    L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(m,1);L2(2,2)=W(m,2);L2(2,3)=-W(m,1)*E(1);L2(2,4)=-W(m,2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(m,1)*E(1);L2(4,2)=W(m,2)*E(2);L2(4,3)=-W(m,1);L2(4,4)=-W(m,2);
                    L{m} = L1/L2;
                end
                Gamma = AlphaSH(m)*Wavenumber*LayerThicknesses(m);
                LSH{m} = [cos(Gamma) 1i*sin(Gamma)/DSH(m);1i*DSH(m)*sin(Gamma) cos(Gamma)];
            end
            M = L{1};
            MSH = LSH{1};
            if  UseTMM
                for m = 2:SuperLayerSize
                    M = M*L{m};
                    MSH = MSH*LSH{m};
                end
                if  Repetitions > 1
                    M = M^Repetitions;
                    MSH = MSH^Repetitions;
                end
                if  SymmetricSystem
                    M2 = L{end};
                    for m = SuperLayerSize-1:-1:1
                        M2 = M2*L{m};
                    end
                    M = M*M2^Repetitions;
                end
                Y(i) = real(det(M(3:4,1:2)));
            else
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                    M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                    MSH = MSH*LSH{m};
                end
                if  Repetitions > 1
                    MSH = MSH^Repetitions;
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
                Y(i) = real(det(MM{end}));
            end
            if  SuperLayerSize == 1
                YSSH(i) = sin(Gamma/2);
                YASH(i) = cos(Gamma/2);
            else
                YSSH(i) = MSH(2,1);
                YASH(i) = MSH(2,2);
            end
            if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
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
                    Wavenumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        if  UseTMM
                            D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                            % R = [1 1 1 1;W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:)];
                            R(1,1)=1;R(1,2)=1;R(1,3)=1;R(1,4)=1;R(2,1)=W(m,1);R(2,2)=W(m,2);R(2,3)=-W(m,1);R(2,4)=-W(m,2);R(3,1)=D3(m,1);R(3,2)=D3(m,2);R(3,3)=D3(m,1);R(3,4)=D3(m,2);R(4,1)=D5(m,1);R(4,2)=D5(m,2);R(4,3)=-D5(m,1);R(4,4)=-D5(m,2);
                            L{m} = R*D/R;
                        else
                            E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                            % L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                            % L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                            L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,1)*E(1);L1(1,4)=D3(m,2)*E(2);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=-D5(m,1)*E(1);L1(2,4)=-D5(m,2)*E(2);L1(3,1)=D3(m,1)*E(1);L1(3,2)=D3(m,2)*E(2);L1(3,3)=D3(m,1);L1(3,4)=D3(m,2);L1(4,1)=D5(m,1)*E(1);L1(4,2)=D5(m,2)*E(2);L1(4,3)=-D5(m,1);L1(4,4)=-D5(m,2);
                            L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(m,1);L2(2,2)=W(m,2);L2(2,3)=-W(m,1)*E(1);L2(2,4)=-W(m,2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(m,1)*E(1);L2(4,2)=W(m,2)*E(2);L2(4,3)=-W(m,1);L2(4,4)=-W(m,2);
                            L{m} = L1/L2;
                        end
                    end
                    M = L{1};                    
                    if  UseTMM
                        for m = 2:SuperLayerSize
                            M = M*L{m};
                        end
                        if  Repetitions > 1
                            M = M^Repetitions;
                        end
                        if  SymmetricSystem
                            M2 = L{end};
                            for m = SuperLayerSize-1:-1:1
                                M2 = M2*L{m};
                            end
                            M = M*M2^Repetitions;
                        end
                        Y(i) = real(det(M(3:4,1:2)));
                    else
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                            M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
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
                        Y(i) = real(det(MM{end}));
                    end
                    if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
                        if  o == Bisections && sign(Y(i-2)) ~= sign(Y(i))
                            XFine(end+1) = Frequency(i-1);
                        end
                        Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                        break
                    end
                end
            end
        end
        for p = 1:length(XFine)
            Wavenumber = 2*pi*XFine(p)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                if  SuperLayerSize > 1
                    L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                    L{m} = L1/L2;
                end
            end
            if  SuperLayerSize == 1
                U = L1(2:4,2:4)\-L1(2:4,1);
                u1(1) = 1+U(1)+U(2)*exp(1i*Wavenumber*Alpha(1)*LayerThicknesses)+U(3)*exp(1i*Wavenumber*Alpha(2)*LayerThicknesses);
                u1(2) = exp(1i*Wavenumber*Alpha(1)*LayerThicknesses)+U(1)*exp(1i*Wavenumber*Alpha(2)*LayerThicknesses)+U(2)+U(3);
            else
                M = L{1};
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                    M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
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
                    Wavenumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        Gamma = AlphaSH(m)*Wavenumber*LayerThicknesses(m);
                        L{m} = [cos(Gamma) 1i*sin(Gamma)/DSH(m);1i*DSH(m)*sin(Gamma) cos(Gamma)];
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
                    Wavenumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        Gamma = AlphaSH(m)*Wavenumber*LayerThicknesses(m);
                        L{m} = [cos(Gamma) 1i*sin(Gamma)/DSH(m);1i*DSH(m)*sin(Gamma) cos(Gamma)];
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
            Wavenumber = 2*pi*SweepRange(i)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                if  UseTMM
                    D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                    % R = [1 1 1 1 1 1;V(m,:) V(m,:);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:);D4(m,:) -D4(m,:)];
                    R(1,1)=1;R(1,2)=1;R(1,3)=1;R(1,4)=1;R(1,5)=1;R(1,6)=1;R(2,1)=V(m,1);R(2,2)=V(m,2);R(2,3)=V(m,3);R(2,4)=V(m,1);R(2,5)=V(m,2);R(2,6)=V(m,3);R(3,1)=W(m,1);R(3,2)=W(m,2);R(3,3)=W(m,3);R(3,4)=-W(m,1);R(3,5)=-W(m,2);R(3,6)=-W(m,3);R(4,1)=D3(m,1);R(4,2)=D3(m,2);R(4,3)=D3(m,3);R(4,4)=D3(m,1);R(4,5)=D3(m,2);R(4,6)=D3(m,3);R(5,1)=D5(m,1);R(5,2)=D5(m,2);R(5,3)=D5(m,3);R(5,4)=-D5(m,1);R(5,5)=-D5(m,2);R(5,6)=-D5(m,3);R(6,1)=D4(m,1);R(6,2)=D4(m,2);R(6,3)=D4(m,3);R(6,4)=-D4(m,1);R(6,5)=-D4(m,2);R(6,6)=-D4(m,3);
                    L{m} = R*D/R;
                else
                    E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                    % L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                    % L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                    L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,3);L1(1,4)=D3(m,1)*E(1);L1(1,5)=D3(m,2)*E(2);L1(1,6)=D3(m,3)*E(3);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=D5(m,3);L1(2,4)=-D5(m,1)*E(1);L1(2,5)=-D5(m,2)*E(2);L1(2,6)=-D5(m,3)*E(3);L1(3,1)=D4(m,1);L1(3,2)=D4(m,2);L1(3,3)=D4(m,3);L1(3,4)=-D4(m,1)*E(1);L1(3,5)=-D4(m,2)*E(2);L1(3,6)=-D4(m,3)*E(3);L1(4,1)=D3(m,1)*E(1);L1(4,2)=D3(m,2)*E(2);L1(4,3)=D3(m,3)*E(3);L1(4,4)=D3(m,1);L1(4,5)=D3(m,2);L1(4,6)=D3(m,3);L1(5,1)=D5(m,1)*E(1);L1(5,2)=D5(m,2)*E(2);L1(5,3)=D5(m,3)*E(3);L1(5,4)=-D5(m,1);L1(5,5)=-D5(m,2);L1(5,6)=-D5(m,3);L1(6,1)=D4(m,1)*E(1);L1(6,2)=D4(m,2)*E(2);L1(6,3)=D4(m,3)*E(3);L1(6,4)=-D4(m,1);L1(6,5)=-D4(m,2);L1(6,6)=-D4(m,3);
                    L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(m,1);L2(2,2)=V(m,2);L2(2,3)=V(m,3);L2(2,4)=V(m,1)*E(1);L2(2,5)=V(m,2)*E(2);L2(2,6)=V(m,3)*E(3);L2(3,1)=W(m,1);L2(3,2)=W(m,2);L2(3,3)=W(m,3);L2(3,4)=-W(m,1)*E(1);L2(3,5)=-W(m,2)*E(2);L2(3,6)=-W(m,3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(m,1)*E(1);L2(5,2)=V(m,2)*E(2);L2(5,3)=V(m,3)*E(3);L2(5,4)=V(m,1);L2(5,5)=V(m,2);L2(5,6)=V(m,3);L2(6,1)=W(m,1)*E(1);L2(6,2)=W(m,2)*E(2);L2(6,3)=W(m,3)*E(3);L2(6,4)=-W(m,1);L2(6,5)=-W(m,2);L2(6,6)=-W(m,3);
                    L{m} = L1/L2;
                end
            end
            M = L{1};
            if  UseTMM
                for m = 2:SuperLayerSize
                    M = M*L{m};
                end
                if  Repetitions > 1
                    M = M^Repetitions;
                end
                Y(i) = imag(det(M(4:6,1:3)));
            else
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                    M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
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
                Y(i) = real(det(MM{end}));
            end
            if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
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
                    Wavenumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        if  UseTMM
                            D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                            % R = [1 1 1 1 1 1;V(m,:) V(m,:);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:);D4(m,:) -D4(m,:)];
                            R(1,1)=1;R(1,2)=1;R(1,3)=1;R(1,4)=1;R(1,5)=1;R(1,6)=1;R(2,1)=V(m,1);R(2,2)=V(m,2);R(2,3)=V(m,3);R(2,4)=V(m,1);R(2,5)=V(m,2);R(2,6)=V(m,3);R(3,1)=W(m,1);R(3,2)=W(m,2);R(3,3)=W(m,3);R(3,4)=-W(m,1);R(3,5)=-W(m,2);R(3,6)=-W(m,3);R(4,1)=D3(m,1);R(4,2)=D3(m,2);R(4,3)=D3(m,3);R(4,4)=D3(m,1);R(4,5)=D3(m,2);R(4,6)=D3(m,3);R(5,1)=D5(m,1);R(5,2)=D5(m,2);R(5,3)=D5(m,3);R(5,4)=-D5(m,1);R(5,5)=-D5(m,2);R(5,6)=-D5(m,3);R(6,1)=D4(m,1);R(6,2)=D4(m,2);R(6,3)=D4(m,3);R(6,4)=-D4(m,1);R(6,5)=-D4(m,2);R(6,6)=-D4(m,3);
                            L{m} = R*D/R;
                        else
                            E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                            % L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                            % L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                            L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,3);L1(1,4)=D3(m,1)*E(1);L1(1,5)=D3(m,2)*E(2);L1(1,6)=D3(m,3)*E(3);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=D5(m,3);L1(2,4)=-D5(m,1)*E(1);L1(2,5)=-D5(m,2)*E(2);L1(2,6)=-D5(m,3)*E(3);L1(3,1)=D4(m,1);L1(3,2)=D4(m,2);L1(3,3)=D4(m,3);L1(3,4)=-D4(m,1)*E(1);L1(3,5)=-D4(m,2)*E(2);L1(3,6)=-D4(m,3)*E(3);L1(4,1)=D3(m,1)*E(1);L1(4,2)=D3(m,2)*E(2);L1(4,3)=D3(m,3)*E(3);L1(4,4)=D3(m,1);L1(4,5)=D3(m,2);L1(4,6)=D3(m,3);L1(5,1)=D5(m,1)*E(1);L1(5,2)=D5(m,2)*E(2);L1(5,3)=D5(m,3)*E(3);L1(5,4)=-D5(m,1);L1(5,5)=-D5(m,2);L1(5,6)=-D5(m,3);L1(6,1)=D4(m,1)*E(1);L1(6,2)=D4(m,2)*E(2);L1(6,3)=D4(m,3)*E(3);L1(6,4)=-D4(m,1);L1(6,5)=-D4(m,2);L1(6,6)=-D4(m,3);
                            L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(m,1);L2(2,2)=V(m,2);L2(2,3)=V(m,3);L2(2,4)=V(m,1)*E(1);L2(2,5)=V(m,2)*E(2);L2(2,6)=V(m,3)*E(3);L2(3,1)=W(m,1);L2(3,2)=W(m,2);L2(3,3)=W(m,3);L2(3,4)=-W(m,1)*E(1);L2(3,5)=-W(m,2)*E(2);L2(3,6)=-W(m,3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(m,1)*E(1);L2(5,2)=V(m,2)*E(2);L2(5,3)=V(m,3)*E(3);L2(5,4)=V(m,1);L2(5,5)=V(m,2);L2(5,6)=V(m,3);L2(6,1)=W(m,1)*E(1);L2(6,2)=W(m,2)*E(2);L2(6,3)=W(m,3)*E(3);L2(6,4)=-W(m,1);L2(6,5)=-W(m,2);L2(6,6)=-W(m,3);
                            L{m} = L1/L2;
                        end
                    end
                    M = L{1};
                    if  UseTMM
                        for m = 2:SuperLayerSize
                            M = M*L{m};
                        end
                        if  Repetitions > 1
                            M = M^Repetitions;
                        end
                        Y(i) = imag(det(M(4:6,1:3)));
                    else
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                            M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
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
                        Y(i) = real(det(MM{end}));
                    end
                    if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
                        if  o == Bisections && sign(Y(i-2)) ~= sign(Y(i))
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
            Wavenumber = 2*pi*SweepRange(i)*1e3/PhaseVelocityLimit;
            for m = 1:SuperLayerSize
                if  UseTMM
                    D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                    % R = [1 1 1 1;W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:)];
                    R(1,1)=1;R(1,2)=1;R(1,3)=1;R(1,4)=1;R(2,1)=W(m,1);R(2,2)=W(m,2);R(2,3)=-W(m,1);R(2,4)=-W(m,2);R(3,1)=D3(m,1);R(3,2)=D3(m,2);R(3,3)=D3(m,1);R(3,4)=D3(m,2);R(4,1)=D5(m,1);R(4,2)=D5(m,2);R(4,3)=-D5(m,1);R(4,4)=-D5(m,2);
                    L{m} = R*D/R;
                else
                    E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                    % L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                    % L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                    L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,1)*E(1);L1(1,4)=D3(m,2)*E(2);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=-D5(m,1)*E(1);L1(2,4)=-D5(m,2)*E(2);L1(3,1)=D3(m,1)*E(1);L1(3,2)=D3(m,2)*E(2);L1(3,3)=D3(m,1);L1(3,4)=D3(m,2);L1(4,1)=D5(m,1)*E(1);L1(4,2)=D5(m,2)*E(2);L1(4,3)=-D5(m,1);L1(4,4)=-D5(m,2);
                    L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(m,1);L2(2,2)=W(m,2);L2(2,3)=-W(m,1)*E(1);L2(2,4)=-W(m,2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(m,1)*E(1);L2(4,2)=W(m,2)*E(2);L2(4,3)=-W(m,1);L2(4,4)=-W(m,2);
                    L{m} = L1/L2;
                end
                Gamma = AlphaSH(m)*Wavenumber*LayerThicknesses(m);
                LSH{m} = [cos(Gamma) 1i*sin(Gamma)/DSH(m);1i*DSH(m)*sin(Gamma) cos(Gamma)];
            end
            M = L{1};
            MSH = LSH{1};
            if  UseTMM
                for m = 2:SuperLayerSize
                    M = M*L{m};
                    MSH = MSH*LSH{m};
                end
                if  Repetitions > 1
                    M = M^Repetitions;
                    MSH = MSH^Repetitions;
                end
                if  SymmetricSystem
                    M2 = L{end};
                    for m = SuperLayerSize-1:-1:1
                        M2 = M2*L{m};
                    end
                    M = M*M2^Repetitions;
                end
                Y(i) = real(det(M(3:4,1:2)));
            else
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                    M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                    MSH = MSH*LSH{m};
                end
                if  Repetitions > 1
                    MSH = MSH^Repetitions;
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
                Y(i) = real(det(MM{end}));
            end
            YSH(i) = MSH(2,1);
            if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))    
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
                    Wavenumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        if  UseTMM
                            D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                            % R = [1 1 1 1;W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:)];
                            R(1,1)=1;R(1,2)=1;R(1,3)=1;R(1,4)=1;R(2,1)=W(m,1);R(2,2)=W(m,2);R(2,3)=-W(m,1);R(2,4)=-W(m,2);R(3,1)=D3(m,1);R(3,2)=D3(m,2);R(3,3)=D3(m,1);R(3,4)=D3(m,2);R(4,1)=D5(m,1);R(4,2)=D5(m,2);R(4,3)=-D5(m,1);R(4,4)=-D5(m,2);
                            L{m} = R*D/R;
                        else
                            E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                            % L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                            % L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                            L1(1,1)=D3(m,1);L1(1,2)=D3(m,2);L1(1,3)=D3(m,1)*E(1);L1(1,4)=D3(m,2)*E(2);L1(2,1)=D5(m,1);L1(2,2)=D5(m,2);L1(2,3)=-D5(m,1)*E(1);L1(2,4)=-D5(m,2)*E(2);L1(3,1)=D3(m,1)*E(1);L1(3,2)=D3(m,2)*E(2);L1(3,3)=D3(m,1);L1(3,4)=D3(m,2);L1(4,1)=D5(m,1)*E(1);L1(4,2)=D5(m,2)*E(2);L1(4,3)=-D5(m,1);L1(4,4)=-D5(m,2);
                            L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(m,1);L2(2,2)=W(m,2);L2(2,3)=-W(m,1)*E(1);L2(2,4)=-W(m,2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(m,1)*E(1);L2(4,2)=W(m,2)*E(2);L2(4,3)=-W(m,1);L2(4,4)=-W(m,2);
                            L{m} = L1/L2;
                        end
                    end
                    M = L{1};
                    if  UseTMM
                        for m = 2:SuperLayerSize
                            M = M*L{m};
                        end
                        if  Repetitions > 1
                            M = M^Repetitions;
                        end
                        if  SymmetricSystem
                            M2 = L{end};
                            for m = SuperLayerSize-1:-1:1
                                M2 = M2*L{m};
                            end
                            M = M*M2^Repetitions;
                        end
                        Y(i) = real(det(M(3:4,1:2)));
                    else
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                            M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
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
                        Y(i) = real(det(MM{end}));
                    end
                    if  i > 2 && abs(Y(i-1)) < abs(Y(i-2)) && abs(Y(i-1)) < abs(Y(i))
                        if  o == Bisections && sign(Y(i-2)) ~= sign(Y(i))
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
                    Wavenumber = 2*pi*Frequency(i)*1e3/PhaseVelocityLimit;
                    for m = 1:SuperLayerSize
                        Gamma = AlphaSH(m)*Wavenumber*LayerThicknesses(m);
                        L{m} = [cos(Gamma) 1i*sin(Gamma)/DSH(m);1i*DSH(m)*sin(Gamma) cos(Gamma)];
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