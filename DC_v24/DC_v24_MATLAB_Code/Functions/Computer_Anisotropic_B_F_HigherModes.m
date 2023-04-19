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
function [B,X1,MissingModes] = Computer_Anisotropic_B_F_HigherModes(FluidLoading,FluidVelocity,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRange,PlateThickness,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,SymmetricSystem,I,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,B,X1,MissingModes,X0,T,p,g,g1,a11,a12,a21,a22,a23,a31,a32,a33,a34)   
%#ok<*AGROW> 
%#ok<*MINV>
%#ok<*GVMIS> 
global Stop 
Stop = 0;
for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
    if  Stop == 1
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    for m = 1:SuperLayerSize
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = Material{m}.Density^2*AngularFrequency^4;
        r3w6(m) = Material{m}.Density^3*AngularFrequency^6;
        b12(m) = a12(m)*rw2(m);
        b22(m) = a22(m)*rw2(m);
        b23(m) = a23(m)*r2w4(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = a33(m)*r2w4(m);
        b34(m) = a34(m)*r3w6(m);
    end
    X(i,1) = 0;
    Neighbors = [];
    for j = 1:length(B)
        if  i <= height(B{j})
            Neighbors(j,:) = [B{j}(i,8) B{j}(i,9)];
        end
    end
    NeighborsNumber = height(Neighbors);
    for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
        if  all(X(:,1) == 0) 
            SweepRangeReal = [1.1*PhaseVelocityLimit X0(2)];
        elseif numel(find(X(:,1) ~= 0)) == 1
            SweepRangeReal = [PhaseVelocityLimit X0(2)];
        elseif numel(find(X(:,1) ~= 0)) == 2
            SweepRangeReal = [1.1*X(end-1,3) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
        else
            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
        end
        if  SweepRangeReal(1) == SweepRangeReal(2)
            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
        end
        if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
            SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
            SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
        end
        if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < X0(2)
            SweepRangeReal(2) = X0(2);
        elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
            SweepRangeReal(2) = 0;
        end
        if  all(X(:,1) == 0)
            SweepRangeImag = [-1000*T 0];
        elseif numel(find(X(:,1) ~= 0)) == 1
            SweepRangeImag = [4*SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; 
        else
            SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
        end
        if  all(SweepRangeImag == [0 0])
            SweepRangeImag = [-20*T 0];
        elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
            SweepRangeImag(2) = 0;
        end
        if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
            SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
            SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
        end
        if  SweepRangeImag(2) > 0
            SweepRangeImag(2) = 0;
        end
        for o = 1:SearchAreaSections+1% increase search resolution
            if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                break
            end 
            for k = 1:1e2 % search minimum in characteristic equation and converge upon it 
                if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                        SweepRangeReal = SweepRangeReal(1):.25^o/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        SweepRangeReal = SweepRangeReal(1):.25/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    end
                else
                    if  length(SweepRangeReal) == 2
                        SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                    end
                    if  length(SweepRangeImag) == 2
                        SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    end
                end
                if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                    break
                end
                WaveNumber = AngularFrequency./(repmat(SweepRangeReal',1,length(SweepRangeImag))+repmat(SweepRangeImag,length(SweepRangeReal),1)*1i);
                WaveNumber2 = WaveNumber.^2;
                WaveNumber4 = WaveNumber.^4;
                WaveNumber6 = WaveNumber.^6;
                n1 = 1:height(WaveNumber);
                n2 = 1:width(WaveNumber);
                if  SuperLayerSize > 1
                    for m = 1:SuperLayerSize
                        A1 = a11(m)*WaveNumber2+b12(m);
                        A2 = a21(m)*WaveNumber4+b22(m)*WaveNumber2+b23(m);
                        A3 = a31(m)*WaveNumber6+b32(m)*WaveNumber4+b33(m)*WaveNumber2+b34(m);
                        k3a = A2/3-A1.^2/9;
                        k3b = A1.^3/27-A1.*A2/6+A3/2;
                        k3c = (sqrt(k3b.^2+k3a.^3)-k3b).^(1/3);
                        k3d = k3a./(2*k3c)-k3c/2;
                        k3e = k3a./k3c;
                        k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
                        k31 = sqrt(k3d-k3f-A1/3);
                        k32 = sqrt(k3d+k3f-A1/3);
                        k33 = -sqrt(k3c-k3e-A1/3);
                        k312 = k31.^2;
                        k322 = k32.^2;
                        k332 = k33.^2;
                        V1(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31));
                        V2(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32));
                        V3(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33));
                        W1(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k312).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31)));
                        W2(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k322).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32)));
                        W3(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k332).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33)));
                        D31(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
                        D32(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
                        D33(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
                        D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+WaveNumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
                        D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+WaveNumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
                        D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+WaveNumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
                        D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+WaveNumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
                        D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+WaveNumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
                        D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+WaveNumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
                        E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
                        E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
                        E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
                    end
                end
                if  numel(find(X(:,1) ~= 0)) <= 3
                    for l = 2:width(WaveNumber)-1 % remove solutions of previously found lower modes
                        for j = 2:height(WaveNumber)-1
                            for n = 1:NeighborsNumber
                                if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1)    
                                    WaveNumber(j,l) = NaN;
                                end
                            end
                        end
                    end
                else
                    for l = 2:width(WaveNumber)-1 % remove solutions of previously found lower modes
                        for j = 2:height(WaveNumber)-1
                            for n = 1:NeighborsNumber
                                if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1) && SweepRangeImag(l-1) < Neighbors(n,2) && SweepRangeImag(l+1) > Neighbors(n,2)
                                    WaveNumber(j,l) = NaN;
                                end
                            end
                        end
                    end
                end
                Y = 0;
                for l = n2
                    for j = n1
                        if  isnan(WaveNumber(j,l))
                            Y(j,l) = NaN;
                        else
                            if  SuperLayerSize == 1
                                A1 = a11*WaveNumber2(j,l)+b12;
                                A2 = a21*WaveNumber4(j,l)+b22*WaveNumber2(j,l)+b23;
                                A3 = a31*WaveNumber6(j,l)+b32*WaveNumber4(j,l)+b33*WaveNumber2(j,l)+b34;
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
                                m11 = c{1}(1,1)*WaveNumber2(j,l)-rw2+c{1}(5,5)*k32;
                                m12 = c{1}(1,6)*WaveNumber2(j,l)+c{1}(4,5)*k32;
                                m13 = (c{1}(1,3)+c{1}(5,5))*WaveNumber(j,l)*k3;
                                m22 = c{1}(6,6)*WaveNumber2(j,l)-rw2+c{1}(4,4)*k32;
                                m23 = (c{1}(3,6)+c{1}(4,5))*WaveNumber(j,l)*k3;
                                V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                                W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                                D3 = 1i*(c{1}(1,3)*WaveNumber(j,l)+c{1}(3,6)*WaveNumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                                D4 = 1i*(c{1}(4,5)*(k3+WaveNumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                                D5 = 1i*(c{1}(5,5)*(k3+WaveNumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                                E = exp(1i*k3*LayerThicknesses);
                                L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                                L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                                L{1} = L1/L2;
                            else
                                for m = 1:SuperLayerSize
                                    L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
                                    L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
                                    L{m} = L1/L2;
                                end
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
                            if  FluidLoading
                                G = inv(MM{end});
                                if  ToggleUpperFluid && ToggleLowerFluid
                                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                                    k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                                    WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                                    WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l); % in the lower fluid
                                    Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                                elseif ToggleUpperFluid && ~ToggleLowerFluid
                                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                                    WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                                    Y(j,l) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                                elseif ~ToggleUpperFluid && ToggleLowerFluid
                                    k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                                    WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                                    Y(j,l) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                                end
                            else
                                Y(j,l) = abs(det(MM{end}));
                            end
                        end
                    end
                end
                if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  p==3&&i>=149
% if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))
% else
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
% close(f)
% end
                Min = zeros(size(Y,1),size(Y,2));
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            Min(j,l) = 1;
                        end
                    end
                end
                [b1,b2] = find(Min);
                if  ~isempty(b1) % one or multiple minima are found
                    if  length(b1) == 1
                        MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                    else
                        delta = [];
                        for l = 1:length(b1) % calculate which minimum lies closest to last solution
                            cp(l) = AngularFrequency/real(WaveNumber(b1(l),b2(l)));
                            delta(l) = abs(cp(l)-X(end-1,1));
                        end
                        if  all(X(:,1) == 0)
                            [~,l] = max(delta);
                        else
                            [~,l] = min(delta);
                        end
                        MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                    end
                else
                    if  numel(find(X(:,1) ~= 0)) <= 3 && k <= 2
                        if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                            SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                            SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                            if  SweepRangeReal(1) > 1.1*PhaseVelocityLimit
                                SweepRangeReal(1) = 1.1*PhaseVelocityLimit;
                            end
                            if  SweepRangeReal(2) < X0(2)
                                SweepRangeReal(2) = X0(2);
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        MIN = 0;
                        break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                    else
                        Min = zeros(size(Y,1),size(Y,2)); % find border minima
                        for j = 2:size(Y,1)-1
                            if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                Min(j,1) = 1;
                            end
                            if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
                                Min(j,size(Y,2)) = 1;
                            end
                        end
                        for j = 2:size(Y,2)-1
                            if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                Min(1,j) = 1;
                            end
                            if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
                                Min(size(Y,1),j) = 1;
                            end
                        end
                        [b1,b2] = find(Min);
                        if  ~isempty(b1) % one or multiple BORDER minima are found
                            if  length(b1) == 1
                                MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            else
                                Value = [];
                                for l = 1:length(b1)
                                    Value(l) = Y(b1(l),b2(l));
                                end
                                [~,l] = min(Value);
                                MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            end
                        else
                            if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                if  SweepRangeReal(2) < 0
                                    SweepRangeReal(2) = 0;
                                end
                                if  SweepRangeImag(2) > 0
                                    SweepRangeImag(2) = 0;
                                end
                            end
                            MIN = 0;
                            break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                        end
                    end
                end 
                if  k == 100 || (Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                    break
                end
                if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                    if  MIN(1) == 1
                        if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                        SweepRangeReal = [SweepRangeReal(1)+4*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
                    elseif MIN(2) == 1
                        if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                    elseif MIN(1) == size(Y,1)
                        if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                        SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                    elseif MIN(2) == size(Y,2)
                        if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                    end
                else
                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                        if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                    end
                end
                if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(1) > 1.1*PhaseVelocityLimit
                    SweepRangeReal(1) = 1.1*PhaseVelocityLimit;
                end
                if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < X0(2)
                    SweepRangeReal(2) = X0(2);
                elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                    SweepRangeReal(2) = 0;
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
            end
            if  any(MIN)
                if  (numel(find(X(:,1) ~= 0)) <= 3 && i <= height(B{p+2}) && AngularFrequency/real(WaveNumber(MIN(1),MIN(2))) < B{p+2}(i,4)*1e3) ||...
                    (numel(find(X(:,1) ~= 0)) >= 1 && numel(find(X(:,1) ~= 0)) <= 3 && AngularFrequency/real(WaveNumber(MIN(1),MIN(2))) > X(end-1,1)) ||...
                    SweepRangeReal(MIN(1)) == 0
                    Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                    NeighborsNumber = height(Neighbors);
                else
                    if  numel(find(X(:,1) ~= 0)) <= 5
                        Outlier = 0;
                    else
                        z1 = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        z2 = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)
                            Outlier = 1;
                        else
                            Outlier = 0;
                        end
                    end
                    if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                        X(i,1) = AngularFrequency/real(WaveNumber(MIN(1),MIN(2))); % phase velocity (m/s)
                        X(i,2) = imag(WaveNumber(MIN(1),MIN(2))); % attenuation (Np/m)
                        X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                        X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                        Misses(i) = 0;
                        BelowCutoff(i) = 0;
                        break
                    end
                end
                if  all(X(:,1) == 0) 
                    SweepRangeReal = [1.1*PhaseVelocityLimit X0(2)];
                elseif numel(find(X(:,1) ~= 0)) == 1
                    SweepRangeReal = [PhaseVelocityLimit X0(2)];
                elseif numel(find(X(:,1) ~= 0)) == 2
                    SweepRangeReal = [1.1*X(end-1,3) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                else
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                end
                if  SweepRangeReal(1) == SweepRangeReal(2)
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                end
                if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                    SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                    SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                end
                if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < X0(2)
                    SweepRangeReal(2) = X0(2);
                elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                    SweepRangeReal(2) = 0;
                end
                if  all(X(:,1) == 0)
                    SweepRangeImag = [-1000*T 0];
                elseif numel(find(X(:,1) ~= 0)) == 1
                    SweepRangeImag = [4*SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; 
                else
                    SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                end
                if  all(SweepRangeImag == [0 0])
                    SweepRangeImag = [-20*T 0];
                elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                    SweepRangeImag(2) = 0;
                end
                if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                    SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                    SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
            end
            if  numel(find(X(:,1) ~= 0)) > 20 && o == SearchAreaSections
                break
            end
        end
        if  X(i,1) > 0 % stop q-loop if minimum has been found
            break
        end
    end
    if  X(i,1) == 0 && any(X(:,1)) % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
        if  numel(find(X(:,1) > 0)) == 1
            X(i-1,:) = 0;
            Misses(i) = 0;
            BelowCutoff(i) = 1;
        else
            Smooth1 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth2 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth3 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth4 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
            Fit1 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth1,'cubicspline');
            Fit2 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth2,'cubicspline');
            Fit3 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth3,'cubicspline');
            Fit4 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth4,'cubicspline');
            X(i,1) = Fit1(FrequencyRange(i));
            X(i,2) = Fit2(FrequencyRange(i));
            X(i,3) = Fit3(FrequencyRange(i));
            X(i,4) = Fit4(FrequencyRange(i));
            if  X(i,2) < 0 % negative attenuation is impossible
                X(i,[2 4]) = 0;
            end
            Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
            BelowCutoff(i) = 0;
        end
    elseif all(X(:,1) == 0) % if we are scanning below the cut-off frequency of the damped mode
        Misses(i) = 0;
        BelowCutoff(i) = 1;
    end
    if  FrequencyResolution*length(BelowCutoff(BelowCutoff == 1)) > BelowCutoffWidth/PlateThickness/1e3 ||... % if we exceed the allowed scanning width (in kHz*mm) below the cut-off frequency of a damped mode without finding it; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one
        length(Misses) == ceil(H(p)/FrequencyResolution)+6 && numel(find(Misses == 1)) >= 2
        MissingModes(p+3) = 1;
        break
    end           
    if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
        X(end-MissingSamples:end,:) = [];
        Misses(end-MissingSamples:end) = 0;
        break
    end
    if  FluidLoading && any(X(:,1)) && X(end,3) < FluidVelocity && abs(X(end,4)) < Resolution % avoid jumping to Scholte modes
        X(end,:) = [];
        break
    end
    if  BelowCutoff(i) == 0 && X(i,1) > 0
        addpoints(g(p+3),FrequencyRange(i),X(i)/1e3);
        drawnow limitrate
    end
% if  p == 3
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRangeF(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
% end
end
if  all(X(:,1) == 0)
    MissingModes(p+3) = 1;
end
if  MissingModes(p+3) == 0
    z = find(X(:,1) == max(X(:,1))); % find non-zero data below the cut-off frequency
    X(1:z-1,:) = 0; % remove them
    Misses(1:z-1) = 0; % remove also misses
    if  X(z,2) == 0 % if frequency cut-off is at zero damping, remove it also
        X(z,:) = 0;
    end
    X(Misses(1:height(X)) == 1,:) = NaN;
    B{p+3}(:,1) = FrequencyRange(1:height(X));
    B{p+3}(:,2) = FrequencyRange(1:height(X))/1e3;
    B{p+3}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
    B{p+3}(:,4) = fillmissing(X(:,1),'spline')/1e3;
    B{p+3}(:,7) = fillmissing(X(:,2),'spline');
    B{p+3}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
    B{p+3}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
    B{p+3}(B{p+3}(:,7) < 0,7) = 0; % negative attenuation is impossible
else
    B{p+3} = B{p+2};
    X1{p+3}(1,1) = 0;
    return
end
if  max(B{p+3}(:,4))*1e3 > PhaseVelocityLimit
    X1{p+3}(1,1) = 0;
else
    [Max,MaxInd] = max(B{p+3}(:,8));
    PhaseVelocityRange = Max+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
    for i = 1:length(PhaseVelocityRange)
        X1{p+3}(i,1) = 0;
        Neighbors = [];
        NeighborsNumber = height(Neighbors);
        for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
            if  i == 1
                SweepRangeFrq = [B{p+3}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 B{p+3}(MaxInd,1)+4/PlateThickness/1e3];
                SweepRangeImag = [SearchWidthImag(1)*B{p+3}(MaxInd,9) SearchWidthImag(2)*B{p+3}(MaxInd,9)];
            else
                SweepRangeFrq = [X1{p+3}(end-1,1)-FrequencyOffset/PlateThickness/1e3 X1{p+3}(end-1,1)+4/PlateThickness/1e3];
                SweepRangeImag = [SearchWidthImag(1)*X1{p+3}(end-1,3) SearchWidthImag(2)*X1{p+3}(end-1,3)];
            end
            if  all(SweepRangeImag == [0 0])
                SweepRangeImag = [-20*T 0];
            end
            if  SweepRangeImag(2) > 0
                SweepRangeImag(2) = 0;
            end
            for o = 1:SearchAreaSections % increase search resolution
                if  isempty(SweepRangeFrq) || isempty(SweepRangeImag)
                    break
                end
                for k = 1:1e2 % search minimum in characteristic equation and converge upon it 
                    if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                        if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                            SweepRangeFrq = SweepRangeFrq(1):.25^o/q*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                            SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                        elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeFrq(1)-SweepRangeFrq(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeFrq = SweepRangeFrq(1):.25/q/o*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                            SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                        elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeFrq = SweepRangeFrq(1):.25/q*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                            SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                        end
                    else
                        if  length(SweepRangeFrq) == 2
                            SweepRangeFrq = SweepRangeFrq(1):.25*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                        end
                        if  length(SweepRangeImag) == 2
                            SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                        end
                    end
                    if  isempty(SweepRangeFrq) || isempty(SweepRangeImag)
                        break
                    end
                    AngularFrequency = 2*pi*SweepRangeFrq*1e3;
                    Omega = repmat(AngularFrequency',1',length(SweepRangeImag));
                    if  SuperLayerSize > 1
                        Omega2 = Omega.^2;
                        Omega4 = Omega.^4;
                        Omega6 = Omega.^6;
                    end
                    WaveNumber = Omega./(repmat((PhaseVelocityRange(i)+SweepRangeImag*1i),length(SweepRangeFrq),1));
                    WaveNumber2 = WaveNumber.^2;
                    WaveNumber4 = WaveNumber.^4;
                    WaveNumber6 = WaveNumber.^6;
                    n1 = 1:height(WaveNumber);
                    n2 = 1:width(WaveNumber);
                    if  SuperLayerSize > 1
                        for m = 1:SuperLayerSize
                            rw2 = Material{m}.Density*Omega2;
                            r2w4 = Material{m}.Density^2*Omega4;
                            A1 = a11(m)*WaveNumber2+a12(m)*rw2;
                            A2 = a21(m)*WaveNumber4+a22(m)*rw2.*WaveNumber2+a23(m)*r2w4;
                            A3 = a31(m)*WaveNumber6+a32(m)*rw2.*WaveNumber4+a33(m)*r2w4.*WaveNumber2+a34(m)*Material{m}.Density^3*Omega6;
                            k3a = A2/3-A1.^2/9;
                            k3b = A1.^3/27-A1.*A2/6+A3/2;
                            k3c = (sqrt(k3b.^2+k3a.^3)-k3b).^(1/3);
                            k3d = k3a./(2*k3c)-k3c/2;
                            k3e = k3a./k3c;
                            k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
                            k31 = sqrt(k3d-k3f-A1/3);
                            k32 = sqrt(k3d+k3f-A1/3);
                            k33 = -sqrt(k3c-k3e-A1/3);
                            k312 = k31.^2;
                            k322 = k32.^2;
                            k332 = k33.^2;
                            V1(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31).*(c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k312)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31));
                            V2(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32).*(c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k322)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32));
                            V3(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33).*(c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k332)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33));
                            W1(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2+c{m}(5,5)*k312).*(c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k312)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31)-((c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31)));
                            W2(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2+c{m}(5,5)*k322).*(c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k322)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32)-((c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32)));
                            W3(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2+c{m}(5,5)*k332).*(c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k332)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33)-((c{m}(6,6)*WaveNumber2-rw2+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33)));
                            D31(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
                            D32(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
                            D33(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
                            D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+WaveNumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
                            D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+WaveNumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
                            D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+WaveNumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
                            D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+WaveNumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
                            D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+WaveNumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
                            D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+WaveNumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
                            E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
                            E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
                            E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
                        end
                    end
                    if  ~isempty(Neighbors)
                        for l = 2:width(WaveNumber)-1 % remove solutions of previously found lower modes
                            for j = 2:height(WaveNumber)-1
                                for n = 1:NeighborsNumber
                                    if  SweepRangeFrq(j-1) > Neighbors(n,1) && SweepRangeFrq(j+1) < Neighbors(n,1) && SweepRangeImag(l-1) < Neighbors(n,2) && SweepRangeImag(l+1) > Neighbors(n,2)   
                                        WaveNumber(j,l) = NaN;
                                    end
                                end
                            end
                        end
                    end
                    Y = 0;
                    for l = n2
                        for j = n1
                            if  isnan(WaveNumber(j,l))
                                Y(j,l) = NaN;
                            else
                                if  SuperLayerSize == 1
                                    rw2 = Material{1}.Density*AngularFrequency(j)^2;
                                    r2w4 = Material{1}.Density^2*AngularFrequency(j)^4;
                                    A1 = a11*WaveNumber2(j,l)+a12*rw2;
                                    A2 = a21*WaveNumber4(j,l)+a22*rw2*WaveNumber2(j,l)+a23*r2w4;
                                    A3 = a31*WaveNumber6(j,l)+a32*rw2*WaveNumber4(j,l)+a33*r2w4*WaveNumber2(j,l)+a34*Material{1}.Density^3*AngularFrequency(j)^6;
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
                                    m11 = c{1}(1,1)*WaveNumber2(j,l)-rw2+c{1}(5,5)*k32;
                                    m12 = c{1}(1,6)*WaveNumber2(j,l)+c{1}(4,5)*k32;
                                    m13 = (c{1}(1,3)+c{1}(5,5))*WaveNumber(j,l)*k3;
                                    m22 = c{1}(6,6)*WaveNumber2(j,l)-rw2+c{1}(4,4)*k32;
                                    m23 = (c{1}(3,6)+c{1}(4,5))*WaveNumber(j,l)*k3;
                                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                                    D3 = 1i*(c{1}(1,3)*WaveNumber(j,l)+c{1}(3,6)*WaveNumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                                    D4 = 1i*(c{1}(4,5)*(k3+WaveNumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                                    D5 = 1i*(c{1}(5,5)*(k3+WaveNumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                                    E = exp(1i*k3*LayerThicknesses);
                                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                                    L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                                    L{1} = L1/L2;
                                else
                                    for m = 1:SuperLayerSize
                                        L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
                                        L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
                                        L{m} = L1/L2;
                                    end
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
                                if  FluidLoading
                                    G = inv(MM{end});
                                    if  ToggleUpperFluid && ToggleLowerFluid
                                        k3UpperFluid = sqrt(AngularFrequency(j)^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                                        k3LowerFluid = sqrt(AngularFrequency(j)^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                                        WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                                        WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                                        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency(j)^2/WaveNumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                                        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency(j)^2/WaveNumber(j,l); % in the lower fluid
                                        Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                                    elseif ToggleUpperFluid && ~ToggleLowerFluid
                                        k3UpperFluid = sqrt(AngularFrequency(j)^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                                        WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                                        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency(j)^2/WaveNumber(j,l);
                                        Y(j,l) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                                    elseif ~ToggleUpperFluid && ToggleLowerFluid
                                        k3LowerFluid = sqrt(AngularFrequency(j)^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                                        WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                                        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency(j)^2/WaveNumber(j,l);
                                        Y(j,l) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                                    end
                                else
                                    Y(j,l) = abs(det(MM{end}));
                                end
                            end
                        end
                    end
                    if  abs(SweepRangeImag(end)) < 1e-3 % in case the fluid has low density, this allows finding solutions at zero attenuation for A0
                        Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                    end
% if  p==4
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeFrq,20*log10(Y))
% else
% f = figure;surf(SweepRangeImag,SweepRangeFrq,20*log10(Y))
% end
% close(f)
% end
                    Min = zeros(size(Y,1),size(Y,2));
                    for l = 2:size(Y,2)-1
                        for j = 2:size(Y,1)-1
                            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                Min(j,l) = 1;
                            end
                        end
                    end
                    [b1,b2] = find(Min);
                    if  ~isempty(b1) % one or multiple minima are found
                        if  length(b1) == 1
                            MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                        else
                            delta = [];
                            for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                frq(l) = AngularFrequency(b1(l))/(2*pi)/1e3;
                                if  i == 1
                                    delta(l) = abs(frq(l)-B{p+3}(MaxInd,1));
                                else
                                    delta(l) = abs(frq(l)-X1{p+3}(end-1,1));
                                end
                            end
                            [~,l] = min(delta);
                            MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                        end
                    else
                        Min = zeros(size(Y,1),size(Y,2)); % find border minima
                        for j = 2:size(Y,1)-1
                            if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                Min(j,1) = 1;
                            end
                            if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
                                Min(j,size(Y,2)) = 1;
                            end
                        end
                        for j = 2:size(Y,2)-1
                            if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                Min(1,j) = 1;
                            end
                            if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
                                Min(size(Y,1),j) = 1;
                            end
                        end
                        [b1,b2] = find(Min);
                        if  ~isempty(b1) % one or multiple BORDER minima are found
                            if  length(b1) == 1
                                MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                            else
                                Value = [];
                                for l = 1:length(b1)
                                    Value(l) = Y(b1(l),b2(l));
                                end
                                [~,l] = min(Value);
                                MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                            end
                        else
                            if  q > 1% && o == 1 % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                                SweepRangeFrq = [SweepRangeFrq(1)+(q-1)*abs(SweepRangeFrq(1)-SweepRangeFrq(end)) SweepRangeFrq(end)-(q-1)*abs(SweepRangeFrq(1)-SweepRangeFrq(end))];
                                SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                if  SweepRangeImag(2) > 0
                                    SweepRangeImag(2) = 0;
                                end
                            end
                            MIN = 0;
                            break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                        end
                    end 
                    if  k == 100 || (Resolution > abs(SweepRangeFrq(1)-SweepRangeFrq(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                        break
                    end
                    if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                        if  MIN(1) == 1
                            if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeFrq(1)-SweepRangeFrq(end))
                                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                            end
                            SweepRangeFrq = [SweepRangeFrq(1)+4*abs(SweepRangeFrq(1)-SweepRangeFrq(end)) SweepRangeFrq(end)];
                        elseif MIN(2) == 1
                            if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                            end
                            SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                        elseif MIN(1) == size(Y,1)
                            if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeFrq(1)-SweepRangeFrq(end))
                                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                            end
                            SweepRangeFrq = [SweepRangeFrq(1) SweepRangeFrq(end)-4*abs(SweepRangeFrq(1)-SweepRangeFrq(end))];
                        elseif MIN(2) == size(Y,2)
                            if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                            end
                            SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                        end
                    else
                        if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                            if  Resolution < abs(SweepRangeFrq(1)-SweepRangeFrq(2))
                                SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                            end
                        elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeFrq(1)-SweepRangeFrq(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            if  Resolution < abs(SweepRangeFrq(1)-SweepRangeFrq(2))
                                SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                            end
                            if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                            end
                        elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                            end
                        end
                    end
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                end
                if  any(MIN)
                    z = isoutlier(vertcat(X1{p+3}(1:end-1,1),AngularFrequency(MIN(1))/(2*pi)/1e3),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                    if  ~z(end) || all(X1{p+3}(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                        X1{p+3}(i,1) = AngularFrequency(MIN(1))/(2*pi)/1e3; % frequency (kHz)
                        X1{p+3}(i,2) = imag(WaveNumber(MIN(1),MIN(2))); % attenuation (Np/m)
                        X1{p+3}(i,3) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                        X1{p+3}(i,4) = AngularFrequency(MIN(1))/real(WaveNumber(MIN(1),MIN(2))); % phase velocity (m/s)
                        break
                    end
                    if  i == 1
                        SweepRangeFrq = [B{p+3}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 B{p+3}(MaxInd,1)+4/PlateThickness/1e3];
                        SweepRangeImag = [SearchWidthImag(1)*B{p+3}(MaxInd,9) SearchWidthImag(2)*B{p+3}(MaxInd,9)];
                    else
                        SweepRangeFrq = [X1{p+3}(end-1,1)-FrequencyOffset/PlateThickness/1e3 X1{p+3}(end-1,1)+4/PlateThickness/1e3];
                        SweepRangeImag = [SearchWidthImag(1)*X1{p+3}(end-1,3) SearchWidthImag(2)*X1{p+3}(end-1,3)];
                    end
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-20*T 0];
                    end
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                end
            end
            if  X1{p+3}(i,1) > 0 % stop q-loop if minimum has been found
                break
            end
        end
        if  X1{p+3}(i,1) == 0 || X1{p+3}(i,4) > PhaseVelocityLimit
            break
        end            
        if  X1{p+3}(i,2) < 0 % negative attenuation is impossible
            X1{p+3}(i,2:3) = 0;
        end
        addpoints(g1(p+3),X1{p+3}(i,1),X1{p+3}(i,4)/1e3);
        drawnow limitrate
% disp(['c = ',num2str(PhaseVelocityRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)]);
    end
    if  length(X1) < p+3
        X1{p+3}(1,1) = 0;
    end
    if  X1{p+3}(1,1) > 0 
        X1{p+3}(:,7) = X1{p+3}(:,2);
        X1{p+3}(:,2) = X1{p+3}(:,1)/1e3;
        X1{p+3}(:,3) = X1{p+3}(:,1)*PlateThickness;
        X1{p+3}(:,4) = X1{p+3}(:,4)/1e3;
        X1{p+3}(X1{p+3}(:,1) == 0,:) = []; 
        X1{p+3} = flipud(X1{p+3});
    end
end