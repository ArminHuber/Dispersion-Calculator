% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function A = Computer_Anisotropic_Plate_Lamb(Multithreading,Q,ax,ModeFamily,Decoupled,DC,Material,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0,a1,a11,a12,a21,a22,a23,a31,a32,a33,a34,H,LineColor)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};B={};
if  ~Multithreading
    g = [animatedline(ax,'color',LineColor) animatedline(ax,'color',LineColor)];
end
Y = zeros(5^(PhaseVelocitySections+1)+1,3);
Neighbors0 = [];
if  ModeFamily == 1
    H = [zeros(1,length(XS0)) H];
elseif ModeFamily == 2
    H = [0 H];
elseif ModeFamily == 3
    H = [0 zeros(1,length(XS0)) H];
end
for p = 1:length(H)
    X = 0;
    Misses = false;
    BelowCutoff = false;
    A0 = ModeFamily > 1 && p == 1;
    S0 = (ModeFamily == 1 && p <= length(XS0)) || (ModeFamily == 3 && p > 1 && p <= 1+length(XS0));
    i1 = ceil(H(p)/FrequencyResolution)+1;
    for i = i1:length(FrequencyRange)
        if  Stop
            return
        end
        if  (A0 || S0) && i == 1
            if  A0
                X = XA0;
            elseif S0
                if  ModeFamily == 1
                    X = XS0(p);
                elseif ModeFamily == 3
                    X = XS0(p-1);
                end
            end
            if  Multithreading
                send(Q(1),[FrequencyRange(1) X/1e3])
            else
                addpoints(g(1),FrequencyRange(1),X/1e3);
                drawnow limitrate
            end
            continue
        end
        AngularFrequency = pi*FrequencyRange(i)*2e3;
        X(i,1) = 0;
        Misses(i) = false;
        BelowCutoff(i) = false;
        Neighbors = Neighbors0;
        for j = 1:length(A)
            if  i <= size(A{j},1)
                Neighbors(end+1) = A{j}(i,4);
            end
        end
        for q = 1:2
            if  q == 1
                if  A0
                    if  i == 2
                        SweepRange = .1:10:XS0(1)/2;
                    else
                        SweepRange = X(i-1)+LambPhaseVelocitySweepRange2*abs(X(i-2)-X(i-1))*[-2 2];
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                        SweepRange = SweepRange+Accuracy*[-1 1];
                    end
                    if  SweepRange(1) < 0
                        SweepRange(1) = 0;
                    end
                    if  ModeFamily == 3 && SweepRange(2) > .99*XS0(1)
                        SweepRange(2) = .99*XS0(1);
                    end
                elseif S0
                    if  i < 6
                        SweepRange = X(i-1)+20:-2:X(i-1)-20;
                    else
                        dy = abs(X(i-2)-X(i-1));
                        if  dy > abs(X(i-3)-X(i-2))
                            Factor = 2*LambPhaseVelocitySweepRange1;
                        else
                            Factor = 2*LambPhaseVelocitySweepRange2;
                        end
                        Top = 4*dy;
                        if  Top > 100
                            Top = 100;
                        end
                        SweepRange = X(i-1)+[Top -Factor*dy];
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                        SweepRange = SweepRange+Accuracy*[1 -1];
                    end
                    if  SweepRange(2) < 0
                        SweepRange(2) = 0;
                    end
                else
                    if  ~any(X)
                        SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                    elseif isscalar(find(X))
                        SweepRange = [X(i-1) max(Neighbors(Neighbors < X(i-1)))+1];
                    else
                        dy = abs(X(i-2)-X(i-1));
                        if  dy > abs(X(i-3)-X(i-2))
                            Factor = LambPhaseVelocitySweepRange1;
                        else
                            Factor = LambPhaseVelocitySweepRange2;
                        end
                        Top = 4*dy;
                        if  Top > 50
                            Top = 50;
                        end
                        SweepRange = X(i-1)+[Top -Factor*dy];
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                        SweepRange = SweepRange+Accuracy*[1 -1];
                    end
                    if  SweepRange(2) < max(Neighbors(Neighbors < X(i-1)))+1
                        SweepRange(2) = max(Neighbors(Neighbors < X(i-1)))+1;
                    end
                end
            else
                if  A0 || S0
                    continue
                end
                if  ~any(X)
                    SweepRange = [PhaseVelocityLimit min(Neighbors)+1];
                elseif isscalar(find(X))
                    SweepRange = [X(i-1) min(Neighbors)+1];
                else
                    SweepRange = X(i-1)+.1*abs(X(i-2)-X(i-1))+[0 50];
                end
            end
            if  ~H(p)
                if  MatrixMethods == 1
                    if  min(SweepRange) > MatrixMethodLimit(i)
                        MatrixMethod = 1; % TMM
                    else
                        MatrixMethod = 2; % SMM
                    end
                else
                    MatrixMethod = 2; % SMM
                end
            end
            for o = 0:PhaseVelocitySections+1
                if  (A0 || S0) && i > 2
                    if  o == 0
                        continue
                    else
                        SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                else
                    if  q == 2 && o == 0 % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                    if  q == 1 && o >= PhaseVelocitySections-1 && numel(find(X)) > 1
                        SweepRange(end) = X(i-1)-Factor*abs(X(i-2)-X(i-1));
                        if  SweepRange(end) < min(Neighbors)+1
                            SweepRange(end) = min(Neighbors)+1;
                        end
                    end
                    if  q == 1 && numel(find(X)) <= 1 && o > 0
                        if  o < PhaseVelocitySections+1
                            SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                        else
                            if  ~any(X)
                                SweepRange = PhaseVelocityLimit:-.2^(o+1)*(SweepRange(1)-SweepRange(end)):PhaseVelocityLimit-.2*(SweepRange(1)-SweepRange(end));
                            elseif isscalar(find(X))
                                SweepRange = X(i-1):-.2^(o+1)*(SweepRange(1)-SweepRange(end)):X(i-1)-.2*(SweepRange(1)-SweepRange(end));
                            end
                        end
                    elseif ((q == 1 && numel(find(X)) > 1) || q > 1) && o > 0
                        SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                    if  q == 1 && o == 0 && any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                end
                if  H(p)
                    if  MatrixMethods == 1
                        if  numel(find(X > 0)) < 3 || min(SweepRange) > MatrixMethodLimit(i)
                            MatrixMethod = 1; % TMM
                        else
                            MatrixMethod = 2; % SMM
                        end
                    else
                        MatrixMethod = 2; % SMM
                    end
                end
                Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Accuracy));
                if  Bisections < 1
                    Bisections = 1;
                end
                for j = 1:length(SweepRange)-1
                    for t = 1:length(Neighbors)
                        if  (SweepRange(j) > Neighbors(t) && SweepRange(j+1) < Neighbors(t)) || (SweepRange(j) < Neighbors(t) && SweepRange(j+1) > Neighbors(t))
                            if  j == 1
                                SweepRange(2) = NaN;
                            else
                                SweepRange(j) = NaN;
                            end
                        end
                    end
                end
                for j = 1:length(SweepRange)-1
                    if  j == 1
                        Idx = 1:3;
                    else
                        Idx = 2:3;
                    end
                    PhaseVelocity = [SweepRange(j) (SweepRange(j)+SweepRange(j+1))/2 SweepRange(j+1)];
                    for k = 1:Bisections
                        if  k > 1
                            Idx = 2;
                        end
                        for l = Idx(1):Idx(end)
                            if  isnan(PhaseVelocity(l))
                                Y(j,l) = NaN;
                            else
                                Wavenumber = AngularFrequency/PhaseVelocity(l);
                                PhaseVelocity2 = PhaseVelocity(l)^2;
                                for m = 1:SuperLayerSize
                                    rc2 = Material{m}.Density*PhaseVelocity2;
                                    if  Decoupled
                                        A2 = a21(m)+a22(m)*rc2;
                                        A3 = a31(m)+a32(m)*rc2+rc2^2;
                                        d1 = sqrt(A2^2-2*a1(m)*A3);
                                        Alpha2 = [d1-A2 -d1-A2]/a1(m);
                                        Alpha = sqrt(Alpha2);
                                        W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
                                        D3 = c{m}(1,3)+c{m}(3,3)*Alpha.*W;
                                        D5 = c{m}(5,5)*(Alpha+W);
                                        if  SuperLayerSize == 1
                                            Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                        else
                                            Phi = 1i*Wavenumber*Alpha*LayerThicknesses(m);
                                        end
                                        E = exp(Phi);
                                        if  MatrixMethod == 1
                                            E_ = exp(-Phi);
                                            L1 = [E E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_];
                                            L2 = [ones(1,4);W -W;D3 D3;D5 -D5];
                                        elseif MatrixMethod == 2
                                            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                                            L2 = [ones(1,2) E;W -W.*E;E ones(1,2);W.*E -W];
                                        end
                                    else
                                        if  DC(m)
                                            A2 = a21(m)+a22(m)*rc2;
                                            A3 = a31(m)+a32(m)*rc2+rc2^2;
                                            d1 = sqrt(A2^2-2*a1(m)*A3);
                                            Alpha2 = [d1-A2 -d1-A2]/a1(m);
                                            Alpha = sqrt(Alpha2);
                                            AlphaSH = sqrt((rc2-c{m}(6,6))/c{m}(4,4));
                                            W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
                                            D3 = c{m}(1,3)+c{m}(3,3)*Alpha.*W;
                                            D4 = c{m}(4,4)*AlphaSH;
                                            D5 = c{m}(5,5)*(Alpha+W);
                                            if  SuperLayerSize == 1
                                                Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                                PhiSH = .5i*Wavenumber*AlphaSH*LayerThicknesses;
                                            else
                                                Phi = 1i*Wavenumber*Alpha*LayerThicknesses(m);
                                                PhiSH = 1i*Wavenumber*AlphaSH*LayerThicknesses(m);
                                            end
                                            E = exp(Phi);
                                            ESH = exp(PhiSH);
                                            if  MatrixMethod == 1
                                                E_ = exp(-Phi);
                                                ESH_ = exp(-PhiSH);
                                                L1 = [E 0 E_ 0;zeros(1,2) ESH zeros(1,2) ESH_;W.*E 0 -W.*E_ 0;D3.*E 0 D3.*E_ 0;D5.*E 0 -D5.*E_ 0;zeros(1,2) D4*ESH zeros(1,2) -D4*ESH_];
                                                L2 = [ones(1,2) 0 ones(1,2) 0;zeros(1,2) 1 zeros(1,2) 1;W 0 -W 0;D3 0 D3 0;D5 0 -D5 0;zeros(1,2) D4 zeros(1,2) -D4];
                                            elseif MatrixMethod == 2
                                                L1 = [D3 0 D3.*E 0;D5 0 -D5.*E 0;zeros(1,2) D4 zeros(1,2) -D4*ESH;D3.*E 0 D3 0;D5.*E 0 -D5 0;zeros(1,2) D4*ESH zeros(1,2) -D4];
                                                L2 = [ones(1,2) 0 E 0;zeros(1,2) 1 zeros(1,2) ESH;W 0 -W.*E 0;E 0 ones(1,2) 0;zeros(1,2) ESH zeros(1,2) 1;W.*E 0 -W 0];
                                            end
                                        else
                                            r2c4 = rc2^2;
                                            A1 = a11(m)+a12(m)*rc2;
                                            A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
                                            A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*rc2^3;
                                            d1 = A1/3;
                                            d2 = A2/3-d1^2;
                                            d3 = d1^3-d1*A2/2+A3/2;
                                            d4 = (sqrt(d2^3+d3^2)-d3)^(1/3);
                                            d5 = d2/d4;
                                            d6 = (d5-d4)/2-d1;
                                            d7 = (d5+d4)/2i*sqrt(3);
                                            Alpha2 = [d6+d7 d6-d7 d4-d5-d1];
                                            Alpha = sqrt(Alpha2);
                                            m11 = c{m}(1,1)+c{m}(5,5)*Alpha2-rc2;
                                            m22 = c{m}(6,6)+c{m}(4,4)*Alpha2-rc2;
                                            m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
                                            m13 = (c{m}(1,3)+c{m}(5,5))*Alpha;
                                            m23 = (c{m}(3,6)+c{m}(4,5))*Alpha;
                                            m1 = m13.*m22-m12.*m23;
                                            V = (m11.*m23-m13.*m12)./m1;
                                            W = (m11.*m22-m12.^2)./-m1;
                                            e1 = Alpha+W;
                                            e2 = Alpha.*V;
                                            D3 = c{m}(1,3)+c{m}(3,6)*V+c{m}(3,3)*Alpha.*W;
                                            D4 = c{m}(4,5)*e1+c{m}(4,4)*e2;
                                            D5 = c{m}(5,5)*e1+c{m}(4,5)*e2;
                                            if  SuperLayerSize == 1
                                                Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                            else
                                                Phi = 1i*Wavenumber*Alpha*LayerThicknesses(m);
                                            end
                                            E = exp(Phi);
                                            if  MatrixMethod == 1
                                                E_ = exp(-Phi);
                                                L1 = [E E_;V.*E V.*E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_;D4.*E -D4.*E_];
                                                L2 = [ones(1,6);V V;W -W;D3 D3;D5 -D5;D4 -D4];
                                            elseif MatrixMethod == 2
                                                L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                                                L2 = [ones(1,3) E;V V.*E;W -W.*E;E ones(1,3);V.*E V;W.*E -W];
                                            end
                                        end
                                    end
                                    L{m} = L1/L2;
                                end
                                M{1} = L{1};
                                if  MatrixMethod == 1
                                    for m = 2:SuperLayerSize
                                        M{1} = M{1}*L{m};
                                    end
                                    for m = 1:length(Pattern)
                                        M{m+1} = M{m}*M{Pattern(m)};
                                    end
                                    if  Decoupled
                                        if  ModeFamily == 1
                                            Y(j,l) = imag(det(M{end}(3:4,[1 3])));
                                        elseif ModeFamily == 2
                                            Y(j,l) = imag(det(M{end}(3:4,[2 4])));
                                        elseif ModeFamily == 3
                                            Y(j,l) = real(det(M{end}(3:4,1:2)));
                                        end
                                    else
                                        if  ModeFamily == 1
                                            Y(j,l) = real(det(M{end}(4:6,[1 2 4])));
                                        elseif ModeFamily == 2
                                            Y(j,l) = imag(det(M{end}(4:6,[3 5 6])));
                                        elseif ModeFamily == 3
                                            Y(j,l) = imag(det(M{end}(4:6,1:3)));
                                        end
                                    end
                                elseif MatrixMethod == 2
                                    if  Decoupled
                                        for m = 2:SuperLayerSize
                                            M0 = L{m}(1:2,1:2)-M{1}(3:4,3:4);
                                            M1 = M{1}(1:2,3:4)/M0;
                                            M2 = L{m}(3:4,1:2)/M0;
                                            M{1} = [M{1}(1:2,1:2)+M1*M{1}(3:4,1:2) -M1*L{m}(1:2,3:4);M2*M{1}(3:4,1:2) L{m}(3:4,3:4)-M2*L{m}(1:2,3:4)];
                                        end
                                        for m = 1:length(Pattern)
                                            M0 = M{Pattern(m)}(1:2,1:2)-M{m}(3:4,3:4);
                                            M1 = M{m}(1:2,3:4)/M0;
                                            M2 = M{Pattern(m)}(3:4,1:2)/M0;
                                            M{m+1} = [M{m}(1:2,1:2)+M1*M{m}(3:4,1:2) -M1*M{Pattern(m)}(1:2,3:4);M2*M{m}(3:4,1:2) M{Pattern(m)}(3:4,3:4)-M2*M{Pattern(m)}(1:2,3:4)];
                                        end
                                        if  ModeFamily == 1
                                            Y(j,l) = imag(det(M{end}([1 2 4],1:3)));
                                        elseif ModeFamily == 2
                                            Y(j,l) = imag(det(M{end}(1:3,[1 2 4])));
                                        elseif ModeFamily == 3
                                            Y(j,l) = real(det(M{end}));
                                        end
                                    else
                                        for m = 2:SuperLayerSize
                                            M0 = L{m}(1:3,1:3)-M{1}(4:6,4:6);
                                            M1 = M{1}(1:3,4:6)/M0;
                                            M2 = L{m}(4:6,1:3)/M0;
                                            M{1} = [M{1}(1:3,1:3)+M1*M{1}(4:6,1:3) -M1*L{m}(1:3,4:6);M2*M{1}(4:6,1:3) L{m}(4:6,4:6)-M2*L{m}(1:3,4:6)];
                                        end
                                        for m = 1:length(Pattern)
                                            M0 = M{Pattern(m)}(1:3,1:3)-M{m}(4:6,4:6);
                                            M1 = M{m}(1:3,4:6)/M0;
                                            M2 = M{Pattern(m)}(4:6,1:3)/M0;
                                            M{m+1} = [M{m}(1:3,1:3)+M1*M{m}(4:6,1:3) -M1*M{Pattern(m)}(1:3,4:6);M2*M{m}(4:6,1:3) M{Pattern(m)}(4:6,4:6)-M2*M{Pattern(m)}(1:3,4:6)];
                                        end
                                        if  ModeFamily == 1
                                            Y(j,l) = imag(det(M{end}([1:3 5 6],1:5)));
                                        elseif ModeFamily == 2
                                            Y(j,l) = real(det(M{end}(1:4,[1:3 6])));
                                        elseif ModeFamily == 3
                                            Y(j,l) = real(det(M{end}));
                                        end
                                    end
                                end
                            end
                        end
% if  p==2&&i>=101
% if  j == 1 || k > 1
%     PhaseVelocity_ = PhaseVelocity;Y_ = Y(j,:);
% else
%     PhaseVelocity_(2*j:2*j+1) = PhaseVelocity(2:3);Y_(2*j:2*j+1) = Y(j,2:3);
% end 
% f1 = figure;plot(PhaseVelocity_,20*log10(abs(Y_)));%f2 = figure;plot(PhaseVelocity_,Y_),yline(0);
% % 1
% close([f1])
% end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  sign(Y(j,1)) ~= sign(Y(j,2)) && (MatrixMethod == 1 || (MatrixMethod == 2 && ((j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))))
                            PhaseVelocity = [PhaseVelocity(1) (PhaseVelocity(1)+PhaseVelocity(2))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif sign(Y(j,2)) ~= sign(Y(j,3)) && (MatrixMethod == 1 || (MatrixMethod == 2 && ((j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))))
                            PhaseVelocity = [PhaseVelocity(2) (PhaseVelocity(2)+PhaseVelocity(3))/2 PhaseVelocity(3)];
                            Y(j,1) = Y(j,2);
                        else
                            PhaseVelocity(2) = 0;
                            break
                        end
                    end
                    if  PhaseVelocity(2) > 0
                        Outlier = false;
                        if  numel(find(X)) > 4
                            Idx = i1:i-1;
                            if  length(Idx) > 50
                                Idx(1:end-50) = [];
                            end
                            z = isoutlier([X(Idx);PhaseVelocity(2)],'movmedian',5,'ThresholdFactor',9);
                            if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
                                Neighbors(end+1) = PhaseVelocity(2);
                                Outlier = true;
                            end
                        end
                        if  ~Outlier
                            X(i) = PhaseVelocity(2);
                            break
                        end
                    end
                end
                if  X(i) || (q == 1 && abs(SweepRange(1)-SweepRange(end)) < 1000 && o == PhaseVelocitySections) || (q == 2 && o == PhaseVelocitySections-1)
                    break
                end
            end
            if  X(i)
                break
            end
        end
        if  ~X(i) && any(X)
            if  isscalar(find(X))
                if  ~(A0 || S0)
                    X(i-1) = 0;
                    BelowCutoff(i-1:i) = true;
                end
            else
                Idx = i1:i-1;
                Idx(BelowCutoff(Idx)) = [];
                if  length(Idx) > 50
                    Idx(1:end-50) = [];
                end
                Fit = fit(FrequencyRange(Idx)',X(Idx),'cubicspline');
                X(i) = Fit(FrequencyRange(i));
                Misses(i) = true; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
% figure,line(FrequencyRange(Idx),X(Idx),'LineWidth',4,'color','r'),line(FrequencyRange([Idx i]),Fit(FrequencyRange([Idx i])),'LineWidth',1.5,'color','g')
            end
        elseif ~any(X)
            BelowCutoff(i) = true;
        end
        if  FrequencyResolution*length(BelowCutoff(BelowCutoff)) > BelowCutoffWidth/PlateThickness/1e3
            break
        elseif i > MissingSamples && all(Misses(i-MissingSamples+1:i))
            X(i-MissingSamples+1:i) = [];
            Misses(i-MissingSamples+1:i) = [];
            break
        end
        if  ~BelowCutoff(i) && ~Misses(i)
            if  Multithreading
                send(Q(1),[FrequencyRange(i) X(i)/1e3])
            else
                addpoints(g(1),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
        end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k),' M = ',num2str(MatrixMethod)];
% if  Misses(i)
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i)
%     disp(['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k),' M = ',num2str(MatrixMethod),' Miss'])
% end
    end
    if  numel(find(X)) < 2
        continue
    end
    X(Misses) = NaN;
    A{end+1} = [FrequencyRange(1:length(X))'*[1 1e-3 PlateThickness] fillmissing(X,'spline')];
    z = A{end}(:,4) > 0;
    if  Multithreading
        send(Q(1),A{end}(z,[1 4]))
    else
        line(ax,A{end}(z),A{end}(z,4)/1e3,'color',LineColor)
        drawnow limitrate
        clearpoints(g(1))
    end
    if  ~any(H) || (~HigherOrderModes && H(p+1))
        break
    elseif A0 || S0
        continue
    end
    B{end+1} = 0;
    Idx0 = find(A{end}(:,4),1);
    PhaseVelocityRange = PhaseVelocityStep*ceil(A{end}(Idx0,4)/PhaseVelocityStep):PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
    PhaseVelocityRange2 = PhaseVelocityRange.^2;
    if  Multithreading
        send(Q(2),[A{end}(Idx0) A{end}(Idx0,4)/1e3])
    else
        addpoints(g(2),A{end}(Idx0),A{end}(Idx0,4)/1e3);
        drawnow limitrate
    end
    for i = 1:length(PhaseVelocityRange)
        if  Stop
            return
        end
        if  MatrixMethods == 1
            if  PhaseVelocityRange(i) > MatrixMethodLimit(Idx0)
                MatrixMethod = 1; % TMM
            else
                MatrixMethod = 2; % SMM
            end
        else
            MatrixMethod = 2; % SMM
        end
        for m = 1:SuperLayerSize
            rc2 = Material{m}.Density*PhaseVelocityRange2(i);
            if  Decoupled || (~Decoupled && DC(m))
                A2 = a21(m)+a22(m)*rc2;
                A3 = a31(m)+a32(m)*rc2+rc2^2;
                d1 = sqrt(A2^2-2*a1(m)*A3);
                Alpha2 = [d1-A2 -d1-A2]/a1(m);
                if  ~Decoupled
                    Alpha2 = [Alpha2 0];
                end
                Alpha(m,1:length(Alpha2)) = sqrt(Alpha2);
                AlphaSH(m,1) = sqrt((rc2-c{m}(6,6))/c{m}(4,4));
                W(m,1:length(Alpha2)) = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha2)./((c{m}(1,3)+c{m}(5,5))*Alpha(m,:));
                D3(m,1:length(Alpha2)) = c{m}(1,3)+c{m}(3,3)*Alpha(m,:).*W(m,:);
                D4(m,1) = c{m}(4,4)*AlphaSH(m);
                D5(m,1:length(Alpha2)) = c{m}(5,5)*(Alpha(m,:)+W(m,:));
            else
                r2c4 = rc2^2;
                A1 = a11(m)+a12(m)*rc2;
                A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
                A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*rc2^3;
                d1 = A1/3;
                d2 = A2/3-d1^2;
                d3 = d1^3-d1*A2/2+A3/2;
                d4 = (sqrt(d2^3+d3^2)-d3)^(1/3);
                d5 = d2/d4;
                d6 = (d5-d4)/2-d1;
                d7 = (d5+d4)/2i*sqrt(3);
                Alpha2 = [d6+d7 d6-d7 d4-d5-d1];
                Alpha(m,1:3) = sqrt(Alpha2);
                m11 = c{m}(1,1)+c{m}(5,5)*Alpha2-rc2;
                m22 = c{m}(6,6)+c{m}(4,4)*Alpha2-rc2;
                m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
                m13 = (c{m}(1,3)+c{m}(5,5))*Alpha(m,:);
                m23 = (c{m}(3,6)+c{m}(4,5))*Alpha(m,:);
                m1 = m13.*m22-m12.*m23;
                V(m,1:3) = (m11.*m23-m13.*m12)./m1;
                W(m,1:3) = (m11.*m22-m12.^2)./-m1;
                e1 = Alpha(m,:)+W(m,:);
                e2 = Alpha(m,:).*V(m,:);
                D3(m,1:3) = c{m}(1,3)+c{m}(3,6)*V(m,:)+c{m}(3,3)*Alpha(m,:).*W(m,:);
                D4(m,1:3) = c{m}(4,5)*e1+c{m}(4,4)*e2;
                D5(m,1:3) = c{m}(5,5)*e1+c{m}(4,5)*e2;
            end
        end
        B{end}(i,1) = 0;
        if  i == 1
            X1 = A{end}(Idx0);
        else
            X1 = B{end}(i-1);
        end
        if  p == 2
            SweepRange = X1+[-2*FrequencyOffset 10]/PlateThickness/1e3;
        else
            SweepRange = X1+[-FrequencyOffset 10]/PlateThickness/1e3;
        end
        if  isscalar(B) && SweepRange(1) < FrequencyRange(1)
            SweepRange(1) = FrequencyRange(1);
        elseif length(B) > 1 && size(A{end-1},1) >= Idx0 && A{end}(Idx0,4) > A{end-1}(Idx0,4)
            if  PhaseVelocityRange(i) <= max(A{end-1}(:,4))
                z = find(PhaseVelocityRange(i) < A{end-1}(:,4),1,'last');
                X2 = max(A{end-1}(z:z+1));
                if  SweepRange(1) < X2
                    SweepRange(1) = X2;
                end
            elseif PhaseVelocityRange(i) >= B{end-1}(1,4) && PhaseVelocityRange(i) <= B{end-1}(end,4)
                z = find(PhaseVelocityRange(i) == B{end-1}(:,4),1);
                X2 = B{end-1}(z)+diff(SweepRange)/100;
                if  any(X2) && SweepRange(1) < X2
                    SweepRange(1) = X2;
                elseif ~any(X2)
                    SweepRange(1) = max(B{end-1}(:,1));
                end
            end
        end
        if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
            SweepRange = SweepRange+Accuracy*[-1 1];
        end
        for o = 0:FrequencySections
            if  o > 0
                SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
            end
            Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Accuracy));
            if  Bisections < 1
                Bisections = 1;
            end
            for j = 1:length(SweepRange)-1
                if  j == 1
                    Idx = 1:3;
                else
                    Idx = 2:3;
                end
                Frequency = [SweepRange(j) (SweepRange(j)+SweepRange(j+1))/2 SweepRange(j+1)];
                for k = 1:Bisections
                    if  k > 1
                        Idx = 2;
                    end
                    for l = Idx(1):Idx(end)
                        Wavenumber = pi*Frequency(l)*2e3/PhaseVelocityRange(i);
                        for m = 1:SuperLayerSize
                            if  Decoupled
                                if  SuperLayerSize == 1
                                    Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                else
                                    Phi = 1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m);
                                end
                                E = exp(Phi);
                                if  MatrixMethod == 1
                                    E_ = exp(-Phi);
                                    L1 = [E E_;W(m,:).*E -W(m,:).*E_;D3(m,:).*E D3(m,:).*E_;D5(m,:).*E -D5(m,:).*E_];
                                    L2 = [ones(1,4);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:)];
                                elseif MatrixMethod == 2
                                    L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                                    L2 = [ones(1,2) E;W(m,:) -W(m,:).*E;E ones(1,2);W(m,:).*E -W(m,:)];
                                end
                            else
                                if  DC(m)
                                    if  SuperLayerSize == 1
                                        Phi = .5i*Wavenumber*Alpha(m,1:2)*LayerThicknesses(m);
                                        PhiSH = .5i*Wavenumber*AlphaSH(m)*LayerThicknesses(m);
                                    else
                                        Phi = 1i*Wavenumber*Alpha(m,1:2)*LayerThicknesses(m);
                                        PhiSH = 1i*Wavenumber*AlphaSH(m)*LayerThicknesses(m);
                                    end
                                    E = exp(Phi);
                                    ESH = exp(PhiSH);
                                    if  MatrixMethod == 1
                                        E_ = exp(-Phi);
                                        ESH_ = exp(-PhiSH);
                                        L1 = [E 0 E_ 0;zeros(1,2) ESH zeros(1,2) ESH_;W(m,1:2).*E 0 -W(m,1:2).*E_ 0;D3(m,1:2).*E 0 D3(m,1:2).*E_ 0;D5(m,1:2).*E 0 -D5(m,1:2).*E_ 0;zeros(1,2) D4(m)*ESH zeros(1,2) -D4(m)*ESH_];
                                        L2 = [ones(1,2) 0 ones(1,2) 0;zeros(1,2) 1 zeros(1,2) 1;W(m,1:2) 0 -W(m,1:2) 0;D3(m,1:2) 0 D3(m,1:2) 0;D5(m,1:2) 0 -D5(m,1:2) 0;zeros(1,2) D4(m) zeros(1,2) -D4(m)];
                                    elseif MatrixMethod == 2
                                        L1 = [D3(m,1:2) 0 D3(m,1:2).*E 0;D5(m,1:2) 0 -D5(m,1:2).*E 0;zeros(1,2) D4(m) zeros(1,2) -D4(m)*ESH;D3(m,1:2).*E 0 D3(m,1:2) 0;D5(m,1:2).*E 0 -D5(m,1:2) 0;zeros(1,2) D4(m)*ESH zeros(1,2) -D4(m)];
                                        L2 = [ones(1,2) 0 E 0;zeros(1,2) 1 zeros(1,2) ESH;W(m,1:2) 0 -W(m,1:2).*E 0;E 0 ones(1,2) 0;zeros(1,2) ESH zeros(1,2) 1;W(m,1:2).*E 0 -W(m,1:2) 0];
                                    end
                                else
                                    if  SuperLayerSize == 1
                                        Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                    else
                                        Phi = 1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m);
                                    end
                                    E = exp(Phi);
                                    if  MatrixMethod == 1
                                        E_ = exp(-Phi);
                                        L1 = [E E_;V(m,:).*E V(m,:).*E_;W(m,:).*E -W(m,:).*E_;D3(m,:).*E D3(m,:).*E_;D5(m,:).*E -D5(m,:).*E_;D4(m,:).*E -D4(m,:).*E_];
                                        L2 = [ones(1,6);V(m,:) V(m,:);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:);D4(m,:) -D4(m,:)];
                                    elseif MatrixMethod == 2
                                        L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                                        L2 = [ones(1,3) E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E ones(1,3);V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                                    end
                                end
                            end
                            L{m} = L1/L2;
                        end
                        M{1} = L{1};
                        if  MatrixMethod == 1
                            for m = 2:SuperLayerSize
                                M{1} = M{1}*L{m};
                            end
                            for m = 1:length(Pattern)
                                M{m+1} = M{m}*M{Pattern(m)};
                            end
                            if  Decoupled
                                if  ModeFamily == 1
                                    Y(j,l) = imag(det(M{end}(3:4,[1 3])));
                                elseif ModeFamily == 2
                                    Y(j,l) = imag(det(M{end}(3:4,[2 4])));
                                elseif ModeFamily == 3
                                    Y(j,l) = real(det(M{end}(3:4,1:2)));
                                end
                            else
                                if  ModeFamily == 1
                                    Y(j,l) = real(det(M{end}(4:6,[1 2 4])));
                                elseif ModeFamily == 2
                                    Y(j,l) = imag(det(M{end}(4:6,[3 5 6])));
                                elseif ModeFamily == 3
                                    Y(j,l) = imag(det(M{end}(4:6,1:3)));
                                end
                            end
                        elseif MatrixMethod == 2
                            if  Decoupled
                                for m = 2:SuperLayerSize
                                    M0 = L{m}(1:2,1:2)-M{1}(3:4,3:4);
                                    M1 = M{1}(1:2,3:4)/M0;
                                    M2 = L{m}(3:4,1:2)/M0;
                                    M{1} = [M{1}(1:2,1:2)+M1*M{1}(3:4,1:2) -M1*L{m}(1:2,3:4);M2*M{1}(3:4,1:2) L{m}(3:4,3:4)-M2*L{m}(1:2,3:4)];
                                end
                                for m = 1:length(Pattern)
                                    M0 = M{Pattern(m)}(1:2,1:2)-M{m}(3:4,3:4);
                                    M1 = M{m}(1:2,3:4)/M0;
                                    M2 = M{Pattern(m)}(3:4,1:2)/M0;
                                    M{m+1} = [M{m}(1:2,1:2)+M1*M{m}(3:4,1:2) -M1*M{Pattern(m)}(1:2,3:4);M2*M{m}(3:4,1:2) M{Pattern(m)}(3:4,3:4)-M2*M{Pattern(m)}(1:2,3:4)];
                                end
                                if  ModeFamily == 1
                                    Y(j,l) = imag(det(M{end}([1 2 4],1:3)));
                                elseif ModeFamily == 2
                                    Y(j,l) = imag(det(M{end}(1:3,[1 2 4])));
                                elseif ModeFamily == 3
                                    Y(j,l) = real(det(M{end}));
                                end
                            else
                                for m = 2:SuperLayerSize
                                    M0 = L{m}(1:3,1:3)-M{1}(4:6,4:6);
                                    M1 = M{1}(1:3,4:6)/M0;
                                    M2 = L{m}(4:6,1:3)/M0;
                                    M{1} = [M{1}(1:3,1:3)+M1*M{1}(4:6,1:3) -M1*L{m}(1:3,4:6);M2*M{1}(4:6,1:3) L{m}(4:6,4:6)-M2*L{m}(1:3,4:6)];
                                end
                                for m = 1:length(Pattern)
                                    M0 = M{Pattern(m)}(1:3,1:3)-M{m}(4:6,4:6);
                                    M1 = M{m}(1:3,4:6)/M0;
                                    M2 = M{Pattern(m)}(4:6,1:3)/M0;
                                    M{m+1} = [M{m}(1:3,1:3)+M1*M{m}(4:6,1:3) -M1*M{Pattern(m)}(1:3,4:6);M2*M{m}(4:6,1:3) M{Pattern(m)}(4:6,4:6)-M2*M{Pattern(m)}(1:3,4:6)];
                                end
                                if  ModeFamily == 1
                                    Y(j,l) = imag(det(M{end}([1:3 5 6],1:5)));
                                elseif ModeFamily == 2
                                    Y(j,l) = real(det(M{end}(1:4,[1:3 6])));
                                elseif ModeFamily == 3
                                    Y(j,l) = real(det(M{end}));
                                end
                            end
                        end
                    end
                    if  k == 1
                        Y(j+1,1) = Y(j,3);
                    end
                    if  sign(Y(j,1)) ~= sign(Y(j,2)) && (MatrixMethod == 1 || (MatrixMethod == 2 && ((j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))))
                        Frequency = [Frequency(1) (Frequency(1)+Frequency(2))/2 Frequency(2)];
                        Y(j,3) = Y(j,2);
                    elseif sign(Y(j,2)) ~= sign(Y(j,3)) && (MatrixMethod == 1 || (MatrixMethod == 2 && ((j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))))
                        Frequency = [Frequency(2) (Frequency(2)+Frequency(3))/2 Frequency(3)];
                        Y(j,1) = Y(j,2);
                    else
                        Frequency(2) = 0;
                        break
                    end
                end
                if  Frequency(2) > 0
                    B{end}(i) = Frequency(2);
                    break
                end
            end
            if  B{end}(i)
                break
            end
        end
        if  ~B{end}(i)
            break
        end
        if  Multithreading
            send(Q(2),[B{end}(i) PhaseVelocityRange(i)/1e3])
        else
            addpoints(g(2),B{end}(i),PhaseVelocityRange(i)/1e3);
            drawnow limitrate
        end
    end
    B{end}(:,2:4) = [B{end}(:,1)*[1e-3 PlateThickness] PhaseVelocityRange(1:size(B{end},1))'];
    if  ~B{end}(i)
        B{end}(i,:) = [H(p)*[1 1e-3 PlateThickness] PhaseVelocityLimit];
    end
    B{end}(B{end}(:,1) <= 0,:) = [];
    if  Multithreading
        send(Q(2),[A{end}(Idx0,[1 4]);B{end}(:,[1 4])])
    else
        line(ax,[A{end}(Idx0);B{end}(:,1)],[A{end}(Idx0,4);B{end}(:,4)]/1e3,'color',LineColor)
        drawnow limitrate
        clearpoints(g(2))
    end
end
for p = 1:length(A) % post processing
    A{p}(A{p}(:,4) <= 0,:) = [];
    if  ModeFamily == 1 && p > length(XS0)
        A{p} = [flipud(B{p-length(XS0)});A{p}];
    elseif ModeFamily == 2 && p > 1
        A{p} = [flipud(B{p-1});A{p}];
    elseif ModeFamily == 3 && p > length(XS0)+1
        A{p} = [flipud(B{p-length(XS0)-1});A{p}];
    end
    A{p}(:,4) = A{p}(:,4)/1e3;
    A{p}(:,8:9) = 1;
% A{p}(:,5) = smooth(((A{p}(:,4)).^2)./(A{p}(:,4)-A{p}(:,1).*differentiate(fit(A{p}(:,1),A{p}(:,4),'cubicspline'),A{p}(:,1))));
end
% figure,hold on
% for p = 1:length(A)
%     plot(A{p}(:,1),A{p}(:,4),'k');
% end
% line(FrequencyRange,MatrixMethodLimit/1e3,'color','r')