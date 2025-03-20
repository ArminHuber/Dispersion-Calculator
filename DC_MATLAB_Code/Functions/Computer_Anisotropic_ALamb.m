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
function ALamb = Computer_Anisotropic_ALamb(Multithreading,Q1,Q2,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
ALamb{1} = [];
N = 1;
if  ~Multithreading
    for p = 1:length(H)+N
        g(p) = animatedline(ax,'color','b');
        g1(p) = animatedline(ax,'color','b');
    end
end
for m = 1:SuperLayerSize
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
X = XA0;
if  Multithreading
    send(Q1,[FrequencyRange(1),XA0/1e3,1])
else
    addpoints(g(1),FrequencyRange(1),XA0/1e3);
    drawnow limitrate
end
for i = 2:length(FrequencyRange)
    if  Stop
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    X(i,1) = 0;
    if  i > 2
        if  X(i-2) == X(i-1)
            x = abs(X-X(i-1));
            x(x == 0) = [];
            delta = min(x);
        else
            delta = abs(X(i-2)-X(i-1));
        end
    end
    if  i == 2
        SweepRange = .1:10:.5*XS0(1);
    else
        SweepRange = [X(i-1)-2*LambPhaseVelocitySweepRange2*delta X(i-1)+2*LambPhaseVelocitySweepRange2*delta];
    end
    if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
        if  delta > PhaseVelocityResolution
            SweepRange(1) = SweepRange(1)-PhaseVelocityResolution;
            SweepRange(2) = SweepRange(2)+PhaseVelocityResolution;
        else
            SweepRange(1) = SweepRange(1)-2*delta;
            SweepRange(2) = SweepRange(2)+2*delta;
        end
    end
    if  SweepRange(1) < 0
        SweepRange(1) = 0;
    end
    if  MatrixMethods == 1
        if  min(SweepRange) > MatrixMethodLimit(i)
            MatrixMethod = 1; % TMM
        else
            MatrixMethod = 2; % SMM
        end
    else
        MatrixMethod = 2; % SMM
    end
    for o = 0:PhaseVelocitySections
        if  o > 0
            SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
        end
        if  i > 2 && delta < PhaseVelocityResolution
            Bisections = ceil(log2(delta/PhaseVelocityResolution)/log2(.5));
        else
            Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
        end
        if  i < 6 && Bisections < 20
            Bisections = 20;
        elseif i >= 6
            if  o == 0 && Bisections < 5
                Bisections = 5;
            elseif o > 0 && Bisections < 1
                Bisections = 1;
            end
        end
        for j = 1:length(SweepRange)-1
            if  j == 1
                PhaseVelocityIndices = [1 2 3];
            else
                PhaseVelocityIndices = [2 3];
            end
            PhaseVelocity = [SweepRange(j) SweepRange(j)+(SweepRange(j+1)-SweepRange(j))/2 SweepRange(j+1)];
            for k = 1:Bisections
                if  k > 1
                    PhaseVelocityIndices = 2;
                end
                for l = PhaseVelocityIndices(1):PhaseVelocityIndices(end)
                    Wavenumber = AngularFrequency/PhaseVelocity(l);
                    PhaseVelocity2 = PhaseVelocity(l)^2;
                    for m = 1:SuperLayerSize
                        rc2 = Material{m}.Density*PhaseVelocity2;
                        if  ~Decoupled
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
                            Alpha(1) = sqrt(d6+d7);
                            Alpha(2) = sqrt(d6-d7);
                            Alpha(3) = sqrt(d4-d5-d1);
                            Alpha2 = Alpha.^2;
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
                        else
                            A2 = a21(m)+a22(m)*rc2;
                            A3 = a31(m)+a32(m)*rc2+rc2^2;
                            d1 = sqrt(A2^2-2*A1(m)*A3);
                            Alpha(1) = sqrt((-A2+d1)/A1(m));
                            Alpha(2) = sqrt((-A2-d1)/A1(m));
                            W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha.^2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
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
                        if  ~Decoupled
                            Y(j,l) = imag(det(M{end}(4:6,[3 5:6])));
                        else
                            Y(j,l) = imag(det(M{end}(3:4,[2 4])));
                        end
                    elseif MatrixMethod == 2
                        if  ~Decoupled
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
                            Y(j,l) = real(det(M{end}(1:4,[1:3 6])));
                        else
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
                            Y(j,l) = imag(det(M{end}(1:3,[1:2 4])));
                        end
                    end
                end
% if  i>=2
% f = figure;line(PhaseVelocity,(Y(j,:))),yline(0),xline(PhaseVelocity(2))
% close(f)
% end
                if  k == 1
                    Y(j+1,1) = Y(j,3);
                end
                if  MatrixMethod == 1
                    if  sign(Y(j,1)) ~= sign(Y(j,2))
                        PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                        Y(j,3) = Y(j,2);
                    elseif sign(Y(j,2)) ~= sign(Y(j,3))
                        PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                        Y(j,1) = Y(j,2);
                    else
                        PhaseVelocity(2) = 0;
                        break
                    end
                elseif MatrixMethod == 2
                    if  (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))
                        PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                        Y(j,3) = Y(j,2);
                    elseif (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))
                        PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                        Y(j,1) = Y(j,2);
                    else
                        PhaseVelocity(2) = 0;
                        break
                    end
                end
            end
            if  PhaseVelocity(2) > 0
                if  i < 4
                    Outlier = 0;
                else
                    z = isoutlier(vertcat(X(1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
                    if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
                        Outlier = 1;
                    else
                        Outlier = 0;
                    end
                end
                if  ~Outlier || all(X == 0)
                    X(i,1) = PhaseVelocity(2);
                    Misses(i) = 0;
                    break
                end
            end
        end
        if  X(i) > 0
            break
        end
    end
    if  i > 2 && X(i) == 0
        Smooth = filloutliers(X(1:i-1),'spline','movmedian',5,'ThresholdFactor',1);
        Fit = fit((1:i-1)',Smooth,'cubicspline');
        X(i,1) = Fit(i);
        Misses(i) = 1;
    end
    if  i > MissingSamples && all(Misses(end-MissingSamples:end))
        X(end-MissingSamples:end) = [];
        Misses(end-MissingSamples:end) = 0;
        break
    end
    if  Multithreading
        send(Q1,[FrequencyRange(i),X(i)/1e3,1])
    else
        addpoints(g(1),FrequencyRange(i),X(i)/1e3);
        drawnow limitrate
    end
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k),' M = ',num2str(MatrixMethod)];
% if  Misses(i) == 1
%     String = append(String,' Miss');
% end
% disp(String)
end
X(Misses(1:length(X)) == 1) = NaN;
ALamb{1}(:,1) = FrequencyRange(1:length(X));
ALamb{1}(:,2) = FrequencyRange(1:length(X))/1e3;
ALamb{1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
ALamb{1}(:,4) = fillmissing(X,'spline')/1e3;
ALamb{1}(:,7) = 0;
if  HigherOrderModes && any(H)
    MissingModes(length(H)+N) = 0;
    for p = 1:length(H)
        X = [];
        Misses = 0;
        BelowCutoff = 0;
        for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
            if  Stop
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            X(i,1) = 0;
            Neighbors = [];
            for j = 1:length(ALamb)
                if  i <= height(ALamb{j})
                    Neighbors(j) = ALamb{j}(i,4)*1e3;
                end
            end
            for q = 1:2
                if  q == 1
                    if  all(X == 0)
                        SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                    elseif isscalar(find(X > 0))
                        SweepRange = [X(i-1) max(Neighbors(Neighbors < X(i-1)))+1];
                    else
                        if  abs(X(i-2)-X(i-1)) > abs(X(i-3)-X(i-2))
                            Factor = LambPhaseVelocitySweepRange1;
                        else
                            Factor = LambPhaseVelocitySweepRange2;
                        end
                        if  ~Hybrid && (strcmp(Material{1}.Class,'Transversely isotropic') || strcmp(Material{1}.Class,'Isotropic'))
                            SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
                        elseif (~Hybrid && (strcmp(Material{1}.Class,'Orthotropic') || strcmp(Material{1}.Class,'Cubic'))) || Hybrid
                            if  4*abs(X(i-2)-X(i-1)) > 50
                                Top = 50;
                            else
                                Top = 4*abs(X(i-2)-X(i-1));
                            end
                            SweepRange = [X(i-1)+Top X(i-1)-Factor*abs(X(i-2)-X(i-1))];
                        end
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
                        SweepRange(1) = SweepRange(1)+PhaseVelocityResolution;
                        SweepRange(2) = SweepRange(2)-PhaseVelocityResolution;
                    end
                    if  SweepRange(end) < max(Neighbors(Neighbors < X(i-1)))+1
                        SweepRange(end) = max(Neighbors(Neighbors < X(i-1)))+1;
                    end
                else
                    if  all(X == 0)
                        SweepRange = [PhaseVelocityLimit min(Neighbors)+1];
                    elseif isscalar(find(X > 0))
                        SweepRange = [X(i-1) min(Neighbors)+1];
                    else
                        SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)+.1*abs(X(i-2)-X(i-1))+50];
                    end
                end
                for o = 0:PhaseVelocitySections+1
                    if  q == 2 && o == 0 % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                    if  q == 1 && o >= PhaseVelocitySections-1 && numel(find(X > 0)) > 1
                        SweepRange(end) = X(i-1)-Factor*abs(X(i-2)-X(i-1));
                        if  SweepRange(end) < min(Neighbors)+1
                            SweepRange(end) = min(Neighbors)+1;
                        end
                    end
                    if  q == 1 && numel(find(X > 0)) <= 1 && o > 0
                        if  o < PhaseVelocitySections+1
                            SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
                        else
                            if  all(X == 0)
                                SweepRange = PhaseVelocityLimit:-.2^(o+1)*(SweepRange(1)-SweepRange(end)):PhaseVelocityLimit-.2*(SweepRange(1)-SweepRange(end));
                            elseif isscalar(find(X > 0))
                                SweepRange = X(i-1):-.2^(o+1)*(SweepRange(1)-SweepRange(end)):X(i-1)-.2*(SweepRange(1)-SweepRange(end));
                            end
                        end
                    elseif ((q == 1 && numel(find(X > 0)) > 1) || q > 1) && o > 0
                        SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                    if  MatrixMethods == 1
                        if  numel(find(X > 0)) < 3 || min(SweepRange) > MatrixMethodLimit(i)
                            MatrixMethod = 1; % TMM
                        else
                            MatrixMethod = 2; % SMM
                        end
                    else
                        MatrixMethod = 2; % SMM
                    end
                    if  PhaseVelocityResolution < abs(SweepRange(1)-SweepRange(2))
                        Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
                    else
                        Bisections = 1;
                    end
                    if  q == 1 && o == 0 && any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                    for j = 1:length(SweepRange)-1
                        for t = 1:length(Neighbors)
                            if  SweepRange(j) > Neighbors(t) && SweepRange(j+1) < Neighbors(t)
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
                            PhaseVelocityIndices = [1 2 3];
                        else
                            PhaseVelocityIndices = [2 3];
                        end
                        PhaseVelocity = [SweepRange(j) SweepRange(j)+(SweepRange(j+1)-SweepRange(j))/2 SweepRange(j+1)];
                        for k = 1:Bisections
                            if  k > 1
                                PhaseVelocityIndices = 2;
                            end
                            for l = PhaseVelocityIndices(1):PhaseVelocityIndices(end)
                                if  isnan(PhaseVelocity(l))
                                    Y(j,l) = NaN;
                                else
                                    Wavenumber = AngularFrequency/PhaseVelocity(l);
                                    PhaseVelocity2 = PhaseVelocity(l)^2;
                                    for m = 1:SuperLayerSize
                                        rc2 = Material{m}.Density*PhaseVelocity2;
                                        if  ~Decoupled
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
                                            Alpha(1) = sqrt(d6+d7);
                                            Alpha(2) = sqrt(d6-d7);
                                            Alpha(3) = sqrt(d4-d5-d1);
                                            Alpha2 = Alpha.^2;
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
                                        else
                                            A2 = a21(m)+a22(m)*rc2;
                                            A3 = a31(m)+a32(m)*rc2+rc2^2;
                                            d1 = sqrt(A2^2-2*A1(m)*A3);
                                            Alpha(1) = sqrt((-A2+d1)/A1(m));
                                            Alpha(2) = sqrt((-A2-d1)/A1(m));
                                            W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha.^2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
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
                                        if  ~Decoupled
                                            Y(j,l) = imag(det(M{end}(4:6,[3 5:6])));
                                        else
                                            Y(j,l) = imag(det(M{end}(3:4,[2 4])));
                                        end
                                    elseif MatrixMethod == 2
                                        if  ~Decoupled
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
                                            Y(j,l) = real(det(M{end}(1:4,[1:3 6])));
                                        else
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
                                            Y(j,l) = imag(det(M{end}(1:3,[1:2 4])));
                                        end
                                    end
                                end
                            end
                            if  k == 1
                                Y(j+1,1) = Y(j,3);
                            end
                            if  MatrixMethod == 1
                                if  sign(Y(j,1)) ~= sign(Y(j,2))
                                    PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                                    Y(j,3) = Y(j,2);
                                elseif sign(Y(j,2)) ~= sign(Y(j,3))
                                    PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                                    Y(j,1) = Y(j,2);
                                else
                                    PhaseVelocity(2) = 0;
                                    break
                                end
                            elseif MatrixMethod == 2
                                if  (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))
                                    PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                                    Y(j,3) = Y(j,2);
                                elseif (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))
                                    PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                                    Y(j,1) = Y(j,2);
                                else
                                    PhaseVelocity(2) = 0;
                                    break
                                end
                            end
                        end
                        if  PhaseVelocity(2) > 0
                            if  numel(find(X > 0)) <= 5
                                Outlier = 0;
                            else
                                z = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
                                if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
                                    Outlier = 1;
                                else
                                    Outlier = 0;
                                end
                            end
                            if  ~Outlier || all(X == 0)
                                X(i,1) = PhaseVelocity(2);
                                Misses(i) = 0;
                                BelowCutoff(i) = 0;
                                break
                            end
                        end
                    end
                    if  X(i) > 0 || (q == 1 && abs(SweepRange(1)-SweepRange(end)) < 1000 && o == PhaseVelocitySections) || (q == 2 && o == PhaseVelocitySections-1)
                        break
                    end
                end
                if  X(i) > 0
                    break
                end
            end
            if  X(i) == 0 && any(X)
                if  isscalar(find(X > 0))
                    X(i-1,1) = 0;
                    Misses(i) = 0;
                    BelowCutoff(i) = 1;
                else
                    Smooth = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1),'spline','movmedian',5,'ThresholdFactor',1);
                    Fit = fit((ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth,'cubicspline');
                    X(i,1) = Fit(i);
                    Misses(i) = 1;
                    BelowCutoff(i) = 0;
                end
            elseif all(X == 0)
                Misses(i) = 0;
                BelowCutoff(i) = 1;
            end
            if  FrequencyResolution*length(BelowCutoff(BelowCutoff == 1)) > BelowCutoffWidth/PlateThickness/1e3
                MissingModes(p+N) = 1;
                break
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end))
                X(end-MissingSamples:end) = [];
                Misses(end-MissingSamples:end) = 0;
                break
            end
            if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                send(Q1,[FrequencyRange(i),X(i)/1e3,p+N])
            elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                addpoints(g(p+N),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k),' M = ',num2str(MatrixMethod)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
        end
        if  all(X == 0)
            MissingModes(p+N) = 1;
        end
        if  ~MissingModes(p+N)
            X(Misses(1:length(X)) == 1) = NaN;
            ALamb{p+N}(:,1) = FrequencyRange(1:length(X));
            ALamb{p+N}(:,2) = FrequencyRange(1:length(X))/1e3;
            ALamb{p+N}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
            ALamb{p+N}(:,4) = fillmissing(X,'spline')/1e3;
        else
            ALamb{p+N} = ALamb{p+N-1};
            X1{p+N}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(ALamb{p+N}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        PhaseVelocityRange2 = PhaseVelocityRange.^2;
        for i = 1:length(PhaseVelocityRange)
            if  Stop
                return
            end
            if  MatrixMethods == 1
                if  PhaseVelocityRange(i) > MatrixMethodLimit(MaxInd)
                    MatrixMethod = 1; % TMM
                else
                    MatrixMethod = 2; % SMM
                end
            else
                MatrixMethod = 2; % SMM
            end
            for m = 1:SuperLayerSize
                rc2 = Material{m}.Density*PhaseVelocityRange2(i);
                if  ~Decoupled
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
                    Alpha(m,1) = sqrt(d6+d7);
                    Alpha(m,2) = sqrt(d6-d7);
                    Alpha(m,3) = sqrt(d4-d5-d1);
                    Alpha2 = Alpha(m,:).^2;
                    m11 = c{m}(1,1)+c{m}(5,5)*Alpha2-rc2;
                    m22 = c{m}(6,6)+c{m}(4,4)*Alpha2-rc2;
                    m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
                    m13 = (c{m}(1,3)+c{m}(5,5))*Alpha(m,:);
                    m23 = (c{m}(3,6)+c{m}(4,5))*Alpha(m,:);
                    m1 = m13.*m22-m12.*m23;
                    V(m,:) = (m11.*m23-m13.*m12)./m1;
                    W(m,:) = (m11.*m22-m12.^2)./-m1;
                    e1 = Alpha(m,:)+W(m,:);
                    e2 = Alpha(m,:).*V(m,:);
                    D3(m,:) = c{m}(1,3)+c{m}(3,6)*V(m,:)+c{m}(3,3)*Alpha(m,:).*W(m,:);
                    D4(m,:) = c{m}(4,5)*e1+c{m}(4,4)*e2;
                    D5(m,:) = c{m}(5,5)*e1+c{m}(4,5)*e2;
                else
                    A2 = a21(m)+a22(m)*rc2;
                    A3 = a31(m)+a32(m)*rc2+rc2^2;
                    d1 = sqrt(A2^2-2*A1(m)*A3);
                    Alpha(m,1) = sqrt((-A2+d1)/A1(m));
                    Alpha(m,2) = sqrt((-A2-d1)/A1(m));
                    W(m,:) = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha(m,:).^2)./((c{m}(1,3)+c{m}(5,5))*Alpha(m,:));
                    D3(m,:) = c{m}(1,3)+c{m}(3,3)*Alpha(m,:).*W(m,:);
                    D5(m,:) = c{m}(5,5)*(Alpha(m,:)+W(m,:));
                end
            end
            X1{p+N}(i,1) = 0;
            if  i == 1
                if  p == 1
                    SweepRange = [ALamb{p+N}(MaxInd,1)-2*FrequencyOffset/PlateThickness/1e3 ALamb{p+N}(MaxInd,1)+2/PlateThickness/1e3];
                else
                    SweepRange = [ALamb{p+N}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 ALamb{p+N}(MaxInd,1)+2/PlateThickness/1e3];
                end
            elseif i == 2
                if  p == 1
                    SweepRange = [X1{p+N}(i-1)-2*FrequencyOffset/PlateThickness/1e3 X1{p+N}(i-1)+2/PlateThickness/1e3];
                else
                    SweepRange = [X1{p+N}(i-1)-FrequencyOffset/PlateThickness/1e3 X1{p+N}(i-1)+2/PlateThickness/1e3];
                end
            else
                SweepRange = X1{p+N}(i-1)-2*abs(X1{p+N}(i-2)-X1{p+N}(i-1));
                if  X1{p+N}(i-1) < X1{p+N}(i-2)
                    SweepRange(2) = X1{p+N}(i-1)+2/PlateThickness/1e3;
                else
                    SweepRange(2) = X1{p+N}(i-1)+5*abs(X1{p+N}(i-2)-X1{p+N}(i-1));
                end
            end 
            if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
                SweepRange(1) = SweepRange(1)-PhaseVelocityResolution;
                SweepRange(2) = SweepRange(2)+PhaseVelocityResolution;
            end
            for o = 0:FrequencySections
                if  o > 0
                    SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
                end
                if  PhaseVelocityResolution < abs(SweepRange(1)-SweepRange(2))
                    Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
                else
                    Bisections = 1;
                end
                for j = 1:length(SweepRange)-1
                    if  j == 1
                        FrequencyIndices = [1 2 3];
                    else
                        FrequencyIndices = [2 3];
                    end
                    Frequency = [SweepRange(j) SweepRange(j)+(SweepRange(j+1)-SweepRange(j))/2 SweepRange(j+1)];
                    for k = 1:Bisections
                        if  k > 1
                            FrequencyIndices = 2;
                        end
                        for l = FrequencyIndices(1):FrequencyIndices(end)
                            Wavenumber = 2*pi*Frequency(l)*1e3/PhaseVelocityRange(i);
                            for m = 1:SuperLayerSize
                                if  SuperLayerSize == 1
                                    Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                else
                                    Phi = 1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m);
                                end
                                E = exp(Phi);
                                if  ~Decoupled
                                    if  MatrixMethod == 1
                                        E_ = exp(-Phi);
                                        L1 = [E E_;V(m,:).*E V(m,:).*E_;W(m,:).*E -W(m,:).*E_;D3(m,:).*E D3(m,:).*E_;D5(m,:).*E -D5(m,:).*E_;D4(m,:).*E -D4(m,:).*E_];
                                        L2 = [ones(1,6);V(m,:) V(m,:);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:);D4(m,:) -D4(m,:)];
                                    elseif MatrixMethod == 2
                                        L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                                        L2 = [ones(1,3) E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E ones(1,3);V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                                    end
                                else
                                    if  MatrixMethod == 1
                                        E_ = exp(-Phi);
                                        L1 = [E E_;W(m,:).*E -W(m,:).*E_;D3(m,:).*E D3(m,:).*E_;D5(m,:).*E -D5(m,:).*E_];
                                        L2 = [ones(1,4);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:)];
                                    elseif MatrixMethod == 2
                                        L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                                        L2 = [ones(1,2) E;W(m,:) -W(m,:).*E;E ones(1,2);W(m,:).*E -W(m,:)];
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
                                if  ~Decoupled
                                    Y(j,l) = imag(det(M{end}(4:6,[3 5:6])));
                                else
                                    Y(j,l) = imag(det(M{end}(3:4,[2 4])));
                                end
                            elseif MatrixMethod == 2
                                if  ~Decoupled
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
                                    Y(j,l) = real(det(M{end}(1:4,[1:3 6])));
                                else
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
                                    Y(j,l) = imag(det(M{end}(1:3,[1:2 4])));
                                end
                            end
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  MatrixMethod == 1
                            if  sign(Y(j,1)) ~= sign(Y(j,2))
                                Frequency = [Frequency(1) Frequency(1)+(Frequency(2)-Frequency(1))/2 Frequency(2)];
                                Y(j,3) = Y(j,2);
                            elseif sign(Y(j,2)) ~= sign(Y(j,3))
                                Frequency = [Frequency(2) Frequency(2)+(Frequency(3)-Frequency(2))/2 Frequency(3)];
                                Y(j,1) = Y(j,2);
                            else
                                Frequency(2) = 0;
                                break
                            end
                        elseif MatrixMethod == 2
                            if  (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))
                                Frequency = [Frequency(1) Frequency(1)+(Frequency(2)-Frequency(1))/2 Frequency(2)];
                                Y(j,3) = Y(j,2);
                            elseif (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))
                                Frequency = [Frequency(2) Frequency(2)+(Frequency(3)-Frequency(2))/2 Frequency(3)];
                                Y(j,1) = Y(j,2);
                            else
                                Frequency(2) = 0;
                                break
                            end
                        end
                    end
                    if  Frequency(2) > 0
                        if  i < 4
                            Outlier = 0;
                        else
                            z = isoutlier(vertcat(X1{p+N}(1:i-1),Frequency(2)),'movmedian',5,'ThresholdFactor',9);
                            if  z(end) && abs(X1{p+N}(i-1)-Frequency(2)) > 1
                                Outlier = 1;
                            else
                                Outlier = 0;
                            end
                        end
                        if  ~Outlier
                            X1{p+N}(i,1) = Frequency(2);
                            break
                        end
                    end
                end
                if  X1{p+N}(i) > 0
                    break
                end
            end
            if  X1{p+N}(i) == 0
                break
            end
            if  Multithreading
                send(Q2,[X1{p+N}(i),PhaseVelocityRange(i)/1e3,p+N])
            else
                addpoints(g1(p+N),X1{p+N}(i),PhaseVelocityRange(i)/1e3);
                drawnow limitrate
            end
        end
        X1{p+N}(:,2) = X1{p+N}/1e3;
        X1{p+N}(:,3) = X1{p+N}(:,1)*PlateThickness;
        X1{p+N}(:,4) = PhaseVelocityRange(1:height(X1{p+N}))/1e3;
        if  X1{p+N}(end,1) == 0
            X1{p+N}(end,1) = H(p);
            X1{p+N}(end,2) = H(p)/1e3;
            X1{p+N}(end,3) = H(p)*PlateThickness;
            X1{p+N}(end,4) = PhaseVelocityLimit/1e3;
        end
        X1{p+N}(X1{p+N}(:,1) == 0,:) = [];
        X1{p+N} = flipud(X1{p+N});
        Alpha = 0;
        V = 0;
        W = 0;
        D3 = 0;
        D4 = 0;
        D5 = 0;
    end
    X1(MissingModes == 1) = [];
    ALamb(MissingModes == 1) = [];
    for p = N+1:length(ALamb)
        ALamb{p}(ALamb{p}(:,4) == 0,:) = [];
        ALamb{p} = vertcat(X1{p},ALamb{p});
        ALamb{p}(:,7) = 0;
    end
end
% figure,hold on
% for p = 1:length(ALamb)
%     plot(ALamb{p}(:,1),ALamb{p}(:,4),'k');
% end
% line(FrequencyRange,MatrixMethodLimit/1e3,'color','r')