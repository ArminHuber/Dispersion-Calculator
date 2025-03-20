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
function F = Computer_Isotropic_Pipe_Fn(Multithreading,Q1,Q2,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,n,H,LineColor,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
F{1} = [];
if  ~Multithreading
    for p = 1:length(H)
        g(p) = animatedline(ax,'color',LineColor);
        g1(p) = animatedline(ax,'color',LineColor);
    end
end
Ro2 = Ro^2;
Ri2 = Ri^2;
Y = zeros(5^PhaseVelocitySections+1,3);
MissingModes(length(H)) = 0;
for p = 1:length(H)
    X = [];
    Misses = 0;
    BelowCutoff = 0;
    if  p == 1
        x = FrequencyRange(ceil(H(p)/FrequencyResolution)):FrequencyResolution/5:FrequencyRange(ceil(H(p)/FrequencyResolution)+1);
        f = x(find(x > H(p),1));
        if  isempty(f)
            f = FrequencyRange(ceil(H(p)/FrequencyResolution)+1);
        end
        FrequencyRange_Edit = [FrequencyRange(1:ceil(H(p)/FrequencyResolution)) f:FrequencyResolution/5:f+20*FrequencyResolution/5 FrequencyRange(find(FrequencyRange > f+20*FrequencyResolution/5))];
    else
        FrequencyRange_Edit = FrequencyRange;
    end
    for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange_Edit)
        if  Stop
            return
        end
        AngularFrequency = 2*pi*FrequencyRange_Edit(i)*1e3;
        kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
        X(i,1) = 0;
        Neighbors = [];
        if  p > 1
            for j = 1:length(F)
                if  FrequencyRange_Edit(i) <= F{j}(end,1)
                    z = find(abs(FrequencyRange_Edit(i)-F{j}(:,1)) == min(abs(FrequencyRange_Edit(i)-F{j}(:,1))));
                    Neighbors(j) = F{j}(z,4)*1e3; %#ok<*FNDSB> 
                end
            end
        end
        if  isempty(Neighbors)
            Neighbors = 0;
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
                    if  p == 1 && X(i-1) > min(X(X > 0))
                        SweepRange = [X(i-1)+10*abs(X(i-2)-X(i-1)) X(i-1)-10*abs(X(i-2)-X(i-1))];
                        if  SweepRange(1) > Material.TransverseVelocity
                            SweepRange(1) = Material.TransverseVelocity;
                        end
                    else
                        SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
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
                if  PhaseVelocityResolution < abs(SweepRange(1)-SweepRange(2))
                    Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
                else
                    Bisections = 1;
                end
                if  q == 1 && o == 0 && (any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) || (SweepRange(1) > Material.TransverseVelocity && SweepRange(2) < Material.TransverseVelocity)) % Neighbors cannot be excluded because SweepRange has only 2 elements
                    continue
                end
                for j = 1:length(SweepRange)-1
                    for m = 1:length(Neighbors)
                        if  SweepRange(j) > Neighbors(m) && SweepRange(j+1) < Neighbors(m)
                            if  j == 1
                                SweepRange(2) = NaN;
                            else
                                SweepRange(j) = NaN;
                            end
                        end
                    end
                    if  SweepRange(j) > Material.TransverseVelocity && SweepRange(j+1) < Material.TransverseVelocity
                        if  j == 1
                            SweepRange(2) = NaN;
                        else
                            SweepRange(j) = NaN;
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
                                k1 = AngularFrequency/PhaseVelocity(l);
                                k2 = k1^2;
                                x = sqrt(kL2-k2);
                                y = sqrt(kT2-k2);
                                y2 = y^2;
                                if  SweepRange(j) > Material.LongitudinalVelocity
                                    Lambda1 = 1;
                                    Lambda2 = 1;
                                    xRo = x*Ro;
                                    yRo = y*Ro;
                                    xRi = x*Ri;
                                    yRi = y*Ri;
                                    Zn0xi = besselj(n,xRi);
                                    Zn0xo = besselj(n,xRo);
                                    Zn0yi = besselj(n,yRi);
                                    Zn0yo = besselj(n,yRo);
                                    Zn1xi = besselj(n+1,xRi);
                                    Zn1xo = besselj(n+1,xRo);
                                    Zn1yi = besselj(n+1,yRi);
                                    Zn1yo = besselj(n+1,yRo);
                                    Wn0xi = bessely(n,xRi);
                                    Wn0xo = bessely(n,xRo);
                                    Wn0yi = bessely(n,yRi);
                                    Wn0yo = bessely(n,yRo);
                                    Wn1xi = bessely(n+1,xRi);
                                    Wn1xo = bessely(n+1,xRo);
                                    Wn1yi = bessely(n+1,yRi);
                                    Wn1yo = bessely(n+1,yRo);
                                elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                    Lambda1 = -1;
                                    Lambda2 = 1;
                                    xRo = abs(x*Ro);
                                    yRo = y*Ro;
                                    xRi = abs(x*Ri);
                                    yRi = y*Ri;
                                    Zn0xi = besseli(n,xRi);
                                    Zn0xo = besseli(n,xRo);
                                    Zn0yi = besselj(n,yRi);
                                    Zn0yo = besselj(n,yRo);
                                    Zn1xi = besseli(n+1,xRi);
                                    Zn1xo = besseli(n+1,xRo);
                                    Zn1yi = besselj(n+1,yRi);
                                    Zn1yo = besselj(n+1,yRo);
                                    Wn0xi = besselk(n,xRi);
                                    Wn0xo = besselk(n,xRo);
                                    Wn0yi = bessely(n,yRi);
                                    Wn0yo = bessely(n,yRo);
                                    Wn1xi = besselk(n+1,xRi);
                                    Wn1xo = besselk(n+1,xRo);
                                    Wn1yi = bessely(n+1,yRi);
                                    Wn1yo = bessely(n+1,yRo);
                                elseif SweepRange(j) < Material.TransverseVelocity
                                    Lambda1 = -1;
                                    Lambda2 = -1;
                                    xRo = abs(x*Ro);
                                    yRo = abs(y*Ro);
                                    xRi = abs(x*Ri);
                                    yRi = abs(y*Ri);
                                    Zn0xi = besseli(n,xRi);
                                    Zn0xo = besseli(n,xRo);
                                    Zn0yi = besseli(n,yRi);
                                    Zn0yo = besseli(n,yRo);
                                    Zn1xi = besseli(n+1,xRi);
                                    Zn1xo = besseli(n+1,xRo);
                                    Zn1yi = besseli(n+1,yRi);
                                    Zn1yo = besseli(n+1,yRo);
                                    Wn0xi = besselk(n,xRi);
                                    Wn0xo = besselk(n,xRo);
                                    Wn0yi = besselk(n,yRi);
                                    Wn0yo = besselk(n,yRo);
                                    Wn1xi = besselk(n+1,xRi);
                                    Wn1xo = besselk(n+1,xRo);
                                    Wn1yi = besselk(n+1,yRi);
                                    Wn1yo = besselk(n+1,yRo);
                                end
                                M(1,1) = (2*n*(n-1)-(y2-k2)*Ri2)*Zn0xi+2*Lambda1*xRi*Zn1xi;
                                M(1,2) = 2*k1*yRi*Ri*Zn0yi-2*k1*Ri*(n+1)*Zn1yi;
                                M(1,3) = -2*n*(n-1)*Zn0yi+2*Lambda2*n*yRi*Zn1yi;
                                M(1,4) = (2*n*(n-1)-(y2-k2)*Ri2)*Wn0xi+2*xRi*Wn1xi;
                                M(1,5) = 2*Lambda2*k1*yRi*Ri*Wn0yi-2*(n+1)*k1*Ri*Wn1yi;
                                M(1,6) = -2*n*(n-1)*Wn0yi+2*n*yRi*Wn1yi;
                                M(2,1) = 2*n*(n-1)*Zn0xi-2*Lambda1*n*xRi*Zn1xi;
                                M(2,2) = -k1*yRi*Ri*Zn0yi+2*k1*Ri*(n+1)*Zn1yi;
                                M(2,3) = -(2*n*(n-1)-y2*Ri2)*Zn0yi-2*Lambda2*yRi*Zn1yi;
                                M(2,4) = 2*n*(n-1)*Wn0xi-2*n*xRi*Wn1xi;
                                M(2,5) = -Lambda2*k1*yRi*Ri*Wn0yi+2*k1*Ri*(n+1)*Wn1yi;
                                M(2,6) = -(2*n*(n-1)-y2*Ri2)*Wn0yi-2*yRi*Wn1yi;
                                M(3,1) = 2*n*k1*Ri*Zn0xi-2*Lambda1*k1*xRi*Ri*Zn1xi;
                                M(3,2) = n*yRi*Zn0yi-(y2-k2)*Ri2*Zn1yi;
                                M(3,3) = -n*k1*Ri*Zn0yi;
                                M(3,4) = 2*n*k1*Ri*Wn0xi-2*k1*xRi*Ri*Wn1xi;
                                M(3,5) = Lambda2*n*yRi*Wn0yi-(y2-k2)*Ri2*Wn1yi;
                                M(3,6) = -n*k1*Ri*Wn0yi;
                                M(4,1) = (2*n*(n-1)-(y2-k2)*Ro2)*Zn0xo+2*Lambda1*xRo*Zn1xo;
                                M(4,2) = 2*k1*yRo*Ro*Zn0yo-2*k1*Ro*(n+1)*Zn1yo;
                                M(4,3) = -2*n*(n-1)*Zn0yo+2*Lambda2*n*yRo*Zn1yo;
                                M(4,4) = (2*n*(n-1)-(y2-k2)*Ro2)*Wn0xo+2*xRo*Wn1xo;
                                M(4,5) = 2*Lambda2*k1*yRo*Ro*Wn0yo-2*(n+1)*k1*Ro*Wn1yo;
                                M(4,6) = -2*n*(n-1)*Wn0yo+2*n*yRo*Wn1yo;
                                M(5,1) = 2*n*(n-1)*Zn0xo-2*Lambda1*n*xRo*Zn1xo;
                                M(5,2) = -k1*yRo*Ro*Zn0yo+2*k1*Ro*(n+1)*Zn1yo;
                                M(5,3) = -(2*n*(n-1)-y2*Ro2)*Zn0yo-2*Lambda2*yRo*Zn1yo;
                                M(5,4) = 2*n*(n-1)*Wn0xo-2*n*xRo*Wn1xo;
                                M(5,5) = -Lambda2*k1*yRo*Ro*Wn0yo+2*k1*Ro*(n+1)*Wn1yo;
                                M(5,6) = -(2*n*(n-1)-y2*Ro2)*Wn0yo-2*yRo*Wn1yo;
                                M(6,1) = 2*n*k1*Ro*Zn0xo-2*Lambda1*k1*xRo*Ro*Zn1xo;
                                M(6,2) = n*yRo*Zn0yo-(y2-k2)*Ro2*Zn1yo;
                                M(6,3) = -n*k1*Ro*Zn0yo;
                                M(6,4) = 2*n*k1*Ro*Wn0xo-2*k1*xRo*Ro*Wn1xo;
                                M(6,5) = Lambda2*n*yRo*Wn0yo-(y2-k2)*Ro2*Wn1yo;
                                M(6,6) = -n*k1*Ro*Wn0yo;
                                Y(j,l) = real(det(M));
                            end
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  (j == 1 && abs(Y(j,2)) < 1e2) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
                            PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif (j == 1 && abs(Y(j,2)) < 1e2) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
                            PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                            Y(j,1) = Y(j,2);
                        else
                            PhaseVelocity(2) = 0;
                            break
                        end
                    end
                    if  PhaseVelocity(2) > 0
%                             if  numel(find(X > 0)) <= 5
%                                 Outlier = 0;
%                             else
%                                 z = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
%                                 if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
%                                     Outlier = 1;
%                                 else
%                                     Outlier = 0;
%                                 end
%                             end
%                             if  ~Outlier || all(X == 0)
                            X(i,1) = PhaseVelocity(2);
                            Misses(i) = 0;
                            BelowCutoff(i) = 0;
                            break
%                             end
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
        if  FrequencyResolution*length(BelowCutoff(BelowCutoff == 1)) > BelowCutoffWidth/(Ro-Ri)/1e3
            MissingModes(p) = 1;
            break
        end
        if  i > MissingSamples && all(Misses(end-MissingSamples:end))
            X(end-MissingSamples:end) = [];
            Misses(end-MissingSamples:end) = 0;
            break
        end
        if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
            send(Q1,[FrequencyRange_Edit(i),X(i)/1e3,p])
        elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
            addpoints(g(p),FrequencyRange_Edit(i),X(i)/1e3);
            drawnow limitrate
        end
    end
    if  all(X == 0)
        MissingModes(p) = 1;
    end
    if  ~MissingModes(p)
        X(Misses(1:length(X)) == 1) = NaN;
        F{p}(:,1) = FrequencyRange_Edit(1:length(X));
        F{p}(:,2) = FrequencyRange_Edit(1:length(X))/1e3;
        F{p}(:,3) = FrequencyRange_Edit(1:length(X))*(Ro-Ri);
        F{p}(:,4) = fillmissing(X,'spline')/1e3;
    else
        if  p == 1
            F{1}(1,1) = 0;
        else
            F{p} = F{p-1};
        end
        X1{p}(1) = 0;
        continue
    end
    [Max,MaxInd] = max(F{p}(:,4));
    PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
    for i = 1:length(PhaseVelocityRange)
        if  Stop
            return
        end
        X1{p}(i,1) = 0;
        if  i == 1
            SweepRange = [F{p}(MaxInd,1)-FrequencyOffset/(Ro-Ri)/1e3 F{p}(MaxInd,1)+10/(Ro-Ri)/1e3];
        elseif i == 2
            SweepRange = [X1{p}(i-1)-FrequencyOffset/(Ro-Ri)/1e3 X1{p}(i-1)+10/(Ro-Ri)/1e3];
        else
            SweepRange = X1{p}(i-1)-2*abs(X1{p}(i-2)-X1{p}(i-1));
            if  X1{p}(i-1) < X1{p}(i-2)
                SweepRange(2) = X1{p}(i-1)+2/(Ro-Ri)/1e3;
            else
                SweepRange(2) = X1{p}(i-1)+5*abs(X1{p}(i-2)-X1{p}(i-1));
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
                        AngularFrequency = 2*pi*Frequency(l)*1e3;
                        k1 = AngularFrequency/PhaseVelocityRange(i);
                        k2 = k1^2;
                        x = sqrt((AngularFrequency/Material.LongitudinalVelocity)^2-k2);
                        y = sqrt((AngularFrequency/Material.TransverseVelocity)^2-k2);
                        y2 = y^2;
                        if  PhaseVelocityRange(i) > Material.LongitudinalVelocity
                            Lambda1 = 1;
                            Lambda2 = 1;
                            xRo = x*Ro;
                            yRo = y*Ro;
                            xRi = x*Ri;
                            yRi = y*Ri;
                            Zn0xi = besselj(n,xRi);
                            Zn0xo = besselj(n,xRo);
                            Zn0yi = besselj(n,yRi);
                            Zn0yo = besselj(n,yRo);
                            Zn1xi = besselj(n+1,xRi);
                            Zn1xo = besselj(n+1,xRo);
                            Zn1yi = besselj(n+1,yRi);
                            Zn1yo = besselj(n+1,yRo);
                            Wn0xi = bessely(n,xRi);
                            Wn0xo = bessely(n,xRo);
                            Wn0yi = bessely(n,yRi);
                            Wn0yo = bessely(n,yRo);
                            Wn1xi = bessely(n+1,xRi);
                            Wn1xo = bessely(n+1,xRo);
                            Wn1yi = bessely(n+1,yRi);
                            Wn1yo = bessely(n+1,yRo);
                        elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity && PhaseVelocityRange(i) > Material.TransverseVelocity
                            Lambda1 = -1;
                            Lambda2 = 1;
                            xRo = abs(x*Ro);
                            yRo = y*Ro;
                            xRi = abs(x*Ri);
                            yRi = y*Ri;
                            Zn0xi = besseli(n,xRi);
                            Zn0xo = besseli(n,xRo);
                            Zn0yi = besselj(n,yRi);
                            Zn0yo = besselj(n,yRo);
                            Zn1xi = besseli(n+1,xRi);
                            Zn1xo = besseli(n+1,xRo);
                            Zn1yi = besselj(n+1,yRi);
                            Zn1yo = besselj(n+1,yRo);
                            Wn0xi = besselk(n,xRi);
                            Wn0xo = besselk(n,xRo);
                            Wn0yi = bessely(n,yRi);
                            Wn0yo = bessely(n,yRo);
                            Wn1xi = besselk(n+1,xRi);
                            Wn1xo = besselk(n+1,xRo);
                            Wn1yi = bessely(n+1,yRi);
                            Wn1yo = bessely(n+1,yRo);
                        elseif PhaseVelocityRange(i) < Material.TransverseVelocity
                            Lambda1 = -1;
                            Lambda2 = -1;
                            xRo = abs(x*Ro);
                            yRo = abs(y*Ro);
                            xRi = abs(x*Ri);
                            yRi = abs(y*Ri);
                            Zn0xi = besseli(n,xRi);
                            Zn0xo = besseli(n,xRo);
                            Zn0yi = besseli(n,yRi);
                            Zn0yo = besseli(n,yRo);
                            Zn1xi = besseli(n+1,xRi);
                            Zn1xo = besseli(n+1,xRo);
                            Zn1yi = besseli(n+1,yRi);
                            Zn1yo = besseli(n+1,yRo);
                            Wn0xi = besselk(n,xRi);
                            Wn0xo = besselk(n,xRo);
                            Wn0yi = besselk(n,yRi);
                            Wn0yo = besselk(n,yRo);
                            Wn1xi = besselk(n+1,xRi);
                            Wn1xo = besselk(n+1,xRo);
                            Wn1yi = besselk(n+1,yRi);
                            Wn1yo = besselk(n+1,yRo);
                        end
                        M(1,1) = (2*n*(n-1)-(y2-k2)*Ri2)*Zn0xi+2*Lambda1*xRi*Zn1xi;
                        M(1,2) = 2*k1*yRi*Ri*Zn0yi-2*k1*Ri*(n+1)*Zn1yi;
                        M(1,3) = -2*n*(n-1)*Zn0yi+2*Lambda2*n*yRi*Zn1yi;
                        M(1,4) = (2*n*(n-1)-(y2-k2)*Ri2)*Wn0xi+2*xRi*Wn1xi;
                        M(1,5) = 2*Lambda2*k1*yRi*Ri*Wn0yi-2*(n+1)*k1*Ri*Wn1yi;
                        M(1,6) = -2*n*(n-1)*Wn0yi+2*n*yRi*Wn1yi;
                        M(2,1) = 2*n*(n-1)*Zn0xi-2*Lambda1*n*xRi*Zn1xi;
                        M(2,2) = -k1*yRi*Ri*Zn0yi+2*k1*Ri*(n+1)*Zn1yi;
                        M(2,3) = -(2*n*(n-1)-y2*Ri2)*Zn0yi-2*Lambda2*yRi*Zn1yi;
                        M(2,4) = 2*n*(n-1)*Wn0xi-2*n*xRi*Wn1xi;
                        M(2,5) = -Lambda2*k1*yRi*Ri*Wn0yi+2*k1*Ri*(n+1)*Wn1yi;
                        M(2,6) = -(2*n*(n-1)-y2*Ri2)*Wn0yi-2*yRi*Wn1yi;
                        M(3,1) = 2*n*k1*Ri*Zn0xi-2*Lambda1*k1*xRi*Ri*Zn1xi;
                        M(3,2) = n*yRi*Zn0yi-(y2-k2)*Ri2*Zn1yi;
                        M(3,3) = -n*k1*Ri*Zn0yi;
                        M(3,4) = 2*n*k1*Ri*Wn0xi-2*k1*xRi*Ri*Wn1xi;
                        M(3,5) = Lambda2*n*yRi*Wn0yi-(y2-k2)*Ri2*Wn1yi;
                        M(3,6) = -n*k1*Ri*Wn0yi;
                        M(4,1) = (2*n*(n-1)-(y2-k2)*Ro2)*Zn0xo+2*Lambda1*xRo*Zn1xo;
                        M(4,2) = 2*k1*yRo*Ro*Zn0yo-2*k1*Ro*(n+1)*Zn1yo;
                        M(4,3) = -2*n*(n-1)*Zn0yo+2*Lambda2*n*yRo*Zn1yo;
                        M(4,4) = (2*n*(n-1)-(y2-k2)*Ro2)*Wn0xo+2*xRo*Wn1xo;
                        M(4,5) = 2*Lambda2*k1*yRo*Ro*Wn0yo-2*(n+1)*k1*Ro*Wn1yo;
                        M(4,6) = -2*n*(n-1)*Wn0yo+2*n*yRo*Wn1yo;
                        M(5,1) = 2*n*(n-1)*Zn0xo-2*Lambda1*n*xRo*Zn1xo;
                        M(5,2) = -k1*yRo*Ro*Zn0yo+2*k1*Ro*(n+1)*Zn1yo;
                        M(5,3) = -(2*n*(n-1)-y2*Ro2)*Zn0yo-2*Lambda2*yRo*Zn1yo;
                        M(5,4) = 2*n*(n-1)*Wn0xo-2*n*xRo*Wn1xo;
                        M(5,5) = -Lambda2*k1*yRo*Ro*Wn0yo+2*k1*Ro*(n+1)*Wn1yo;
                        M(5,6) = -(2*n*(n-1)-y2*Ro2)*Wn0yo-2*yRo*Wn1yo;
                        M(6,1) = 2*n*k1*Ro*Zn0xo-2*Lambda1*k1*xRo*Ro*Zn1xo;
                        M(6,2) = n*yRo*Zn0yo-(y2-k2)*Ro2*Zn1yo;
                        M(6,3) = -n*k1*Ro*Zn0yo;
                        M(6,4) = 2*n*k1*Ro*Wn0xo-2*k1*xRo*Ro*Wn1xo;
                        M(6,5) = Lambda2*n*yRo*Wn0yo-(y2-k2)*Ro2*Wn1yo;
                        M(6,6) = -n*k1*Ro*Wn0yo;
                        Y(j,l) = det(M);
                    end
                    if  k == 1
                        Y(j+1,1) = Y(j,3);
                    end
                    if  (j == 1 && abs(Y(j,2)) < 1e5) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
                        Frequency = [Frequency(1) Frequency(1)+(Frequency(2)-Frequency(1))/2 Frequency(2)];
                        Y(j,3) = Y(j,2);
                    elseif (j == 1 && abs(Y(j,2)) < 1e5) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
                        Frequency = [Frequency(2) Frequency(2)+(Frequency(3)-Frequency(2))/2 Frequency(3)];
                        Y(j,1) = Y(j,2);
                    else
                        Frequency(2) = 0;
                        break
                    end
                end
                if  Frequency(2) > 0
%                         if  i < 4
%                             Outlier = 0;
%                         else
%                             z = isoutlier(vertcat(X1{p}(1:i-1),Frequency(2)),'movmedian',5,'ThresholdFactor',9);
%                             if  z(end) && abs(X1{p}(i-1)-Frequency(2)) > 1
%                                 Outlier = 1;
%                             else
%                                 Outlier = 0;
%                             end
%                         end
%                         if  ~Outlier
                        X1{p}(i,1) = Frequency(2);
                        break
%                         end
                end
            end
            if  X1{p}(i) > 0
                break
            end
        end
        if  X1{p}(i) == 0
            break
        end
        if  Multithreading
            send(Q2,[X1{p}(i),PhaseVelocityRange(i)/1e3,p])
        else
            addpoints(g1(p),X1{p}(i),PhaseVelocityRange(i)/1e3);
            drawnow limitrate
        end
    end
    X1{p}(:,2) = X1{p}(:,1)/1e3;
    X1{p}(:,3) = X1{p}(:,1)*(Ro-Ri);
    X1{p}(:,4) = PhaseVelocityRange(1:height(X1{p}))/1e3;
    if  X1{p}(end,1) == 0
        X1{p}(end,1) = H(p);
        X1{p}(end,2) = H(p)/1e3;
        X1{p}(end,3) = H(p)*(Ro-Ri);
        X1{p}(end,4) = PhaseVelocityLimit/1e3;
    end
    X1{p}(X1{p}(:,1) == 0,:) = [];
    X1{p} = flipud(X1{p});
end
X1(MissingModes == 1) = [];
F(MissingModes == 1) = [];
for p = 1:length(F)
    F{p}(F{p}(:,4) == 0,:) = [];
    F{p} = vertcat(X1{p},F{p});
    F{p}(:,6) = 0;
% F{p}(:,5) = smooth(((F{p}(:,4)).^2)./(F{p}(:,4)-F{p}(:,1).*differentiate(fit(F{p}(:,1),F{p}(:,4),'cubicspline'),F{p}(:,1))));
end