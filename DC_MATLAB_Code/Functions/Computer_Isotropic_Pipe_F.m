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
function F = Computer_Isotropic_Pipe_F(Multithreading,Q1,Q2,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
F{1} = [];
if  ~Multithreading
    if  HigherOrderModes && any(H)
        for p = 1:length(H)+1
            g(p) = animatedline(ax,'color','b');
            g1(p) = animatedline(ax,'color','b');
        end
    else
        g = animatedline(ax,'color','b');
    end
end
Ro2 = Ro^2;
Ri2 = Ri^2;
Y = zeros(5^PhaseVelocitySections+1,3);
for i = 1:length(FrequencyRange)
    if  Stop
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
    kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
    X(i,1) = 0;
    if  i == 1
        SweepRange = 1e-2:2e2;
    elseif i == 2
        SweepRange = X(i-1):2e3;
    else
        SweepRange = [X(i-1)-10*abs((X(i-2)-X(i-1))) X(i-1)+10*abs((X(i-2)-X(i-1)))];
        if  SweepRange(end) > 1.01*Material.RayleighVelocity
            SweepRange(end) = 1.01*Material.RayleighVelocity;
        end
    end
    for o = 0:PhaseVelocitySections
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
                    k1 = AngularFrequency/PhaseVelocity(l);
                    k2 = k1^2;
                    x = sqrt(kL2-k2);
                    y = sqrt(kT2-k2);
                    y2 = y^2;
                    Lambda1 = -1;
                    Lambda2 = -1;
                    xRo = abs(x*Ro);
                    yRo = abs(y*Ro);
                    xRi = abs(x*Ri);
                    yRi = abs(y*Ri);
                    Z1xi = besseli(1,xRi);
                    Z1xo = besseli(1,xRo);
                    Z1yi = besseli(1,yRi);
                    Z1yo = besseli(1,yRo);
                    Z2xi = besseli(2,xRi);
                    Z2xo = besseli(2,xRo);
                    Z2yi = besseli(2,yRi);
                    Z2yo = besseli(2,yRo);
                    W1xi = besselk(1,xRi);
                    W1xo = besselk(1,xRo);
                    W1yi = besselk(1,yRi);
                    W1yo = besselk(1,yRo);
                    W2xi = besselk(2,xRi);
                    W2xo = besselk(2,xRo);
                    W2yi = besselk(2,yRi);
                    W2yo = besselk(2,yRo);
                    M(1,1) = -(y2-k2)*Ri2*Z1xi+2*Lambda1*xRi*Z2xi;
                    M(1,2) = 2*k1*yRi*Ri*Z1yi-4*k1*Ri*Z2yi;
                    M(1,3) = 2*Lambda2*yRi*Z2yi;
                    M(1,4) = -(y2-k2)*Ri2*W1xi+2*xRi*W2xi;
                    M(1,5) = 2*Lambda2*k1*yRi*Ri*W1yi-4*k1*Ri*W2yi;
                    M(1,6) = 2*yRi*W2yi;
                    M(2,1) = -2*Lambda1*xRi*Z2xi;
                    M(2,2) = -k1*yRi*Ri*Z1yi+4*k1*Ri*Z2yi;
                    M(2,3) = y2*Ri2*Z1yi-2*Lambda2*yRi*Z2yi;
                    M(2,4) = -2*xRi*W2xi;
                    M(2,5) = -Lambda2*k1*yRi*Ri*W1yi+4*k1*Ri*W2yi;
                    M(2,6) = y2*Ri2*W1yi-2*yRi*W2yi;
                    M(3,1) = 2*k1*Ri*Z1xi-2*Lambda1*k1*xRi*Ri*Z2xi;
                    M(3,2) = yRi*Z1yi-(y2-k2)*Ri2*Z2yi;
                    M(3,3) = -k1*Ri*Z1yi;
                    M(3,4) = 2*k1*Ri*W1xi-2*k1*xRi*Ri*W2xi;
                    M(3,5) = Lambda2*yRi*W1yi-(y2-k2)*Ri2*W2yi;
                    M(3,6) = -k1*Ri*W1yi;
                    M(4,1) = -(y2-k2)*Ro2*Z1xo+2*Lambda1*xRo*Z2xo;
                    M(4,2) = 2*k1*yRo*Ro*Z1yo-4*k1*Ro*Z2yo;
                    M(4,3) = 2*Lambda2*yRo*Z2yo;
                    M(4,4) = -(y2-k2)*Ro2*W1xo+2*xRo*W2xo;
                    M(4,5) = 2*Lambda2*k1*yRo*Ro*W1yo-4*k1*Ro*W2yo;
                    M(4,6) = 2*yRo*W2yo;
                    M(5,1) = -2*Lambda1*xRo*Z2xo;
                    M(5,2) = -k1*yRo*Ro*Z1yo+4*k1*Ro*Z2yo;
                    M(5,3) = y2*Ro2*Z1yo-2*Lambda2*yRo*Z2yo;
                    M(5,4) = -2*xRo*W2xo;
                    M(5,5) = -Lambda2*k1*yRo*Ro*W1yo+4*k1*Ro*W2yo;
                    M(5,6) = y2*Ro2*W1yo-2*yRo*W2yo;
                    M(6,1) = 2*k1*Ro*Z1xo-2*Lambda1*k1*xRo*Ro*Z2xo;
                    M(6,2) = yRo*Z1yo-(y2-k2)*Ro2*Z2yo;
                    M(6,3) = -k1*Ro*Z1yo;
                    M(6,4) = 2*k1*Ro*W1xo-2*k1*xRo*Ro*W2xo;
                    M(6,5) = Lambda2*yRo*W1yo-(y2-k2)*Ro2*W2yo;
                    M(6,6) = -k1*Ro*W1yo;
                    Y(j,l) = real(det(M));
                end
                if  k == 1
                    Y(j+1,1) = Y(j,3);
                end
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
            end
            if  PhaseVelocity(2) > 0
%                 if  i < 4
%                     Outlier = 0;
%                 else
%                     z = isoutlier(vertcat(X(1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
%                     if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
%                         Outlier = 1;
%                     else
%                         Outlier = 0;
%                     end
%                 end
%                 if  ~Outlier || all(X == 0)
                    X(i,1) = PhaseVelocity(2);
                    Misses(i) = 0;
                    break
%                 end
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
end
X(Misses(1:length(X)) == 1) = NaN;
F{1}(:,1) = FrequencyRange(1:length(X));
F{1}(:,2) = FrequencyRange(1:length(X))/1e3;
F{1}(:,3) = FrequencyRange(1:length(X))*(Ro-Ri);
F{1}(:,4) = fillmissing(X,'spline')/1e3;
F{1}(:,6) = 0;
% F{1}(:,5) = smooth(((F{1}(:,4)).^2)./(F{1}(:,4)-F{1}(:,1).*differentiate(fit(F{1}(:,1),F{1}(:,4),'cubicspline'),F{1}(:,1))));
if  HigherOrderModes && any(H)
    MissingModes(length(H)+1) = 0;
    for p = 1:length(H)
        X = [];
        Misses = 0;
        BelowCutoff = 0;
        for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
            if  Stop
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
            kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
            X(i,1) = 0;
            Neighbors = [];
            for j = 1:length(F)
                if  i <= height(F{j})
                    Neighbors(j) = F{j}(i,4)*1e3;
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
                        SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
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
                        for n = 1:length(Neighbors)
                            if  SweepRange(j) > Neighbors(n) && SweepRange(j+1) < Neighbors(n)
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
                                        Z1xi = besselj(1,xRi);
                                        Z1xo = besselj(1,xRo);
                                        Z1yi = besselj(1,yRi);
                                        Z1yo = besselj(1,yRo);
                                        Z2xi = besselj(2,xRi);
                                        Z2xo = besselj(2,xRo);
                                        Z2yi = besselj(2,yRi);
                                        Z2yo = besselj(2,yRo);
                                        W1xi = bessely(1,xRi);
                                        W1xo = bessely(1,xRo);
                                        W1yi = bessely(1,yRi);
                                        W1yo = bessely(1,yRo);
                                        W2xi = bessely(2,xRi);
                                        W2xo = bessely(2,xRo);
                                        W2yi = bessely(2,yRi);
                                        W2yo = bessely(2,yRo);
                                    elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                        Lambda1 = -1;
                                        Lambda2 = 1;
                                        xRo = abs(x*Ro);
                                        yRo = y*Ro;
                                        xRi = abs(x*Ri);
                                        yRi = y*Ri;
                                        Z1xi = besseli(1,xRi);
                                        Z1xo = besseli(1,xRo);
                                        Z1yi = besselj(1,yRi);
                                        Z1yo = besselj(1,yRo);
                                        Z2xi = besseli(2,xRi);
                                        Z2xo = besseli(2,xRo);
                                        Z2yi = besselj(2,yRi);
                                        Z2yo = besselj(2,yRo);
                                        W1xi = besselk(1,xRi);
                                        W1xo = besselk(1,xRo);
                                        W1yi = bessely(1,yRi);
                                        W1yo = bessely(1,yRo);
                                        W2xi = besselk(2,xRi);
                                        W2xo = besselk(2,xRo);
                                        W2yi = bessely(2,yRi);
                                        W2yo = bessely(2,yRo);
                                    elseif SweepRange(j) < Material.TransverseVelocity
                                        Lambda1 = -1;
                                        Lambda2 = -1;
                                        xRo = abs(x*Ro);
                                        yRo = abs(y*Ro);
                                        xRi = abs(x*Ri);
                                        yRi = abs(y*Ri);
                                        Z1xi = besseli(1,xRi);
                                        Z1xo = besseli(1,xRo);
                                        Z1yi = besseli(1,yRi);
                                        Z1yo = besseli(1,yRo);
                                        Z2xi = besseli(2,xRi);
                                        Z2xo = besseli(2,xRo);
                                        Z2yi = besseli(2,yRi);
                                        Z2yo = besseli(2,yRo);
                                        W1xi = besselk(1,xRi);
                                        W1xo = besselk(1,xRo);
                                        W1yi = besselk(1,yRi);
                                        W1yo = besselk(1,yRo);
                                        W2xi = besselk(2,xRi);
                                        W2xo = besselk(2,xRo);
                                        W2yi = besselk(2,yRi);
                                        W2yo = besselk(2,yRo);
                                    end
                                    M(1,1) = -(y2-k2)*Ri2*Z1xi+2*Lambda1*xRi*Z2xi;
                                    M(1,2) = 2*k1*yRi*Ri*Z1yi-4*k1*Ri*Z2yi;
                                    M(1,3) = 2*Lambda2*yRi*Z2yi;
                                    M(1,4) = -(y2-k2)*Ri2*W1xi+2*xRi*W2xi;
                                    M(1,5) = 2*Lambda2*k1*yRi*Ri*W1yi-4*k1*Ri*W2yi;
                                    M(1,6) = 2*yRi*W2yi;
                                    M(2,1) = -2*Lambda1*xRi*Z2xi;
                                    M(2,2) = -k1*yRi*Ri*Z1yi+4*k1*Ri*Z2yi;
                                    M(2,3) = y2*Ri2*Z1yi-2*Lambda2*yRi*Z2yi;
                                    M(2,4) = -2*xRi*W2xi;
                                    M(2,5) = -Lambda2*k1*yRi*Ri*W1yi+4*k1*Ri*W2yi;
                                    M(2,6) = y2*Ri2*W1yi-2*yRi*W2yi;
                                    M(3,1) = 2*k1*Ri*Z1xi-2*Lambda1*k1*xRi*Ri*Z2xi;
                                    M(3,2) = yRi*Z1yi-(y2-k2)*Ri2*Z2yi;
                                    M(3,3) = -k1*Ri*Z1yi;
                                    M(3,4) = 2*k1*Ri*W1xi-2*k1*xRi*Ri*W2xi;
                                    M(3,5) = Lambda2*yRi*W1yi-(y2-k2)*Ri2*W2yi;
                                    M(3,6) = -k1*Ri*W1yi;
                                    M(4,1) = -(y2-k2)*Ro2*Z1xo+2*Lambda1*xRo*Z2xo;
                                    M(4,2) = 2*k1*yRo*Ro*Z1yo-4*k1*Ro*Z2yo;
                                    M(4,3) = 2*Lambda2*yRo*Z2yo;
                                    M(4,4) = -(y2-k2)*Ro2*W1xo+2*xRo*W2xo;
                                    M(4,5) = 2*Lambda2*k1*yRo*Ro*W1yo-4*k1*Ro*W2yo;
                                    M(4,6) = 2*yRo*W2yo;
                                    M(5,1) = -2*Lambda1*xRo*Z2xo;
                                    M(5,2) = -k1*yRo*Ro*Z1yo+4*k1*Ro*Z2yo;
                                    M(5,3) = y2*Ro2*Z1yo-2*Lambda2*yRo*Z2yo;
                                    M(5,4) = -2*xRo*W2xo;
                                    M(5,5) = -Lambda2*k1*yRo*Ro*W1yo+4*k1*Ro*W2yo;
                                    M(5,6) = y2*Ro2*W1yo-2*yRo*W2yo;
                                    M(6,1) = 2*k1*Ro*Z1xo-2*Lambda1*k1*xRo*Ro*Z2xo;
                                    M(6,2) = yRo*Z1yo-(y2-k2)*Ro2*Z2yo;
                                    M(6,3) = -k1*Ro*Z1yo;
                                    M(6,4) = 2*k1*Ro*W1xo-2*k1*xRo*Ro*W2xo;
                                    M(6,5) = Lambda2*yRo*W1yo-(y2-k2)*Ro2*W2yo;
                                    M(6,6) = -k1*Ro*W1yo;
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
                MissingModes(p+1) = 1;
                break
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end))
                X(end-MissingSamples:end) = [];
                Misses(end-MissingSamples:end) = 0;
                break
            end
            if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                send(Q1,[FrequencyRange(i),X(i)/1e3,p+1])
            elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                addpoints(g(p+1),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
        end
        if  all(X == 0)
            MissingModes(p+1) = 1;
        end
        if  ~MissingModes(p+1)
            X(Misses(1:length(X)) == 1) = NaN;
            F{p+1}(:,1) = FrequencyRange(1:length(X));
            F{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
            F{p+1}(:,3) = FrequencyRange(1:length(X))*(Ro-Ri);
            F{p+1}(:,4) = fillmissing(X,'spline')/1e3;
        else
            F{p+1} = F{p};
            X1{p+1}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(F{p+1}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop
                return
            end
            X1{p+1}(i,1) = 0;
            if  i == 1
                if  p == 1
                    SweepRange = [F{p+1}(MaxInd,1)-2*FrequencyOffset/(Ro-Ri)/1e3 F{p+1}(MaxInd,1)+10/(Ro-Ri)/1e3];
                else
                    SweepRange = [F{p+1}(MaxInd,1)-FrequencyOffset/(Ro-Ri)/1e3 F{p+1}(MaxInd,1)+10/(Ro-Ri)/1e3];
                end
            elseif i == 2
                if  p == 1
                    SweepRange = [X1{p+1}(i-1)-2*FrequencyOffset/(Ro-Ri)/1e3 X1{p+1}(i-1)+10/(Ro-Ri)/1e3];
                else
                    SweepRange = [X1{p+1}(i-1)-FrequencyOffset/(Ro-Ri)/1e3 X1{p+1}(i-1)+10/(Ro-Ri)/1e3];
                end
            else
                SweepRange = X1{p+1}(i-1)-2*abs(X1{p+1}(i-2)-X1{p+1}(i-1));
                if  X1{p+1}(i-1) < X1{p+1}(i-2)
                    SweepRange(2) = X1{p+1}(i-1)+2/(Ro-Ri)/1e3;
                else
                    SweepRange(2) = X1{p+1}(i-1)+5*abs(X1{p+1}(i-2)-X1{p+1}(i-1));
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
                                Z1xi = besselj(1,xRi);
                                Z1xo = besselj(1,xRo);
                                Z1yi = besselj(1,yRi);
                                Z1yo = besselj(1,yRo);
                                Z2xi = besselj(2,xRi);
                                Z2xo = besselj(2,xRo);
                                Z2yi = besselj(2,yRi);
                                Z2yo = besselj(2,yRo);
                                W1xi = bessely(1,xRi);
                                W1xo = bessely(1,xRo);
                                W1yi = bessely(1,yRi);
                                W1yo = bessely(1,yRo);
                                W2xi = bessely(2,xRi);
                                W2xo = bessely(2,xRo);
                                W2yi = bessely(2,yRi);
                                W2yo = bessely(2,yRo);
                            elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity
                                Lambda1 = -1;
                                Lambda2 = 1;
                                xRo = abs(x*Ro);
                                yRo = y*Ro;
                                xRi = abs(x*Ri);
                                yRi = y*Ri;
                                Z1xi = besseli(1,xRi);
                                Z1xo = besseli(1,xRo);
                                Z1yi = besselj(1,yRi);
                                Z1yo = besselj(1,yRo);
                                Z2xi = besseli(2,xRi);
                                Z2xo = besseli(2,xRo);
                                Z2yi = besselj(2,yRi);
                                Z2yo = besselj(2,yRo);
                                W1xi = besselk(1,xRi);
                                W1xo = besselk(1,xRo);
                                W1yi = bessely(1,yRi);
                                W1yo = bessely(1,yRo);
                                W2xi = besselk(2,xRi);
                                W2xo = besselk(2,xRo);
                                W2yi = bessely(2,yRi);
                                W2yo = bessely(2,yRo);
                            end
                            M(1,1) = -(y2-k2)*Ri2*Z1xi+2*Lambda1*xRi*Z2xi;
                            M(1,2) = 2*k1*yRi*Ri*Z1yi-4*k1*Ri*Z2yi;
                            M(1,3) = 2*Lambda2*yRi*Z2yi;
                            M(1,4) = -(y2-k2)*Ri2*W1xi+2*xRi*W2xi;
                            M(1,5) = 2*Lambda2*k1*yRi*Ri*W1yi-4*k1*Ri*W2yi;
                            M(1,6) = 2*yRi*W2yi;
                            M(2,1) = -2*Lambda1*xRi*Z2xi;
                            M(2,2) = -k1*yRi*Ri*Z1yi+4*k1*Ri*Z2yi;
                            M(2,3) = y2*Ri2*Z1yi-2*Lambda2*yRi*Z2yi;
                            M(2,4) = -2*xRi*W2xi;
                            M(2,5) = -Lambda2*k1*yRi*Ri*W1yi+4*k1*Ri*W2yi;
                            M(2,6) = y2*Ri2*W1yi-2*yRi*W2yi;
                            M(3,1) = 2*k1*Ri*Z1xi-2*Lambda1*k1*xRi*Ri*Z2xi;
                            M(3,2) = yRi*Z1yi-(y2-k2)*Ri2*Z2yi;
                            M(3,3) = -k1*Ri*Z1yi;
                            M(3,4) = 2*k1*Ri*W1xi-2*k1*xRi*Ri*W2xi;
                            M(3,5) = Lambda2*yRi*W1yi-(y2-k2)*Ri2*W2yi;
                            M(3,6) = -k1*Ri*W1yi;
                            M(4,1) = -(y2-k2)*Ro2*Z1xo+2*Lambda1*xRo*Z2xo;
                            M(4,2) = 2*k1*yRo*Ro*Z1yo-4*k1*Ro*Z2yo;
                            M(4,3) = 2*Lambda2*yRo*Z2yo;
                            M(4,4) = -(y2-k2)*Ro2*W1xo+2*xRo*W2xo;
                            M(4,5) = 2*Lambda2*k1*yRo*Ro*W1yo-4*k1*Ro*W2yo;
                            M(4,6) = 2*yRo*W2yo;
                            M(5,1) = -2*Lambda1*xRo*Z2xo;
                            M(5,2) = -k1*yRo*Ro*Z1yo+4*k1*Ro*Z2yo;
                            M(5,3) = y2*Ro2*Z1yo-2*Lambda2*yRo*Z2yo;
                            M(5,4) = -2*xRo*W2xo;
                            M(5,5) = -Lambda2*k1*yRo*Ro*W1yo+4*k1*Ro*W2yo;
                            M(5,6) = y2*Ro2*W1yo-2*yRo*W2yo;
                            M(6,1) = 2*k1*Ro*Z1xo-2*Lambda1*k1*xRo*Ro*Z2xo;
                            M(6,2) = yRo*Z1yo-(y2-k2)*Ro2*Z2yo;
                            M(6,3) = -k1*Ro*Z1yo;
                            M(6,4) = 2*k1*Ro*W1xo-2*k1*xRo*Ro*W2xo;
                            M(6,5) = Lambda2*yRo*W1yo-(y2-k2)*Ro2*W2yo;
                            M(6,6) = -k1*Ro*W1yo;
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
%                             z = isoutlier(vertcat(X1{p+1}(1:i-1),Frequency(2)),'movmedian',5,'ThresholdFactor',9);
%                             if  z(end) && abs(X1{p+1}(i-1)-Frequency(2)) > 1
%                                 Outlier = 1;
%                             else
%                                 Outlier = 0;
%                             end
%                         end
%                         if  ~Outlier
                            X1{p+1}(i,1) = Frequency(2);
                            break
%                         end
                    end
                end
                if  X1{p+1}(i) > 0
                    break
                end
            end
            if  X1{p+1}(i) == 0
                break
            end
            if  Multithreading
                send(Q2,[X1{p+1}(i),PhaseVelocityRange(i)/1e3,p+1])
            else
                addpoints(g1(p+1),X1{p+1}(i),PhaseVelocityRange(i)/1e3);
                drawnow limitrate
            end
        end
        X1{p+1}(:,2) = X1{p+1}(:,1)/1e3;
        X1{p+1}(:,3) = X1{p+1}(:,1)*(Ro-Ri);
        X1{p+1}(:,4) = PhaseVelocityRange(1:height(X1{p+1}))/1e3;
        if  X1{p+1}(end,1) == 0
            X1{p+1}(end,1) = H(p);
            X1{p+1}(end,2) = H(p)/1e3;
            X1{p+1}(end,3) = H(p)*(Ro-Ri);
            X1{p+1}(end,4) = PhaseVelocityLimit/1e3;
        end
        X1{p+1}(X1{p+1}(:,1) == 0,:) = [];
        X1{p+1} = flipud(X1{p+1});
    end
    X1(MissingModes == 1) = [];
    F(MissingModes == 1) = [];
    for p = 2:length(F)
        F{p}(F{p}(:,4) == 0,:) = [];
        F{p} = vertcat(X1{p},F{p});
        F{p}(:,6) = 0;
% F{p}(:,5) = smooth(((F{p}(:,4)).^2)./(F{p}(:,4)-F{p}(:,1).*differentiate(fit(F{p}(:,1),F{p}(:,4),'cubicspline'),F{p}(:,1))));
    end
end