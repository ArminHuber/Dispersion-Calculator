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
function L = Computer_Isotropic_Pipe_L(Multithreading,Q1,Q2,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
L{1} = [];
if  ~Multithreading
    if  HigherOrderModes && any(H)
        for p = 1:length(H)+1
            g(p) = animatedline(ax,'color','r');
            g1(p) = animatedline(ax,'color','r');
        end
    else
        g = animatedline(ax,'color','r');
    end
end
Ro2 = Ro^2;
Ri2 = Ri^2;
Y = zeros(5^PhaseVelocitySections+1,3);
X = Material.CylinderVelocity;
if  Multithreading
    send(Q1,[FrequencyRange(1),X(1)/1e3,1])
else
    addpoints(g(1),FrequencyRange(1),X(1)/1e3);
    drawnow limitrate
end
for i = 2:length(FrequencyRange)
    if  Stop
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
    kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
    X(i,1) = 0;
    if  i == 2
        SweepRange = X(i-1)+5:-1:X(i-1)-5;
    elseif i == 3
        SweepRange = [X(i-1)+1*abs(X(i-2)-X(i-1)) X(i-1)-20*abs((X(i-2)-X(i-1)))];
    else
        SweepRange = [X(i-1)+10*abs(X(i-2)-X(i-1)) X(i-1)-10*abs((X(i-2)-X(i-1)))];  
    end
    if  SweepRange(1) > Material.CylinderVelocity
        SweepRange(1) = Material.CylinderVelocity;
    end
    if  SweepRange(end) < 0
        SweepRange(end) = 0;
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
                    if  SweepRange(j) > Material.TransverseVelocity
                        Lambda1 = -1;
                        Lambda2 = 1;
                        xRo = abs(x*Ro);
                        yRo = y*Ro;
                        xRi = abs(x*Ri);
                        yRi = y*Ri;
                        Z0xi = besseli(0,xRi);
                        Z0xo = besseli(0,xRo);
                        Z0yi = besselj(0,yRi);
                        Z0yo = besselj(0,yRo);
                        Z1xi = besseli(1,xRi);
                        Z1xo = besseli(1,xRo);
                        Z1yi = besselj(1,yRi);
                        Z1yo = besselj(1,yRo);
                        W0xi = besselk(0,xRi);
                        W0xo = besselk(0,xRo);
                        W0yi = bessely(0,yRi);
                        W0yo = bessely(0,yRo);
                        W1xi = besselk(1,xRi);
                        W1xo = besselk(1,xRo);
                        W1yi = bessely(1,yRi);
                        W1yo = bessely(1,yRo);
                    elseif SweepRange(j) < Material.TransverseVelocity
                        Lambda1 = -1;
                        Lambda2 = -1;
                        xRo = abs(x*Ro);
                        yRo = abs(y*Ro);
                        xRi = abs(x*Ri);
                        yRi = abs(y*Ri);
                        Z0xi = besseli(0,xRi);
                        Z0xo = besseli(0,xRo);
                        Z0yi = besseli(0,yRi);
                        Z0yo = besseli(0,yRo);
                        Z1xi = besseli(1,xRi);
                        Z1xo = besseli(1,xRo);
                        Z1yi = besseli(1,yRi);
                        Z1yo = besseli(1,yRo);
                        W0xi = besselk(0,xRi);
                        W0xo = besselk(0,xRo);
                        W0yi = besselk(0,yRi);
                        W0yo = besselk(0,yRo);
                        W1xi = besselk(1,xRi);
                        W1xo = besselk(1,xRo);
                        W1yi = besselk(1,yRi);
                        W1yo = besselk(1,yRo);
                    end
                    M(1,1) = -(y2-k2)*Ri2*Z0xi+2*Lambda1*xRi*Z1xi;
                    M(1,2) = 2*k1*Ri*(yRi*Z0yi-Z1yi);
                    M(1,3) = -(y2-k2)*Ri2*W0xi+2*xRi*W1xi;
                    M(1,4) = 2*k1*Ri*(Lambda2*yRi*W0yi-W1yi);
                    M(2,1) = -2*k1*xRi*Lambda1*Ri*Z1xi;
                    M(2,2) = -(y2-k2)*Ri2*Z1yi;
                    M(2,3) = -2*k1*Ri*xRi*W1xi;
                    M(2,4) = -(y2-k2)*Ri2*W1yi;
                    M(3,1) = -(y2-k2)*Ro2*Z0xo+2*Lambda1*xRo*Z1xo;
                    M(3,2) = 2*k1*Ro*(yRo*Z0yo-Z1yo);
                    M(3,3) = -(y2-k2)*Ro2*W0xo+2*xRo*W1xo;
                    M(3,4) = 2*k1*Ro*(Lambda2*yRo*W0yo-W1yo);
                    M(4,1) = -2*k1*xRo*Lambda1*Ro*Z1xo;
                    M(4,2) = -(y2-k2)*Ro2*Z1yo;
                    M(4,3) = -2*k1*Ro*xRo*W1xo;
                    M(4,4) = -(y2-k2)*Ro2*W1yo;
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
L{1}(:,1) = FrequencyRange(1:length(X));
L{1}(:,2) = FrequencyRange(1:length(X))/1e3;
L{1}(:,3) = FrequencyRange(1:length(X))*(Ro-Ri);
L{1}(:,4) = fillmissing(X,'spline')/1e3;
L{1}(:,6) = 0;
% L{1}(:,5) = smooth(((L{1}(:,4)).^2)./(L{1}(:,4)-L{1}(:,1).*differentiate(fit(L{1}(:,1),L{1}(:,4),'cubicspline'),L{1}(:,1))));
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
            for j = 1:length(L)
                if  i <= height(L{j})
                    Neighbors(j) = L{j}(i,4)*1e3;
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
                    if  q == 1 && o == 0 && any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) % Neighbors cannot be excluded because SweepRange has only 2 elements
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
                                        Z0xi = besselj(0,xRi);
                                        Z0xo = besselj(0,xRo);
                                        Z0yi = besselj(0,yRi);
                                        Z0yo = besselj(0,yRo);
                                        Z1xi = besselj(1,xRi);
                                        Z1xo = besselj(1,xRo);
                                        Z1yi = besselj(1,yRi);
                                        Z1yo = besselj(1,yRo);
                                        W0xi = bessely(0,xRi);
                                        W0xo = bessely(0,xRo);
                                        W0yi = bessely(0,yRi);
                                        W0yo = bessely(0,yRo);
                                        W1xi = bessely(1,xRi);
                                        W1xo = bessely(1,xRo);
                                        W1yi = bessely(1,yRi);
                                        W1yo = bessely(1,yRo);
                                    elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                        Lambda1 = -1;
                                        Lambda2 = 1;
                                        xRo = abs(x*Ro);
                                        yRo = y*Ro;
                                        xRi = abs(x*Ri);
                                        yRi = y*Ri;
                                        Z0xi = besseli(0,xRi);
                                        Z0xo = besseli(0,xRo);
                                        Z0yi = besselj(0,yRi);
                                        Z0yo = besselj(0,yRo);
                                        Z1xi = besseli(1,xRi);
                                        Z1xo = besseli(1,xRo);
                                        Z1yi = besselj(1,yRi);
                                        Z1yo = besselj(1,yRo);
                                        W0xi = besselk(0,xRi);
                                        W0xo = besselk(0,xRo);
                                        W0yi = bessely(0,yRi);
                                        W0yo = bessely(0,yRo);
                                        W1xi = besselk(1,xRi);
                                        W1xo = besselk(1,xRo);
                                        W1yi = bessely(1,yRi);
                                        W1yo = bessely(1,yRo);
                                    elseif SweepRange(j) < Material.TransverseVelocity
                                        Lambda1 = -1;
                                        Lambda2 = -1;
                                        xRo = abs(x*Ro);
                                        yRo = abs(y*Ro);
                                        xRi = abs(x*Ri);
                                        yRi = abs(y*Ri);
                                        Z0xi = besseli(0,xRi);
                                        Z0xo = besseli(0,xRo);
                                        Z0yi = besseli(0,yRi);
                                        Z0yo = besseli(0,yRo);
                                        Z1xi = besseli(1,xRi);
                                        Z1xo = besseli(1,xRo);
                                        Z1yi = besseli(1,yRi);
                                        Z1yo = besseli(1,yRo);
                                        W0xi = besselk(0,xRi);
                                        W0xo = besselk(0,xRo);
                                        W0yi = besselk(0,yRi);
                                        W0yo = besselk(0,yRo);
                                        W1xi = besselk(1,xRi);
                                        W1xo = besselk(1,xRo);
                                        W1yi = besselk(1,yRi);
                                        W1yo = besselk(1,yRo);
                                    end
                                    M(1,1) = -(y2-k2)*Ri2*Z0xi+2*Lambda1*xRi*Z1xi;
                                    M(1,2) = 2*k1*Ri*(yRi*Z0yi-Z1yi);
                                    M(1,3) = -(y2-k2)*Ri2*W0xi+2*xRi*W1xi;
                                    M(1,4) = 2*k1*Ri*(Lambda2*yRi*W0yi-W1yi);
                                    M(2,1) = -2*k1*xRi*Lambda1*Ri*Z1xi;
                                    M(2,2) = -(y2-k2)*Ri2*Z1yi;
                                    M(2,3) = -2*k1*Ri*xRi*W1xi;
                                    M(2,4) = -(y2-k2)*Ri2*W1yi;
                                    M(3,1) = -(y2-k2)*Ro2*Z0xo+2*Lambda1*xRo*Z1xo;
                                    M(3,2) = 2*k1*Ro*(yRo*Z0yo-Z1yo);
                                    M(3,3) = -(y2-k2)*Ro2*W0xo+2*xRo*W1xo;
                                    M(3,4) = 2*k1*Ro*(Lambda2*yRo*W0yo-W1yo);
                                    M(4,1) = -2*k1*xRo*Lambda1*Ro*Z1xo;
                                    M(4,2) = -(y2-k2)*Ro2*Z1yo;
                                    M(4,3) = -2*k1*Ro*xRo*W1xo;
                                    M(4,4) = -(y2-k2)*Ro2*W1yo;
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
            L{p+1}(:,1) = FrequencyRange(1:length(X));
            L{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
            L{p+1}(:,3) = FrequencyRange(1:length(X))*(Ro-Ri);
            L{p+1}(:,4) = fillmissing(X,'spline')/1e3;
        else
            L{p+1} = L{p};
            X1{p+1}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(L{p+1}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop
                return
            end
            X1{p+1}(i,1) = 0;
            if  i == 1
                if  p == 1
                    SweepRange = [L{p+1}(MaxInd,1)-2*FrequencyOffset/(Ro-Ri)/1e3 L{p+1}(MaxInd,1)+10/(Ro-Ri)/1e3];
                else
                    SweepRange = [L{p+1}(MaxInd,1)-FrequencyOffset/(Ro-Ri)/1e3 L{p+1}(MaxInd,1)+10/(Ro-Ri)/1e3];
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
                                Z0xi = besselj(0,xRi);
                                Z0xo = besselj(0,xRo);
                                Z0yi = besselj(0,yRi);
                                Z0yo = besselj(0,yRo);
                                Z1xi = besselj(1,xRi);
                                Z1xo = besselj(1,xRo);
                                Z1yi = besselj(1,yRi);
                                Z1yo = besselj(1,yRo);
                                W0xi = bessely(0,xRi);
                                W0xo = bessely(0,xRo);
                                W0yi = bessely(0,yRi);
                                W0yo = bessely(0,yRo);
                                W1xi = bessely(1,xRi);
                                W1xo = bessely(1,xRo);
                                W1yi = bessely(1,yRi);
                                W1yo = bessely(1,yRo);
                            elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity
                                Lambda1 = -1;
                                Lambda2 = 1;
                                xRo = abs(x*Ro);
                                yRo = y*Ro;
                                xRi = abs(x*Ri);
                                yRi = y*Ri;
                                Z0xi = besseli(0,xRi);
                                Z0xo = besseli(0,xRo);
                                Z0yi = besselj(0,yRi);
                                Z0yo = besselj(0,yRo);
                                Z1xi = besseli(1,xRi);
                                Z1xo = besseli(1,xRo);
                                Z1yi = besselj(1,yRi);
                                Z1yo = besselj(1,yRo);
                                W0xi = besselk(0,xRi);
                                W0xo = besselk(0,xRo);
                                W0yi = bessely(0,yRi);
                                W0yo = bessely(0,yRo);
                                W1xi = besselk(1,xRi);
                                W1xo = besselk(1,xRo);
                                W1yi = bessely(1,yRi);
                                W1yo = bessely(1,yRo);
                            end
                            M(1,1) = -(y2-k2)*Ri2*Z0xi+2*Lambda1*xRi*Z1xi;
                            M(1,2) = 2*k1*Ri*(yRi*Z0yi-Z1yi);
                            M(1,3) = -(y2-k2)*Ri2*W0xi+2*xRi*W1xi;
                            M(1,4) = 2*k1*Ri*(Lambda2*yRi*W0yi-W1yi);
                            M(2,1) = -2*k1*xRi*Lambda1*Ri*Z1xi;
                            M(2,2) = -(y2-k2)*Ri2*Z1yi;
                            M(2,3) = -2*k1*Ri*xRi*W1xi;
                            M(2,4) = -(y2-k2)*Ri2*W1yi;
                            M(3,1) = -(y2-k2)*Ro2*Z0xo+2*Lambda1*xRo*Z1xo;
                            M(3,2) = 2*k1*Ro*(yRo*Z0yo-Z1yo);
                            M(3,3) = -(y2-k2)*Ro2*W0xo+2*xRo*W1xo;
                            M(3,4) = 2*k1*Ro*(Lambda2*yRo*W0yo-W1yo);
                            M(4,1) = -2*k1*xRo*Lambda1*Ro*Z1xo;
                            M(4,2) = -(y2-k2)*Ro2*Z1yo;
                            M(4,3) = -2*k1*Ro*xRo*W1xo;
                            M(4,4) = -(y2-k2)*Ro2*W1yo;
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
    L(MissingModes == 1) = [];
    for p = 2:length(L)
        L{p}(L{p}(:,4) == 0,:) = [];
        L{p} = vertcat(X1{p},L{p});
        L{p}(:,6) = 0;
% L{p}(:,5) = smooth(((L{p}(:,4)).^2)./(L{p}(:,4)-L{p}(:,1).*differentiate(fit(L{p}(:,1),L{p}(:,4),'cubicspline'),L{p}(:,1))));
    end
end