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
function F = Computer_Isotropic_Pipe_F_UD(Multithreading,Q1,Q2,ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth)        
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
Densityi = InnerFluid.Density/Material.Density;
Y = zeros(5^PhaseVelocitySections+1,3);
for i = 1:length(FrequencyRange)
    if  Stop
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
    kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
    kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
    X(i,1) = 0;
    if  i == 1
        SweepRange = 1e-2:50;
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
        for j = 2:length(SweepRange)-1
            if  SweepRange(j-1) < InnerFluid.Velocity && SweepRange(j+1) > InnerFluid.Velocity
                SweepRange(j) = NaN;
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
                    k1 = AngularFrequency/PhaseVelocity(l);
                    k2 = k1^2;
                    x = sqrt(k2-kL2);
                    y = sqrt(k2-kT2);
                    y2 = kT2-k2;
                    xRo = x*Ro;
                    yRo = y*Ro;
                    xRi = x*Ri;
                    yRi = y*Ri;
                    Z1xi = besseli(1,xRi);
                    Z1xo = besseli(1,xRo);
                    Z1yi = besseli(1,yRi);
                    Z1yo = besseli(1,yRo);
                    W1xi = besselk(1,xRi);
                    W1xo = besselk(1,xRo);
                    W1yi = besselk(1,yRi);
                    W1yo = besselk(1,yRo);
                    dZ1xiRi = Z1xi+xRi*besseli(2,xRi);
                    dZ1xoRo = Z1xo+xRo*besseli(2,xRo);
                    dZ1yiRi = Z1yi+yRi*besseli(2,yRi);
                    dZ1yoRo = Z1yo+yRo*besseli(2,yRo);
                    dW1xiRi = W1xi-xRi*besselk(2,xRi);
                    dW1xoRo = W1xo-xRo*besselk(2,xRo);
                    dW1yiRi = W1yi-yRi*besselk(2,yRi);
                    dW1yoRo = W1yo-yRo*besselk(2,yRo);
                    zRi = sqrt(kInnerFluid2-k2)*Ri;
                    Z1zi = besselj(1,zRi);
                    dZ1ziRi = Z1zi-zRi*besselj(2,zRi);
                    M(1,1) = dZ1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*Z1xi;
                    M(1,2) = dW1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*W1xi;
                    M(1,3) = -k1*(dZ1yiRi+(y2*Ri2-1)*Z1yi);
                    M(1,4) = -k1*(dW1yiRi+(y2*Ri2-1)*W1yi);
                    M(1,5) = dZ1yiRi-Z1yi;
                    M(1,6) = dW1yiRi-W1yi;
                    M(1,7) = .5*kT2*Ri2*Densityi*Z1zi;
                    M(2,1) = 2*(dZ1xiRi-Z1xi);
                    M(2,2) = 2*(dW1xiRi-W1xi);
                    M(2,3) = 2*k1*(Z1yi-dZ1yiRi);
                    M(2,4) = 2*k1*(W1yi-dW1yiRi);
                    M(2,5) = 2*dZ1yiRi+(y2*Ri2-2)*Z1yi;
                    M(2,6) = 2*dW1yiRi+(y2*Ri2-2)*W1yi;
                    M(3,1) = 2*k1*dZ1xiRi;
                    M(3,2) = 2*k1*dW1xiRi;
                    M(3,3) = (y2-k2)*dZ1yiRi;
                    M(3,4) = (y2-k2)*dW1yiRi;
                    M(3,5) = -k1*Z1yi;
                    M(3,6) = -k1*W1yi;
                    M(4,1) = dZ1xiRi;
                    M(4,2) = dW1xiRi;
                    M(4,3) = -k1*dZ1yiRi;
                    M(4,4) = -k1*dW1yiRi;
                    M(4,5) = -Z1yi;
                    M(4,6) = -W1yi;
                    M(4,7) = dZ1ziRi;
                    M(5,1) = dZ1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*Z1xo;
                    M(5,2) = dW1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*W1xo;
                    M(5,3) = -k1*(dZ1yoRo+(y2*Ro2-1)*Z1yo);
                    M(5,4) = -k1*(dW1yoRo+(y2*Ro2-1)*W1yo);
                    M(5,5) = dZ1yoRo-Z1yo;
                    M(5,6) = dW1yoRo-W1yo;
                    M(6,1) = 2*(dZ1xoRo-Z1xo);
                    M(6,2) = 2*(dW1xoRo-W1xo);
                    M(6,3) = 2*k1*(Z1yo-dZ1yoRo);
                    M(6,4) = 2*k1*(W1yo-dW1yoRo);
                    M(6,5) = 2*dZ1yoRo+(y2*Ro2-2)*Z1yo;
                    M(6,6) = 2*dW1yoRo+(y2*Ro2-2)*W1yo;
                    M(7,1) = 2*k1*dZ1xoRo;
                    M(7,2) = 2*k1*dW1xoRo;
                    M(7,3) = (y2-k2)*dZ1yoRo;
                    M(7,4) = (y2-k2)*dW1yoRo;
                    M(7,5) = -k1*Z1yo;
                    M(7,6) = -k1*W1yo;
                    Y(j,l) = real(det(M/Z1zi));
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
            kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
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
                        SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)+10*abs(X(i-2)-X(i-1))];
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
                    if  q == 1 && o == 0 && (any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) || (SweepRange(1) > InnerFluid.Velocity && SweepRange(2) < InnerFluid.Velocity) || (SweepRange(1) > Material.LongitudinalVelocity && SweepRange(2) < Material.LongitudinalVelocity)) % Neighbors cannot be excluded because SweepRange has only 2 elements
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
                        if  (SweepRange(j) > InnerFluid.Velocity && SweepRange(j+1) < InnerFluid.Velocity) || (SweepRange(j) > Material.LongitudinalVelocity && SweepRange(j+1) < Material.LongitudinalVelocity)
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
                                    if  SweepRange(j) > Material.LongitudinalVelocity
                                        x = sqrt(kL2-k2);
                                        y = sqrt(kT2-k2);
                                    elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                        x = sqrt(k2-kL2);
                                        y = sqrt(kT2-k2);
                                    elseif SweepRange(j) < Material.TransverseVelocity
                                        x = sqrt(k2-kL2);
                                        y = sqrt(k2-kT2);
                                    end     
                                    y2 = kT2-k2;
                                    xRo = x*Ro;
                                    yRo = y*Ro;
                                    xRi = x*Ri;
                                    yRi = y*Ri;
                                    if  SweepRange(j) > Material.LongitudinalVelocity
                                        Z1xi = besselj(1,xRi);
                                        Z1xo = besselj(1,xRo);
                                        Z1yi = besselj(1,yRi);
                                        Z1yo = besselj(1,yRo);
                                        W1xi = bessely(1,xRi);
                                        W1xo = bessely(1,xRo);
                                        W1yi = bessely(1,yRi);
                                        W1yo = bessely(1,yRo);
                                        dZ1xiRi = Z1xi-xRi*besselj(2,xRi);
                                        dZ1xoRo = Z1xo-xRo*besselj(2,xRo);
                                        dZ1yiRi = Z1yi-yRi*besselj(2,yRi);
                                        dZ1yoRo = Z1yo-yRo*besselj(2,yRo);
                                        dW1xiRi = W1xi-xRi*bessely(2,xRi);
                                        dW1xoRo = W1xo-xRo*bessely(2,xRo);
                                        dW1yiRi = W1yi-yRi*bessely(2,yRi);
                                        dW1yoRo = W1yo-yRo*bessely(2,yRo);
                                    elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                        Z1xi = besseli(1,xRi);
                                        Z1xo = besseli(1,xRo);
                                        Z1yi = besselj(1,yRi);
                                        Z1yo = besselj(1,yRo);
                                        W1xi = besselk(1,xRi);
                                        W1xo = besselk(1,xRo);
                                        W1yi = bessely(1,yRi);
                                        W1yo = bessely(1,yRo);
                                        dZ1xiRi = Z1xi+xRi*besseli(2,xRi);
                                        dZ1xoRo = Z1xo+xRo*besseli(2,xRo);
                                        dZ1yiRi = Z1yi-yRi*besselj(2,yRi);
                                        dZ1yoRo = Z1yo-yRo*besselj(2,yRo);
                                        dW1xiRi = W1xi-xRi*besselk(2,xRi);
                                        dW1xoRo = W1xo-xRo*besselk(2,xRo);
                                        dW1yiRi = W1yi-yRi*bessely(2,yRi);
                                        dW1yoRo = W1yo-yRo*bessely(2,yRo);
                                    elseif SweepRange(j) < Material.TransverseVelocity
                                        Z1xi = besseli(1,xRi);
                                        Z1xo = besseli(1,xRo);
                                        Z1yi = besseli(1,yRi);
                                        Z1yo = besseli(1,yRo);
                                        W1xi = besselk(1,xRi);
                                        W1xo = besselk(1,xRo);
                                        W1yi = besselk(1,yRi);
                                        W1yo = besselk(1,yRo);
                                        dZ1xiRi = Z1xi+xRi*besseli(2,xRi);
                                        dZ1xoRo = Z1xo+xRo*besseli(2,xRo);
                                        dZ1yiRi = Z1yi+yRi*besseli(2,yRi);
                                        dZ1yoRo = Z1yo+yRo*besseli(2,yRo);
                                        dW1xiRi = W1xi-xRi*besselk(2,xRi);
                                        dW1xoRo = W1xo-xRo*besselk(2,xRo);
                                        dW1yiRi = W1yi-yRi*besselk(2,yRi);
                                        dW1yoRo = W1yo-yRo*besselk(2,yRo);
                                    end
                                    zRi = sqrt(kInnerFluid2-k2)*Ri;
                                    Z1zi = besselj(1,zRi);
                                    dZ1ziRi = Z1zi-zRi*besselj(2,zRi);
                                    M(1,1) = dZ1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*Z1xi;
                                    M(1,2) = dW1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*W1xi;
                                    M(1,3) = -k1*(dZ1yiRi+(y2*Ri2-1)*Z1yi);
                                    M(1,4) = -k1*(dW1yiRi+(y2*Ri2-1)*W1yi);
                                    M(1,5) = dZ1yiRi-Z1yi;
                                    M(1,6) = dW1yiRi-W1yi;
                                    M(1,7) = .5*kT2*Ri2*Densityi*Z1zi;
                                    M(2,1) = 2*(dZ1xiRi-Z1xi);
                                    M(2,2) = 2*(dW1xiRi-W1xi);
                                    M(2,3) = 2*k1*(Z1yi-dZ1yiRi);
                                    M(2,4) = 2*k1*(W1yi-dW1yiRi);
                                    M(2,5) = 2*dZ1yiRi+(y2*Ri2-2)*Z1yi;
                                    M(2,6) = 2*dW1yiRi+(y2*Ri2-2)*W1yi;
                                    M(3,1) = 2*k1*dZ1xiRi;
                                    M(3,2) = 2*k1*dW1xiRi;
                                    M(3,3) = (y2-k2)*dZ1yiRi;
                                    M(3,4) = (y2-k2)*dW1yiRi;
                                    M(3,5) = -k1*Z1yi;
                                    M(3,6) = -k1*W1yi;
                                    M(4,1) = dZ1xiRi;
                                    M(4,2) = dW1xiRi;
                                    M(4,3) = -k1*dZ1yiRi;
                                    M(4,4) = -k1*dW1yiRi;
                                    M(4,5) = -Z1yi;
                                    M(4,6) = -W1yi;
                                    M(4,7) = dZ1ziRi;
                                    M(5,1) = dZ1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*Z1xo;
                                    M(5,2) = dW1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*W1xo;
                                    M(5,3) = -k1*(dZ1yoRo+(y2*Ro2-1)*Z1yo);
                                    M(5,4) = -k1*(dW1yoRo+(y2*Ro2-1)*W1yo);
                                    M(5,5) = dZ1yoRo-Z1yo;
                                    M(5,6) = dW1yoRo-W1yo;
                                    M(6,1) = 2*(dZ1xoRo-Z1xo);
                                    M(6,2) = 2*(dW1xoRo-W1xo);
                                    M(6,3) = 2*k1*(Z1yo-dZ1yoRo);
                                    M(6,4) = 2*k1*(W1yo-dW1yoRo);
                                    M(6,5) = 2*dZ1yoRo+(y2*Ro2-2)*Z1yo;
                                    M(6,6) = 2*dW1yoRo+(y2*Ro2-2)*W1yo;
                                    M(7,1) = 2*k1*dZ1xoRo;
                                    M(7,2) = 2*k1*dW1xoRo;
                                    M(7,3) = (y2-k2)*dZ1yoRo;
                                    M(7,4) = (y2-k2)*dW1yoRo;
                                    M(7,5) = -k1*Z1yo;
                                    M(7,6) = -k1*W1yo;
                                    Y(j,l) = real(det(M)/Z1zi);
                                end
                            end
                            if  k == 1
                                Y(j+1,1) = Y(j,3);
                            end
                            if  (j == 1 && abs(Y(j,2)) < 1e22) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
                                PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                                Y(j,3) = Y(j,2);
                            elseif (j == 1 && abs(Y(j,2)) < 1e22) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
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
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
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
                    SweepRange = [F{p+1}(MaxInd,1)-2/4*FrequencyOffset/(Ro-Ri)/1e3 F{p+1}(MaxInd,1)+10/4/(Ro-Ri)/1e3];
                else
                    SweepRange = [F{p+1}(MaxInd,1)-1/4*FrequencyOffset/(Ro-Ri)/1e3 F{p+1}(MaxInd,1)+10/4/(Ro-Ri)/1e3];
                end
            elseif i == 2
                if  p == 1
                    SweepRange = [X1{p+1}(i-1)-2/4*FrequencyOffset/(Ro-Ri)/1e3 X1{p+1}(i-1)+10/4/(Ro-Ri)/1e3];
                else
                    SweepRange = [X1{p+1}(i-1)-1/4*FrequencyOffset/(Ro-Ri)/1e3 X1{p+1}(i-1)+10/4/(Ro-Ri)/1e3];
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
                            kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
                            kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
                            if  PhaseVelocityRange(i) > Material.LongitudinalVelocity
                                x = sqrt(kL2-k2);
                                y = sqrt(kT2-k2);
                            elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity
                                x = sqrt(k2-kL2);
                                y = sqrt(kT2-k2);
                            end
                            y2 = kT2-k2;
                            xRo = x*Ro;
                            yRo = y*Ro;
                            xRi = x*Ri;
                            yRi = y*Ri;
                            if  PhaseVelocityRange(i) > Material.LongitudinalVelocity
                                Z1xi = besselj(1,xRi);
                                Z1xo = besselj(1,xRo);
                                Z1yi = besselj(1,yRi);
                                Z1yo = besselj(1,yRo);
                                W1xi = bessely(1,xRi);
                                W1xo = bessely(1,xRo);
                                W1yi = bessely(1,yRi);
                                W1yo = bessely(1,yRo);
                                dZ1xiRi = Z1xi-xRi*besselj(2,xRi);
                                dZ1xoRo = Z1xo-xRo*besselj(2,xRo);
                                dZ1yiRi = Z1yi-yRi*besselj(2,yRi);
                                dZ1yoRo = Z1yo-yRo*besselj(2,yRo);
                                dW1xiRi = W1xi-xRi*bessely(2,xRi);
                                dW1xoRo = W1xo-xRo*bessely(2,xRo);
                                dW1yiRi = W1yi-yRi*bessely(2,yRi);
                                dW1yoRo = W1yo-yRo*bessely(2,yRo);
                            elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity
                                Z1xi = besseli(1,xRi);
                                Z1xo = besseli(1,xRo);
                                Z1yi = besselj(1,yRi);
                                Z1yo = besselj(1,yRo);
                                W1xi = besselk(1,xRi);
                                W1xo = besselk(1,xRo);
                                W1yi = bessely(1,yRi);
                                W1yo = bessely(1,yRo);
                                dZ1xiRi = Z1xi+xRi*besseli(2,xRi);
                                dZ1xoRo = Z1xo+xRo*besseli(2,xRo);
                                dZ1yiRi = Z1yi-yRi*besselj(2,yRi);
                                dZ1yoRo = Z1yo-yRo*besselj(2,yRo);
                                dW1xiRi = W1xi-xRi*besselk(2,xRi);
                                dW1xoRo = W1xo-xRo*besselk(2,xRo);
                                dW1yiRi = W1yi-yRi*bessely(2,yRi);
                                dW1yoRo = W1yo-yRo*bessely(2,yRo);
                            end
                            zRi = sqrt((AngularFrequency/InnerFluid.Velocity)^2-k2)*Ri;
                            Z1zi = besselj(1,zRi);
                            dZ1ziRi = Z1zi-zRi*besselj(2,zRi);
                            M(1,1) = dZ1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*Z1xi;
                            M(1,2) = dW1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*W1xi;
                            M(1,3) = -k1*(dZ1yiRi+(y2*Ri2-1)*Z1yi);
                            M(1,4) = -k1*(dW1yiRi+(y2*Ri2-1)*W1yi);
                            M(1,5) = dZ1yiRi-Z1yi;
                            M(1,6) = dW1yiRi-W1yi;
                            M(1,7) = .5*kT2*Ri2*Densityi*Z1zi;
                            M(2,1) = 2*(dZ1xiRi-Z1xi);
                            M(2,2) = 2*(dW1xiRi-W1xi);
                            M(2,3) = 2*k1*(Z1yi-dZ1yiRi);
                            M(2,4) = 2*k1*(W1yi-dW1yiRi);
                            M(2,5) = 2*dZ1yiRi+(y2*Ri2-2)*Z1yi;
                            M(2,6) = 2*dW1yiRi+(y2*Ri2-2)*W1yi;
                            M(3,1) = 2*k1*dZ1xiRi;
                            M(3,2) = 2*k1*dW1xiRi;
                            M(3,3) = (y2-k2)*dZ1yiRi;
                            M(3,4) = (y2-k2)*dW1yiRi;
                            M(3,5) = -k1*Z1yi;
                            M(3,6) = -k1*W1yi;
                            M(4,1) = dZ1xiRi;
                            M(4,2) = dW1xiRi;
                            M(4,3) = -k1*dZ1yiRi;
                            M(4,4) = -k1*dW1yiRi;
                            M(4,5) = -Z1yi;
                            M(4,6) = -W1yi;
                            M(4,7) = dZ1ziRi;
                            M(5,1) = dZ1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*Z1xo;
                            M(5,2) = dW1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*W1xo;
                            M(5,3) = -k1*(dZ1yoRo+(y2*Ro2-1)*Z1yo);
                            M(5,4) = -k1*(dW1yoRo+(y2*Ro2-1)*W1yo);
                            M(5,5) = dZ1yoRo-Z1yo;
                            M(5,6) = dW1yoRo-W1yo;
                            M(6,1) = 2*(dZ1xoRo-Z1xo);
                            M(6,2) = 2*(dW1xoRo-W1xo);
                            M(6,3) = 2*k1*(Z1yo-dZ1yoRo);
                            M(6,4) = 2*k1*(W1yo-dW1yoRo);
                            M(6,5) = 2*dZ1yoRo+(y2*Ro2-2)*Z1yo;
                            M(6,6) = 2*dW1yoRo+(y2*Ro2-2)*W1yo;
                            M(7,1) = 2*k1*dZ1xoRo;
                            M(7,2) = 2*k1*dW1xoRo;
                            M(7,3) = (y2-k2)*dZ1yoRo;
                            M(7,4) = (y2-k2)*dW1yoRo;
                            M(7,5) = -k1*Z1yo;
                            M(7,6) = -k1*W1yo;
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