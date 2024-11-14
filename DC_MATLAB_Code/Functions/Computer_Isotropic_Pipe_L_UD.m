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
function L = Computer_Isotropic_Pipe_L_UD(Multithreading,Q1,Q2,ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
L{1} = [];
if  ~Multithreading
    if  HigherOrderModes && any(H)
        for p = 1:length(H)+2
            g(p) = animatedline(ax,'color','r');
            g1(p) = animatedline(ax,'color','r');
        end
    else
        g = animatedline(ax,'color','r');
        g(2) = animatedline(ax,'color','r');
    end
end
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityi = InnerFluid.Density/Material.Density;
AngularFrequency = 2*pi*FrequencyRange(1)*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
SweepRange = 50:10:Material.CylinderVelocity+500;
XRough = [];
for i = 1:length(SweepRange)
    k1 = AngularFrequency/SweepRange(i);
    k2 = k1^2;
    if  SweepRange(i) > Material.TransverseVelocity
        x = sqrt(k2-kL2);
        y = sqrt(kT2-k2);
    elseif SweepRange(i) < Material.TransverseVelocity
        x = sqrt(k2-kL2);
        y = sqrt(k2-kT2);
    end     
    y2 = kT2-k2;
    xRo = x*Ro;
    yRo = y*Ro;
    xRi = x*Ri;
    yRi = y*Ri;
    if  SweepRange(i) > Material.TransverseVelocity
        Z0xi = besseli(0,xRi);
        Z0xo = besseli(0,xRo);
        Z0yi = besselj(0,yRi);
        Z0yo = besselj(0,yRo);
        W0xi = besselk(0,xRi);
        W0xo = besselk(0,xRo);
        W0yi = bessely(0,yRi);
        W0yo = bessely(0,yRo);
        Z1xiRi = -xRi*besseli(1,xRi);
        Z1xoRo = -xRo*besseli(1,xRo);
        Z1yiRi = yRi*besselj(1,yRi);
        Z1yoRo = yRo*besselj(1,yRo);
        W1xiRi = xRi*besselk(1,xRi);
        W1xoRo = xRo*besselk(1,xRo);
        W1yiRi = yRi*bessely(1,yRi);
        W1yoRo = yRo*bessely(1,yRo);
    elseif SweepRange(i) < Material.TransverseVelocity
        Z0xi = besseli(0,xRi);
        Z0xo = besseli(0,xRo);
        Z0yi = besseli(0,yRi);
        Z0yo = besseli(0,yRo);
        W0xi = besselk(0,xRi);
        W0xo = besselk(0,xRo);
        W0yi = besselk(0,yRi);
        W0yo = besselk(0,yRo);
        Z1xiRi = -xRi*besseli(1,xRi);
        Z1xoRo = -xRo*besseli(1,xRo);
        Z1yiRi = -yRi*besseli(1,yRi);
        Z1yoRo = -yRo*besseli(1,yRo);
        W1xiRi = xRi*besselk(1,xRi);
        W1xoRo = xRo*besselk(1,xRo);
        W1yiRi = yRi*besselk(1,yRi);
        W1yoRo = yRo*besselk(1,yRo);
    end
    zRi = sqrt(kInnerFluid2-k2)*Ri;
    Z0zi = besselj(0,zRi);
    Z1ziRi = -zRi*besselj(1,zRi);
    M(1,1) = Ri2*(.5*kT2-k2)*Z0xi-Z1xiRi;
    M(1,2) = Ri2*(.5*kT2-k2)*W0xi-W1xiRi;
    M(1,3) = k1*(Z1yiRi-y2*Ri2*Z0yi);
    M(1,4) = k1*(W1yiRi-y2*Ri2*W0yi);
    M(1,5) = .5*kT2*Ri2*Densityi*Z0zi;
    M(2,1) = -2*k1*Z1xiRi;
    M(2,2) = -2*k1*W1xiRi;
    M(2,3) = (k2-y2)*Z1yiRi;
    M(2,4) = (k2-y2)*W1yiRi;
    M(3,1) = -Z1xiRi;
    M(3,2) = -W1xiRi;
    M(3,3) = k1*Z1yiRi;
    M(3,4) = k1*W1yiRi;
    M(3,5) = Z1ziRi;
    M(4,1) = Ro2*(.5*kT2-k2)*Z0xo-Z1xoRo;
    M(4,2) = Ro2*(.5*kT2-k2)*W0xo-W1xoRo;
    M(4,3) = k1*(Z1yoRo-y2*Ro2*Z0yo);
    M(4,4) = k1*(W1yoRo-y2*Ro2*W0yo);
    M(5,1) = -2*k1*Z1xoRo;
    M(5,2) = -2*k1*W1xoRo;
    M(5,3) = (k2-y2)*Z1yoRo;
    M(5,4) = (k2-y2)*W1yoRo;
    Y(i) = abs(det(M));
    if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i) && ~(SweepRange(i-2) < Material.TransverseVelocity && SweepRange(i) > Material.TransverseVelocity)
        XRough(end+1) = SweepRange(i-1);
    end
    if  length(XRough) == 2
        break
    end
end
Bisections = ceil(log2(1e-6/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
for j = 1:length(XRough)
    PhaseVelocity = [XRough(j)-(SweepRange(2)-SweepRange(1)) XRough(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
        for i = 1:length(PhaseVelocity)
            k1 = AngularFrequency/PhaseVelocity(i);
            k2 = k1^2;
            if  PhaseVelocity(i) > Material.TransverseVelocity
                x = sqrt(k2-kL2);
                y = sqrt(kT2-k2);
            elseif PhaseVelocity(i) < Material.TransverseVelocity
                x = sqrt(k2-kL2);
                y = sqrt(k2-kT2);
            end
            y2 = kT2-k2;
            xRo = x*Ro;
            yRo = y*Ro;
            xRi = x*Ri;
            yRi = y*Ri;
            if  PhaseVelocity(i) > Material.TransverseVelocity
                Z0xi = besseli(0,xRi);
                Z0xo = besseli(0,xRo);
                Z0yi = besselj(0,yRi);
                Z0yo = besselj(0,yRo);
                W0xi = besselk(0,xRi);
                W0xo = besselk(0,xRo);
                W0yi = bessely(0,yRi);
                W0yo = bessely(0,yRo);
                Z1xiRi = -xRi*besseli(1,xRi);
                Z1xoRo = -xRo*besseli(1,xRo);
                Z1yiRi = yRi*besselj(1,yRi);
                Z1yoRo = yRo*besselj(1,yRo);
                W1xiRi = xRi*besselk(1,xRi);
                W1xoRo = xRo*besselk(1,xRo);
                W1yiRi = yRi*bessely(1,yRi);
                W1yoRo = yRo*bessely(1,yRo);
            elseif PhaseVelocity(i) < Material.TransverseVelocity
                Z0xi = besseli(0,xRi);
                Z0xo = besseli(0,xRo);
                Z0yi = besseli(0,yRi);
                Z0yo = besseli(0,yRo);
                W0xi = besselk(0,xRi);
                W0xo = besselk(0,xRo);
                W0yi = besselk(0,yRi);
                W0yo = besselk(0,yRo);
                Z1xiRi = -xRi*besseli(1,xRi);
                Z1xoRo = -xRo*besseli(1,xRo);
                Z1yiRi = -yRi*besseli(1,yRi);
                Z1yoRo = -yRo*besseli(1,yRo);
                W1xiRi = xRi*besselk(1,xRi);
                W1xoRo = xRo*besselk(1,xRo);
                W1yiRi = yRi*besselk(1,yRi);
                W1yoRo = yRo*besselk(1,yRo);
            end
            zRi = sqrt(kInnerFluid2-k2)*Ri;
            Z0zi = besselj(0,zRi);
            Z1ziRi = -zRi*besselj(1,zRi);
            M(1,1) = Ri2*(.5*kT2-k2)*Z0xi-Z1xiRi;
            M(1,2) = Ri2*(.5*kT2-k2)*W0xi-W1xiRi;
            M(1,3) = k1*(Z1yiRi-y2*Ri2*Z0yi);
            M(1,4) = k1*(W1yiRi-y2*Ri2*W0yi);
            M(1,5) = .5*kT2*Ri2*Densityi*Z0zi;
            M(2,1) = -2*k1*Z1xiRi;
            M(2,2) = -2*k1*W1xiRi;
            M(2,3) = (k2-y2)*Z1yiRi;
            M(2,4) = (k2-y2)*W1yiRi;
            M(3,1) = -Z1xiRi;
            M(3,2) = -W1xiRi;
            M(3,3) = k1*Z1yiRi;
            M(3,4) = k1*W1yiRi;
            M(3,5) = Z1ziRi;
            M(4,1) = Ro2*(.5*kT2-k2)*Z0xo-Z1xoRo;
            M(4,2) = Ro2*(.5*kT2-k2)*W0xo-W1xoRo;
            M(4,3) = k1*(Z1yoRo-y2*Ro2*Z0yo);
            M(4,4) = k1*(W1yoRo-y2*Ro2*W0yo);
            M(5,1) = -2*k1*Z1xoRo;
            M(5,2) = -2*k1*W1xoRo;
            M(5,3) = (k2-y2)*Z1yoRo;
            M(5,4) = (k2-y2)*W1yoRo;
            Y(i) = abs(det(M));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    X0(j) = PhaseVelocity(i-1);
                end
                PhaseVelocity = [PhaseVelocity(i-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(i-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
                break
            end
        end
    end
end
Y = zeros(5^PhaseVelocitySections+1,3);
for p = 1:2
    X = X0(p);
    if  Multithreading
        send(Q1,[FrequencyRange(1),X0(p)/1e3,p])
    else
        addpoints(g(p),FrequencyRange(1),X0(p)/1e3);
        drawnow limitrate
    end
    Misses = 0;
    for i = 2:length(FrequencyRange)
        if  Stop
            return
        end
        AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
        kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
        kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
        X(i,1) = 0;
        if  p == 2
            Neighbors = [];
            if  i <= height(L{1})
                Neighbors = L{1}(i,4)*1e3;
            end
        end
        if  i == 2
            SweepRange = X(i-1)+5:-1:X(i-1)-5;
        elseif i == 3
            SweepRange = [X(i-1)+1*abs(X(i-2)-X(i-1)) X(i-1)-20*abs((X(i-2)-X(i-1)))];
        else
            SweepRange = [X(i-1)+10*abs(X(i-2)-X(i-1)) X(i-1)-10*abs((X(i-2)-X(i-1)))];  
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
            if  p == 2 && ~isempty(Neighbors)
                for j = 2:length(SweepRange)-1
                    if  SweepRange(j-1) > Neighbors && SweepRange(j+1) < Neighbors
                        SweepRange(j) = NaN;
                    end
                end
            end
            for j = 2:length(SweepRange)-1
                if  (SweepRange(j-1) > Material.TransverseVelocity && SweepRange(j+1) < Material.TransverseVelocity)% || (SweepRange(j-1) > InnerFluid.Velocity && SweepRange(j+1) < InnerFluid.Velocity)
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
                        if  isnan(PhaseVelocity(l))
                            Y(j,l) = NaN;
                        else
                            k1 = AngularFrequency/PhaseVelocity(l);
                            k2 = k1^2;
                            if  SweepRange(j) > Material.TransverseVelocity
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
                            if  SweepRange(j) > Material.TransverseVelocity
                                Z0xi = besseli(0,xRi);
                                Z0xo = besseli(0,xRo);
                                Z0yi = besselj(0,yRi);
                                Z0yo = besselj(0,yRo);
                                W0xi = besselk(0,xRi);
                                W0xo = besselk(0,xRo);
                                W0yi = bessely(0,yRi);
                                W0yo = bessely(0,yRo);
                                Z1xiRi = -xRi*besseli(1,xRi);
                                Z1xoRo = -xRo*besseli(1,xRo);
                                Z1yiRi = yRi*besselj(1,yRi);
                                Z1yoRo = yRo*besselj(1,yRo);
                                W1xiRi = xRi*besselk(1,xRi);
                                W1xoRo = xRo*besselk(1,xRo);
                                W1yiRi = yRi*bessely(1,yRi);
                                W1yoRo = yRo*bessely(1,yRo);
                            elseif SweepRange(j) < Material.TransverseVelocity
                                Z0xi = besseli(0,xRi);
                                Z0xo = besseli(0,xRo);
                                Z0yi = besseli(0,yRi);
                                Z0yo = besseli(0,yRo);
                                W0xi = besselk(0,xRi);
                                W0xo = besselk(0,xRo);
                                W0yi = besselk(0,yRi);
                                W0yo = besselk(0,yRo);
                                Z1xiRi = -xRi*besseli(1,xRi);
                                Z1xoRo = -xRo*besseli(1,xRo);
                                Z1yiRi = -yRi*besseli(1,yRi);
                                Z1yoRo = -yRo*besseli(1,yRo);
                                W1xiRi = xRi*besselk(1,xRi);
                                W1xoRo = xRo*besselk(1,xRo);
                                W1yiRi = yRi*besselk(1,yRi);
                                W1yoRo = yRo*besselk(1,yRo);
                            end
                            zRi = sqrt(kInnerFluid2-k2)*Ri;
                            Z0zi = besselj(0,zRi);
                            Z1ziRi = -zRi*besselj(1,zRi);
                            M(1,1) = Ri2*(.5*kT2-k2)*Z0xi-Z1xiRi;
                            M(1,2) = Ri2*(.5*kT2-k2)*W0xi-W1xiRi;
                            M(1,3) = k1*(Z1yiRi-y2*Ri2*Z0yi);
                            M(1,4) = k1*(W1yiRi-y2*Ri2*W0yi);
                            M(1,5) = .5*kT2*Ri2*Densityi*Z0zi;
                            M(2,1) = -2*k1*Z1xiRi;
                            M(2,2) = -2*k1*W1xiRi;
                            M(2,3) = (k2-y2)*Z1yiRi;
                            M(2,4) = (k2-y2)*W1yiRi;
                            M(3,1) = -Z1xiRi;
                            M(3,2) = -W1xiRi;
                            M(3,3) = k1*Z1yiRi;
                            M(3,4) = k1*W1yiRi;
                            M(3,5) = Z1ziRi;
                            M(4,1) = Ro2*(.5*kT2-k2)*Z0xo-Z1xoRo;
                            M(4,2) = Ro2*(.5*kT2-k2)*W0xo-W1xoRo;
                            M(4,3) = k1*(Z1yoRo-y2*Ro2*Z0yo);
                            M(4,4) = k1*(W1yoRo-y2*Ro2*W0yo);
                            M(5,1) = -2*k1*Z1xoRo;
                            M(5,2) = -2*k1*W1xoRo;
                            M(5,3) = (k2-y2)*Z1yoRo;
                            M(5,4) = (k2-y2)*W1yoRo;
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
            send(Q1,[FrequencyRange(i),X(i)/1e3,p])
        else
            addpoints(g(p),FrequencyRange(i),X(i)/1e3);
            drawnow limitrate
        end
    end
    X(Misses(1:length(X)) == 1) = NaN;
    L{p}(:,1) = FrequencyRange(1:length(X));
    L{p}(:,2) = FrequencyRange(1:length(X))/1e3;
    L{p}(:,3) = FrequencyRange(1:length(X))*(Ro-Ri);
    L{p}(:,4) = fillmissing(X,'spline')/1e3;
    L{p}(:,6) = 0;
% L{p}(:,5) = smooth(((L{p}(:,4)).^2)./(L{p}(:,4)-L{p}(:,1).*differentiate(fit(L{p}(:,1),L{p}(:,4),'cubicspline'),L{p}(:,1))));
end
if  HigherOrderModes && any(H)
    MissingModes(length(H)+2) = 0;
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
                    if  q == 1 && o == 0 && (any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) || (SweepRange(1) > InnerFluid.Velocity && SweepRange(2) < InnerFluid.Velocity) || (SweepRange(1) > Material.TransverseVelocity && SweepRange(2) < Material.TransverseVelocity)) % Neighbors cannot be excluded because SweepRange has only 2 elements
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
                        if  (SweepRange(j) > InnerFluid.Velocity && SweepRange(j+1) < InnerFluid.Velocity) || (SweepRange(j) > Material.TransverseVelocity && SweepRange(j+1) < Material.TransverseVelocity)
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
                                        Z0xi = besselj(0,xRi);
                                        Z0xo = besselj(0,xRo);
                                        Z0yi = besselj(0,yRi);
                                        Z0yo = besselj(0,yRo);
                                        W0xi = bessely(0,xRi);
                                        W0xo = bessely(0,xRo);
                                        W0yi = bessely(0,yRi);
                                        W0yo = bessely(0,yRo);
                                        Z1xiRi = xRi*besselj(1,xRi);
                                        Z1xoRo = xRo*besselj(1,xRo);
                                        Z1yiRi = yRi*besselj(1,yRi);
                                        Z1yoRo = yRo*besselj(1,yRo);
                                        W1xiRi = xRi*bessely(1,xRi);
                                        W1xoRo = xRo*bessely(1,xRo);
                                        W1yiRi = yRi*bessely(1,yRi);
                                        W1yoRo = yRo*bessely(1,yRo);
                                    elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                        Z0xi = besseli(0,xRi);
                                        Z0xo = besseli(0,xRo);
                                        Z0yi = besselj(0,yRi);
                                        Z0yo = besselj(0,yRo);
                                        W0xi = besselk(0,xRi);
                                        W0xo = besselk(0,xRo);
                                        W0yi = bessely(0,yRi);
                                        W0yo = bessely(0,yRo);
                                        Z1xiRi = -xRi*besseli(1,xRi);
                                        Z1xoRo = -xRo*besseli(1,xRo);
                                        Z1yiRi = yRi*besselj(1,yRi);
                                        Z1yoRo = yRo*besselj(1,yRo);
                                        W1xiRi = xRi*besselk(1,xRi);
                                        W1xoRo = xRo*besselk(1,xRo);
                                        W1yiRi = yRi*bessely(1,yRi);
                                        W1yoRo = yRo*bessely(1,yRo);
                                    elseif SweepRange(j) < Material.TransverseVelocity
                                        Z0xi = besseli(0,xRi);
                                        Z0xo = besseli(0,xRo);
                                        Z0yi = besseli(0,yRi);
                                        Z0yo = besseli(0,yRo);
                                        W0xi = besselk(0,xRi);
                                        W0xo = besselk(0,xRo);
                                        W0yi = besselk(0,yRi);
                                        W0yo = besselk(0,yRo);
                                        Z1xiRi = -xRi*besseli(1,xRi);
                                        Z1xoRo = -xRo*besseli(1,xRo);
                                        Z1yiRi = -yRi*besseli(1,yRi);
                                        Z1yoRo = -yRo*besseli(1,yRo);
                                        W1xiRi = xRi*besselk(1,xRi);
                                        W1xoRo = xRo*besselk(1,xRo);
                                        W1yiRi = yRi*besselk(1,yRi);
                                        W1yoRo = yRo*besselk(1,yRo);
                                    end
                                    zRi = sqrt(kInnerFluid2-k2)*Ri;
                                    Z0zi = besselj(0,zRi);
                                    Z1ziRi = -zRi*besselj(1,zRi);
                                    M(1,1) = Ri2*(.5*kT2-k2)*Z0xi-Z1xiRi;
                                    M(1,2) = Ri2*(.5*kT2-k2)*W0xi-W1xiRi;
                                    M(1,3) = k1*(Z1yiRi-y2*Ri2*Z0yi);
                                    M(1,4) = k1*(W1yiRi-y2*Ri2*W0yi);
                                    M(1,5) = .5*kT2*Ri2*Densityi*Z0zi;
                                    M(2,1) = -2*k1*Z1xiRi;
                                    M(2,2) = -2*k1*W1xiRi;
                                    M(2,3) = (k2-y2)*Z1yiRi;
                                    M(2,4) = (k2-y2)*W1yiRi;
                                    M(3,1) = -Z1xiRi;
                                    M(3,2) = -W1xiRi;
                                    M(3,3) = k1*Z1yiRi;
                                    M(3,4) = k1*W1yiRi;
                                    M(3,5) = Z1ziRi;
                                    M(4,1) = Ro2*(.5*kT2-k2)*Z0xo-Z1xoRo;
                                    M(4,2) = Ro2*(.5*kT2-k2)*W0xo-W1xoRo;
                                    M(4,3) = k1*(Z1yoRo-y2*Ro2*Z0yo);
                                    M(4,4) = k1*(W1yoRo-y2*Ro2*W0yo);
                                    M(5,1) = -2*k1*Z1xoRo;
                                    M(5,2) = -2*k1*W1xoRo;
                                    M(5,3) = (k2-y2)*Z1yoRo;
                                    M(5,4) = (k2-y2)*W1yoRo;
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
                MissingModes(p+2) = 1;
                break
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end))
                X(end-MissingSamples:end) = [];
                Misses(end-MissingSamples:end) = 0;
                break
            end
            if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                send(Q1,[FrequencyRange(i),X(i)/1e3,p+2])
            elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                addpoints(g(p+2),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
        end
        if  all(X == 0)
            MissingModes(p+2) = 1;
        end
        if  ~MissingModes(p+2)
            X(Misses(1:length(X)) == 1) = NaN;
            L{p+2}(:,1) = FrequencyRange(1:length(X));
            L{p+2}(:,2) = FrequencyRange(1:length(X))/1e3;
            L{p+2}(:,3) = FrequencyRange(1:length(X))*(Ro-Ri);
            L{p+2}(:,4) = fillmissing(X,'spline')/1e3;
        else
            L{p+2} = L{p+1};
            X1{p+2}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(L{p+2}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop
                return
            end
            X1{p+2}(i,1) = 0;
            if  i == 1
                if  p == 1
                    SweepRange = [L{p+2}(MaxInd,1)-2/4*FrequencyOffset/(Ro-Ri)/1e3 L{p+2}(MaxInd,1)+10/4/(Ro-Ri)/1e3];
                else
                    SweepRange = [L{p+2}(MaxInd,1)-1/4*FrequencyOffset/(Ro-Ri)/1e3 L{p+2}(MaxInd,1)+10/4/(Ro-Ri)/1e3];
                end
            elseif i == 2
                if  p == 1
                    SweepRange = [X1{p+2}(i-1)-2/4*FrequencyOffset/(Ro-Ri)/1e3 X1{p+2}(i-1)+10/4/(Ro-Ri)/1e3];
                else
                    SweepRange = [X1{p+2}(i-1)-1/4*FrequencyOffset/(Ro-Ri)/1e3 X1{p+2}(i-1)+10/4/(Ro-Ri)/1e3];
                end
            else
                SweepRange = X1{p+2}(i-1)-2*abs(X1{p+2}(i-2)-X1{p+2}(i-1));
                if  X1{p+2}(i-1) < X1{p+2}(i-2)
                    SweepRange(2) = X1{p+2}(i-1)+2/(Ro-Ri)/1e3;
                else
                    SweepRange(2) = X1{p+2}(i-1)+5*abs(X1{p+2}(i-2)-X1{p+2}(i-1));
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
                                Z0xi = besselj(0,xRi);
                                Z0xo = besselj(0,xRo);
                                Z0yi = besselj(0,yRi);
                                Z0yo = besselj(0,yRo);
                                W0xi = bessely(0,xRi);
                                W0xo = bessely(0,xRo);
                                W0yi = bessely(0,yRi);
                                W0yo = bessely(0,yRo);
                                Z1xiRi = xRi*besselj(1,xRi);
                                Z1xoRo = xRo*besselj(1,xRo);
                                Z1yiRi = yRi*besselj(1,yRi);
                                Z1yoRo = yRo*besselj(1,yRo);
                                W1xiRi = xRi*bessely(1,xRi);
                                W1xoRo = xRo*bessely(1,xRo);
                                W1yiRi = yRi*bessely(1,yRi);
                                W1yoRo = yRo*bessely(1,yRo);
                            elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity
                                Z0xi = besseli(0,xRi);
                                Z0xo = besseli(0,xRo);
                                Z0yi = besselj(0,yRi);
                                Z0yo = besselj(0,yRo);
                                W0xi = besselk(0,xRi);
                                W0xo = besselk(0,xRo);
                                W0yi = bessely(0,yRi);
                                W0yo = bessely(0,yRo);
                                Z1xiRi = -xRi*besseli(1,xRi);
                                Z1xoRo = -xRo*besseli(1,xRo);
                                Z1yiRi = yRi*besselj(1,yRi);
                                Z1yoRo = yRo*besselj(1,yRo);
                                W1xiRi = xRi*besselk(1,xRi);
                                W1xoRo = xRo*besselk(1,xRo);
                                W1yiRi = yRi*bessely(1,yRi);
                                W1yoRo = yRo*bessely(1,yRo);
                            end
                            zRi = sqrt((AngularFrequency/InnerFluid.Velocity)^2-k2)*Ri;
                            Z0zi = besselj(0,zRi);
                            Z1ziRi = -zRi*besselj(1,zRi);
                            M(1,1) = Ri2*(.5*kT2-k2)*Z0xi-Z1xiRi;
                            M(1,2) = Ri2*(.5*kT2-k2)*W0xi-W1xiRi;
                            M(1,3) = k1*(Z1yiRi-y2*Ri2*Z0yi);
                            M(1,4) = k1*(W1yiRi-y2*Ri2*W0yi);
                            M(1,5) = .5*kT2*Ri2*Densityi*Z0zi;
                            M(2,1) = -2*k1*Z1xiRi;
                            M(2,2) = -2*k1*W1xiRi;
                            M(2,3) = (k2-y2)*Z1yiRi;
                            M(2,4) = (k2-y2)*W1yiRi;
                            M(3,1) = -Z1xiRi;
                            M(3,2) = -W1xiRi;
                            M(3,3) = k1*Z1yiRi;
                            M(3,4) = k1*W1yiRi;
                            M(3,5) = Z1ziRi;
                            M(4,1) = Ro2*(.5*kT2-k2)*Z0xo-Z1xoRo;
                            M(4,2) = Ro2*(.5*kT2-k2)*W0xo-W1xoRo;
                            M(4,3) = k1*(Z1yoRo-y2*Ro2*Z0yo);
                            M(4,4) = k1*(W1yoRo-y2*Ro2*W0yo);
                            M(5,1) = -2*k1*Z1xoRo;
                            M(5,2) = -2*k1*W1xoRo;
                            M(5,3) = (k2-y2)*Z1yoRo;
                            M(5,4) = (k2-y2)*W1yoRo;
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
%                             z = isoutlier(vertcat(X1{p+2}(1:i-1),Frequency(2)),'movmedian',5,'ThresholdFactor',9);
%                             if  z(end) && abs(X1{p+2}(i-1)-Frequency(2)) > 1
%                                 Outlier = 1;
%                             else
%                                 Outlier = 0;
%                             end
%                         end
%                         if  ~Outlier
                            X1{p+2}(i,1) = Frequency(2);
                            break
%                         end
                    end
                end
                if  X1{p+2}(i) > 0
                    break
                end
            end
            if  X1{p+2}(i) == 0
                break
            end
            if  Multithreading
                send(Q2,[X1{p+2}(i),PhaseVelocityRange(i)/1e3,p+2])
            else
                addpoints(g1(p+2),X1{p+2}(i),PhaseVelocityRange(i)/1e3);
                drawnow limitrate
            end
        end
        X1{p+2}(:,2) = X1{p+2}(:,1)/1e3;
        X1{p+2}(:,3) = X1{p+2}(:,1)*(Ro-Ri);
        X1{p+2}(:,4) = PhaseVelocityRange(1:height(X1{p+2}))/1e3;
        if  X1{p+2}(end,1) == 0
            X1{p+2}(end,1) = H(p);
            X1{p+2}(end,2) = H(p)/1e3;
            X1{p+2}(end,3) = H(p)*(Ro-Ri);
            X1{p+2}(end,4) = PhaseVelocityLimit/1e3;
        end
        X1{p+2}(X1{p+2}(:,1) == 0,:) = [];
        X1{p+2} = flipud(X1{p+2});
    end
    X1(MissingModes == 1) = [];
    L(MissingModes == 1) = [];
    for p = 3:length(L)
        L{p}(L{p}(:,4) == 0,:) = [];
        L{p} = vertcat(X1{p},L{p});
        L{p}(:,6) = 0;
% L{p}(:,5) = smooth(((L{p}(:,4)).^2)./(L{p}(:,4)-L{p}(:,1).*differentiate(fit(L{p}(:,1),L{p}(:,4),'cubicspline'),L{p}(:,1))));
    end
end