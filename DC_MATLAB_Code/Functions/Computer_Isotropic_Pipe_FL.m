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
function A = Computer_Isotropic_Pipe_FL(Multithreading,Q,ax,ToggleInnerFluid,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,H,LineColor,n)        
MissingSamples = MissingSamples+2;

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};B={};
n2 = n^2;
Ro2 = Ro^2;
Ri2 = Ri^2;
cL2 = Material.LongitudinalVelocity^2;
cT2 = Material.TransverseVelocity^2;
cFi2 = InnerFluid.Velocity^2;
Densityi = InnerFluid.Density/Material.Density;
if  n == 0
    if  ToggleInnerFluid
        X0 = [];
        XRough = [];
        AngularFrequency = pi*FrequencyRange(1)*2e3;
        AngularFrequency2 = AngularFrequency^2;
        kL2 = AngularFrequency2/cL2;
        kT2 = AngularFrequency2/cT2;
        kFi2 = AngularFrequency2/cFi2;
        SweepRange = 50:Material.CylinderVelocity+500;
        Y = Computer(AngularFrequency./SweepRange,Ro,Ri,Ro2,Ri2,kL2,kT2,kFi2,Densityi,length(SweepRange));
        for i = 2:length(SweepRange)-1
            if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
                XRough(end+1) = SweepRange(i);
            end
        end
        SweepRangeStep = abs(SweepRange(1)-SweepRange(2));
        XRough(XRough > Material.TransverseVelocity-SweepRangeStep/2 & XRough < Material.TransverseVelocity+SweepRangeStep/2) = [];
        Bisections = ceil(log2(SweepRangeStep/1e-6));
        for j = 1:length(XRough)
            PhaseVelocity = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
            for o = 1:Bisections
                PhaseVelocity = PhaseVelocity(1):(PhaseVelocity(end)-PhaseVelocity(1))/4:PhaseVelocity(end);
                [Y,Yc] = Computer(AngularFrequency./PhaseVelocity,Ro,Ri,Ro2,Ri2,kL2,kT2,kFi2,Densityi,length(PhaseVelocity));
                for i = 2:length(PhaseVelocity)-1
                    if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
                        if  o == Bisections && sign(Yc(i-1)) ~= sign(Yc(i+1))
                            X0(end+1) = PhaseVelocity(i);
                        end
                        PhaseVelocity = [PhaseVelocity(i-1) PhaseVelocity(i+1)];
                        break
                    end
                end
            end
        end
    else
        X0 = Material.CylinderVelocity;
    end
end
if  ~Multithreading
    g = [animatedline(ax,'color',LineColor) animatedline(ax,'color',LineColor)];
end
Y = zeros(5^(PhaseVelocitySections+1)+1,3);
if  n == 0
    Neighbors0 = Material.TransverseVelocity;
    H = [zeros(1,length(X0)) H];
else
    Neighbors0 = [];
    if  n == 1
        H = [0 H];
    end
end
for p = 1:length(H)
    X = 0;
    Misses = false;
    BelowCutoff = false;
    i1 = ceil(H(p)/FrequencyResolution)+1;
    if  p == 1
        Step = FrequencyResolution/5;
        if  i1 == 1
            x = FrequencyRange(1);
        else
            x = FrequencyRange(i1-1):Step:FrequencyRange(i1);
        end
        f = x(find(x > H(p),1));
        if  isempty(f)
            f = FrequencyRange(i1);
        end
        Range = f+40*Step;
        FrequencyRangeEdit = [FrequencyRange(1:i1-1) f:Step:Range FrequencyRange(FrequencyRange > Range)];
    else
        FrequencyRangeEdit = FrequencyRange;
    end
    for i = i1:length(FrequencyRangeEdit)
        if  Stop
            return
        end
        AngularFrequency = pi*FrequencyRangeEdit(i)*2e3;
        AngularFrequency2 = AngularFrequency^2;
        kL2 = AngularFrequency2/cL2;
        kT2 = AngularFrequency2/cT2;
        kFi2 = AngularFrequency2/cFi2;
        if  n == 0
            a1 = kT2/2;
        else
            a1 = kT2/2*Ri2;
            a2 = kT2/2*Ro2;
        end
        X(i,1) = 0;
        Misses(i) = false;
        BelowCutoff(i) = false;
        Neighbors = Neighbors0;
        for j = 1:length(A)
            if  FrequencyRangeEdit(i) <= A{j}(end,1)
                [~,z] = min(abs(FrequencyRangeEdit(i)-A{j}(:,1)));
                Neighbors(end+1) = A{j}(z,4);
            end
        end
        if  ~any(Neighbors)
            Neighbors = 0;
        end
        for q = 1:2
            if  q == 1
                if  n == 0 && p <= length(X0)
                    if  i == 1
                        X = X0(p);
                        if  Multithreading
                            send(Q(1),[FrequencyRange(1) X/1e3])
                        else
                            addpoints(g(1),FrequencyRange(1),X/1e3);
                            drawnow limitrate
                        end
                        continue
                    elseif i == 2
                        SweepRange = X(i-1)+5:-1:X(i-1)-5;
                    elseif i == 3
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[1 -20];
                    else
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[10 -10];
                    end
                    if  p == length(X0) && SweepRange(1) > X0(p)
                        SweepRange(1) = X0(p);
                    end
                    if  SweepRange(2) < 0
                        SweepRange(2) = 0;
                    end
                elseif n == 1 && p == 1
                    if  i == 1
                        SweepRange = 1e-2:2e2;
                    elseif i == 2
                        SweepRange = X(i-1):4e3;
                    else
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[-10 10];
                        if  SweepRange(2) > 1.01*Material.RayleighVelocity
                            SweepRange(2) = 1.01*Material.RayleighVelocity;
                        end
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
                        if  n > 1 && p == 1 && X(i-1) > min(X(X > 0))
                            SweepRange = X(i-1)+dy*[10 -10];
                            if  SweepRange(1) > Material.TransverseVelocity
                                SweepRange(1) = Material.TransverseVelocity;
                            end
                        else
                            SweepRange = X(i-1)+dy*[.1 -Factor];
                        end
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                        SweepRange = SweepRange+Accuracy*[1 -1];
                    end
                    if  SweepRange(2) < max(Neighbors(Neighbors < X(i-1)))+1
                        SweepRange(2) = max(Neighbors(Neighbors < X(i-1)))+1;
                    end
                end
            else
                if  (n == 0 && p <= length(X0)) || (n == 1 && p == 1)
                    continue
                end
                if  ~any(X)
                    SweepRange = [PhaseVelocityLimit min(Neighbors)+1];
                elseif isscalar(find(X))
                    SweepRange = [X(i-1) min(Neighbors)+1];
                else
                    if  n < 2 && ToggleInnerFluid
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[.1 10];
                    else
                        SweepRange = X(i-1)+.1*abs(X(i-2)-X(i-1))+[0 50];
                    end
                end
            end
            for o = 0:PhaseVelocitySections+1
                if  ((n == 0 && p <= length(X0)) || (n == 1 && p == 1)) && i > 2
                    if  o == 0
                        continue
                    else
                        SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
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
                            SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
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
                    if  q == 1 && o == 0 && (any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) || (n > 1 && ToggleInnerFluid && SweepRange(1) > InnerFluid.Velocity && SweepRange(2) < InnerFluid.Velocity)) % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                end
                dy = abs(SweepRange(1)-SweepRange(2));
                if  Accuracy < dy
                    Bisections = ceil(log2(dy/Accuracy));
                else
                    Bisections = 1;
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
                    if  n > 1 && ToggleInnerFluid && SweepRange(j) > InnerFluid.Velocity && SweepRange(j+1) < InnerFluid.Velocity
                        if  j == 1
                            SweepRange(2) = NaN;
                        else
                            SweepRange(j) = NaN;
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
                                k1 = AngularFrequency/PhaseVelocity(l);
                                k2 = k1^2;
                                y2 = kT2-k2;
                                if  SweepRange(j) > Material.LongitudinalVelocity
                                    x = sqrt(kL2-k2);
                                    y = sqrt(y2);
                                elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                    x = sqrt(k2-kL2);
                                    y = sqrt(y2);
                                elseif SweepRange(j) < Material.TransverseVelocity
                                    x = sqrt(k2-kL2);
                                    y = sqrt(-y2);
                                end
                                xR = x*[Ri Ro];
                                yR = y*[Ri Ro];
                                if  n == 0
                                    a2 = a1-k2;
                                    a3 = k2-y2;
                                    if  SweepRange(j) > Material.LongitudinalVelocity
                                        Z0xR2 = [Ri2 Ro2].*besselj(0,xR);
                                        Z0yR2 = [Ri2 Ro2].*besselj(0,yR);
                                        W0xR2 = [Ri2 Ro2].*bessely(0,xR);
                                        W0yR2 = [Ri2 Ro2].*bessely(0,yR);
                                        Z1xR = xR.*besselj(1,xR);
                                        Z1yR = yR.*besselj(1,yR);
                                        W1xR = xR.*bessely(1,xR);
                                        W1yR = yR.*bessely(1,yR);
                                    elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                        Z0xR2 = [Ri2 Ro2].*besseli(0,xR);
                                        Z0yR2 = [Ri2 Ro2].*besselj(0,yR);
                                        W0xR2 = [Ri2 Ro2].*besselk(0,xR);
                                        W0yR2 = [Ri2 Ro2].*bessely(0,yR);
                                        Z1xR = -xR.*besseli(1,xR);
                                        Z1yR = yR.*besselj(1,yR);
                                        W1xR = xR.*besselk(1,xR);
                                        W1yR = yR.*bessely(1,yR);
                                    elseif SweepRange(j) < Material.TransverseVelocity
                                        Z0xR2 = [Ri2 Ro2].*besseli(0,xR);
                                        Z0yR2 = [Ri2 Ro2].*besseli(0,yR);
                                        W0xR2 = [Ri2 Ro2].*besselk(0,xR);
                                        W0yR2 = [Ri2 Ro2].*besselk(0,yR);
                                        Z1xR = -xR.*besseli(1,xR);
                                        Z1yR = -yR.*besseli(1,yR);
                                        W1xR = xR.*besselk(1,xR);
                                        W1yR = yR.*besselk(1,yR);
                                    end
                                    M(1,1) = a2*Z0xR2(1)-Z1xR(1);
                                    M(1,2) = a2*W0xR2(1)-W1xR(1);
                                    M(1,3) = k1*(Z1yR(1)-y2*Z0yR2(1));
                                    M(1,4) = k1*(W1yR(1)-y2*W0yR2(1));
                                    M(2,1) = -2*k1*Z1xR(1);
                                    M(2,2) = -2*k1*W1xR(1);
                                    M(2,3) = a3*Z1yR(1);
                                    M(2,4) = a3*W1yR(1);
                                    if  ToggleInnerFluid
                                        zRi = sqrt(kFi2-k2)*Ri;
                                        M(1,5) = a1*Ri2*Densityi*besselj(0,zRi);
                                        M(3,1) = -Z1xR(1);
                                        M(3,2) = -W1xR(1);
                                        M(3,3) = k1*Z1yR(1);
                                        M(3,4) = k1*W1yR(1);
                                        M(3,5) = -zRi*besselj(1,zRi);
                                        M(4,1) = a2*Z0xR2(2)-Z1xR(2);
                                        M(4,2) = a2*W0xR2(2)-W1xR(2);
                                        M(4,3) = k1*(Z1yR(2)-y2*Z0yR2(2));
                                        M(4,4) = k1*(W1yR(2)-y2*W0yR2(2));
                                        M(5,1) = -2*k1*Z1xR(2);
                                        M(5,2) = -2*k1*W1xR(2);
                                        M(5,3) = a3*Z1yR(2);
                                        M(5,4) = a3*W1yR(2);
                                    else
                                        M(3,1) = a2*Z0xR2(2)-Z1xR(2);
                                        M(3,2) = a2*W0xR2(2)-W1xR(2);
                                        M(3,3) = k1*(Z1yR(2)-y2*Z0yR2(2));
                                        M(3,4) = k1*(W1yR(2)-y2*W0yR2(2));
                                        M(4,1) = -2*k1*Z1xR(2);
                                        M(4,2) = -2*k1*W1xR(2);
                                        M(4,3) = a3*Z1yR(2);
                                        M(4,4) = a3*W1yR(2);
                                    end
                                else
                                    a3 = y2-k2;
                                    if  SweepRange(j) > Material.LongitudinalVelocity
                                        Znx = besselj(n,xR);
                                        Zny = besselj(n,yR);
                                        Wnx = bessely(n,xR);
                                        Wny = bessely(n,yR);
                                        dZnxR = n*Znx-xR.*besselj(n+1,xR);
                                        dZnyR = n*Zny-yR.*besselj(n+1,yR);
                                        dWnxR = n*Wnx-xR.*bessely(n+1,xR);
                                        dWnyR = n*Wny-yR.*bessely(n+1,yR);
                                    elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                        Znx = besseli(n,xR);
                                        Zny = besselj(n,yR);
                                        Wnx = besselk(n,xR);
                                        Wny = bessely(n,yR);
                                        dZnxR = n*Znx+xR.*besseli(n+1,xR);
                                        dZnyR = n*Zny-yR.*besselj(n+1,yR);
                                        dWnxR = n*Wnx-xR.*besselk(n+1,xR);
                                        dWnyR = n*Wny-yR.*bessely(n+1,yR);
                                    elseif SweepRange(j) < Material.TransverseVelocity
                                        Znx = besseli(n,xR);
                                        Zny = besseli(n,yR);
                                        Wnx = besselk(n,xR);
                                        Wny = besselk(n,yR);
                                        dZnxR = n*Znx+xR.*besseli(n+1,xR);
                                        dZnyR = n*Zny+yR.*besseli(n+1,yR);
                                        dWnxR = n*Wnx-xR.*besselk(n+1,xR);
                                        dWnyR = n*Wny-yR.*besselk(n+1,yR);
                                    end
                                    M(1,1) = dZnxR(1)+(a1-k2*Ri2-n2)*Znx(1);
                                    M(1,2) = dWnxR(1)+(a1-k2*Ri2-n2)*Wnx(1);
                                    M(1,3) = -k1*(dZnyR(1)+(y2*Ri2-n2)*Zny(1));
                                    M(1,4) = -k1*(dWnyR(1)+(y2*Ri2-n2)*Wny(1));
                                    M(1,5) = n*(dZnyR(1)-Zny(1));
                                    M(1,6) = n*(dWnyR(1)-Wny(1));
                                    M(2,1) = 2*n*(dZnxR(1)-Znx(1));
                                    M(2,2) = 2*n*(dWnxR(1)-Wnx(1));
                                    M(2,3) = 2*k1*n*(Zny(1)-dZnyR(1));
                                    M(2,4) = 2*k1*n*(Wny(1)-dWnyR(1));
                                    M(2,5) = 2*dZnyR(1)+(y2*Ri2-2*n2)*Zny(1);
                                    M(2,6) = 2*dWnyR(1)+(y2*Ri2-2*n2)*Wny(1);
                                    M(3,1) = 2*k1*dZnxR(1);
                                    M(3,2) = 2*k1*dWnxR(1);
                                    M(3,3) = a3*dZnyR(1);
                                    M(3,4) = a3*dWnyR(1);
                                    M(3,5) = -k1*n*Zny(1);
                                    M(3,6) = -k1*n*Wny(1);
                                    if  ToggleInnerFluid
                                        zRi = sqrt(kFi2-k2)*Ri;
                                        Znzi = besselj(n,zRi);
                                        M(1,7) = a1*Densityi*Znzi;
                                        M(4,1) = dZnxR(1);
                                        M(4,2) = dWnxR(1);
                                        M(4,3) = -k1*dZnyR(1);
                                        M(4,4) = -k1*dWnyR(1);
                                        M(4,5) = -n*Zny(1);
                                        M(4,6) = -n*Wny(1);
                                        M(4,7) = n*Znzi-zRi*besselj(n+1,zRi);
                                        M(5,1) = dZnxR(2)+(a2-k2*Ro2-n2)*Znx(2);
                                        M(5,2) = dWnxR(2)+(a2-k2*Ro2-n2)*Wnx(2);
                                        M(5,3) = -k1*(dZnyR(2)+(y2*Ro2-n2)*Zny(2));
                                        M(5,4) = -k1*(dWnyR(2)+(y2*Ro2-n2)*Wny(2));
                                        M(5,5) = n*(dZnyR(2)-Zny(2));
                                        M(5,6) = n*(dWnyR(2)-Wny(2));
                                        M(6,1) = 2*n*(dZnxR(2)-Znx(2));
                                        M(6,2) = 2*n*(dWnxR(2)-Wnx(2));
                                        M(6,3) = 2*k1*n*(Zny(2)-dZnyR(2));
                                        M(6,4) = 2*k1*n*(Wny(2)-dWnyR(2));
                                        M(6,5) = 2*dZnyR(2)+(y2*Ro2-2*n2)*Zny(2);
                                        M(6,6) = 2*dWnyR(2)+(y2*Ro2-2*n2)*Wny(2);
                                        M(7,1) = 2*k1*dZnxR(2);
                                        M(7,2) = 2*k1*dWnxR(2);
                                        M(7,3) = a3*dZnyR(2);
                                        M(7,4) = a3*dWnyR(2);
                                        M(7,5) = -k1*n*Zny(2);
                                        M(7,6) = -k1*n*Wny(2);
                                    else
                                        M(4,1) = dZnxR(2)+(a2-k2*Ro2-n2)*Znx(2);
                                        M(4,2) = dWnxR(2)+(a2-k2*Ro2-n2)*Wnx(2);
                                        M(4,3) = -k1*(dZnyR(2)+(y2*Ro2-n2)*Zny(2));
                                        M(4,4) = -k1*(dWnyR(2)+(y2*Ro2-n2)*Wny(2));
                                        M(4,5) = n*(dZnyR(2)-Zny(2));
                                        M(4,6) = n*(dWnyR(2)-Wny(2));
                                        M(5,1) = 2*n*(dZnxR(2)-Znx(2));
                                        M(5,2) = 2*n*(dWnxR(2)-Wnx(2));
                                        M(5,3) = 2*k1*n*(Zny(2)-dZnyR(2));
                                        M(5,4) = 2*k1*n*(Wny(2)-dWnyR(2));
                                        M(5,5) = 2*dZnyR(2)+(y2*Ro2-2*n2)*Zny(2);
                                        M(5,6) = 2*dWnyR(2)+(y2*Ro2-2*n2)*Wny(2);
                                        M(6,1) = 2*k1*dZnxR(2);
                                        M(6,2) = 2*k1*dWnxR(2);
                                        M(6,3) = a3*dZnyR(2);
                                        M(6,4) = a3*dWnyR(2);
                                        M(6,5) = -k1*n*Zny(2);
                                        M(6,6) = -k1*n*Wny(2);
                                    end
                                end
                                if  ToggleInnerFluid && mod(n,2)
                                    Y(j,l) = real(det(M)/Znzi);
                                else
                                    Y(j,l) = real(det(M));
                                end
                            end
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  sign(Y(j,1)) ~= sign(Y(j,2)) && ((j == 1 && abs(Y(j,2)) < abs(Y(j,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))
                            PhaseVelocity = [PhaseVelocity(1) (PhaseVelocity(1)+PhaseVelocity(2))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif sign(Y(j,2)) ~= sign(Y(j,3)) && ((j == 1 && abs(Y(j,2)) < abs(Y(j,1))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))
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
                if  ~((n == 0 && p <= length(X0)) || (n == 1 && p == 1))
                    X(i-1) = 0;
                    BelowCutoff(i-1:i) = true;
                end
            else
                Idx = i1:i-1;
                Idx(BelowCutoff(Idx)) = [];
                if  length(Idx) > 50
                    Idx(1:end-50) = [];
                end
                Fit = fit(FrequencyRangeEdit(Idx)',X(Idx),'cubicspline');
                X(i) = Fit(FrequencyRangeEdit(i));
                Misses(i) = true; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
% figure,line(FrequencyRangeEdit(Idx),X(Idx),'LineWidth',4,'color','r'),line(FrequencyRangeEdit([Idx i]),Fit(FrequencyRangeEdit([Idx i])),'LineWidth',1.5,'color','g')
            end
        elseif ~any(X)
            if  n == 1 && p == 1
                X(i,1) = 1e-2;
            else
                BelowCutoff(i) = true;
            end
        end
        if  FrequencyResolution*length(BelowCutoff(BelowCutoff)) > BelowCutoffWidth/(Ro-Ri)/1e3
            break
        elseif i > MissingSamples && all(Misses(i-MissingSamples+1:i))
            X(i-MissingSamples+1:i) = [];
            Misses(i-MissingSamples+1:i) = [];
            break
        end
        if  ~BelowCutoff(i) && ~Misses(i)
            if  Multithreading
                send(Q(1),[FrequencyRangeEdit(i) X(i)/1e3])
            else
                addpoints(g(1),FrequencyRangeEdit(i),X(i)/1e3);
                drawnow limitrate
            end
        end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRangeEdit(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k)];
% if  Misses(i)
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i)
%     disp(['p = ',num2str(p),' f = ',num2str(FrequencyRangeEdit(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k),' Miss'])
% end
    end
    if  numel(find(X)) < 2
        continue
    end
    X(Misses) = NaN;
    A{end+1} = [FrequencyRangeEdit(1:length(X))'*[1 1e-3 Ro-Ri] fillmissing(X,'spline')];
    z = A{end}(:,4) > 0;
    if  Multithreading
        send(Q(1),[A{end}(z,[1 4])])
    else
        line(ax,A{end}(z),A{end}(z,4)/1e3,'color',LineColor)
        drawnow limitrate
        clearpoints(g(1))
    end
    if  ~HigherOrderModes || ~any(H)
        break
    elseif (n == 0 && p <= length(X0)) || (n == 1 && p == 1)
        continue
    end
    B{end+1} = 0;
    Idx0 = find(A{end}(:,4),1);
    PhaseVelocityRange = PhaseVelocityStep*ceil(A{end}(Idx0,4)/PhaseVelocityStep):PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
    PhaseVelocityRange(PhaseVelocityRange == Material.TransverseVelocity | PhaseVelocityRange == Material.LongitudinalVelocity) = [];
    if  ToggleInnerFluid
        PhaseVelocityRange(PhaseVelocityRange == InnerFluid.Velocity) = [];
    end
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
        B{end}(i,1) = 0;
        if  ToggleInnerFluid
            Frac = 4;
        else
            Frac = 1;
        end
        if  i == 1
            X1 = A{end}(Idx0);
        else
            X1 = B{end}(i-1);
        end
        if  (n == 0 && p == length(X0)+1) || (n == 1 && p == 2) || (n > 1 && p == 1)
            SweepRange = X1+[-2*FrequencyOffset 10]/Frac/(Ro-Ri)/1e3;
        else
            SweepRange = X1+[-FrequencyOffset 10]/Frac/(Ro-Ri)/1e3;
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
                SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
            end
            dy = abs(SweepRange(1)-SweepRange(2));
            if  Accuracy < dy
                Bisections = ceil(log2(dy/Accuracy));
            else
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
                        AngularFrequency = pi*Frequency(l)*2e3;
                        AngularFrequency2 = AngularFrequency^2;
                        k1 = AngularFrequency/PhaseVelocityRange(i);
                        k2 = k1^2;
                        kT2 = AngularFrequency2/cT2;
                        y2 = kT2-k2;
                        if  PhaseVelocityRange(i) > Material.LongitudinalVelocity
                            x = sqrt(AngularFrequency2/cL2-k2);
                            y = sqrt(y2);
                        elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity && PhaseVelocityRange(i) > Material.TransverseVelocity
                            x = sqrt(k2-AngularFrequency2/cL2);
                            y = sqrt(y2);
                        elseif PhaseVelocityRange(i) < Material.TransverseVelocity
                            x = sqrt(k2-AngularFrequency2/cL2);
                            y = sqrt(-y2);
                        end
                        xR = x*[Ri Ro];
                        yR = y*[Ri Ro];
                        if  n == 0
                            a1 = kT2/2;
                            a2 = a1-k2;
                            a3 = k2-y2;
                            if  PhaseVelocityRange(i) > Material.LongitudinalVelocity
                                Z0xR2 = [Ri2 Ro2].*besselj(0,xR);
                                Z0yR2 = [Ri2 Ro2].*besselj(0,yR);
                                W0xR2 = [Ri2 Ro2].*bessely(0,xR);
                                W0yR2 = [Ri2 Ro2].*bessely(0,yR);
                                Z1xR = xR.*besselj(1,xR);
                                Z1yR = yR.*besselj(1,yR);
                                W1xR = xR.*bessely(1,xR);
                                W1yR = yR.*bessely(1,yR);
                            elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity && PhaseVelocityRange(i) > Material.TransverseVelocity
                                Z0xR2 = [Ri2 Ro2].*besseli(0,xR);
                                Z0yR2 = [Ri2 Ro2].*besselj(0,yR);
                                W0xR2 = [Ri2 Ro2].*besselk(0,xR);
                                W0yR2 = [Ri2 Ro2].*bessely(0,yR);
                                Z1xR = -xR.*besseli(1,xR);
                                Z1yR = yR.*besselj(1,yR);
                                W1xR = xR.*besselk(1,xR);
                                W1yR = yR.*bessely(1,yR);
                            elseif PhaseVelocityRange(i) < Material.TransverseVelocity
                                Z0xR2 = [Ri2 Ro2].*besseli(0,xR);
                                Z0yR2 = [Ri2 Ro2].*besseli(0,yR);
                                W0xR2 = [Ri2 Ro2].*besselk(0,xR);
                                W0yR2 = [Ri2 Ro2].*besselk(0,yR);
                                Z1xR = -xR.*besseli(1,xR);
                                Z1yR = -yR.*besseli(1,yR);
                                W1xR = xR.*besselk(1,xR);
                                W1yR = yR.*besselk(1,yR);
                            end
                            M(1,1) = a2*Z0xR2(1)-Z1xR(1);
                            M(1,2) = a2*W0xR2(1)-W1xR(1);
                            M(1,3) = k1*(Z1yR(1)-y2*Z0yR2(1));
                            M(1,4) = k1*(W1yR(1)-y2*W0yR2(1));
                            M(2,1) = -2*k1*Z1xR(1);
                            M(2,2) = -2*k1*W1xR(1);
                            M(2,3) = a3*Z1yR(1);
                            M(2,4) = a3*W1yR(1);
                            if  ToggleInnerFluid
                                zRi = sqrt(AngularFrequency2/cFi2-k2)*Ri;
                                M(1,5) = a1*Ri2*Densityi*besselj(0,zRi);
                                M(3,1) = -Z1xR(1);
                                M(3,2) = -W1xR(1);
                                M(3,3) = k1*Z1yR(1);
                                M(3,4) = k1*W1yR(1);
                                M(3,5) = -zRi*besselj(1,zRi);
                                M(4,1) = a2*Z0xR2(2)-Z1xR(2);
                                M(4,2) = a2*W0xR2(2)-W1xR(2);
                                M(4,3) = k1*(Z1yR(2)-y2*Z0yR2(2));
                                M(4,4) = k1*(W1yR(2)-y2*W0yR2(2));
                                M(5,1) = -2*k1*Z1xR(2);
                                M(5,2) = -2*k1*W1xR(2);
                                M(5,3) = a3*Z1yR(2);
                                M(5,4) = a3*W1yR(2);
                            else
                                M(3,1) = a2*Z0xR2(2)-Z1xR(2);
                                M(3,2) = a2*W0xR2(2)-W1xR(2);
                                M(3,3) = k1*(Z1yR(2)-y2*Z0yR2(2));
                                M(3,4) = k1*(W1yR(2)-y2*W0yR2(2));
                                M(4,1) = -2*k1*Z1xR(2);
                                M(4,2) = -2*k1*W1xR(2);
                                M(4,3) = a3*Z1yR(2);
                                M(4,4) = a3*W1yR(2);
                            end
                        else
                            a1 = kT2/2*Ri2;
                            a2 = kT2/2*Ro2;
                            a3 = y2-k2;
                            if  PhaseVelocityRange(i) > Material.LongitudinalVelocity
                                Znx = besselj(n,xR);
                                Zny = besselj(n,yR);
                                Wnx = bessely(n,xR);
                                Wny = bessely(n,yR);
                                dZnxR = n*Znx-xR.*besselj(n+1,xR);
                                dZnyR = n*Zny-yR.*besselj(n+1,yR);
                                dWnxR = n*Wnx-xR.*bessely(n+1,xR);
                                dWnyR = n*Wny-yR.*bessely(n+1,yR);
                            elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity && PhaseVelocityRange(i) > Material.TransverseVelocity
                                Znx = besseli(n,xR);
                                Zny = besselj(n,yR);
                                Wnx = besselk(n,xR);
                                Wny = bessely(n,yR);
                                dZnxR = n*Znx+xR.*besseli(n+1,xR);
                                dZnyR = n*Zny-yR.*besselj(n+1,yR);
                                dWnxR = n*Wnx-xR.*besselk(n+1,xR);
                                dWnyR = n*Wny-yR.*bessely(n+1,yR);
                            elseif PhaseVelocityRange(i) < Material.TransverseVelocity
                                Znx = besseli(n,xR);
                                Zny = besseli(n,yR);
                                Wnx = besselk(n,xR);
                                Wny = besselk(n,yR);
                                dZnxR = n*Znx+xR.*besseli(n+1,xR);
                                dZnyR = n*Zny+yR.*besseli(n+1,yR);
                                dWnxR = n*Wnx-xR.*besselk(n+1,xR);
                                dWnyR = n*Wny-yR.*besselk(n+1,yR);
                            end
                            M(1,1) = dZnxR(1)+(a1-k2*Ri2-n2)*Znx(1);
                            M(1,2) = dWnxR(1)+(a1-k2*Ri2-n2)*Wnx(1);
                            M(1,3) = -k1*(dZnyR(1)+(y2*Ri2-n2)*Zny(1));
                            M(1,4) = -k1*(dWnyR(1)+(y2*Ri2-n2)*Wny(1));
                            M(1,5) = n*(dZnyR(1)-Zny(1));
                            M(1,6) = n*(dWnyR(1)-Wny(1));
                            M(2,1) = 2*n*(dZnxR(1)-Znx(1));
                            M(2,2) = 2*n*(dWnxR(1)-Wnx(1));
                            M(2,3) = 2*k1*n*(Zny(1)-dZnyR(1));
                            M(2,4) = 2*k1*n*(Wny(1)-dWnyR(1));
                            M(2,5) = 2*dZnyR(1)+(y2*Ri2-2*n2)*Zny(1);
                            M(2,6) = 2*dWnyR(1)+(y2*Ri2-2*n2)*Wny(1);
                            M(3,1) = 2*k1*dZnxR(1);
                            M(3,2) = 2*k1*dWnxR(1);
                            M(3,3) = a3*dZnyR(1);
                            M(3,4) = a3*dWnyR(1);
                            M(3,5) = -k1*n*Zny(1);
                            M(3,6) = -k1*n*Wny(1);
                            if  ToggleInnerFluid
                                zRi = sqrt(AngularFrequency2/cFi2-k2)*Ri;
                                Znzi = besselj(n,zRi);
                                M(1,7) = a1*Densityi*Znzi;
                                M(4,1) = dZnxR(1);
                                M(4,2) = dWnxR(1);
                                M(4,3) = -k1*dZnyR(1);
                                M(4,4) = -k1*dWnyR(1);
                                M(4,5) = -n*Zny(1);
                                M(4,6) = -n*Wny(1);
                                M(4,7) = n*Znzi-zRi*besselj(n+1,zRi);
                                M(5,1) = dZnxR(2)+(a2-k2*Ro2-n2)*Znx(2);
                                M(5,2) = dWnxR(2)+(a2-k2*Ro2-n2)*Wnx(2);
                                M(5,3) = -k1*(dZnyR(2)+(y2*Ro2-n2)*Zny(2));
                                M(5,4) = -k1*(dWnyR(2)+(y2*Ro2-n2)*Wny(2));
                                M(5,5) = n*(dZnyR(2)-Zny(2));
                                M(5,6) = n*(dWnyR(2)-Wny(2));
                                M(6,1) = 2*n*(dZnxR(2)-Znx(2));
                                M(6,2) = 2*n*(dWnxR(2)-Wnx(2));
                                M(6,3) = 2*k1*n*(Zny(2)-dZnyR(2));
                                M(6,4) = 2*k1*n*(Wny(2)-dWnyR(2));
                                M(6,5) = 2*dZnyR(2)+(y2*Ro2-2*n2)*Zny(2);
                                M(6,6) = 2*dWnyR(2)+(y2*Ro2-2*n2)*Wny(2);
                                M(7,1) = 2*k1*dZnxR(2);
                                M(7,2) = 2*k1*dWnxR(2);
                                M(7,3) = a3*dZnyR(2);
                                M(7,4) = a3*dWnyR(2);
                                M(7,5) = -k1*n*Zny(2);
                                M(7,6) = -k1*n*Wny(2);
                            else
                                M(4,1) = dZnxR(2)+(a2-k2*Ro2-n2)*Znx(2);
                                M(4,2) = dWnxR(2)+(a2-k2*Ro2-n2)*Wnx(2);
                                M(4,3) = -k1*(dZnyR(2)+(y2*Ro2-n2)*Zny(2));
                                M(4,4) = -k1*(dWnyR(2)+(y2*Ro2-n2)*Wny(2));
                                M(4,5) = n*(dZnyR(2)-Zny(2));
                                M(4,6) = n*(dWnyR(2)-Wny(2));
                                M(5,1) = 2*n*(dZnxR(2)-Znx(2));
                                M(5,2) = 2*n*(dWnxR(2)-Wnx(2));
                                M(5,3) = 2*k1*n*(Zny(2)-dZnyR(2));
                                M(5,4) = 2*k1*n*(Wny(2)-dWnyR(2));
                                M(5,5) = 2*dZnyR(2)+(y2*Ro2-2*n2)*Zny(2);
                                M(5,6) = 2*dWnyR(2)+(y2*Ro2-2*n2)*Wny(2);
                                M(6,1) = 2*k1*dZnxR(2);
                                M(6,2) = 2*k1*dWnxR(2);
                                M(6,3) = a3*dZnyR(2);
                                M(6,4) = a3*dWnyR(2);
                                M(6,5) = -k1*n*Zny(2);
                                M(6,6) = -k1*n*Wny(2);
                            end
                        end
                        Y(j,l) = det(M);
                    end
                    if  k == 1
                        Y(j+1,1) = Y(j,3);
                    end
                    if  sign(Y(j,1)) ~= sign(Y(j,2)) && ((j == 1 && abs(Y(j,2)) < abs(Y(j,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))
                        Frequency = [Frequency(1) (Frequency(1)+Frequency(2))/2 Frequency(2)];
                        Y(j,3) = Y(j,2);
                    elseif sign(Y(j,2)) ~= sign(Y(j,3)) && ((j == 1 && abs(Y(j,2)) < abs(Y(j,1))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))
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
    B{end}(:,2:4) = [B{end}(:,1)*[1e-3 Ro-Ri] PhaseVelocityRange(1:size(B{end},1))'];
    if  ~B{end}(i)
        B{end}(i,:) = [H(p)*[1 1e-3 Ro-Ri] PhaseVelocityLimit];
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
    if  n == 0 && p > length(X0)
        A{p} = [flipud(B{p-length(X0)});A{p}];
    elseif n == 1 && p > 1
        A{p} = [flipud(B{p-1});A{p}];
    elseif n > 1
        A{p} = [flipud(B{p});A{p}];
    end
    A{p}(:,4) = A{p}(:,4)/1e3;
    A{p}(:,8:9) = 1;
% A{p}(:,5) = smooth(((A{p}(:,4)).^2)./(A{p}(:,4)-A{p}(:,1).*differentiate(fit(A{p}(:,1),A{p}(:,4),'cubicspline'),A{p}(:,1))));
end
end
function [Y,Yc] = Computer(k,Ro,Ri,Ro2,Ri2,kL2,kT2,kFi2,Densityi,YSize)
    k2 = k.^2;
    y2 = kT2-k2;
    x = sqrt(k2-kL2);
    y = sqrt(-y2);
    xRo = x*Ro;
    yRo = y*Ro;
    xRi = x*Ri;
    yRi = y*Ri;
    Z0xiRi2 = Ri2*besseli(0,xRi);
    Z0xoRo2 = Ro2*besseli(0,xRo);
    Z0yiRi2 = Ri2*besseli(0,yRi);
    Z0yoRo2 = Ro2*besseli(0,yRo);
    W0xiRi2 = Ri2*besselk(0,xRi);
    W0xoRo2 = Ro2*besselk(0,xRo);
    W0yiRi2 = Ri2*besselk(0,yRi);
    W0yoRo2 = Ro2*besselk(0,yRo);
    Z1xiRi = -xRi.*besseli(1,xRi);
    Z1xoRo = -xRo.*besseli(1,xRo);
    Z1yiRi = -yRi.*besseli(1,yRi);
    Z1yoRo = -yRo.*besseli(1,yRo);
    W1xiRi = xRi.*besselk(1,xRi);
    W1xoRo = xRo.*besselk(1,xRo);
    W1yiRi = yRi.*besselk(1,yRi);
    W1yoRo = yRo.*besselk(1,yRo);
    zRi = sqrt(kFi2-k2)*Ri;
    a1 = kT2/2;
    a2 = a1-k2;
    a3 = k2-y2;
    a4 = a1*Ri2*Densityi*besselj(0,zRi);
    Z1ziRi = -zRi.*besselj(1,zRi);
    Yc = NaN(1,YSize);
    for j = 1:YSize
        M(1,1) = a2(j)*Z0xiRi2(j)-Z1xiRi(j);
        M(1,2) = a2(j)*W0xiRi2(j)-W1xiRi(j);
        M(1,3) = k(j)*(Z1yiRi(j)-y2(j)*Z0yiRi2(j));
        M(1,4) = k(j)*(W1yiRi(j)-y2(j)*W0yiRi2(j));
        M(1,5) = a4(j);
        M(2,1) = -2*k(j)*Z1xiRi(j);
        M(2,2) = -2*k(j)*W1xiRi(j);
        M(2,3) = a3(j)*Z1yiRi(j);
        M(2,4) = a3(j)*W1yiRi(j);
        M(3,1) = -Z1xiRi(j);
        M(3,2) = -W1xiRi(j);
        M(3,3) = k(j)*Z1yiRi(j);
        M(3,4) = k(j)*W1yiRi(j);
        M(3,5) = Z1ziRi(j);
        M(4,1) = a2(j)*Z0xoRo2(j)-Z1xoRo(j);
        M(4,2) = a2(j)*W0xoRo2(j)-W1xoRo(j);
        M(4,3) = k(j)*(Z1yoRo(j)-y2(j)*Z0yoRo2(j));
        M(4,4) = k(j)*(W1yoRo(j)-y2(j)*W0yoRo2(j));
        M(5,1) = -2*k(j)*Z1xoRo(j);
        M(5,2) = -2*k(j)*W1xoRo(j);
        M(5,3) = a3(j)*Z1yoRo(j);
        M(5,4) = a3(j)*W1yoRo(j);
        Yc(j) = det(M);
    end
    Y = abs(Yc);
end