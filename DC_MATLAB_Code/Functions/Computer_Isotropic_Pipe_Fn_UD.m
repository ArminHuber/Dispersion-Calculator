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
function F = Computer_Isotropic_Pipe_Fn_UD(Multithreading,Q1,Q2,ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,n,H,LineColor,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth)        
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
n2 = n^2;
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityi = InnerFluid.Density/Material.Density;
Y = zeros(5^PhaseVelocitySections+1,3);
MissingModes(length(H)) = 0;
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
        if  p == 1
            Neighbors = 0;
        else
            for j = 1:length(F)
                if  FrequencyRange(i) <= F{j}(end,1)
                    z = find(abs(FrequencyRange(i)-F{j}(:,1)) == min(abs(FrequencyRange(i)-F{j}(:,1))));
                    Neighbors(j) = F{j}(z,4)*1e3; %#ok<*FNDSB>
                else
                    Neighbors(j) = 0;
                end
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
                if  q == 1 && o == 0 && (any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) || (SweepRange(1) > InnerFluid.Velocity && SweepRange(2) < InnerFluid.Velocity) || (SweepRange(1) > Material.TransverseVelocity && SweepRange(2) < Material.TransverseVelocity)) % Neighbors cannot be excluded because SweepRange has only 2 elements
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
                                    Znxi = besselj(n,xRi);
                                    Znxo = besselj(n,xRo);
                                    Znyi = besselj(n,yRi);
                                    Znyo = besselj(n,yRo);
                                    Wnxi = bessely(n,xRi);
                                    Wnxo = bessely(n,xRo);
                                    Wnyi = bessely(n,yRi);
                                    Wnyo = bessely(n,yRo);
                                    dZnxiRi = n*Znxi-xRi*besselj(n+1,xRi);
                                    dZnxoRo = n*Znxo-xRo*besselj(n+1,xRo);
                                    dZnyiRi = n*Znyi-yRi*besselj(n+1,yRi);
                                    dZnyoRo = n*Znyo-yRo*besselj(n+1,yRo);
                                    dWnxiRi = n*Wnxi-xRi*bessely(n+1,xRi);
                                    dWnxoRo = n*Wnxo-xRo*bessely(n+1,xRo);
                                    dWnyiRi = n*Wnyi-yRi*bessely(n+1,yRi);
                                    dWnyoRo = n*Wnyo-yRo*bessely(n+1,yRo);
                                elseif SweepRange(j) < Material.LongitudinalVelocity && SweepRange(j) > Material.TransverseVelocity
                                    Znxi = besseli(n,xRi);
                                    Znxo = besseli(n,xRo);
                                    Znyi = besselj(n,yRi);
                                    Znyo = besselj(n,yRo);
                                    Wnxi = besselk(n,xRi);
                                    Wnxo = besselk(n,xRo);
                                    Wnyi = bessely(n,yRi);
                                    Wnyo = bessely(n,yRo);
                                    dZnxiRi = n*Znxi+xRi*besseli(n+1,xRi);
                                    dZnxoRo = n*Znxo+xRo*besseli(n+1,xRo);
                                    dZnyiRi = n*Znyi-yRi*besselj(n+1,yRi);
                                    dZnyoRo = n*Znyo-yRo*besselj(n+1,yRo);
                                    dWnxiRi = n*Wnxi-xRi*besselk(n+1,xRi);
                                    dWnxoRo = n*Wnxo-xRo*besselk(n+1,xRo);
                                    dWnyiRi = n*Wnyi-yRi*bessely(n+1,yRi);
                                    dWnyoRo = n*Wnyo-yRo*bessely(n+1,yRo);
                                elseif SweepRange(j) < Material.TransverseVelocity
                                    Znxi = besseli(n,xRi);
                                    Znxo = besseli(n,xRo);
                                    Znyi = besseli(n,yRi);
                                    Znyo = besseli(n,yRo);
                                    Wnxi = besselk(n,xRi);
                                    Wnxo = besselk(n,xRo);
                                    Wnyi = besselk(n,yRi);
                                    Wnyo = besselk(n,yRo);
                                    dZnxiRi = n*Znxi+xRi*besseli(n+1,xRi);
                                    dZnxoRo = n*Znxo+xRo*besseli(n+1,xRo);
                                    dZnyiRi = n*Znyi+yRi*besseli(n+1,yRi);
                                    dZnyoRo = n*Znyo+yRo*besseli(n+1,yRo);
                                    dWnxiRi = n*Wnxi-xRi*besselk(n+1,xRi);
                                    dWnxoRo = n*Wnxo-xRo*besselk(n+1,xRo);
                                    dWnyiRi = n*Wnyi-yRi*besselk(n+1,yRi);
                                    dWnyoRo = n*Wnyo-yRo*besselk(n+1,yRo);
                                end
                                zRi = sqrt(kInnerFluid2-k2)*Ri;
                                Znzi = besselj(n,zRi);
                                dZnziRi = n*Znzi-zRi*besselj(n+1,zRi);
                                M(1,1) = dZnxiRi+(.5*kT2*Ri2-k2*Ri2-n2)*Znxi;
                                M(1,2) = dWnxiRi+(.5*kT2*Ri2-k2*Ri2-n2)*Wnxi;
                                M(1,3) = -k1*(dZnyiRi+(y2*Ri2-n2)*Znyi);
                                M(1,4) = -k1*(dWnyiRi+(y2*Ri2-n2)*Wnyi);
                                M(1,5) = n*(dZnyiRi-Znyi);
                                M(1,6) = n*(dWnyiRi-Wnyi);
                                M(1,7) = .5*kT2*Ri2*Densityi*Znzi;
                                M(2,1) = 2*n*(dZnxiRi-Znxi);
                                M(2,2) = 2*n*(dWnxiRi-Wnxi);
                                M(2,3) = 2*k1*n*(Znyi-dZnyiRi);
                                M(2,4) = 2*k1*n*(Wnyi-dWnyiRi);
                                M(2,5) = 2*dZnyiRi+(y2*Ri2-2*n2)*Znyi;
                                M(2,6) = 2*dWnyiRi+(y2*Ri2-2*n2)*Wnyi;
                                M(3,1) = 2*k1*dZnxiRi;
                                M(3,2) = 2*k1*dWnxiRi;
                                M(3,3) = (y2-k2)*dZnyiRi;
                                M(3,4) = (y2-k2)*dWnyiRi;
                                M(3,5) = -k1*n*Znyi;
                                M(3,6) = -k1*n*Wnyi;
                                M(4,1) = dZnxiRi;
                                M(4,2) = dWnxiRi;
                                M(4,3) = -k1*dZnyiRi;
                                M(4,4) = -k1*dWnyiRi;
                                M(4,5) = -n*Znyi;
                                M(4,6) = -n*Wnyi;
                                M(4,7) = dZnziRi;
                                M(5,1) = dZnxoRo+(.5*kT2*Ro2-k2*Ro2-n2)*Znxo;
                                M(5,2) = dWnxoRo+(.5*kT2*Ro2-k2*Ro2-n2)*Wnxo;
                                M(5,3) = -k1*(dZnyoRo+(y2*Ro2-n2)*Znyo);
                                M(5,4) = -k1*(dWnyoRo+(y2*Ro2-n2)*Wnyo);
                                M(5,5) = n*(dZnyoRo-Znyo);
                                M(5,6) = n*(dWnyoRo-Wnyo);
                                M(6,1) = 2*n*(dZnxoRo-Znxo);
                                M(6,2) = 2*n*(dWnxoRo-Wnxo);
                                M(6,3) = 2*k1*n*(Znyo-dZnyoRo);
                                M(6,4) = 2*k1*n*(Wnyo-dWnyoRo);
                                M(6,5) = 2*dZnyoRo+(y2*Ro2-2*n2)*Znyo;
                                M(6,6) = 2*dWnyoRo+(y2*Ro2-2*n2)*Wnyo;
                                M(7,1) = 2*k1*dZnxoRo;
                                M(7,2) = 2*k1*dWnxoRo;
                                M(7,3) = (y2-k2)*dZnyoRo;
                                M(7,4) = (y2-k2)*dWnyoRo;
                                M(7,5) = -k1*n*Znyo;
                                M(7,6) = -k1*n*Wnyo;
                                if  ~mod(n,2) 
                                    Y(j,l) = real(det(M));
                                else
                                    Y(j,l) = real(det(M)/Znzi);
                                end
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
            MissingModes(p) = 1;
            break
        end
        if  i > MissingSamples && all(Misses(end-MissingSamples:end))
            X(end-MissingSamples:end) = [];
            Misses(end-MissingSamples:end) = 0;
            break
        end
        if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
            send(Q1,[FrequencyRange(i),X(i)/1e3,p])
        elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
            addpoints(g(p),FrequencyRange(i),X(i)/1e3);
            drawnow limitrate
        end
% String = ['n = ',num2str(n),' p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
    end
    if  all(X == 0)
        MissingModes(p) = 1;
    end
    if  ~MissingModes(p)
        X(Misses(1:length(X)) == 1) = NaN;
        F{p}(:,1) = FrequencyRange(1:length(X));
        F{p}(:,2) = FrequencyRange(1:length(X))/1e3;
        F{p}(:,3) = FrequencyRange(1:length(X))*(Ro-Ri);
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
            if  p == 1
                SweepRange = [F{p}(MaxInd,1)-2/4*FrequencyOffset/(Ro-Ri)/1e3 F{p}(MaxInd,1)+10/4/(Ro-Ri)/1e3];
            else
                SweepRange = [F{p}(MaxInd,1)-1/4*FrequencyOffset/(Ro-Ri)/1e3 F{p}(MaxInd,1)+10/4/(Ro-Ri)/1e3];
            end
        elseif i == 2
            if  p == 1
                SweepRange = [X1{p}(i-1)-2/4*FrequencyOffset/(Ro-Ri)/1e3 X1{p}(i-1)+10/4/(Ro-Ri)/1e3];
            else
                SweepRange = [X1{p}(i-1)-1/4*FrequencyOffset/(Ro-Ri)/1e3 X1{p}(i-1)+10/4/(Ro-Ri)/1e3];
            end
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
                        kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
                        kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
                        if  PhaseVelocityRange(i) > Material.LongitudinalVelocity
                            x = sqrt(kL2-k2);
                            y = sqrt(kT2-k2);
                        elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity && PhaseVelocityRange(i) > Material.TransverseVelocity
                            x = sqrt(k2-kL2);
                            y = sqrt(kT2-k2);
                        elseif PhaseVelocityRange(i) < Material.TransverseVelocity
                            x = sqrt(k2-kL2);
                            y = sqrt(k2-kT2);
                        end
                        y2 = kT2-k2;
                        xRo = x*Ro;
                        yRo = y*Ro;
                        xRi = x*Ri;
                        yRi = y*Ri;
                        if  PhaseVelocityRange(i) > Material.LongitudinalVelocity
                            Znxi = besselj(n,xRi);
                            Znxo = besselj(n,xRo);
                            Znyi = besselj(n,yRi);
                            Znyo = besselj(n,yRo);
                            Wnxi = bessely(n,xRi);
                            Wnxo = bessely(n,xRo);
                            Wnyi = bessely(n,yRi);
                            Wnyo = bessely(n,yRo);
                            dZnxiRi = n*Znxi-xRi*besselj(n+1,xRi);
                            dZnxoRo = n*Znxo-xRo*besselj(n+1,xRo);
                            dZnyiRi = n*Znyi-yRi*besselj(n+1,yRi);
                            dZnyoRo = n*Znyo-yRo*besselj(n+1,yRo);
                            dWnxiRi = n*Wnxi-xRi*bessely(n+1,xRi);
                            dWnxoRo = n*Wnxo-xRo*bessely(n+1,xRo);
                            dWnyiRi = n*Wnyi-yRi*bessely(n+1,yRi);
                            dWnyoRo = n*Wnyo-yRo*bessely(n+1,yRo);
                        elseif PhaseVelocityRange(i) < Material.LongitudinalVelocity && PhaseVelocityRange(i) > Material.TransverseVelocity
                            Znxi = besseli(n,xRi);
                            Znxo = besseli(n,xRo);
                            Znyi = besselj(n,yRi);
                            Znyo = besselj(n,yRo);
                            Wnxi = besselk(n,xRi);
                            Wnxo = besselk(n,xRo);
                            Wnyi = bessely(n,yRi);
                            Wnyo = bessely(n,yRo);
                            dZnxiRi = n*Znxi+xRi*besseli(n+1,xRi);
                            dZnxoRo = n*Znxo+xRo*besseli(n+1,xRo);
                            dZnyiRi = n*Znyi-yRi*besselj(n+1,yRi);
                            dZnyoRo = n*Znyo-yRo*besselj(n+1,yRo);
                            dWnxiRi = n*Wnxi-xRi*besselk(n+1,xRi);
                            dWnxoRo = n*Wnxo-xRo*besselk(n+1,xRo);
                            dWnyiRi = n*Wnyi-yRi*bessely(n+1,yRi);
                            dWnyoRo = n*Wnyo-yRo*bessely(n+1,yRo);
                        elseif PhaseVelocityRange(i) < Material.TransverseVelocity
                            Znxi = besseli(n,xRi);
                            Znxo = besseli(n,xRo);
                            Znyi = besseli(n,yRi);
                            Znyo = besseli(n,yRo);
                            Wnxi = besselk(n,xRi);
                            Wnxo = besselk(n,xRo);
                            Wnyi = besselk(n,yRi);
                            Wnyo = besselk(n,yRo);
                            dZnxiRi = n*Znxi+xRi*besseli(n+1,xRi);
                            dZnxoRo = n*Znxo+xRo*besseli(n+1,xRo);
                            dZnyiRi = n*Znyi+yRi*besseli(n+1,yRi);
                            dZnyoRo = n*Znyo+yRo*besseli(n+1,yRo);
                            dWnxiRi = n*Wnxi-xRi*besselk(n+1,xRi);
                            dWnxoRo = n*Wnxo-xRo*besselk(n+1,xRo);
                            dWnyiRi = n*Wnyi-yRi*besselk(n+1,yRi);
                            dWnyoRo = n*Wnyo-yRo*besselk(n+1,yRo);
                        end
                        zRi = sqrt((AngularFrequency/InnerFluid.Velocity)^2-k2)*Ri;
                        Znzi = besselj(n,zRi);
                        dZnziRi = n*Znzi-zRi*besselj(n+1,zRi);
                        M(1,1) = dZnxiRi+(.5*kT2*Ri2-k2*Ri2-n2)*Znxi;
                        M(1,2) = dWnxiRi+(.5*kT2*Ri2-k2*Ri2-n2)*Wnxi;
                        M(1,3) = -k1*(dZnyiRi+(y2*Ri2-n2)*Znyi);
                        M(1,4) = -k1*(dWnyiRi+(y2*Ri2-n2)*Wnyi);
                        M(1,5) = n*(dZnyiRi-Znyi);
                        M(1,6) = n*(dWnyiRi-Wnyi);
                        M(1,7) = .5*kT2*Ri2*Densityi*Znzi;
                        M(2,1) = 2*n*(dZnxiRi-Znxi);
                        M(2,2) = 2*n*(dWnxiRi-Wnxi);
                        M(2,3) = 2*k1*n*(Znyi-dZnyiRi);
                        M(2,4) = 2*k1*n*(Wnyi-dWnyiRi);
                        M(2,5) = 2*dZnyiRi+(y2*Ri2-2*n2)*Znyi;
                        M(2,6) = 2*dWnyiRi+(y2*Ri2-2*n2)*Wnyi;
                        M(3,1) = 2*k1*dZnxiRi;
                        M(3,2) = 2*k1*dWnxiRi;
                        M(3,3) = (y2-k2)*dZnyiRi;
                        M(3,4) = (y2-k2)*dWnyiRi;
                        M(3,5) = -k1*n*Znyi;
                        M(3,6) = -k1*n*Wnyi;
                        M(4,1) = dZnxiRi;
                        M(4,2) = dWnxiRi;
                        M(4,3) = -k1*dZnyiRi;
                        M(4,4) = -k1*dWnyiRi;
                        M(4,5) = -n*Znyi;
                        M(4,6) = -n*Wnyi;
                        M(4,7) = dZnziRi;
                        M(5,1) = dZnxoRo+(.5*kT2*Ro2-k2*Ro2-n2)*Znxo;
                        M(5,2) = dWnxoRo+(.5*kT2*Ro2-k2*Ro2-n2)*Wnxo;
                        M(5,3) = -k1*(dZnyoRo+(y2*Ro2-n2)*Znyo);
                        M(5,4) = -k1*(dWnyoRo+(y2*Ro2-n2)*Wnyo);
                        M(5,5) = n*(dZnyoRo-Znyo);
                        M(5,6) = n*(dWnyoRo-Wnyo);
                        M(6,1) = 2*n*(dZnxoRo-Znxo);
                        M(6,2) = 2*n*(dWnxoRo-Wnxo);
                        M(6,3) = 2*k1*n*(Znyo-dZnyoRo);
                        M(6,4) = 2*k1*n*(Wnyo-dWnyoRo);
                        M(6,5) = 2*dZnyoRo+(y2*Ro2-2*n2)*Znyo;
                        M(6,6) = 2*dWnyoRo+(y2*Ro2-2*n2)*Wnyo;
                        M(7,1) = 2*k1*dZnxoRo;
                        M(7,2) = 2*k1*dWnxoRo;
                        M(7,3) = (y2-k2)*dZnyoRo;
                        M(7,4) = (y2-k2)*dWnyoRo;
                        M(7,5) = -k1*n*Znyo;
                        M(7,6) = -k1*n*Wnyo;
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
% F{p}(:,5) = smooth(((F{p}(:,4)).^2)./(F{p}(:,4)-F{p}(:,1).*differentiate(fit(F{p}(:,1),F{p}(:,4),'cubicspline'),F{p}(:,1))));
    F{p}(:,6) = 0;
end