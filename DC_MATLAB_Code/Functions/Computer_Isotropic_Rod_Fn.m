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
function F = Computer_Isotropic_Rod_Fn(Multithreading,Q1,Q2,ax,Material,FrequencyRange,PhaseVelocitySections,R,n,H,LineColor,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth) 
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
F{1} = [];
n2 = n^2;
R2 = R^2;
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
        X(i,1) = 0;
        Neighbors = [];
        for j = 1:length(F)
            if  i <= height(F{j})
                Neighbors(j) = F{j}(i,4)*1e3;
            end
        end
        if  isempty(Neighbors)
            Neighbors = Material.TransverseVelocity;
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
                if  p > 1 && any(Neighbors < X(i-1)) && SweepRange(end) < max(Neighbors(Neighbors < X(i-1)))+1
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
                    for m = 1:length(Neighbors)
                        if  SweepRange(j) > Neighbors(m) && SweepRange(j+1) < Neighbors(m)
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
                                y2 = kT2-k2;
                                xR = sqrt(kL2-k2)*R;
                                yR = sqrt(y2)*R;
                                Jnx = besselj(n,xR);
                                Jny = besselj(n,yR);
                                dJnxR = n*Jnx-xR*besselj(n+1,xR);
                                dJnyR = n*Jny-yR*besselj(n+1,yR);
                                M(1,1) = dJnxR+(.5*kT2*R2-k2*R2-n2)*Jnx;
                                M(1,2) = -k1*(dJnyR+(y2*R2-n2)*Jny);
                                M(1,3) = n*(dJnyR-Jny);
                                M(2,1) = 2*n*(dJnxR-Jnx);
                                M(2,2) = 2*k1*n*(Jny-dJnyR);
                                M(2,3) = 2*dJnyR+(y2*R2-2*n2)*Jny;
                                M(3,1) = 2*k1*dJnxR;
                                M(3,2) = (y2-k2)*dJnyR;
                                M(3,3) = -k1*n*Jny;
                                Y(j,l) = det(M/Jnx/Jny);
                            end
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  (j == 1 && abs(Y(j,2)) < 1e-5*10^-n) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
                            PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif (j == 1 && abs(Y(j,2)) < 1e-5*10^-n) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
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
        if  FrequencyResolution*length(BelowCutoff(BelowCutoff == 1)) > BelowCutoffWidth/R/2e3
            MissingModes(p) = 1;
            break
        end
        if  i > MissingSamples && all(Misses(end-MissingSamples:end))
            X(end-MissingSamples:end) = [];
            Misses(end-MissingSamples:end) = 0;
            break
        end
    end
    if  all(X == 0)
        MissingModes(p) = 1;
    end
    if  ~MissingModes(p)
        X(Misses(1:length(X)) == 1) = NaN;
        F{p}(:,1) = FrequencyRange(1:length(X));
        F{p}(:,2) = FrequencyRange(1:length(X))/1e3;
        F{p}(:,3) = FrequencyRange(1:length(X))*2*R;
        F{p}(:,4) = fillmissing(X,'spline')/1e3;
    else
        F{p} = F{p};
        X1{p}(1) = 0;
        continue
    end
    if  Multithreading
        send(Q1,{F{p},'-',LineColor})
    else
        line(ax,F{p}((F{p}(:,4) ~= 0),1),F{p}((F{p}(:,4) ~= 0),4),'color',LineColor)
        drawnow limitrate
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
                SweepRange = [F{p}(MaxInd,1)-2*FrequencyOffset/R/2e3 F{p}(MaxInd,1)+10/R/2e3];
            else
                SweepRange = [F{p}(MaxInd,1)-FrequencyOffset/R/2e3 F{p}(MaxInd,1)+10/R/2e3];
            end
        elseif i == 2
            if  p == 1
                SweepRange = [X1{p}(i-1)-2*FrequencyOffset/R/2e3 X1{p}(i-1)+10/R/2e3];
            else
                SweepRange = [X1{p}(i-1)-FrequencyOffset/R/2e3 X1{p}(i-1)+10/R/2e3];
            end
        else
            SweepRange = X1{p}(i-1)-2*abs(X1{p}(i-2)-X1{p}(i-1));
            if  X1{p}(i-1) < X1{p}(i-2)
                SweepRange(2) = X1{p}(i-1)+2/R/2e3;
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
                        kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
                        y2 = kT2-k2;
                        xR = sqrt((AngularFrequency/Material.LongitudinalVelocity)^2-k2)*R;
                        yR = sqrt(y2)*R;
                        Jnx = besselj(n,xR);
                        Jny = besselj(n,yR);
                        dJnxR = n*Jnx-xR*besselj(n+1,xR);
                        dJnyR = n*Jny-yR*besselj(n+1,yR);
                        M(1,1) = dJnxR+(.5*kT2*R2-k2*R2-n2)*Jnx;
                        M(1,2) = -k1*(dJnyR+(y2*R2-n2)*Jny);
                        M(1,3) = n*(dJnyR-Jny);
                        M(2,1) = 2*n*(dJnxR-Jnx);
                        M(2,2) = 2*k1*n*(Jny-dJnyR);
                        M(2,3) = 2*dJnyR+(y2*R2-2*n2)*Jny;
                        M(3,1) = 2*k1*dJnxR;
                        M(3,2) = (y2-k2)*dJnyR;
                        M(3,3) = -k1*n*Jny;
                        Y(j,l) = det(M/Jnx/Jny);
                    end
                    if  k == 1
                        Y(j+1,1) = Y(j,3);
                    end
                    if  (j == 1 && abs(Y(j,2)) < 1e-5) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
                        Frequency = [Frequency(1) Frequency(1)+(Frequency(2)-Frequency(1))/2 Frequency(2)];
                        Y(j,3) = Y(j,2);
                    elseif (j == 1 && abs(Y(j,2)) < 1e-5) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
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
    end
    X1{p}(:,2) = X1{p}(:,1)/1e3;
    X1{p}(:,3) = X1{p}(:,1)*2*R;
    X1{p}(:,4) = PhaseVelocityRange(1:height(X1{p}))/1e3;
    if  X1{p}(end,1) == 0
        X1{p}(end,1) = H(p);
        X1{p}(end,2) = H(p)/1e3;
        X1{p}(end,3) = H(p)*2*R;
        X1{p}(end,4) = PhaseVelocityLimit/1e3;
    end
    X1{p}(X1{p}(:,1) == 0,:) = [];
    X1{p} = flipud(X1{p});
    if  Multithreading
        send(Q2,{X1{p},'-',LineColor})
    else
        line(ax,X1{p}(:,1),X1{p}(:,4),'color',LineColor)
        drawnow limitrate
    end
end
X1(MissingModes == 1) = [];
F(MissingModes == 1) = [];
for p = 1:length(F)
    F{p}(F{p}(:,4) == 0,:) = [];
    F{p} = vertcat(X1{p},F{p});
    F{p}(:,6) = 0;
end