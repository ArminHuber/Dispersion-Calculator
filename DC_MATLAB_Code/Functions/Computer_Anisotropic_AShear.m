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
function AShear = Computer_Anisotropic_AShear(Multithreading,Q1,Q2,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth)        
Bisections = 24;

%#ok<*AGROW>
%#ok<*GVMIS> 
global Stop 
Stop = 0;
AShear{1} = [];
if  SuperLayerSize == 1
    for p = 1:length(H)
        i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange);
        if  isempty(i)
            continue
        end
        FrqRange = FrequencyRange(i)*1e3;
        AShear{p}(i,4) = sqrt(c{1}(6,6)./(Material{1}.Density-c{1}(4,4)*((2*p-1)./(2*FrqRange*PlateThickness)).^2))/1e3;
        AShear{p}(:,1) = FrequencyRange(1:height(AShear{p}));
        AShear{p}(:,2) = FrequencyRange(1:height(AShear{p}))/1e3;
        AShear{p}(:,3) = FrequencyRange(1:height(AShear{p}))*PlateThickness;
        PhaseVelocityRange = max(AShear{p}(:,4))*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        X1{p}(:,1) = (2*p-1)*PhaseVelocityRange.*sqrt(c{1}(4,4)./(Material{1}.Density*PhaseVelocityRange.^2-c{1}(6,6)))/PlateThickness/2e3;
        X1{p}(:,2) = X1{p}/1e3;
        X1{p}(:,3) = X1{p}(:,1)*PlateThickness;
        X1{p}(:,4) = PhaseVelocityRange(1:height(X1{p}))/1e3;
        X1{p}(X1{p}(:,1) == 0,:) = [];
        X1{p} = flipud(X1{p});
    end
else
    if  ~Multithreading
        for p = 1:length(H)
            g(p) = animatedline(ax,'LineStyle','--','color','b');
            g1(p) = animatedline(ax,'LineStyle','--','color','b');
        end
    end
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
            X(i,1) = 0;
            Neighbors = [];
            if  p > 1
                for j = 1:length(AShear)
                    if  i <= height(AShear{j})
                        Neighbors(j) = AShear{j}(i,4)*1e3;
                    end
                end
            end
            for q = 1:2
                if  q == 1
                    if  all(X == 0)
                        if  p == 1
                            SweepRange = [PhaseVelocityLimit 2e3];
                        else
                            SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                        end
                    elseif isscalar(find(X > 0))
                        if  p == 1
                            if  ~Hybrid
                                SweepRange = [X(i-1) sqrt(c{1}(6,6)/Material{1}.Density)+1];
                            else
                                SweepRange = [X(i-1) 2e3];
                            end
                        else
                            SweepRange = [X(i-1) max(Neighbors(Neighbors < X(i-1)))+1];
                        end
                    else
                        if  ~Hybrid
                            Factor = 2;
                        else
                            Factor = ShearPhaseVelocitySweepRange;
                        end
                        SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
                        if  p == 1
                            if  ~Hybrid && SweepRange(end) < sqrt(c{1}(6,6)/Material{1}.Density)
                                SweepRange(end) = sqrt(c{1}(6,6)/Material{1}.Density);
                            end
                        else
                            if  SweepRange(end) < max(Neighbors(Neighbors < X(i-1)))+1
                                SweepRange(end) = max(Neighbors(Neighbors < X(i-1)))+1;
                            end
                        end
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
                        SweepRange(end) = X(i-1)-ShearPhaseVelocitySweepRange*abs(X(i-2)-X(i-1));
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
                                    Wavenumber = AngularFrequency/PhaseVelocity(l);
                                    for m = 1:SuperLayerSize
                                        Alpha = sqrt((Material{m}.Density*PhaseVelocity(l)^2-c{m}(6,6))/c{m}(4,4));
                                        D = Alpha*c{m}(4,4);
                                        Gamma = Alpha*Wavenumber*LayerThicknesses(m);
                                        L{m} = [cos(Gamma) 1i*sin(Gamma)/D;1i*D*sin(Gamma) cos(Gamma)];
                                    end
                                    M = L{1};
                                    for m = 2:SuperLayerSize
                                        M = M*L{m};
                                    end
                                    if  Repetitions > 1
                                        M = M^Repetitions;
                                    end
                                    Y(j,l) = M(2,2);
                                end
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
                            X(i,1) = PhaseVelocity(2);
                            Misses(i) = 0;
                            BelowCutoff(i) = 0;
                            break
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
        end
        if  all(X == 0)
            MissingModes(p) = 1;
        end
        if  ~MissingModes(p)
            X(Misses(1:length(X)) == 1) = NaN;
            AShear{p}(:,1) = FrequencyRange(1:length(X));
            AShear{p}(:,2) = FrequencyRange(1:length(X))/1e3;
            AShear{p}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
            AShear{p}(:,4) = fillmissing(X,'spline')/1e3;
        else
            if  p == 1
                AShear{1}(1,1) = 0;
            else
                AShear{p} = AShear{p-1};
            end
            X1{p}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(AShear{p}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop
                return
            end
            X1{p}(i,1) = 0;
            if  i == 1
                SweepRange = [AShear{p}(MaxInd,1)-10/PlateThickness/1e3 AShear{p}(MaxInd,1)+2/PlateThickness/1e3];
            else
                SweepRange = [X1{p}(i-1)-10/PlateThickness/1e3 X1{p}(i-1)+2/PlateThickness/1e3];
            end
            for o = 0:FrequencySections
                if  o > 0
                    SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
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
                                Alpha = sqrt((Material{m}.Density*PhaseVelocityRange(i)^2-c{m}(6,6))/c{m}(4,4));
                                D = Alpha*c{m}(4,4);
                                Gamma = Alpha*Wavenumber*LayerThicknesses(m);
                                L{m} = [cos(Gamma) 1i*sin(Gamma)/D;1i*D*sin(Gamma) cos(Gamma)];
                            end
                            M = L{1};
                            for m = 2:SuperLayerSize
                                M = M*L{m};
                            end
                            if  Repetitions > 1
                                M = M^Repetitions;
                            end
                            Y(j,l) = M(2,2);
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
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
                    end
                    if  Frequency(2) > 0
                        X1{p}(i,1) = Frequency(2);
                        break
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
        X1{p}(:,2) = X1{p}/1e3;
        X1{p}(:,3) = X1{p}(:,1)*PlateThickness;
        X1{p}(:,4) = PhaseVelocityRange(1:height(X1{p}))/1e3;
        if  X1{p}(end,1) == 0
            X1{p}(end,1) = H(p);
            X1{p}(end,2) = H(p)/1e3;
            X1{p}(end,3) = H(p)*PlateThickness;
            X1{p}(end,4) = PhaseVelocityLimit/1e3;
        end
        X1{p}(X1{p}(:,1) == 0,:) = [];
        X1{p} = flipud(X1{p});
    end
    X1(MissingModes == 1) = [];
    AShear(MissingModes == 1) = [];
end
for p = 1:length(AShear)
    AShear{p}(AShear{p}(:,4) == 0,:) = [];
    AShear{p} = vertcat(X1{p},AShear{p});
    AShear{p}(:,6) = 0;
end