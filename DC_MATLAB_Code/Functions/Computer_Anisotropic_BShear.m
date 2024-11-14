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
function BShear = Computer_Anisotropic_BShear(Multithreading,Q1,Q2,Q3,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth)        
Bisections = 24;

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
BShear{1} = [];
if  ~Multithreading
    if  HigherOrderModes && any(H)
        for p = 1:length(H)+1
            g(p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
            g1(p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
        end
    else
        g = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
    end
end
if  ~Hybrid
    BShear{1}(:,1) = FrequencyRange;
    BShear{1}(:,2) = FrequencyRange/1e3;
    BShear{1}(:,3) = FrequencyRange*PlateThickness;
    BShear{1}(:,4) = sqrt(c{1}(6,6)/Material{1}.Density)/1e3;
    BShear{1}(:,6) = 0;
    line(ax,BShear{1}(:,1),BShear{1}(:,4),'LineStyle','--','color',[.5 0 1])
    if  Multithreading
        send(Q3,{BShear{1},'--',[.5 0 1]})
    else
        line(ax,BShear{1}(:,1),BShear{1}(:,4),'LineStyle','--','color',[.5 0 1])
    end
else
    for i = 1:length(FrequencyRange)
        if  Stop
            return
        end
        AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
        X(i,1) = 0;
        if  i == 1
            SweepRange = 10e3:-10:100;
        elseif i == 2
            SweepRange = X(i-1)+10:-1:X(i-1)-10;
        elseif i >= 3 && i <= 4
            SweepRange = [X(i-1)+1*abs(X(i-2)-X(i-1)) X(i-1)-100*abs(X(i-2)-X(i-1))];
        else
            SweepRange = [X(i-1)+1*abs(X(i-2)-X(i-1)) X(i-1)-2*ShearPhaseVelocitySweepRange*abs(X(i-2)-X(i-1))];
            if  SweepRange(end) < 0
                SweepRange(end) = 0;
            end
        end
        for o = 0:PhaseVelocitySections
            if  o > 0
                SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
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
                        Y(j,l) = M(2,1);
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
                    break
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
    BShear{1}(:,1) = FrequencyRange;
    BShear{1}(:,2) = FrequencyRange/1e3;
    BShear{1}(:,3) = FrequencyRange*PlateThickness;
    BShear{1}(:,4) = fillmissing(X,'spline')/1e3;
    BShear{1}(:,6) = 0;
end
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
            X(i,1) = 0;
            Neighbors = [];
            for j = 1:length(BShear)
                if  i <= height(BShear{j})
                    Neighbors(j) = BShear{j}(i,4)*1e3;
                end
            end
            for q = 1:2
                if  q == 1
                    if  all(X == 0)
                        SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                    elseif isscalar(find(X > 0))
                        SweepRange = [X(i-1) max(Neighbors(Neighbors < X(i-1)))+1];
                    else
                        if  ~Hybrid
                            Factor = 2;
                        else
                            Factor = ShearPhaseVelocitySweepRange;
                        end
                        SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
                        if  SweepRange(end) < max(Neighbors(Neighbors < X(i-1)))+1
                            SweepRange(end) = max(Neighbors(Neighbors < X(i-1)))+1;
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
                                    Y(j,l) = M(2,1);
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
            BShear{p+1}(:,1) = FrequencyRange(1:length(X));
            BShear{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
            BShear{p+1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
            BShear{p+1}(:,4) = fillmissing(X,'spline')/1e3;
        else
            BShear{p+1} = BShear{p};
            X1{p+1}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(BShear{p+1}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop
                return
            end
            X1{p+1}(i,1) = 0;
            if  i == 1
                SweepRange = [BShear{p+1}(MaxInd,1)-10/PlateThickness/1e3 BShear{p+1}(MaxInd,1)+2/PlateThickness/1e3];
            elseif i > 1
                SweepRange = [X1{p+1}(i-1)-10/PlateThickness/1e3 X1{p+1}(i-1)+2/PlateThickness/1e3];
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
                            Y(j,l) = M(2,1);
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
                        X1{p+1}(i,1) = Frequency(2);
                        break
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
        X1{p+1}(:,2) = X1{p+1}/1e3;
        X1{p+1}(:,3) = X1{p+1}(:,1)*PlateThickness;
        X1{p+1}(:,4) = PhaseVelocityRange(1:height(X1{p+1}))/1e3;
        if  X1{p+1}(end,1) == 0
            X1{p+1}(end,1) = H(p);
            X1{p+1}(end,2) = H(p)/1e3;
            X1{p+1}(end,3) = H(p)*PlateThickness;
            X1{p+1}(end,4) = PhaseVelocityLimit/1e3;
        end
        X1{p+1}(X1{p+1}(:,1) == 0,:) = []; 
        X1{p+1} = flipud(X1{p+1});
    end
    X1(MissingModes == 1) = [];
    BShear(MissingModes == 1) = [];
    for p = 2:length(BShear)
        BShear{p}(BShear{p}(:,4) == 0,:) = [];
        BShear{p} = vertcat(X1{p},BShear{p});
        BShear{p}(:,6) = 0;
    end
end