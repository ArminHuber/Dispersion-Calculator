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
function SH = Computer_Anisotropic_SH(Multithreading,Q1,Q2,Q3,ax,ModeFamily,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0)        
Bisections = 24;

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
SH{1} = [];
if  ModeFamily < 3
    N = 1; % S/B
elseif ModeFamily == 3
    N = 0; % A
end
if  ModeFamily == 1
    Color = 'r';
elseif ModeFamily == 2
    Color = [.5 0 1];
elseif ModeFamily == 3
    Color = 'b';
end
if  SuperLayerSize == 1
    if  ModeFamily == 1
        SH{1}(:,1) = FrequencyRange;
        SH{1}(:,2) = FrequencyRange/1e3;
        SH{1}(:,3) = FrequencyRange*PlateThickness;
        SH{1}(:,4) = XSH0/1e3;
        SH{1}(:,7) = 0;
    end
    if  HigherOrderModes && any(H)
        for p = 1:length(H)
            i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange);
            if  isempty(i)
                continue
            end
            FrqRange = FrequencyRange(i)*1e3;
            if  ModeFamily == 1
                SH{p+N}(i,4) = sqrt(c{1}(6,6)./(Material{1}.Density-c{1}(4,4)*(p./(FrqRange*PlateThickness)).^2))/1e3;
            elseif ModeFamily == 3
                SH{p+N}(i,4) = sqrt(c{1}(6,6)./(Material{1}.Density-c{1}(4,4)*((2*p-1)./(2*FrqRange*PlateThickness)).^2))/1e3;
            end
            SH{p+N}(:,1) = FrequencyRange(1:height(SH{p+N}));
            SH{p+N}(:,2) = FrequencyRange(1:height(SH{p+N}))/1e3;
            SH{p+N}(:,3) = FrequencyRange(1:height(SH{p+N}))*PlateThickness;
            PhaseVelocityRange = max(SH{p+N}(:,4))*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
            if  ModeFamily == 1
                X1{p+N}(:,1) = p*PhaseVelocityRange.*sqrt(c{1}(4,4)./(Material{1}.Density*PhaseVelocityRange.^2-c{1}(6,6)))/PlateThickness/1e3;
            elseif ModeFamily == 3
                X1{p+N}(:,1) = (2*p-1)*PhaseVelocityRange.*sqrt(c{1}(4,4)./(Material{1}.Density*PhaseVelocityRange.^2-c{1}(6,6)))/PlateThickness/2e3;
            end
            X1{p+N}(:,2) = X1{p+N}/1e3;
            X1{p+N}(:,3) = X1{p+N}(:,1)*PlateThickness;
            X1{p+N}(:,4) = PhaseVelocityRange(1:height(X1{p+N}))/1e3;
            X1{p+N}(X1{p+N}(:,1) == 0,:) = [];
            X1{p+N} = flipud(X1{p+N});
        end
    end
else
    if  ~Multithreading
        for p = 1:length(H)+N
            g(p) = animatedline(ax,'LineStyle','--','color',Color);
            g1(p) = animatedline(ax,'LineStyle','--','color',Color);
        end
    end
    if  ModeFamily < 3
        if  ~Hybrid
            SH{1}(:,1) = FrequencyRange;
            SH{1}(:,2) = FrequencyRange/1e3;
            SH{1}(:,3) = FrequencyRange*PlateThickness;
            SH{1}(:,4) = XSH0/1e3;
            SH{1}(:,7) = 0;
            if  Multithreading
                send(Q3,{SH{1},'--',Color})
            else
                line(ax,SH{1}(:,1),SH{1}(:,4),'LineStyle','--','color',Color)
            end
        else
            X = XSH0;
            if  Multithreading
                send(Q1,[FrequencyRange(1),XSH0/1e3,1])
            else
                addpoints(g(1),FrequencyRange(1),XSH0/1e3);
                drawnow limitrate
            end
            for i = 2:length(FrequencyRange)
                if  Stop
                    return
                end
                AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
                X(i,1) = 0;
                if  i == 2
                    SweepRange = X(i-1)+10:-1:X(i-1)-10;
                elseif i == 3 || i == 4
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
                                PhaseVelocity2 = PhaseVelocity(l)^2;
                                for m = 1:SuperLayerSize
                                    Alpha = sqrt((Material{m}.Density*PhaseVelocity2-c{m}(6,6))/c{m}(4,4));
                                    D = Alpha*c{m}(4,4);
                                    G = Wavenumber*Alpha*LayerThicknesses(m);
                                    CosG = cos(G);
                                    SinG = sin(G);
                                    L{m} = [CosG 1i*SinG./D;1i*SinG.*D CosG];
                                end
                                M{1} = L{1};
                                for m = 2:SuperLayerSize
                                    M{1} = M{1}*L{m};
                                end
                                for m = 1:length(Pattern)
                                    M{m+1} = M{m}*M{Pattern(m)};
                                end
                                if  ModeFamily < 3
                                    Y(j,l) = M{end}(2,1);
                                elseif ModeFamily == 3
                                    Y(j,l) = M{end}(2,2);
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
                            break
                        end
                    end
                    if  X(i) > 0
                        break
                    end
                end
                if  X(i) == 0
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
            SH{1}(:,1) = FrequencyRange(1:length(X));
            SH{1}(:,2) = FrequencyRange(1:length(X))/1e3;
            SH{1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
            SH{1}(:,4) = fillmissing(X,'spline')/1e3;
            SH{1}(:,7) = 0;
        end
    end
    if  HigherOrderModes && any(H)
        MissingModes(length(H)+N) = 0;
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
                for j = 1:length(SH)
                    if  i <= height(SH{j})
                        Neighbors(j) = SH{j}(i,4)*1e3;
                    end
                end
                for q = 1:2
                    if  q == 1
                        if  all(X == 0)
                            if  ModeFamily == 3 && p == 1
                                SweepRange = [PhaseVelocityLimit XSH0];
                            else
                                SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                            end
                        elseif isscalar(find(X > 0))
                            if  ModeFamily == 3 && p == 1
                                SweepRange = [X(i-1) XSH0];
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
                            if  ModeFamily == 3 && p == 1
                                if  SweepRange(end) < XSH0
                                    SweepRange(end) = XSH0;
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
                            for t = 1:length(Neighbors)
                                if  SweepRange(j) > Neighbors(t) && SweepRange(j+1) < Neighbors(t)
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
                                        PhaseVelocity2 = PhaseVelocity(l)^2;
                                        for m = 1:SuperLayerSize
                                            Alpha = sqrt((Material{m}.Density*PhaseVelocity2-c{m}(6,6))/c{m}(4,4));
                                            D = Alpha*c{m}(4,4);
                                            G = Wavenumber*Alpha*LayerThicknesses(m);
                                            CosG = cos(G);
                                            SinG = sin(G);
                                            L{m} = [CosG 1i*SinG./D;1i*SinG.*D CosG];
                                        end
                                        M{1} = L{1};
                                        for m = 2:SuperLayerSize
                                            M{1} = M{1}*L{m};
                                        end
                                        for m = 1:length(Pattern)
                                            M{m+1} = M{m}*M{Pattern(m)};
                                        end
                                        if  ModeFamily < 3
                                            Y(j,l) = M{end}(2,1);
                                        elseif ModeFamily == 3
                                            Y(j,l) = M{end}(2,2);
                                        end
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
                    MissingModes(p+N) = 1;
                    break
                end
                if  i > MissingSamples && all(Misses(end-MissingSamples:end))
                    X(end-MissingSamples:end) = [];
                    Misses(end-MissingSamples:end) = 0;
                    break
                end
                if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                    send(Q1,[FrequencyRange(i),X(i)/1e3,p+N])
                elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                    addpoints(g(p+N),FrequencyRange(i),X(i)/1e3);
                    drawnow limitrate
                end
            end
            if  all(X == 0)
                MissingModes(p+N) = 1;
            end
            if  ~MissingModes(p+N)
                X(Misses(1:length(X)) == 1) = NaN;
                SH{p+N}(:,1) = FrequencyRange(1:length(X));
                SH{p+N}(:,2) = FrequencyRange(1:length(X))/1e3;
                SH{p+N}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
                SH{p+N}(:,4) = fillmissing(X,'spline')/1e3;
            else
                if  ModeFamily == 3 && p == 1
                    SH{1}(1,1) = 0;
                else
                    SH{p+N} = SH{p+N-1};
                end
                X1{p+N}(1) = 0;
                continue
            end
            [Max,MaxInd] = max(SH{p+N}(:,4));
            PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
            PhaseVelocityRange2 = PhaseVelocityRange.^2;
            for i = 1:length(PhaseVelocityRange)
                if  Stop
                    return
                end
                for m = 1:SuperLayerSize
                    Alpha(m) = sqrt((Material{m}.Density*PhaseVelocityRange2(i)-c{m}(6,6))/c{m}(4,4));
                    D(m) = Alpha(m)*c{m}(4,4);
                end
                X1{p+N}(i,1) = 0;
                if  i == 1
                    SweepRange = [SH{p+N}(MaxInd,1)-10/PlateThickness/1e3 SH{p+N}(MaxInd,1)+2/PlateThickness/1e3];
                else
                    SweepRange = [X1{p+N}(i-1)-10/PlateThickness/1e3 X1{p+N}(i-1)+2/PlateThickness/1e3];
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
                                    G = Wavenumber*Alpha(m)*LayerThicknesses(m);
                                    CosG = cos(G);
                                    SinG = sin(G);
                                    L{m} = [CosG 1i*SinG./D(m);1i*SinG.*D(m) CosG];
                                end
                                M{1} = L{1};
                                for m = 2:SuperLayerSize
                                    M{1} = M{1}*L{m};
                                end
                                for m = 1:length(Pattern)
                                    M{m+1} = M{m}*M{Pattern(m)};
                                end
                                if  ModeFamily < 3
                                    Y(j,l) = M{end}(2,1);
                                elseif ModeFamily == 3
                                    Y(j,l) = M{end}(2,2);
                                end
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
                            X1{p+N}(i,1) = Frequency(2);
                            break
                        end
                    end
                    if  X1{p+N}(i) > 0
                        break
                    end
                end
                if  X1{p+N}(i) == 0
                    break
                end
                if  Multithreading
                    send(Q2,[X1{p+N}(i),PhaseVelocityRange(i)/1e3,p+N])
                else
                    addpoints(g1(p+N),X1{p+N}(i),PhaseVelocityRange(i)/1e3);
                    drawnow limitrate
                end
            end
            X1{p+N}(:,2) = X1{p+N}/1e3;
            X1{p+N}(:,3) = X1{p+N}(:,1)*PlateThickness;
            X1{p+N}(:,4) = PhaseVelocityRange(1:height(X1{p+N}))/1e3;
            if  X1{p+N}(end,1) == 0
                X1{p+N}(end,1) = H(p);
                X1{p+N}(end,2) = H(p)/1e3;
                X1{p+N}(end,3) = H(p)*PlateThickness;
                X1{p+N}(end,4) = PhaseVelocityLimit/1e3;
            end
            X1{p+N}(X1{p+N}(:,1) == 0,:) = [];
            X1{p+N} = flipud(X1{p+N});
        end
        X1(MissingModes == 1) = [];
        SH(MissingModes == 1) = [];
    end
end
for p = N+1:length(SH)
    SH{p}(SH{p}(:,4) == 0,:) = [];
    SH{p} = vertcat(X1{p},SH{p});
    SH{p}(:,7) = 0;
end