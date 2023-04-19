% =========================================================================
% Dispersion Calculator
% Copyright (C) 2018-2023 DLR
% Created by Armin Huber
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function SShear = Computer_Anisotropic_SShear(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth)        
Bisections = 24;

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop  
Stop = 0;
if  SuperLayerSize == 1
    SShear{1}(:,1) = FrequencyRange;
    SShear{1}(:,2) = FrequencyRange/1e3;
    SShear{1}(:,3) = FrequencyRange*PlateThickness;
    SShear{1}(:,4) = sqrt(c{1}(6,6)/Material{1}.Density)/1e3;
    SShear{1}(:,6) = 0;
    if  HigherOrderModes && any(H)
        for p = 1:length(H)
            i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange);
            if  isempty(i)
                continue
            end
            FrqRange = FrequencyRange(i)*1e3;
            SShear{p+1}(i,4) = sqrt(c{1}(6,6)./(Material{1}.Density-c{1}(4,4)*(p./(FrqRange*PlateThickness)).^2))/1e3;
            SShear{p+1}(:,1) = FrequencyRange(1:height(SShear{p+1}));
            SShear{p+1}(:,2) = FrequencyRange(1:height(SShear{p+1}))/1e3;
            SShear{p+1}(:,3) = FrequencyRange(1:height(SShear{p+1}))*PlateThickness;
            PhaseVelocityRange = max(SShear{p+1}(:,4))*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
            X1{p+1}(:,1) = p*PhaseVelocityRange.*sqrt(c{1}(4,4)./(Material{1}.Density*PhaseVelocityRange.^2-c{1}(6,6)))/PlateThickness/1e3;
            X1{p+1}(:,2) = X1{p+1}/1e3;
            X1{p+1}(:,3) = X1{p+1}(:,1)*PlateThickness;
            X1{p+1}(:,4) = PhaseVelocityRange(1:height(X1{p+1}))/1e3;
            X1{p+1}(X1{p+1}(:,1) == 0,:) = [];
            X1{p+1} = flipud(X1{p+1});
        end
    end
else
    if  HigherOrderModes && any(H)
        for p = 1:length(H)+1
            g(p) = animatedline(ax,'LineStyle','--','color','r');
            g1(p) = animatedline(ax,'LineStyle','--','color','r');
        end
    else
        g = animatedline(ax,'LineStyle','--','color','r');
    end
    if  ~Hybrid
        SShear{1}(:,1) = FrequencyRange;
        SShear{1}(:,2) = FrequencyRange/1e3;
        SShear{1}(:,3) = FrequencyRange*PlateThickness;
        SShear{1}(:,4) = sqrt(c{1}(6,6)/Material{1}.Density)/1e3;
        SShear{1}(:,6) = 0;
        line(ax,SShear{1}(:,1),SShear{1}(:,4),'LineStyle','--','color','r')
    else
        for i = 1:length(FrequencyRange)
            if  Stop == 1
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
                            WaveNumber = AngularFrequency/PhaseVelocity(l);
                            for m = 1:SuperLayerSize
                                Alpha = sqrt((Material{m}.Density*PhaseVelocity(l)^2-c{m}(6,6))/c{m}(4,4));
                                D = Alpha*c{m}(4,4);
                                Gamma = Alpha*WaveNumber*LayerThicknesses(m);
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
            addpoints(g(1),FrequencyRange(i),X(i)/1e3);
            drawnow limitrate
        end
        X(Misses(1:length(X)) == 1) = NaN;
        SShear{1}(:,1) = FrequencyRange;
        SShear{1}(:,2) = FrequencyRange/1e3;
        SShear{1}(:,3) = FrequencyRange*PlateThickness;
        SShear{1}(:,4) = fillmissing(X,'spline')/1e3;    
        SShear{1}(:,6) = 0;
    end
    if  HigherOrderModes && any(H)
        MissingModes(length(H)+1) = 0;
        for p = 1:length(H)
            X = [];
            Misses = 0;
            BelowCutoff = 0;
            for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
                if  Stop == 1
                    return
                end
                AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
                X(i,1) = 0;
                Neighbors = [];
                for j = 1:length(SShear)
                    if  i <= height(SShear{j})
                        Neighbors(j) = SShear{j}(i,4)*1e3;
                    end
                end
                if  all(X == 0)
                    SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                elseif numel(find(X > 0)) == 1
                    SweepRange = [X(i-1) max(Neighbors)+1];
                else
                    SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-ShearPhaseVelocitySweepRange*abs(X(i-2)-X(i-1))];
                    if  SweepRange(end) < max(Neighbors)+1
                        SweepRange(end) = max(Neighbors)+1;
                    end
                end
                for o = 0:PhaseVelocitySections+1
                    if  numel(find(X > 0)) <= 1 && o > 0
                        if  o < PhaseVelocitySections+1
                            SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
                        elseif o == PhaseVelocitySections+1
                            if  all(X == 0)
                                SweepRange = PhaseVelocityLimit:-.2^(o+1)*(SweepRange(1)-SweepRange(end)):PhaseVelocityLimit-.2*(SweepRange(1)-SweepRange(end));
                            elseif numel(find(X > 0)) == 1
                                SweepRange = X(i-1):-.2^(o+1)*(SweepRange(1)-SweepRange(end)):X(i-1)-.2*(SweepRange(1)-SweepRange(end));
                            end
                        end
                    elseif numel(find(X > 0)) > 1 && o > 0
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
                                WaveNumber = AngularFrequency/PhaseVelocity(l);
                                for m = 1:SuperLayerSize
                                    Alpha = sqrt((Material{m}.Density*PhaseVelocity(l)^2-c{m}(6,6))/c{m}(4,4));
                                    D = Alpha*c{m}(4,4);
                                    Gamma = Alpha*WaveNumber*LayerThicknesses(m);
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
                            BelowCutoff(i) = 0;
                            break
                        end
                    end
                    if  X(i) > 0 || (abs(SweepRange(1)-SweepRange(end)) < 1000 && o == PhaseVelocitySections)
                        break
                    end
                end
                if  X(i) == 0 && any(X)
                    if  numel(find(X > 0)) == 1
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
                if  X(i) > 0 && BelowCutoff(i) == 0
                    addpoints(g(p+1),FrequencyRange(i),X(i)/1e3);
                    drawnow limitrate
                end
            end
            if  all(X == 0)
                MissingModes(p+1) = 1;
            end
            if  ~MissingModes(p+1)
                X(Misses(1:length(X)) == 1) = NaN;
                SShear{p+1}(:,1) = FrequencyRange(1:length(X));
                SShear{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
                SShear{p+1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
                SShear{p+1}(:,4) = fillmissing(X,'spline')/1e3;
            else
                SShear{p+1} = SShear{p};
                X1{p+1}(1) = 0;
                continue
            end
            [Max,MaxInd] = max(SShear{p+1}(:,4));
            PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
            for i = 1:length(PhaseVelocityRange)
                if  Stop == 1
                    return
                end
                X1{p+1}(i,1) = 0;
                if  i == 1
                    SweepRange = [SShear{p+1}(MaxInd,1)-10/PlateThickness/1e3 SShear{p+1}(MaxInd,1)+2/PlateThickness/1e3];
                else
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
                                WaveNumber = 2*pi*Frequency(l)*1e3/PhaseVelocityRange(i);
                                for m = 1:SuperLayerSize
                                    Alpha = sqrt((Material{m}.Density*PhaseVelocityRange(i)^2-c{m}(6,6))/c{m}(4,4));
                                    D = Alpha*c{m}(4,4);
                                    Gamma = Alpha*WaveNumber*LayerThicknesses(m);
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
                addpoints(g1(p+1),X1{p+1}(i),PhaseVelocityRange(i)/1e3);
                drawnow limitrate
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
        SShear(MissingModes == 1) = [];
    end
end
for p = 2:length(SShear)
    SShear{p}(SShear{p}(:,4) == 0,:) = [];
    SShear{p} = vertcat(X1{p},SShear{p});
    SShear{p}(:,6) = 0;
end