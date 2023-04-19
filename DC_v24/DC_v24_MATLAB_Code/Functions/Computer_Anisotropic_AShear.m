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
function AShear = Computer_Anisotropic_AShear(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth)        
Bisections = 24;

%#ok<*AGROW>
%#ok<*GVMIS> 
global Stop 
Stop = 0;
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
    for p = 1:length(H)
        g(p) = animatedline(ax,'LineStyle','--','color','b');
        g1(p) = animatedline(ax,'LineStyle','--','color','b');
    end
    MissingModes(length(H)) = 0;
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
            if  p > 1
                for j = 1:length(AShear)
                    if  i <= height(AShear{j})
                        Neighbors(j) = AShear{j}(i,4)*1e3;
                    end
                end
            end
            if  all(X == 0)
                if  p == 1
                    SweepRange = [PhaseVelocityLimit 2e3];
                else
                    SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                end
            elseif numel(find(X > 0)) == 1
                if  p == 1
                    if  ~Hybrid
                        SweepRange = [X(i-1) sqrt(c{1}(6,6)/Material{1}.Density)+1];
                    else
                        SweepRange = [X(i-1) 2e3];
                    end
                else
                    SweepRange = [X(i-1) max(Neighbors)+1];
                end
            else
                SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-ShearPhaseVelocitySweepRange*abs((X(i-2)-X(i-1)))];
                if  p == 1
                    if  ~Hybrid && SweepRange(end) < sqrt(c{1}(6,6)/Material{1}.Density)
                        SweepRange(end) = sqrt(c{1}(6,6)/Material{1}.Density);
                    end
                else
                    if  SweepRange(end) < max(Neighbors)+1
                        SweepRange(end) = max(Neighbors)+1;
                    end
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
                            Y(j,l) = M(2,2);
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
                MissingModes(p) = 1;
                break
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end))
                X(end-MissingSamples:end) = [];
                Misses(end-MissingSamples:end) = 0;
                break
            end
            if  X(i) > 0 && BelowCutoff(i) == 0
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
            AShear{p} = AShear{p};
            X1{p}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(AShear{p}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop == 1
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
            addpoints(g1(p),X1{p}(i),PhaseVelocityRange(i)/1e3);
            drawnow limitrate
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