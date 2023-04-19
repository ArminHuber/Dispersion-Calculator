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
function ALamb = Computer_Anisotropic_ALamb(Multithreading,Q1,Q2,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*MINV>
%#ok<*GVMIS>
global Stop
Stop = 0;
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
for m = 1:SuperLayerSize
    A1(m) = 2*c{m}(3,3)*c{m}(5,5);
    a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
    a22(m) = -c{m}(3,3)-c{m}(5,5);
    a31(m) = c{m}(1,1)*c{m}(5,5);
    a32(m) = -c{m}(1,1)-c{m}(5,5);
end
for i = 1:length(FrequencyRange)
    if  Stop == 1
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    X(i,1) = 0;
    if  i > 2
        if  X(i-2) == X(i-1)
            x = abs(X-X(i-1));
            x(x == 0) = [];
            delta = min(x);
        else
            delta = abs(X(i-2)-X(i-1));
        end
    end
    if  i == 1
        SweepRange = .1:50;
    elseif i == 2
        SweepRange = .1:1e3;
    else
        SweepRange = [X(i-1)-2*LambPhaseVelocitySweepRange2*delta X(i-1)+2*LambPhaseVelocitySweepRange2*delta];
    end
    if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
        if  delta > PhaseVelocityResolution
            SweepRange(1) = SweepRange(1)-PhaseVelocityResolution;
            SweepRange(2) = SweepRange(2)+PhaseVelocityResolution;
        else
            SweepRange(1) = SweepRange(1)-2*delta;
            SweepRange(2) = SweepRange(2)+2*delta;
        end
    end
    if  SweepRange(1) < 0
        SweepRange(1) = 0;
    end
    for o = 0:PhaseVelocitySections
        if  o > 0
            SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
        end
        if  i > 2 && delta < PhaseVelocityResolution
            Bisections = ceil(log2(delta/PhaseVelocityResolution)/log2(.5));
        else
            Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
        end
        if  i < 6 && Bisections < 20
            Bisections = 20;
        elseif i >= 6
            if  o == 0 && Bisections < 5
                Bisections = 5;
            elseif o > 0 && Bisections < 1
                Bisections = 1;
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
                    WaveNumber = AngularFrequency/PhaseVelocity(l);
                    for m = 1:SuperLayerSize
                        rc2 = Material{m}.Density*PhaseVelocity(l)^2;
                        A2 = a21(m)+a22(m)*rc2;
                        A3 = a31(m)+a32(m)*rc2+Material{m}.Density^2*PhaseVelocity(l)^4;
                        Alpha(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                        Alpha(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                        W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha.^2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
                        D3 = c{m}(1,3)+c{m}(3,3)*Alpha.*W;
                        D5 = c{m}(5,5)*(Alpha+W);
                        if  SuperLayerSize == 1
                            E = exp(.5i*WaveNumber*Alpha*LayerThicknesses);
                        else
                            E = exp(1i*WaveNumber*Alpha*LayerThicknesses(m));
                        end
                        L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                        L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                        L{m} = L1/L2;
                    end
                    M = L{1};
                    for m = 2:SuperLayerSize
                        N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                        M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                    end
                    MM{1} = M;
                    for m = 2:Repetitions
                        N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                        MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
                    end
                    Y(j,l) = det(MM{end}(1:3,[1:2,4]));
                end
                if  k == 1
                    Y(j+1,1) = Y(j,3);
                end
                if  (j == 1 && abs(imag(Y(1,2))) < abs(imag(Y(1,1))) && abs(imag(Y(1,2))) < abs(imag(Y(1,3)))) | (j > 1 && abs(imag(Y(j,2))) < abs(imag(Y(j-1,2)))) && sign(imag(Y(j,1))) ~= sign(imag(Y(j,2)))
                    PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                    Y(j,3) = Y(j,2);
                elseif (j == 1 && abs(imag(Y(1,2))) < abs(imag(Y(1,1))) && abs(imag(Y(1,2))) < abs(imag(Y(1,3)))) | (j > 1 && abs(imag(Y(j,2))) < abs(imag(Y(j-1,2)))) && sign(imag(Y(j,2))) ~= sign(imag(Y(j,3)))
                    PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                    Y(j,1) = Y(j,2);
                else
                    PhaseVelocity(2) = 0;
                    break
                end
            end
            if  PhaseVelocity(2) > 0
                if  i < 4
                    Outlier = 0;
                else
                    z = isoutlier(vertcat(X(1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
                    if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
                        Outlier = 1;
                    else
                        Outlier = 0;
                    end
                end
                if  ~Outlier || all(X == 0)
                    X(i,1) = PhaseVelocity(2);
                    Misses(i) = 0;
                    break
                end
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
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i) == 1
%     disp(['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
% end
end
X(Misses(1:length(X)) == 1) = NaN;
ALamb{1}(:,1) = FrequencyRange(1:length(X));
ALamb{1}(:,2) = FrequencyRange(1:length(X))/1e3;
ALamb{1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
ALamb{1}(:,4) = fillmissing(X,'spline')/1e3;
ALamb{1}(:,6) = 0;
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
            for j = 1:length(ALamb)
                if  i <= height(ALamb{j})
                    Neighbors(j) = ALamb{j}(i,4)*1e3;
                end
            end
            if  all(X == 0)
                SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
            elseif numel(find(X > 0)) == 1
                SweepRange = [X(i-1) max(Neighbors)+1];
            else
                if  abs(X(i-2)-X(i-1)) > abs(X(i-3)-X(i-2))
                    Factor = LambPhaseVelocitySweepRange1;
                else
                    Factor = LambPhaseVelocitySweepRange2;
                end
                if  ~Hybrid && (strcmp(Material{1}.Class,'Transversely isotropic') || strcmp(Material{1}.Class,'Isotropic'))
                    SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
                elseif (~Hybrid && (strcmp(Material{1}.Class,'Orthotropic') || strcmp(Material{1}.Class,'Cubic'))) || Hybrid
                    if  4*abs(X(i-2)-X(i-1)) > 50
                        Top = 50;
                    else
                        Top = 4*abs(X(i-2)-X(i-1));
                    end
                    SweepRange = [X(i-1)+Top X(i-1)-Factor*abs(X(i-2)-X(i-1))];
                end
            end
            if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
                SweepRange(1) = SweepRange(1)+PhaseVelocityResolution;
                SweepRange(2) = SweepRange(2)-PhaseVelocityResolution;
            end
            if  SweepRange(end) < max(Neighbors)+1
                SweepRange(end) = max(Neighbors)+1;
            end
            for o = 0:PhaseVelocitySections+1
                if  numel(find(X > 0)) <= 1 && o > 0
                    if  o < PhaseVelocitySections+1
                        SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
                    else
                        if  all(X == 0)
                            SweepRange = PhaseVelocityLimit:-.2^(o+1)*(SweepRange(1)-SweepRange(end)):PhaseVelocityLimit-.2*(SweepRange(1)-SweepRange(end));
                        elseif numel(find(X > 0)) == 1
                            SweepRange = X(i-1):-.2^(o+1)*(SweepRange(1)-SweepRange(end)):X(i-1)-.2*(SweepRange(1)-SweepRange(end));
                        end
                    end
                elseif numel(find(X > 0)) > 1 && o > 0
                    SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
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
                                rc2 = Material{m}.Density*PhaseVelocity(l)^2;
                                A2 = a21(m)+a22(m)*rc2;
                                A3 = a31(m)+a32(m)*rc2+Material{m}.Density^2*PhaseVelocity(l)^4;
                                Alpha(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                                Alpha(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                                W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha.^2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
                                D3 = c{m}(1,3)+c{m}(3,3)*Alpha.*W;
                                D5 = c{m}(5,5)*(Alpha+W);
                                if  SuperLayerSize == 1
                                    E = exp(.5i*WaveNumber*Alpha*LayerThicknesses);
                                else
                                    E = exp(1i*WaveNumber*Alpha*LayerThicknesses(m));
                                end
                                L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                                L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                                L{m} = L1/L2;
                            end
                            M = L{1};
                            for m = 2:SuperLayerSize
                                N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                                M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                            end
                            MM{1} = M;
                            for m = 2:Repetitions
                                N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                                MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
                            end
                            Y(j,l) = det(MM{end}(1:3,[1:2,4]));
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  (j == 1 && abs(imag(Y(1,2))) < abs(imag(Y(1,1))) && abs(imag(Y(1,2))) < abs(imag(Y(1,3)))) | (j > 1 && abs(imag(Y(j,2))) < abs(imag(Y(j-1,2)))) && sign(imag(Y(j,1))) ~= sign(imag(Y(j,2)))
                            PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif (j == 1 && abs(imag(Y(1,2))) < abs(imag(Y(1,1))) && abs(imag(Y(1,2))) < abs(imag(Y(1,3)))) | (j > 1 && abs(imag(Y(j,2))) < abs(imag(Y(j-1,2)))) && sign(imag(Y(j,2))) ~= sign(imag(Y(j,3)))
                            PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                            Y(j,1) = Y(j,2);
                        else
                            PhaseVelocity(2) = 0;
                            break
                        end
                    end
                    if  PhaseVelocity(2) > 0
                        if  numel(find(X > 0)) <= 5
                            Outlier = 0;
                        else
                            z = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
                            if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
                                Outlier = 1;
                            else
                                Outlier = 0;
                            end
                        end
                        if  ~Outlier || all(X == 0)
                            X(i,1) = PhaseVelocity(2);
                            Misses(i) = 0;
                            BelowCutoff(i) = 0;
                            break
                        end
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
            if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                send(Q1,[FrequencyRange(i),X(i)/1e3,p+1])
            elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                addpoints(g(p+1),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i) == 1
%     disp(['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
% end
        end
        if  all(X == 0)
            MissingModes(p+1) = 1;
        end
        if  ~MissingModes(p+1)
            X(Misses(1:length(X)) == 1) = NaN;
            ALamb{p+1}(:,1) = FrequencyRange(1:length(X));
            ALamb{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
            ALamb{p+1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
            ALamb{p+1}(:,4) = fillmissing(X,'spline')/1e3;
        else
            ALamb{p+1} = ALamb{p};
            X1{p+1}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(ALamb{p+1}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop == 1
                return
            end
            for m = 1:SuperLayerSize
                rc2 = Material{m}.Density*PhaseVelocityRange(i)^2;
                A2 = a21(m)+a22(m)*rc2;
                A3 = a31(m)+a32(m)*rc2+Material{m}.Density^2*PhaseVelocityRange(i)^4;
                Alpha(m,1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                Alpha(m,2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                W(m,:) = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha(m,:).^2)./((c{m}(1,3)+c{m}(5,5))*Alpha(m,:));
                D3(m,:) = c{m}(1,3)+c{m}(3,3)*Alpha(m,:).*W(m,:);
                D5(m,:) = c{m}(5,5)*(Alpha(m,:)+W(m,:));
            end
            X1{p+1}(i,1) = 0;
            if  i == 1
                if  p == 1
                    SweepRange = [ALamb{p+1}(MaxInd,1)-2*FrequencyOffset/PlateThickness/1e3 ALamb{p+1}(MaxInd,1)+2/PlateThickness/1e3];
                else
                    SweepRange = [ALamb{p+1}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 ALamb{p+1}(MaxInd,1)+2/PlateThickness/1e3];
                end
            else
                if  p == 1
                    SweepRange = [X1{p+1}(i-1)-2*FrequencyOffset/PlateThickness/1e3 X1{p+1}(i-1)+2/PlateThickness/1e3];
                else
                    SweepRange = [X1{p+1}(i-1)-FrequencyOffset/PlateThickness/1e3 X1{p+1}(i-1)+2/PlateThickness/1e3];
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
                            WaveNumber = 2*pi*Frequency(l)*1e3/PhaseVelocityRange(i);
                            for m = 1:SuperLayerSize
                                if  SuperLayerSize == 1
                                    E = exp(.5i*WaveNumber*Alpha*LayerThicknesses);
                                else
                                    E = exp(1i*WaveNumber*Alpha(m,:)*LayerThicknesses(m));
                                end
                                L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
                                L2 = [1 1 E;W(m,:) -W(m,:).*E;E 1 1;W(m,:).*E -W(m,:)];
                                L{m} = L1/L2;
                            end
                            M = L{1};
                            for m = 2:SuperLayerSize
                                N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                                M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                            end
                            MM{1} = M;
                            for m = 2:Repetitions
                                N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                                MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
                            end
                            Y(j,l) = det(MM{end}(1:3,[1:2,4]));
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  (j == 1 && abs(imag(Y(1,2))) < abs(imag(Y(1,1))) && abs(imag(Y(1,2))) < abs(imag(Y(1,3)))) | (j > 1 && abs(imag(Y(j,2))) < abs(imag(Y(j-1,2)))) && sign(imag(Y(j,1))) ~= sign(imag(Y(j,2)))
                            Frequency = [Frequency(1) Frequency(1)+(Frequency(2)-Frequency(1))/2 Frequency(2)];
                            Y(j,3) = Y(j,2);
                        elseif (j == 1 && abs(imag(Y(1,2))) < abs(imag(Y(1,1))) && abs(imag(Y(1,2))) < abs(imag(Y(1,3)))) | (j > 1 && abs(imag(Y(j,2))) < abs(imag(Y(j-1,2)))) && sign(imag(Y(j,2))) ~= sign(imag(Y(j,3)))
                            Frequency = [Frequency(2) Frequency(2)+(Frequency(3)-Frequency(2))/2 Frequency(3)];
                            Y(j,1) = Y(j,2);
                        else
                            Frequency(2) = 0;
                            break
                        end
                    end
                    if  Frequency(2) > 0
                        if  i < 4
                            Outlier = 0;
                        else
                            z = isoutlier(vertcat(X1{p+1}(1:i-1),Frequency(2)),'movmedian',5,'ThresholdFactor',9);
                            if  z(end) && abs(X1{p+1}(i-1)-Frequency(2)) > 1
                                Outlier = 1;
                            else
                                Outlier = 0;
                            end
                        end
                        if  ~Outlier
                            X1{p+1}(i,1) = Frequency(2);
                            break
                        end
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
        Alpha = 0;
        W = 0;
        D3 = 0;
        D5 = 0;
    end
    X1(MissingModes == 1) = [];
    ALamb(MissingModes == 1) = [];
    for p = 2:length(ALamb)
        ALamb{p}(ALamb{p}(:,4) == 0,:) = [];
        ALamb{p} = vertcat(X1{p},ALamb{p});
        ALamb{p}(:,6) = 0;
    end
end