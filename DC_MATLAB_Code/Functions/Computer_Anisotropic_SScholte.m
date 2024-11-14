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
function SScholte = Computer_Anisotropic_SScholte(Multithreading,Q,ax,Fluid,Material,H,FrequencyRange,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*MINV>
%#ok<*GVMIS> 
global Stop 
Stop = 0;
SScholte{1} = [];
if  ~Multithreading
    for p = 1:length(H)+1
        g(p) = animatedline(ax,'LineStyle','-.','color','r');
    end
end
for m = 1:SuperLayerSize
    A1(m) = 2*c{m}(3,3)*c{m}(5,5);
    a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
    a22(m) = -c{m}(3,3)-c{m}(5,5);
    a31(m) = c{m}(1,1)*c{m}(5,5);
    a32(m) = -c{m}(1,1)-c{m}(5,5);
end
X = Fluid.Velocity; % phase velocity (m/ms)
if  Multithreading
    send(Q,[FrequencyRange(1),X(1)/1e3,1])
else
    addpoints(g(1),FrequencyRange(1),X(1)/1e3);
    drawnow limitrate
end
for i = 2:length(FrequencyRange)
    if  Stop
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    kFluid2 = (AngularFrequency/Fluid.Velocity)^2;
    for m = 1:SuperLayerSize
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = rw2(m)^2;
        b22(m) = a22(m)*rw2(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = r2w4(m);
    end
    X(i,1) = 0;
    if  i <= 3
        SweepRange = [1.001*Fluid.Velocity .99*Fluid.Velocity];
    else
        if  abs(X(i-2)-X(i-1)) > abs(X(i-3)-X(i-2))
            Factor = LambPhaseVelocitySweepRange1;
        else
            Factor = LambPhaseVelocitySweepRange2;
        end        
        SweepRange = [X(i-1)+3*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
    end
    if  SweepRange(1) == SweepRange(2)
        if  i < 4 || SweepRange(2) > .99*Fluid.Velocity
            SweepRange = [1.001*Fluid.Velocity .99*Fluid.Velocity];
        else
            SweepRange = [X(i-1)+3*abs(X(i-3)-X(i-2)) X(i-1)-Factor*abs(X(i-3)-X(i-2))];
        end
    end
    if  SweepRange(1) > 1.001*Fluid.Velocity
        SweepRange(1) = 1.001*Fluid.Velocity;
    end
    if  SweepRange(2) < 0
        SweepRange(2) = 0;
    end
    for o = 1:PhaseVelocitySections % increase search resolution
        if  isempty(SweepRange)
            break
        end
        if  1e-6 < abs(SweepRange(1)-SweepRange(2)) % how many bisections have to performed to reach the desired PhaseVelocityResolution
            Bisections = ceil(log2(1e-6/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25^o)); % 2* because each new SweepRange spans two increments of the preceding SweepRange
        else
            Bisections = 1;
        end
        for k = 1:Bisections % search minimum in characteristic equation and converge upon it
            if  k == 1 % divide SweepRange into sections where the characteristic equation is evaluated
                SweepRange = SweepRange(1):.25^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
            else
                SweepRange = SweepRange(1):.25*(SweepRange(end)-SweepRange(1)):SweepRange(end);
            end
            MIN = 0;
            Y = 0;
            Wavenumber = AngularFrequency./SweepRange;
            Wavenumber2 = Wavenumber.^2;
            Wavenumber4 = Wavenumber2.^2;
            for j = 1:length(Wavenumber)
                for m = 1:SuperLayerSize
                    A2 = a21(m)*Wavenumber2(j)+b22(m);
                    A3 = a31(m)*Wavenumber4(j)+b32(m)*Wavenumber2(j)+b33(m);
                    k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    W = (rw2(m)-c{m}(1,1)*Wavenumber2(j)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j)*k3);
                    D3 = 1i*(c{m}(1,3)*Wavenumber(j)+c{m}(3,3)*k3.*W); % sigma33
                    D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j)*W)); % sigma13
                    if  SuperLayerSize == 1
                        E = exp(.5i*k3*LayerThicknesses);
                    else
                        E = exp(1i*k3*LayerThicknesses(m));
                    end
                    L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                    L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                    L{m} = L2/L1;
                end
                M = L{1};
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                    M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                end
                MM{1} = M;
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                end
                G = MM{end};
                k3Fluid = sqrt(kFluid2-Wavenumber2(j));
                MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency^2);
                Y(j) = abs(MFluid+G(2,1)-G(2,3)*G(4,1)/G(4,3));
            end
            Min = zeros(1,length(Y));
            for j = 2:length(Y)-1
                if  Y(j) < Y(j-1) && Y(j) < Y(j+1)
                    Min(j) = 1;
                end
            end
            b = find(Min);
            if  ~isempty(b) % one or multiple minima are found
                if  isscalar(b)
                    MIN = b; % SweepRange-index
                else
                    cp = SweepRange(b);
                    delta = abs(cp-X(i-1));
                    [~,j] = min(delta);
                    MIN = b(j); % SweepRange-index
                end
            end
            if  MIN == 0
                break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
            end
            if  k < Bisections % set the new search area around the found minimum
                SweepRange = [SweepRange(MIN)-(SweepRange(2)-SweepRange(1)) SweepRange(MIN)+(SweepRange(2)-SweepRange(1))];
            end
        end
        if  MIN > 0
            z = isoutlier(vertcat(X(1:end-1),SweepRange(MIN)),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
            if  (i < 4 && z(end)) || (i >= 4 && z(end) && abs(X(i-1)-SweepRange(MIN)) > 1)
                Outlier = 1;
            else
                Outlier = 0;
            end
            if  ~Outlier || all(X == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                X(i,1) = SweepRange(MIN); % phase velocity (m/s)
                Misses(i) = 0;
                break
            end
        end
    end
    if  i == 2 && X(2) == 0
        X(1:length(FrequencyRange),1) = 0;
        break
    elseif i > 2 && X(i) == 0 % fit phase velocity where we missed the solution to obtain useful sweep range for the next frequency step
        Fit = fit(FrequencyRange(1:i-1)',X(1:i-1),'cubicspline');
        X(i,1) = Fit(FrequencyRange(i));
        Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
    end
    if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops
        X(end-MissingSamples:end) = [];
        Misses(end-MissingSamples:end) = 0;
        break
    end
    if  Multithreading
        send(Q,[FrequencyRange(i),X(i)/1e3,1])
    else
        addpoints(g(1),FrequencyRange(i),X(i)/1e3);
        drawnow limitrate
    end
end
if  any(X)
    X(Misses(1:length(X)) == 1) = NaN;
end
SScholte{1}(:,1) = FrequencyRange(1:length(X));
SScholte{1}(:,2) = FrequencyRange(1:length(X))/1e3;
SScholte{1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
SScholte{1}(:,4) = fillmissing(X,'spline')/1e3;
SScholte{1}(:,6) = 0;
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
            kFluid2 = (AngularFrequency/Fluid.Velocity)^2;
            for m = 1:SuperLayerSize
                rw2(m) = Material{m}.Density*AngularFrequency^2;
                r2w4(m) = rw2(m)^2;
                b22(m) = a22(m)*rw2(m);
                b32(m) = a32(m)*rw2(m);
                b33(m) = r2w4(m);
            end
            X(i,1) = 0;
            Neighbors = [];
            for j = 1:length(SScholte)
                if  i <= height(SScholte{j})
                    Neighbors(j) = SScholte{j}(i,4)*1e3;
                end
            end
            if  numel(find(X > 0)) <= 5
                SweepRange = [1.001*Fluid.Velocity .99*Fluid.Velocity];
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
            if  SweepRange(1) > 1.001*Fluid.Velocity
                SweepRange(1) = 1.001*Fluid.Velocity;
            end
            if  SweepRange(2) < max(Neighbors)+1
                SweepRange(2) = max(Neighbors)+1;
            end
            for o = 1:PhaseVelocitySections % increase search resolution
                if  isempty(SweepRange)
                    break
                end
                if  PhaseVelocityResolution < abs(SweepRange(1)-SweepRange(2)) % how many bisections have to performed to reach the desired PhaseVelocityResolution
                    Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25^o)); % 2* because each new SweepRange spans two increments of the preceding SweepRange
                else
                    Bisections = 1;
                end
                for k = 1:Bisections % search minimum in characteristic equation and converge upon it
                    if  k == 1 % divide SweepRange into sections where the characteristic equation is evaluated
                        SweepRange = SweepRange(1):.25^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    else
                        SweepRange = SweepRange(1):.25*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                    MIN = 0;
                    Y = 0;
                    Wavenumber = AngularFrequency./SweepRange;
                    Wavenumber2 = Wavenumber.^2;
                    Wavenumber4 = Wavenumber2.^2;
                    for j = 1:length(Wavenumber)
                        for m = 1:SuperLayerSize
                            A2 = a21(m)*Wavenumber2(j)+b22(m);
                            A3 = a31(m)*Wavenumber4(j)+b32(m)*Wavenumber2(j)+b33(m);
                            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                            W = (rw2(m)-c{m}(1,1)*Wavenumber2(j)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j)*k3);
                            D3 = 1i*(c{m}(1,3)*Wavenumber(j)+c{m}(3,3)*k3.*W); % sigma33
                            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j)*W)); % sigma13
                            if  SuperLayerSize == 1
                                E = exp(.5i*k3*LayerThicknesses);
                            else
                                E = exp(1i*k3*LayerThicknesses(m));
                            end
                            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                            L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                            L{m} = L2/L1;
                        end
                        M = L{1};
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                            M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                        end
                        MM{1} = M;
                        for m = 2:log2(Repetitions)+1
                            N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                        end
                        for m = m+1:length(Pattern)
                            N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                        end
                        G = MM{end};
                        k3Fluid = sqrt(kFluid2-Wavenumber2(j));
                        MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency^2);
                        Y(j) = abs(MFluid+G(2,1)-G(2,3)*G(4,1)/G(4,3));
                    end
                    for j = 1:length(Y)-1
                        for n = 1:length(Neighbors)
                            if  SweepRange(j) > Neighbors(n) && SweepRange(j+1) < Neighbors(n)
                                if  j == 1
                                    Y(2) = NaN;
                                else
                                    Y(j) = NaN;
                                end
                            end
                        end
                    end
                    Min = zeros(1,length(Y));
                    for j = 2:length(Y)-1
                        if  Y(j) < Y(j-1) && Y(j) < Y(j+1)
                            Min(j) = 1;
                        end
                    end
                    b = find(Min);
                    if  ~isempty(b) % one or multiple minima are found
                        if  isscalar(b)
                            MIN = b; % SweepRange-index
                        else
                            cp = SweepRange(b);
                            delta = abs(cp-X(i-1));
                            [~,j] = min(delta);
                            MIN = b(j); % SweepRange-index
                        end
                    end
                    if  MIN == 0
                        break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                    end
                    if  k < Bisections % set the new search area around the found minimum
                        SweepRange = [SweepRange(MIN)-(SweepRange(2)-SweepRange(1)) SweepRange(MIN)+(SweepRange(2)-SweepRange(1))];
                    end
                end
                if  MIN > 0
                    if  numel(find(X > 0)) <= 5
                        Outlier = 0;
                    else
                        z = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:i-1),SweepRange(MIN)),'movmedian',5,'ThresholdFactor',9);
                        if  z(end) && abs(X(i-1)-SweepRange(MIN)) > 1
                            Outlier = 1;
                        else
                            Outlier = 0;
                        end
                    end
                    if  ~Outlier || all(X == 0)
                        X(i,1) = SweepRange(MIN);
                        Misses(i) = 0;
                        BelowCutoff(i) = 0;
                        break
                    end
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
            if  X(i) > 0 && BelowCutoff(i) == 0
                if  Multithreading
                    send(Q,[FrequencyRange(i),X(i)/1e3,p+1])
                else
                    addpoints(g(p+1),FrequencyRange(i),X(i)/1e3);
                    drawnow limitrate
                end
            end
        end
        if  all(X == 0)
            MissingModes(p+1) = 1;
        end
        if  ~MissingModes(p+1)
            if  any(X)
                X(Misses(1:length(X)) == 1) = NaN;
            end
            SScholte{p+1}(:,1) = FrequencyRange(1:length(X));
            SScholte{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
            SScholte{p+1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
            SScholte{p+1}(:,4) = fillmissing(X,'spline')/1e3;
        else
            SScholte{p+1} = SScholte{p};
        end
    end
    SScholte(MissingModes == 1) = [];
    for p = 2:length(SScholte)
        SScholte{p}(SScholte{p}(:,4) == 0,:) = [];
        SScholte{p}(:,6) = 0;
    end    
end