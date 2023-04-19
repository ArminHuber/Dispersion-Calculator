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
function BScholte = Computer_Anisotropic_BScholte_Coupled(ax,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FScholte,H,FrequencyRange,PhaseVelocitySections,PlateThickness,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,SymmetricSystem,I,Delta,MissingSamples)
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*MINV>
%#ok<*GVMIS>
global Stop
Stop = 0;
if  ToggleUpperFluid && ToggleLowerFluid
    if  strcmp(UpperFluid.Name,LowerFluid.Name)
        FluidVelocity = [UpperFluid.Velocity 0];
        SlowFluidVelocity = UpperFluid.Velocity;
        FastFluidVelocity = UpperFluid.Velocity;
    else
        if  UpperFluid.Velocity > LowerFluid.Velocity
            FluidVelocity = [LowerFluid.Velocity UpperFluid.Velocity]; % slower and faster velocity
        else
            FluidVelocity = [UpperFluid.Velocity LowerFluid.Velocity];
        end
        SlowFluidVelocity = FluidVelocity(1);
        FastFluidVelocity = FluidVelocity(2);
    end
elseif ToggleUpperFluid && ~ToggleLowerFluid
    FluidVelocity = [UpperFluid.Velocity 0];
    SlowFluidVelocity = UpperFluid.Velocity;
    FastFluidVelocity = UpperFluid.Velocity;
elseif ~ToggleUpperFluid && ToggleLowerFluid
    FluidVelocity = [LowerFluid.Velocity 0];
    SlowFluidVelocity = LowerFluid.Velocity;
    FastFluidVelocity = LowerFluid.Velocity;
end
FModes = numel(find(FluidVelocity > 0));
for p = 1:length(H)+1+FModes
    g(p) = animatedline(ax,'LineStyle','-.','color',[.5 0 1]);
end
for m = 1:SuperLayerSize
    a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
    a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
    a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
    a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
    a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
    a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
    a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
    a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
    a34(m) = -1/Delta(m);
end
X = FScholte(1);
addpoints(g(1),FrequencyRange(1),X/1e3);
drawnow limitrate
for i = 2:length(FrequencyRange)
    if  Stop == 1
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    for m = 1:SuperLayerSize
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = Material{m}.Density^2*AngularFrequency^4;
        r3w6(m) = Material{m}.Density^3*AngularFrequency^6;
        b12(m) = a12(m)*rw2(m);
        b22(m) = a22(m)*rw2(m);
        b23(m) = a23(m)*r2w4(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = a33(m)*r2w4(m);
        b34(m) = a34(m)*r3w6(m);
    end
    X(i,1) = 0;
    if  i <= 3
        SweepRange = [6*X(i-1) X(i-1)];
    elseif i > 3
        if  abs((X(i-2)-X(i-1))) > abs((X(i-3)-X(i-2)))
            Factor = LambPhaseVelocitySweepRange1;
        else
            Factor = LambPhaseVelocitySweepRange2;
        end
        SweepRange = [X(i-1)+3*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs((X(i-2)-X(i-1)))];        
    end
    if  SweepRange(1) == SweepRange(2) && i > 3
        SweepRange = [X(i-1)+3*abs(X(i-3)-X(i-2)) X(i-1)-Factor*abs((X(i-3)-X(i-2)))];
    end
    if  i <= 10 && SweepRange(1) > SlowFluidVelocity
        SweepRange(1) = SlowFluidVelocity;
    elseif i > 10 && SweepRange(1) > 1.001*FastFluidVelocity
        SweepRange(1) = 1.001*FastFluidVelocity;
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
            WaveNumber = AngularFrequency./SweepRange;
            WaveNumber2 = WaveNumber.^2;
            WaveNumber4 = WaveNumber.^4;
            WaveNumber6 = WaveNumber.^6;
            for j = 1:length(WaveNumber)
                for m = 1:SuperLayerSize
                    A1 = a11(m)*WaveNumber2(j)+b12(m);
                    A2 = a21(m)*WaveNumber4(j)+b22(m)*WaveNumber2(j)+b23(m);
                    A3 = a31(m)*WaveNumber6(j)+b32(m)*WaveNumber4(j)+b33(m)*WaveNumber2(j)+b34(m);
                    k3a = A2/3-A1^2/9;
                    k3b = A1^3/27-A1*A2/6+A3/2;
                    k3c = (sqrt(k3b^2+k3a^3)-k3b)^(1/3);
                    k3d = k3a/(2*k3c)-k3c/2;
                    k3e = k3a/k3c;
                    k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
                    k3(1) = sqrt(k3d-k3f-A1/3);
                    k3(2) = sqrt(k3d+k3f-A1/3);
                    k3(3) = -sqrt(k3c-k3e-A1/3);
                    k32 = k3.^2;
                    m11 = c{m}(1,1)*WaveNumber2(j)-rw2(m)+c{m}(5,5)*k32;
                    m12 = c{m}(1,6)*WaveNumber2(j)+c{m}(4,5)*k32;
                    m13 = (c{m}(1,3)+c{m}(5,5))*WaveNumber(j)*k3;
                    m22 = c{m}(6,6)*WaveNumber2(j)-rw2(m)+c{m}(4,4)*k32;
                    m23 = (c{m}(3,6)+c{m}(4,5))*WaveNumber(j)*k3;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D3 = 1i*(c{m}(1,3)*WaveNumber(j)+c{m}(3,6)*WaveNumber(j)*V+c{m}(3,3)*k3.*W); % sigma33
                    D4 = 1i*(c{m}(4,5)*(k3+WaveNumber(j)*W)+c{m}(4,4)*k3.*V); % sigma23
                    D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(j)*W)+c{m}(4,5)*k3.*V); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                    L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                    L{m} = L1/L2;
                end
                M = L{1};
                for m = 2:SuperLayerSize
                    N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                    M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
                end
                MM{1} = M;
                for m = 2:Repetitions
                    N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                    MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                end
                G = inv(MM{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j));
                    k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j));
                    WUpperFluid = k3UpperFluid/WaveNumber(j);
                    WLowerFluid = k3LowerFluid/WaveNumber(j);
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j); % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j); % in the lower fluid
                    Y(j) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j));
                    WUpperFluid = k3UpperFluid/WaveNumber(j);
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j);
                    Y(j) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j));
                    WLowerFluid = k3LowerFluid/WaveNumber(j);
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j);
                    Y(j) = abs(WLowerFluid-G(6,4)*DLowerFluid);
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
                if  length(b) == 1
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
            z = isoutlier(vertcat(X(1:i-1),SweepRange(MIN)),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
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
    addpoints(g(1),FrequencyRange(i),X(i)/1e3);
    drawnow limitrate
end
if  any(X)
    X(Misses(1:length(X)) == 1) = NaN;
end
BScholte{1}(:,1) = FrequencyRange(1:length(X));
BScholte{1}(:,2) = FrequencyRange(1:length(X))/1e3;
BScholte{1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
BScholte{1}(:,4) = fillmissing(X,'spline')/1e3;
BScholte{1}(:,7) = 0;
for p = 1:FModes 
    X = FluidVelocity(p);
    addpoints(g(p+1),FrequencyRange(1),X/1e3);
    drawnow limitrate
    Misses = 0;
    for i = 2:length(FrequencyRange)
        if  Stop == 1
            return
        end
        AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
        for m = 1:SuperLayerSize
            rw2(m) = Material{m}.Density*AngularFrequency^2;
            r2w4(m) = Material{m}.Density^2*AngularFrequency^4;
            r3w6(m) = Material{m}.Density^3*AngularFrequency^6;
            b12(m) = a12(m)*rw2(m);
            b22(m) = a22(m)*rw2(m);
            b23(m) = a23(m)*r2w4(m);
            b32(m) = a32(m)*rw2(m);
            b33(m) = a33(m)*r2w4(m);
            b34(m) = a34(m)*r3w6(m);
        end
        X(i,1) = 0;
        if  i <= 3
            SweepRange = [1.001*FluidVelocity(p) .95*FluidVelocity(p)];
        elseif i > 3
            if  abs((X(i-2)-X(i-1))) > abs((X(i-3)-X(i-2)))
                Factor = LambPhaseVelocitySweepRange1;
            else
                Factor = LambPhaseVelocitySweepRange2;
            end        
            SweepRange = [X(i-1)+3*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs((X(i-2)-X(i-1)))];
        end
        if  SweepRange(1) == SweepRange(2)
            if  i < 4 || SweepRange(2) > .95*FluidVelocity(p)
                SweepRange = [1.001*FluidVelocity(p) .95*FluidVelocity(p)];
            else
                SweepRange = [X(i-1)+3*abs(X(i-3)-X(i-2)) X(i-1)-Factor*abs((X(i-3)-X(i-2)))];
            end
        end
        if  SweepRange(1) > 1.001*FluidVelocity(p)
            SweepRange(1) = 1.001*FluidVelocity(p);
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
                WaveNumber = AngularFrequency./SweepRange;
                WaveNumber2 = WaveNumber.^2;
                WaveNumber4 = WaveNumber.^4;
                WaveNumber6 = WaveNumber.^6;
                for j = 1:length(WaveNumber)
                    for m = 1:SuperLayerSize
                        A1 = a11(m)*WaveNumber2(j)+b12(m);
                        A2 = a21(m)*WaveNumber4(j)+b22(m)*WaveNumber2(j)+b23(m);
                        A3 = a31(m)*WaveNumber6(j)+b32(m)*WaveNumber4(j)+b33(m)*WaveNumber2(j)+b34(m);
                        k3a = A2/3-A1^2/9;
                        k3b = A1^3/27-A1*A2/6+A3/2;
                        k3c = (sqrt(k3b^2+k3a^3)-k3b)^(1/3);
                        k3d = k3a/(2*k3c)-k3c/2;
                        k3e = k3a/k3c;
                        k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
                        k3(1) = sqrt(k3d-k3f-A1/3);
                        k3(2) = sqrt(k3d+k3f-A1/3);
                        k3(3) = -sqrt(k3c-k3e-A1/3);
                        k32 = k3.^2;
                        m11 = c{m}(1,1)*WaveNumber2(j)-rw2(m)+c{m}(5,5)*k32;
                        m12 = c{m}(1,6)*WaveNumber2(j)+c{m}(4,5)*k32;
                        m13 = (c{m}(1,3)+c{m}(5,5))*WaveNumber(j)*k3;
                        m22 = c{m}(6,6)*WaveNumber2(j)-rw2(m)+c{m}(4,4)*k32;
                        m23 = (c{m}(3,6)+c{m}(4,5))*WaveNumber(j)*k3;
                        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                        D3 = 1i*(c{m}(1,3)*WaveNumber(j)+c{m}(3,6)*WaveNumber(j)*V+c{m}(3,3)*k3.*W); % sigma33
                        D4 = 1i*(c{m}(4,5)*(k3+WaveNumber(j)*W)+c{m}(4,4)*k3.*V); % sigma23
                        D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(j)*W)+c{m}(4,5)*k3.*V); % sigma13
                        E = exp(1i*k3*LayerThicknesses(m));
                        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                        L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                        L{m} = L1/L2;
                    end
                    M = L{1};
                    for m = 2:SuperLayerSize
                        N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                        M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
                    end
                    MM{1} = M;
                    for m = 2:Repetitions
                        N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                        MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
                    end
                    if  SymmetricSystem
                        N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                        MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                    end
                    G = inv(MM{end});
                    if  ToggleUpperFluid && ToggleLowerFluid
                        k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j));
                        k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j));
                        WUpperFluid = k3UpperFluid/WaveNumber(j);
                        WLowerFluid = k3LowerFluid/WaveNumber(j);
                        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j); % sigma11, sigma22, sigma33 in the upper fluid
                        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j); % in the lower fluid
                        Y(j) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                    elseif ToggleUpperFluid && ~ToggleLowerFluid
                        k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j));
                        WUpperFluid = k3UpperFluid/WaveNumber(j);
                        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j);
                        Y(j) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                    elseif ~ToggleUpperFluid && ToggleLowerFluid
                        k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j));
                        WLowerFluid = k3LowerFluid/WaveNumber(j);
                        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j);
                        Y(j) = abs(WLowerFluid-G(6,4)*DLowerFluid);
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
                    if  length(b) == 1
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
                z = isoutlier(vertcat(X(1:i-1),SweepRange(MIN)),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
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
        addpoints(g(p+1),FrequencyRange(i),X(i)/1e3);
        drawnow limitrate
    end
    if  any(X)
        X(Misses(1:length(X)) == 1) = NaN;
    end
    BScholte{p+1}(:,1) = FrequencyRange(1:length(X));
    BScholte{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
    BScholte{p+1}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
    BScholte{p+1}(:,4) = fillmissing(X,'spline')/1e3;    
    BScholte{p+1}(:,7) = 0;
    if  p == 2 && ~isempty(H) && BScholte{p+1}(end,1) > H(1)    
        H(1) = [];
    end
end
for p = 1:length(BScholte)
    if  all(BScholte{p}(:,4) == 0)
        x(p) = 1;
    else
        x(p) = 0;
    end
end
BScholte(x == 1) = [];