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
function SScholte = Computer_Isotropic_SScholte(Multithreading,Q,ax,Fluid,Material,H,FrequencyRange,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS> 
global Stop
Stop = 0;
if  ~Multithreading
    for p = 1:length(H)+1
        g(p) = animatedline(ax,'LineStyle','-.','color','r');
    end
end
if  Material.PlateVelocity > Fluid.Velocity
    X = Fluid.Velocity;
else
    X = Material.PlateVelocity;
end
if  Multithreading
    send(Q,[FrequencyRange(1),X(1)/1e3,1])
else
    addpoints(g(1),FrequencyRange(1),X(1)/1e3);
    drawnow limitrate
end
for i = 2:length(FrequencyRange)
    if  Stop == 1
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
    kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
    kF2 = (AngularFrequency/Fluid.Velocity)^2;
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
    if  i <= 3
        if  Material.PlateVelocity > Fluid.Velocity
            SweepRange = [1.001*Fluid.Velocity .99*Fluid.Velocity];
        else
            SweepRange = [1.001*Material.PlateVelocity .99*Material.PlateVelocity];
        end
    else
        if  delta > abs(X(i-3)-X(i-2))
            Factor = LambPhaseVelocitySweepRange1;
        else
            Factor = LambPhaseVelocitySweepRange2;
        end
        SweepRange = [X(i-1)+3*delta X(i-1)-Factor*abs((X(i-2)-X(i-1)))];
    end
    if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
        if  delta > PhaseVelocityResolution
            SweepRange(1) = SweepRange(1)+PhaseVelocityResolution;
            SweepRange(2) = SweepRange(2)-PhaseVelocityResolution;
        else
            SweepRange(1) = SweepRange(1)+2*delta;
            SweepRange(2) = SweepRange(2)-2*delta;
        end
    end
    if  Material.PlateVelocity > Fluid.Velocity
        if  SweepRange(1) > 1.001*Fluid.Velocity
            SweepRange(1) = 1.001*Fluid.Velocity;
        end
    else
        if  SweepRange(1) > 1.001*Material.PlateVelocity
            SweepRange(1) = 1.001*Material.PlateVelocity;
        end
    end
    if  SweepRange(2) < 0
        SweepRange(2) = 0;
    end
    for o = 1:PhaseVelocitySections % increase search resolution
        if  isempty(SweepRange)
            break
        end
        if  i > 2 && delta < 1e-10
            Bisections = ceil(log2(delta/1e-10)/log2(.5));
        else
            Bisections = ceil(log2(1e-10/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
        end
        if  i <= 3 && Bisections < 10
            Bisections = 10;
        elseif i > 3 && Bisections < 1
            Bisections = 1;
        end
        for k = 1:Bisections % search minimum in characteristic equation and converge upon it
            if  k == 1 % divide SweepRange into sections where the characteristic equation is evaluated
                SweepRange = SweepRange(1):.25^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
            else
                SweepRange = SweepRange(1):.25*(SweepRange(end)-SweepRange(1)):SweepRange(end);
            end
            MIN = 0;
            for j = 1:length(SweepRange)
                k2 = (AngularFrequency/SweepRange(j))^2;
                x = sqrt(kL2-k2);
                y = sqrt(kT2-k2);
                Y(j) = abs((y^2-k2)^2/(y*tan(x*Half))+4*k2*x/tan(y*Half)-1i*Fluid.Density*kT2^2*x/(y*Material.Density*sqrt(kF2-k2)));
                if  j > 2 && Y(j-1) < Y(j-2) && Y(j-1) < Y(j)
                    MIN = j-1;
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
            if  ~z(end) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
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
SScholte{1}(:,3) = FrequencyRange(1:length(X))*2*Half;
SScholte{1}(:,4) = fillmissing(X,'spline')/1e3;
SScholte{1}(:,6) = 0;
if  HigherOrderModes && any(H)
    for p = 1:length(H) 
        X = [];
        Misses = 0;
        for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
            if  Stop == 1
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            kL2 = (AngularFrequency/Material.LongitudinalVelocity)^2;
            kT2 = (AngularFrequency/Material.TransverseVelocity)^2;
            kF2 = (AngularFrequency/Fluid.Velocity)^2;
            X(i,1) = 0;
            if  i <= ceil(H(p)/FrequencyResolution)+5
                SweepRange = [1.001*Fluid.Velocity .99*Fluid.Velocity];
            elseif i > ceil(H(p)/FrequencyResolution)+5
                if  abs((X(i-2)-X(i-1))) > abs((X(i-3)-X(i-2)))
                    Factor = LambPhaseVelocitySweepRange1;
                else
                    Factor = LambPhaseVelocitySweepRange2;
                end 
                SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs((X(i-2)-X(i-1)))];
            end
            if  SweepRange(1) == SweepRange(2)
                if  i < ceil(H(p)/FrequencyResolution)+4
                    SweepRange = [1.001*Fluid.Velocity .99*Fluid.Velocity];
                else
                    SweepRange = [X(i-1)+.1*abs(X(i-3)-X(i-2)) X(i-1)-Factor*abs((X(i-3)-X(end-2)))];
                end
            end
            if  SweepRange(1) > 1.001*Fluid.Velocity
                SweepRange(1) = 1.001*Fluid.Velocity;
            end
            if  SweepRange(2) < SScholte{p}(i,4)*1e3+.1
                SweepRange(2) = SScholte{p}(i,4)*1e3+.1;
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
                    for j = 1:length(SweepRange)
                        k2 = (AngularFrequency/SweepRange(j))^2;
                        x = sqrt(kL2-k2);
                        y = sqrt(kT2-k2);
                        Y(j) = abs((y^2-k2)^2/(y*tan(x*Half))+4*k2*x/tan(y*Half)-1i*Fluid.Density*kT2^2*x/(y*Material.Density*sqrt(kF2-k2)));
                        if  j > 2 && Y(j-1) < Y(j-2) && Y(j-1) < Y(j)
                            MIN = j-1;
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
                    z = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:i-1),SweepRange(MIN)),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                    if  ~z(end) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                        X(i,1) = SweepRange(MIN); % phase velocity (m/s)
                        Misses(i) = 0;
                        break
                    end
                end
            end
            if  X(i) == 0 % fit phase velocity where we missed the solution to obtain useful sweep range for the next frequency step
                Fit = fit(FrequencyRange(1:i-1)',X(1:i-1),'cubicspline');
                X(i,1) = Fit(FrequencyRange(i));
                Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
            end
            if  Multithreading
                send(Q,[FrequencyRange(i),X(i)/1e3,p+1])
            else
                addpoints(g(p+1),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
        end
        if  any(X)
            X(Misses(1:length(X)) == 1) = NaN;
        end
        SScholte{p+1}(:,1) = FrequencyRange(1:length(X));
        SScholte{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
        SScholte{p+1}(:,3) = FrequencyRange(1:length(X))*2*Half;
        SScholte{p+1}(:,4) = fillmissing(X,'spline')/1e3;
    end
    for p = 2:length(SScholte)
        SScholte{p}(SScholte{p}(:,4) == 0,:) = [];
        SScholte{p}(:,6) = 0;
    end
end