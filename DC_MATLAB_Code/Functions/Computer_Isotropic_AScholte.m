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
function AScholte = Computer_Isotropic_AScholte(Multithreading,Q,ax,Fluid,Material,FScholte,H,FrequencyRange,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS> 
global Stop 
Stop = 0;
AScholte{1} = [];
if  ~Multithreading
    for p = 1:length(H)+1
        g(p) = animatedline(ax,'LineStyle','-.','color','b');
    end
end
X = FScholte(1);
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
        SweepRange = [6*X(i-1) X(i-1)];
    else
        if  delta > abs(X(i-3)-X(i-2))
            Factor = LambPhaseVelocitySweepRange1;
        else
            Factor = LambPhaseVelocitySweepRange2;
        end
        SweepRange = [X(i-1)+3*delta X(i-1)-Factor*delta];
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
    if  SweepRange(1) > 1.001*Fluid.Velocity
        SweepRange(1) = 1.001*Fluid.Velocity;
    end
    if  SweepRange(2) < 0
        SweepRange(2) = 0;
    end
%     if  Fluid.Density < 20 % for gasesous media (density < 20 kg/m^3), we use search for minima, while for liquid media we use search for sign change
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
                    Y(j) = abs((y^2-k2)^2/y*tan(x*Half)+4*k2*x*tan(y*Half)+1i*Fluid.Density*kT2^2*x/(y*Material.Density*sqrt(kF2-k2)));
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
%                 z = isoutlier(vertcat(X(1:i-1),SweepRange(MIN)),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
%                 if  ~z(end) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                    X(i,1) = SweepRange(MIN); % phase velocity (m/s)
                    Misses(i) = 0;
                    break
%                 end
            end
        end
%     else
%         for o = 0:PhaseVelocitySections
%             if  o > 0
%                 SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
%             end
%             if  i > 2 && delta < 1e-10
%                 Bisections = ceil(log2(delta/1e-10)/log2(.5));
%             else
%                 Bisections = ceil(log2(1e-10/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
%             end
%             if  i <= 3 && Bisections < 10
%                 Bisections = 10;
%             elseif i > 3 && Bisections < 1
%                 Bisections = 1;
%             end
%             for j = 1:length(SweepRange)-1
%                 if  j == 1
%                     PhaseVelocityIndices = [1 2 3];
%                 else
%                     PhaseVelocityIndices = [2 3];
%                 end
%                 PhaseVelocity = [SweepRange(j) SweepRange(j)+(SweepRange(j+1)-SweepRange(j))/2 SweepRange(j+1)];
%                 for k = 1:Bisections
%                     if  k > 1
%                         PhaseVelocityIndices = 2;
%                     end
%                     for l = PhaseVelocityIndices(1):PhaseVelocityIndices(end)
%                         k2 = (AngularFrequency/PhaseVelocity(l))^2;
%                         x = sqrt(kL2-k2);
%                         y = sqrt(kT2-k2);
%                         Y(j,l) = (y^2-k2)^2/y*tan(x*Half)+4*k2*x*tan(y*Half)+1i*Fluid.Density*kT2^2*x/(y*Material.Density*sqrt(kF2-k2));
%                     end
%                     if  k == 1
%                         Y(j+1,1) = Y(j,3);
%                     end
%                     if  (j == 1 && abs(Y(j,2)) < 1e1) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
%                         PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
%                         Y(j,3) = Y(j,2);
%                     elseif (j == 1 && abs(Y(j,2)) < 1e1) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
%                         PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
%                         Y(j,1) = Y(j,2);
%                     else
%                         PhaseVelocity(2) = 0;
%                         break
%                     end
%                 end
%                 if  PhaseVelocity(2) > 0
% %                     if  i < 4
% %                         Outlier = 0;
% %                     else
% %                         z = isoutlier(vertcat(X(1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
% %                         if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
% %                             Outlier = 1;
% %                         else
% %                             Outlier = 0;
% %                         end
% %                     end
% %                     if  ~Outlier || all(X == 0)
%                         X(i,1) = PhaseVelocity(2);
%                         Misses(i) = 0;
%                         break
% %                     end
%                 end
%             end
%             if  PhaseVelocity(2) > 0
%                 break
%             end
%         end
%     end
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
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
end
if  any(X)
    X(Misses(1:length(X)) == 1) = NaN;
end
AScholte{1}(:,1) = FrequencyRange(1:length(X));
AScholte{1}(:,2) = FrequencyRange(1:length(X))/1e3;
AScholte{1}(:,3) = FrequencyRange(1:length(X))*2*Half;
AScholte{1}(:,4) = fillmissing(X,'spline')/1e3;
AScholte{1}(:,6) = 0;
if  HigherOrderModes && any(H)
    for p = 1:length(H)
        X = [];
        Misses = 0;
        for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
            if  Stop
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
                if  abs(X(i-2)-X(i-1)) > abs(X(i-3)-X(i-2))
                    Factor = LambPhaseVelocitySweepRange1;
                else
                    Factor = LambPhaseVelocitySweepRange2;
                end 
                SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
            end
            if  SweepRange(1) == SweepRange(2)
                if  i < ceil(H(p)/FrequencyResolution)+4
                    SweepRange = [1.001*Fluid.Velocity .99*Fluid.Velocity];
                else
                    SweepRange = [X(i-1)+.1*abs(X(i-3)-X(i-2)) X(i-1)-Factor*abs(X(i-3)-X(i-2))];
                end
            end
            if  SweepRange(1) > 1.001*Fluid.Velocity
                SweepRange(1) = 1.001*Fluid.Velocity;
            end
            if  i <= height(AScholte{p}) && SweepRange(2) < AScholte{p}(i,4)*1e3+.1
                SweepRange(2) = AScholte{p}(i,4)*1e3+.1;
            end
            for o = 0:PhaseVelocitySections
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
                            k2 = (AngularFrequency/PhaseVelocity(l))^2;
                            x = sqrt(kL2-k2);
                            y = sqrt(kT2-k2);
                            Y(j,l) = (y^2-k2)^2/y*tan(x*Half)+4*k2*x*tan(y*Half)+1i*Fluid.Density*kT2^2*x/(y*Material.Density*sqrt(kF2-k2));
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  (j == 1 && abs(Y(j,2)) < 1e-1) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
                            PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif (j == 1 && abs(Y(j,2)) < 1e-1) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
                            PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                            Y(j,1) = Y(j,2);
                        else
                            PhaseVelocity(2) = 0;
                            break
                        end
                    end
                    if  PhaseVelocity(2) > 0
%                         if  numel(find(X > 0)) <= 5
%                             Outlier = 0;
%                         else
%                             z = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
%                             if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
%                                 Outlier = 1;
%                             else
%                                 Outlier = 0;
%                             end
%                         end
%                         if  ~Outlier || all(X == 0)
                            X(i,1) = PhaseVelocity(2);
                            Misses(i) = 0;
                            break
%                         end
                    end
                end
                if  PhaseVelocity(2) > 0
                    break
                end
            end
            if  X(i) == 0 % fit phase velocity where we missed the solution to obtain useful sweep range for the next frequency step
                Fit = fit(FrequencyRange(1:i-1)',X(1:i-1),'cubicspline');
                X(i,1) = Fit(FrequencyRange(i));
                Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end))
                X(end-MissingSamples:end) = [];
                Misses(end-MissingSamples:end) = 0;
                break
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
        AScholte{p+1}(:,1) = FrequencyRange(1:length(X));
        AScholte{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
        AScholte{p+1}(:,3) = FrequencyRange(1:length(X))*2*Half;
        AScholte{p+1}(:,4) = fillmissing(X,'spline')/1e3;
    end
    for p = 2:length(AScholte)
        AScholte{p}(AScholte{p}(:,4) == 0,:) = [];
        AScholte{p}(:,6) = 0;
    end
end