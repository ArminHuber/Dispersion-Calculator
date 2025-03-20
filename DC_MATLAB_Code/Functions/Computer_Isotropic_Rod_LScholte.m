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
function LScholte = Computer_Isotropic_Rod_LScholte(Multithreading,Q,ax,Fluid,Material,H,FrequencyRange,PhaseVelocityResolution,PhaseVelocitySections,R,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS> 
global Stop
Stop = 0;
LScholte{1} = [];
R2 = R^2;
Density = Fluid.Density/Material.Density;
if  Material.CylinderVelocity <= Fluid.Velocity
    CylinderWaveCase = 1;
    H = [FrequencyRange(1) H];
else
    CylinderWaveCase = 0;
end
for p = 1:length(H)
    X = [];
    Misses = 0;
    if  CylinderWaveCase && p == 1
        X = Material.CylinderVelocity;
    end
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
            if  ~CylinderWaveCase || (CylinderWaveCase && p > 1)
                SweepRange = [Fluid.Velocity .99*Fluid.Velocity];
            elseif CylinderWaveCase && p == 1
                SweepRange = [1.001*Material.CylinderVelocity .99*Material.CylinderVelocity];
            end
        else
            if  abs(X(i-2)-X(i-1)) > abs(X(i-3)-X(i-2))
                Factor = LambPhaseVelocitySweepRange1;
            else
                Factor = LambPhaseVelocitySweepRange2;
            end 
            SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
        end
        if  SweepRange(1) == SweepRange(2)
            if  i < ceil(H(p)/FrequencyResolution)+4
                if  ~CylinderWaveCase || (CylinderWaveCase && p > 1)
                    SweepRange = [Fluid.Velocity .99*Fluid.Velocity];
                elseif CylinderWaveCase && p == 1
                    SweepRange = [1.001*Material.CylinderVelocity .99*Material.CylinderVelocity];
                end
            else
                SweepRange = [X(i-1)+.1*abs(X(i-3)-X(i-2)) X(i-1)-10*abs(X(i-3)-X(i-2))];
            end
        end
        if  (~CylinderWaveCase || (CylinderWaveCase && p > 1)) && SweepRange(1) > Fluid.Velocity
            SweepRange(1) = Fluid.Velocity;
        elseif CylinderWaveCase && p == 1 && SweepRange(1) > 1.001*Material.CylinderVelocity
            SweepRange(1) = 1.001*Material.CylinderVelocity;
        end
        if  p > 1 && i <= height(LScholte{p-1}) && SweepRange(2) < LScholte{p-1}(i,4)*1e3+.1
            SweepRange(2) = LScholte{p-1}(i,4)*1e3+.1;
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
            if  o == 0 && SweepRange(1) > Material.TransverseVelocity && SweepRange(2) < Material.TransverseVelocity % Neighbors cannot be excluded because SweepRange has only 2 elements
                continue
            end
            for j = 1:length(SweepRange)-1
                if  SweepRange(j) > Material.TransverseVelocity && SweepRange(j+1) < Material.TransverseVelocity
                    if  j == 1
                        SweepRange(2) = NaN;
                    else
                        SweepRange(j) = NaN;
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
                            k1 = AngularFrequency/PhaseVelocity(l);
                            k2 = k1^2;
                            y2 = kT2-k2;
                            xR = sqrt(kL2-k2)*R;
                            yR = sqrt(kT2-k2)*R;
                            J0x = besselj(0,xR);
                            J0y = besselj(0,yR);
                            J1xR = xR*besselj(1,xR);
                            J1yR = yR*besselj(1,yR);
                            zR = sqrt(kF2-k2)*R;
                            H0z = besselh(0,zR);
                            H1zR = -zR*besselh(1,zR);
                            M(1,1) = R2*(.5*kT2-k2)*J0x-J1xR;
                            M(1,2) = k1*(J1yR-y2*R2*J0y);
                            M(1,3) = .5*kT2*R2*Density*H0z;
                            M(2,1) = -2*k1*J1xR;
                            M(2,2) = (k2-y2)*J1yR;
                            M(3,1) = -J1xR;
                            M(3,2) = k1*J1yR;
                            M(3,3) = H1zR;
                            Y(j,l) = det(M);
                        end
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
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
    end
    if  any(X)
        X(Misses(1:length(X)) == 1) = NaN;
    end
    LScholte{p}(:,1) = FrequencyRange(1:length(X));
    LScholte{p}(:,2) = FrequencyRange(1:length(X))/1e3;
    LScholte{p}(:,3) = FrequencyRange(1:length(X))*2*R;
    LScholte{p}(:,4) = fillmissing(X,'spline')/1e3;
    if  Multithreading
        send(Q,{LScholte{p},'-.','r'})
    else
        line(ax,LScholte{p}((LScholte{p}(:,4) ~= 0),1),LScholte{p}((LScholte{p}(:,4) ~= 0),4),'LineStyle','-.','color','r')
        drawnow limitrate
    end
end
for p = 1:length(LScholte)
    LScholte{p}(LScholte{p}(:,4) == 0,:) = [];
    LScholte{p}(:,6) = 0; 
end