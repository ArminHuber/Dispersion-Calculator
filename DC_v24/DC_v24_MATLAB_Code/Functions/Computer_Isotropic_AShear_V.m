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
function AShear = Computer_Isotropic_AShear_V(Material,FrequencyRange,Thickness,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,ASteps,AOffset)
%#ok<*AGROW>
T = Material.LongitudinalAttenuation+Material.TransverseAttenuation;
for p = 1:length(H)
    i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange);
    if  isempty(i)
        continue
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    WaveNumber = sqrt((AngularFrequency/Material.TransverseVelocity_complex).^2-((2*p-1)*pi/Thickness)^2);
    AShear{p}(i,4) = AngularFrequency./real(WaveNumber)/1e3;
    AShear{p}(i,6) = abs(imag(WaveNumber));
    AShear{p}(:,1) = FrequencyRange(1:height(AShear{p}));
    AShear{p}(:,2) = FrequencyRange(1:height(AShear{p}))/1e3;
    AShear{p}(:,3) = FrequencyRange(1:height(AShear{p}))*Thickness;
    if  max(AShear{p}(:,4))*1e3 > PhaseVelocityLimit
        X1{p}(1,1) = 0;
    else
        [~,MaxInd] = max(AShear{p}(:,4));
        AngularFrequency = 2*pi*AShear{p}(MaxInd,1)*1e3;
        WaveNumber = AngularFrequency/(AShear{p}(MaxInd,4)*1e3)+1i*AShear{p}(MaxInd,6);
        PhaseVelocityComplex = AngularFrequency/WaveNumber;
        PhaseVelocityRange = real(PhaseVelocityComplex)+ASteps:ASteps:PhaseVelocityLimit+ASteps;
        for i = 1:length(PhaseVelocityRange)
            X1{p}(i,1) = 0;
            for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                if  i == 1
                    SweepRangeFrq = [AShear{p}(MaxInd,1)-AOffset/Thickness/1e3 AShear{p}(MaxInd,1)+4/Thickness/1e3];
                    SweepRangeImag = [SearchWidthImag(1)*imag(PhaseVelocityComplex) SearchWidthImag(2)*imag(PhaseVelocityComplex)];
                else
                    SweepRangeFrq = [X1{p}(end-1,1)-AOffset/Thickness/1e3 X1{p}(end-1,1)+4/Thickness/1e3];
                    SweepRangeImag = [SearchWidthImag(1)*X1{p}(end-1,3) SearchWidthImag(2)*X1{p}(end-1,3)];
                end
                if  all(SweepRangeImag == [0 0])
                    SweepRangeImag = [-20*T 0];
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
                for o = 1:SearchAreaSections % increase search resolution
                    if  isempty(SweepRangeFrq) || isempty(SweepRangeImag)
                        break
                    end
                    for k = 1:1e2 % search minimum in characteristic equation and converge upon it
                        if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                            if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                SweepRangeFrq = SweepRangeFrq(1):.25^o/q*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeFrq(1)-SweepRangeFrq(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                SweepRangeFrq = SweepRangeFrq(1):.25/q/o*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                SweepRangeFrq = SweepRangeFrq(1):.25/q*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            end
                        else
                            if  length(SweepRangeFrq) == 2
                                SweepRangeFrq = SweepRangeFrq(1):.25*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                            end
                            if  length(SweepRangeImag) == 2
                                SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            end
                        end
                        if  isempty(SweepRangeFrq) || isempty(SweepRangeImag)
                            break
                        end
                        AngularFrequency = 2*pi*SweepRangeFrq*1e3;
                        WaveNumber = repmat(AngularFrequency',1',length(SweepRangeImag))./(repmat((PhaseVelocityRange(i)+SweepRangeImag*1i),length(SweepRangeFrq),1));
                        y = sqrt((AngularFrequency'/Material.TransverseVelocity_complex).^2-WaveNumber.^2);
                        Y = abs(cos(y*Thickness/2));
                        if  abs(SweepRangeImag(end)) < 1e-3 % in case the fluid has low density, this allows finding solutions at zero attenuation for A0
                            Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                        end
% if  i>1
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeFrq,20*log10(Y))
% else
% f = figure;surf(SweepRangeImag,SweepRangeFrq,20*log10(Y))
% end
% close(f)
% end
                        Min = zeros(size(Y,1),size(Y,2));
                        for l = 2:size(Y,2)-1
                            for j = 2:size(Y,1)-1
                                if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                    Min(j,l) = 1;
                                end
                            end
                        end
                        [b1,b2] = find(Min);
                        if  ~isempty(b1) % one or multiple minima are found
                            if  length(b1) == 1
                                MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                            else
                                delta = [];
                                for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                    frq(l) = AngularFrequency(b1(l))/(2*pi)/1e3;
                                    if  i == 1
                                        delta(l) = abs(frq(l)-AShear{p}(MaxInd,1));
                                    else
                                        delta(l) = abs(frq(l)-X1{p}(end-1,1));
                                    end
                                end
                                [~,l] = min(delta);
                                MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                            end
                        else
                            Min = zeros(size(Y,1),size(Y,2)); % find border minima
                            for j = 2:size(Y,1)-1
                                if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                    Min(j,1) = 1;
                                end
                                if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
                                    Min(j,size(Y,2)) = 1;
                                end
                            end
                            for j = 2:size(Y,2)-1
                                if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                    Min(1,j) = 1;
                                end
                                if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
                                    Min(size(Y,1),j) = 1;
                                end
                            end
                            [b1,b2] = find(Min);
                            if  ~isempty(b1) % one or multiple BORDER minima are found
                                if  length(b1) == 1
                                    MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                else
                                    Value = [];
                                    for l = 1:length(b1)
                                        Value(l) = Y(b1(l),b2(l));
                                    end
                                    [~,l] = min(Value);
                                    MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                end
                            else
                                if  q > 1% && o == 1 % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                                    SweepRangeFrq = [SweepRangeFrq(1)+(q-1)*abs(SweepRangeFrq(1)-SweepRangeFrq(end)) SweepRangeFrq(end)-(q-1)*abs(SweepRangeFrq(1)-SweepRangeFrq(end))];
                                    SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                    if  SweepRangeImag(2) > 0
                                        SweepRangeImag(2) = 0;
                                    end
                                end
                                MIN = 0;
                                break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                            end
                        end 
                        if  k == 100 || (Resolution > abs(SweepRangeFrq(1)-SweepRangeFrq(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                            break
                        end
                        if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                            if  MIN(1) == 1
                                if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeFrq(1)-SweepRangeFrq(end))
                                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                end
                                SweepRangeFrq = [SweepRangeFrq(1)+4*abs(SweepRangeFrq(1)-SweepRangeFrq(end)) SweepRangeFrq(end)];
                            elseif MIN(2) == 1
                                if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                end
                                SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                            elseif MIN(1) == size(Y,1)
                                if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeFrq(1)-SweepRangeFrq(end))
                                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                end
                                SweepRangeFrq = [SweepRangeFrq(1) SweepRangeFrq(end)-4*abs(SweepRangeFrq(1)-SweepRangeFrq(end))];
                            elseif MIN(2) == size(Y,2)
                                if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                end
                                SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                            end
                        else
                            if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                if  Resolution < abs(SweepRangeFrq(1)-SweepRangeFrq(2))
                                    SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                end
                            elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeFrq(1)-SweepRangeFrq(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                if  Resolution < abs(SweepRangeFrq(1)-SweepRangeFrq(2))
                                    SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                end
                                if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                end
                            elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                end
                            end
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                    end
                    if  any(MIN)
                        z = isoutlier(vertcat(X1{p}(1:end-1,1),AngularFrequency(MIN(1))/(2*pi)/1e3),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                        if  ~z(end) || all(X1{p}(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                            X1{p}(i,4) = AngularFrequency(MIN(1))/real(WaveNumber(MIN(1),MIN(2))); % phase velocity (m/s)
                            X1{p}(i,1) = AngularFrequency(MIN(1))/(2*pi)/1e3; % frequency (kHz)
                            X1{p}(i,2) = imag(WaveNumber(MIN(1),MIN(2))); % attenuation (Np/m)
                            X1{p}(i,3) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                            break
                        end
                        if  i == 1
                            SweepRangeFrq = [AShear{p}(MaxInd,1)-AOffset/Thickness/1e3 AShear{p}(MaxInd,1)+4/Thickness/1e3];
                            SweepRangeImag = [SearchWidthImag(1)*imag(PhaseVelocityComplex) SearchWidthImag(2)*imag(PhaseVelocityComplex)];
                        else
                            SweepRangeFrq = [X1{p}(end-1,1)-AOffset/Thickness/1e3 X1{p}(end-1,1)+4/Thickness/1e3];
                            SweepRangeImag = [SearchWidthImag(1)*X1{p}(end-1,3) SearchWidthImag(2)*X1{p}(end-1,3)];
                        end
                        if  all(SweepRangeImag == [0 0])
                            SweepRangeImag = [-20*T 0];
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                    end
                end
                if  X1{p}(i) > 0 % stop q-loop if minimum has been found
                    break
                end
            end
            if  X1{p}(i) == 0
                break
            end
            if  X1{p}(i,2) < 0 % negative attenuation is impossible
                X1{p}(i,2:3) = 0;
            end
        end
        if  length(X1) < p
            X1{p}(1,1) = 0;
        end
        if  X1{p}(1,1) > 0
            X1{p}(:,6) = X1{p}(:,2);
            X1{p}(:,2) = X1{p}(:,1)/1e3;
            X1{p}(:,3) = X1{p}(:,1)*Thickness;
            X1{p}(:,4) = X1{p}(:,4)/1e3;           
            X1{p}(X1{p}(:,1) == 0,:) = []; 
            X1{p} = flipud(X1{p});
        end
    end
end
if  ~exist('AShear','var')
    AShear{1} = 0;
else
    for p = 1:length(AShear)
        AShear{p}(AShear{p}(:,4) == 0,:) = [];
        if  X1{p}(1,1) > 0
            AShear{p} = vertcat(X1{p},AShear{p});
        end
    end
end