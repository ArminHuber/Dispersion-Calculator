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
function SH = Computer_Anisotropic_SH_D(Multithreading,Q1,Q2,ax,ModeFamily,Material,FrequencyRange,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0)        
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
if  SuperLayerSize > 1
    if  ~Multithreading
        for p = 1:length(H)+N
            if  ModeFamily == 1
                g(p) = animatedline(ax,'LineStyle','--','color','r');
                g1(p) = animatedline(ax,'LineStyle','--','color','r');
            elseif ModeFamily == 2
                g(p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
                g1(p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
            elseif ModeFamily == 3
                g(p) = animatedline(ax,'LineStyle','--','color','b');
                g1(p) = animatedline(ax,'LineStyle','--','color','b');
            end
        end
    end
    for i = 1:length(Material)
        if  ~isreal(Material{i}.C)
            T = 1e2*pi*(imag(Material{i}.C(1,1))/real(Material{i}.C(1,1))+imag(Material{i}.C(6,6))/real(Material{i}.C(6,6)));
            break
        end
    end
end
if  ModeFamily < 3
    if  SuperLayerSize == 1
        AngularFrequency = 2*pi*FrequencyRange*1e3;
        Wavenumber = AngularFrequency*sqrt(Material{1}.Density/c{1}(6,6));
        SH{1}(:,1) = FrequencyRange;
        SH{1}(:,2) = FrequencyRange/1e3;
        SH{1}(:,3) = FrequencyRange*PlateThickness;
        SH{1}(:,4) = AngularFrequency./real(Wavenumber)/1e3;
        SH{1}(:,7) = abs(imag(Wavenumber));
        SH{1}(end,9) = 0;
    else
        for i = 1:length(FrequencyRange)
            if  Stop
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            AngularFrequency2 = AngularFrequency^2;
            for m = 1:SuperLayerSize
                rw2(m) = Material{m}.Density*AngularFrequency2;
            end
            X(i,1) = 0;
            for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                if  i < 4
                    SweepRangeReal = [1.01*XSH0 .99*XSH0];
                else
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                end
                if  SweepRangeReal(1) == SweepRangeReal(2)
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                end
                if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                    SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                    SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                end
                if  SweepRangeReal(2) < 0
                    SweepRangeReal(2) = 0;
                end
                if  i == 1
                    SweepRangeImag = [-10*T 0]; % the imaginary phase velocity corresponds to the attenuation
                elseif i == 2 || i == 3
                    SweepRangeImag = [-20*T 0]; % the imaginary phase velocity corresponds to the attenuation
                else
                    SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-20*T 0];
                    elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                        SweepRangeImag(2) = 0;
                    end
                end
                if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                    SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                    SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
                for o = 1:SearchAreaSections % increase search resolution
                    if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                        break
                    end
                    for k = 1:1e2
                        if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                            if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                SweepRangeReal = SweepRangeReal(1):.25^o/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                SweepRangeReal = SweepRangeReal(1):.25/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            end
                        else
                            if  length(SweepRangeReal) == 2
                                SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                            end
                            if  length(SweepRangeImag) == 2
                                SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            end
                        end
                        if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                            break
                        end
                        Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
                        Y = Computer(1,ModeFamily,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,c,rw2);
                        if  abs(SweepRangeImag(end)) < 1e-3
                            Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
                        end
% if  i>=371
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))    
% else
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
% close(f)
% end
                        Min = zeros(size(Y));
                        for l = 2:size(Y,2)-1
                            for j = 2:size(Y,1)-1
                                if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                    Min(j,l) = 1;
                                end
                            end
                        end
                        [b1,b2] = find(Min);
                        if  ~isempty(b1) % one or multiple minima are found
                            if  isscalar(b1)
                                MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            else
                                delta = [];
                                for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                    cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                                    if  i == 1
                                        delta(l) = abs(cp(l)-XSH0);
                                    else
                                        delta(l) = abs(cp(l)-X(end-1,1));
                                    end
                                end
                                [~,l] = min(delta);
                                MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            end
                        else
                            Min = zeros(size(Y)); % find border minima
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
                                if  isscalar(b1)
                                    MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                else
                                    Value = [];
                                    for l = 1:length(b1)
                                        Value(l) = Y(b1(l),b2(l));
                                    end
                                    [~,l] = min(Value);
                                    MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                end
                            else
                                if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal x SweepRangeImag)
                                    SweepRangeReal = [SweepRangeReal(1)+(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                    SweepRangeImag = [SweepRangeImag(1)-(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                    if  SweepRangeReal(2) < 0
                                        SweepRangeReal(2) = 0;
                                    end
                                    if  SweepRangeImag(2) > 0
                                        SweepRangeImag(2) = 0;
                                    end
                                end
                                MIN = 0;
                                break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                            end
                        end
                        if  k == 100 || (Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                            break
                        end
                        if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                            if  MIN(1) == 1
                                if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                end
                                SweepRangeReal = [SweepRangeReal(1)+4*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
                            elseif MIN(2) == 1
                                if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                end
                                SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                            elseif MIN(1) == size(Y,1)
                                if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                end
                                SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                            elseif MIN(2) == size(Y,2)
                                if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                end
                                SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                            end
                        else
                            if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                                    SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                end
                            elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                                    SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                end
                                if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                end
                            elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                end
                            end
                        end
                        if  SweepRangeReal(2) < 0
                            SweepRangeReal(2) = 0;
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                    end
                    if  any(MIN)
                        if  i < 4
                            Outlier = 0;
                        else
                            z1 = isoutlier(vertcat(X(1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                            z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                            if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)
                                Outlier = 1;
                            else
                                Outlier = 0;
                            end
                        end
                        if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                            X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                            X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                            X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                            X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                            Misses(i) = 0;
                            break
                        end
                        if  i < 4
                            SweepRangeReal = [1.01*XSH0 .99*XSH0];
                        else
                            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                        end
                        if  SweepRangeReal(1) == SweepRangeReal(2)
                            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                        end
                        if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                            SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                            SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                        end
                        if  SweepRangeReal(2) < 0
                            SweepRangeReal(2) = 0;
                        end
                        if  i == 1
                            SweepRangeImag = [-10*T 0]; % the imaginary phase velocity corresponds to the attenuation
                        elseif i == 2 || i == 3
                            SweepRangeImag = [-20*T 0]; % the imaginary phase velocity corresponds to the attenuation
                        else
                            SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                            if  all(SweepRangeImag == [0 0])
                                SweepRangeImag = [-20*T 0];
                            elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                            SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                            SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                    end
                end
                if  X(i,1) > 0 % stop q-loop if minimum has been found
                    break
                end
            end
            if  i == 1 && isinf(X(1)) || X(1) == 0
                X = [XSH0 0 XSH0 0];
            end
            if  X(i,1) == 0 % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
                Smooth1 = filloutliers(X(1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
                Smooth2 = filloutliers(X(1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
                Smooth3 = filloutliers(X(1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
                Smooth4 = filloutliers(X(1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
                Fit1 = fit(FrequencyRange(1:i-1)',Smooth1,'cubicspline');
                Fit2 = fit(FrequencyRange(1:i-1)',Smooth2,'cubicspline');
                Fit3 = fit(FrequencyRange(1:i-1)',Smooth3,'cubicspline');
                Fit4 = fit(FrequencyRange(1:i-1)',Smooth4,'cubicspline');
                X(i,1) = Fit1(FrequencyRange(i));
                X(i,2) = Fit2(FrequencyRange(i));
                X(i,3) = Fit3(FrequencyRange(i));
                X(i,4) = Fit4(FrequencyRange(i));
                if  X(i,2) < 0 % negative attenuation is impossible
                    X(i,[2 4]) = 0;
                end
                Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
                X(end-MissingSamples:end,:) = [];
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
        X(Misses(1:height(X)) == 1,:) = NaN;
        SH{1}(:,1) = FrequencyRange(1:height(X));
        SH{1}(:,2) = FrequencyRange(1:height(X))/1e3;
        SH{1}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
        SH{1}(:,4) = fillmissing(X(:,1),'spline')/1e3;
        SH{1}(:,7) = fillmissing(X(:,2),'spline');
        SH{1}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
        SH{1}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
        SH{1}(SH{1}(:,7) < 0,7) = 0; % negative attenuation is impossible
    end
end
if  HigherOrderModes && any(H)
    MissingModes(length(H)+N) = 0;
    for p = 1:length(H)
        if  SuperLayerSize == 1
            i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange);
            if  isempty(i)
                continue
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            if  ModeFamily == 1
                Wavenumber = sqrt((Material{1}.Density*AngularFrequency.^2-4*c{1}(4,4)*(p*pi/PlateThickness)^2)/c{1}(6,6));
            elseif ModeFamily == 3
                Wavenumber = sqrt((Material{1}.Density*AngularFrequency.^2-c{1}(4,4)*((2*p-1)*pi/PlateThickness)^2)/c{1}(6,6));
            end
            c_complex = AngularFrequency./Wavenumber;
            SH{p+N}(i,4) = AngularFrequency./real(Wavenumber)/1e3;
            SH{p+N}(i,7) = abs(imag(Wavenumber));
            SH{p+N}(:,1) = FrequencyRange(1:height(SH{p+N}));
            SH{p+N}(:,2) = FrequencyRange(1:height(SH{p+N}))/1e3;
            SH{p+N}(:,3) = FrequencyRange(1:height(SH{p+N}))*PlateThickness;
            SH{p+N}(end-length(c_complex)+1:end,8) = real(c_complex);
            SH{p+N}(end-length(c_complex)+1:end,9) = imag(c_complex);
            rw2 = 0;
        else
            X = [];
            Misses = 0;
            BelowCutoff = 0;
            for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
                if  Stop
                    return
                end
                AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                for m = 1:SuperLayerSize
                    rw2(m) = Material{m}.Density*AngularFrequency2;
                end
                X(i,1) = 0;
                Neighbors = [];
                for j = 1:length(SH)
                    if  i <= height(SH{j})
                        Neighbors(j,:) = [SH{j}(i,8) SH{j}(i,9)];
                    end
                end
                NeighborsNumber = height(Neighbors);
                for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                    if  all(X(:,1) == 0) 
                        SweepRangeReal = [1.1*PhaseVelocityLimit XSH0];
                    elseif isscalar(find(X(:,1) ~= 0))
                        SweepRangeReal = [PhaseVelocityLimit XSH0];
                    elseif numel(find(X(:,1) ~= 0)) == 2
                        SweepRangeReal = [1.1*X(end-1,3) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                    else
                        SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                    end
                    if  SweepRangeReal(1) == SweepRangeReal(2)
                        SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                    end
                    if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                        SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                        SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                    end
                    if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < XSH0
                        SweepRangeReal(2) = XSH0;
                    elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                        SweepRangeReal(2) = 0;
                    end
                    if  all(X(:,1) == 0)
                        SweepRangeImag = [-1000*T 0];
                    elseif isscalar(find(X(:,1) ~= 0))
                        SweepRangeImag = [4*SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; 
                    else
                        SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                    end
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-20*T 0];
                    elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                        SweepRangeImag(2) = 0;
                    end
                    if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                        SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                        SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                    end
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                    for o = 1:SearchAreaSections+1 % increase search resolution
                        if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                            break
                        end 
                        for k = 1:1e2 % search minimum in characteristic equation and converge upon it 
                            if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                                if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                    SweepRangeReal = SweepRangeReal(1):.25^o/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeReal = SweepRangeReal(1):.25/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                end
                            else
                                if  length(SweepRangeReal) == 2
                                    SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                end
                                if  length(SweepRangeImag) == 2
                                    SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                end
                            end
                            if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                                break
                            end
                            Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
                            Y = Computer(1,ModeFamily,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,c,rw2);
                            if  numel(find(X(:,1) ~= 0)) <= 3
                                for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                                    for j = 2:height(Wavenumber)-1
                                        for t = 1:NeighborsNumber
                                            if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1)
                                                Y(j,l) = NaN;
                                            end
                                        end
                                    end
                                end
                            else
                                for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                                    for j = 2:height(Wavenumber)-1
                                        for t = 1:NeighborsNumber
                                            if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1) && SweepRangeImag(l-1) < Neighbors(t,2) && SweepRangeImag(l+1) > Neighbors(t,2)
                                                Y(j,l) = NaN;
                                            end
                                        end
                                    end
                                end
                            end
                            if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
                                Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
                            end
% if  i>=ceil(HSShear(p)/FrequencyResolution)+10
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))
% else
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
% close(f)
% end
                            Min = zeros(size(Y));
                            for l = 2:size(Y,2)-1
                                for j = 2:size(Y,1)-1
                                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                        Min(j,l) = 1;
                                    end
                                end
                            end
                            [b1,b2] = find(Min);
                            if  ~isempty(b1) % one or multiple minima are found
                                if  isscalar(b1)
                                    MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                else
                                    delta = [];
                                    for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                        cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                                        delta(l) = abs(cp(l)-X(end-1,1));
                                    end
                                    if  all(X(:,1) == 0)
                                        [~,l] = max(delta);
                                    else
                                        [~,l] = min(delta);
                                    end
                                    MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                end
                            else
                                if  numel(find(X(:,1) ~= 0)) <= 3 && k <= 2
                                    if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                        SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                        SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                        if  SweepRangeReal(1) > 1.1*PhaseVelocityLimit
                                            SweepRangeReal(1) = 1.1*PhaseVelocityLimit;
                                        end
                                        if  SweepRangeReal(2) < XSH0
                                            SweepRangeReal(2) = XSH0;
                                        end
                                        if  SweepRangeImag(2) > 0
                                            SweepRangeImag(2) = 0;
                                        end
                                    end
                                    MIN = 0;
                                    break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                                else
                                    Min = zeros(size(Y)); % find border minima
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
                                        if  isscalar(b1)
                                            MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                        else
                                            Value = [];
                                            for l = 1:length(b1)
                                                Value(l) = Y(b1(l),b2(l));
                                            end
                                            [~,l] = min(Value);
                                            MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                        end
                                    else
                                        if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                            SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                            SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                            if  SweepRangeReal(2) < 0
                                                SweepRangeReal(2) = 0;
                                            end
                                            if  SweepRangeImag(2) > 0
                                                SweepRangeImag(2) = 0;
                                            end
                                        end
                                        MIN = 0;
                                        break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                                    end
                                end
                            end
                            if  k == 100 || (Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                                break
                            end
                            if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                                if  MIN(1) == 1
                                    if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                    SweepRangeReal = [SweepRangeReal(1)+4*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
                                elseif MIN(2) == 1
                                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                    end
                                    SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                                elseif MIN(1) == size(Y,1)
                                    if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                    SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                elseif MIN(2) == size(Y,2)
                                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                    end
                                    SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                end
                            else
                                if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                    if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                                        SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                    end
                                elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                                        SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                    end
                                    if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                end
                            end
                            if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(1) > 1.1*PhaseVelocityLimit
                                SweepRangeReal(1) = 1.1*PhaseVelocityLimit;
                            end
                            if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < XSH0
                                SweepRangeReal(2) = XSH0;
                            elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                                SweepRangeReal(2) = 0;
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  any(MIN)
                            if  (N > 0 && (numel(find(X(:,1) ~= 0)) <= 3 && i <= height(SH{p+N-1}) && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < SH{p+N-1}(i,4)*1e3)) ||...
                                (numel(find(X(:,1) ~= 0)) >= 1 && numel(find(X(:,1) ~= 0)) <= 3 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) > X(end-1,1)) ||...
                                SweepRangeReal(MIN(1)) == 0
                                Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                                NeighborsNumber = height(Neighbors);
                            else
                                if  numel(find(X(:,1) ~= 0)) <= 5
                                    Outlier = 0;
                                else
                                    z1 = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                                    z2 = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                                    if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)
                                        Outlier = 1;
                                        Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                                        NeighborsNumber = height(Neighbors);
                                    else
                                        Outlier = 0;
                                    end
                                end
                                if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                                    X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                                    X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                                    X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                                    X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                                    Misses(i) = 0;
                                    BelowCutoff(i) = 0;
                                    break
                                end
                            end
                            if  all(X(:,1) == 0) 
                                SweepRangeReal = [1.1*PhaseVelocityLimit XSH0];
                            elseif isscalar(find(X(:,1) ~= 0))
                                SweepRangeReal = [PhaseVelocityLimit XSH0];
                            elseif numel(find(X(:,1) ~= 0)) == 2
                                SweepRangeReal = [1.1*X(end-1,3) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                            else
                                SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                            end
                            if  SweepRangeReal(1) == SweepRangeReal(2)
                                SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                            end
                            if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                                SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                                SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                            end
                            if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < XSH0
                                SweepRangeReal(2) = XSH0;
                            elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                                SweepRangeReal(2) = 0;
                            end
                            if  all(X(:,1) == 0)
                                SweepRangeImag = [-1000*T 0];
                            elseif isscalar(find(X(:,1) ~= 0))
                                SweepRangeImag = [4*SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; 
                            else
                                SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                            end
                            if  all(SweepRangeImag == [0 0])
                                SweepRangeImag = [-20*T 0];
                            elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                                SweepRangeImag(2) = 0;
                            end
                            if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                                SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                                SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  numel(find(X(:,1) ~= 0)) > 20 && o == SearchAreaSections
                            break
                        end
                    end
                    if  X(i,1) > 0 % stop q-loop if minimum has been found
                        break
                    end
                end
                if  X(i,1) == 0 && any(X(:,1)) % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
                    if  isscalar(find(X(:,1) > 0))
                        X(i-1,:) = 0;
                        Misses(i) = 0;
                        BelowCutoff(i) = 1;
                    else
                        Smooth1 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
                        Smooth2 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
                        Smooth3 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
                        Smooth4 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
                        Fit1 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth1,'cubicspline');
                        Fit2 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth2,'cubicspline');
                        Fit3 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth3,'cubicspline');
                        Fit4 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth4,'cubicspline');
                        X(i,1) = Fit1(FrequencyRange(i));
                        X(i,2) = Fit2(FrequencyRange(i));
                        X(i,3) = Fit3(FrequencyRange(i));
                        X(i,4) = Fit4(FrequencyRange(i));
                        if  X(i,2) < 0 % negative attenuation is impossible
                            X(i,[2 4]) = 0;
                        end
                        Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
                        BelowCutoff(i) = 0;
                    end
                elseif all(X(:,1) == 0) % if we are scanning below the cut-off frequency of the damped mode
                    Misses(i) = 0;
                    BelowCutoff(i) = 1;
                end
                if  FrequencyResolution*length(BelowCutoff(BelowCutoff == 1)) > BelowCutoffWidth/PlateThickness/1e3 ||... % if we exceed the allowed scanning width (in kHz*mm) below the cut-off frequency of a damped mode without finding it; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one
                    length(Misses) == ceil(H(p)/FrequencyResolution)+6 && numel(find(Misses == 1)) >= 2
                    MissingModes(p+N) = 1;
                    break
                end
                if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
                    X(end-MissingSamples:end,:) = [];
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
            if  all(X(:,1) == 0)
                MissingModes(p+N) = 1;
            end
            if  ~MissingModes(p+N)
                [~,z] = max(X(:,1)); % find non-zero data below the cut-off frequency
                X(1:z-1,:) = 0; % remove them
                Misses(1:z-1) = 0; % remove also misses
                if  X(z,2) == 0 % if frequency cut-off is at zero damping, remove it also
                    X(z,:) = 0;
                end
                X(Misses(1:height(X)) == 1,:) = NaN;
                SH{p+N}(:,1) = FrequencyRange(1:height(X));
                SH{p+N}(:,2) = FrequencyRange(1:height(X))/1e3;
                SH{p+N}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
                SH{p+N}(:,4) = fillmissing(X(:,1),'spline')/1e3;
                SH{p+N}(:,7) = fillmissing(X(:,2),'spline');
                SH{p+N}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
                SH{p+N}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
                SH{p+N}(SH{p+N}(:,7) < 0,7) = 0; % negative attenuation is impossible
            else
                if  ModeFamily == 3 && p == 1
                    SH{1}(1,1) = 0;
                else
                    SH{p+N} = SH{p+N-1};
                end
                X1{p+N}(1) = 0;
                continue
            end
        end
        if  max(SH{p+N}(:,4))*1e3 > PhaseVelocityLimit
            X1{p+N}(1,1) = 0;
        else
            [Max,MaxInd] = max(SH{p+N}(:,8));
            PhaseVelocityRange = Max+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
            for i = 1:length(PhaseVelocityRange)
                if  Stop
                    return
                end
                X1{p+N}(i,1) = 0;
                for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                    if  i == 1
                        SweepRangeFrq = [SH{p+N}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 SH{p+N}(MaxInd,1)+4/PlateThickness/1e3];
                        SweepRangeImag = [SearchWidthImag(1)*SH{p+N}(MaxInd,9) SearchWidthImag(2)*SH{p+N}(MaxInd,9)];
                    else
                        SweepRangeFrq = [X1{p+N}(end-1,1)-FrequencyOffset/PlateThickness/1e3 X1{p+N}(end-1,1)+4/PlateThickness/1e3];
                        SweepRangeImag = [SearchWidthImag(1)*X1{p+N}(end-1,3) SearchWidthImag(2)*X1{p+N}(end-1,3)];
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
                            AngularFrequency2 = AngularFrequency'.^2;
                            Wavenumber = AngularFrequency'./(PhaseVelocityRange(i)+SweepRangeImag*1i);
                            if  SuperLayerSize == 1
                                Y = Computer_SingleLayer(ModeFamily,AngularFrequency2,Wavenumber,Material,LayerThicknesses,c);
                            else
                                Y = Computer(2,ModeFamily,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,c,rw2);
                            end
                            if  abs(SweepRangeImag(end)) < 1e-3 % in case the fluid has low density, this allows finding solutions at zero attenuation for A0
                                Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
                            end
% if  p==4
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeFrq,20*log10(Y))
% else
% f = figure;surf(SweepRangeImag,SweepRangeFrq,20*log10(Y))
% end
% close(f)
% end
                            Min = zeros(size(Y));
                            for l = 2:size(Y,2)-1
                                for j = 2:size(Y,1)-1
                                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                        Min(j,l) = 1;
                                    end
                                end
                            end
                            [b1,b2] = find(Min);
                            if  ~isempty(b1) % one or multiple minima are found
                                if  isscalar(b1)
                                    MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                else
                                    delta = [];
                                    for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                        frq(l) = AngularFrequency(b1(l))/(2*pi)/1e3;
                                        if  i == 1
                                            delta(l) = abs(frq(l)-SH{p+N}(MaxInd,1));
                                        else
                                            delta(l) = abs(frq(l)-X1{p+N}(end-1,1));
                                        end
                                    end
                                    [~,l] = min(delta);
                                    MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                end
                            else
                                Min = zeros(size(Y)); % find border minima
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
                                    if  isscalar(b1)
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
                            z = isoutlier(vertcat(X1{p+N}(1:end-1,1),AngularFrequency(MIN(1))/(2*pi)/1e3),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                            if  ~z(end) || all(X1{p+N}(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                                X1{p+N}(i,1) = AngularFrequency(MIN(1))/(2*pi)/1e3; % frequency (kHz)
                                X1{p+N}(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                                X1{p+N}(i,3) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                                X1{p+N}(i,4) = AngularFrequency(MIN(1))/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                                break
                            end
                            if  i == 1
                                SweepRangeFrq = [SH{p+N}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 SH{p+N}(MaxInd,1)+4/PlateThickness/1e3];
                                SweepRangeImag = [SearchWidthImag(1)*SH{p+N}(MaxInd,9) SearchWidthImag(2)*SH{p+N}(MaxInd,9)];
                            else
                                SweepRangeFrq = [X1{p+N}(end-1,1)-FrequencyOffset/PlateThickness/1e3 X1{p+N}(end-1,1)+4/PlateThickness/1e3];
                                SweepRangeImag = [SearchWidthImag(1)*X1{p+N}(end-1,3) SearchWidthImag(2)*X1{p+N}(end-1,3)];
                            end
                            if  all(SweepRangeImag == [0 0])
                                SweepRangeImag = [-20*T 0];
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                    end
                    if  X1{p+N}(i,1) > 0 % stop q-loop if minimum has been found
                        break
                    end
                end
                if  X1{p+N}(i,1) == 0 || X1{p+N}(i,4) > PhaseVelocityLimit
                    break
                end
                if  X1{p+N}(i,2) < 0 % negative attenuation is impossible
                    X1{p+N}(i,2:3) = 0;
                end
                if  SuperLayerSize > 1
                    if  Multithreading
                        send(Q2,[X1{p+N}(i,1),X1{p+N}(i,4)/1e3,p+N])
                    else
                        addpoints(g1(p+N),X1{p+N}(i,1),X1{p+N}(i,4)/1e3);
                        drawnow limitrate
                    end
                end
            end
            if  length(X1) < p+N
                X1{p+N}(1,1) = 0;
            end
            if  X1{p+N}(1,1) > 0
                X1{p+N}(:,7) = X1{p+N}(:,2);
                X1{p+N}(:,2) = X1{p+N}(:,1)/1e3;
                X1{p+N}(:,3) = X1{p+N}(:,1)*PlateThickness;
                X1{p+N}(:,4) = X1{p+N}(:,4)/1e3;
                X1{p+N}(X1{p+N}(:,1) == 0,:) = []; 
                X1{p+N} = flipud(X1{p+N});
            end
        end
    end
    X1(MissingModes == 1) = [];
    SH(MissingModes == 1) = [];
    for p = 1:N
        SH{p}(:,8:9) = [];
    end
    for p = N+1:length(SH)
        SH{p}(SH{p}(:,4) == 0,:) = [];
        SH{p}(:,8:9) = [];
        if  X1{p}(1,1) > 0
            SH{p} = vertcat(X1{p},SH{p});
        end
    end
end
end
function Y = Computer(Mode,ModeFamily,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,c,rw2)
    if  Mode == 2
        AngularFrequency2 = reshape(repmat(AngularFrequency2,width(Wavenumber),1),1,1,[]);
    end
    Wavenumber2 = reshape(Wavenumber,1,1,[]).^2;
    for m = 1:SuperLayerSize
        if  Mode == 1
            k3 = sqrt((rw2(m)-Wavenumber2*c{m}(6,6))/c{m}(4,4));
        elseif Mode == 2
            k3 = sqrt((Material{m}.Density*AngularFrequency2-Wavenumber2*c{m}(6,6))/c{m}(4,4));
        end
        D = k3*c{m}(4,4);
        G = k3*LayerThicknesses(m);
        CosG = cos(G);
        SinG = sin(G);
        L{m} = [CosG 1i*SinG./D;1i*SinG.*D CosG];
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M{1} = pagemtimes(M{1},L{m});
    end
    for m = 1:length(Pattern)
        M{m+1} = pagemtimes(M{m},M{Pattern(m)});
    end
    if  ModeFamily < 3
        Y = reshape(abs(M{end}(2,1,:)),size(Wavenumber));
    elseif ModeFamily == 3
        Y = reshape(abs(M{end}(2,2,:)),size(Wavenumber));
    end
end
function Y = Computer_SingleLayer(ModeFamily,AngularFrequency2,Wavenumber,Material,LayerThicknesses,c)
    AngularFrequency2 = reshape(repmat(AngularFrequency2,width(Wavenumber),1),1,1,[]);
    Wavenumber2 = reshape(Wavenumber,1,1,[]).^2;
    k3 = sqrt((Material{1}.Density*AngularFrequency2-Wavenumber2*c{1}(6,6))/c{1}(4,4));
    G = k3*LayerThicknesses/2;
    if  ModeFamily == 1
        Y = reshape(abs(sin(G)),size(Wavenumber));
    elseif ModeFamily == 3
        Y = reshape(abs(cos(G)),size(Wavenumber));
    end
end