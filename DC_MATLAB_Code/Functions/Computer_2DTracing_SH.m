% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function A = Computer_2DTracing_SH(Multithreading,Q,ax,Geometry,Material,Ro,Ri,FrequencyRange,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,h,LineColor,c,SuperLayerSize,LayerThicknesses,Pattern,XS0)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};Fit={};
if  isstruct(Material)
    MaterialType = 1;
    Velocity1 = Material.LongitudinalVelocity;
    cT2 = Material.TransverseVelocity_complex^2;
else
    MaterialType = 2;
    Velocity1 = 1.1*XS0(end);
end
if  strcmp(Geometry,'Plate')
    Geometry = 1;
    ModeFamily = Ri;
    ScaleFactor = Ro;
elseif strcmp(Geometry,'Rod')
    Geometry = 2;
    ScaleFactor = 2*Ro;
elseif strcmp(Geometry,'Pipe')
    Geometry = 3;
    ScaleFactor = Ro-Ri;
end
H = FrequencyRange(2);
for i = 1:2*Sweeps
    H(end+1,1) = FrequencyRange(find(FrequencyRange >= i/(2*Sweeps+1)*FrequencyRange(end),1));
end
H(end+1,1) = FrequencyRange(end-10);
H(:,6) = 1;
H = [h' zeros(length(h),5);H];
Neighbors0 = zeros(0,2);
if  ~Multithreading
    g = [animatedline(ax,'LineStyle','--','color',LineColor) animatedline(ax,'LineStyle','--','color',LineColor)];
end
if  MaterialType == 1
    AngularFrequency = 2*pi*FrequencyRange*1e3;
    Wavenumber = AngularFrequency/Material.TransverseVelocity_complex;
    SweepRangeComplex = AngularFrequency./Wavenumber;
    A{1}(:,[1:4 7:9]) = [FrequencyRange'*[1 1e-3 ScaleFactor] (AngularFrequency./real(Wavenumber))' imag(Wavenumber)' real(SweepRangeComplex)' imag(SweepRangeComplex)']; % add the dispersion curve to the output data container
    Fit{1,1} = fit(A{1}(:,1),A{1}(:,8),'cubicspline'); % fit the dispersion curve to query interpolated data to be ignored in the tracing of the next dispersion curves
    Fit{2,1} = fit(A{1}(:,1),A{1}(:,9),'cubicspline');
    if  Multithreading
        send(Q(1),A{1}(:,[1 4]))
    else
        line(ax,A{1}(:,1),A{1}(:,4)/1e3,'LineStyle','--','color',LineColor)
    end
end
p = 0;
while p < size(H,1) % step through the modes
    p = p+1;
    X = 0;
    Misses = false;
    BelowCutoff = false;
    if  H(p,2) || H(p,6)
        Frequency = H(p);
    else
        Frequency = FrequencyRange(ceil(H(p)/FrequencyResolution)+1);
    end
    Gear = 1;
    for r = 1:2 % 1: trace forward 2: trace backward
        if  r == 1
            Direction = 1;
            i = 0;
        else
            Direction = -1;
            X = flipud(X(~BelowCutoff,:));
            Misses = fliplr(Misses(~BelowCutoff));
            Frequency = fliplr(Frequency(~BelowCutoff));
            Gear = fliplr(Gear(~BelowCutoff));
            i = length(Frequency);
            Safety = i+5;
            if  Multithreading
                send(Q(r),[Frequency(i) X(i)/1e3])
            else
                addpoints(g(r),Frequency(i),X(i)/1e3);
                drawnow limitrate
            end
        end
        while (r == 1 && Frequency(end) < FrequencyRange(end)) || (r == 2 && Frequency(end) > FrequencyRange(1)) % step through the frequency range
            if  Stop
                return
            end
            i = i+1;
            X(i,4) = 0;
            Misses(i) = false;
            BelowCutoff(i) = false;
            if  i > 1 && ~any(X(:,1)) % determine the next frequency step depending on the slope of the dispersion curve; a higher slope needs a finer stepping; a finer stepping sets in only if "CriticalSlope" is exceeded or in some special cases
                Gear(i) = 1;
                Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution;
            elseif i > 1 && numel(find(X(:,1))) < 6
                Gear(i) = 3;
                Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution/10;
            elseif numel(find(X(:,1))) > 5
                Slope = abs((X(i-1)-X(i-2))/(Frequency(i-1)-Frequency(i-2)))*FrequencyResolution;
                for j = 1:length(GearSlopeRanges)-1
                    if  Slope >= GearSlopeRanges(j) && Slope < GearSlopeRanges(j+1)
                        Gear(i) = j;
                        break
                    end
                end
                if  Gear(i)-Gear(i-1) > 1
                    Gear(i) = Gear(i-1)+1;
                elseif Gear(i)-Gear(i-1) < -1
                    Gear(i) = Gear(i-1)-1;
                end
                if  r == 2 && i <= Safety && Gear(i) < 3
                    Gear(i) = 3;
                end
                if  Gear(i) <= length(GearFrequencyResolution)
                    Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution/GearFrequencyResolution(Gear(i));
                    if  Gear(i) == 1 && ~any(isapprox(FrequencyRange,Frequency(i),'tight'))
                        if  r == 1 && FrequencyRange(end) > Frequency(i)
                            Frequency(i) = FrequencyRange(find(FrequencyRange > Frequency(i-1),1));
                        elseif r == 2 && FrequencyRange(1) < Frequency(i)
                            Frequency(i) = FrequencyRange(find(FrequencyRange < Frequency(i-1),1,'last'));
                        end
                    end
                else
                    if  sqrt(CriticalSlope/Slope) < 1/GearFrequencyResolution(end)
                        Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution*sqrt(CriticalSlope/Slope);
                    else
                        Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution/GearFrequencyResolution(end);
                    end
                end
            end
            if  Frequency(i) > FrequencyRange(end)
                Frequency(i) = FrequencyRange(end);
            elseif Frequency(i) < FrequencyRange(1)
                Frequency(i) = FrequencyRange(1);
            end
            if  isapprox(Frequency(i),FrequencyRange(end),'tight')
                Frequency(i) = FrequencyRange(end);
            elseif isapprox(Frequency(i),FrequencyRange(1),'tight')
                Frequency(i) = FrequencyRange(1);
            end
            Neighbors = Neighbors0;
            for j = 1:length(A) % add previously traced solutions to be ignored at the current frequency
                if  Frequency(i) >= A{j}(1) && Frequency(i) <= A{j}(end,1)
                    [~,z] = min(abs(Frequency(i)-A{j}(:,1)));
                    if  A{j}(z) == Frequency(i)
                        Neighbors(end+1,:) = A{j}(z,8:9);
                    else
                        Neighbors(end+1,:) = [Fit{1,j}(Frequency(i)) Fit{2,j}(Frequency(i))]; % if we are in between frequency steps, we use interpolation
                    end
                end
            end
            AngularFrequency = 2*pi*Frequency(i)*1e3;
            AngularFrequency2 = AngularFrequency^2;
            if  MaterialType == 1
                kT2 = AngularFrequency2/cT2;
            else
                for m = 1:SuperLayerSize
                    rw2(m) = Material{m}.Density*AngularFrequency2;
                end
            end
            for q = 1:SearchAreaExtensions % q > 1 iterations are with extended search area (SweepRangeReal x SweepRangeImag)
                if  numel(find(X(:,1))) < 2
                    if  ~H(p,2) || isscalar(find(X(:,1)))
                        if  ~H(p,6)
                            if  any(X(:,1))
                                SweepRangeReal = X(i-1,3)+PhaseVelocityLimit*[.25 -.25];
                                SweepRangeImag = X(i-1,4)+PhaseVelocityLimit*[-.25 .25];
                            else
                                SweepRangeReal = [PhaseVelocityLimit 1];
                                SweepRangeImag = PhaseVelocityLimit*[-1 1];
                            end
                        else
                            if  any(X(:,1))
                                SweepRangeReal = X(i-1,3)+Velocity1*[.25 -.25];
                                SweepRangeImag = X(i-1,4)+Velocity1*[-.25 .25];
                            else
                                SweepRangeReal = [Velocity1 1];
                                SweepRangeImag = Velocity1*[-1 1];
                            end
                        end
                    else
                        SweepRangeReal = H(p,2)+H(p,4)*[2 -2];
                        SweepRangeImag = H(p,3)+H(p,5)*[-2 2];
                    end
                else
                    dx = X(i-1,3:4)-X(i-2,3:4);
                    df = abs(diff(Frequency(i-2:i)));
                    SweepRangeReal = X(i-1,3)+dx(1)/df(1)*df(2)+abs(dx(1))*SearchWidth;
                    SweepRangeImag = X(i-1,4)+dx(2)/df(1)*df(2)+abs(dx(2))*-SearchWidth;
                end
                if  SweepRangeReal(1)-SweepRangeReal(2) < Accuracy
                    SweepRangeReal = SweepRangeReal+Accuracy*[1 -1];
                end
                if  SweepRangeImag(2)-SweepRangeImag(1) < Accuracy
                    SweepRangeImag = SweepRangeImag+Accuracy*[-1 1];
                end
                if  SweepRangeReal(2) < 0
                    SweepRangeReal(2) = 0;
                end
                SweepRangeReal0 = SweepRangeReal;
                SweepRangeImag0 = SweepRangeImag;
                for o = 1:SearchAreaSections % increase search resolution
                    for k = 1:100 % search minimum in characteristic equation and converge upon it (bisections)
                        if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                            if  (~any(X(:,1)) && ~H(p,2)) || isscalar(find(X(:,1)))
                                if  any(X(:,1))
                                    if  ~H(p,6)
                                        SweepRangeReal = SweepRangeReal(1):-PhaseVelocityLimit/SweepSections:SweepRangeReal(end);
                                        SweepRangeImag = SweepRangeImag(1):PhaseVelocityLimit/SweepSections:SweepRangeImag(end);
                                    else
                                        SweepRangeReal = SweepRangeReal(1):-Velocity1/SweepSections:SweepRangeReal(end);
                                        SweepRangeImag = SweepRangeImag(1):Velocity1/SweepSections:SweepRangeImag(end);
                                    end
                                else
                                    SweepRangeReal = SweepRangeReal(1):(SweepRangeReal(end)-SweepRangeReal(1))/SweepSections:SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):(SweepRangeImag(end)-SweepRangeImag(1))/SweepSections/2:SweepRangeImag(end);
                                end
                            else
                                NeighborsInside = false;
                                for j = 1:size(Neighbors,1)
                                    if  Neighbors(j,1) < SweepRangeReal(1) && Neighbors(j,1) > SweepRangeReal(end) && Neighbors(j,2) > SweepRangeImag(1) && Neighbors(j,2) < SweepRangeImag(end)
                                        NeighborsInside = true;
                                        break
                                    end
                                end
                                if  q == 1 && NeighborsInside
                                    Quotient = 10;
                                else
                                    Quotient = 1;
                                end
                                if  SweepRangeReal(1)-SweepRangeReal(end) > 100*(SweepRangeImag(end)-SweepRangeImag(1)) % set the new search area around the found minimum
                                    SweepRangeReal = SweepRangeReal(1):.25^o/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25/q/Quotient*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                elseif SweepRangeReal(1)-SweepRangeReal(end) <= 100*(SweepRangeImag(end)-SweepRangeImag(1)) && SweepRangeReal(1)-SweepRangeReal(end) >= .01*(SweepRangeImag(end)-SweepRangeImag(1))
                                    SweepRangeReal = SweepRangeReal(1):.25/q/o/Quotient*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25/q/o/Quotient*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                elseif SweepRangeReal(1)-SweepRangeReal(end) < .01*(SweepRangeImag(end)-SweepRangeImag(1))
                                    SweepRangeReal = SweepRangeReal(1):.25/q/Quotient*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                end
                            end
                        else
                            if  length(SweepRangeReal) == 2
                                SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                            end
                            if  length(SweepRangeImag) == 2
                                SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                            end
                        end
                        SweepRangeRealRange = SweepRangeReal(1)-SweepRangeReal(end);
                        SweepRangeImagRange = SweepRangeImag(end)-SweepRangeImag(1);
                        SweepRangeRealStep = SweepRangeReal(1)-SweepRangeReal(2);
                        SweepRangeImagStep = SweepRangeImag(2)-SweepRangeImag(1);
                        Wavenumber = AngularFrequency./(SweepRangeReal'+1i*SweepRangeImag);
                        YSize = size(Wavenumber);
                        if  k == 1 && ~any(X(:,1)) && ~H(p,2)
                            Wavenumber(abs(imag(Wavenumber)) > AttenuationLimit*real(Wavenumber)) = NaN; % exclude strongly damped modes
                        end
                        for j = 1:size(Neighbors,1) % ignore previously traced solutions
                            Wavenumber(SweepRangeReal > Neighbors(j,1)-SweepRangeRealStep/2 & SweepRangeReal < Neighbors(j,1)+SweepRangeRealStep/2,SweepRangeImag > Neighbors(j,2)-SweepRangeImagStep/2 & SweepRangeImag < Neighbors(j,2)+SweepRangeImagStep/2) = NaN;
                        end
                        if  Geometry == 1
                            for m = 1:SuperLayerSize
                                k3 = sqrt((rw2(m)-reshape(Wavenumber,1,1,[]).^2*c{m}(6,6))/c{m}(4,4));
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
                            if  ModeFamily == 2
                                Yc = reshape(M{end}(2,2,:),YSize);
                            else
                                Yc = reshape(M{end}(2,1,:),YSize);
                            end
                        elseif Geometry == 2
                            yRo = sqrt(kT2-Wavenumber.^2)*Ro;
                            Yc = yRo-2*besselj(1,yRo)./besselj(0,yRo);
                        elseif Geometry == 3
                            y = sqrt(Wavenumber.^2-kT2);
                            yRo = y*Ro;
                            yRi = y*Ri;
                            Yc = besseli(2,yRi).*besselk(2,yRo)-besseli(2,yRo).*besselk(2,yRi);
                        end
                        Y = abs(Yc);
                        Y(isinf(Y)) = NaN;
% if numel(find(X(:,1))) < 3 && k < 3
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y));rotate3d,view(-60,60)
% 1
% close(f)
% end
                        Min = [];
                        if  k == 1 && ((~any(X(:,1)) && ~H(p,2)) || isscalar(find(X(:,1))))
                            for l = 2:YSize(2)-1
                                for j = 2:YSize(1)-1
                                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1) && Y(j,l) < Y(j+1,l-1) && Y(j,l) < Y(j+1,l+1) && Y(j,l) < Y(j-1,l-1) && Y(j,l) < Y(j-1,l+1)
                                        Min(end+1,:) = [j l];
                                    end
                                end
                            end
                        else
                            for l = 2:YSize(2)-1
                                for j = 2:YSize(1)-1
                                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                        Min(end+1,:) = [j l];
                                    end
                                end
                            end
                        end
                        if  any(Min) % one or multiple minima are found
                            if  size(Min,1) == 1
                                MIN = Min; % SweepRangeReal-index, SweepRangeImag-index
                            else
                                cp = zeros(size(Min,1),1);
                                for l = 1:length(cp)
                                    if  i == 1
                                        cp(l) = AngularFrequency/real(Wavenumber(Min(l,1),Min(l,2)));
                                    else
                                        cp(l) = abs(AngularFrequency/real(Wavenumber(Min(l,1),Min(l,2)))-X(i-1));
                                    end
                                end
                                if  i == 1
                                    [~,l] = min(cp);
                                else
                                    dm = sign(X(i-1,4)) == sign(SweepRangeImag(Min(:,2)));
                                    if  all(dm) || all(~dm)
                                        [~,l] = min(cp);
                                    else
                                        l = find(dm);
                                        if  length(l) > 1
                                            cp(~dm) = Inf;
                                            [~,l] = min(cp);
                                        end
                                    end
                                end
                                MIN = Min(l,:);
                                if  k == 1 && ~any(X(:,1)) && ~H(p,2)
                                    for j = 1:size(Min,1)
                                        if  j ~= l
                                            H = [H(1:p,:);Frequency(i) SweepRangeReal(Min(j,1)) SweepRangeImag(Min(j,2)) SweepRangeRealStep SweepRangeImagStep H(p,6);H(p+1:end,:)];
                                        end
                                    end
                                end
                            end
                        else
                            if  k == 1 && ((~any(X(:,1)) && ~H(p,2)) || isscalar(find(X(:,1))))
                                MIN = [];
                                break % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                            else
                                for j = 2:YSize(1)-1 % find border minima
                                    if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                        Min(end+1,:) = [j 1];
                                    end
                                    if  Y(j,YSize(2)) < Y(j-1,YSize(2)) && Y(j,YSize(2)) < Y(j+1,YSize(2)) && Y(j,YSize(2)) < Y(j,YSize(2)-1)
                                        Min(end+1,:) = [j YSize(2)];
                                    end
                                end
                                for j = 2:YSize(2)-1
                                    if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                        Min(end+1,:) = [1 j];
                                    end
                                    if  Y(YSize(1),j) < Y(YSize(1),j-1) && Y(YSize(1),j) < Y(YSize(1),j+1) && Y(YSize(1),j) < Y(YSize(1)-1,j) && SweepRangeReal(end) > 0 % ignore solutions with zero real part
                                        Min(end+1,:) = [YSize(1) j];
                                    end
                                end
                                if  any(Min) % one or multiple border minima are found
                                    if  size(Min,1) == 1
                                        MIN = Min;
                                    else
                                        Value = zeros(size(Min,1),1);
                                        for l = 1:length(Value)
                                            Value(l) = Y(Min(l,1),Min(l,2));
                                        end
                                        [~,l] = min(Value);
                                        MIN = Min(l,:);
                                    end
                                else
                                    if  q > 1 % extend search area (SweepRangeReal x SweepRangeImag)
                                        SweepRangeReal = (q-1)*SweepRangeRealRange*[.5 -.5]+[SweepRangeReal(1) SweepRangeReal(end)];
                                        SweepRangeImag = (q-1)*SweepRangeImagRange*[-.5 .5]+[SweepRangeImag(1) SweepRangeImag(end)];
                                        if  SweepRangeReal(2) < 0
                                            SweepRangeReal(2) = 0;
                                        end
                                    end
                                    MIN = [];
                                    break % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                                end
                            end
                        end
                        if  k == 100 || (Accuracy > SweepRangeRealStep && Accuracy > SweepRangeImagStep)
                            break
                        end
                        if  MIN(1) == 1 % set the new search area around the found minimum
                            if  SweepRangeImagRange > 10*SweepRangeRealRange
                                SweepRangeImag =  SweepRangeImagStep*[-1 1]+SweepRangeImag(MIN(2));
                            end
                            SweepRangeReal = [SweepRangeReal(1)+4*SweepRangeRealRange SweepRangeReal(end)];
                        elseif MIN(2) == 1
                            if  SweepRangeRealRange > 10*SweepRangeImagRange
                                SweepRangeReal = SweepRangeRealStep*[1 -1]+SweepRangeReal(MIN(1));
                            end
                            SweepRangeImag = [SweepRangeImag(1)-4*SweepRangeImagRange SweepRangeImag(end)];
                        elseif MIN(1) == YSize(1)
                            if  SweepRangeImagRange > 10*SweepRangeRealRange
                                SweepRangeImag =  SweepRangeImagStep*[-1 1]+SweepRangeImag(MIN(2));
                            end
                            SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*SweepRangeRealRange];
                        elseif MIN(2) == YSize(2)
                            if  SweepRangeRealRange > 10*SweepRangeImagRange
                                SweepRangeReal = SweepRangeRealStep*[1 -1]+SweepRangeReal(MIN(1));
                            end
                            SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*SweepRangeImagRange];
                        else
                            if  SweepRangeRealRange > 10*SweepRangeImagRange
                                if  Accuracy < SweepRangeRealStep
                                    SweepRangeReal = SweepRangeRealStep*[1 -1]+SweepRangeReal(MIN(1));
                                end
                            elseif SweepRangeRealRange <= 10*SweepRangeImagRange && SweepRangeRealRange >= .1*SweepRangeImagRange
                                if  Accuracy < SweepRangeRealStep
                                    SweepRangeReal = SweepRangeRealStep*[1 -1]+SweepRangeReal(MIN(1));
                                end
                                if  Accuracy < SweepRangeImagStep
                                    SweepRangeImag =  SweepRangeImagStep*[-1 1]+SweepRangeImag(MIN(2));
                                end
                            elseif SweepRangeRealRange < .1*SweepRangeImagRange
                                if  Accuracy < SweepRangeImagStep
                                    SweepRangeImag =  SweepRangeImagStep*[-1 1]+SweepRangeImag(MIN(2));
                                end
                            end
                        end
                        if  SweepRangeReal(2) < 0
                            SweepRangeReal(2) = 0;
                        end
                    end
                    if  any(MIN)
                        cp = AngularFrequency/real(Wavenumber(MIN(1),MIN(2)));
                        if  numel(find(X(:,1))) > 3 && ((cp < Velocity1 && X(i-1) > 3*Velocity1) || (cp > 3*Velocity1 && X(i-1) < Velocity1)) % determine if the solution is an outlier; if yes, ignore it and continue searching 
                            Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))]; % add the outlier to the solutions to be ignored
                            SweepRangeReal = SweepRangeReal0; % reset to initial search range
                            SweepRangeImag = SweepRangeImag0;
% disp(['p = ',num2str(p),' r = ',num2str(r),' f = ',num2str(Frequency(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Outlier'])
                        else
                            if  MaterialType == 1 && size(Min,1) == 1 && MIN(1) > 1 && MIN(2) > 1 && MIN(1) < YSize(1) && MIN(2) < YSize(2) && all(abs(diff(angle([Yc(MIN(1)-1,MIN(2)) Yc(MIN(1),MIN(2)-1);Yc(MIN(1)+1,MIN(2)) Yc(MIN(1),MIN(2)+1)]))) < pi/2) % check that there is a sufficient change in phase; only then is the minimum a zero and thereby a valid modal solution
                                Misses(i) = true;
                            end
                            X(i,:) = [cp imag(Wavenumber(MIN(1),MIN(2))) SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                            break
                        end
                    end
                    if  numel(find(X(:,1))) < 2
                        break
                    end
                end
                if  X(i) || numel(find(X(:,1))) < 2 % stop q-loop if minimum has been found
                    break
                end
            end
            if  ~X(i) && any(X(:,1)) % fit phase velocity, attenuation, real part, and imaginary part velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
                if  isscalar(find(X(:,1)))
                    X(i-1,:) = 0;
                    BelowCutoff(i-1:i) = true;
                else
                    Idx = 1:i-1;
                    if  r == 1
                        Idx(BelowCutoff(Idx)) = [];
                    end
                    if  length(Idx) > 50
                        Idx(1:end-50) = [];
                    end
                    Fit1 = fit(Frequency(Idx)',X(Idx,1),'cubicspline');
                    Fit2 = fit(Frequency(Idx)',X(Idx,2),'cubicspline');
                    Fit3 = fit(Frequency(Idx)',X(Idx,3),'cubicspline');
                    Fit4 = fit(Frequency(Idx)',X(Idx,4),'cubicspline');
                    X(i,:) = [Fit1(Frequency(i)) Fit2(Frequency(i)) Fit3(Frequency(i)) Fit4(Frequency(i))];
                    if  X(i,3) < 0
                        X(i,3) = 0;
                    end
                    Misses(i) = true; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
% figure,line(Frequency(Idx),X(Idx,1),'LineWidth',4,'color','r'),line(Frequency([Idx i]),Fit1(Frequency([Idx i])),'LineWidth',1.5,'color','g')
% figure,line(Frequency(Idx),X(Idx,2),'LineWidth',4,'color','r'),line(Frequency([Idx i]),Fit2(Frequency([Idx i])),'LineWidth',1.5,'color','g')
% figure,line(Frequency(Idx),X(Idx,3),'LineWidth',4,'color','r'),line(Frequency([Idx i]),Fit3(Frequency([Idx i])),'LineWidth',1.5,'color','g')
% figure,line(Frequency(Idx),X(Idx,4),'LineWidth',4,'color','r'),line(Frequency([Idx i]),Fit4(Frequency([Idx i])),'LineWidth',1.5,'color','g')
                end
            elseif ~any(X(:,1)) % if we are scanning below the cut-off frequency
                BelowCutoff(i) = true;
            end
            if  FrequencyResolution*length(BelowCutoff(BelowCutoff)) > BelowCutoffWidth/ScaleFactor/1e3 ||... % if we exceed the allowed scanning width (in kHz*mm) below the cut-off frequency without finding it
                (H(p,6) && BelowCutoff(i)) ||...
                X(i) > PhaseVelocityLimit2
                break
            elseif i > MissingSamples && all(Misses(i-MissingSamples+1:i)) % if more than 'MissingSamples' sample points are missing, the tracing stops
                X(i-MissingSamples+1:i,:) = [];
                Misses(i-MissingSamples+1:i) = [];
                Frequency(i-MissingSamples+1:i) = [];
                Gear(i-MissingSamples+1:i) = [];
                BelowCutoff(i-MissingSamples+1:i) = [];
                break
            end
            if  ~BelowCutoff(i) && ~Misses(i)
                if  Multithreading
                    send(Q(r),[Frequency(i) X(i)/1e3])
                else
                    addpoints(g(r),Frequency(i),X(i)/1e3);
                    drawnow limitrate
                end
            end
% String = ['p = ',num2str(p),' r = ',num2str(r),' f = ',num2str(Frequency(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i)
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i)
%     disp(['p = ',num2str(p),' r = ',num2str(r),' f = ',num2str(Frequency(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
% end
        end
        if  ~any(X(:,1)) || all(Misses)
            break
        end
    end
    if  any(X(:,1)) && (size(X,1) >= 20 || (size(X,1) > 1 && size(X,1) < 20 && ~any(Misses)))
        X(Misses,1:4) = NaN;
        if  r == 2
            X = flipud(X);
            Frequency = fliplr(Frequency);
        end
        A{end+1}(:,[1:4 7:9]) = [Frequency'*[1 1e-3 ScaleFactor] fillmissing(X(:,1:4),'spline')]; % add the dispersion curve to the output data container
        Fit{1,end+1} = fit(A{end}(:,1),A{end}(:,8),'cubicspline'); % fit the dispersion curve to query interpolated data to be ignored in the tracing of the next dispersion curves
        Fit{2,end} = fit(A{end}(:,1),A{end}(:,9),'cubicspline');
        if  Multithreading
            send(Q(r),A{end}(:,[1 4]))
        else
            line(ax,A{end}(:,1),A{end}(:,4)/1e3,'LineStyle','--','color',LineColor)
            clearpoints(g(1))
            clearpoints(g(2))
        end
    end
end
for p = 1:length(A) % post processing
    A{p}(A{p}(:,4) < 0,:) = [];
    A{p}(:,4) = A{p}(:,4)/1e3;
    A{p}(:,8:9) = [];
end