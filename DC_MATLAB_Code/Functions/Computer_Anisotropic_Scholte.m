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
function Scholte = Computer_Anisotropic_Scholte(Multithreading,Q1,Q2,ax,ModeFamily,Decoupled,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRange,PlateThickness,HigherOrderModes,FAScholte,H,H2,FrequencyResolution,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Symmetric,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth)        
TracingMode = 3; % 1: forward tracing only 2: backward tracing only 3: forward and backward tracing
BelowCutoffWidth = 5*BelowCutoffWidth;
SweepRangeImagStart = .1; % x*SlowFluidVelocity;
PhaseVelocityOffset = 1; % (m/s)
x = 20*FrequencyResolution/5; % limit up to which a 5 times smaller frequency step is used (kHz/mm)

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
Scholte{1} = [];
Scholte1{1} = [];
Scholte2{1} = [];
for m = 1:SuperLayerSize
    if  ~Decoupled
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
    else
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
    end
end
if  TracingMode ~=2
    if  ModeFamily == 1
        N = 0; % S
    elseif ModeFamily ~= 1
        N = 1; % A/B
    end
    if  ~Multithreading
        for p = 1:length(H)+N+1
            if  ModeFamily == 1
                g(p) = animatedline(ax,'LineStyle','-.','color','r');
            elseif ModeFamily == 2
                g(p) = animatedline(ax,'LineStyle','-.','color','b');
            elseif ModeFamily == 0
                g(p) = animatedline(ax,'LineStyle','-.','color',[.5 0 1]);
            end
        end
    end
    if  Symmetric
        SlowFluidVelocity = UpperFluid.Velocity;
    else
        if  ToggleUpperFluid && ToggleLowerFluid
            if  UpperFluid.Velocity > LowerFluid.Velocity
                SlowFluidVelocity = LowerFluid.Velocity;
            else
                SlowFluidVelocity = UpperFluid.Velocity;
            end
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            SlowFluidVelocity = UpperFluid.Velocity;
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            SlowFluidVelocity = LowerFluid.Velocity;
        end
    end
    if  ModeFamily ~= 1
        X = FAScholte(1,:);
        if  Multithreading
            send(Q1,[FrequencyRange(1),X(1)/1e3,1])
        else
            addpoints(g(1),FrequencyRange(1),X(1)/1e3);
            drawnow limitrate
        end
        for i = 2:length(FrequencyRange)
            if  Stop
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            AngularFrequency2 = AngularFrequency^2;
            kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
            kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
            gUpperFluid = 1i*UpperFluid.Density*AngularFrequency2;
            gLowerFluid = 1i*LowerFluid.Density*AngularFrequency2;
            for m = 1:SuperLayerSize
                rw2(m) = Material{m}.Density*AngularFrequency2;
                r2w4(m) = rw2(m)^2;
                if  ~Decoupled
                    r3w6(m) = rw2(m)^3;
                    b12(m) = a12(m)*rw2(m);
                    b22(m) = a22(m)*rw2(m);
                    b23(m) = a23(m)*r2w4(m);
                    b32(m) = a32(m)*rw2(m);
                    b33(m) = a33(m)*r2w4(m);
                    b34(m) = a34(m)*r3w6(m);
                else
                    b22(m) = a22(m)*rw2(m);
                    b32(m) = a32(m)*rw2(m);
                    b33(m) = r2w4(m);
                end
            end
            X(i,1) = 0;
            Neighbors = [];
            NeighborsNumber = 0;
            for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                if  i == 2
                    SweepRangeReal = [5*X(end-1,3) X(end-1,3)];
                else
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                end
                if  SweepRangeReal(1) == SweepRangeReal(2)
                    if  SweepRangeReal(2) > .99*SlowFluidVelocity
                        SweepRangeReal = [1.001*SlowFluidVelocity .99*SlowFluidVelocity];
                    else
                        SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                    end
                end
                if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                    SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                    SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                end
                if  SweepRangeReal(1) > 1.001*SlowFluidVelocity
                    SweepRangeReal(1) = 1.001*SlowFluidVelocity;
                end
                if  SweepRangeReal(2) < 0
                    SweepRangeReal(2) = 0;
                end
                if  i == 2
                    SweepRangeImag = [5*X(end-1,4) X(end-1,4)];
                else
                    SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                    if  -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                        SweepRangeImag(2) = -1e-10;
                    end
                end
                if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                    SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                    SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                end
                if  SweepRangeImag(2) > -1e-10
                    SweepRangeImag(2) = -1e-10;
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
                        if  ~Decoupled
                            Y = Computer_Coupled(1,ModeFamily,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,X,SlowFluidVelocity,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
                        else
                            Y = Computer_Decoupled(1,ModeFamily,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,X,SlowFluidVelocity,c,rw2,A1,a21,a31,b22,b32,b33);
                        end
                        if  ~isempty(Neighbors) % remove solutions of previously found lower modes
                            if  Viscoelastic
                                for l = 2:width(Wavenumber)-1
                                    for j = 2:height(Wavenumber)-1
                                        for t = 1:NeighborsNumber
                                            if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1) && SweepRangeImag(l-1) < Neighbors(t,2) && SweepRangeImag(l+1) > Neighbors(t,2)
                                                Y(j,l) = NaN;
                                            end
                                        end
                                    end
                                end
                            else
                                for j = 2:height(Wavenumber)-1
                                    for t = 1:NeighborsNumber
                                        if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1)
                                            Y(j,:) = NaN;
                                        end
                                    end
                                end
                            end
                        end
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
                                for l = 1:length(b1) % calculate which minimum lies closest to previous solution
                                    cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                                    delta(l) = abs(cp(l)-X(end-1,1));
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
                                    if  SweepRangeReal(1) > 1.001*SlowFluidVelocity
                                        SweepRangeReal(1) = 1.001*SlowFluidVelocity;
                                    end
                                    if  SweepRangeReal(2) < 0
                                        SweepRangeReal(2) = 0;
                                    end
                                    if  SweepRangeImag(2) > -1e-10
                                        SweepRangeImag(2) = -1e-10;
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
                        if  SweepRangeReal(1) > 1.001*SlowFluidVelocity
                            SweepRangeReal(1) = 1.001*SlowFluidVelocity;
                        end
                        if  SweepRangeReal(2) < 0
                            SweepRangeReal(2) = 0;
                        end
                        if  SweepRangeImag(2) > -1e-10
                            SweepRangeImag(2) = -1e-10;
                        end
                    end
                    if  any(MIN)
                        if  Viscoelastic && SweepRangeImag(MIN(2)) == -1e-10
                            Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                            NeighborsNumber = height(Neighbors);
                        else
                            if  i < 4
                                Outlier = 0;
                            else
                                z1 = isoutlier(vertcat(X(1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                                z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
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
                                break
                            end
                        end
                        if  i == 2
                            SweepRangeReal = [5*X(end-1,3) X(end-1,3)];
                        else
                            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                        end
                        if  SweepRangeReal(1) == SweepRangeReal(2)
                            if  SweepRangeReal(2) > .99*SlowFluidVelocity
                                SweepRangeReal = [1.001*SlowFluidVelocity .99*SlowFluidVelocity];
                            else
                                SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                            end
                        end
                        if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                            SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                            SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                        end
                        if  SweepRangeReal(1) > 1.001*SlowFluidVelocity
                            SweepRangeReal(1) = 1.001*SlowFluidVelocity;
                        end
                        if  SweepRangeReal(2) < -1e-10
                            SweepRangeReal(2) = -1e-10;
                        end
                        if  i == 2
                            SweepRangeImag = [5*X(end-1,4) X(end-1,4)];
                        else
                            SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                            if  -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                                SweepRangeImag(2) = -1e-10;
                            end
                        end
                        if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                            SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                            SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                        end
                        if  SweepRangeImag(2) > -1e-10
                            SweepRangeImag(2) = -1e-10;
                        end
                    end
                end
                if  X(i,1) > 0 % stop q-loop if minimum has been found
                    break
                end
            end
            if  X(i,1) == 0 % fit phase velocity, attenuation, and imaginary velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
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
        Scholte1{1}(:,1) = FrequencyRange(1:height(X));
        Scholte1{1}(:,2) = FrequencyRange(1:height(X))/1e3;
        Scholte1{1}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
        Scholte1{1}(:,4) = fillmissing(X(:,1),'spline')/1e3;
        Scholte1{1}(:,7) = fillmissing(X(:,2),'spline');
        Scholte1{1}(:,8) = fillmissing(X(:,3),'spline'); % real velocity (m/s)
        Scholte1{1}(:,9) = fillmissing(X(:,4),'spline'); % imaginary velocity (m/s)
        Scholte1{1}(Scholte1{1}(:,7) < 0,7) = 0; % negative attenuation is impossible
    end
    if  ToggleUpperFluid && ToggleLowerFluid && strcmp(UpperFluid.Name,LowerFluid.Name)
        H = [FrequencyRange(1) H];
    end
    if  HigherOrderModes && any(H)
        MissingModes(length(H)+N) = 0;
        for p = 1:length(H)
            X = [];
            Misses = 0;
            BelowCutoff = 0;
            for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
                if  Stop
                    return
                end
                AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
                kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
                gUpperFluid = 1i*UpperFluid.Density*AngularFrequency2;
                gLowerFluid = 1i*LowerFluid.Density*AngularFrequency2;
                for m = 1:SuperLayerSize
                    rw2(m) = Material{m}.Density*AngularFrequency2;
                    r2w4(m) = rw2(m)^2;
                    if  ~Decoupled
                        r3w6(m) = rw2(m)^3;
                        b12(m) = a12(m)*rw2(m);
                        b22(m) = a22(m)*rw2(m);
                        b23(m) = a23(m)*r2w4(m);
                        b32(m) = a32(m)*rw2(m);
                        b33(m) = a33(m)*r2w4(m);
                        b34(m) = a34(m)*r3w6(m);
                    else
                        b22(m) = a22(m)*rw2(m);
                        b32(m) = a32(m)*rw2(m);
                        b33(m) = r2w4(m);
                    end
                end
                X(i,1) = 0;
                Neighbors = [];
                for j = 1:length(Scholte1)
                    if  i <= height(Scholte1{j})
                        Neighbors(end+1,:) = [Scholte1{j}(i,8) Scholte1{j}(i,9)];
                    end
                end
                NeighborsNumber = height(Neighbors);
                for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                    if  numel(find(X(:,1) ~= 0)) < 4
                        SweepRangeReal = [1.001*SlowFluidVelocity .95*SlowFluidVelocity];
                    else
                        SweepRangeReal = [X(end-1,3)+.1*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                    end
                    if  SweepRangeReal(1) == SweepRangeReal(2)
                        if  numel(find(X(:,1) ~= 0)) < 4
                            SweepRangeReal = [1.001*SlowFluidVelocity .95*SlowFluidVelocity];
                        else
                            SweepRangeReal = [X(end-1,3)+.1*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                        end
                    end
                    if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                        SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                        SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                    end
                    if  SweepRangeReal(1) > 1.001*SlowFluidVelocity
                        SweepRangeReal(1) = 1.001*SlowFluidVelocity;
                    end
                    if  SweepRangeReal(2) < 0
                        SweepRangeReal(2) = 0;
                    end
                    if  numel(find(X(:,1) ~= 0)) < 4
                        SweepRangeImag = [-SweepRangeImagStart*SlowFluidVelocity -1e-10];
                    else
                        SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary velocity corresponds to the attenuation
                        if  -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                            SweepRangeImag(2) = -1e-10;
                        end
                    end
                    if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                        SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                        SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                    end
                    if  SweepRangeImag(2) > -1e-10
                        SweepRangeImag(2) = -1e-10;
                    end
                    for o = 1:SearchAreaSections % increase search resolution
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
                            if  ~Decoupled
                                Y = Computer_Coupled(1,ModeFamily,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,X,SlowFluidVelocity,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
                            else
                                Y = Computer_Decoupled(1,ModeFamily,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,X,SlowFluidVelocity,c,rw2,A1,a21,a31,b22,b32,b33);
                            end
                            if  Viscoelastic % remove solutions of previously found lower modes
                                if  numel(find(X(:,1) ~= 0)) < 4
                                    for l = 2:width(Wavenumber)-1
                                        for j = 2:height(Wavenumber)-1
                                            for t = 1:NeighborsNumber
                                                if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1)    
                                                    Y(j,l) = NaN;
                                                end
                                            end
                                        end
                                    end
                                else
                                    for l = 2:width(Wavenumber)-1
                                        for j = 2:height(Wavenumber)-1
                                            for t = 1:NeighborsNumber
                                                if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1) && SweepRangeImag(l-1) < Neighbors(t,2) && SweepRangeImag(l+1) > Neighbors(t,2)
                                                    Y(j,l) = NaN;
                                                end
                                            end
                                        end
                                    end
                                end
                            else
                                for j = 2:height(Wavenumber)-1
                                    for t = 1:NeighborsNumber
                                        if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1)
                                            Y(j,:) = NaN;
                                        end
                                    end
                                end
                            end
                            if  abs(SweepRangeImag(end)) < 1e-3
                                Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
                            end
% if  p==4&& i >= 235%q==1 && o == 1&&k==1
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
                                if  numel(find(X(:,1) ~= 0)) < 4 && k <= 2
                                    if  q > 1
                                        SweepRangeReal = [1.001*SlowFluidVelocity .95*SlowFluidVelocity];
                                        SweepRangeImag = [-SweepRangeImagStart*SlowFluidVelocity -1e-10];
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
                                            if  SweepRangeReal(1) > 1.001*SlowFluidVelocity
                                                SweepRangeReal(1) = 1.001*SlowFluidVelocity;
                                            end
                                            if  SweepRangeReal(2) < 0
                                                SweepRangeReal(2) = 0;
                                            end
                                            if  SweepRangeImag(2) > -1e-10
                                                SweepRangeImag(2) = -1e-10;
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
                            if  SweepRangeReal(1) > 1.001*SlowFluidVelocity
                                SweepRangeReal(1) = 1.001*SlowFluidVelocity;
                            end
                            if  SweepRangeReal(2) < 0
                                SweepRangeReal(2) = 0;
                            end
                            if  SweepRangeImag(2) > -1e-10
                                SweepRangeImag(2) = -1e-10;
                            end
                        end
                        if  any(MIN)
                            if  (N > 0 && (numel(find(X(:,1) ~= 0)) < 4 && i <= height(Scholte1{p+N-1}) && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < Scholte1{p+N-1}(i,4)*1e3)) ||...
                                (numel(find(X(:,1) ~= 0)) >= 1 && numel(find(X(:,1) ~= 0)) < 4 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) > X(end-1,1)) ||...
                                SweepRangeReal(MIN(1)) == 0 ||...
                                (Viscoelastic && SweepRangeImag(MIN(2)) == -1e-10)
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
                            if  numel(find(X(:,1) ~= 0)) < 4
                                SweepRangeReal = [1.001*SlowFluidVelocity .95*SlowFluidVelocity];
                            else
                                SweepRangeReal = [X(end-1,3)+.1*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                            end
                            if  SweepRangeReal(1) == SweepRangeReal(2)
                                if  numel(find(X(:,1) ~= 0)) < 4
                                    SweepRangeReal = [1.001*SlowFluidVelocity .95*SlowFluidVelocity];
                                else
                                    SweepRangeReal = [X(end-1,3)+.1*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                                end
                            end
                            if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                                SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                                SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                            end
                            if  SweepRangeReal(1) > 1.001*SlowFluidVelocity
                                SweepRangeReal(1) = 1.001*SlowFluidVelocity;
                            end
                            if  SweepRangeReal(2) < 0
                                SweepRangeReal(2) = 0;
                            end
                            if  numel(find(X(:,1) ~= 0)) < 4
                                SweepRangeImag = [-SweepRangeImagStart*SlowFluidVelocity -1e-10];
                            else
                                SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary velocity corresponds to the attenuation
                                if  -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                                    SweepRangeImag(2) = -1e-10;
                                end
                            end
                            if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                                SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                                SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                            end
                            if  SweepRangeImag(2) > -1e-10
                                SweepRangeImag(2) = -1e-10;
                            end
                        end
                    end
                    if  X(i,1) > 0 || (all(X(:,1) == 0) && q == SearchAreaExtensions-1) % stop q-loop if minimum has been found
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
                if  X(i,1) > 0
                    if  Multithreading
                        send(Q1,[FrequencyRange(i),X(i)/1e3,p+N])
                    else
                        addpoints(g(p+N),FrequencyRange(i),X(i)/1e3);
                        drawnow limitrate
                    end
                end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
            end
            if  all(X(:,1) == 0)
                MissingModes(p+N) = 1;
            end
            if  ~MissingModes(p+N)
                % [~,z] = max(X(:,1)); % find non-zero data below the cut-off frequency
                % X(1:z-1,:) = 0; % remove them
                % Misses(1:z-1) = 0; % remove also misses
                % if  X(z,2) == 0 % if frequency cut-off is at zero damping, remove it also
                %     X(z,:) = 0;
                % end
                X(Misses(1:height(X)) == 1,:) = NaN;
                Scholte1{p+N}(:,1) = FrequencyRange(1:height(X));
                Scholte1{p+N}(:,2) = FrequencyRange(1:height(X))/1e3;
                Scholte1{p+N}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
                Scholte1{p+N}(:,4) = fillmissing(X(:,1),'spline')/1e3;
                Scholte1{p+N}(:,7) = fillmissing(X(:,2),'spline');
                Scholte1{p+N}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
                Scholte1{p+N}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
                Scholte1{p+N}(Scholte1{p+N}(:,7) < 0,7) = 0; % negative attenuation is impossible
            else
                if  ModeFamily == 1 && p == 1
                    Scholte1{1}(1,1) = 0;
                else
                    Scholte1{p+N} = Scholte1{p+N-1};
                end
            end
        end
        Scholte1(MissingModes == 1) = [];
    end
end
if  TracingMode > 1
    if  ~isempty(H2) && ~isempty(Scholte1{1})
        zi = [];
        for p = 1:length(Scholte1)
            Scholte1{p} = flipud(Scholte1{p});
            z = find(isapprox(Scholte1{p}(1,8),H2(:,3),'veryloose'));
            if  ~isempty(z)
                zi(end+1) = z(1);
            end
        end
        H2(zi,:) = [];
    end
    if  ~isempty(H2)
        FrequencyRange_ = [FrequencyRange(1) 2*FrequencyRange(1):FrequencyResolution/5:x x+FrequencyResolution:FrequencyResolution:FrequencyRange(end)];
        FrequencyRange_(FrequencyRange_ > FrequencyRange(end)) = [];
        if  FrequencyRange_(end) ~= FrequencyRange(end)
            FrequencyRange_(end) = FrequencyRange(end);
        end
        FrequencyRange = fliplr(FrequencyRange_);
        if  ~Multithreading
            for p = 1:height(H2)
                if  ModeFamily == 1
                    g2(p) = animatedline(ax,'LineStyle','-.','color','r');
                elseif ModeFamily == 2
                    g2(p) = animatedline(ax,'LineStyle','-.','color','b');
                elseif ModeFamily == 0
                    g2(p) = animatedline(ax,'LineStyle','-.','color',[.5 0 1]);
                end
            end
        end
        if  ToggleUpperFluid && ToggleLowerFluid && ~strcmp(UpperFluid.Name,LowerFluid.Name)
            if  UpperFluid.Velocity > LowerFluid.Velocity
                FastFluidVelocity = UpperFluid.Velocity;
                SlowFluidVelocity = LowerFluid.Velocity;
            else
                FastFluidVelocity = LowerFluid.Velocity;
                SlowFluidVelocity = UpperFluid.Velocity;
            end
        elseif ToggleUpperFluid && (~ToggleLowerFluid || (ToggleLowerFluid && strcmp(UpperFluid.Name,LowerFluid.Name)))
            FastFluidVelocity = UpperFluid.Velocity;
            SlowFluidVelocity = UpperFluid.Velocity;
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            FastFluidVelocity = LowerFluid.Velocity;
            SlowFluidVelocity = LowerFluid.Velocity;
        end
        MissingModes = zeros(1,height(H2));
        Diff = diff(H2(:,3));
        Hdiff = min([[Diff;Inf] [Inf;Diff]],[],2);
        for p = 1:height(H2)
            X = H2(p,:);
            if  Multithreading
                send(Q2,[FrequencyRange(1),X(1)/1e3,p])
            else
                addpoints(g2(p),FrequencyRange(1),X(1)/1e3);
                drawnow limitrate
            end
            Misses = 0;
            for i = 2:length(FrequencyRange)
                if  Stop
                    return
                end
                AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
                kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
                gUpperFluid = 1i*UpperFluid.Density*AngularFrequency2;
                gLowerFluid = 1i*LowerFluid.Density*AngularFrequency2;
                for m = 1:SuperLayerSize
                    rw2(m) = Material{m}.Density*AngularFrequency2;
                    r2w4(m) = rw2(m)^2;
                    if  ~Decoupled
                        r3w6(m) = rw2(m)^3;
                        b12(m) = a12(m)*rw2(m);
                        b22(m) = a22(m)*rw2(m);
                        b23(m) = a23(m)*r2w4(m);
                        b32(m) = a32(m)*rw2(m);
                        b33(m) = a33(m)*r2w4(m);
                        b34(m) = a34(m)*r3w6(m);
                    else
                        b22(m) = a22(m)*rw2(m);
                        b32(m) = a32(m)*rw2(m);
                        b33(m) = r2w4(m);
                    end
                end
                X(i,1) = 0;
                Neighbors = [];
                for j = 1:length(Scholte1)
                    if  i <= height(Scholte1{j})
                        Neighbors(end+1,:) = [Scholte1{j}(i,8) Scholte1{j}(i,9)];
                    end
                end
                if  p > 1
                    for j = 1:length(Scholte2)
                        if  i <= height(Scholte2{j})
                            Neighbors(end+1,:) = [Scholte2{j}(i,8) Scholte2{j}(i,9)];
                        end
                    end
                end
                NeighborsNumber = height(Neighbors);
                for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal x SweepRangeImag)
                    if  numel(find(X(:,1) ~= 0)) < 4
                        if  Hdiff(p) > .001*X(end-1,3)
                            SweepRangeReal = [1.001*X(end-1,3) .999*X(end-1,3)];
                        else
                            SweepRangeReal = [X(end-1,3)+Hdiff(p) X(end-1,3)-Hdiff(p)];
                        end
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
                    if  numel(find(X(:,1) ~= 0)) < 4
                        SweepRangeImag = [1.1*X(end-1,4) .9*X(end-1,4)];
                    else
                        SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary velocity corresponds to the attenuation
                        if  -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                            SweepRangeImag(2) = -1e-10;
                        end
                    end
                    if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                        SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                        SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                    end
                    if  SweepRangeImag(2) > -1e-10
                        SweepRangeImag(2) = -1e-10;
                    end
                    for o = 1:SearchAreaSections % increase search resolution
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
                            if  ~Decoupled
                                Y = Computer_Coupled(2,ModeFamily,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,X,SlowFluidVelocity,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34);
                            else
                                Y = Computer_Decoupled(2,ModeFamily,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,X,SlowFluidVelocity,c,rw2,A1,a21,a31,b22,b32,b33);
                            end
                            if  ~isempty(Neighbors) % remove solutions of previously found lower modes
                                if  Viscoelastic
                                    for l = 2:width(Wavenumber)-1
                                        for j = 2:height(Wavenumber)-1
                                            for t = 1:NeighborsNumber
                                                if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1) && SweepRangeImag(l-1) < Neighbors(t,2) && SweepRangeImag(l+1) > Neighbors(t,2)
                                                    Y(j,l) = NaN;
                                                end
                                            end
                                        end
                                    end
                                else
                                    for j = 2:height(Wavenumber)-1
                                        for t = 1:NeighborsNumber
                                            if  SweepRangeReal(j-1) > Neighbors(t,1) && SweepRangeReal(j+1) < Neighbors(t,1)
                                                Y(j,:) = NaN;
                                            end
                                        end
                                    end
                                end
                            end
                            if  abs(SweepRangeImag(end)) < 1e-3
                                Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
                            end
% if  p==2%&&k==1
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
                                    for l = 1:length(b1) % calculate which minimum lies closest to previous solution
                                        cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                                        delta(l) = abs(cp(l)-X(end-1,1));
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
                                        if  SweepRangeImag(2) > -1e-10
                                            SweepRangeImag(2) = -1e-10;
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
                            if  SweepRangeImag(2) > -1e-10
                                SweepRangeImag(2) = -1e-10;
                            end
                        end
                        if  any(MIN)
                            if  numel(find(X(:,1) ~= 0)) < 4
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
                                X(i,3) = SweepRangeReal(MIN(1)); % real velocity (m/s)
                                X(i,4) = SweepRangeImag(MIN(2)); % imaginary velocity (m/s)
                                Misses(i) = 0;
                                break
                            end
                            if  numel(find(X(:,1) ~= 0)) < 4
                                if  Hdiff(p) > .001*X(end-1,3)
                                    SweepRangeReal = [1.001*X(end-1,3) .999*X(end-1,3)];
                                else
                                    SweepRangeReal = [X(end-1,3)+Hdiff(p) X(end-1,3)-Hdiff(p)];
                                end
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
                            if  numel(find(X(:,1) ~= 0)) < 4
                                SweepRangeImag = [1.1*X(end-1,4) .9*X(end-1,4)];
                            else
                                SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary velocity corresponds to the attenuation
                                if  -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                                    SweepRangeImag(2) = -1e-10;
                                end
                            end
                            if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                                SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                                SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                            end
                            if  SweepRangeImag(2) > -1e-10
                                SweepRangeImag(2) = -1e-10;
                            end
                        end
                    end
                    if  X(i,1) > 0 % stop q-loop if minimum has been found
                        break
                    end
                end
                if  X(i,1) == 0 % fit phase velocity, attenuation, and imaginary velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
                    if  isscalar(find(X(:,1) > 0))
                        MissingModes(p) = 1;
                        break
                    else
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
                        if  X(i,2) < 0 % exclude negative attenuation
                            X(i,[2 4]) = 0;
                        end
                        Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
                    end
                end
                if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % stop tracing the current dispersion curve if more than 'MissingSamples' have been missed
                    X(end-MissingSamples:end,:) = [];
                    Misses(end-MissingSamples:end) = 0;
                    break
                end
                if  X(end,1) > FastFluidVelocity-PhaseVelocityOffset || (X(1) < SlowFluidVelocity && X(end,1) > SlowFluidVelocity-PhaseVelocityOffset) || (X(1) > SlowFluidVelocity && X(end,1) < SlowFluidVelocity+PhaseVelocityOffset)
                    X(end,:) = [];
                    if  height(X) <= 10
                        MissingModes(p) = 1;
                    end
                    break
                end
                if  Multithreading
                    send(Q2,[FrequencyRange(i),X(i)/1e3,p])
                else
                    addpoints(g2(p),FrequencyRange(i),X(i)/1e3);
                    drawnow limitrate
                end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
            end
            if  ~MissingModes(p)
                X(Misses(1:height(X)) == 1,:) = NaN;
                Scholte2{p}(:,1) = FrequencyRange(1:height(X));
                Scholte2{p}(:,2) = FrequencyRange(1:height(X))/1e3;
                Scholte2{p}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
                Scholte2{p}(:,4) = fillmissing(X(:,1),'spline')/1e3;
                Scholte2{p}(:,7) = fillmissing(X(:,2),'spline');
                Scholte2{p}(:,8) = fillmissing(X(:,3),'spline'); % real velocity (m/s)
                Scholte2{p}(:,9) = fillmissing(X(:,4),'spline'); % imaginary velocity (m/s)
                Scholte2{p}(Scholte2{p}(:,7) < 0,7) = 0; % exclude negative attenuation
            else
                if  p == 1
                    Scholte2{1}(1,1) = 0;
                else
                    Scholte2{p} = Scholte2{p-1};
                end
            end
        end
        Scholte2(MissingModes == 1) = [];
    end
end
StartSorted = [];
if  ~isempty(Scholte1{1})
    for p = 1:length(Scholte1)
        Scholte1{p}(Scholte1{p}(:,4) == 0,:) = [];
        StartSorted(end+1,:) = [Scholte1{p}(end,1) 1 p];
    end
end
if  ~isempty(Scholte2) && ~isempty(Scholte2{1})
    for p = 1:length(Scholte2)
        StartSorted(end+1,:) = [Scholte2{p}(end,1) 2 p];
    end
end
if  ~isempty(StartSorted)
    StartSorted = sortrows(StartSorted,1);
    for p = 1:height(StartSorted)
        if  StartSorted(p,2) == 1
            Scholte{p} = flipud(Scholte1{StartSorted(p,3)}(:,1:7));
        elseif StartSorted(p,2) == 2
            Scholte{p} = flipud(Scholte2{StartSorted(p,3)}(:,1:7));
        end
    end
end
end
function Y = Computer_Coupled(Mode,ModeFamily,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,X,SlowFluidVelocity,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34)
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Wavenumber6 = Wavenumber2.^3;
    Length = length(Wavenumber);
    for m = 1:SuperLayerSize
        A1 = a11(m)*Wavenumber2+b12(m);
        A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
        A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
        d1 = A1/3;
        d2 = A2/3-d1.^2;
        d3 = d1.^3-d1.*A2/2+A3/2;
        d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
        d5 = d2./d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        k3(1,1,:) = sqrt(d6+d7);
        k3(1,2,:) = sqrt(d6-d7);
        k3(1,3,:) = sqrt(d4-d5-d1);
        k32 = k3.^2;
        k3k = k3.*Wavenumber;
        m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2(m);
        m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2(m);
        m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
        m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
        m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
        m1 = m13.*m22-m12.*m23;
        V = (m11.*m23-m13.*m12)./m1;
        W = (m11.*m22-m12.^2)./-m1;
        e1 = Wavenumber.*W+k3;
        e2 = k3.*V;
        D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V).*Wavenumber+c{m}(3,3)*k3.*W);
        D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
        D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
        if  Symmetric && SuperLayerSize == 1
            E = exp(.5i*k3*LayerThicknesses);
        else
            E = exp(1i*k3*LayerThicknesses(m));
        end
        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
        L2 = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W];
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M0 = L{m}(1:3,1:3,:)-M{1}(4:6,4:6,:);
        M1 = pagemrdivide(M{1}(1:3,4:6,:),M0);
        M2 = pagemrdivide(L{m}(4:6,1:3,:),M0);
        M{1} = [M{1}(1:3,1:3,:)+pagemtimes(M1,M{1}(4:6,1:3,:)) -pagemtimes(M1,L{m}(1:3,4:6,:));pagemtimes(M2,M{1}(4:6,1:3,:)) L{m}(4:6,4:6,:)-pagemtimes(M2,L{m}(1:3,4:6,:))];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1:3,1:3,:)-M{m}(4:6,4:6,:);
        M1 = pagemrdivide(M{m}(1:3,4:6,:),M0);
        M2 = pagemrdivide(M{Pattern(m)}(4:6,1:3,:),M0);
        M{m+1} = [M{m}(1:3,1:3,:)+pagemtimes(M1,M{m}(4:6,1:3,:)) -pagemtimes(M1,M{Pattern(m)}(1:3,4:6,:));pagemtimes(M2,M{m}(4:6,1:3,:)) M{Pattern(m)}(4:6,4:6,:)-pagemtimes(M2,M{Pattern(m)}(1:3,4:6,:))];
    end
    if  Symmetric
        M{end} = pageinv(M{end});
        if  ModeFamily == 1
            Y = reshape(abs(-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(3,1,:)-M{end}(3,4,:).*M{end}(6,1,:)./M{end}(6,4,:)),Size);
        elseif ModeFamily == 2
            Y = reshape(abs(-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(3,1,:)-M{end}(3,5,:).*M{end}(5,1,:)./M{end}(5,5,:)-(M{end}(3,5,:).*M{end}(5,6,:)./M{end}(5,5,:)-M{end}(3,6,:)).*(M{end}(4,1,:).*M{end}(5,5,:)-M{end}(4,5,:).*M{end}(5,1,:))./(M{end}(4,5,:).*M{end}(5,6,:)-M{end}(4,6,:).*M{end}(5,5,:))),Size);
        end
    else
        if  SymmetricSystem
            M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
            M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
            M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
            M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
        end
        M{end} = pageinv(M{end});
        if  ToggleUpperFluid && ToggleLowerFluid
            if  Mode == 1 || X(end-1,1) < SlowFluidVelocity
                k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
                k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
            else
                if  UpperFluid.Velocity > LowerFluid.Velocity
                    k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                else
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
                end
            end
            WUpperFluid = k3UpperFluid./Wavenumber;
            WLowerFluid = k3LowerFluid./Wavenumber;
            DUpperFluid = gUpperFluid./Wavenumber;
            DLowerFluid = gLowerFluid./Wavenumber;
            Y = reshape(abs((WLowerFluid-M{end}(6,4,:).*DLowerFluid).*(WUpperFluid+M{end}(3,1,:).*DUpperFluid)+M{end}(3,4,:).*M{end}(6,1,:).*DUpperFluid.*DLowerFluid),Size); 
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            Y = reshape(abs(-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(3,1,:)),Size);
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            Y = reshape(abs(-sqrt(kLowerFluid2-Wavenumber2)/gLowerFluid-M{end}(6,4,:)),Size);
        end
    end
end
function Y = Computer_Decoupled(Mode,ModeFamily,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,X,SlowFluidVelocity,c,rw2,A1,a21,a31,b22,b32,b33)
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Length = length(Wavenumber);
    for m = 1:SuperLayerSize
        A2 = a21(m)*Wavenumber2+b22(m);
        A3 = a31(m)*Wavenumber4+b32(m)*Wavenumber2+b33(m);
        d1 = sqrt(A2.^2-2*A1(m)*A3);
        k3(1,1,:) = sqrt((-A2+d1)/A1(m));
        k3(1,2,:) = sqrt((-A2-d1)/A1(m));
        W = (rw2(m)-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k3);
        D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W);
        D5 = 1i*(c{m}(5,5)*(k3+Wavenumber.*W));
        if  Symmetric && SuperLayerSize == 1
            E = exp(.5i*k3*LayerThicknesses);
        else
            E = exp(1i*k3*LayerThicknesses(m));
        end
        L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
        L2 = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W];
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M0 = L{m}(1:2,1:2,:)-M{1}(3:4,3:4,:);
        M1 = pagemrdivide(M{1}(1:2,3:4,:),M0);
        M2 = pagemrdivide(L{m}(3:4,1:2,:),M0);
        M{1} = [M{1}(1:2,1:2,:)+pagemtimes(M1,M{1}(3:4,1:2,:)) -pagemtimes(M1,L{m}(1:2,3:4,:));pagemtimes(M2,M{1}(3:4,1:2,:)) L{m}(3:4,3:4,:)-pagemtimes(M2,L{m}(1:2,3:4,:))];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1:2,1:2,:)-M{m}(3:4,3:4,:);
        M1 = pagemrdivide(M{m}(1:2,3:4,:),M0);
        M2 = pagemrdivide(M{Pattern(m)}(3:4,1:2,:),M0);
        M{m+1} = [M{m}(1:2,1:2,:)+pagemtimes(M1,M{m}(3:4,1:2,:)) -pagemtimes(M1,M{Pattern(m)}(1:2,3:4,:));pagemtimes(M2,M{m}(3:4,1:2,:)) M{Pattern(m)}(3:4,3:4,:)-pagemtimes(M2,M{Pattern(m)}(1:2,3:4,:))];
    end
    if  Symmetric
        M{end} = pageinv(M{end});
        if  ModeFamily == 1
            Y = reshape(abs(-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(2,1,:)-M{end}(2,3,:).*M{end}(4,1,:)./M{end}(4,3,:)),Size);
        elseif ModeFamily == 2
            Y = reshape(abs(-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(2,1,:)-M{end}(2,4,:).*M{end}(3,1,:)./M{end}(3,4,:)),Size);
        end
    else
        if  SymmetricSystem
            M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
            M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
            M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
            M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
        end
        M{end} = pageinv(M{end});
        if  ToggleUpperFluid && ToggleLowerFluid
            if  Mode == 1 || X(end-1,1) < SlowFluidVelocity
                k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
                k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
            else
                if  UpperFluid.Velocity > LowerFluid.Velocity
                    k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                else
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2);
                end
            end
            WUpperFluid = k3UpperFluid./Wavenumber;
            WLowerFluid = k3LowerFluid./Wavenumber;
            DUpperFluid = gUpperFluid./Wavenumber;
            DLowerFluid = gLowerFluid./Wavenumber;
            Y = reshape(abs((WLowerFluid-M{end}(4,3,:).*DLowerFluid).*(WUpperFluid+M{end}(2,1,:).*DUpperFluid)+M{end}(2,3,:).*M{end}(4,1,:).*DUpperFluid.*DLowerFluid),Size);
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            Y = reshape(abs(-sqrt(kUpperFluid2-Wavenumber2)/gUpperFluid+M{end}(2,1,:)),Size);
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            Y = reshape(abs(-sqrt(kLowerFluid2-Wavenumber2)/gLowerFluid-M{end}(4,3,:)),Size);
        end
    end
end