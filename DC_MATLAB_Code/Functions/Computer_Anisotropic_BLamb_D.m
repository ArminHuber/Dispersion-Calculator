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
function BLamb = Computer_Anisotropic_BLamb_D(Multithreading,Q1,Q2,ax,Decoupled,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,FrequencyRange,PlateThickness,HigherOrderModes,FSLamb,FALamb,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,MatrixMethods,MatrixMethodLimit,XS0)        
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
BLamb{1} = [];
if  ~Decoupled
    N = 3;
else
    N = 2;
end
if  ~Multithreading
    for p = 1:length(H)+N
        g(p) = animatedline(ax,'color',[.5 0 1]);
        g1(p) = animatedline(ax,'color',[.5 0 1]);
    end
end
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
if  strcmp(Material{1}.Class,'Isotropic')
    UpperMaterialVelocity = Material{1}.PlateVelocity;
else
    UpperMaterialVelocity = Material{1}.LongitudinalVelocity_1;
end
if  strcmp(Material{end}.Class,'Isotropic')
    LowerMaterialVelocity = Material{end}.PlateVelocity;
else
    LowerMaterialVelocity = Material{end}.LongitudinalVelocity_1;
end
MaterialVelocity = .5*(UpperMaterialVelocity+LowerMaterialVelocity);
MaterialDensity = .5*(Material{1}.Density+Material{end}.Density);
if  ToggleUpperFluid && ToggleLowerFluid
    FluidVelocity = .5*(UpperFluid.Velocity+LowerFluid.Velocity);
    FluidDensity = .5*(UpperFluid.Density+LowerFluid.Density);
    if  UpperFluid.Velocity > LowerFluid.Velocity
        FastFluidVelocity = UpperFluid.Velocity;
    else
        FastFluidVelocity = LowerFluid.Velocity;
    end
elseif ToggleUpperFluid && ~ToggleLowerFluid
    FluidVelocity = .5*UpperFluid.Velocity;
    FluidDensity = .5*UpperFluid.Density;
    FastFluidVelocity = UpperFluid.Velocity;
elseif ~ToggleUpperFluid && ToggleLowerFluid
    FluidVelocity = .5*LowerFluid.Velocity;
    FluidDensity = .5*LowerFluid.Density;
    FastFluidVelocity = LowerFluid.Velocity;
end
if  Viscoelastic
    for i = 1:length(Material)
        if  ~isreal(Material{i}.C)
            TV = 1e2*pi*(imag(Material{i}.C(1,1))/real(Material{i}.C(1,1))+imag(Material{i}.C(6,6))/real(Material{i}.C(6,6)));
            break
        end
    end
end
if  FluidLoading && Viscoelastic
    T = FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity)+TV;
elseif FluidLoading && ~Viscoelastic
    T = FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity);
elseif ~FluidLoading && Viscoelastic
    T = TV;
end
X = FALamb(1,:);
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
    if  MatrixMethods == 1
        if  X(i-1) > MatrixMethodLimit(i-1)
            MatrixMethod = 1; % TMM
        else
            MatrixMethod = 2; % SMM
        end
    else
        MatrixMethod = 2; % SMM
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
    for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
        if  i == 2
            SweepRangeReal = [5*X(end-1,3) X(end-1,3)];
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
        if  i == 2
            SweepRangeImag = [5*X(end-1,4) X(end-1,4)];
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
        if  FluidLoading && i > 2 && SweepRangeImag(1) > min(X(:,4)) && Scholte(i-1)
            SweepRangeImag(1) = min(X(:,4));
        end
        for o = 1:SearchAreaSections % increase search resolution
            if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                break
            end
            for k = 1:1e2
                if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                    if  i <= 3
                        SweepRangeReal = SweepRangeReal(1):.25^o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);                    
                    else
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
                    Y = Computer_Coupled(1,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,a12,a22,a23,a32,a33,a34,b12,b22,b23,b32,b33,b34);
                else
                    Y = Computer_Decoupled(1,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,a22,a32,b22,b32,b33);
                end
                if  FluidLoading
                    for j = 2:height(Wavenumber)-1 % avoid jumping to the Scholte mode
                        if  ToggleUpperFluid
                            if  SweepRangeReal(j-1) > UpperFluid.Velocity && SweepRangeReal(j+1) < UpperFluid.Velocity
                                if  UpperFluid.Density < FluidDensityThreshold
                                    Y(j,:) = NaN;
                                elseif UpperFluid.Density >= FluidDensityThreshold && SweepRangeImag(end) == 0
                                    Y(j,end) = NaN;
                                end
                            end
                        end
                        if  ToggleLowerFluid
                            if  SweepRangeReal(j-1) > LowerFluid.Velocity && SweepRangeReal(j+1) < LowerFluid.Velocity
                                if  LowerFluid.Density < FluidDensityThreshold
                                    Y(j,:) = NaN;
                                elseif LowerFluid.Density >= FluidDensityThreshold && SweepRangeImag(end) == 0
                                    Y(j,end) = NaN;
                                end
                            end
                        end
                    end
                end
                if  ((FluidLoading && (FluidDensity < FluidDensityThreshold || i > 100)) || ~FluidLoading) && abs(SweepRangeImag(end)) < 1e-3 % in case the fluid has low density, this allows finding solutions at zero attenuation for A0
                    Y(:,end+1) = max(Y,[],'all'); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  i>=5
% if  (FluidDensity < FluidDensityThreshold || i > 100) && abs(SweepRangeImag(end)) < 1e-3
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
                        [~,l] = min(delta);
                        MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                    end
                else
                    if  i <= 3 && k == 1
                        if  q < SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                            SweepRangeImag = [SweepRangeImag(1)-q*o*10*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                        else
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
                    else
                        Min = zeros(size(Y)); % find border minima
                        for j = 2:size(Y,1)-1
                            if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                Min(j,1) = 1;
                            end
                            if  SweepRangeImag(end) < 0 && Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
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
                    if  AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) > XS0(1)
                        Outlier = 1;
                    else
                        Outlier = 0;
                    end
                else
                    z1 = isoutlier(vertcat(X(1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
%                     z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                    if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1)% || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)
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
                    if  FluidLoading && SweepRangeImag(end) == 0
                        Scholte(i) = 1;
                    else
                        Scholte(i) = 0;
                    end
                    break
                end
                if  i == 2
                    SweepRangeReal = [5*X(end-1,3) X(end-1,3)];
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
                if  i == 2
                    SweepRangeImag = [5*X(end-1,4) X(end-1,4)];
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
                if  FluidLoading && i > 2 && SweepRangeImag(1) > min(X(:,4)) && Scholte(i-1)
                    SweepRangeImag(1) = min(X(:,4));
                end
            end
        end
        if  X(i,1) > 0 % stop q-loop if minimum has been found
            break
        end
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
        if  FluidLoading
            Scholte(i) = 0;
        end
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
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' M = ',num2str(MatrixMethod)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
end
X(Misses(1:height(X)) == 1,:) = NaN;
BLamb{1}(:,1) = FrequencyRange(1:height(X));
BLamb{1}(:,2) = FrequencyRange(1:height(X))/1e3;
BLamb{1}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
BLamb{1}(:,4) = fillmissing(X(:,1),'spline')/1e3;
BLamb{1}(:,7) = fillmissing(X(:,2),'spline');
BLamb{1}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
BLamb{1}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
BLamb{1}(BLamb{1}(:,7) < 0,7) = 0; % negative attenuation is impossible
for p = 2:N
    X = FSLamb(p-1,:);
    if  Multithreading
        send(Q1,[FrequencyRange(1),X(1)/1e3,p])
    else
        addpoints(g(p),FrequencyRange(1),X(1)/1e3);
        drawnow limitrate
    end
    Misses = 0;
    for i = 2:length(FrequencyRange)
        if  Stop
            return
        end
        if  MatrixMethods == 1
            if  X(i-1) > MatrixMethodLimit(i-1)
                MatrixMethod = 1; % TMM
            else
                MatrixMethod = 2; % SMM
            end
        else
            MatrixMethod = 2; % SMM
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
        for j = 1:length(BLamb)
            if  i <= height(BLamb{j})
                Neighbors(j,:) = [BLamb{j}(i,8) BLamb{j}(i,9)];
            end
        end
        NeighborsNumber = height(Neighbors);
        for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
            if  numel(find(X(:,1) ~= 0)) <= 3
                SweepRangeReal = [1.01*X(end-1,3) .99*X(end-1,3)];
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
            if  numel(find(X(:,1) ~= 0)) <= 3
                SweepRangeImag = [1.01*X(end-1,4) .99*X(end-1,4)];
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
                    if  ~Decoupled
                        Y = Computer_Coupled(1,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,a12,a22,a23,a32,a33,a34,b12,b22,b23,b32,b33,b34);
                    else
                        Y = Computer_Decoupled(1,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,a22,a32,b22,b32,b33);
                    end
                    if  ~isempty(Neighbors)
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
                            if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
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
                    if  numel(find(X(:,1) ~= 0)) <= 3
                        if  SweepRangeReal(1) > 1.01*XS0(p-1)
                            SweepRangeReal(1) = 1.01*XS0(p-1);
                        end
                        if  SweepRangeReal(2) < .99*XS0(p-1)
                            SweepRangeReal(2) = .99*XS0(p-1);
                        end
                    else
                        if  SweepRangeReal(2) < 0
                            SweepRangeReal(2) = 0;
                        end
                    end
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                end
                if  any(MIN)
                    if  (numel(find(X(:,1) ~= 0)) <= 3 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < .9*XS0(p-1)) ||...
                        SweepRangeReal(MIN(1)) == 0
                        Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                        NeighborsNumber = height(Neighbors);
                    else
                        if  numel(find(X(:,1) ~= 0)) <= 3
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
                    if  numel(find(X(:,1) ~= 0)) <= 3
                        SweepRangeReal = [1.01*X(end-1,3) .99*X(end-1,3)];
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
                    if  numel(find(X(:,1) ~= 0)) <= 3
                        SweepRangeImag = [1.01*X(end-1,4) .99*X(end-1,4)];
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
        if  FluidLoading && any(X(:,1)) && X(end,1) < FastFluidVelocity % avoid jumping to Scholte modes
            X(end,:) = [];
            break
        end
        if  Multithreading
            send(Q1,[FrequencyRange(i),X(i)/1e3,p])
        else
            addpoints(g(p),FrequencyRange(i),X(i)/1e3);
            drawnow limitrate
        end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' M = ',num2str(MatrixMethod)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
    end
    X(Misses(1:height(X)) == 1,:) = NaN;
    BLamb{p}(:,1) = FrequencyRange(1:height(X));
    BLamb{p}(:,2) = FrequencyRange(1:height(X))/1e3;
    BLamb{p}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
    BLamb{p}(:,4) = fillmissing(X(:,1),'spline')/1e3;
    BLamb{p}(:,7) = fillmissing(X(:,2),'spline');
    BLamb{p}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
    BLamb{p}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
    BLamb{p}(BLamb{p}(:,7) < 0,7) = 0; % negative attenuation is impossible
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
            X(i,1) = 0;
            if  MatrixMethods == 1
                if  all(X(:,1) == 0) || X(i-1) > MatrixMethodLimit(i-1)
                    MatrixMethod = 1; % TMM
                else
                    MatrixMethod = 2; % SMM
                end
            else
                MatrixMethod = 2; % SMM
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
            Neighbors = [];
            for j = 1:length(BLamb)
                if  i <= height(BLamb{j})
                    Neighbors(j,:) = [BLamb{j}(i,8) BLamb{j}(i,9)];
                end
            end
            NeighborsNumber = height(Neighbors);
            for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                if  all(X(:,1) == 0) 
                    SweepRangeReal = [1.1*PhaseVelocityLimit XS0(end)];
                elseif isscalar(find(X(:,1) ~= 0))
                    SweepRangeReal = [PhaseVelocityLimit XS0(end)];
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
                if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < XS0(end)
                    SweepRangeReal(2) = XS0(end);
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
                        if  ~Decoupled
                            Y = Computer_Coupled(1,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,a12,a22,a23,a32,a33,a34,b12,b22,b23,b32,b33,b34);
                        else
                            Y = Computer_Decoupled(1,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,a22,a32,b22,b32,b33);
                        end
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
% if  p==3&&i>=149
% if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
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
                                    if  SweepRangeReal(2) < XS0(end)
                                        SweepRangeReal(2) = XS0(end);
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
                        if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < XS0(end)
                            SweepRangeReal(2) = XS0(end);
                        elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                            SweepRangeReal(2) = 0;
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                    end
                    if  any(MIN)
                        if  (numel(find(X(:,1) ~= 0)) <= 3 && i <= height(BLamb{p+N-1}) && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < BLamb{p+N-1}(i,4)*1e3) ||...
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
                            SweepRangeReal = [1.1*PhaseVelocityLimit XS0(end)];
                        elseif isscalar(find(X(:,1) ~= 0))
                            SweepRangeReal = [PhaseVelocityLimit XS0(end)];
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
                        if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < XS0(end)
                            SweepRangeReal(2) = XS0(end);
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
            if  FluidLoading && any(X(:,1)) && X(end,1) < FastFluidVelocity % avoid jumping to Scholte modes
                X(end,:) = [];
                break
            end
            if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                send(Q1,[FrequencyRange(i),X(i)/1e3,p+N])
            elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                addpoints(g(p+N),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' M = ',num2str(MatrixMethod)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
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
            BLamb{p+N}(:,1) = FrequencyRange(1:height(X));
            BLamb{p+N}(:,2) = FrequencyRange(1:height(X))/1e3;
            BLamb{p+N}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
            BLamb{p+N}(:,4) = fillmissing(X(:,1),'spline')/1e3;
            BLamb{p+N}(:,7) = fillmissing(X(:,2),'spline');
            BLamb{p+N}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
            BLamb{p+N}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
            BLamb{p+N}(BLamb{p+N}(:,7) < 0,7) = 0; % negative attenuation is impossible
        else
            BLamb{p+N} = BLamb{p+N-1};
            X1{p+N}(1,1) = 0;
            continue
        end
        if  max(BLamb{p+N}(:,4))*1e3 > PhaseVelocityLimit
            X1{p+N}(1,1) = 0;
        else
            [Max,MaxInd] = max(BLamb{p+N}(:,8));
            PhaseVelocityRange = Max+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
            for i = 1:length(PhaseVelocityRange)
                if  Stop
                    return
                end
                if  MatrixMethods == 1
                    if  PhaseVelocityRange(i) > MatrixMethodLimit(MaxInd)
                        MatrixMethod = 1; % TMM
                    else
                        MatrixMethod = 2; % SMM
                    end
                else
                    MatrixMethod = 2; % SMM
                end
                X1{p+N}(i,1) = 0;
                for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                    if  i == 1
                        SweepRangeFrq = [BLamb{p+N}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 BLamb{p+N}(MaxInd,1)+4/PlateThickness/1e3];
                        SweepRangeImag = [SearchWidthImag(1)*BLamb{p+N}(MaxInd,9) SearchWidthImag(2)*BLamb{p+N}(MaxInd,9)];
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
                            kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
                            kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
                            gUpperFluid = 1i*UpperFluid.Density*AngularFrequency2;
                            gLowerFluid = 1i*LowerFluid.Density*AngularFrequency2;
                            Wavenumber = AngularFrequency'./(PhaseVelocityRange(i)+SweepRangeImag*1i);
                            if  ~Decoupled
                                Y = Computer_Coupled(2,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,a12,a22,a23,a32,a33,a34,b12,b22,b23,b32,b33,b34);
                            else
                                Y = Computer_Decoupled(2,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,a22,a32,b22,b32,b33);
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
                                            delta(l) = abs(frq(l)-BLamb{p+N}(MaxInd,1));
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
                                SweepRangeFrq = [BLamb{p+N}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 BLamb{p+N}(MaxInd,1)+4/PlateThickness/1e3];
                                SweepRangeImag = [SearchWidthImag(1)*BLamb{p+N}(MaxInd,9) SearchWidthImag(2)*BLamb{p+N}(MaxInd,9)];
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
                if  Multithreading
                    send(Q2,[X1{p+N}(i,1),X1{p+N}(i,4)/1e3,p+N])
                else
                    addpoints(g1(p+N),X1{p+N}(i,1),X1{p+N}(i,4)/1e3);
                    drawnow limitrate
                end
% disp(['c = ',num2str(PhaseVelocityRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' M = ',num2str(MatrixMethod)]);
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
    BLamb(MissingModes == 1) = [];
    for p = 1:N
        BLamb{p}(:,8:9) = [];
    end
    for p = N+1:length(BLamb)
        BLamb{p}(BLamb{p}(:,4) == 0,:) = [];
        BLamb{p}(:,8:9) = [];
        if  X1{p}(1,1) > 0
            BLamb{p} = vertcat(X1{p},BLamb{p});
        end
    end
end
% figure,hold on
% for p = 1:length(BLamb)
%     plot(BLamb{p}(:,1),BLamb{p}(:,4),'k');
% end
% line(FrequencyRange,MatrixMethodLimit/1e3,'color','r')
end
function Y = Computer_Coupled(Mode,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,a11,a21,a31,a12,a22,a23,a32,a33,a34,b12,b22,b23,b32,b33,b34)
    % MatrixMethod = 2;

    if  Mode == 2
        AngularFrequency2 = reshape(repmat(AngularFrequency2,width(Wavenumber),1),1,1,[]);
        kUpperFluid2 = reshape(repmat(kUpperFluid2,width(Wavenumber),1),1,1,[]);
        kLowerFluid2 = reshape(repmat(kLowerFluid2,width(Wavenumber),1),1,1,[]);
        gUpperFluid = reshape(repmat(gUpperFluid,width(Wavenumber),1),1,1,[]);
        gLowerFluid = reshape(repmat(gLowerFluid,width(Wavenumber),1),1,1,[]);
    end
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Wavenumber6 = Wavenumber2.^3;
    Length = length(Wavenumber);
    Y = NaN(size(Wavenumber)); % to be discarded once pagedet exists
    for m = 1:SuperLayerSize
        if  Mode == 1
            A1 = a11(m)*Wavenumber2+b12(m);
            A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
            A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
        elseif Mode == 2
            rw2 = Material{m}.Density*AngularFrequency2;
            r2w4 = rw2.^2;
            A1 = a11(m)*Wavenumber2+a12(m)*rw2;
            A2 = a21(m)*Wavenumber4+a22(m)*rw2.*Wavenumber2+a23(m)*r2w4;
            A3 = a31(m)*Wavenumber6+a32(m)*rw2.*Wavenumber4+a33(m)*r2w4.*Wavenumber2+a34(m)*rw2.^3;
        end
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
        if  Mode == 1
            m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2(m);
            m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2(m);
        elseif Mode == 2
            m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2;
            m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2;
        end
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
        Phi = 1i*k3*LayerThicknesses(m);
        E = exp(Phi);
        if  MatrixMethod == 1
            E_ = exp(-Phi);
            L1 = [E E_;V.*E V.*E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_;D4.*E -D4.*E_];
            L2 = [ones(1,6,Length);V V;W -W;D3 D3;D5 -D5;D4 -D4];
        elseif MatrixMethod == 2
            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
            L2 = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W];
        end
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    if  MatrixMethod == 1
        for m = 2:SuperLayerSize
            M{1} = pagemtimes(M{1},L{m});
        end
        for m = 1:length(Pattern)
            M{m+1} = pagemtimes(M{m},M{Pattern(m)});
        end
        if  SymmetricSystem
            M2{1} = L{end};
            for m = SuperLayerSize-1:-1:1
                M2{1} = pagemtimes(M2{1},L{m});
            end
            for m = 1:length(Pattern)
                M2{m+1} = pagemtimes(M2{m},M2{Pattern(m)});
            end
            M{end} = pagemtimes(M{end},M2{end});
        end
        if  FluidLoading
            if  ToggleUpperFluid && ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                QUpperFluid = gUpperFluid./k3UpperFluid;
                QLowerFluid = gLowerFluid./k3LowerFluid;
                for j = 1:length(Wavenumber)
                    M11(1,1,j) = det(M{end}([3 5:6],1:3,j));
                    M12(1,1,j) = det(M{end}([3 5:6],[1:2 4],j));
                    M21(1,1,j) = det(M{end}(4:6,1:3,j));
                    M22(1,1,j) = det(M{end}(4:6,[1:2 4],j));
                end
                Y = abs(M21-QLowerFluid.*M22-QUpperFluid.*(M11-QLowerFluid.*M12));
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                for j = 1:length(Wavenumber)
                    M11(1,1,j) = det(M{end}([3 5:6],1:3,j));
                    M21(1,1,j) = det(M{end}(4:6,1:3,j));
                end
                Y = abs(M21-gUpperFluid./k3UpperFluid.*M11);
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                for j = 1:length(Wavenumber)
                    M21(1,1,j) = det(M{end}(4:6,1:3,j));
                    M22(1,1,j) = det(M{end}(4:6,[1:2 4],j));
                end
                Y = abs(M21-gLowerFluid./k3LowerFluid.*M22);
            end
        else
            for j = 1:length(Wavenumber)
                Y(j) = abs(det(M{end}(4:6,1:3,j)));
            end
        end
    elseif MatrixMethod == 2
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
        if  SymmetricSystem
            M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
            M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
            M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
            M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
        end
        if  FluidLoading
            M{end} = pageinv(M{end});
            if  ToggleUpperFluid && ToggleLowerFluid
                WUpperFluid = sqrt(kUpperFluid2-Wavenumber2)./Wavenumber;
                WLowerFluid = sqrt(kLowerFluid2-Wavenumber2)./Wavenumber;
                DUpperFluid = gUpperFluid./Wavenumber;
                DLowerFluid = gLowerFluid./Wavenumber;
                Y = abs((WLowerFluid-M{end}(6,4,:).*DLowerFluid).*(WUpperFluid+M{end}(3,1,:).*DUpperFluid)+M{end}(3,4,:).*M{end}(6,1,:).*DUpperFluid.*DLowerFluid); 
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                Y = abs(sqrt(kUpperFluid2-Wavenumber2)./gUpperFluid+M{end}(3,1,:));
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                Y = abs(sqrt(kLowerFluid2-Wavenumber2)./gLowerFluid-M{end}(6,4,:));
            end
        else
            for j = 1:length(Wavenumber)
                Y(j) = abs(det(M{end}(:,:,j)));
            end
        end
    end
    Y = reshape(Y,Size);
end
function Y = Computer_Decoupled(Mode,MatrixMethod,AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,rw2,A1,a21,a31,a22,a32,b22,b32,b33)
    % MatrixMethod = 2;

    if  Mode == 2
        AngularFrequency2 = reshape(repmat(AngularFrequency2,width(Wavenumber),1),1,1,[]);
        kUpperFluid2 = reshape(repmat(kUpperFluid2,width(Wavenumber),1),1,1,[]);
        kLowerFluid2 = reshape(repmat(kLowerFluid2,width(Wavenumber),1),1,1,[]);
        gUpperFluid = reshape(repmat(gUpperFluid,width(Wavenumber),1),1,1,[]);
        gLowerFluid = reshape(repmat(gLowerFluid,width(Wavenumber),1),1,1,[]);
    end
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Length = length(Wavenumber);
    Y = NaN(size(Wavenumber)); % to be discarded once pagedet exists
    for m = 1:SuperLayerSize
        if  Mode == 1
            A2 = a21(m)*Wavenumber2+b22(m);
            A3 = a31(m)*Wavenumber4+b32(m)*Wavenumber2+b33(m);
        elseif Mode == 2
            rw2 = Material{m}.Density*AngularFrequency2;
            A2 = a21(m)*Wavenumber2+a22(m)*rw2;
            A3 = a31(m)*Wavenumber4+a32(m)*rw2.*Wavenumber2+rw2.^2;
        end
        d1 = sqrt(A2.^2-2*A1(m)*A3);
        k3(1,1,:) = sqrt((-A2+d1)/A1(m));
        k3(1,2,:) = sqrt((-A2-d1)/A1(m));
        if  Mode == 1
            W = (rw2(m)-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k3);
        elseif Mode == 2
            W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k3);
        end
        D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W);
        D5 = 1i*(c{m}(5,5)*(k3+Wavenumber.*W));
        Phi = 1i*k3*LayerThicknesses(m);
        E = exp(Phi);
        if  MatrixMethod == 1
            E_ = exp(-Phi);
            L1 = [E E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_];
            L2 = [ones(1,4,Length);W -W;D3 D3;D5 -D5];
        elseif MatrixMethod == 2
            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
            L2 = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W];
        end
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    if  MatrixMethod == 1
        for m = 2:SuperLayerSize
            M{1} = pagemtimes(M{1},L{m});
        end
        for m = 1:length(Pattern)
            M{m+1} = pagemtimes(M{m},M{Pattern(m)});
        end
        if  SymmetricSystem
            M2{1} = L{end};
            for m = SuperLayerSize-1:-1:1
                M2{1} = pagemtimes(M2{1},L{m});
            end
            for m = 1:length(Pattern)
                M2{m+1} = pagemtimes(M2{m},M2{Pattern(m)});
            end
            M{end} = pagemtimes(M{end},M2{end});
        end
        if  FluidLoading
            if  ToggleUpperFluid && ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                QUpperFluid = gUpperFluid./k3UpperFluid;
                QLowerFluid = gLowerFluid./k3LowerFluid;
                for j = 1:length(Wavenumber)
                    M11(1,1,j) = det(M{end}([2 4],1:2,j));
                    M12(1,1,j) = det(M{end}([2 4],[1 3],j));
                    M21(1,1,j) = det(M{end}(3:4,1:2,j));
                    M22(1,1,j) = det(M{end}(3:4,[1 3],j));
                end
                Y = abs(M21-QLowerFluid.*M22-QUpperFluid.*(M11-QLowerFluid.*M12));
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                for j = 1:length(Wavenumber)
                    M11(1,1,j) = det(M{end}([2 4],1:2,j));
                    M21(1,1,j) = det(M{end}(3:4,1:2,j));
                end
                Y = abs(M21-gUpperFluid./k3UpperFluid.*M11);
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                for j = 1:length(Wavenumber)
                    M21(1,1,j) = det(M{end}(3:4,1:2,j));
                    M22(1,1,j) = det(M{end}(3:4,[1 3],j));
                end
                Y = abs(M21-gLowerFluid./k3LowerFluid.*M22);
            end
        else
            for j = 1:length(Wavenumber)
                Y(j) = abs(det(M{end}(3:4,1:2,j)));
            end
        end
    elseif MatrixMethod == 2
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
        if  SymmetricSystem
            M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
            M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
            M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
            M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
        end
        if  FluidLoading
            M{end} = pageinv(M{end});
            if  ToggleUpperFluid && ToggleLowerFluid
                WUpperFluid = sqrt(kUpperFluid2-Wavenumber2)./Wavenumber;
                WLowerFluid = sqrt(kLowerFluid2-Wavenumber2)./Wavenumber;
                DUpperFluid = gUpperFluid./Wavenumber;
                DLowerFluid = gLowerFluid./Wavenumber;
                Y = abs((WLowerFluid-M{end}(4,3,:).*DLowerFluid).*(WUpperFluid+M{end}(2,1,:).*DUpperFluid)+M{end}(2,3,:).*M{end}(4,1,:).*DUpperFluid.*DLowerFluid); 
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                Y = abs(sqrt(kUpperFluid2-Wavenumber2)./gUpperFluid+M{end}(2,1,:));
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                Y = abs(sqrt(kLowerFluid2-Wavenumber2)./gLowerFluid-M{end}(4,3,:));
            end
        else
            for j = 1:length(Wavenumber)
                Y(j) = abs(det(M{end}(:,:,j)));
            end
        end
    end
    Y = reshape(Y,Size);
end