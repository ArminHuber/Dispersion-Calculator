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
function A = Computer_2DTracing(Multithreading,Q,ax,Geometry,Material,Ro,Ri,FluidLoading,ToggleFluid1,ToggleFluid2,Fluid1,Fluid2,Sink,Symmetric,FrequencyRange,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,h,LineColor,n,Decoupled,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};Fit={};
if  isstruct(Material)
    MaterialType = 1;
    Velocity1 = Material.LongitudinalVelocity;
    Density1 = Fluid1.Density/Material.Density;
    Density2 = Fluid2.Density/Material.Density;
    cL2 = Material.LongitudinalVelocity_complex^2;
    cT2 = Material.TransverseVelocity_complex^2;
else
    MaterialType = 2;
    Velocity1 = 1.1*XS0(end);
end
if  strcmp(Geometry,'Plate')
    Geometry = 1;
    ModeFamily = Ri;
    if  MaterialType == 1
        Lambda = conj(Material.Lambda_complex);
        Mu = conj(Material.Mu_complex);
        Half = Ro;
        ScaleFactor = 2*Half;
    else
        ScaleFactor = Ro;
    end
elseif strcmp(Geometry,'Rod')
    Geometry = 2;
    ScaleFactor = 2*Ro;
elseif strcmp(Geometry,'Pipe')
    Geometry = 3;
    ScaleFactor = Ro-Ri;
end
n2 = n^2;
Ro2 = Ro^2;
Ri2 = Ri^2;
cF12 = Fluid1.Velocity^2;
cF22 = Fluid2.Velocity^2;
H0 = [FrequencyRange(2);FrequencyRange(10)];
for i = 1:2*Sweeps
    H0(end+1,1) = FrequencyRange(find(FrequencyRange >= i/(2*Sweeps+1)*FrequencyRange(end),1));
end
H0(end+1,1) = FrequencyRange(end-10);
H0(:,8) = 1;
if  FluidLoading
    if  Geometry < 3 || (Geometry == 3 && ~(~ToggleFluid1 && ToggleFluid2 && ~Sink))
        H0 = [H0;H0(2:end,:)];
        H0(end-floor(size(H0,1)/2)+1:end,8) = 2;
    end
    SolutionTypeRange = 1:2;
    if  Geometry == 2
        ContinuousAcrossFastFluidVelocity = false;
        ContinuousAcrossSlowFluidVelocity = true;
        FastFluidVelocity = Fluid1.Velocity;
    else
        if  Geometry == 1
            ContinuousAcrossFastFluidVelocity = false;
            if  ToggleFluid1 && ToggleFluid2 && ~Symmetric
                ContinuousAcrossSlowFluidVelocity = false;
            else
                ContinuousAcrossSlowFluidVelocity = true;
            end
        elseif Geometry == 3
            if  ToggleFluid2 && ~Sink && (~ToggleFluid1 || (ToggleFluid1 && Fluid1.Velocity < Fluid2.Velocity))
                ContinuousAcrossFastFluidVelocity = true;
            else
                ContinuousAcrossFastFluidVelocity = false;
            end
            if  ToggleFluid1 && ToggleFluid2 && ~Symmetric && (Sink || Fluid1.Velocity < Fluid2.Velocity)
                ContinuousAcrossSlowFluidVelocity = false;
            else
                ContinuousAcrossSlowFluidVelocity = true;
            end
        end
        if  ~Symmetric && ToggleFluid1 && ToggleFluid2
            if  Fluid1.Velocity > Fluid2.Velocity
                FastFluidVelocity = Fluid1.Velocity;
                SlowFluidVelocity = Fluid2.Velocity;
            else
                FastFluidVelocity = Fluid2.Velocity;
                SlowFluidVelocity = Fluid1.Velocity;
            end
        elseif Symmetric || (ToggleFluid1 && ~ToggleFluid2)
            FastFluidVelocity = Fluid1.Velocity;
            SlowFluidVelocity = Fluid1.Velocity;
        elseif ~ToggleFluid1 && ToggleFluid2
            FastFluidVelocity = Fluid2.Velocity;
            SlowFluidVelocity = Fluid2.Velocity;
        end
    end
else
    SolutionTypeRange = 1;
    SignFluid1 = 1;
    SignFluid2 = 1;
    if  MaterialType == 2
        SignFluid11 = 1;
        SignFluid21 = 1;
        SignFluid12 = 1;
        SignFluid22 = 1;
    end
end
if  Geometry == 3
    if  any(h)
        FineFrequencyLimit = 1.5*h(1);
    else
        FineFrequencyLimit = FrequencyRange(end);
    end
    if  ToggleFluid2 && ~Sink
        FineFrequencyStep = 2;
    else
        FineFrequencyStep = 10;
    end
else
    FineFrequencyLimit = FrequencyRange(10);
    FineFrequencyStep = 5;
end
if  MaterialType == 1
    if  Geometry == 1
        if  ModeFamily == 1
            Neighbors0 = zeros(0,2); % initiate solutions to be ignored during dispersion curve tracing
        elseif ModeFamily == 2
            Neighbors0 = [real(Material.LongitudinalVelocity_complex) imag(Material.LongitudinalVelocity_complex)];
        else
            Neighbors0 = [real(Material.TransverseVelocity_complex) imag(Material.TransverseVelocity_complex);real(Material.LongitudinalVelocity_complex) imag(Material.LongitudinalVelocity_complex)];
        end
    elseif Geometry == 2
        Neighbors0 = zeros(0,2);
        if  ~(n == 0 && ~ToggleFluid1)
            Neighbors0(end+1,:) = [real(Material.TransverseVelocity_complex) imag(Material.TransverseVelocity_complex)];
        end
        if  (n == 1 && ToggleFluid1) || n > 1
            Neighbors0(end+1,:) = [real(Material.LongitudinalVelocity_complex) imag(Material.LongitudinalVelocity_complex)];
        end
    elseif Geometry == 3
        Neighbors0 = [real(Material.TransverseVelocity_complex) imag(Material.TransverseVelocity_complex)];
        if  n > 0 && ToggleFluid2 && ~Sink && (~ToggleFluid1 || (ToggleFluid1 && ~Symmetric))
            Neighbors0(end+1,:) = [Fluid2.Velocity 0];
        end
    end
else
    Neighbors0 = zeros(0,2);
    if  ToggleFluid1
        Neighbors0(end+1,:) = [Fluid1.Velocity 0];
    end
    if  ToggleFluid2
        Neighbors0(end+1,:) = [Fluid2.Velocity 0];
    end
    Neighbors0 = unique(Neighbors0,'rows');
end
if  ~Multithreading
    g = [animatedline(ax,'color',LineColor) animatedline(ax,'color',LineColor)];
end
for SolutionType = SolutionTypeRange % 1: trace forward leaky and backward Scholte modes 2: trace backward leaky and forward Scholte modes 
    H = H0;
    if  SolutionType == 1
        H = [h' zeros(length(h),7);H];
    end
    p = 0;
    while p < size(H,1) % step through the modes
        p = p+1;
        X = 0;
        Misses = false;
        BelowCutoff = false;
        if  H(p,2) || H(p,8)
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
                X(i,6) = 0;
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
                        if  Frequency(i-1) < FineFrequencyLimit
                            Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution/GearFrequencyResolution(Gear(i))/FineFrequencyStep;
                            if  Gear(i) == 1
                                if  r == 1 && FrequencyRange(end) > Frequency(i)
                                    z = find(FrequencyRange > Frequency(i-1),1);
                                    if  z == length(FrequencyRange)
                                        f = FrequencyRange(z-1):FrequencyResolution/FineFrequencyStep:FrequencyRange(z);
                                    else
                                        f = FrequencyRange(z-1):FrequencyResolution/FineFrequencyStep:FrequencyRange(z+1);
                                    end
                                    if  ~any(isapprox(f,Frequency(i),'tight'))
                                        Frequency(i) = f(find(f > Frequency(i-1),1));
                                    end
                                elseif r == 2 && FrequencyRange(1) < Frequency(i)
                                    z = find(FrequencyRange < Frequency(i-1),1,'last');
                                    if  z == 1
                                        f = FrequencyRange(z):FrequencyResolution/FineFrequencyStep:FrequencyRange(z+1);
                                    else
                                        f = FrequencyRange(z-1):FrequencyResolution/FineFrequencyStep:FrequencyRange(z+1);
                                    end
                                    if  ~any(isapprox(f,Frequency(i),'tight'))
                                        Frequency(i) = f(find(f < Frequency(i-1),1,'last'));
                                    end
                                end
                            end
                        else
                            Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution/GearFrequencyResolution(Gear(i));
                            if  Gear(i) == 1 && ~any(isapprox(FrequencyRange,Frequency(i),'tight'))
                                if  r == 1 && FrequencyRange(end) > Frequency(i)
                                    Frequency(i) = FrequencyRange(find(FrequencyRange > Frequency(i-1),1));
                                elseif r == 2 && FrequencyRange(1) < Frequency(i)
                                    Frequency(i) = FrequencyRange(find(FrequencyRange < Frequency(i-1),1,'last'));
                                end
                            end
                        end
                    else
                        if  Frequency(i-1) < FineFrequencyLimit
                            if  sqrt(CriticalSlope/Slope) < 1/GearFrequencyResolution(end)
                                Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution*sqrt(CriticalSlope/Slope)/FineFrequencyStep;
                            else
                                Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution/GearFrequencyResolution(end)/FineFrequencyStep;
                            end
                        else
                            if  sqrt(CriticalSlope/Slope) < 1/GearFrequencyResolution(end)
                                Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution*sqrt(CriticalSlope/Slope);
                            else
                                Frequency(i) = Frequency(i-1)+Direction*FrequencyResolution/GearFrequencyResolution(end);
                            end
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
                kF12 = AngularFrequency2/cF12;
                kF22 = AngularFrequency2/cF22;
                if  MaterialType == 1
                    kL2 = AngularFrequency2/cL2;
                    kT2 = AngularFrequency2/cT2;
                    gF1 = Fluid1.Density*AngularFrequency2;
                    gF2 = Fluid2.Density*AngularFrequency2;
                else
                    if  MatrixMethods == 1
                        [~,z] = min(abs(Frequency(i)-FrequencyRange));
                        Limit = MatrixMethodLimit(z);
                        if  ~any(X(:,1)) || X(i-1) > Limit
                            MatrixMethod = 1; % TMM
                        else
                            MatrixMethod = 2; % SMM
                        end
                    else
                        MatrixMethod = 2; % SMM
                    end
                    if  Frequency(i) < FrequencyRange(3) % SMM becomes unstable at very low frequency-thicknesses while TMM becomes unstable at high f*d
                        MatrixMethod = 1; % TMM
                    end
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
                    gF1 = 1i*Fluid1.Density*AngularFrequency2;
                    gF2 = 1i*Fluid2.Density*AngularFrequency2; 
                end
                for q = 1:SearchAreaExtensions % q > 1 iterations are with extended search area (SweepRangeReal x SweepRangeImag)
                    if  numel(find(X(:,1))) < 2
                        if  ~H(p,2) || isscalar(find(X(:,1)))
                            if  ~H(p,8)
                                if  any(X(:,1))
                                    SweepRangeReal = X(i-1,3)+PhaseVelocityLimit*[.25 -.25];
                                    SweepRangeImag = X(i-1,4)+PhaseVelocityLimit*[-.25 .25];
                                else
                                    SweepRangeReal = [PhaseVelocityLimit 1];
                                    SweepRangeImag = PhaseVelocityLimit*[-1 1];
                                end
                            elseif H(p,8) == 1
                                if  any(X(:,1))
                                    SweepRangeReal = X(i-1,3)+Velocity1*[.25 -.25];
                                    SweepRangeImag = X(i-1,4)+Velocity1*[-.25 .25];
                                else
                                    SweepRangeReal = [Velocity1 1];
                                    SweepRangeImag = Velocity1*[-1 1];
                                end
                            elseif H(p,8) == 2
                                if  any(X(:,1))
                                    SweepRangeReal = X(i-1,3)+FastFluidVelocity*[.25 -.25];
                                    SweepRangeImag = X(i-1,4)+FastFluidVelocity*[-.25 .25];
                                else
                                    SweepRangeReal = [FastFluidVelocity 1];
                                    SweepRangeImag = FastFluidVelocity*[-1 1];
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
                                        if  ~H(p,8)
                                            SweepRangeReal = SweepRangeReal(1):-PhaseVelocityLimit/SweepSections:SweepRangeReal(end);
                                            SweepRangeImag = SweepRangeImag(1):PhaseVelocityLimit/SweepSections:SweepRangeImag(end);
                                        elseif H(p,8) == 1
                                            SweepRangeReal = SweepRangeReal(1):-Velocity1/SweepSections:SweepRangeReal(end);
                                            SweepRangeImag = SweepRangeImag(1):Velocity1/SweepSections:SweepRangeImag(end);
                                        elseif H(p,8) == 2
                                            SweepRangeReal = SweepRangeReal(1):-FastFluidVelocity/SweepSections:SweepRangeReal(end);
                                            SweepRangeImag = SweepRangeImag(1):FastFluidVelocity/SweepSections:SweepRangeImag(end);
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
                            if  FluidLoading
                                if  (any(SweepRangeImag >= 0) && any(SweepRangeImag < 0)) || (numel(find(X(:,1))) > 1 && all(sign(SweepRangeImag) ~= sign(X(i-1,4))))
                                    cp = AngularFrequency./real(Wavenumber);
                                    if  SolutionType == 1
                                        Idx3 = SweepRangeImag >= 0;
                                    else
                                        Idx3 = SweepRangeImag < 0;
                                    end
                                    SignFluid1 = ones(YSize);
                                    if  Geometry == 2
                                        SignFluid1(cp <= FastFluidVelocity & Idx3) = -1;
                                    else
                                        Idx1 = cp <= SlowFluidVelocity;
                                        Idx2 = cp <= FastFluidVelocity & ~Idx1;
                                        SignFluid2 = ones(YSize);
                                        SignFluid1(Idx1 & Idx3) = -1;
                                        SignFluid2(Idx1 & Idx3) = -1;
                                        if  Fluid1.Velocity > Fluid2.Velocity
                                            SignFluid1(Idx2 & Idx3) = -1;
                                        elseif Fluid1.Velocity < Fluid2.Velocity
                                            SignFluid2(Idx2 & Idx3) = -1;
                                        end
                                    end
                                    if  k == 1 && ~any(X(:,1)) && ~H(p,2)
                                        if  H(p,8) < 2
                                            if  SolutionType == 1
                                                Idx = SweepRangeImag > 1.5*SweepRangeImagStep;
                                            else
                                                Idx = Idx3;
                                            end
                                        else
                                            Wavenumber(cp > FastFluidVelocity) = NaN;
                                            if  ~ContinuousAcrossSlowFluidVelocity
                                                Border = false(YSize);
                                                dm = cp < SlowFluidVelocity;
                                                for j = 1:size(Border,1)
                                                    if  SolutionType == 1
                                                        l = find(dm(j,:),1,'last');
                                                        if  l ~= YSize(2)
                                                            Border(j,l) = true;
                                                        end
                                                    else
                                                        l = find(dm(j,:),1);
                                                        if  l ~= 1
                                                            Border(j,l) = true;
                                                        end
                                                    end
                                                end
                                                for j = 1:size(Border,2)
                                                    l = find(dm(:,j),1);
                                                    if  l ~= 1
                                                        Border(l,j) = true;
                                                    end
                                                    l = find(dm(:,j),1,'last');
                                                    if  l ~= YSize(1)
                                                        Border(l,j) = true;
                                                    end
                                                end
                                                Wavenumber(Border & Idx3) = NaN;
                                            end
                                            if  SolutionType == 1
                                                Wavenumber(imag(Wavenumber) < -real(Wavenumber)/2) = NaN; % exclude strongly damped backward Scholte modes
                                                Idx = ~Idx3 | SweepRangeImag > SweepRangeImag(ceil(.69*length(SweepRangeImag)));
                                            else
                                                Wavenumber(imag(Wavenumber) > real(Wavenumber)/2) = NaN; % exclude strongly damped forward Scholte modes
                                                Idx = SweepRangeImag > 1.5*SweepRangeImagStep | SweepRangeImag < SweepRangeImag(ceil(.31*length(SweepRangeImag)));
                                            end
                                        end
                                        Wavenumber(:,Idx) = []; % cut off excessive NaNs
                                        SweepRangeImag(Idx) = [];
                                        SignFluid1(:,Idx) = [];
                                        if  Geometry ~= 2
                                            SignFluid2(:,Idx) = [];
                                        end
                                        YSize = size(Wavenumber);
                                    end
                                else
                                    if  ~any(X(:,1))
                                        if  H(p,2)
                                            SignFluid1 = H(p,6)*ones(YSize);
                                            if  Geometry ~= 2
                                                SignFluid2 = H(p,7)*ones(YSize);
                                            end
                                        else
                                            SignFluid1 = SignFluid1(MIN(1),MIN(2))*ones(YSize);
                                            if  Geometry ~= 2
                                                SignFluid2 = SignFluid2(MIN(1),MIN(2))*ones(YSize);
                                            end
                                        end
                                    else
                                        SignFluid1 = X(i-1,5)*ones(YSize);
                                        if  Geometry ~= 2
                                            SignFluid2 = X(i-1,6)*ones(YSize);
                                        end
                                    end
                                end
                            end
                            if  k == 1 && ~any(X(:,1)) && ~H(p,2)
                                Wavenumber(abs(imag(Wavenumber)) > AttenuationLimit*real(Wavenumber)) = NaN; % exclude strongly damped modes
                            end
                            for j = 1:size(Neighbors,1) % ignore previously traced solutions
                                Wavenumber(SweepRangeReal > Neighbors(j,1)-SweepRangeRealStep/2 & SweepRangeReal < Neighbors(j,1)+SweepRangeRealStep/2,SweepRangeImag > Neighbors(j,2)-SweepRangeImagStep/2 & SweepRangeImag < Neighbors(j,2)+SweepRangeImagStep/2) = NaN;
                            end
                            if  Geometry == 1
                                if  MaterialType == 1
                                    [Y,Yc] = Computer_Plate(ModeFamily,Wavenumber,Lambda,Mu,Half,kL2,kT2,kF12,kF22,gF1,gF2,SignFluid1,SignFluid2,Density1,ToggleFluid1,ToggleFluid2,YSize);
                                else
                                    if  k == 1 && ~any(X(:,1)) && ~H(p,2) && MatrixMethods == 1
                                        Idx = SweepRangeReal > Limit;
                                        Num = numel(find(Idx));
                                        if  FluidLoading
                                            SignFluid11 = SignFluid1(Idx,:);
                                            SignFluid21 = SignFluid2(Idx,:);
                                            SignFluid12 = SignFluid1(~Idx,:);
                                            SignFluid22 = SignFluid2(~Idx,:);
                                        end
                                        if  ~Decoupled
                                            if  Num > 0
                                                Y1 = Computer_Coupled(ModeFamily,1,Wavenumber(Idx,:),LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleFluid1,ToggleFluid2,kF12,kF22,gF1,gF2,SignFluid11,SignFluid21,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34,[Num YSize(2)]);
                                            end
                                            Y2 = Computer_Coupled(ModeFamily,2,Wavenumber(~Idx,:),LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleFluid1,ToggleFluid2,kF12,kF22,gF1,gF2,SignFluid12,SignFluid22,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34,[YSize(1)-Num YSize(2)]);
                                        else
                                            if  Num > 0
                                                Y1 = Computer_Decoupled(ModeFamily,1,Wavenumber(Idx,:),LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleFluid1,ToggleFluid2,kF12,kF22,gF1,gF2,SignFluid11,SignFluid21,c,rw2,A1,a21,a31,b22,b32,b33,[Num YSize(2)]);
                                            end
                                            Y2 = Computer_Decoupled(ModeFamily,2,Wavenumber(~Idx,:),LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleFluid1,ToggleFluid2,kF12,kF22,gF1,gF2,SignFluid12,SignFluid22,c,rw2,A1,a21,a31,b22,b32,b33,[YSize(1)-Num YSize(2)]);
                                        end
                                        if  Num > 0
                                            Y = [Y1;Y2];
                                            Y(Num,:) = NaN;
                                        else
                                            Y = Y2;
                                        end
                                    else
                                        if  ~any(X(:,1)) && MatrixMethods == 1
                                            if  AngularFrequency/real(Wavenumber(3,3)) > Limit
                                                MatrixMethod = 1; % TMM
                                            else
                                                MatrixMethod = 2; % SMM
                                            end
                                        end
                                        if  ~Decoupled
                                            Y = Computer_Coupled(ModeFamily,MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleFluid1,ToggleFluid2,kF12,kF22,gF1,gF2,SignFluid1,SignFluid2,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34,YSize);
                                        else
                                            Y = Computer_Decoupled(ModeFamily,MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleFluid1,ToggleFluid2,kF12,kF22,gF1,gF2,SignFluid1,SignFluid2,c,rw2,A1,a21,a31,b22,b32,b33,YSize);
                                        end
                                    end
                                end
                            elseif Geometry == 2
                                if  n == 0
                                    [Y,Yc] = Computer_Rod_L(Wavenumber,Ro,Ro2,kL2,kT2,kF12,SignFluid1,Density1,ToggleFluid1);
                                else
                                    [Y,Yc] = Computer_Rod_F(Wavenumber,Ro,Ro2,kL2,kT2,kF12,SignFluid1,Density1,ToggleFluid1,YSize,n,n2);
                                end
                            elseif Geometry == 3
                                if  n == 0
                                    [Y,Yc] = Computer_Pipe_L(Wavenumber,Ro,Ri,Ro2,Ri2,kL2,kT2,kF12,kF22,SignFluid1,SignFluid2,Density1,Density2,ToggleFluid1,ToggleFluid2,Sink,YSize);
                                else
                                    [Y,Yc] = Computer_Pipe_F(Wavenumber,Ro,Ri,Ro2,Ri2,kL2,kT2,kF12,kF22,SignFluid1,SignFluid2,Density1,Density2,ToggleFluid1,ToggleFluid2,Sink,YSize,n,n2);
                                end
                            end
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
                                                if  Geometry == 2
                                                    if  FluidLoading
                                                        H = [H(1:p,:);Frequency(i) SweepRangeReal(Min(j,1)) SweepRangeImag(Min(j,2)) SweepRangeRealStep SweepRangeImagStep SignFluid1(Min(j,1),Min(j,2)) 0 H(p,8);H(p+1:end,:)];
                                                    else
                                                        H = [H(1:p,:);Frequency(i) SweepRangeReal(Min(j,1)) SweepRangeImag(Min(j,2)) SweepRangeRealStep SweepRangeImagStep SignFluid1 0 H(p,8);H(p+1:end,:)];
                                                    end
                                                else
                                                    if  FluidLoading
                                                        H = [H(1:p,:);Frequency(i) SweepRangeReal(Min(j,1)) SweepRangeImag(Min(j,2)) SweepRangeRealStep SweepRangeImagStep SignFluid1(Min(j,1),Min(j,2)) SignFluid2(Min(j,1),Min(j,2)) H(p,8);H(p+1:end,:)];
                                                    else
                                                        H = [H(1:p,:);Frequency(i) SweepRangeReal(Min(j,1)) SweepRangeImag(Min(j,2)) SweepRangeRealStep SweepRangeImagStep SignFluid1 SignFluid2 H(p,8);H(p+1:end,:)];
                                                    end
                                                end
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
% disp(['Type = ',num2str(SolutionType),' p = ',num2str(p),' r = ',num2str(r),' f = ',num2str(Frequency(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Outlier'])
                            else
                                if  MaterialType == 1 && size(Min,1) == 1 && MIN(1) > 1 && MIN(2) > 1 && MIN(1) < YSize(1) && MIN(2) < YSize(2) && all(abs(diff(angle([Yc(MIN(1)-1,MIN(2)) Yc(MIN(1),MIN(2)-1);Yc(MIN(1)+1,MIN(2)) Yc(MIN(1),MIN(2)+1)]))) < pi/2) % check that there is a sufficient change in phase; only then is the minimum a zero and thereby a valid modal solution
                                    Misses(i) = true;
                                end
                                if  Geometry == 2
                                    if  FluidLoading
                                        X(i,:) = [cp imag(Wavenumber(MIN(1),MIN(2))) SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2)) SignFluid1(MIN(1),MIN(2)) 0]; % phase velocity (m/s), attenuation (Np/m), real velocity (m/s), imaginary velocity (m/s), growing (leaky wave)/decaying (Scholte wave) out-of-plane wavenumber component in the outer fluid (+1/-1), ...in the inner fluid                              
                                    else
                                        X(i,:) = [cp imag(Wavenumber(MIN(1),MIN(2))) SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2)) SignFluid1 0];
                                    end
                                else
                                    if  FluidLoading
                                        X(i,:) = [cp imag(Wavenumber(MIN(1),MIN(2))) SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2)) SignFluid1(MIN(1),MIN(2)) SignFluid2(MIN(1),MIN(2))];
                                    else
                                        X(i,:) = [cp imag(Wavenumber(MIN(1),MIN(2))) SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2)) SignFluid1 SignFluid2];
                                    end
                                end
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
                        X(i,:) = [Fit1(Frequency(i)) Fit2(Frequency(i)) Fit3(Frequency(i)) Fit4(Frequency(i)) X(i-1,5:6)];
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
                    (H(p,8) && BelowCutoff(i)) ||...
                    X(i) > PhaseVelocityLimit2
                    break
                elseif i > 1 && any(X(i-1,5:6) == -1) && (...
                    (~ContinuousAcrossSlowFluidVelocity && ((X(i-1) < SlowFluidVelocity && X(i) > SlowFluidVelocity) || (X(i-1) > SlowFluidVelocity && X(i) < SlowFluidVelocity))) ||...
                    (~ContinuousAcrossFastFluidVelocity && X(i) > FastFluidVelocity))
                    X(i,:) = [];
                    Misses(i) = [];
                    Frequency(i) = [];
                    Gear(i) = [];
                    BelowCutoff(i) = [];
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
% String = ['Type = ',num2str(SolutionType),' p = ',num2str(p),' r = ',num2str(r),' f = ',num2str(Frequency(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i)
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i)
%     disp(['Type = ',num2str(SolutionType),' p = ',num2str(p),' r = ',num2str(r),' f = ',num2str(Frequency(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
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
            A{end+1}(:,[1:4 7:11]) = [Frequency'*[1 1e-3 ScaleFactor] fillmissing(X(:,1:4),'spline') X(:,5:6)]; % add the dispersion curve to the output data container
            Fit{1,end+1} = fit(A{end}(:,1),A{end}(:,8),'cubicspline'); % fit the dispersion curve to query interpolated data to be ignored in the tracing of the next dispersion curves
            Fit{2,end} = fit(A{end}(:,1),A{end}(:,9),'cubicspline');
            if  Multithreading
                send(Q(r),A{end}(:,[1 4]))
            else
                line(ax,A{end}(:,1),A{end}(:,4)/1e3,'color',LineColor)
                clearpoints(g(1))
                clearpoints(g(2))
            end
        end
    end
end
for p = 1:length(A) % post processing
    A{p}(A{p}(:,4) < 0,:) = [];
    A{p}(:,4) = A{p}(:,4)/1e3;
    A{p}(:,8:9) = [];
end
% figure,hold on
% for p = 1:length(A)
%     plot(A{p}(:,1),A{p}(:,4),'k');
% end
% line(FrequencyRange,MatrixMethodLimit/1e3,'color','r')
end
function [Y,Yc] = Computer_Plate(ModeFamily,k,Lambda,Mu,Half,kL2,kT2,kFu2,kFl2,gFu,gFl,SignUpperFluid,SignLowerFluid,Densityu,ToggleUpperFluid,ToggleLowerFluid,YSize)
    k2 = k.^2;
    x2 = kL2-k2;
    y2 = kT2-k2;
    x = sqrt(x2);
    y = sqrt(y2);
    xH = x*Half;
    yH = y*Half;
    if  ModeFamily == 1
        if  ToggleUpperFluid
            Yc = (y2-k2).^2./y./tan(xH)+4*k2.*x./tan(yH)-1i*SignUpperFluid*Densityu*kT2^2.*x./y./sqrt(kFu2-k2);
        else
            Yc = (y2-k2).^2./y./tan(xH)+4*k2.*x./tan(yH);
        end
    elseif ModeFamily == 2
        if  ToggleUpperFluid
            Yc = (y2-k2).^2./y.*tan(xH)+4*k2.*x.*tan(yH)+1i*SignUpperFluid*Densityu*kT2^2.*x./y./sqrt(kFu2-k2);
        else
            Yc = (y2-k2).^2./y.*tan(xH)+4*k2.*x.*tan(yH);
        end
    else
        SinxH = sin(xH);
        SinyH = sin(yH);
        CosxH = cos(xH);
        CosyH = cos(yH);
        Sin_xH = sin(-xH);
        Sin_yH = sin(-yH);
        Cos_xH = cos(-xH);
        Cos_yH = cos(-yH);
        a0 = 2i*Mu*k;
        a1 = -Lambda*kL2-2*Mu*x2;
        a2 = Mu*(k2-y2);
        a3 = a0.*x;
        a4 = a0.*y;
        a5 = -1i*k;
        Yc = NaN(YSize)+1i*NaN;
        if  ToggleUpperFluid && ToggleLowerFluid
            zu = 1i*SignUpperFluid.*sqrt(kFu2-k2);
            zl = 1i*SignLowerFluid.*sqrt(kFl2-k2);
            E = exp(zu*Half);
            E_ = exp(-zl*Half);
            for l = 1:YSize(2)
                for j = 1:YSize(1)
                    if  ~isnan(k(j,l))
                        M(1,1) = a1(j,l)*SinxH(j,l);
                        M(1,2) = a1(j,l)*CosxH(j,l);
                        M(1,3) = -a4(j,l)*CosyH(j,l);
                        M(1,4) = a4(j,l)*SinyH(j,l);
                        M(1,5) = gFu*E(j,l);
                        M(2,1) = a3(j,l)*CosxH(j,l);
                        M(2,2) = -a3(j,l)*SinxH(j,l);
                        M(2,3) = a2(j,l)*SinyH(j,l);
                        M(2,4) = a2(j,l)*CosyH(j,l);
                        M(3,1) = x(j,l)*CosxH(j,l);
                        M(3,2) = -x(j,l)*SinxH(j,l);
                        M(3,3) = a5(j,l)*SinyH(j,l);
                        M(3,4) = a5(j,l)*CosyH(j,l);
                        M(3,5) = -zu(j,l)*E(j,l);
                        M(4,1) = a1(j,l)*Sin_xH(j,l);
                        M(4,2) = a1(j,l)*Cos_xH(j,l);
                        M(4,3) = -a4(j,l)*Cos_yH(j,l);
                        M(4,4) = a4(j,l)*Sin_yH(j,l);
                        M(4,6) = gFl*E_(j,l);
                        M(5,1) = a3(j,l)*Cos_xH(j,l);
                        M(5,2) = -a3(j,l)*Sin_xH(j,l);
                        M(5,3) = a2(j,l)*Sin_yH(j,l);
                        M(5,4) = a2(j,l)*Cos_yH(j,l);
                        M(6,1) = x(j,l)*Cos_xH(j,l);
                        M(6,2) = -x(j,l)*Sin_xH(j,l);
                        M(6,3) = a5(j,l)*Sin_yH(j,l);
                        M(6,4) = a5(j,l)*Cos_yH(j,l);
                        M(6,6) = zl(j,l)*E_(j,l);
                        Yc(j,l) = det(M);
                    end
                end
            end
        else
            if  ToggleUpperFluid
                z = 1i*SignUpperFluid.*sqrt(kFu2-k2);
                gF = gFu;
            else
                z = 1i*SignLowerFluid.*sqrt(kFl2-k2);
                gF = gFl;
            end
            E = exp(z*Half);
            for l = 1:YSize(2)
                for j = 1:YSize(1)
                    if  ~isnan(k(j,l))
                        M(1,1) = a1(j,l)*SinxH(j,l);
                        M(1,2) = a1(j,l)*CosxH(j,l);
                        M(1,3) = -a4(j,l)*CosyH(j,l);
                        M(1,4) = a4(j,l)*SinyH(j,l);
                        M(1,5) = gF*E(j,l);
                        M(2,1) = a3(j,l)*CosxH(j,l);
                        M(2,2) = -a3(j,l)*SinxH(j,l);
                        M(2,3) = a2(j,l)*SinyH(j,l);
                        M(2,4) = a2(j,l)*CosyH(j,l);
                        M(3,1) = x(j,l)*CosxH(j,l);
                        M(3,2) = -x(j,l)*SinxH(j,l);
                        M(3,3) = a5(j,l)*SinyH(j,l);
                        M(3,4) = a5(j,l)*CosyH(j,l);
                        M(3,5) = -z(j,l)*E(j,l);
                        M(4,1) = a1(j,l)*Sin_xH(j,l);
                        M(4,2) = a1(j,l)*Cos_xH(j,l);
                        M(4,3) = -a4(j,l)*Cos_yH(j,l);
                        M(4,4) = a4(j,l)*Sin_yH(j,l);
                        M(5,1) = a3(j,l)*Cos_xH(j,l);
                        M(5,2) = -a3(j,l)*Sin_xH(j,l);
                        M(5,3) = a2(j,l)*Sin_yH(j,l);
                        M(5,4) = a2(j,l)*Cos_yH(j,l);
                        Yc(j,l) = det(M);
                    end
                end
            end
        end
    end
    Y = abs(Yc);
end
function [Y,Yc] = Computer_Rod_F(k,R,R2,kL2,kT2,kF2,SignFluid,Density,ToggleFluid,YSize,n,n2)
    if  n == 1 && ~ToggleFluid
        k2 = k.^2*R2;
        k4 = k2.^2;
        x = sqrt(kL2*R2-k2);
        y2 = kT2*R2-k2;
        y = sqrt(y2);
        y4 = y2.^2;
        y6 = y2.^3;
        Zx = x.*besselj(0,x)./besselj(1,x);
        Zy = y.*besselj(0,y)./besselj(1,y);
        Zy2 = Zy.^2;
        a1 = y2.*k2;
        a2 = y4.*k2;
        a3 = y2.*k4;
        A1 = 2*(y2-k2).^2;
        A2 = 2*y4+10*a1;
        A3 = y6-10*y4-2*a2+2*a1+a3-4*k4;
        A4 = 4*a2-2*y4-18*a1;
        A5 = -y6+8*y4-2*a2+8*a1-a3;
        Yc = A1+A2.*Zx./Zy+A3./Zy+A4.*Zx./Zy2+A5./Zy2;
    else
        k2 = k.^2;
        y2 = kT2-k2;
        xR = sqrt(kL2-k2)*R;
        yR = sqrt(y2)*R;
        Jnx = besselj(n,xR);
        Jny = besselj(n,yR);
        dJnxR = n*Jnx-xR.*besselj(n+1,xR);
        dJnyR = n*Jny-yR.*besselj(n+1,yR);
        a1 = kT2/2*R2;
        if  ToggleFluid
            zR = SignFluid.*sqrt(kF2-k2)*R;
            Hnz = besselh(n,zR);
            dHnzR = n*Hnz-zR.*besselh(n+1,zR);
            a2 = a1*Density*Hnz;
        end
        Yc = NaN(YSize)+1i*NaN;
        for l = 1:YSize(2)
            for j = 1:YSize(1)
                if  ~isnan(k(j,l))
                    if  ToggleFluid
                        M(1,1) = dJnxR(j,l)+(a1-k2(j,l)*R2-n2)*Jnx(j,l);
                        M(1,2) = -k(j,l)*(dJnyR(j,l)+(y2(j,l)*R2-n2)*Jny(j,l));
                        M(1,3) = n*(dJnyR(j,l)-Jny(j,l));
                        M(1,4) = a2(j,l);
                        M(2,1) = 2*n*(dJnxR(j,l)-Jnx(j,l));
                        M(2,2) = 2*k(j,l)*n*(Jny(j,l)-dJnyR(j,l));
                        M(2,3) = 2*dJnyR(j,l)+(y2(j,l)*R2-2*n2)*Jny(j,l);
                        M(3,1) = 2*k(j,l)*dJnxR(j,l);
                        M(3,2) = (y2(j,l)-k2(j,l))*dJnyR(j,l);
                        M(3,3) = -k(j,l)*n*Jny(j,l);
                        M(4,1) = dJnxR(j,l);
                        M(4,2) = -k(j,l)*dJnyR(j,l);
                        M(4,3) = -n*Jny(j,l);
                        M(4,4) = dHnzR(j,l);
                    else
                        M(1,1) = dJnxR(j,l)+(a1-k2(j,l)*R2-n2)*Jnx(j,l);
                        M(1,2) = -k(j,l)*(dJnyR(j,l)+(y2(j,l)*R2-n2)*Jny(j,l));
                        M(1,3) = n*(dJnyR(j,l)-Jny(j,l));
                        M(2,1) = 2*n*(dJnxR(j,l)-Jnx(j,l));
                        M(2,2) = 2*k(j,l)*n*(Jny(j,l)-dJnyR(j,l));
                        M(2,3) = 2*dJnyR(j,l)+(y2(j,l)*R2-2*n2)*Jny(j,l);
                        M(3,1) = 2*k(j,l)*dJnxR(j,l);
                        M(3,2) = (y2(j,l)-k2(j,l))*dJnyR(j,l);
                        M(3,3) = -k(j,l)*n*Jny(j,l);
                    end
                    Yc(j,l) = det(M);
                end
            end
        end
    end
    Y = abs(Yc);
end
function [Y,Yc] = Computer_Rod_L(k,R,R2,kL2,kT2,kF2,SignFluid,Density,ToggleFluid)
    k2 = k.^2;
    y2 = kT2-k2;
    x = sqrt(kL2-k2);
    y = sqrt(y2);
    xR = x*R;
    yR = y*R;
    if  ToggleFluid
        zR = SignFluid.*sqrt(kF2-k2)*R;
        J1xR = xR.*besselj(1,xR);
        J1yR = yR.*besselj(1,yR);
        H1zR = zR.*besselh(1,zR);
        Yc = H1zR*R2.*((k2-y2).^2/2.*J1yR.*besselj(0,xR)+2*k2.*y2.*J1xR.*besselj(0,yR))-(H1zR+kT2/2*R2*Density*besselh(0,zR))*kT2.*J1xR.*J1yR;
    else
        Yc = 2*x/R.*kT2-(y2-k2).^2.*besselj(0,xR)./besselj(1,xR)-4*k2.*x.*y.*besselj(0,yR)./besselj(1,yR);
    end
    Y = abs(Yc);
end
function [Y,Yc] = Computer_Pipe_F(k,Ro,Ri,Ro2,Ri2,kL2,kT2,kFo2,kFi2,SignOuterFluid,SignInnerFluid,Densityo,Densityi,ToggleOuterFluid,ToggleInnerFluid,Sink,YSize,n,n2)
    k2 = k.^2;
    y2 = kT2-k2;
    x = sqrt(k2-kL2);
    y = sqrt(-y2);
    xRo = x*Ro;
    yRo = y*Ro;
    xRi = x*Ri;
    yRi = y*Ri;
    Znxi = besseli(n,xRi);
    Znxo = besseli(n,xRo);
    Znyi = besseli(n,yRi);
    Znyo = besseli(n,yRo);
    Wnxi = besselk(n,xRi);
    Wnxo = besselk(n,xRo);
    Wnyi = besselk(n,yRi);
    Wnyo = besselk(n,yRo);
    dZnxiRi = n*Znxi+xRi.*besseli(n+1,xRi);
    dZnxoRo = n*Znxo+xRo.*besseli(n+1,xRo);
    dZnyiRi = n*Znyi+yRi.*besseli(n+1,yRi);
    dZnyoRo = n*Znyo+yRo.*besseli(n+1,yRo);
    dWnxiRi = n*Wnxi-xRi.*besselk(n+1,xRi);
    dWnxoRo = n*Wnxo-xRo.*besselk(n+1,xRo);
    dWnyiRi = n*Wnyi-yRi.*besselk(n+1,yRi);
    dWnyoRo = n*Wnyo-yRo.*besselk(n+1,yRo);
    a1 = kT2/2*Ri2;
    a2 = kT2/2*Ro2;
    if  ToggleInnerFluid
        zRi = SignInnerFluid.*sqrt(kFi2-k2)*Ri;
        if  Sink
            Znzi = besselh(n,2,zRi);
            dZnziRi = n*Znzi-zRi.*besselh(n+1,2,zRi);
        else
            Znzi = besselj(n,zRi);
            dZnziRi = n*Znzi-zRi.*besselj(n+1,zRi);
        end
        a3 = a1*Densityi*Znzi;
    end
    if  ToggleOuterFluid
        zRo = SignOuterFluid.*sqrt(kFo2-k2)*Ro;
        Hnzo = besselh(n,zRo);
        dHnzoRo = n*Hnzo-zRo.*besselh(n+1,zRo);
        a4 = a2*Densityo*Hnzo;
    end
    Yc = NaN(YSize)+1i*NaN;
    for l = 1:YSize(2)
        for j = 1:YSize(1)
            if  ~isnan(k(j,l))
                if  ToggleInnerFluid && ToggleOuterFluid
                    M(1,1) = dZnxiRi(j,l)+(a1-k2(j,l)*Ri2-n2)*Znxi(j,l);
                    M(1,2) = dWnxiRi(j,l)+(a1-k2(j,l)*Ri2-n2)*Wnxi(j,l);
                    M(1,3) = -k(j,l)*(dZnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Znyi(j,l));
                    M(1,4) = -k(j,l)*(dWnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Wnyi(j,l));
                    M(1,5) = n*(dZnyiRi(j,l)-Znyi(j,l));
                    M(1,6) = n*(dWnyiRi(j,l)-Wnyi(j,l));
                    M(1,7) = a3(j,l);
                    M(2,1) = 2*n*(dZnxiRi(j,l)-Znxi(j,l));
                    M(2,2) = 2*n*(dWnxiRi(j,l)-Wnxi(j,l));
                    M(2,3) = 2*k(j,l)*n*(Znyi(j,l)-dZnyiRi(j,l));
                    M(2,4) = 2*k(j,l)*n*(Wnyi(j,l)-dWnyiRi(j,l));
                    M(2,5) = 2*dZnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Znyi(j,l);
                    M(2,6) = 2*dWnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Wnyi(j,l);
                    M(3,1) = 2*k(j,l)*dZnxiRi(j,l);
                    M(3,2) = 2*k(j,l)*dWnxiRi(j,l);
                    M(3,3) = (y2(j,l)-k2(j,l))*dZnyiRi(j,l);
                    M(3,4) = (y2(j,l)-k2(j,l))*dWnyiRi(j,l);
                    M(3,5) = -k(j,l)*n*Znyi(j,l);
                    M(3,6) = -k(j,l)*n*Wnyi(j,l);
                    M(4,1) = dZnxiRi(j,l);
                    M(4,2) = dWnxiRi(j,l);
                    M(4,3) = -k(j,l)*dZnyiRi(j,l);
                    M(4,4) = -k(j,l)*dWnyiRi(j,l);
                    M(4,5) = -n*Znyi(j,l);
                    M(4,6) = -n*Wnyi(j,l);
                    M(4,7) = dZnziRi(j,l);
                    M(5,1) = dZnxoRo(j,l)+(a2-k2(j,l)*Ro2-n2)*Znxo(j,l);
                    M(5,2) = dWnxoRo(j,l)+(a2-k2(j,l)*Ro2-n2)*Wnxo(j,l);
                    M(5,3) = -k(j,l)*(dZnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Znyo(j,l));
                    M(5,4) = -k(j,l)*(dWnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Wnyo(j,l));
                    M(5,5) = n*(dZnyoRo(j,l)-Znyo(j,l));
                    M(5,6) = n*(dWnyoRo(j,l)-Wnyo(j,l));
                    M(5,8) = a4(j,l);
                    M(6,1) = 2*n*(dZnxoRo(j,l)-Znxo(j,l));
                    M(6,2) = 2*n*(dWnxoRo(j,l)-Wnxo(j,l));
                    M(6,3) = 2*k(j,l)*n*(Znyo(j,l)-dZnyoRo(j,l));
                    M(6,4) = 2*k(j,l)*n*(Wnyo(j,l)-dWnyoRo(j,l));
                    M(6,5) = 2*dZnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Znyo(j,l);
                    M(6,6) = 2*dWnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Wnyo(j,l);
                    M(7,1) = 2*k(j,l)*dZnxoRo(j,l);
                    M(7,2) = 2*k(j,l)*dWnxoRo(j,l);
                    M(7,3) = (y2(j,l)-k2(j,l))*dZnyoRo(j,l);
                    M(7,4) = (y2(j,l)-k2(j,l))*dWnyoRo(j,l);
                    M(7,5) = -k(j,l)*n*Znyo(j,l);
                    M(7,6) = -k(j,l)*n*Wnyo(j,l);
                    M(8,1) = dZnxoRo(j,l);
                    M(8,2) = dWnxoRo(j,l);
                    M(8,3) = -k(j,l)*dZnyoRo(j,l);
                    M(8,4) = -k(j,l)*dWnyoRo(j,l);
                    M(8,5) = -n*Znyo(j,l);
                    M(8,6) = -n*Wnyo(j,l);
                    M(8,8) = dHnzoRo(j,l);
                elseif ToggleInnerFluid && ~ToggleOuterFluid
                    M(1,1) = dZnxiRi(j,l)+(a1-k2(j,l)*Ri2-n2)*Znxi(j,l);
                    M(1,2) = dWnxiRi(j,l)+(a1-k2(j,l)*Ri2-n2)*Wnxi(j,l);
                    M(1,3) = -k(j,l)*(dZnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Znyi(j,l));
                    M(1,4) = -k(j,l)*(dWnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Wnyi(j,l));
                    M(1,5) = n*(dZnyiRi(j,l)-Znyi(j,l));
                    M(1,6) = n*(dWnyiRi(j,l)-Wnyi(j,l));
                    M(1,7) = a3(j,l);
                    M(2,1) = 2*n*(dZnxiRi(j,l)-Znxi(j,l));
                    M(2,2) = 2*n*(dWnxiRi(j,l)-Wnxi(j,l));
                    M(2,3) = 2*k(j,l)*n*(Znyi(j,l)-dZnyiRi(j,l));
                    M(2,4) = 2*k(j,l)*n*(Wnyi(j,l)-dWnyiRi(j,l));
                    M(2,5) = 2*dZnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Znyi(j,l);
                    M(2,6) = 2*dWnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Wnyi(j,l);
                    M(3,1) = 2*k(j,l)*dZnxiRi(j,l);
                    M(3,2) = 2*k(j,l)*dWnxiRi(j,l);
                    M(3,3) = (y2(j,l)-k2(j,l))*dZnyiRi(j,l);
                    M(3,4) = (y2(j,l)-k2(j,l))*dWnyiRi(j,l);
                    M(3,5) = -k(j,l)*n*Znyi(j,l);
                    M(3,6) = -k(j,l)*n*Wnyi(j,l);
                    M(4,1) = dZnxiRi(j,l);
                    M(4,2) = dWnxiRi(j,l);
                    M(4,3) = -k(j,l)*dZnyiRi(j,l);
                    M(4,4) = -k(j,l)*dWnyiRi(j,l);
                    M(4,5) = -n*Znyi(j,l);
                    M(4,6) = -n*Wnyi(j,l);
                    M(4,7) = dZnziRi(j,l);
                    M(5,1) = dZnxoRo(j,l)+(a2-k2(j,l)*Ro2-n2)*Znxo(j,l);
                    M(5,2) = dWnxoRo(j,l)+(a2-k2(j,l)*Ro2-n2)*Wnxo(j,l);
                    M(5,3) = -k(j,l)*(dZnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Znyo(j,l));
                    M(5,4) = -k(j,l)*(dWnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Wnyo(j,l));
                    M(5,5) = n*(dZnyoRo(j,l)-Znyo(j,l));
                    M(5,6) = n*(dWnyoRo(j,l)-Wnyo(j,l));
                    M(6,1) = 2*n*(dZnxoRo(j,l)-Znxo(j,l));
                    M(6,2) = 2*n*(dWnxoRo(j,l)-Wnxo(j,l));
                    M(6,3) = 2*k(j,l)*n*(Znyo(j,l)-dZnyoRo(j,l));
                    M(6,4) = 2*k(j,l)*n*(Wnyo(j,l)-dWnyoRo(j,l));
                    M(6,5) = 2*dZnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Znyo(j,l);
                    M(6,6) = 2*dWnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Wnyo(j,l);
                    M(7,1) = 2*k(j,l)*dZnxoRo(j,l);
                    M(7,2) = 2*k(j,l)*dWnxoRo(j,l);
                    M(7,3) = (y2(j,l)-k2(j,l))*dZnyoRo(j,l);
                    M(7,4) = (y2(j,l)-k2(j,l))*dWnyoRo(j,l);
                    M(7,5) = -k(j,l)*n*Znyo(j,l);
                    M(7,6) = -k(j,l)*n*Wnyo(j,l);
                elseif ~ToggleInnerFluid && ToggleOuterFluid
                    M(1,1) = dZnxiRi(j,l)+(a1-k2(j,l)*Ri2-n2)*Znxi(j,l);
                    M(1,2) = dWnxiRi(j,l)+(a1-k2(j,l)*Ri2-n2)*Wnxi(j,l);
                    M(1,3) = -k(j,l)*(dZnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Znyi(j,l));
                    M(1,4) = -k(j,l)*(dWnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Wnyi(j,l));
                    M(1,5) = n*(dZnyiRi(j,l)-Znyi(j,l));
                    M(1,6) = n*(dWnyiRi(j,l)-Wnyi(j,l));
                    M(2,1) = 2*n*(dZnxiRi(j,l)-Znxi(j,l));
                    M(2,2) = 2*n*(dWnxiRi(j,l)-Wnxi(j,l));
                    M(2,3) = 2*k(j,l)*n*(Znyi(j,l)-dZnyiRi(j,l));
                    M(2,4) = 2*k(j,l)*n*(Wnyi(j,l)-dWnyiRi(j,l));
                    M(2,5) = 2*dZnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Znyi(j,l);
                    M(2,6) = 2*dWnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Wnyi(j,l);
                    M(3,1) = 2*k(j,l)*dZnxiRi(j,l);
                    M(3,2) = 2*k(j,l)*dWnxiRi(j,l);
                    M(3,3) = (y2(j,l)-k2(j,l))*dZnyiRi(j,l);
                    M(3,4) = (y2(j,l)-k2(j,l))*dWnyiRi(j,l);
                    M(3,5) = -k(j,l)*n*Znyi(j,l);
                    M(3,6) = -k(j,l)*n*Wnyi(j,l);
                    M(4,1) = dZnxoRo(j,l)+(a2-k2(j,l)*Ro2-n2)*Znxo(j,l);
                    M(4,2) = dWnxoRo(j,l)+(a2-k2(j,l)*Ro2-n2)*Wnxo(j,l);
                    M(4,3) = -k(j,l)*(dZnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Znyo(j,l));
                    M(4,4) = -k(j,l)*(dWnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Wnyo(j,l));
                    M(4,5) = n*(dZnyoRo(j,l)-Znyo(j,l));
                    M(4,6) = n*(dWnyoRo(j,l)-Wnyo(j,l));
                    M(4,7) = a4(j,l);
                    M(5,1) = 2*n*(dZnxoRo(j,l)-Znxo(j,l));
                    M(5,2) = 2*n*(dWnxoRo(j,l)-Wnxo(j,l));
                    M(5,3) = 2*k(j,l)*n*(Znyo(j,l)-dZnyoRo(j,l));
                    M(5,4) = 2*k(j,l)*n*(Wnyo(j,l)-dWnyoRo(j,l));
                    M(5,5) = 2*dZnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Znyo(j,l);
                    M(5,6) = 2*dWnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Wnyo(j,l);
                    M(6,1) = 2*k(j,l)*dZnxoRo(j,l);
                    M(6,2) = 2*k(j,l)*dWnxoRo(j,l);
                    M(6,3) = (y2(j,l)-k2(j,l))*dZnyoRo(j,l);
                    M(6,4) = (y2(j,l)-k2(j,l))*dWnyoRo(j,l);
                    M(6,5) = -k(j,l)*n*Znyo(j,l);
                    M(6,6) = -k(j,l)*n*Wnyo(j,l);
                    M(7,1) = dZnxoRo(j,l);
                    M(7,2) = dWnxoRo(j,l);
                    M(7,3) = -k(j,l)*dZnyoRo(j,l);
                    M(7,4) = -k(j,l)*dWnyoRo(j,l);
                    M(7,5) = -n*Znyo(j,l);
                    M(7,6) = -n*Wnyo(j,l);
                    M(7,7) = dHnzoRo(j,l);
                elseif ~ToggleInnerFluid && ~ToggleOuterFluid
                    M(1,1) = dZnxiRi(j,l)+(a1-k2(j,l)*Ri2-n2)*Znxi(j,l);
                    M(1,2) = dWnxiRi(j,l)+(a1-k2(j,l)*Ri2-n2)*Wnxi(j,l);
                    M(1,3) = -k(j,l)*(dZnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Znyi(j,l));
                    M(1,4) = -k(j,l)*(dWnyiRi(j,l)+(y2(j,l)*Ri2-n2)*Wnyi(j,l));
                    M(1,5) = n*(dZnyiRi(j,l)-Znyi(j,l));
                    M(1,6) = n*(dWnyiRi(j,l)-Wnyi(j,l));
                    M(2,1) = 2*n*(dZnxiRi(j,l)-Znxi(j,l));
                    M(2,2) = 2*n*(dWnxiRi(j,l)-Wnxi(j,l));
                    M(2,3) = 2*k(j,l)*n*(Znyi(j,l)-dZnyiRi(j,l));
                    M(2,4) = 2*k(j,l)*n*(Wnyi(j,l)-dWnyiRi(j,l));
                    M(2,5) = 2*dZnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Znyi(j,l);
                    M(2,6) = 2*dWnyiRi(j,l)+(y2(j,l)*Ri2-2*n2)*Wnyi(j,l);
                    M(3,1) = 2*k(j,l)*dZnxiRi(j,l);
                    M(3,2) = 2*k(j,l)*dWnxiRi(j,l);
                    M(3,3) = (y2(j,l)-k2(j,l))*dZnyiRi(j,l);
                    M(3,4) = (y2(j,l)-k2(j,l))*dWnyiRi(j,l);
                    M(3,5) = -k(j,l)*n*Znyi(j,l);
                    M(3,6) = -k(j,l)*n*Wnyi(j,l);
                    M(4,1) = dZnxoRo(j,l)+(a2-k2(j,l)*Ro2-n2)*Znxo(j,l);
                    M(4,2) = dWnxoRo(j,l)+(a2-k2(j,l)*Ro2-n2)*Wnxo(j,l);
                    M(4,3) = -k(j,l)*(dZnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Znyo(j,l));
                    M(4,4) = -k(j,l)*(dWnyoRo(j,l)+(y2(j,l)*Ro2-n2)*Wnyo(j,l));
                    M(4,5) = n*(dZnyoRo(j,l)-Znyo(j,l));
                    M(4,6) = n*(dWnyoRo(j,l)-Wnyo(j,l));
                    M(5,1) = 2*n*(dZnxoRo(j,l)-Znxo(j,l));
                    M(5,2) = 2*n*(dWnxoRo(j,l)-Wnxo(j,l));
                    M(5,3) = 2*k(j,l)*n*(Znyo(j,l)-dZnyoRo(j,l));
                    M(5,4) = 2*k(j,l)*n*(Wnyo(j,l)-dWnyoRo(j,l));
                    M(5,5) = 2*dZnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Znyo(j,l);
                    M(5,6) = 2*dWnyoRo(j,l)+(y2(j,l)*Ro2-2*n2)*Wnyo(j,l);
                    M(6,1) = 2*k(j,l)*dZnxoRo(j,l);
                    M(6,2) = 2*k(j,l)*dWnxoRo(j,l);
                    M(6,3) = (y2(j,l)-k2(j,l))*dZnyoRo(j,l);
                    M(6,4) = (y2(j,l)-k2(j,l))*dWnyoRo(j,l);
                    M(6,5) = -k(j,l)*n*Znyo(j,l);
                    M(6,6) = -k(j,l)*n*Wnyo(j,l);
                end
                Yc(j,l) = det(M);
            end
        end
    end
    Y = abs(Yc);
end
function [Y,Yc] = Computer_Pipe_L(k,Ro,Ri,Ro2,Ri2,kL2,kT2,kFo2,kFi2,SignOuterFluid,SignInnerFluid,Densityo,Densityi,ToggleOuterFluid,ToggleInnerFluid,Sink,YSize)
    k2 = k.^2;
    y2 = kT2-k2;
    x = sqrt(k2-kL2);
    y = sqrt(-y2);
    xRo = x*Ro;
    yRo = y*Ro;
    xRi = x*Ri;
    yRi = y*Ri;
    Z0xiRi2 = Ri2*besseli(0,xRi);
    Z0xoRo2 = Ro2*besseli(0,xRo);
    Z0yiRi2 = Ri2*besseli(0,yRi);
    Z0yoRo2 = Ro2*besseli(0,yRo);
    W0xiRi2 = Ri2*besselk(0,xRi);
    W0xoRo2 = Ro2*besselk(0,xRo);
    W0yiRi2 = Ri2*besselk(0,yRi);
    W0yoRo2 = Ro2*besselk(0,yRo);
    Z1xiRi = -xRi.*besseli(1,xRi);
    Z1xoRo = -xRo.*besseli(1,xRo);
    Z1yiRi = -yRi.*besseli(1,yRi);
    Z1yoRo = -yRo.*besseli(1,yRo);
    W1xiRi = xRi.*besselk(1,xRi);
    W1xoRo = xRo.*besselk(1,xRo);
    W1yiRi = yRi.*besselk(1,yRi);
    W1yoRo = yRo.*besselk(1,yRo);
    a1 = kT2/2;
    if  ToggleInnerFluid
        zRi = SignInnerFluid.*sqrt(kFi2-k2)*Ri;
        if  Sink
            a2 = a1*Ri2*Densityi*besselh(0,2,zRi);
            Z1ziRi = -zRi.*besselh(1,2,zRi);
        else
            a2 = a1*Ri2*Densityi*besselj(0,zRi);
            Z1ziRi = -zRi.*besselj(1,zRi);
        end
    end
    if  ToggleOuterFluid
        zRo = SignOuterFluid.*sqrt(kFo2-k2)*Ro;
        a3 = a1*Ro2*Densityo*besselh(0,zRo);
        H1zoRo = -zRo.*besselh(1,zRo);
    end
    Yc = NaN(YSize)+1i*NaN;
    for l = 1:YSize(2)
        for j = 1:YSize(1)
            if  ~isnan(k(j,l))
                if  ToggleInnerFluid && ToggleOuterFluid
                    M(1,1) = (a1-k2(j,l))*Z0xiRi2(j,l)-Z1xiRi(j,l);
                    M(1,2) = (a1-k2(j,l))*W0xiRi2(j,l)-W1xiRi(j,l);
                    M(1,3) = k(j,l)*(Z1yiRi(j,l)-y2(j,l)*Z0yiRi2(j,l));
                    M(1,4) = k(j,l)*(W1yiRi(j,l)-y2(j,l)*W0yiRi2(j,l));
                    M(1,5) = a2(j,l);
                    M(2,1) = -2*k(j,l)*Z1xiRi(j,l);
                    M(2,2) = -2*k(j,l)*W1xiRi(j,l);
                    M(2,3) = (k2(j,l)-y2(j,l))*Z1yiRi(j,l);
                    M(2,4) = (k2(j,l)-y2(j,l))*W1yiRi(j,l);
                    M(3,1) = -Z1xiRi(j,l);
                    M(3,2) = -W1xiRi(j,l);
                    M(3,3) = k(j,l)*Z1yiRi(j,l);
                    M(3,4) = k(j,l)*W1yiRi(j,l);
                    M(3,5) = Z1ziRi(j,l);
                    M(4,1) = (a1-k2(j,l))*Z0xoRo2(j,l)-Z1xoRo(j,l);
                    M(4,2) = (a1-k2(j,l))*W0xoRo2(j,l)-W1xoRo(j,l);
                    M(4,3) = k(j,l)*(Z1yoRo(j,l)-y2(j,l)*Z0yoRo2(j,l));
                    M(4,4) = k(j,l)*(W1yoRo(j,l)-y2(j,l)*W0yoRo2(j,l));
                    M(4,6) = a3(j,l);
                    M(5,1) = -2*k(j,l)*Z1xoRo(j,l);
                    M(5,2) = -2*k(j,l)*W1xoRo(j,l);
                    M(5,3) = (k2(j,l)-y2(j,l))*Z1yoRo(j,l);
                    M(5,4) = (k2(j,l)-y2(j,l))*W1yoRo(j,l);
                    M(6,1) = -Z1xoRo(j,l);
                    M(6,2) = -W1xoRo(j,l);
                    M(6,3) = k(j,l)*Z1yoRo(j,l);
                    M(6,4) = k(j,l)*W1yoRo(j,l);
                    M(6,6) = H1zoRo(j,l);
                elseif ToggleInnerFluid && ~ToggleOuterFluid
                    M(1,1) = (a1-k2(j,l))*Z0xiRi2(j,l)-Z1xiRi(j,l);
                    M(1,2) = (a1-k2(j,l))*W0xiRi2(j,l)-W1xiRi(j,l);
                    M(1,3) = k(j,l)*(Z1yiRi(j,l)-y2(j,l)*Z0yiRi2(j,l));
                    M(1,4) = k(j,l)*(W1yiRi(j,l)-y2(j,l)*W0yiRi2(j,l));
                    M(1,5) = a2(j,l);
                    M(2,1) = -2*k(j,l)*Z1xiRi(j,l);
                    M(2,2) = -2*k(j,l)*W1xiRi(j,l);
                    M(2,3) = (k2(j,l)-y2(j,l))*Z1yiRi(j,l);
                    M(2,4) = (k2(j,l)-y2(j,l))*W1yiRi(j,l);
                    M(3,1) = -Z1xiRi(j,l);
                    M(3,2) = -W1xiRi(j,l);
                    M(3,3) = k(j,l)*Z1yiRi(j,l);
                    M(3,4) = k(j,l)*W1yiRi(j,l);
                    M(3,5) = Z1ziRi(j,l);
                    M(4,1) = (a1-k2(j,l))*Z0xoRo2(j,l)-Z1xoRo(j,l);
                    M(4,2) = (a1-k2(j,l))*W0xoRo2(j,l)-W1xoRo(j,l);
                    M(4,3) = k(j,l)*(Z1yoRo(j,l)-y2(j,l)*Z0yoRo2(j,l));
                    M(4,4) = k(j,l)*(W1yoRo(j,l)-y2(j,l)*W0yoRo2(j,l));
                    M(5,1) = -2*k(j,l)*Z1xoRo(j,l);
                    M(5,2) = -2*k(j,l)*W1xoRo(j,l);
                    M(5,3) = (k2(j,l)-y2(j,l))*Z1yoRo(j,l);
                    M(5,4) = (k2(j,l)-y2(j,l))*W1yoRo(j,l);
                elseif ~ToggleInnerFluid && ToggleOuterFluid
                    M(1,1) = (a1-k2(j,l))*Z0xiRi2(j,l)-Z1xiRi(j,l);
                    M(1,2) = (a1-k2(j,l))*W0xiRi2(j,l)-W1xiRi(j,l);
                    M(1,3) = k(j,l)*(Z1yiRi(j,l)-y2(j,l)*Z0yiRi2(j,l));
                    M(1,4) = k(j,l)*(W1yiRi(j,l)-y2(j,l)*W0yiRi2(j,l));
                    M(2,1) = -2*k(j,l)*Z1xiRi(j,l);
                    M(2,2) = -2*k(j,l)*W1xiRi(j,l);
                    M(2,3) = (k2(j,l)-y2(j,l))*Z1yiRi(j,l);
                    M(2,4) = (k2(j,l)-y2(j,l))*W1yiRi(j,l);
                    M(3,1) = (a1-k2(j,l))*Z0xoRo2(j,l)-Z1xoRo(j,l);
                    M(3,2) = (a1-k2(j,l))*W0xoRo2(j,l)-W1xoRo(j,l);
                    M(3,3) = k(j,l)*(Z1yoRo(j,l)-y2(j,l)*Z0yoRo2(j,l));
                    M(3,4) = k(j,l)*(W1yoRo(j,l)-y2(j,l)*W0yoRo2(j,l));
                    M(3,5) = a3(j,l);
                    M(4,1) = -2*k(j,l)*Z1xoRo(j,l);
                    M(4,2) = -2*k(j,l)*W1xoRo(j,l);
                    M(4,3) = (k2(j,l)-y2(j,l))*Z1yoRo(j,l);
                    M(4,4) = (k2(j,l)-y2(j,l))*W1yoRo(j,l);
                    M(5,1) = -Z1xoRo(j,l);
                    M(5,2) = -W1xoRo(j,l);
                    M(5,3) = k(j,l)*Z1yoRo(j,l);
                    M(5,4) = k(j,l)*W1yoRo(j,l);
                    M(5,5) = H1zoRo(j,l);
                elseif ~ToggleInnerFluid && ~ToggleOuterFluid
                    M(1,1) = (a1-k2(j,l))*Z0xiRi2(j,l)-Z1xiRi(j,l);
                    M(1,2) = (a1-k2(j,l))*W0xiRi2(j,l)-W1xiRi(j,l);
                    M(1,3) = k(j,l)*(Z1yiRi(j,l)-y2(j,l)*Z0yiRi2(j,l));
                    M(1,4) = k(j,l)*(W1yiRi(j,l)-y2(j,l)*W0yiRi2(j,l));
                    M(2,1) = -2*k(j,l)*Z1xiRi(j,l);
                    M(2,2) = -2*k(j,l)*W1xiRi(j,l);
                    M(2,3) = (k2(j,l)-y2(j,l))*Z1yiRi(j,l);
                    M(2,4) = (k2(j,l)-y2(j,l))*W1yiRi(j,l);
                    M(3,1) = (a1-k2(j,l))*Z0xoRo2(j,l)-Z1xoRo(j,l);
                    M(3,2) = (a1-k2(j,l))*W0xoRo2(j,l)-W1xoRo(j,l);
                    M(3,3) = k(j,l)*(Z1yoRo(j,l)-y2(j,l)*Z0yoRo2(j,l));
                    M(3,4) = k(j,l)*(W1yoRo(j,l)-y2(j,l)*W0yoRo2(j,l));
                    M(4,1) = -2*k(j,l)*Z1xoRo(j,l);
                    M(4,2) = -2*k(j,l)*W1xoRo(j,l);
                    M(4,3) = (k2(j,l)-y2(j,l))*Z1yoRo(j,l);
                    M(4,4) = (k2(j,l)-y2(j,l))*W1yoRo(j,l);
                end
                Yc(j,l) = det(M);
            end
        end
    end
    Y = abs(Yc);
end
function Y = Computer_Coupled(ModeFamily,MatrixMethod,k,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kFu2,kFl2,gFu,gFl,SignUpperFluid,SignLowerFluid,c,rw2,a11,a21,a31,b12,b22,b23,b32,b33,b34,YSize)
    k2 = k.^2;
    if  ToggleUpperFluid
        k3Fu = SignUpperFluid.*sqrt(kFu2-k2);
        k3Fu = reshape(k3Fu,1,1,[]);
    end
    if  ToggleLowerFluid
        k3Fl = SignLowerFluid.*sqrt(kFl2-k2);
        k3Fl = reshape(k3Fl,1,1,[]);
    end
    k = reshape(k,1,1,[]);
    k2 = reshape(k2,1,1,[]);
    k4 = k2.^2;
    k6 = k2.^3;
    Length = length(k);
    Yc = NaN(YSize)+1i*NaN;
    for m = 1:SuperLayerSize
        A1 = a11(m)*k2+b12(m);
        A2 = a21(m)*k4+b22(m)*k2+b23(m);
        A3 = a31(m)*k6+b32(m)*k4+b33(m)*k2+b34(m);
        d1 = A1/3;
        d2 = A2/3-d1.^2;
        d3 = d1.^3-d1.*A2/2+A3/2;
        d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
        d5 = d2./d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        k32 = [d6+d7 d6-d7 d4-d5-d1];
        k3 = sqrt(k32);
        k3k = k3.*k;
        m11 = c{m}(1,1)*k2+c{m}(5,5)*k32-rw2(m);
        m22 = c{m}(6,6)*k2+c{m}(4,4)*k32-rw2(m);
        m12 = c{m}(1,6)*k2+c{m}(4,5)*k32;
        m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
        m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
        m1 = m13.*m22-m12.*m23;
        V = (m11.*m23-m13.*m12)./m1;
        W = (m11.*m22-m12.^2)./-m1;
        e1 = k.*W+k3;
        e2 = k3.*V;
        D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V).*k+c{m}(3,3)*k3.*W);
        D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
        D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
        if  Symmetric && SuperLayerSize == 1
            Phi = .5i*k3*LayerThicknesses;
        else
            Phi = 1i*k3*LayerThicknesses(m);
        end
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
        if  Symmetric
            if  ModeFamily == 1
                if  FluidLoading
                    M{end} = [M{end}(3:6,[1 2 4],:) [k3Fu./k;gFu./k;zeros(2,1,Length)]];
                    for j = 1:Length
                        Yc(j) = det(M{end}(:,:,j));
                    end
                else
                    for j = 1:Length
                        Yc(j) = det(M{end}(4:6,[1 2 4],j));
                    end
                end
            elseif ModeFamily == 2
                if  FluidLoading
                    M{end} = [M{end}(3:6,[3 5 6],:) [k3Fu./k;gFu./k;zeros(2,1,Length)]];
                    for j = 1:Length
                        Yc(j) = det(M{end}(:,:,j));
                    end
                else
                    for j = 1:Length
                        Yc(j) = det(M{end}(4:6,[3 5 6],j));
                    end
                end
            end
        else
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
                    QFu = gFu./k3Fu;
                    QFl = gFl./k3Fl;
                    for j = 1:Length
                        M11(1,1,j) = det(M{end}([3 5 6],1:3,j));
                        M12(1,1,j) = det(M{end}([3 5 6],[1 2 4],j));
                        M21(1,1,j) = det(M{end}(4:6,1:3,j));
                        M22(1,1,j) = det(M{end}(4:6,[1 2 4],j));
                    end
                    Yc = M21-QFl.*M22-QFu.*(M11-QFl.*M12);
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    for j = 1:Length
                        M11(1,1,j) = det(M{end}([3 5 6],1:3,j));
                        M21(1,1,j) = det(M{end}(4:6,1:3,j));
                    end
                    Yc = M21-gFu./k3Fu.*M11;
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    for j = 1:Length
                        M21(1,1,j) = det(M{end}(4:6,1:3,j));
                        M22(1,1,j) = det(M{end}(4:6,[1 2 4],j));
                    end
                    Yc = M21-gFl./k3Fl.*M22;
                end
            else
                for j = 1:Length
                    Yc(j) = det(M{end}(4:6,1:3,j));
                end
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
        if  Symmetric
            if  ModeFamily == 1
                if  FluidLoading
                    M{end} = pageinv(M{end});
                    Yc = k3Fu./gFu+M{end}(3,1,:)-M{end}(3,4,:).*M{end}(6,1,:)./M{end}(6,4,:);
                else
                    for j = 1:Length
                        Yc(j) = det(M{end}([1:3 5 6],1:5,j));
                    end
                end
            elseif ModeFamily == 2
                if  FluidLoading
                    M{end} = pageinv(M{end});
                    Yc = k3Fu./gFu+M{end}(3,1,:)-M{end}(3,5,:).*M{end}(5,1,:)./M{end}(5,5,:)-(M{end}(3,5,:).*M{end}(5,6,:)./M{end}(5,5,:)-M{end}(3,6,:)).*(M{end}(4,1,:).*M{end}(5,5,:)-M{end}(4,5,:).*M{end}(5,1,:))./(M{end}(4,5,:).*M{end}(5,6,:)-M{end}(4,6,:).*M{end}(5,5,:));
                else
                    for j = 1:Length
                        Yc(j) = det(M{end}(1:4,[1:3 6],j));
                    end
                end
            end
        else
            if  SymmetricSystem
                M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
                M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
                M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
                M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
            end
            if  FluidLoading
                M{end} = pageinv(M{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    WFu = k3Fu./k;
                    WFl = k3Fl./k;
                    DFu = gFu./k;
                    DFl = gFl./k;
                    Yc = (WFl-M{end}(6,4,:).*DFl).*(WFu+M{end}(3,1,:).*DFu)+M{end}(3,4,:).*M{end}(6,1,:).*DFu.*DFl; 
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    Yc = k3Fu./gFu+M{end}(3,1,:);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    Yc = k3Fl./gFl-M{end}(6,4,:);
                end
            else
                for j = 1:Length
                    Yc(j) = det(M{end}(:,:,j));
                end
            end
        end
    end
    Yc = reshape(Yc,YSize);
    Y = abs(Yc);
end
function Y = Computer_Decoupled(ModeFamily,MatrixMethod,k,LayerThicknesses,SuperLayerSize,Pattern,Symmetric,SymmetricSystem,I1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,kFu2,kFl2,gFu,gFl,SignUpperFluid,SignLowerFluid,c,rw2,A1,a21,a31,b22,b32,b33,YSize)
    k2 = k.^2;
    if  ToggleUpperFluid
        k3Fu = SignUpperFluid.*sqrt(kFu2-k2);
        k3Fu = reshape(k3Fu,1,1,[]);
    end
    if  ToggleLowerFluid
        k3Fl = SignLowerFluid.*sqrt(kFl2-k2);
        k3Fl = reshape(k3Fl,1,1,[]);
    end
    k = reshape(k,1,1,[]);
    k2 = reshape(k2,1,1,[]);
    k4 = k2.^2;
    Length = length(k);
    Yc = NaN(YSize)+1i*NaN;
    for m = 1:SuperLayerSize
        A2 = a21(m)*k2+b22(m);
        A3 = a31(m)*k4+b32(m)*k2+b33(m);
        d1 = sqrt(A2.^2-2*A1(m)*A3);
        k32 = [d1-A2 -d1-A2]/A1(m);
        k3 = sqrt(k32);
        W = (rw2(m)-c{m}(1,1)*k2-c{m}(5,5)*k32)./((c{m}(1,3)+c{m}(5,5))*k.*k3);
        D3 = 1i*(c{m}(1,3)*k+c{m}(3,3)*k3.*W);
        D5 = 1i*(c{m}(5,5)*(k3+k.*W));
        if  Symmetric && SuperLayerSize == 1
            Phi = .5i*k3*LayerThicknesses;
        else
            Phi = 1i*k3*LayerThicknesses(m);
        end
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
        if  Symmetric
            if  ModeFamily == 1
                if  FluidLoading
                    M{end} = [M{end}(2:4,[1 3],:) [k3Fu./k;gFu./k;zeros(1,1,Length)]];
                    for j = 1:Length
                        Yc(j) = det(M{end}(:,:,j));
                    end
                else
                    for j = 1:Length
                        Yc(j) = det(M{end}(3:4,[1 3],j));
                    end
                end
            elseif ModeFamily == 2
                if  FluidLoading
                    M{end} = [M{end}(2:4,[2 4],:) [k3Fu./k;gFu./k;zeros(1,1,Length)]];
                    for j = 1:Length
                        Yc(j) = det(M{end}(:,:,j));
                    end
                else
                    for j = 1:Length
                        Yc(j) = det(M{end}(3:4,[2 4],j));
                    end
                end
            end
        else
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
                    QFu = gFu./k3Fu;
                    QFl = gFl./k3Fl;
                    for j = 1:Length
                        M11(1,1,j) = det(M{end}([2 4],1:2,j));
                        M12(1,1,j) = det(M{end}([2 4],[1 3],j));
                        M21(1,1,j) = det(M{end}(3:4,1:2,j));
                        M22(1,1,j) = det(M{end}(3:4,[1 3],j));
                    end
                    Yc = M21-QFl.*M22-QFu.*(M11-QFl.*M12);
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    for j = 1:Length
                        M11(1,1,j) = det(M{end}([2 4],1:2,j));
                        M21(1,1,j) = det(M{end}(3:4,1:2,j));
                    end
                    Yc = M21-gFu./k3Fu.*M11;
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    for j = 1:Length
                        M21(1,1,j) = det(M{end}(3:4,1:2,j));
                        M22(1,1,j) = det(M{end}(3:4,[1 3],j));
                    end
                    Yc = M21-gFl./k3Fl.*M22;
                end
            else
                for j = 1:Length
                    Yc(j) = det(M{end}(3:4,1:2,j));
                end
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
        if  Symmetric
            if  ModeFamily == 1
                if  FluidLoading
                    M{end} = pageinv(M{end});
                    Yc = k3Fu./gFu+M{end}(2,1,:)-M{end}(2,3,:).*M{end}(4,1,:)./M{end}(4,3,:);
                else
                    for j = 1:Length
                        Yc(j) = det(M{end}([1 2 4],1:3,j));
                    end
                end
            elseif ModeFamily == 2
                if  FluidLoading
                    M{end} = pageinv(M{end});
                    Yc = k3Fu./gFu+M{end}(2,1,:)-M{end}(2,4,:).*M{end}(3,1,:)./M{end}(3,4,:);
                else
                    for j = 1:Length
                        Yc(j) = det(M{end}(1:3,[1 2 4],j));
                    end
                end
            end
        else
            if  SymmetricSystem
                M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
                M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
                M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
                M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
            end
            if  FluidLoading
                M{end} = pageinv(M{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    WFu = k3Fu./k;
                    WFl = k3Fl./k;
                    DFu = gFu./k;
                    DFl = gFl./k;
                    Yc = (WFl-M{end}(4,3,:).*DFl).*(WFu+M{end}(2,1,:).*DFu)+M{end}(2,3,:).*M{end}(4,1,:).*DFu.*DFl;
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    Yc = k3Fu./gFu+M{end}(2,1,:);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    Yc = k3Fl./gFl-M{end}(4,3,:);
                end
            else
                for j = 1:Length
                    Yc(j) = det(M{end}(:,:,j));
                end
            end
        end
    end
    Yc = reshape(Yc,YSize);
    Y = abs(Yc);
end