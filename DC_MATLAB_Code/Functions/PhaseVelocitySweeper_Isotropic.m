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
function [FLambF,FScholte] = PhaseVelocitySweeper_Isotropic(Material,Half,Symmetric,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Frequency)
Resolution = 1e-6; % (m/s)
DensityThreshold = .01; % rho_fluid/rho_solid

%#ok<*AGROW>
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
if  Symmetric
    Fluid = UpperFluid;
    kF2 = (AngularFrequency/Fluid.Velocity)^2;
    Density_ = Fluid.Density/Material.Density;
    kUpperFluid2 = NaN;
    kLowerFluid2 = NaN;
else
    if  ToggleUpperFluid && ToggleLowerFluid
        kUpperFluid2 = (AngularFrequency/UpperFluid.Velocity)^2;
        kLowerFluid2 = (AngularFrequency/LowerFluid.Velocity)^2;
        Density_ = .5*(UpperFluid.Density+LowerFluid.Density)/Material.Density;
        Fluid = NaN;
        kF2 = NaN;
    else
        if  ToggleUpperFluid
            Fluid = UpperFluid;
        elseif ToggleLowerFluid
            Fluid = LowerFluid;
        end
        kF2 = (AngularFrequency/Fluid.Velocity)^2;
        Density_ = .5*Fluid.Density/Material.Density;
        kUpperFluid2 = NaN;
        kLowerFluid2 = NaN;
    end
end

Range1 = .02; % for the weakly damped solution near the real axis
Steps1 = 2e2;
if  Density_ > DensityThreshold % for the strongly damped solutions when fluid-loading is present
    Range2 = 1.8;
else
    Range2 = .1;
end
Steps2 = 2e2;

% rough search for the undamped A0/Scholte mode on the real axis
SweepRange = 1:1000;
Wavenumber = AngularFrequency./SweepRange;
k2 = Wavenumber.^2;
x = sqrt(kL2-k2);
y = sqrt(kT2-k2);
xH = x*Half;
yH = y*Half;
if  Symmetric
    R = (y.^2-k2).^2./y.*tan(xH)+4*k2.*x.*tan(yH);
    if  FluidLoading
        a = Fluid.Density*kT2^2*x./(y*Material.Density.*sqrt(kF2-k2));
        Y = abs(R+1i*a);
    else
        Y = abs(R);
    end
    for i = 2:length(Y)-1
        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
            XRough = SweepRange(i);
            break
        end
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
    a1 = -(Lambda*kL2+2*Mu*x.^2);
    a2 = Mu*(k2-y.^2);
    a3 = 2i*Mu*Wavenumber.*x;
    a4 = 2i*Mu*Wavenumber.*y;
    a5 = -1i*Wavenumber;
    if  ToggleUpperFluid && ToggleLowerFluid
        k3UpperFluid = sqrt(kUpperFluid2-k2);
        k3LowerFluid = sqrt(kLowerFluid2-k2);
        e = exp(1i*k3UpperFluid*Half);
        e_ = exp(-1i*k3LowerFluid*Half); 
        a6 = -UpperFluid.Velocity^2*UpperFluid.Density*kUpperFluid2;
        a7 = LowerFluid.Velocity^2*LowerFluid.Density*kLowerFluid2;
        a8 = 1i*k3UpperFluid;
        a9 = 1i*k3LowerFluid;
        for i = 1:length(SweepRange)
            M(1,1) = a1(i)*SinxH(i);
            M(1,2) = a1(i)*CosxH(i);
            M(1,3) = -a4(i)*CosyH(i);
            M(1,4) = a4(i)*SinyH(i);
            M(1,5) = a6*e(i);
            M(2,1) = a3(i)*CosxH(i);
            M(2,2) = -a3(i)*SinxH(i);
            M(2,3) = a2(i)*SinyH(i);
            M(2,4) = a2(i)*CosyH(i);
            M(3,1) = x(i)*CosxH(i);
            M(3,2) = -x(i)*SinxH(i);
            M(3,3) = a5(i)*SinyH(i);
            M(3,4) = a5(i)*CosyH(i);
            M(3,5) = a8(i)*e(i);
            M(4,1) = a1(i)*Sin_xH(i);
            M(4,2) = a1(i)*Cos_xH(i);
            M(4,3) = -a4(i)*Cos_yH(i);
            M(4,4) = a4(i)*Sin_yH(i);
            M(4,6) = a7*e_(i);
            M(5,1) = a3(i)*Cos_xH(i);
            M(5,2) = -a3(i)*Sin_xH(i);
            M(5,3) = a2(i)*Sin_yH(i);
            M(5,4) = a2(i)*Cos_yH(i);
            M(6,1) = x(i)*Cos_xH(i);
            M(6,2) = -x(i)*Sin_xH(i);
            M(6,3) = a5(i)*Sin_yH(i);
            M(6,4) = a5(i)*Cos_yH(i);
            M(6,6) = a9(i)*e_(i);
            Y(i) = abs(det(M));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                XRough = SweepRange(i-1);
                break
            end
        end
    else
        k3Fluid = sqrt(kF2-k2);
        e = exp(1i*k3Fluid*Half);
        a6 = -Fluid.Velocity^2*Fluid.Density*kF2;
        a8 = 1i*k3Fluid;
        for i = 1:length(SweepRange)
            M(1,1) = a1(i)*SinxH(i);
            M(1,2) = a1(i)*CosxH(i);
            M(1,3) = -a4(i)*CosyH(i);
            M(1,4) = a4(i)*SinyH(i);
            M(1,5) = a6*e(i);
            M(2,1) = a3(i)*CosxH(i);
            M(2,2) = -a3(i)*SinxH(i);
            M(2,3) = a2(i)*SinyH(i);
            M(2,4) = a2(i)*CosyH(i);
            M(3,1) = x(i)*CosxH(i);
            M(3,2) = -x(i)*SinxH(i);
            M(3,3) = a5(i)*SinyH(i);
            M(3,4) = a5(i)*CosyH(i);
            M(3,5) = a8(i)*e(i);
            M(4,1) = a1(i)*Sin_xH(i);
            M(4,2) = a1(i)*Cos_xH(i);
            M(4,3) = -a4(i)*Cos_yH(i);
            M(4,4) = a4(i)*Sin_yH(i);
            M(5,1) = a3(i)*Cos_xH(i);
            M(5,2) = -a3(i)*Sin_xH(i);
            M(5,3) = a2(i)*Sin_yH(i);
            M(5,4) = a2(i)*Cos_yH(i);
            Y(i) = abs(det(M));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                XRough = SweepRange(i-1);
                break
            end
        end
    end
end
% figure,plot(20*log10(Y))

% fine search for the undamped A0/Scholte mode on the real axis 
Bisections = ceil(log2(Resolution/(SweepRange(2)-SweepRange(1)))/log2(2*.25));
PhaseVelocity = [XRough-(SweepRange(2)-SweepRange(1)) XRough+(SweepRange(2)-SweepRange(1))];
for o = 1:Bisections
    PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
    for i = 1:length(PhaseVelocity)
        Wavenumber = AngularFrequency/PhaseVelocity(i);
        k2 = Wavenumber^2;
        x = sqrt(kL2-k2);
        y = sqrt(kT2-k2);
        xH = x*Half;
        yH = y*Half;
        if  Symmetric
            R = (y^2-k2)^2/y*tan(xH)+4*k2*x*tan(yH);
            if  FluidLoading
                a = Fluid.Density*kT2^2*x/(y*Material.Density*sqrt(kF2-k2));
                Y(i) = abs(R+1i*a);
            else
                Y(i) = abs(R);
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
            a1 = -(Lambda*kL2+2*Mu*x^2);
            a2 = Mu*(k2-y^2);
            a3 = 2i*Mu*Wavenumber*x;
            a4 = 2i*Mu*Wavenumber*y;
            a5 = -1i*Wavenumber;
            if  ToggleUpperFluid && ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-k2);
                k3LowerFluid = sqrt(kLowerFluid2-k2);
                e = exp(1i*k3UpperFluid*Half);
                e_ = exp(-1i*k3LowerFluid*Half); 
                a6 = -UpperFluid.Velocity^2*UpperFluid.Density*kUpperFluid2;
                a7 = LowerFluid.Velocity^2*LowerFluid.Density*kLowerFluid2;
                a8 = 1i*k3UpperFluid;
                a9 = 1i*k3LowerFluid;
                M(1,1) = a1*SinxH;
                M(1,2) = a1*CosxH;
                M(1,3) = -a4*CosyH;
                M(1,4) = a4*SinyH;
                M(1,5) = a6*e;
                M(2,1) = a3*CosxH;
                M(2,2) = -a3*SinxH;
                M(2,3) = a2*SinyH;
                M(2,4) = a2*CosyH;
                M(3,1) = x*CosxH;
                M(3,2) = -x*SinxH;
                M(3,3) = a5*SinyH;
                M(3,4) = a5*CosyH;
                M(3,5) = a8*e;
                M(4,1) = a1*Sin_xH;
                M(4,2) = a1*Cos_xH;
                M(4,3) = -a4*Cos_yH;
                M(4,4) = a4*Sin_yH;
                M(4,6) = a7*e_;
                M(5,1) = a3*Cos_xH;
                M(5,2) = -a3*Sin_xH;
                M(5,3) = a2*Sin_yH;
                M(5,4) = a2*Cos_yH;
                M(6,1) = x*Cos_xH;
                M(6,2) = -x*Sin_xH;
                M(6,3) = a5*Sin_yH;
                M(6,4) = a5*Cos_yH;
                M(6,6) = a9*e_;
                Y(i) = abs(det(M));
            else
                k3Fluid = sqrt(kF2-k2);
                e = exp(1i*k3Fluid*Half);
                a6 = -Fluid.Velocity^2*Fluid.Density*kF2;
                a8 = 1i*k3Fluid;
                M(1,1) = a1*SinxH;
                M(1,2) = a1*CosxH;
                M(1,3) = -a4*CosyH;
                M(1,4) = a4*SinyH;
                M(1,5) = a6*e;
                M(2,1) = a3*CosxH;
                M(2,2) = -a3*SinxH;
                M(2,3) = a2*SinyH;
                M(2,4) = a2*CosyH;
                M(3,1) = x*CosxH;
                M(3,2) = -x*SinxH;
                M(3,3) = a5*SinyH;
                M(3,4) = a5*CosyH;
                M(3,5) = a8*e;
                M(4,1) = a1*Sin_xH;
                M(4,2) = a1*Cos_xH;
                M(4,3) = -a4*Cos_yH;
                M(4,4) = a4*Sin_yH;
                M(5,1) = a3*Cos_xH;
                M(5,2) = -a3*Sin_xH;
                M(5,3) = a2*Sin_yH;
                M(5,4) = a2*Cos_yH;
                Y(i) = abs(det(M));
            end
        end
        if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
            if  o == Bisections
                FLamb = PhaseVelocity(i-1);
            end
            PhaseVelocity = [PhaseVelocity(i-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(i-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
            break
        end
    end
end

% search for the weakly damped Scholte mode near the real axis
if  Viscoelastic

    % rough search for the weakly damped Scholte mode near the real axis
    XRough = [];
    SweepRangeReal = (1-.5*Range1)*FLamb:Range1*FLamb/Steps1:(1+.5*Range1)*FLamb;
    SweepRangeImag = -Range1*FLamb:Range1*FLamb/Steps1:0;
    Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
    Y = Computer(-1,Wavenumber,kL2,kT2,Half,Material,Lambda,Mu,Symmetric,FluidLoading,ToggleLowerFluid,ToggleUpperFluid,Fluid,LowerFluid,UpperFluid,kF2,kLowerFluid2,kUpperFluid2);
    if  ~FluidLoading
        Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
    end
% if  ~FluidLoading
% figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
% else
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
% end
    for l = 2:size(Y,2)-1
        for j = 2:size(Y,1)-1
            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                XRough = [SweepRangeReal(j) SweepRangeImag(l)]; % real velocity (m/s), imaginary velocity (m/s)
            end
        end
    end
    
    % fine search for the weakly damped Scholte mode near the real axis
    if  ~isempty(XRough)
        Bisections = ceil(log2(Resolution/(2*(SweepRangeReal(2)-SweepRangeReal(1))))/log2(2*.25));
        RangeReal = [XRough(1)-2*(SweepRangeReal(2)-SweepRangeReal(1)) XRough(1)+2*(SweepRangeReal(2)-SweepRangeReal(1))];
        RangeImag = [XRough(2)-2*(SweepRangeImag(2)-SweepRangeImag(1)) XRough(2)+2*(SweepRangeImag(2)-SweepRangeImag(1))];
        for k = 1:Bisections % search minimum in characteristic equation and converge upon it
            if  k == 1
                RangeReal = RangeReal(1):.125*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                RangeImag = RangeImag(1):.125*(RangeImag(end)-RangeImag(1)):RangeImag(end);
            else
                RangeReal = RangeReal(1):.25*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                RangeImag = RangeImag(1):.25*(RangeImag(end)-RangeImag(1)):RangeImag(end);
            end
            Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
            Y = Computer(-1,Wavenumber,kL2,kT2,Half,Material,Lambda,Mu,Symmetric,FluidLoading,ToggleLowerFluid,ToggleUpperFluid,Fluid,LowerFluid,UpperFluid,kF2,kLowerFluid2,kUpperFluid2);
            if  abs(RangeImag(end)) < 1e-3
                Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
            end
% if  abs(RangeImag(end)) < 1e-3
% f = figure;surf(horzcat(RangeImag,-RangeImag(end-1)),RangeReal,20*log10(Y))    
% else
% f = figure;surf(RangeImag,RangeReal,20*log10(Y))
% end
% close(f)
            for l = 2:size(Y,2)-1
                for j = 2:size(Y,1)-1
                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                        MIN = [j l];
                    end
                end
            end
            if  k < Bisections % set the new search area around the found minimum
                RangeReal = [RangeReal(MIN(1))-(RangeReal(2)-RangeReal(1)) RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))];
                RangeImag = [RangeImag(MIN(2))-(RangeImag(2)-RangeImag(1)) RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))];
            end
        end
        XFine(1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
        XFine(2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
        XFine(3) = RangeReal(MIN(1)); % real velocity (m/s)
        XFine(4) = RangeImag(MIN(2)); % imaginary velocity (m/s)
    else
        XFine = [FLamb 0 FLamb 0];
    end
else
    XFine = [FLamb 0 FLamb 0];
end

% search for the strongly damped A0 when fluid-loading is present
if  FluidLoading

    % rough search for the strongly damped A0 when fluid-loading is present
    XRough = [];
    for o = [1 2 5]
        SweepRangeReal = (1-.5*Range2)*FLamb:Range2*FLamb/Steps2/o:(1+.5*Range2)*FLamb;
        SweepRangeImag = -Range2*FLamb:Range2*FLamb/Steps2/o:0;
        Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
        Y = Computer(1,Wavenumber,kL2,kT2,Half,Material,Lambda,Mu,Symmetric,FluidLoading,ToggleLowerFluid,ToggleUpperFluid,Fluid,LowerFluid,UpperFluid,kF2,kLowerFluid2,kUpperFluid2);
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
        for l = 2:size(Y,2)-1
            for j = 2:size(Y,1)-1
                if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                    XRough(end+1,:) = [SweepRangeReal(j) SweepRangeImag(l) imag(Wavenumber(j,l))]; % real velocity (m/s), imaginary velocity (m/s)
                end
            end
        end
        if  ~isempty(XRough)
            break
        end
    end

    % fine search for the strongly damped A0 when fluid-loading is present
    if  ~isempty(XRough)
        Bisections = ceil(log2(Resolution/(2*(SweepRangeReal(2)-SweepRangeReal(1))))/log2(2*.25));
        RangeReal = [XRough(1)-2*(SweepRangeReal(2)-SweepRangeReal(1)) XRough(1)+2*(SweepRangeReal(2)-SweepRangeReal(1))];
        RangeImag = [XRough(2)-2*(SweepRangeImag(2)-SweepRangeImag(1)) XRough(2)+2*(SweepRangeImag(2)-SweepRangeImag(1))];
        for k = 1:Bisections % search minimum in characteristic equation and converge upon it
            if  k == 1
                RangeReal = RangeReal(1):.125*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                RangeImag = RangeImag(1):.125*(RangeImag(end)-RangeImag(1)):RangeImag(end);
            else
                RangeReal = RangeReal(1):.25*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                RangeImag = RangeImag(1):.25*(RangeImag(end)-RangeImag(1)):RangeImag(end);
            end
            Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
            Y = Computer(1,Wavenumber,kL2,kT2,Half,Material,Lambda,Mu,Symmetric,FluidLoading,ToggleLowerFluid,ToggleUpperFluid,Fluid,LowerFluid,UpperFluid,kF2,kLowerFluid2,kUpperFluid2);
            if  abs(RangeImag(end)) < 1e-3
                Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
            end
% if  abs(RangeImag(end)) < 1e-3
% f = figure;surf(horzcat(RangeImag,-RangeImag(end-1)),RangeReal,20*log10(Y))
% else
% f = figure;surf(RangeImag,RangeReal,20*log10(Y))
% end
% close(f)
            for l = 2:size(Y,2)-1
                for j = 2:size(Y,1)-1
                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                        MIN = [j l];
                    end
                end
            end
            if  k < Bisections % set the new search area around the found minimum
                RangeReal = [RangeReal(MIN(1))-(RangeReal(2)-RangeReal(1)) RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))];
                RangeImag = [RangeImag(MIN(2))-(RangeImag(2)-RangeImag(1)) RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))];
            end
        end
        FLambF(1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
        FLambF(2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
        FLambF(3) = RangeReal(MIN(1)); % real phase velocity (m/s)
        FLambF(4) = RangeImag(MIN(2)); % imaginary phase velocity (m/s)
    else
        FLambF = XFine; 
    end
    FScholte = XFine;
else
    FLambF = XFine;
    FScholte = [0 0 0 0];
end
% disp(['Undamped  ',num2str(FLamb),newline...
%       'A0Scholte ',num2str(FScholte(1)),', ',num2str(FScholte(2)),', ',num2str(FScholte(3)),', ',num2str(FScholte(4)),newline...;    
%       'A0        ',num2str(FLambF(1)),', ',num2str(FLambF(2)),', ',num2str(FLambF(3)),', ',num2str(FLambF(4)),newline,'-----------------------'])
end
function Y = Computer(Sign,Wavenumber,kL2,kT2,Half,Material,Lambda,Mu,Symmetric,FluidLoading,ToggleLowerFluid,ToggleUpperFluid,Fluid,LowerFluid,UpperFluid,kF2,kLowerFluid2,kUpperFluid2)
    k2 = Wavenumber.^2;
    x = sqrt(kL2-k2);
    y = sqrt(kT2-k2);
    xH = x*Half;
    yH = y*Half;
    if  Symmetric
        R = (y.^2-k2).^2./y.*tan(xH)+4*k2.*x.*tan(yH);
        if  FluidLoading
            a = Sign*Fluid.Density*kT2^2*x./(y*Material.Density.*sqrt(kF2-k2));
            Y = abs(R+1i*a);
        else
            Y = abs(R);
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
        a1 = -(Lambda*kL2+2*Mu*x.^2);
        a2 = Mu*(k2-y.^2);
        a3 = 2i*Mu*Wavenumber.*x;
        a4 = 2i*Mu*Wavenumber.*y;
        a5 = -1i*Wavenumber;
        Y = NaN(size(Wavenumber));
        if  ToggleUpperFluid && ToggleLowerFluid
            k3UpperFluid = Sign*sqrt(kUpperFluid2-k2);
            k3LowerFluid = Sign*sqrt(kLowerFluid2-k2);
            e = exp(1i*k3UpperFluid*Half);
            e_ = exp(-1i*k3LowerFluid*Half); 
            a6 = -UpperFluid.Velocity^2*UpperFluid.Density*kUpperFluid2;
            a7 = LowerFluid.Velocity^2*LowerFluid.Density*kLowerFluid2;
            a8 = 1i*k3UpperFluid;
            a9 = 1i*k3LowerFluid;
            for l = 1:width(Wavenumber)
                for j = 1:height(Wavenumber)
                    M(1,1) = a1(j,l)*SinxH(j,l);
                    M(1,2) = a1(j,l)*CosxH(j,l);
                    M(1,3) = -a4(j,l)*CosyH(j,l);
                    M(1,4) = a4(j,l)*SinyH(j,l);
                    M(1,5) = a6*e(j,l);
                    M(2,1) = a3(j,l)*CosxH(j,l);
                    M(2,2) = -a3(j,l)*SinxH(j,l);
                    M(2,3) = a2(j,l)*SinyH(j,l);
                    M(2,4) = a2(j,l)*CosyH(j,l);
                    M(3,1) = x(j,l)*CosxH(j,l);
                    M(3,2) = -x(j,l)*SinxH(j,l);
                    M(3,3) = a5(j,l)*SinyH(j,l);
                    M(3,4) = a5(j,l)*CosyH(j,l);
                    M(3,5) = a8(j,l)*e(j,l);
                    M(4,1) = a1(j,l)*Sin_xH(j,l);
                    M(4,2) = a1(j,l)*Cos_xH(j,l);
                    M(4,3) = -a4(j,l)*Cos_yH(j,l);
                    M(4,4) = a4(j,l)*Sin_yH(j,l);
                    M(4,6) = a7*e_(j,l);
                    M(5,1) = a3(j,l)*Cos_xH(j,l);
                    M(5,2) = -a3(j,l)*Sin_xH(j,l);
                    M(5,3) = a2(j,l)*Sin_yH(j,l);
                    M(5,4) = a2(j,l)*Cos_yH(j,l);
                    M(6,1) = x(j,l)*Cos_xH(j,l);
                    M(6,2) = -x(j,l)*Sin_xH(j,l);
                    M(6,3) = a5(j,l)*Sin_yH(j,l);
                    M(6,4) = a5(j,l)*Cos_yH(j,l);
                    M(6,6) = a9(j,l)*e_(j,l);
                    Y(j,l) = abs(det(M));
                end
            end
        else
            k3Fluid = Sign*sqrt(kF2-k2);
            e = exp(1i*k3Fluid*Half);
            a6 = -Fluid.Velocity^2*Fluid.Density*kF2;
            a8 = 1i*k3Fluid;
            for l = 1:width(Wavenumber)
                for j = 1:height(Wavenumber)
                    M(1,1) = a1(j,l)*SinxH(j,l);
                    M(1,2) = a1(j,l)*CosxH(j,l);
                    M(1,3) = -a4(j,l)*CosyH(j,l);
                    M(1,4) = a4(j,l)*SinyH(j,l);
                    M(1,5) = a6*e(j,l);
                    M(2,1) = a3(j,l)*CosxH(j,l);
                    M(2,2) = -a3(j,l)*SinxH(j,l);
                    M(2,3) = a2(j,l)*SinyH(j,l);
                    M(2,4) = a2(j,l)*CosyH(j,l);
                    M(3,1) = x(j,l)*CosxH(j,l);
                    M(3,2) = -x(j,l)*SinxH(j,l);
                    M(3,3) = a5(j,l)*SinyH(j,l);
                    M(3,4) = a5(j,l)*CosyH(j,l);
                    M(3,5) = a8(j,l)*e(j,l);
                    M(4,1) = a1(j,l)*Sin_xH(j,l);
                    M(4,2) = a1(j,l)*Cos_xH(j,l);
                    M(4,3) = -a4(j,l)*Cos_yH(j,l);
                    M(4,4) = a4(j,l)*Sin_yH(j,l);
                    M(5,1) = a3(j,l)*Cos_xH(j,l);
                    M(5,2) = -a3(j,l)*Sin_xH(j,l);
                    M(5,3) = a2(j,l)*Sin_yH(j,l);
                    M(5,4) = a2(j,l)*Cos_yH(j,l);
                    Y(j,l) = abs(det(M));
                end
            end
        end
    end   
end