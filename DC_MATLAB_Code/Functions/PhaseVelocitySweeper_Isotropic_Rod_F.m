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
function [FLambF,FScholte] = PhaseVelocitySweeper_Isotropic_Rod_F(FLamb,Material,Viscoelastic,R,Fluid,FluidLoading,Frequency)
Resolution = 1e-6; % (m/s)
Range1 = .02; % for the weakly damped solution near the real axis
Steps1 = 2e2;
Range2 = .3; % for the strongly damped solution when fluid-loading is present
Steps2 = 2e2;

%#ok<*AGROW>
R2 = R^2;
Density = Fluid.Density/Material.Density;
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
kF2 = (AngularFrequency/Fluid.Velocity)^2;

% search for the weakly damped Scholte mode near the real axis
if  Viscoelastic

    % rough search for the weakly damped Scholte mode near the real axis
    XRough = [];
    SweepRangeReal = (1-.5*Range1)*FLamb:Range1*FLamb/Steps1:(1+.5*Range1)*FLamb;
    SweepRangeImag = -Range1*FLamb:Range1*FLamb/Steps1:0;
    Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
    Y = Computer(-1,Wavenumber,kL2,kT2,R,R2,Density,FluidLoading,kF2);
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
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
            Y = Computer(-1,Wavenumber,kL2,kT2,R,R2,Density,FluidLoading,kF2);
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

% search for the strongly damped F(1,1) when fluid-loading is present
if  FluidLoading

    % rough search for the strongly damped F(1,1) when fluid-loading is present
    XRough = [];
    for o = [1 2 5]
        SweepRangeReal = 0:Range2*FLamb/Steps2/o:Range2*FLamb;
        SweepRangeImag = -SweepRangeReal;
        Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
        Y = Computer(1,Wavenumber,kL2,kT2,R,R2,Density,FluidLoading,kF2);
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
    if  height(XRough) > 1
        [~,z] = min(abs(XRough(:,3)));
        XRough = XRough(z,:);
    end

    % fine search for the strongly damped F(1,1) when fluid-loading is present
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
            Y = Computer(1,Wavenumber,kL2,kT2,R,R2,Density,FluidLoading,kF2);
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
% disp(['Undamped      ',num2str(FLamb),newline...
%       'F(1,1)Scholte ',num2str(FScholte(1)),', ',num2str(FScholte(2)),', ',num2str(FScholte(3)),', ',num2str(FScholte(4)),newline...    
%       'F(1,1)        ',num2str(FLambF(1)),', ',num2str(FLambF(2)),', ',num2str(FLambF(3)),', ',num2str(FLambF(4)),newline,'-----------------------------------------']);
end
function Y = Computer(Sign,Wavenumber,kL2,kT2,R,R2,Density,FluidLoading,kF2)
    k2 = Wavenumber.^2;
    y2 = kT2-k2;
    xR = sqrt(kL2-k2)*R;
    yR = sqrt(kT2-k2)*R;
    J1x = besselj(1,xR);
    J1y = besselj(1,yR);
    dJ1xR = J1x-xR.*besselj(2,xR);
    dJ1yR = J1y-yR.*besselj(2,yR);
    if  FluidLoading
        zR = Sign*sqrt(kF2-k2)*R;
        H1z = besselh(1,zR);
        dH1zR = H1z-zR.*besselh(2,zR);
    end
    Y = NaN(size(Wavenumber));
    for l = 1:width(Wavenumber)
        for j = 1:height(Wavenumber)
            if  FluidLoading
                M(1,1) = dJ1xR(j,l)+(.5*kT2*R2-k2(j,l)*R2-1)*J1x(j,l);
                M(1,2) = -Wavenumber(j,l)*(dJ1yR(j,l)+(y2(j,l)*R2-1)*J1y(j,l));
                M(1,3) = dJ1yR(j,l)-J1y(j,l);
                M(1,4) = .5*kT2*R2*Density*H1z(j,l);
                M(2,1) = 2*(dJ1xR(j,l)-J1x(j,l));
                M(2,2) = 2*Wavenumber(j,l)*(J1y(j,l)-dJ1yR(j,l));
                M(2,3) = 2*dJ1yR(j,l)+(y2(j,l)*R2-2)*J1y(j,l);
                M(3,1) = 2*Wavenumber(j,l)*dJ1xR(j,l);
                M(3,2) = (y2(j,l)-k2(j,l))*dJ1yR(j,l);
                M(3,3) = -Wavenumber(j,l)*J1y(j,l);
                M(4,1) = dJ1xR(j,l);
                M(4,2) = -Wavenumber(j,l)*dJ1yR(j,l);
                M(4,3) = -J1y(j,l);
                M(4,4) = dH1zR(j,l);
            else
                M(1,1) = dJ1xR(j,l)+(.5*kT2*R2-k2(j,l)*R2-1)*J1x(j,l);
                M(1,2) = -Wavenumber(j,l)*(dJ1yR(j,l)+(y2(j,l)*R2-1)*J1y(j,l));
                M(1,3) = dJ1yR(j,l)-J1y(j,l);
                M(2,1) = 2*(dJ1xR(j,l)-J1x(j,l));
                M(2,2) = 2*Wavenumber(j,l)*(J1y(j,l)-dJ1yR(j,l));
                M(2,3) = 2*dJ1yR(j,l)+(y2(j,l)*R2-2)*J1y(j,l);
                M(3,1) = 2*Wavenumber(j,l)*dJ1xR(j,l);
                M(3,2) = (y2(j,l)-k2(j,l))*dJ1yR(j,l);
                M(3,3) = -Wavenumber(j,l)*J1y(j,l);
            end
            Y(j,l) = abs(det(M));
        end
    end
end