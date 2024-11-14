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
function FLamb = PhaseVelocitySweeper_Isotropic_Rod_F_UD(Material,R,Fluid,FluidLoading,Frequency)
%#ok<*AGROW>
Resolution = 1e-6; % (m/s)
R2 = R^2;
Density = Fluid.Density/Material.Density;
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
kF2 = (AngularFrequency/Fluid.Velocity)^2;

% rough search for the undamped F(1,1)/FScholte(1,1) mode on the real axis
SweepRange = 1-1e-10:1000;
Wavenumber = AngularFrequency./SweepRange;
k2 = Wavenumber.^2;
y2 = kT2-k2;
xR = sqrt(kL2-k2)*R;
yR = sqrt(kT2-k2)*R;
J1x = besselj(1,xR);
J1y = besselj(1,yR);
dJ1xR = J1x-xR.*besselj(2,xR);
dJ1yR = J1y-yR.*besselj(2,yR);
if  FluidLoading
    zR = sqrt(kF2-k2)*R;
    H1z = besselh(1,zR);
    dH1zR = H1z-zR.*besselh(2,zR);
end
for i = 1:length(SweepRange)
    if  FluidLoading
        M(1,1) = dJ1xR(i)+(.5*kT2*R2-k2(i)*R2-1)*J1x(i);
        M(1,2) = -Wavenumber(i)*(dJ1yR(i)+(y2(i)*R2-1)*J1y(i));
        M(1,3) = dJ1yR(i)-J1y(i);
        M(1,4) = .5*kT2*R2*Density*H1z(i);
        M(2,1) = 2*(dJ1xR(i)-J1x(i));
        M(2,2) = 2*Wavenumber(i)*(J1y(i)-dJ1yR(i));
        M(2,3) = 2*dJ1yR(i)+(y2(i)*R2-2)*J1y(i);
        M(3,1) = 2*Wavenumber(i)*dJ1xR(i);
        M(3,2) = (y2(i)-k2(i))*dJ1yR(i);
        M(3,3) = -Wavenumber(i)*J1y(i);
        M(4,1) = dJ1xR(i);
        M(4,2) = -Wavenumber(i)*dJ1yR(i);
        M(4,3) = -J1y(i);
        M(4,4) = dH1zR(i);
    else
        M(1,1) = dJ1xR(i)+(.5*kT2*R2-k2(i)*R2-1)*J1x(i);
        M(1,2) = -Wavenumber(i)*(dJ1yR(i)+(y2(i)*R2-1)*J1y(i));
        M(1,3) = dJ1yR(i)-J1y(i);
        M(2,1) = 2*(dJ1xR(i)-J1x(i));
        M(2,2) = 2*Wavenumber(i)*(J1y(i)-dJ1yR(i));
        M(2,3) = 2*dJ1yR(i)+(y2(i)*R2-2)*J1y(i);
        M(3,1) = 2*Wavenumber(i)*dJ1xR(i);
        M(3,2) = (y2(i)-k2(i))*dJ1yR(i);
        M(3,3) = -Wavenumber(i)*J1y(i);
    end
    Y(i) = abs(det(M));
    if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
        XRough = SweepRange(i-1);
        break
    end
end
% figure,plot(20*log10(Y))

% fine search for the undamped F(1,1)/FScholte(1,1) mode on the real axis
Bisections = ceil(log2(Resolution/(SweepRange(2)-SweepRange(1)))/log2(2*.25));
PhaseVelocity = [XRough-(SweepRange(2)-SweepRange(1)) XRough+(SweepRange(2)-SweepRange(1))];
for o = 1:Bisections
    PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
    for i = 1:length(PhaseVelocity)
        Wavenumber = AngularFrequency/PhaseVelocity(i);
        k2 = Wavenumber^2;
        y2 = kT2-k2;
        xR = sqrt(kL2-k2)*R;
        yR = sqrt(kT2-k2)*R;
        J1x = besselj(1,xR);
        J1y = besselj(1,yR);
        dJ1xR = J1x-xR*besselj(2,xR);
        dJ1yR = J1y-yR*besselj(2,yR);
        if  FluidLoading
            zR = sqrt(kF2-k2)*R;
            H1z = besselh(1,zR);
            dH1zR = H1z-zR*besselh(2,zR);
            M(1,1) = dJ1xR+(.5*kT2*R2-k2*R2-1)*J1x;
            M(1,2) = -Wavenumber*(dJ1yR+(y2*R2-1)*J1y);
            M(1,3) = dJ1yR-J1y;
            M(1,4) = .5*kT2*R2*Density*H1z;
            M(2,1) = 2*(dJ1xR-J1x);
            M(2,2) = 2*Wavenumber*(J1y-dJ1yR);
            M(2,3) = 2*dJ1yR+(y2*R2-2)*J1y;
            M(3,1) = 2*Wavenumber*dJ1xR;
            M(3,2) = (y2-k2)*dJ1yR;
            M(3,3) = -Wavenumber*J1y;
            M(4,1) = dJ1xR;
            M(4,2) = -Wavenumber*dJ1yR;
            M(4,3) = -J1y;
            M(4,4) = dH1zR;
        else
            M(1,1) = dJ1xR+(.5*kT2*R2-k2*R2-1)*J1x;
            M(1,2) = -Wavenumber*(dJ1yR+(y2*R2-1)*J1y);
            M(1,3) = dJ1yR-J1y;
            M(2,1) = 2*(dJ1xR-J1x);
            M(2,2) = 2*Wavenumber*(J1y-dJ1yR);
            M(2,3) = 2*dJ1yR+(y2*R2-2)*J1y;
            M(3,1) = 2*Wavenumber*dJ1xR;
            M(3,2) = (y2-k2)*dJ1yR;
            M(3,3) = -Wavenumber*J1y;
        end
        Y(i) = abs(det(M));
        if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
            if  o == Bisections
                FLamb = PhaseVelocity(i-1);
            end
            PhaseVelocity = [PhaseVelocity(i-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(i-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
            break
        end
    end
end