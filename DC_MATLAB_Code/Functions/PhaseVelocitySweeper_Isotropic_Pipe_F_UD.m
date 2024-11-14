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
function FLamb = PhaseVelocitySweeper_Isotropic_Pipe_F_UD(Material,Ro,Ri,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Frequency)
%#ok<*AGROW>
Resolution = 1e-6; % (m/s)
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityi = InnerFluid.Density/Material.Density;
Densityo = OuterFluid.Density/Material.Density;
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
kInnerFluid2 = (AngularFrequency/InnerFluid.Velocity)^2;
kOuterFluid2 = (AngularFrequency/OuterFluid.Velocity)^2;

% rough search for the undamped F(1,1)/FScholte(1,1) mode on the real axis
SweepRange = 1-1e-10:Material.TransverseVelocity+2;
Wavenumber = AngularFrequency./SweepRange;
k2 = Wavenumber.^2;
x = sqrt(k2-kL2);
y = sqrt(k2-kT2);
y2 = kT2-k2;
xRo = x*Ro;
yRo = y*Ro;
xRi = x*Ri;
yRi = y*Ri;
Z1xi = besseli(1,xRi);
Z1xo = besseli(1,xRo);
Z1yi = besseli(1,yRi);
Z1yo = besseli(1,yRo);
W1xi = besselk(1,xRi);
W1xo = besselk(1,xRo);
W1yi = besselk(1,yRi);
W1yo = besselk(1,yRo);
dZ1xiRi = Z1xi+xRi.*besseli(2,xRi);
dZ1xoRo = Z1xo+xRo.*besseli(2,xRo);
dZ1yiRi = Z1yi+yRi.*besseli(2,yRi);
dZ1yoRo = Z1yo+yRo.*besseli(2,yRo);
dW1xiRi = W1xi-xRi.*besselk(2,xRi);
dW1xoRo = W1xo-xRo.*besselk(2,xRo);
dW1yiRi = W1yi-yRi.*besselk(2,yRi);
dW1yoRo = W1yo-yRo.*besselk(2,yRo);
if  ToggleInnerFluid
    zRi = sqrt(kInnerFluid2-k2)*Ri;
    if  Sink
        Z1zi = besselh(1,2,zRi);
        dZ1ziRi = Z1zi-zRi.*besselh(2,2,zRi);
    else
        Z1zi = besselj(1,zRi);
        dZ1ziRi = Z1zi-zRi.*besselj(2,zRi);
    end
end
if  ToggleOuterFluid
    zRo = sqrt(kOuterFluid2-k2)*Ro;
    H1zo = besselh(1,zRo);
    dH1zoRo = H1zo-zRo.*besselh(2,zRo);
end
for i = 1:length(SweepRange)
    if  ToggleInnerFluid && ToggleOuterFluid
        M(1,1) = dZ1xiRi(i)+(.5*kT2*Ri2-k2(i)*Ri2-1)*Z1xi(i);
        M(1,2) = dW1xiRi(i)+(.5*kT2*Ri2-k2(i)*Ri2-1)*W1xi(i);
        M(1,3) = -Wavenumber(i)*(dZ1yiRi(i)+(y2(i)*Ri2-1)*Z1yi(i));
        M(1,4) = -Wavenumber(i)*(dW1yiRi(i)+(y2(i)*Ri2-1)*W1yi(i));
        M(1,5) = dZ1yiRi(i)-Z1yi(i);
        M(1,6) = dW1yiRi(i)-W1yi(i);
        M(1,7) = .5*kT2*Ri2*Densityi*Z1zi(i);
        M(2,1) = 2*(dZ1xiRi(i)-Z1xi(i));
        M(2,2) = 2*(dW1xiRi(i)-W1xi(i));
        M(2,3) = 2*Wavenumber(i)*(Z1yi(i)-dZ1yiRi(i));
        M(2,4) = 2*Wavenumber(i)*(W1yi(i)-dW1yiRi(i));
        M(2,5) = 2*dZ1yiRi(i)+(y2(i)*Ri2-2)*Z1yi(i);
        M(2,6) = 2*dW1yiRi(i)+(y2(i)*Ri2-2)*W1yi(i);
        M(3,1) = 2*Wavenumber(i)*dZ1xiRi(i);
        M(3,2) = 2*Wavenumber(i)*dW1xiRi(i);
        M(3,3) = (y2(i)-k2(i))*dZ1yiRi(i);
        M(3,4) = (y2(i)-k2(i))*dW1yiRi(i);
        M(3,5) = -Wavenumber(i)*Z1yi(i);
        M(3,6) = -Wavenumber(i)*W1yi(i);
        M(4,1) = dZ1xiRi(i);
        M(4,2) = dW1xiRi(i);
        M(4,3) = -Wavenumber(i)*dZ1yiRi(i);
        M(4,4) = -Wavenumber(i)*dW1yiRi(i);
        M(4,5) = -Z1yi(i);
        M(4,6) = -W1yi(i);
        M(4,7) = dZ1ziRi(i);
        M(5,1) = dZ1xoRo(i)+(.5*kT2*Ro2-k2(i)*Ro2-1)*Z1xo(i);
        M(5,2) = dW1xoRo(i)+(.5*kT2*Ro2-k2(i)*Ro2-1)*W1xo(i);
        M(5,3) = -Wavenumber(i)*(dZ1yoRo(i)+(y2(i)*Ro2-1)*Z1yo(i));
        M(5,4) = -Wavenumber(i)*(dW1yoRo(i)+(y2(i)*Ro2-1)*W1yo(i));
        M(5,5) = dZ1yoRo(i)-Z1yo(i);
        M(5,6) = dW1yoRo(i)-W1yo(i);
        M(5,8) = .5*kT2*Ro2*Densityo*H1zo(i);
        M(6,1) = 2*(dZ1xoRo(i)-Z1xo(i));
        M(6,2) = 2*(dW1xoRo(i)-W1xo(i));
        M(6,3) = 2*Wavenumber(i)*(Z1yo(i)-dZ1yoRo(i));
        M(6,4) = 2*Wavenumber(i)*(W1yo(i)-dW1yoRo(i));
        M(6,5) = 2*dZ1yoRo(i)+(y2(i)*Ro2-2)*Z1yo(i);
        M(6,6) = 2*dW1yoRo(i)+(y2(i)*Ro2-2)*W1yo(i);
        M(7,1) = 2*Wavenumber(i)*dZ1xoRo(i);
        M(7,2) = 2*Wavenumber(i)*dW1xoRo(i);
        M(7,3) = (y2(i)-k2(i))*dZ1yoRo(i);
        M(7,4) = (y2(i)-k2(i))*dW1yoRo(i);
        M(7,5) = -Wavenumber(i)*Z1yo(i);
        M(7,6) = -Wavenumber(i)*W1yo(i);
        M(8,1) = dZ1xoRo(i);
        M(8,2) = dW1xoRo(i);
        M(8,3) = -Wavenumber(i)*dZ1yoRo(i);
        M(8,4) = -Wavenumber(i)*dW1yoRo(i);
        M(8,5) = -Z1yo(i);
        M(8,6) = -W1yo(i);
        M(8,8) = dH1zoRo(i);
    elseif ToggleInnerFluid && ~ToggleOuterFluid
        M(1,1) = dZ1xiRi(i)+(.5*kT2*Ri2-k2(i)*Ri2-1)*Z1xi(i);
        M(1,2) = dW1xiRi(i)+(.5*kT2*Ri2-k2(i)*Ri2-1)*W1xi(i);
        M(1,3) = -Wavenumber(i)*(dZ1yiRi(i)+(y2(i)*Ri2-1)*Z1yi(i));
        M(1,4) = -Wavenumber(i)*(dW1yiRi(i)+(y2(i)*Ri2-1)*W1yi(i));
        M(1,5) = dZ1yiRi(i)-Z1yi(i);
        M(1,6) = dW1yiRi(i)-W1yi(i);
        M(1,7) = .5*kT2*Ri2*Densityi*Z1zi(i);
        M(2,1) = 2*(dZ1xiRi(i)-Z1xi(i));
        M(2,2) = 2*(dW1xiRi(i)-W1xi(i));
        M(2,3) = 2*Wavenumber(i)*(Z1yi(i)-dZ1yiRi(i));
        M(2,4) = 2*Wavenumber(i)*(W1yi(i)-dW1yiRi(i));
        M(2,5) = 2*dZ1yiRi(i)+(y2(i)*Ri2-2)*Z1yi(i);
        M(2,6) = 2*dW1yiRi(i)+(y2(i)*Ri2-2)*W1yi(i);
        M(3,1) = 2*Wavenumber(i)*dZ1xiRi(i);
        M(3,2) = 2*Wavenumber(i)*dW1xiRi(i);
        M(3,3) = (y2(i)-k2(i))*dZ1yiRi(i);
        M(3,4) = (y2(i)-k2(i))*dW1yiRi(i);
        M(3,5) = -Wavenumber(i)*Z1yi(i);
        M(3,6) = -Wavenumber(i)*W1yi(i);
        M(4,1) = dZ1xiRi(i);
        M(4,2) = dW1xiRi(i);
        M(4,3) = -Wavenumber(i)*dZ1yiRi(i);
        M(4,4) = -Wavenumber(i)*dW1yiRi(i);
        M(4,5) = -Z1yi(i);
        M(4,6) = -W1yi(i);
        M(4,7) = dZ1ziRi(i);
        M(5,1) = dZ1xoRo(i)+(.5*kT2*Ro2-k2(i)*Ro2-1)*Z1xo(i);
        M(5,2) = dW1xoRo(i)+(.5*kT2*Ro2-k2(i)*Ro2-1)*W1xo(i);
        M(5,3) = -Wavenumber(i)*(dZ1yoRo(i)+(y2(i)*Ro2-1)*Z1yo(i));
        M(5,4) = -Wavenumber(i)*(dW1yoRo(i)+(y2(i)*Ro2-1)*W1yo(i));
        M(5,5) = dZ1yoRo(i)-Z1yo(i);
        M(5,6) = dW1yoRo(i)-W1yo(i);
        M(6,1) = 2*(dZ1xoRo(i)-Z1xo(i));
        M(6,2) = 2*(dW1xoRo(i)-W1xo(i));
        M(6,3) = 2*Wavenumber(i)*(Z1yo(i)-dZ1yoRo(i));
        M(6,4) = 2*Wavenumber(i)*(W1yo(i)-dW1yoRo(i));
        M(6,5) = 2*dZ1yoRo(i)+(y2(i)*Ro2-2)*Z1yo(i);
        M(6,6) = 2*dW1yoRo(i)+(y2(i)*Ro2-2)*W1yo(i);
        M(7,1) = 2*Wavenumber(i)*dZ1xoRo(i);
        M(7,2) = 2*Wavenumber(i)*dW1xoRo(i);
        M(7,3) = (y2(i)-k2(i))*dZ1yoRo(i);
        M(7,4) = (y2(i)-k2(i))*dW1yoRo(i);
        M(7,5) = -Wavenumber(i)*Z1yo(i);
        M(7,6) = -Wavenumber(i)*W1yo(i);
    elseif ~ToggleInnerFluid && ToggleOuterFluid
        M(1,1) = dZ1xiRi(i)+(.5*kT2*Ri2-k2(i)*Ri2-1)*Z1xi(i);
        M(1,2) = dW1xiRi(i)+(.5*kT2*Ri2-k2(i)*Ri2-1)*W1xi(i);
        M(1,3) = -Wavenumber(i)*(dZ1yiRi(i)+(y2(i)*Ri2-1)*Z1yi(i));
        M(1,4) = -Wavenumber(i)*(dW1yiRi(i)+(y2(i)*Ri2-1)*W1yi(i));
        M(1,5) = dZ1yiRi(i)-Z1yi(i);
        M(1,6) = dW1yiRi(i)-W1yi(i);
        M(2,1) = 2*(dZ1xiRi(i)-Z1xi(i));
        M(2,2) = 2*(dW1xiRi(i)-W1xi(i));
        M(2,3) = 2*Wavenumber(i)*(Z1yi(i)-dZ1yiRi(i));
        M(2,4) = 2*Wavenumber(i)*(W1yi(i)-dW1yiRi(i));
        M(2,5) = 2*dZ1yiRi(i)+(y2(i)*Ri2-2)*Z1yi(i);
        M(2,6) = 2*dW1yiRi(i)+(y2(i)*Ri2-2)*W1yi(i);
        M(3,1) = 2*Wavenumber(i)*dZ1xiRi(i);
        M(3,2) = 2*Wavenumber(i)*dW1xiRi(i);
        M(3,3) = (y2(i)-k2(i))*dZ1yiRi(i);
        M(3,4) = (y2(i)-k2(i))*dW1yiRi(i);
        M(3,5) = -Wavenumber(i)*Z1yi(i);
        M(3,6) = -Wavenumber(i)*W1yi(i);
        M(4,1) = dZ1xoRo(i)+(.5*kT2*Ro2-k2(i)*Ro2-1)*Z1xo(i);
        M(4,2) = dW1xoRo(i)+(.5*kT2*Ro2-k2(i)*Ro2-1)*W1xo(i);
        M(4,3) = -Wavenumber(i)*(dZ1yoRo(i)+(y2(i)*Ro2-1)*Z1yo(i));
        M(4,4) = -Wavenumber(i)*(dW1yoRo(i)+(y2(i)*Ro2-1)*W1yo(i));
        M(4,5) = dZ1yoRo(i)-Z1yo(i);
        M(4,6) = dW1yoRo(i)-W1yo(i);
        M(4,7) = .5*kT2*Ro2*Densityo*H1zo(i);
        M(5,1) = 2*(dZ1xoRo(i)-Z1xo(i));
        M(5,2) = 2*(dW1xoRo(i)-W1xo(i));
        M(5,3) = 2*Wavenumber(i)*(Z1yo(i)-dZ1yoRo(i));
        M(5,4) = 2*Wavenumber(i)*(W1yo(i)-dW1yoRo(i));
        M(5,5) = 2*dZ1yoRo(i)+(y2(i)*Ro2-2)*Z1yo(i);
        M(5,6) = 2*dW1yoRo(i)+(y2(i)*Ro2-2)*W1yo(i);
        M(6,1) = 2*Wavenumber(i)*dZ1xoRo(i);
        M(6,2) = 2*Wavenumber(i)*dW1xoRo(i);
        M(6,3) = (y2(i)-k2(i))*dZ1yoRo(i);
        M(6,4) = (y2(i)-k2(i))*dW1yoRo(i);
        M(6,5) = -Wavenumber(i)*Z1yo(i);
        M(6,6) = -Wavenumber(i)*W1yo(i);
        M(7,1) = dZ1xoRo(i);
        M(7,2) = dW1xoRo(i);
        M(7,3) = -Wavenumber(i)*dZ1yoRo(i);
        M(7,4) = -Wavenumber(i)*dW1yoRo(i);
        M(7,5) = -Z1yo(i);
        M(7,6) = -W1yo(i);
        M(7,7) = dH1zoRo(i);
    elseif ~ToggleInnerFluid && ~ToggleOuterFluid
        M(1,1) = dZ1xiRi(i)+(.5*kT2*Ri2-k2(i)*Ri2-1)*Z1xi(i);
        M(1,2) = dW1xiRi(i)+(.5*kT2*Ri2-k2(i)*Ri2-1)*W1xi(i);
        M(1,3) = -Wavenumber(i)*(dZ1yiRi(i)+(y2(i)*Ri2-1)*Z1yi(i));
        M(1,4) = -Wavenumber(i)*(dW1yiRi(i)+(y2(i)*Ri2-1)*W1yi(i));
        M(1,5) = dZ1yiRi(i)-Z1yi(i);
        M(1,6) = dW1yiRi(i)-W1yi(i);
        M(2,1) = 2*(dZ1xiRi(i)-Z1xi(i));
        M(2,2) = 2*(dW1xiRi(i)-W1xi(i));
        M(2,3) = 2*Wavenumber(i)*(Z1yi(i)-dZ1yiRi(i));
        M(2,4) = 2*Wavenumber(i)*(W1yi(i)-dW1yiRi(i));
        M(2,5) = 2*dZ1yiRi(i)+(y2(i)*Ri2-2)*Z1yi(i);
        M(2,6) = 2*dW1yiRi(i)+(y2(i)*Ri2-2)*W1yi(i);
        M(3,1) = 2*Wavenumber(i)*dZ1xiRi(i);
        M(3,2) = 2*Wavenumber(i)*dW1xiRi(i);
        M(3,3) = (y2(i)-k2(i))*dZ1yiRi(i);
        M(3,4) = (y2(i)-k2(i))*dW1yiRi(i);
        M(3,5) = -Wavenumber(i)*Z1yi(i);
        M(3,6) = -Wavenumber(i)*W1yi(i);
        M(4,1) = dZ1xoRo(i)+(.5*kT2*Ro2-k2(i)*Ro2-1)*Z1xo(i);
        M(4,2) = dW1xoRo(i)+(.5*kT2*Ro2-k2(i)*Ro2-1)*W1xo(i);
        M(4,3) = -Wavenumber(i)*(dZ1yoRo(i)+(y2(i)*Ro2-1)*Z1yo(i));
        M(4,4) = -Wavenumber(i)*(dW1yoRo(i)+(y2(i)*Ro2-1)*W1yo(i));
        M(4,5) = dZ1yoRo(i)-Z1yo(i);
        M(4,6) = dW1yoRo(i)-W1yo(i);
        M(5,1) = 2*(dZ1xoRo(i)-Z1xo(i));
        M(5,2) = 2*(dW1xoRo(i)-W1xo(i));
        M(5,3) = 2*Wavenumber(i)*(Z1yo(i)-dZ1yoRo(i));
        M(5,4) = 2*Wavenumber(i)*(W1yo(i)-dW1yoRo(i));
        M(5,5) = 2*dZ1yoRo(i)+(y2(i)*Ro2-2)*Z1yo(i);
        M(5,6) = 2*dW1yoRo(i)+(y2(i)*Ro2-2)*W1yo(i);
        M(6,1) = 2*Wavenumber(i)*dZ1xoRo(i);
        M(6,2) = 2*Wavenumber(i)*dW1xoRo(i);
        M(6,3) = (y2(i)-k2(i))*dZ1yoRo(i);
        M(6,4) = (y2(i)-k2(i))*dW1yoRo(i);
        M(6,5) = -Wavenumber(i)*Z1yo(i);
        M(6,6) = -Wavenumber(i)*W1yo(i);
    end
    Y(i) = abs(det(M));
    if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
        XRough = SweepRange(i-1);
        break
    end
end
% figure,plot(20*log10(Y))

% rough search for the undamped F(1,1)/FScholte(1,1) mode on the real axis
Bisections = ceil(log2(Resolution/(SweepRange(2)-SweepRange(1)))/log2(2*.25));
PhaseVelocity = [XRough-(SweepRange(2)-SweepRange(1)) XRough+(SweepRange(2)-SweepRange(1))];
for o = 1:Bisections
    PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
    for i = 1:length(PhaseVelocity)
        Wavenumber = AngularFrequency/PhaseVelocity(i);
        k2 = Wavenumber^2;
        x = sqrt(k2-kL2);
        y = sqrt(k2-kT2);
        y2 = kT2-k2;
        xRo = x*Ro;
        yRo = y*Ro;
        xRi = x*Ri;
        yRi = y*Ri;
        Z1xi = besseli(1,xRi);
        Z1xo = besseli(1,xRo);
        Z1yi = besseli(1,yRi);
        Z1yo = besseli(1,yRo);
        W1xi = besselk(1,xRi);
        W1xo = besselk(1,xRo);
        W1yi = besselk(1,yRi);
        W1yo = besselk(1,yRo);
        dZ1xiRi = Z1xi+xRi*besseli(2,xRi);
        dZ1xoRo = Z1xo+xRo*besseli(2,xRo);
        dZ1yiRi = Z1yi+yRi*besseli(2,yRi);
        dZ1yoRo = Z1yo+yRo*besseli(2,yRo);
        dW1xiRi = W1xi-xRi*besselk(2,xRi);
        dW1xoRo = W1xo-xRo*besselk(2,xRo);
        dW1yiRi = W1yi-yRi*besselk(2,yRi);
        dW1yoRo = W1yo-yRo*besselk(2,yRo);
        if  ToggleInnerFluid
            zRi = sqrt(kInnerFluid2-k2)*Ri;
            if  Sink
                Z1zi = besselh(1,2,zRi);
                dZ1ziRi = Z1zi-zRi*besselh(2,2,zRi);
            else
                Z1zi = besselj(1,zRi);
                dZ1ziRi = Z1zi-zRi*besselj(2,zRi);
            end
        end
        if  ToggleOuterFluid
            zRo = sqrt(kOuterFluid2-k2)*Ro;
            H1zo = besselh(1,zRo);
            dH1zoRo = H1zo-zRo*besselh(2,zRo);
        end
        if  ToggleInnerFluid && ToggleOuterFluid
            M(1,1) = dZ1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*Z1xi;
            M(1,2) = dW1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*W1xi;
            M(1,3) = -Wavenumber*(dZ1yiRi+(y2*Ri2-1)*Z1yi);
            M(1,4) = -Wavenumber*(dW1yiRi+(y2*Ri2-1)*W1yi);
            M(1,5) = dZ1yiRi-Z1yi;
            M(1,6) = dW1yiRi-W1yi;
            M(1,7) = .5*kT2*Ri2*Densityi*Z1zi;
            M(2,1) = 2*(dZ1xiRi-Z1xi);
            M(2,2) = 2*(dW1xiRi-W1xi);
            M(2,3) = 2*Wavenumber*(Z1yi-dZ1yiRi);
            M(2,4) = 2*Wavenumber*(W1yi-dW1yiRi);
            M(2,5) = 2*dZ1yiRi+(y2*Ri2-2)*Z1yi;
            M(2,6) = 2*dW1yiRi+(y2*Ri2-2)*W1yi;
            M(3,1) = 2*Wavenumber*dZ1xiRi;
            M(3,2) = 2*Wavenumber*dW1xiRi;
            M(3,3) = (y2-k2)*dZ1yiRi;
            M(3,4) = (y2-k2)*dW1yiRi;
            M(3,5) = -Wavenumber*Z1yi;
            M(3,6) = -Wavenumber*W1yi;
            M(4,1) = dZ1xiRi;
            M(4,2) = dW1xiRi;
            M(4,3) = -Wavenumber*dZ1yiRi;
            M(4,4) = -Wavenumber*dW1yiRi;
            M(4,5) = -Z1yi;
            M(4,6) = -W1yi;
            M(4,7) = dZ1ziRi;
            M(5,1) = dZ1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*Z1xo;
            M(5,2) = dW1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*W1xo;
            M(5,3) = -Wavenumber*(dZ1yoRo+(y2*Ro2-1)*Z1yo);
            M(5,4) = -Wavenumber*(dW1yoRo+(y2*Ro2-1)*W1yo);
            M(5,5) = dZ1yoRo-Z1yo;
            M(5,6) = dW1yoRo-W1yo;
            M(5,8) = .5*kT2*Ro2*Densityo*H1zo;
            M(6,1) = 2*(dZ1xoRo-Z1xo);
            M(6,2) = 2*(dW1xoRo-W1xo);
            M(6,3) = 2*Wavenumber*(Z1yo-dZ1yoRo);
            M(6,4) = 2*Wavenumber*(W1yo-dW1yoRo);
            M(6,5) = 2*dZ1yoRo+(y2*Ro2-2)*Z1yo;
            M(6,6) = 2*dW1yoRo+(y2*Ro2-2)*W1yo;
            M(7,1) = 2*Wavenumber*dZ1xoRo;
            M(7,2) = 2*Wavenumber*dW1xoRo;
            M(7,3) = (y2-k2)*dZ1yoRo;
            M(7,4) = (y2-k2)*dW1yoRo;
            M(7,5) = -Wavenumber*Z1yo;
            M(7,6) = -Wavenumber*W1yo;
            M(8,1) = dZ1xoRo;
            M(8,2) = dW1xoRo;
            M(8,3) = -Wavenumber*dZ1yoRo;
            M(8,4) = -Wavenumber*dW1yoRo;
            M(8,5) = -Z1yo;
            M(8,6) = -W1yo;
            M(8,8) = dH1zoRo;
        elseif ToggleInnerFluid && ~ToggleOuterFluid
            M(1,1) = dZ1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*Z1xi;
            M(1,2) = dW1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*W1xi;
            M(1,3) = -Wavenumber*(dZ1yiRi+(y2*Ri2-1)*Z1yi);
            M(1,4) = -Wavenumber*(dW1yiRi+(y2*Ri2-1)*W1yi);
            M(1,5) = dZ1yiRi-Z1yi;
            M(1,6) = dW1yiRi-W1yi;
            M(1,7) = .5*kT2*Ri2*Densityi*Z1zi;
            M(2,1) = 2*(dZ1xiRi-Z1xi);
            M(2,2) = 2*(dW1xiRi-W1xi);
            M(2,3) = 2*Wavenumber*(Z1yi-dZ1yiRi);
            M(2,4) = 2*Wavenumber*(W1yi-dW1yiRi);
            M(2,5) = 2*dZ1yiRi+(y2*Ri2-2)*Z1yi;
            M(2,6) = 2*dW1yiRi+(y2*Ri2-2)*W1yi;
            M(3,1) = 2*Wavenumber*dZ1xiRi;
            M(3,2) = 2*Wavenumber*dW1xiRi;
            M(3,3) = (y2-k2)*dZ1yiRi;
            M(3,4) = (y2-k2)*dW1yiRi;
            M(3,5) = -Wavenumber*Z1yi;
            M(3,6) = -Wavenumber*W1yi;
            M(4,1) = dZ1xiRi;
            M(4,2) = dW1xiRi;
            M(4,3) = -Wavenumber*dZ1yiRi;
            M(4,4) = -Wavenumber*dW1yiRi;
            M(4,5) = -Z1yi;
            M(4,6) = -W1yi;
            M(4,7) = dZ1ziRi;
            M(5,1) = dZ1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*Z1xo;
            M(5,2) = dW1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*W1xo;
            M(5,3) = -Wavenumber*(dZ1yoRo+(y2*Ro2-1)*Z1yo);
            M(5,4) = -Wavenumber*(dW1yoRo+(y2*Ro2-1)*W1yo);
            M(5,5) = dZ1yoRo-Z1yo;
            M(5,6) = dW1yoRo-W1yo;
            M(6,1) = 2*(dZ1xoRo-Z1xo);
            M(6,2) = 2*(dW1xoRo-W1xo);
            M(6,3) = 2*Wavenumber*(Z1yo-dZ1yoRo);
            M(6,4) = 2*Wavenumber*(W1yo-dW1yoRo);
            M(6,5) = 2*dZ1yoRo+(y2*Ro2-2)*Z1yo;
            M(6,6) = 2*dW1yoRo+(y2*Ro2-2)*W1yo;
            M(7,1) = 2*Wavenumber*dZ1xoRo;
            M(7,2) = 2*Wavenumber*dW1xoRo;
            M(7,3) = (y2-k2)*dZ1yoRo;
            M(7,4) = (y2-k2)*dW1yoRo;
            M(7,5) = -Wavenumber*Z1yo;
            M(7,6) = -Wavenumber*W1yo;
        elseif ~ToggleInnerFluid && ToggleOuterFluid
            M(1,1) = dZ1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*Z1xi;
            M(1,2) = dW1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*W1xi;
            M(1,3) = -Wavenumber*(dZ1yiRi+(y2*Ri2-1)*Z1yi);
            M(1,4) = -Wavenumber*(dW1yiRi+(y2*Ri2-1)*W1yi);
            M(1,5) = dZ1yiRi-Z1yi;
            M(1,6) = dW1yiRi-W1yi;
            M(2,1) = 2*(dZ1xiRi-Z1xi);
            M(2,2) = 2*(dW1xiRi-W1xi);
            M(2,3) = 2*Wavenumber*(Z1yi-dZ1yiRi);
            M(2,4) = 2*Wavenumber*(W1yi-dW1yiRi);
            M(2,5) = 2*dZ1yiRi+(y2*Ri2-2)*Z1yi;
            M(2,6) = 2*dW1yiRi+(y2*Ri2-2)*W1yi;
            M(3,1) = 2*Wavenumber*dZ1xiRi;
            M(3,2) = 2*Wavenumber*dW1xiRi;
            M(3,3) = (y2-k2)*dZ1yiRi;
            M(3,4) = (y2-k2)*dW1yiRi;
            M(3,5) = -Wavenumber*Z1yi;
            M(3,6) = -Wavenumber*W1yi;
            M(4,1) = dZ1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*Z1xo;
            M(4,2) = dW1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*W1xo;
            M(4,3) = -Wavenumber*(dZ1yoRo+(y2*Ro2-1)*Z1yo);
            M(4,4) = -Wavenumber*(dW1yoRo+(y2*Ro2-1)*W1yo);
            M(4,5) = dZ1yoRo-Z1yo;
            M(4,6) = dW1yoRo-W1yo;
            M(4,7) = .5*kT2*Ro2*Densityo*H1zo;
            M(5,1) = 2*(dZ1xoRo-Z1xo);
            M(5,2) = 2*(dW1xoRo-W1xo);
            M(5,3) = 2*Wavenumber*(Z1yo-dZ1yoRo);
            M(5,4) = 2*Wavenumber*(W1yo-dW1yoRo);
            M(5,5) = 2*dZ1yoRo+(y2*Ro2-2)*Z1yo;
            M(5,6) = 2*dW1yoRo+(y2*Ro2-2)*W1yo;
            M(6,1) = 2*Wavenumber*dZ1xoRo;
            M(6,2) = 2*Wavenumber*dW1xoRo;
            M(6,3) = (y2-k2)*dZ1yoRo;
            M(6,4) = (y2-k2)*dW1yoRo;
            M(6,5) = -Wavenumber*Z1yo;
            M(6,6) = -Wavenumber*W1yo;
            M(7,1) = dZ1xoRo;
            M(7,2) = dW1xoRo;
            M(7,3) = -Wavenumber*dZ1yoRo;
            M(7,4) = -Wavenumber*dW1yoRo;
            M(7,5) = -Z1yo;
            M(7,6) = -W1yo;
            M(7,7) = dH1zoRo;
        elseif ~ToggleInnerFluid && ~ToggleOuterFluid
            M(1,1) = dZ1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*Z1xi;
            M(1,2) = dW1xiRi+(.5*kT2*Ri2-k2*Ri2-1)*W1xi;
            M(1,3) = -Wavenumber*(dZ1yiRi+(y2*Ri2-1)*Z1yi);
            M(1,4) = -Wavenumber*(dW1yiRi+(y2*Ri2-1)*W1yi);
            M(1,5) = dZ1yiRi-Z1yi;
            M(1,6) = dW1yiRi-W1yi;
            M(2,1) = 2*(dZ1xiRi-Z1xi);
            M(2,2) = 2*(dW1xiRi-W1xi);
            M(2,3) = 2*Wavenumber*(Z1yi-dZ1yiRi);
            M(2,4) = 2*Wavenumber*(W1yi-dW1yiRi);
            M(2,5) = 2*dZ1yiRi+(y2*Ri2-2)*Z1yi;
            M(2,6) = 2*dW1yiRi+(y2*Ri2-2)*W1yi;
            M(3,1) = 2*Wavenumber*dZ1xiRi;
            M(3,2) = 2*Wavenumber*dW1xiRi;
            M(3,3) = (y2-k2)*dZ1yiRi;
            M(3,4) = (y2-k2)*dW1yiRi;
            M(3,5) = -Wavenumber*Z1yi;
            M(3,6) = -Wavenumber*W1yi;
            M(4,1) = dZ1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*Z1xo;
            M(4,2) = dW1xoRo+(.5*kT2*Ro2-k2*Ro2-1)*W1xo;
            M(4,3) = -Wavenumber*(dZ1yoRo+(y2*Ro2-1)*Z1yo);
            M(4,4) = -Wavenumber*(dW1yoRo+(y2*Ro2-1)*W1yo);
            M(4,5) = dZ1yoRo-Z1yo;
            M(4,6) = dW1yoRo-W1yo;
            M(5,1) = 2*(dZ1xoRo-Z1xo);
            M(5,2) = 2*(dW1xoRo-W1xo);
            M(5,3) = 2*Wavenumber*(Z1yo-dZ1yoRo);
            M(5,4) = 2*Wavenumber*(W1yo-dW1yoRo);
            M(5,5) = 2*dZ1yoRo+(y2*Ro2-2)*Z1yo;
            M(5,6) = 2*dW1yoRo+(y2*Ro2-2)*W1yo;
            M(6,1) = 2*Wavenumber*dZ1xoRo;
            M(6,2) = 2*Wavenumber*dW1xoRo;
            M(6,3) = (y2-k2)*dZ1yoRo;
            M(6,4) = (y2-k2)*dW1yoRo;
            M(6,5) = -Wavenumber*Z1yo;
            M(6,6) = -Wavenumber*W1yo;
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