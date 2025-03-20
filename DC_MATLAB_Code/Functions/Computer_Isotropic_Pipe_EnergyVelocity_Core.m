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
function [X,Counter] = Computer_Isotropic_Pipe_EnergyVelocity_Core(Multithreading,X,ModeType,n,ModeTotal,Counter,h,Material,OuterFluid,InnerFluid,Sink,Ro,Ri,SamplesR)
%#ok<*AGROW>
%#ok<*PFBNS>
%#ok<*PFTIN>
%#ok<*PFTUSW>
%#ok<*GVMIS>
global Stop 
Stop = 0;
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
Ro2 = Ro^2;
Ri2 = Ri^2;
r = (Ri:(Ro-Ri)/SamplesR:Ro)';
if  Sink
    rInnerFluid = (Ri/SamplesR:Ri/SamplesR:Ri)';
else
    rInnerFluid = (0:Ri/SamplesR:Ri)';
    rInnerFluid(1) = 1e-10;
end
rTotal = [rInnerFluid;r];
r2 = r.^2;
rInnerFluid2 = rInnerFluid.^2;
rTotal2 = rTotal.^2;
Diff = diff(r2)';
DiffTotal = diff(rTotal2)';
for p = 1:length(X)
    if  contains(ModeType,'F')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            Attenuation = X{p}(:,6).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                AngularFrequency2 = AngularFrequency(q)^2;
                k = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                k2 = k^2;
                kL2 = (AngularFrequency(q)/Material.LongitudinalVelocity_complex)^2; 
                kT2 = (AngularFrequency(q)/Material.TransverseVelocity_complex)^2;
                if  PhaseVelocity(q) > Material.LongitudinalVelocity
                    x = sqrt(kL2-k2); 
                    y = sqrt(kT2-k2);
                elseif PhaseVelocity(q) <= Material.LongitudinalVelocity && PhaseVelocity(q) > Material.TransverseVelocity
                    x = sqrt(k2-kL2);
                    y = sqrt(kT2-k2);
                elseif PhaseVelocity(q) <= Material.TransverseVelocity
                    x = sqrt(k2-kL2);
                    y = sqrt(k2-kT2);
                end
                zi = sqrt((AngularFrequency(q)/InnerFluid.Velocity)^2-k2);
                zo = sqrt((AngularFrequency(q)/OuterFluid.Velocity)^2-k2);
                if  strcmp(ModeType,'FScholte') && Attenuation(q) ~= 0
                    if  PhaseVelocity(q) < InnerFluid.Velocity
                        zi = -zi;
                    end
                    if  PhaseVelocity(q) < OuterFluid.Velocity
                        zo = -zo;
                    end
                end
                x2 = kL2-k2;
                y2 = kT2-k2;
                xr = x*r;
                yr = y*r;
                zri = zi*rInnerFluid;
                zRo = zo*Ro;
                if  PhaseVelocity(q) > Material.LongitudinalVelocity
                    Znx = besselj(n,xr);
                    Zny = besselj(n,yr);
                    Wnx = bessely(n,xr);
                    Wny = bessely(n,yr);
                    Zn1x = x*besselj(n+1,xr);
                    Zn1y = y*besselj(n+1,yr);
                    Wn1x = x*bessely(n+1,xr);
                    Wn1y = y*bessely(n+1,yr);
                    dZnx = n*Znx./r-Zn1x;
                    dZny = n*Zny./r-Zn1y;
                    dWnx = n*Wnx./r-Wn1x;
                    dWny = n*Wny./r-Wn1y;
                    dg1Zn1y = -(n+1)*Zn1y./r+y2*Zny;
                    dg1Wn1y = -(n+1)*Wn1y./r+y2*Wny;
                    d2Znx = (n*(n-1)./r2-x2).*Znx+Zn1x./r;
                    d2Zny = (n*(n-1)./r2-y2).*Zny+Zn1y./r;
                    d2Wnx = (n*(n-1)./r2-x2).*Wnx+Wn1x./r;
                    d2Wny = (n*(n-1)./r2-y2).*Wny+Wn1y./r;
                elseif PhaseVelocity(q) <= Material.LongitudinalVelocity && PhaseVelocity(q) > Material.TransverseVelocity
                    Znx = besseli(n,xr);
                    Zny = besselj(n,yr);
                    Wnx = besselk(n,xr);
                    Wny = bessely(n,yr);
                    Zn1x = x*besseli(n+1,xr);
                    Zn1y = y*besselj(n+1,yr);
                    Wn1x = x*besselk(n+1,xr);
                    Wn1y = y*bessely(n+1,yr);
                    dZnx = n*Znx./r+Zn1x;
                    dZny = n*Zny./r-Zn1y;
                    dWnx = n*Wnx./r-Wn1x;
                    dWny = n*Wny./r-Wn1y;
                    dg1Zn1y = -(n+1)*Zn1y./r+y2*Zny;
                    dg1Wn1y = -(n+1)*Wn1y./r+y2*Wny;
                    d2Znx = (n*(n-1)./r2-x2).*Znx-Zn1x./r;
                    d2Zny = (n*(n-1)./r2-y2).*Zny+Zn1y./r;
                    d2Wnx = (n*(n-1)./r2-x2).*Wnx+Wn1x./r;
                    d2Wny = (n*(n-1)./r2-y2).*Wny+Wn1y./r;
                elseif PhaseVelocity(q) <= Material.TransverseVelocity
                    Znx = besseli(n,xr);
                    Zny = besseli(n,yr);
                    Wnx = besselk(n,xr);
                    Wny = besselk(n,yr);
                    Zn1x = x*besseli(n+1,xr);
                    Zn1y = y*besseli(n+1,yr);
                    Wn1x = x*besselk(n+1,xr);
                    Wn1y = y*besselk(n+1,yr);
                    dZnx = n*Znx./r+Zn1x;
                    dZny = n*Zny./r+Zn1y;
                    dWnx = n*Wnx./r-Wn1x;
                    dWny = n*Wny./r-Wn1y;
                    dg1Zn1y = -(n+1)*Zn1y./r-y2*Zny;
                    dg1Wn1y = -(n+1)*Wn1y./r+y2*Wny;
                    d2Znx = (n*(n-1)./r2-x2).*Znx-Zn1x./r;
                    d2Zny = (n*(n-1)./r2-y2).*Zny-Zn1y./r;
                    d2Wnx = (n*(n-1)./r2-x2).*Wnx+Wn1x./r;
                    d2Wny = (n*(n-1)./r2-y2).*Wny+Wn1y./r;
                end
                if  Sink
                    Znzi = besselh(n,2,zri);
                    Zn1zi = zi*besselh(n+1,2,zri);
                else
                    Znzi = besselj(n,zri);
                    Zn1zi = zi*besselj(n+1,zri);
                end
                dZnzi = n*Znzi./rInnerFluid-Zn1zi;
                Hnzo = besselh(n,zRo);
                dHnzo = n*Hnzo/Ro-zo*besselh(n+1,zRo);
                Z1 = [Mu*2i*n*(-Wnx(1)/Ri2+dWnx(1)/Ri) Mu*1i*k*((n+1)*Zn1y(1)/Ri-dg1Zn1y(1)) Mu*1i*k*((n+1)*Wn1y(1)/Ri-dg1Wn1y(1)) Mu*1i*(2*d2Zny(1)+y2*Zny(1)) Mu*1i*(2*d2Wny(1)+y2*Wny(1)) 0 0;Mu*2i*k*dWnx(1) Mu*1i*((n*(n+1)/Ri2+k2-y2)*Zn1y(1)+n*dg1Zn1y(1)/Ri) Mu*1i*((n*(n+1)/Ri2+k2-y2)*Wn1y(1)+n*dg1Wn1y(1)/Ri) Mu*1i*n*k*Zny(1)/Ri Mu*1i*n*k*Wny(1)/Ri 0 0;dWnx(1) k*Zn1y(1) k*Wn1y(1) n*Zny(1)/Ri n*Wny(1)/Ri -dZnzi(end) 0;-Lambda*kL2*Wnx(end)+2*Mu*d2Wnx(end) 2*Mu*k*dg1Zn1y(end) 2*Mu*k*dg1Wn1y(end) 2*Mu*n*(-Zny(end)/Ro2+dZny(end)/Ro) 2*Mu*n*(-Wny(end)/Ro2+dWny(end)/Ro) 0 OuterFluid.Density*AngularFrequency2*Hnzo;Mu*2i*n*(-Wnx(end)/Ro2+dWnx(end)/Ro) Mu*1i*k*((n+1)*Zn1y(end)/Ro-dg1Zn1y(end)) Mu*1i*k*((n+1)*Wn1y(end)/Ro-dg1Wn1y(end)) Mu*1i*(2*d2Zny(end)+y2*Zny(end)) Mu*1i*(2*d2Wny(end)+y2*Wny(end)) 0 0;Mu*2i*k*dWnx(end) Mu*1i*((n*(n+1)/Ro2+k2-y2)*Zn1y(end)+n*dg1Zn1y(end)/Ro) Mu*1i*((n*(n+1)/Ro2+k2-y2)*Wn1y(end)+n*dg1Wn1y(end)/Ro) Mu*1i*n*k*Zny(end)/Ro Mu*1i*n*k*Wny(end)/Ro 0 0;dWnx(end) k*Zn1y(end) k*Wn1y(end) n*Zny(end)/Ro n*Wny(end)/Ro 0 -dHnzo];
                Z2 = [Mu*2i*n*(-Znx(1)/Ri2+dZnx(1)/Ri);Mu*2i*k*dZnx(1);dZnx(1);-Lambda*kL2*Znx(end)+2*Mu*d2Znx(end);Mu*2i*n*(-Znx(end)/Ro2+dZnx(end)/Ro);Mu*2i*k*dZnx(end);dZnx(end)];
                U = Z1\-Z2;
                f = Znx+U(1)*Wnx; % L_in + L_out
                df = dZnx+U(1)*dWnx;
                d2f = d2Znx+U(1)*d2Wnx;
                g1 = U(2)*Zn1y+U(3)*Wn1y; % SV_in + SV_out
                dg1 = U(2)*dg1Zn1y+U(3)*dg1Wn1y;
                g3 = U(4)*Zny+U(5)*Wny; % SH_in + SH_out
                dg3 = U(4)*dZny+U(5)*dWny;
                d2g3 = U(4)*d2Zny+U(5)*d2Wny;
                Delta = -Lambda*kL2*f;
                u1 = 1i*(k*f+(n+1)*g1./r+dg1);
                u2 = 1i*(n*f./r-k*g1+dg3);
                u3 = df+k*g1+n*g3./r;
                v1 = -1i*AngularFrequency(q)*u1;
                v2 = -1i*AngularFrequency(q)*u2;
                v3 = -1i*AngularFrequency(q)*u3;
                epsilon1 = 1i*k*u1;
                epsilon2 = (1i*n*u2+u3)./r;
                epsilon3 = d2f+k*dg1-n*g3./r2+n*dg3./r;
                epsilon4 = 1i*(-2*n*f./r2+2*n*df./r+(n+1)*k*g1./r-k*dg1+2*d2g3+y2*g3);
                epsilon5 = 1i*(2*k*df+(n*(n+1)./r2+k2-y2).*g1+n*dg1./r+n*k*g3./r);
                epsilon6 = 1i*(n*u1./r+k*u2);
                sigma1 = Delta+2*Mu*epsilon1;
                sigma2 = Delta+2*Mu*epsilon2;
                sigma3 = Delta+2*Mu*epsilon3;
                sigma4 = Mu*epsilon4;
                sigma5 = Mu*epsilon5;
                sigma6 = Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4)+epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2);
                PowerFlowDensity = -.5*real(sigma1.*conj(v1)+sigma6.*conj(v2)+sigma5.*conj(v3));
                f = U(6)*Znzi;
                df = U(6)*dZnzi;
                d2f = U(6)*((n*(n-1)./rInnerFluid2-zi^2).*Znzi+Zn1zi./rInnerFluid);
                u1 = 1i*k*f;
                u2 = 1i*n*f./rInnerFluid;
                u3 = df;
                v1 = -1i*AngularFrequency(q)*u1;
                v2 = -1i*AngularFrequency(q)*u2;
                v3 = -1i*AngularFrequency(q)*u3;
                epsilon1 = 1i*k*u1;
                epsilon2 = (1i*n*u2+u3)./rInnerFluid;
                epsilon3 = d2f;
                sigma = -InnerFluid.Density*AngularFrequency2*f;
                StrainEnergyDensityInnerFluid = .5*real(epsilon1.*conj(sigma)+epsilon2.*conj(sigma)+epsilon3.*conj(sigma));
                KineticEnergyDensityInnerFluid = .5*InnerFluid.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2);
                PowerFlowDensityInnerFluid = -.5*real(sigma.*conj(v1));
                StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
                KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
                PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
                PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
                TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
                EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % cez (m/ms)    
            end
            X{p}(:,5) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                PhaseVelocity = X{p}(q,4)*1e3;
                Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
                k2 = k^2;
                kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
                kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
                if  PhaseVelocity > Material.LongitudinalVelocity
                    x = sqrt(kL2-k2);
                    y = sqrt(kT2-k2);
                elseif PhaseVelocity <= Material.LongitudinalVelocity && PhaseVelocity > Material.TransverseVelocity
                    x = sqrt(k2-kL2);
                    y = sqrt(kT2-k2);
                elseif PhaseVelocity <= Material.TransverseVelocity
                    x = sqrt(k2-kL2);
                    y = sqrt(k2-kT2);
                end
                zi = sqrt((AngularFrequency/InnerFluid.Velocity)^2-k2);
                zo = sqrt((AngularFrequency/OuterFluid.Velocity)^2-k2);
                if  strcmp(ModeType,'FScholte') && Attenuation ~= 0
                    if  PhaseVelocity < InnerFluid.Velocity
                        zi = -zi;
                    end
                    if  PhaseVelocity < OuterFluid.Velocity
                        zo = -zo;
                    end
                end
                x2 = kL2-k2;
                y2 = kT2-k2;
                xr = x*r;
                yr = y*r;
                zri = zi*rInnerFluid;
                zRo = zo*Ro;
                if  PhaseVelocity > Material.LongitudinalVelocity
                    Znx = besselj(n,xr);
                    Zny = besselj(n,yr);
                    Wnx = bessely(n,xr);
                    Wny = bessely(n,yr);
                    Zn1x = x*besselj(n+1,xr);
                    Zn1y = y*besselj(n+1,yr);
                    Wn1x = x*bessely(n+1,xr);
                    Wn1y = y*bessely(n+1,yr);
                    dZnx = n*Znx./r-Zn1x;
                    dZny = n*Zny./r-Zn1y;
                    dWnx = n*Wnx./r-Wn1x;
                    dWny = n*Wny./r-Wn1y;
                    dg1Zn1y = -(n+1)*Zn1y./r+y2*Zny;
                    dg1Wn1y = -(n+1)*Wn1y./r+y2*Wny;
                    d2Znx = (n*(n-1)./r2-x2).*Znx+Zn1x./r;
                    d2Zny = (n*(n-1)./r2-y2).*Zny+Zn1y./r;
                    d2Wnx = (n*(n-1)./r2-x2).*Wnx+Wn1x./r;
                    d2Wny = (n*(n-1)./r2-y2).*Wny+Wn1y./r;
                elseif PhaseVelocity <= Material.LongitudinalVelocity && PhaseVelocity > Material.TransverseVelocity
                    Znx = besseli(n,xr);
                    Zny = besselj(n,yr);
                    Wnx = besselk(n,xr);
                    Wny = bessely(n,yr);
                    Zn1x = x*besseli(n+1,xr);
                    Zn1y = y*besselj(n+1,yr);
                    Wn1x = x*besselk(n+1,xr);
                    Wn1y = y*bessely(n+1,yr);
                    dZnx = n*Znx./r+Zn1x;
                    dZny = n*Zny./r-Zn1y;
                    dWnx = n*Wnx./r-Wn1x;
                    dWny = n*Wny./r-Wn1y;
                    dg1Zn1y = -(n+1)*Zn1y./r+y2*Zny;
                    dg1Wn1y = -(n+1)*Wn1y./r+y2*Wny;
                    d2Znx = (n*(n-1)./r2-x2).*Znx-Zn1x./r;
                    d2Zny = (n*(n-1)./r2-y2).*Zny+Zn1y./r;
                    d2Wnx = (n*(n-1)./r2-x2).*Wnx+Wn1x./r;
                    d2Wny = (n*(n-1)./r2-y2).*Wny+Wn1y./r;
                elseif PhaseVelocity <= Material.TransverseVelocity
                    Znx = besseli(n,xr);
                    Zny = besseli(n,yr);
                    Wnx = besselk(n,xr);
                    Wny = besselk(n,yr);
                    Zn1x = x*besseli(n+1,xr);
                    Zn1y = y*besseli(n+1,yr);
                    Wn1x = x*besselk(n+1,xr);
                    Wn1y = y*besselk(n+1,yr);
                    dZnx = n*Znx./r+Zn1x;
                    dZny = n*Zny./r+Zn1y;
                    dWnx = n*Wnx./r-Wn1x;
                    dWny = n*Wny./r-Wn1y;
                    dg1Zn1y = -(n+1)*Zn1y./r-y2*Zny;
                    dg1Wn1y = -(n+1)*Wn1y./r+y2*Wny;
                    d2Znx = (n*(n-1)./r2-x2).*Znx-Zn1x./r;
                    d2Zny = (n*(n-1)./r2-y2).*Zny-Zn1y./r;
                    d2Wnx = (n*(n-1)./r2-x2).*Wnx+Wn1x./r;
                    d2Wny = (n*(n-1)./r2-y2).*Wny+Wn1y./r;
                end
                if  Sink
                    Znzi = besselh(n,2,zri);
                    Zn1zi = zi*besselh(n+1,2,zri);
                else
                    Znzi = besselj(n,zri);
                    Zn1zi = zi*besselj(n+1,zri);
                end
                dZnzi = n*Znzi./rInnerFluid-Zn1zi;
                Hnzo = besselh(n,zRo);
                dHnzo = n*Hnzo/Ro-zo*besselh(n+1,zRo);
                Z1(1,1) = Mu*2i*n*(-Wnx(1)/Ri2+dWnx(1)/Ri);
                Z1(1,2) = Mu*1i*k*((n+1)*Zn1y(1)/Ri-dg1Zn1y(1));
                Z1(1,3) = Mu*1i*k*((n+1)*Wn1y(1)/Ri-dg1Wn1y(1));
                Z1(1,4) = Mu*1i*(2*d2Zny(1)+y2*Zny(1));
                Z1(1,5) = Mu*1i*(2*d2Wny(1)+y2*Wny(1));
                Z1(2,1) = Mu*2i*k*dWnx(1);
                Z1(2,2) = Mu*1i*((n*(n+1)/Ri2+k2-y2)*Zn1y(1)+n*dg1Zn1y(1)/Ri);
                Z1(2,3) = Mu*1i*((n*(n+1)/Ri2+k2-y2)*Wn1y(1)+n*dg1Wn1y(1)/Ri);
                Z1(2,4) = Mu*1i*n*k*Zny(1)/Ri;
                Z1(2,5) = Mu*1i*n*k*Wny(1)/Ri;
                Z1(3,1) = dWnx(1);
                Z1(3,2) = k*Zn1y(1);
                Z1(3,3) = k*Wn1y(1);
                Z1(3,4) = n*Zny(1)/Ri;
                Z1(3,5) = n*Wny(1)/Ri;
                Z1(3,6) = -dZnzi(end);
                Z1(4,1) = -Lambda*kL2*Wnx(end)+2*Mu*d2Wnx(end);
                Z1(4,2) = 2*Mu*k*dg1Zn1y(end);
                Z1(4,3) = 2*Mu*k*dg1Wn1y(end);
                Z1(4,4) = 2*Mu*n*(-Zny(end)/Ro2+dZny(end)/Ro);
                Z1(4,5) = 2*Mu*n*(-Wny(end)/Ro2+dWny(end)/Ro);
                Z1(4,7) = OuterFluid.Density*AngularFrequency2*Hnzo;
                Z1(5,1) = Mu*2i*n*(-Wnx(end)/Ro2+dWnx(end)/Ro);
                Z1(5,2) = Mu*1i*k*((n+1)*Zn1y(end)/Ro-dg1Zn1y(end));
                Z1(5,3) = Mu*1i*k*((n+1)*Wn1y(end)/Ro-dg1Wn1y(end));
                Z1(5,4) = Mu*1i*(2*d2Zny(end)+y2*Zny(end));
                Z1(5,5) = Mu*1i*(2*d2Wny(end)+y2*Wny(end));
                Z1(6,1) = Mu*2i*k*dWnx(end);
                Z1(6,2) = Mu*1i*((n*(n+1)/Ro2+k2-y2)*Zn1y(end)+n*dg1Zn1y(end)/Ro);
                Z1(6,3) = Mu*1i*((n*(n+1)/Ro2+k2-y2)*Wn1y(end)+n*dg1Wn1y(end)/Ro);
                Z1(6,4) = Mu*1i*n*k*Zny(end)/Ro;
                Z1(6,5) = Mu*1i*n*k*Wny(end)/Ro;
                Z1(7,1) = dWnx(end);
                Z1(7,2) = k*Zn1y(end);
                Z1(7,3) = k*Wn1y(end);
                Z1(7,4) = n*Zny(end)/Ro;
                Z1(7,5) = n*Wny(end)/Ro;
                Z1(7,7) = -dHnzo;
                Z2(1,1) = Mu*2i*n*(-Znx(1)/Ri2+dZnx(1)/Ri);
                Z2(2,1) = Mu*2i*k*dZnx(1);
                Z2(3,1) = dZnx(1);
                Z2(4,1) = -Lambda*kL2*Znx(end)+2*Mu*d2Znx(end);
                Z2(5,1) = Mu*2i*n*(-Znx(end)/Ro2+dZnx(end)/Ro);
                Z2(6,1) = Mu*2i*k*dZnx(end);
                Z2(7,1) = dZnx(end);
                U = Z1\-Z2;
                f = Znx+U(1)*Wnx; % L_in + L_out
                df = dZnx+U(1)*dWnx;
                d2f = d2Znx+U(1)*d2Wnx;
                g1 = U(2)*Zn1y+U(3)*Wn1y; % SV_in + SV_out
                dg1 = U(2)*dg1Zn1y+U(3)*dg1Wn1y;
                g3 = U(4)*Zny+U(5)*Wny; % SH_in + SH_out
                dg3 = U(4)*dZny+U(5)*dWny;
                d2g3 = U(4)*d2Zny+U(5)*d2Wny;
                Delta = -Lambda*kL2*f;
                u1 = 1i*(k*f+(n+1)*g1./r+dg1);
                u2 = 1i*(n*f./r-k*g1+dg3);
                u3 = df+k*g1+n*g3./r;
                v1 = -1i*AngularFrequency*u1;
                v2 = -1i*AngularFrequency*u2;
                v3 = -1i*AngularFrequency*u3;
                epsilon1 = 1i*k*u1;
                epsilon2 = (1i*n*u2+u3)./r;
                epsilon3 = d2f+k*dg1-n*g3./r2+n*dg3./r;
                epsilon4 = 1i*(-2*n*f./r2+2*n*df./r+(n+1)*k*g1./r-k*dg1+2*d2g3+y2*g3);
                epsilon5 = 1i*(2*k*df+(n*(n+1)./r2+k2-y2).*g1+n*dg1./r+n*k*g3./r);
                epsilon6 = 1i*(n*u1./r+k*u2);
                sigma1 = Delta+2*Mu*epsilon1;
                sigma2 = Delta+2*Mu*epsilon2;
                sigma3 = Delta+2*Mu*epsilon3;
                sigma4 = Mu*epsilon4;
                sigma5 = Mu*epsilon5;
                sigma6 = Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4)+epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2);
                PowerFlowDensity = -.5*real(sigma1.*conj(v1)+sigma6.*conj(v2)+sigma5.*conj(v3));
                f = U(6)*Znzi;
                df = U(6)*dZnzi;
                d2f = U(6)*((n*(n-1)./rInnerFluid2-zi^2).*Znzi+Zn1zi./rInnerFluid);
                u1 = 1i*k*f;
                u2 = 1i*n*f./rInnerFluid;
                u3 = df;
                v1 = -1i*AngularFrequency*u1;
                v2 = -1i*AngularFrequency*u2;
                v3 = -1i*AngularFrequency*u3;
                epsilon1 = 1i*k*u1;
                epsilon2 = (1i*n*u2+u3)./rInnerFluid;
                epsilon3 = d2f;
                sigma = -InnerFluid.Density*AngularFrequency2*f;
                StrainEnergyDensityInnerFluid = .5*real(epsilon1.*conj(sigma)+epsilon2.*conj(sigma)+epsilon3.*conj(sigma));
                KineticEnergyDensityInnerFluid = .5*InnerFluid.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2);
                PowerFlowDensityInnerFluid = -.5*real(sigma.*conj(v1));
                StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
                KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
                PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
                PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
                TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % cez (m/ms)
            end
        end
    elseif contains(ModeType,'L')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            Attenuation = X{p}(:,6).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                AngularFrequency2 = AngularFrequency(q)^2;
                k = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                k2 = k^2;
                kL2 = (AngularFrequency(q)/Material.LongitudinalVelocity_complex)^2;
                kT2 = (AngularFrequency(q)/Material.TransverseVelocity_complex)^2;
                if  PhaseVelocity(q) > Material.LongitudinalVelocity
                    x = sqrt(kL2-k2);
                    y = sqrt(kT2-k2);
                elseif PhaseVelocity(q) <= Material.LongitudinalVelocity && PhaseVelocity(q) > Material.TransverseVelocity
                    x = sqrt(k2-kL2);
                    y = sqrt(kT2-k2);
                elseif PhaseVelocity(q) <= Material.TransverseVelocity
                    x = sqrt(k2-kL2);
                    y = sqrt(k2-kT2);
                end
                zi = sqrt((AngularFrequency(q)/InnerFluid.Velocity)^2-k2);
                zo = sqrt((AngularFrequency(q)/OuterFluid.Velocity)^2-k2);
                if  strcmp(ModeType,'LScholte') && Attenuation(q) ~= 0
                    if  PhaseVelocity(q) < InnerFluid.Velocity
                        zi = -zi;
                    end
                    if  PhaseVelocity(q) < OuterFluid.Velocity
                        zo = -zo;
                    end
                end
                x2 = kL2-k2;
                y2 = kT2-k2;
                xr = x*r;
                yr = y*r;
                zri = zi*rInnerFluid;
                zRo = zo*Ro;
                if  PhaseVelocity(q) > Material.LongitudinalVelocity
                    Z0x = besselj(0,xr);
                    Z0y = besselj(0,yr);
                    W0x = bessely(0,xr);
                    W0y = bessely(0,yr);
                    Z1x = x*besselj(1,xr);
                    Z1y = y*besselj(1,yr);
                    W1x = x*bessely(1,xr);
                    W1y = y*bessely(1,yr);
                elseif PhaseVelocity(q) <= Material.LongitudinalVelocity && PhaseVelocity(q) > Material.TransverseVelocity
                    Z0x = besseli(0,xr);
                    Z0y = besselj(0,yr);
                    W0x = besselk(0,xr);
                    W0y = bessely(0,yr);
                    Z1x = -x*besseli(1,xr);
                    Z1y = y*besselj(1,yr);
                    W1x = x*besselk(1,xr);
                    W1y = y*bessely(1,yr);
                elseif PhaseVelocity(q) <= Material.TransverseVelocity
                    Z0x = besseli(0,xr);
                    Z0y = besseli(0,yr);
                    W0x = besselk(0,xr);
                    W0y = besselk(0,yr);
                    Z1x = -x*besseli(1,xr);
                    Z1y = -y*besseli(1,yr);
                    W1x = x*besselk(1,xr);
                    W1y = y*besselk(1,yr);
                end
                if  Sink
                    Z0zi = besselh(0,2,zri);
                    Z1zi = zi*besselh(1,2,zri);
                else
                    Z0zi = besselj(0,zri);
                    Z1zi = zi*besselj(1,zri);
                end
                H0zo = besselh(0,zRo);
                H1zo = zo*besselh(1,zRo);
                Z1 = [-Mu*2i*k*W1x(1) Mu*1i*(k2-y2)*Z1y(1) Mu*1i*(k2-y2)*W1y(1) 0 0;-W1x(1) k*Z1y(1) k*W1y(1) Z1zi(end) 0;-Lambda*kL2*W0x(end)+2*Mu*(W1x(end)/Ro-x2*W0x(end)) 2*Mu*k*(-Z1y(end)/Ro+y2*Z0y(end)) 2*Mu*k*(-W1y(end)/Ro+y2*W0y(end)) 0 OuterFluid.Density*AngularFrequency2*H0zo;-Mu*2i*k*W1x(end) Mu*1i*(k2-y2)*Z1y(end) Mu*1i*(k2-y2)*W1y(end) 0 0;-W1x(end) k*Z1y(end) k*W1y(end) 0 H1zo];
                Z2 = [Mu*2i*k*Z1x(1);Z1x(1);Lambda*kL2*Z0x(end)-2*Mu*(Z1x(end)/Ro-x2*Z0x(end));Mu*2i*k*Z1x(end);Z1x(end)];
                U = Z1\Z2;
                f = Z0x+U(1)*W0x; % L_in + L_out
                df = -Z1x-U(1)*W1x;
                d2f = Z1x./r-x2*Z0x+U(1)*(W1x./r-x2*W0x);
                g = U(2)*Z1y+U(3)*W1y; % SV_in + SV_out
                dg = U(2)*(-Z1y./r+y2*Z0y)+U(3)*(-W1y./r+y2*W0y);
                Delta = -Lambda*kL2*f;
                u1 = 1i*(k*f+g./r+dg);
                u3 = df+k*g;
                v1 = -1i*AngularFrequency(q)*u1;
                v3 = -1i*AngularFrequency(q)*u3;
                epsilon1 = 1i*k*u1;
                epsilon2 = u3./r;
                epsilon3 = d2f+k*dg;
                epsilon5 = 1i*(2*k*df+(k2-y2)*g);
                sigma1 = Delta+2*Mu*epsilon1;
                sigma2 = Delta+2*Mu*epsilon2;
                sigma3 = Delta+2*Mu*epsilon3;
                sigma5 = Mu*epsilon5;
                StrainEnergyDensity = .5*real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon5.*conj(sigma5));
                KineticEnergyDensity = .5*Material.Density*(abs(v1).^2+abs(v3).^2);
                PowerFlowDensity = -.5*real(sigma1.*conj(v1)+sigma5.*conj(v3));
                f = U(4)*Z0zi;
                df = -U(4)*Z1zi;
                d2f = U(4)*(Z1zi./rInnerFluid-zi^2*Z0zi);
                u1 = 1i*k*f;
                u3 = df;
                v1 = -1i*AngularFrequency(q)*u1;
                v3 = -1i*AngularFrequency(q)*u3;
                epsilon1 = 1i*k*u1;
                epsilon2 = u3./rInnerFluid;
                epsilon3 = d2f;
                sigma = -InnerFluid.Density*AngularFrequency2*f;
                StrainEnergyDensityInnerFluid = .5*real(epsilon1.*conj(sigma)+epsilon2.*conj(sigma)+epsilon3.*conj(sigma));
                KineticEnergyDensityInnerFluid = .5*InnerFluid.Density*(abs(v1).^2+abs(v3).^2);
                PowerFlowDensityInnerFluid = -.5*real(sigma.*conj(v1));
                StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
                KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
                PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
                PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
                TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
                EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % cez (m/ms)
            end
            X{p}(:,5) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                PhaseVelocity = X{p}(q,4)*1e3;
                Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                AngularFrequency2 = AngularFrequency^2;
                k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
                k2 = k^2;
                kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
                kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
                if  PhaseVelocity > Material.LongitudinalVelocity
                    x = sqrt(kL2-k2);
                    y = sqrt(kT2-k2);
                elseif PhaseVelocity <= Material.LongitudinalVelocity && PhaseVelocity > Material.TransverseVelocity
                    x = sqrt(k2-kL2);
                    y = sqrt(kT2-k2);
                elseif PhaseVelocity <= Material.TransverseVelocity
                    x = sqrt(k2-kL2);
                    y = sqrt(k2-kT2);
                end
                zi = sqrt((AngularFrequency/InnerFluid.Velocity)^2-k2);
                zo = sqrt((AngularFrequency/OuterFluid.Velocity)^2-k2);
                if  strcmp(ModeType,'LScholte') && Attenuation ~= 0
                    if  PhaseVelocity < InnerFluid.Velocity
                        zi = -zi;
                    end
                    if  PhaseVelocity < OuterFluid.Velocity
                        zo = -zo;
                    end
                end
                x2 = kL2-k2;
                y2 = kT2-k2;
                xr = x*r;
                yr = y*r;
                zri = zi*rInnerFluid;
                zRo = zo*Ro;
                if  PhaseVelocity > Material.LongitudinalVelocity
                    Z0x = besselj(0,xr);
                    Z0y = besselj(0,yr);
                    W0x = bessely(0,xr);
                    W0y = bessely(0,yr);
                    Z1x = x*besselj(1,xr);
                    Z1y = y*besselj(1,yr);
                    W1x = x*bessely(1,xr);
                    W1y = y*bessely(1,yr);
                elseif PhaseVelocity <= Material.LongitudinalVelocity && PhaseVelocity > Material.TransverseVelocity
                    Z0x = besseli(0,xr);
                    Z0y = besselj(0,yr);
                    W0x = besselk(0,xr);
                    W0y = bessely(0,yr);
                    Z1x = -x*besseli(1,xr);
                    Z1y = y*besselj(1,yr);
                    W1x = x*besselk(1,xr);
                    W1y = y*bessely(1,yr);
                elseif PhaseVelocity <= Material.TransverseVelocity
                    Z0x = besseli(0,xr);
                    Z0y = besseli(0,yr);
                    W0x = besselk(0,xr);
                    W0y = besselk(0,yr);
                    Z1x = -x*besseli(1,xr);
                    Z1y = -y*besseli(1,yr);
                    W1x = x*besselk(1,xr);
                    W1y = y*besselk(1,yr);
                end
                if  Sink
                    Z0zi = besselh(0,2,zri);
                    Z1zi = zi*besselh(1,2,zri);
                else
                    Z0zi = besselj(0,zri);
                    Z1zi = zi*besselj(1,zri);
                end
                H0zo = besselh(0,zRo);
                H1zo = zo*besselh(1,zRo);
                Z1(1,1) = -Mu*2i*k*W1x(1);
                Z1(1,2) = Mu*1i*(k2-y2)*Z1y(1);
                Z1(1,3) = Mu*1i*(k2-y2)*W1y(1);
                Z1(2,1) = -W1x(1);
                Z1(2,2) = k*Z1y(1);
                Z1(2,3) = k*W1y(1);
                Z1(2,4) = Z1zi(end);
                Z1(3,1) = -Lambda*kL2*W0x(end)+2*Mu*(W1x(end)/Ro-x2*W0x(end));
                Z1(3,2) = 2*Mu*k*(-Z1y(end)/Ro+y2*Z0y(end));
                Z1(3,3) = 2*Mu*k*(-W1y(end)/Ro+y2*W0y(end));
                Z1(3,5) = OuterFluid.Density*AngularFrequency2*H0zo;
                Z1(4,1) = -Mu*2i*k*W1x(end);
                Z1(4,2) = Mu*1i*(k2-y2)*Z1y(end);
                Z1(4,3) = Mu*1i*(k2-y2)*W1y(end);
                Z1(5,1) = -W1x(end);
                Z1(5,2) = k*Z1y(end);
                Z1(5,3) = k*W1y(end);
                Z1(5,5) = H1zo;
                Z2(1,1) = Mu*2i*k*Z1x(1);
                Z2(2,1) = Z1x(1);
                Z2(3,1) = Lambda*kL2*Z0x(end)-2*Mu*(Z1x(end)/Ro-x2*Z0x(end));
                Z2(4,1) = Mu*2i*k*Z1x(end);
                Z2(5,1) = Z1x(end);
                U = Z1\Z2;
                f = Z0x+U(1)*W0x; % L_in + L_out
                df = -Z1x-U(1)*W1x;
                d2f = Z1x./r-x2*Z0x+U(1)*(W1x./r-x2*W0x);
                g = U(2)*Z1y+U(3)*W1y; % SV_in + SV_out
                dg = U(2)*(-Z1y./r+y2*Z0y)+U(3)*(-W1y./r+y2*W0y);
                Delta = -Lambda*kL2*f;
                u1 = 1i*(k*f+g./r+dg);
                u3 = df+k*g;
                v1 = -1i*AngularFrequency*u1;
                v3 = -1i*AngularFrequency*u3;
                epsilon1 = 1i*k*u1;
                epsilon2 = u3./r;
                epsilon3 = d2f+k*dg;
                epsilon5 = 1i*(2*k*df+(k2-y2)*g);
                sigma1 = Delta+2*Mu*epsilon1;
                sigma2 = Delta+2*Mu*epsilon2;
                sigma3 = Delta+2*Mu*epsilon3;
                sigma5 = Mu*epsilon5;
                StrainEnergyDensity = .5*real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon5.*conj(sigma5));
                KineticEnergyDensity = .5*Material.Density*(abs(v1).^2+abs(v3).^2);
                PowerFlowDensity = -.5*real(sigma1.*conj(v1)+sigma5.*conj(v3));
                f = U(4)*Z0zi;
                df = -U(4)*Z1zi;
                d2f = U(4)*(Z1zi./rInnerFluid-zi^2*Z0zi);
                u1 = 1i*k*f;
                u3 = df;
                v1 = -1i*AngularFrequency*u1;
                v3 = -1i*AngularFrequency*u3;
                epsilon1 = 1i*k*u1;
                epsilon2 = u3./rInnerFluid;
                epsilon3 = d2f;
                sigma = -InnerFluid.Density*AngularFrequency2*f;
                StrainEnergyDensityInnerFluid = .5*real(epsilon1.*conj(sigma)+epsilon2.*conj(sigma)+epsilon3.*conj(sigma));
                KineticEnergyDensityInnerFluid = .5*InnerFluid.Density*(abs(v1).^2+abs(v3).^2);
                PowerFlowDensityInnerFluid = -.5*real(sigma.*conj(v1));
                StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
                KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
                PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
                PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
                TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % cez (m/ms)
            end
        end
    elseif strcmp(ModeType,'T')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            Attenuation = X{p}(:,6).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                k = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                k2 = k^2;
                kT2 = (AngularFrequency(q)/Material.TransverseVelocity_complex)^2;
                y = sqrt(kT2-k2);
                y2 = kT2-k2;
                yr = y*r;
                U = bessely(2,yr(end))\-besselj(2,yr(end));
                J0 = besselj(0,yr);
                Y0 = bessely(0,yr);
                J1 = y*besselj(1,yr);
                Y1 = y*bessely(1,yr);
                g = J0+U*Y0; % SH_in + SH_out
                dg = -J1-U*Y1;
                d2g = J1./r-y2*J0+U*(Y1./r-y2*Y0);
                u2 = 1i*dg;
                v2 = -1i*AngularFrequency(q)*u2;
                epsilon4 = 1i*(2*d2g+y2*g);
                epsilon6 = 1i*k*u2;
                sigma4 = Mu*epsilon4;
                sigma6 = Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon4.*conj(sigma4)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*abs(v2).^2;
                PowerFlowDensity = -.5*real(sigma6.*conj(v2));
                PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
                TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
                EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % cez (m/ms)
            end
            X{p}(:,5) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                PhaseVelocity = X{p}(q,4)*1e3;
                Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
                k2 = k^2;
                kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
                y = sqrt(kT2-k2);
                y2 = kT2-k2;
                yr = y*r;
                U = bessely(2,yr(end))\-besselj(2,yr(end));
                J0 = besselj(0,yr);
                Y0 = bessely(0,yr);
                J1 = y*besselj(1,yr);
                Y1 = y*bessely(1,yr);
                g = J0+U*Y0; % SH_in + SH_out
                dg = -J1-U*Y1;
                d2g = J1./r-y2*J0+U*(Y1./r-y2*Y0);
                u2 = 1i*dg;
                v2 = -1i*AngularFrequency*u2;
                epsilon4 = 1i*(2*d2g+y2*g);
                epsilon6 = 1i*k*u2;
                sigma4 = Mu*epsilon4;
                sigma6 = Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon4.*conj(sigma4)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*abs(v2).^2;
                PowerFlowDensity = -.5*real(sigma6.*conj(v2));
                PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
                TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % cez (m/ms)
            end
        end
    end
    X{p}(:,5) = filloutliers(X{p}(:,5),'spline','movmedian',5,'ThresholdFactor',1);
    Counter = Counter+1;
    if  ModeTotal > 0
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
    end
    if  Stop
        return
    end
end