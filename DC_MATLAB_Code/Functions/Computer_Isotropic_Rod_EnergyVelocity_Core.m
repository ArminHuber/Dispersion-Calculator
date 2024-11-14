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
function [X,Counter] = Computer_Isotropic_Rod_EnergyVelocity_Core(Multithreading,X,ModeType,n,ModeTotal,Counter,h,Material,Fluid,Viscoelastic,FluidLoading,R,SamplesR)
%#ok<*PFBNS>
%#ok<*GVMIS>
global Stop 
Stop = 0;
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
if  ~FluidLoading
    Fluid.Velocity = 1e-10;
    Fluid.Density = 1e-10;
end
R2 = R^2;
r = (0:R/SamplesR:R)';
r(1) = 1e-10;
r2 = r.^2;
Diff = diff(r2)';
if  strcmp(ModeType,'L') && ~Viscoelastic && ~FluidLoading
    cL2 = Material.LongitudinalVelocity^2;
    cT2 = Material.TransverseVelocity^2;
end
for p = 1:length(X)
    if  contains(ModeType,'F')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            Attenuation = X{p}(:,6).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                k = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                k2 = k^2;
                kL2 = (AngularFrequency(q)/Material.LongitudinalVelocity_complex)^2;
                kT2 = (AngularFrequency(q)/Material.TransverseVelocity_complex)^2;
                x = sqrt(kL2-k2);
                y = sqrt(kT2-k2);
                z = sqrt((AngularFrequency(q)/Fluid.Velocity)^2-k2);
                x2 = kL2-k2;
                y2 = kT2-k2;
                xr = x*r;
                yr = y*r;
                zR = z*R;
                Jnx = besselj(n,xr);
                Jny = besselj(n,yr);
                Jn1x = x*besselj(n+1,xr);
                Jn1y = y*besselj(n+1,yr);
                dJnx = n*Jnx./r-Jn1x;
                dJny = n*Jny./r-Jn1y;
                dg1Jn1y = -(n+1)*Jn1y./r+y2*Jny;
                d2Jnx = (n*(n-1)./r2-x2).*Jnx+Jn1x./r;
                d2Jny = (n*(n-1)./r2-y2).*Jny+Jn1y./r;
                dHnz = -n*besselh(n,zR)/R+z.*besselh(n+1,zR);
                Z1 = [Mu*1i*k*((n+1)*Jn1y(end)/R-dg1Jn1y(end)) Mu*1i*(2*d2Jny(end)+y2*Jny(end)) 0;Mu*1i*((n*(n+1)/R2+k2-y2)*Jn1y(end)+n*dg1Jn1y(end)/R) Mu*1i*n*k*Jny(end)/R 0;k*Jn1y(end) n*Jny(end)/R dHnz];
                Z2 = [Mu*2i*n*(-Jnx(end)/R2+dJnx(end)/R);Mu*2i*k*dJnx(end);dJnx(end)];
                U = Z1\-Z2;
                f = Jnx;
                df = dJnx;
                d2f = d2Jnx;
                g1 = U(1)*Jn1y;
                dg1 = U(1)*dg1Jn1y;
                g3 = U(2)*Jny;
                dg3 = U(2)*dJny;
                d2g3 = U(2)*d2Jny;
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
                epsilon6 = 1i*n*u1./r+1i*k*u2;
                sigma1 = Delta+2*Mu*epsilon1;
                sigma2 = Delta+2*Mu*epsilon2;
                sigma3 = Delta+2*Mu*epsilon3;
                sigma4 = Mu*epsilon4;
                sigma5 = Mu*epsilon5;
                sigma6 = Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4)+epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2);
                PowerFlowDensity = -.5*real(sigma1.*conj(v1)+sigma6.*conj(v2)+sigma5.*conj(v3));
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
                kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
                kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
                x = sqrt(kL2-k2);
                y = sqrt(kT2-k2);
                z = sqrt((AngularFrequency/Fluid.Velocity)^2-k2);
                x2 = kL2-k2;
                y2 = kT2-k2;
                xr = x*r;
                yr = y*r;
                zR = z*R;
                Jnx = besselj(n,xr);
                Jny = besselj(n,yr);
                Jn1x = x*besselj(n+1,xr);
                Jn1y = y*besselj(n+1,yr);
                dJnx = n*Jnx./r-Jn1x;
                dJny = n*Jny./r-Jn1y;
                dg1Jn1y = -(n+1)*Jn1y./r+y2*Jny;
                d2Jnx = (n*(n-1)./r2-x2).*Jnx+Jn1x./r;
                d2Jny = (n*(n-1)./r2-y2).*Jny+Jn1y./r;
                dHnz = -n*besselh(n,zR)/R+z.*besselh(n+1,zR);
                Z1(1,1) = Mu*1i*k*((n+1)*Jn1y(end)/R-dg1Jn1y(end));
                Z1(1,2) = Mu*1i*(2*d2Jny(end)+y2*Jny(end));
                Z1(2,1) = Mu*1i*((n*(n+1)/R2+k2-y2)*Jn1y(end)+n*dg1Jn1y(end)/R);
                Z1(2,2) = Mu*1i*n*k*Jny(end)/R;
                Z1(3,1) = k*Jn1y(end);
                Z1(3,2) = n*Jny(end)/R;
                Z1(3,3) = dHnz;
                Z2(1,1) = Mu*2i*n*(-Jnx(end)/R2+dJnx(end)/R);
                Z2(2,1) = Mu*2i*k*dJnx(end);
                Z2(3,1) = dJnx(end);
                U = Z1\-Z2;
                f = Jnx;
                df = dJnx;
                d2f = d2Jnx;
                g1 = U(1)*Jn1y;
                dg1 = U(1)*dg1Jn1y;
                g3 = U(2)*Jny;
                dg3 = U(2)*dJny;
                d2g3 = U(2)*d2Jny;
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
                epsilon6 = 1i*n*u1./r+1i*k*u2;
                sigma1 = Delta+2*Mu*epsilon1;
                sigma2 = Delta+2*Mu*epsilon2;
                sigma3 = Delta+2*Mu*epsilon3;
                sigma4 = Mu*epsilon4;
                sigma5 = Mu*epsilon5;
                sigma6 = Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4)+epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2);
                PowerFlowDensity = -.5*real(sigma1.*conj(v1)+sigma6.*conj(v2)+sigma5.*conj(v3));
                PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
                TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % cez (m/ms)
            end
        end
    elseif contains(ModeType,'L')
        if  strcmp(ModeType,'L') && ~Viscoelastic && ~FluidLoading
            PhaseVelocity = X{p}(:,4)*1e3;
            AngularFrequency2 = (2*pi*X{p}(:,1)*1e3).^2;
            k2 = AngularFrequency2./PhaseVelocity.^2;
            x = sqrt(AngularFrequency2/cL2-k2);
            y = sqrt(AngularFrequency2/cT2-k2);
            y2 = y.^2;
            xR = x*R;
            yR = y*R;
            J0x = besselj(0,xR);
            J0y = besselj(0,yR);
            J2x = besselj(2,xR);
            J2y = besselj(2,yR);
            J1x = xR/2.*(J0x+J2x);
            J1y = yR/2.*(J0y+J2y);
            a1 = J0x./J1x;
            a2 = J0y./J1y;
            a3 = 2./xR.*(y2+k2);
            a4 = R./x.*(y2-k2).^2.*(1+J0x.*(J0x-J2x)./J1x.^2/2);
            a5 = R.*x.*k2.*(1+J0y.*(J0y-J2y)./J1y.^2/2);
            A1 = a3+(8*(k2-y2).*a1+a4)+4*(2*x.*y.*a2-k2./x.*y.*a2-k2.*x./y.*a2+a5);
            A2 = a3/cL2+4*x/R/cT2-(4*(y2-k2).*a1/cT2-a4/cL2)-4*(k2./x.*y.*a2/cL2+k2.*x./y.*a2/cT2-a5/cT2);
            X{p}(:,5) = A1./A2./PhaseVelocity/1e3; % cez (m/ms)
        else
            if  Multithreading
                EnergyVelocity = 0;
                PhaseVelocity = X{p}(:,4)*1e3;
                Attenuation = X{p}(:,6).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
                AngularFrequency = 2*pi*X{p}(:,1)*1e3;
                parfor q = 1:height(X{p})
                    k = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                    k2 = k^2;
                    kL2 = (AngularFrequency(q)/Material.LongitudinalVelocity_complex)^2;
                    kT2 = (AngularFrequency(q)/Material.TransverseVelocity_complex)^2;
                    x = sqrt(kL2-k2);
                    y = sqrt(kT2-k2);
                    z = sqrt((AngularFrequency(q)/Fluid.Velocity)^2-k2);
                    x2 = kL2-k2;
                    y2 = kT2-k2;
                    xr = x*r;
                    yr = y*r;
                    zR = z*R;
                    J0x = besselj(0,xr);
                    J0y = besselj(0,yr);
                    J1x = x*besselj(1,xr);
                    J1y = y*besselj(1,yr);
                    H1z = z*besselh(1,zR);
                    Z1 = [Mu*1i*(k2-y2)*J1y(end) 0;k*J1y(end) H1z];
                    Z2 = [Mu*2i*k*J1x(end);J1x(end)];
                    U = Z1\Z2;
                    f = J0x;
                    df = -J1x;
                    d2f = J1x./r-x2*J0x;
                    g = U(1)*J1y;
                    dg = U(1)*(-J1y./r+y2*J0y);
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
                    kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
                    kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
                    x = sqrt(kL2-k2);
                    y = sqrt(kT2-k2);
                    z = sqrt((AngularFrequency/Fluid.Velocity)^2-k2);
                    x2 = kL2-k2;
                    y2 = kT2-k2;
                    xr = x*r;
                    yr = y*r;
                    zR = z*R;
                    J0x = besselj(0,xr);
                    J0y = besselj(0,yr);
                    J1x = x*besselj(1,xr);
                    J1y = y*besselj(1,yr);
                    H1z = z*besselh(1,zR);
                    Z1(1,1) = Mu*1i*(k2-y2)*J1y(end);
                    Z1(2,1) = k*J1y(end);
                    Z1(2,2) = H1z;
                    Z2(1,1) = Mu*2i*k*J1x(end);
                    Z2(2,1) = J1x(end);
                    U = Z1\Z2;
                    f = J0x;
                    df = -J1x;
                    d2f = J1x./r-x2*J0x;
                    g = U(1)*J1y;
                    dg = U(1)*(-J1y./r+y2*J0y);
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
                    PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
                    TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
                    X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % cez (m/ms)
                end
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
                J0 = besselj(0,yr);
                J1 = y*besselj(1,yr);
                g = J0;
                dg = -J1;
                d2g = J1./r-y2*J0;
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
                J0 = besselj(0,yr);
                J1 = y*besselj(1,yr);
                g = J0;
                dg = -J1;
                d2g = J1./r-y2*J0;
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