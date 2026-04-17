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
function [u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,r,n,p,Frequency] = ModeShapeLinesComputer_Isotropic_Rod_Pipe(Geometry,Material,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,F,L,T,Frequency,Mode,Ro,Ri,SamplesR,ShowHalfSpace,HalfSpaces,Phase)
uPhase=0;epsilonPhase=0;sigmaPhase=0;
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
p = str2double(regexp(Mode,'\d*','match'));
n = p(1);
if  contains(Mode,'F')
    [PhaseVelocity,Attenuation,SignFluids,Direction,Frequency] = ExtractData(Geometry,F{n}{p(2)},Frequency,Mode);
elseif contains(Mode,'L')
    [PhaseVelocity,Attenuation,SignFluids,Direction,Frequency] = ExtractData(Geometry,L{p(2)},Frequency,Mode);
elseif contains(Mode,'T')
    [PhaseVelocity,Attenuation,SignFluids,Direction,Frequency] = ExtractData(Geometry,T{p(2)},Frequency,Mode);
end
AngularFrequency = Frequency*pi*2e3;
AngularFrequency2 = AngularFrequency^2;
k = Direction*(AngularFrequency/PhaseVelocity+1i*Attenuation);
k2 = k^2;
y2 = AngularFrequency2/Material.TransverseVelocity_complex^2-k2;
if  strcmp(Geometry,'Pipe')
    if  ~ToggleOuterFluid
        OuterFluid.Velocity = 1e-10;
        OuterFluid.Density = 1e-10;
    end
    if  ~ToggleInnerFluid
        InnerFluid.Velocity = 1e-10;
        InnerFluid.Density = 1e-10;
    end
    r = (Ri:(Ro-Ri)/SamplesR:Ro)';
    rInnerFluid = (Ri/SamplesR:Ri/SamplesR:Ri)';
    rOuterFluid = (Ro:Ri/SamplesR:Ro+HalfSpaces*Ri)';
    r2 = r.^2;
    Ro2 = Ro^2;
    Ri2 = Ri^2;
    if  contains(Mode,'F')
        kL2 = AngularFrequency2/Material.LongitudinalVelocity_complex^2;
        x2 = kL2-k2;
        x = sqrt(-x2);
        y = sqrt(-y2);
        zo = SignFluids(1)*sqrt(AngularFrequency2/OuterFluid.Velocity^2-k2);
        zi = SignFluids(2)*sqrt(AngularFrequency2/InnerFluid.Velocity^2-k2);
        xr = x*r;
        yr = y*r;
        zri = zi*rInnerFluid;
        zro = zo*rOuterFluid;
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
        if  Sink
            Znzi = besselh(n,2,zri);
            Zn1zi = zi*besselh(n+1,2,zri);
        else
            Znzi = besselj(n,zri);
            Zn1zi = zi*besselj(n+1,zri);
        end
        dZnzi = n*Znzi./rInnerFluid-Zn1zi;
        Hnzo = besselh(n,zro);
        Hn1zo = zo*besselh(n+1,zro);
        dHnzo = n*Hnzo./rOuterFluid-Hn1zo;
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
        Z1(4,7) = OuterFluid.Density*AngularFrequency2*Hnzo(1);
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
        Z1(7,7) = -dHnzo(1);
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
        u(:,1) = 1i*(k*f+(n+1)*g1./r+dg1);
        u(:,2) = 1i*(n*f./r-k*g1+dg3);
        u(:,3) = df+k*g1+n*g3./r;
        v = -1i*AngularFrequency*u;
        epsilon(:,1) = 1i*k*u(:,1);
        epsilon(:,2) = (1i*n*u(:,2)+u(:,3))./r;
        epsilon(:,3) = d2f+k*dg1-n*g3./r2+n*dg3./r;
        epsilon(:,4) = 1i*(-2*n*f./r2+2*n*df./r+(n+1)*k*g1./r-k*dg1+2*d2g3+y2*g3);
        epsilon(:,5) = 1i*(2*k*df+(n*(n+1)./r2+k2-y2).*g1+n*dg1./r+n*k*g3./r);
        epsilon(:,6) = 1i*(n*u(:,1)./r+k*u(:,2));
        sigma(:,1:3) = -Lambda*kL2*f+2*Mu*epsilon(:,1:3);
        sigma(:,4:6) = Mu*epsilon(:,4:6);
        StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
        KineticEnergyDensity = Material.Density*sum(abs(v).^2,2)/2;
        PowerFlowDensity(:,1) = -real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)))/2;
        PowerFlowDensity(:,2) = -real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)))/2;
        PowerFlowDensity(:,3) = -real(sigma(:,5).*conj(v(:,1))+sigma(:,4).*conj(v(:,2))+sigma(:,3).*conj(v(:,3)))/2;
        if  Phase
            sigmaPhase456 = rad2deg(angle(sigma(:,4:6)*exp(-1i*angle(u(1,1)))));
            epsilonPhase456 = rad2deg(angle(epsilon(:,4:6)*exp(-1i*angle(u(1,1)))));
        end
        epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4))); % zero in the fluid, so we make phase shift now
        epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
        epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
        sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
        sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
        sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
        if  ToggleInnerFluid
            f = U(6)*Znzi;
            df = U(6)*dZnzi;
            d2f = U(6)*((n*(n-1)./rInnerFluid.^2-zi^2).*Znzi+Zn1zi./rInnerFluid);
            uInnerFluid(:,1) = 1i*k*f;
            uInnerFluid(:,2) = 1i*n*f./rInnerFluid;
            uInnerFluid(:,3) = df;
            vInnerFluid = -1i*AngularFrequency*uInnerFluid;
            epsilonInnerFluid(:,1) = 1i*k*uInnerFluid(:,1);
            epsilonInnerFluid(:,2) = (1i*n*uInnerFluid(:,2)+uInnerFluid(:,3))./rInnerFluid;
            epsilonInnerFluid(:,3) = d2f;
            epsilonInnerFluid(:,6) = 0;
            sigmaInnerFluid(:,1) = -InnerFluid.Density*AngularFrequency2*f;
            sigmaInnerFluid(:,2) = sigmaInnerFluid(:,1);
            sigmaInnerFluid(:,3) = sigmaInnerFluid(:,1);
            sigmaInnerFluid(:,6) = 0;
            StrainEnergyDensityInnerFluid = real(sum(epsilonInnerFluid.*conj(sigmaInnerFluid),2))/2;
            KineticEnergyDensityInnerFluid = InnerFluid.Density*sum(abs(vInnerFluid).^2,2)/2;
            PowerFlowDensityInnerFluid(:,1) = -real(sigmaInnerFluid(:,1).*conj(vInnerFluid(:,1)))/2;
            PowerFlowDensityInnerFluid(:,2) = -real(sigmaInnerFluid(:,2).*conj(vInnerFluid(:,2)))/2;
            PowerFlowDensityInnerFluid(:,3) = -real(sigmaInnerFluid(:,3).*conj(vInnerFluid(:,3)))/2;
            r = [rInnerFluid;r];
            u = [uInnerFluid;u];
            epsilon = [epsilonInnerFluid;epsilon];
            sigma = [sigmaInnerFluid;sigma];
            StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
            KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
            PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
        end
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(2*PowerFlow/pi/trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
        if  ToggleOuterFluid && ShowHalfSpace
            f = U(7)*Hnzo;
            df = U(7)*dHnzo;
            d2f = U(7)*((n*(n-1)./rOuterFluid.^2-zo^2).*Hnzo+Hn1zo./rOuterFluid);
            uOuterFluid(:,1) = 1i*k*f;
            uOuterFluid(:,2) = 1i*n*f./rOuterFluid;
            uOuterFluid(:,3) = df;
            vOuterFluid = -1i*AngularFrequency*uOuterFluid;
            epsilonOuterFluid(:,1) = 1i*k*uOuterFluid(:,1);
            epsilonOuterFluid(:,2) = (1i*n*uOuterFluid(:,2)+uOuterFluid(:,3))./rOuterFluid;
            epsilonOuterFluid(:,3) = d2f;
            epsilonOuterFluid(:,6) = 0;
            sigmaOuterFluid(:,1) = -OuterFluid.Density*AngularFrequency2*f;
            sigmaOuterFluid(:,2) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,3) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,6) = 0;
            StrainEnergyDensityOuterFluid = real(sum(epsilonOuterFluid.*conj(sigmaOuterFluid),2))/2;
            KineticEnergyDensityOuterFluid = OuterFluid.Density*sum(abs(vOuterFluid).^2,2)/2;
            PowerFlowDensityOuterFluid(:,1) = -real(sigmaOuterFluid(:,1).*conj(vOuterFluid(:,1)))/2;
            PowerFlowDensityOuterFluid(:,2) = -real(sigmaOuterFluid(:,2).*conj(vOuterFluid(:,2)))/2;
            PowerFlowDensityOuterFluid(:,3) = -real(sigmaOuterFluid(:,3).*conj(vOuterFluid(:,3)))/2;
            r = [r;rOuterFluid];
            u = [u;uOuterFluid];
            epsilon = [epsilon;epsilonOuterFluid];
            sigma = [sigma;sigmaOuterFluid];
            StrainEnergyDensity = [StrainEnergyDensity;StrainEnergyDensityOuterFluid];
            KineticEnergyDensity = [KineticEnergyDensity;KineticEnergyDensityOuterFluid];
            PowerFlowDensity = [PowerFlowDensity;PowerFlowDensityOuterFluid];
        end
        if  Phase
            if  FluidLoading
                if  ToggleInnerFluid
                    uPhase = rad2deg(angle(u*exp(-1i*angle(u(length(rInnerFluid)+1,1)))));
                    sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(length(rInnerFluid)+1,1)))));
                    epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(length(rInnerFluid)+1,1)))));
                    sigmaPhase(length(rInnerFluid)+1:length(rInnerFluid)+length(sigmaPhase456),4:6) = sigmaPhase456;
                    epsilonPhase(length(rInnerFluid)+1:length(rInnerFluid)+length(epsilonPhase456),4:6) = epsilonPhase456;
                else
                    uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
                    sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
                    epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
                    sigmaPhase(1:length(sigmaPhase456),4:6) = sigmaPhase456;
                    epsilonPhase(1:length(epsilonPhase456),4:6) = epsilonPhase456;
                end
            else
                uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
                sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
                epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
                sigmaPhase(:,4:6) = sigmaPhase456;
                epsilonPhase(:,4:6) = epsilonPhase456;
            end
        end
        u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
        u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
        u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
        epsilon(:,1) = epsilon(:,1)*exp(-1i*angle(epsilon(2,1)));
        epsilon(:,2) = epsilon(:,2)*exp(-1i*angle(epsilon(2,2)));
        epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
        sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
        sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
        sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
    elseif contains(Mode,'L')
        kL2 = AngularFrequency2/Material.LongitudinalVelocity_complex^2;
        x2 = kL2-k2;
        x = sqrt(-x2);
        y = sqrt(-y2);
        zo = SignFluids(1)*sqrt(AngularFrequency2/OuterFluid.Velocity^2-k2);
        zi = SignFluids(2)*sqrt(AngularFrequency2/InnerFluid.Velocity^2-k2);
        xr = x*r;
        yr = y*r;
        zri = zi*rInnerFluid;
        zro = zo*rOuterFluid;
        Z0x = besseli(0,xr);
        Z0y = besseli(0,yr);
        W0x = besselk(0,xr);
        W0y = besselk(0,yr);
        Z1x = -x*besseli(1,xr);
        Z1y = -y*besseli(1,yr);
        W1x = x*besselk(1,xr);
        W1y = y*besselk(1,yr);
        if  Sink
            Z0zi = besselh(0,2,zri);
            Z1zi = zi*besselh(1,2,zri);
        else
            Z0zi = besselj(0,zri);
            Z1zi = zi*besselj(1,zri);
        end
        H0zo = besselh(0,zro);
        H1zo = zo*besselh(1,zro);
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
        Z1(3,5) = OuterFluid.Density*AngularFrequency2*H0zo(1);
        Z1(4,1) = -Mu*2i*k*W1x(end);
        Z1(4,2) = Mu*1i*(k2-y2)*Z1y(end);
        Z1(4,3) = Mu*1i*(k2-y2)*W1y(end);
        Z1(5,1) = -W1x(end);
        Z1(5,2) = k*Z1y(end);
        Z1(5,3) = k*W1y(end);
        Z1(5,5) = H1zo(1);
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
        u(:,1) = 1i*(k*f+g./r+dg);
        u(:,3) = df+k*g;
        v = -1i*AngularFrequency*u;
        epsilon(:,1) = 1i*k*u(:,1);
        epsilon(:,2) = u(:,3)./r;
        epsilon(:,3) = d2f+k*dg;
        epsilon(:,5) = 1i*(2*k*df+(k2-y2)*g);
        epsilon(:,6) = 0;
        sigma(:,1:3) = -Lambda*kL2*f+2*Mu*epsilon(:,1:3);
        sigma(:,5) = Mu*epsilon(:,5);
        sigma(:,6) = 0;
        StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
        KineticEnergyDensity = Material.Density*sum(abs(v).^2,2)/2;
        PowerFlowDensity(:,1) = -real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)))/2;
        PowerFlowDensity(:,3) = -real(sigma(:,5).*conj(v(:,1))+sigma(:,3).*conj(v(:,3)))/2;
        if  Phase
            sigmaPhase5 = rad2deg(angle(sigma(:,5)*exp(-1i*angle(u(1,1)))));
            epsilonPhase5 = rad2deg(angle(epsilon(:,5)*exp(-1i*angle(u(1,1)))));
        end
        epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5))); % zero in the fluid, so we make phase shift now 
        sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
        if  ToggleInnerFluid
            f = U(4)*Z0zi;
            df = -U(4)*Z1zi;
            d2f = U(4)*(Z1zi./rInnerFluid-zi^2*Z0zi);
            uInnerFluid(:,1) = 1i*k*f;
            uInnerFluid(:,3) = df;
            vInnerFluid = -1i*AngularFrequency*uInnerFluid;
            epsilonInnerFluid(:,1) = 1i*k*uInnerFluid(:,1);
            epsilonInnerFluid(:,2) = uInnerFluid(:,3)./rInnerFluid;
            epsilonInnerFluid(:,3) = d2f;
            epsilonInnerFluid(:,6) = 0;
            sigmaInnerFluid(:,1) = -InnerFluid.Density*AngularFrequency2*f;
            sigmaInnerFluid(:,2) = sigmaInnerFluid(:,1);
            sigmaInnerFluid(:,3) = sigmaInnerFluid(:,1);
            sigmaInnerFluid(:,6) = 0;
            StrainEnergyDensityInnerFluid = real(sum(epsilonInnerFluid.*conj(sigmaInnerFluid),2))/2;
            KineticEnergyDensityInnerFluid = InnerFluid.Density*sum(abs(vInnerFluid).^2,2)/2;
            PowerFlowDensityInnerFluid(:,1) = -real(sigmaInnerFluid(:,1).*conj(vInnerFluid(:,1)))/2;
            PowerFlowDensityInnerFluid(:,3) = -real(sigmaInnerFluid(:,3).*conj(vInnerFluid(:,3)))/2;
            r = [rInnerFluid;r];
            u = [uInnerFluid;u];
            epsilon = [epsilonInnerFluid;epsilon];
            sigma = [sigmaInnerFluid;sigma];
            StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
            KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
            PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
        end
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(2*PowerFlow/pi/trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
        if  ToggleOuterFluid && ShowHalfSpace
            f = U(5)*H0zo;
            df = -U(5)*H1zo;
            d2f = U(5)*(H1zo./rOuterFluid-zo^2*H0zo);
            uOuterFluid(:,1) = 1i*k*f;
            uOuterFluid(:,3) = df;
            vOuterFluid = -1i*AngularFrequency*uOuterFluid;
            epsilonOuterFluid(:,1) = 1i*k*uOuterFluid(:,1);
            epsilonOuterFluid(:,2) = uOuterFluid(:,3)./rOuterFluid;
            epsilonOuterFluid(:,3) = d2f;
            epsilonOuterFluid(:,6) = 0;
            sigmaOuterFluid(:,1) = -OuterFluid.Density*AngularFrequency2*f;
            sigmaOuterFluid(:,2) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,3) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,6) = 0;
            StrainEnergyDensityOuterFluid = real(sum(epsilonOuterFluid.*conj(sigmaOuterFluid),2))/2;
            KineticEnergyDensityOuterFluid = OuterFluid.Density*sum(abs(vOuterFluid).^2,2)/2;
            PowerFlowDensityOuterFluid(:,1) = -real(sigmaOuterFluid(:,1).*conj(vOuterFluid(:,1)))/2;
            PowerFlowDensityOuterFluid(:,3) = -real(sigmaOuterFluid(:,3).*conj(vOuterFluid(:,3)))/2;
            r = [r;rOuterFluid];
            u = [u;uOuterFluid];
            epsilon = [epsilon;epsilonOuterFluid];
            sigma = [sigma;sigmaOuterFluid];
            StrainEnergyDensity = [StrainEnergyDensity;StrainEnergyDensityOuterFluid];
            KineticEnergyDensity = [KineticEnergyDensity;KineticEnergyDensityOuterFluid];
            PowerFlowDensity = [PowerFlowDensity;PowerFlowDensityOuterFluid];
        end
        if  Phase
            if  FluidLoading
                if  ToggleInnerFluid
                    uPhase = rad2deg(angle(u*exp(-1i*angle(u(length(rInnerFluid)+1,1)))));
                    sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(length(rInnerFluid)+1,1)))));
                    epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(length(rInnerFluid)+1,1)))));
                    sigmaPhase(length(rInnerFluid)+1:length(rInnerFluid)+length(sigmaPhase5),5) = sigmaPhase5;
                    epsilonPhase(length(rInnerFluid)+1:length(rInnerFluid)+length(epsilonPhase5),5) = epsilonPhase5;
                else
                    uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
                    sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
                    epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
                    sigmaPhase(1:length(sigmaPhase5),5) = sigmaPhase5;
                    epsilonPhase(1:length(epsilonPhase5),5) = epsilonPhase5;
                end
            else
                uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
                sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
                epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
                sigmaPhase(:,5) = sigmaPhase5;
                epsilonPhase(:,5) = epsilonPhase5;
            end
        end
        u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
        u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
        epsilon(:,1) = epsilon(:,1)*exp(-1i*angle(epsilon(2,1)));
        epsilon(:,2) = epsilon(:,2)*exp(-1i*angle(epsilon(2,2)));
        epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
        sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
        sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
        sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
    elseif contains(Mode,'T')
        y = sqrt(y2);
        if  y == 0
            y = 1e-10;
            y2 = y^2;
        end
        yr = y*r;
        U = bessely(2,yr(end))\-besselj(2,yr(end));
        J0 = besselj(0,yr);
        Y0 = bessely(0,yr);
        J1 = y*besselj(1,yr);
        Y1 = y*bessely(1,yr);
        g = J0+U*Y0; % SH_in + SH_out
        dg = -J1-U*Y1;
        d2g = J1./r-y2*J0+U*(Y1./r-y2*Y0);
        u(:,2) = 1i*dg;
        u(:,3) = 0;
        v = -1i*AngularFrequency*u;
        epsilon(:,4) = 1i*(2*d2g+y2*g);
        epsilon(:,6) = 1i*k*u(:,2);
        sigma(:,4:6) = Mu*epsilon(:,4:6);
        StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
        KineticEnergyDensity = Material.Density*abs(v(:,2)).^2/2;
        PowerFlowDensity(:,1) = -real(sigma(:,6).*conj(v(:,2)))/2;
        PowerFlowDensity(:,3) = -real(sigma(:,4).*conj(v(:,2)))/2;
        if  Phase
            uPhase = rad2deg(angle(u));
            sigmaPhase = rad2deg(angle(sigma));
            epsilonPhase = rad2deg(angle(epsilon));
        end
        u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
        epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
        epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
        sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
        sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
        if  ToggleInnerFluid
            uInnerFluid(length(rInnerFluid),3) = 0;
            epsilonInnerFluid(length(rInnerFluid),6) = 0;
            sigmaInnerFluid(length(rInnerFluid),6) = 0;
            StrainEnergyDensityInnerFluid(length(rInnerFluid),1) = 0;
            KineticEnergyDensityInnerFluid(length(rInnerFluid),1) = 0;
            PowerFlowDensityInnerFluid(length(rInnerFluid),3) = 0;
            r = [rInnerFluid;r];
            u = [uInnerFluid;u];
            epsilon = [epsilonInnerFluid;epsilon];
            sigma = [sigmaInnerFluid;sigma];
            StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
            KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
            PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
            if  Phase
                uPhase = [zeros(length(rInnerFluid),3);uPhase];
                sigmaPhase = [zeros(length(rInnerFluid),6);sigmaPhase];
                epsilonPhase = [zeros(length(rInnerFluid),6);epsilonPhase];
            end
        end
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(2*PowerFlow/pi/trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
        if  ToggleOuterFluid && ShowHalfSpace
            uOuterFluid(length(rOuterFluid),3) = 0;
            epsilonOuterFluid(length(rOuterFluid),6) = 0;
            sigmaOuterFluid(length(rOuterFluid),6) = 0;
            StrainEnergyDensityOuterFluid(length(rOuterFluid),1) = 0;
            KineticEnergyDensityOuterFluid(length(rOuterFluid),1) = 0;
            PowerFlowDensityOuterFluid(length(rOuterFluid),3) = 0;
            r = [r;rOuterFluid];
            u = [u;uOuterFluid];
            epsilon = [epsilon;epsilonOuterFluid];
            sigma = [sigma;sigmaOuterFluid];
            StrainEnergyDensity = [StrainEnergyDensity;StrainEnergyDensityOuterFluid];
            KineticEnergyDensity = [KineticEnergyDensity;KineticEnergyDensityOuterFluid];
            PowerFlowDensity = [PowerFlowDensity;PowerFlowDensityOuterFluid];
            if  Phase
                uPhase = [uPhase;zeros(length(rOuterFluid),3)];
                sigmaPhase = [sigmaPhase;zeros(length(rOuterFluid),6)];
                epsilonPhase = [epsilonPhase;zeros(length(rOuterFluid),6)];
            end
        end
    end
elseif strcmp(Geometry,'Rod')
    if  ~ToggleOuterFluid
        OuterFluid.Velocity = 1e-10;
        OuterFluid.Density = 1e-10;
    end
    R = Ro;
    r = (0:R/SamplesR:R)';
    r(1) = 1e-10;
    rOuterFluid = (R:R/SamplesR:R+HalfSpaces*R)';
    r2 = r.^2;
    R2 = R^2;
    if  contains(Mode,'F')
        kL2 = AngularFrequency2/Material.LongitudinalVelocity_complex^2;
        x2 = kL2-k2;
        x = sqrt(x2);
        y = sqrt(y2);
        z = SignFluids*sqrt(AngularFrequency2/OuterFluid.Velocity^2-k2);
        xr = x*r;
        yr = y*r;
        zr = z*rOuterFluid;
        Jnx = besselj(n,xr);
        Jny = besselj(n,yr);
        Jn1x = x*besselj(n+1,xr);
        Jn1y = y*besselj(n+1,yr);
        dJnx = n*Jnx./r-Jn1x;
        dJny = n*Jny./r-Jn1y;
        dg1Jn1y = -(n+1)*Jn1y./r+y2*Jny;
        d2Jnx = (n*(n-1)./r2-x2).*Jnx+Jn1x./r;
        d2Jny = (n*(n-1)./r2-y2).*Jny+Jn1y./r;
        Hnz = besselh(n,zr);
        Hn1z = z*besselh(n+1,zr);
        dHnz = n*Hnz./rOuterFluid-Hn1z;
        Z1(1,1) = Mu*1i*k*((n+1)*Jn1y(end)/R-dg1Jn1y(end));
        Z1(1,2) = Mu*1i*(2*d2Jny(end)+y2*Jny(end));
        Z1(2,1) = Mu*1i*((n*(n+1)/R2+k2-y2)*Jn1y(end)+n*dg1Jn1y(end)/R);
        Z1(2,2) = Mu*1i*n*k*Jny(end)/R;
        Z1(3,1) = k*Jn1y(end);
        Z1(3,2) = n*Jny(end)/R;
        Z1(3,3) = -dHnz(1);
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
        u(:,1) = 1i*(k*f+(n+1)*g1./r+dg1);
        u(:,2) = 1i*(n*f./r-k*g1+dg3);
        u(:,3) = df+k*g1+n*g3./r;
        v = -1i*AngularFrequency*u;
        epsilon(:,1) = 1i*k*u(:,1);
        epsilon(:,2) = (1i*n*u(:,2)+u(:,3))./r;
        epsilon(:,3) = d2f+k*dg1-n*g3./r2+n*dg3./r;
        epsilon(:,4) = 1i*(-2*n*f./r2+2*n*df./r+(n+1)*k*g1./r-k*dg1+2*d2g3+y2*g3);
        epsilon(:,5) = 1i*(2*k*df+(n*(n+1)./r2+k2-y2).*g1+n*dg1./r+n*k*g3./r);
        epsilon(:,6) = 1i*(n*u(:,1)./r+k*u(:,2));
        sigma(:,1:3) = -Lambda*kL2*f+2*Mu*epsilon(:,1:3);
        sigma(:,4:6) = Mu*epsilon(:,4:6);
        StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
        KineticEnergyDensity = Material.Density*sum(abs(v).^2,2)/2;
        PowerFlowDensity(:,1) = -real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)))/2;
        PowerFlowDensity(:,2) = -real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)))/2;
        PowerFlowDensity(:,3) = -real(sigma(:,5).*conj(v(:,1))+sigma(:,4).*conj(v(:,2))+sigma(:,3).*conj(v(:,3)))/2;
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(2*PowerFlow/pi/trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
        if  ToggleOuterFluid && ShowHalfSpace
            f = U(3)*Hnz;
            df = U(3)*dHnz;
            d2f = U(3)*((n*(n-1)./rOuterFluid.^2-z^2).*Hnz+Hn1z./rOuterFluid);
            uOuterFluid(:,1) = 1i*k*f;
            uOuterFluid(:,2) = 1i*n*f./rOuterFluid;
            uOuterFluid(:,3) = df;
            vOuterFluid = -1i*AngularFrequency*uOuterFluid;
            epsilonOuterFluid(:,1) = 1i*k*uOuterFluid(:,1);
            epsilonOuterFluid(:,2) = (1i*n*uOuterFluid(:,2)+uOuterFluid(:,3))./rOuterFluid;
            epsilonOuterFluid(:,3) = d2f;
            epsilonOuterFluid(:,6) = 0;
            sigmaOuterFluid(:,1) = -OuterFluid.Density*AngularFrequency2*f;
            sigmaOuterFluid(:,2) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,3) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,6) = 0;
            StrainEnergyDensityOuterFluid = real(sum(epsilonOuterFluid.*conj(sigmaOuterFluid),2))/2;
            KineticEnergyDensityOuterFluid = OuterFluid.Density*sum(abs(vOuterFluid).^2,2)/2;
            PowerFlowDensityOuterFluid(:,1) = -real(sigmaOuterFluid(:,1).*conj(vOuterFluid(:,1)))/2;
            PowerFlowDensityOuterFluid(:,2) = -real(sigmaOuterFluid(:,2).*conj(vOuterFluid(:,2)))/2;
            PowerFlowDensityOuterFluid(:,3) = -real(sigmaOuterFluid(:,3).*conj(vOuterFluid(:,3)))/2;
            r = [r;rOuterFluid];
            u = [u;uOuterFluid];
            epsilon = [epsilon;epsilonOuterFluid];
            sigma = [sigma;sigmaOuterFluid];
            StrainEnergyDensity = [StrainEnergyDensity;StrainEnergyDensityOuterFluid];
            KineticEnergyDensity = [KineticEnergyDensity;KineticEnergyDensityOuterFluid];
            PowerFlowDensity = [PowerFlowDensity;PowerFlowDensityOuterFluid];
        end
        if  Phase
            uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
            sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
            epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
        end
        u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
        u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
        u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
        epsilon(:,1) = epsilon(:,1)*exp(-1i*angle(epsilon(2,1)));
        epsilon(:,2) = epsilon(:,2)*exp(-1i*angle(epsilon(2,2)));
        epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
        epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
        epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
        epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
        sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
        sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
        sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
        sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
        sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
        sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
    elseif contains(Mode,'L')
        kL2 = AngularFrequency2/Material.LongitudinalVelocity_complex^2;
        x2 = kL2-k2;
        x = sqrt(x2);
        y = sqrt(y2);
        z = SignFluids*sqrt(AngularFrequency2/OuterFluid.Velocity^2-k2);
        xr = x*r;
        yr = y*r;
        zr = z*rOuterFluid;
        J0x = besselj(0,xr);
        J0y = besselj(0,yr);
        J1x = x*besselj(1,xr);
        J1y = y*besselj(1,yr);
        H0z = besselh(0,zr);
        H1z = z*besselh(1,zr);
        Z1(1,1) = Mu*1i*(k2-y2)*J1y(end);
        Z1(2,1) = k*J1y(end);
        Z1(2,2) = H1z(1);
        Z2(1,1) = Mu*2i*k*J1x(end);
        Z2(2,1) = J1x(end);
        U = Z1\Z2;
        f = J0x;
        df = -J1x;
        d2f = J1x./r-x2*J0x;
        g = U(1)*J1y;
        dg = U(1)*(-J1y./r+y2*J0y);
        u(:,1) = 1i*(k*f+g./r+dg);
        u(:,3) = df+k*g;
        v = -1i*AngularFrequency*u;
        epsilon(:,1) = 1i*k*u(:,1);
        epsilon(:,2) = u(:,3)./r;
        epsilon(:,3) = d2f+k*dg;
        epsilon(:,5) = 1i*(2*k*df+(k2-y2)*g);
        epsilon(:,6) = 0;
        sigma(:,1:3) = -Lambda*kL2*f+2*Mu*epsilon(:,1:3);
        sigma(:,5) = Mu*epsilon(:,5);
        sigma(:,6) = 0;
        StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
        KineticEnergyDensity = Material.Density*sum(abs(v).^2,2)/2;
        PowerFlowDensity(:,1) = -real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)))/2;
        PowerFlowDensity(:,3) = -real(sigma(:,5).*conj(v(:,1))+sigma(:,3).*conj(v(:,3)))/2;
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(2*PowerFlow/pi/trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
        if  ToggleOuterFluid && ShowHalfSpace
            f = U(2)*H0z;
            df = -U(2)*H1z;
            d2f = U(2)*(H1z./rOuterFluid-z^2*H0z);
            uOuterFluid(:,1) = 1i*k*f;
            uOuterFluid(:,3) = df;
            vOuterFluid = -1i*AngularFrequency*uOuterFluid;
            epsilonOuterFluid(:,1) = 1i*k*uOuterFluid(:,1);
            epsilonOuterFluid(:,2) = uOuterFluid(:,3)./rOuterFluid;
            epsilonOuterFluid(:,3) = d2f;
            epsilonOuterFluid(:,6) = 0;
            sigmaOuterFluid(:,1) = -OuterFluid.Density*AngularFrequency2*f;
            sigmaOuterFluid(:,2) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,3) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,6) = 0;
            StrainEnergyDensityOuterFluid = real(sum(epsilonOuterFluid.*conj(sigmaOuterFluid),2))/2;
            KineticEnergyDensityOuterFluid = OuterFluid.Density*sum(abs(vOuterFluid).^2,2)/2;
            PowerFlowDensityOuterFluid(:,1) = -real(sigmaOuterFluid(:,1).*conj(vOuterFluid(:,1)))/2;
            PowerFlowDensityOuterFluid(:,3) = -real(sigmaOuterFluid(:,3).*conj(vOuterFluid(:,3)))/2;
            r = [r;rOuterFluid];
            u = [u;uOuterFluid];
            epsilon = [epsilon;epsilonOuterFluid];
            sigma = [sigma;sigmaOuterFluid];
            StrainEnergyDensity = [StrainEnergyDensity;StrainEnergyDensityOuterFluid];
            KineticEnergyDensity = [KineticEnergyDensity;KineticEnergyDensityOuterFluid];
            PowerFlowDensity = [PowerFlowDensity;PowerFlowDensityOuterFluid];
        end
        if  Phase
            uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
            sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
            epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
        end
        u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
        u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
        epsilon(:,1) = epsilon(:,1)*exp(-1i*angle(epsilon(2,1)));
        epsilon(:,2) = epsilon(:,2)*exp(-1i*angle(epsilon(2,2)));
        epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
        epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
        sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
        sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
        sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
        sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
    elseif contains(Mode,'T')
        y = sqrt(y2);
        if  y == 0
            y = 1e-10;
            y2 = y^2;
        end
        yr = y*r;
        J0 = besselj(0,yr);
        J1 = y*besselj(1,yr);
        g = J0;
        dg = -J1;
        d2g = J1./r-y2*J0;
        u(:,2) = 1i*dg;
        u(:,3) = 0;
        v = -1i*AngularFrequency*u;
        epsilon(:,4) = 1i*(2*d2g+y2*g);
        epsilon(:,6) = 1i*k*u(:,2);
        sigma(:,4:6) = Mu*epsilon(:,4:6);
        StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
        KineticEnergyDensity = Material.Density*abs(v(:,2)).^2/2;
        PowerFlowDensity(:,1) = -real(sigma(:,6).*conj(v(:,2)))/2;
        PowerFlowDensity(:,3) = -real(sigma(:,4).*conj(v(:,2)))/2;
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(2*PowerFlow/pi/trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
        if  Phase
            uPhase = rad2deg(angle(u));
            sigmaPhase = rad2deg(angle(sigma));
            epsilonPhase = rad2deg(angle(epsilon));
        end
        u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
        epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
        epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
        sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
        sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
        if  ToggleOuterFluid && ShowHalfSpace
            uOuterFluid(length(rOuterFluid),3) = 0;
            epsilonOuterFluid(length(rOuterFluid),6) = 0;
            sigmaOuterFluid(length(rOuterFluid),6) = 0;
            StrainEnergyDensityOuterFluid(length(rOuterFluid),1) = 0;
            KineticEnergyDensityOuterFluid(length(rOuterFluid),1) = 0;
            PowerFlowDensityOuterFluid(length(rOuterFluid),3) = 0;
            r = [r;rOuterFluid];
            u = [u;uOuterFluid];
            epsilon = [epsilon;epsilonOuterFluid];
            sigma = [sigma;sigmaOuterFluid];
            StrainEnergyDensity = [StrainEnergyDensity;StrainEnergyDensityOuterFluid];
            KineticEnergyDensity = [KineticEnergyDensity;KineticEnergyDensityOuterFluid];
            PowerFlowDensity = [PowerFlowDensity;PowerFlowDensityOuterFluid];
            if  Phase
                uPhase = [uPhase;zeros(length(rOuterFluid),3)];
                sigmaPhase = [sigmaPhase;zeros(length(rOuterFluid),6)];
                epsilonPhase = [epsilonPhase;zeros(length(rOuterFluid),6)];
            end
        end
    end
end
u = u/sqrt(PowerFlow);
epsilon = epsilon/sqrt(PowerFlow);
sigma = sigma/sqrt(PowerFlow);
if  real(u(2,1)) < 0
    u = -u;
end
if  real(epsilon(2,1)) < 0
    epsilon = -epsilon;
end
if  real(sigma(2,1)) < 0
    sigma = -sigma;
end
if  Phase
    uPhase(round(uPhase) == -180) = 180;
    sigmaPhase(round(sigmaPhase) == -180) = 180;
    epsilonPhase(round(epsilonPhase) == -180) = 180;
end
StrainEnergyDensity = StrainEnergyDensity/PowerFlow/2;
KineticEnergyDensity = KineticEnergyDensity/PowerFlow/2;
TotalEnergyDensity = StrainEnergyDensity+KineticEnergyDensity;
PowerFlowDensity = Direction*PowerFlowDensity/PowerFlow;
end
function [PhaseVelocity,Attenuation,SignFluids,Direction,Frequency] = ExtractData(Geometry,X,Frequency,Mode)
    [Min,Max] = bounds(X(:,1));
    if  Frequency < Min || Frequency > Max
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(Min),' and ',num2str(Max),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(X(:,1)-Frequency));
        Frequency = X(q);
    end
    PhaseVelocity = X(q,4)*1e3;
    Attenuation = X(q,7);
    if  contains(Mode,'T')
        SignFluids = 0;
    else
        SignFluids = X(q,8);
        if  strcmp(Geometry,'Pipe') 
            SignFluids(2) = X(q,9);
        end
    end
    if  X(q,5) > 0
        Direction = 1;
    else
        Direction = -1; % reverse propagation direction for backward propagating modes
    end
end