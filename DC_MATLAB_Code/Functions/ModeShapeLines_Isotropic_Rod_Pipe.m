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
function ModeShapeLines_Isotropic_Rod_Pipe(FunctionMode,MakeFigure,ExportData,XSLX,TXT,MAT,Geometry,Plot,PNGresolution,Material,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Color1,Color2,Color3,Color4,Color5,Color6,F,L,T,FScholte_,LScholte,BoxLineWidth,Directory,FileName,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeLegend,Frequency,HeadLine,LegendLocation,LineWidth,Mode,PDF,PNG,Ro,Ri,SamplesR,ShowHalfSpace,HalfSpaces,Phase)
%#ok<*FNDSB>
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
p = str2double(regexp(Mode,'\d*','match'));
n = p(1);
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
        if  ~contains(Mode,'Scholte')
            q = find(F{n}{p(2)}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(F{n}{p(2)}(:,1))) || Frequency < min(F{n}{p(2)}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(F{n}{p(2)}(:,1))),' and ',num2str(ceil(max(F{n}{p(2)}(:,1)))),' kHz.'],'Error');
                    return
                else
                    [~,q] = min(abs(F{n}{p(2)}(:,1)-Frequency));
                    Frequency = F{n}{p(2)}(q,1);
                end
            end
            PhaseVelocity = F{n}{p(2)}(q,4)*1e3;
            Attenuation = F{n}{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        else
            q = find(FScholte_{n}{p(2)}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(FScholte_{n}{p(2)}(:,1))) || Frequency < min(FScholte_{n}{p(2)}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(FScholte_{n}{p(2)}(:,1))),' and ',num2str(ceil(max(FScholte_{n}{p(2)}(:,1)))),' kHz.'],'Error');
                    return
                else
                    [~,q] = min(abs(FScholte_{n}{p(2)}(:,1)-Frequency));
                    Frequency = FScholte_{n}{p(2)}(q,1);
                end
            end
            PhaseVelocity = FScholte_{n}{p(2)}(q,4)*1e3;
            Attenuation = FScholte_{n}{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        end
        AngularFrequency = 2*pi*Frequency*1e3;
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
        if  contains(Mode,'Scholte') && Attenuation ~= 0
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
        zro = zo*rOuterFluid;
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
        StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,2).*conj(sigma(:,2))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,4).*conj(sigma(:,4))+epsilon(:,5).*conj(sigma(:,5))+epsilon(:,6).*conj(sigma(:,6)));
        KineticEnergyDensity = .5*Material.Density*(abs(v(:,1)).^2+abs(v(:,2)).^2+abs(v(:,3)).^2);
        PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)));
        PowerFlowDensity(:,2) = -.5*real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)));
        PowerFlowDensity(:,3) = -.5*real(sigma(:,5).*conj(v(:,1))+sigma(:,4).*conj(v(:,2))+sigma(:,3).*conj(v(:,3)));
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
            StrainEnergyDensityInnerFluid = .5*real(epsilonInnerFluid(:,1).*conj(sigmaInnerFluid(:,1))+epsilonInnerFluid(:,2).*conj(sigmaInnerFluid(:,2))+epsilonInnerFluid(:,3).*conj(sigmaInnerFluid(:,3)));
            KineticEnergyDensityInnerFluid = .5*InnerFluid.Density*(abs(vInnerFluid(:,1)).^2+abs(vInnerFluid(:,2)).^2+abs(vInnerFluid(:,3)).^2);
            PowerFlowDensityInnerFluid(:,1) = -.5*real(sigmaInnerFluid(:,1).*conj(vInnerFluid(:,1)));
            PowerFlowDensityInnerFluid(:,2) = -.5*real(sigmaInnerFluid(:,2).*conj(vInnerFluid(:,2)));
            PowerFlowDensityInnerFluid(:,3) = -.5*real(sigmaInnerFluid(:,3).*conj(vInnerFluid(:,3)));
            r = [rInnerFluid;r];
            u = [uInnerFluid;u];
            epsilon = [epsilonInnerFluid;epsilon];
            sigma = [sigmaInnerFluid;sigma];
            StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
            KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
            PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
        end
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(PowerFlow/(.5*pi*trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity))),' m/s'])
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
            StrainEnergyDensityOuterFluid = .5*real(epsilonOuterFluid(:,1).*conj(sigmaOuterFluid(:,1))+epsilonOuterFluid(:,2).*conj(sigmaOuterFluid(:,2))+epsilonOuterFluid(:,3).*conj(sigmaOuterFluid(:,3)));
            KineticEnergyDensityOuterFluid = .5*OuterFluid.Density*(abs(vOuterFluid(:,1)).^2+abs(vOuterFluid(:,2)).^2+abs(vOuterFluid(:,3)).^2);
            PowerFlowDensityOuterFluid(:,1) = -.5*real(sigmaOuterFluid(:,1).*conj(vOuterFluid(:,1)));
            PowerFlowDensityOuterFluid(:,2) = -.5*real(sigmaOuterFluid(:,2).*conj(vOuterFluid(:,2)));
            PowerFlowDensityOuterFluid(:,3) = -.5*real(sigmaOuterFluid(:,3).*conj(vOuterFluid(:,3)));
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
        if  ~contains(Mode,'Scholte')
            q = find(L{p(2)}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(L{p(2)}(:,1))) || Frequency < min(L{p(2)}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(L{p(2)}(:,1))),' and ',num2str(ceil(max(L{p(2)}(:,1)))),' kHz.'],'Error');
                    return
                else
                    [~,q] = min(abs(L{p(2)}(:,1)-Frequency));
                    Frequency = L{p(2)}(q,1);
                end
            end
            PhaseVelocity = L{p(2)}(q,4)*1e3;
            Attenuation = L{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        else
            q = find(LScholte{p(2)}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(LScholte{p(2)}(:,1))) || Frequency < min(LScholte{p(2)}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(LScholte{p(2)}(:,1))),' and ',num2str(ceil(max(LScholte{p(2)}(:,1)))),' kHz.'],'Error');
                    return
                else
                    [~,q] = min(abs(LScholte{p(2)}(:,1)-Frequency));
                    Frequency = LScholte{p(2)}(q,1);
                end
            end
            PhaseVelocity = LScholte{p(2)}(q,4)*1e3;
            Attenuation = LScholte{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        end
        AngularFrequency = 2*pi*Frequency*1e3;
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
        if  contains(Mode,'Scholte') && Attenuation ~= 0
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
        zro = zo*rOuterFluid;
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
        StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,2).*conj(sigma(:,2))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,5).*conj(sigma(:,5)));
        KineticEnergyDensity = .5*Material.Density*(abs(v(:,1)).^2+abs(v(:,3)).^2);
        PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
        PowerFlowDensity(:,3) = -.5*real(sigma(:,5).*conj(v(:,1))+sigma(:,3).*conj(v(:,3)));
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
            StrainEnergyDensityInnerFluid = .5*real(epsilonInnerFluid(:,1).*conj(sigmaInnerFluid(:,1))+epsilonInnerFluid(:,2).*conj(sigmaInnerFluid(:,2))+epsilonInnerFluid(:,3).*conj(sigmaInnerFluid(:,3)));
            KineticEnergyDensityInnerFluid = .5*InnerFluid.Density*(abs(vInnerFluid(:,1)).^2+abs(vInnerFluid(:,3)).^2);
            PowerFlowDensityInnerFluid(:,1) = -.5*real(sigmaInnerFluid(:,1).*conj(vInnerFluid(:,1)));
            PowerFlowDensityInnerFluid(:,3) = -.5*real(sigmaInnerFluid(:,3).*conj(vInnerFluid(:,3)));
            r = [rInnerFluid;r];
            u = [uInnerFluid;u];
            epsilon = [epsilonInnerFluid;epsilon];
            sigma = [sigmaInnerFluid;sigma];
            StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
            KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
            PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
        end
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(PowerFlow/(.5*pi*trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity))),' m/s'])
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
            StrainEnergyDensityOuterFluid = .5*real(epsilonOuterFluid(:,1).*conj(sigmaOuterFluid(:,1))+epsilonOuterFluid(:,2).*conj(sigmaOuterFluid(:,2))+epsilonOuterFluid(:,3).*conj(sigmaOuterFluid(:,3)));
            KineticEnergyDensityOuterFluid = .5*OuterFluid.Density*(abs(vOuterFluid(:,1)).^2+abs(vOuterFluid(:,3)).^2);
            PowerFlowDensityOuterFluid(:,1) = -.5*real(sigmaOuterFluid(:,1).*conj(vOuterFluid(:,1)));
            PowerFlowDensityOuterFluid(:,3) = -.5*real(sigmaOuterFluid(:,3).*conj(vOuterFluid(:,3)));
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
        q = find(T{p(2)}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(T{p(2)}(:,1))) || Frequency < min(T{p(2)}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(T{p(2)}(:,1))),' and ',num2str(ceil(max(T{p(2)}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(T{p(2)}(:,1)-Frequency));
                Frequency = T{p(2)}(q,1);
            end
        end
        PhaseVelocity = T{p(2)}(q,4)*1e3;
        Attenuation = T{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k2 = k^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
        y = sqrt(kT2-k2);
        if  y == 0
            y = 1e-10;
        end
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
        u(:,2) = 1i*dg;
        u(:,3) = 0;
        v = -1i*AngularFrequency*u;
        epsilon(:,4) = 1i*(2*d2g+y2*g);
        epsilon(:,6) = 1i*k*u(:,2);
        sigma(:,4:6) = Mu*epsilon(:,4:6);
        StrainEnergyDensity = .5*real(epsilon(:,4).*conj(sigma(:,4))+epsilon(:,6).*conj(sigma(:,6)));
        KineticEnergyDensity = .5*Material.Density*abs(v(:,2)).^2;
        PowerFlowDensity(:,1) = -.5*real(sigma(:,6).*conj(v(:,2)));
        PowerFlowDensity(:,3) = -.5*real(sigma(:,4).*conj(v(:,2)));
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
        % disp(['ce: ',num2str(PowerFlow/(.5*pi*trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity))),' m/s'])
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
        if  ~contains(Mode,'Scholte')
            q = find(F{n}{p(2)}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(F{n}{p(2)}(:,1))) || Frequency < min(F{n}{p(2)}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(F{n}{p(2)}(:,1))),' and ',num2str(ceil(max(F{n}{p(2)}(:,1)))),' kHz.'],'Error');
                    return
                else
                    [~,q] = min(abs(F{n}{p(2)}(:,1)-Frequency));
                    Frequency = F{n}{p(2)}(q,1);
                end
            end
            PhaseVelocity = F{n}{p(2)}(q,4)*1e3;
            Attenuation = F{n}{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        else
            q = find(FScholte_{n}{p(2)}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(FScholte_{n}{p(2)}(:,1))) || Frequency < min(FScholte_{n}{p(2)}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(FScholte_{n}{p(2)}(:,1))),' and ',num2str(ceil(max(FScholte_{n}{p(2)}(:,1)))),' kHz.'],'Error');
                    return
                else
                    [~,q] = min(abs(FScholte_{n}{p(2)}(:,1)-Frequency));
                    Frequency = FScholte_{n}{p(2)}(q,1);
                end
            end
            PhaseVelocity = FScholte_{n}{p(2)}(q,4)*1e3;
            Attenuation = FScholte_{n}{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        end
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k2 = k^2;
        kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
        x = sqrt(kL2-k2);
        y = sqrt(kT2-k2);
        z = sqrt((AngularFrequency/OuterFluid.Velocity)^2-k2);
        if  contains(Mode,'Scholte') && Attenuation ~= 0 && PhaseVelocity < OuterFluid.Velocity
            z = -z;
        end
        x2 = kL2-k2;
        y2 = kT2-k2;
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
        StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,2).*conj(sigma(:,2))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,4).*conj(sigma(:,4))+epsilon(:,5).*conj(sigma(:,5))+epsilon(:,6).*conj(sigma(:,6)));
        KineticEnergyDensity = .5*Material.Density*(abs(v(:,1)).^2+abs(v(:,2)).^2+abs(v(:,3)).^2);
        PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,6).*conj(v(:,2))+sigma(:,5).*conj(v(:,3)));
        PowerFlowDensity(:,2) = -.5*real(sigma(:,6).*conj(v(:,1))+sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)));
        PowerFlowDensity(:,3) = -.5*real(sigma(:,5).*conj(v(:,1))+sigma(:,4).*conj(v(:,2))+sigma(:,3).*conj(v(:,3)));
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(PowerFlow/(.5*pi*trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity))),' m/s'])
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
            sigmaOuterFluid(:,1) = -OuterFluid.Density*AngularFrequency^2*f;
            sigmaOuterFluid(:,2) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,3) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,6) = 0;
            StrainEnergyDensityOuterFluid = .5*real(epsilonOuterFluid(:,1).*conj(sigmaOuterFluid(:,1))+epsilonOuterFluid(:,2).*conj(sigmaOuterFluid(:,2))+epsilonOuterFluid(:,3).*conj(sigmaOuterFluid(:,3)));
            KineticEnergyDensityOuterFluid = .5*OuterFluid.Density*(abs(vOuterFluid(:,1)).^2+abs(vOuterFluid(:,2)).^2+abs(vOuterFluid(:,3)).^2);
            PowerFlowDensityOuterFluid(:,1) = -.5*real(sigmaOuterFluid(:,1).*conj(vOuterFluid(:,1)));
            PowerFlowDensityOuterFluid(:,2) = -.5*real(sigmaOuterFluid(:,2).*conj(vOuterFluid(:,2)));
            PowerFlowDensityOuterFluid(:,3) = -.5*real(sigmaOuterFluid(:,3).*conj(vOuterFluid(:,3)));
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
        if  ~contains(Mode,'Scholte')
            q = find(L{p(2)}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(L{p(2)}(:,1))) || Frequency < min(L{p(2)}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(L{p(2)}(:,1))),' and ',num2str(ceil(max(L{p(2)}(:,1)))),' kHz.'],'Error');
                    return
                else
                    [~,q] = min(abs(L{p(2)}(:,1)-Frequency));
                    Frequency = L{p(2)}(q,1);
                end
            end
            PhaseVelocity = L{p(2)}(q,4)*1e3;
            Attenuation = L{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        else
            q = find(LScholte{p(2)}(:,1) == Frequency);
            if  isempty(q)
                if  Frequency > ceil(max(LScholte{p(2)}(:,1))) || Frequency < min(LScholte{p(2)}(:,1))
                    errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(LScholte{p(2)}(:,1))),' and ',num2str(ceil(max(LScholte{p(2)}(:,1)))),' kHz.'],'Error');
                    return
                else
                    [~,q] = min(abs(LScholte{p(2)}(:,1)-Frequency));
                    Frequency = LScholte{p(2)}(q,1);
                end
            end
            PhaseVelocity = LScholte{p(2)}(q,4)*1e3;
            Attenuation = LScholte{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        end
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k2 = k^2;
        kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
        x = sqrt(kL2-k2);
        y = sqrt(kT2-k2);
        z = sqrt((AngularFrequency/OuterFluid.Velocity)^2-k2);
        if  contains(Mode,'Scholte') && Attenuation ~= 0 && PhaseVelocity < OuterFluid.Velocity
            z = -z;
        end
        x2 = kL2-k2;
        y2 = kT2-k2;
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
        StrainEnergyDensity = .5*real(epsilon(:,1).*conj(sigma(:,1))+epsilon(:,2).*conj(sigma(:,2))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,5).*conj(sigma(:,5)));
        KineticEnergyDensity = .5*Material.Density*(abs(v(:,1)).^2+abs(v(:,3)).^2);
        PowerFlowDensity(:,1) = -.5*real(sigma(:,1).*conj(v(:,1))+sigma(:,5).*conj(v(:,3)));
        PowerFlowDensity(:,3) = -.5*real(sigma(:,5).*conj(v(:,1))+sigma(:,3).*conj(v(:,3)));
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(PowerFlow/(.5*pi*trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity))),' m/s'])
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
            sigmaOuterFluid(:,1) = -OuterFluid.Density*AngularFrequency^2*f;
            sigmaOuterFluid(:,2) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,3) = sigmaOuterFluid(:,1);
            sigmaOuterFluid(:,6) = 0;
            StrainEnergyDensityOuterFluid = .5*real(epsilonOuterFluid(:,1).*conj(sigmaOuterFluid(:,1))+epsilonOuterFluid(:,2).*conj(sigmaOuterFluid(:,2))+epsilonOuterFluid(:,3).*conj(sigmaOuterFluid(:,3)));
            KineticEnergyDensityOuterFluid = .5*OuterFluid.Density*(abs(vOuterFluid(:,1)).^2+abs(vOuterFluid(:,3)).^2);
            PowerFlowDensityOuterFluid(:,1) = -.5*real(sigmaOuterFluid(:,1).*conj(vOuterFluid(:,1)));
            PowerFlowDensityOuterFluid(:,3) = -.5*real(sigmaOuterFluid(:,3).*conj(vOuterFluid(:,3)));
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
        q = find(T{p(2)}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(T{p(2)}(:,1))) || Frequency < min(T{p(2)}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(T{p(2)}(:,1))),' and ',num2str(ceil(max(T{p(2)}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(T{p(2)}(:,1)-Frequency));
                Frequency = T{p(2)}(q,1);
            end
        end
        PhaseVelocity = T{p(2)}(q,4)*1e3;
        Attenuation = T{p(2)}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k2 = k^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
        y = sqrt(kT2-k2);
        if  y == 0
            y = 1e-10;
        end
        y2 = kT2-k2;
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
        StrainEnergyDensity = .5*real(epsilon(:,4).*conj(sigma(:,4))+epsilon(:,6).*conj(sigma(:,6)));
        KineticEnergyDensity = .5*Material.Density*abs(v(:,2)).^2;
        PowerFlowDensity(:,1) = -.5*real(sigma(:,6).*conj(v(:,2)));
        PowerFlowDensity(:,3) = -.5*real(sigma(:,4).*conj(v(:,2)));
        PowerFlow = trapz(r.^2,PowerFlowDensity(:,1))*pi;
        % disp(['ce: ',num2str(PowerFlow/(.5*pi*trapz(r.^2,StrainEnergyDensity+KineticEnergyDensity))),' m/s'])
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
    uPhase(find(round(uPhase) == -180)) = 180;
    sigmaPhase(find(round(sigmaPhase) == -180)) = 180;
    epsilonPhase(find(round(epsilonPhase) == -180)) = 180;
end
StrainEnergyDensity = .5*StrainEnergyDensity/PowerFlow;
KineticEnergyDensity = .5*KineticEnergyDensity/PowerFlow;
TotalEnergyDensity = StrainEnergyDensity+KineticEnergyDensity;
PowerFlowDensity = PowerFlowDensity/PowerFlow;
if  ExportData
    if  strcmp(Geometry,'Rod')
        Product = num2str(Frequency*2*Ro);
    else
        Product = num2str(Frequency*(Ro-Ri));
    end
    if  FunctionMode == 1
        if  ~Phase
            Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'r (mm)','uz (nm)','uq (nm)','ur (nm)'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
        else
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','uz (nm)','uq (nm)','ur (nm)','Phase uz (deg)','Phase uq (deg)','Phase ur (deg)'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
            Table(:,5) = num2cell(uPhase(:,1));
            Table(:,6) = num2cell(uPhase(:,2));
            Table(:,7) = num2cell(uPhase(:,3));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',Product,'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',Product,'MHzmm.txt']));
            end
            if  MAT
                Displacement = Table;
                save(fullfile(Directory,[FileName,'_Displacement_',Mode,'@',Product,'MHzmm.mat']),'Displacement')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end  
    elseif FunctionMode == 2
        if  ~Phase
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','sigmazz (kPa)','sigmaqq (kPa)','sigmarr (kPa)','sigmaqr (kPa)','sigmazr (kPa)','sigmazq (kPa)'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
        else
            Table = table('Size',[length(r) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','sigmazz (kPa)','sigmaqq (kPa)','sigmarr (kPa)','sigmaqr (kPa)','sigmazr (kPa)','sigmazq (kPa)','Phase sigmazz (deg)','Phase sigmaqq (deg)','Phase sigmarr (deg)','Phase sigmaqr (deg)','Phase sigmazr (deg)','Phase sigmazq (deg)'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
            Table(:,8) = num2cell(sigmaPhase(:,1));
            Table(:,9) = num2cell(sigmaPhase(:,2));
            Table(:,10) = num2cell(sigmaPhase(:,3));
            Table(:,11) = num2cell(sigmaPhase(:,4));
            Table(:,12) = num2cell(sigmaPhase(:,5));
            Table(:,13) = num2cell(sigmaPhase(:,6));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',Product,'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',Product,'MHzmm.txt']));
            end
            if  MAT
                Stress = Table;
                save(fullfile(Directory,[FileName,'_Stress_',Mode,'@',Product,'MHzmm.mat']),'Stress')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end   
    elseif FunctionMode == 3
        if  ~Phase
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','epsilonzz','epsilonqq','epsilonrr','epsilonqr','epsilonzr','epsilonzq'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
        else
            Table = table('Size',[length(r) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','epsilonzz','epsilonqq','epsilonrr','epsilonqr','epsilonzr','epsilonzq','Phase epsilonzz (deg)','Phase epsilonqq (deg)','Phase epsilonrr (deg)','Phase epsilonqr (deg)','Phase epsilonzr (deg)','Phase epsilonzq (deg)'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
            Table(:,8) = num2cell(epsilonPhase(:,1));
            Table(:,9) = num2cell(epsilonPhase(:,2));
            Table(:,10) = num2cell(epsilonPhase(:,3));
            Table(:,11) = num2cell(epsilonPhase(:,4));
            Table(:,12) = num2cell(epsilonPhase(:,5));
            Table(:,13) = num2cell(epsilonPhase(:,6));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',Product,'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',Product,'MHzmm.txt']));
            end
            if  MAT
                Strain = Table;
                save(fullfile(Directory,[FileName,'_Strain_',Mode,'@',Product,'MHzmm.mat']),'Strain')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end   
    elseif FunctionMode == 4
        Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'r (mm)','Estrain (J/m2)','Ekin (J/m2)','Etotal (J/m2)'});
        Table(:,1) = num2cell(1e3*r);
        Table(:,2) = num2cell(StrainEnergyDensity);
        Table(:,3) = num2cell(KineticEnergyDensity);
        Table(:,4) = num2cell(TotalEnergyDensity);
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',Product,'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',Product,'MHzmm.txt']));
            end
            if  MAT
                EnergyDensity = Table;
                save(fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',Product,'MHzmm.mat']),'EnergyDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end        
    elseif FunctionMode == 5
        Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'r (mm)','Pz (W/m)','Pq (W/m)','Pr (W/m)'});
        Table(:,1) = num2cell(1e3*r);
        Table(:,2) = num2cell(PowerFlowDensity(:,1));
        Table(:,3) = num2cell(PowerFlowDensity(:,2));
        Table(:,4) = num2cell(PowerFlowDensity(:,3));
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',Product,'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',Product,'MHzmm.txt']));
            end
            if  MAT
                PowerFlowDensity = Table;
                save(fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',Product,'MHzmm.mat']),'PowerFlowDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    end
end
if  MakeFigure
    if  Phase
        if  FunctionMode == 1 
            f = figure('Name','Displacement phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        elseif FunctionMode == 2
            f = figure('Name','Stress phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        elseif FunctionMode == 3
            f = figure('Name','Strain phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        end   
        datacursormode on
        jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
        jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));    
        x = xline(0,'Color',[.6 .6 .6]); 
        hasbehavior(x,'legend',false);
        hold on
        if  FunctionMode == 1
            if  Plot(3)
                plot(uPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),r*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),r*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),r*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),r*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),r*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),r*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end
        end    
        ax = gca;
        ax.Box = 'on';
        ax.LineWidth = BoxLineWidth;
        ax.FontSize = FontSizeAxes;
        ax.Title.Interpreter = 'latex';
        ax.Title.FontSize = FontSizeHeadLine;
        if  HeadLine
            if  ~contains(Mode,'Scholte')
                ModeName = [Mode(1),'(',num2str(n),',',num2str(p(2)),')'];
            else
                ModeName = [Mode(1),'(',num2str(n),',',num2str(p(2)),')$^\mathrm{Scholte}$'];
            end
            if  strcmp(Geometry,'Rod')
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
            elseif strcmp(Geometry,'Pipe')
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,$\times$\,',num2str((Ro-Ri)*1e3),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
            end
            if  FluidLoading
                if  strcmp(Geometry,'Rod')
                    String = append(String,' in ',replace(OuterFluid.Name,'_','\_'));
                elseif strcmp(Geometry,'Pipe')
                    if  ToggleOuterFluid && ToggleInnerFluid
                        String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/',replace(InnerFluid.Name,'_','\_'));
                    elseif ToggleOuterFluid && ~ToggleInnerFluid
                        String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/vacuum');
                    elseif ~ToggleOuterFluid && ToggleInnerFluid
                        String = append(String,' in vacuum/',replace(InnerFluid.Name,'_','\_'));
                    end
                end
            end
            ax.Title.String = String;
        end
        ax.XLabel.Interpreter = 'latex';
        ax.XLabel.FontSize = FontSizeAxesLabels;
        ax.XLim = max(abs(ax.XLim))*[-1 1];
        ax.YLabel.Interpreter = 'latex';
        ax.YLabel.FontSize = FontSizeAxesLabels;
        ax.YLabel.String = '$r$ (mm)';
        ax.TickLabelInterpreter = 'latex';
        if  strcmp(Geometry,'Rod') || (strcmp(Geometry,'pipe') && ToggleInnerFluid)
            ax.YLim = [0 r(end)*1e3];
        else
            ax.YLim = [r(1) r(end)]*1e3;
        end
        if  strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
            if  strcmp(OuterFluid.Name,InnerFluid.Name)
                OuterFaceAlpha = .2;
                InnerFaceAlpha = .2;
            else
                if  OuterFluid.Density*OuterFluid.Velocity > InnerFluid.Density*InnerFluid.Velocity
                    OuterFaceAlpha = .2;
                    InnerFaceAlpha = .1;
                else
                    OuterFaceAlpha = .1;
                    InnerFaceAlpha = .2;
                end
            end
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[Ro*1e3 Ro*1e3 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',OuterFaceAlpha,'EdgeColor','none')
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 Ri*1e3 Ri*1e3],'b','FaceAlpha',InnerFaceAlpha,'EdgeColor','none')
        elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[Ro*1e3 Ro*1e3 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none')
        elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 Ri*1e3 Ri*1e3],'b','FaceAlpha',.2,'EdgeColor','none')
        end
        if  FunctionMode == 1
            if  Plot(3)
                plot(uPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),r*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),r*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),r*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),r*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),r*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),r*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end
        end
        if  FunctionMode == 1
            ax.XLabel.String = 'Displacement phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$u_r$','$u_z$','$u_\theta$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
            else
                LegendNames = {'Radial ($u_r$)','Axial ($u_z$)','Circumferential ($u_\theta$)'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')
            end
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        elseif FunctionMode == 2
            ax.XLabel.String = 'Stress phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$\sigma_{rr}$','$\sigma_{zz}$','$\sigma_{\theta\theta}$','$\sigma_{zr}$','$\sigma_{\theta r}$','$\sigma_{z\theta}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Radial ($\sigma_{rr}$)','Axial ($\sigma_{zz}$)','Circumferential ($\sigma_{\theta\theta}$)','Shear ($\sigma_{zr}$)','Shear ($\sigma_{\theta r}$)','Shear ($\sigma_{z\theta}$)'}; 
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
            end
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        elseif FunctionMode == 3
            ax.XLabel.String = 'Strain phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$\varepsilon_{rr}$','$\varepsilon_{zz}$','$\varepsilon_{\theta\theta}$','$\varepsilon_{zr}$','$\varepsilon_{\theta r}$','$\varepsilon_{z\theta}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Radial ($\varepsilon_{rr}$)','Axial ($\varepsilon_{zz}$)','Circumferential ($\varepsilon_{\theta\theta}$)','Shear ($\varepsilon_{zr}$)','Shear ($\varepsilon_{\theta r}$)','Shear ($\varepsilon_{z\theta}$)'}; 
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
            end 
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        end
        tb = axtoolbar('default');
        tb.Visible = 'on';
        d = datacursormode(f);
        d.Interpreter = 'latex';
        d.UpdateFcn = @CursorPhase;
    end
    if  FunctionMode == 1
        f = figure('Name','Displacement','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 2
        f = figure('Name','Stress','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 3
        f = figure('Name','Strain','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 4
        f = figure('Name','Energy density','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 5 
        f = figure('Name','Power flow density','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    end
    datacursormode on
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));    
    x = xline(0,'Color',[.6 .6 .6]); 
    hasbehavior(x,'legend',false);
    hold on
    if  FunctionMode == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),r*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),r*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),r*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
    end    
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.Interpreter = 'latex';
    ax.Title.FontSize = FontSizeHeadLine;
    if  HeadLine
        if  ~contains(Mode,'Scholte')
            ModeName = [Mode(1),'(',num2str(n),',',num2str(p(2)),')'];
        else
            ModeName = [Mode(1),'(',num2str(n),',',num2str(p(2)),')$^\mathrm{Scholte}$'];
        end
        if  strcmp(Geometry,'Rod')
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
        elseif strcmp(Geometry,'Pipe')
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,$\times$\,',num2str((Ro-Ri)*1e3),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
        end
        if  FluidLoading
            if  strcmp(Geometry,'Rod')
                String = append(String,' in ',replace(OuterFluid.Name,'_','\_'));
            elseif strcmp(Geometry,'Pipe')
                if  ToggleOuterFluid && ToggleInnerFluid
                    String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/',replace(InnerFluid.Name,'_','\_'));
                elseif ToggleOuterFluid && ~ToggleInnerFluid
                    String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/vacuum');
                elseif ~ToggleOuterFluid && ToggleInnerFluid
                    String = append(String,' in vacuum/',replace(InnerFluid.Name,'_','\_'));
                end
            end
        end
        ax.Title.String = String;
    end
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.XLim = max(abs(ax.XLim))*[-1 1];
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.String = '$r$ (mm)';    
    ax.TickLabelInterpreter = 'latex';    
    if  strcmp(Geometry,'Rod') || (strcmp(Geometry,'Pipe') && ToggleInnerFluid)
        ax.YLim = [0 r(end)*1e3];
    else
        ax.YLim = [r(1) r(end)]*1e3;
    end
    if  strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
        if  strcmp(OuterFluid.Name,InnerFluid.Name)
            OuterFaceAlpha = .2;
            InnerFaceAlpha = .2;
        else
            if  OuterFluid.Density*OuterFluid.Velocity > InnerFluid.Density*InnerFluid.Velocity
                OuterFaceAlpha = .2;
                InnerFaceAlpha = .1;
            else
                OuterFaceAlpha = .1;
                InnerFaceAlpha = .2;
            end
        end
        patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[Ro*1e3 Ro*1e3 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',OuterFaceAlpha,'EdgeColor','none')
        patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 Ri*1e3 Ri*1e3],'b','FaceAlpha',InnerFaceAlpha,'EdgeColor','none')
    elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
        patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[Ro*1e3 Ro*1e3 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none')
    elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
        patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 Ri*1e3 Ri*1e3],'b','FaceAlpha',.2,'EdgeColor','none')
    end
    if  FunctionMode == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),r*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),r*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),r*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
    end
    if  FunctionMode == 1
        ax.XLabel.String = 'Displacement (nm)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$u_r$','$u_z$','$u_\theta$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
        else
            LegendNames = {'Radial ($u_r$)','Axial ($u_z$)','Circumferential ($u_\theta$)'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')
        end
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 2
        ax.XLabel.String = 'Stress (kPa)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$\sigma_{rr}$','$\sigma_{zz}$','$\sigma_{\theta\theta}$','$\sigma_{zr}$','$\sigma_{\theta r}$','$\sigma_{z\theta}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Radial ($\sigma_{rr}$)','Axial ($\sigma_{zz}$)','Circumferential ($\sigma_{\theta\theta}$)','Shear ($\sigma_{zr}$)','Shear ($\sigma_{\theta r}$)','Shear ($\sigma_{z\theta}$)'}; 
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
        end
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 3
        ax.XLabel.String = 'Strain';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$\varepsilon_{rr}$','$\varepsilon_{zz}$','$\varepsilon_{\theta\theta}$','$\varepsilon_{zr}$','$\varepsilon_{\theta r}$','$\varepsilon_{z\theta}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Radial ($\varepsilon_{rr}$)','Axial ($\varepsilon_{zz}$)','Circumferential ($\varepsilon_{\theta\theta}$)','Shear ($\varepsilon_{zr}$)','Shear ($\varepsilon_{\theta r}$)','Shear ($\varepsilon_{z\theta}$)'}; 
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
        end 
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    elseif FunctionMode == 4
        ax.XLabel.String = 'Energy density (J/m$^2$)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$E_\mathrm{strain}$','$E_\mathrm{kin}$','$E_\mathrm{total}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Strain energy density','Kinetic energy density','Total energy density'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')        
        end       
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end 
    elseif FunctionMode == 5
        ax.XLabel.String = 'Power flow density (W/m)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$p_r$','$p_z$','$p_\theta$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Radial ($p_r$)','Axial ($p_z$)','Circumferential ($p_\theta$)'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')        
        end
        if  Export
            try
                if  PDF
                    exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
                end
                if  PNG
                    exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                return
            end
        end
    end
    tb = axtoolbar('default');
    tb.Visible = 'on';
    d = datacursormode(f);
    d.Interpreter = 'latex';
    d.UpdateFcn = @Cursor;
end
function output_txt = Cursor(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = replace(InnerFluid.Name,'_','\_');
        else
            output_txt = replace(OuterFluid.Name,'_','\_');
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$u_z$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$u_\theta$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$u_r$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\sigma_{rr}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\sigma_{zr}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\sigma_{\theta r}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\sigma_{zz}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\sigma_{\theta\theta}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\sigma_{z\theta}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varepsilon_{rr}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varepsilon_{zr}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varepsilon_{\theta r}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varepsilon_{zz}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varepsilon_{\theta\theta}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varepsilon_{z\theta}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 4
            if  event_obj.Target.Color == Color1
                output_txt = {['$E_\mathrm{strain}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$E_\mathrm{kin}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$E_\mathrm{total}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 5 
            if  event_obj.Target.Color == Color1
                output_txt = {['$p_z$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$p_\theta$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$p_r$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        end
    end
end
function output_txt = CursorPhase(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = replace(InnerFluid.Name,'_','\_');
        else
            output_txt = replace(OuterFluid.Name,'_','\_');
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$\varphi(u_z)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(u_\theta)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$\varphi(u_r)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\sigma_{rr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\sigma_{zr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\sigma_{\theta r})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\sigma_{zz})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\sigma_{\theta\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\sigma_{z\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\varepsilon_{rr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\varepsilon_{zr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\varepsilon_{\theta r})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\varepsilon_{zz})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\varepsilon_{\theta\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\varepsilon_{z\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        end
    end
end
end