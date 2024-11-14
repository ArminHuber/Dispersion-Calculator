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
function [Time,u,ua,ub,x1,r,Thetaa,Thetab,p] = ModeShapeGridAnimationComputer_Isotropic_Rod_Pipe(Geometry,Plane,Material,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,F,L,T,FScholte_,LScholte,CycleDuration,FrameRate,Frequency,Length,Mode,Ro,Ri,SamplesZ,SamplesR,ShowHalfSpace,HalfSpaces)
SamplesTheta = 4*SamplesR;
ThetaRange = [.3 -.8]*pi;

%#ok<*AGROW>
u = {};
ua = {};
ub = {};
Thetaa = (2*pi*(-.5:1/SamplesTheta:.5))'; % for the cross section
Thetab = Thetaa(Thetaa < ThetaRange(1) & Thetaa > ThetaRange(2)); % for the visible outer cylinder surface
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
Time = (0:1/(Frequency*1e3*CycleDuration*FrameRate):1/(Frequency*1e3))';
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
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesZ*Frequency):Length*PhaseVelocity/(1e3*Frequency);
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
        g1 = U(2)*Zn1y+U(3)*Wn1y; % SV_in + SV_out
        dg1 = U(2)*dg1Zn1y+U(3)*dg1Wn1y;
        g3 = U(4)*Zny+U(5)*Wny; % SH_in + SH_out
        dg3 = U(4)*dZny+U(5)*dWny;
        u1 = 1i*(k*f+(n+1)*g1./r+dg1);
        u2 = 1i*(n*f./r-k*g1+dg3);
        u3 = df+k*g1+n*g3./r;
        if  strcmp(Plane,'r-z')
            E = exp(1i*(k*x1-AngularFrequency*Time));
            for i = 1:length(x1)
                for j = 1:length(Time)
                    u{j,i}(:,1) = u1*E(j,i);
                    u{j,i}(:,3) = u3*E(j,i);
                end
            end
            if  ToggleInnerFluid
                f = U(6)*Znzi;
                df = U(6)*dZnzi;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    for j = 1:length(Time)
                        uInnerFluid{j,i}(:,1) = u1*E(j,i);
                        uInnerFluid{j,i}(:,3) = df*E(j,i);
                        u{j,i} = [uInnerFluid{j,i};u{j,i}];
                    end
                end
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(7)*Hnzo;
                df = U(7)*dHnzo;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    for j = 1:length(Time)
                        uOuterFluid{j,i}(:,1) = u1*E(j,i);
                        uOuterFluid{j,i}(:,3) = df*E(j,i);
                        u{j,i} = [u{j,i};uOuterFluid{j,i}];
                    end
                end
                r = [r;rOuterFluid];
            end
        else
            Ea = exp(1i*(n*Thetaa'-AngularFrequency*Time));
            for i = 1:length(Thetaa)
                for j = 1:length(Time)
                    ua{j,i}(:,1) = u1*Ea(j,i);
                    ua{j,i}(:,2) = u2*Ea(j,i);
                    ua{j,i}(:,3) = u3*Ea(j,i);
                end
            end
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                for i = 1:length(Thetab)
                    Eb = exp(1i*(k*x1'+n*Thetab(i)-AngularFrequency*Time'));
                    for j = 1:length(Time)
                        ub{j,i}(:,1) = u1(end)*Eb(:,j);
                        ub{j,i}(:,2) = u2(end)*Eb(:,j);
                        ub{j,i}(:,3) = u3(end)*Eb(:,j);
                    end
                end
            end
            if  ToggleInnerFluid
                f = U(6)*Znzi;
                df = U(6)*dZnzi;
                u1 = 1i*k*f;
                u2 = 1i*n*f./rInnerFluid;
                for i = 1:length(Thetaa)
                    for j = 1:length(Time)
                        uInnerFluida{j,i}(:,1) = u1*Ea(j,i);
                        uInnerFluida{j,i}(:,2) = u2*Ea(j,i);
                        uInnerFluida{j,i}(:,3) = df*Ea(j,i);
                        ua{j,i} = [uInnerFluida{j,i};ua{j,i}];
                    end
                end
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(7)*Hnzo;
                df = U(7)*dHnzo;
                u1 = 1i*k*f;
                u2 = 1i*n*f./rOuterFluid;
                for i = 1:length(Thetaa)
                    for j = 1:length(Time)
                        uOuterFluida{j,i}(:,1) = u1*Ea(j,i);
                        uOuterFluida{j,i}(:,2) = u2*Ea(j,i);
                        uOuterFluida{j,i}(:,3) = df*Ea(j,i);
                        ua{j,i} = [ua{j,i};uOuterFluida{j,i}];
                    end
                end
                for i = 1:length(Thetab)
                    Eb = exp(1i*(k*x1'+n*Thetab(i)-AngularFrequency*Time'));
                    for j = 1:length(Time)
                        ub{j,i}(:,1) = u1(end)*Eb(:,j);
                        ub{j,i}(:,2) = u2(end)*Eb(:,j);
                        ub{j,i}(:,3) = df(end)*Eb(:,j);
                    end
                end
                r = [r;rOuterFluid];
            end
        end
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
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesZ*Frequency):Length*PhaseVelocity/(1e3*Frequency);
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
        g = U(2)*Z1y+U(3)*W1y; % SV_in + SV_out
        dg = U(2)*(-Z1y./r+y2*Z0y)+U(3)*(-W1y./r+y2*W0y);
        u1 = 1i*(k*f+g./r+dg);
        u3 = df+k*g;
        E = exp(1i*(k*x1-AngularFrequency*Time));
        if  strcmp(Plane,'r-z')
            for i = 1:length(x1)
                for j = 1:length(Time)
                    u{j,i}(:,1) = u1*E(j,i);
                    u{j,i}(:,3) = u3*E(j,i);
                end
            end
            if  ToggleInnerFluid
                f = U(4)*Z0zi;
                df = -U(4)*Z1zi;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    for j = 1:length(Time)
                        uInnerFluid{j,i}(:,1) = u1*E(j,i);
                        uInnerFluid{j,i}(:,3) = df*E(j,i);
                        u{j,i} = [uInnerFluid{j,i};u{j,i}];
                    end
                end
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(5)*H0zo;
                df = -U(5)*H1zo;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    for j = 1:length(Time)
                        uOuterFluid{j,i}(:,1) = u1*E(j,i);
                        uOuterFluid{j,i}(:,3) = df*E(j,i);
                        u{j,i} = [u{j,i};uOuterFluid{j,i}];
                    end
                end
                r = [r;rOuterFluid];
            end
        else
            for j = 1:length(Time)
                ua{j,1}(:,1) = u1*E(j,1);
                ua{j,1}(:,3) = u3*E(j,1);
            end
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                for j = 1:length(Time)
                    ub{j,1}(:,1) = u1(end)*E(j,:);
                    ub{j,1}(:,3) = u3(end)*E(j,:);
                end
            end
            if  ToggleInnerFluid
                f = U(4)*Z0zi;
                df = -U(4)*Z1zi;
                u1 = 1i*k*f;
                for j = 1:length(Time)
                    uInnerFluida{j,1}(:,1) = u1*E(j,1);
                    uInnerFluida{j,1}(:,3) = df*E(j,1);
                    ua{j,1} = [uInnerFluida{j};ua{j}];
                end
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(5)*H0zo;
                df = -U(5)*H1zo;
                u1 = 1i*k*f;
                for j = 1:length(Time)
                    uOuterFluida{j,1}(:,1) = u1*E(j,1);
                    uOuterFluida{j,1}(:,3) = df*E(j,1);
                    ua{j,1} = [ua{j};uOuterFluida{j}];
                    ub{j,1}(:,1) = u1(end)*E(j,:);
                    ub{j,1}(:,3) = df(end)*E(j,:);
                end
                r = [r;rOuterFluid];
            end
        end
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
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesZ*Frequency):Length*PhaseVelocity/(1e3*Frequency);
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k2 = k^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
        y = sqrt(kT2-k2);
        if  y == 0
            y = 1e-10;
        end
        yr = y*r;
        U = bessely(2,yr(end))\-besselj(2,yr(end));
        J1 = y*besselj(1,yr);
        Y1 = y*bessely(1,yr);
        dg = -J1-U*Y1; % SH_in + SH_out
        if  strcmp(Plane,'r-z')
            if  ToggleInnerFluid
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                r = [r;rOuterFluid];
            end
        else
            E = exp(1i*(k*x1-AngularFrequency*Time));
            u2 = 1i*dg*exp(.5i*pi);
            for j = 1:length(Time)
                ua{j,1}(:,2) = u2*E(j,1);
                ua{j,1}(:,3) = 0;
            end
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                for j = 1:length(Time)
                    ub{j,1}(:,2) = u2(end)*E(j,:);
                    ub{j,1}(:,3) = 0;
                end
            end
            if  ToggleInnerFluid
                for j = 1:length(Time)
                    ua{j,1} = [zeros(length(rInnerFluid),3);ua{j,1}];
                end
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                for j = 1:length(Time)
                    ua{j,1} = [ua{j,1};zeros(length(rOuterFluid),3)];
                    ub{j,1}(length(x1),3) = 0;
                end
                r = [r;rOuterFluid];
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
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesZ*Frequency):Length*PhaseVelocity/(1e3*Frequency);
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
        g1 = U(1)*Jn1y;
        dg1 = U(1)*dg1Jn1y;
        g3 = U(2)*Jny;
        dg3 = U(2)*dJny;
        u1 = 1i*(k*f+(n+1)*g1./r+dg1);
        u2 = 1i*(n*f./r-k*g1+dg3);
        u3 = df+k*g1+n*g3./r;
        if  strcmp(Plane,'r-z')
            E = exp(1i*(k*x1-AngularFrequency*Time));
            for i = 1:length(x1)
                for j = 1:length(Time)
                    u{j,i}(:,1) = u1*E(j,i);
                    u{j,i}(:,3) = u3*E(j,i);
                end
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(3)*Hnz;
                df = U(3)*dHnz;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    for j = 1:length(Time)
                        uOuterFluid{j,i}(:,1) = u1*E(j,i);
                        uOuterFluid{j,i}(:,3) = df*E(j,i);
                        u{j,i} = [u{j,i};uOuterFluid{j,i}];
                    end
                end
                r = [r;rOuterFluid];
            end
        else
            Ea = exp(1i*(n*Thetaa'-AngularFrequency*Time));
            for i = 1:length(Thetaa)
                for j = 1:length(Time)
                    ua{j,i}(:,1) = u1*Ea(j,i);
                    ua{j,i}(:,2) = u2*Ea(j,i);
                    ua{j,i}(:,3) = u3*Ea(j,i);
                end
            end
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                for i = 1:length(Thetab)
                    Eb = exp(1i*(k*x1'+n*Thetab(i)-AngularFrequency*Time'));
                    for j = 1:length(Time)
                        ub{j,i}(:,1) = u1(end)*Eb(:,j);
                        ub{j,i}(:,2) = u2(end)*Eb(:,j);
                        ub{j,i}(:,3) = u3(end)*Eb(:,j);
                    end
                end
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(3)*Hnz;
                df = U(3)*dHnz;
                u1 = 1i*k*f;
                u2 = 1i*n*f./rOuterFluid;
                for i = 1:length(Thetaa)
                    for j = 1:length(Time)
                        uOuterFluida{j,i}(:,1) = u1*Ea(j,i);
                        uOuterFluida{j,i}(:,2) = u2*Ea(j,i);
                        uOuterFluida{j,i}(:,3) = df*Ea(j,i);
                        ua{j,i} = [ua{j,i};uOuterFluida{j,i}];
                    end
                end
                for i = 1:length(Thetab)
                    Eb = exp(1i*(k*x1'+n*Thetab(i)-AngularFrequency*Time'));
                    for j = 1:length(Time)
                        ub{j,i}(:,1) = u1(end)*Eb(:,j);
                        ub{j,i}(:,2) = u2(end)*Eb(:,j);
                        ub{j,i}(:,3) = df(end)*Eb(:,j);
                    end
                end
                r = [r;rOuterFluid];
            end
        end
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
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesZ*Frequency):Length*PhaseVelocity/(1e3*Frequency);
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
        g = U(1)*J1y;
        dg = U(1)*(-J1y./r+y2*J0y);
        u1 = 1i*(k*f+g./r+dg);
        u3 = df+k*g;
        E = exp(1i*(k*x1-AngularFrequency*Time));
        if  strcmp(Plane,'r-z')
            for i = 1:length(x1)
                for j = 1:length(Time)
                    u{j,i}(:,1) = u1*E(j,i);
                    u{j,i}(:,3) = u3*E(j,i);
                end
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(2)*H0z;
                df = -U(2)*H1z;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    for j = 1:length(Time)
                        uOuterFluid{j,i}(:,1) = u1*E(j,i);
                        uOuterFluid{j,i}(:,3) = df*E(j,i);
                        u{j,i} = [u{j,i};uOuterFluid{j,i}];
                    end
                end
                r = [r;rOuterFluid];
            end
        else
            for j = 1:length(Time)
                ua{j,1}(:,1) = u1*E(j,1);
                ua{j,1}(:,3) = u3*E(j,1);
            end
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                for j = 1:length(Time)
                    ub{j,1}(:,1) = u1(end)*E(j,:);
                    ub{j,1}(:,3) = u3(end)*E(j,:);
                end
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(2)*H0z;
                df = -U(2)*H1z;
                u1 = 1i*k*f;
                for j = 1:length(Time)
                    uOuterFluida{j,1}(:,1) = u1*E(j,1);
                    uOuterFluida{j,1}(:,3) = df*E(j,1);
                    ua{j,1} = [ua{j};uOuterFluida{j}];
                    ub{j,1}(:,1) = u1(end)*E(j,:);
                    ub{j,1}(:,3) = df(end)*E(j,:);
                end
                r = [r;rOuterFluid];
            end
        end
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
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesZ*Frequency):Length*PhaseVelocity/(1e3*Frequency);
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k2 = k^2;
        kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
        y = sqrt(kT2-k2);
        if  y == 0
            y = 1e-10;
        end
        yr = y*r;
        J1 = y*besselj(1,yr);
        dg = -J1;
        if  strcmp(Plane,'r-z')
            if  ToggleOuterFluid && ShowHalfSpace
                r = [r;rOuterFluid];
            end
        else
            E = exp(1i*(k*x1-AngularFrequency*Time));
            u2 = 1i*dg*exp(.5i*pi);
            for j = 1:length(Time)
                ua{j,1}(:,2) = u2*E(j,1);
                ua{j,1}(:,3) = 0;
            end
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                for j = 1:length(Time)
                    ub{j,1}(:,2) = u2(end)*E(j,:);
                    ub{j,1}(:,3) = 0;
                end
            end
            if  ToggleOuterFluid && ShowHalfSpace
                for j = 1:length(Time)
                    ua{j,1} = [ua{j,1};zeros(length(rOuterFluid),3)];
                    ub{j,1}(length(x1),3) = 0;
                end
                r = [r;rOuterFluid];
            end
        end
    end
end