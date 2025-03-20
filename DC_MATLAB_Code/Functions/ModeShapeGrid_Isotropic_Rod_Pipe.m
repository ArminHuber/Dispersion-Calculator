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
function ModeShapeGrid_Isotropic_Rod_Pipe(Geometry,Plane,PNGresolution,Material,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,F,L,T,FScholte_,LScholte,Directory,Export,FontSizeHeadLine,FontSizeAxesLabels,Frequency,GridLine,HeadLine,Length,LineWidth,Mode,FileName,PDF,PNG,Ro,Ri,SamplesZ,SamplesR,Gain,Undistorted,ShowHalfSpace,HalfSpaces)
SamplesTheta = 4*SamplesR;
ThetaRange = [.3 -.8]*pi;
Azimuth = -82;
Elevation = 8;

%#ok<*AGROW>
if  strcmp(Plane,'r-theta')
    Thetaa = (2*pi*(-.5:1/SamplesTheta:.5))'; % for the cross section
    Thetab = Thetaa(Thetaa < ThetaRange(1) & Thetaa > ThetaRange(2)); % for the visible outer cylinder surface
end
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
            E = exp(1i*k*x1);
            for i = 1:length(x1)
                u{i,1}(:,1) = u1*E(i);
                u{i,1}(:,3) = u3*E(i);
            end
            if  ToggleInnerFluid
                f = U(6)*Znzi;
                df = U(6)*dZnzi;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    uInnerFluid{i,1}(:,1) = u1*E(i);
                    uInnerFluid{i,1}(:,3) = df*E(i);
                    u{i,1} = [uInnerFluid{i};u{i}];
                end
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(7)*Hnzo;
                df = U(7)*dHnzo;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    uOuterFluid{i,1}(:,1) = u1*E(i);
                    uOuterFluid{i,1}(:,3) = df*E(i);
                    u{i,1} = [u{i};uOuterFluid{i}];
                end
                r = [r;rOuterFluid];
            end
        else
            Ea = exp(1i*n*Thetaa);
            Eb = exp(1i*(k*x1+n*Thetab));
            for i = 1:length(Thetaa)
                ua{i,1}(:,1) = u1*Ea(i);
                ua{i,1}(:,2) = u2*Ea(i);
                ua{i,1}(:,3) = u3*Ea(i);
            end
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                for i = 1:length(Thetab)
                    ub{i,1}(:,1) = u1(end)*Eb(i,:);
                    ub{i,1}(:,2) = u2(end)*Eb(i,:);
                    ub{i,1}(:,3) = u3(end)*Eb(i,:);
                end
            end
            if  ToggleInnerFluid
                f = U(6)*Znzi;
                df = U(6)*dZnzi;
                u1 = 1i*k*f;
                u2 = 1i*n*f./rInnerFluid;
                for i = 1:length(Thetaa)
                    uInnerFluida{i,1}(:,1) = u1*Ea(i);
                    uInnerFluida{i,1}(:,2) = u2*Ea(i);
                    uInnerFluida{i,1}(:,3) = df*Ea(i);
                    ua{i,1} = [uInnerFluida{i};ua{i}];
                end
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(7)*Hnzo;
                df = U(7)*dHnzo;
                u1 = 1i*k*f;
                u2 = 1i*n*f./rOuterFluid;
                for i = 1:length(Thetaa)
                    uOuterFluida{i,1}(:,1) = u1*Ea(i);
                    uOuterFluida{i,1}(:,2) = u2*Ea(i);
                    uOuterFluida{i,1}(:,3) = df*Ea(i);
                    ua{i,1} = [ua{i};uOuterFluida{i}];
                end
                for i = 1:length(Thetab)
                    ub{i,1}(:,1) = u1(end)*Eb(i,:);
                    ub{i,1}(:,2) = u2(end)*Eb(i,:);
                    ub{i,1}(:,3) = df(end)*Eb(i,:);
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
        E = exp(1i*k*x1);
        if  strcmp(Plane,'r-z')
            for i = 1:length(x1)
                u{i,1}(:,1) = u1*E(i);
                u{i,1}(:,3) = u3*E(i);
            end
            if  ToggleInnerFluid
                f = U(4)*Z0zi;
                df = -U(4)*Z1zi;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    uInnerFluid{i,1}(:,1) = u1*E(i);
                    uInnerFluid{i,1}(:,3) = df*E(i);
                    u{i,1} = [uInnerFluid{i};u{i}];
                end
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(5)*H0zo;
                df = -U(5)*H1zo;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    uOuterFluid{i,1}(:,1) = u1*E(i);
                    uOuterFluid{i,1}(:,3) = df*E(i);
                    u{i,1} = [u{i};uOuterFluid{i}];
                end
                r = [r;rOuterFluid];
            end
        else
            ua(:,1) = u1;
            ua(:,3) = u3;
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                ub(:,1) = u1(end)*E;
                ub(:,3) = u3(end)*E;
            end
            if  ToggleInnerFluid
                f = U(4)*Z0zi;
                df = -U(4)*Z1zi;
                uInnerFluida(:,1) = 1i*k*f;
                uInnerFluida(:,3) = df;
                ua = [uInnerFluida;ua];
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(5)*H0zo;
                df = -U(5)*H1zo;
                uOuterFluida(:,1) = 1i*k*f;
                uOuterFluida(:,3) = df;
                ua = [ua;uOuterFluida];
                ub(:,1) = 1i*k*f(end)*E;
                ub(:,3) = df(end)*E;
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
            ua(:,2) = 1i*dg*exp(.5i*pi);
            ua(:,3) = 0;
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                ub(:,2) = ua(end,2)*exp(1i*k*x1);
                ub(:,3) = 0;
            end
            if  ToggleInnerFluid
                ua = [zeros(length(rInnerFluid),3);ua];
                r = [rInnerFluid;r];
            end
            if  ToggleOuterFluid && ShowHalfSpace
                ua = [ua;zeros(length(rOuterFluid),3)];
                ub(length(x1),3) = 0;
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
            E = exp(1i*k*x1);
            for i = 1:length(x1)
                u{i,1}(:,1) = u1*E(i);
                u{i,1}(:,3) = u3*E(i);
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(3)*Hnz;
                df = U(3)*dHnz;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    uOuterFluid{i,1}(:,1) = u1*E(i);
                    uOuterFluid{i,1}(:,3) = df*E(i);
                    u{i,1} = [u{i};uOuterFluid{i}];
                end
                r = [r;rOuterFluid];
            end
        else
            Ea = exp(1i*n*Thetaa);
            Eb = exp(1i*(k*x1+n*Thetab));
            for i = 1:length(Thetaa)
                ua{i,1}(:,1) = u1*Ea(i);
                ua{i,1}(:,2) = u2*Ea(i);
                ua{i,1}(:,3) = u3*Ea(i);
            end
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                for i = 1:length(Thetab)
                    ub{i,1}(:,1) = u1(end)*Eb(i,:);
                    ub{i,1}(:,2) = u2(end)*Eb(i,:);
                    ub{i,1}(:,3) = u3(end)*Eb(i,:);
                end
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(3)*Hnz;
                df = U(3)*dHnz;
                u1 = 1i*k*f;
                u2 = 1i*n*f./rOuterFluid;
                for i = 1:length(Thetaa)
                    uOuterFluida{i,1}(:,1) = u1*Ea(i);
                    uOuterFluida{i,1}(:,2) = u2*Ea(i);
                    uOuterFluida{i,1}(:,3) = df*Ea(i);
                    ua{i,1} = [ua{i};uOuterFluida{i}];
                end
                for i = 1:length(Thetab)
                    ub{i,1}(:,1) = u1(end)*Eb(i,:);
                    ub{i,1}(:,2) = u2(end)*Eb(i,:);
                    ub{i,1}(:,3) = df(end)*Eb(i,:);
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
        E = exp(1i*k*x1); 
        if  strcmp(Plane,'r-z')
            for i = 1:length(x1)
                u{i,1}(:,1) = u1*E(i);
                u{i,1}(:,3) = u3*E(i);
            end
            if  ToggleOuterFluid && ShowHalfSpace
                f = U(2)*H0z;
                df = -U(2)*H1z;
                u1 = 1i*k*f;
                for i = 1:length(x1)
                    uOuterFluid{i,1}(:,1) = u1*E(i);
                    uOuterFluid{i,1}(:,3) = df*E(i);
                    u{i,1} = [u{i};uOuterFluid{i}];
                end
                r = [r;rOuterFluid];
            end
        else
            ua(:,1) = u1;
            ua(:,3) = u3;
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                ub(:,1) = u1(end)*E;
                ub(:,3) = u3(end)*E;
            else
                f = U(2)*H0z;
                df = -U(2)*H1z;
                uOuterFluida(:,1) = 1i*k*f;
                uOuterFluida(:,3) = df;
                ua = [ua;uOuterFluida];
                ub(:,1) = 1i*k*f(end)*E;
                ub(:,3) = df(end)*E;
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
            ua(:,2) = 1i*dg*exp(.5i*pi);
            ua(:,3) = 0;
            if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
                ub(:,2) = ua(end,2)*exp(1i*k*x1);
                ub(:,3) = 0;
            else
                ua = [ua;zeros(length(rOuterFluid),3)];
                ub(length(x1),3) = 0;
                r = [r;rOuterFluid];
            end
        end
    end
end
if  strcmp(Plane,'r-z')
    [X1,X3] = meshgrid(x1,r);
    if  ~contains(Mode,'T')
        Ratio = (r(end)-r(1))/x1(end);
        u1Max = max(abs(real(u{1}(:,1))));
        u3Max = max(abs(real(u{1}(:,3))));
        if  ToggleOuterFluid && ShowHalfSpace
            k = 1+HalfSpaces;
        else
            k = 1;
        end
        if  u1Max > u3Max
            Compensation = Gain*x1(end)/40/k/u1Max;
        else
            Compensation = Gain*(r(end)-r(1))/20/k/u3Max;
        end
        for i = 1:size(X1,2)
            X1Distorted(:,i) = X1(:,i)+real(u{i}(:,1))*Compensation;
            X3Distorted(:,i) = X3(:,i)+real(u{i}(:,3))*Compensation*Ratio;
        end
    else
        X1Distorted = X1;
        X3Distorted = X3;
    end
else
    X1a = repmat(r,[1,length(Thetaa)]);
    X2a = repmat(Thetaa',[length(r),1]);
    X3a = zeros(size(X1a));
    X1b = repmat(r(end),[length(x1),length(Thetab)]);
    X2b = repmat(Thetab',[length(x1),1]);
    X3b = repmat(x1',[1,length(Thetab)]);
    if  contains(Mode,'F')
        for i = 1:length(ua)
            uMax(i) = max(abs(real(ua{i})),[],'all');
        end
        uMax = max(uMax);
        Compensation = Gain*r(end)/20/uMax;
        for i = 1:size(X1a,2)
            X1Distorteda(:,i) = X1a(:,i)+real(ua{i}(:,3))*Compensation; % r
            X2Distorteda(:,i) = X2a(:,i)+real(ua{i}(:,2))./r*Compensation; % theta[rad] = b[m]/r[m]
            X3Distorteda(:,i) = real(ua{i}(:,1))*Compensation; % z
        end
        for i = 1:size(X1b,2)
            X1Distortedb(:,i) = X1b(:,i)+real(ub{i}(:,3))*Compensation;
            X2Distortedb(:,i) = X2b(:,i)+real(ub{i}(:,2))/r(end)*Compensation;
            X3Distortedb(:,i) = X3b(:,i)+real(ub{i}(:,1))*Compensation;
        end
    else
        uMax = max(abs(real(ua)),[],'all');
        Compensation = Gain*r(end)/20/uMax;
        X1Distorteda = repmat(X1a(:,1)+real(ua(:,3))*Compensation,[1,size(X1a,2)]);
        X2Distorteda = X2a(1,:)+real(ua(:,2))./r*Compensation;
        X3Distorteda = repmat(real(ua(:,1))*Compensation,[1,size(X1a,2)]);
        X1Distortedb = repmat(X1b(:,1)+real(ub(:,3))*Compensation,[1,size(X1b,2)]);
        X2Distortedb = X2b(1,:)+real(ub(:,2))/r(end)*Compensation;
        X3Distortedb = repmat(X3b(:,1)+real(ub(:,1))*Compensation,[1,size(X1b,2)]);
    end
    [X1a,X2a] = pol2cart(X2a,X1a);
    [X1b,X2b] = pol2cart(X2b,X1b);
    [X1Distorteda,X2Distorteda] = pol2cart(X2Distorteda,X1Distorteda);
    [X1Distortedb,X2Distortedb] = pol2cart(X2Distortedb,X1Distortedb);
end
f = figure('Name','Mode shape','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
ax = gca;
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
if  strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
    k = HalfSpaces*SamplesR;
    if  strcmp(OuterFluid.Name,InnerFluid.Name)
        OuterFluidColor = [0 .5 1];
        InnerFluidColor = [0 .5 1];
    else
        if  OuterFluid.Density*OuterFluid.Velocity > InnerFluid.Density*InnerFluid.Velocity
            OuterFluidColor = [0 .5 1];
            InnerFluidColor = [0 .7 1];
        else
            OuterFluidColor = [0 .7 1];
            InnerFluidColor = [0 .5 1];
        end
    end
end
if  strcmp(Plane,'r-z')
    axis off
    if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
        if  Undistorted
            line(X1(:,1:GridLine:end),X3(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(1:GridLine:end,:)',X3(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
        line(X1Distorted(:,1:GridLine:end),X3Distorted(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(1:GridLine:end,:)',X3Distorted(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
    elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
        if  Undistorted
            line(X1(1:SamplesR,1:GridLine:end),X3(1:SamplesR,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(1:GridLine:SamplesR,:)',X3(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(end-k:end,1:GridLine:end),X3(end-k:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(end-k:GridLine:end,:)',X3(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(SamplesR+1:end-k-1,1:GridLine:end),X3(SamplesR+1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(SamplesR+1:GridLine:end-k-1,:)',X3(SamplesR+1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
        line(X1Distorted(1:SamplesR,1:GridLine:end),X3Distorted(1:SamplesR,1:GridLine:end),'LineWidth',LineWidth,'Color',InnerFluidColor)
        line(X1Distorted(1:GridLine:SamplesR,:)',X3Distorted(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',InnerFluidColor)
        line(X1Distorted(end-k:end,1:GridLine:end),X3Distorted(end-k:end,1:GridLine:end),'LineWidth',LineWidth,'Color',OuterFluidColor)
        line(X1Distorted(end-k:GridLine:end,:)',X3Distorted(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',OuterFluidColor)
        line(X1Distorted(SamplesR+1:end-k-1,1:GridLine:end),X3Distorted(SamplesR+1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(SamplesR+1:GridLine:end-k-1,:)',X3Distorted(SamplesR+1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color','k')
    elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
        if  Undistorted
            line(X1(SamplesR+2:end,1:GridLine:end),X3(SamplesR+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(SamplesR+2:GridLine:end,:)',X3(SamplesR+2:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(1:SamplesR+1,1:GridLine:end),X3(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(1:GridLine:SamplesR+1,:)',X3(1:GridLine:SamplesR+1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
        line(X1Distorted(SamplesR+2:end,1:GridLine:end),X3Distorted(SamplesR+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(SamplesR+2:GridLine:end,:)',X3Distorted(SamplesR+2:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(1:SamplesR+1,1:GridLine:end),X3Distorted(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(1:GridLine:SamplesR+1,:)',X3Distorted(1:GridLine:SamplesR+1,:)','LineWidth',LineWidth,'Color','k')
    elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
        if  Undistorted
            line(X1(1:SamplesR,1:GridLine:end),X3(1:SamplesR,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(1:GridLine:SamplesR,:)',X3(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(SamplesR+1:end,1:GridLine:end),X3(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X1(SamplesR+1:GridLine:end,:)',X3(SamplesR+1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
        line(X1Distorted(1:SamplesR,1:GridLine:end),X3Distorted(1:SamplesR,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(1:GridLine:SamplesR,:)',X3Distorted(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(SamplesR+1:end,1:GridLine:end),X3Distorted(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(SamplesR+1:GridLine:end,:)',X3Distorted(SamplesR+1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
    end
    text(-.02,.5-.05*FontSizeAxesLabels/30,'$r$','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels,'rotation',90)
    text(.5-.155*FontSizeAxesLabels/30,-.01,'Propagation direction ($z$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels)
    ax.XLim = [-.075*x1(end) x1(end)+.075*x1(end)];
else
    rotate3d on
    axis equal tight off
    view(Azimuth,Elevation)
    if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
        if  Undistorted
            line(X3a(:,1:GridLine:end),X2a(:,1:GridLine:end),X1a(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(1:GridLine:end,:)',X2a(1:GridLine:end,:)',X1a(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3b(:,GridLine:GridLine:end),X2b(:,GridLine:GridLine:end),X1b(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3b(1:GridLine:end,:)',X2b(1:GridLine:end,:)',X1b(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
        line(X3Distorteda(:,1:GridLine:end),X2Distorteda(:,1:GridLine:end),X1Distorteda(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X3Distorteda(1:GridLine:end,:)',X2Distorteda(1:GridLine:end,:)',X1Distorteda(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
        line(X3Distortedb(:,GridLine:GridLine:end),X2Distortedb(:,GridLine:GridLine:end),X1Distortedb(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X3Distortedb(1:GridLine:end,:)',X2Distortedb(1:GridLine:end,:)',X1Distortedb(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
    elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
        if  Undistorted
            line(X3a(1:SamplesR+1,1:GridLine:end),X2a(1:SamplesR+1,1:GridLine:end),X1a(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(1:GridLine:SamplesR,:)',X2a(1:GridLine:SamplesR,:)',X1a(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(end-k-1:end,1:GridLine:end),X2a(end-k-1:end,1:GridLine:end),X1a(end-k-1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(end-k:GridLine:end,:)',X2a(end-k:GridLine:end,:)',X1a(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(SamplesR+1:end-k-1,1:GridLine:end),X2a(SamplesR+1:end-k-1,1:GridLine:end),X1a(SamplesR+1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(SamplesR+1:GridLine:end-k-1,:)',X2a(SamplesR+1:GridLine:end-k-1,:)',X1a(SamplesR+1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3b(:,GridLine:GridLine:end),X2b(:,GridLine:GridLine:end),X1b(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3b(1:GridLine:end,:)',X2b(1:GridLine:end,:)',X1b(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
        line(X3Distorteda(1:SamplesR+1,1:GridLine:end),X2Distorteda(1:SamplesR+1,1:GridLine:end),X1Distorteda(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',InnerFluidColor)
        line(X3Distorteda(1:GridLine:SamplesR,:)',X2Distorteda(1:GridLine:SamplesR,:)',X1Distorteda(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',InnerFluidColor)
        line(X3Distorteda(end-k-1:end,1:GridLine:end),X2Distorteda(end-k-1:end,1:GridLine:end),X1Distorteda(end-k-1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',OuterFluidColor)
        line(X3Distorteda(end-k:GridLine:end,:)',X2Distorteda(end-k:GridLine:end,:)',X1Distorteda(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',OuterFluidColor)
        line(X3Distorteda(SamplesR+1:end-k-1,1:GridLine:end),X2Distorteda(SamplesR+1:end-k-1,1:GridLine:end),X1Distorteda(SamplesR+1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X3Distorteda(SamplesR+1:GridLine:end-k-1,:)',X2Distorteda(SamplesR+1:GridLine:end-k-1,:)',X1Distorteda(SamplesR+1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color','k')
        line(X3Distortedb(:,GridLine:GridLine:end),X2Distortedb(:,GridLine:GridLine:end),X1Distortedb(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',OuterFluidColor)
        line(X3Distortedb(1:GridLine:end,:)',X2Distortedb(1:GridLine:end,:)',X1Distortedb(1:GridLine:end,:)','LineWidth',LineWidth,'Color',OuterFluidColor)
    elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
        if  Undistorted
            line(X3a(SamplesR+1:end,1:GridLine:end),X2a(SamplesR+1:end,1:GridLine:end),X1a(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(SamplesR+2:GridLine:end,:)',X2a(SamplesR+2:GridLine:end,:)',X1a(SamplesR+2:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(1:SamplesR+1,1:GridLine:end),X2a(1:SamplesR+1,1:GridLine:end),X1a(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(1:GridLine:SamplesR+1,:)',X2a(1:GridLine:SamplesR+1,:)',X1a(1:GridLine:SamplesR+1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3b(:,GridLine:GridLine:end),X2b(:,GridLine:GridLine:end),X1b(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3b(1:GridLine:end,:)',X2b(1:GridLine:end,:)',X1b(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
        line(X3Distorteda(SamplesR+1:end,1:GridLine:end),X2Distorteda(SamplesR+1:end,1:GridLine:end),X1Distorteda(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X3Distorteda(SamplesR+2:GridLine:end,:)',X2Distorteda(SamplesR+2:GridLine:end,:)',X1Distorteda(SamplesR+2:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1])
        line(X3Distorteda(1:SamplesR+1,1:GridLine:end),X2Distorteda(1:SamplesR+1,1:GridLine:end),X1Distorteda(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X3Distorteda(1:GridLine:SamplesR+1,:)',X2Distorteda(1:GridLine:SamplesR+1,:)',X1Distorteda(1:GridLine:SamplesR+1,:)','LineWidth',LineWidth,'Color','k')
        line(X3Distortedb(:,GridLine:GridLine:end),X2Distortedb(:,GridLine:GridLine:end),X1Distortedb(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X3Distortedb(1:GridLine:end,:)',X2Distortedb(1:GridLine:end,:)',X1Distortedb(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1])
    elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
        if  Undistorted
            line(X3a(1:SamplesR+1,1:GridLine:end),X2a(1:SamplesR+1,1:GridLine:end),X1a(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(1:GridLine:SamplesR,:)',X2a(1:GridLine:SamplesR,:)',X1a(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(SamplesR+1:end,1:GridLine:end),X2a(SamplesR+1:end,1:GridLine:end),X1a(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3a(SamplesR+1:GridLine:end,:)',X2a(SamplesR+1:GridLine:end,:)',X1a(SamplesR+1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3b(:,GridLine:GridLine:end),X2b(:,GridLine:GridLine:end),X1b(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(X3b(1:GridLine:end,:)',X2b(1:GridLine:end,:)',X1b(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
        line(X3Distorteda(1:SamplesR+1,1:GridLine:end),X2Distorteda(1:SamplesR+1,1:GridLine:end),X1Distorteda(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X3Distorteda(1:GridLine:SamplesR,:)',X2Distorteda(1:GridLine:SamplesR,:)',X1Distorteda(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[0 .5 1])
        line(X3Distorteda(SamplesR+1:end,1:GridLine:end),X2Distorteda(SamplesR+1:end,1:GridLine:end),X1Distorteda(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X3Distorteda(SamplesR+1:GridLine:end,:)',X2Distorteda(SamplesR+1:GridLine:end,:)',X1Distorteda(SamplesR+1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
        line(X3Distortedb(:,GridLine:GridLine:end),X2Distortedb(:,GridLine:GridLine:end),X1Distortedb(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X3Distortedb(1:GridLine:end,:)',X2Distortedb(1:GridLine:end,:)',X1Distortedb(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
    end
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
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,~)
    output_txt = {[]};
end
end