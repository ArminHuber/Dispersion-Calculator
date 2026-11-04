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
function [u,ua,ub,x1,rTotal,Thetaa,Thetab,p] = ModeShapeGridComputer_Rod_Pipe(Geometry,Plane,Material,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,F,L,T,Frequency,Length,Mode,Layers,Radius,SamplesZ,SamplesR,ShowHalfSpace,HalfSpaces,CycleDuration,FrameRate)
SamplesTheta = 4*SamplesR;
ThetaRange = [.3 -.8]*pi;

%#ok<*AGROW>
u=0;ua=0;ub=0;
Thetaa = (2*pi*(-.5:1/SamplesTheta:.5)); % for the cross section
Thetab = Thetaa(Thetaa < ThetaRange(1) & Thetaa > ThetaRange(2)); % for the visible outer cylinder surface
p = str2double(regexp(Mode,'\d*','match'));
n = p(1);
if  contains(Mode,'F')
    [PhaseVelocity,Attenuation,SignFluids,Direction,Frequency] = ExtractData(Geometry,F{n}{p(2)},Frequency,Mode);
elseif contains(Mode,'L')
    [PhaseVelocity,Attenuation,SignFluids,Direction,Frequency] = ExtractData(Geometry,L{p(2)},Frequency,Mode);
elseif contains(Mode,'T')
    [PhaseVelocity,Attenuation,SignFluids,Direction,Frequency] = ExtractData(Geometry,T{p(2)},Frequency,Mode);
end
for m = 1:Layers
    Lambda(m) = conj(Material{m}.Lambda_complex);
    Mu(m) = conj(Material{m}.Mu_complex);
    cL2(m) = Material{m}.LongitudinalVelocity_complex^2;
    cT2(m) = Material{m}.TransverseVelocity_complex^2;
end
if  strcmp(Geometry,'Rod')
    rFo = (Radius(end):Radius(end)/SamplesR:Radius(end)+HalfSpaces*Radius(end))';
elseif strcmp(Geometry,'Pipe')
    rFi = (Radius(1)/SamplesR:Radius(1)/SamplesR:Radius(1))';
    rFo = (Radius(end):Radius(1)/SamplesR:Radius(end)+HalfSpaces*Radius(1))';
    if  ~ToggleInnerFluid
        InnerFluid.Velocity = 1e-10;
        InnerFluid.Density = 1e-10;
    end
end
if  ~ToggleOuterFluid
    OuterFluid.Velocity = 1e-10;
    OuterFluid.Density = 1e-10;
end
Time = reshape(Direction/Frequency/1e3*(0:1/CycleDuration/FrameRate:1),1,1,[]);
x1 = Length*PhaseVelocity/Frequency/1e3*(0:1/SamplesZ:1);
rTotal = [];
AngularFrequency = Frequency*pi*2e3;
AngularFrequency2 = AngularFrequency^2;
k = Direction*(AngularFrequency/PhaseVelocity+1i*Attenuation);
k2 = k^2;
wt = AngularFrequency*Time;
if  strcmp(Plane,'r-z')
    E = exp(1i*(k*x1-wt));
else
    Ea = exp(1i*(n*Thetaa-wt));
    Eb = exp(1i*(k*x1'+n*Thetab-wt));
end
if  contains(Mode,'F')
    for m = 1:Layers
        kL2 = AngularFrequency2/cL2(m);
        kT2 = AngularFrequency2/cT2(m);
        x2 = kL2-k2;
        y2 = kT2-k2;
        x = sqrt(-x2);
        y = sqrt(-y2);
        if  m == 1 && strcmp(Geometry,'Rod')
            R = Radius(2);
            R2 = R^2;
            xR = x*R;
            yR = y*R;
            ZWnx = besseli(n,xR);
            ZWny = besseli(n,yR);
            ZWn1x = x*besseli(n+1,xR);
            ZWn1y = y*besseli(n+1,yR);
            dZWnx = n*ZWnx/R+ZWn1x;
            dZWny = n*ZWny/R+ZWn1y;
            dg1ZWn1y = -(n+1)*ZWn1y/R-y2*ZWny;
            d2ZWnx = (n*(n-1)/R2-x2)*ZWnx-ZWn1x/R;
            d2ZWny = (n*(n-1)/R2-y2)*ZWny-ZWn1y/R;
            L1(4,4) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2*ZWnx; % sigma33
            L1(4,5) = 2*Mu(m)*k*dg1ZWn1y;
            L1(4,6) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
            L1(5,4) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
            L1(5,5) = 1i*Mu(m)*k*((n+1)*ZWn1y/R-dg1ZWn1y);
            L1(5,6) = 1i*Mu(m)*(2*d2ZWny+y2*ZWny);
            L1(6,4) = 2i*Mu(m)*k*dZWnx; % sigma13
            L1(6,5) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2)*ZWn1y+n*dg1ZWn1y/R);
            L1(6,6) = 1i*Mu(m)*n*k*ZWny/R;
            A{m,3}(4,4) = 1i*k*ZWnx; % u1
            A{m,3}(4,5) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
            A{m,3}(5,4) = 1i*n*ZWnx/R; % u2
            A{m,3}(5,5) = -1i*k*ZWn1y;
            A{m,3}(5,6) = 1i*dZWny;
            A{m,3}(6,4) = dZWnx; % u3
            A{m,3}(6,5) = k*ZWn1y;
            A{m,3}(6,6) = n*ZWny/R;
            A{m,1}(4:6,4:6) = L1(4:6,4:6)/A{m,3}(4:6,4:6);
        else
            R = Radius(m:m+1);
            R2 = R.^2;
            xR = x*R;
            yR = y*R;
            ZWnx = [besseli(n,xR) besselk(n,xR)];
            ZWny = [besseli(n,yR) besselk(n,yR)];
            ZWn1x = x*[besseli(n+1,xR) besselk(n+1,xR)];
            ZWn1y = y*[besseli(n+1,yR) besselk(n+1,yR)];
            dZWnx = n*ZWnx./R+[ZWn1x(:,1) -ZWn1x(:,2)];
            dZWny = n*ZWny./R+[ZWn1y(:,1) -ZWn1y(:,2)];
            dg1ZWn1y = -(n+1)*ZWn1y./R+y2*[-ZWny(:,1) ZWny(:,2)];
            d2ZWnx = (n*(n-1)./R2-x2).*ZWnx+[-ZWn1x(:,1) ZWn1x(:,2)]./R;
            d2ZWny = (n*(n-1)./R2-y2).*ZWny+[-ZWn1y(:,1) ZWn1y(:,2)]./R;
            L1([1 4],[1 4]) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2*ZWnx; % sigma33
            L1([1 4],[2 5]) = 2*Mu(m)*k*dg1ZWn1y;
            L1([1 4],[3 6]) = 2*Mu(m)*n*(dZWny-ZWny./R)./R;
            L1([2 5],[1 4]) = 2i*Mu(m)*n*(dZWnx-ZWnx./R)./R; % sigma23
            L1([2 5],[2 5]) = 1i*Mu(m)*k*((n+1)*ZWn1y./R-dg1ZWn1y);
            L1([2 5],[3 6]) = 1i*Mu(m)*(2*d2ZWny+y2*ZWny);
            L1([3 6],[1 4]) = 2i*Mu(m)*k*dZWnx; % sigma13
            L1([3 6],[2 5]) = 1i*Mu(m)*((n*(n+1)./R2+k2-y2).*ZWn1y+n*dg1ZWn1y./R);
            L1([3 6],[3 6]) = 1i*Mu(m)*n*k*ZWny./R;
            A{m,3}([1 4],[1 4]) = 1i*k*ZWnx; % u1
            A{m,3}([1 4],[2 5]) = 1i*((n+1)*ZWn1y./R+dg1ZWn1y);
            A{m,3}([2 5],[1 4]) = 1i*n*ZWnx./R; % u2
            A{m,3}([2 5],[2 5]) = -1i*k*ZWn1y;
            A{m,3}([2 5],[3 6]) = 1i*dZWny;
            A{m,3}([3 6],[1 4]) = dZWnx; % u3
            A{m,3}([3 6],[2 5]) = k*ZWn1y;
            A{m,3}([3 6],[3 6]) = n*ZWny./R;
            A{m,1} = L1/A{m,3};
        end
    end
    A{1,2} = A{1};
    for m = 2:Layers
        M0 = A{m}(1:3,1:3)-A{m-1,2}(4:6,4:6);
        M1 = A{m-1,2}(1:3,4:6)/M0;
        M2 = A{m}(4:6,1:3)/M0;
        A{m,2} = [A{m-1,2}(1:3,1:3)+M1*A{m-1,2}(4:6,1:3) -M1*A{m}(1:3,4:6);M2*A{m-1,2}(4:6,1:3) A{m}(4:6,4:6)-M2*A{m}(1:3,4:6)];
    end
    zo = SignFluids(1)*sqrt(AngularFrequency2/OuterFluid.Velocity^2-k2);
    zRo = zo*Radius(end);
    Hnzo = besselh(n,zRo);
    DFo = -OuterFluid.Density*AngularFrequency2*Hnzo;
    if  strcmp(Geometry,'Rod')
        A{end,2} = inv(A{end,2}(4:6,4:6));
        uInterfaces(:,Layers+1) = A{end,2}(1:3,1)*DFo;
        for m = Layers:-1:2
            uInterfaces(:,m) = (A{m-1,2}(4:6,4:6)-A{m}(1:3,1:3))\A{m}(1:3,4:6)*uInterfaces(:,m+1);
        end
    elseif strcmp(Geometry,'Pipe')
        zi = SignFluids(2)*sqrt(AngularFrequency2/InnerFluid.Velocity^2-k2);
        zRi = zi*Radius(1);
        if  Sink
            Znzi = besselh(n,2,zRi);
            Zn1zi = zi*besselh(n+1,2,zRi);
        else
            Znzi = besselj(n,zRi);
            Zn1zi = zi*besselj(n+1,zRi);
        end
        WFi = Zn1zi-n*Znzi/Radius(1);
        DFi = -InnerFluid.Density*AngularFrequency2*Znzi;
        A{end,2} = inv(A{end,2});
        UFi = -A{end,2}(3,4)*DFo/(WFi+A{end,2}(3,1)*DFi);
        uInterfaces = A{end,2}(1:3,1)*DFi*UFi+A{end,2}(1:3,4)*DFo;
        uInterfaces(:,Layers+1) = A{end,2}(4:6,1)*DFi*UFi+A{end,2}(4:6,4)*DFo;
        for m = Layers:-1:2
            M0 = A{m}(1:3,1:3)-A{m-1,2}(4:6,4:6);
            uInterfaces(:,m) = M0\A{m-1,2}(4:6,1:3)*uInterfaces(:,1)-M0\A{m}(1:3,4:6)*uInterfaces(:,m+1);
        end
    end
    for m = 1:Layers
        if  m == 1 && strcmp(Geometry,'Rod')
            r = (0:Radius(2)/SamplesR:Radius(2))';
            r(1) = 1e-10;
        else
            r = (Radius(m):(Radius(m+1)-Radius(m))/SamplesR:Radius(m+1))';
        end
        rTotal = [rTotal;r];
        s = (m-1)*length(r)+1:m*length(r);
        U = A{m,3}\[uInterfaces(:,m);uInterfaces(:,m+1)];
        kL2 = AngularFrequency2/cL2(m);
        kT2 = AngularFrequency2/cT2(m);
        x2 = kL2-k2;
        y2 = kT2-k2;
        x = sqrt(-x2);
        y = sqrt(-y2);
        xr = x*r;
        yr = y*r;
        Znx = besseli(n,xr);
        Zny = besseli(n,yr);
        Zn1x = x*besseli(n+1,xr);
        Zn1y = y*besseli(n+1,yr);
        dZnx = n*Znx./r+Zn1x;
        dZny = n*Zny./r+Zn1y;
        dg1Zn1y = -(n+1)*Zn1y./r-y2*Zny;
        if  m == 1 && strcmp(Geometry,'Rod')
            f = U(4)*Znx; % L_in
            df = U(4)*dZnx;
            g1 = U(5)*Zn1y; % SV_in
            dg1 = U(5)*dg1Zn1y;
            g3 = U(6)*Zny; % SH_in
            dg3 = U(6)*dZny;
        else
            Wnx = besselk(n,xr);
            Wny = besselk(n,yr);
            Wn1x = x*besselk(n+1,xr);
            Wn1y = y*besselk(n+1,yr);
            dWnx = n*Wnx./r-Wn1x;
            dWny = n*Wny./r-Wn1y;
            dg1Wn1y = -(n+1)*Wn1y./r+y2*Wny;
            f = U(1)*Znx+U(4)*Wnx; % L_in + L_out
            df = U(1)*dZnx+U(4)*dWnx;
            g1 = U(2)*Zn1y+U(5)*Wn1y; % SV_in + SV_out
            dg1 = U(2)*dg1Zn1y+U(5)*dg1Wn1y;
            g3 = U(3)*Zny+U(6)*Wny; % SH_in + SH_out
            dg3 = U(3)*dZny+U(6)*dWny;
        end
        u1(s,1) = 1i*(k*f+(n+1)*g1./r+dg1); % uz
        u2(s,1) = 1i*(n*f./r-k*g1+dg3); % uq
        u3(s,1) = df+k*g1+n*g3./r; % ur
    end
    if  Layers > 1
        [rTotal,ia] = unique(rTotal);
        u1 = u1(ia);
        u2 = u2(ia);
        u3 = u3(ia);
    end
    if  strcmp(Plane,'r-z')
        u = real(u1.*E); % uz
        u(:,:,:,3) = real(u3.*E); % ur
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            zri = zi*rFi;
            if  Sink
                Znzi = besselh(n,2,zri);
                Zn1zi = zi*besselh(n+1,2,zri);
            else
                Znzi = besselj(n,zri);
                Zn1zi = zi*besselj(n+1,zri);
            end
            f = UFi*Znzi;
            df = UFi*(n*Znzi./rFi-Zn1zi);
            uFi(:,:,:,1) = real(1i*k*f.*E);
            uFi(:,:,:,3) = real(df.*E);
            rTotal = [rFi;rTotal];
            u = [uFi;u];
        end
        if  ToggleOuterFluid && ShowHalfSpace
            zro = zo*rFo;
            Hnzo = besselh(n,zro);
            Hn1zo = zo*besselh(n+1,zro);
            f = Hnzo;
            df = n*Hnzo./rFo-Hn1zo;
            uFo(:,:,:,1) = real(1i*k*f.*E);
            uFo(:,:,:,3) = real(df.*E);
            rTotal = [rTotal;rFo];
            u = [u;uFo];
        end
    else
        ua = real(u1.*Ea); % uz
        ua(:,:,:,2) = real(u2.*Ea); % uq
        ua(:,:,:,3) = real(u3.*Ea); % ur
        if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
            ub = real(u1(end).*Eb);
            ub(:,:,:,2) = real(u2(end).*Eb);
            ub(:,:,:,3) = real(u3(end).*Eb);
        end
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            zri = zi*rFi;
            if  Sink
                Znzi = besselh(n,2,zri);
                Zn1zi = zi*besselh(n+1,2,zri);
            else
                Znzi = besselj(n,zri);
                Zn1zi = zi*besselj(n+1,zri);
            end
            f = UFi*Znzi;
            df = UFi*(n*Znzi./rFi-Zn1zi);
            uFi(:,:,:,1) = real(1i*k*f.*Ea);
            uFi(:,:,:,2) = real(1i*n*f./rFi.*Ea);
            uFi(:,:,:,3) = real(df.*Ea);
            rTotal = [rFi;rTotal];
            ua = [uFi;ua];
        end
        if  ToggleOuterFluid && ShowHalfSpace
            zro = zo*rFo;
            Hnzo = besselh(n,zro);
            Hn1zo = zo*besselh(n+1,zro);
            f = Hnzo;
            df = n*Hnzo./rFo-Hn1zo;
            u1 = 1i*k*f;
            u2 = 1i*n*f./rFo;
            uFo(:,:,:,1) = real(u1.*Ea);
            uFo(:,:,:,2) = real(u2.*Ea);
            uFo(:,:,:,3) = real(df.*Ea);
            ub = real(u1(end).*Eb);
            ub(:,:,:,2) = real(u2(end).*Eb);
            ub(:,:,:,3) = real(df(end).*Eb);
            rTotal = [rTotal;rFo];
            ua = [ua;uFo];
        end
    end
elseif contains(Mode,'L')
    for m = 1:Layers
        kL2 = AngularFrequency2/cL2(m);
        kT2 = AngularFrequency2/cT2(m);
        x2 = kL2-k2;
        y2 = kT2-k2;
        x = sqrt(-x2);
        y = sqrt(-y2);
        if  m == 1 && strcmp(Geometry,'Rod')
            R = Radius(2);
            xR = x*R;
            yR = y*R;
            Z0x = besseli(0,xR);
            Z0y = besseli(0,yR);
            Z1x = x*besseli(1,xR);
            Z1y = y*besseli(1,yR);
            d2Z0y = -y2*Z0y-Z1y/R;
            d2Z0x = -x2*Z0x-Z1x/R;
            L1(3,3) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2*Z0x; % sigma33
            L1(3,4) = 2*Mu(m)*k*d2Z0y;
            L1(4,3) = 2i*Mu(m)*k*Z1x; % sigma13
            L1(4,4) = 1i*Mu(m)*(k2-y2)*Z1y;
            A{m,3}(3,3) = 1i*k*Z0x; % u1
            A{m,3}(3,4) = 1i*(Z1y/R+d2Z0y);
            A{m,3}(4,3) = Z1x; % u3
            A{m,3}(4,4) = k*Z1y;
            A{m,1}(3:4,3:4) = L1(3:4,3:4)/A{m,3}(3:4,3:4);
        else
            R = Radius(m:m+1);
            xR = x*R;
            yR = y*R;
            ZW0x = [besseli(0,xR) besselk(0,xR)];
            ZW0y = [besseli(0,yR) besselk(0,yR)];
            ZW1x = x*[besseli(1,xR) besselk(1,xR)];
            ZW1y = y*[besseli(1,yR) besselk(1,yR)];
            dZW0x = [ZW1x(:,1) -ZW1x(:,2)];
            dg1ZW1y = y2*[-ZW0y(:,1) ZW0y(:,2)]-ZW1y./R;
            d2ZW0x = -x2*ZW0x-dZW0x./R;
            L1([1 3],[1 3]) = 2*Mu(m)*d2ZW0x-Lambda(m)*kL2*ZW0x; % sigma33
            L1([1 3],[2 4]) = 2*Mu(m)*k*dg1ZW1y;
            L1([2 4],[1 3]) = 2i*Mu(m)*k*dZW0x; % sigma13
            L1([2 4],[2 4]) = 1i*Mu(m)*(k2-y2)*ZW1y;
            A{m,3}([1 3],[1 3]) = 1i*k*ZW0x; % u1
            A{m,3}([1 3],[2 4]) = 1i*(ZW1y./R+dg1ZW1y);
            A{m,3}([2 4],[1 3]) = dZW0x; % u3
            A{m,3}([2 4],[2 4]) = k*ZW1y;
            A{m,1} = L1/A{m,3};
        end
    end
    A{1,2} = A{1};
    for m = 2:Layers
        M0 = A{m}(1:2,1:2)-A{m-1,2}(3:4,3:4);
        M1 = A{m-1,2}(1:2,3:4)/M0;
        M2 = A{m}(3:4,1:2)/M0;
        A{m,2} = [A{m-1,2}(1:2,1:2)+M1*A{m-1,2}(3:4,1:2) -M1*A{m}(1:2,3:4);M2*A{m-1,2}(3:4,1:2) A{m}(3:4,3:4)-M2*A{m}(1:2,3:4)];
    end
    zo = SignFluids(1)*sqrt(AngularFrequency2/OuterFluid.Velocity^2-k2);
    zRo = zo*Radius(end);
    H0zo = besselh(0,zRo);
    DFo = -OuterFluid.Density*AngularFrequency2*H0zo;
    if  strcmp(Geometry,'Rod')
        A{end,2} = inv(A{end,2}(3:4,3:4));
        uInterfaces(:,Layers+1) = A{end,2}(1:2,1)*DFo;
        for m = Layers:-1:2
            uInterfaces(:,m) = (A{m-1,2}(3:4,3:4)-A{m}(1:2,1:2))\A{m}(1:2,3:4)*uInterfaces(:,m+1);
        end
    elseif strcmp(Geometry,'Pipe')
        zi = SignFluids(2)*sqrt(AngularFrequency2/InnerFluid.Velocity^2-k2);
        zRi = zi*Radius(1);
        if  Sink
            Z0zi = besselh(0,2,zRi);
            Z1zi = zi*besselh(1,2,zRi);
        else
            Z0zi = besselj(0,zRi);
            Z1zi = zi*besselj(1,zRi);
        end
        WFi = Z1zi;
        DFi = -InnerFluid.Density*AngularFrequency2*Z0zi;
        A{end,2} = inv(A{end,2});
        UFi = -A{end,2}(2,3)*DFo/(WFi+A{end,2}(2,1)*DFi);
        uInterfaces = A{end,2}(1:2,1)*DFi*UFi+A{end,2}(1:2,3)*DFo;
        uInterfaces(:,Layers+1) = A{end,2}(3:4,1)*DFi*UFi+A{end,2}(3:4,3)*DFo;
        for m = Layers:-1:2
            M0 = A{m}(1:2,1:2)-A{m-1,2}(3:4,3:4);
            uInterfaces(:,m) = M0\A{m-1,2}(3:4,1:2)*uInterfaces(:,1)-M0\A{m}(1:2,3:4)*uInterfaces(:,m+1);
        end
    end
    for m = 1:Layers
        if  m == 1 && strcmp(Geometry,'Rod')
            r = (0:Radius(2)/SamplesR:Radius(2))';
            r(1) = 1e-10;
        else
            r = (Radius(m):(Radius(m+1)-Radius(m))/SamplesR:Radius(m+1))';
        end
        rTotal = [rTotal;r];
        s = (m-1)*length(r)+1:m*length(r);
        U = A{m,3}\[uInterfaces(:,m);uInterfaces(:,m+1)];
        kL2 = AngularFrequency2/cL2(m);
        kT2 = AngularFrequency2/cT2(m);
        x2 = kL2-k2;
        y2 = kT2-k2;
        x = sqrt(-x2);
        y = sqrt(-y2);
        xr = x*r;
        yr = y*r;
        Z0x = besseli(0,xr);
        Z0y = besseli(0,yr);
        Z1x = x*besseli(1,xr);
        Z1y = y*besseli(1,yr);
        dg1Z1y = -Z1y./r-y2*Z0y;
        if  m == 1 && strcmp(Geometry,'Rod')
            f = U(3)*Z0x; % L_in
            df = U(3)*Z1x;
            g = U(4)*Z1y; % SV_in
            dg = U(4)*dg1Z1y;
        else
            Wnx = besselk(n,xr);
            Wny = besselk(n,yr);
            Wn1x = x*besselk(n+1,xr);
            Wn1y = y*besselk(n+1,yr);
            dg1Wn1y = -Wn1y./r+y2*Wny;
            f = U(1)*Z0x+U(3)*Wnx; % L_in + L_out
            df = U(1)*Z1x-U(3)*Wn1x;
            g = U(2)*Z1y+U(4)*Wn1y; % SV_in + SV_out
            dg = U(2)*dg1Z1y+U(4)*dg1Wn1y;
        end
        u1(s,1) = 1i*(k*f+g./r+dg); % uz
        u3(s,1) = df+k*g; % ur
    end
    if  Layers > 1
        [rTotal,ia] = unique(rTotal);
        u1 = u1(ia);
        u3 = u3(ia);
    end
    if  strcmp(Plane,'r-z')
        u = real(u1.*E); % uz
        u(:,:,:,3) = real(u3.*E); % ur
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            zri = zi*rFi;
            if  Sink
                Z0zi = besselh(0,2,zri);
                Z1zi = zi*besselh(1,2,zri);
            else
                Z0zi = besselj(0,zri);
                Z1zi = zi*besselj(1,zri);
            end
            f = UFi*Z0zi;
            df = -UFi*Z1zi;
            uFi(:,:,:,1) = real(1i*k*f.*E);
            uFi(:,:,:,3) = real(df.*E);
            rTotal = [rFi;rTotal];
            u = [uFi;u];
        end
        if  ToggleOuterFluid && ShowHalfSpace
            zro = zo*rFo;
            Z0zo = besselh(0,zro);
            Z1zo = zo*besselh(1,zro);
            f = Z0zo;
            df = -Z1zo;
            uFo(:,:,:,1) = real(1i*k*f.*E);
            uFo(:,:,:,3) = real(df.*E);
            rTotal = [rTotal;rFo];
            u = [u;uFo];
        end
    else
        ua = real(u1.*Ea); % uz
        ua(:,:,:,3) = real(u3.*Ea); % ur
        if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
            ub = real(u1(end).*Eb);
            ub(:,:,:,3) = real(u3(end).*Eb);
        end
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            zri = zi*rFi;
            if  Sink
                Z0zi = besselh(0,2,zri);
                Z1zi = zi*besselh(1,2,zri);
            else
                Z0zi = besselj(0,zri);
                Z1zi = zi*besselj(1,zri);
            end
            f = UFi*Z0zi;
            df = -UFi*Z1zi;
            uFi(:,:,:,1) = real(1i*k*f.*Ea);
            uFi(:,:,:,3) = real(df.*Ea);
            rTotal = [rFi;rTotal];
            ua = [uFi;ua];
        end
        if  ToggleOuterFluid && ShowHalfSpace
            zro = zo*rFo;
            Z0zo = besselh(0,zro);
            Z1zo = zo*besselh(1,zro);
            f = Z0zo;
            df = -Z1zo;
            u1 = 1i*k*f;
            uFo(:,:,:,1) = real(u1.*Ea);
            uFo(:,:,:,3) = real(df.*Ea);
            ub = real(u1(end).*Eb);
            ub(:,:,:,3) = real(df(end).*Eb);
            rTotal = [rTotal;rFo];
            ua = [ua;uFo];
        end
    end
elseif contains(Mode,'T')
    for m = 1:Layers
        kT2 = AngularFrequency2/cT2(m);
        y2 = kT2-k2;
        if  y2 == 0
            y2 = 1e-10;
        end
        y = sqrt(-y2);
        if  m == 1 && strcmp(Geometry,'Rod')
            R = Radius(2);
            yR = y*R;
            Z0y = besseli(0,yR);
            Z1y = y*besseli(1,yR);
            A{m,3} = 1i*Z1y; % u2;
            A{m,1}(2,2) = -1i*Mu(m)*(y2*Z0y+2*Z1y/R)/A{m,3};
        else
            R = Radius(m:m+1);
            yR = y*R;
            ZW0y = [besseli(0,yR) besselk(0,yR)];
            ZW1y = y*[besseli(1,yR) besselk(1,yR)];
            dZW0y = [ZW1y(:,1) -ZW1y(:,2)];
            d2ZW0y = -y2*ZW0y-dZW0y./R;
            A{m,3} = 1i*dZW0y;
            A{m,1} = 1i*Mu(m)*(2*d2ZW0y+y2*ZW0y)/A{m,3};
        end
    end
    A{1,2} = A{1};
    for m = 2:Layers
        M0 = A{m}(1,1)-A{m-1,2}(2,2);
        M1 = A{m-1,2}(1,2)/M0;
        M2 = A{m}(2,1)/M0;
        A{m,2} = [A{m-1,2}(1,1)+M1*A{m-1,2}(2,1) -M1*A{m}(1,2);M2*A{m-1,2}(2,1) A{m}(2,2)-M2*A{m}(1,2)];
    end
    if  strcmp(Geometry,'Rod')
        uInterfaces(:,Layers+1) = 1;
        for m = Layers:-1:2
            uInterfaces(:,m) = (A{m-1,2}(2,2)-A{m}(1,1))\A{m}(1,2)*uInterfaces(:,m+1);
        end
    elseif strcmp(Geometry,'Pipe')
        In0 = A{1,3}(1,1);
        Out0 = A{1,3}(1,2);
        OutN = A{end,3}(2,2);
        Z1 = [A{end,2}(1,1)*In0 -A{end,2}(1,2)*OutN;A{end,2}(2,1)*In0 A{end,2}(2,2)*OutN]; % changed sign in Z1(1,1)!
        Z2 = [A{end,2}(1,1)*Out0;-A{end,2}(2,1)*Out0];
        RT = Z1\Z2;
        uInterfaces = Out0+In0*RT(1);
        uInterfaces(:,Layers+1) = OutN*RT(2);
        for m = Layers:-1:2
            M0 = A{m}(1,1)-A{m-1,2}(2,2);
            uInterfaces(:,m) = M0\A{m-1,2}(2,1)*uInterfaces(:,1)-M0\A{m}(1,2)*uInterfaces(:,m+1);
        end
    end
    for m = 1:Layers
        if  m == 1 && strcmp(Geometry,'Rod')
            r = (0:Radius(2)/SamplesR:Radius(2))';
            r(1) = 1e-10;
        else
            r = (Radius(m):(Radius(m+1)-Radius(m))/SamplesR:Radius(m+1))';
        end
        rTotal = [rTotal;r];
        s = (m-1)*length(r)+1:m*length(r);
        U = A{m,3}\[uInterfaces(:,m);uInterfaces(:,m+1)];
        kT2 = AngularFrequency2/cT2(m);
        y2 = kT2-k2;
        if  y2 == 0
            y2 = 1e-10;
        end
        y = sqrt(-y2);
        yr = y*r;
        Z1y = y*besseli(1,yr);
        if  m == 1 && strcmp(Geometry,'Rod')
            dg = U(2)*Z1y; % SH_in
        else
            W1y = y*besselk(1,yr);
            dg = U(1)*Z1y-U(2)*W1y; % SH_in + SH_out
        end
        u2(s,1) = 1i*dg; % uq
    end
    if  Layers > 1
        [rTotal,ia] = unique(rTotal);
        u2 = u2(ia);
    end
    if  strcmp(Plane,'r-z')
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            rTotal = [rFi;rTotal];
        end
        if  ToggleOuterFluid && ShowHalfSpace
            rTotal = [rTotal;rFo];
        end
    else
        ua = zeros([length(rTotal) size(Ea,2:3) 3]);
        ua(:,:,:,2) = real(u2.*Ea); % uq
        if  ~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)
            ub = zeros([size(Eb,1:3) 3]);
            ub(:,:,:,2) = real(u2(end).*Eb);
        end
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            uFi = zeros([length(rFi) size(Ea,2:3) 3]);
            rTotal = [rFi;rTotal];
            ua = [uFi;ua];
        end
        if  ToggleOuterFluid && ShowHalfSpace
            uFo = zeros([length(rFo) size(Ea,2:3) 3]);
            ub = zeros([size(Eb,1:3) 3]);
            rTotal = [rTotal;rFo];
            ua = [ua;uFo];
        end
    end
end
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