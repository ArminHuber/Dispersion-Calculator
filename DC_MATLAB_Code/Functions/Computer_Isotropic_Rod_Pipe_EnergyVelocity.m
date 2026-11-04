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
function [F,L,T] = Computer_Isotropic_Rod_Pipe_EnergyVelocity(Geometry,F,L,T,Material,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Layers,Radius,SamplesR)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
ModeTotal = 0;
if  ~isempty(F{1})
    for n = 1:length(F)
        ModeTotal = ModeTotal+length(F{n});
    end
end
if  ~isempty(L{1})
    ModeTotal = ModeTotal+length(L);
end
if  ~isempty(T{1})
    ModeTotal = ModeTotal+length(T);
end
if  ModeTotal == 0
    return
end
SamplesR = ceil(SamplesR/Layers); % samples per layer
Height = Layers*(SamplesR+1);
if  ~ToggleOuterFluid
    OuterFluid.Velocity = 1e-10;
    OuterFluid.Density = 1e-10;
end
if  ~ToggleInnerFluid
    InnerFluid.Velocity = 1e-10;
    InnerFluid.Density = 1e-10;
end
cFo2 = OuterFluid.Velocity^2;
cFi2 = InnerFluid.Velocity^2;
for m = 1:Layers
    Lambda(m) = conj(Material{m}.Lambda_complex);
    Mu(m) = conj(Material{m}.Mu_complex);
    cL2(m) = Material{m}.LongitudinalVelocity_complex^2;
    cT2(m) = Material{m}.TransverseVelocity_complex^2;
end
if  Sink
    rFi = (Radius(1)/SamplesR:Radius(1)/SamplesR:Radius(1))';
else
    rFi = (0:Radius(1)/SamplesR:Radius(1))';
    rFi(1) = 1e-10;
end
rFi2 = rFi.^2;
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(F{1})
    for n = 1:length(F)
        [F{n},Counter] = Computer_Isotropic_Multi_Rod_Pipe_EnergyVelocity_Core(Geometry,F{n},'F',n,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFo2,cFi2,ToggleInnerFluid,OuterFluid,InnerFluid,Sink,Layers,Radius,rFi,rFi2,SamplesR,Height);
        if  Stop
            return
        end
    end
end
if  ~isempty(L{1})
    [L,Counter] = Computer_Isotropic_Multi_Rod_Pipe_EnergyVelocity_Core(Geometry,L,'L',0,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFo2,cFi2,ToggleInnerFluid,OuterFluid,InnerFluid,Sink,Layers,Radius,rFi,rFi2,SamplesR,Height);
    if  Stop
        return
    end
end
if  ~isempty(T{1})
    [T,~] = Computer_Isotropic_Multi_Rod_Pipe_EnergyVelocity_Core(Geometry,T,'T',0,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFo2,cFi2,ToggleInnerFluid,OuterFluid,InnerFluid,Sink,Layers,Radius,rFi,rFi2,SamplesR,Height);
    if  Stop
        return
    end
end
close(h)
end
function [X,Counter] = Computer_Isotropic_Multi_Rod_Pipe_EnergyVelocity_Core(Geometry,X,ModeType,n,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFo2,cFi2,ToggleInnerFluid,OuterFluid,InnerFluid,Sink,Layers,Radius,rFi,rFi2,SamplesR,Height)
    global Stop 
    Stop = 0;
    for p = 1:length(X)
        PhaseVelocity = reshape(X{p}(:,4),1,1,[])*1e3;
        Attenuation = reshape(X{p}(:,7),1,1,[]); % Np/m
        AngularFrequency = reshape(X{p}(:,1),1,1,[])*pi*2e3;
        AngularFrequency2 = AngularFrequency.^2;
        k = AngularFrequency./PhaseVelocity+1i*Attenuation;
        k2 = k.^2;
        Length = length(k);
        KineticEnergyDensity = zeros(Height,1,Length);
        A = cell(0);
        rTotal = [];
        if  strcmp(ModeType,'F')
            L1 = zeros(6,6,Length);
            uInterfaces = zeros(3,Layers+1,Length);
            u = zeros(Height,3,Length);
            v = zeros(Height,3,Length);
            epsilon = zeros(Height,6,Length);
            sigma = zeros(Height,6,Length);
            uFi = zeros(length(rFi),3,Length);
            epsilonFi = zeros(length(rFi),3,Length);
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
                    ZWn1x = x.*besseli(n+1,xR);
                    ZWn1y = y.*besseli(n+1,yR);
                    dZWnx = n*ZWnx/R+ZWn1x;
                    dZWny = n*ZWny/R+ZWn1y;
                    dg1ZWn1y = -(n+1)*ZWn1y/R-y2.*ZWny;
                    d2ZWnx = (n*(n-1)/R2-x2).*ZWnx-ZWn1x/R;
                    d2ZWny = (n*(n-1)/R2-y2).*ZWny-ZWn1y/R;
                    L1(4,4,:) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2.*ZWnx; % sigma33
                    L1(4,5,:) = 2*Mu(m)*k.*dg1ZWn1y;
                    L1(4,6,:) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
                    L1(5,4,:) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
                    L1(5,5,:) = 1i*Mu(m)*k.*((n+1)*ZWn1y/R-dg1ZWn1y);
                    L1(5,6,:) = 1i*Mu(m)*(2*d2ZWny+y2.*ZWny);
                    L1(6,4,:) = 2i*Mu(m)*k.*dZWnx; % sigma13
                    L1(6,5,:) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2).*ZWn1y+n*dg1ZWn1y/R);
                    L1(6,6,:) = 1i*Mu(m)*n*k.*ZWny/R;
                    A{m,3}(4,4,:) = 1i*k.*ZWnx; % u1
                    A{m,3}(4,5,:) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
                    A{m,3}(5,4,:) = 1i*n*ZWnx/R; % u2
                    A{m,3}(5,5,:) = -1i*k.*ZWn1y;
                    A{m,3}(5,6,:) = 1i*dZWny;
                    A{m,3}(6,4,:) = dZWnx; % u3
                    A{m,3}(6,5,:) = k.*ZWn1y;
                    A{m,3}(6,6,:) = n*ZWny/R;
                    A{m,1}(4:6,4:6,:) = pagemrdivide(L1(4:6,4:6,:),A{m,3}(4:6,4:6,:));
                else
                    R = Radius(m:m+1);
                    R2 = R.^2;
                    xR = x.*R;
                    yR = y.*R;
                    ZWnx = [besseli(n,xR) besselk(n,xR)];
                    ZWny = [besseli(n,yR) besselk(n,yR)];
                    ZWn1x = x.*[besseli(n+1,xR) besselk(n+1,xR)];
                    ZWn1y = y.*[besseli(n+1,yR) besselk(n+1,yR)];
                    dZWnx = n*ZWnx./R+[ZWn1x(:,1,:) -ZWn1x(:,2,:)];
                    dZWny = n*ZWny./R+[ZWn1y(:,1,:) -ZWn1y(:,2,:)];
                    dg1ZWn1y = -(n+1)*ZWn1y./R+y2.*[-ZWny(:,1,:) ZWny(:,2,:)];
                    d2ZWnx = (n*(n-1)./R2-x2).*ZWnx+[-ZWn1x(:,1,:) ZWn1x(:,2,:)]./R;
                    d2ZWny = (n*(n-1)./R2-y2).*ZWny+[-ZWn1y(:,1,:) ZWn1y(:,2,:)]./R;
                    L1([1 4],[1 4],:) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2.*ZWnx; % sigma33
                    L1([1 4],[2 5],:) = 2*Mu(m)*k.*dg1ZWn1y;
                    L1([1 4],[3 6],:) = 2*Mu(m)*n*(dZWny-ZWny./R)./R;
                    L1([2 5],[1 4],:) = 2i*Mu(m)*n*(dZWnx-ZWnx./R)./R; % sigma23
                    L1([2 5],[2 5],:) = 1i*Mu(m)*k.*((n+1)*ZWn1y./R-dg1ZWn1y);
                    L1([2 5],[3 6],:) = 1i*Mu(m)*(2*d2ZWny+y2.*ZWny);
                    L1([3 6],[1 4],:) = 2i*Mu(m)*k.*dZWnx; % sigma13
                    L1([3 6],[2 5],:) = 1i*Mu(m)*((n*(n+1)./R2+k2-y2).*ZWn1y+n*dg1ZWn1y./R);
                    L1([3 6],[3 6],:) = 1i*Mu(m)*n*k.*ZWny./R;
                    A{m,3}([1 4],[1 4],:) = 1i*k.*ZWnx; % u1
                    A{m,3}([1 4],[2 5],:) = 1i*((n+1)*ZWn1y./R+dg1ZWn1y);
                    A{m,3}([2 5],[1 4],:) = 1i*n*ZWnx./R; % u2
                    A{m,3}([2 5],[2 5],:) = -1i*k.*ZWn1y;
                    A{m,3}([2 5],[3 6],:) = 1i*dZWny;
                    A{m,3}([3 6],[1 4],:) = dZWnx; % u3
                    A{m,3}([3 6],[2 5],:) = k.*ZWn1y;
                    A{m,3}([3 6],[3 6],:) = n*ZWny./R;
                    A{m,1} = pagemrdivide(L1,A{m,3});
                end
            end
            A{1,2} = A{1};
            for m = 2:Layers
                M0 = A{m}(1:3,1:3,:)-A{m-1,2}(4:6,4:6,:);
                M1 = pagemrdivide(A{m-1,2}(1:3,4:6,:),M0);
                M2 = pagemrdivide(A{m}(4:6,1:3,:),M0);
                A{m,2} = [A{m-1,2}(1:3,1:3,:)+pagemtimes(M1,A{m-1,2}(4:6,1:3,:)) -pagemtimes(M1,A{m}(1:3,4:6,:));pagemtimes(M2,A{m-1,2}(4:6,1:3,:)) A{m}(4:6,4:6,:)-pagemtimes(M2,A{m}(1:3,4:6,:))];
            end
            zo = reshape(X{p}(:,8),1,1,[]).*sqrt(AngularFrequency2/cFo2-k2);
            zRo = zo*Radius(end);
            Hnzo = besselh(n,zRo);
            DFo = -OuterFluid.Density*AngularFrequency2.*Hnzo;
            if  strcmp(Geometry,'Rod')
                A{end,2} = pageinv(A{end,2}(4:6,4:6,:));
                uInterfaces(:,Layers+1,:) = A{end,2}(1:3,1,:).*DFo;
                for m = Layers:-1:2
                    uInterfaces(:,m,:) = pagemtimes(pagemldivide(A{m-1,2}(4:6,4:6,:)-A{m}(1:3,1:3,:),A{m}(1:3,4:6,:)),uInterfaces(:,m+1,:));
                end
            elseif strcmp(Geometry,'Pipe')
                zi = reshape(X{p}(:,9),1,1,[]).*sqrt(AngularFrequency2/cFi2-k2);
                zRi = zi*Radius(1);
                if  Sink
                    Znzi = besselh(n,2,zRi);
                    Zn1zi = zi.*besselh(n+1,2,zRi);
                else
                    Znzi = besselj(n,zRi);
                    Zn1zi = zi.*besselj(n+1,zRi);
                end
                WFi = Zn1zi-n*Znzi./Radius(1);
                DFi = -InnerFluid.Density*AngularFrequency2.*Znzi;
                A{end,2} = pageinv(A{end,2});
                UFi = -A{end,2}(3,4,:).*DFo./(WFi+A{end,2}(3,1,:).*DFi);
                uInterfaces = A{end,2}(1:3,1,:).*DFi.*UFi+A{end,2}(1:3,4,:).*DFo;
                uInterfaces(:,Layers+1,:) = A{end,2}(4:6,1,:).*DFi.*UFi+A{end,2}(4:6,4,:).*DFo;
                for m = Layers:-1:2
                    M0 = A{m}(1:3,1:3,:)-A{m-1,2}(4:6,4:6,:);
                    uInterfaces(:,m,:) = pagemtimes(pagemldivide(M0,A{m-1,2}(4:6,1:3,:)),uInterfaces(:,1,:))-pagemtimes(pagemldivide(M0,A{m}(1:3,4:6,:)),uInterfaces(:,m+1,:));
                end
            end
            for m = 1:Layers
                if  m == 1 && strcmp(Geometry,'Rod')
                    r = (0:Radius(2)/SamplesR:Radius(2))';
                    r(1) = 1e-10;
                else
                    r = (Radius(m):(Radius(m+1)-Radius(m))/SamplesR:Radius(m+1))';
                end
                r2 = r.^2;
                rTotal = [rTotal;r];
                s = (m-1)*length(r)+1:m*length(r);
                U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                kL2 = AngularFrequency2/cL2(m);
                kT2 = AngularFrequency2/cT2(m);
                x2 = kL2-k2;
                y2 = kT2-k2;
                x = sqrt(-x2);
                y = sqrt(-y2);
                xr = x.*r;
                yr = y.*r;
                Znx = besseli(n,xr);
                Zny = besseli(n,yr);
                Zn1x = x.*besseli(n+1,xr);
                Zn1y = y.*besseli(n+1,yr);
                dZnx = n*Znx./r+Zn1x;
                dZny = n*Zny./r+Zn1y;
                dg1Zn1y = -(n+1)*Zn1y./r-y2.*Zny;
                d2Znx = (n*(n-1)./r2-x2).*Znx-Zn1x./r;
                d2Zny = (n*(n-1)./r2-y2).*Zny-Zn1y./r;
                if  m == 1 && strcmp(Geometry,'Rod')
                    f = U(4,1,:).*Znx; % L_in
                    df = U(4,1,:).*dZnx;
                    d2f = U(4,1,:).*d2Znx;
                    g1 = U(5,1,:).*Zn1y; % SV_in
                    dg1 = U(5,1,:).*dg1Zn1y;
                    g3 = U(6,1,:).*Zny; % SH_in
                    dg3 = U(6,1,:).*dZny;
                    d2g3 = U(6,1,:).*d2Zny;
                else
                    Wnx = besselk(n,xr);
                    Wny = besselk(n,yr);
                    Wn1x = x.*besselk(n+1,xr);
                    Wn1y = y.*besselk(n+1,yr);
                    dWnx = n*Wnx./r-Wn1x;
                    dWny = n*Wny./r-Wn1y;
                    dg1Wn1y = -(n+1)*Wn1y./r+y2.*Wny;
                    d2Wnx = (n*(n-1)./r2-x2).*Wnx+Wn1x./r;
                    d2Wny = (n*(n-1)./r2-y2).*Wny+Wn1y./r;
                    f = U(1,1,:).*Znx+U(4,1,:).*Wnx; % L_in + L_out
                    df = U(1,1,:).*dZnx+U(4,1,:).*dWnx;
                    d2f = U(1,1,:).*d2Znx+U(4,1,:).*d2Wnx;
                    g1 = U(2,1,:).*Zn1y+U(5,1,:).*Wn1y; % SV_in + SV_out
                    dg1 = U(2,1,:).*dg1Zn1y+U(5,1,:).*dg1Wn1y;
                    g3 = U(3,1,:).*Zny+U(6,1,:).*Wny; % SH_in + SH_out
                    dg3 = U(3,1,:).*dZny+U(6,1,:).*dWny;
                    d2g3 = U(3,1,:).*d2Zny+U(6,1,:).*d2Wny;
                end
                u(s,1,:) = 1i*(k.*f+(n+1)*g1./r+dg1);
                u(s,2,:) = 1i*(n*f./r-k.*g1+dg3);
                u(s,3,:) = df+k.*g1+n*g3./r;
                v(s,:,:) = -1i*AngularFrequency.*u(s,:,:);
                epsilon(s,1,:) = 1i*k.*u(s,1,:);
                epsilon(s,2,:) = (1i*n*u(s,2,:)+u(s,3,:))./r;
                epsilon(s,3,:) = d2f+k.*dg1-n*g3./r2+n*dg3./r;
                epsilon(s,4,:) = 1i*(-2*n*f./r2+2*n*df./r+(n+1)*k.*g1./r-k.*dg1+2*d2g3+y2.*g3);
                epsilon(s,5,:) = 1i*(2*k.*df+(n*(n+1)./r2+k2-y2).*g1+n*dg1./r+n*k.*g3./r);
                epsilon(s,6,:) = 1i*(n*u(s,1,:)./r+k.*u(s,2,:));
                sigma(s,1:3,:) = -Lambda(m)*kL2.*f+2*Mu(m)*epsilon(s,1:3,:);
                sigma(s,4:6,:) = Mu(m)*epsilon(s,4:6,:);
                KineticEnergyDensity(s,1,:) = Material{m}.Density*sum(abs(v(s,:,:)).^2,2)/2;
            end
            PowerFlowDensity = -real(sigma(:,1,:).*conj(v(:,1,:))+sigma(:,6,:).*conj(v(:,2,:))+sigma(:,5,:).*conj(v(:,3,:)))/2;
            if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
                zri = zi.*rFi;
                if  Sink
                    Znzi = besselh(n,2,zri);
                    Zn1zi = zi.*besselh(n+1,2,zri);
                else
                    Znzi = besselj(n,zri);
                    Zn1zi = zi.*besselj(n+1,zri);
                end
                f = UFi.*Znzi;
                df = UFi.*(n*Znzi./rFi-Zn1zi);
                d2f = UFi.*((n*(n-1)./rFi2-zi.^2).*Znzi+Zn1zi./rFi);
                uFi(:,2,:) = 1i*n*f./rFi;
                uFi(:,3,:) = df;
                epsilonFi(:,2,:) = (1i*n*uFi(:,2,:)+uFi(:,3,:))./rFi;
            end
        elseif strcmp(ModeType,'L')
            L1 = zeros(4,4,Length);
            uInterfaces = zeros(2,Layers+1,Length);
            u = zeros(Height,2,Length);
            v = zeros(Height,2,Length);
            epsilon = zeros(Height,4,Length);
            sigma = zeros(Height,4,Length);
            uFi = zeros(length(rFi),2,Length);
            epsilonFi = zeros(length(rFi),2,Length);
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
                    Z1x = x.*besseli(1,xR);
                    Z1y = y.*besseli(1,yR);
                    d2Z0y = -y2.*Z0y-Z1y/R;
                    d2Z0x = -x2.*Z0x-Z1x/R;
                    L1(3,3,:) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2.*Z0x; % sigma33
                    L1(3,4,:) = 2*Mu(m)*k.*d2Z0y;
                    L1(4,3,:) = 2i*Mu(m)*k.*Z1x; % sigma13
                    L1(4,4,:) = 1i*Mu(m)*(k2-y2).*Z1y;
                    A{m,3}(3,3,:) = 1i*k.*Z0x; % u1
                    A{m,3}(3,4,:) = 1i*(Z1y/R+d2Z0y);
                    A{m,3}(4,3,:) = Z1x; % u3
                    A{m,3}(4,4,:) = k.*Z1y;
                    A{m,1}(3:4,3:4,:) = pagemrdivide(L1(3:4,3:4,:),A{m,3}(3:4,3:4,:));
                else
                    R = Radius(m:m+1);
                    xR = x.*R;
                    yR = y.*R;
                    ZW0x = [besseli(0,xR) besselk(0,xR)];
                    ZW0y = [besseli(0,yR) besselk(0,yR)];
                    ZW1x = x.*[besseli(1,xR) besselk(1,xR)];
                    ZW1y = y.*[besseli(1,yR) besselk(1,yR)];
                    dZW0x = [ZW1x(:,1,:) -ZW1x(:,2,:)];
                    dg1ZW1y = y2.*[-ZW0y(:,1,:) ZW0y(:,2,:)]-ZW1y./R;
                    d2ZW0x = -x2.*ZW0x-dZW0x./R;
                    L1([1 3],[1 3],:) = 2*Mu(m)*d2ZW0x-Lambda(m)*kL2.*ZW0x; % sigma33
                    L1([1 3],[2 4],:) = 2*Mu(m)*k.*dg1ZW1y;
                    L1([2 4],[1 3],:) = 2i*Mu(m)*k.*dZW0x; % sigma13
                    L1([2 4],[2 4],:) = 1i*Mu(m)*(k2-y2).*ZW1y;
                    A{m,3}([1 3],[1 3],:) = 1i*k.*ZW0x; % u1
                    A{m,3}([1 3],[2 4],:) = 1i*(ZW1y./R+dg1ZW1y);
                    A{m,3}([2 4],[1 3],:) = dZW0x; % u3
                    A{m,3}([2 4],[2 4],:) = k.*ZW1y;
                    A{m,1} = pagemrdivide(L1,A{m,3});
                end
            end
            A{1,2} = A{1};
            for m = 2:Layers
                M0 = A{m}(1:2,1:2,:)-A{m-1,2}(3:4,3:4,:);
                M1 = pagemrdivide(A{m-1,2}(1:2,3:4,:),M0);
                M2 = pagemrdivide(A{m}(3:4,1:2,:),M0);
                A{m,2} = [A{m-1,2}(1:2,1:2,:)+pagemtimes(M1,A{m-1,2}(3:4,1:2,:)) -pagemtimes(M1,A{m}(1:2,3:4,:));pagemtimes(M2,A{m-1,2}(3:4,1:2,:)) A{m}(3:4,3:4,:)-pagemtimes(M2,A{m}(1:2,3:4,:))];
            end
            zo = reshape(X{p}(:,8),1,1,[]).*sqrt(AngularFrequency2/cFo2-k2);
            zRo = zo*Radius(end);
            H0zo = besselh(0,zRo);
            DFo = -OuterFluid.Density*AngularFrequency2.*H0zo;
            if  strcmp(Geometry,'Rod')
                A{end,2} = pageinv(A{end,2}(3:4,3:4,:));
                uInterfaces(:,Layers+1,:) = A{end,2}(1:2,1,:).*DFo;
                for m = Layers:-1:2
                    uInterfaces(:,m,:) = pagemtimes(pagemldivide(A{m-1,2}(3:4,3:4,:)-A{m}(1:2,1:2,:),A{m}(1:2,3:4,:)),uInterfaces(:,m+1,:));
                end
            elseif strcmp(Geometry,'Pipe')
                zi = reshape(X{p}(:,9),1,1,[]).*sqrt(AngularFrequency2/cFi2-k2);
                zRi = zi*Radius(1);
                if  Sink
                    Z0zi = besselh(0,2,zRi);
                    Z1zi = zi.*besselh(1,2,zRi);
                else
                    Z0zi = besselj(0,zRi);
                    Z1zi = zi.*besselj(1,zRi);
                end
                WFi = Z1zi;
                DFi = -InnerFluid.Density*AngularFrequency2.*Z0zi;
                A{end,2} = pageinv(A{end,2});
                UFi = -A{end,2}(2,3,:).*DFo./(WFi+A{end,2}(2,1,:).*DFi);
                uInterfaces = A{end,2}(1:2,1,:).*DFi.*UFi+A{end,2}(1:2,3,:).*DFo;
                uInterfaces(:,Layers+1,:) = A{end,2}(3:4,1,:).*DFi.*UFi+A{end,2}(3:4,3,:).*DFo;
                for m = Layers:-1:2
                    M0 = A{m}(1:2,1:2,:)-A{m-1,2}(3:4,3:4,:);
                    uInterfaces(:,m,:) = pagemtimes(pagemldivide(M0,A{m-1,2}(3:4,1:2,:)),uInterfaces(:,1,:))-pagemtimes(pagemldivide(M0,A{m}(1:2,3:4,:)),uInterfaces(:,m+1,:));
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
                U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                kL2 = AngularFrequency2/cL2(m);
                kT2 = AngularFrequency2/cT2(m);
                x2 = kL2-k2;
                y2 = kT2-k2;
                x = sqrt(-x2);
                y = sqrt(-y2);
                xr = x.*r;
                yr = y.*r;
                Z0x = besseli(0,xr);
                Z0y = besseli(0,yr);
                Z1x = x.*besseli(1,xr);
                Z1y = y.*besseli(1,yr);
                dg1Z1y = -Z1y./r-y2.*Z0y;
                d2Z0x = -x2.*Z0x-Z1x./r;
                if  m == 1 && strcmp(Geometry,'Rod')
                    f = U(3,1,:).*Z0x; % L_in
                    df = U(3,1,:).*Z1x;
                    d2f = U(3,1,:).*d2Z0x;
                    g = U(4,1,:).*Z1y; % SV_in
                    dg = U(4,1,:).*dg1Z1y;
                else
                    Wnx = besselk(n,xr);
                    Wny = besselk(n,yr);
                    Wn1x = x.*besselk(n+1,xr);
                    Wn1y = y.*besselk(n+1,yr);
                    dg1Wn1y = -Wn1y./r+y2.*Wny;
                    d2Wnx = -x2.*Wnx+Wn1x./r;
                    f = U(1,1,:).*Z0x+U(3,1,:).*Wnx; % L_in + L_out
                    df = U(1,1,:).*Z1x-U(3,1,:).*Wn1x;
                    d2f = U(1,1,:).*d2Z0x+U(3,1,:).*d2Wnx;
                    g = U(2,1,:).*Z1y+U(4,1,:).*Wn1y; % SV_in + SV_out
                    dg = U(2,1,:).*dg1Z1y+U(4,1,:).*dg1Wn1y;
                end
                u(s,1,:) = 1i*(k.*f+g./r+dg); % u1
                u(s,2,:) = df+k.*g; % u3
                v(s,:,:) = -1i*AngularFrequency.*u(s,:,:);
                epsilon(s,1,:) = 1i*k.*u(s,1,:); % epsilon11
                epsilon(s,2,:) = u(s,2,:)./r; % epsilon22
                epsilon(s,3,:) = d2f+k.*dg; % epsilon33
                epsilon(s,4,:) = 1i*(2*k.*df+(k2-y2).*g); % epsilon13
                sigma(s,1:3,:) = -Lambda(m)*kL2.*f+2*Mu(m)*epsilon(s,1:3,:);
                sigma(s,4,:) = Mu(m)*epsilon(s,4,:); % sigma13
                KineticEnergyDensity(s,1,:) = Material{m}.Density*sum(abs(v(s,:,:)).^2,2)/2;
            end
            PowerFlowDensity = -real(sigma(:,1,:).*conj(v(:,1,:))+sigma(:,4,:).*conj(v(:,2,:)))/2;
            if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
                zri = zi.*rFi;
                if  Sink
                    Z0zi = besselh(0,2,zri);
                    Z1zi = zi.*besselh(1,2,zri);
                else
                    Z0zi = besselj(0,zri);
                    Z1zi = zi.*besselj(1,zri);
                end
                f = UFi.*Z0zi;
                df = -UFi.*Z1zi;
                d2f = UFi.*(Z1zi./rFi-zi.^2.*Z0zi);
                uFi(:,2,:) = df;
                epsilonFi(:,2,:) = uFi(:,2,:)./rFi;
            end
        elseif strcmp(ModeType,'T')
            uInterfaces = ones(1,Layers+1,Length);
            u = zeros(Height,1,Length);
            v = zeros(Height,1,Length);
            epsilon = zeros(Height,2,Length);
            sigma = zeros(Height,2,Length);
            for m = 1:Layers
                kT2 = AngularFrequency2/cT2(m);
                y2 = kT2-k2;
                y = sqrt(-y2);
                if  m == 1 && strcmp(Geometry,'Rod')
                    R = Radius(2);
                    yR = y*R;
                    Z0y = besseli(0,yR);
                    Z1y = y.*besseli(1,yR);
                    A{m,3} = 1i*Z1y; % u2;
                    A{m,1}(2,2,:) = -1i*Mu(m)*(y2.*Z0y+2*Z1y/R)./A{m,3};
                else
                    R = Radius(m:m+1);
                    yR = y.*R;
                    ZW0y = [besseli(0,yR) besselk(0,yR)];
                    ZW1y = y.*[besseli(1,yR) besselk(1,yR)];
                    dZW0y = [ZW1y(:,1,:) -ZW1y(:,2,:)];
                    d2ZW0y = -y2.*ZW0y-dZW0y./R;
                    A{m,3} = 1i*dZW0y;
                    A{m,1} = pagemrdivide(1i*Mu(m)*(2*d2ZW0y+y2.*ZW0y),A{m,3});
                end
            end
            A{1,2} = A{1};
            for m = 2:Layers
                M0 = A{m}(1,1,:)-A{m-1,2}(2,2,:);
                M1 = A{m-1,2}(1,2,:)./M0;
                M2 = A{m}(2,1,:)./M0;
                A{m,2} = [A{m-1,2}(1,1,:)+M1.*A{m-1,2}(2,1,:) -M1.*A{m}(1,2,:);M2.*A{m-1,2}(2,1,:) A{m}(2,2,:)-M2.*A{m}(1,2,:)];
            end
            if  strcmp(Geometry,'Rod')
                for m = Layers:-1:2
                    uInterfaces(:,m,:) = (A{m-1,2}(2,2,:)-A{m}(1,1,:)).\A{m}(1,2,:).*uInterfaces(:,m+1,:);
                end
            elseif strcmp(Geometry,'Pipe')
                In0 = A{1,3}(1,1,:);
                Out0 = A{1,3}(1,2,:);
                OutN = A{end,3}(2,2,:);
                Z1 = [A{end,2}(1,1,:).*In0 -A{end,2}(1,2,:).*OutN;A{end,2}(2,1,:).*In0 A{end,2}(2,2,:).*OutN]; % changed sign in Z1(1,1)!
                Z2 = [A{end,2}(1,1,:).*Out0;-A{end,2}(2,1,:).*Out0];
                RT = pagemldivide(Z1,Z2);
                uInterfaces = Out0+In0.*RT(1,1,:);
                uInterfaces(:,Layers+1,:) = OutN.*RT(2,1,:);
                for m = Layers:-1:2
                    M0 = A{m}(1,1,:)-A{m-1,2}(2,2,:);
                    uInterfaces(:,m,:) = M0.\A{m-1,2}(2,1,:).*uInterfaces(:,1,:)-M0.\A{m}(1,2,:).*uInterfaces(:,m+1,:);
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
                U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                kT2 = AngularFrequency2/cT2(m);
                y2 = kT2-k2;
                y = sqrt(-y2);
                yr = y.*r;
                Z0y = besseli(0,yr);
                Z1y = y.*besseli(1,yr);
                d2Z0y = -y2.*Z0y-Z1y./r;
                if  m == 1 && strcmp(Geometry,'Rod')
                    g = U(2,1,:).*Z0y; % SH_in
                    dg = U(2,1,:).*Z1y;
                    d2g = U(2,1,:).*d2Z0y;
                else
                    W0y = besselk(0,yr);
                    W1y = y.*besselk(1,yr);
                    d2W0y = -y2.*W0y+W1y./r;
                    g = U(1,1,:).*Z0y+U(2,1,:).*W0y; % SH_in + SH_out
                    dg = U(1,1,:).*Z1y-U(2,1,:).*W1y;
                    d2g = U(1,1,:).*d2Z0y+U(2,1,:).*d2W0y;
                end
                u(s,1,:) = 1i*dg; % u2
                v(s,1,:) = -1i*AngularFrequency.*u(s,1,:);
                epsilon(s,1,:) = 1i*(2*d2g+y2.*g); % epsilon23
                epsilon(s,2,:) = 1i*k.*u(s,1,:); % epsilon12
                sigma(s,:,:) = Mu(m)*epsilon(s,:,:);
                KineticEnergyDensity(s,1,:) = Material{m}.Density*abs(v(s,1,:)).^2/2;
            end
            PowerFlowDensity = -real(sigma(:,2,:).*conj(v(:,1,:)))/2;
        end
        StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid && ~strcmp(ModeType,'T')
            uFi(:,1,:) = 1i*k.*f;
            vFi = -1i*AngularFrequency.*uFi;
            epsilonFi(:,1,:) = 1i*k.*uFi(:,1,:);
            epsilonFi(:,3,:) = d2f;
            sigmaFi = -InnerFluid.Density*AngularFrequency2.*f;
            StrainEnergyDensity = [real(sum(epsilonFi.*conj(sigmaFi),2))/2;StrainEnergyDensity];
            KineticEnergyDensity = [InnerFluid.Density*sum(abs(vFi).^2,2)/2;KineticEnergyDensity];
            PowerFlowDensity = [-real(sigmaFi.*conj(vFi(:,1,:)))/2;PowerFlowDensity];
            rTotal2pi = [rFi;rTotal].^2*pi;
        else
            rTotal2pi = rTotal.^2*pi;
        end
        PowerFlow = trapz(rTotal2pi,squeeze(PowerFlowDensity));
        TotalEnergy = trapz(rTotal2pi,squeeze(StrainEnergyDensity+KineticEnergyDensity))/2;
        X{p}(:,5) = fillmissing(filloutliers((PowerFlow./TotalEnergy)'/1e3,'spline','movmedian',5,'ThresholdFactor',1),'spline'); % cez (m/ms)
        Counter = Counter+1;
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
        if  Stop
            close(h)
            return
        end
    end
end