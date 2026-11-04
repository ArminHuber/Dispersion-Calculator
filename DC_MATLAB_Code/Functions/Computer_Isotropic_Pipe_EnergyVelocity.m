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
function [F,L,T] = Computer_Isotropic_Pipe_EnergyVelocity(F,L,T,Material,Force2DTracing,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Ro,Ri,SamplesR)
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
if  ~isempty(T{1}) && (Force2DTracing || Viscoelastic)
    ModeTotal = ModeTotal+length(T);
end
if  ModeTotal == 0
    return
end
if  ~ToggleOuterFluid
    OuterFluid.Velocity = 1e-10;
    OuterFluid.Density = 1e-10;
end
if  ~ToggleInnerFluid
    InnerFluid.Velocity = 1e-10;
    InnerFluid.Density = 1e-10;
end
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
cL2 = Material.LongitudinalVelocity_complex^2;
cT2 = Material.TransverseVelocity_complex^2;
cFo2 = OuterFluid.Velocity^2;
cFi2 = InnerFluid.Velocity^2;
Ro2 = Ro^2;
Ri2 = Ri^2;
r = (Ri:(Ro-Ri)/SamplesR:Ro)';
if  Sink
    rFi = (Ri/SamplesR:Ri/SamplesR:Ri)';
else
    rFi = (0:Ri/SamplesR:Ri)';
    rFi(1) = 1e-10;
end
rTotal = [rFi;r];
r2 = r.^2;
rFi2 = rFi.^2;
rTotal2 = rTotal.^2;
Diff = diff(r2)';
DiffTotal = diff(rTotal2)';
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(F{1})
    for n = 1:length(F)
        [F{n},Counter] = Computer_Isotropic_Pipe_EnergyVelocity_Core(F{n},'F',n,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFo2,cFi2,OuterFluid,InnerFluid,Sink,r,r2,rFi,rFi2,Ro,Ro2,Ri,Ri2,DiffTotal);
        if  Stop
            return
        end
    end
end
if  ~isempty(L{1})
    [L,Counter] = Computer_Isotropic_Pipe_EnergyVelocity_Core(L,'L',0,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFo2,cFi2,OuterFluid,InnerFluid,Sink,r,r2,rFi,rFi2,Ro,Ro2,Ri,Ri2,DiffTotal);
    if  Stop
        return
    end
end
if  ~isempty(T{1}) && (Force2DTracing || Viscoelastic)
    [T,~] = Computer_Isotropic_Pipe_EnergyVelocity_Core(T,'T',0,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFo2,cFi2,OuterFluid,InnerFluid,Sink,r,r2,rFi,rFi2,Ro,Ro2,Ri,Ri2,Diff);
    if  Stop
        return
    end
end
close(h)
end
function [X,Counter] = Computer_Isotropic_Pipe_EnergyVelocity_Core(X,ModeType,n,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFo2,cFi2,OuterFluid,InnerFluid,Sink,r,r2,rFi,rFi2,Ro,Ro2,Ri,Ri2,DiffTotal)
    global Stop 
    Stop = 0;
    for p = 1:length(X)
        PhaseVelocity = reshape(X{p}(:,4),1,1,[])*1e3;
        Attenuation = reshape(X{p}(:,7),1,1,[]); % Np/m
        AngularFrequency = reshape(X{p}(:,1),1,1,[])*pi*2e3;
        AngularFrequency2 = AngularFrequency.^2;
        k = AngularFrequency./PhaseVelocity+1i*Attenuation;
        k2 = k.^2;
        y2 = AngularFrequency2/cT2-k2;
        y = sqrt(-y2);
        yr = y.*r;
        Length = length(k);
        if  strcmp(ModeType,'F')
            kL2 = AngularFrequency2/cL2;
            x2 = kL2-k2;
            x = sqrt(-x2);
            zo = reshape(X{p}(:,8),1,1,[]).*sqrt(AngularFrequency2/cFo2-k2);
            zi = reshape(X{p}(:,9),1,1,[]).*sqrt(AngularFrequency2/cFi2-k2);
            xr = x.*r;
            zri = zi.*rFi;
            zRo = zo*Ro;
            Znx = besseli(n,xr);
            Zny = besseli(n,yr);
            Wnx = besselk(n,xr);
            Wny = besselk(n,yr);
            Zn1x = x.*besseli(n+1,xr);
            Zn1y = y.*besseli(n+1,yr);
            Wn1x = x.*besselk(n+1,xr);
            Wn1y = y.*besselk(n+1,yr);
            dZnx = n*Znx./r+Zn1x;
            dZny = n*Zny./r+Zn1y;
            dWnx = n*Wnx./r-Wn1x;
            dWny = n*Wny./r-Wn1y;
            dg1Zn1y = -(n+1)*Zn1y./r-y2.*Zny;
            dg1Wn1y = -(n+1)*Wn1y./r+y2.*Wny;
            d2Znx = (n*(n-1)./r2-x2).*Znx-Zn1x./r;
            d2Zny = (n*(n-1)./r2-y2).*Zny-Zn1y./r;
            d2Wnx = (n*(n-1)./r2-x2).*Wnx+Wn1x./r;
            d2Wny = (n*(n-1)./r2-y2).*Wny+Wn1y./r;
            if  Sink
                Znzi = besselh(n,2,zri);
                Zn1zi = zi.*besselh(n+1,2,zri);
            else
                Znzi = besselj(n,zri);
                Zn1zi = zi.*besselj(n+1,zri);
            end
            dZnzi = n*Znzi./rFi-Zn1zi;
            Hnzo = besselh(n,zRo);
            dHnzo = n*Hnzo/Ro-zo.*besselh(n+1,zRo);
            Z1 = [Mu*2i*n*(-Wnx(1,1,:)/Ri2+dWnx(1,1,:)/Ri) Mu*1i*k.*((n+1)*Zn1y(1,1,:)/Ri-dg1Zn1y(1,1,:)) Mu*1i*k.*((n+1)*Wn1y(1,1,:)/Ri-dg1Wn1y(1,1,:)) Mu*1i*(2*d2Zny(1,1,:)+y2.*Zny(1,1,:)) Mu*1i*(2*d2Wny(1,1,:)+y2.*Wny(1,1,:)) zeros(1,2,Length);Mu*2i*k.*dWnx(1,1,:) Mu*1i*((n*(n+1)/Ri2+k2-y2).*Zn1y(1,1,:)+n*dg1Zn1y(1,1,:)/Ri) Mu*1i*((n*(n+1)/Ri2+k2-y2).*Wn1y(1,1,:)+n*dg1Wn1y(1,1,:)/Ri) Mu*1i*n*k.*Zny(1,1,:)/Ri Mu*1i*n*k.*Wny(1,1,:)/Ri zeros(1,2,Length);dWnx(1,1,:) k.*Zn1y(1,1,:) k.*Wn1y(1,1,:) n*Zny(1,1,:)/Ri n*Wny(1,1,:)/Ri -dZnzi(end,1,:) zeros(1,1,Length);-Lambda*kL2.*Wnx(end,1,:)+2*Mu*d2Wnx(end,1,:) 2*Mu*k.*dg1Zn1y(end,1,:) 2*Mu*k.*dg1Wn1y(end,1,:) 2*Mu*n*(-Zny(end,1,:)/Ro2+dZny(end,1,:)/Ro) 2*Mu*n*(-Wny(end,1,:)/Ro2+dWny(end,1,:)/Ro) zeros(1,1,Length) OuterFluid.Density*AngularFrequency2.*Hnzo;Mu*2i*n*(-Wnx(end,1,:)/Ro2+dWnx(end,1,:)/Ro) Mu*1i*k.*((n+1)*Zn1y(end,1,:)/Ro-dg1Zn1y(end,1,:)) Mu*1i*k.*((n+1)*Wn1y(end,1,:)/Ro-dg1Wn1y(end,1,:)) Mu*1i*(2*d2Zny(end,1,:)+y2.*Zny(end,1,:)) Mu*1i*(2*d2Wny(end,1,:)+y2.*Wny(end,1,:)) zeros(1,2,Length);Mu*2i*k.*dWnx(end,1,:) Mu*1i*((n*(n+1)/Ro2+k2-y2).*Zn1y(end,1,:)+n*dg1Zn1y(end,1,:)/Ro) Mu*1i*((n*(n+1)/Ro2+k2-y2).*Wn1y(end,1,:)+n*dg1Wn1y(end,1,:)/Ro) Mu*1i*n*k.*Zny(end,1,:)/Ro Mu*1i*n*k.*Wny(end,1,:)/Ro zeros(1,2,Length);dWnx(end,1,:) k.*Zn1y(end,1,:) k.*Wn1y(end,1,:) n*Zny(end,1,:)/Ro n*Wny(end,1,:)/Ro zeros(1,1,Length) -dHnzo];
            Z2 = [Mu*2i*n*(-Znx(1,1,:)/Ri2+dZnx(1,1,:)/Ri);Mu*2i*k.*dZnx(1,1,:);dZnx(1,1,:);-Lambda*kL2.*Znx(end,1,:)+2*Mu*d2Znx(end,1,:);Mu*2i*n*(-Znx(end,1,:)/Ro2+dZnx(end,1,:)/Ro);Mu*2i*k.*dZnx(end,1,:);dZnx(end,1,:)];
            U = -pagemldivide(Z1,Z2);
            f = Znx+U(1,1,:).*Wnx; % L_in + L_out
            df = dZnx+U(1,1,:).*dWnx;
            d2f = d2Znx+U(1,1,:).*d2Wnx;
            g1 = U(2,1,:).*Zn1y+U(3,1,:).*Wn1y; % SV_in + SV_out
            dg1 = U(2,1,:).*dg1Zn1y+U(3,1,:).*dg1Wn1y;
            g3 = U(4,1,:).*Zny+U(5,1,:).*Wny; % SH_in + SH_out
            dg3 = U(4,1,:).*dZny+U(5,1,:).*dWny;
            d2g3 = U(4,1,:).*d2Zny+U(5,1,:).*d2Wny;
            Delta = -Lambda*kL2.*f;
            u1 = 1i*(k.*f+(n+1)*g1./r+dg1);
            u2 = 1i*(n*f./r-k.*g1+dg3);
            u3 = df+k.*g1+n*g3./r;
            v1 = -1i*AngularFrequency.*u1;
            v2 = -1i*AngularFrequency.*u2;
            v3 = -1i*AngularFrequency.*u3;
            epsilon1 = 1i*k.*u1;
            epsilon2 = (1i*n*u2+u3)./r;
            epsilon3 = d2f+k.*dg1-n*g3./r2+n*dg3./r;
            epsilon4 = 1i*(-2*n*f./r2+2*n*df./r+(n+1)*k.*g1./r-k.*dg1+2*d2g3+y2.*g3);
            epsilon5 = 1i*(2*k.*df+(n*(n+1)./r2+k2-y2).*g1+n*dg1./r+n*k.*g3./r);
            epsilon6 = 1i*(n*u1./r+k.*u2);
            sigma1 = Delta+2*Mu*epsilon1;
            sigma2 = Delta+2*Mu*epsilon2;
            sigma3 = Delta+2*Mu*epsilon3;
            sigma4 = Mu*epsilon4;
            sigma5 = Mu*epsilon5;
            sigma6 = Mu*epsilon6;
            StrainEnergyDensity = real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4)+epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6))/2;
            KineticEnergyDensity = Material.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2)/2;
            PowerFlowDensity = reshape(-real(sigma1.*conj(v1)+sigma6.*conj(v2)+sigma5.*conj(v3))/2,length(r),Length);
            f = U(6,1,:).*Znzi;
            df = U(6,1,:).*dZnzi;
            d2f = U(6,1,:).*((n*(n-1)./rFi2-zi.^2).*Znzi+Zn1zi./rFi);
            u1 = 1i*k.*f;
            u2 = 1i*n*f./rFi;
            u3 = df;
            v1 = -1i*AngularFrequency.*u1;
            v2 = -1i*AngularFrequency.*u2;
            v3 = -1i*AngularFrequency.*u3;
            epsilon1 = 1i*k.*u1;
            epsilon2 = (1i*n*u2+u3)./rFi;
            epsilon3 = d2f;
            sigma = -InnerFluid.Density*AngularFrequency2.*f;
            StrainEnergyDensityInnerFluid = real(epsilon1.*conj(sigma)+epsilon2.*conj(sigma)+epsilon3.*conj(sigma))/2;
            KineticEnergyDensityInnerFluid = InnerFluid.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2)/2;
            PowerFlowDensityInnerFluid = reshape(-real(sigma.*conj(v1))/2,length(rFi),Length);
            StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
            KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
            PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
        elseif strcmp(ModeType,'L')
            kL2 = AngularFrequency2/cL2;
            x2 = kL2-k2;
            x = sqrt(-x2);
            zo = reshape(X{p}(:,8),1,1,[]).*sqrt(AngularFrequency2/cFo2-k2);
            zi = reshape(X{p}(:,9),1,1,[]).*sqrt(AngularFrequency2/cFi2-k2);
            xr = x.*r;
            zri = zi.*rFi;
            zRo = zo*Ro;
            Z0x = besseli(0,xr);
            Z0y = besseli(0,yr);
            W0x = besselk(0,xr);
            W0y = besselk(0,yr);
            Z1x = -x.*besseli(1,xr);
            Z1y = -y.*besseli(1,yr);
            W1x = x.*besselk(1,xr);
            W1y = y.*besselk(1,yr);
            if  Sink
                Z0zi = besselh(0,2,zri);
                Z1zi = zi.*besselh(1,2,zri);
            else
                Z0zi = besselj(0,zri);
                Z1zi = zi.*besselj(1,zri);
            end
            H0zo = besselh(0,zRo);
            H1zo = zo.*besselh(1,zRo);
            Z1 = [-Mu*2i*k.*W1x(1,1,:) Mu*1i*(k2-y2).*Z1y(1,1,:) Mu*1i*(k2-y2).*W1y(1,1,:) zeros(1,2,Length);-W1x(1,1,:) k.*Z1y(1,1,:) k.*W1y(1,1,:) Z1zi(end,1,:) zeros(1,1,Length);-Lambda*kL2.*W0x(end,1,:)+2*Mu*(W1x(end,1,:)/Ro-x2.*W0x(end,1,:)) 2*Mu*k.*(-Z1y(end,1,:)/Ro+y2.*Z0y(end,1,:)) 2*Mu*k.*(-W1y(end,1,:)/Ro+y2.*W0y(end,1,:)) zeros(1,1,Length) OuterFluid.Density*AngularFrequency2.*H0zo;-Mu*2i*k.*W1x(end,1,:) Mu*1i*(k2-y2).*Z1y(end,1,:) Mu*1i*(k2-y2).*W1y(end,1,:) zeros(1,2,Length);-W1x(end,1,:) k.*Z1y(end,1,:) k.*W1y(end,1,:) zeros(1,1,Length) H1zo];
            Z2 = [Mu*2i*k.*Z1x(1,1,:);Z1x(1,1,:);Lambda*kL2.*Z0x(end,1,:)-2*Mu*(Z1x(end,1,:)/Ro-x2.*Z0x(end,1,:));Mu*2i*k.*Z1x(end,1,:);Z1x(end,1,:)];
            U = pagemldivide(Z1,Z2);
            f = Z0x+U(1,1,:).*W0x; % L_in + L_out
            df = -Z1x-U(1,1,:).*W1x;
            d2f = Z1x./r-x2.*Z0x+U(1,1,:).*(W1x./r-x2.*W0x);
            g = U(2,1,:).*Z1y+U(3,1,:).*W1y; % SV_in + SV_out
            dg = U(2,1,:).*(-Z1y./r+y2.*Z0y)+U(3,1,:).*(-W1y./r+y2.*W0y);
            Delta = -Lambda*kL2.*f;
            u1 = 1i*(k.*f+g./r+dg);
            u3 = df+k.*g;
            v1 = -1i*AngularFrequency.*u1;
            v3 = -1i*AngularFrequency.*u3;
            epsilon1 = 1i*k.*u1;
            epsilon2 = u3./r;
            epsilon3 = d2f+k.*dg;
            epsilon5 = 1i*(2*k.*df+(k2-y2).*g);
            sigma1 = Delta+2*Mu*epsilon1;
            sigma2 = Delta+2*Mu*epsilon2;
            sigma3 = Delta+2*Mu*epsilon3;
            sigma5 = Mu*epsilon5;
            StrainEnergyDensity = real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon5.*conj(sigma5))/2;
            KineticEnergyDensity = Material.Density*(abs(v1).^2+abs(v3).^2)/2;
            PowerFlowDensity = reshape(-real(sigma1.*conj(v1)+sigma5.*conj(v3))/2,length(r),Length);
            f = U(4,1,:).*Z0zi;
            df = -U(4,1,:).*Z1zi;
            d2f = U(4,1,:).*(Z1zi./rFi-zi.^2.*Z0zi);
            u1 = 1i*k.*f;
            u3 = df;
            v1 = -1i*AngularFrequency.*u1;
            v3 = -1i*AngularFrequency.*u3;
            epsilon1 = 1i*k.*u1;
            epsilon2 = u3./rFi;
            epsilon3 = d2f;
            sigma = -InnerFluid.Density*AngularFrequency2.*f;
            StrainEnergyDensityInnerFluid = real(epsilon1.*conj(sigma)+epsilon2.*conj(sigma)+epsilon3.*conj(sigma))/2;
            KineticEnergyDensityInnerFluid = InnerFluid.Density*(abs(v1).^2+abs(v3).^2)/2;
            PowerFlowDensityInnerFluid = reshape(-real(sigma.*conj(v1))/2,length(rFi),Length);
            StrainEnergyDensity = [StrainEnergyDensityInnerFluid;StrainEnergyDensity];
            KineticEnergyDensity = [KineticEnergyDensityInnerFluid;KineticEnergyDensity];
            PowerFlowDensity = [PowerFlowDensityInnerFluid;PowerFlowDensity];
        elseif strcmp(ModeType,'T')
            U = -besselk(2,yr(end,1,:)).\besseli(2,yr(end,1,:));
            J0 = besseli(0,yr);
            Y0 = besselk(0,yr);
            J1 = -y.*besseli(1,yr);
            Y1 = y.*besselk(1,yr);
            g = J0+U.*Y0; % SH_in + SH_out
            dg = -J1-U.*Y1;
            d2g = J1./r-y2.*J0+U.*(Y1./r-y2.*Y0);
            u2 = 1i*dg;
            v2 = -1i*AngularFrequency.*u2;
            epsilon4 = 1i*(2*d2g+y2.*g);
            epsilon6 = 1i*k.*u2;
            sigma4 = Mu*epsilon4;
            sigma6 = Mu*epsilon6;
            StrainEnergyDensity = real(epsilon4.*conj(sigma4)+epsilon6.*conj(sigma6))/2;
            KineticEnergyDensity = Material.Density*abs(v2).^2/2;
            PowerFlowDensity = reshape(-real(sigma6.*conj(v2))/2,length(r),Length);
        end
        PowerFlow = DiffTotal*(PowerFlowDensity(1:end-1,:)+PowerFlowDensity(2:end,:)); % PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
        TotalEnergy = DiffTotal*(StrainEnergyDensity(1:end-1,:)+StrainEnergyDensity(2:end,:)+KineticEnergyDensity(1:end-1,:)+KineticEnergyDensity(2:end,:))/2; % TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
        X{p}(:,5) = fillmissing(filloutliers((PowerFlow./TotalEnergy)'/1e3,'spline','movmedian',5,'ThresholdFactor',1),'spline'); % cez (m/ms)
        Counter = Counter+1;
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
        if  Stop
            close(h)
            return
        end
    end
end