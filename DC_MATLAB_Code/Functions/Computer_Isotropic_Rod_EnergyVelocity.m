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
function [F,L,T] = Computer_Isotropic_Rod_EnergyVelocity(F,L,T,Material,Force2DTracing,Viscoelastic,FluidLoading,Fluid,R,SamplesR)
%#ok<*GVMIS>
global Stop 
Stop = 0;
ModeTotal = 0;
if  ~isempty(F{1})
    for n = 1:length(F)
        ModeTotal = ModeTotal+length(F{n});
    end
end
if  ~isempty(L{1}) && (Force2DTracing || Viscoelastic || FluidLoading)
    ModeTotal = ModeTotal+length(L);
end
if  ~isempty(T{1}) && (Force2DTracing || Viscoelastic)
    ModeTotal = ModeTotal+length(T);
end
if  ModeTotal == 0
    return
end
if  ~FluidLoading
    Fluid.Velocity = 1e-10;
    Fluid.Density = 1e-10;
end
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
cL2 = Material.LongitudinalVelocity_complex^2;
cT2 = Material.TransverseVelocity_complex^2;
cF2 = Fluid.Velocity^2;
R2 = R^2;
r = (0:R/SamplesR:R)';
r(1) = 1e-10;
r2 = r.^2;
Diff = diff(r2)';
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(F{1})
    for n = 1:length(F)
        [F{n},Counter] = Computer_Isotropic_Rod_EnergyVelocity_Core(F{n},'F',n,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cF2,r,r2,R,R2,Diff);
        if  Stop
            return
        end
    end
end
if  ~isempty(L{1}) && (Force2DTracing || Viscoelastic || FluidLoading)
    [L,Counter] = Computer_Isotropic_Rod_EnergyVelocity_Core(L,'L',0,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cF2,r,r2,R,R2,Diff);
    if  Stop
        return
    end
end
if  ~isempty(T{1}) && (Force2DTracing || Viscoelastic)
    [T,~] = Computer_Isotropic_Rod_EnergyVelocity_Core(T,'T',0,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cF2,r,r2,R,R2,Diff);
    if  Stop
        return
    end
end
close(h)
end
function [X,Counter] = Computer_Isotropic_Rod_EnergyVelocity_Core(X,ModeType,n,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cF2,r,r2,R,R2,Diff)
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
        y = sqrt(y2);
        yr = y.*r;
        Length = length(k);
        if  strcmp(ModeType,'F')
            kL2 = AngularFrequency2/cL2;
            x2 = kL2-k2;
            x = sqrt(x2);
            z = reshape(X{p}(:,8),1,1,[]).*sqrt(AngularFrequency2/cF2-k2);
            xr = x.*r;
            zR = z*R;
            Jnx = besselj(n,xr);
            Jny = besselj(n,yr);
            Jn1x = x.*besselj(n+1,xr);
            Jn1y = y.*besselj(n+1,yr);
            dJnx = n*Jnx./r-Jn1x;
            dJny = n*Jny./r-Jn1y;
            dg1Jn1y = -(n+1)*Jn1y./r+y2.*Jny;
            d2Jnx = (n*(n-1)./r2-x2).*Jnx+Jn1x./r;
            d2Jny = (n*(n-1)./r2-y2).*Jny+Jn1y./r;
            dHnz = -n*besselh(n,zR)/R+z.*besselh(n+1,zR);
            Z1 = [Mu*1i*k.*((n+1)*Jn1y(end,1,:)/R-dg1Jn1y(end,1,:)) Mu*1i*(2*d2Jny(end,1,:)+y2.*Jny(end,1,:)) zeros(1,1,Length);Mu*1i*((n*(n+1)/R2+k2-y2).*Jn1y(end,1,:)+n*dg1Jn1y(end,1,:)/R) Mu*1i*n*k.*Jny(end,1,:)/R zeros(1,1,Length);k.*Jn1y(end,1,:) n*Jny(end,1,:)/R dHnz];
            Z2 = [Mu*2i*n*(-Jnx(end,1,:)/R2+dJnx(end,1,:)/R);Mu*2i*k.*dJnx(end,1,:);dJnx(end,1,:)];
            U = -pagemldivide(Z1,Z2);
            f = Jnx;
            df = dJnx;
            d2f = d2Jnx;
            g1 = U(1,1,:).*Jn1y;
            dg1 = U(1,1,:).*dg1Jn1y;
            g3 = U(2,1,:).*Jny;
            dg3 = U(2,1,:).*dJny;
            d2g3 = U(2,1,:).*d2Jny;
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
            epsilon6 = 1i*n*u1./r+1i*k.*u2;
            sigma1 = Delta+2*Mu*epsilon1;
            sigma2 = Delta+2*Mu*epsilon2;
            sigma3 = Delta+2*Mu*epsilon3;
            sigma4 = Mu*epsilon4;
            sigma5 = Mu*epsilon5;
            sigma6 = Mu*epsilon6;
            StrainEnergyDensity = real(epsilon1.*conj(sigma1)+epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4)+epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6))/2;
            KineticEnergyDensity = Material.Density*(abs(v1).^2+abs(v2).^2+abs(v3).^2)/2;
            PowerFlowDensity = reshape(-real(sigma1.*conj(v1)+sigma6.*conj(v2)+sigma5.*conj(v3))/2,length(r),Length);
        elseif strcmp(ModeType,'L')
            kL2 = AngularFrequency2/cL2;
            x2 = kL2-k2;
            x = sqrt(x2);
            z = reshape(X{p}(:,8),1,1,[]).*sqrt(AngularFrequency2/cF2-k2);
            xr = x.*r;
            zR = z*R;
            J0x = besselj(0,xr);
            J0y = besselj(0,yr);
            J1x = x.*besselj(1,xr);
            J1y = y.*besselj(1,yr);
            H1z = z.*besselh(1,zR);
            Z1 = [Mu*1i*(k2-y2).*J1y(end,1,:) zeros(1,1,Length);k.*J1y(end) H1z];
            Z2 = [Mu*2i*k.*J1x(end,1,:);J1x(end,1,:)];
            U = pagemldivide(Z1,Z2);
            f = J0x;
            df = -J1x;
            d2f = J1x./r-x2.*J0x;
            g = U(1,1,:).*J1y;
            dg = U(1,1,:).*(-J1y./r+y2.*J0y);
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
        elseif strcmp(ModeType,'T')
            J0 = besselj(0,yr);
            J1 = y.*besselj(1,yr);
            g = J0;
            dg = -J1;
            d2g = J1./r-y2.*J0;
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
        PowerFlow = Diff*(PowerFlowDensity(1:end-1,:)+PowerFlowDensity(2:end,:)); % PowerFlow = Diff*(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))*pi/2;
        TotalEnergy = Diff*(StrainEnergyDensity(1:end-1,:)+StrainEnergyDensity(2:end,:)+KineticEnergyDensity(1:end-1,:)+KineticEnergyDensity(2:end,:))/2; % TotalEnergy = Diff*(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))*pi/4;
        X{p}(:,5) = fillmissing(filloutliers((PowerFlow./TotalEnergy)'/1e3,'spline','movmedian',5,'ThresholdFactor',1),'spline'); % cez (m/ms)
        Counter = Counter+1;
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
        if  Stop
            close(h)
            return
        end
    end
end