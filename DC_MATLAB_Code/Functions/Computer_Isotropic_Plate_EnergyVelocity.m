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
function [SLamb,ALamb,BLamb,SShear,AShear] = Computer_Isotropic_Plate_EnergyVelocity(SLamb,ALamb,BLamb,SShear,AShear,Material,Force2DTracing,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3)
%#ok<*GVMIS>
global Stop
Stop = 0;
ModeTotal = 0;
if  ~isempty(SLamb{1}) && (Force2DTracing || FluidLoading || Viscoelastic)
    ModeTotal = ModeTotal+length(SLamb);
end
if  ~isempty(ALamb{1}) && (Force2DTracing || FluidLoading || Viscoelastic)
    ModeTotal = ModeTotal+length(ALamb);
end
if  ~isempty(BLamb{1})
    ModeTotal = ModeTotal+length(BLamb);
end
if  ~isempty(SShear{1}) && Viscoelastic
    ModeTotal = ModeTotal+length(SShear);
end
if  ~isempty(AShear{1}) && Viscoelastic
    ModeTotal = ModeTotal+length(AShear);
end
if  ModeTotal == 0
    return
end
if  ~ToggleUpperFluid
    UpperFluid.Velocity = 1e-10;
    UpperFluid.Density = 1e-10;
end
if  ~ToggleLowerFluid
    LowerFluid.Velocity = 1e-10;
    LowerFluid.Density = 1e-10;
end
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
cL2 = Material.LongitudinalVelocity_complex^2;
cT2 = Material.TransverseVelocity_complex^2;
cFu2 = UpperFluid.Velocity^2;
cFl2 = LowerFluid.Velocity^2;
T = (-Half:2*Half/SamplesX3:Half)';
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(SLamb{1}) && (Force2DTracing || FluidLoading || Viscoelastic)
    [SLamb,Counter] = Computer(SLamb,'Lamb',ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFu2,cFl2,LowerFluid,Half,T);
    if  Stop
        return
    end
end
if  ~isempty(ALamb{1}) && (Force2DTracing || FluidLoading || Viscoelastic)
    [ALamb,Counter] = Computer(ALamb,'Lamb',ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFu2,cFl2,LowerFluid,Half,T);
    if  Stop
        return
    end
end
if  ~isempty(BLamb{1})
    [BLamb,Counter] = Computer(BLamb,'Lamb',ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFu2,cFl2,LowerFluid,Half,T);
    if  Stop
        return
    end
end
if  ~isempty(SShear{1}) && Viscoelastic
    [SShear,Counter] = Computer(SShear,'SShear',ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFu2,cFl2,LowerFluid,Half,T);
    if  Stop
        return
    end
end
if  ~isempty(AShear{1}) && Viscoelastic
    AShear = Computer(AShear,'AShear',ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFu2,cFl2,LowerFluid,Half,T);
    if  Stop
        return
    end
end
close(h)
end
function [X,Counter] = Computer(X,ModeType,ModeTotal,Counter,h,Material,Lambda,Mu,cL2,cT2,cFu2,cFl2,LowerFluid,Half,T)
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
        if  strcmp(ModeType,'Lamb')
            kL2 = AngularFrequency2/cL2;
            kT2 = AngularFrequency2/cT2;
            x2 = kL2-k2;
            y2 = kT2-k2;
            x = sqrt(x2);
            y = sqrt(y2);
            zu = 1i*reshape(X{p}(:,8),1,1,[]).*sqrt(AngularFrequency2/cFu2-k2);
            zl = 1i*reshape(X{p}(:,9),1,1,[]).*sqrt(AngularFrequency2/cFl2-k2);
            xT = x.*T;
            yT = y.*T;
            SinxT = sin(xT);
            SinyT = sin(yT);
            CosxT = cos(xT);
            CosyT = cos(yT);
            E_ = exp(zl*Half);
            a0 = 2i*Mu*k;
            a1 = -Lambda*kL2-2*Mu*x2;
            a2 = Mu*(k2-y2);
            a3 = a0.*x;
            a4 = a0.*y;
            a5 = -1i*k;
            Z1 = [-a3.*SinxT(end,1,:) a2.*SinyT(end,1,:) a2.*CosyT(end,1,:) zeros(1,2,Length);-x.*SinxT(end,1,:) a5.*SinyT(end,1,:) a5.*CosyT(end,1,:) -zu.*exp(zu*Half) zeros(1,1,Length);a1.*CosxT(1,1,:) -a4.*CosyT(1,1,:) a4.*SinyT(1,1,:) zeros(1,1,Length) LowerFluid.Density*AngularFrequency2.*E_;-a3.*SinxT(1,1,:) a2.*SinyT(1,1,:) a2.*CosyT(1,1,:) zeros(1,2,Length);-x.*SinxT(1,1,:) a5.*SinyT(1,1,:) a5.*CosyT(1,1,:) zeros(1,1,Length) zl.*E_];
            Z2 = [a3.*CosxT(end,1,:);x.*CosxT(end,1,:);a1.*SinxT(1,1,:);a3.*CosxT(1,1,:);x.*CosxT(1,1,:)];
            U = -pagemldivide(Z1,Z2);
            f = SinxT+U(1,1,:).*CosxT; % L- + L+
            df = x.*(CosxT-U(1,1,:).*SinxT);
            d2f = -x2.*f;
            g = U(2,1,:).*SinyT+U(3,1,:).*CosyT; % SV- + SV+
            dg = y.*(U(2,1,:).*CosyT-U(3,1,:).*SinyT);
            d2g = -y2.*g;
            Delta = -Lambda*kL2.*f;
            u1 = 1i*k.*f+dg;
            u3 = df-1i*k.*g;
            v1 = -1i*AngularFrequency.*u1;
            v3 = -1i*AngularFrequency.*u3;
            epsilon1 = 1i*k.*u1;
            epsilon3 = d2f-1i*k.*dg;
            epsilon5 = 2i*k.*df+k2.*g+d2g;
            sigma1 = Delta+2*Mu*epsilon1;
            sigma3 = Delta+2*Mu*epsilon3;
            sigma5 = Mu*epsilon5;
            StrainEnergyDensity = real(epsilon1.*conj(sigma1)+epsilon3.*conj(sigma3)+epsilon5.*conj(sigma5))/2;
            KineticEnergyDensity = Material.Density*(abs(v1).^2+abs(v3).^2)/2;
            PowerFlowDensity = reshape(-real(sigma1.*conj(v1)+sigma5.*conj(v3))/2,length(T),Length);
        elseif contains(ModeType,'Shear')
            y = sqrt(AngularFrequency2/cT2-k2);
            yT = y.*T;
            if  strcmp(ModeType,'SShear')
                g = cos(yT);
                dg = -y.*sin(yT);
            elseif strcmp(ModeType,'AShear')
                g = sin(yT);
                dg = y.*cos(yT);
            end
            v2 = -1i*AngularFrequency.*g;
            epsilon4 = dg;
            epsilon6 = 1i*k.*g;
            sigma4 = Mu*epsilon4;
            sigma6 = Mu*epsilon6;
            StrainEnergyDensity = real(epsilon4.*conj(sigma4)+epsilon6.*conj(sigma6))/2;
            KineticEnergyDensity = Material.Density*abs(v2).^2/2;
            PowerFlowDensity = reshape(-real(sigma6.*conj(v2))/2,length(T),Length);
        end
        PowerFlow = sum(PowerFlowDensity(1:end-1,:)+PowerFlowDensity(2:end,:)); % PowerFlow = deltaT*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
        TotalEnergy = sum(StrainEnergyDensity(1:end-1,:)+StrainEnergyDensity(2:end,:)+KineticEnergyDensity(1:end-1,:)+KineticEnergyDensity(2:end,:))/2; % TotalEnergy = deltaT*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
        X{p}(:,5) = fillmissing(filloutliers((PowerFlow./TotalEnergy)'/1e3,'spline','movmedian',5,'ThresholdFactor',1),'spline'); % ce1 (m/ms)
        Counter = Counter+1;
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
        if  Stop
            close(h)
            return
        end
    end
end