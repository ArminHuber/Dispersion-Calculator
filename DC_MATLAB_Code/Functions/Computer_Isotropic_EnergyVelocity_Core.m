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
function [X,Counter] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,X,ModeType,ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3)
%#ok<*PFBNS>
%#ok<*PFTIN>
%#ok<*PFTUSW>
%#ok<*GVMIS>
global Stop 
Stop = 0;
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
if  ~ToggleUpperFluid
    UpperFluid.Velocity = 1e-10;
    UpperFluid.Density = 1e-10;
end
if  ~ToggleLowerFluid
    LowerFluid.Velocity = 1e-10;
    LowerFluid.Density = 1e-10;
end
deltaT = 2*Half/SamplesX3;
T = (-Half:deltaT:Half)';
if  (strcmp(ModeType,'SLamb') || strcmp(ModeType,'ALamb')) && ~Viscoelastic && ~FluidLoading
    cL2 = Material.LongitudinalVelocity^2;
    cT2 = Material.TransverseVelocity^2;
end
for p = 1:length(X)
    if  contains(ModeType,'Lamb') || contains(ModeType,'Scholte')
        if  contains(ModeType,'Lamb') && ~Viscoelastic && ~FluidLoading
            PhaseVelocity = X{p}(:,4)*1e3;
            AngularFrequency2 = (2*pi*X{p}(:,1)*1e3).^2;
            k2 = AngularFrequency2./PhaseVelocity.^2;
            x = sqrt(AngularFrequency2/cL2-k2);
            y = sqrt(AngularFrequency2/cT2-k2);
            xH = x*Half;
            yH = y*Half;
            a1 = 2*cT2*k2;
            a2 = AngularFrequency2-a1;
            if  strcmp(ModeType,'SLamb')
                a3 = a2.^2*Half./x./sin(xH).^2;
                a4 = 2*a1.*xH./sin(yH).^2;
                a5 = 2*a1*cT2.*cot(yH);
                a6 = 4*a2.*cot(xH);
                A1 = 2*a6*cT2+(y./x+x./y-2*x.*y./k2).*a5-a3-a4*cT2;
                A2 = a6+(y./x/cL2+x./y/cT2).*a5-a3/cL2-a4;
            elseif strcmp(ModeType,'ALamb')
                a3 = a2.^2*Half./x./cos(xH).^2;
                a4 = 2*a1.*xH./cos(yH).^2;
                a5 = 2*a1*cT2.*tan(yH);
                a6 = 4*a2.*tan(xH);
                A1 = 2*a6*cT2+(y./x+x./y-2*x.*y./k2).*a5+a3+a4*cT2;
                A2 = a6+(y./x/cL2+x./y/cT2).*a5+a3/cL2+a4;
            end
            X{p}(:,5) = A1./A2./PhaseVelocity/1e3; % ce1 (m/ms)
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
                    kUpperFluid2 = (AngularFrequency(q)/UpperFluid.Velocity)^2;
                    kLowerFluid2 = (AngularFrequency(q)/LowerFluid.Velocity)^2;
                    x = sqrt(kL2-k2);
                    y = sqrt(kT2-k2);
                    zu = sqrt(kUpperFluid2-k2);
                    zl = sqrt(kLowerFluid2-k2);
                    if  contains(ModeType,'Scholte') && Attenuation(q) ~= 0
                        if  PhaseVelocity(q) < UpperFluid.Velocity
                            zu = -zu;
                        end
                        if  PhaseVelocity(q) < LowerFluid.Velocity
                            zl = -zl;
                        end
                    end
                    x2 = kL2-k2;
                    y2 = kT2-k2;
                    xT = x*T;
                    yT = y*T;
                    SinxT = sin(xT);
                    SinyT = sin(yT);
                    CosxT = cos(xT);
                    CosyT = cos(yT);
                    e_ = exp(-1i*zl*Half);
                    a1 = -(Lambda*kL2+2*Mu*x2);
                    a2 = Mu*(k2-y2);
                    a3 = 2i*Mu*k*x;
                    a4 = 2i*Mu*k*y;
                    a5 = -1i*k;
                    Z1 = [-a3*SinxT(end) a2*SinyT(end) a2*CosyT(end) 0 0;-x*SinxT(end) a5*SinyT(end) a5*CosyT(end) 1i*zu*exp(1i*zu*Half) 0;a1*CosxT(1) -a4*CosyT(1) a4*SinyT(1) 0 LowerFluid.Velocity^2*LowerFluid.Density*kLowerFluid2*e_;-a3*SinxT(1) a2*SinyT(1) a2*CosyT(1) 0 0;-x*SinxT(1) a5*SinyT(1) a5*CosyT(1) 0 1i*zl*e_];
                    Z2 = [a3*CosxT(end);x*CosxT(end);a1*SinxT(1);a3*CosxT(1);x*CosxT(1)];
                    U = -Z1\Z2;
                    f = SinxT+U(1)*CosxT; % L- + L+
                    df = x*(CosxT-U(1)*SinxT);
                    d2f = -x2*f;
                    g = U(2)*SinyT+U(3)*CosyT; % SV- + SV+
                    dg = y*(U(2)*CosyT-U(3)*SinyT);
                    d2g = -y2*g;
                    Delta = -Lambda*kL2*f;
                    u1 = 1i*k*f+dg;
                    u3 = df-1i*k*g;
                    v1 = -1i*AngularFrequency(q)*u1;
                    v3 = -1i*AngularFrequency(q)*u3;
                    epsilon1 = 1i*k*u1;
                    epsilon3 = d2f-1i*k*dg;
                    epsilon5 = 2i*k*df+k2*g+d2g;
                    sigma1 = Delta+2*Mu*epsilon1;
                    sigma3 = Delta+2*Mu*epsilon3;
                    sigma5 = Mu*epsilon5;
                    StrainEnergyDensity = .5*real(epsilon1.*conj(sigma1)+epsilon3.*conj(sigma3)+epsilon5.*conj(sigma5));
                    KineticEnergyDensity = .5*Material.Density*(abs(v1).^2+abs(v3).^2);
                    PowerFlowDensity = -.5*real(sigma1.*conj(v1)+sigma5.*conj(v3));
                    PowerFlow = sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = deltaT*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                    TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = deltaT*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                    EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % ce1 (m/ms)
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
                    kUpperFluid2 = (AngularFrequency/UpperFluid.Velocity)^2;
                    kLowerFluid2 = (AngularFrequency/LowerFluid.Velocity)^2;
                    x = sqrt(kL2-k2);
                    y = sqrt(kT2-k2);
                    zu = sqrt(kUpperFluid2-k2);
                    zl = sqrt(kLowerFluid2-k2);
                    if  contains(ModeType,'Scholte') && Attenuation ~= 0
                        if  PhaseVelocity < UpperFluid.Velocity
                            zu = -zu;
                        end
                        if  PhaseVelocity < LowerFluid.Velocity
                            zl = -zl;
                        end
                    end
                    x2 = kL2-k2;
                    y2 = kT2-k2;
                    xT = x*T;
                    yT = y*T;
                    SinxT = sin(xT);
                    SinyT = sin(yT);
                    CosxT = cos(xT);
                    CosyT = cos(yT);
                    e_ = exp(-1i*zl*Half);
                    a1 = -(Lambda*kL2+2*Mu*x2);
                    a2 = Mu*(k2-y2);
                    a3 = 2i*Mu*k*x;
                    a4 = 2i*Mu*k*y;
                    a5 = -1i*k;
                    Z1(1,1) = -a3*SinxT(end);
                    Z1(1,2) = a2*SinyT(end);
                    Z1(1,3) = a2*CosyT(end);
                    Z1(2,1) = -x*SinxT(end);
                    Z1(2,2) = a5*SinyT(end);
                    Z1(2,3) = a5*CosyT(end);
                    Z1(2,4) = 1i*zu*exp(1i*zu*Half);
                    Z1(3,1) = a1*CosxT(1);
                    Z1(3,2) = -a4*CosyT(1);
                    Z1(3,3) = a4*SinyT(1);
                    Z1(3,5) = LowerFluid.Velocity^2*LowerFluid.Density*kLowerFluid2*e_;
                    Z1(4,1) = -a3*SinxT(1);
                    Z1(4,2) = a2*SinyT(1);
                    Z1(4,3) = a2*CosyT(1);
                    Z1(5,1) = -x*SinxT(1);
                    Z1(5,2) = a5*SinyT(1);
                    Z1(5,3) = a5*CosyT(1);
                    Z1(5,5) = 1i*zl*e_;
                    Z2(1,1) = a3*CosxT(end);
                    Z2(2,1) = x*CosxT(end);
                    Z2(3,1) = a1*SinxT(1);
                    Z2(4,1) = a3*CosxT(1);
                    Z2(5,1) = x*CosxT(1);
                    U = -Z1\Z2;
                    f = SinxT+U(1)*CosxT; % L- + L+
                    df = x*(CosxT-U(1)*SinxT);
                    d2f = -x2*f;
                    g = U(2)*SinyT+U(3)*CosyT; % SV- + SV+
                    dg = y*(U(2)*CosyT-U(3)*SinyT);
                    d2g = -y2*g;
                    Delta = -Lambda*kL2*f;
                    u1 = 1i*k*f+dg;
                    u3 = df-1i*k*g;
                    v1 = -1i*AngularFrequency*u1;
                    v3 = -1i*AngularFrequency*u3;
                    epsilon1 = 1i*k*u1;
                    epsilon3 = d2f-1i*k*dg;
                    epsilon5 = 2i*k*df+k2*g+d2g;
                    sigma1 = Delta+2*Mu*epsilon1;
                    sigma3 = Delta+2*Mu*epsilon3;
                    sigma5 = Mu*epsilon5;
                    StrainEnergyDensity = .5*real(epsilon1.*conj(sigma1)+epsilon3.*conj(sigma3)+epsilon5.*conj(sigma5));
                    KineticEnergyDensity = .5*Material.Density*(abs(v1).^2+abs(v3).^2);
                    PowerFlowDensity = -.5*real(sigma1.*conj(v1)+sigma5.*conj(v3));
                    PowerFlow = sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = deltaT*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                    TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = deltaT*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                    X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % ce1 (m/ms)
                end
            end
        end
    elseif contains(ModeType,'Shear')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            Attenuation = X{p}(:,6).*PhaseVelocity./(X{p}(:,1)*1e3); % Np/wavelength
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                k = AngularFrequency(q)/PhaseVelocity(q)*(1+1i*Attenuation(q)/(2*pi));
                y = sqrt(AngularFrequency(q)^2/Material.TransverseVelocity_complex^2-k^2);
                yT = y*T;
                SinyT = sin(yT);
                CosyT = cos(yT);
                if  strcmp(ModeType,'SShear')
                    g = CosyT;
                    dg = -y*SinyT;
                elseif strcmp(ModeType,'AShear')
                    g = SinyT;
                    dg = y*CosyT;
                end
                v2 = -1i*AngularFrequency(q)*g;
                epsilon4 = dg;
                epsilon6 = 1i*k*g;
                sigma4 = Mu*epsilon4;
                sigma6 = Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon4.*conj(sigma4)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*abs(v2).^2;
                PowerFlowDensity = -.5*real(sigma6.*conj(v2));
                PowerFlow = sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = deltaT*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = deltaT*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % ce1 (m/ms)
            end
            X{p}(:,5) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                PhaseVelocity = X{p}(q,4)*1e3;
                Attenuation = X{p}(q,6)*PhaseVelocity/(X{p}(q,1)*1e3); % Np/wavelength
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
                y = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
                yT = y*T;
                SinyT = sin(yT);
                CosyT = cos(yT);
                if  strcmp(ModeType,'SShear')
                    g = CosyT;
                    dg = -y*SinyT;
                elseif strcmp(ModeType,'AShear')
                    g = SinyT;
                    dg = y*CosyT;
                end
                v2 = -1i*AngularFrequency*g;
                epsilon4 = dg;
                epsilon6 = 1i*k*g;
                sigma4 = Mu*epsilon4;
                sigma6 = Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon4.*conj(sigma4)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*abs(v2).^2;
                PowerFlowDensity = -.5*real(sigma6.*conj(v2));
                PowerFlow = sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = deltaT*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = deltaT*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % ce1 (m/ms)
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