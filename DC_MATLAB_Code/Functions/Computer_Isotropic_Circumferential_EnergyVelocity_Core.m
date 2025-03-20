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
function [X,Counter] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(Multithreading,X,ModeType,ModeTotal,Counter,h,Material,Ro,Ri,SamplesR)
%#ok<*PFBNS>
%#ok<*GVMIS>
global Stop 
Stop = 0;
Ro2 = Ro^2;
Ri2 = Ri^2;
deltar = (Ro-Ri)/SamplesR;
r = (Ri:deltar:Ro)';
r2 = r.^2;
for p = 1:length(X)
    if  strcmp(ModeType,'Lamb')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                k = AngularFrequency(q)/PhaseVelocity(q);
                kL = AngularFrequency(q)/Material.LongitudinalVelocity; 
                kT = AngularFrequency(q)/Material.TransverseVelocity;
                kL2 = kL^2;
                kT2 = kT^2;
                kLr = kL*r;
                kTr = kT*r;
                kRo = k*Ro;
                kRo2 = kRo^2;
                kRo1 = 1/(kRo+1);
                kRo_1 = 1/(kRo-1);
                JkLr = besselj(kRo,kLr);
                JkTr = besselj(kRo,kTr);
                YkLr = bessely(kRo,kLr);
                YkTr = bessely(kRo,kTr);
                J_2kLr = besselj(kRo-2,kLr);
                J_2kTr = besselj(kRo-2,kTr);
                Y_2kLr = bessely(kRo-2,kLr);
                Y_2kTr = bessely(kRo-2,kTr);
                J2kLr = besselj(kRo+2,kLr);
                J2kTr = besselj(kRo+2,kTr);
                Y2kLr = bessely(kRo+2,kLr);
                Y2kTr = bessely(kRo+2,kTr);
                dJkLr = kL2*r/4.*(kRo_1*(J_2kLr+JkLr)-kRo1*(JkLr+J2kLr));
                dJkTr = kT2*r/4.*(kRo_1*(J_2kTr+JkTr)-kRo1*(JkTr+J2kTr));
                dYkLr = kL2*r/4.*(kRo_1*(Y_2kLr+YkLr)-kRo1*(YkLr+Y2kLr));
                dYkTr = kT2*r/4.*(kRo_1*(Y_2kTr+YkTr)-kRo1*(YkTr+Y2kTr));
                d2JkLr = kL2/4*(J_2kLr-2*JkLr+J2kLr);
                d2JkTr = kT2/4*(J_2kTr-2*JkTr+J2kTr);
                d2YkLr = kL2/4*(Y_2kLr-2*YkLr+Y2kLr);
                d2YkTr = kT2/4*(Y_2kTr-2*YkTr+Y2kTr);
                Z1 = [2i*Material.Mu*kRo*(dYkLr(1)/Ri-YkLr(1)/Ri2) -Material.Mu*(d2JkTr(1)-dJkTr(1)/Ri+kRo2*JkTr(1)/Ri2) -Material.Mu*(d2YkTr(1)-dYkTr(1)/Ri+kRo2*YkTr(1)/Ri2);(Material.Lambda+2*Material.Mu)*d2YkLr(end)+Material.Lambda*(dYkLr(end)/Ro-kRo2*YkLr(end)/Ro2) 2i*Material.Mu*kRo*(dJkTr(end)/Ro-JkTr(end)/Ro2) 2i*Material.Mu*kRo*(dYkTr(end)/Ro-YkTr(end)/Ro2);2i*Material.Mu*kRo*(dYkLr(end)/Ro-YkLr(end)/Ro2) -Material.Mu*(d2JkTr(end)-dJkTr(end)/Ro+kRo2*JkTr(end)/Ro2) -Material.Mu*(d2YkTr(end)-dYkTr(end)/Ro+kRo2*YkTr(end)/Ro2)];
                Z2 = [2i*Material.Mu*kRo*(dJkLr(1)/Ri-JkLr(1)/Ri2);(Material.Lambda+2*Material.Mu)*d2JkLr(end)+Material.Lambda*(dJkLr(end)/Ro-kRo2*JkLr(end)/Ro2);2i*Material.Mu*kRo*(dJkLr(end)/Ro-JkLr(end)/Ro2)];
                U = Z1\-Z2;
                f = JkLr+U(1)*YkLr; % L_in + L_out
                df = dJkLr+U(1)*dYkLr;
                d2f = d2JkLr+U(1)*d2YkLr;
                g = U(2)*JkTr+U(3)*YkTr; % SV_in + SV_out
                dg = U(2)*dJkTr+U(3)*dYkTr;
                d2g = U(2)*d2JkTr+U(3)*d2YkTr;
                Delta = Material.Lambda*(d2f+df./r-kRo2*f./r2);
                u2 = 1i*kRo*f./r-dg;
                u3 = df+1i*kRo*g./r;
                v2 = -1i*AngularFrequency(q)*u2;
                v3 = -1i*AngularFrequency(q)*u3;
                epsilon2 = (1i*kRo*u2+u3)./r;
                epsilon3 = d2f+1i*kRo*(dg./r-g./r2);
                epsilon4 = 1i*kRo*(df./r-f./r2+u3./r)-d2g-u2./r;
                sigma2 = Delta+2*Material.Mu*epsilon2;
                sigma3 = Delta+2*Material.Mu*epsilon3;
                sigma4 = Material.Mu*epsilon4;
                StrainEnergyDensity = .5*real(epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4));
                KineticEnergyDensity = .5*Material.Density*(abs(v2).^2+abs(v3).^2);
                PowerFlowDensity = -.5*real(sigma2.*conj(v2)+sigma4.*conj(v3));
                PowerFlow = sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = deltar*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = deltar*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % cetheta (m/ms)
            end
            X{p}(:,5) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                PhaseVelocity = X{p}(q,4)*1e3;
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                k = AngularFrequency/PhaseVelocity;
                kL = AngularFrequency/Material.LongitudinalVelocity;
                kT = AngularFrequency/Material.TransverseVelocity;
                kL2 = kL^2;
                kT2 = kT^2;
                kLr = kL*r;
                kTr = kT*r;
                kRo = k*Ro;
                kRo2 = kRo^2;
                kRo1 = 1/(kRo+1);
                kRo_1 = 1/(kRo-1);
                JkLr = besselj(kRo,kLr);
                JkTr = besselj(kRo,kTr);
                YkLr = bessely(kRo,kLr);
                YkTr = bessely(kRo,kTr);
                J_2kLr = besselj(kRo-2,kLr);
                J_2kTr = besselj(kRo-2,kTr);
                Y_2kLr = bessely(kRo-2,kLr);
                Y_2kTr = bessely(kRo-2,kTr);
                J2kLr = besselj(kRo+2,kLr);
                J2kTr = besselj(kRo+2,kTr);
                Y2kLr = bessely(kRo+2,kLr);
                Y2kTr = bessely(kRo+2,kTr);
                dJkLr = kL2*r/4.*(kRo_1*(J_2kLr+JkLr)-kRo1*(JkLr+J2kLr));
                dJkTr = kT2*r/4.*(kRo_1*(J_2kTr+JkTr)-kRo1*(JkTr+J2kTr));
                dYkLr = kL2*r/4.*(kRo_1*(Y_2kLr+YkLr)-kRo1*(YkLr+Y2kLr));
                dYkTr = kT2*r/4.*(kRo_1*(Y_2kTr+YkTr)-kRo1*(YkTr+Y2kTr));
                d2JkLr = kL2/4*(J_2kLr-2*JkLr+J2kLr);
                d2JkTr = kT2/4*(J_2kTr-2*JkTr+J2kTr);
                d2YkLr = kL2/4*(Y_2kLr-2*YkLr+Y2kLr);
                d2YkTr = kT2/4*(Y_2kTr-2*YkTr+Y2kTr);
                Z1(1,1) = 2i*Material.Mu*kRo*(dYkLr(1)/Ri-YkLr(1)/Ri2);
                Z1(1,2) = -Material.Mu*(d2JkTr(1)-dJkTr(1)/Ri+kRo2*JkTr(1)/Ri2);
                Z1(1,3) = -Material.Mu*(d2YkTr(1)-dYkTr(1)/Ri+kRo2*YkTr(1)/Ri2);
                Z1(2,1) = (Material.Lambda+2*Material.Mu)*d2YkLr(end)+Material.Lambda*(dYkLr(end)/Ro-kRo2*YkLr(end)/Ro2);
                Z1(2,2) = 2i*Material.Mu*kRo*(dJkTr(end)/Ro-JkTr(end)/Ro2);
                Z1(2,3) = 2i*Material.Mu*kRo*(dYkTr(end)/Ro-YkTr(end)/Ro2);
                Z1(3,1) = 2i*Material.Mu*kRo*(dYkLr(end)/Ro-YkLr(end)/Ro2);
                Z1(3,2) = -Material.Mu*(d2JkTr(end)-dJkTr(end)/Ro+kRo2*JkTr(end)/Ro2);
                Z1(3,3) = -Material.Mu*(d2YkTr(end)-dYkTr(end)/Ro+kRo2*YkTr(end)/Ro2);
                Z2(1,1) = 2i*Material.Mu*kRo*(dJkLr(1)/Ri-JkLr(1)/Ri2);
                Z2(2,1) = (Material.Lambda+2*Material.Mu)*d2JkLr(end)+Material.Lambda*(dJkLr(end)/Ro-kRo2*JkLr(end)/Ro2);
                Z2(3,1) = 2i*Material.Mu*kRo*(dJkLr(end)/Ro-JkLr(end)/Ro2);
                U = Z1\-Z2;
                f = JkLr+U(1)*YkLr; % L_in + L_out
                df = dJkLr+U(1)*dYkLr;
                d2f = d2JkLr+U(1)*d2YkLr;
                g = U(2)*JkTr+U(3)*YkTr; % SV_in + SV_out
                dg = U(2)*dJkTr+U(3)*dYkTr;
                d2g = U(2)*d2JkTr+U(3)*d2YkTr;
                Delta = Material.Lambda*(d2f+df./r-kRo2*f./r2);
                u2 = 1i*kRo*f./r-dg;
                u3 = df+1i*kRo*g./r;
                v2 = -1i*AngularFrequency*u2;
                v3 = -1i*AngularFrequency*u3;
                epsilon2 = (1i*kRo*u2+u3)./r;
                epsilon3 = d2f+1i*kRo*(dg./r-g./r2);
                epsilon4 = 1i*kRo*(df./r-f./r2+u3./r)-d2g-u2./r;
                sigma2 = Delta+2*Material.Mu*epsilon2;
                sigma3 = Delta+2*Material.Mu*epsilon3;
                sigma4 = Material.Mu*epsilon4;
                StrainEnergyDensity = .5*real(epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4));
                KineticEnergyDensity = .5*Material.Density*(abs(v2).^2+abs(v3).^2);
                PowerFlowDensity = -.5*real(sigma2.*conj(v2)+sigma4.*conj(v3));
                PowerFlow = sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = deltar*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = deltar*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % cetheta (m/ms)
            end
        end
    elseif contains(ModeType,'Shear')
        if  Multithreading
            EnergyVelocity = 0;
            PhaseVelocity = X{p}(:,4)*1e3;
            AngularFrequency = 2*pi*X{p}(:,1)*1e3;
            parfor q = 1:height(X{p})
                k = AngularFrequency(q)/PhaseVelocity(q);
                kT = AngularFrequency(q)/Material.TransverseVelocity;
                kTr = kT*r;
                kRo = k*Ro;
                J = besselj(kRo,kTr);
                Y = bessely(kRo,kTr);
                J_1 = besselj(kRo-1,kTr);
                Y_1 = bessely(kRo-1,kTr);
                J1 = besselj(kRo+1,kTr);
                Y1 = bessely(kRo+1,kTr);
                U = (Y_1(end)-Y1(end))\-(J_1(end)-J1(end));
                u1 = J+U*Y;
                v1 = -1i*AngularFrequency(q)*u1;
                epsilon5 = kT/2*(J_1-J1+U*(Y_1-Y1));
                epsilon6 = 1i*kRo./r.*u1;
                sigma5 = Material.Mu*epsilon5;
                sigma6 = Material.Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*abs(v1).^2;
                PowerFlowDensity = -.5*real(sigma6.*conj(v1));
                PowerFlow = sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = deltar*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = deltar*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                EnergyVelocity(q) = PowerFlow/TotalEnergy/1e3; % cetheta (m/ms)
            end
            X{p}(:,5) = EnergyVelocity;
        else
            for q = 1:height(X{p})
                PhaseVelocity = X{p}(q,4)*1e3;
                AngularFrequency = 2*pi*X{p}(q,1)*1e3;
                k = AngularFrequency/PhaseVelocity;
                kT = AngularFrequency/Material.TransverseVelocity;
                kTr = kT*r;
                kRo = k*Ro;
                J = besselj(kRo,kTr);
                Y = bessely(kRo,kTr);
                J_1 = besselj(kRo-1,kTr);
                Y_1 = bessely(kRo-1,kTr);
                J1 = besselj(kRo+1,kTr);
                Y1 = bessely(kRo+1,kTr);
                U = (Y_1(end)-Y1(end))\-(J_1(end)-J1(end));
                u1 = J+U*Y;
                v1 = -1i*AngularFrequency*u1;
                epsilon5 = kT/2*(J_1-J1+U*(Y_1-Y1));
                epsilon6 = 1i*kRo./r.*u1;
                sigma5 = Material.Mu*epsilon5;
                sigma6 = Material.Mu*epsilon6;
                StrainEnergyDensity = .5*real(epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6));
                KineticEnergyDensity = .5*Material.Density*abs(v1).^2;
                PowerFlowDensity = -.5*real(sigma6.*conj(v1));
                PowerFlow = sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end)); % PowerFlow = deltar*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
                TotalEnergy = sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/2; % TotalEnergy = deltar*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
                X{p}(q,5) = PowerFlow/TotalEnergy/1e3; % cetheta (m/ms)
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