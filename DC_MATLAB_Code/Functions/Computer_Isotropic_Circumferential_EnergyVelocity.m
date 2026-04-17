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
function [CLamb,CShear] = Computer_Isotropic_Circumferential_EnergyVelocity(CLamb,CShear,Material,Ro,Ri,SamplesR)
%#ok<*GVMIS>
global Stop
Stop = 0;
ModeTotal = 0;
if  ~isempty(CLamb{1})
    ModeTotal = ModeTotal+length(CLamb);
end
if  ~isempty(CShear{1})
    ModeTotal = ModeTotal+length(CShear);
end
if  ModeTotal == 0
    return
end
Ro2 = Ro^2;
Ri2 = Ri^2;
deltar = (Ro-Ri)/SamplesR;
r = (Ri:deltar:Ro)';
r2 = r.^2;
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(CLamb{1})
    [CLamb,Counter] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(CLamb,'Lamb',ModeTotal,Counter,h,Material,r,r2,Ro,Ro2,Ri,Ri2);
    if  Stop
        return
    end
end
if  ~isempty(CShear{1})
    [CShear,~] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(CShear,'Shear',ModeTotal,Counter,h,Material,r,r2,Ro,Ro2,Ri,Ri2);
    if  Stop
        return
    end
end
close(h)
end
function [X,Counter] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(X,ModeType,ModeTotal,Counter,h,Material,r,r2,Ro,Ro2,Ri,Ri2)
    global Stop 
    Stop = 0;
    for p = 1:length(X)
        PhaseVelocity = reshape(X{p}(:,4),1,1,[])*1e3;
        AngularFrequency = reshape(X{p}(:,1),1,1,[])*pi*2e3;
        k = AngularFrequency./PhaseVelocity;
        kT = AngularFrequency/Material.TransverseVelocity;
        kTr = kT.*r;
        kRo = k*Ro;
        kRoRep = repmat(kRo,length(r),1,1);
        if  strcmp(ModeType,'Lamb')
            kL = AngularFrequency/Material.LongitudinalVelocity;
            kL2 = kL.^2;
            kT2 = kT.^2;
            kLr = kL.*r;
            kRo2 = kRo.^2;
            kRo1 = 1/(kRo+1);
            kRo_1 = 1/(kRo-1);
            JkLr = besselj(kRoRep,kLr);
            JkTr = besselj(kRoRep,kTr);
            YkLr = bessely(kRoRep,kLr);
            YkTr = bessely(kRoRep,kTr);
            J_2kLr = besselj(kRoRep-2,kLr);
            J_2kTr = besselj(kRoRep-2,kTr);
            Y_2kLr = bessely(kRoRep-2,kLr);
            Y_2kTr = bessely(kRoRep-2,kTr);
            J2kLr = besselj(kRoRep+2,kLr);
            J2kTr = besselj(kRoRep+2,kTr);
            Y2kLr = bessely(kRoRep+2,kLr);
            Y2kTr = bessely(kRoRep+2,kTr);
            dJkLr = kL2.*r/4.*(kRo_1.*(J_2kLr+JkLr)-kRo1.*(JkLr+J2kLr));
            dJkTr = kT2.*r/4.*(kRo_1.*(J_2kTr+JkTr)-kRo1.*(JkTr+J2kTr));
            dYkLr = kL2.*r/4.*(kRo_1.*(Y_2kLr+YkLr)-kRo1.*(YkLr+Y2kLr));
            dYkTr = kT2.*r/4.*(kRo_1.*(Y_2kTr+YkTr)-kRo1.*(YkTr+Y2kTr));
            d2JkLr = kL2/4.*(J_2kLr-2*JkLr+J2kLr);
            d2JkTr = kT2/4.*(J_2kTr-2*JkTr+J2kTr);
            d2YkLr = kL2/4.*(Y_2kLr-2*YkLr+Y2kLr);
            d2YkTr = kT2/4.*(Y_2kTr-2*YkTr+Y2kTr);
            Z1 = [2i*Material.Mu*kRo.*(dYkLr(1,1,:)/Ri-YkLr(1,1,:)/Ri2) -Material.Mu*(d2JkTr(1,1,:)-dJkTr(1,1,:)/Ri+kRo2.*JkTr(1,1,:)/Ri2) -Material.Mu*(d2YkTr(1,1,:)-dYkTr(1,1,:)/Ri+kRo2.*YkTr(1,1,:)/Ri2);(Material.Lambda+2*Material.Mu)*d2YkLr(end,1,:)+Material.Lambda*(dYkLr(end,1,:)/Ro-kRo2.*YkLr(end,1,:)/Ro2) 2i*Material.Mu*kRo.*(dJkTr(end,1,:)/Ro-JkTr(end,1,:)/Ro2) 2i*Material.Mu*kRo.*(dYkTr(end,1,:)/Ro-YkTr(end,1,:)/Ro2);2i*Material.Mu*kRo.*(dYkLr(end,1,:)/Ro-YkLr(end,1,:)/Ro2) -Material.Mu*(d2JkTr(end,1,:)-dJkTr(end,1,:)/Ro+kRo2.*JkTr(end,1,:)/Ro2) -Material.Mu*(d2YkTr(end,1,:)-dYkTr(end,1,:)/Ro+kRo2.*YkTr(end,1,:)/Ro2)];
            Z2 = [2i*Material.Mu*kRo.*(dJkLr(1,1,:)/Ri-JkLr(1,1,:)/Ri2);(Material.Lambda+2*Material.Mu)*d2JkLr(end,1,:)+Material.Lambda*(dJkLr(end,1,:)/Ro-kRo2.*JkLr(end,1,:)/Ro2);2i*Material.Mu*kRo.*(dJkLr(end,1,:)/Ro-JkLr(end,1,:)/Ro2)];
            U = -pagemldivide(Z1,Z2);
            f = JkLr+U(1,1,:).*YkLr; % L_in + L_out
            df = dJkLr+U(1,1,:).*dYkLr;
            d2f = d2JkLr+U(1,1,:).*d2YkLr;
            g = U(2,1,:).*JkTr+U(3,1,:).*YkTr; % SV_in + SV_out
            dg = U(2,1,:).*dJkTr+U(3,1,:).*dYkTr;
            d2g = U(2,1,:).*d2JkTr+U(3,1,:).*d2YkTr;
            Delta = Material.Lambda*(d2f+df./r-kRo2.*f./r2);
            u2 = 1i*kRo.*f./r-dg;
            u3 = df+1i*kRo.*g./r;
            v2 = -1i*AngularFrequency.*u2;
            v3 = -1i*AngularFrequency.*u3;
            epsilon2 = (1i*kRo.*u2+u3)./r;
            epsilon3 = d2f+1i*kRo.*(dg./r-g./r2);
            epsilon4 = 1i*kRo.*(df./r-f./r2+u3./r)-d2g-u2./r;
            sigma2 = Delta+2*Material.Mu*epsilon2;
            sigma3 = Delta+2*Material.Mu*epsilon3;
            sigma4 = Material.Mu*epsilon4;
            StrainEnergyDensity = real(epsilon2.*conj(sigma2)+epsilon3.*conj(sigma3)+epsilon4.*conj(sigma4))/2;
            KineticEnergyDensity = Material.Density*(abs(v2).^2+abs(v3).^2)/2;
            PowerFlowDensity = reshape(-real(sigma2.*conj(v2)+sigma4.*conj(v3))/2,length(r),length(k));
        elseif strcmp(ModeType,'Shear')
            J_1 = besselj(kRoRep-1,kTr);
            Y_1 = bessely(kRoRep-1,kTr);
            J1 = besselj(kRoRep+1,kTr);
            Y1 = bessely(kRoRep+1,kTr);
            J = kTr./kRo/2.*(J_1+J1);
            Y = kTr./kRo/2.*(Y_1+Y1);
            U = (Y1(end,1,:)-Y_1(end,1,:)).\(J_1(end,1,:)-J1(end,1,:));
            u1 = J+U.*Y;
            v1 = -1i*AngularFrequency.*u1;
            epsilon5 = kT/2.*(J_1-J1+U.*(Y_1-Y1));
            epsilon6 = 1i*kRo./r.*u1;
            sigma5 = Material.Mu*epsilon5;
            sigma6 = Material.Mu*epsilon6;
            StrainEnergyDensity = real(epsilon5.*conj(sigma5)+epsilon6.*conj(sigma6))/2;
            KineticEnergyDensity = Material.Density*abs(v1).^2/2;
            PowerFlowDensity = reshape(-real(sigma6.*conj(v1))/2,length(r),length(k));
        end
        PowerFlow = sum(PowerFlowDensity(1:end-1,:)+PowerFlowDensity(2:end,:)); % PowerFlow = deltar*sum(PowerFlowDensity(1:end-1)+PowerFlowDensity(2:end))/2;
        TotalEnergy = sum(StrainEnergyDensity(1:end-1,:)+StrainEnergyDensity(2:end,:)+KineticEnergyDensity(1:end-1,:)+KineticEnergyDensity(2:end,:))/2; % TotalEnergy = deltar*sum(StrainEnergyDensity(1:end-1)+StrainEnergyDensity(2:end)+KineticEnergyDensity(1:end-1)+KineticEnergyDensity(2:end))/4;
        X{p}(:,5) = fillmissing(filloutliers((PowerFlow./TotalEnergy)'/1e3,'spline','movmedian',5,'ThresholdFactor',1),'spline'); % cetheta (m/ms)
        Counter = Counter+1;
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
        if  Stop
            close(h)
            return
        end
    end
end