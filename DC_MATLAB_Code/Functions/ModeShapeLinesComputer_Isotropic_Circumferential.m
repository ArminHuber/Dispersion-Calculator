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
function [u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,r,p,Frequency] = ModeShapeLinesComputer_Isotropic_Circumferential(Material,CLamb,CShear,Frequency,Mode,Ro,Ri,SamplesR,Phase)
uPhase=0;epsilonPhase=0;sigmaPhase=0;
r = (Ri:(Ro-Ri)/SamplesR:Ro)';
r2 = r.^2;
Ro2 = Ro^2;
Ri2 = Ri^2;
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'SH')
    [PhaseVelocity,Direction,Frequency] = ExtractData(CLamb{p},Frequency);
elseif contains(Mode,'SH')
    [PhaseVelocity,Direction,Frequency] = ExtractData(CShear{p},Frequency);
end
AngularFrequency = Frequency*pi*2e3;
kT = AngularFrequency/Material.TransverseVelocity;
kTr = kT*r;
kRo = Direction*AngularFrequency/PhaseVelocity*Ro;
if  ~contains(Mode,'SH')
    kL = AngularFrequency/Material.LongitudinalVelocity;
    kL2 = kL^2;
    kT2 = kT^2;
    kLr = kL*r;
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
    u(:,2) = 1i*kRo*f./r-dg;
    u(:,3) = df+1i*kRo*g./r;
    v = -1i*AngularFrequency*u;
    epsilon(:,2) = (1i*kRo*u(:,2)+u(:,3))./r;
    epsilon(:,3) = d2f+1i*kRo*(dg./r-g./r2);
    epsilon(:,4) = 1i*kRo*(df./r-f./r2+u(:,3)./r)-d2g-u(:,2)./r;
    epsilon(:,6) = 0;
    sigma(:,1:3) = Material.Lambda*(d2f+df./r-kRo2*f./r2)+2*Material.Mu*epsilon(:,1:3);
    sigma(:,4) = Material.Mu*epsilon(:,4);
    sigma(:,6) = 0;
    StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
    KineticEnergyDensity = Material.Density*sum(abs(v).^2,2)/2;
    PowerFlowDensity(:,2) = -real(sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)))/2;
    if  Phase
        uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
        sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
        epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
    end
    u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
    u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
    epsilon(:,2) = epsilon(:,2)*exp(-1i*angle(epsilon(2,2)));
    epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
    epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
    sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
    sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
    sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
    sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
elseif contains(Mode,'SH')
    J_1 = besselj(kRo-1,kTr);
    Y_1 = bessely(kRo-1,kTr);
    J1 = besselj(kRo+1,kTr);
    Y1 = bessely(kRo+1,kTr);
    J = kTr/kRo/2.*(J_1+J1);
    Y = kTr/kRo/2.*(Y_1+Y1);
    U = (Y1(end)-Y_1(end))\(J_1(end)-J1(end));
    u(:,1) = J+U*Y;
    u(:,3) = 0;
    v = -1i*AngularFrequency*u;
    epsilon(:,5) = kT/2*(J_1-J1+U*(Y_1-Y1));
    epsilon(:,6) = 1i*kRo./r.*u(:,1);
    sigma(:,5:6) = Material.Mu*epsilon(:,5:6);
    StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
    KineticEnergyDensity = Material.Density*abs(v(:,1)).^2/2;
    PowerFlowDensity(:,2) = -real(sigma(:,6).*conj(v(:,1)))/2;
    if  Phase
        uPhase = rad2deg(angle(u));
        sigmaPhase = rad2deg(angle(sigma));
        epsilonPhase = rad2deg(angle(epsilon));
    end
    u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
    epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
    epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
    sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
    sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
end
PowerFlowDensity(:,3) = 0;
PowerFlow = trapz(r,PowerFlowDensity(:,2));
% disp(['ce: ',num2str(2*PowerFlow/trapz(r,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
u = u/sqrt(PowerFlow);
epsilon = epsilon/sqrt(PowerFlow);
sigma = sigma/sqrt(PowerFlow);
if  real(u(2,1)) < 0
    u = -u;
end
if  real(epsilon(2,1)) < 0
    epsilon = -epsilon;
end
if  real(sigma(2,1)) < 0
    sigma = -sigma;
end
if  Phase
    uPhase(round(uPhase) == -180) = 180;
    sigmaPhase(round(sigmaPhase) == -180) = 180;
    epsilonPhase(round(epsilonPhase) == -180) = 180;
end
StrainEnergyDensity = StrainEnergyDensity/PowerFlow/2;
KineticEnergyDensity = KineticEnergyDensity/PowerFlow/2;
TotalEnergyDensity = StrainEnergyDensity+KineticEnergyDensity;
PowerFlowDensity = Direction*PowerFlowDensity/PowerFlow;
end
function [PhaseVelocity,Direction,Frequency] = ExtractData(X,Frequency)
    [Min,Max] = bounds(X(:,1));
    if  Frequency < Min || Frequency > Max
        errordlg(['Selected frequency outside frequency range! Select between ',num2str(Min),' and ',num2str(Max),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(X(:,1)-Frequency));
        Frequency = X(q);
    end
    PhaseVelocity = X(q,4)*1e3;
    if  X(q,5) > 0
        Direction = 1;
    else
        Direction = -1; % reverse propagation direction for backward propagating modes
    end
end