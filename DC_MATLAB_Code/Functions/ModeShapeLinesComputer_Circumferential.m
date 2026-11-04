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
function [u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,rTotal,p,Frequency] = ModeShapeLinesComputer_Circumferential(Material,CLamb,CShear,Frequency,Mode,Layers,Radius,SamplesR,Phase)
%#ok<*AGROW>
uPhase=0;epsilonPhase=0;sigmaPhase=0;
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'SH')
    [PhaseVelocity,Direction,Frequency] = ExtractData(CLamb{p},Frequency);
elseif contains(Mode,'SH')
    [PhaseVelocity,Direction,Frequency] = ExtractData(CShear{p},Frequency);
end
for m = 1:Layers
    Lambda(m) = Material{m}.Lambda;
    Mu(m) = Material{m}.Mu;
    cL(m) = Material{m}.LongitudinalVelocity;
    cT(m) = Material{m}.TransverseVelocity;
end
rTotal = [];
AngularFrequency = Frequency*pi*2e3;
kRo = Direction*AngularFrequency/PhaseVelocity*Radius(end);
kRoRep = repmat(kRo,2,1,1);
if  ~contains(Mode,'SH')
    kRo2 = kRo^2;
    kRo1 = 1/(kRo+1);
    kRo_1 = 1/(kRo-1);
    for m = 1:Layers
        R = Radius(m:m+1);
        R2 = R.^2;
        kL = AngularFrequency/cL(m);
        kT = AngularFrequency/cT(m);
        kL2 = kL^2;
        kT2 = kT^2;
        kLr = kL*R;
        kTr = kT*R;
        JYkLr = [besselj(kRoRep,kLr) bessely(kRoRep,kLr)];
        JYkTr = [besselj(kRoRep,kTr) bessely(kRoRep,kTr)];
        JY_2kLr = [besselj(kRoRep-2,kLr) bessely(kRoRep-2,kLr)];
        JY_2kTr = [besselj(kRoRep-2,kTr) bessely(kRoRep-2,kTr)];
        JY2kLr = [besselj(kRoRep+2,kLr) bessely(kRoRep+2,kLr)];
        JY2kTr = [besselj(kRoRep+2,kTr) bessely(kRoRep+2,kTr)];
        dJYkLr = kL2*R/4.*[kRo_1*(JY_2kLr(:,1)+JYkLr(:,1))-kRo1*(JYkLr(:,1)+JY2kLr(:,1)) kRo_1*(JY_2kLr(:,2)+JYkLr(:,2))-kRo1*(JYkLr(:,2)+JY2kLr(:,2))];
        dJYkTr = kT2*R/4.*[kRo_1*(JY_2kTr(:,1)+JYkTr(:,1))-kRo1*(JYkTr(:,1)+JY2kTr(:,1)) kRo_1*(JY_2kTr(:,2)+JYkTr(:,2))-kRo1*(JYkTr(:,2)+JY2kTr(:,2))];
        d2JYkLr = kL2/4*[JY_2kLr(:,1)-2*JYkLr(:,1)+JY2kLr(:,1) JY_2kLr(:,2)-2*JYkLr(:,2)+JY2kLr(:,2)];
        d2JYkTr = kT2/4*[JY_2kTr(:,1)-2*JYkTr(:,1)+JY2kTr(:,1) JY_2kTr(:,2)-2*JYkTr(:,2)+JY2kTr(:,2)];
        L1([1 3],[1 3]) = 2*Mu(m)*d2JYkLr+Lambda(m)*(d2JYkLr+dJYkLr./R-kRo2*JYkLr./R2); % sigma33
        L1([1 3],[2 4]) = 2i*Mu(m)*kRo*(dJYkTr-JYkTr./R)./R;
        L1([2 4],[1 3]) = 2i*Mu(m)*kRo*(dJYkLr-JYkLr./R)./R; % sigma23
        L1([2 4],[2 4]) = Mu(m)*(dJYkTr./R-kRo2*JYkTr./R2-d2JYkTr);
        A{m,3}([1 3],[1 3]) = 1i*kRo*JYkLr./R; % u2
        A{m,3}([1 3],[2 4]) = -dJYkTr;
        A{m,3}([2 4],[1 3]) = dJYkLr; % u3
        A{m,3}([2 4],[2 4]) = 1i*kRo*JYkTr./R;
        A{m,1} = L1/A{m,3};
    end
    A{1,2} = A{1};
    for m = 2:Layers
        M0 = A{m}(1:2,1:2)-A{m-1,2}(3:4,3:4);
        M1 = A{m-1,2}(1:2,3:4)/M0;
        M2 = A{m}(3:4,1:2)/M0;
        A{m,2} = [A{m-1,2}(1:2,1:2)+M1*A{m-1,2}(3:4,1:2) -M1*A{m}(1:2,3:4);M2*A{m-1,2}(3:4,1:2) A{m}(3:4,3:4)-M2*A{m}(1:2,3:4)];
    end
    In0 = A{1,3}(1:2,1:2);
    Out0 = A{1,3}(1:2,3:4);
    OutN = A{end,3}(3:4,3:4);
    Z1 = [-A{end,2}(1:2,1:2)*In0 -A{end,2}(1:2,3:4)*OutN;A{end,2}(3:4,1:2)*In0 A{end,2}(3:4,3:4)*OutN];
    Z2 = [A{end,2}(1:2,1:2)*Out0;-A{end,2}(3:4,1:2)*Out0];
    RT = Z1\Z2(:,2);
    uInterfaces = Out0(:,2)+In0*RT(1:2);
    uInterfaces(:,Layers+1) = OutN*RT(3:4);
    for m = Layers:-1:2
        M0 = A{m}(1:2,1:2)-A{m-1,2}(3:4,3:4);
        uInterfaces(:,m) = M0\A{m-1,2}(3:4,1:2)*uInterfaces(:,1)-M0\A{m}(1:2,3:4)*uInterfaces(:,m+1);
    end
    for m = 1:Layers
        r = (Radius(m):(Radius(m+1)-Radius(m))/SamplesR:Radius(m+1))';
        r2 = r.^2;
        rTotal = [rTotal;r];
        s = (m-1)*length(r)+1:m*length(r);
        U = A{m,3}\[uInterfaces(:,m);uInterfaces(:,m+1)];
        kL = AngularFrequency/cL(m);
        kT = AngularFrequency/cT(m);
        kL2 = kL^2;
        kT2 = kT^2;
        kLr = kL*r;
        kTr = kT*r;
        kRoRep = repmat(kRo,length(r),1,1);
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
        dJkLr = kL2*r/4.*(kRo_1*(J_2kLr+JkLr)-kRo1*(JkLr+J2kLr));
        dJkTr = kT2*r/4.*(kRo_1*(J_2kTr+JkTr)-kRo1*(JkTr+J2kTr));
        dYkLr = kL2*r/4.*(kRo_1*(Y_2kLr+YkLr)-kRo1*(YkLr+Y2kLr));
        dYkTr = kT2*r/4.*(kRo_1*(Y_2kTr+YkTr)-kRo1*(YkTr+Y2kTr));
        d2JkLr = kL2/4*(J_2kLr-2*JkLr+J2kLr);
        d2JkTr = kT2/4*(J_2kTr-2*JkTr+J2kTr);
        d2YkLr = kL2/4*(Y_2kLr-2*YkLr+Y2kLr);
        d2YkTr = kT2/4*(Y_2kTr-2*YkTr+Y2kTr);
        f = U(1)*JkLr+U(3)*YkLr; % L_in + L_out
        df = U(1)*dJkLr+U(3)*dYkLr;
        d2f = U(1)*d2JkLr+U(3)*d2YkLr;
        g = U(2)*JkTr+U(4)*YkTr; % SV_in + SV_out
        dg = U(2)*dJkTr+U(4)*dYkTr;
        d2g = U(2)*d2JkTr+U(4)*d2YkTr;
        u(s,2) = 1i*kRo*f./r-dg; % u2
        u(s,3) = df+1i*kRo*g./r; % u3
        v(s,:) = -1i*AngularFrequency*u(s,:);
        epsilon(s,2) = (1i*kRo*u(s,2)+u(s,3))./r; % epsilon22
        epsilon(s,3) = d2f+1i*kRo*(dg./r-g./r2); % epsilon33
        epsilon(s,4) = 1i*kRo*(df./r-f./r2+u(s,3)./r)-d2g-u(s,2)./r; % epsilon23
        sigma(s,1:3) = Lambda(m)*(d2f+df./r-kRo2*f./r2)+2*Mu(m)*epsilon(s,1:3);
        sigma(s,4) = Mu(m)*epsilon(s,4); % sigma23
        KineticEnergyDensity(s,1) = Material{m}.Density*sum(abs(v(s,:)).^2,2)/2;
    end
    epsilon(:,6) = 0;
    sigma(:,6) = 0;
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
    for m = 1:Layers
        kT = AngularFrequency/cT(m);
        kTr = kT*Radius(m:m+1);
        JY_1 = [besselj(kRoRep-1,kTr) bessely(kRoRep-1,kTr)];
        JY1 = [besselj(kRoRep+1,kTr) bessely(kRoRep+1,kTr)];
        A{m,3} = kTr/kRo.*(JY_1+JY1); % multiplied by two
        A{m,1} = Mu(m)*kT*(JY_1-JY1)/A{m,3};
    end
    A{1,2} = A{1};
    for m = 2:Layers
        M0 = A{m}(1,1)-A{m-1,2}(2,2);
        M1 = A{m-1,2}(1,2)/M0;
        M2 = A{m}(2,1)/M0;
        A{m,2} = [A{m-1,2}(1,1)+M1*A{m-1,2}(2,1) -M1*A{m}(1,2);M2*A{m-1,2}(2,1) A{m}(2,2)-M2*A{m}(1,2)];
    end
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
    for m = 1:Layers
        r = (Radius(m):(Radius(m+1)-Radius(m))/SamplesR:Radius(m+1))';
        rTotal = [rTotal;r];
        s = (m-1)*length(r)+1:m*length(r);
        U = A{m,3}\[uInterfaces(:,m);uInterfaces(:,m+1)];
        kT = AngularFrequency/cT(m);
        kTr = kT*r;
        kRoRep = repmat(kRo,length(r),1,1);
        J_1 = besselj(kRoRep-1,kTr);
        Y_1 = bessely(kRoRep-1,kTr);
        J1 = besselj(kRoRep+1,kTr);
        Y1 = bessely(kRoRep+1,kTr);
        J = kTr/kRo/2.*(J_1+J1);
        Y = kTr/kRo/2.*(Y_1+Y1);
        u(s,1) = U(1)*J+U(2)*Y;
        v(s,1) = -1i*AngularFrequency*u(s,1);
        epsilon(s,5) = kT/2*(U(1)*(J_1-J1)+U(2)*(Y_1-Y1)); % epsilon13
        epsilon(s,6) = 1i*kRo./r.*u(s,1); % epsilon12
        sigma(s,5:6) = Mu(m)*epsilon(s,5:6);
        KineticEnergyDensity(s,1) = Material{m}.Density*abs(v(s,1)).^2/2;
    end
    u(:,3) = 0;
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
StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
PowerFlowDensity(:,3) = 0;
PowerFlow = trapz(rTotal,PowerFlowDensity(:,2));
% disp(['ce: ',num2str(2*PowerFlow/trapz(rTotal,StrainEnergyDensity+KineticEnergyDensity)),' m/s'])
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