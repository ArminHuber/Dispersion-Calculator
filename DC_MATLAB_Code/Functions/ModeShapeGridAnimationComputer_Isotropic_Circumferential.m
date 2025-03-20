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
function [Time,u,x1,r,p] = ModeShapeGridAnimationComputer_Isotropic_Circumferential(Material,CLamb,CShear,CycleDuration,FrameRate,Frequency,Length,Mode,Ro,Ri,SamplesX1,SamplesR)
%#ok<*AGROW>
r = (Ri:(Ro-Ri)/SamplesR:Ro)';
Ro2 = Ro^2;
Ri2 = Ri^2;
Time = (0:1/(Frequency*1e3*CycleDuration*FrameRate):1/(Frequency*1e3))';
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'SH')
    q = find(CLamb{p}(:,1) == Frequency);
    if  isempty(q)
        if  Frequency > ceil(max(CLamb{p}(:,1))) || Frequency < min(CLamb{p}(:,1))
            errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(CLamb{p}(:,1))),' and ',num2str(ceil(max(CLamb{p}(:,1)))),' kHz.'],'Error');
            return
        else
            [~,q] = min(abs(CLamb{p}(:,1)-Frequency));
            Frequency = CLamb{p}(q,1);
        end
    end
    PhaseVelocity = CLamb{p}(q,4)*1e3;
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
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
    g = U(2)*JkTr+U(3)*YkTr; % SV_in + SV_out
    dg = U(2)*dJkTr+U(3)*dYkTr;
    u2 = 1i*kRo*f./r-dg; % uq
    u3 = df+1i*kRo*g./r; % ur
    E = exp(1i*(k*x1-AngularFrequency*Time));
    for i = 1:length(x1)
        for j = 1:length(Time)
            u{j,i}(:,1) = u2*E(j,i);
            u{j,i}(:,2) = u3*E(j,i);
        end
    end
elseif contains(Mode,'SH')
    q = find(CShear{p}(:,1) == Frequency);
    if  isempty(q)
        if  Frequency > ceil(max(CShear{p}(:,1))) || Frequency < min(CShear{p}(:,1))
            errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(CShear{p}(:,1))),' and ',num2str(ceil(max(CShear{p}(:,1)))),' kHz.'],'Error');
            return
        else
            [~,q] = min(abs(CShear{p}(:,1)-Frequency));
            Frequency = CShear{p}(q,1);
        end
    end
    PhaseVelocity = CShear{p}(q,4)*1e3;
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
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
    u1 = J+U*Y; % uz
    E = exp(1i*(k*x1-AngularFrequency*Time));
    for i = 1:length(x1)
        for j = 1:length(Time)
            u{j,i}(:,1) = u1*E(j,i);
            u{j,i}(:,2) = 0;
        end
    end
end