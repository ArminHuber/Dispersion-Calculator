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
function [CLamb,CShear] = Computer_Isotropic_Circumferential_EnergyVelocity_Multi(CLamb,CShear,Material,Layers,Radius,SamplesR)
%#ok<*AGROW>
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
SamplesR = ceil(SamplesR/Layers); % samples per layer
Height = Layers*(SamplesR+1);
for m = 1:Layers
    Lambda(m) = Material{m}.Lambda;
    Mu(m) = Material{m}.Mu;
    cL(m) = Material{m}.LongitudinalVelocity;
    cT(m) = Material{m}.TransverseVelocity;
end
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(CLamb{1})
    [CLamb,Counter] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(CLamb,'Lamb',ModeTotal,Counter,h,Material,Lambda,Mu,cL,cT,Layers,Radius,SamplesR,Height);
    if  Stop
        return
    end
end
if  ~isempty(CShear{1})
    [CShear,~] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(CShear,'Shear',ModeTotal,Counter,h,Material,Lambda,Mu,cL,cT,Layers,Radius,SamplesR,Height);
    if  Stop
        return
    end
end
close(h)
end
function [X,Counter] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(X,ModeType,ModeTotal,Counter,h,Material,Lambda,Mu,cL,cT,Layers,Radius,SamplesR,Height)
    global Stop 
    Stop = 0;
    for p = 1:length(X)
        PhaseVelocity = reshape(X{p}(:,4),1,1,[])*1e3;
        AngularFrequency = reshape(X{p}(:,1),1,1,[])*pi*2e3;
        k = AngularFrequency./PhaseVelocity;
        kRo = k*Radius(end);
        kRoRep = repmat(kRo,2,1,1);
        Length = length(k);
        KineticEnergyDensity = zeros(Height,1,Length);
        A = cell(0);
        rTotal = [];
        if  strcmp(ModeType,'Lamb')
            L1 = zeros(4,4,Length);
            u = zeros(Height,2,Length);
            v = zeros(Height,2,Length);
            epsilon = zeros(Height,3,Length);
            sigma = zeros(Height,3,Length);
            kRo2 = kRo.^2;
            kRo1 = 1/(kRo+1);
            kRo_1 = 1/(kRo-1);
            for m = 1:Layers
                R = Radius(m:m+1);
                R2 = R.^2;
                kL = AngularFrequency/cL(m);
                kT = AngularFrequency/cT(m);
                kL2 = kL.^2;
                kT2 = kT.^2;
                kLr = kL.*R;
                kTr = kT.*R;
                JYkLr = [besselj(kRoRep,kLr) bessely(kRoRep,kLr)];
                JYkTr = [besselj(kRoRep,kTr) bessely(kRoRep,kTr)];
                JY_2kLr = [besselj(kRoRep-2,kLr) bessely(kRoRep-2,kLr)];
                JY_2kTr = [besselj(kRoRep-2,kTr) bessely(kRoRep-2,kTr)];
                JY2kLr = [besselj(kRoRep+2,kLr) bessely(kRoRep+2,kLr)];
                JY2kTr = [besselj(kRoRep+2,kTr) bessely(kRoRep+2,kTr)];
                dJYkLr = kL2.*R/4.*[kRo_1.*(JY_2kLr(:,1,:)+JYkLr(:,1,:))-kRo1.*(JYkLr(:,1,:)+JY2kLr(:,1,:)) kRo_1.*(JY_2kLr(:,2,:)+JYkLr(:,2,:))-kRo1.*(JYkLr(:,2,:)+JY2kLr(:,2,:))];
                dJYkTr = kT2.*R/4.*[kRo_1.*(JY_2kTr(:,1,:)+JYkTr(:,1,:))-kRo1.*(JYkTr(:,1,:)+JY2kTr(:,1,:)) kRo_1.*(JY_2kTr(:,2,:)+JYkTr(:,2,:))-kRo1.*(JYkTr(:,2,:)+JY2kTr(:,2,:))];
                d2JYkLr = kL2/4.*[JY_2kLr(:,1,:)-2*JYkLr(:,1,:)+JY2kLr(:,1,:) JY_2kLr(:,2,:)-2*JYkLr(:,2,:)+JY2kLr(:,2,:)];
                d2JYkTr = kT2/4.*[JY_2kTr(:,1,:)-2*JYkTr(:,1,:)+JY2kTr(:,1,:) JY_2kTr(:,2,:)-2*JYkTr(:,2,:)+JY2kTr(:,2,:)];
                L1([1 3],[1 3],:) = 2*Mu(m)*d2JYkLr+Lambda(m)*(d2JYkLr+dJYkLr./R-kRo2.*JYkLr./R2); % sigma33
                L1([1 3],[2 4],:) = 2i*Mu(m)*kRo.*(dJYkTr-JYkTr./R)./R;
                L1([2 4],[1 3],:) = 2i*Mu(m)*kRo.*(dJYkLr-JYkLr./R)./R; % sigma23
                L1([2 4],[2 4],:) = Mu(m)*(dJYkTr./R-kRo2.*JYkTr./R2-d2JYkTr);
                A{m,3}([1 3],[1 3],:) = 1i*kRo.*JYkLr./R; % u2
                A{m,3}([1 3],[2 4],:) = -dJYkTr;
                A{m,3}([2 4],[1 3],:) = dJYkLr; % u3
                A{m,3}([2 4],[2 4],:) = 1i*kRo.*JYkTr./R;
                A{m,1} = pagemrdivide(L1,A{m,3});
            end
            A{1,2} = A{1};
            for m = 2:Layers
                M0 = A{m}(1:2,1:2,:)-A{m-1,2}(3:4,3:4,:);
                M1 = pagemrdivide(A{m-1,2}(1:2,3:4,:),M0);
                M2 = pagemrdivide(A{m}(3:4,1:2,:),M0);
                A{m,2} = [A{m-1,2}(1:2,1:2,:)+pagemtimes(M1,A{m-1,2}(3:4,1:2,:)) -pagemtimes(M1,A{m}(1:2,3:4,:));pagemtimes(M2,A{m-1,2}(3:4,1:2,:)) A{m}(3:4,3:4,:)-pagemtimes(M2,A{m}(1:2,3:4,:))];
            end
            In0 = A{1,3}(1:2,1:2,:);
            Out0 = A{1,3}(1:2,3:4,:);
            OutN = A{end,3}(3:4,3:4,:);
            Z1 = [pagemtimes(-A{end,2}(1:2,1:2,:),In0) pagemtimes(-A{end,2}(1:2,3:4,:),OutN);pagemtimes(A{end,2}(3:4,1:2,:),In0) pagemtimes(A{end,2}(3:4,3:4,:),OutN)];
            Z2 = [pagemtimes(A{end,2}(1:2,1:2,:),Out0);pagemtimes(-A{end,2}(3:4,1:2,:),Out0)];
            RT = pagemldivide(Z1,Z2(:,2,:));
            uInterfaces = Out0(:,2,:)+pagemtimes(In0,RT(1:2,1,:));
            uInterfaces(:,Layers+1,:) = pagemtimes(OutN,RT(3:4,1,:));
            for m = Layers:-1:2
                M0 = A{m}(1:2,1:2,:)-A{m-1,2}(3:4,3:4,:);
                uInterfaces(:,m,:) = pagemtimes(pagemldivide(M0,A{m-1,2}(3:4,1:2,:)),uInterfaces(:,1,:))-pagemtimes(pagemldivide(M0,A{m}(1:2,3:4,:)),uInterfaces(:,m+1,:));
            end
            for m = 1:Layers
                r = (Radius(m):(Radius(m+1)-Radius(m))/SamplesR:Radius(m+1))';
                r2 = r.^2;
                rTotal = [rTotal;r];
                s = (m-1)*length(r)+1:m*length(r);
                U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                kL = AngularFrequency/cL(m);
                kT = AngularFrequency/cT(m);
                kL2 = kL.^2;
                kT2 = kT.^2;
                kLr = kL.*r;
                kTr = kT.*r;
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
                dJkLr = kL2.*r/4.*(kRo_1.*(J_2kLr+JkLr)-kRo1.*(JkLr+J2kLr));
                dJkTr = kT2.*r/4.*(kRo_1.*(J_2kTr+JkTr)-kRo1.*(JkTr+J2kTr));
                dYkLr = kL2.*r/4.*(kRo_1.*(Y_2kLr+YkLr)-kRo1.*(YkLr+Y2kLr));
                dYkTr = kT2.*r/4.*(kRo_1.*(Y_2kTr+YkTr)-kRo1.*(YkTr+Y2kTr));
                d2JkLr = kL2/4.*(J_2kLr-2*JkLr+J2kLr);
                d2JkTr = kT2/4.*(J_2kTr-2*JkTr+J2kTr);
                d2YkLr = kL2/4.*(Y_2kLr-2*YkLr+Y2kLr);
                d2YkTr = kT2/4.*(Y_2kTr-2*YkTr+Y2kTr);
                f = U(1,1,:).*JkLr+U(3,1,:).*YkLr; % L_in + L_out
                df = U(1,1,:).*dJkLr+U(3,1,:).*dYkLr;
                d2f = U(1,1,:).*d2JkLr+U(3,1,:).*d2YkLr;
                g = U(2,1,:).*JkTr+U(4,1,:).*YkTr; % SV_in + SV_out
                dg = U(2,1,:).*dJkTr+U(4,1,:).*dYkTr;
                d2g = U(2,1,:).*d2JkTr+U(4,1,:).*d2YkTr;
                u(s,1,:) = 1i*kRo.*f./r-dg; % u2
                u(s,2,:) = df+1i*kRo.*g./r; % u3
                v(s,:,:) = -1i*AngularFrequency.*u(s,:,:);
                epsilon(s,1,:) = (1i*kRo.*u(s,1,:)+u(s,2,:))./r; % epsilon22
                epsilon(s,2,:) = d2f+1i*kRo.*(dg./r-g./r2); % epsilon33
                epsilon(s,3,:) = 1i*kRo.*(df./r-f./r2+u(s,2,:)./r)-d2g-u(s,1,:)./r; % epsilon23
                sigma(s,1:2,:) = Lambda(m)*(d2f+df./r-kRo2.*f./r2)+2*Mu(m)*epsilon(s,1:2,:);
                sigma(s,3,:) = Mu(m)*epsilon(s,3,:); % sigma23
                KineticEnergyDensity(s,1,:) = Material{m}.Density*sum(abs(v(s,:,:)).^2,2)/2;
            end
            PowerFlowDensity = -real(sigma(:,1,:).*conj(v(:,1,:))+sigma(:,3,:).*conj(v(:,2,:)))/2;
        elseif strcmp(ModeType,'Shear')
            u = zeros(Height,1,Length);
            v = zeros(Height,1,Length);
            epsilon = zeros(Height,2,Length);
            sigma = zeros(Height,2,Length);
            for m = 1:Layers
                kT = AngularFrequency/cT(m);
                kTr = kT.*Radius(m:m+1);
                JY_1 = [besselj(kRoRep-1,kTr) bessely(kRoRep-1,kTr)];
                JY1 = [besselj(kRoRep+1,kTr) bessely(kRoRep+1,kTr)];
                A{m,3} = kTr./kRo.*(JY_1+JY1); % multiplied by two
                A{m,1} = pagemrdivide(Mu(m)*kT.*(JY_1-JY1),A{m,3});
            end
            A{1,2} = A{1};
            for m = 2:Layers
                M0 = A{m}(1,1,:)-A{m-1,2}(2,2,:);
                M1 = A{m-1,2}(1,2,:)./M0;
                M2 = A{m}(2,1,:)./M0;
                A{m,2} = [A{m-1,2}(1,1,:)+M1.*A{m-1,2}(2,1,:) -M1.*A{m}(1,2,:);M2.*A{m-1,2}(2,1,:) A{m}(2,2,:)-M2.*A{m}(1,2,:)];
            end
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
            for m = 1:Layers
                r = (Radius(m):(Radius(m+1)-Radius(m))/SamplesR:Radius(m+1))';
                rTotal = [rTotal;r];
                s = (m-1)*length(r)+1:m*length(r);
                U = pagemldivide(A{m,3},[uInterfaces(:,m,:);uInterfaces(:,m+1,:)]);
                kT = AngularFrequency/cT(m);
                kTr = kT.*r;
                kRoRep = repmat(kRo,length(r),1,1);
                J_1 = besselj(kRoRep-1,kTr);
                Y_1 = bessely(kRoRep-1,kTr);
                J1 = besselj(kRoRep+1,kTr);
                Y1 = bessely(kRoRep+1,kTr);
                J = kTr./kRo/2.*(J_1+J1);
                Y = kTr./kRo/2.*(Y_1+Y1);
                u(s,1,:) = U(1,1,:).*J+U(2,1,:).*Y;
                v(s,1,:) = -1i*AngularFrequency.*u(s,1,:);
                epsilon(s,1,:) = kT/2.*(U(1,1,:).*(J_1-J1)+U(2,1,:).*(Y_1-Y1)); % epsilon13
                epsilon(s,2,:) = 1i*kRo./r.*u(s,1,:); % epsilon12
                sigma(s,:,:) = Mu(m)*epsilon(s,:,:);
                KineticEnergyDensity(s,1,:) = Material{m}.Density*abs(v(s,1,:)).^2/2;
            end
            PowerFlowDensity = -real(sigma(:,2,:).*conj(v(:,1,:)))/2;
        end
        StrainEnergyDensity = real(sum(epsilon.*conj(sigma),2))/2;
        PowerFlow = trapz(rTotal,squeeze(PowerFlowDensity));
        TotalEnergy = trapz(rTotal,squeeze(StrainEnergyDensity+KineticEnergyDensity))/2;
        X{p}(:,5) = fillmissing(filloutliers((PowerFlow./TotalEnergy)'/1e3,'spline','movmedian',5,'ThresholdFactor',1),'spline'); % cetheta (m/ms)
        Counter = Counter+1;
        waitbar(Counter/ModeTotal,h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',Counter,ModeTotal,100*Counter/ModeTotal,toc))
        if  Stop
            close(h)
            return
        end
    end
end
% InN = A{end,3}(3:4,1:2,:);
% OutN = A{end,3}(3:4,3:4,:);
% In0 = A{1,3}(1:2,1:2,:);
% Z1 = [pagemtimes(-A{end,2}(1:2,1:2,:),OutN) pagemtimes(-A{end,2}(1:2,3:4,:),In0);pagemtimes(A{end,2}(3:4,1:2,:),OutN) pagemtimes(A{end,2}(3:4,3:4,:),In0)];
% Z2 = [pagemtimes(A{end,2}(1:2,1:2,:),InN);pagemtimes(-A{end,2}(3:4,1:2,:),InN)];
% RT = pagemldivide(Z1,Z2(:,2,:));
% uInterfaces = InN(:,2,:)+pagemtimes(OutN,RT(1:2,1,:));
% uInterfaces(:,Layers+1,:) = pagemtimes(In0,RT(3:4,1,:));

% InN = A{end,3}(2,1,:);
% OutN = A{end,3}(2,2,:);
% In0 = A{1,3}(1,1,:);
% Z1 = [A{end,2}(1,1,:).*OutN -A{end,2}(1,2,:).*In0;A{end,2}(2,1,:).*OutN A{end,2}(2,2,:).*In0]; % changed sign in Z1(1,1)!
% Z2 = [A{end,2}(1,1,:).*InN;-A{end,2}(2,1,:).*InN];
% RT = pagemldivide(Z1,Z2);
% uInterfaces = InN+OutN.*RT(1,1,:);
% uInterfaces(:,Layers+1,:) = In0.*RT(2,1,:);