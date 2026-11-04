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
function [MatrixMethodLimit,EvanescenceLimit,XC0,XCSH0] = MatrixMethodLimitComputer_Circumferential(Material,Layers,Radius,FrequencyRange)
Resolution = 1e-6; % (m/s)
InstabilityLimitFactor = 3; % multiples of the phase velocity limit at which TMM becomes unstable; below this multiple, SMM is used instead of TMM
%#ok<*AGROW>

% find C0 at low frequency
AngularFrequency = pi*FrequencyRange(1)*2e3;
for m = 1:Layers
    Lambda(m) = Material{m}.Lambda;
    Mu(m) = Material{m}.Mu;
    kL(m) = AngularFrequency/Material{m}.LongitudinalVelocity;
    kT(m) = AngularFrequency/Material{m}.TransverseVelocity;
    LongitudinalVelocity(m) = Material{m}.LongitudinalVelocity;
end
kL2 = kL.^2;
kT2 = kT.^2;
EvanescenceLimit = max(LongitudinalVelocity);
SweepRange = .1:.1:200;
Wavenumber = AngularFrequency./SweepRange;
Y = Computer_C(Wavenumber,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius);
XRough = [];
for i = 2:length(SweepRange)-1
    if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
        XRough = SweepRange(i);
        break
    end
end
% figure,plot(SweepRange,20*log10(abs(Y)))
if  isempty(XRough)
    XC0 = 0;
else
    XC0 = Converger('Lamb',XRough,Resolution,AngularFrequency,SweepRange,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius);
    if  isempty(XC0)
        XC0 = 0;
    end
end

% find CSH0 at low frequency
SweepRange = XC0:10:EvanescenceLimit;
Wavenumber = AngularFrequency./SweepRange;
Y = Computer_CSH(Wavenumber,Mu,kT,Layers,Radius);
for i = 2:length(SweepRange)-1
    if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
        XRough = SweepRange(i);
        break
    end
end
XCSH0 = Converger('SH',XRough,Resolution,AngularFrequency,SweepRange,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius);

% find instability limit at high frequency
AngularFrequency = pi*FrequencyRange(end)*2e3;
for m = 1:Layers
    kL(m) = AngularFrequency/Material{m}.LongitudinalVelocity;
    kT(m) = AngularFrequency/Material{m}.TransverseVelocity;
end
kL2 = kL.^2;
kT2 = kT.^2;
SweepRange = 0:EvanescenceLimit;
Wavenumber = AngularFrequency./SweepRange;
Y = Computer_C(Wavenumber,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius);
% figure,plot(SweepRange,20*log10(abs(Y)))

% SweepRange = 0:1e4;
% Wavenumber = AngularFrequency./SweepRange;
% YC = Computer_C(Wavenumber,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius);
% YCSH = Computer_CSH(Wavenumber,Mu,kT,Layers,Radius);
% figure('name','Lamb'),plot(SweepRange,20*log10(abs(YC)))
% figure('name','SH'),plot(SweepRange,20*log10(abs(YCSH)))

XRough = [];
for i = length(SweepRange)-1:-1:2
    if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1))
        XRough(end+1) = SweepRange(i);
    end
end
Diff = abs(diff(XRough));
InstabilityLimit = EvanescenceLimit/4;
for i = 3:length(Diff)
    if  all(Diff(i-2:i) <= 10*(SweepRange(2)-SweepRange(1)))
        InstabilityLimit = XRough(i-2);
        break
    end
end
MatrixMethodLimit = InstabilityLimitFactor*InstabilityLimit*FrequencyRange/FrequencyRange(end);
MatrixMethodLimit(MatrixMethodLimit > EvanescenceLimit) = EvanescenceLimit;
end
function X = Converger(ModeFamily,XRough,Resolution,AngularFrequency,SweepRange,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius)
    X = [];
    Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Resolution));
    for j = 1:length(XRough)
        PhaseVelocity = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            PhaseVelocity = PhaseVelocity(1):(PhaseVelocity(end)-PhaseVelocity(1))/4:PhaseVelocity(end);
            Wavenumber = AngularFrequency./PhaseVelocity;
            if  strcmp(ModeFamily,'Lamb')
                Y = Computer_C(Wavenumber,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius);
            elseif strcmp(ModeFamily,'SH')
                Y = Computer_CSH(Wavenumber,Mu,kT,Layers,Radius);
            end
            for i = 2:length(PhaseVelocity)-1
                if  sign(Y(i-1)) ~= sign(Y(i+1))
                    if  o == Bisections
                        X(end+1) = PhaseVelocity(i);
                    end
                    PhaseVelocity = [PhaseVelocity(i-1) PhaseVelocity(i+1)];
                    break
                end
            end
        end
    end
end
function Y = Computer_C(k,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius)
    Size = size(k);
    k = reshape(k,1,1,[]);
    kRo = k*Radius(end);
    kRoRep = repmat(kRo,2,1,1);
    kRo2 = kRo.^2;
    kRo1 = 1/(kRo+1);
    kRo_1 = 1/(kRo-1);
    Length = length(k);
    Y = NaN(Size);
    for m = 1:Layers
        R = Radius(m:m+1);
        R2 = R.^2;
        kLrRep = repmat(kL(m)*R,1,1,Length);
        kTrRep = repmat(kT(m)*R,1,1,Length);
        JYkLr = [besselj(kRoRep,kLrRep) bessely(kRoRep,kLrRep)];
        JYkTr = [besselj(kRoRep,kTrRep) bessely(kRoRep,kTrRep)];
        JY_2kLr = [besselj(kRoRep-2,kLrRep) bessely(kRoRep-2,kLrRep)];
        JY_2kTr = [besselj(kRoRep-2,kTrRep) bessely(kRoRep-2,kTrRep)];
        JY2kLr = [besselj(kRoRep+2,kLrRep) bessely(kRoRep+2,kLrRep)];
        JY2kTr = [besselj(kRoRep+2,kTrRep) bessely(kRoRep+2,kTrRep)];
        dJYkLr = kL2(m)*R/4.*[kRo_1.*(JY_2kLr(:,1,:)+JYkLr(:,1,:))-kRo1.*(JYkLr(:,1,:)+JY2kLr(:,1,:)) kRo_1.*(JY_2kLr(:,2,:)+JYkLr(:,2,:))-kRo1.*(JYkLr(:,2,:)+JY2kLr(:,2,:))];
        dJYkTr = kT2(m)*R/4.*[kRo_1.*(JY_2kTr(:,1,:)+JYkTr(:,1,:))-kRo1.*(JYkTr(:,1,:)+JY2kTr(:,1,:)) kRo_1.*(JY_2kTr(:,2,:)+JYkTr(:,2,:))-kRo1.*(JYkTr(:,2,:)+JY2kTr(:,2,:))];
        d2JYkLr = kL2(m)/4*[JY_2kLr(:,1,:)-2*JYkLr(:,1,:)+JY2kLr(:,1,:) JY_2kLr(:,2,:)-2*JYkLr(:,2,:)+JY2kLr(:,2,:)];
        d2JYkTr = kT2(m)/4*[JY_2kTr(:,1,:)-2*JYkTr(:,1,:)+JY2kTr(:,1,:) JY_2kTr(:,2,:)-2*JYkTr(:,2,:)+JY2kTr(:,2,:)];
        L1(1,[1 3],:) = 1i*kRo.*JYkLr(2,:,:)./R(2); % u2
        L1(1,[2 4],:) = -dJYkTr(2,:,:);
        L1(2,[1 3],:) = dJYkLr(2,:,:); % u3
        L1(2,[2 4],:) = 1i*kRo.*JYkTr(2,:,:)./R(2);
        L1(3,[1 3],:) = 2*Mu(m)*d2JYkLr(2,:,:)+Lambda(m)*(d2JYkLr(2,:,:)+dJYkLr(2,:,:)./R(2)-kRo2.*JYkLr(2,:,:)./R2(2)); % sigma33
        L1(3,[2 4],:) = 2i*Mu(m)*kRo.*(dJYkTr(2,:,:)-JYkTr(2,:,:)./R(2))./R(2);
        L1(4,[1 3],:) = 2i*Mu(m)*kRo.*(dJYkLr(2,:,:)-JYkLr(2,:,:)./R(2))./R(2); % sigma23
        L1(4,[2 4],:) = Mu(m)*(dJYkTr(2,:,:)./R(2)-kRo2.*JYkTr(2,:,:)./R2(2)-d2JYkTr(2,:,:));
        L2(1,[1 3],:) = 1i*kRo.*JYkLr(1,:,:)./R(1); % u2
        L2(1,[2 4],:) = -dJYkTr(1,:,:);
        L2(2,[1 3],:) = dJYkLr(1,:,:); % u3
        L2(2,[2 4],:) = 1i*kRo.*JYkTr(1,:,:)./R(1);
        L2(3,[1 3],:) = 2*Mu(m)*d2JYkLr(1,:,:)+Lambda(m)*(d2JYkLr(1,:,:)+dJYkLr(1,:,:)./R(1)-kRo2.*JYkLr(1,:,:)./R2(1)); % sigma33
        L2(3,[2 4],:) = 2i*Mu(m)*kRo.*(dJYkTr(1,:,:)-JYkTr(1,:,:)./R(1))./R(1);
        L2(4,[1 3],:) = 2i*Mu(m)*kRo.*(dJYkLr(1,:,:)-JYkLr(1,:,:)./R(1))./R(1); % sigma23
        L2(4,[2 4],:) = Mu(m)*(dJYkTr(1,:,:)./R(1)-kRo2.*JYkTr(1,:,:)./R2(1)-d2JYkTr(1,:,:));
        L{m} = pagemrdivide(L1,L2);
    end
    M = L{end};
    for m = Layers-1:-1:1
        M = pagemtimes(M,L{m});
    end
    for j = 1:Length
        Y(j) = det(M(3:4,1:2,j));
    end
    Y = reshape(real(Y),Size);
end
function Y = Computer_CSH(k,Mu,kT,Layers,Radius)
    Size = size(k);
    k = reshape(k,1,1,[]);
    kRo = k*Radius(end);
    kRoRep = repmat(kRo,2,1,1);
    Length = length(k);
    for m = 1:Layers
        kTr = kT(m)*Radius(m:m+1);
        kTrRep = repmat(kTr,1,1,Length);
        JY_1 = [besselj(kRoRep-1,kTrRep) bessely(kRoRep-1,kTrRep)];
        JY1 = [besselj(kRoRep+1,kTrRep) bessely(kRoRep+1,kTrRep)];
        JY = kTr./kRo.*[JY_1(:,1,:)+JY1(:,1,:) JY_1(:,2,:)+JY1(:,2,:)]; % multiplied by two
        L1(1,:,:) = JY(2,:,:); % u1
        L1(2,:,:) = Mu(m)*kT(m)*(JY_1(2,:,:)-JY1(2,:,:)); % sigma13
        L2(1,:,:) = JY(1,:,:); % u1
        L2(2,:,:) = Mu(m)*kT(m)*(JY_1(1,:,:)-JY1(1,:,:)); % sigma13
        L{m} = pagemrdivide(L1,L2);
    end
    M = L{end};
    for m = Layers-1:-1:1
        M = pagemtimes(M,L{m});
    end
    Y = reshape(real(M(2,1,:)),Size);
end