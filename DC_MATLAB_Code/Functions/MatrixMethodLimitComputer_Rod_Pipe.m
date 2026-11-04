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
function [MatrixMethodLimit,EvanescenceLimit,XL01,XF11,XT01] = MatrixMethodLimitComputer_Rod_Pipe(Geometry,Material,Layers,Radius,FrequencyRange,InnerFluid,ToggleInnerFluid)
Resolution = 1e-6; % (m/s)
InstabilityLimitFactor = 3; % multiples of the phase velocity limit at which TMM becomes unstable; below this multiple, SMM is used instead of TMM
%#ok<*AGROW>

% find L(0,1) at low frequency
cFi2 = InnerFluid.Velocity^2;
AngularFrequency = pi*FrequencyRange(1)*2e3;
AngularFrequency2 = AngularFrequency^2;
kFi2 = AngularFrequency2/cFi2;
gFi = InnerFluid.Density*AngularFrequency2;
for m = 1:Layers
    Lambda(m) = Material{m}.Lambda;
    Mu(m) = Material{m}.Mu;
    kL2(m) = AngularFrequency2/Material{m}.LongitudinalVelocity^2;
    kT2(m) = AngularFrequency2/Material{m}.TransverseVelocity^2;
    LongitudinalVelocity(m) = Material{m}.LongitudinalVelocity;
end
EvanescenceLimit = max(LongitudinalVelocity);
if  Layers == 1 && ~ToggleInnerFluid
    XL01 = Material{1}.CylinderVelocity;
else
    SweepRangeStep = 10;
    SweepRange = 50.1:SweepRangeStep:25e3;
    if  strcmp(Geometry,'Rod')
        SweepRange(SweepRange > Material{1}.TransverseVelocity-SweepRangeStep/2 & SweepRange < Material{1}.TransverseVelocity+SweepRangeStep/2) = NaN;
    end
    Wavenumber = AngularFrequency./SweepRange;
    Y = Computer_L(Geometry,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
    XRough = [];
    for i = 2:length(SweepRange)-1
        if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
            XRough(end+1) = SweepRange(i);
            if  (ToggleInnerFluid && length(XRough) == 2) || (~ToggleInnerFluid && isscalar(XRough))
                break
            end
        end
    end
% figure,plot(SweepRange,20*log10(abs(Y)))
    XL01 = Converger('L',Geometry,XRough,Resolution,AngularFrequency,SweepRange,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
end

% find F(1,1) at low frequency
SweepRange = .1:.1:200;
Wavenumber = AngularFrequency./SweepRange;
Y = Computer_F(Geometry,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
XRough = [];
for i = 2:length(SweepRange)-1
    if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
        XRough = SweepRange(i);
        break
    end
end
% figure,plot(SweepRange,20*log10(abs(Y)))
if  isempty(XRough)
    XF11 = 0;
else
    XF11 = Converger('F',Geometry,XRough,Resolution,AngularFrequency,SweepRange,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
    if  isempty(XF11)
        XF11 = 0;
    end
end

% find T(0,1) at low frequency
if  Layers == 1
    XT01 = Material{1}.TransverseVelocity;
else
    SweepRange = XF11:10:EvanescenceLimit;
    if  strcmp(Geometry,'Rod')
        SweepRange(SweepRange > Material{1}.TransverseVelocity-SweepRangeStep/2 & SweepRange < Material{1}.TransverseVelocity+SweepRangeStep/2) = NaN;
    end
    Wavenumber = AngularFrequency./SweepRange;
    Y = Computer_T(Geometry,Wavenumber,Mu,kT2,Layers,Radius);
    for i = 2:length(SweepRange)-1
        if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
            XRough = SweepRange(i);
            break
        end
    end
    XT01 = Converger('T',Geometry,XRough,Resolution,AngularFrequency,SweepRange,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
end

% find instability limit at high frequency
AngularFrequency = pi*FrequencyRange(end)*2e3;
AngularFrequency2 = AngularFrequency^2;
kFi2 = AngularFrequency2/cFi2;
gFi = InnerFluid.Density*AngularFrequency2;
for m = 1:Layers
    kL2(m) = AngularFrequency2/Material{m}.LongitudinalVelocity^2;
    kT2(m) = AngularFrequency2/Material{m}.TransverseVelocity^2;
end
SweepRange = 0:EvanescenceLimit;
Wavenumber = AngularFrequency./SweepRange;
Y = Computer_L(Geometry,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
% figure,plot(SweepRange,20*log10(abs(Y)))

% SweepRange = 0:1e4;
% Wavenumber = AngularFrequency./SweepRange;
% YF = Computer_F(Geometry,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
% YL = Computer_L(Geometry,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
% YT = Computer_T(Geometry,Wavenumber,Mu,kT2,Layers,Radius);
% figure('name','F'),plot(SweepRange,20*log10(abs(YF)))
% figure('name','L'),plot(SweepRange,20*log10(abs(YL)))
% figure('name','T'),plot(SweepRange,20*log10(abs(YT)))

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
function X = Converger(ModeFamily,Geometry,XRough,Resolution,AngularFrequency,SweepRange,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi)
    X = [];
    Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Resolution));
    for j = 1:length(XRough)
        PhaseVelocity = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            PhaseVelocity = PhaseVelocity(1):(PhaseVelocity(end)-PhaseVelocity(1))/4:PhaseVelocity(end);
            Wavenumber = AngularFrequency./PhaseVelocity;
            if  strcmp(ModeFamily,'F')
                Y = Computer_F(Geometry,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
            elseif strcmp(ModeFamily,'L')
                Y = Computer_L(Geometry,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
            elseif strcmp(ModeFamily,'T')
                Y = Computer_T(Geometry,Wavenumber,Mu,kT2,Layers,Radius);
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
function Y = Computer_F(Geometry,k,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi)
    n = 1;
    Size = size(k);
    k = reshape(k,1,1,[]);
    k2 = k.^2;
    Length = length(k);
    Y = NaN(Size);
    for m = 1:Layers
        x2 = kL2(m)-k2;
        y2 = kT2(m)-k2;
        x = sqrt(-x2);
        y = sqrt(-y2);
        if  m == 1 && strcmp(Geometry,'Rod')
            R = Radius(2);
            R2 = R^2;
            xR = x*R;
            yR = y*R;
            ZWnx = besseli(n,xR);
            ZWny = besseli(n,yR);
            ZWn1x = x.*besseli(n+1,xR);
            ZWn1y = y.*besseli(n+1,yR);
            dZWnx = n*ZWnx/R+ZWn1x;
            dZWny = n*ZWny/R+ZWn1y;
            dg1ZWn1y = -(n+1)*ZWn1y/R-y2.*ZWny;
            d2ZWnx = (n*(n-1)/R2-x2).*ZWnx-ZWn1x/R;
            d2ZWny = (n*(n-1)/R2-y2).*ZWny-ZWn1y/R;
            L{m} = zeros(6,6,Length);
            L{m}(1,1,:) = 1i*k.*ZWnx; % u1
            L{m}(1,2,:) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
            L{m}(2,1,:) = 1i*n*ZWnx/R; % u2
            L{m}(2,2,:) = -1i*k.*ZWn1y;
            L{m}(2,3,:) = 1i*dZWny;
            L{m}(3,1,:) = dZWnx; % u3
            L{m}(3,2,:) = k.*ZWn1y;
            L{m}(3,3,:) = n*ZWny/R;
            L{m}(4,1,:) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2(m)*ZWnx; % sigma33
            L{m}(4,2,:) = 2*Mu(m)*k.*dg1ZWn1y;
            L{m}(4,3,:) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
            L{m}(5,1,:) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
            L{m}(5,2,:) = 1i*Mu(m)*k.*((n+1)*ZWn1y/R-dg1ZWn1y);
            L{m}(5,3,:) = 1i*Mu(m)*(2*d2ZWny+y2.*ZWny);
            L{m}(6,1,:) = 2i*Mu(m)*k.*dZWnx; % sigma13
            L{m}(6,2,:) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2).*ZWn1y+n*dg1ZWn1y/R);
            L{m}(6,3,:) = 1i*Mu(m)*n*k.*ZWny/R;
        else
            R = Radius(m:m+1);
            R2 = R.^2;
            xR = x.*R;
            yR = y.*R;
            ZWnx = [besseli(n,xR) besselk(n,xR)];
            ZWny = [besseli(n,yR) besselk(n,yR)];
            ZWn1x = x.*[besseli(n+1,xR) besselk(n+1,xR)];
            ZWn1y = y.*[besseli(n+1,yR) besselk(n+1,yR)];
            dZWnx = n*ZWnx./R+[ZWn1x(:,1,:) -ZWn1x(:,2,:)];
            dZWny = n*ZWny./R+[ZWn1y(:,1,:) -ZWn1y(:,2,:)];
            dg1ZWn1y = -(n+1)*ZWn1y./R+y2.*[-ZWny(:,1,:) ZWny(:,2,:)];
            d2ZWnx = (n*(n-1)./R2-x2).*ZWnx+[-ZWn1x(:,1,:) ZWn1x(:,2,:)]./R;
            d2ZWny = (n*(n-1)./R2-y2).*ZWny+[-ZWn1y(:,1,:) ZWn1y(:,2,:)]./R;
            L1(1,[1 4],:) = 1i*k.*ZWnx(2,:,:); % u1
            L1(1,[2 5],:) = 1i*((n+1)*ZWn1y(2,:,:)/R(2)+dg1ZWn1y(2,:,:));
            L1(2,[1 4],:) = 1i*n*ZWnx(2,:,:)/R(2); % u2
            L1(2,[2 5],:) = -1i*k.*ZWn1y(2,:,:);
            L1(2,[3 6],:) = 1i*dZWny(2,:,:);
            L1(3,[1 4],:) = dZWnx(2,:,:); % u3
            L1(3,[2 5],:) = k.*ZWn1y(2,:,:);
            L1(3,[3 6],:) = n*ZWny(2,:,:)/R(2);
            L1(4,[1 4],:) = 2*Mu(m)*d2ZWnx(2,:,:)-Lambda(m)*kL2(m)*ZWnx(2,:,:); % sigma33
            L1(4,[2 5],:) = 2*Mu(m)*k.*dg1ZWn1y(2,:,:);
            L1(4,[3 6],:) = 2*Mu(m)*n*(dZWny(2,:,:)-ZWny(2,:,:)/R(2))/R(2);
            L1(5,[1 4],:) = 2i*Mu(m)*n*(dZWnx(2,:,:)-ZWnx(2,:,:)/R(2))/R(2); % sigma23
            L1(5,[2 5],:) = 1i*Mu(m)*k.*((n+1)*ZWn1y(2,:,:)/R(2)-dg1ZWn1y(2,:,:));
            L1(5,[3 6],:) = 1i*Mu(m)*(2*d2ZWny(2,:,:)+y2.*ZWny(2,:,:));
            L1(6,[1 4],:) = 2i*Mu(m)*k.*dZWnx(2,:,:); % sigma13
            L1(6,[2 5],:) = 1i*Mu(m)*((n*(n+1)/R2(2)+k2-y2).*ZWn1y(2,:,:)+n*dg1ZWn1y(2,:,:)/R(2));
            L1(6,[3 6],:) = 1i*Mu(m)*n*k.*ZWny(2,:,:)/R(2);
            L2(1,[1 4],:) = 1i*k.*ZWnx(1,:,:); % u1
            L2(1,[2 5],:) = 1i*((n+1)*ZWn1y(1,:,:)/R(1)+dg1ZWn1y(1,:,:));
            L2(2,[1 4],:) = 1i*n*ZWnx(1,:,:)/R(1); % u2
            L2(2,[2 5],:) = -1i*k.*ZWn1y(1,:,:);
            L2(2,[3 6],:) = 1i*dZWny(1,:,:);
            L2(3,[1 4],:) = dZWnx(1,:,:); % u3
            L2(3,[2 5],:) = k.*ZWn1y(1,:,:);
            L2(3,[3 6],:) = n*ZWny(1,:,:)/R(1);
            L2(4,[1 4],:) = 2*Mu(m)*d2ZWnx(1,:,:)-Lambda(m)*kL2(m)*ZWnx(1,:,:); % sigma33
            L2(4,[2 5],:) = 2*Mu(m)*k.*dg1ZWn1y(1,:,:);
            L2(4,[3 6],:) = 2*Mu(m)*n*(dZWny(1,:,:)-ZWny(1,:,:)/R(1))/R(1);
            L2(5,[1 4],:) = 2i*Mu(m)*n*(dZWnx(1,:,:)-ZWnx(1,:,:)/R(1))/R(1); % sigma23
            L2(5,[2 5],:) = 1i*Mu(m)*k.*((n+1)*ZWn1y(1,:,:)/R(1)-dg1ZWn1y(1,:,:));
            L2(5,[3 6],:) = 1i*Mu(m)*(2*d2ZWny(1,:,:)+y2.*ZWny(1,:,:));
            L2(6,[1 4],:) = 2i*Mu(m)*k.*dZWnx(1,:,:); % sigma13
            L2(6,[2 5],:) = 1i*Mu(m)*((n*(n+1)/R2(1)+k2-y2).*ZWn1y(1,:,:)+n*dg1ZWn1y(1,:,:)/R(1));
            L2(6,[3 6],:) = 1i*Mu(m)*n*k.*ZWny(1,:,:)/R(1);
            L{m} = pagemrdivide(L1,L2);
        end
    end
    M = L{end};
    for m = Layers-1:-1:1
        M = pagemtimes(M,L{m});
    end
    if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
        zi = sqrt(kFi2-k2);
        zRi = zi*Radius(1);
        Znzi = besselj(n,zRi);
        WFi = n*Znzi./Radius(1)-zi.*besselj(n+1,zRi);
        DFi = gFi*Znzi;
        for j = 1:Length
            M21(1,1,j) = det(M(4:6,1:3,j));
            M22(1,1,j) = det(M(4:6,[1 2 4],j));
        end
        Y = M21-DFi./WFi.*M22;
    else
        for j = 1:Length
            Y(j) = det(M(4:6,1:3,j));
        end
    end
    Y = reshape(real(Y),Size);
end
function Y = Computer_L(Geometry,k,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi)
    Size = size(k);
    k = reshape(k,1,1,[]);
    k2 = k.^2;
    Length = length(k);
    Y = NaN(Size);
    for m = 1:Layers
        x2 = kL2(m)-k2;
        y2 = kT2(m)-k2;
        x = sqrt(-x2);
        y = sqrt(-y2);
        if  m == 1 && strcmp(Geometry,'Rod')
            R = Radius(2);
            xR = x*R;
            yR = y*R;
            Z0x = besseli(0,xR);
            Z0y = besseli(0,yR);
            Z1x = x.*besseli(1,xR);
            Z1y = y.*besseli(1,yR);
            d2Z0y = -y2.*Z0y-Z1y/R;
            d2Z0x = -x2.*Z0x-Z1x/R;
            L{m} = zeros(4,4,Length);
            L{m}(1,1,:) = 1i*k.*Z0x; % u1
            L{m}(1,2,:) = 1i*(Z1y/R+d2Z0y);
            L{m}(2,1,:) = Z1x; % u3
            L{m}(2,2,:) = k.*Z1y;
            L{m}(3,1,:) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2(m)*Z0x; % sigma33
            L{m}(3,2,:) = 2*Mu(m)*k.*d2Z0y;
            L{m}(4,1,:) = 2i*Mu(m)*k.*Z1x; % sigma13
            L{m}(4,2,:) = 1i*Mu(m)*(k2-y2).*Z1y;
        else
            R = Radius(m:m+1);
            xR = x.*R;
            yR = y.*R;
            ZW0x = [besseli(0,xR) besselk(0,xR)];
            ZW0y = [besseli(0,yR) besselk(0,yR)];
            ZW1x = x.*[besseli(1,xR) besselk(1,xR)];
            ZW1y = y.*[besseli(1,yR) besselk(1,yR)];
            dZW0x = [ZW1x(:,1,:) -ZW1x(:,2,:)];
            dg1ZW1y = y2.*[-ZW0y(:,1,:) ZW0y(:,2,:)]-ZW1y./R;
            d2ZW0x = -x2.*ZW0x-dZW0x./R;
            L1(1,[1 3],:) = 1i*k.*ZW0x(2,:,:); % u1
            L1(1,[2 4],:) = 1i*(ZW1y(2,:,:)/R(2)+dg1ZW1y(2,:,:));
            L1(2,[1 3],:) = dZW0x(2,:,:); % u3
            L1(2,[2 4],:) = k.*ZW1y(2,:,:);
            L1(3,[1 3],:) = 2*Mu(m)*d2ZW0x(2,:,:)-Lambda(m)*kL2(m)*ZW0x(2,:,:); % sigma33
            L1(3,[2 4],:) = 2*Mu(m)*k.*dg1ZW1y(2,:,:);
            L1(4,[1 3],:) = 2i*Mu(m)*k.*dZW0x(2,:,:); % sigma13
            L1(4,[2 4],:) = 1i*Mu(m)*(k2-y2).*ZW1y(2,:,:);
            L2(1,[1 3],:) = 1i*k.*ZW0x(1,:,:); % u1
            L2(1,[2 4],:) = 1i*(ZW1y(1,:,:)/R(1)+dg1ZW1y(1,:,:));
            L2(2,[1 3],:) = dZW0x(1,:,:); % u3
            L2(2,[2 4],:) = k.*ZW1y(1,:,:);
            L2(3,[1 3],:) = 2*Mu(m)*d2ZW0x(1,:,:)-Lambda(m)*kL2(m)*ZW0x(1,:,:); % sigma33
            L2(3,[2 4],:) = 2*Mu(m)*k.*dg1ZW1y(1,:,:);
            L2(4,[1 3],:) = 2i*Mu(m)*k.*dZW0x(1,:,:); % sigma13
            L2(4,[2 4],:) = 1i*Mu(m)*(k2-y2).*ZW1y(1,:,:);
            L{m} = pagemrdivide(L1,L2);
        end
    end
    M = L{end};
    for m = Layers-1:-1:1
        M = pagemtimes(M,L{m});
    end
    if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
        zi = sqrt(kFi2-k2);
        zRi = zi*Radius(1);
        Z1zi = zi.*besselj(1,zRi);
        WFi = -Z1zi;
        DFi = gFi*besselj(0,zRi);
        for j = 1:Length
            M21(1,1,j) = det(M(3:4,1:2,j));
            M22(1,1,j) = det(M(3:4,[1 3],j));
        end
        Y = M21-DFi./WFi.*M22;
    else
        for j = 1:Length
            Y(j) = det(M(3:4,1:2,j));
        end
    end
    if  strcmp(Geometry,'Rod')
        Y = reshape(imag(Y),Size);
    elseif strcmp(Geometry,'Pipe')
        Y = reshape(real(Y),Size);
    end
end
function Y = Computer_T(Geometry,k,Mu,kT2,Layers,Radius)
    Size = size(k);
    k = reshape(k,1,1,[]);
    k2 = k.^2;
    Length = length(k);
    for m = 1:Layers
        y2 = kT2(m)-k2;
        y = sqrt(-y2);
        if  m == 1 && strcmp(Geometry,'Rod')
            R = Radius(2);
            yR = y*R;
            Z0y = besseli(0,yR);
            Z1y = y.*besseli(1,yR);
            L{m} = [[1i*Z1y;-1i*Mu(m)*(y2.*Z0y+2*Z1y/R)] zeros(2,1,Length)]; % u2; sigma23
        else
            R = Radius(m:m+1);
            yR = y.*R;
            ZW0y = [besseli(0,yR) besselk(0,yR)];
            ZW1y = y.*[besseli(1,yR) besselk(1,yR)];
            dZW0y = [ZW1y(:,1,:) -ZW1y(:,2,:)];
            d2ZW0y = -y2.*ZW0y-dZW0y./R;
            L1(1,:,:) = 1i*dZW0y(2,:,:); % u2
            L1(2,:,:) = 1i*Mu(m)*(2*d2ZW0y(2,:,:)+y2.*ZW0y(2,:,:)); % sigma23
            L2(1,:,:) = 1i*dZW0y(1,:,:); % u2
            L2(2,:,:) = 1i*Mu(m)*(2*d2ZW0y(1,:,:)+y2.*ZW0y(1,:,:)); % sigma23
            L{m} = pagemrdivide(L1,L2);
        end
    end
    M = L{end};
    for m = Layers-1:-1:1
        M = pagemtimes(M,L{m});
    end
    if  strcmp(Geometry,'Rod')
        Y = reshape(imag(M(2,1,:)),Size);
    elseif strcmp(Geometry,'Pipe')
        Y = reshape(real(M(2,1,:)),Size);
    end
end