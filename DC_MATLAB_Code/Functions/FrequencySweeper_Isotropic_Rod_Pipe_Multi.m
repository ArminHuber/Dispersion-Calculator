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
function [HF,HL,HT] = FrequencySweeper_Isotropic_Rod_Pipe_Multi(Geometry,PhaseVelocity,SweepRange,FlexuralModeOrders,Material,Layers,Radius,InnerFluid,ToggleInnerFluid,EvanescenceLimit,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Bisections = 17; % ceil(log2(1e5))

%#ok<*AGROW>
HF{FlexuralModeOrders}=[];HL=[];HT=[];XRoughF{FlexuralModeOrders}=[];XRoughL=[];XRoughT=[];
if  PhaseVelocity > EvanescenceLimit
    MatrixMethod = 1; % TMM
else
    MatrixMethod = 2; % SMM
end
cFi2 = InnerFluid.Velocity^2;
AngularFrequency = pi*SweepRange*2e3;
AngularFrequency2 = AngularFrequency.^2;
kFi2 = AngularFrequency2/cFi2;
gFi = InnerFluid.Density*AngularFrequency2;
for m = 1:Layers
    Lambda(m) = Material{m}.Lambda;
    Mu(m) = Material{m}.Mu;
    kL2(m,1,:) = AngularFrequency2/Material{m}.LongitudinalVelocity^2;
    kT2(m,1,:) = AngularFrequency2/Material{m}.TransverseVelocity^2;
end
CriticalVelocity = Material{1}.LongitudinalVelocity; % for F-modes in a rod computed with TMM
Wavenumber = AngularFrequency/PhaseVelocity;
for n = 1:FlexuralModeOrders
    Y = Computer_F(n,Geometry,MatrixMethod,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi,PhaseVelocity,CriticalVelocity);
    for i = 2:length(SweepRange)-1
        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
            XRoughF{n}(end+1) = SweepRange(i);
        end
    end
    if  isempty(XRoughF{n})
        XRoughF(n:end) = [];
        break
    end
    % figure('name','F'),plot(SweepRange,20*log10(Y))
end
YL = Computer_L(Geometry,MatrixMethod,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
YT = Computer_T(Geometry,Wavenumber,Mu,kT2,Layers,Radius);
% figure('name','L'),plot(SweepRange,20*log10(YL))
% figure('name','T'),plot(SweepRange,20*log10(YT))
for i = 2:length(SweepRange)-1
    if  YL(i) < YL(i-1) && YL(i) < YL(i+1)
        XRoughL(end+1) = SweepRange(i);
    end
    if  YT(i) < YT(i-1) && YT(i) < YT(i+1)
        XRoughT(end+1) = SweepRange(i);
    end  
end
if  isempty(XRoughL) && isempty(XRoughF)
    disp('No higher order modes found!')
    return
end
for n = 1:length(XRoughF)
    HF{n} = Converger('F',n,Geometry,XRoughF{n},Bisections,PhaseVelocity,SweepRange,Material,MatrixMethod,Lambda,Mu,Layers,Radius,InnerFluid,ToggleInnerFluid,cFi2,CriticalVelocity);
end
HL = Converger('L',n,Geometry,XRoughL,Bisections,PhaseVelocity,SweepRange,Material,MatrixMethod,Lambda,Mu,Layers,Radius,InnerFluid,ToggleInnerFluid,cFi2,CriticalVelocity);
HT = Converger('T',n,Geometry,XRoughT,Bisections,PhaseVelocity,SweepRange,Material,MatrixMethod,Lambda,Mu,Layers,Radius,InnerFluid,ToggleInnerFluid,cFi2,CriticalVelocity);
ModeNameLength = 1+7;
String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode',pad('Frq.(kHz)',ModeNameLength+8,'left')];
if  any(HF{1})
    NumFmodes = 0;
    for n = 1:length(HF)
        for i = 1:length(HF{n})
            String = append(String,newline,pad(['F(',num2str(n),',',num2str(i+1),')'],ModeNameLength),pad(num2str(HF{n}(i),'%.3f'),11,'left'));
        end
        NumFmodes = NumFmodes+i;
    end
end
OutputWindow1aUI1.String = String;
disp(String)
String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode',pad('Frq.(kHz)',ModeNameLength+8,'left')];
if  any(HL)
    for i = 1:length(HL)
        String = append(String,newline,pad(['L(0,',num2str(i+1),')'],ModeNameLength),pad(num2str(HL(i),'%.3f'),11,'left'));
    end
end
String = append(String,newline,' ');
if  any(HT)
    for i = 1:length(HT)
        String = append(String,newline,pad(['T(0,',num2str(i+1),')'],ModeNameLength),pad(num2str(HT(i),'%.3f'),11,'left'));
    end
end
OutputWindow1bUI1.String = String;
disp(extractAfter(String,')'))
disp(' ')
String = '';
if  any(HF{1})
    String = append(String,'F: ',num2str(NumFmodes));
end
OutputWindow2aUI1.String = String;
disp(String)
String = '';
if  any(HL)
    String = append(String,'L: ',num2str(length(HL)));
end
if  any(HT)
    String = append(String,newline,'T: ',num2str(length(HT)));
end
OutputWindow2bUI1.String = String;
disp([String,newline,'----------------'])
end
function X = Converger(ModeFamily,n,Geometry,XRough,Bisections,PhaseVelocity,SweepRange,Material,MatrixMethod,Lambda,Mu,Layers,Radius,InnerFluid,ToggleInnerFluid,cFi2,CriticalVelocity)
    X = [];
    for j = 1:length(XRough)
        Frequency = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            AngularFrequency = pi*Frequency*2e3;
            AngularFrequency2 = AngularFrequency.^2;
            kFi2 = AngularFrequency2/cFi2;
            gFi = InnerFluid.Density*AngularFrequency2;
            for m = 1:Layers
                kL2(m,1,:) = AngularFrequency2/Material{m}.LongitudinalVelocity^2;
                kT2(m,1,:) = AngularFrequency2/Material{m}.TransverseVelocity^2;
            end
            Wavenumber = AngularFrequency/PhaseVelocity;
            if  strcmp(ModeFamily,'F')
                [Y,Yc] = Computer_F(n,Geometry,MatrixMethod,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi,PhaseVelocity,CriticalVelocity);
            elseif strcmp(ModeFamily,'L')
                [Y,Yc] = Computer_L(Geometry,MatrixMethod,Wavenumber,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi);
            elseif strcmp(ModeFamily,'T')
                [Y,Yc] = Computer_T(Geometry,Wavenumber,Mu,kT2,Layers,Radius);
            end
            for i = 2:length(Frequency)-1
                if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
                    if  o == Bisections && sign(Yc(i-1)) ~= sign(Yc(i+1))
                        X(end+1) = Frequency(i);
                    end
                    Frequency = [Frequency(i-1) Frequency(i+1)];
                    break
                end
            end
        end
    end
end
function [Y,Yc] = Computer_F(n,Geometry,MatrixMethod,k,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi,PhaseVelocity,CriticalVelocity)
    Size = size(k);
    k = reshape(k,1,1,[]);
    k2 = k.^2;
    Length = length(k);
    Yc = NaN(Size);
    for m = 1:Layers
        x2 = kL2(m,1,:)-k2;
        y2 = kT2(m,1,:)-k2;
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
            if  MatrixMethod == 1
                L{m}(1,1,:) = 1i*k.*ZWnx; % u1
                L{m}(1,2,:) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
                L{m}(2,1,:) = 1i*n*ZWnx/R; % u2
                L{m}(2,2,:) = -1i*k.*ZWn1y;
                L{m}(2,3,:) = 1i*dZWny;
                L{m}(3,1,:) = dZWnx; % u3
                L{m}(3,2,:) = k.*ZWn1y;
                L{m}(3,3,:) = n*ZWny/R;
                L{m}(4,1,:) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2(m,1,:).*ZWnx; % sigma33
                L{m}(4,2,:) = 2*Mu(m)*k.*dg1ZWn1y;
                L{m}(4,3,:) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
                L{m}(5,1,:) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
                L{m}(5,2,:) = 1i*Mu(m)*k.*((n+1)*ZWn1y/R-dg1ZWn1y);
                L{m}(5,3,:) = 1i*Mu(m)*(2*d2ZWny+y2.*ZWny);
                L{m}(6,1,:) = 2i*Mu(m)*k.*dZWnx; % sigma13
                L{m}(6,2,:) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2).*ZWn1y+n*dg1ZWn1y/R);
                L{m}(6,3,:) = 1i*Mu(m)*n*k.*ZWny/R;
            elseif MatrixMethod == 2
                L1(1,1,:) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2(m,1,:).*ZWnx; % sigma33
                L1(1,2,:) = 2*Mu(m)*k.*dg1ZWn1y;
                L1(1,3,:) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
                L1(2,1,:) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
                L1(2,2,:) = 1i*Mu(m)*k.*((n+1)*ZWn1y/R-dg1ZWn1y);
                L1(2,3,:) = 1i*Mu(m)*(2*d2ZWny+y2.*ZWny);
                L1(3,1,:) = 2i*Mu(m)*k.*dZWnx; % sigma13
                L1(3,2,:) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2).*ZWn1y+n*dg1ZWn1y/R);
                L1(3,3,:) = 1i*Mu(m)*n*k.*ZWny/R;
                L2(1,1,:) = 1i*k.*ZWnx; % u1
                L2(1,2,:) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
                L2(2,1,:) = 1i*n*ZWnx/R; % u2
                L2(2,2,:) = -1i*k.*ZWn1y;
                L2(2,3,:) = 1i*dZWny;
                L2(3,1,:) = dZWnx; % u3
                L2(3,2,:) = k.*ZWn1y;
                L2(3,3,:) = n*ZWny/R;
                L{m}(4:6,4:6,:) = pagemrdivide(L1,L2);
            end
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
            if  MatrixMethod == 1
                L1(1,[1 4],:) = 1i*k.*ZWnx(2,:,:); % u1
                L1(1,[2 5],:) = 1i*((n+1)*ZWn1y(2,:,:)/R(2)+dg1ZWn1y(2,:,:));
                L1(2,[1 4],:) = 1i*n*ZWnx(2,:,:)/R(2); % u2
                L1(2,[2 5],:) = -1i*k.*ZWn1y(2,:,:);
                L1(2,[3 6],:) = 1i*dZWny(2,:,:);
                L1(3,[1 4],:) = dZWnx(2,:,:); % u3
                L1(3,[2 5],:) = k.*ZWn1y(2,:,:);
                L1(3,[3 6],:) = n*ZWny(2,:,:)/R(2);
                L1(4,[1 4],:) = 2*Mu(m)*d2ZWnx(2,:,:)-Lambda(m)*kL2(m,1,:).*ZWnx(2,:,:); % sigma33
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
                L2(4,[1 4],:) = 2*Mu(m)*d2ZWnx(1,:,:)-Lambda(m)*kL2(m,1,:).*ZWnx(1,:,:); % sigma33
                L2(4,[2 5],:) = 2*Mu(m)*k.*dg1ZWn1y(1,:,:);
                L2(4,[3 6],:) = 2*Mu(m)*n*(dZWny(1,:,:)-ZWny(1,:,:)/R(1))/R(1);
                L2(5,[1 4],:) = 2i*Mu(m)*n*(dZWnx(1,:,:)-ZWnx(1,:,:)/R(1))/R(1); % sigma23
                L2(5,[2 5],:) = 1i*Mu(m)*k.*((n+1)*ZWn1y(1,:,:)/R(1)-dg1ZWn1y(1,:,:));
                L2(5,[3 6],:) = 1i*Mu(m)*(2*d2ZWny(1,:,:)+y2.*ZWny(1,:,:));
                L2(6,[1 4],:) = 2i*Mu(m)*k.*dZWnx(1,:,:); % sigma13
                L2(6,[2 5],:) = 1i*Mu(m)*((n*(n+1)/R2(1)+k2-y2).*ZWn1y(1,:,:)+n*dg1ZWn1y(1,:,:)/R(1));
                L2(6,[3 6],:) = 1i*Mu(m)*n*k.*ZWny(1,:,:)/R(1);
            elseif MatrixMethod == 2
                L1([1 4],[1 4],:) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2(m,1,:).*ZWnx; % sigma33
                L1([1 4],[2 5],:) = 2*Mu(m)*k.*dg1ZWn1y;
                L1([1 4],[3 6],:) = 2*Mu(m)*n*(dZWny-ZWny./R)./R;
                L1([2 5],[1 4],:) = 2i*Mu(m)*n*(dZWnx-ZWnx./R)./R; % sigma23
                L1([2 5],[2 5],:) = 1i*Mu(m)*k.*((n+1)*ZWn1y./R-dg1ZWn1y);
                L1([2 5],[3 6],:) = 1i*Mu(m)*(2*d2ZWny+y2.*ZWny);
                L1([3 6],[1 4],:) = 2i*Mu(m)*k.*dZWnx; % sigma13
                L1([3 6],[2 5],:) = 1i*Mu(m)*((n*(n+1)./R2+k2-y2).*ZWn1y+n*dg1ZWn1y./R);
                L1([3 6],[3 6],:) = 1i*Mu(m)*n*k.*ZWny./R;
                L2([1 4],[1 4],:) = 1i*k.*ZWnx; % u1
                L2([1 4],[2 5],:) = 1i*((n+1)*ZWn1y./R+dg1ZWn1y);
                L2([2 5],[1 4],:) = 1i*n*ZWnx./R; % u2
                L2([2 5],[2 5],:) = -1i*k.*ZWn1y;
                L2([2 5],[3 6],:) = 1i*dZWny;
                L2([3 6],[1 4],:) = dZWnx; % u3
                L2([3 6],[2 5],:) = k.*ZWn1y;
                L2([3 6],[3 6],:) = n*ZWny./R;
            end
            L{m} = pagemrdivide(L1,L2);
        end
    end
    if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
        zi = sqrt(reshape(kFi2,1,1,[])-k2);
        zRi = zi*Radius(1);
        Znzi = besselj(n,zRi);
        WFi = n*Znzi./Radius(1)-zi.*besselj(n+1,zRi);
        DFi = reshape(gFi,1,1,[]).*Znzi;
    end
    if  MatrixMethod == 1
        M = L{end};
        for m = Layers-1:-1:1
            M = pagemtimes(M,L{m});
        end
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            for j = 1:Length
                M21(1,1,j) = det(M(4:6,1:3,j));
                M22(1,1,j) = det(M(4:6,[1 2 4],j));
            end
            Yc = real(M21-DFi./WFi.*M22);
        else
            for j = 1:Length
                Yc(j) = det(M(4:6,1:3,j));
            end
            if  strcmp(Geometry,'Rod')
                if  PhaseVelocity > CriticalVelocity
                    if  mod(n,2)
                        Yc = imag(Yc);
                    else
                        Yc = real(Yc);
                    end
                else
                    Yc = real(Yc);
                end
            elseif strcmp(Geometry,'Pipe')
                Yc = real(Yc);
            end
        end
    elseif MatrixMethod == 2
        M = L{1};
        for m = 2:Layers
            M0 = L{m}(1:3,1:3,:)-M(4:6,4:6,:);
            M1 = pagemrdivide(M(1:3,4:6,:),M0);
            M2 = pagemrdivide(L{m}(4:6,1:3,:),M0);
            M = [M(1:3,1:3,:)+pagemtimes(M1,M(4:6,1:3,:)) -pagemtimes(M1,L{m}(1:3,4:6,:));pagemtimes(M2,M(4:6,1:3,:)) L{m}(4:6,4:6,:)-pagemtimes(M2,L{m}(1:3,4:6,:))];
        end
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            M = pageinv(M);
            Yc = real(WFi./DFi+M(3,1,:));
        else
            if  strcmp(Geometry,'Rod')
                for j = 1:Length
                    Yc(j) = real(det(M(4:6,4:6,j)));
                end
            elseif strcmp(Geometry,'Pipe')
                for j = 1:Length
                    Yc(j) = real(det(M(:,:,j)));
                end
            end
        end
    end
    Yc = reshape(Yc,Size);
    Y = abs(Yc);
end
function [Y,Yc] = Computer_L(Geometry,MatrixMethod,k,Lambda,Mu,kL2,kT2,Layers,Radius,ToggleInnerFluid,kFi2,gFi)
    Size = size(k);
    k = reshape(k,1,1,[]);
    k2 = k.^2;
    Length = length(k);
    Yc = NaN(Size);
    for m = 1:Layers
        x2 = kL2(m,1,:)-k2;
        y2 = kT2(m,1,:)-k2;
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
            if  MatrixMethod == 1
                L{m}(1,1,:) = 1i*k.*Z0x; % u1
                L{m}(1,2,:) = 1i*(Z1y/R+d2Z0y);
                L{m}(2,1,:) = Z1x; % u3
                L{m}(2,2,:) = k.*Z1y;
                L{m}(3,1,:) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2(m,1,:).*Z0x; % sigma33
                L{m}(3,2,:) = 2*Mu(m)*k.*d2Z0y;
                L{m}(4,1,:) = 2i*Mu(m)*k.*Z1x; % sigma13
                L{m}(4,2,:) = 1i*Mu(m)*(k2-y2).*Z1y;
            elseif MatrixMethod == 2
                L1(1,1,:) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2(m,1,:).*Z0x; % sigma33
                L1(1,2,:) = 2*Mu(m)*k.*d2Z0y;
                L1(2,1,:) = 2i*Mu(m)*k.*Z1x; % sigma13
                L1(2,2,:) = 1i*Mu(m)*(k2-y2).*Z1y;
                L2(1,1,:) = 1i*k.*Z0x; % u1
                L2(1,2,:) = 1i*(Z1y/R+d2Z0y);
                L2(2,1,:) = Z1x; % u3
                L2(2,2,:) = k.*Z1y;
                L{m}(3:4,3:4,:) = pagemrdivide(L1,L2);
            end
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
            if  MatrixMethod == 1
                L1(1,[1 3],:) = 1i*k.*ZW0x(2,:,:); % u1
                L1(1,[2 4],:) = 1i*(ZW1y(2,:,:)/R(2)+dg1ZW1y(2,:,:));
                L1(2,[1 3],:) = dZW0x(2,:,:); % u3
                L1(2,[2 4],:) = k.*ZW1y(2,:,:);
                L1(3,[1 3],:) = 2*Mu(m)*d2ZW0x(2,:,:)-Lambda(m)*kL2(m,1,:).*ZW0x(2,:,:); % sigma33
                L1(3,[2 4],:) = 2*Mu(m)*k.*dg1ZW1y(2,:,:);
                L1(4,[1 3],:) = 2i*Mu(m)*k.*dZW0x(2,:,:); % sigma13
                L1(4,[2 4],:) = 1i*Mu(m)*(k2-y2).*ZW1y(2,:,:);
                L2(1,[1 3],:) = 1i*k.*ZW0x(1,:,:); % u1
                L2(1,[2 4],:) = 1i*(ZW1y(1,:,:)/R(1)+dg1ZW1y(1,:,:));
                L2(2,[1 3],:) = dZW0x(1,:,:); % u3
                L2(2,[2 4],:) = k.*ZW1y(1,:,:);
                L2(3,[1 3],:) = 2*Mu(m)*d2ZW0x(1,:,:)-Lambda(m)*kL2(m,1,:).*ZW0x(1,:,:); % sigma33
                L2(3,[2 4],:) = 2*Mu(m)*k.*dg1ZW1y(1,:,:);
                L2(4,[1 3],:) = 2i*Mu(m)*k.*dZW0x(1,:,:); % sigma13
                L2(4,[2 4],:) = 1i*Mu(m)*(k2-y2).*ZW1y(1,:,:);
            elseif MatrixMethod == 2
                L1([1 3],[1 3],:) = 2*Mu(m)*d2ZW0x-Lambda(m)*kL2(m,1,:).*ZW0x; % sigma33
                L1([1 3],[2 4],:) = 2*Mu(m)*k.*dg1ZW1y;
                L1([2 4],[1 3],:) = 2i*Mu(m)*k.*dZW0x; % sigma13
                L1([2 4],[2 4],:) = 1i*Mu(m)*(k2-y2).*ZW1y;
                L2([1 3],[1 3],:) = 1i*k.*ZW0x; % u1
                L2([1 3],[2 4],:) = 1i*(ZW1y./R+dg1ZW1y);
                L2([2 4],[1 3],:) = dZW0x; % u3
                L2([2 4],[2 4],:) = k.*ZW1y;
            end
            L{m} = pagemrdivide(L1,L2);
        end
    end
    if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
        zi = sqrt(reshape(kFi2,1,1,[])-k2);
        zRi = zi*Radius(1);
        WFi = -zi.*besselj(1,zRi);
        DFi = reshape(gFi,1,1,[]).*besselj(0,zRi);
    end
    if  MatrixMethod == 1
        M = L{end};
        for m = Layers-1:-1:1
            M = pagemtimes(M,L{m});
        end
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            for j = 1:Length
                M21(1,1,j) = det(M(3:4,1:2,j));
                M22(1,1,j) = det(M(3:4,[1 3],j));
            end
            Yc = real(M21-DFi./WFi.*M22);
        else
            for j = 1:Length
                Yc(j) = det(M(3:4,1:2,j));
            end
            if  strcmp(Geometry,'Rod')
                Yc = imag(Yc);
            elseif strcmp(Geometry,'Pipe')
                Yc = real(Yc);
            end
        end
    elseif MatrixMethod == 2
        M = L{1};
        for m = 2:Layers
            M0 = L{m}(1:2,1:2,:)-M(3:4,3:4,:);
            M1 = pagemrdivide(M(1:2,3:4,:),M0);
            M2 = pagemrdivide(L{m}(3:4,1:2,:),M0);
            M = [M(1:2,1:2,:)+pagemtimes(M1,M(3:4,1:2,:)) -pagemtimes(M1,L{m}(1:2,3:4,:));pagemtimes(M2,M(3:4,1:2,:)) L{m}(3:4,3:4,:)-pagemtimes(M2,L{m}(1:2,3:4,:))];
        end
        if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
            M = pageinv(M);
            Yc = real(WFi./DFi+M(2,1,:));
        else
            if  strcmp(Geometry,'Rod')
                for j = 1:Length
                    Yc(j) = real(det(M(3:4,3:4,j)));
                end
            elseif strcmp(Geometry,'Pipe')
                for j = 1:Length
                    Yc(j) = real(det(M(:,:,j)));
                end
            end
        end
    end
    Yc = reshape(Yc,Size);
    Y = abs(Yc);
end
function [Y,Yc] = Computer_T(Geometry,k,Mu,kT2,Layers,Radius)
    Size = size(k);
    k = reshape(k,1,1,[]);
    k2 = k.^2;
    Length = length(k);
    for m = 1:Layers
        y2 = kT2(m,1,:)-k2;
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
        Yc = reshape(imag(M(2,1,:)),Size);
    elseif strcmp(Geometry,'Pipe')
        Yc = reshape(real(M(2,1,:)),Size);
    end
    Y = abs(Yc);
end