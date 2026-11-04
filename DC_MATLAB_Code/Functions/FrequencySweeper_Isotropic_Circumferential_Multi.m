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
function [HCLamb,HCShear] = FrequencySweeper_Isotropic_Circumferential_Multi(PhaseVelocity,SweepRange,Material,Layers,Radius,EvanescenceLimit,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Bisections = 17; % ceil(log2(1e5))

%#ok<*AGROW>
HCLamb=[];HCShear=[];XRough=[];XSHRough=[];
if  PhaseVelocity > EvanescenceLimit
    MatrixMethod = 1; % TMM
else
    MatrixMethod = 2; % SMM
end
AngularFrequency = pi*SweepRange*2e3;
for m = 1:Layers
    Lambda(m) = Material{m}.Lambda;
    Mu(m) = Material{m}.Mu;
    kL(m,1,:) = AngularFrequency/Material{m}.LongitudinalVelocity;
    kT(m,1,:) = AngularFrequency/Material{m}.TransverseVelocity;
end
kL2 = kL.^2;
kT2 = kT.^2;
Wavenumber = AngularFrequency/PhaseVelocity;
Y = Computer_C(MatrixMethod,Wavenumber,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius);
YSH = Computer_CSH(Wavenumber,Mu,kT,Layers,Radius);
% figure('name','Lamb'),plot(SweepRange,20*log10(Y))
% figure('name','SH'),plot(SweepRange,20*log10(YSH))
for i = 2:length(SweepRange)-1
    if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
        XRough(end+1) = SweepRange(i);
    end
    if  YSH(i) < YSH(i-1) && YSH(i) < YSH(i+1)
        XSHRough(end+1) = SweepRange(i);
    end
end
if  isempty(XRough) && isempty(XSHRough)
    String = 'No higher order modes found!';
    OutputWindow1aUI1.String = String;
    OutputWindow1bUI1.String = '';
    OutputWindow2aUI1.String = '';
    OutputWindow2bUI1.String = '';
    disp(String)
    return
end
HCLamb = Converger('Lamb',XRough,Bisections,PhaseVelocity,SweepRange,Material,MatrixMethod,Lambda,Mu,Layers,Radius);
HCShear = Converger('SH',XSHRough,Bisections,PhaseVelocity,SweepRange,Material,MatrixMethod,Lambda,Mu,Layers,Radius);
ModeNameLength = 1+5;
String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
if  any(HCLamb)
    for i = 1:length(HCLamb)
        String = append(String,newline,pad(['C',num2str(i)],ModeNameLength),pad(num2str(HCLamb(i),'%.3f'),11,'left'));
    end
end
OutputWindow1aUI1.String = String;
disp(String)
String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
if  any(HCShear)
    for i = 1:length(HCShear)
        String = append(String,newline,pad(['CSH',num2str(i)],ModeNameLength),pad(num2str(HCShear(i),'%.3f'),11,'left'));
    end
end
OutputWindow1bUI1.String = String;
disp(extractAfter(String,')'))
disp(' ')
String = '';
if  any(HCLamb)
    String = append(String,'C:   ',num2str(length(HCLamb)));
end
OutputWindow2aUI1.String = String;
disp(String)
String = '';
if  any(HCShear)
    String = append(String,'CSH: ',num2str(length(HCShear)));
end
OutputWindow2bUI1.String = String;
disp([String,newline,'----------------'])
end
function X = Converger(ModeFamily,XRough,Bisections,PhaseVelocity,SweepRange,Material,MatrixMethod,Lambda,Mu,Layers,Radius)
    X = [];
    for j = 1:length(XRough)
        Frequency = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            AngularFrequency = pi*Frequency*2e3;
            for m = 1:Layers
                kL(m,1,:) = AngularFrequency/Material{m}.LongitudinalVelocity;
                kT(m,1,:) = AngularFrequency/Material{m}.TransverseVelocity;
            end
            kL2 = kL.^2;
            kT2 = kT.^2;
            Wavenumber = AngularFrequency/PhaseVelocity;
            if  strcmp(ModeFamily,'Lamb')
                [Y,Yc] = Computer_C(MatrixMethod,Wavenumber,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius);
            elseif strcmp(ModeFamily,'SH')
                [Y,Yc] = Computer_CSH(Wavenumber,Mu,kT,Layers,Radius);
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
function [Y,Yc] = Computer_C(MatrixMethod,k,Lambda,Mu,kL,kT,kL2,kT2,Layers,Radius)
    Size = size(k);
    k = reshape(k,1,1,[]);
    kRo = k*Radius(end);
    kRoRep = repmat(kRo,2,1,1);
    kRo2 = kRo.^2;
    kRo1 = 1/(kRo+1);
    kRo_1 = 1/(kRo-1);
    Length = length(k);
    Yc = NaN(Size);
    for m = 1:Layers
        R = Radius(m:m+1);
        R2 = R.^2;
        kLr = kL(m,1,:).*R;
        kTr = kT(m,1,:).*R;
        JYkLr = [besselj(kRoRep,kLr) bessely(kRoRep,kLr)];
        JYkTr = [besselj(kRoRep,kTr) bessely(kRoRep,kTr)];
        JY_2kLr = [besselj(kRoRep-2,kLr) bessely(kRoRep-2,kLr)];
        JY_2kTr = [besselj(kRoRep-2,kTr) bessely(kRoRep-2,kTr)];
        JY2kLr = [besselj(kRoRep+2,kLr) bessely(kRoRep+2,kLr)];
        JY2kTr = [besselj(kRoRep+2,kTr) bessely(kRoRep+2,kTr)];
        dJYkLr = kL2(m,1,:).*R/4.*[kRo_1.*(JY_2kLr(:,1,:)+JYkLr(:,1,:))-kRo1.*(JYkLr(:,1,:)+JY2kLr(:,1,:)) kRo_1.*(JY_2kLr(:,2,:)+JYkLr(:,2,:))-kRo1.*(JYkLr(:,2,:)+JY2kLr(:,2,:))];
        dJYkTr = kT2(m,1,:).*R/4.*[kRo_1.*(JY_2kTr(:,1,:)+JYkTr(:,1,:))-kRo1.*(JYkTr(:,1,:)+JY2kTr(:,1,:)) kRo_1.*(JY_2kTr(:,2,:)+JYkTr(:,2,:))-kRo1.*(JYkTr(:,2,:)+JY2kTr(:,2,:))];
        d2JYkLr = kL2(m,1,:)/4.*[JY_2kLr(:,1,:)-2*JYkLr(:,1,:)+JY2kLr(:,1,:) JY_2kLr(:,2,:)-2*JYkLr(:,2,:)+JY2kLr(:,2,:)];
        d2JYkTr = kT2(m,1,:)/4.*[JY_2kTr(:,1,:)-2*JYkTr(:,1,:)+JY2kTr(:,1,:) JY_2kTr(:,2,:)-2*JYkTr(:,2,:)+JY2kTr(:,2,:)];
        if  MatrixMethod == 1
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
        elseif MatrixMethod == 2
            L1([1 3],[1 3],:) = 2*Mu(m)*d2JYkLr+Lambda(m)*(d2JYkLr+dJYkLr./R-kRo2.*JYkLr./R2); % sigma33
            L1([1 3],[2 4],:) = 2i*Mu(m)*kRo.*(dJYkTr-JYkTr./R)./R;
            L1([2 4],[1 3],:) = 2i*Mu(m)*kRo.*(dJYkLr-JYkLr./R)./R; % sigma23
            L1([2 4],[2 4],:) = Mu(m)*(dJYkTr./R-kRo2.*JYkTr./R2-d2JYkTr);
            L2([1 3],[1 3],:) = 1i*kRo.*JYkLr./R; % u2
            L2([1 3],[2 4],:) = -dJYkTr;
            L2([2 4],[1 3],:) = dJYkLr; % u3
            L2([2 4],[2 4],:) = 1i*kRo.*JYkTr./R;
        end
        L{m} = pagemrdivide(L1,L2);
    end
    if  MatrixMethod == 1
        M = L{end};
        for m = Layers-1:-1:1
            M = pagemtimes(M,L{m});
        end
        for j = 1:Length
            Yc(j) = det(M(3:4,1:2,j));
        end
    elseif MatrixMethod == 2
        M = L{1};
        for m = 2:Layers
            M0 = L{m}(1:2,1:2,:)-M(3:4,3:4,:);
            M1 = pagemrdivide(M(1:2,3:4,:),M0);
            M2 = pagemrdivide(L{m}(3:4,1:2,:),M0);
            M = [M(1:2,1:2,:)+pagemtimes(M1,M(3:4,1:2,:)) -pagemtimes(M1,L{m}(1:2,3:4,:));pagemtimes(M2,M(3:4,1:2,:)) L{m}(3:4,3:4,:)-pagemtimes(M2,L{m}(1:2,3:4,:))];
        end
        for j = 1:Length
            Yc(j) = det(M(:,:,j));
        end
    end
    Yc = reshape(real(Yc),Size);
    Y = abs(Yc);
end
function [Y,Yc] = Computer_CSH(k,Mu,kT,Layers,Radius)
    Size = size(k);
    k = reshape(k,1,1,[]);
    kRo = k*Radius(end);
    kRoRep = repmat(kRo,2,1,1);
    for m = 1:Layers
        kTr = kT(m,1,:).*Radius(m:m+1);
        JY_1 = [besselj(kRoRep-1,kTr) bessely(kRoRep-1,kTr)];
        JY1 = [besselj(kRoRep+1,kTr) bessely(kRoRep+1,kTr)];
        JY = kTr./kRo.*[JY_1(:,1,:)+JY1(:,1,:) JY_1(:,2,:)+JY1(:,2,:)]; % multiplied by two
        L1(1,:,:) = JY(2,:,:); % u1
        L1(2,:,:) = Mu(m)*kT(m,1,:).*(JY_1(2,:,:)-JY1(2,:,:)); % sigma13
        L2(1,:,:) = JY(1,:,:); % u1
        L2(2,:,:) = Mu(m)*kT(m,1,:).*(JY_1(1,:,:)-JY1(1,:,:)); % sigma13
        L{m} = pagemrdivide(L1,L2);
    end
    M = L{end};
    for m = Layers-1:-1:1
        M = pagemtimes(M,L{m});
    end
    Yc = reshape(real(M(2,1,:)),Size);
    Y = abs(Yc);
end