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
function [HSScholte,HAScholte,HBScholte] = FrequencySweeper_Anisotropic_Scholte(c,Delta,Material,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,SweepRange,I,I1,Pattern,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled)
Resolution = 1e-5; % (kHz)
PhaseVelocityOffset = 1-1e-4; % (m/s)

%#ok<*AGROW>
HSScholte=[];HAScholte=[];HBScholte=[];
if  length(SweepRange) < 2
    return
end
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Resolution));
else
    Bisections = 1;
end
if  Symmetric
    PhaseVelocity = PhaseVelocityOffset*UpperFluid.Velocity;
else
    if  ToggleUpperFluid && ToggleLowerFluid
        if  UpperFluid.Velocity > LowerFluid.Velocity
            PhaseVelocity = PhaseVelocityOffset*LowerFluid.Velocity;
        else
            PhaseVelocity = PhaseVelocityOffset*UpperFluid.Velocity;
        end
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        PhaseVelocity = PhaseVelocityOffset*UpperFluid.Velocity;
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        PhaseVelocity = PhaseVelocityOffset*LowerFluid.Velocity;
    end
end
AngularFrequency = 2*pi*SweepRange*1e3;
AngularFrequency2 = AngularFrequency.^2;
kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
gUpperFluid = 1i*UpperFluid.Density*AngularFrequency2;
gLowerFluid = 1i*LowerFluid.Density*AngularFrequency2;
for m = 1:SuperLayerSize
    c{m} = real(c{m});
    if  ~Decoupled
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
        A1=0;
    else
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
        a11=0;a12=0;a23=0;a33=0;a34=0;
    end
end
XRough = [];
Wavenumber = AngularFrequency/PhaseVelocity;
if  ~Decoupled
    Y = Computer_Coupled(AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,a11,a21,a31,a12,a22,a23,a32,a33,a34);
else
    Y = Computer_Decoupled(AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,A1,a21,a31,a22,a32);
end
for i = 2:length(SweepRange)-1
    if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
        XRough(end+1) = SweepRange(i);
    end
end
if  isempty(XRough)
    disp('No higher order Scholte modes found!')
    return
end
X = Converger(Decoupled,XRough,PhaseVelocity,SweepRange,Bisections,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,c,A1,a11,a21,a31,a12,a22,a23,a32,a33,a34);
if  Symmetric
    if  ~Decoupled
        [HSScholte,HAScholte] = SymmetryChecker_Coupled(X,PhaseVelocity,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,UpperFluid,c,a11,a21,a31,a12,a22,a23,a32,a33,a34,HSScholte,HAScholte);
    else
        [HSScholte,HAScholte] = SymmetryChecker_Decoupled(X,PhaseVelocity,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,UpperFluid,c,A1,a21,a31,a22,a32,HSScholte,HAScholte);
    end
else
    HBScholte = X;
end
if  Symmetric
    if  isempty(HAScholte)
        HSScholte = [];
    else
        HSScholte(HSScholte < HAScholte(1)) = [];
    end
    if  isempty(HSScholte) && isempty(HAScholte)
        disp('No higher order Scholte modes found!')
        return
    end
    String = ['Frq. @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode         Frq.(kHz)'];
    if  any(HAScholte)
        for i = 1:length(HAScholte)
            if  i < 10
                if  HAScholte(i) < 1e2
                    String = append(String,newline,'AScholte',num2str(i),'      ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e2 && HAScholte(i) < 1e3
                    String = append(String,newline,'AScholte',num2str(i),'     ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e3 && HAScholte(i) < 1e4
                    String = append(String,newline,'AScholte',num2str(i),'    ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e4
                    String = append(String,newline,'AScholte',num2str(i),'   ',num2str(HAScholte(i),'%.3f'));
                end
            else
                if  HAScholte(i) < 1e2
                    String = append(String,newline,'AScholte',num2str(i),'     ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e2 && HAScholte(i) < 1e3
                    String = append(String,newline,'AScholte',num2str(i),'    ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e3 && HAScholte(i) < 1e4
                    String = append(String,newline,'AScholte',num2str(i),'   ',num2str(HAScholte(i),'%.3f'));
                elseif HAScholte(i) >= 1e4
                    String = append(String,newline,'AScholte',num2str(i),'  ',num2str(HAScholte(i),'%.3f'));
                end
            end
        end
    end
    String = append(String,newline);
    if  any(HSScholte)
        for i = 1:length(HSScholte)
            if  i < 10
                if  HSScholte(i) < 1e2
                    String = append(String,newline,'SScholte',num2str(i),'      ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e2 && HSScholte(i) < 1e3
                    String = append(String,newline,'SScholte',num2str(i),'     ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e3 && HSScholte(i) < 1e4
                    String = append(String,newline,'SScholte',num2str(i),'    ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e4
                    String = append(String,newline,'SScholte',num2str(i),'   ',num2str(HSScholte(i),'%.3f'));
                end
            else
                if  HSScholte(i) < 1e2
                    String = append(String,newline,'SScholte',num2str(i),'     ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e2 && HSScholte(i) < 1e3
                    String = append(String,newline,'SScholte',num2str(i),'    ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e3 && HSScholte(i) < 1e4
                    String = append(String,newline,'SScholte',num2str(i),'   ',num2str(HSScholte(i),'%.3f'));
                elseif HSScholte(i) >= 1e4
                    String = append(String,newline,'SScholte',num2str(i),'  ',num2str(HSScholte(i),'%.3f'));
                end
            end
        end
    end
    String = append(String,newline,newline);
    if  any(HAScholte)
        String = append(String,'AScholte: ',num2str(length(HAScholte)));
    end
    String = append(String,newline);
    if  any(HSScholte)
        String = append(String,'SScholte: ',num2str(length(HSScholte)));
    end
else
    String = ['Frq. @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode         Frq.(kHz)'];
    if  any(HBScholte)
        for i = 1:length(HBScholte)
            if  i < 10
                if  HBScholte(i) < 1e2
                    String = append(String,newline,'BScholte',num2str(i),'      ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e2 && HBScholte(i) < 1e3
                    String = append(String,newline,'BScholte',num2str(i),'     ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e3 && HBScholte(i) < 1e4
                    String = append(String,newline,'BScholte',num2str(i),'    ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e4
                    String = append(String,newline,'BScholte',num2str(i),'   ',num2str(HBScholte(i),'%.3f'));
                end
            else
                if  HBScholte(i) < 1e2
                    String = append(String,newline,'BScholte',num2str(i),'     ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e2 && HBScholte(i) < 1e3
                    String = append(String,newline,'BScholte',num2str(i),'    ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e3 && HBScholte(i) < 1e4
                    String = append(String,newline,'BScholte',num2str(i),'   ',num2str(HBScholte(i),'%.3f'));
                elseif HBScholte(i) >= 1e4
                    String = append(String,newline,'BScholte',num2str(i),'  ',num2str(HBScholte(i),'%.3f'));
                end
            end
        end
    end
    String = append(String,newline,newline);
    if  any(HBScholte)
        String = append(String,'BScholte: ',num2str(length(HBScholte)));
    end
end
disp([String,newline,'----------------'])
end
function X = Converger(Decoupled,XRough,PhaseVelocity,SweepRange,Bisections,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,I1,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,c,A1,a11,a21,a31,a12,a22,a23,a32,a33,a34)
    X = [];
    for j = 1:length(XRough)
        Frequency = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            AngularFrequency = 2*pi*Frequency*1e3;
            AngularFrequency2 = AngularFrequency.^2;
            kUpperFluid2 = AngularFrequency2/UpperFluid.Velocity^2;
            kLowerFluid2 = AngularFrequency2/LowerFluid.Velocity^2;
            gUpperFluid = 1i*UpperFluid.Density*AngularFrequency2;
            gLowerFluid = 1i*LowerFluid.Density*AngularFrequency2;
            Wavenumber = AngularFrequency/PhaseVelocity;
            if  ~Decoupled
                Y = Computer_Coupled(AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,a11,a21,a31,a12,a22,a23,a32,a33,a34);
            else
                Y = Computer_Decoupled(AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,A1,a21,a31,a22,a32);
            end
            for i = 2:length(Frequency)-1
                if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1))
                    if  o == Bisections
                        X(end+1) = Frequency(i);
                    end
                    Frequency = [Frequency(i-1) Frequency(i+1)];
                    break
                end
            end
        end
    end
end
function Y = Computer_Coupled(AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,a11,a21,a31,a12,a22,a23,a32,a33,a34)
    AngularFrequency2 = reshape(AngularFrequency2,1,1,[]);
    kUpperFluid2 = reshape(kUpperFluid2,1,1,[]);
    kLowerFluid2 = reshape(kLowerFluid2,1,1,[]);
    gUpperFluid = reshape(gUpperFluid,1,1,[]);
    gLowerFluid = reshape(gLowerFluid,1,1,[]);
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Wavenumber6 = Wavenumber2.^3;
    Length = length(Wavenumber);
    for m = 1:SuperLayerSize
        rw2 = Material{m}.Density*AngularFrequency2;
        r2w4 = rw2.^2;
        A1 = a11(m)*Wavenumber2+a12(m)*rw2;
        A2 = a21(m)*Wavenumber4+a22(m)*rw2.*Wavenumber2+a23(m)*r2w4;
        A3 = a31(m)*Wavenumber6+a32(m)*rw2.*Wavenumber4+a33(m)*r2w4.*Wavenumber2+a34(m)*rw2.^3;
        d1 = A1/3;
        d2 = A2/3-d1.^2;
        d3 = d1.^3-d1.*A2/2+A3/2;
        d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
        d5 = d2./d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        k3(1,1,:) = sqrt(d6+d7);
        k3(1,2,:) = sqrt(d6-d7);
        k3(1,3,:) = sqrt(d4-d5-d1);
        k32 = k3.^2;
        k3k = k3.*Wavenumber;
        m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2;
        m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2;
        m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
        m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
        m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
        m1 = m13.*m22-m12.*m23;
        V = (m11.*m23-m13.*m12)./m1;
        W = (m11.*m22-m12.^2)./-m1;
        e1 = Wavenumber.*W+k3;
        e2 = k3.*V;
        D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V).*Wavenumber+c{m}(3,3)*k3.*W);
        D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
        D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
        E = exp(1i*k3*LayerThicknesses(m));
        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
        L2 = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W];
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M0 = L{m}(1:3,1:3,:)-M{1}(4:6,4:6,:);
        M1 = pagemrdivide(M{1}(1:3,4:6,:),M0);
        M2 = pagemrdivide(L{m}(4:6,1:3,:),M0);
        M{1} = [M{1}(1:3,1:3,:)+pagemtimes(M1,M{1}(4:6,1:3,:)) -pagemtimes(M1,L{m}(1:3,4:6,:));pagemtimes(M2,M{1}(4:6,1:3,:)) L{m}(4:6,4:6,:)-pagemtimes(M2,L{m}(1:3,4:6,:))];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1:3,1:3,:)-M{m}(4:6,4:6,:);
        M1 = pagemrdivide(M{m}(1:3,4:6,:),M0);
        M2 = pagemrdivide(M{Pattern(m)}(4:6,1:3,:),M0);
        M{m+1} = [M{m}(1:3,1:3,:)+pagemtimes(M1,M{m}(4:6,1:3,:)) -pagemtimes(M1,M{Pattern(m)}(1:3,4:6,:));pagemtimes(M2,M{m}(4:6,1:3,:)) M{Pattern(m)}(4:6,4:6,:)-pagemtimes(M2,M{Pattern(m)}(1:3,4:6,:))];
    end
    if  SymmetricSystem
        M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
        M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
        M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
        M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
    end
    M{end} = pageinv(M{end});
    if  ToggleUpperFluid && ToggleLowerFluid
        k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
        k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
        WUpperFluid = k3UpperFluid./Wavenumber;
        WLowerFluid = k3LowerFluid./Wavenumber;
        DUpperFluid = gUpperFluid./Wavenumber;
        DLowerFluid = gLowerFluid./Wavenumber;
        Y = reshape(real((WLowerFluid-M{end}(6,4,:).*DLowerFluid).*(WUpperFluid+M{end}(3,1,:).*DUpperFluid)+M{end}(3,4,:).*M{end}(6,1,:).*DUpperFluid.*DLowerFluid),Size); 
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
        Y = reshape(real(k3UpperFluid./gUpperFluid+M{end}(3,1,:)),Size);
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
        Y = reshape(real(k3LowerFluid./gLowerFluid-M{end}(6,4,:)),Size);
    end
end
function Y = Computer_Decoupled(AngularFrequency2,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,ToggleUpperFluid,ToggleLowerFluid,kUpperFluid2,kLowerFluid2,gUpperFluid,gLowerFluid,c,A1,a21,a31,a22,a32)
    AngularFrequency2 = reshape(AngularFrequency2,1,1,[]);
    kUpperFluid2 = reshape(kUpperFluid2,1,1,[]);
    kLowerFluid2 = reshape(kLowerFluid2,1,1,[]);
    gUpperFluid = reshape(gUpperFluid,1,1,[]);
    gLowerFluid = reshape(gLowerFluid,1,1,[]);
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Wavenumber2 = Wavenumber.^2;
    Wavenumber4 = Wavenumber2.^2;
    Length = length(Wavenumber);
    for m = 1:SuperLayerSize
        rw2 = Material{m}.Density*AngularFrequency2;
        A2 = a21(m)*Wavenumber2+a22(m)*rw2;
        A3 = a31(m)*Wavenumber4+a32(m)*rw2.*Wavenumber2+rw2.^2;
        d1 = sqrt(A2.^2-2*A1(m)*A3);
        k3(1,1,:) = sqrt((-A2+d1)/A1(m));
        k3(1,2,:) = sqrt((-A2-d1)/A1(m));
        W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k3);
        D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W);
        D5 = 1i*(c{m}(5,5)*(k3+Wavenumber.*W));
        E = exp(1i*k3*LayerThicknesses(m));
        L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
        L2 = [ones(1,2,Length) E;W -W.*E;E ones(1,2,Length);W.*E -W];
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M0 = L{m}(1:2,1:2,:)-M{1}(3:4,3:4,:);
        M1 = pagemrdivide(M{1}(1:2,3:4,:),M0);
        M2 = pagemrdivide(L{m}(3:4,1:2,:),M0);
        M{1} = [M{1}(1:2,1:2,:)+pagemtimes(M1,M{1}(3:4,1:2,:)) -pagemtimes(M1,L{m}(1:2,3:4,:));pagemtimes(M2,M{1}(3:4,1:2,:)) L{m}(3:4,3:4,:)-pagemtimes(M2,L{m}(1:2,3:4,:))];
    end
    for m = 1:length(Pattern)
        M0 = M{Pattern(m)}(1:2,1:2,:)-M{m}(3:4,3:4,:);
        M1 = pagemrdivide(M{m}(1:2,3:4,:),M0);
        M2 = pagemrdivide(M{Pattern(m)}(3:4,1:2,:),M0);
        M{m+1} = [M{m}(1:2,1:2,:)+pagemtimes(M1,M{m}(3:4,1:2,:)) -pagemtimes(M1,M{Pattern(m)}(1:2,3:4,:));pagemtimes(M2,M{m}(3:4,1:2,:)) M{Pattern(m)}(3:4,3:4,:)-pagemtimes(M2,M{Pattern(m)}(1:2,3:4,:))];
    end
    if  SymmetricSystem
        M0 = M{end}(3:4,3:4,:).*I1-M{end}(3:4,3:4,:);
        M1 = pagemrdivide(M{end}(1:2,3:4,:),M0);
        M2 = pagemrdivide(M{end}(1:2,3:4,:).*I1,M0);
        M{end} = [M{end}(1:2,1:2,:)+pagemtimes(M1,M{end}(3:4,1:2,:)) -pagemtimes(M1,M{end}(3:4,1:2,:).*I1);pagemtimes(M2,M{end}(3:4,1:2,:)) M{end}(1:2,1:2,:).*I1-pagemtimes(M2,M{end}(3:4,1:2,:).*I1)];
    end
    M{end} = pageinv(M{end});
    if  ToggleUpperFluid && ToggleLowerFluid
        k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
        k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
        WUpperFluid = k3UpperFluid./Wavenumber;
        WLowerFluid = k3LowerFluid./Wavenumber;
        DUpperFluid = gUpperFluid./Wavenumber;
        DLowerFluid = gLowerFluid./Wavenumber;
        Y = reshape(real((WLowerFluid-M{end}(4,3,:).*DLowerFluid).*(WUpperFluid+M{end}(2,1,:).*DUpperFluid)+M{end}(2,3,:).*M{end}(4,1,:).*DUpperFluid.*DLowerFluid),Size); 
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
        Y = reshape(real(k3UpperFluid./gUpperFluid+M{end}(2,1,:)),Size);
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
        Y = reshape(real(k3LowerFluid./gLowerFluid-M{end}(4,3,:)),Size);
    end
end
function [HSScholte,HAScholte] = SymmetryChecker_Coupled(X,PhaseVelocity,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,Fluid,c,a11,a21,a31,a12,a22,a23,a32,a33,a34,HSScholte,HAScholte)
    for p = 1:length(X)
        AngularFrequency = 2*pi*X(p)*1e3;
        AngularFrequency2 = AngularFrequency^2;
        kFluid2 = AngularFrequency2/Fluid.Velocity^2;
        gFluid = 1i*Fluid.Density*AngularFrequency2;
        Wavenumber = AngularFrequency/PhaseVelocity;
        Wavenumber2 = Wavenumber^2;
        Wavenumber4 = Wavenumber2^2;
        Wavenumber6 = Wavenumber2^3;
        for m = 1:SuperLayerSize
            rw2 = Material{m}.Density*AngularFrequency2;
            r2w4 = rw2^2;
            A1 = a11(m)*Wavenumber2+a12(m)*rw2;
            A2 = a21(m)*Wavenumber4+a22(m)*rw2*Wavenumber2+a23(m)*r2w4;
            A3 = a31(m)*Wavenumber6+a32(m)*rw2*Wavenumber4+a33(m)*r2w4*Wavenumber2+a34(m)*rw2^3;
            d1 = A1/3;
            d2 = A2/3-d1^2;
            d3 = d1^3-d1*A2/2+A3/2;
            d4 = (sqrt(d2^3+d3^2)-d3)^(1/3);
            d5 = d2/d4;
            d6 = (d5-d4)/2-d1;
            d7 = (d5+d4)/2i*sqrt(3);
            k3(1) = sqrt(d6+d7);
            k3(2) = sqrt(d6-d7);
            k3(3) = sqrt(d4-d5-d1);
            k32 = k3.^2;
            k3k = k3*Wavenumber;
            m11 = c{m}(1,1)*Wavenumber2+c{m}(5,5)*k32-rw2;
            m22 = c{m}(6,6)*Wavenumber2+c{m}(4,4)*k32-rw2;
            m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
            m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
            m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
            m1 = m13.*m22-m12.*m23;
            V = (m11.*m23-m13.*m12)./m1;
            W = (m11.*m22-m12.^2)./-m1;
            e1 = Wavenumber*W+k3;
            e2 = k3.*V;
            D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V)*Wavenumber+c{m}(3,3)*k3.*W);
            D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2);
            D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2);
            E = exp(1i*k3*LayerThicknesses(m));
            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
            L2 = [ones(1,3) E;V V.*E;W -W.*E;E ones(1,3);V.*E V;W.*E -W];
            L{m} = L1/L2;
        end
        M{1} = L{1};
        for m = 2:SuperLayerSize
            M0 = L{m}(1:3,1:3)-M{1}(4:6,4:6);
            M1 = M{1}(1:3,4:6)/M0;
            M2 = L{m}(4:6,1:3)/M0;
            M{1} = [M{1}(1:3,1:3)+M1*M{1}(4:6,1:3) -M1*L{m}(1:3,4:6);M2*M{1}(4:6,1:3) L{m}(4:6,4:6)-M2*L{m}(1:3,4:6)];
        end
        for m = 1:length(Pattern)
            M0 = M{Pattern(m)}(1:3,1:3)-M{m}(4:6,4:6);
            M1 = M{m}(1:3,4:6)/M0;
            M2 = M{Pattern(m)}(4:6,1:3)/M0;
            M{m+1} = [M{m}(1:3,1:3)+M1*M{m}(4:6,1:3) -M1*M{Pattern(m)}(1:3,4:6);M2*M{m}(4:6,1:3) M{Pattern(m)}(4:6,4:6)-M2*M{Pattern(m)}(1:3,4:6)];
        end
        if  SymmetricSystem
            M0 = M{end}(4:6,4:6).*I-M{end}(4:6,4:6);
            M1 = M{end}(1:3,4:6)/M0;
            M2 = M{end}(1:3,4:6).*I/M0;
            M{end} = [M{end}(1:3,1:3,:)+M1*M{end}(4:6,1:3) -M1*(M{end}(4:6,1:3).*I);M2*M{end}(4:6,1:3,:) M{end}(1:3,1:3,:).*I-M2*(M{end}(4:6,1:3,:).*I)];
        end
        M{end} = inv(M{end});
        WFluid = sqrt(kFluid2-Wavenumber2)/Wavenumber;
        DFluid = gFluid/Wavenumber;
        UFluid = -(WFluid+M{end}(3,1)*DFluid)/(M{end}(3,4)*DFluid);
        if  sign(real(UFluid)) == 1 % this works in the nonattenuated case only
            HSScholte(end+1) = X(p);
        else
            HAScholte(end+1) = X(p);
        end
    end
end
function [HSScholte,HAScholte] = SymmetryChecker_Decoupled(X,PhaseVelocity,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,Fluid,c,A1,a21,a31,a22,a32,HSScholte,HAScholte)
    for p = 1:length(X)
        AngularFrequency = 2*pi*X(p)*1e3;
        AngularFrequency2 = AngularFrequency^2;
        kFluid2 = AngularFrequency2/Fluid.Velocity^2;
        gFluid = 1i*Fluid.Density*AngularFrequency2;
        Wavenumber = AngularFrequency/PhaseVelocity;
        Wavenumber2 = Wavenumber^2;
        Wavenumber4 = Wavenumber2^2;
        for m = 1:SuperLayerSize
            rw2 = Material{m}.Density*AngularFrequency2;
            A2 = a21(m)*Wavenumber2+a22(m)*rw2;
            A3 = a31(m)*Wavenumber4+a32(m)*rw2*Wavenumber2+rw2^2;
            d1 = sqrt(A2^2-2*A1(m)*A3);
            k3(1) = sqrt((-A2+d1)/A1(m));
            k3(2) = sqrt((-A2-d1)/A1(m));
            W = (rw2-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber*k3);
            D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W);
            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W));
            E = exp(1i*k3*LayerThicknesses(m));
            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
            L2 = [ones(1,2) E;W -W.*E;E ones(1,2);W.*E -W];
            L{m} = L1/L2;
        end
        M{1} = L{1};
        for m = 2:SuperLayerSize
            M0 = L{m}(1:2,1:2)-M{1}(3:4,3:4);
            M1 = M{1}(1:2,3:4)/M0;
            M2 = L{m}(3:4,1:2)/M0;
            M{1} = [M{1}(1:2,1:2)+M1*M{1}(3:4,1:2) -M1*L{m}(1:2,3:4);M2*M{1}(3:4,1:2) L{m}(3:4,3:4)-M2*L{m}(1:2,3:4)];
        end
        for m = 1:length(Pattern)
            M0 = M{Pattern(m)}(1:2,1:2)-M{m}(3:4,3:4);
            M1 = M{m}(1:2,3:4)/M0;
            M2 = M{Pattern(m)}(3:4,1:2)/M0;
            M{m+1} = [M{m}(1:2,1:2)+M1*M{m}(3:4,1:2) -M1*M{Pattern(m)}(1:2,3:4);M2*M{m}(3:4,1:2) M{Pattern(m)}(3:4,3:4)-M2*M{Pattern(m)}(1:2,3:4)];
        end
        if  SymmetricSystem
            M0 = M{end}(3:4,3:4).*I1-M{end}(3:4,3:4);
            M1 = M{end}(1:2,3:4)/M0;
            M2 = M{end}(1:2,3:4).*I1/M0;
            M{end} = [M{end}(1:2,1:2,:)+M1*M{end}(3:4,1:2) -M1*(M{end}(3:4,1:2).*I1);M2*M{end}(3:4,1:2,:) M{end}(1:2,1:2,:).*I1-M2*(M{end}(3:4,1:2,:).*I1)];
        end
        M{end} = inv(M{end});
        WFluid = sqrt(kFluid2-Wavenumber2)/Wavenumber;
        DFluid = gFluid/Wavenumber;
        UFluid = -(WFluid+M{end}(2,1)*DFluid)/(M{end}(2,3)*DFluid);
        if  sign(real(UFluid)) == 1 % this works in the nonattenuated case only
            HSScholte(end+1) = X(p);
        else
            HAScholte(end+1) = X(p);
        end
    end
end