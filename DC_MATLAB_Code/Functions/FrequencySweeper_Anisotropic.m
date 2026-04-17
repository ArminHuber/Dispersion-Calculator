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
function [HSLamb,HSShear,HALamb,HAShear,HBLamb,HBShear] = FrequencySweeper_Anisotropic(c,Delta,Material,SweepRange,PhaseVelocity,I,I1,Pattern,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled,EvanescenceLimit,OutputWindow1aUI2,OutputWindow1bUI2,OutputWindow2aUI2,OutputWindow2bUI2)
Bisections = 17; % ceil(log2(1e5))

%#ok<*AGROW>
HSLamb=[];HSShear=[];HALamb=[];HAShear=[];HBLamb=[];HBShear=[];XRough=[];XSHRough=[];
if  length(SweepRange) < 2
    return
end
if  PhaseVelocity > EvanescenceLimit
    MatrixMethod = 1; % TMM
else
    MatrixMethod = 2; % SMM
end
PhaseVelocity2 = PhaseVelocity^2;
for m = 1:SuperLayerSize
    c{m} = real(c{m});
    rc2 = Material{m}.Density*PhaseVelocity2;
    if  ~Decoupled
        a11 = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12 = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21 = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22 = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23 = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31 = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32 = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33 = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34 = -1/Delta(m);
        r2c4 = rc2^2;
        A1 = a11+a12*rc2;
        A2 = a21+a22*rc2+a23*r2c4;
        A3 = a31+a32*rc2+a33*r2c4+a34*rc2^3;
        d1 = A1/3;
        d2 = A2/3-d1^2;
        d3 = d1^3-d1*A2/2+A3/2;
        d4 = (sqrt(d2^3+d3^2)-d3)^(1/3);
        d5 = d2/d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        Alpha2 = [d6+d7 d6-d7 d4-d5-d1];
        Alpha(m,:) = sqrt(Alpha2);
        m11 = c{m}(1,1)+c{m}(5,5)*Alpha2-rc2;
        m22 = c{m}(6,6)+c{m}(4,4)*Alpha2-rc2;
        m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
        m13 = (c{m}(1,3)+c{m}(5,5))*Alpha(m,:);
        m23 = (c{m}(3,6)+c{m}(4,5))*Alpha(m,:);
        m1 = m13.*m22-m12.*m23;
        V(m,:) = (m11.*m23-m13.*m12)./m1;
        W(m,:) = (m11.*m22-m12.^2)./-m1;
        e1 = Alpha(m,:)+W(m,:);
        e2 = Alpha(m,:).*V(m,:);
        D3(m,:) = c{m}(1,3)+c{m}(3,6)*V(m,:)+c{m}(3,3)*Alpha(m,:).*W(m,:);
        D4(m,:) = c{m}(4,5)*e1+c{m}(4,4)*e2;
        D5(m,:) = c{m}(5,5)*e1+c{m}(4,5)*e2;
        AlphaSH=0;DSH=0;
    else
        a21 = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22 = -c{m}(3,3)-c{m}(5,5);
        a31 = c{m}(1,1)*c{m}(5,5);
        a32 = -c{m}(1,1)-c{m}(5,5);
        A1 = 2*c{m}(3,3)*c{m}(5,5);
        A2 = a21+a22*rc2;
        A3 = a31+a32*rc2+rc2^2;
        d1 = sqrt(A2^2-2*A1*A3);
        Alpha2 = [d1-A2 -d1-A2]/A1;
        Alpha(m,:) = sqrt(Alpha2);
        W(m,:) = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha2)./((c{m}(1,3)+c{m}(5,5))*Alpha(m,:));
        D3(m,:) = c{m}(1,3)+c{m}(3,3)*Alpha(m,:).*W(m,:);
        D5(m,:) = c{m}(5,5)*(Alpha(m,:)+W(m,:));
        AlphaSH(m,1) = sqrt((rc2-c{m}(6,6))/c{m}(4,4));
        DSH(m,1) = AlphaSH(m)*c{m}(4,4);
        V=0;D4=0;
    end
end
Wavenumber = 2*pi*SweepRange*1e3/PhaseVelocity;
if  ~Decoupled
    Y = Computer_Coupled(MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,Alpha,V,W,D3,D4,D5);
    for i = 2:length(SweepRange)-1
        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
            XRough(end+1) = SweepRange(i);
        end
    end
else
    Y = Computer_Decoupled(MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,Alpha,W,D3,D5);
    YSH = Computer_Decoupled_SH(Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,AlphaSH,DSH);
    for i = 2:length(SweepRange)-1
        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
            XRough(end+1) = SweepRange(i);
        end
        if  YSH(i) < YSH(i-1) && YSH(i) < YSH(i+1)
            XSHRough(end+1) = SweepRange(i);
        end
    end
end
% figure,plot(SweepRange,20*log10(abs(Y)))
if  isempty(XRough) && isempty(XSHRough)
    String = 'No higher order modes found!';
    OutputWindow1aUI2.String = String;
    OutputWindow1bUI2.String = '';
    OutputWindow2aUI2.String = '';
    OutputWindow2bUI2.String = '';
    disp(String)
    return
end
X = Converger(MatrixMethod,Decoupled,0,XRough,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,I1,Bisections,PhaseVelocity,SweepRange,Alpha,V,W,D3,D4,D5,AlphaSH,DSH);
if  Decoupled
    XSH = Converger(MatrixMethod,Decoupled,1,XSHRough,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,I1,Bisections,PhaseVelocity,SweepRange,Alpha,V,W,D3,D4,D5,AlphaSH,DSH);
end
if  ~Decoupled
    if  SuperLayerSize == 1 || SymmetricSystem
        [HSLamb,HALamb] = SymmetryChecker_Coupled(X,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,PhaseVelocity,Alpha,V,W,D3,D4,D5,HSLamb,HALamb);
    else
        HBLamb = X;
    end
else
    if  SuperLayerSize == 1 || SymmetricSystem
        [HSLamb,HALamb] = SymmetryChecker_Decoupled(X,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,PhaseVelocity,Alpha,W,D3,D5,HSLamb,HALamb);
        [HSShear,HAShear] = SymmetryChecker_Decoupled_SH(XSH,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,PhaseVelocity,AlphaSH,DSH,HSShear,HAShear);
    else
        HBLamb = X;
        HBShear = XSH;
    end
end
ModeNameLength = 1+5;
if  SuperLayerSize == 1 || SymmetricSystem
    if  ~Decoupled
        String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
        if  Symmetric
            if  any(HALamb)
                for i = 1:length(HALamb)
                    String = append(String,newline,pad(['A',num2str(i)],ModeNameLength),pad(num2str(HALamb(i),'%.3f'),11,'left'));
                end
            end
        else
            HBLamb = sort(horzcat(HALamb,HSLamb));
            if  any(HBLamb)
                for i = 1:length(HBLamb)
                    String = append(String,newline,pad(['B',num2str(i+2)],ModeNameLength),pad(num2str(HBLamb(i),'%.3f'),11,'left'));
                end
            end
        end
        OutputWindow1aUI2.String = String;
        disp(String)
        String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
        if  any(HSLamb) && Symmetric
            for i = 1:length(HSLamb)
                String = append(String,newline,pad(['S',num2str(i+1)],ModeNameLength),pad(num2str(HSLamb(i),'%.3f'),11,'left'));
            end
        end
        OutputWindow1bUI2.String = String;
        disp(extractAfter(String,')'))
        disp(' ')
        if  Symmetric
            if  any(HALamb)
                disp(['A: ',num2str(length(HALamb))])
                OutputWindow2aUI2.String = ['A: ',num2str(length(HALamb))];
            else
                OutputWindow2aUI2.String = '';
            end
        else
            if  any(HBLamb)
                disp(['B: ',num2str(length(HBLamb))])
                OutputWindow2aUI2.String = ['B: ',num2str(length(HBLamb))];
            else
                OutputWindow2aUI2.String = '';
            end
        end
        if  any(HSLamb) && Symmetric
            disp(['S: ',num2str(length(HSLamb))])
            OutputWindow2bUI2.String = ['S: ',num2str(length(HSLamb))];
        else
            OutputWindow2bUI2.String = '';
        end
    else
        String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
        if  Symmetric
            if  any(HALamb)
                for i = 1:length(HALamb)
                    String = append(String,newline,pad(['A',num2str(i)],ModeNameLength),pad(num2str(HALamb(i),'%.3f'),11,'left'));
                end
            end
        else
            HBLamb = sort(horzcat(HALamb,HSLamb));
            if  any(HBLamb)
                for i = 1:length(HBLamb)
                    String = append(String,newline,pad(['B',num2str(i+1)],ModeNameLength),pad(num2str(HBLamb(i),'%.3f'),11,'left'));
                end
            end
        end
        String = append(String,newline,' ');
        if  any(HAShear)
            for i = 1:length(HAShear)
                String = append(String,newline,pad(['ASH',num2str(i)],ModeNameLength),pad(num2str(HAShear(i),'%.3f'),11,'left'));
            end
        end
        OutputWindow1aUI2.String = String;
        disp(String)
        String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
        if  any(HSLamb) && Symmetric
            for i = 1:length(HSLamb)
                String = append(String,newline,pad(['S',num2str(i)],ModeNameLength),pad(num2str(HSLamb(i),'%.3f'),11,'left'));
            end
        end
        String = append(String,newline,' ');
        if  any(HSShear)
            for i = 1:length(HSShear)
                String = append(String,newline,pad(['SSH',num2str(i)],ModeNameLength),pad(num2str(HSShear(i),'%.3f'),11,'left'));
            end
        end
        OutputWindow1bUI2.String = String;
        disp(extractAfter(String,')'))
        disp(' ')
        String = '';
        if  Symmetric
            if  any(HALamb)
                String = append(String,'A:   ',num2str(length(HALamb)));
            end
        else
            if  any(HBLamb)
                String = append(String,'B:   ',num2str(length(HBLamb)));
            end
        end
        if  any(HAShear)
            String = append(String,newline,'ASH: ',num2str(length(HAShear)));
        end
        OutputWindow2aUI2.String = String;
        disp(String)
        String = '';
        if  any(HSLamb) && Symmetric
            String = append(String,'S:   ',num2str(length(HSLamb)));
        end
        if  any(HSShear)
            String = append(String,newline,'SSH: ',num2str(length(HSShear)));
        end
        OutputWindow2bUI2.String = String;
        disp(String)
    end
else
    if  ~Decoupled
        String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
        if  any(HBLamb)
            for i = 1:length(HBLamb)
                String = append(String,newline,pad(['B',num2str(i+2)],ModeNameLength),pad(num2str(HBLamb(i),'%.3f'),11,'left'));
            end
        end
        OutputWindow1aUI2.String = String;
        disp(String)
        disp(' ')
        if  any(HBLamb)
            disp(['B: ',num2str(length(HBLamb))])
            OutputWindow2aUI2.String = ['B: ',num2str(length(HBLamb))];
        end
        OutputWindow1bUI2.String = '';
        OutputWindow2bUI2.String = '';
    else
        String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
        if  any(HBLamb)
            for i = 1:length(HBLamb)
                String = append(String,newline,pad(['B',num2str(i+1)],ModeNameLength),pad(num2str(HBLamb(i),'%.3f'),11,'left'));
            end
        end
        OutputWindow1aUI2.String = String;
        disp(String)
        String = ['Modes @ ',num2str(PhaseVelocity/1e3),' m/ms:',newline,'Mode     Frq.(kHz)'];
        if  any(HBShear)
            for i = 1:length(HBShear)
                String = append(String,newline,pad(['BSH',num2str(i)],ModeNameLength),pad(num2str(HBShear(i),'%.3f'),11,'left'));
            end
        end
        OutputWindow1bUI2.String = String;
        disp(extractAfter(String,')'))
        disp(' ')
        if  any(HBLamb)
            disp(['B:   ',num2str(length(HBLamb))])
            OutputWindow2aUI2.String = ['B:   ',num2str(length(HBLamb))];
        else
            OutputWindow2aUI2.String = '';
        end
        if  any(HBShear)
            disp(['BSH: ',num2str(length(HBShear))])
            OutputWindow2bUI2.String = ['BSH: ',num2str(length(HBShear))];
        else
            OutputWindow2bUI2.String = '';
        end
    end
end
disp('----------------')
end
function X = Converger(MatrixMethod,Decoupled,SHMode,XRough,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,I1,Bisections,PhaseVelocity,SweepRange,Alpha,V,W,D3,D4,D5,AlphaSH,DSH)
    X = [];
    for j = 1:length(XRough)
        Frequency = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            Wavenumber = 2*pi*Frequency*1e3/PhaseVelocity;
            if  ~Decoupled
                [Y,Yc] = Computer_Coupled(MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,Alpha,V,W,D3,D4,D5);
            else
                if  ~SHMode
                    [Y,Yc] = Computer_Decoupled(MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,Alpha,W,D3,D5);
                else
                    [Y,Yc] = Computer_Decoupled_SH(Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,AlphaSH,DSH);
                end
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
function [Y,Yc] = Computer_Coupled(MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,Alpha,V,W,D3,D4,D5)
    Length = length(Wavenumber);
    V = repmat(V,1,1,Length);
    W = repmat(W,1,1,Length);
    D3 = repmat(D3,1,1,Length);
    D4 = repmat(D4,1,1,Length);
    D5 = repmat(D5,1,1,Length);
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Yc = NaN(Size);
    for m = 1:SuperLayerSize
        Phi = 1i*Wavenumber.*Alpha(m,:)*LayerThicknesses(m);
        E = exp(Phi);
        if  MatrixMethod == 1
            E_ = exp(-Phi);
            L1 = [E E_;V(m,:,:).*E V(m,:,:).*E_;W(m,:,:).*E -W(m,:,:).*E_;D3(m,:,:).*E D3(m,:,:).*E_;D5(m,:,:).*E -D5(m,:,:).*E_;D4(m,:,:).*E -D4(m,:,:).*E_];
            L2 = [ones(1,6,Length);V(m,:,:) V(m,:,:);W(m,:,:) -W(m,:,:);D3(m,:,:) D3(m,:,:);D5(m,:,:) -D5(m,:,:);D4(m,:,:) -D4(m,:,:)];
        elseif MatrixMethod == 2
            L1 = [D3(m,:,:) D3(m,:,:).*E;D5(m,:,:) -D5(m,:,:).*E;D4(m,:,:) -D4(m,:,:).*E;D3(m,:,:).*E D3(m,:,:);D5(m,:,:).*E -D5(m,:,:);D4(m,:,:).*E -D4(m,:,:)];
            L2 = [ones(1,3,Length) E;V(m,:,:) V(m,:,:).*E;W(m,:,:) -W(m,:,:).*E;E ones(1,3,Length);V(m,:,:).*E V(m,:,:);W(m,:,:).*E -W(m,:,:)];
        end
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    if  MatrixMethod == 1
        for m = 2:SuperLayerSize
            M{1} = pagemtimes(M{1},L{m});
        end
        for m = 1:length(Pattern)
            M{m+1} = pagemtimes(M{m},M{Pattern(m)});
        end
        if  SymmetricSystem
            M2{1} = L{end};
            for m = SuperLayerSize-1:-1:1
                M2{1} = pagemtimes(M2{1},L{m});
            end
            for m = 1:length(Pattern)
                M2{m+1} = pagemtimes(M2{m},M2{Pattern(m)});
            end
            M{end} = pagemtimes(M{end},M2{end});
        end
        for j = 1:Length
            Yc(j) = imag(det(M{end}(4:6,1:3,j)));
        end
    elseif MatrixMethod == 2
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
        for j = 1:Length
            Yc(j) = real(det(M{end}(:,:,j)));
        end
    end
    Yc = reshape(Yc,Size);
    Y = abs(Yc);
end
function [Y,Yc] = Computer_Decoupled(MatrixMethod,Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,Alpha,W,D3,D5)
    Length = length(Wavenumber);    
    W = repmat(W,1,1,Length);
    D3 = repmat(D3,1,1,Length);
    D5 = repmat(D5,1,1,Length);
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Yc = NaN(Size);
    for m = 1:SuperLayerSize
        Phi = 1i*Wavenumber.*Alpha(m,:)*LayerThicknesses(m);
        E = exp(Phi);
        if  MatrixMethod == 1
            E_ = exp(-Phi);
            L1 = [E E_;W(m,:,:).*E -W(m,:,:).*E_;D3(m,:,:).*E D3(m,:,:).*E_;D5(m,:,:).*E -D5(m,:,:).*E_];
            L2 = [ones(1,4,Length);W(m,:,:) -W(m,:,:);D3(m,:,:) D3(m,:,:);D5(m,:,:) -D5(m,:,:)];
        elseif MatrixMethod == 2
            L1 = [D3(m,:,:) D3(m,:,:).*E;D5(m,:,:) -D5(m,:,:).*E;D3(m,:,:).*E D3(m,:,:);D5(m,:,:).*E -D5(m,:,:)];
            L2 = [ones(1,2,Length) E;W(m,:,:) -W(m,:,:).*E;E ones(1,2,Length);W(m,:,:).*E -W(m,:,:)];
        end
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    if  MatrixMethod == 1
        for m = 2:SuperLayerSize
            M{1} = pagemtimes(M{1},L{m});
        end
        for m = 1:length(Pattern)
            M{m+1} = pagemtimes(M{m},M{Pattern(m)});
        end
        if  SymmetricSystem
            M2{1} = L{end};
            for m = SuperLayerSize-1:-1:1
                M2{1} = pagemtimes(M2{1},L{m});
            end
            for m = 1:length(Pattern)
                M2{m+1} = pagemtimes(M2{m},M2{Pattern(m)});
            end
            M{end} = pagemtimes(M{end},M2{end});
        end
        for j = 1:Length
            Yc(j) = real(det(M{end}(3:4,1:2,j)));
        end
    elseif MatrixMethod == 2
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
        for j = 1:Length
            Yc(j) = real(det(M{end}(:,:,j)));
        end
    end
    Yc = reshape(Yc,Size);
    Y = abs(Yc);
end
function [Y,Yc] = Computer_Decoupled_SH(Wavenumber,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,Alpha,D)
    D = repmat(D,1,1,length(Wavenumber));
    Size = size(Wavenumber);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    for m = 1:SuperLayerSize
        G = Wavenumber*Alpha(m)*LayerThicknesses(m);
        CosG = cos(G);
        SinG = sin(G);
        L{m} = [CosG 1i*SinG./D(m,1,:);1i*SinG.*D(m,1,:) CosG];
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M{1} = pagemtimes(M{1},L{m});
    end
    for m = 1:length(Pattern)
        M{m+1} = pagemtimes(M{m},M{Pattern(m)});
    end
    if  SymmetricSystem
        M2{1} = L{end};
        for m = SuperLayerSize-1:-1:1
            M2{1} = pagemtimes(M2{1},L{m});
        end
        for m = 1:length(Pattern)
            M2{m+1} = pagemtimes(M2{m},M2{Pattern(m)});
        end
        M{end} = pagemtimes(M{end},M2{end});
    end
    Yc = reshape(M{end}(2,1,:),Size);
    Y = abs(Yc);
end
function [HSLamb,HALamb] = SymmetryChecker_Coupled(X,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I,PhaseVelocity,Alpha,V,W,D3,D4,D5,HSLamb,HALamb)
    for p = 1:length(X)
        Wavenumber = 2*pi*X(p)*1e3/PhaseVelocity;
        for m = 1:SuperLayerSize
            E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
            L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
            L2 = [ones(1,3) E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E ones(1,3);V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
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
        Z1 = [-M{end}(1:3,1:3)*[1 1 1;V(1,:);-W(1,:)] -M{end}(1:3,4:6)*[1 1 1;V(1,:);W(1,:)];M{end}(4:6,1:3)*[1 1 1;V(1,:);-W(1,:)] M{end}(4:6,4:6)*[1 1 1;V(1,:);W(1,:)]];
        Z2 = [M{end}(1:3,1:3)*[1 1 1;V(1,:);W(1,:)];-M{end}(4:6,1:3)*[1 1 1;V(1,:);W(1,:)]];
        RT = Z1\Z2(:,1);
        u1 = [1+sum(RT(1:3)) sum(RT(4:6))];
        u1 = u1*exp(-1i*angle(u1(1)));
        if  sign(real(u1(1))) == sign(real(u1(2)))
            HSLamb(end+1) = X(p);
        else
            HALamb(end+1) = X(p);
        end
    end
end
function [HSLamb,HALamb] = SymmetryChecker_Decoupled(X,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,I1,PhaseVelocity,Alpha,W,D3,D5,HSLamb,HALamb)
    for p = 1:length(X)
        Wavenumber = 2*pi*X(p)*1e3/PhaseVelocity;
        for m = 1:SuperLayerSize
            E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
            L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:)];
            L2 = [ones(1,2) E;W(m,:) -W(m,:).*E;E ones(1,2);W(m,:).*E -W(m,:)];
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
        Z1 = [-M{end}(1:2,1:2)*[1 1;-W(1,:)] -M{end}(1:2,3:4)*[1 1;W(1,:)];M{end}(3:4,1:2)*[1 1;-W(1,:)] M{end}(3:4,3:4)*[1 1;W(1,:)]];
        Z2 = [M{end}(1:2,1:2)*[1 1;W(1,:)];-M{end}(3:4,1:2)*[1 1;W(1,:)]];
        RT = Z1\Z2(:,1);
        u1 = [1+sum(RT(1:2)) sum(RT(3:4))];
        u1 = u1*exp(-1i*angle(u1(1)));
        if  sign(real(u1(1))) == sign(real(u1(2)))
            HSLamb(end+1) = X(p);
        else
            HALamb(end+1) = X(p);
        end
    end
end
function [HSShear,HAShear] = SymmetryChecker_Decoupled_SH(X,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,PhaseVelocity,Alpha,D,HSShear,HAShear)
    for p = 1:length(X)
        Wavenumber = 2*pi*X(p)*1e3/PhaseVelocity;
        for m = 1:SuperLayerSize
            E = exp(1i*Wavenumber*Alpha(m)*LayerThicknesses(m));
            E2 = E^2;
            L{m} = Wavenumber*D(m)/(E2-1)*[-1-E2 2*E;-2*E 1+E2];
        end
        M{1} = L{1};
        for m = 2:SuperLayerSize
            M0 = L{m}(1,1)-M{1}(2,2);
            M1 = M{1}(1,2)/M0;
            M2 = L{m}(2,1)/M0;
            M{1} = [M{1}(1,1)+M1*M{1}(2,1) -M1*L{m}(1,2);M2*M{1}(2,1) L{m}(2,2)-M2*L{m}(1,2)];
        end
        for m = 1:length(Pattern)
            M0 = M{Pattern(m)}(1,1)-M{m}(2,2);
            M1 = M{m}(1,2)/M0;
            M2 = M{Pattern(m)}(2,1)/M0;
            M{m+1} = [M{m}(1,1)+M1*M{m}(2,1) -M1*M{Pattern(m)}(1,2);M2*M{m}(2,1) M{Pattern(m)}(2,2)-M2*M{Pattern(m)}(1,2)];
        end
        if  SymmetricSystem
            M1 = -M{end}(1,2)/(2*M{end}(2,2));
            M{end} = [M{end}(1,1)+M1*M{end}(2,1) M1*M{end}(2,1);-M1*M{end}(2,1) -M{end}(1,1)-M1*M{end}(2,1)];
        end
        Z1 = [M{end}(1,1) -M{end}(1,2);M{end}(2,:)]; % changed sign in Z1(1,1)!
        Z2 = [M{end}(1,1);-M{end}(2,1)];
        RT = Z1\Z2;
        u1 = [1+RT(1) RT(2)];
        u1 = u1*exp(-1i*angle(u1(1)));
        if  sign(real(u1(1))) == sign(real(u1(2)))
            HSShear(end+1) = X(p);
        else
            HAShear(end+1) = X(p);
        end
    end
end