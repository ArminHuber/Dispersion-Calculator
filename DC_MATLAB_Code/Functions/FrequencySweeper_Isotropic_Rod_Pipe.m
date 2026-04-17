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
function [HF,HL,HT] = FrequencySweeper_Isotropic_Rod_Pipe(Geometry,Material,Ro,Ri,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,SweepRange,FlexuralModeOrders,PhaseVelocity,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Bisections = 17; % ceil(log2(1e5))

%#ok<*AGROW>
HF{FlexuralModeOrders}=[];HL=[];HT=[];XRoughF{FlexuralModeOrders}=[];XRoughL=[];XRoughT=[];
if  PhaseVelocity > Material.LongitudinalVelocity
    Region = 1;
elseif PhaseVelocity <= Material.LongitudinalVelocity && PhaseVelocity >= Material.TransverseVelocity
    Region = 2;
elseif PhaseVelocity < Material.TransverseVelocity
    Region = 3;
end
cL2 = Material.LongitudinalVelocity^2;
cT2 = Material.TransverseVelocity^2;
cFi2 = InnerFluid.Velocity^2;
cFo2 = OuterFluid.Velocity^2;
Ro2 = Ro^2;
Ri2 = Ri^2;
Densityi = InnerFluid.Density/Material.Density;
Densityo = OuterFluid.Density/Material.Density;
AngularFrequency = 2*pi*SweepRange*1e3;
AngularFrequency2 = AngularFrequency.^2;
PhaseVelocity2 = PhaseVelocity^2;
Wavenumber = AngularFrequency/PhaseVelocity;
k2 = Wavenumber.^2;
kT2 = AngularFrequency2/cT2;
y2 = kT2-k2;
if  strcmp(Geometry,'Rod')
    kL2 = AngularFrequency2/cL2;
    x = sqrt(kL2-k2);
    y = sqrt(y2);
    xRo = x*Ro;
    yRo = y*Ro;
    zRo = sqrt(AngularFrequency2/cFo2-k2)*Ro;
    Y = abs(yRo-2*besselj(1,yRo)./besselj(0,yRo));
    for i = 2:length(SweepRange)-1
        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
            XRoughT(end+1) = SweepRange(i);
        end
    end
    for j = 1:length(XRoughT)
        Frequency = XRoughT(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            yRoT = sqrt(1/cT2-1/PhaseVelocity2).*2*pi*Frequency*1e3*Ro;
            Yc = yRoT-2*besselj(1,yRoT)./besselj(0,yRoT);
            Y = abs(Yc);
            for i = 2:length(Frequency)-1
                if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
                    if  o == Bisections && sign(Yc(i-1)) ~= sign(Yc(i+1))
                        HT(end+1) = Frequency(i);
                    end
                    Frequency = [Frequency(i-1) Frequency(i+1)];
                    break
                end
            end
        end
    end
    Y = abs(2*x/Ro.*kT2-(y2-k2).^2.*besselj(0,xRo)./besselj(1,xRo)-4*k2.*x.*y.*besselj(0,yRo)./besselj(1,yRo));
elseif strcmp(Geometry,'Pipe')
    if  Region == 1
        x = sqrt(AngularFrequency2/cL2-k2);
        y = sqrt(y2);
    elseif Region == 2
        x = sqrt(k2-AngularFrequency2/cL2);
        y = sqrt(y2);
    elseif Region == 3
        x = sqrt(k2-AngularFrequency2/cL2);
        y = sqrt(-y2);
    end
    xRo = x*Ro;
    xRi = x*Ri;
    yRo = y*Ro;
    yRi = y*Ri;
    zRo = sqrt(AngularFrequency2/cFo2-k2)*Ro;
    zRi = sqrt(AngularFrequency2/cFi2-k2)*Ri;
    Y = abs(besselj(2,yRi).*bessely(2,yRo)-besselj(2,yRo).*bessely(2,yRi));
    for i = 2:length(SweepRange)-1
        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
            XRoughT(end+1) = SweepRange(i);
        end
    end
    for j = 1:length(XRoughT)
        Frequency = XRoughT(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            y = sqrt(1/cT2-1/PhaseVelocity2).*2*pi*Frequency*1e3;
            yRoT = y*Ro;
            yRiT = y*Ri;
            Yc = besselj(2,yRiT).*bessely(2,yRoT)-besselj(2,yRoT).*bessely(2,yRiT);
            Y = abs(Yc);
            for i = 2:length(Frequency)-1
                if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
                    if  o == Bisections && sign(Yc(i-1)) ~= sign(Yc(i+1))
                        HT(end+1) = Frequency(i);
                    end
                    Frequency = [Frequency(i-1) Frequency(i+1)];
                    break
                end
            end
        end
    end
    Y = Computer_Pipe_L(Wavenumber,Region,k2,kT2,Ro2,Ri2,y2,xRo,yRo,zRo,xRi,yRi,zRi,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink);
end
for i = 2:length(SweepRange)-1
    if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
        XRoughL(end+1) = SweepRange(i);
    end
end
% figure,plot(SweepRange,20*log10(Y))
for n = 1:FlexuralModeOrders
    if  strcmp(Geometry,'Rod')
        Y = Computer_Rod_F(n,Wavenumber,k2,kT2,Ro2,y2,xRo,yRo);
    elseif strcmp(Geometry,'Pipe')
        Y = Computer_Pipe_F(n,Wavenumber,Region,k2,kT2,Ro2,Ri2,y2,xRo,yRo,zRo,xRi,yRi,zRi,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink);
    end
    for i = 2:length(SweepRange)-1
        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
            XRoughF{n}(end+1) = SweepRange(i);
        end
    end
    if  isempty(XRoughF{n})
        XRoughF(n:end) = [];
        break
    end
% figure,plot(SweepRange,20*log10(Y))
end
if  isempty(XRoughL) && isempty(XRoughF)
    disp('No higher order modes found!')
    return
end
HL = Converger(Geometry,XRoughL,Bisections,SweepRange,cL2,cT2,cFo2,cFi2,0,Region,PhaseVelocity,Ro,Ri,Ro2,Ri2,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink);
for n = 1:length(XRoughF)
    HF{n} = Converger(Geometry,XRoughF{n},Bisections,SweepRange,cL2,cT2,cFo2,cFi2,n,Region,PhaseVelocity,Ro,Ri,Ro2,Ri2,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink);
end
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
function X = Converger(Geometry,XRough,Bisections,SweepRange,cL2,cT2,cFo2,cFi2,n,Region,PhaseVelocity,Ro,Ri,Ro2,Ri2,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink)
    X = [];
    for j = 1:length(XRough)
        Frequency = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            AngularFrequency = 2*pi*Frequency*1e3;
            AngularFrequency2 = AngularFrequency.^2;
            Wavenumber = AngularFrequency/PhaseVelocity;
            k2 = Wavenumber.^2;
            kT2 = AngularFrequency2/cT2;
            y2 = kT2-k2;
            if  strcmp(Geometry,'Rod')
                kL2 = AngularFrequency2/cL2;
                x = sqrt(kL2-k2);
                y = sqrt(y2);
                xRo = x*Ro;
                yRo = y*Ro;
                if  n == 0
                    Yc = 2*x/Ro.*kT2-(y2-k2).^2.*besselj(0,xRo)./besselj(1,xRo)-4*k2.*x.*y.*besselj(0,yRo)./besselj(1,yRo);
                    Y = abs(Yc);
                else
                    [Y,Yc] = Computer_Rod_F(n,Wavenumber,k2,kT2,Ro2,y2,xRo,yRo);
                end
            elseif strcmp(Geometry,'Pipe')
                if  Region == 1
                    x = sqrt(AngularFrequency2/cL2-k2);
                    y = sqrt(y2);
                elseif Region == 2
                    x = sqrt(k2-AngularFrequency2/cL2);
                    y = sqrt(y2);
                elseif Region == 3
                    x = sqrt(k2-AngularFrequency2/cL2);
                    y = sqrt(-y2);
                end
                xRo = x*Ro;
                xRi = x*Ri;
                yRo = y*Ro;
                yRi = y*Ri;
                zRo = sqrt(AngularFrequency2/cFo2-k2)*Ro;
                zRi = sqrt(AngularFrequency2/cFi2-k2)*Ri;
                if  n == 0
                    [Y,Yc] = Computer_Pipe_L(Wavenumber,Region,k2,kT2,Ro2,Ri2,y2,xRo,yRo,zRo,xRi,yRi,zRi,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink);
                else
                    [Y,Yc] = Computer_Pipe_F(n,Wavenumber,Region,k2,kT2,Ro2,Ri2,y2,xRo,yRo,zRo,xRi,yRi,zRi,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink);
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
function [Y,Yc] = Computer_Pipe_F(n,Wavenumber,Region,k2,kT2,Ro2,Ri2,y2,xRo,yRo,zRo,xRi,yRi,zRi,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink)
    n2 = n^2;
    if  Region == 1
        Znxi = besselj(n,xRi);
        Znxo = besselj(n,xRo);
        Znyi = besselj(n,yRi);
        Znyo = besselj(n,yRo);
        Wnxi = bessely(n,xRi);
        Wnxo = bessely(n,xRo);
        Wnyi = bessely(n,yRi);
        Wnyo = bessely(n,yRo);
        dZnxiRi = n*Znxi-xRi.*besselj(n+1,xRi);
        dZnxoRo = n*Znxo-xRo.*besselj(n+1,xRo);
        dZnyiRi = n*Znyi-yRi.*besselj(n+1,yRi);
        dZnyoRo = n*Znyo-yRo.*besselj(n+1,yRo);
        dWnxiRi = n*Wnxi-xRi.*bessely(n+1,xRi);
        dWnxoRo = n*Wnxo-xRo.*bessely(n+1,xRo);
        dWnyiRi = n*Wnyi-yRi.*bessely(n+1,yRi);
        dWnyoRo = n*Wnyo-yRo.*bessely(n+1,yRo);
    elseif Region == 2
        Znxi = besseli(n,xRi);
        Znxo = besseli(n,xRo);
        Znyi = besselj(n,yRi);
        Znyo = besselj(n,yRo);
        Wnxi = besselk(n,xRi);
        Wnxo = besselk(n,xRo);
        Wnyi = bessely(n,yRi);
        Wnyo = bessely(n,yRo);
        dZnxiRi = n*Znxi+xRi.*besseli(n+1,xRi);
        dZnxoRo = n*Znxo+xRo.*besseli(n+1,xRo);
        dZnyiRi = n*Znyi-yRi.*besselj(n+1,yRi);
        dZnyoRo = n*Znyo-yRo.*besselj(n+1,yRo);
        dWnxiRi = n*Wnxi-xRi.*besselk(n+1,xRi);
        dWnxoRo = n*Wnxo-xRo.*besselk(n+1,xRo);
        dWnyiRi = n*Wnyi-yRi.*bessely(n+1,yRi);
        dWnyoRo = n*Wnyo-yRo.*bessely(n+1,yRo);
    elseif Region == 3
        Znxi = besseli(n,xRi);
        Znxo = besseli(n,xRo);
        Znyi = besseli(n,yRi);
        Znyo = besseli(n,yRo);
        Wnxi = besselk(n,xRi);
        Wnxo = besselk(n,xRo);
        Wnyi = besselk(n,yRi);
        Wnyo = besselk(n,yRo);
        dZnxiRi = n*Znxi+xRi.*besseli(n+1,xRi);
        dZnxoRo = n*Znxo+xRo.*besseli(n+1,xRo);
        dZnyiRi = n*Znyi+yRi.*besseli(n+1,yRi);
        dZnyoRo = n*Znyo+yRo.*besseli(n+1,yRo);
        dWnxiRi = n*Wnxi-xRi.*besselk(n+1,xRi);
        dWnxoRo = n*Wnxo-xRo.*besselk(n+1,xRo);
        dWnyiRi = n*Wnyi-yRi.*besselk(n+1,yRi);
        dWnyoRo = n*Wnyo-yRo.*besselk(n+1,yRo);
    end
    if  ToggleInnerFluid
        if  Sink
            Znzi = besselh(n,2,zRi);
            dZnziRi = n*Znzi-zRi.*besselh(n+1,2,zRi);
        else
            Znzi = besselj(n,zRi);
            dZnziRi = n*Znzi-zRi.*besselj(n+1,zRi);
        end
    end
    if  ToggleOuterFluid
        Hnzo = besselh(n,zRo);
        dHnzoRo = n*Hnzo-zRo.*besselh(n+1,zRo);
    end
    Yc = NaN(size(Wavenumber));
    for i = 1:length(Wavenumber)
        if  ToggleInnerFluid && ToggleOuterFluid
            M(1,1) = dZnxiRi(i)+(kT2(i)/2*Ri2-k2(i)*Ri2-n2)*Znxi(i);
            M(1,2) = dWnxiRi(i)+(kT2(i)/2*Ri2-k2(i)*Ri2-n2)*Wnxi(i);
            M(1,3) = -Wavenumber(i)*(dZnyiRi(i)+(y2(i)*Ri2-n2)*Znyi(i));
            M(1,4) = -Wavenumber(i)*(dWnyiRi(i)+(y2(i)*Ri2-n2)*Wnyi(i));
            M(1,5) = n*(dZnyiRi(i)-Znyi(i));
            M(1,6) = n*(dWnyiRi(i)-Wnyi(i));
            M(1,7) = kT2(i)/2*Ri2*Densityi*Znzi(i);
            M(2,1) = 2*n*(dZnxiRi(i)-Znxi(i));
            M(2,2) = 2*n*(dWnxiRi(i)-Wnxi(i));
            M(2,3) = 2*Wavenumber(i)*n*(Znyi(i)-dZnyiRi(i));
            M(2,4) = 2*Wavenumber(i)*n*(Wnyi(i)-dWnyiRi(i));
            M(2,5) = 2*dZnyiRi(i)+(y2(i)*Ri2-2*n2)*Znyi(i);
            M(2,6) = 2*dWnyiRi(i)+(y2(i)*Ri2-2*n2)*Wnyi(i);
            M(3,1) = 2*Wavenumber(i)*dZnxiRi(i);
            M(3,2) = 2*Wavenumber(i)*dWnxiRi(i);
            M(3,3) = (y2(i)-k2(i))*dZnyiRi(i);
            M(3,4) = (y2(i)-k2(i))*dWnyiRi(i);
            M(3,5) = -Wavenumber(i)*n*Znyi(i);
            M(3,6) = -Wavenumber(i)*n*Wnyi(i);
            M(4,1) = dZnxiRi(i);
            M(4,2) = dWnxiRi(i);
            M(4,3) = -Wavenumber(i)*dZnyiRi(i);
            M(4,4) = -Wavenumber(i)*dWnyiRi(i);
            M(4,5) = -n*Znyi(i);
            M(4,6) = -n*Wnyi(i);
            M(4,7) = dZnziRi(i);
            M(5,1) = dZnxoRo(i)+(kT2(i)/2*Ro2-k2(i)*Ro2-n2)*Znxo(i);
            M(5,2) = dWnxoRo(i)+(kT2(i)/2*Ro2-k2(i)*Ro2-n2)*Wnxo(i);
            M(5,3) = -Wavenumber(i)*(dZnyoRo(i)+(y2(i)*Ro2-n2)*Znyo(i));
            M(5,4) = -Wavenumber(i)*(dWnyoRo(i)+(y2(i)*Ro2-n2)*Wnyo(i));
            M(5,5) = n*(dZnyoRo(i)-Znyo(i));
            M(5,6) = n*(dWnyoRo(i)-Wnyo(i));
            M(5,8) = kT2(i)/2*Ro2*Densityo*Hnzo(i);
            M(6,1) = 2*n*(dZnxoRo(i)-Znxo(i));
            M(6,2) = 2*n*(dWnxoRo(i)-Wnxo(i));
            M(6,3) = 2*Wavenumber(i)*n*(Znyo(i)-dZnyoRo(i));
            M(6,4) = 2*Wavenumber(i)*n*(Wnyo(i)-dWnyoRo(i));
            M(6,5) = 2*dZnyoRo(i)+(y2(i)*Ro2-2*n2)*Znyo(i);
            M(6,6) = 2*dWnyoRo(i)+(y2(i)*Ro2-2*n2)*Wnyo(i);
            M(7,1) = 2*Wavenumber(i)*dZnxoRo(i);
            M(7,2) = 2*Wavenumber(i)*dWnxoRo(i);
            M(7,3) = (y2(i)-k2(i))*dZnyoRo(i);
            M(7,4) = (y2(i)-k2(i))*dWnyoRo(i);
            M(7,5) = -Wavenumber(i)*n*Znyo(i);
            M(7,6) = -Wavenumber(i)*n*Wnyo(i);
            M(8,1) = dZnxoRo(i);
            M(8,2) = dWnxoRo(i);
            M(8,3) = -Wavenumber(i)*dZnyoRo(i);
            M(8,4) = -Wavenumber(i)*dWnyoRo(i);
            M(8,5) = -n*Znyo(i);
            M(8,6) = -n*Wnyo(i);
            M(8,8) = dHnzoRo(i);
        elseif ToggleInnerFluid && ~ToggleOuterFluid
            M(1,1) = dZnxiRi(i)+(kT2(i)/2*Ri2-k2(i)*Ri2-n2)*Znxi(i);
            M(1,2) = dWnxiRi(i)+(kT2(i)/2*Ri2-k2(i)*Ri2-n2)*Wnxi(i);
            M(1,3) = -Wavenumber(i)*(dZnyiRi(i)+(y2(i)*Ri2-n2)*Znyi(i));
            M(1,4) = -Wavenumber(i)*(dWnyiRi(i)+(y2(i)*Ri2-n2)*Wnyi(i));
            M(1,5) = n*(dZnyiRi(i)-Znyi(i));
            M(1,6) = n*(dWnyiRi(i)-Wnyi(i));
            M(1,7) = kT2(i)/2*Ri2*Densityi*Znzi(i);
            M(2,1) = 2*n*(dZnxiRi(i)-Znxi(i));
            M(2,2) = 2*n*(dWnxiRi(i)-Wnxi(i));
            M(2,3) = 2*Wavenumber(i)*n*(Znyi(i)-dZnyiRi(i));
            M(2,4) = 2*Wavenumber(i)*n*(Wnyi(i)-dWnyiRi(i));
            M(2,5) = 2*dZnyiRi(i)+(y2(i)*Ri2-2*n2)*Znyi(i);
            M(2,6) = 2*dWnyiRi(i)+(y2(i)*Ri2-2*n2)*Wnyi(i);
            M(3,1) = 2*Wavenumber(i)*dZnxiRi(i);
            M(3,2) = 2*Wavenumber(i)*dWnxiRi(i);
            M(3,3) = (y2(i)-k2(i))*dZnyiRi(i);
            M(3,4) = (y2(i)-k2(i))*dWnyiRi(i);
            M(3,5) = -Wavenumber(i)*n*Znyi(i);
            M(3,6) = -Wavenumber(i)*n*Wnyi(i);
            M(4,1) = dZnxiRi(i);
            M(4,2) = dWnxiRi(i);
            M(4,3) = -Wavenumber(i)*dZnyiRi(i);
            M(4,4) = -Wavenumber(i)*dWnyiRi(i);
            M(4,5) = -n*Znyi(i);
            M(4,6) = -n*Wnyi(i);
            M(4,7) = dZnziRi(i);
            M(5,1) = dZnxoRo(i)+(kT2(i)/2*Ro2-k2(i)*Ro2-n2)*Znxo(i);
            M(5,2) = dWnxoRo(i)+(kT2(i)/2*Ro2-k2(i)*Ro2-n2)*Wnxo(i);
            M(5,3) = -Wavenumber(i)*(dZnyoRo(i)+(y2(i)*Ro2-n2)*Znyo(i));
            M(5,4) = -Wavenumber(i)*(dWnyoRo(i)+(y2(i)*Ro2-n2)*Wnyo(i));
            M(5,5) = n*(dZnyoRo(i)-Znyo(i));
            M(5,6) = n*(dWnyoRo(i)-Wnyo(i));
            M(6,1) = 2*n*(dZnxoRo(i)-Znxo(i));
            M(6,2) = 2*n*(dWnxoRo(i)-Wnxo(i));
            M(6,3) = 2*Wavenumber(i)*n*(Znyo(i)-dZnyoRo(i));
            M(6,4) = 2*Wavenumber(i)*n*(Wnyo(i)-dWnyoRo(i));
            M(6,5) = 2*dZnyoRo(i)+(y2(i)*Ro2-2*n2)*Znyo(i);
            M(6,6) = 2*dWnyoRo(i)+(y2(i)*Ro2-2*n2)*Wnyo(i);
            M(7,1) = 2*Wavenumber(i)*dZnxoRo(i);
            M(7,2) = 2*Wavenumber(i)*dWnxoRo(i);
            M(7,3) = (y2(i)-k2(i))*dZnyoRo(i);
            M(7,4) = (y2(i)-k2(i))*dWnyoRo(i);
            M(7,5) = -Wavenumber(i)*n*Znyo(i);
            M(7,6) = -Wavenumber(i)*n*Wnyo(i);
        elseif ~ToggleInnerFluid && ToggleOuterFluid
            M(1,1) = dZnxiRi(i)+(kT2(i)/2*Ri2-k2(i)*Ri2-n2)*Znxi(i);
            M(1,2) = dWnxiRi(i)+(kT2(i)/2*Ri2-k2(i)*Ri2-n2)*Wnxi(i);
            M(1,3) = -Wavenumber(i)*(dZnyiRi(i)+(y2(i)*Ri2-n2)*Znyi(i));
            M(1,4) = -Wavenumber(i)*(dWnyiRi(i)+(y2(i)*Ri2-n2)*Wnyi(i));
            M(1,5) = n*(dZnyiRi(i)-Znyi(i));
            M(1,6) = n*(dWnyiRi(i)-Wnyi(i));
            M(2,1) = 2*n*(dZnxiRi(i)-Znxi(i));
            M(2,2) = 2*n*(dWnxiRi(i)-Wnxi(i));
            M(2,3) = 2*Wavenumber(i)*n*(Znyi(i)-dZnyiRi(i));
            M(2,4) = 2*Wavenumber(i)*n*(Wnyi(i)-dWnyiRi(i));
            M(2,5) = 2*dZnyiRi(i)+(y2(i)*Ri2-2*n2)*Znyi(i);
            M(2,6) = 2*dWnyiRi(i)+(y2(i)*Ri2-2*n2)*Wnyi(i);
            M(3,1) = 2*Wavenumber(i)*dZnxiRi(i);
            M(3,2) = 2*Wavenumber(i)*dWnxiRi(i);
            M(3,3) = (y2(i)-k2(i))*dZnyiRi(i);
            M(3,4) = (y2(i)-k2(i))*dWnyiRi(i);
            M(3,5) = -Wavenumber(i)*n*Znyi(i);
            M(3,6) = -Wavenumber(i)*n*Wnyi(i);
            M(4,1) = dZnxoRo(i)+(kT2(i)/2*Ro2-k2(i)*Ro2-n2)*Znxo(i);
            M(4,2) = dWnxoRo(i)+(kT2(i)/2*Ro2-k2(i)*Ro2-n2)*Wnxo(i);
            M(4,3) = -Wavenumber(i)*(dZnyoRo(i)+(y2(i)*Ro2-n2)*Znyo(i));
            M(4,4) = -Wavenumber(i)*(dWnyoRo(i)+(y2(i)*Ro2-n2)*Wnyo(i));
            M(4,5) = n*(dZnyoRo(i)-Znyo(i));
            M(4,6) = n*(dWnyoRo(i)-Wnyo(i));
            M(4,7) = kT2(i)/2*Ro2*Densityo*Hnzo(i);
            M(5,1) = 2*n*(dZnxoRo(i)-Znxo(i));
            M(5,2) = 2*n*(dWnxoRo(i)-Wnxo(i));
            M(5,3) = 2*Wavenumber(i)*n*(Znyo(i)-dZnyoRo(i));
            M(5,4) = 2*Wavenumber(i)*n*(Wnyo(i)-dWnyoRo(i));
            M(5,5) = 2*dZnyoRo(i)+(y2(i)*Ro2-2*n2)*Znyo(i);
            M(5,6) = 2*dWnyoRo(i)+(y2(i)*Ro2-2*n2)*Wnyo(i);
            M(6,1) = 2*Wavenumber(i)*dZnxoRo(i);
            M(6,2) = 2*Wavenumber(i)*dWnxoRo(i);
            M(6,3) = (y2(i)-k2(i))*dZnyoRo(i);
            M(6,4) = (y2(i)-k2(i))*dWnyoRo(i);
            M(6,5) = -Wavenumber(i)*n*Znyo(i);
            M(6,6) = -Wavenumber(i)*n*Wnyo(i);
            M(7,1) = dZnxoRo(i);
            M(7,2) = dWnxoRo(i);
            M(7,3) = -Wavenumber(i)*dZnyoRo(i);
            M(7,4) = -Wavenumber(i)*dWnyoRo(i);
            M(7,5) = -n*Znyo(i);
            M(7,6) = -n*Wnyo(i);
            M(7,7) = dHnzoRo(i);
        elseif ~ToggleInnerFluid && ~ToggleOuterFluid
            M(1,1) = dZnxiRi(i)+(kT2(i)/2*Ri2-k2(i)*Ri2-n2)*Znxi(i);
            M(1,2) = dWnxiRi(i)+(kT2(i)/2*Ri2-k2(i)*Ri2-n2)*Wnxi(i);
            M(1,3) = -Wavenumber(i)*(dZnyiRi(i)+(y2(i)*Ri2-n2)*Znyi(i));
            M(1,4) = -Wavenumber(i)*(dWnyiRi(i)+(y2(i)*Ri2-n2)*Wnyi(i));
            M(1,5) = n*(dZnyiRi(i)-Znyi(i));
            M(1,6) = n*(dWnyiRi(i)-Wnyi(i));
            M(2,1) = 2*n*(dZnxiRi(i)-Znxi(i));
            M(2,2) = 2*n*(dWnxiRi(i)-Wnxi(i));
            M(2,3) = 2*Wavenumber(i)*n*(Znyi(i)-dZnyiRi(i));
            M(2,4) = 2*Wavenumber(i)*n*(Wnyi(i)-dWnyiRi(i));
            M(2,5) = 2*dZnyiRi(i)+(y2(i)*Ri2-2*n2)*Znyi(i);
            M(2,6) = 2*dWnyiRi(i)+(y2(i)*Ri2-2*n2)*Wnyi(i);
            M(3,1) = 2*Wavenumber(i)*dZnxiRi(i);
            M(3,2) = 2*Wavenumber(i)*dWnxiRi(i);
            M(3,3) = (y2(i)-k2(i))*dZnyiRi(i);
            M(3,4) = (y2(i)-k2(i))*dWnyiRi(i);
            M(3,5) = -Wavenumber(i)*n*Znyi(i);
            M(3,6) = -Wavenumber(i)*n*Wnyi(i);
            M(4,1) = dZnxoRo(i)+(kT2(i)/2*Ro2-k2(i)*Ro2-n2)*Znxo(i);
            M(4,2) = dWnxoRo(i)+(kT2(i)/2*Ro2-k2(i)*Ro2-n2)*Wnxo(i);
            M(4,3) = -Wavenumber(i)*(dZnyoRo(i)+(y2(i)*Ro2-n2)*Znyo(i));
            M(4,4) = -Wavenumber(i)*(dWnyoRo(i)+(y2(i)*Ro2-n2)*Wnyo(i));
            M(4,5) = n*(dZnyoRo(i)-Znyo(i));
            M(4,6) = n*(dWnyoRo(i)-Wnyo(i));
            M(5,1) = 2*n*(dZnxoRo(i)-Znxo(i));
            M(5,2) = 2*n*(dWnxoRo(i)-Wnxo(i));
            M(5,3) = 2*Wavenumber(i)*n*(Znyo(i)-dZnyoRo(i));
            M(5,4) = 2*Wavenumber(i)*n*(Wnyo(i)-dWnyoRo(i));
            M(5,5) = 2*dZnyoRo(i)+(y2(i)*Ro2-2*n2)*Znyo(i);
            M(5,6) = 2*dWnyoRo(i)+(y2(i)*Ro2-2*n2)*Wnyo(i);
            M(6,1) = 2*Wavenumber(i)*dZnxoRo(i);
            M(6,2) = 2*Wavenumber(i)*dWnxoRo(i);
            M(6,3) = (y2(i)-k2(i))*dZnyoRo(i);
            M(6,4) = (y2(i)-k2(i))*dWnyoRo(i);
            M(6,5) = -Wavenumber(i)*n*Znyo(i);
            M(6,6) = -Wavenumber(i)*n*Wnyo(i);
        end
        Yc(i) = det(M);
    end
    Y = abs(Yc);
end
function [Y,Yc] = Computer_Pipe_L(Wavenumber,Region,k2,kT2,Ro2,Ri2,y2,xRo,yRo,zRo,xRi,yRi,zRi,Densityo,Densityi,ToggleInnerFluid,ToggleOuterFluid,Sink)
    if  Region == 1
        Z0xi = besselj(0,xRi);
        Z0xo = besselj(0,xRo);
        Z0yi = besselj(0,yRi);
        Z0yo = besselj(0,yRo);
        W0xi = bessely(0,xRi);
        W0xo = bessely(0,xRo);
        W0yi = bessely(0,yRi);
        W0yo = bessely(0,yRo);
        Z1xiRi = xRi.*besselj(1,xRi);
        Z1xoRo = xRo.*besselj(1,xRo);
        Z1yiRi = yRi.*besselj(1,yRi);
        Z1yoRo = yRo.*besselj(1,yRo);
        W1xiRi = xRi.*bessely(1,xRi);
        W1xoRo = xRo.*bessely(1,xRo);
        W1yiRi = yRi.*bessely(1,yRi);
        W1yoRo = yRo.*bessely(1,yRo);
    elseif Region == 2
        Z0xi = besseli(0,xRi);
        Z0xo = besseli(0,xRo);
        Z0yi = besselj(0,yRi);
        Z0yo = besselj(0,yRo);
        W0xi = besselk(0,xRi);
        W0xo = besselk(0,xRo);
        W0yi = bessely(0,yRi);
        W0yo = bessely(0,yRo);
        Z1xiRi = -xRi.*besseli(1,xRi);
        Z1xoRo = -xRo.*besseli(1,xRo);
        Z1yiRi = yRi.*besselj(1,yRi);
        Z1yoRo = yRo.*besselj(1,yRo);
        W1xiRi = xRi.*besselk(1,xRi);
        W1xoRo = xRo.*besselk(1,xRo);
        W1yiRi = yRi.*bessely(1,yRi);
        W1yoRo = yRo.*bessely(1,yRo);
    elseif Region == 3
        Z0xi = besseli(0,xRi);
        Z0xo = besseli(0,xRo);
        Z0yi = besseli(0,yRi);
        Z0yo = besseli(0,yRo);
        W0xi = besselk(0,xRi);
        W0xo = besselk(0,xRo);
        W0yi = besselk(0,yRi);
        W0yo = besselk(0,yRo);
        Z1xiRi = -xRi.*besseli(1,xRi);
        Z1xoRo = -xRo.*besseli(1,xRo);
        Z1yiRi = -yRi.*besseli(1,yRi);
        Z1yoRo = -yRo.*besseli(1,yRo);
        W1xiRi = xRi.*besselk(1,xRi);
        W1xoRo = xRo.*besselk(1,xRo);
        W1yiRi = yRi.*besselk(1,yRi);
        W1yoRo = yRo.*besselk(1,yRo);
    end
    if  ToggleInnerFluid
        if  Sink
            Z0zi = besselh(0,2,zRi);
            Z1ziRi = -zRi.*besselh(1,2,zRi);
        else
            Z0zi = besselj(0,zRi);
            Z1ziRi = -zRi.*besselj(1,zRi);
        end
    end
    if  ToggleOuterFluid
        H0zo = besselh(0,zRo);
        H1zoRo = -zRo.*besselh(1,zRo);
    end
    Yc = NaN(size(Wavenumber));
    for i = 1:length(Wavenumber)
        if  ToggleInnerFluid && ToggleOuterFluid
            M(1,1) = Ri2*(kT2(i)/2-k2(i))*Z0xi(i)-Z1xiRi(i);
            M(1,2) = Ri2*(kT2(i)/2-k2(i))*W0xi(i)-W1xiRi(i);
            M(1,3) = Wavenumber(i)*(Z1yiRi(i)-y2(i)*Ri2*Z0yi(i));
            M(1,4) = Wavenumber(i)*(W1yiRi(i)-y2(i)*Ri2*W0yi(i));
            M(1,5) = kT2(i)/2*Ri2*Densityi*Z0zi(i);
            M(2,1) = -2*Wavenumber(i)*Z1xiRi(i);
            M(2,2) = -2*Wavenumber(i)*W1xiRi(i);
            M(2,3) = (k2(i)-y2(i))*Z1yiRi(i);
            M(2,4) = (k2(i)-y2(i))*W1yiRi(i);
            M(3,1) = -Z1xiRi(i);
            M(3,2) = -W1xiRi(i);
            M(3,3) = Wavenumber(i)*Z1yiRi(i);
            M(3,4) = Wavenumber(i)*W1yiRi(i);
            M(3,5) = Z1ziRi(i);
            M(4,1) = Ro2*(kT2(i)/2-k2(i))*Z0xo(i)-Z1xoRo(i);
            M(4,2) = Ro2*(kT2(i)/2-k2(i))*W0xo(i)-W1xoRo(i);
            M(4,3) = Wavenumber(i)*(Z1yoRo(i)-y2(i)*Ro2*Z0yo(i));
            M(4,4) = Wavenumber(i)*(W1yoRo(i)-y2(i)*Ro2*W0yo(i));
            M(4,6) = kT2(i)/2*Ro2*Densityo*H0zo(i);
            M(5,1) = -2*Wavenumber(i)*Z1xoRo(i);
            M(5,2) = -2*Wavenumber(i)*W1xoRo(i);
            M(5,3) = (k2(i)-y2(i))*Z1yoRo(i);
            M(5,4) = (k2(i)-y2(i))*W1yoRo(i);
            M(6,1) = -Z1xoRo(i);
            M(6,2) = -W1xoRo(i);
            M(6,3) = Wavenumber(i)*Z1yoRo(i);
            M(6,4) = Wavenumber(i)*W1yoRo(i);
            M(6,6) = H1zoRo(i);
        elseif ToggleInnerFluid && ~ToggleOuterFluid
            M(1,1) = Ri2*(kT2(i)/2-k2(i))*Z0xi(i)-Z1xiRi(i);
            M(1,2) = Ri2*(kT2(i)/2-k2(i))*W0xi(i)-W1xiRi(i);
            M(1,3) = Wavenumber(i)*(Z1yiRi(i)-y2(i)*Ri2*Z0yi(i));
            M(1,4) = Wavenumber(i)*(W1yiRi(i)-y2(i)*Ri2*W0yi(i));
            M(1,5) = kT2(i)/2*Ri2*Densityi*Z0zi(i);
            M(2,1) = -2*Wavenumber(i)*Z1xiRi(i);
            M(2,2) = -2*Wavenumber(i)*W1xiRi(i);
            M(2,3) = (k2(i)-y2(i))*Z1yiRi(i);
            M(2,4) = (k2(i)-y2(i))*W1yiRi(i);
            M(3,1) = -Z1xiRi(i);
            M(3,2) = -W1xiRi(i);
            M(3,3) = Wavenumber(i)*Z1yiRi(i);
            M(3,4) = Wavenumber(i)*W1yiRi(i);
            M(3,5) = Z1ziRi(i);
            M(4,1) = Ro2*(kT2(i)/2-k2(i))*Z0xo(i)-Z1xoRo(i);
            M(4,2) = Ro2*(kT2(i)/2-k2(i))*W0xo(i)-W1xoRo(i);
            M(4,3) = Wavenumber(i)*(Z1yoRo(i)-y2(i)*Ro2*Z0yo(i));
            M(4,4) = Wavenumber(i)*(W1yoRo(i)-y2(i)*Ro2*W0yo(i));
            M(5,1) = -2*Wavenumber(i)*Z1xoRo(i);
            M(5,2) = -2*Wavenumber(i)*W1xoRo(i);
            M(5,3) = (k2(i)-y2(i))*Z1yoRo(i);
            M(5,4) = (k2(i)-y2(i))*W1yoRo(i);
        elseif ~ToggleInnerFluid && ToggleOuterFluid
            M(1,1) = Ri2*(kT2(i)/2-k2(i))*Z0xi(i)-Z1xiRi(i);
            M(1,2) = Ri2*(kT2(i)/2-k2(i))*W0xi(i)-W1xiRi(i);
            M(1,3) = Wavenumber(i)*(Z1yiRi(i)-y2(i)*Ri2*Z0yi(i));
            M(1,4) = Wavenumber(i)*(W1yiRi(i)-y2(i)*Ri2*W0yi(i));
            M(2,1) = -2*Wavenumber(i)*Z1xiRi(i);
            M(2,2) = -2*Wavenumber(i)*W1xiRi(i);
            M(2,3) = (k2(i)-y2(i))*Z1yiRi(i);
            M(2,4) = (k2(i)-y2(i))*W1yiRi(i);
            M(3,1) = Ro2*(kT2(i)/2-k2(i))*Z0xo(i)-Z1xoRo(i);
            M(3,2) = Ro2*(kT2(i)/2-k2(i))*W0xo(i)-W1xoRo(i);
            M(3,3) = Wavenumber(i)*(Z1yoRo(i)-y2(i)*Ro2*Z0yo(i));
            M(3,4) = Wavenumber(i)*(W1yoRo(i)-y2(i)*Ro2*W0yo(i));
            M(3,5) = kT2(i)/2*Ro2*Densityo*H0zo(i);
            M(4,1) = -2*Wavenumber(i)*Z1xoRo(i);
            M(4,2) = -2*Wavenumber(i)*W1xoRo(i);
            M(4,3) = (k2(i)-y2(i))*Z1yoRo(i);
            M(4,4) = (k2(i)-y2(i))*W1yoRo(i);
            M(5,1) = -Z1xoRo(i);
            M(5,2) = -W1xoRo(i);
            M(5,3) = Wavenumber(i)*Z1yoRo(i);
            M(5,4) = Wavenumber(i)*W1yoRo(i);
            M(5,5) = H1zoRo(i);
        elseif ~ToggleInnerFluid && ~ToggleOuterFluid
            M(1,1) = Ri2*(kT2(i)/2-k2(i))*Z0xi(i)-Z1xiRi(i);
            M(1,2) = Ri2*(kT2(i)/2-k2(i))*W0xi(i)-W1xiRi(i);
            M(1,3) = Wavenumber(i)*(Z1yiRi(i)-y2(i)*Ri2*Z0yi(i));
            M(1,4) = Wavenumber(i)*(W1yiRi(i)-y2(i)*Ri2*W0yi(i));
            M(2,1) = -2*Wavenumber(i)*Z1xiRi(i);
            M(2,2) = -2*Wavenumber(i)*W1xiRi(i);
            M(2,3) = (k2(i)-y2(i))*Z1yiRi(i);
            M(2,4) = (k2(i)-y2(i))*W1yiRi(i);
            M(3,1) = Ro2*(kT2(i)/2-k2(i))*Z0xo(i)-Z1xoRo(i);
            M(3,2) = Ro2*(kT2(i)/2-k2(i))*W0xo(i)-W1xoRo(i);
            M(3,3) = Wavenumber(i)*(Z1yoRo(i)-y2(i)*Ro2*Z0yo(i));
            M(3,4) = Wavenumber(i)*(W1yoRo(i)-y2(i)*Ro2*W0yo(i));
            M(4,1) = -2*Wavenumber(i)*Z1xoRo(i);
            M(4,2) = -2*Wavenumber(i)*W1xoRo(i);
            M(4,3) = (k2(i)-y2(i))*Z1yoRo(i);
            M(4,4) = (k2(i)-y2(i))*W1yoRo(i);  
        end
        Yc(i) = det(M);
    end
    Y = abs(Yc);
end
function [Y,Yc] = Computer_Rod_F(n,Wavenumber,k2,kT2,R2,y2,xR,yR)
    n2 = n^2;
    Jnx = besselj(n,xR);
    Jny = besselj(n,yR);
    dJnxR = n*Jnx-xR.*besselj(n+1,xR);
    dJnyR = n*Jny-yR.*besselj(n+1,yR);
    a1 = kT2/2*R2;
    Yc = NaN(size(Wavenumber));
    for i = 1:length(Wavenumber)
        M(1,1) = dJnxR(i)+(a1(i)-k2(i)*R2-n2)*Jnx(i);
        M(1,2) = -Wavenumber(i)*(dJnyR(i)+(y2(i)*R2-n2)*Jny(i));
        M(1,3) = n*(dJnyR(i)-Jny(i));
        M(2,1) = 2*n*(dJnxR(i)-Jnx(i));
        M(2,2) = 2*Wavenumber(i)*n*(Jny(i)-dJnyR(i));
        M(2,3) = 2*dJnyR(i)+(y2(i)*R2-2*n2)*Jny(i);
        M(3,1) = 2*Wavenumber(i)*dJnxR(i);
        M(3,2) = (y2(i)-k2(i))*dJnyR(i);
        M(3,3) = -Wavenumber(i)*n*Jny(i);
        Yc(i) = det(M);
    end
    Y = abs(Yc);
end