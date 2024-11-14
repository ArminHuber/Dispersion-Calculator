% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2024 DLR
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
function [HL,HF,HT] = FrequencySweeper_Isotropic_Pipe(SweepRange,Material,PhaseVelocityLimit,Ro,Ri,FlexuralModeOrders,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Resolution = 1e-5; % (kHz)

%#ok<*AGROW>
HL = [];
HT = [];
HF = [];
XRoughL = [];
XRoughT = [];
XRoughF = zeros(FlexuralModeOrders,1);
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
else
    Bisections = 1;
end
Ro2 = Ro^2;
Ri2 = Ri^2;
AngularFrequency = 2*pi*SweepRange*1e3;
k = AngularFrequency/PhaseVelocityLimit;
k2 = k.^2;
kT2 = (AngularFrequency/Material.TransverseVelocity).^2;
x = sqrt((AngularFrequency/Material.LongitudinalVelocity).^2-k2);
y = sqrt(kT2-k2);
y2 = kT2-k2;
xRo = x*Ro;
yRo = y*Ro;
xRi = x*Ri;
yRi = y*Ri;
YT = abs(besselj(2,yRi).*bessely(2,yRo)-besselj(2,yRo).*bessely(2,yRi));
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
for i = 1:length(SweepRange)
    M(1,1) = Ri2*(.5*kT2(i)-k2(i))*Z0xi(i)-Z1xiRi(i);
    M(1,2) = Ri2*(.5*kT2(i)-k2(i))*W0xi(i)-W1xiRi(i);
    M(1,3) = k(i)*(Z1yiRi(i)-y2(i)*Ri2*Z0yi(i));
    M(1,4) = k(i)*(W1yiRi(i)-y2(i)*Ri2*W0yi(i));
    M(2,1) = -2*k(i)*Z1xiRi(i);
    M(2,2) = -2*k(i)*W1xiRi(i);
    M(2,3) = (k2(i)-y2(i))*Z1yiRi(i);
    M(2,4) = (k2(i)-y2(i))*W1yiRi(i);
    M(3,1) = Ro2*(.5*kT2(i)-k2(i))*Z0xo(i)-Z1xoRo(i);
    M(3,2) = Ro2*(.5*kT2(i)-k2(i))*W0xo(i)-W1xoRo(i);
    M(3,3) = k(i)*(Z1yoRo(i)-y2(i)*Ro2*Z0yo(i));
    M(3,4) = k(i)*(W1yoRo(i)-y2(i)*Ro2*W0yo(i));
    M(4,1) = -2*k(i)*Z1xoRo(i);
    M(4,2) = -2*k(i)*W1xoRo(i);
    M(4,3) = (k2(i)-y2(i))*Z1yoRo(i);
    M(4,4) = (k2(i)-y2(i))*W1yoRo(i);
    YL(i) = abs(det(M));
    if  i > 2 && YL(i-1) < YL(i-2) && YL(i-1) < YL(i)
        XRoughL(end+1) = SweepRange(i-1);
    end
end
for n = 1:FlexuralModeOrders
    n2 = n^2;
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
    for i = 1:length(SweepRange)
        M(1,1) = dZnxiRi(i)+(.5*kT2(i)*Ri2-k2(i)*Ri2-n2)*Znxi(i);
        M(1,2) = dWnxiRi(i)+(.5*kT2(i)*Ri2-k2(i)*Ri2-n2)*Wnxi(i);
        M(1,3) = -k(i)*(dZnyiRi(i)+(y2(i)*Ri2-n2)*Znyi(i));
        M(1,4) = -k(i)*(dWnyiRi(i)+(y2(i)*Ri2-n2)*Wnyi(i));
        M(1,5) = n*(dZnyiRi(i)-Znyi(i));
        M(1,6) = n*(dWnyiRi(i)-Wnyi(i));
        M(2,1) = 2*n*(dZnxiRi(i)-Znxi(i));
        M(2,2) = 2*n*(dWnxiRi(i)-Wnxi(i));
        M(2,3) = 2*k(i)*n*(Znyi(i)-dZnyiRi(i));
        M(2,4) = 2*k(i)*n*(Wnyi(i)-dWnyiRi(i));
        M(2,5) = 2*dZnyiRi(i)+(y2(i)*Ri2-2*n2)*Znyi(i);
        M(2,6) = 2*dWnyiRi(i)+(y2(i)*Ri2-2*n2)*Wnyi(i);
        M(3,1) = 2*k(i)*dZnxiRi(i);
        M(3,2) = 2*k(i)*dWnxiRi(i);
        M(3,3) = (y2(i)-k2(i))*dZnyiRi(i);
        M(3,4) = (y2(i)-k2(i))*dWnyiRi(i);
        M(3,5) = -k(i)*n*Znyi(i);
        M(3,6) = -k(i)*n*Wnyi(i);
        M(4,1) = dZnxoRo(i)+(.5*kT2(i)*Ro2-k2(i)*Ro2-n2)*Znxo(i);
        M(4,2) = dWnxoRo(i)+(.5*kT2(i)*Ro2-k2(i)*Ro2-n2)*Wnxo(i);
        M(4,3) = -k(i)*(dZnyoRo(i)+(y2(i)*Ro2-n2)*Znyo(i));
        M(4,4) = -k(i)*(dWnyoRo(i)+(y2(i)*Ro2-n2)*Wnyo(i));
        M(4,5) = n*(dZnyoRo(i)-Znyo(i));
        M(4,6) = n*(dWnyoRo(i)-Wnyo(i));
        M(5,1) = 2*n*(dZnxoRo(i)-Znxo(i));
        M(5,2) = 2*n*(dWnxoRo(i)-Wnxo(i));
        M(5,3) = 2*k(i)*n*(Znyo(i)-dZnyoRo(i));
        M(5,4) = 2*k(i)*n*(Wnyo(i)-dWnyoRo(i));
        M(5,5) = 2*dZnyoRo(i)+(y2(i)*Ro2-2*n2)*Znyo(i);
        M(5,6) = 2*dWnyoRo(i)+(y2(i)*Ro2-2*n2)*Wnyo(i);
        M(6,1) = 2*k(i)*dZnxoRo(i);
        M(6,2) = 2*k(i)*dWnxoRo(i);
        M(6,3) = (y2(i)-k2(i))*dZnyoRo(i);
        M(6,4) = (y2(i)-k2(i))*dWnyoRo(i);
        M(6,5) = -k(i)*n*Znyo(i);
        M(6,6) = -k(i)*n*Wnyo(i);
        YF(n,i) = abs(det(M));
        if  i > 2 && YF(n,i-1) < YF(n,i-2) && YF(n,i-1) < YF(n,i)
            XRoughF(n,1+numel(find(XRoughF(n,:) > 0))) = SweepRange(i-1);
        end
    end
    if  all(XRoughF(n,:) == 0)
        XRoughF(n:end,:) = [];
        break
    end
end
for i = 2:length(SweepRange)-1
    if  YT(i) < YT(i-1) && YT(i) < YT(i+1)
        XRoughT(end+1) = SweepRange(i);
    end
end
if  isempty(XRoughL) && all(all(XRoughF == 0)) && isempty(XRoughT)
    String = 'No higher order modes found!';
    OutputWindow1aUI1.String = String;
    OutputWindow1bUI1.String = '';
    OutputWindow2aUI1.String = '';
    OutputWindow2bUI1.String = '';
    disp(String)
    return
end
M = 0;
for j = 1:length(XRoughL)
    Frequency = [XRoughL(j)-(SweepRange(2)-SweepRange(1)) XRoughL(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            k = AngularFrequency(i)/PhaseVelocityLimit;
            k2 = k^2;
            kT2 = (AngularFrequency(i)/Material.TransverseVelocity)^2;
            x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
            y = sqrt(kT2-k2);
            y2 = kT2-k2;
            xRo = x*Ro;
            yRo = y*Ro;
            xRi = x*Ri;
            yRi = y*Ri;
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
            M(1,1) = Ri2*(.5*kT2-k2)*Z0xi-Z1xiRi;
            M(1,2) = Ri2*(.5*kT2-k2)*W0xi-W1xiRi;
            M(1,3) = k*(Z1yiRi-y2*Ri2*Z0yi);
            M(1,4) = k*(W1yiRi-y2*Ri2*W0yi);
            M(2,1) = -2*k*Z1xiRi;
            M(2,2) = -2*k*W1xiRi;
            M(2,3) = (k2-y2)*Z1yiRi;
            M(2,4) = (k2-y2)*W1yiRi;
            M(3,1) = Ro2*(.5*kT2-k2)*Z0xo-Z1xoRo;
            M(3,2) = Ro2*(.5*kT2-k2)*W0xo-W1xoRo;
            M(3,3) = k*(Z1yoRo-y2*Ro2*Z0yo);
            M(3,4) = k*(W1yoRo-y2*Ro2*W0yo);
            M(4,1) = -2*k*Z1xoRo;
            M(4,2) = -2*k*W1xoRo;
            M(4,3) = (k2-y2)*Z1yoRo;
            M(4,4) = (k2-y2)*W1yoRo;
            Y(i) = abs(det(M));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HL(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
end
HF = zeros(size(XRoughF));
for n = 1:height(XRoughF)
    n2 = n^2;
    for j = 1:numel(find(XRoughF(n,:) > 0))
        Frequency = [XRoughF(n,j)-(SweepRange(2)-SweepRange(1)) XRoughF(n,j)+(SweepRange(2)-SweepRange(1))];
        for o = 1:Bisections
            Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
            AngularFrequency = 2*pi*Frequency*1e3;
            for i = 1:length(Frequency)
                k = AngularFrequency(i)/PhaseVelocityLimit;
                k2 = k^2;
                kT2 = (AngularFrequency(i)/Material.TransverseVelocity)^2;
                x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
                y = sqrt(kT2-k2);
                y2 = kT2-k2;
                xRo = x*Ro;
                yRo = y*Ro;
                xRi = x*Ri;
                yRi = y*Ri;
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
                M(1,1) = dZnxiRi+(.5*kT2*Ri2-k2*Ri2-n2)*Znxi;
                M(1,2) = dWnxiRi+(.5*kT2*Ri2-k2*Ri2-n2)*Wnxi;
                M(1,3) = -k*(dZnyiRi+(y2*Ri2-n2)*Znyi);
                M(1,4) = -k*(dWnyiRi+(y2*Ri2-n2)*Wnyi);
                M(1,5) = n*(dZnyiRi-Znyi);
                M(1,6) = n*(dWnyiRi-Wnyi);
                M(2,1) = 2*n*(dZnxiRi-Znxi);
                M(2,2) = 2*n*(dWnxiRi-Wnxi);
                M(2,3) = 2*k*n*(Znyi-dZnyiRi);
                M(2,4) = 2*k*n*(Wnyi-dWnyiRi);
                M(2,5) = 2*dZnyiRi+(y2*Ri2-2*n2)*Znyi;
                M(2,6) = 2*dWnyiRi+(y2*Ri2-2*n2)*Wnyi;
                M(3,1) = 2*k*dZnxiRi;
                M(3,2) = 2*k*dWnxiRi;
                M(3,3) = (y2-k2)*dZnyiRi;
                M(3,4) = (y2-k2)*dWnyiRi;
                M(3,5) = -k*n*Znyi;
                M(3,6) = -k*n*Wnyi;
                M(4,1) = dZnxoRo+(.5*kT2*Ro2-k2*Ro2-n2)*Znxo;
                M(4,2) = dWnxoRo+(.5*kT2*Ro2-k2*Ro2-n2)*Wnxo;
                M(4,3) = -k*(dZnyoRo+(y2*Ro2-n2)*Znyo);
                M(4,4) = -k*(dWnyoRo+(y2*Ro2-n2)*Wnyo);
                M(4,5) = n*(dZnyoRo-Znyo);
                M(4,6) = n*(dWnyoRo-Wnyo);
                M(5,1) = 2*n*(dZnxoRo-Znxo);
                M(5,2) = 2*n*(dWnxoRo-Wnxo);
                M(5,3) = 2*k*n*(Znyo-dZnyoRo);
                M(5,4) = 2*k*n*(Wnyo-dWnyoRo);
                M(5,5) = 2*dZnyoRo+(y2*Ro2-2*n2)*Znyo;
                M(5,6) = 2*dWnyoRo+(y2*Ro2-2*n2)*Wnyo;
                M(6,1) = 2*k*dZnxoRo;
                M(6,2) = 2*k*dWnxoRo;
                M(6,3) = (y2-k2)*dZnyoRo;
                M(6,4) = (y2-k2)*dWnyoRo;
                M(6,5) = -k*n*Znyo;
                M(6,6) = -k*n*Wnyo;
                Y(i) = abs(det(M));
                if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                    if  o == Bisections
                        HF(n,j) = Frequency(i-1);
                    end
                    Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                    break
                end
            end
        end
    end
end
for j = 1:length(XRoughT)
    Frequency = [XRoughT(j)-(SweepRange(2)-SweepRange(1)) XRoughT(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            y = sqrt(1/Material.TransverseVelocity^2-1/PhaseVelocityLimit^2)*AngularFrequency(i);
            yRo = y*Ro;
            yRi = y*Ri;
            Y(i) = abs(besselj(2,yRi)*bessely(2,yRo)-besselj(2,yRo)*bessely(2,yRi));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HT(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
end
String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode       Frq.(kHz)'];
if  any(HF)
    for i = 1:numel(find(HF(1,:) > 0))
        if  i < 9
            if  HF(1,i) < 1e2
                String = append(String,newline,'F(1,',num2str(i+1),')       ',num2str(HF(1,i),'%.3f'));
            elseif HF(1,i) >= 1e2 && HF(1,i) < 1e3
                String = append(String,newline,'F(1,',num2str(i+1),')      ',num2str(HF(1,i),'%.3f'));
            elseif HF(1,i) >= 1e3 && HF(1,i) < 1e4
                String = append(String,newline,'F(1,',num2str(i+1),')     ',num2str(HF(1,i),'%.3f'));
            elseif HF(1,i) >= 1e4
                String = append(String,newline,'F(1,',num2str(i+1),')    ',num2str(HF(1,i),'%.3f'));
            end
        else
            if  HF(1,i) < 1e2
                String = append(String,newline,'F(1,',num2str(i+1),')      ',num2str(HF(1,i),'%.3f'));
            elseif HF(1,i) >= 1e2 && HF(1,i) < 1e3
                String = append(String,newline,'F(1,',num2str(i+1),')     ',num2str(HF(1,i),'%.3f'));
            elseif HF(1,i) >= 1e3 && HF(1,i) < 1e4
                String = append(String,newline,'F(1,',num2str(i+1),')    ',num2str(HF(1,i),'%.3f'));
            elseif HF(1,i) >= 1e4
                String = append(String,newline,'F(1,',num2str(i+1),')   ',num2str(HF(1,i),'%.3f'));
            end
        end
    end    
    for n = 2:height(HF)
        for i = 1:numel(find(HF(n,:) > 0))
            if  i < 10
                if  HF(n,i) < 1e2
                    String = append(String,newline,'F(',num2str(n),',',num2str(i),')       ',num2str(HF(n,i),'%.3f'));
                elseif HF(n,i) >= 1e2 && HF(n,i) < 1e3
                    String = append(String,newline,'F(',num2str(n),',',num2str(i),')      ',num2str(HF(n,i),'%.3f'));
                elseif HF(n,i) >= 1e3 && HF(n,i) < 1e4
                    String = append(String,newline,'F(',num2str(n),',',num2str(i),')     ',num2str(HF(n,i),'%.3f'));
                elseif HF(n,i) >= 1e4
                    String = append(String,newline,'F(',num2str(n),',',num2str(i),')    ',num2str(HF(n,i),'%.3f'));
                end
            else
                if  HF(n,i) < 1e2
                    String = append(String,newline,'F(',num2str(n),',',num2str(i),')      ',num2str(HF(n,i),'%.3f'));
                elseif HF(n,i) >= 1e2 && HF(n,i) < 1e3
                    String = append(String,newline,'F(',num2str(n),',',num2str(i),')     ',num2str(HF(n,i),'%.3f'));
                elseif HF(n,i) >= 1e3 && HF(n,i) < 1e4
                    String = append(String,newline,'F(',num2str(n),',',num2str(i),')    ',num2str(HF(n,i),'%.3f'));
                elseif HF(n,i) >= 1e4
                    String = append(String,newline,'F(',num2str(n),',',num2str(i),')   ',num2str(HF(n,i),'%.3f'));
                end
            end
        end
    end
end
OutputWindow1aUI1.String = String;
disp(String)
String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode       Frq.(kHz)'];
if  any(HL)
    for i = 1:length(HL)
        if  i < 9
            if  HL(i) < 1e2
                String = append(String,newline,'L(0,',num2str(i+1),')       ',num2str(HL(i),'%.3f'));
            elseif HL(i) >= 1e2 && HL(i) < 1e3
                String = append(String,newline,'L(0,',num2str(i+1),')      ',num2str(HL(i),'%.3f'));
            elseif HL(i) >= 1e3 && HL(i) < 1e4
                String = append(String,newline,'L(0,',num2str(i+1),')     ',num2str(HL(i),'%.3f'));
            elseif HL(i) >= 1e4
                String = append(String,newline,'L(0,',num2str(i+1),')    ',num2str(HL(i),'%.3f'));
            end
        else
            if  HL(i) < 1e2
                String = append(String,newline,'L(0,',num2str(i+1),')      ',num2str(HL(i),'%.3f'));
            elseif HL(i) >= 1e2 && HL(i) < 1e3
                String = append(String,newline,'L(0,',num2str(i+1),')     ',num2str(HL(i),'%.3f'));
            elseif HL(i) >= 1e3 && HL(i) < 1e4
                String = append(String,newline,'L(0,',num2str(i+1),')    ',num2str(HL(i),'%.3f'));
            elseif HL(i) >= 1e4
                String = append(String,newline,'L(0,',num2str(i+1),')   ',num2str(HL(i),'%.3f'));
            end
        end
    end
end
String = append(String,newline,' ');
if  any(HT)
    for i = 1:length(HT)
        if  i < 9
            if  HT(i) < 1e2
                String = append(String,newline,'T(0,',num2str(i+1),')       ',num2str(HT(i),'%.3f'));
            elseif HT(i) >= 1e2 && HT(i) < 1e3
                String = append(String,newline,'T(0,',num2str(i+1),')      ',num2str(HT(i),'%.3f'));
            elseif HT(i) >= 1e3 && HT(i) < 1e4
                String = append(String,newline,'T(0,',num2str(i+1),')     ',num2str(HT(i),'%.3f'));
            elseif HT(i) >= 1e4
                String = append(String,newline,'T(0,',num2str(i+1),')    ',num2str(HT(i),'%.3f'));
            end
        else
            if  HT(i) < 1e2
                String = append(String,newline,'T(0,',num2str(i+1),')      ',num2str(HT(i),'%.3f'));
            elseif HT(i) >= 1e2 && HT(i) < 1e3
                String = append(String,newline,'T(0,',num2str(i+1),')     ',num2str(HT(i),'%.3f'));
            elseif HT(i) >= 1e3 && HT(i) < 1e4
                String = append(String,newline,'T(0,',num2str(i+1),')    ',num2str(HT(i),'%.3f'));
            elseif HT(i) >= 1e4
                String = append(String,newline,'T(0,',num2str(i+1),')   ',num2str(HT(i),'%.3f'));
            end
        end
    end
end
OutputWindow1bUI1.String = String;
disp(extractAfter(String,')'))
disp(' ')
String = '';
if  any(HF)
    String = append(String,'F: ',num2str(numel(find(HF > 0))));
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