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
function [HCLamb,HCShear] = FrequencySweeper_Isotropic_Circumferential(SweepRange,Material,PhaseVelocity,Ro,Ri,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Bisections = 17; % ceil(log2(1e5))

%#ok<*AGROW>
HCLamb=[];HCShear=[];XRough=[];XSHRough=[];
Xi = Material.LongitudinalVelocity/Material.TransverseVelocity;
Xi2 = Xi^2;
b = Ri/Ro;
b2 = b^2;
bXi2 = b2/Xi2;
AngularFrequency = 2*pi*SweepRange*1e3;
kRo = AngularFrequency/PhaseVelocity*Ro;
kT = AngularFrequency/Material.TransverseVelocity;
kTRo = kT*Ro;
kTRi = kT*Ri;
Y = Computer(kRo,kTRo,b,b2,Xi,Xi2,bXi2);
YSH = abs((besselj(kRo-1,kTRi)-besselj(kRo+1,kTRi)).*(bessely(kRo-1,kTRo)-bessely(kRo+1,kTRo))-(besselj(kRo-1,kTRo)-besselj(kRo+1,kTRo)).*(bessely(kRo-1,kTRi)-bessely(kRo+1,kTRi)));
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
HCLamb = Converger('Lamb',XRough,Bisections,PhaseVelocity,SweepRange,Material,Ro,Ri,b,b2,Xi,Xi2,bXi2);
HCShear = Converger('SH',XSHRough,Bisections,PhaseVelocity,SweepRange,Material,Ro,Ri,b,b2,Xi,Xi2,bXi2);
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
function X = Converger(ModeFamily,XRough,Bisections,PhaseVelocity,SweepRange,Material,Ro,Ri,b,b2,Xi,Xi2,bXi2)
    X = [];
    for j = 1:length(XRough)
        Frequency = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            Frequency = Frequency(1):(Frequency(end)-Frequency(1))/4:Frequency(end);
            AngularFrequency = 2*pi*Frequency*1e3;
            kRo = AngularFrequency/PhaseVelocity*Ro;
            kT = AngularFrequency/Material.TransverseVelocity;
            kTRo = kT*Ro;
            if  strcmp(ModeFamily,'Lamb')
                [Y,Yc] = Computer(kRo,kTRo,b,b2,Xi,Xi2,bXi2);
            elseif strcmp(ModeFamily,'SH')
                kTRi = kT*Ri;
                Yc = (besselj(kRo-1,kTRi)-besselj(kRo+1,kTRi)).*(bessely(kRo-1,kTRo)-bessely(kRo+1,kTRo))-(besselj(kRo-1,kTRo)-besselj(kRo+1,kTRo)).*(bessely(kRo-1,kTRi)-bessely(kRo+1,kTRi));
                Y = abs(Yc);
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
function [Y,Yc] = Computer(kRo,e,b,b2,Xi,Xi2,bXi2)
    be = b*e;
    eXi = e/Xi;
    beXi = b*eXi;
    J_2eX = besselj(kRo-2,eXi);
    J2eX = besselj(kRo+2,eXi);
    JeX = besselj(kRo,eXi);
    J_2e = besselj(kRo-2,e);
    J2e = besselj(kRo+2,e);
    Y_2eX = bessely(kRo-2,eXi);
    Y2eX = bessely(kRo+2,eXi);
    YeX = bessely(kRo,eXi);
    Y_2e = bessely(kRo-2,e);
    Y2e = bessely(kRo+2,e);
    J_2beX = besselj(kRo-2,beXi);
    J2beX = besselj(kRo+2,beXi);
    JbeX = besselj(kRo,beXi);
    J_2be = besselj(kRo-2,be);
    J2be = besselj(kRo+2,be);
    Y_2beX = bessely(kRo-2,beXi);
    Y2beX = bessely(kRo+2,beXi);
    YbeX = bessely(kRo,beXi);
    Y_2be = bessely(kRo-2,be);
    Y2be = bessely(kRo+2,be);
    Yc = NaN(size(kRo));
    for i = 1:length(kRo)
        M(1,1) = (J_2eX(i)+J2eX(i)-2*(Xi2-1)*JeX(i))/Xi2;
        M(1,2) = 1i*(J_2e(i)-J2e(i));
        M(1,3) = (Y_2eX(i)+Y2eX(i)-2*(Xi2-1)*YeX(i))/Xi2;
        M(1,4) = 1i*(Y_2e(i)-Y2e(i));
        M(2,1) = 1i*(J_2eX(i)-J2eX(i))/Xi2;
        M(2,2) = -(J_2e(i)+J2e(i));
        M(2,3) = 1i*(Y_2eX(i)-Y2eX(i))/Xi2;
        M(2,4) = -(Y_2e(i)+Y2e(i));
        M(3,1) = (J_2beX(i)+J2beX(i)-2*(Xi2-1)*JbeX(i))*bXi2;
        M(3,2) = 1i*(J_2be(i)-J2be(i))*b2;
        M(3,3) = (Y_2beX(i)+Y2beX(i)-2*(Xi2-1)*YbeX(i))*bXi2;
        M(3,4) = 1i*(Y_2be(i)-Y2be(i))*b2;
        M(4,1) = 1i*(J_2beX(i)-J2beX(i))*bXi2;
        M(4,2) = -(J_2be(i)+J2be(i))*b2;
        M(4,3) = 1i*(Y_2beX(i)-Y2beX(i))*bXi2;
        M(4,4) = -(Y_2be(i)+Y2be(i))*b2;
        Yc(i) = det(M);
    end
    Y = abs(Yc);
end