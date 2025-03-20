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
function [HCLamb,HCShear] = FrequencySweeper_Isotropic_Circumferential(SweepRange,Material,PhaseVelocityLimit,Ro,Ri,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Resolution = 1e-5; % (kHz)

%#ok<*AGROW>
HCLamb = [];
HCShear = [];
XRoughLamb = [];
XRoughShear = [];
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
else
    Bisections = 1;
end
Xi = Material.LongitudinalVelocity/Material.TransverseVelocity;
Xi2 = Xi^2;
g = Ri/Ro;
g2 = g^2;
gXi2 = g2/Xi2;
AngularFrequency = 2*pi*SweepRange*1e3;
kRo = AngularFrequency/PhaseVelocityLimit*Ro;
kTRo = AngularFrequency/Material.TransverseVelocity*Ro;
kTRi = AngularFrequency/Material.TransverseVelocity*Ri;
YShear = abs((besselj(kRo-1,kTRi)-besselj(kRo+1,kTRi)).*(bessely(kRo-1,kTRo)-bessely(kRo+1,kTRo))-(besselj(kRo-1,kTRo)-besselj(kRo+1,kTRo)).*(bessely(kRo-1,kTRi)-bessely(kRo+1,kTRi)));
e = AngularFrequency/Material.TransverseVelocity*Ro;
ge = g*e;
eXi = e/Xi;
geXi = g*eXi;
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
J_2geX = besselj(kRo-2,geXi);
J2geX = besselj(kRo+2,geXi);
JgeX = besselj(kRo,geXi);
J_2ge = besselj(kRo-2,ge);
J2ge = besselj(kRo+2,ge);
Y_2geX = bessely(kRo-2,geXi);
Y2geX = bessely(kRo+2,geXi);
YgeX = bessely(kRo,geXi);
Y_2ge = bessely(kRo-2,ge);
Y2ge = bessely(kRo+2,ge);
for i = 1:length(SweepRange)
    M(1,1) = (J_2eX(i)+J2eX(i)-2*(Xi2-1)*JeX(i))/Xi2;
    M(1,2) = 1i*(J_2e(i)-J2e(i));
    M(1,3) = (Y_2eX(i)+Y2eX(i)-2*(Xi2-1)*YeX(i))/Xi2;
    M(1,4) = 1i*(Y_2e(i)-Y2e(i));
    M(2,1) = 1i*(J_2eX(i)-J2eX(i))/Xi2;
    M(2,2) = -(J_2e(i)+J2e(i));
    M(2,3) = 1i*(Y_2eX(i)-Y2eX(i))/Xi2;
    M(2,4) = -(Y_2e(i)+Y2e(i));
    M(3,1) = (J_2geX(i)+J2geX(i)-2*(Xi2-1)*JgeX(i))*gXi2;
    M(3,2) = 1i*(J_2ge(i)-J2ge(i))*g2;
    M(3,3) = (Y_2geX(i)+Y2geX(i)-2*(Xi2-1)*YgeX(i))*gXi2;
    M(3,4) = 1i*(Y_2ge(i)-Y2ge(i))*g2;
    M(4,1) = 1i*(J_2geX(i)-J2geX(i))*gXi2;
    M(4,2) = -(J_2ge(i)+J2ge(i))*g2;
    M(4,3) = 1i*(Y_2geX(i)-Y2geX(i))*gXi2;
    M(4,4) = -(Y_2ge(i)+Y2ge(i))*g2;
    YLamb(i) = abs(det(M));
    if  i > 2 && YLamb(i-1) < YLamb(i-2) && YLamb(i-1) < YLamb(i)
        XRoughLamb(end+1) = SweepRange(i-1);
    end
end
for i = 2:length(SweepRange)-1
    if  YShear(i) < YShear(i-1) && YShear(i) < YShear(i+1)
        XRoughShear(end+1) = SweepRange(i);
    end
end
if  isempty(XRoughLamb) && isempty(XRoughShear)
    String = 'No higher order modes found!';
    OutputWindow1aUI1.String = String;
    OutputWindow1bUI1.String = '';
    OutputWindow2aUI1.String = '';
    OutputWindow2bUI1.String = '';
    disp(String)
    return
end
for j = 1:length(XRoughLamb)
    Frequency = [XRoughLamb(j)-(SweepRange(2)-SweepRange(1)) XRoughLamb(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            kRo = AngularFrequency(i)/PhaseVelocityLimit*Ro;
            e = AngularFrequency(i)/Material.TransverseVelocity*Ro;
            ge = g*e;
            eXi = e/Xi;
            geXi = g*eXi;
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
            J_2geX = besselj(kRo-2,geXi);
            J2geX = besselj(kRo+2,geXi);
            JgeX = besselj(kRo,geXi);
            J_2ge = besselj(kRo-2,ge);
            J2ge = besselj(kRo+2,ge);
            Y_2geX = bessely(kRo-2,geXi);
            Y2geX = bessely(kRo+2,geXi);
            YgeX = bessely(kRo,geXi);
            Y_2ge = bessely(kRo-2,ge);
            Y2ge = bessely(kRo+2,ge);
            M(1,1) = (J_2eX+J2eX-2*(Xi2-1)*JeX)/Xi2;
            M(1,2) = 1i*(J_2e-J2e);
            M(1,3) = (Y_2eX+Y2eX-2*(Xi2-1)*YeX)/Xi2;
            M(1,4) = 1i*(Y_2e-Y2e);
            M(2,1) = 1i*(J_2eX-J2eX)/Xi2;
            M(2,2) = -(J_2e+J2e);
            M(2,3) = 1i*(Y_2eX-Y2eX)/Xi2;
            M(2,4) = -(Y_2e+Y2e);
            M(3,1) = (J_2geX+J2geX-2*(Xi2-1)*JgeX)*gXi2;
            M(3,2) = 1i*(J_2ge-J2ge)*g2;
            M(3,3) = (Y_2geX+Y2geX-2*(Xi2-1)*YgeX)*gXi2;
            M(3,4) = 1i*(Y_2ge-Y2ge)*g2;
            M(4,1) = 1i*(J_2geX-J2geX)*gXi2;
            M(4,2) = -(J_2ge+J2ge)*g2;
            M(4,3) = 1i*(Y_2geX-Y2geX)*gXi2;
            M(4,4) = -(Y_2ge+Y2ge)*g2;
            Y(i) = abs(det(M));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HCLamb(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
end
for j = 1:length(XRoughShear)
    Frequency = [XRoughShear(j)-(SweepRange(2)-SweepRange(1)) XRoughShear(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            kRo = AngularFrequency(i)/PhaseVelocityLimit*Ro;
            kTRo = AngularFrequency(i)/Material.TransverseVelocity*Ro;
            kTRi = AngularFrequency(i)/Material.TransverseVelocity*Ri;
            Y(i) = abs((besselj(kRo-1,kTRi)-besselj(kRo+1,kTRi))*(bessely(kRo-1,kTRo)-bessely(kRo+1,kTRo))-(besselj(kRo-1,kTRo)-besselj(kRo+1,kTRo))*(bessely(kRo-1,kTRi)-bessely(kRo+1,kTRi)));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HCShear(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
end
String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
if  any(HCLamb)
    for i = 1:length(HCLamb)
        if  i < 10
            if  HCLamb(i) < 1e2
                String = append(String,newline,'C',num2str(i),'       ',num2str(HCLamb(i),'%.3f'));
            elseif HCLamb(i) >= 1e2 && HCLamb(i) < 1e3
                String = append(String,newline,'C',num2str(i),'      ',num2str(HCLamb(i),'%.3f'));
            elseif HCLamb(i) >= 1e3 && HCLamb(i) < 1e4
                String = append(String,newline,'C',num2str(i),'     ',num2str(HCLamb(i),'%.3f'));
            elseif HCLamb(i) >= 1e4
                String = append(String,newline,'C',num2str(i),'    ',num2str(HCLamb(i),'%.3f'));
            end
        else
            if  HCLamb(i) < 1e2
                String = append(String,newline,'C',num2str(i),'      ',num2str(HCLamb(i),'%.3f'));
            elseif HCLamb(i) >= 1e2 && HCLamb(i) < 1e3
                String = append(String,newline,'C',num2str(i),'     ',num2str(HCLamb(i),'%.3f'));
            elseif HCLamb(i) >= 1e3 && HCLamb(i) < 1e4
                String = append(String,newline,'C',num2str(i),'    ',num2str(HCLamb(i),'%.3f'));
            elseif HCLamb(i) >= 1e4
                String = append(String,newline,'C',num2str(i),'   ',num2str(HCLamb(i),'%.3f'));
            end
        end
    end
end
OutputWindow1aUI1.String = String;
disp(String)
String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
if  any(HCShear)
    for i = 1:length(HCShear)
        if  i < 10
            if  HCShear(i) < 1e2
                String = append(String,newline,'CSH',num2str(i),'     ',num2str(HCShear(i),'%.3f'));
            elseif HCShear(i) >= 1e2 && HCShear(i) < 1e3
                String = append(String,newline,'CSH',num2str(i),'    ',num2str(HCShear(i),'%.3f'));
            elseif HCShear(i) >= 1e3 && HCShear(i) < 1e4
                String = append(String,newline,'CSH',num2str(i),'   ',num2str(HCShear(i),'%.3f'));
            elseif HCShear(i) >= 1e4
                String = append(String,newline,'CSH',num2str(i),'  ',num2str(HCShear(i),'%.3f'));
            end
        else
            if  HCShear(i) < 1e2
                String = append(String,newline,'CSH',num2str(i),'    ',num2str(HCShear(i),'%.3f'));
            elseif HCShear(i) >= 1e2 && HCShear(i) < 1e3
                String = append(String,newline,'CSH',num2str(i),'   ',num2str(HCShear(i),'%.3f'));
            elseif HCShear(i) >= 1e3 && HCShear(i) < 1e4
                String = append(String,newline,'CSH',num2str(i),'  ',num2str(HCShear(i),'%.3f'));
            elseif HCShear(i) >= 1e4
                String = append(String,newline,'CSH',num2str(i),' ',num2str(HCShear(i),'%.3f'));
            end
        end
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