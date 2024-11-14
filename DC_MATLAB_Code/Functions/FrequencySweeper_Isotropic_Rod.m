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
function [HL,HF,HT] = FrequencySweeper_Isotropic_Rod(SweepRange,Material,PhaseVelocityLimit,R,FlexuralModeOrders,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
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
AngularFrequency = 2*pi*SweepRange*1e3;
k = AngularFrequency/PhaseVelocityLimit;
k2 = k.^2;
kT2 = (AngularFrequency/Material.TransverseVelocity).^2;
R2 = R^2;
kR2 = k2*R2;
kR4 = kR2.^2;
x = sqrt((AngularFrequency/Material.LongitudinalVelocity).^2-k2);
y = sqrt(kT2-k2);
y2 = y.^2;
xR = x*R;
yR = y*R;
yR2 = yR.^2;
yR4 = yR2.^2;
J0x = besselj(0,xR);
J0y = besselj(0,yR);
J1x = besselj(1,xR);
J1y = besselj(1,yR);
Zx = xR.*J0x./J1x;
Zy = yR.*J0y./J1y;
Zy2 = Zy.^2;
F1 = 2*(yR2-kR2).^2;
F2 = 2*yR2.*(yR2+5*kR2);
F3 = yR2.^3-10*yR4-2*yR4.*kR2+2*yR2.*kR2+yR2.*kR4-4*kR4;
F4 = 2*yR2.*(2*yR2.*kR2-yR2-9*kR2);
F5 = yR2.*(-yR4+8*yR2-2*yR2.*kR2+8*kR2-kR4);
YL = 2*x/R.*(y2+k2)-(y2-k2).^2.*J0x./J1x-4*k2.*x.*y.*J0y./J1y;
YF = F1+F2.*Zx./Zy+F3./Zy+F4.*Zx./Zy2+F5./Zy2;
YT = abs(Zy-2);
for n = 2:FlexuralModeOrders
    n2 = n^2;
    Jnx = besselj(n,xR);
    Jny = besselj(n,yR);
    dJnxR = n*Jnx-xR.*besselj(n+1,xR);
    dJnyR = n*Jny-yR.*besselj(n+1,yR);
    for i = 1:length(SweepRange)
        M(1,1) = dJnxR(i)+(.5*kT2(i)*R2-k2(i)*R2-n2)*Jnx(i);
        M(1,2) = -k(i)*(dJnyR(i)+(y2(i)*R2-n2)*Jny(i));
        M(1,3) = n*(dJnyR(i)-Jny(i));
        M(2,1) = 2*n*(dJnxR(i)-Jnx(i));
        M(2,2) = 2*k(i)*n*(Jny(i)-dJnyR(i));
        M(2,3) = 2*dJnyR(i)+(y2(i)*R2-2*n2)*Jny(i);
        M(3,1) = 2*k(i)*dJnxR(i);
        M(3,2) = (y2(i)-k2(i))*dJnyR(i);
        M(3,3) = -k(i)*n*Jny(i);
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
    if  abs(YL(1,i)) < abs(YL(1,i-1)) && abs(YL(1,i)) < abs(YL(1,i+1)) && sign(YL(1,i-1)) ~= sign(YL(1,i+1))
        XRoughL(end+1) = SweepRange(i);
    end
    if  abs(YF(1,i)) < abs(YF(1,i-1)) && abs(YF(1,i)) < abs(YF(1,i+1)) && sign(YF(1,i-1)) ~= sign(YF(1,i+1))
        XRoughF(1,1+numel(find(XRoughF(1,:) > 0))) = SweepRange(i);
    end
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
for j = 1:length(XRoughL)
    Frequency = [XRoughL(j)-(SweepRange(2)-SweepRange(1)) XRoughL(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            k2 = (AngularFrequency(i)/PhaseVelocityLimit)^2;
            x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
            y = sqrt((AngularFrequency(i)/Material.TransverseVelocity)^2-k2);
            y2 = y^2;
            xR = x*R;
            yR = y*R;
            Y(i) = abs(2*x/R*(y2+k2)-(y2-k2)^2*besselj(0,xR)/besselj(1,xR)-4*k2*x*y*besselj(0,yR)/besselj(1,yR));
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
                y2 = kT2-k2;
                xR = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2)*R;
                yR = sqrt(y2)*R;
                Jnx = besselj(n,xR);
                Jny = besselj(n,yR);
                dJnxR = n*Jnx-xR*besselj(n+1,xR);
                dJnyR = n*Jny-yR*besselj(n+1,yR);
                M(1,1) = dJnxR+(.5*kT2*R2-k2*R2-n2)*Jnx;
                M(1,2) = -k*(dJnyR+(y2*R2-n2)*Jny);
                M(1,3) = n*(dJnyR-Jny);
                M(2,1) = 2*n*(dJnxR-Jnx);
                M(2,2) = 2*k*n*(Jny-dJnyR);
                M(2,3) = 2*dJnyR+(y2*R2-2*n2)*Jny;
                M(3,1) = 2*k*dJnxR;
                M(3,2) = (y2-k2)*dJnyR;
                M(3,3) = -k*n*Jny;
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
            y = sqrt(1/Material.TransverseVelocity^2-1/PhaseVelocityLimit^2)*AngularFrequency(i)*R;
            Y(i) = abs(y-2*besselj(1,y)/besselj(0,y));
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