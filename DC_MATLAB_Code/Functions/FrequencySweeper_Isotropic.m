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
function [HSLamb,HALamb,HBLamb,HSShear,HAShear] = FrequencySweeper_Isotropic(SweepRange,Material,PhaseVelocityLimit,Half,Symmetric,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Resolution = 1e-5; % (kHz)

%#ok<*AGROW>
HSLamb = [];
HALamb = [];
HBLamb = [];
HSShear = [];
HAShear = [];
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
else
    Bisections = 1;
end
AngularFrequency = 2*pi*SweepRange*1e3;
k2 = (AngularFrequency/PhaseVelocityLimit).^2;
x = sqrt((AngularFrequency/Material.LongitudinalVelocity).^2-k2);
y = sqrt((AngularFrequency/Material.TransverseVelocity).^2-k2);
a1 = (y.^2-k2).^2./y;
a2 = 4*k2.*x;
a3 = tan(x*Half);
a4 = tan(y*Half);
YSLamb = abs(a1./a3+a2./a4);
YALamb = abs(a1.*a3+a2.*a4);
XRoughSLamb = [];
XRoughALamb = [];
for i = 2:length(SweepRange)-1
    if  YSLamb(i) < YSLamb(i-1) && YSLamb(i) < YSLamb(i+1)
        XRoughSLamb(end+1) = SweepRange(i);
    end
    if  YALamb(i) < YALamb(i-1) && YALamb(i) < YALamb(i+1)
        XRoughALamb(end+1)  = SweepRange(i);
    end
end
for n = 1:1e5
    X = n*PhaseVelocityLimit*Material.TransverseVelocity/sqrt(PhaseVelocityLimit^2-Material.TransverseVelocity^2)/Half/4e3;
    if  X > SweepRange(end)
        break
    end
    if  mod(n,2) == 0
        HSShear(end+1) = X;
    else
        HAShear(end+1) = X;
    end
end
if  isempty(XRoughSLamb) && isempty(XRoughALamb) && isempty(HSShear) && isempty(HAShear)
    String = 'No higher order modes found!';
    OutputWindow1aUI1.String = String;
    OutputWindow1bUI1.String = '';
    OutputWindow2aUI1.String = '';
    OutputWindow2bUI1.String = '';
    disp(String)
    return
end    
for j = 1:length(XRoughSLamb)
    Frequency = [XRoughSLamb(j)-(SweepRange(2)-SweepRange(1)) XRoughSLamb(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            k2 = (AngularFrequency(i)/PhaseVelocityLimit)^2;
            x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
            y = sqrt((AngularFrequency(i)/Material.TransverseVelocity)^2-k2);
            Y(i) = abs(tan(y*Half)/y+4*k2*x*tan(x*Half)/(y^2-k2)^2);
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HSLamb(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
end
for j = 1:length(XRoughALamb)
    Frequency = [XRoughALamb(j)-(SweepRange(2)-SweepRange(1)) XRoughALamb(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            k2 = (AngularFrequency(i)/PhaseVelocityLimit)^2;
            x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
            y = sqrt((AngularFrequency(i)/Material.TransverseVelocity)^2-k2);
            Y(i) = abs(y*tan(y*Half)+(y^2-k2)^2*tan(x*Half)/(4*k2*x));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HALamb(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
end
String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
if  Symmetric
    if  any(HALamb)
        for i = 1:length(HALamb)
            if  i < 10
                if  HALamb(i) < 1e2
                    String = append(String,newline,'A',num2str(i),'       ',num2str(HALamb(i),'%.3f'));
                elseif HALamb(i) >= 1e2 && HALamb(i) < 1e3
                    String = append(String,newline,'A',num2str(i),'      ',num2str(HALamb(i),'%.3f'));
                elseif HALamb(i) >= 1e3 && HALamb(i) < 1e4
                    String = append(String,newline,'A',num2str(i),'     ',num2str(HALamb(i),'%.3f'));
                elseif HALamb(i) >= 1e4
                    String = append(String,newline,'A',num2str(i),'    ',num2str(HALamb(i),'%.3f'));
                end
            else
                if  HALamb(i) < 1e2
                    String = append(String,newline,'A',num2str(i),'      ',num2str(HALamb(i),'%.3f'));
                elseif HALamb(i) >= 1e2 && HALamb(i) < 1e3
                    String = append(String,newline,'A',num2str(i),'     ',num2str(HALamb(i),'%.3f'));
                elseif HALamb(i) >= 1e3 && HALamb(i) < 1e4
                    String = append(String,newline,'A',num2str(i),'    ',num2str(HALamb(i),'%.3f'));
                elseif HALamb(i) >= 1e4
                    String = append(String,newline,'A',num2str(i),'   ',num2str(HALamb(i),'%.3f'));
                end
            end
        end
    end
else
    HBLamb = sort(horzcat(HALamb,HSLamb));
    if  any(HBLamb)
        for i = 1:length(HBLamb)
            if  i < 9
                if  HBLamb(i) < 1e2
                    String = append(String,newline,'B',num2str(i+1),'       ',num2str(HBLamb(i),'%.3f'));
                elseif HBLamb(i) >= 1e2 && HBLamb(i) < 1e3
                    String = append(String,newline,'B',num2str(i+1),'      ',num2str(HBLamb(i),'%.3f'));
                elseif HBLamb(i) >= 1e3 && HBLamb(i) < 1e4
                    String = append(String,newline,'B',num2str(i+1),'     ',num2str(HBLamb(i),'%.3f'));
                elseif HBLamb(i) >= 1e4
                    String = append(String,newline,'B',num2str(i+1),'    ',num2str(HBLamb(i),'%.3f'));
                end
            else
                if  HBLamb(i) < 1e2
                    String = append(String,newline,'B',num2str(i+1),'      ',num2str(HBLamb(i),'%.3f'));
                elseif HBLamb(i) >= 1e2 && HBLamb(i) < 1e3
                    String = append(String,newline,'B',num2str(i+1),'     ',num2str(HBLamb(i),'%.3f'));
                elseif HBLamb(i) >= 1e3 && HBLamb(i) < 1e4
                    String = append(String,newline,'B',num2str(i+1),'    ',num2str(HBLamb(i),'%.3f'));
                elseif HBLamb(i) >= 1e4
                    String = append(String,newline,'B',num2str(i+1),'   ',num2str(HBLamb(i),'%.3f'));
                end
            end
        end
    end
end
String = append(String,newline,' ');
if  any(HAShear)
    for i = 1:length(HAShear)
        if  i < 10
            if  HAShear(i) < 1e2
                String = append(String,newline,'ASH',num2str(i),'     ',num2str(HAShear(i),'%.3f'));
            elseif HAShear(i) >= 1e2 && HAShear(i) < 1e3
                String = append(String,newline,'ASH',num2str(i),'    ',num2str(HAShear(i),'%.3f'));
            elseif HAShear(i) >= 1e3 && HAShear(i) < 1e4
                String = append(String,newline,'ASH',num2str(i),'   ',num2str(HAShear(i),'%.3f'));
            elseif HAShear(i) >= 1e4
                String = append(String,newline,'ASH',num2str(i),'  ',num2str(HAShear(i),'%.3f'));
            end
        else
            if  HAShear(i) < 1e2
                String = append(String,newline,'ASH',num2str(i),'    ',num2str(HAShear(i),'%.3f'));
            elseif HAShear(i) >= 1e2 && HAShear(i) < 1e3
                String = append(String,newline,'ASH',num2str(i),'   ',num2str(HAShear(i),'%.3f'));
            elseif HAShear(i) >= 1e3 && HAShear(i) < 1e4
                String = append(String,newline,'ASH',num2str(i),'  ',num2str(HAShear(i),'%.3f'));
            elseif HAShear(i) >= 1e4
                String = append(String,newline,'ASH',num2str(i),' ',num2str(HAShear(i),'%.3f'));
            end
        end
    end
end
OutputWindow1aUI1.String = String;
disp(String)
String = ['Frq. @ ',num2str(PhaseVelocityLimit/1e3),' m/ms:',newline,'Mode   Frq.(kHz)'];
if  any(HSLamb) && Symmetric
    for i = 1:length(HSLamb)
        if  i < 10
            if  HSLamb(i) < 1e2
                String = append(String,newline,'S',num2str(i),'       ',num2str(HSLamb(i),'%.3f'));
            elseif HSLamb(i) >= 1e2 && HSLamb(i) < 1e3
                String = append(String,newline,'S',num2str(i),'      ',num2str(HSLamb(i),'%.3f'));
            elseif HSLamb(i) >= 1e3 && HSLamb(i) < 1e4
                String = append(String,newline,'S',num2str(i),'     ',num2str(HSLamb(i),'%.3f'));
            elseif HSLamb(i) >= 1e4
                String = append(String,newline,'S',num2str(i),'    ',num2str(HSLamb(i),'%.3f'));
            end
        else
            if  HSLamb(i) < 1e2
                String = append(String,newline,'S',num2str(i),'      ',num2str(HSLamb(i),'%.3f'));
            elseif HSLamb(i) >= 1e2 && HSLamb(i) < 1e3
                String = append(String,newline,'S',num2str(i),'     ',num2str(HSLamb(i),'%.3f'));
            elseif HSLamb(i) >= 1e3 && HSLamb(i) < 1e4
                String = append(String,newline,'S',num2str(i),'    ',num2str(HSLamb(i),'%.3f'));
            elseif HSLamb(i) >= 1e4
                String = append(String,newline,'S',num2str(i),'   ',num2str(HSLamb(i),'%.3f'));
            end
        end
    end
end
String = append(String,newline,' ');
if  any(HSShear)
    for i = 1:length(HSShear)
        if  i < 10
            if  HSShear(i) < 1e2
                String = append(String,newline,'SSH',num2str(i),'     ',num2str(HSShear(i),'%.3f'));
            elseif HSShear(i) >= 1e2 && HSShear(i) < 1e3
                String = append(String,newline,'SSH',num2str(i),'    ',num2str(HSShear(i),'%.3f'));
            elseif HSShear(i) >= 1e3 && HSShear(i) < 1e4
                String = append(String,newline,'SSH',num2str(i),'   ',num2str(HSShear(i),'%.3f'));
            elseif HSShear(i) >= 1e4
                String = append(String,newline,'SSH',num2str(i),'  ',num2str(HSShear(i),'%.3f'));
            end
        else
            if  HSShear(i) < 1e2
                String = append(String,newline,'SSH',num2str(i),'    ',num2str(HSShear(i),'%.3f'));
            elseif HSShear(i) >= 1e2 && HSShear(i) < 1e3
                String = append(String,newline,'SSH',num2str(i),'   ',num2str(HSShear(i),'%.3f'));
            elseif HSShear(i) >= 1e3 && HSShear(i) < 1e4
                String = append(String,newline,'SSH',num2str(i),'  ',num2str(HSShear(i),'%.3f'));
            elseif HSShear(i) >= 1e4
                String = append(String,newline,'SSH',num2str(i),' ',num2str(HSShear(i),'%.3f'));
            end
        end
    end
end
OutputWindow1bUI1.String = String;
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
OutputWindow2aUI1.String = String;
disp(String)
String = '';
if  any(HSLamb) && Symmetric
    String = append(String,'S:   ',num2str(length(HSLamb)));
end
if  any(HSShear)
    String = append(String,newline,'SSH: ',num2str(length(HSShear)));
end
OutputWindow2bUI1.String = String;
disp([String,newline,'----------------'])