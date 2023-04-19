% =========================================================================
% Dispersion Calculator
% Copyright (C) 2018-2023 DLR
% Created by Armin Huber
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function [HSLamb,HALamb,HSShear,HAShear] = FrequencySweeper_Isotropic(SweepRange,Material,PhaseVelocityLimit,Half,OutputWindow1aUI1,OutputWindow1bUI1,OutputWindow2aUI1,OutputWindow2bUI1)
Resolution = 1e-5; % (kHz)

%#ok<*AGROW>
HSLamb = [];
HALamb = [];
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
a1 = (y.^2-k2).^2;
a2 = tan(x*Half);
a3 = tan(y*Half);
YSLamb = abs(real(a1./(y.*a2)+4*k2.*x./a3));
YALamb = abs(real(a1./y.*a2+4*k2.*x.*a3));
XRoughSLamb = [];
XRoughALamb = [];
for i = 2:length(YSLamb)-1
    if  YSLamb(i) < YSLamb(i-1) && YSLamb(i) < YSLamb(i+1)
        XRoughSLamb(end+1) = SweepRange(i);
    end
end
for i = 2:length(YALamb)-1
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
            Y(i) = abs(real(tan(y*Half)/y+4*k2*x*tan(x*Half)/(y^2-k2)^2));
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
            Y(i) = abs(real(y*tan(y*Half)+(y^2-k2)^2*tan(x*Half)/(4*k2*x)));
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
if  any(HSLamb)
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
if  any(HALamb)
    String = append(String,'A:   ',num2str(length(HALamb)));
end
if  any(HAShear)
    String = append(String,newline,'ASH: ',num2str(length(HAShear)));
end
OutputWindow2aUI1.String = String;
disp(String)
String = '';
if  any(HSLamb)
    String = append(String,'S:   ',num2str(length(HSLamb)));
end
if  any(HSShear)
    String = append(String,newline,'SSH: ',num2str(length(HSShear)));
end
OutputWindow2bUI1.String = String;
disp([String,newline,'----------------'])