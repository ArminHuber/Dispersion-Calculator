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
function [HSScholte,HAScholte] = FrequencySweeper_Isotropic_Scholte(Material,Half,Fluid,SweepRange)
Resolution = 1e-5; % (kHz)
PhaseVelocity = (1-1e-4)*Fluid.Velocity;

%#ok<*AGROW>
HSScholte = [];
HAScholte = [];
if  Resolution < abs(SweepRange(1)-SweepRange(2))
    Bisections = ceil(log2(Resolution/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
else
    Bisections = 1;
end
AngularFrequency = 2*pi*SweepRange*1e3;
k2 = (AngularFrequency/PhaseVelocity).^2;
kL2 = (AngularFrequency/Material.LongitudinalVelocity).^2;
kT2 = (AngularFrequency/Material.TransverseVelocity).^2;
kF2 = (AngularFrequency/Fluid.Velocity).^2;
x = sqrt(kL2-k2);
y = sqrt(kT2-k2);
a1 = (y.^2-k2).^2;
a2 = tan(x*Half);
a3 = tan(y*Half);
F = Fluid.Density*kT2.^2.*x./(y*Material.Density.*sqrt(kF2-k2));
YSScholte = abs(a1./(y.*a2)+4*k2.*x./a3-1i*F);
YAScholte = abs(a1./y.*a2+4*k2.*x.*a3+1i*F);
XRoughSScholte = [];
XRoughAScholte = [];
for i = 2:length(YSScholte)-1
    if  YSScholte(i) < YSScholte(i-1) && YSScholte(i) < YSScholte(i+1)
        XRoughSScholte(end+1) = SweepRange(i);
    end
end
for i = 2:length(YAScholte)-1
    if  YAScholte(i) < YAScholte(i-1) && YAScholte(i) < YAScholte(i+1)
        XRoughAScholte(end+1)  = SweepRange(i);
    end
end
if  isempty(XRoughAScholte)
    XRoughSScholte = [];
else
    XRoughSScholte(XRoughSScholte < XRoughAScholte(1)) = [];
end
if  isempty(XRoughSScholte) && isempty(XRoughAScholte)
    disp('No higher order Scholte modes found!')
    return
end
for j = 1:length(XRoughSScholte)
    Frequency = [XRoughSScholte(j)-(SweepRange(2)-SweepRange(1)) XRoughSScholte(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            k2 = (AngularFrequency(i)/PhaseVelocity)^2;
            kT = AngularFrequency(i)/Material.TransverseVelocity;
            kF = AngularFrequency(i)/Fluid.Velocity;
            x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
            y = sqrt(kT^2-k2);
            Y(i) = abs((y^2-k2)^2/(y*tan(x*Half))+4*k2*x/tan(y*Half)-1i*Fluid.Density*kT^4*x/(y*Material.Density*sqrt(kF^2-k2)));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HSScholte(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
end
for j = 1:length(XRoughAScholte)
    Frequency = [XRoughAScholte(j)-(SweepRange(2)-SweepRange(1)) XRoughAScholte(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        Frequency = Frequency(1):.25*(Frequency(end)-Frequency(1)):Frequency(end);
        AngularFrequency = 2*pi*Frequency*1e3;
        for i = 1:length(Frequency)
            k2 = (AngularFrequency(i)/PhaseVelocity)^2;
            kT = AngularFrequency(i)/Material.TransverseVelocity;
            kF = AngularFrequency(i)/Fluid.Velocity;
            x = sqrt((AngularFrequency(i)/Material.LongitudinalVelocity)^2-k2);
            y = sqrt(kT^2-k2);
            Y(i) = abs((y^2-k2)^2/y*tan(x*Half)+4*k2*x*tan(y*Half)+1i*Fluid.Density*kT^4*x/(y*Material.Density*sqrt(kF^2-k2)));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    HAScholte(end+1) = Frequency(i-1);
                end
                Frequency = [Frequency(i-1)-(Frequency(2)-Frequency(1)) Frequency(i-1)+(Frequency(2)-Frequency(1))];
                break
            end
        end
    end
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
disp([String,newline,'----------------'])