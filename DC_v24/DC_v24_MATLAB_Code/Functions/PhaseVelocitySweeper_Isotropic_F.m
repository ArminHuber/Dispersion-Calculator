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
function [FALambF,FAScholte] = PhaseVelocitySweeper_Isotropic_F(Material,Half,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Frequency)
RangeReal = 2;
RangeImag = -10;
StepsReal = 6e2;
StepsImag = 6e2;

%#ok<*AGROW>
FALambF = [];
FAScholte = [];
if  FluidLoading && Viscoelastic
    if  Fluid.Density < FluidDensityThreshold 
        T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material.Density*Material.PlateVelocity);
    else
        T = 1/(Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material.Density*Material.PlateVelocity));
    end
    T = T+Material.LongitudinalAttenuation+Material.TransverseAttenuation;
elseif FluidLoading && ~Viscoelastic
    if  Fluid.Density < FluidDensityThreshold
        T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material.Density*Material.PlateVelocity);
    else
        T = 1/(Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material.Density*Material.PlateVelocity));
    end
elseif ~FluidLoading && Viscoelastic
    T = Material.LongitudinalAttenuation+Material.TransverseAttenuation;
end
AngularFrequency = 2*pi*Frequency*1e3;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
SweepRangeReal = 1:500;
k2 = (AngularFrequency./SweepRangeReal).^2;
x = sqrt(kL2-k2);
y = sqrt(kT2-k2);
Y = abs(y.*tan(y*Half)+(y.^2-k2).^2.*tan(x*Half)./(4*k2.*x));
for i = 2:length(Y)-1
    if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
        XRough = SweepRangeReal(i);
        break
    end
end
PhaseVelocity = [XRough-(SweepRangeReal(2)-SweepRangeReal(1)) XRough+(SweepRangeReal(2)-SweepRangeReal(1))];
Bisections = ceil(log2(1e-6/(abs(SweepRangeReal(1)-SweepRangeReal(2))))/log2(2*.25));
for o = 1:Bisections
    PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
    k2 = (AngularFrequency./PhaseVelocity).^2;
    for j = 1:length(PhaseVelocity)
        x = sqrt(kL2-k2(j));
        y = sqrt(kT2-k2(j));
        Y(j) = abs(y*tan(y*Half)+(y^2-k2(j))^2*tan(x*Half)/(4*k2(j)*x));
        if  j > 2 && Y(j-1) < Y(j-2) && Y(j-1) < Y(j)
            if  o == Bisections
                FALamb = PhaseVelocity(j-1);
            end
            PhaseVelocity = [PhaseVelocity(j-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(j-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
            break
        end
    end
end
SweepRangeReal = 1:RangeReal*FALamb/StepsReal:RangeReal*FALamb;
SweepRangeImag = 0:RangeImag*T/StepsImag:RangeImag*T;
WaveNumberComplex = AngularFrequency./(repmat(SweepRangeReal',1,length(SweepRangeImag))+repmat(SweepRangeImag,length(SweepRangeReal),1)*1i);
k2 = WaveNumberComplex.^2;
x = sqrt(kL2-k2);
y = sqrt(kT2-k2);
if  FluidLoading
    kF2 = (AngularFrequency/Fluid.Velocity)^2;
    Y = abs((y.^2-k2).^2./y.*tan(x*Half)+4*k2.*x.*tan(y*Half)+1i*Fluid.Density*kT2.^2.*x./(y.*Material.Density.*sqrt(kF2-k2)));
else
    Y = abs((y.^2-k2).^2./y.*tan(x*Half)+4*k2.*x.*tan(y*Half));
end
Y = horzcat(zeros(size(Y,1),1),Y);
Y(:,1) = max(max(Y));
% figure;surf(horzcat(-SweepRangeImag(2),SweepRangeImag(1:size(Y,2)-1)),SweepRangeReal(1:size(Y,1)),20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
Min = zeros(size(Y,1),size(Y,2));
for l = 2:size(Y,2)-1
    for j = 2:size(Y,1)-1
        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
            Min(j,l) = 1;
        end
    end
end
[b1,b2] = find(Min);
if  ~isempty(b1)
    for j = 1:length(b1)
        XFRough(j,1) = SweepRangeReal(b1(j)); % real phase velocity (m/s)    
        XFRough(j,2) = SweepRangeImag(b2(j)); % imaginary phase velocity (m/s)
    end
    Bisections = ceil(log2(1e-6/(RangeReal*FALamb/StepsReal))/log2(2*.25));
    for p = 1:size(XFRough,1)
        SweepRangeReal = [XFRough(p,1)-RangeReal*FALamb/StepsReal XFRough(p,1)+RangeReal*FALamb/StepsReal];
        SweepRangeImag = [XFRough(p,2)+RangeImag*T/StepsImag XFRough(p,2)-RangeImag*T/StepsImag];
        if  SweepRangeImag(2) > 0
            SweepRangeImag(2) = 0;
        end
        for k = 1:Bisections % search minimum in characteristic equation and converge upon it
            SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
            SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
            WaveNumberComplex = AngularFrequency./(repmat(SweepRangeReal',1,length(SweepRangeImag))+repmat(SweepRangeImag,length(SweepRangeReal),1)*1i);
            k2 = WaveNumberComplex.^2;
            x = sqrt(kL2-k2);
            y = sqrt(kT2-k2);
            if  FluidLoading
                kF2 = (AngularFrequency/Fluid.Velocity)^2;
                Y = abs((y.^2-k2).^2./y.*tan(x*Half)+4*k2.*x.*tan(y*Half)+1i*Fluid.Density*kT2.^2.*x./(y.*Material.Density.*sqrt(kF2-k2)));
            else
                Y = abs((y.^2-k2).^2./y.*tan(x*Half)+4*k2.*x.*tan(y*Half));
            end
            Y(:,end+1) = max(max(Y));
% if k==10
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))
% close(f)
% end
            Min = zeros(size(Y,1),size(Y,2));
            for l = 2:size(Y,2)-1
                for j = 2:size(Y,1)-1
                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                        Min(j,l) = 1;
                    end
                end
            end
            [b1,b2] = find(Min);
            if  isempty(b1)
                break
            end
            MIN = [b1(1) b2(1)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
            if  k < Bisections % set the new search area around the found minimum
                SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
            end
        end
        X(p,1) = AngularFrequency/real(WaveNumberComplex(MIN(1),MIN(2))); % phase velocity (m/s)
        X(p,2) = 2*pi*imag(WaveNumberComplex(MIN(1),MIN(2)))/real(WaveNumberComplex(MIN(1),MIN(2)))/(X(p,1)/(Frequency*1e3)); % attenuation (Np/m)
        X(p,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
        X(p,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
    end
    FAScholte = X(X(:,2) == 0,:);
    FALambF = X(X(:,2) ~= 0,:);
    if  isempty(FAScholte)
        FAScholte = [FALamb 0 FALamb 0];
    end
    if  isempty(FALambF)
        FALambF = [FAScholte 0 FAScholte 0];
    end
end
% String = ['Fluid density: ',num2str(Fluid.Density),newline,'Fluid velocity: ',num2str(Fluid.Velocity),newline,'A0 undamped: ',num2str(FALamb),newline];
% for i = 1:size(FALambF,1)
%     String = append(String,'A0_',num2str(i-1),': ',num2str(FALambF(i,1)),', ',num2str(FALambF(i,2)),', ',num2str(FALambF(i,3)),', ',num2str(FALambF(i,4)),newline);    
% end
% disp([String,'A0_Scholte: ',num2str(FAScholte(1)),newline,'-----------------------'])