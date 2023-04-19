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
function [BulkWaves,X,Y] = Computer_Anisotropic_SnellsLaw(Fluid,Solid,Phi,Theta,OutputWindowUI8)
%#ok<*AGROW>
ThetaRange = 0:90;
n(:,1) = cosd(ThetaRange)'.*cosd(Phi)'; % propagation direction vector matrix for the bulk waves in the solid
n(:,2) = cosd(ThetaRange)'.*sind(Phi)';
n(:,3) = sind(ThetaRange)';
C = real(Solid.C);
a21 = (C(1,1)*C(2,2)+C(1,1)*C(4,4)+C(2,2)*C(5,5)-2*C(1,2)*C(6,6)+C(4,4)*C(6,6)+C(5,5)*C(6,6)-C(1,2)^2)/Solid.Density^2;
a22 = (C(1,1)*C(3,3)+C(1,1)*C(4,4)-2*C(1,3)*C(5,5)+C(3,3)*C(6,6)+C(4,4)*C(5,5)+C(5,5)*C(6,6)-C(1,3)^2)/Solid.Density^2;
a23 = (C(2,2)*C(3,3)-2*C(2,3)*C(4,4)+C(2,2)*C(5,5)+C(3,3)*C(6,6)+C(4,4)*C(5,5)+C(4,4)*C(6,6)-C(2,3)^2)/Solid.Density^2;
a24 = (C(1,1)*C(5,5)+C(1,1)*C(6,6)+C(5,5)*C(6,6))/Solid.Density^2;
a25 = (C(2,2)*C(4,4)+C(2,2)*C(6,6)+C(4,4)*C(6,6))/Solid.Density^2;
a26 = (C(3,3)*C(4,4)+C(3,3)*C(5,5)+C(4,4)*C(5,5))/Solid.Density^2;
a31 = (C(1,1)*C(2,3)^2+C(1,3)^2*C(2,2)+C(1,2)^2*C(3,3)-2*C(1,2)*C(1,3)*C(2,3)-C(1,1)*C(2,2)*C(3,3)-2*C(1,2)*C(1,3)*C(4,4)+2*C(1,1)*C(2,3)*C(4,4)-2*C(1,2)*C(2,3)*C(5,5)+2*C(1,3)*C(2,2)*C(5,5)-2*C(1,3)*C(2,3)*C(6,6)+2*C(1,2)*C(3,3)*C(6,6)-2*C(1,2)*C(4,4)*C(5,5)-2*C(1,3)*C(4,4)*C(6,6)-2*C(2,3)*C(5,5)*C(6,6)-4*C(4,4)*C(5,5)*C(6,6))/Solid.Density^3;
a32 = (C(1,2)^2*C(4,4)-C(1,1)*C(2,2)*C(4,4)+2*C(1,2)*C(4,4)*C(6,6)-C(2,2)*C(5,5)*C(6,6))/Solid.Density^3;
a33 = (C(1,3)^2*C(4,4)-C(1,1)*C(3,3)*C(4,4)+2*C(1,3)*C(4,4)*C(5,5)-C(3,3)*C(5,5)*C(6,6))/Solid.Density^3;
a34 = (C(1,2)^2*C(5,5)-C(1,1)*C(2,2)*C(5,5)-C(1,1)*C(4,4)*C(6,6)+2*C(1,2)*C(5,5)*C(6,6))/Solid.Density^3;
a35 = (C(1,3)^2*C(6,6)-C(1,1)*C(3,3)*C(6,6)-C(1,1)*C(4,4)*C(5,5)+2*C(1,3)*C(5,5)*C(6,6))/Solid.Density^3;
a36 = (C(2,3)^2*C(5,5)-C(2,2)*C(3,3)*C(5,5)+2*C(2,3)*C(4,4)*C(5,5)-C(3,3)*C(4,4)*C(6,6))/Solid.Density^3;
a37 = (C(2,3)^2*C(6,6)-C(2,2)*C(3,3)*C(6,6)-C(2,2)*C(4,4)*C(5,5)+2*C(2,3)*C(4,4)*C(6,6))/Solid.Density^3;
a38 = -C(1,1)*C(5,5)*C(6,6)/Solid.Density^3;
a39 = -C(2,2)*C(4,4)*C(6,6)/Solid.Density^3;
a30 = -C(3,3)*C(4,4)*C(5,5)/Solid.Density^3;
A1 = -(C(1,1)*n(:,1).^2+C(2,2)*n(:,2).^2+C(3,3)*n(:,3).^2+C(4,4)*n(:,2).^2+C(4,4)*n(:,3).^2+C(5,5)*n(:,1).^2+C(5,5)*n(:,3).^2+C(6,6)*n(:,1).^2+C(6,6)*n(:,2).^2)/Solid.Density;
A2 = a21*n(:,1).^2.*n(:,2).^2+a22*n(:,1).^2.*n(:,3).^2+a23*n(:,2).^2.*n(:,3).^2+a24*n(:,1).^4+a25*n(:,2).^4+a26*n(:,3).^4;
A3 = a31*n(:,1).^2.*n(:,2).^2.*n(:,3).^2+a32*n(:,1).^2.*n(:,2).^4+a33*n(:,1).^2.*n(:,3).^4+a34*n(:,1).^4.*n(:,2).^2+a35*n(:,1).^4.*n(:,3).^2+a36*n(:,2).^2.*n(:,3).^4+a37*n(:,2).^4.*n(:,3).^2+a38*n(:,1).^6+a39*n(:,2).^6+a30*n(:,3).^6;
Xa = A2/3-A1.^2/9;
Xb = A1.^3/27-A1.*A2/6+A3/2;
Xc = (sqrt(Xb.^2+Xa.^3)-Xb).^(1/3);
Xd = Xa./(2*Xc)-Xc/2;
Xe = Xa./Xc;
Xf = (sqrt(3)*(Xc+Xe)*1i)/2;
X(:,2) = abs(real(sqrt(Xd-Xf-A1/3)))/1e3; % phase velocity in the solid (m/ms)
X(:,3) = abs(real(sqrt(Xd+Xf-A1/3)))/1e3;
X(:,4) = abs(real(-sqrt(Xc-Xe-A1/3)))/1e3;
ThetaCritical_S_fast = asind(Fluid.Velocity/1e3/X(1,2));
ThetaCritical_S_slow = asind(Fluid.Velocity/1e3/X(1,3));
ThetaCritical_L = asind(Fluid.Velocity/1e3/X(1,4));
X(:,5:7) = 1./X(:,2:4); % slowness in the solid (ms/m)
X = vertcat(X(1:end,:),flipud(X(1:end-1,:))); % extending the data from 180 deg to 360 deg
X(:,1) = pi:2*pi/(4*(length(n)-1)):2*pi; % Theta (rad)
Y(1:size(X,1),2) = 1e3/Fluid.Velocity; % slowness in the fluid (ms/m)
Y(:,1) = 0:2*pi/(4*(length(n)-1)):pi; % Theta (rad)
Phi2 = Phi;
if  Phi2 == 0
    Phi2 = 0.001;
elseif Phi2 == 45 && strcmp(Solid.Class,'Cubic')
    Phi2 = 45.001;
elseif Phi2 == 90
    Phi2 = 89.999;
end
Theta2 = Theta;
if  Theta2 == 0
    Theta2 = 0.001;
elseif Theta2 == 90
    Theta2 = 89.999;
end
s = sind(Phi2);
g = cosd(Phi2);
c(1,1) = C(1,1)*g^4+C(2,2)*s^4+2*(C(1,2)+2*C(6,6))*s^2*g^2;
c(1,2) = (C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2+C(1,2);
c(1,3) = C(1,3)*g^2+C(2,3)*s^2;
c(1,6) = (C(1,2)+2*C(6,6)-C(1,1))*s*g^3+(C(2,2)-C(1,2)-2*C(6,6))*g*s^3;
c(2,2) = C(1,1)*s^4+C(2,2)*g^4+2*(C(1,2)+2*C(6,6))*s^2*g^2;
c(2,3) = C(2,3)*g^2+C(1,3)*s^2;
c(2,6) = (C(1,2)+2*C(6,6)-C(1,1))*g*s^3+(C(2,2)-C(1,2)-2*C(6,6))*s*g^3;
c(3,3) = C(3,3);
c(3,6) = (C(2,3)-C(1,3))*s*g;
c(4,4) = C(4,4)*g^2+C(5,5)*s^2;
c(4,5) = (C(4,4)-C(5,5))*s*g;
c(5,5) = C(5,5)*g^2+C(4,4)*s^2;
c(6,6) = C(6,6)+(C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2;
Delta = c(3,3)*c(4,4)*c(5,5)-c(3,3)*c(4,5)^2;
PhaseVelocity = Fluid.Velocity/sind(Theta2); % phase velocity component along x1 shared among all participating bulk waves, as required by Snell's law
A1 = (2*(c(1,3)*c(4,5)^2+c(1,3)*c(3,6)*c(4,5)-c(1,6)*c(3,3)*c(4,5)-c(1,3)*c(4,4)*c(5,5))+c(1,1)*c(3,3)*c(4,4)-c(1,3)^2*c(4,4)-c(3,6)^2*c(5,5)+c(3,3)*c(5,5)*c(6,6)-(c(3,3)*c(4,4)+c(3,3)*c(5,5)+c(4,4)*c(5,5)-c(4,5)^2)*Solid.Density*PhaseVelocity^2)/Delta;
A2 = (2*(c(1,3)*c(1,6)*c(3,6)+c(1,3)*c(1,6)*c(4,5)+c(1,6)*c(3,6)*c(5,5)-c(1,3)*c(5,5)*c(6,6)-c(1,1)*c(3,6)*c(4,5))+c(1,1)*c(3,3)*c(6,6)+c(1,1)*c(4,4)*c(5,5)-c(1,1)*c(3,6)^2-c(1,6)^2*c(3,3)-c(1,1)*c(4,5)^2-c(1,3)^2*c(6,6)+(2*(c(1,6)*c(4,5)+c(1,3)*c(5,5)+c(3,6)*c(4,5))+c(1,3)^2+c(3,6)^2+c(4,5)^2-c(1,1)*c(3,3)-c(1,1)*c(4,4)-c(3,3)*c(6,6)-c(4,4)*c(5,5)-c(5,5)*c(6,6))*Solid.Density*PhaseVelocity^2+(c(3,3)+c(4,4)+c(5,5))*Solid.Density^2*PhaseVelocity^4)/Delta;
A3 = (c(1,1)*c(5,5)*c(6,6)-c(1,6)^2*c(5,5)-(c(1,1)*c(5,5)+c(1,1)*c(6,6)+c(5,5)*c(6,6)-c(1,6)^2)*Solid.Density*PhaseVelocity^2+(c(1,1)+c(5,5)+c(6,6))*Solid.Density^2*PhaseVelocity^4-Solid.Density^3*PhaseVelocity^6)/Delta;
Alphaa = A2/3-A1^2/9;
Alphab = A1^3/27-A1*A2/6+A3/2;
Alphac = (sqrt(Alphab^2+Alphaa^3)-Alphab)^(1/3);
Alphad = Alphaa/(2*Alphac)-Alphac/2;
Alphae = Alphaa/Alphac;
Alphaf = (sqrt(3)*(Alphac+Alphae)*1i)/2;
Alpha(1) = sqrt(Alphad-Alphaf-A1/3);
Alpha(3) = sqrt(Alphad+Alphaf-A1/3);
Alpha(2) = sqrt(Alphac-Alphae-A1/3);
BulkWaves = [atand(real(Alpha));1e3/PhaseVelocity*sqrt(1+real(Alpha).^2)];
BulkWaves(3,:) = cosd(BulkWaves(1,:))'.*cosd(Phi)'; % propagation direction vector matrix for the bulk waves in the solid
BulkWaves(4,:) = cosd(BulkWaves(1,:))'.*sind(Phi)';
BulkWaves(5,:) = sind(BulkWaves(1,:))';
BulkWaves(1,:) = 360-BulkWaves(1,:); % shift the bulk waves from the 3rd quadrant to the 4th
BulkWaves(9,1) = 0;
for i = 1:3
    if  BulkWaves(1,i) > 359.9
        BulkWaves(1,i) = 360;
    end
    if  BulkWaves(1,i) < 270.1
        BulkWaves(1,i) = 270;
    end
end
BulkWaves(BulkWaves < 1e-3) = 0;
for i = 1:3 % calculating their polarizations, i=1: S_fast, 1=2: S_slow, i=3: L
    if  BulkWaves(1,i) ~= 360 && Theta ~= 0
        B = [C(1,1)*BulkWaves(3,i)^2+C(6,6)*BulkWaves(4,i)^2+C(5,5)*BulkWaves(5,i)^2 (C(1,2)+C(6,6))*BulkWaves(3,i)*BulkWaves(4,i) (C(1,3)+C(5,5))*BulkWaves(3,i)*BulkWaves(5,i);0 C(6,6)*BulkWaves(3,i)^2+C(2,2)*BulkWaves(4,i)^2+C(4,4)*BulkWaves(5,i)^2 (C(2,3)+C(4,4))*BulkWaves(4,i)*BulkWaves(5,i);0 0 C(5,5)*BulkWaves(3,i)^2+C(4,4)*BulkWaves(4,i)^2+C(3,3)*BulkWaves(5,i)^2]/Solid.Density;
        if  Phi == 0 && ~strcmp(Solid.Class,'Cubic') && i == 2 % S_slow is polarized only along x2; U2 = 1; [1] p. 35 
            BulkWaves(6,i) = 0; % bulk wave amplitude U1, i.e., polarization component along x1
            BulkWaves(7,i) = 1; % bulk wave amplitude U2, i.e., polarization component along x2
            BulkWaves(8,i) = 0; % bulk wave amplitude U3, i.e., polarization component along x3
        elseif Phi == 0 && strcmp(Solid.Class,'Cubic') && i == 1 % S_fast
            BulkWaves(6,i) = 0;
            BulkWaves(7,i) = 1;
            BulkWaves(8,i) = 0;
        elseif Phi == 90 && i == 1 % S_fast
            BulkWaves(6,i) = 1;
            BulkWaves(7,i) = 0;
            BulkWaves(8,i) = 0;
        elseif Phi == 90 && strcmp(Solid.Class,'Cubic') && i == 2 % S_slow
            BulkWaves(6,i) = 0;
            BulkWaves(7,i) = 1;
            BulkWaves(8,i) = -(B(2,2)-(1/BulkWaves(2,i)*1e3)^2)/B(2,3);
        elseif Phi == 90 && i == 3 % L
            BulkWaves(6,i) = 0;
            BulkWaves(7,i) = 1;
            BulkWaves(8,i) = -(B(2,2)-(1/BulkWaves(2,i)*1e3)^2)/B(2,3); % bulk wave amplitude ratio U3/U2, i.e., polarization component along x3
        else
            if  strcmp(Solid.Class,'Cubic') || strcmp(Solid.Class,'Transversely isotropic') && i ~= 2 || strcmp(Solid.Class,'Orthotropic') && Phi ~= 90
                BulkWaves(6,i) = 1; % bulk wave amplitude U1, i.e., polarization component along x1; U1 = 1
                BulkWaves(7,i) = (B(2,3)*(B(1,1)-(1/BulkWaves(2,i)*1e3)^2)-B(1,3)*B(1,2))/(B(1,3)*(B(2,2)-(1/BulkWaves(2,i)*1e3)^2)-B(1,2)*B(2,3)); % bulk wave amplitude ratio U2/U1, i.e., polarization component along x2
                BulkWaves(8,i) = ((B(1,1)-(1/BulkWaves(2,i)*1e3)^2)*(B(2,2)-(1/BulkWaves(2,i)*1e3)^2)-B(1,2)^2)/(B(1,2)*B(2,3)-B(1,3)*(B(2,2)-(1/BulkWaves(2,i)*1e3)^2)); % bulk wave amplitude ratio U3/U1, i.e., polarization component along x3
            else
                BulkWaves(6,i) = 0;
                BulkWaves(7,i) = 1; % S_slow is always polarized normal to the fibers, i.e., in the x2-x3-plane only
                BulkWaves(8,i) = -(B(2,2)-(1/BulkWaves(2,i)*1e3)^2)/B(2,3); % bulk wave amplitude ratio U3/U2, i.e., polarization component along x3
            end
        end
        if  Phi == 0 && i == 2
            BulkWaves(9,i) = 90;
        elseif Phi == 90 && i == 3 && strcmp(Solid.Class,'Transversely isotropic')
            BulkWaves(9,i) = 0;
        else
            BulkWaves(9,i) = acosd(dot([BulkWaves(6,i) BulkWaves(7,i) BulkWaves(8,i)],[BulkWaves(3,i) BulkWaves(4,i) BulkWaves(5,i)])/(sqrt(BulkWaves(6,i)^2+BulkWaves(7,i)^2+BulkWaves(8,i)^2))); % polarization angles
        end
    elseif BulkWaves(1,i) ~= 360 && Theta == 0
        BulkWaves(6,1) = 1; % S_fast
        BulkWaves(7,1) = 0;
        BulkWaves(8,1) = 0;
        BulkWaves(9,1) = 90;
        BulkWaves(6,2) = 0; % S_slow
        BulkWaves(7,2) = 1;
        BulkWaves(8,2) = 0;
        BulkWaves(9,2) = 90;
        BulkWaves(6,3) = 0; % L
        BulkWaves(7,3) = 0;
        BulkWaves(8,3) = 1;
        BulkWaves(9,3) = 0;
        break
    end
    if  BulkWaves(1,i) == 360 % remove those which are evanescent
        BulkWaves(2,i) = sind(Theta)*1e3/Fluid.Velocity;
    end
end
BulkWaves = abs(BulkWaves);
% content of "BulkWaves": 
% columns 1: S_fast 2: S_slow 3: L
% lines: 
% 1: polar angle
% 2: slowness
% 3: propagation direction component n1 
% 4: propagation direction component n2 
% 5: propagation direction component n3
% 6: polarization component along x1
% 7: polarization component along x2
% 8: polarization component along x3
% 9: angle between propagation direction and polarization
AlphaFluid = 1/Fluid.Velocity*PhaseVelocity*cosd(Theta2);
if  angle(Alpha(1)) < 0
    Alpha(1) = conj(Alpha(1)); % T(S_fast)
end
if  angle(Alpha(2)) < 0
    Alpha(2) = conj(Alpha(2)); % T(S_slow)
end
if  angle(Alpha(3)) < 0
    Alpha(3) = conj(Alpha(3)); % T(L)
end
if  Theta2 > ThetaCritical_S_slow
    Alpha(2) = -Alpha(2); % T(S_slow)
end
if  strcmp(Solid.Class,'Cubic') && Theta2 > ThetaCritical_S_fast
    Alpha(3) = -conj(Alpha(3)); % T(L)
end
m11 = c(1,1)-Solid.Density*PhaseVelocity^2+c(5,5)*Alpha.^2;
m12 = c(1,6)+c(4,5)*Alpha.^2;
m13 = (c(1,3)+c(5,5))*Alpha;
m22 = c(6,6)-Solid.Density*PhaseVelocity^2+c(4,4)*Alpha.^2;
m23 = (c(3,6)+c(4,5))*Alpha;
V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
D = 1i/PhaseVelocity*[0 0 0;0 0 0;c(1,3)+c(3,6)*V+c(3,3)*Alpha.*W;c(4,5)*(Alpha+W)+c(4,4)*Alpha.*V;c(5,5)*(Alpha+W)+c(4,5)*Alpha.*V]; % sigma33 sigma23 sigma13
DFluid = 1i*Fluid.Density*PhaseVelocity; % sigmaFluid33
Z1 = [W AlphaFluid;D(3,:) -DFluid;D(5,:) 0;D(4,:) 0];
Z2 = [AlphaFluid;DFluid;0;0];
U = Z1\Z2; % UT(S_fast) UT(S_slow) UT(L) UR(L)
v = -1i*[1 0 AlphaFluid;U(1) V(1)*U(1) W(1)*U(1);U(2) V(2)*U(2) W(2)*U(2);U(3) V(3)*U(3) W(3)*U(3);U(4) 0 AlphaFluid*U(4)]; % I(L) T(S_fast) T(S_slow) T(L) R(L)
sigma = [DFluid 0 0;D(3,1)*U(1) D(5,1)*U(1) D(4,1)*U(1);D(3,2)*U(2) D(5,2)*U(2) D(4,2)*U(2);D(3,3)*U(3) D(5,3)*U(3) D(4,3)*U(3);DFluid*U(4) 0 0]; % I(L) T(S_fast) T(S_slow) T(L) R(L)
P = -.5*[real(sigma(1,1)*conj(v(1,3)));real(sigma(2,1)*conj(v(2,3))+sigma(2,2)*conj(v(2,1))+sigma(2,3)*conj(v(2,2)));real(sigma(3,1)*conj(v(3,3))+sigma(3,2)*conj(v(3,1))+sigma(3,3)*conj(v(3,2)));real(sigma(4,1)*conj(v(4,3))+sigma(4,2)*conj(v(4,1))+sigma(4,3)*conj(v(4,2)));real(sigma(5,1)*conj(v(5,3)))]; % I(L) T(S_fast) T(S_slow) T(L) R(L)
E(1:4) = P(2:5)/P(1); % T(S_fast) T(S_slow) T(L) R(L) Sum
if  E(4) > 1.001
    E = [0 0 0 1 0];
end
E(5) = sum(E,2);
if  abs(E(1)) < 1e-7
    E(1) = 0;
end
if  abs(E(2)) < 1e-7
    E(2) = 0;
end
if  abs(E(3)) < 1e-7
    E(3) = 0;
end
if  abs(E(4)) < 1e-7
    E(4) = 0;
end
format long
String = ['Fluid: ',Fluid.Name,newline,'Solid: ',Solid.Name,newline,newline,'Critical angles @ Phi = ',num2str(Phi),' deg:'];
if  isreal(ThetaCritical_L)
    String = append(String,newline,'  Theta_L = ',num2str(ThetaCritical_L),' deg');
else
    disp('  Theta_L = -')
end
if  isreal(ThetaCritical_S_fast)
    String = append(String,newline,'  Theta_Sfast = ',num2str(ThetaCritical_S_fast),' deg');
else
    String = append(String,newline,'  Theta_Sfast = -');
end
if  isreal(ThetaCritical_S_slow)
    String = append(String,newline,'  Theta_Sslow = ',num2str(ThetaCritical_S_slow),' deg');
else
    String = append(String,newline,'  Theta_Sslow = -');
end
String = append(String,newline,newline,...
'Plane wave incidence:',newline,...
'  Phi = ',num2str(Phi),' deg',newline,...
'  Theta = ',num2str(Theta),' deg',newline,...
newline,...
'(QUASI) LONGITUDINAL WAVE (L)',newline,...
'Refraction angle = ',num2str(BulkWaves(1,3)-270),' deg',newline,...
'Phase velocity = ',num2str(1/BulkWaves(2,3)),' m/ms',newline,...
'Phase velocity x''1 = ',num2str(cosd(BulkWaves(1,3))*1/BulkWaves(2,3)*cosd(Phi)),' m/ms',newline,...
'Phase velocity x''2 = ',num2str(cosd(BulkWaves(1,3))*1/BulkWaves(2,3)*sind(Phi)),' m/ms',newline,...
'Phase velocity x''3 = ',num2str(-sind(BulkWaves(1,3))*1/BulkWaves(2,3)),' m/ms');
if  BulkWaves(1,3) ~= 360 
    String = append(String,newline,...
    'Polarization x''1 = ',num2str(BulkWaves(6,3)),newline,...
    'Polarization x''2 = ',num2str(BulkWaves(7,3)),newline,...
    'Polarization x''3 = ',num2str(BulkWaves(8,3)),newline,...
    'Polarization angle = ',num2str(BulkWaves(9,3)),' deg');
else
    String = append(String,newline,'Evanescent surface wave');
end
String = append(String,newline,newline,...
'FAST (QUASI) SHEAR WAVE (Sfast)',newline,...
'Refraction angle = ',num2str(BulkWaves(1,1)-270),' deg',newline,...
'Phase velocity = ',num2str(1/BulkWaves(2,1)),' m/ms',newline,...
'Phase velocity x''1 = ',num2str(cosd(BulkWaves(1,1))*1/BulkWaves(2,1)*cosd(Phi)),' m/ms',newline,...
'Phase velocity x''2 = ',num2str(cosd(BulkWaves(1,1))*1/BulkWaves(2,1)*sind(Phi)),' m/ms',newline,...
'Phase velocity x''3 = ',num2str(-sind(BulkWaves(1,1))*1/BulkWaves(2,1)),' m/ms');
if  BulkWaves(1,1) ~= 360
    String = append(String,newline,...
    'Polarization x''1 = ',num2str(BulkWaves(6,1)),newline,...
    'Polarization x''2 = ',num2str(BulkWaves(7,1)),newline,...
    'Polarization x''3 = ',num2str(BulkWaves(8,1)),newline,...
    'Polarization angle = ',num2str(180-BulkWaves(9,1)),' deg');
else
    String = append(String,newline,'Evanescent surface wave');
end
if  strcmp(Solid.Class,'Transversely isotropic')
    String = append(String,newline,newline,...
    'SLOW SHEAR WAVE (Sslow)');
else
    String = append(String,newline,newline,...
    'SLOW (QUASI) SHEAR WAVE (Sslow)');
 end
String = append(String,newline,...
'Refraction angle = ',num2str(BulkWaves(1,2)-270),' deg',newline,...
'Phase velocity = ',num2str(1/BulkWaves(2,2)),' m/ms',newline,...
'Phase velocity x''1 = ',num2str(cosd(BulkWaves(1,2))*1/BulkWaves(2,2)*cosd(Phi)),' m/ms',newline,...
'Phase velocity x''2 = ',num2str(cosd(BulkWaves(1,2))*1/BulkWaves(2,2)*sind(Phi)),' m/ms',newline,...
'Phase velocity x''3 = ',num2str(-sind(BulkWaves(1,2))*1/BulkWaves(2,2)),' m/ms');
if  BulkWaves(1,2) ~= 360
    String = append(String,newline,...
    'Polarization x''1 = ',num2str(BulkWaves(6,2)),newline,...
    'Polarization x''2 = ',num2str(BulkWaves(7,2)),newline,...
    'Polarization x''3 = ',num2str(BulkWaves(8,2)),newline,...
    'Polarization angle = ',num2str(BulkWaves(9,2)),' deg');
else
    String = append(String,newline,'Evanescent surface wave');
end
String = append(String,newline,newline,...
'Energy scattering coefficients:',newline,...
'  R(L) = ',num2str(E(4)),newline,...
'  T(L) = ',num2str(E(3)),newline,...
'  T(Sfast) = ',num2str(E(1)),newline,...
'  T(Sslow) = ',num2str(E(2)),newline,...
'  Sum = ',num2str(E(5)));
OutputWindowUI8.String = String;
disp([String,newline,'---------------------------------'])