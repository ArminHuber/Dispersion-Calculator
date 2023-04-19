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
function BulkWaves3D(Crop,FontSizeModeLabels,Material,FontSizeAxesLabels,BoxLineWidth,ViewTheta,ViewPhi,Export,PDF,PNG,FileName,Directory,HeadLine,FontSizeHeadLine,FontSizeAxes,WaveVectorLineWidth,PNGresolution,Phi,Theta)
LArrowTheta = 2;
LArrowLength = .1;

%#ok<*AGROW>
C = real(Material.C);
n = [cosd(Theta)*cosd(Phi) cosd(Theta)*sind(Phi) sind(Theta)]; % propagation direction vector of the bulk waves
A1 = -(C(1,1)*n(1)^2+C(2,2)*n(2)^2+C(3,3)*n(3)^2+C(4,4)*n(2)^2+C(4,4)*n(3)^2+C(5,5)*n(1)^2+C(5,5)*n(3)^2+C(6,6)*n(1)^2+C(6,6)*n(2)^2)/Material.Density;
A2 = ((C(1,1)*C(2,2)+C(1,1)*C(4,4)+C(2,2)*C(5,5)-2*C(1,2)*C(6,6)+C(4,4)*C(6,6)+C(5,5)*C(6,6)-C(1,2)^2)*n(1)^2*n(2)^2+(C(1,1)*C(3,3)+C(1,1)*C(4,4)-2*C(1,3)*C(5,5)+C(3,3)*C(6,6)+C(4,4)*C(5,5)+C(5,5)*C(6,6)-C(1,3)^2)*n(1)^2*n(3)^2+(C(2,2)*C(3,3)-2*C(2,3)*C(4,4)+C(2,2)*C(5,5)+C(3,3)*C(6,6)+C(4,4)*C(5,5)+C(4,4)*C(6,6)-C(2,3)^2)*n(2)^2*n(3)^2+(C(1,1)*C(5,5)+C(1,1)*C(6,6)+C(5,5)*C(6,6))*n(1)^4+(C(2,2)*C(4,4)+C(2,2)*C(6,6)+C(4,4)*C(6,6))*n(2)^4+(C(3,3)*C(4,4)+C(3,3)*C(5,5)+C(4,4)*C(5,5))*n(3)^4)/Material.Density^2;
A3 = ((C(1,1)*C(2,3)^2+C(1,3)^2*C(2,2)+C(1,2)^2*C(3,3)-2*C(1,2)*C(1,3)*C(2,3)-C(1,1)*C(2,2)*C(3,3)-2*C(1,2)*C(1,3)*C(4,4)+2*C(1,1)*C(2,3)*C(4,4)-2*C(1,2)*C(2,3)*C(5,5)+2*C(1,3)*C(2,2)*C(5,5)-2*C(1,3)*C(2,3)*C(6,6)+2*C(1,2)*C(3,3)*C(6,6)-2*C(1,2)*C(4,4)*C(5,5)-2*C(1,3)*C(4,4)*C(6,6)-2*C(2,3)*C(5,5)*C(6,6)-4*C(4,4)*C(5,5)*C(6,6))*n(1)^2*n(2)^2*n(3)^2+(C(1,2)^2*C(4,4)-C(1,1)*C(2,2)*C(4,4)+2*C(1,2)*C(4,4)*C(6,6)-C(2,2)*C(5,5)*C(6,6))*n(1)^2*n(2)^4+(C(1,3)^2*C(4,4)-C(1,1)*C(3,3)*C(4,4)+2*C(1,3)*C(4,4)*C(5,5)-C(3,3)*C(5,5)*C(6,6))*n(1)^2*n(3)^4+(C(1,2)^2*C(5,5)-C(1,1)*C(2,2)*C(5,5)-C(1,1)*C(4,4)*C(6,6)+2*C(1,2)*C(5,5)*C(6,6))*n(1)^4*n(2)^2+(C(1,3)^2*C(6,6)-C(1,1)*C(3,3)*C(6,6)-C(1,1)*C(4,4)*C(5,5)+2*C(1,3)*C(5,5)*C(6,6))*n(1)^4*n(3)^2+(C(2,3)^2*C(5,5)-C(2,2)*C(3,3)*C(5,5)+2*C(2,3)*C(4,4)*C(5,5)-C(3,3)*C(4,4)*C(6,6))*n(2)^2*n(3)^4+(C(2,3)^2*C(6,6)-C(2,2)*C(3,3)*C(6,6)-C(2,2)*C(4,4)*C(5,5)+2*C(2,3)*C(4,4)*C(6,6))*n(2)^4*n(3)^2-C(1,1)*C(5,5)*C(6,6)*n(1)^6-C(2,2)*C(4,4)*C(6,6)*n(2)^6-C(3,3)*C(4,4)*C(5,5)*n(3)^6)/Material.Density^3;
Xa = A2/3-A1^2/9;
Xb = A1^3/27-A1*A2/6+A3/2;
Xc = (sqrt(Xb^2+Xa^3)-Xb)^(1/3);
Xd = Xa/(2*Xc)-Xc/2;
Xe = Xa/Xc;
Xf = (sqrt(3)*(Xc+Xe)*1i)/2;
X(1) = abs(real(sqrt(Xd-Xf-A1/3)));
X(2) = abs(real(sqrt(Xd+Xf-A1/3)));
X(3) = abs(real(-sqrt(Xc-Xe-A1/3)));
B = [C(1,1)*n(1)^2+C(6,6)*n(2)^2+C(5,5)*n(3)^2 (C(1,2)+C(6,6))*n(1)*n(2) (C(1,3)+C(5,5))*n(1)*n(3);0 C(6,6)*n(1)^2+C(2,2)*n(2)^2+C(4,4)*n(3)^2 (C(2,3)+C(4,4))*n(2)*n(3);0 0 C(5,5)*n(1)^2+C(4,4)*n(2)^2+C(3,3)*n(3)^2]/Material.Density;   
for i = 1:3 % calculating their polarizations, i=1: S_fast, i=2: S_slow, i=3: L
    if  Theta ~= 0 && Theta ~= 90
        if  Phi == 0 && ~strcmp(Material.Class,'Cubic') && i == 2 % S_slow is polarized only along x2; U2 = 1; [1] p. 35 
            p(1,i) = 0; % bulk wave amplitude U1, i.e., polarization component along x1
            p(2,i) = 1; % bulk wave amplitude U2, i.e., polarization component along x2
            p(3,i) = 0; % bulk wave amplitude U3, i.e., polarization component along x3
        elseif Phi == 0 && strcmp(Material.Class,'Cubic') && i == 1 % S_fast
            p(1,i) = 0;
            p(2,i) = 1;
            p(3,i) = 0;
        elseif Phi == 90 && i == 1 % S_fast
            p(1,i) = 1;
            p(2,i) = 0;
            p(3,i) = 0;
        elseif Phi == 90 && strcmp(Material.Class,'Cubic') && i == 2 % S_slow
            p(1,i) = 0;
            p(2,i) = 1;
            p(3,i) = -(B(2,2)-X(i)^2)/B(2,3);
        elseif Phi == 90 && i == 3 % L
            p(1,i) = 0;
            p(2,i) = 1;
            p(3,i) = -(B(2,2)-X(i)^2)/B(2,3); % bulk wave amplitude ratio U3/U2, i.e., polarization component along x3
        else
            if  strcmp(Material.Class,'Isotropic') && i ~= 2 || strcmp(Material.Class,'Cubic') || strcmp(Material.Class,'Transversely isotropic') && i ~= 2 || strcmp(Material.Class,'Orthotropic') && Phi ~= 90
                p(1,i) = 1; % bulk wave amplitude U1, i.e., polarization component along x1; U1 = 1
                p(2,i) = (B(2,3)*(B(1,1)-X(i)^2)-B(1,3)*B(1,2))/(B(1,3)*(B(2,2)-X(i)^2)-B(1,2)*B(2,3)); % bulk wave amplitude ratio U2/U1, i.e., polarization component along x2
                p(3,i)  = ((B(1,1)-X(i)^2)*(B(2,2)-X(i)^2)-B(1,2)^2)/(B(1,2)*B(2,3)-B(1,3)*(B(2,2)-X(i)^2)); % bulk wave amplitude ratio U3/U1, i.e., polarization component along x3
            elseif strcmp(Material.Class,'Transversely isotropic') && i == 2 || strcmp(Material.Class,'Orthotropic') && Phi == 90
                p(1,i) = 0;
                p(2,i) = 1; % S_slow is always polarized normal to the fibers in transversely isotropic media, i.e., in the x2-x3-plane only
                p(3,i) = -(B(2,2)-X(i)^2)/B(2,3); % bulk wave amplitude ratio U3/U2, i.e., polarization component along x3
            elseif strcmp(Material.Class,'Isotropic') && i == 2
                p(:,i) = cross(p(:,1),n)';
            end
        end
        if  Phi == 0 && i == 2
            Beta(i) = 90;
        else
            Beta(i) = acosd(dot(n,p(:,i))/norm(p(:,i))); % polarization angles
        end
    elseif Theta == 0
        if  Phi == 0
            p(1,1) = 0; % S_fast
            p(2,1) = 0;
            p(3,1) = 1;
            Beta(1) = 90;
            p(1,2) = 0; % S_slow
            p(2,2) = 1;
            p(3,2) = 0;
            Beta(2) = 90;
            p(1,3) = 1; % L
            p(2,3) = 0;
            p(3,3) = 0;
            Beta(3) = 0;
        elseif Phi == 90
            p(1,1) = 1; % S_fast
            p(2,1) = 0;
            p(3,1) = 0;
            Beta(1) = 90;
            p(1,2) = 0; % S_slow
            p(2,2) = 0;
            p(3,2) = 1;
            Beta(2) = 90;
            p(1,3) = 0; % L
            p(2,3) = 1;
            p(3,3) = 0;
            Beta(3) = 0;
        else
            if  strcmp(Material.Class,'Cubic')
                p(1,1) = 0; % S_fast
                p(2,1) = 0;
                p(3,1) = 1;
                Beta(1) = acosd(dot(n,p(:,1))/norm(p(:,1)));
                p(1,2) = 1; % S_slow
                p(2,2) = -(B(1,1)-X(2)^2)/B(1,2);
                p(3,2) = 0;
                Beta(1) = acosd(dot(n,p(:,1))/norm(p(:,1)));
            else
                p(1,1) = 1; % S_fast
                p(2,1) = -(B(1,1)-X(1)^2)/B(1,2);
                p(3,1) = 0;
                Beta(1) = acosd(dot(n,p(:,1))/norm(p(:,1)));
                p(1,2) = 0; % S_slow
                p(2,2) = 0;
                p(3,2) = 1;
                Beta(2) = 90;
            end
            p(1,3) = 1; % L
            p(2,3) = -(B(1,1)-X(3)^2)/B(1,2);
            p(3,3) = 0;
            Beta(3) = acosd(dot(n,p(:,3))/norm(p(:,3)));
        end
        break
    elseif Theta == 90
        p(1,1) = 1; % S_fast
        p(2,1) = 0;
        p(3,1) = 0;
        Beta(1) = 90;
        p(1,2) = 0; % S_slow
        p(2,2) = 1;
        p(3,2) = 0;
        Beta(2) = 90;
        p(1,3) = 0; % L
        p(2,3) = 0;
        p(3,3) = 1;
        Beta(3) = 0;
        break
    end
end
% disp(['L,S_fast:      ',num2str(acosd(dot(p(:,3),p(:,1))/(norm(p(:,3))*norm(p(:,1))))),' deg']) % angles between polarizations of S_fast,S_slow,L must be 90 deg
% disp(['L,S_slow:      ',num2str(acosd(dot(p(:,3),p(:,2))/(norm(p(:,3))*norm(p(:,2))))),' deg'])
% disp(['S_fast,S_slow: ',num2str(acosd(dot(p(:,1),p(:,2))/(norm(p(:,1))*norm(p(:,2))))),' deg'])
SfArrowTheta = X(1)/X(3)*LArrowTheta;
SfArrowLength = X(1)/X(3)*LArrowLength;
SsArrowTheta = X(2)/X(3)*LArrowTheta;
SsArrowLength = X(2)/X(3)*LArrowLength;
f = figure('Name','Bulk waves 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
datacursormode on 
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
scatter3(0,0,0,1,'.','MarkerEdgeColor','k')
hold on
if  Phi == 0 && Theta == 0
    line([0 n(1)*1e3/X(2)],[0 n(2)*1e3/X(2)],[0 n(3)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[0 .5 1]); % slowness of S_slow
    line([0 n(1)*1e3/X(3)],[0 n(2)*1e3/X(3)],[0 n(3)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % slowness of L
    line([n(1)*1e3/X(2) (1-SsArrowLength)*cosd(Theta+SsArrowTheta)*cosd(Phi)*1e3/X(2)],[n(2)*1e3/X(2) (1-SsArrowLength)*cosd(Theta+SsArrowTheta)*sind(Phi)*1e3/X(2)],[n(3)*1e3/X(2) (1-SsArrowLength)*sind(Theta+SsArrowTheta)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[0 .5 1]); % upper arrow line of S_slow
    line([n(1)*1e3/X(2) (1-SsArrowLength)*cosd(Theta-SsArrowTheta)*cosd(Phi)*1e3/X(2)],[n(2)*1e3/X(2) (1-SsArrowLength)*cosd(Theta-SsArrowTheta)*sind(Phi)*1e3/X(2)],[n(3)*1e3/X(2) (1-SsArrowLength)*sind(Theta-SsArrowTheta)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[0 .5 1]); % lower arrow line of S_slow
    line([n(1)*1e3/X(3) (1-LArrowLength)*cosd(Theta+LArrowTheta)*cosd(Phi)*1e3/X(3)],[n(2)*1e3/X(3) (1-LArrowLength)*cosd(Theta+LArrowTheta)*sind(Phi)*1e3/X(3)],[n(3)*1e3/X(3) (1-LArrowLength)*sind(Theta+LArrowTheta)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % upper arrow line of L
    line([n(1)*1e3/X(3) (1-LArrowLength)*cosd(Theta-LArrowTheta)*cosd(Phi)*1e3/X(3)],[n(2)*1e3/X(3) (1-LArrowLength)*cosd(Theta-LArrowTheta)*sind(Phi)*1e3/X(3)],[n(3)*1e3/X(3) (1-LArrowLength)*sind(Theta-LArrowTheta)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % lower arrow line of L 
    FactorL = .125*1e3/X(2)/sqrt(p(1,3)^2+p(2,3)^2+p(3,3)^2); % vector length ratio slowness S_slow/polarization of L
    FactorSH = .125*1e3/X(2)/sqrt(p(1,2)^2+p(2,2)^2+p(3,2)^2); % vector length ratio slowness S_slow/polarization of S_slow
    line([.5*n(1)*1e3/X(2) FactorL*p(1,3)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorL*p(2,3)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorL*p(3,3)+.5*n(3)*1e3/X(2)],'color','r','linewidth',2*WaveVectorLineWidth) % polarization of L in positive direction
    line([.5*n(1)*1e3/X(2) -FactorL*p(1,3)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorL*p(2,3)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorL*p(3,3)+.5*n(3)*1e3/X(2)],'color','r','linewidth',2*WaveVectorLineWidth) % polarization of L in negative direction
    line([.5*n(1)*1e3/X(2) FactorSH*p(1,1)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorSH*p(2,1)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorSH*p(3,1)+.5*n(3)*1e3/X(2)],'color',[0 .5 1],'linewidth',2*WaveVectorLineWidth) % polarization of S_fast in positive direction
    line([.5*n(1)*1e3/X(2) -FactorSH*p(1,1)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorSH*p(2,1)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorSH*p(3,1)+.5*n(3)*1e3/X(2)],'color',[0 .5 1],'linewidth',2*WaveVectorLineWidth) % polarization of S_fast in negative direction
    line([.5*n(1)*1e3/X(2) FactorSH*p(1,2)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorSH*p(2,2)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorSH*p(3,2)+.5*n(3)*1e3/X(2)],'color',[0 .5 1],'linewidth',2*WaveVectorLineWidth) % polarization of S_slow in positive direction
    line([.5*n(1)*1e3/X(2) -FactorSH*p(1,2)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorSH*p(2,2)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorSH*p(3,2)+.5*n(3)*1e3/X(2)],'color',[0 .5 1],'linewidth',2*WaveVectorLineWidth) % polarization of S_slow in negative direction
    if  ~strcmp(Material.Class,'Isotropic')
        text(n(1)*1e3/X(2),n(2)*1e3/X(2),n(3)*1e3/X(2),'S$_1$,S$_2$','FontSize',FontSizeModeLabels,'Color','k','Interpreter','latex');
    else
        text(n(1)*1e3/X(2),n(2)*1e3/X(2),n(3)*1e3/X(2),'SV,SH','FontSize',FontSizeModeLabels,'Color','k','Interpreter','latex');
    end
elseif (~Phi == 0 || ~Theta == 0) && ~strcmp(Material.Class,'Isotropic')
    line([0 n(1)*1e3/X(2)],[0 n(2)*1e3/X(2)],[0 n(3)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]); % slowness of S_slow
    line([0 n(1)*1e3/X(1)],[0 n(2)*1e3/X(1)],[0 n(3)*1e3/X(1)],'LineWidth',WaveVectorLineWidth,'Color','b'); % slowness of S_fast
    line([0 n(1)*1e3/X(3)],[0 n(2)*1e3/X(3)],[0 n(3)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % slowness of L
    line([n(1)*1e3/X(2) (1-SsArrowLength)*cosd(Theta+SsArrowTheta)*cosd(Phi)*1e3/X(2)],[n(2)*1e3/X(2) (1-SsArrowLength)*cosd(Theta+SsArrowTheta)*sind(Phi)*1e3/X(2)],[n(3)*1e3/X(2) (1-SsArrowLength)*sind(Theta+SsArrowTheta)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]); % upper arrow line of S_slow
    line([n(1)*1e3/X(2) (1-SsArrowLength)*cosd(Theta-SsArrowTheta)*cosd(Phi)*1e3/X(2)],[n(2)*1e3/X(2) (1-SsArrowLength)*cosd(Theta-SsArrowTheta)*sind(Phi)*1e3/X(2)],[n(3)*1e3/X(2) (1-SsArrowLength)*sind(Theta-SsArrowTheta)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]); % lower arrow line of S_slow
    line([n(1)*1e3/X(1) (1-SfArrowLength)*cosd(Theta+SfArrowTheta)*cosd(Phi)*1e3/X(1)],[n(2)*1e3/X(1) (1-SfArrowLength)*cosd(Theta+SfArrowTheta)*sind(Phi)*1e3/X(1)],[n(3)*1e3/X(1) (1-SfArrowLength)*sind(Theta+SfArrowTheta)*1e3/X(1)],'LineWidth',WaveVectorLineWidth,'Color','b'); % upper arrow line of S_fast
    line([n(1)*1e3/X(1) (1-SfArrowLength)*cosd(Theta-SfArrowTheta)*cosd(Phi)*1e3/X(1)],[n(2)*1e3/X(1) (1-SfArrowLength)*cosd(Theta-SfArrowTheta)*sind(Phi)*1e3/X(1)],[n(3)*1e3/X(1) (1-SfArrowLength)*sind(Theta-SfArrowTheta)*1e3/X(1)],'LineWidth',WaveVectorLineWidth,'Color','b'); % lower arrow line of S_fast
    line([n(1)*1e3/X(3) (1-LArrowLength)*cosd(Theta+LArrowTheta)*cosd(Phi)*1e3/X(3)],[n(2)*1e3/X(3) (1-LArrowLength)*cosd(Theta+LArrowTheta)*sind(Phi)*1e3/X(3)],[n(3)*1e3/X(3) (1-LArrowLength)*sind(Theta+LArrowTheta)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % upper arrow line of L
    line([n(1)*1e3/X(3) (1-LArrowLength)*cosd(Theta-LArrowTheta)*cosd(Phi)*1e3/X(3)],[n(2)*1e3/X(3) (1-LArrowLength)*cosd(Theta-LArrowTheta)*sind(Phi)*1e3/X(3)],[n(3)*1e3/X(3) (1-LArrowLength)*sind(Theta-LArrowTheta)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % lower arrow line of L
    FactorL = .125*1e3/X(2)/sqrt(p(1,3)^2+p(2,3)^2+p(3,3)^2); % vector length ratio slowness S_slow/polarization of L
    FactorSV = .125*1e3/X(2)/sqrt(p(1,1)^2+p(2,1)^2+p(3,1)^2); % vector length ratio slowness S_slow/polarization of S_fast
    FactorSH = .125*1e3/X(2)/sqrt(p(1,2)^2+p(2,2)^2+p(3,2)^2); % vector length ratio slowness S_slow/polarization of S_slow
    line([.5*n(1)*1e3/X(2) FactorL*p(1,3)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorL*p(2,3)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorL*p(3,3)+.5*n(3)*1e3/X(2)],'color','r','linewidth',2*WaveVectorLineWidth) % polarization of L in positive direction
    line([.5*n(1)*1e3/X(2) -FactorL*p(1,3)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorL*p(2,3)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorL*p(3,3)+.5*n(3)*1e3/X(2)],'color','r','linewidth',2*WaveVectorLineWidth) % polarization of L in negative direction
    line([.5*n(1)*1e3/X(2) FactorSV*p(1,1)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorSV*p(2,1)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorSV*p(3,1)+.5*n(3)*1e3/X(2)],'color','b','linewidth',2*WaveVectorLineWidth) % polarization of S_fast in positive direction
    line([.5*n(1)*1e3/X(2) -FactorSV*p(1,1)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorSV*p(2,1)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorSV*p(3,1)+.5*n(3)*1e3/X(2)],'color','b','linewidth',2*WaveVectorLineWidth) % polarization of S_fast in negative direction
    line([.5*n(1)*1e3/X(2) FactorSH*p(1,2)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorSH*p(2,2)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorSH*p(3,2)+.5*n(3)*1e3/X(2)],'color',[.13 .55 .13],'linewidth',2*WaveVectorLineWidth) % polarization of S_slow in positive direction
    line([.5*n(1)*1e3/X(2) -FactorSH*p(1,2)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorSH*p(2,2)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorSH*p(3,2)+.5*n(3)*1e3/X(2)],'color',[.13 .55 .13],'linewidth',2*WaveVectorLineWidth) % polarization of S_slow in negative direction
    text(n(1)*1e3/X(2),n(2)*1e3/X(2),n(3)*1e3/X(2),'S$_\mathrm{slow}$','FontSize',FontSizeModeLabels,'Color','k','Interpreter','latex');
    text((1-2*SfArrowLength)*cosd(Theta-2*SfArrowTheta)*cosd(Phi)*1e3/X(1),(1-2*SfArrowLength)*cosd(Theta-2*SfArrowTheta)*sind(Phi)*1e3/X(1),(1-2*SfArrowLength)*sind(Theta-2*SfArrowTheta)*1e3/X(1),'S$_\mathrm{fast}$','FontSize',FontSizeModeLabels,'Color','k','Interpreter','latex');
elseif (~Phi == 0 || ~Theta == 0) && strcmp(Material.Class,'Isotropic')
    line([0 n(1)*1e3/X(2)],[0 n(2)*1e3/X(2)],[0 n(3)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[0 .5 1]); % slowness of S
    line([0 n(1)*1e3/X(3)],[0 n(2)*1e3/X(3)],[0 n(3)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % slowness of L
    line([n(1)*1e3/X(2) (1-SsArrowLength)*cosd(Theta+SsArrowTheta)*cosd(Phi)*1e3/X(2)],[n(2)*1e3/X(2) (1-SsArrowLength)*cosd(Theta+SsArrowTheta)*sind(Phi)*1e3/X(2)],[n(3)*1e3/X(2) (1-SsArrowLength)*sind(Theta+SsArrowTheta)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[0 .5 1]); % upper arrow line of S
    line([n(1)*1e3/X(2) (1-SsArrowLength)*cosd(Theta-SsArrowTheta)*cosd(Phi)*1e3/X(2)],[n(2)*1e3/X(2) (1-SsArrowLength)*cosd(Theta-SsArrowTheta)*sind(Phi)*1e3/X(2)],[n(3)*1e3/X(2) (1-SsArrowLength)*sind(Theta-SsArrowTheta)*1e3/X(2)],'LineWidth',WaveVectorLineWidth,'Color',[0 .5 1]); % lower arrow line of S
    line([n(1)*1e3/X(3) (1-LArrowLength)*cosd(Theta+LArrowTheta)*cosd(Phi)*1e3/X(3)],[n(2)*1e3/X(3) (1-LArrowLength)*cosd(Theta+LArrowTheta)*sind(Phi)*1e3/X(3)],[n(3)*1e3/X(3) (1-LArrowLength)*sind(Theta+LArrowTheta)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % upper arrow line of L
    line([n(1)*1e3/X(3) (1-LArrowLength)*cosd(Theta-LArrowTheta)*cosd(Phi)*1e3/X(3)],[n(2)*1e3/X(3) (1-LArrowLength)*cosd(Theta-LArrowTheta)*sind(Phi)*1e3/X(3)],[n(3)*1e3/X(3) (1-LArrowLength)*sind(Theta-LArrowTheta)*1e3/X(3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % lower arrow line of L
    FactorL = .125*1e3/X(2)/sqrt(p(1,3)^2+p(2,3)^2+p(3,3)^2); % vector length ratio slowness S_slow/polarization of L
    FactorSV = .125*1e3/X(2)/sqrt(p(1,1)^2+p(2,1)^2+p(3,1)^2); % vector length ratio slowness S_slow/polarization of S_fast
    FactorSH = .125*1e3/X(2)/sqrt(p(1,2)^2+p(2,2)^2+p(3,2)^2); % vector length ratio slowness S_slow/polarization of S_slow
    line([.5*n(1)*1e3/X(2) FactorL*p(1,3)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorL*p(2,3)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorL*p(3,3)+.5*n(3)*1e3/X(2)],'color','r','linewidth',2*WaveVectorLineWidth) % polarization of L in positive direction
    line([.5*n(1)*1e3/X(2) -FactorL*p(1,3)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorL*p(2,3)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorL*p(3,3)+.5*n(3)*1e3/X(2)],'color','r','linewidth',2*WaveVectorLineWidth) % polarization of L in negative direction
    line([.5*n(1)*1e3/X(2) FactorSV*p(1,1)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorSV*p(2,1)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorSV*p(3,1)+.5*n(3)*1e3/X(2)],'color',[0 .5 1],'linewidth',2*WaveVectorLineWidth) % polarization of S1 in positive direction
    line([.5*n(1)*1e3/X(2) -FactorSV*p(1,1)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorSV*p(2,1)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorSV*p(3,1)+.5*n(3)*1e3/X(2)],'color',[0 .5 1],'linewidth',2*WaveVectorLineWidth) % polarization of S1 in negative direction
    line([.5*n(1)*1e3/X(2) FactorSH*p(1,2)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) FactorSH*p(2,2)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) FactorSH*p(3,2)+.5*n(3)*1e3/X(2)],'color',[0 .5 1],'linewidth',2*WaveVectorLineWidth) % polarization of S2 in positive direction
    line([.5*n(1)*1e3/X(2) -FactorSH*p(1,2)+.5*n(1)*1e3/X(2)],[.5*n(2)*1e3/X(2) -FactorSH*p(2,2)+.5*n(2)*1e3/X(2)],[.5*n(3)*1e3/X(2) -FactorSH*p(3,2)+.5*n(3)*1e3/X(2)],'color',[0 .5 1],'linewidth',2*WaveVectorLineWidth) % polarization of S2 in negative direction
    text(n(1)*1e3/X(2),n(2)*1e3/X(2),n(3)*1e3/X(2),'SV,SH','FontSize',FontSizeModeLabels,'Color','k','Interpreter','latex');
end
text((1-2*LArrowLength)*cosd(Theta-2*LArrowTheta)*cosd(Phi)*1e3/X(3),(1-2*LArrowLength)*cosd(Theta-2*LArrowTheta)*sind(Phi)*1e3/X(3),(1-2*LArrowLength)*sind(Theta-2*LArrowTheta)*1e3/X(3),'L','FontSize',FontSizeModeLabels,'Color','k','Interpreter','latex');axis equal
ax = gca;
ax.Clipping = 'off';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
if  HeadLine == 1
    ax.Title.Interpreter = 'latex';
    ax.Title.FontSize = FontSizeHeadLine;
    ax.Title.String = ['Bulk waves propagating @ $\phi$ = ',num2str(Phi),'\,$^{\circ}$, $\theta$ = ',num2str(Theta),'\,$^{\circ}$ in ',char(join(split(Material.Name,'_'),'\_'))];
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
ax.XLabel.String = '$\frac{\zeta^\prime_1}{\omega}$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.YLabel.String = '$\frac{\zeta^\prime_2}{\omega}$';
ax.ZLabel.Interpreter = 'latex';
ax.ZLabel.FontSize = FontSizeAxesLabels;
ax.ZLabel.String = '$\frac{\zeta^\prime_3}{\omega}$ (ms/m)';
ax.TickLabelInterpreter = 'latex';
view(ViewPhi,ViewTheta)
if  Export == 1
    try
        if  PDF == 1
            if  Crop == 0
                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[50 35])
                print(f,fullfile(Directory,FileName),'-dpdf','-vector')
            elseif Crop == 1
                exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
            end
        end
        if  PNG == 1
            if  Crop == 0
                print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
            elseif Crop == 1
                exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
            end
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
        return
    end
end
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj)
    if  event_obj.Target.LineWidth == WaveVectorLineWidth
        output_txt = {['$s^\prime$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'} ms/m'],['$s^\prime_1$: \textbf{',num2str(event_obj.Position(1),6),'} ms/m'],['$s^\prime_2$: \textbf{',num2str(event_obj.Position(2),6),'} ms/m'],['$s^\prime_3$: \textbf{',num2str(event_obj.Position(3),6),'} ms/m']};
    else
        if  all(event_obj.Target.Color == [1 0 0])
            output_txt = {['$\beta$: \textbf{',num2str(Beta(3)),'}\,$^{\circ}$']};
        elseif all(event_obj.Target.Color == [0 0 1])
            output_txt = {['$\beta$: \textbf{',num2str(180-Beta(1)),'}\,$^{\circ}$']};
        else
            output_txt = {['$\beta$: \textbf{',num2str(Beta(2)),'}\,$^{\circ}$']};
        end
    end
end
end