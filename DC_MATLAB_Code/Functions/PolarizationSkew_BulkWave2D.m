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
function PolarizationSkew_BulkWave2D(Material,Plane,Phi,ThetaStep,Export,PDF,PNG,FileName,Directory,HeadLine,LineWidth,FontSizeHeadLine,FontSizeAxesLabels,FontSizeAxes,FontSizeModeLabels,PNGresolution)
%#ok<*AGROW>
if  strcmp(Material.Class,'Isotropic')
    msgbox('The polarization of bulk waves in isotropic media is purely longitudinal for L waves and purely perpendicular for SV and SH for any propgation direction. Therefore, the polarization skew angle is always zero.','Info');
    return
end
Theta = 0:ThetaStep:90;
Theta(1) = .001;
Theta(end) = 89.999;
if  Phi == 0 && strcmp(Material.Class,'Cubic')
    n = [cosd(Theta)'.*cosd(.001)' cosd(Theta)'.*sind(.001)' sind(Theta)']; % propagation direction vector matrix for the bulk waves
elseif Phi == 45 && strcmp(Material.Class,'Cubic')
    n = [cosd(Theta)'.*cosd(45.001)' cosd(Theta)'.*sind(45.001)' sind(Theta)'];
elseif Phi == 90 && strcmp(Material.Class,'Cubic')
    n = [cosd(Theta)'.*cosd(89.999)' cosd(Theta)'.*sind(89.999)' sind(Theta)'];
else
    n = [cosd(Theta)'.*cosd(Phi)' cosd(Theta)'.*sind(Phi)' sind(Theta)'];
end
C = real(Material.C);
for i = 1:length(n) % solving the Christoffel equation for the phase velocities X of the bulk waves in terms of the propagation direction vector matrix n [1], p. 41-42: the determinant is given as polynomial v^6+A1*v^4+A2*v^2+A3 = 0; the roots X are given by analytical expressions                
    if  Plane == 13
        A1 = -(C(1,1)*n(i,1)^2+C(2,2)*n(i,2)^2+C(3,3)*n(i,3)^2+C(4,4)*n(i,2)^2+C(4,4)*n(i,3)^2+C(5,5)*n(i,1)^2+C(5,5)*n(i,3)^2+C(6,6)*n(i,1)^2+C(6,6)*n(i,2)^2)/Material.Density;
        A2 = ((C(1,1)*C(2,2)+C(1,1)*C(4,4)+C(2,2)*C(5,5)-2*C(1,2)*C(6,6)+C(4,4)*C(6,6)+C(5,5)*C(6,6)-C(1,2)^2)*n(i,1)^2*n(i,2)^2+(C(1,1)*C(3,3)+C(1,1)*C(4,4)-2*C(1,3)*C(5,5)+C(3,3)*C(6,6)+C(4,4)*C(5,5)+C(5,5)*C(6,6)-C(1,3)^2)*n(i,1)^2*n(i,3)^2+(C(2,2)*C(3,3)-2*C(2,3)*C(4,4)+C(2,2)*C(5,5)+C(3,3)*C(6,6)+C(4,4)*C(5,5)+C(4,4)*C(6,6)-C(2,3)^2)*n(i,2)^2*n(i,3)^2+(C(1,1)*C(5,5)+C(1,1)*C(6,6)+C(5,5)*C(6,6))*n(i,1)^4+(C(2,2)*C(4,4)+C(2,2)*C(6,6)+C(4,4)*C(6,6))*n(i,2)^4+(C(3,3)*C(4,4)+C(3,3)*C(5,5)+C(4,4)*C(5,5))*n(i,3)^4)/Material.Density^2;
        A3 = ((C(1,1)*C(2,3)^2+C(1,3)^2*C(2,2)+C(1,2)^2*C(3,3)-2*C(1,2)*C(1,3)*C(2,3)-C(1,1)*C(2,2)*C(3,3)-2*C(1,2)*C(1,3)*C(4,4)+2*C(1,1)*C(2,3)*C(4,4)-2*C(1,2)*C(2,3)*C(5,5)+2*C(1,3)*C(2,2)*C(5,5)-2*C(1,3)*C(2,3)*C(6,6)+2*C(1,2)*C(3,3)*C(6,6)-2*C(1,2)*C(4,4)*C(5,5)-2*C(1,3)*C(4,4)*C(6,6)-2*C(2,3)*C(5,5)*C(6,6)-4*C(4,4)*C(5,5)*C(6,6))*n(i,1)^2*n(i,2)^2*n(i,3)^2+(C(1,2)^2*C(4,4)-C(1,1)*C(2,2)*C(4,4)+2*C(1,2)*C(4,4)*C(6,6)-C(2,2)*C(5,5)*C(6,6))*n(i,1)^2*n(i,2)^4+(C(1,3)^2*C(4,4)-C(1,1)*C(3,3)*C(4,4)+2*C(1,3)*C(4,4)*C(5,5)-C(3,3)*C(5,5)*C(6,6))*n(i,1)^2*n(i,3)^4+(C(1,2)^2*C(5,5)-C(1,1)*C(2,2)*C(5,5)-C(1,1)*C(4,4)*C(6,6)+2*C(1,2)*C(5,5)*C(6,6))*n(i,1)^4*n(i,2)^2+(C(1,3)^2*C(6,6)-C(1,1)*C(3,3)*C(6,6)-C(1,1)*C(4,4)*C(5,5)+2*C(1,3)*C(5,5)*C(6,6))*n(i,1)^4*n(i,3)^2+(C(2,3)^2*C(5,5)-C(2,2)*C(3,3)*C(5,5)+2*C(2,3)*C(4,4)*C(5,5)-C(3,3)*C(4,4)*C(6,6))*n(i,2)^2*n(i,3)^4+(C(2,3)^2*C(6,6)-C(2,2)*C(3,3)*C(6,6)-C(2,2)*C(4,4)*C(5,5)+2*C(2,3)*C(4,4)*C(6,6))*n(i,2)^4*n(i,3)^2-C(1,1)*C(5,5)*C(6,6)*n(i,1)^6-C(2,2)*C(4,4)*C(6,6)*n(i,2)^6-C(3,3)*C(4,4)*C(5,5)*n(i,3)^6)/Material.Density^3;
        B = [C(1,1)*n(i,1)^2+C(6,6)*n(i,2)^2+C(5,5)*n(i,3)^2 (C(1,2)+C(6,6))*n(i,1)*n(i,2) (C(1,3)+C(5,5))*n(i,1)*n(i,3);0 C(6,6)*n(i,1)^2+C(2,2)*n(i,2)^2+C(4,4)*n(i,3)^2 (C(2,3)+C(4,4))*n(i,2)*n(i,3);0 0 C(5,5)*n(i,1)^2+C(4,4)*n(i,2)^2+C(3,3)*n(i,3)^2]/Material.Density;
    elseif Plane == 12 % swap C13 - C12, C22 - C33, C55 - C66
        A1 = -(C(1,1)*n(i,1)^2+C(3,3)*n(i,2)^2+C(2,2)*n(i,3)^2+C(4,4)*n(i,2)^2+C(4,4)*n(i,3)^2+C(6,6)*n(i,1)^2+C(6,6)*n(i,3)^2+C(5,5)*n(i,1)^2+C(5,5)*n(i,2)^2)/Material.Density;
        A2 = ((C(1,1)*C(3,3)+C(1,1)*C(4,4)+C(3,3)*C(6,6)-2*C(1,3)*C(5,5)+C(4,4)*C(5,5)+C(6,6)*C(5,5)-C(1,3)^2)*n(i,1)^2*n(i,2)^2+(C(1,1)*C(2,2)+C(1,1)*C(4,4)-2*C(1,2)*C(6,6)+C(2,2)*C(5,5)+C(4,4)*C(6,6)+C(6,6)*C(5,5)-C(1,2)^2)*n(i,1)^2*n(i,3)^2+(C(3,3)*C(2,2)-2*C(2,3)*C(4,4)+C(3,3)*C(6,6)+C(2,2)*C(5,5)+C(4,4)*C(6,6)+C(4,4)*C(5,5)-C(2,3)^2)*n(i,2)^2*n(i,3)^2+(C(1,1)*C(6,6)+C(1,1)*C(5,5)+C(6,6)*C(5,5))*n(i,1)^4+(C(3,3)*C(4,4)+C(3,3)*C(5,5)+C(4,4)*C(5,5))*n(i,2)^4+(C(2,2)*C(4,4)+C(2,2)*C(6,6)+C(4,4)*C(6,6))*n(i,3)^4)/Material.Density^2;
        A3 = ((C(1,1)*C(2,3)^2+C(1,2)^2*C(3,3)+C(1,3)^2*C(2,2)-2*C(1,3)*C(1,2)*C(2,3)-C(1,1)*C(3,3)*C(2,2)-2*C(1,3)*C(1,2)*C(4,4)+2*C(1,1)*C(2,3)*C(4,4)-2*C(1,3)*C(2,3)*C(6,6)+2*C(1,2)*C(3,3)*C(6,6)-2*C(1,2)*C(2,3)*C(5,5)+2*C(1,3)*C(2,2)*C(5,5)-2*C(1,3)*C(4,4)*C(6,6)-2*C(1,2)*C(4,4)*C(5,5)-2*C(2,3)*C(6,6)*C(5,5)-4*C(4,4)*C(6,6)*C(5,5))*n(i,1)^2*n(i,2)^2*n(i,3)^2+(C(1,3)^2*C(4,4)-C(1,1)*C(3,3)*C(4,4)+2*C(1,3)*C(4,4)*C(5,5)-C(3,3)*C(6,6)*C(5,5))*n(i,1)^2*n(i,2)^4+(C(1,2)^2*C(4,4)-C(1,1)*C(2,2)*C(4,4)+2*C(1,2)*C(4,4)*C(6,6)-C(2,2)*C(6,6)*C(5,5))*n(i,1)^2*n(i,3)^4+(C(1,3)^2*C(6,6)-C(1,1)*C(3,3)*C(6,6)-C(1,1)*C(4,4)*C(5,5)+2*C(1,3)*C(6,6)*C(5,5))*n(i,1)^4*n(i,2)^2+(C(1,2)^2*C(5,5)-C(1,1)*C(2,2)*C(5,5)-C(1,1)*C(4,4)*C(6,6)+2*C(1,2)*C(6,6)*C(5,5))*n(i,1)^4*n(i,3)^2+(C(2,3)^2*C(6,6)-C(3,3)*C(2,2)*C(6,6)+2*C(2,3)*C(4,4)*C(6,6)-C(2,2)*C(4,4)*C(5,5))*n(i,2)^2*n(i,3)^4+(C(2,3)^2*C(5,5)-C(3,3)*C(2,2)*C(5,5)-C(3,3)*C(4,4)*C(6,6)+2*C(2,3)*C(4,4)*C(5,5))*n(i,2)^4*n(i,3)^2-C(1,1)*C(6,6)*C(5,5)*n(i,1)^6-C(3,3)*C(4,4)*C(5,5)*n(i,2)^6-C(2,2)*C(4,4)*C(6,6)*n(i,3)^6)/Material.Density^3;
        B = [C(1,1)*n(i,1)^2+C(5,5)*n(i,2)^2+C(6,6)*n(i,3)^2 (C(1,3)+C(5,5))*n(i,1)*n(i,2) (C(1,2)+C(6,6))*n(i,1)*n(i,3);0 C(5,5)*n(i,1)^2+C(3,3)*n(i,2)^2+C(4,4)*n(i,3)^2 (C(2,3)+C(4,4))*n(i,2)*n(i,3);0 0 C(6,6)*n(i,1)^2+C(4,4)*n(i,2)^2+C(2,2)*n(i,3)^2]/Material.Density;
    end
    Xa = A2/3-A1^2/9;
    Xb = A1^3/27-A1*A2/6+A3/2;
    Xc = (sqrt(Xb^2+Xa^3)-Xb)^(1/3);
    Xd = Xa/(2*Xc)-Xc/2;
    Xe = Xa/Xc;
    Xf = (sqrt(3)*(Xc+Xe)*1i)/2;
    X(1) = abs(real(sqrt(Xd-Xf-A1/3)));
    X(2) = abs(real(sqrt(Xd+Xf-A1/3)));
    X(3) = abs(real(-sqrt(Xc-Xe-A1/3)));
    for j = 1:3 % calculating the polarizations, j=1: S_fast, j=2: S_slow, j=3: L
        p(1,j) = 1; % bulk wave amplitude U1, i.e., polarization component along x1; U1 = 1
        p(2,j) = (B(2,3)*(B(1,1)-X(j)^2)-B(1,3)*B(1,2))/(B(1,3)*(B(2,2)-X(j)^2)-B(1,2)*B(2,3)); % bulk wave amplitude ratio U2/U1, i.e., polarization component along x2
        p(3,j)  = ((B(1,1)-X(j)^2)*(B(2,2)-X(j)^2)-B(1,2)^2)/(B(1,2)*B(2,3)-B(1,3)*(B(2,2)-X(j)^2)); % bulk wave amplitude ratio U3/U1, i.e., polarization component along x3
        if  strcmp(Material.Class,'Cubic') || strcmp(Material.Class,'Transversely isotropic') && j ~= 2 || strcmp(Material.Class,'Orthotropic')
            p(1,j) = 1; % bulk wave amplitude U1, i.e., polarization component along x1; U1 = 1
            p(2,j) = (B(2,3)*(B(1,1)-(X(j))^2)-B(1,3)*B(1,2))/(B(1,3)*(B(2,2)-(X(j))^2)-B(1,2)*B(2,3)); % bulk wave amplitude ratio U2/U1, i.e., polarization component along x2
            p(3,j) = ((B(1,1)-(X(j))^2)*(B(2,2)-(X(j))^2)-B(1,2)^2)/(B(1,2)*B(2,3)-B(1,3)*(B(2,2)-(X(j))^2)); % bulk wave amplitude ratio U3/U1, i.e., polarization component along x3
        else
            p(1,j) = 0;
            p(2,j) = 1; % S_slow is always polarized normal to the fibers, i.e., in the x2-x3-plane only
            p(3,j) = -(B(2,2)-X(j)^2)/B(2,3); % bulk wave amplitude ratio U3/U2, i.e., polarization component along x3
        end
        if  Phi == 0 && j == 2 && ~strcmp(Material.Class,'Cubic')
            Delta(i,j+1) = 0;
        elseif Phi == 90 && ~strcmp(Material.Class,'Cubic')
            Delta(i,j+1) = 0;
        else
            if  j < 3
                Delta(i,j+1) = asind(dot(n(i,:),p(:,j)/norm(p(:,j)))); % polarization skew angle of S_fast,S_slow
            else
                Delta(i,j+1) = acosd(dot(n(i,:),p(:,j)/norm(p(:,j)))); % polarization skew angle of L
            end
        end
    end
end
Delta = vertcat(Delta,flipud(Delta(1:end-1,:)),Delta(2:end,:),flipud(Delta(1:end-1,:))); % extending the data to 360 deg
Delta(:,1) = 0:2*pi/(4*(length(n)-1)):2*pi; % Theta (rad)
f = figure('Name','Polarization 2-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 .6 1],'color','w');
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
polarplot(Delta(:,1),Delta(:,4),'LineWidth',LineWidth,'Color','r')
[~,z] = max(Delta(:,4));
text(Delta(z(1)),1.04*Delta(z(1),4),'L','FontSize',FontSizeModeLabels,'Interpreter','latex');
hold on
polarplot(Delta(:,1),abs(Delta(:,2)),'LineWidth',LineWidth,'Color','b')
[~,z] = min(Delta(:,2));
text(pi-Delta(z(1)),1.04*Delta(z(1),2),'S$_\mathrm{fast}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
polarplot(Delta(:,1),abs(Delta(:,3)),'LineWidth',LineWidth,'Color',[.13 .55 .13])
[~,z] = min(Delta(:,3));
text(Delta(z(1)),1.12*Delta(z(1),3),'S$_\mathrm{slow}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
ax = gca;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
if  HeadLine == 0
    ax.Title.FontSize = FontSizeAxesLabels;
    ax.Title.String = 'Polarization skew angle ($^\circ$)';
elseif HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    if  Plane == 13
        ax.Title.String = ['Polarization skew angle @ $\phi$ = ',num2str(Phi),'\,$^{\circ}$ in ',replace(Material.Name,'_','\_')];
    elseif Plane == 12
        ax.Title.String = ['Polarization skew angle @ $\theta$ = ',num2str(Phi),'\,$^{\circ}$ in ',replace(Material.Name,'_','\_')];
    end
end
z = ax.RTick;
z(1) = '';
ax.RTick = z;
ax.ThetaTick = [0 45 90 135 180 225 270 315];
ax.ThetaAxis.TickLabelFormat = '%g';
ax.TickLabelInterpreter = 'latex';
if  HeadLine > 0
    text(.59,1.04,'($^\circ$)','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
end
if  Plane == 13
    text(.41,1.03,'$x_3$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
elseif Plane == 12
    text(.41,1.03,'$x_2$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
end
text(1,.45,'$x_1$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
if  Export
    try
        if  PDF
            exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
        end
        if  PNG
            exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
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
    if  Plane == 13
        output_txt = {['$\theta$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$\delta$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
    elseif Plane == 12
        output_txt = {['$\phi$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$\delta$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};        
    end
end
end