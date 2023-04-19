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
function Slowness_BulkWave2D(Crop,Material,Plane,Phi,ThetaStep,Export,PDF,PNG,FileName,Directory,HeadLine,LineWidth,FontSizeHeadLine,FontSizeAxesLabels,FontSizeAxes,FontSizeModeLabels,PNGresolution)
Theta = 0:ThetaStep:90;
n = [cosd(Theta)'.*cosd(Phi)' cosd(Theta)'.*sind(Phi)' sind(Theta)']; % propagation direction vector matrix for the bulk waves
C = real(Material.C);
X(length(n),4) = 0;
if  ~strcmp(Material.Class,'Isotropic')
    for i = 1:length(n) % solving the Christoffel equation for the phase velocities X of the bulk waves in terms of the propagation direction vector matrix n [1], p. 41-42: the determinant is given as polynomial v^6+A1*v^4+A2*v^2+A3 = 0; the roots X are given by analytical expressions            
        if  Plane == 13
            A1 = -(C(1,1)*n(i,1)^2+C(2,2)*n(i,2)^2+C(3,3)*n(i,3)^2+C(4,4)*n(i,2)^2+C(4,4)*n(i,3)^2+C(5,5)*n(i,1)^2+C(5,5)*n(i,3)^2+C(6,6)*n(i,1)^2+C(6,6)*n(i,2)^2)/Material.Density;
            A2 = ((C(1,1)*C(2,2)+C(1,1)*C(4,4)+C(2,2)*C(5,5)-2*C(1,2)*C(6,6)+C(4,4)*C(6,6)+C(5,5)*C(6,6)-C(1,2)^2)*n(i,1)^2*n(i,2)^2+(C(1,1)*C(3,3)+C(1,1)*C(4,4)-2*C(1,3)*C(5,5)+C(3,3)*C(6,6)+C(4,4)*C(5,5)+C(5,5)*C(6,6)-C(1,3)^2)*n(i,1)^2*n(i,3)^2+(C(2,2)*C(3,3)-2*C(2,3)*C(4,4)+C(2,2)*C(5,5)+C(3,3)*C(6,6)+C(4,4)*C(5,5)+C(4,4)*C(6,6)-C(2,3)^2)*n(i,2)^2*n(i,3)^2+(C(1,1)*C(5,5)+C(1,1)*C(6,6)+C(5,5)*C(6,6))*n(i,1)^4+(C(2,2)*C(4,4)+C(2,2)*C(6,6)+C(4,4)*C(6,6))*n(i,2)^4+(C(3,3)*C(4,4)+C(3,3)*C(5,5)+C(4,4)*C(5,5))*n(i,3)^4)/Material.Density^2;
            A3 = ((C(1,1)*C(2,3)^2+C(1,3)^2*C(2,2)+C(1,2)^2*C(3,3)-2*C(1,2)*C(1,3)*C(2,3)-C(1,1)*C(2,2)*C(3,3)-2*C(1,2)*C(1,3)*C(4,4)+2*C(1,1)*C(2,3)*C(4,4)-2*C(1,2)*C(2,3)*C(5,5)+2*C(1,3)*C(2,2)*C(5,5)-2*C(1,3)*C(2,3)*C(6,6)+2*C(1,2)*C(3,3)*C(6,6)-2*C(1,2)*C(4,4)*C(5,5)-2*C(1,3)*C(4,4)*C(6,6)-2*C(2,3)*C(5,5)*C(6,6)-4*C(4,4)*C(5,5)*C(6,6))*n(i,1)^2*n(i,2)^2*n(i,3)^2+(C(1,2)^2*C(4,4)-C(1,1)*C(2,2)*C(4,4)+2*C(1,2)*C(4,4)*C(6,6)-C(2,2)*C(5,5)*C(6,6))*n(i,1)^2*n(i,2)^4+(C(1,3)^2*C(4,4)-C(1,1)*C(3,3)*C(4,4)+2*C(1,3)*C(4,4)*C(5,5)-C(3,3)*C(5,5)*C(6,6))*n(i,1)^2*n(i,3)^4+(C(1,2)^2*C(5,5)-C(1,1)*C(2,2)*C(5,5)-C(1,1)*C(4,4)*C(6,6)+2*C(1,2)*C(5,5)*C(6,6))*n(i,1)^4*n(i,2)^2+(C(1,3)^2*C(6,6)-C(1,1)*C(3,3)*C(6,6)-C(1,1)*C(4,4)*C(5,5)+2*C(1,3)*C(5,5)*C(6,6))*n(i,1)^4*n(i,3)^2+(C(2,3)^2*C(5,5)-C(2,2)*C(3,3)*C(5,5)+2*C(2,3)*C(4,4)*C(5,5)-C(3,3)*C(4,4)*C(6,6))*n(i,2)^2*n(i,3)^4+(C(2,3)^2*C(6,6)-C(2,2)*C(3,3)*C(6,6)-C(2,2)*C(4,4)*C(5,5)+2*C(2,3)*C(4,4)*C(6,6))*n(i,2)^4*n(i,3)^2-C(1,1)*C(5,5)*C(6,6)*n(i,1)^6-C(2,2)*C(4,4)*C(6,6)*n(i,2)^6-C(3,3)*C(4,4)*C(5,5)*n(i,3)^6)/Material.Density^3;
        elseif Plane == 12 % swap C13 - C12, C22 - C33, C55 - C66
            A1 = -(C(1,1)*n(i,1)^2+C(3,3)*n(i,2)^2+C(2,2)*n(i,3)^2+C(4,4)*n(i,2)^2+C(4,4)*n(i,3)^2+C(6,6)*n(i,1)^2+C(6,6)*n(i,3)^2+C(5,5)*n(i,1)^2+C(5,5)*n(i,2)^2)/Material.Density;
            A2 = ((C(1,1)*C(3,3)+C(1,1)*C(4,4)+C(3,3)*C(6,6)-2*C(1,3)*C(5,5)+C(4,4)*C(5,5)+C(6,6)*C(5,5)-C(1,3)^2)*n(i,1)^2*n(i,2)^2+(C(1,1)*C(2,2)+C(1,1)*C(4,4)-2*C(1,2)*C(6,6)+C(2,2)*C(5,5)+C(4,4)*C(6,6)+C(6,6)*C(5,5)-C(1,2)^2)*n(i,1)^2*n(i,3)^2+(C(3,3)*C(2,2)-2*C(2,3)*C(4,4)+C(3,3)*C(6,6)+C(2,2)*C(5,5)+C(4,4)*C(6,6)+C(4,4)*C(5,5)-C(2,3)^2)*n(i,2)^2*n(i,3)^2+(C(1,1)*C(6,6)+C(1,1)*C(5,5)+C(6,6)*C(5,5))*n(i,1)^4+(C(3,3)*C(4,4)+C(3,3)*C(5,5)+C(4,4)*C(5,5))*n(i,2)^4+(C(2,2)*C(4,4)+C(2,2)*C(6,6)+C(4,4)*C(6,6))*n(i,3)^4)/Material.Density^2;
            A3 = ((C(1,1)*C(2,3)^2+C(1,2)^2*C(3,3)+C(1,3)^2*C(2,2)-2*C(1,3)*C(1,2)*C(2,3)-C(1,1)*C(3,3)*C(2,2)-2*C(1,3)*C(1,2)*C(4,4)+2*C(1,1)*C(2,3)*C(4,4)-2*C(1,3)*C(2,3)*C(6,6)+2*C(1,2)*C(3,3)*C(6,6)-2*C(1,2)*C(2,3)*C(5,5)+2*C(1,3)*C(2,2)*C(5,5)-2*C(1,3)*C(4,4)*C(6,6)-2*C(1,2)*C(4,4)*C(5,5)-2*C(2,3)*C(6,6)*C(5,5)-4*C(4,4)*C(6,6)*C(5,5))*n(i,1)^2*n(i,2)^2*n(i,3)^2+(C(1,3)^2*C(4,4)-C(1,1)*C(3,3)*C(4,4)+2*C(1,3)*C(4,4)*C(5,5)-C(3,3)*C(6,6)*C(5,5))*n(i,1)^2*n(i,2)^4+(C(1,2)^2*C(4,4)-C(1,1)*C(2,2)*C(4,4)+2*C(1,2)*C(4,4)*C(6,6)-C(2,2)*C(6,6)*C(5,5))*n(i,1)^2*n(i,3)^4+(C(1,3)^2*C(6,6)-C(1,1)*C(3,3)*C(6,6)-C(1,1)*C(4,4)*C(5,5)+2*C(1,3)*C(6,6)*C(5,5))*n(i,1)^4*n(i,2)^2+(C(1,2)^2*C(5,5)-C(1,1)*C(2,2)*C(5,5)-C(1,1)*C(4,4)*C(6,6)+2*C(1,2)*C(6,6)*C(5,5))*n(i,1)^4*n(i,3)^2+(C(2,3)^2*C(6,6)-C(3,3)*C(2,2)*C(6,6)+2*C(2,3)*C(4,4)*C(6,6)-C(2,2)*C(4,4)*C(5,5))*n(i,2)^2*n(i,3)^4+(C(2,3)^2*C(5,5)-C(3,3)*C(2,2)*C(5,5)-C(3,3)*C(4,4)*C(6,6)+2*C(2,3)*C(4,4)*C(5,5))*n(i,2)^4*n(i,3)^2-C(1,1)*C(6,6)*C(5,5)*n(i,1)^6-C(3,3)*C(4,4)*C(5,5)*n(i,2)^6-C(2,2)*C(4,4)*C(6,6)*n(i,3)^6)/Material.Density^3;
        end
        Xa = A2/3-A1^2/9;
        Xb = A1^3/27-A1*A2/6+A3/2;
        Xc = (sqrt(Xb^2+Xa^3)-Xb)^(1/3);
        Xd = Xa/(2*Xc)-Xc/2;
        Xe = Xa/Xc;
        Xf = (sqrt(3)*(Xc+Xe)*1i)/2;
        X(i,2) = sqrt(Xd-Xf-A1/3);
        X(i,3) = sqrt(Xd+Xf-A1/3);
        X(i,4) = -sqrt(Xc-Xe-A1/3);
    end
else
    X(:,2:3) = Material.TransverseVelocity;
    X(:,4) = Material.LongitudinalVelocity;
end    
X = abs(real(X)/1e3); % phase velocity (m/ms)
X(:,5:7) = 1./X(:,2:4); % slowness (ms/m)
X = vertcat(X,flipud(X(1:end-1,:)),X(2:end,:),flipud(X(1:end-1,:))); % extending the data to 360 deg
X(:,1) = 0:2*pi/(4*(length(n)-1)):2*pi; % Theta (rad)
f = figure('Name','Slowness 2-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 .6 1],'color','w');
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
polarplot(X(:,1),X(:,5),'LineWidth',LineWidth,'Color','b')
z = find(abs(X(:,1)-pi/4) == min(abs(X(:,1)-pi/4)));
if  ~strcmp(Material.Class,'Isotropic')
    text(pi/4,1.04*X(z(1),5),'S$_\mathrm{fast}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
end
hold on
polarplot(X(:,1),X(:,6),'LineWidth',LineWidth,'Color',[.13 .55 .13])
z = find(abs(X(:,1)-pi/8) == min(abs(X(:,1)-pi/8)));
if  ~strcmp(Material.Class,'Isotropic')
    text(pi/8,1.04*X(z(1),6),'S$_\mathrm{slow}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
else
    text(pi/8,1.04*X(z(1),6),'SV,SH','FontSize',FontSizeModeLabels,'Interpreter','latex');
end
polarplot(X(:,1),X(:,7),'LineWidth',LineWidth,'Color','r')
z = find(abs(X(:,1)+pi) == min(abs(X(:,1)+pi)));
text(pi,.95*X(z(1),7),'L','FontSize',FontSizeModeLabels,'Interpreter','latex');    
ax = gca;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
if  HeadLine == 0
    ax.Title.FontSize = FontSizeAxesLabels;
    ax.Title.String = 'Slowness (ms/m)';
elseif HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    if  Plane == 13
        ax.Title.String = ['Bulk wave slowness @ $\phi$ = ',num2str(Phi),'\,$^{\circ}$ in ',char(join(split(Material.Name,'_'),'\_'))];
    elseif Plane == 12
        ax.Title.String = ['Bulk wave slowness @ $\theta$ = ',num2str(Phi),'\,$^{\circ}$ in ',char(join(split(Material.Name,'_'),'\_'))];
    end
end
z = ax.RTick;
z(1) = '';
ax.RTick = z;
ax.ThetaTick = [0 45 90 135 180 225 270 315];
ax.ThetaAxis.TickLabelFormat = '%g';
ax.TickLabelInterpreter = 'latex';
if  HeadLine > 0
    text(.59,1.04,'(ms/m)','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
end
if  Plane == 13
    text(.41,1.03,'$s_3$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
elseif Plane == 12
    text(.41,1.03,'$s_2$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
end
text(1,.45,'$s_1$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')  
if  Export == 1
    try
        if  PDF == 1
            if  Crop == 0
                set(f,'PaperUnits','centimeters','PaperSize',[30 35])
                print(f,fullfile(Directory,FileName),'-dpdf');
            elseif Crop == 1
                exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
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
    if  Plane == 13
        output_txt = {['$\theta$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$s$: \textbf{',num2str(event_obj.Position(2),6),'} ms/m']};
    elseif Plane == 12
        output_txt = {['$\phi$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$s$: \textbf{',num2str(event_obj.Position(2),6),'} ms/m']};        
    end
end 
end