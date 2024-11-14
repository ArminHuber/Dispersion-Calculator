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
function GroupVelocity_BulkWave2D(Material,Plane,Phi,ThetaStep,Export,PDF,PNG,FileName,Directory,HeadLine,LineWidth,FontSizeHeadLine,FontSizeAxesLabels,FontSizeAxes,FontSizeModeLabels,PNGresolution)
%#ok<*AGROW>
Theta = 0:ThetaStep:90;
Theta(1) = .001;
Theta(end) = 89.999;
if  Phi == 0
    n = [cosd(Theta)'.*cosd(.001)' cosd(Theta)'.*sind(.001)' sind(Theta)']; % propagation direction vector matrix for the bulk waves
elseif Phi == 45 && strcmp(Material.Class,'Cubic')
    n = [cosd(Theta)'.*cosd(45.001)' cosd(Theta)'.*sind(45.001)' sind(Theta)'];
elseif Phi == 90
    n = [cosd(Theta)'.*cosd(89.999)' cosd(Theta)'.*sind(89.999)' sind(Theta)'];
else
    n = [cosd(Theta)'.*cosd(Phi)' cosd(Theta)'.*sind(Phi)' sind(Theta)'];
end
C = real(Material.C);
G(length(n),4) = 0;
if  ~strcmp(Material.Class,'Isotropic')
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
            if  strcmp(Material.Class,'Cubic') || strcmp(Material.Class,'Transversely isotropic') && j ~= 2 || strcmp(Material.Class,'Orthotropic')
                p(1,j) = 1; % bulk wave amplitude U1, i.e., polarization component along x1; U1 = 1
                p(2,j) = (B(2,3)*(B(1,1)-(X(j))^2)-B(1,3)*B(1,2))/(B(1,3)*(B(2,2)-(X(j))^2)-B(1,2)*B(2,3)); % bulk wave amplitude ratio U2/U1, i.e., polarization component along x2
                p(3,j) = ((B(1,1)-(X(j))^2)*(B(2,2)-(X(j))^2)-B(1,2)^2)/(B(1,2)*B(2,3)-B(1,3)*(B(2,2)-(X(j))^2)); % bulk wave amplitude ratio U3/U1, i.e., polarization component along x3
            else
                p(1,j) = 0;
                p(2,j) = 1; % S_slow is always polarized normal to the fibers, i.e., in the x2-x3-plane only
                p(3,j) = -(B(2,2)-X(j)^2)/B(2,3); % bulk wave amplitude ratio U3/U2, i.e., polarization component along x3
            end
            p(:,j) = p(:,j)/norm(p(:,j)); % unit vector
            if  Plane == 13
                g(1,j) = (C(1,1)*n(i,1)*p(1,j)^2+C(1,2)*n(i,2)*p(2,j)*p(1,j)+C(1,3)*n(i,3)*p(3,j)*p(1,j)+C(6,6)*(n(i,1)*p(2,j)^2+n(i,2)*p(1,j)*p(2,j))+C(5,5)*(n(i,1)*p(3,j)^2+n(i,3)*p(1,j)*p(3,j)))/(Material.Density*X(j)*1e3); % group velocity component along x1
                g(2,j) = (C(6,6)*(n(i,1)*p(2,j)*p(1,j)+n(i,2)*p(1,j)^2)+C(1,2)*n(i,1)*p(1,j)*p(2,j)+C(2,2)*n(i,2)*p(2,j)^2+C(2,3)*n(i,3)*p(3,j)*p(2,j)+C(4,4)*(n(i,2)*p(3,j)^2+n(i,3)*p(2,j)*p(3,j)))/(Material.Density*X(j)*1e3); % group velocity component along x2
                g(3,j) = (C(5,5)*(n(i,1)*p(3,j)*p(1,j)+n(i,3)*p(1,j)^2)+C(4,4)*(n(i,2)*p(3,j)*p(2,j)+n(i,3)*p(2,j)^2)+C(1,3)*n(i,1)*p(1,j)*p(3,j)+C(2,3)*n(i,2)*p(2,j)*p(3,j)+C(3,3)*n(i,3)*p(3,j)^2)/(Material.Density*X(j)*1e3); % group velocity component along x3
            elseif Plane == 12
                g(1,j) = (C(1,1)*n(i,1)*p(1,j)^2+C(1,3)*n(i,2)*p(2,j)*p(1,j)+C(1,2)*n(i,3)*p(3,j)*p(1,j)+C(5,5)*(n(i,1)*p(2,j)^2+n(i,2)*p(1,j)*p(2,j))+C(6,6)*(n(i,1)*p(3,j)^2+n(i,3)*p(1,j)*p(3,j)))/(Material.Density*X(j)*1e3); % group velocity component along x1
                g(2,j) = (C(5,5)*(n(i,1)*p(2,j)*p(1,j)+n(i,2)*p(1,j)^2)+C(1,3)*n(i,1)*p(1,j)*p(2,j)+C(3,3)*n(i,2)*p(2,j)^2+C(2,3)*n(i,3)*p(3,j)*p(2,j)+C(4,4)*(n(i,2)*p(3,j)^2+n(i,3)*p(2,j)*p(3,j)))/(Material.Density*X(j)*1e3); % group velocity component along x2
                g(3,j) = (C(6,6)*(n(i,1)*p(3,j)*p(1,j)+n(i,3)*p(1,j)^2)+C(4,4)*(n(i,2)*p(3,j)*p(2,j)+n(i,3)*p(2,j)^2)+C(1,2)*n(i,1)*p(1,j)*p(3,j)+C(2,3)*n(i,2)*p(2,j)*p(3,j)+C(2,2)*n(i,3)*p(3,j)^2)/(Material.Density*X(j)*1e3); % group velocity component along x3                
            end
            [G(i,3*(j-1)+1),G(i,3*(j-1)+2),G(i,3*(j-1)+3)] = cart2sph(g(1,j),g(2,j),g(3,j));
        end
    end
    G = vertcat(G,flipud(G(1:end-1,:)),G(2:end,:),flipud(G(1:end-1,:))); % extending the data to 360 deg
    G(:,[2,5,8]) = vertcat(G(1:length(n),[2,5,8]),pi-flipud(G(1:length(n)-1,[2,5,8])),pi+G(2:length(n),[2,5,8]),2*pi-flipud(G(1:length(n)-1,[2,5,8])));    
else
    G(:,2:3) = Material.TransverseVelocity/1e3;
    G(:,4) = Material.LongitudinalVelocity/1e3;
    G = vertcat(G,flipud(G(1:end-1,:)),G(2:end,:),flipud(G(1:end-1,:))); % extending the data to 360 deg
    G(:,1) = 0:2*pi/(4*(length(n)-1)):2*pi; % Theta (rad)
end
f = figure('Name','Group velocity 2-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 .6 1],'color','w');
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
if  ~strcmp(Material.Class,'Isotropic')
    polarplot(G(:,8),G(:,9),'LineWidth',LineWidth,'Color','r')
    z = find(abs(G(:,8)-pi/8) == min(abs(G(:,8)-pi/8)));
    text(pi/8,1.04*G(z(1),9),'L','FontSize',FontSizeModeLabels,'Interpreter','latex');
    hold on
    polarplot(G(:,2),G(:,3),'LineWidth',LineWidth,'Color','b')
    z = find(abs(G(:,2)-pi/4) == min(abs(G(:,2)-pi/4)));
    text(pi/4,1.04*G(z(1),3),'S$_\mathrm{fast}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
    polarplot(G(:,5),G(:,6),'LineWidth',LineWidth,'Color',[.13 .55 .13])
    z = find(abs(G(:,5)+pi) == min(abs(G(:,5)+pi)));
    text(pi,.95*G(z(1),6),'S$_\mathrm{slow}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
else
    polarplot(G(:,1),G(:,4),'LineWidth',LineWidth,'Color','r')
    z = find(abs(G(:,1)-pi/8) == min(abs(G(:,1)-pi/8)));
    text(pi/8,1.04*G(z(1),4),'L','FontSize',FontSizeModeLabels,'Interpreter','latex');
    hold on
    polarplot(G(:,1),G(:,2),'LineWidth',LineWidth,'Color','b')
    z = find(abs(G(:,1)-pi/4) == min(abs(G(:,1)-pi/4)));
    text(pi/4,1.04*G(z(1),2),'SV,SH','FontSize',FontSizeModeLabels,'Interpreter','latex');
end
ax = gca;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
if  HeadLine == 0
    ax.Title.FontSize = FontSizeAxesLabels;
    ax.Title.String = 'Group velocity (m/ms)';
elseif HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    if  Plane == 13
        ax.Title.String = ['Bulk wave group velocity @ $\phi_\mathrm r$ = ',num2str(Phi),'\,$^{\circ}$ in ',replace(Material.Name,'_','\_')];
    elseif Plane == 12
        ax.Title.String = ['Bulk wave group velocity @ $\theta_\mathrm r$ = ',num2str(Phi),'\,$^{\circ}$ in ',replace(Material.Name,'_','\_')];
    end
end
z = ax.RTick;
z(1) = '';
ax.RTick = z;
ax.ThetaTick = [0 45 90 135 180 225 270 315];
ax.ThetaAxis.TickLabelFormat = '%g';
ax.TickLabelInterpreter = 'latex';
if  HeadLine > 0
    text(.59,1.04,'(m/ms)','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
end
if  Plane == 13
    text(.41,1.03,'$v_{\mathrm g3}$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
elseif Plane == 12
    text(.41,1.03,'$v_{\mathrm g2}$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')    
end
text(1,.45,'$v_{\mathrm g1}$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
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
        output_txt = {['$\theta_\mathrm r$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$v_{\mathrm{g}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif Plane == 12
        output_txt = {['$\phi_\mathrm r$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$v_{\mathrm{g}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    end
end
end