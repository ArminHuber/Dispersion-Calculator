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
function EnergyScattering_Anisotropic(Fluid,Solid,Export,Phi,Theta,ThetaStep,Crop,Mode,PDF,PNG,FileName,Directory,HeadLine,LineWidth,BoxLineWidth,FontSizeHeadLine,FontSizeAxesLabels,FontSizeAxes,PNGresolution)
Phi2 = Phi;
if  Phi2 == 0
    Phi2 = 0.001;
elseif Phi2 == 45 && strcmp(Solid.Class,'Cubic')
    Phi2 = 45.001;
elseif Phi2 == 90
    Phi2 = 89.999;
end
C = real(Solid.C);
n = [cosd(Phi2) sind(Phi2)];
A1 = -(C(1,1)*n(1)^2+C(2,2)*n(2)^2+C(4,4)*n(2)^2+C(5,5)*n(1)^2+C(6,6)*n(1)^2+C(6,6)*n(2)^2)/Solid.Density;
A2 = ((C(1,1)*C(2,2)+C(1,1)*C(4,4)+C(2,2)*C(5,5)-2*C(1,2)*C(6,6)+C(4,4)*C(6,6)+C(5,5)*C(6,6)-C(1,2)^2)*n(1)^2*n(2)^2+(C(1,1)*C(5,5)+C(1,1)*C(6,6)+C(5,5)*C(6,6))*n(1)^4+(C(2,2)*C(4,4)+C(2,2)*C(6,6)+C(4,4)*C(6,6))*n(2)^4)/Solid.Density^2;
A3 = ((C(1,2)^2*C(4,4)-C(1,1)*C(2,2)*C(4,4)+2*C(1,2)*C(4,4)*C(6,6)-C(2,2)*C(5,5)*C(6,6))*n(1)^2*n(2)^4+(C(1,2)^2*C(5,5)-C(1,1)*C(2,2)*C(5,5)-C(1,1)*C(4,4)*C(6,6)+2*C(1,2)*C(5,5)*C(6,6))*n(1)^4*n(2)^2-C(1,1)*C(5,5)*C(6,6)*n(1)^6-C(2,2)*C(4,4)*C(6,6)*n(2)^6)/Solid.Density^3;
Xa = A2/3-A1^2/9;
Xb = A1^3/27-A1*A2/6+A3/2;
Xc = (sqrt(Xb^2+Xa^3)-Xb)^(1/3);
Xd = Xa/(2*Xc)-Xc/2;
Xe = Xa/Xc;
Xf = (sqrt(3)*(Xc+Xe)*1i)/2;
X(1) = abs(real(sqrt(Xd-Xf-A1/3)));
X(2) = abs(real(sqrt(Xd+Xf-A1/3)));
X(3) = abs(real(-sqrt(Xc-Xe-A1/3)));
ThetaCritical_S_fast = asind(Fluid.Velocity/X(1));
ThetaCritical_S_slow = asind(Fluid.Velocity/X(2));
% ThetaCritical_L = asind(Fluid.Velocity/X(3));
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
a11 = (c(1,1)*c(3,3)*c(4,4)+c(3,3)*c(5,5)*c(6,6)-c(3,6)^2*c(5,5)-c(1,3)^2*c(4,4)+2*(c(1,3)*c(3,6)*c(4,5)+c(1,3)*c(4,5)^2-c(1,3)*c(4,4)*c(5,5)-c(1,6)*c(3,3)*c(4,5)))/Delta;
a12 = (c(4,5)^2-c(3,3)*c(4,4)-c(3,3)*c(5,5)-c(4,4)*c(5,5))/Delta;
a21 = (c(1,1)*c(3,3)*c(6,6)+c(1,1)*c(4,4)*c(5,5)-c(1,1)*c(3,6)^2-c(1,1)*c(4,5)^2-c(1,3)^2*c(6,6)-c(1,6)^2*c(3,3)+2*(c(1,6)*c(3,6)*c(5,5)+c(1,3)*c(1,6)*c(3,6)+c(1,3)*c(1,6)*c(4,5)-c(1,1)*c(3,6)*c(4,5)-c(1,3)*c(5,5)*c(6,6)))/Delta;
a22 = (c(1,3)^2+c(4,5)^2+c(3,6)^2-c(1,1)*c(3,3)-c(1,1)*c(4,4)-c(3,3)*c(6,6)-c(5,5)*c(6,6)-c(4,4)*c(5,5)+2*(c(1,3)*c(5,5)+c(1,6)*c(4,5)+c(3,6)*c(4,5)))/Delta;
a23 = (c(4,4)+c(3,3)+c(5,5))/Delta;
a31 = (c(1,1)*c(5,5)*c(6,6)-c(1,6)^2*c(5,5))/Delta;
a32 = (c(1,6)^2-c(5,5)*c(6,6)-c(1,1)*c(5,5)-c(1,1)*c(6,6))/Delta;
a33 = (c(1,1)+c(5,5)+c(6,6))/Delta;
a34 = -1/Delta;
if  Mode == 1
    ThetaRange = Theta;
    if  ThetaRange == 0
        ThetaRange = 0.001;
    elseif ThetaRange == 90
        ThetaRange = 89.999;
    end
elseif Mode == 2
    ThetaRange = 0:ThetaStep:90;
    ThetaRange(1) = .001;
    ThetaRange(end) = 89.999;   
end
E(length(ThetaRange),5) = 0;
for i = 1:length(ThetaRange)
    PhaseVelocity = Fluid.Velocity/sind(ThetaRange(i)); % phase velocity component along x1 shared among all participating bulk waves, as required by Snell's law
    rc2 = Solid.Density*PhaseVelocity^2;
    r2c4 = Solid.Density^2*PhaseVelocity^4;
    A1 = a11+a12*rc2;
    A2 = a21+a22*rc2+a23*r2c4;
    A3 = a31+a32*rc2+a33*r2c4+a34*Solid.Density^3*PhaseVelocity^6;
    Alphaa = A2/3-A1^2/9;
    Alphab = A1^3/27-A1*A2/6+A3/2;
    Alphac = (sqrt(Alphab^2+Alphaa^3)-Alphab)^(1/3);
    Alphad = Alphaa/(2*Alphac)-Alphac/2;
    Alphae = Alphaa/Alphac;
    Alphaf = (sqrt(3)*(Alphac+Alphae)*1i)/2;
    Alpha(1) = sqrt(Alphad-Alphaf-A1/3);
    Alpha(3) = sqrt(Alphad+Alphaf-A1/3);
    Alpha(2) = sqrt(Alphac-Alphae-A1/3);
    AlphaFluid = 1/Fluid.Velocity*PhaseVelocity*cosd(ThetaRange(i));
    if  angle(Alpha(1)) < 0
        Alpha(1) = conj(Alpha(1)); % T(S_fast)
    end
    if  angle(Alpha(2)) < 0
        Alpha(2) = conj(Alpha(2)); % T(S_slow)
    end
    if  angle(Alpha(3)) < 0
        Alpha(3) = conj(Alpha(3)); % T(L)
    end
    if  ThetaRange(i) > ThetaCritical_S_slow
        Alpha(2) = -Alpha(2); % T(S_slow)
    end
    if  strcmp(Solid.Class,'Cubic') && ThetaRange(i) > ThetaCritical_S_fast
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
    DFluid = 1i*Fluid.Density*PhaseVelocity; % sigma11, sigma22, sigma33 in the fluid
    Z1 = [W AlphaFluid;D(3,:) -DFluid;D(5,:) 0;D(4,:) 0];
    Z2 = [AlphaFluid;DFluid;0;0];
    U = Z1\Z2; % UT(S_fast) UT(S_slow) UT(L) UR(L)
    v = -1i*[1 0 AlphaFluid;U(1) V(1)*U(1) W(1)*U(1);U(2) V(2)*U(2) W(2)*U(2);U(3) V(3)*U(3) W(3)*U(3);U(4) 0 AlphaFluid*U(4)]; % I(L) T(S_fast) T(S_slow) T(L) R(L)
    sigma = [DFluid 0 0;D(3,1)*U(1) D(5,1)*U(1) D(4,1)*U(1);D(3,2)*U(2) D(5,2)*U(2) D(4,2)*U(2);D(3,3)*U(3) D(5,3)*U(3) D(4,3)*U(3);DFluid*U(4) 0 0]; % I(L) T(S_fast) T(S_slow) T(L) R(L)
    P = -.5*[real(sigma(1,1)*conj(v(1,3)));real(sigma(2,1)*conj(v(2,3))+sigma(2,2)*conj(v(2,1))+sigma(2,3)*conj(v(2,2)));real(sigma(3,1)*conj(v(3,3))+sigma(3,2)*conj(v(3,1))+sigma(3,3)*conj(v(3,2)));real(sigma(4,1)*conj(v(4,3))+sigma(4,2)*conj(v(4,1))+sigma(4,3)*conj(v(4,2)));real(sigma(5,1)*conj(v(5,3)))]; % I(L) T(S_fast) T(S_slow) T(L) R(L)
    E(i,1:4) = P(2:5)/P(1); % T(S_fast) T(S_slow) T(L) R(L) Sum
    if  E(i,4) > 1.001
        E(i,:) = [0 0 0 1 0];
    end
% AlphaX(i,:) = Alpha;
end
E(:,5) = sum(E,2);
% String = ['Fluid: ',Fluid.Name,newline,'Solid: ',Solid.Name,newline,newline,'Critical angles @ Phi = ',num2str(Phi),' deg:'];
% if  isreal(ThetaCritical_L)
%     String = append(String,newline,'  Theta_L = ',num2str(ThetaCritical_L),' deg');
% else
%     disp('  Theta_L = -')
% end
% if  isreal(ThetaCritical_S_fast)
%     String = append(String,newline,'  Theta_Sfast = ',num2str(ThetaCritical_S_fast),' deg');
% else
%     String = append(String,newline,'  Theta_Sfast = -');
% end
% if  isreal(ThetaCritical_S_slow)
%     String = append(String,newline,'  Theta_Sslow = ',num2str(ThetaCritical_S_slow),' deg');
% else
%     String = append(String,newline,'  Theta_Sslow = -');
% end
if  Mode == 2
%     String = append(String,newline,'--------------------------------');
elseif Mode == 1
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
%     String = append(String,newline,newline,...
%     'Plane wave incidence:',newline,...
%     '  Phi = ',num2str(Phi),' deg',newline,...
%     '  Theta = ',num2str(Theta),' deg',newline,...
%     newline,...
%     'Energy scattering coefficients:',newline,...
%     '  R(L) = ',num2str(E(4)),newline,...
%     '  T(L) = ',num2str(E(3)),newline,...
%     '  T(Sfast) = ',num2str(E(1)),newline,...
%     '  T(Sslow) = ',num2str(E(2)),newline,...
%     '  Sum = ',num2str(E(5)),newline,'--------------------------------');
end
% disp(String)
% figure('name','angle');plot([angle(AlphaX(100:end,1)),angle(AlphaX(100:end,2)),angle(AlphaX(100:end,3))])
% figure('name','real');plot([real(AlphaX(100:end,1)),real(AlphaX(100:end,2)),real(AlphaX(100:end,3))])
% figure('name','imag');plot([imag(AlphaX(100:end,1)),imag(AlphaX(100:end,2)),imag(AlphaX(100:end,3))])
if  Mode == 2
    f = figure('Name','Energy scattering','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    datacursormode on
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))))
    t = tiledlayout(2,2);
    if  HeadLine == 1
        t.Title.FontSize = FontSizeHeadLine;
        t.Title.Interpreter = 'latex';
        t.Title.String = ['Energy scattering for $\phi$ = ',num2str(Phi,'%.0f'),'\,$^{\circ}$ on ',char(join(split(Fluid.Name,'_'),'\_')),'/',char(join(split(Solid.Name,'_'),'\_'))];
    end
    ax1 = nexttile;
    plot(0:ThetaStep:90,E(:,4),'LineWidth',LineWidth,'Color','b')
    ax1.YLim = [0 1.05];
    ax1.Box = 'on';
    ax1.LineWidth = BoxLineWidth;
    ax1.FontSize = .75*FontSizeAxes;
    ax1.Title.FontSize = .75*FontSizeHeadLine;
    ax1.XLabel.FontSize = .75*FontSizeAxesLabels;
    ax1.YLabel.FontSize = .75*FontSizeAxesLabels;
    ax1.Title.Interpreter = 'latex';
    ax1.XLabel.Interpreter = 'latex';
    ax1.YLabel.Interpreter = 'latex';
    ax1.TickLabelInterpreter = 'latex';
    ax1.Title.String = '$R_{\mathrm{L}}$';
    ax1.XLabel.String = 'Incidence angle ($^{\circ}$)';
    ax1.YLabel.String = 'Reflection coefficient';
    ax1.XLim = [0 90];
    ax1.YLim = [0 1.05];
    ax2 = nexttile;
    hold on
    plot(0:ThetaStep:90,E(:,3),'LineWidth',LineWidth,'Color','b')
    plot(0:ThetaStep:90,E(:,5),'LineWidth',LineWidth,'Color','r')
    ax2.YLim = [0 1.05];
    ax2.Box = 'on';
    ax2.LineWidth = BoxLineWidth;
    ax2.FontSize = .75*FontSizeAxes;
    ax2.Title.FontSize = .75*FontSizeHeadLine;
    ax2.XLabel.FontSize = .75*FontSizeAxesLabels;
    ax2.YLabel.FontSize = .75*FontSizeAxesLabels;
    ax2.Title.Interpreter = 'latex';
    ax2.XLabel.Interpreter = 'latex';
    ax2.YLabel.Interpreter = 'latex';
    ax2.TickLabelInterpreter = 'latex';
    ax2.Title.String = '$T_{\mathrm{L}}$';
    ax2.XLabel.String = 'Incidence angle ($^{\circ}$)';
    ax2.YLabel.String = 'Transmission coefficient';
    ax2.XLim = [0 90];
    ax2.YLim = [0 1.05];
    text(.1,.9,'$R_{\mathrm{L}}+T_{\mathrm{L}}+T_{\mathrm{S_{fast}}}+T_{\mathrm{S_{slow}}}$','Units','normalized','FontSize',.75*FontSizeAxes,'Interpreter','latex')
    ax3 = nexttile;
    plot(0:ThetaStep:90,E(:,1),'LineWidth',LineWidth,'Color','b')
    ax3.YLim = [0 1.05];
    ax3.Box = 'on';
    ax3.LineWidth = BoxLineWidth;
    ax3.FontSize = .75*FontSizeAxes;
    ax3.Title.FontSize = .75*FontSizeHeadLine;
    ax3.XLabel.FontSize = .75*FontSizeAxesLabels;
    ax3.YLabel.FontSize = .75*FontSizeAxesLabels;
    ax3.Title.Interpreter = 'latex';
    ax3.XLabel.Interpreter = 'latex';
    ax3.YLabel.Interpreter = 'latex';
    ax3.TickLabelInterpreter = 'latex';
    ax3.Title.String = '$T_{\mathrm{S_{fast}}}$';
    ax3.XLabel.String = 'Incidence angle ($^{\circ}$)';
    ax3.YLabel.String = 'Transmission coefficient';
    ax3.XLim = [0 90];
    ax3.YLim = [0 1.05];
    ax4 = nexttile;
    plot(0:ThetaStep:90,E(:,2),'LineWidth',LineWidth,'Color','b')
    ax4.YLim = [0 1.05];
    ax4.Box = 'on';
    ax4.LineWidth = BoxLineWidth;
    ax4.FontSize = .75*FontSizeAxes;
    ax4.Title.FontSize = .75*FontSizeHeadLine;
    ax4.XLabel.FontSize = .75*FontSizeAxesLabels;
    ax4.YLabel.FontSize = .75*FontSizeAxesLabels;
    ax4.Title.Interpreter = 'latex';
    ax4.XLabel.Interpreter = 'latex';
    ax4.YLabel.Interpreter = 'latex';
    ax4.TickLabelInterpreter = 'latex';
    ax4.Title.String = '$T_{\mathrm{S_{slow}}}$';
    ax4.XLabel.String = 'Incidence angle ($^{\circ}$)';
    ax4.YLabel.String = 'Transmission coefficient';
    ax4.XLim = [0 90];
    ax4.YLim = [0 1.05];
    if  Export == 1
        try
            if  PDF == 1
                if  Crop == 0
                    set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[50 30])
                    print(f,fullfile(Directory,FileName),'-dpdf')
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
end
function OutputTxt = Cursor(~,event_obj)
    if  event_obj.Target.YData(1) == E(1,4)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['$R_{\mathrm{L}}$: \textbf{',num2str(event_obj.Position(2),5),'}']};
    elseif event_obj.Target.YData(1) == E(1,3)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['$T_{\mathrm{L}}$: \textbf{',num2str(event_obj.Position(2),5),'}']};
    elseif event_obj.Target.YData(1) == E(1,1)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['$T_{\mathrm{S_{fast}}}$: \textbf{',num2str(event_obj.Position(2),5),'}']};
    elseif event_obj.Target.YData(1) == E(1,2)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['$T_{\mathrm{S_{slow}}}$: \textbf{',num2str(event_obj.Position(2),5),'}']};
    elseif event_obj.Target.YData(1) == E(1,5)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['Sum: \textbf{',num2str(event_obj.Position(2),5),'}']};
    end
end 
end