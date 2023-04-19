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
function EnergyScattering_Isotropic(Fluid,Solid,Export,Theta,ThetaStep,Crop,Mode,PDF,PNG,FileName,Directory,HeadLine,LineWidth,BoxLineWidth,FontSizeHeadLine,FontSizeAxesLabels,FontSizeAxes,PNGresolution)
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
E(length(ThetaRange),4) = 0;
for i = 1:length(ThetaRange)
    PhaseVelocity = Fluid.Velocity/sind(ThetaRange(i)); % phase velocity component along x1 shared among all participating bulk waves, as required by Snell's law
    Alpha(1) = sqrt(PhaseVelocity^2/Solid.LongitudinalVelocity^2-1);
    Alpha(2) = sqrt(PhaseVelocity^2/Solid.TransverseVelocity^2-1);
    AlphaFluid = 1/Fluid.Velocity*PhaseVelocity*cosd(ThetaRange(i));
    W = (PhaseVelocity^2-Solid.LongitudinalVelocity^2-Solid.TransverseVelocity^2*Alpha.^2)./((Solid.LongitudinalVelocity^2-Solid.TransverseVelocity^2)*Alpha);
    D = 1i/PhaseVelocity*[0 0;0 0;Solid.Lambda+(Solid.Lambda+2*Solid.Mu)*Alpha.*W;0 0;Solid.Mu*(Alpha+W)]; % sigma33 sigma13
    DFluid = 1i*Fluid.Density*PhaseVelocity; % sigma11, sigma22, sigma33 in the fluid
    Z1 = [W AlphaFluid;D(3,:) -DFluid;D(5,:) 0];
    Z2 = [AlphaFluid;DFluid;0];
    U = Z1\Z2; % UT(L) UT(S) UR(L)
    v = -1i*[1 AlphaFluid;U(1) W(1)*U(1);U(2) W(2)*U(2);U(3) AlphaFluid*U(3)]; % I(L) T(L) T(S) R(L)
    sigma = [DFluid 0;D(3,1)*U(1) D(5,1)*U(1);D(3,2)*U(2) D(5,2)*U(2);DFluid*U(3) 0]; % I(L) T(L) T(S) R(L)
    P = -.5*[real(sigma(1,1)*conj(v(1,2)));real(sigma(2,1)*conj(v(2,2))+sigma(2,2)*conj(v(2,1)));real(sigma(3,1)*conj(v(3,2))+sigma(3,2)*conj(v(3,1)));real(sigma(4,1)*conj(v(4,2)))]; % I(L) T(L) T(S) R(L)
    E(i,1:3) = P(2:4)/P(1); % T(L) T(S) R(L) Sum
    if  E(i,3) > 1.001
        E(i,:) = [0 0 1 0];
    end
end
E(:,4) = sum(E,2);
% String = ['Fluid: ',Fluid.Name,newline,'Solid: ',Solid.Name,newline,newline,'Critical angles:'];
% if  isreal(asind(Fluid.Velocity/Solid.LongitudinalVelocity))
%     String = append(String,newline,'  Theta_L = ',num2str(asind(Fluid.Velocity/Solid.LongitudinalVelocity)),' deg');
% else
%     String = append(String,newline,'  Theta_L = -');
% end
% if  isreal(asind(Fluid.Velocity/Solid.TransverseVelocity))
%     String = append(String,newline,'  Theta_SV = ',num2str(asind(Fluid.Velocity/Solid.TransverseVelocity)),' deg');
% else
%     String = append(String,newline,'  Theta_SV = -');
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
%     String = append(String,newline,newline,...
%     'Plane wave incidence:',newline,...
%     '  Theta = ',num2str(Theta),' deg',newline,...
%     newline,...
%     'Energy scattering coefficients:',newline,...
%     '  R(L) = ',num2str(E(3)),newline,...
%     '  T(L) = ',num2str(E(1)),newline,...
%     '  T(SV) = ',num2str(E(2)),newline,...
%     '  Sum = ',num2str(E(4)),newline,'--------------------------------');
end
% disp(String)
if  Mode == 2
    f = figure('Name','Energy scattering','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    datacursormode on
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))))
    t = tiledlayout(2,2);
    if  HeadLine == 1
        t.Title.FontSize = FontSizeHeadLine;
        t.Title.Interpreter = 'latex';
        t.Title.String = ['Energy scattering on ',char(join(split(Fluid.Name,'_'),'\_')),'/',char(join(split(Solid.Name,'_'),'\_'))];
    end
    ax1 = nexttile;
    plot(0:ThetaStep:90,E(:,3),'LineWidth',LineWidth,'Color','b')
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
    plot(0:ThetaStep:90,E(:,1),'LineWidth',LineWidth,'Color','b')
    plot(0:ThetaStep:90,E(:,4),'LineWidth',LineWidth,'Color','r')
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
    text(.1,.9,'$R_{\mathrm{L}}+T_{\mathrm{L}}+T_{\mathrm{S}}$','Units','normalized','FontSize',.75*FontSizeAxes,'Interpreter','latex')
    ax3 = nexttile;
    plot(0:ThetaStep:90,E(:,2),'LineWidth',LineWidth,'Color','b')
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
    ax3.Title.String = '$T_{\mathrm{S}}$';
    ax3.XLabel.String = 'Incidence angle ($^{\circ}$)';
    ax3.YLabel.String = 'Transmission coefficient';
    ax3.XLim = [0 90];
    ax3.YLim = [0 1.05];
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
    if  event_obj.Target.YData(1) == E(1,3)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['$R_{\mathrm{L}}$: \textbf{',num2str(event_obj.Position(2),5),'}']};
    elseif event_obj.Target.YData(1) == E(1,1)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['$T_{\mathrm{L}}$: \textbf{',num2str(event_obj.Position(2),5),'}']};
    elseif event_obj.Target.YData(1) == E(1,2)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['$T_{\mathrm{S}}$: \textbf{',num2str(event_obj.Position(2),5),'}']};
    elseif event_obj.Target.YData(1) == E(1,4)
        OutputTxt = {['$\theta_{\mathrm{I}}$: \textbf{',num2str(event_obj.Position(1)),'}\,$^\circ$'] ['Sum: \textbf{',num2str(event_obj.Position(2),5),'}']};
    end
end
end