% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2025 DLR
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
function PhaseVelocity_Polar(Hybrid,LayupString,PropagationAngleResolution,BulkVelocities,A,A0,Directory,Export,FileName,FontSizeAxes,FontSizeHeadLine,FontSizeModeLabels,Frequency,FrequencyRange,HeadLine,LineWidth,Material,PDF,PlateThickness,PNG,PNGresolution,PropagationAngle,PropagationAngleLimit,Repetitions,S0,SH0,SuperLayerSize,SymmetricSystem)
%#ok<*AGROW>
q = find(FrequencyRange == Frequency);
if  isempty(q)
    if  Frequency > ceil(max(FrequencyRange)) || Frequency < 0
        errordlg(['Selected frequency outside frequency range! Select between 0 and ',num2str(ceil(max(FrequencyRange))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(FrequencyRange-Frequency));
        Frequency = FrequencyRange(q);
    end
end
for n = 1:length(PropagationAngle)
    if  A0
        X(n,1) = A{n,1}(q,1)/1e3;
    end
    if  SH0
        X(n,2) = A{n,2}(q,1)/1e3;
    end
    if  S0
        X(n,3) = A{n,3}(q,1)/1e3;
    end
end
if  SuperLayerSize == 1 && BulkVelocities == 1  
    C = Material{1}.C;
    a21 = (C(1,1)*C(2,2)+C(1,1)*C(4,4)+C(2,2)*C(5,5)+C(4,4)*C(6,6)+C(5,5)*C(6,6)-2*C(1,2)*C(6,6)-C(1,2)^2)/Material{1}.Density^2;
    a24 = (C(1,1)*C(5,5)+C(1,1)*C(6,6)+C(5,5)*C(6,6))/Material{1}.Density^2;
    a25 = (C(2,2)*C(4,4)+C(2,2)*C(6,6)+C(4,4)*C(6,6))/Material{1}.Density^2;
    a32 = (C(1,2)^2*C(4,4)-C(1,1)*C(2,2)*C(4,4)-C(2,2)*C(5,5)*C(6,6)+2*C(1,2)*C(4,4)*C(6,6))/Material{1}.Density^3;
    a34 = (C(1,2)^2*C(5,5)-C(1,1)*C(2,2)*C(5,5)-C(1,1)*C(4,4)*C(6,6)+2*C(1,2)*C(5,5)*C(6,6))/Material{1}.Density^3;
    a38 = -C(1,1)*C(5,5)*C(6,6)/Material{1}.Density^3;
    a39 = -C(2,2)*C(4,4)*C(6,6)/Material{1}.Density^3;
    if  PropagationAngleLimit == 1
        Phi = 0:PropagationAngleResolution:180;
    elseif PropagationAngleLimit == 2
        Phi = 0:PropagationAngleResolution:90;
    end
    n = [cosd(Phi)' sind(Phi)'];
    A1 = -(C(1,1)*n(:,1).^2+C(2,2)*n(:,2).^2+C(4,4)*n(:,2).^2+C(5,5)*n(:,1).^2+C(6,6)*n(:,1).^2+C(6,6)*n(:,2).^2)/Material{1}.Density;
    A2 = a21*n(:,1).^2.*n(:,2).^2+a24*n(:,1).^4+a25*n(:,2).^4;
    A3 = a32*n(:,1).^2.*n(:,2).^4+a34*n(:,1).^4.*n(:,2).^2+a38*n(:,1).^6+a39*n(:,2).^6;
    Xa = A2/3-A1.^2/9;
    Xb = A1.^3/27-A1.*A2/6+A3/2;
    Xc = (sqrt(Xb.^2+Xa.^3)-Xb).^(1/3);
    Xd = Xa./(2*Xc)-Xc/2;
    Xe = Xa./Xc;
    Xf = (sqrt(3)*(Xc+Xe)*1i)/2;
    X(:,5) = abs(real(sqrt(Xd-Xf-A1/3)))/1e3; % phase velocity in the solid (m/ms)
    X(:,6) = abs(real(sqrt(Xd+Xf-A1/3)))/1e3;
    X(:,7) = abs(real(-sqrt(Xc-Xe-A1/3)))/1e3;
end
if  PropagationAngleLimit == 1
    X = vertcat(X,(X(2:end,:)));
    X(:,4) = 0:2*pi/(2*length(PropagationAngle)-2):2*pi;
elseif PropagationAngleLimit == 2
    X = vertcat(X,flipud(X(1:end-1,:)),X(2:end,:),flipud(X(1:end-1,:))); % extending data for polar plot; the wave propagation angle in column 4 is given from 0 to 2*pi
    X(:,4) = 0:2*pi/(4*(length(PropagationAngle)-1)):2*pi;
end
f = figure('Name','Phase velocity profile','Toolbar','none','Units','normalized','OuterPosition',[0 0 .6 1],'color','w');
datacursormode on    
if  SuperLayerSize == 1 || SymmetricSystem
    if  SuperLayerSize == 1 && BulkVelocities == 1
        polarplot(X(:,4),X(:,5),'LineWidth',LineWidth,'Color','c')
        z = find(abs(X(:,4)-pi/16) == min(abs(X(:,4)-pi/16)));
        text(pi/16,1.05*X(z(1),5),'S$_{\mathrm{fast}}$','FontSize',FontSizeModeLabels,'Color','c','Interpreter','latex');
        hold on
        polarplot(X(:,4),X(:,6),'LineWidth',LineWidth,'Color','c')
        z = find(abs(X(:,4)-1.3*pi) == min(abs(X(:,4)-1.3*pi)));
        text(1.3*pi,X(z(1),6),'S$_{\mathrm{slow}}$','FontSize',FontSizeModeLabels,'Color','c','Interpreter','latex');
        polarplot(X(:,4),X(:,7),'LineWidth',LineWidth,'Color','c')
        z = find(abs(X(:,4)-pi/16) == min(abs(X(:,4)-pi/16)));
        text(pi/16,1.04*X(z(1),7),'L','FontSize',FontSizeModeLabels,'Color','c','Interpreter','latex');
    end        
    if  S0
        polarplot(X(:,4),X(:,3),'LineWidth',LineWidth,'Color','r')
        z = find(abs(X(:,4)-pi/8) == min(abs(X(:,4)-pi/8)));
        text(pi/8,1.04*X(z(1),3),'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        hold on
    end
    if  SH0
        polarplot(X(:,4),X(:,2),'LineWidth',LineWidth,'Color',[1 .7 0])
        z = find(abs(X(:,4)-pi/4) == min(abs(X(:,4)-pi/4)));
        text(pi/4,1.05*X(z(1),2),'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        hold on
    end
    if  A0
        polarplot(X(:,4),X(:,1),'LineWidth',LineWidth,'Color','b')
        z = find(abs(X(:,4)-pi) == min(abs(X(:,4)-pi)));
        text(pi,.95*X(z(1),1),'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
    end
else
    if  S0
        polarplot(X(:,4),X(:,3),'LineWidth',LineWidth,'Color',[.5 0 1])
        z = find(abs(X(:,4)-pi/8) == min(abs(X(:,4)-pi/8)));
        text(pi/8,1.04*X(z(1),3),'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        hold on
    end
    if  SH0
        polarplot(X(:,4),X(:,2),'LineWidth',LineWidth,'Color',[1 0 1])
        z = find(abs(X(:,4)-pi/4) == min(abs(X(:,4)-pi/4)));
        text(pi/4,1.05*X(z(1),2),'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        hold on
    end
    if  A0
        polarplot(X(:,4),X(:,1),'LineWidth',LineWidth,'Color',[1 0 .5])
        z = find(abs(X(:,4)-pi) == min(abs(X(:,4)-pi)));
        text(pi,.95*X(z(1),1),'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex'); 
    end
end
ax = gca;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
if  Hybrid
    Material{1}.Name = 'hybrid';
end
if  HeadLine == 0
    ax.Title.FontSize = FontSizeAxes*1.25;
    ax.Title.String = 'Phase velocity (m/ms)';
elseif HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' @ ',num2str(Frequency),'\,kHz'];
elseif ~SymmetricSystem  && HeadLine == 2
    ax.Title.FontSize = FontSizeHeadLine;
    if  Repetitions == 1
        ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,'] @ ',num2str(Frequency),'\,kHz'];
    elseif Repetitions > 1 
        ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$ @ ',num2str(Frequency),'\,kHz'];
    end
elseif SymmetricSystem && HeadLine == 2
    ax.Title.FontSize = FontSizeHeadLine;
    if  Repetitions == 1
        ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$ @ ',num2str(Frequency),'\,kHz'];
    elseif Repetitions > 1
        ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$ @ ',num2str(Frequency),'\,kHz'];
    end
end
z = ax.RTick;
z(1) = '';
ax.RTick = z;
ax.ThetaTick = [0 45 90 135 180 225 270 315];
ax.ThetaAxis.TickLabelFormat = '%g';
ax.TickLabelInterpreter = 'latex';
if  HeadLine > 0
    text(.59,1.04,'(m/ms)','units','normalized','FontSize',FontSizeAxes*1.25,'interpreter','latex')
end
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
    output_txt = {['$\phi$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
end    
end