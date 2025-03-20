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
function PropagationTime_Polar(Mode,Hybrid,LayupString,Distance,A,A0,Directory,Export,FileName,FontSizeAxes,FontSizeHeadLine,FontSizeModeLabels,Frequency,FrequencyRange,HeadLine,LineWidth,Material,PDF,PlateThickness,PNG,PNGresolution,PropagationAngle,PropagationAngleLimit,Repetitions,S0,SH0,SuperLayerSize,SymmetricSystem)
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
    if  Mode == 1
        if  A0
            X(n,1) = Distance/A{n,1}(q,2)*1e3;
        end
        if  SH0
            X(n,2) = Distance/A{n,2}(q,2)*1e3;
        end
        if  S0
            X(n,3) = Distance/A{n,3}(q,2)*1e3;
        end
    elseif Mode == 2
        if  A0
            X(n,1) = Distance/sqrt(A{n,1}(q,2)^2+A{n,1}(q,3)^2)*1e3;
            SkewAngle(n,1) = -atand(A{n,1}(q,3)/A{n,1}(q,2));
        end
        if  SH0
            X(n,2) = Distance/sqrt(A{n,2}(q,2)^2+A{n,2}(q,3)^2)*1e3;
            SkewAngle(n,2) = -atand(A{n,2}(q,3)/A{n,2}(q,2));
        end
        if  S0
            X(n,3) = Distance/sqrt(A{n,3}(q,2)^2+A{n,3}(q,3)^2)*1e3;
            SkewAngle(n,3) = -atand(A{n,3}(q,3)/A{n,3}(q,2));
        end        
    end
end
if  Mode == 2 && abs(SkewAngle(1,1)) < .5
    SkewAngle(1,:) = 0;
    SkewAngle(end,:) = 0;
    if  PropagationAngleLimit == 1
        SkewAngle(ceil(.5*length(SkewAngle)),:) = 0;
    end
end
if  PropagationAngleLimit == 1
    X = vertcat(X,(X(2:end,:)));
    if  Mode == 1
        X(:,4) = 0:2*pi/(2*length(PropagationAngle)-2):2*pi;
    elseif Mode == 2
        if  A0
            RayAngle(:,1) = deg2rad(PropagationAngle'-SkewAngle(:,1));
            X(:,8) = vertcat(RayAngle(:,1),pi+RayAngle(2:end,1));
        end
        if  SH0
            RayAngle(:,2) = deg2rad(PropagationAngle'-SkewAngle(:,2));
            X(:,9) = vertcat(RayAngle(:,2),pi+RayAngle(2:end,2));
        end
        if  S0
            RayAngle(:,3) = deg2rad(PropagationAngle'-SkewAngle(:,3));
            X(:,10) = vertcat(RayAngle(:,3),pi+RayAngle(2:end,3));
        end 
    end
elseif PropagationAngleLimit == 2
    X = vertcat(X,flipud(X(1:end-1,:)),X(2:end,:),flipud(X(1:end-1,:))); % extending data for polar plot; the wave propagation angle in column 4 is given from 0 to 2*pi
    if  Mode == 1
        X(:,4) = 0:2*pi/(4*(length(PropagationAngle)-1)):2*pi;
    elseif Mode == 2
        if  A0
            RayAngle(:,1) = deg2rad(PropagationAngle'-SkewAngle(:,1));
            X(:,8) = vertcat(RayAngle(:,1),pi-flipud(RayAngle(1:end-1,1)),pi+RayAngle(2:end,1),2*pi-flipud(RayAngle(1:end-1,1)));
        end
        if  SH0
            RayAngle(:,2) = deg2rad(PropagationAngle'-SkewAngle(:,2));
            X(:,9) = vertcat(RayAngle(:,2),pi-flipud(RayAngle(1:end-1,2)),pi+RayAngle(2:end,2),2*pi-flipud(RayAngle(1:end-1,2)));
        end
        if  S0
            RayAngle(:,3) = deg2rad(PropagationAngle'-SkewAngle(:,3));
            X(:,10) = vertcat(RayAngle(:,3),pi-flipud(RayAngle(1:end-1,3)),pi+RayAngle(2:end,3),2*pi-flipud(RayAngle(1:end-1,3)));
        end
    end
end
if  Mode == 1
    f = figure('Name','Propagation time profile for energy velocity ce1','Toolbar','none','Units','normalized','OuterPosition',[0 0 .6 1],'color','w');
elseif Mode == 2
    f = figure('Name','Propagation time profile for energy velocity absolute','Toolbar','none','Units','normalized','OuterPosition',[0 0 .6 1],'color','w');    
end
datacursormode on
if  SuperLayerSize == 1 || SymmetricSystem
    if  A0
        if  Mode == 1
            polarplot(X(:,4),X(:,1),'LineWidth',LineWidth,'Color','b')
            z = find(abs(X(:,4)-pi/8) == min(abs(X(:,4)-pi/8)));
            text(pi/8,1.04*X(z(1),1),'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
            hold on
        elseif Mode == 2
            polarplot(X(:,8),X(:,1),'LineWidth',LineWidth,'Color','b')
            z = find(abs(X(:,8)-pi/8) == min(abs(X(:,8)-pi/8)));
            text(pi/8,1.04*X(z(1),1),'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
            hold on
        end
    end
    if  SH0
        if  Mode == 1
            polarplot(X(:,4),X(:,2),'LineWidth',LineWidth,'Color',[1 .7 0])
            z = find(abs(X(:,4)-pi/4) == min(abs(X(:,4)-pi/4)));
            text(pi/4,1.05*X(z(1),2),'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
            hold on
        elseif Mode == 2
            polarplot(X(:,9),X(:,2),'LineWidth',LineWidth,'Color',[1 .7 0])
            z = find(abs(X(:,9)-pi/4) == min(abs(X(:,9)-pi/4)));
            text(pi/4,1.05*X(z(1),2),'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
            hold on            
        end
    end
    if  S0
        if  Mode == 1
            polarplot(X(:,4),X(:,3),'LineWidth',LineWidth,'Color','r')
            z = find(abs(X(:,4)-pi) == min(abs(X(:,4)-pi)));
            text(pi,.95*X(z(1),3),'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        elseif Mode == 2
            polarplot(X(:,10),X(:,3),'LineWidth',LineWidth,'Color','r')
            z = find(abs(X(:,10)-pi) == min(abs(X(:,10)-pi)));
            text(pi,.95*X(z(1),3),'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');            
        end
    end
else
    if  A0
        if  Mode == 1
            polarplot(X(:,4),X(:,1),'LineWidth',LineWidth,'Color',[1 0 .5])
            z = find(abs(X(:,4)-pi/8) == min(abs(X(:,4)-pi/8)));
            text(pi/8,1.04*X(z(1),1),'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
            hold on
        elseif Mode == 2
            polarplot(X(:,8),X(:,1),'LineWidth',LineWidth,'Color',[1 0 .5])
            z = find(abs(X(:,8)-pi/8) == min(abs(X(:,8)-pi/8)));
            text(pi/8,1.04*X(z(1),1),'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
            hold on            
        end
    end
    if  SH0
        if  Mode == 1
            polarplot(X(:,4),X(:,2),'LineWidth',LineWidth,'Color',[1 0 1])
            z = find(abs(X(:,4)-pi/4) == min(abs(X(:,4)-pi/4)));
            text(pi/4,1.05*X(z(1),2),'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
            hold on
        elseif Mode == 2
            polarplot(X(:,9),X(:,2),'LineWidth',LineWidth,'Color',[1 0 1])
            z = find(abs(X(:,9)-pi/4) == min(abs(X(:,9)-pi/4)));
            text(pi/4,1.05*X(z(1),2),'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
            hold on            
        end
    end
    if  S0
        if  Mode == 1
            polarplot(X(:,4),X(:,3),'LineWidth',LineWidth,'Color',[.5 0 1])
            z = find(abs(X(:,4)-pi) == min(abs(X(:,4)-pi)));
            text(pi,.95*X(z(1),3),'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        elseif Mode == 2
            polarplot(X(:,10),X(:,3),'LineWidth',LineWidth,'Color',[.5 0 1])
            z = find(abs(X(:,10)-pi) == min(abs(X(:,10)-pi)));
            text(pi,.95*X(z(1),3),'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');            
        end
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
    if  Mode == 1
        ax.Title.String = '$t_{\mathrm{Prop}}$ for $c_{\mathrm{e}1}$ ($\mu$s)';
    elseif Mode == 2
        ax.Title.String = '$t_{\mathrm{Prop}}$ for $|\vec{c}_{\mathrm e}|$ ($\mu$s)';
    end
elseif HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    if  Mode == 1
        ax.Title.String = ['$t_{\mathrm{Prop}}$ for $c_{\mathrm{e}1}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
    elseif Mode == 2
        ax.Title.String = ['$t_{\mathrm{Prop}}$ for $|\vec{c}_{\mathrm e}|$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
    end
elseif ~SymmetricSystem  && HeadLine == 2
    ax.Title.FontSize = FontSizeHeadLine;
    if  Repetitions == 1
        if  Mode == 1
            ax.Title.String = ['$t_{\mathrm{Prop}}$ for $c_{\mathrm{e}1}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,'] @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
        elseif Mode == 2
            ax.Title.String = ['$t_{\mathrm{Prop}}$ for $|\vec{c}_{\mathrm e}|$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,'] @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
        end
    elseif Repetitions > 1
        if  Mode == 1
            ax.Title.String = ['$t_{\mathrm{Prop}}$ for $c_{\mathrm{e}1}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$ @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
        elseif Mode == 2
            ax.Title.String = ['$t_{\mathrm{Prop}}$ for $|\vec{c}_{\mathrm e}|$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$ @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
        end
    end
elseif SymmetricSystem && HeadLine == 2
    ax.Title.FontSize = FontSizeHeadLine;
    if  Repetitions == 1
        if  Mode == 1
            ax.Title.String = ['$t_{\mathrm{Prop}}$ for $c_{\mathrm{e}1}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$ @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
        elseif Mode == 2
            ax.Title.String = ['$t_{\mathrm{Prop}}$ for $|\vec{c}_{\mathrm e}|$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$ @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
        end
    elseif Repetitions > 1
        if  Mode == 1
            ax.Title.String = ['$t_{\mathrm{Prop}}$ for $c_{\mathrm{e}1}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$ @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
        elseif Mode == 2
            ax.Title.String = ['$t_{\mathrm{Prop}}$ for $|\vec{c}_{\mathrm e}|$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$ @ ',num2str(Frequency),'\,kHz and ',num2str(Distance),'\,mm'];
        end
    end
end
z = ax.RTick;
z(1) = '';
ax.RTick = z;
ax.ThetaTick = [0 45 90 135 180 225 270 315];
ax.ThetaAxis.TickLabelFormat = '%g';
ax.TickLabelInterpreter = 'latex';
if  HeadLine > 0
    text(.59,1.04,'($\mu$s)','units','normalized','FontSize',FontSizeAxes*1.25,'interpreter','latex')
end
if  Export
    try
        if  PDF
            if  Mode == 1
                exportgraphics(f,fullfile(Directory,[FileName,'_tce1.pdf']),'ContentType','vector')
            elseif Mode == 2
                exportgraphics(f,fullfile(Directory,[FileName,'_tceAbs.pdf']),'ContentType','vector')
            end
        end
        if  PNG
            if  Mode == 1
                exportgraphics(f,fullfile(Directory,[FileName,'_tce1.png']),'Resolution',PNGresolution)
            elseif Mode == 2
                exportgraphics(f,fullfile(Directory,[FileName,'_tceAbs.png']),'Resolution',PNGresolution)
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
    if  Mode == 1
        output_txt = {['$\phi$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$t(c_{\mathrm{e}1}$): \textbf{',num2str(event_obj.Position(2),6),'} $\mu$s']};
    elseif Mode == 2
        output_txt = {['$\phi_\mathrm r$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$t(|\vec{c}_{\mathrm e}|$): \textbf{',num2str(event_obj.Position(2),6),'} $\mu$s']};
    end
end
end