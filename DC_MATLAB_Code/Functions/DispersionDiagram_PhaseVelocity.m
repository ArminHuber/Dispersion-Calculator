% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function DispersionDiagram_PhaseVelocity(Geometry,Hybrid,Phi,LayupString,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,BulkVelocities,PNGresolution,SColor,AColor,BColor,ALamb,AntisymmetricModes,AShear,BLamb,BoxLineWidth,BShear,CLamb,CShear,FlexuralModes,FlexuralModeOrders,LongitudinalModes,TorsionalModes,F,L,T,LineColors,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,HeadLine,LambModes,LineWidth,Material,PDF,FileName,Thickness,ThicknessInner,PNG,PropagationAngle,ShearHorizontalModes,SLamb,SShear,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,XAxis,XAxisMode,YAxis,Decoupled)
%#ok<*FXUP>
%#ok<*AGROW>
if  isempty(LayupString)
    MaterialType = 1;
else
    MaterialType = 2;
end
f = figure('Icon',which('DC_Logo16.png'),'Name','Phase velocity','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
hold on
Line = line().empty;
if  strcmp(Geometry,'Plate')
    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
        for i = 1:size(SLamb,2)
            Line(end+1) = plot(SLamb{i}(:,XAxisMode),SLamb{i}(:,4),'LineWidth',LineWidth,'Color',SColor);
            Line(end).UserData.ModeName = ['S$_{',num2str(i-1),'}$'];
            Line(end).UserData.FullData = SLamb{i}(:,[5 7]);
        end
    end
    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
        for i = 1:size(ALamb,2)
            Line(end+1) = plot(ALamb{i}(:,XAxisMode),ALamb{i}(:,4),'LineWidth',LineWidth,'Color',AColor);
            Line(end).UserData.ModeName = ['A$_{',num2str(i-1),'}$'];
            Line(end).UserData.FullData = ALamb{i}(:,[5 7]);
        end
    end
    if  LambModes && ~isempty(BLamb{1})
        for i = 1:size(BLamb,2)
            Line(end+1) = plot(BLamb{i}(:,XAxisMode),BLamb{i}(:,4),'LineWidth',LineWidth,'Color',BColor);
            Line(end).UserData.ModeName = ['B$_{',num2str(i-1),'}$'];
            Line(end).UserData.FullData = BLamb{i}(:,[5 7]);
        end
    end    
    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
        for i = 1:size(SShear,2)
            Line(end+1) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            Line(end).UserData.ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
            Line(end).UserData.FullData = SShear{i}(:,[5 7]);
        end
    end
    if  ShearHorizontalModes && AntisymmetricModes && ~isempty(AShear{1})
        for i = 1:size(AShear,2)
            Line(end+1) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
            Line(end).UserData.ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
            Line(end).UserData.FullData = AShear{i}(:,[5 7]);
        end
    end
    if  ShearHorizontalModes && ~isempty(BShear{1})
        for i = 1:size(BShear,2)
            Line(end+1) = plot(BShear{i}(:,XAxisMode),BShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
            Line(end).UserData.ModeName = ['B$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
            Line(end).UserData.FullData = BShear{i}(:,[5 7]);
        end
    end
elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
    if  LongitudinalModes && ~isempty(L{1})
        for i = 1:length(L)
            Line(end+1) = plot(L{i}(:,XAxisMode),L{i}(:,4),'LineWidth',LineWidth,'Color',SColor);
            Line(end).UserData.ModeName = ['L(0,',num2str(i),')'];
            Line(end).UserData.FullData = L{i}(:,[5 7]);
        end
    end
    if  FlexuralModes && ~isempty(F{1})
        LineColorsIndex = 0;
        if  FlexuralModeOrders > length(F)
            FlexuralModeOrders = length(F);
        end
        for n = 1:FlexuralModeOrders
            LineColorsIndex = LineColorsIndex+1;
            for i = 1:length(F{n})
                Line(end+1) = plot(F{n}{i}(:,XAxisMode),F{n}{i}(:,4),'LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));
                Line(end).UserData.ModeName = ['F(',num2str(n),',',num2str(i),')'];
                Line(end).UserData.FullData = F{n}{i}(:,[5 7]);
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    end
    if  TorsionalModes && ~isempty(T{1})
        for i = 1:length(T)
            Line(end+1) = plot(T{i}(:,XAxisMode),T{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            Line(end).UserData.ModeName = ['T(0,',num2str(i),')'];
            Line(end).UserData.FullData = T{i}(:,[5 7]);
        end
    end
elseif strcmp(Geometry,'Circumferential')
    if  LambModes && ~isempty(CLamb{1})
        for i = 1:size(CLamb,2)
            Line(end+1) = plot(CLamb{i}(:,XAxisMode),CLamb{i}(:,4),'LineWidth',LineWidth,'Color',BColor);
            Line(end).UserData.ModeName = ['C$_{',num2str(i-1),'}$'];
            Line(end).UserData.FullData = CLamb{i}(:,[5 7]);
        end
    end
    if  ShearHorizontalModes && ~isempty(CShear{1})
        for i = 1:size(CShear,2)
            Line(end+1) = plot(CShear{i}(:,XAxisMode),CShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
            Line(end).UserData.ModeName = ['C$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
            Line(end).UserData.FullData = CShear{i}(:,[5 7]);
        end
    end
end
delete(f.Children(end).Children(end))
if  BulkVelocities
    if  MaterialType == 1
        plot(XAxis,[Material.LongitudinalVelocity Material.LongitudinalVelocity]/1e3,'LineWidth',LineWidth,'Color','c');
        plot(XAxis,[Material.TransverseVelocity Material.TransverseVelocity]/1e3,'LineWidth',LineWidth,'Color','c');
    elseif MaterialType == 2 && SuperLayerSize == 1
        n = [cosd(Phi) sind(Phi)];
        C = real(Material{1}.C);
        A1 = (-C(1,1)*n(1)^2-C(2,2)*n(2)^2-C(4,4)*n(2)^2-C(5,5)*n(1)^2-C(6,6)*n(1)^2-C(6,6)*n(2)^2)/Material{1}.Density;
        A2 = (C(1,1)*C(5,5)*n(1)^4+C(2,2)*C(4,4)*n(2)^4+C(1,1)*C(6,6)*n(1)^4+C(2,2)*C(6,6)*n(2)^4+C(4,4)*C(6,6)*n(2)^4+C(5,5)*C(6,6)*n(1)^4-C(1,2)^2*n(1)^2*n(2)^2+C(1,1)*C(2,2)*n(1)^2*n(2)^2+C(1,1)*C(4,4)*n(1)^2*n(2)^2+C(2,2)*C(5,5)*n(1)^2*n(2)^2-2*C(1,2)*C(6,6)*n(1)^2*n(2)^2+C(4,4)*C(6,6)*n(1)^2*n(2)^2+C(5,5)*C(6,6)*n(1)^2*n(2)^2)/Material{1}.Density^2;
        A3 = (C(1,2)^2*C(4,4)*n(1)^2*n(2)^4+C(1,2)^2*C(5,5)*n(1)^4*n(2)^2-C(1,1)*C(5,5)*C(6,6)*n(1)^6-C(2,2)*C(4,4)*C(6,6)*n(2)^6-C(1,1)*C(2,2)*C(4,4)*n(1)^2*n(2)^4-C(1,1)*C(2,2)*C(5,5)*n(1)^4*n(2)^2-C(1,1)*C(4,4)*C(6,6)*n(1)^4*n(2)^2+2*C(1,2)*C(4,4)*C(6,6)*n(1)^2*n(2)^4+2*C(1,2)*C(5,5)*C(6,6)*n(1)^4*n(2)^2-C(2,2)*C(5,5)*C(6,6)*n(1)^2*n(2)^4)/Material{1}.Density^3;
        d1 = A1/3;
        d2 = A2/3-d1^2;
        d3 = d1^3-d1*A2/2+A3/2;
        d4 = (sqrt(d2^3+d3^2)-d3)^(1/3);
        d5 = d2/d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        S_fast = abs(real(sqrt(d6+d7))); % phase velocity in the solid (m/s)
        S_slow = abs(real(sqrt(d6-d7)));
        L = abs(real(sqrt(d4-d5-d1)));
        plot(XAxis,[S_fast S_fast]/1e3,'LineWidth',LineWidth,'Color','c');
        plot(XAxis,[S_slow S_slow]/1e3,'LineWidth',LineWidth,'Color','c');
        plot(XAxis,[L L]/1e3,'LineWidth',LineWidth,'Color','c');
    end
    if  FluidLoading && (MaterialType == 1 || (MaterialType == 2 && SuperLayerSize == 1))
        if  Symmetric || strcmp(Geometry,'Rod')
            plot(XAxis,[UpperFluid.Velocity UpperFluid.Velocity]/1e3,'LineWidth',LineWidth,'Color','c');
        else
            if  ToggleUpperFluid
                plot(XAxis,[UpperFluid.Velocity UpperFluid.Velocity]/1e3,'LineWidth',LineWidth,'Color','c');
            end
            if  ToggleLowerFluid
                plot(XAxis,[LowerFluid.Velocity LowerFluid.Velocity]/1e3,'LineWidth',LineWidth,'Color','c');
            end
        end
    end
end
ax = gca;
ax.Box = 'on';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
if  MaterialType == 1
    if  HeadLine
        if  XAxisMode ~= 3
            if  strcmp(Geometry,'Plate')
                String = ['Dispersion diagram of ',num2str(Thickness*1e3),'\,mm ',replace(Material.Name,'_','\_'),' plate'];
            elseif strcmp(Geometry,'Rod')
                String = ['Dispersion diagram of ',num2str(Thickness*1e3),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
            elseif strcmp(Geometry,'Pipe')
                String = ['Dispersion diagram of ',num2str(Thickness*1e3),'\,$\times$\,',num2str((Thickness-ThicknessInner)*1e3/2),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
            elseif strcmp(Geometry,'Circumferential')
                String = ['Dispersion diagram of ',num2str(Thickness*1e3),'\,$\times$\,',num2str((Thickness-ThicknessInner)*1e3/2),'\,mm ',replace(Material.Name,'_','\_'),' circumference'];
            end
        else
            if  strcmp(Geometry,'Plate')
                String = ['Dispersion diagram of ',replace(Material.Name,'_','\_'),' plate'];
            elseif strcmp(Geometry,'Rod')
                String = ['Dispersion diagram of ',replace(Material.Name,'_','\_'),' rod'];
            elseif strcmp(Geometry,'Pipe')
                String = ['Dispersion diagram of $d_\mathrm{i}$/$d_\mathrm{o}=$\,',num2str(ThicknessInner/Thickness),' ',replace(Material.Name,'_','\_'),' pipe'];
            elseif strcmp(Geometry,'Circumferential')
                String = ['Dispersion diagram of $d_\mathrm{i}$/$d_\mathrm{o}=$\,',num2str(ThicknessInner/Thickness),' ',replace(Material.Name,'_','\_'),' circumference'];
            end
        end
        if  FluidLoading
            if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Pipe')
                if  ToggleUpperFluid && ToggleLowerFluid
                    String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
                end
            elseif strcmp(Geometry,'Rod')
                String = append(String,' in ',replace(UpperFluid.Name,'_','\_'));
            end
        end
        ax.Title.String = String;
    end
else
    if  Hybrid
        Material{1}.Name = 'hybrid';
    end
    if  HeadLine == 1
        if  XAxisMode ~= 3
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
        else
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_')];
        end
    elseif ~SymmetricSystem && HeadLine == 2
        if  Repetitions == 1
            if  XAxisMode ~= 3
                String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
            else
                String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
            end
        elseif Repetitions > 1
            if  XAxisMode ~= 3
                String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
            else
                String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
            end
        end
    elseif SymmetricSystem && HeadLine == 2
        if  Repetitions == 1
            if  XAxisMode ~= 3
                String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
            else
                String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
            end
        elseif Repetitions > 1
            if  XAxisMode ~= 3
                String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
            else
                String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
            end
        end
    end
    if  HeadLine > 0
        if  FluidLoading
            if  ToggleUpperFluid && ToggleLowerFluid
                String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
            end
        end
        ax.Title.String = String;
    end
end
if  ~Decoupled
    ce = '$c_{\mathrm e1}$';
else
    ce = '$c_{\mathrm e}$';
end
if  XAxisMode == 1
    ax.XLabel.String = 'Frequency (kHz)';
    ax.XLim = XAxis;
elseif XAxisMode == 2
    ax.XLabel.String = 'Frequency (MHz)';
    ax.XLim = XAxis/1e3;
elseif XAxisMode == 3
    if  strcmp(Geometry,'Plate')
        ax.XLabel.String = 'Frequency$\cdot$thickness (MHz$\cdot$mm)';
        ax.XLim = XAxis*Thickness;
    elseif strcmp(Geometry,'Rod')
        ax.XLabel.String = 'Frequency$\cdot$diameter (MHz$\cdot$mm)';
        ax.XLim = XAxis*Thickness;
    elseif strcmp(Geometry,'Pipe') || strcmp(Geometry,'Circumferential')
        ax.XLabel.String = 'Frequency$\cdot$wall thickness (MHz$\cdot$mm)';
        ax.XLim = XAxis/2*(Thickness-ThicknessInner);
    end
end
ax.YLabel.Interpreter = 'latex';
if  strcmp(Geometry,'Circumferential')
    ax.YLabel.String = 'Phase velocity @ $d_\mathrm{o}$ (m/ms)';
else
    ax.YLabel.String = 'Phase velocity (m/ms)';
end
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.YLim = YAxis;
ax.TickLabelInterpreter = 'latex';
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
datacursormode on
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj)
    if  length(event_obj.Target.XData) ~= 2
        SecondaryData = event_obj.Target.UserData.FullData(event_obj.Target.Children(1).DataIndex,:);
        if  XAxisMode == 3
            if  strcmp(Geometry,'Plate')
                SecondaryData(2) = SecondaryData(2)*Thickness*1e3;
            elseif strcmp(Geometry,'Circumferential')
                SecondaryData(2) = SecondaryData(2)*(Thickness-ThicknessInner)*1e3/2;
            end
        end
    end
    if  MaterialType == 1
        if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'] [ce,': \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'}\,Np/m']};
        elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'] [ce,': \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'}\,Np/m']};
        elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz$\cdot$mm'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'] [ce,': \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] ['$\alpha d$: \textbf{',num2str(SecondaryData(2),6),'}\,Np$\cdot$mm/m']};
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == Material.LongitudinalVelocity/1e3
            output_txt = {'\textbf{Longitudinal velocity}' ['$v_{\mathrm L}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == Material.TransverseVelocity/1e3
            output_txt = {'\textbf{Transverse velocity}' ['$v_{\mathrm T}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == UpperFluid.Velocity/1e3
            output_txt = {['\textbf{',replace(UpperFluid.Name,'_','\_'),'}'] ['$v$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == LowerFluid.Velocity/1e3
            output_txt = {['\textbf{',replace(LowerFluid.Name,'_','\_'),'}'] ['$v$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
            output_txt = 'Marker';
        end
    else
        if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'] [ce,': \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'}\,Np/m']};
        elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'] [ce,': \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'}\,Np/m']};
        elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz$\cdot$mm'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'] [ce,': \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] ['$\alpha d$: \textbf{',num2str(SecondaryData(2),6),'}\,Np$\cdot$mm/m']};
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == S_fast/1e3
            if  ~Decoupled
                output_txt = {'\textbf{Fast quasi shear velocity}' ['$v_{\mathrm S_{\mathrm{fast}}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
            else
                output_txt = {'\textbf{Fast shear velocity}' ['$v_{\mathrm S_{\mathrm{fast}}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
            end
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == S_slow/1e3
            if  Decoupled && PropagationAngle == 0 && strcmp(Material{1}.Class,'Transversely isotropic')
                output_txt = {'\textbf{Fast \& slow shear velocity}' ['$v_{\mathrm S}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
            elseif Decoupled && ((PropagationAngle == 90 && strcmp(Material{1}.Class,'Transversely isotropic')) || strcmp(Material{1}.Class,'Orthotropic'))
                output_txt = {'\textbf{Slow shear velocity}' ['$v_{\mathrm S_{\mathrm{slow}}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
            else
                output_txt = {'\textbf{Slow quasi shear velocity}' ['$v_{\mathrm S_{\mathrm{slow}}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
            end
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == L/1e3
            if  ~Decoupled
                output_txt = {'\textbf{Quasi longitudinal velocity}' ['$v_{\mathrm L}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
            else
                output_txt = {'\textbf{Longitudinal velocity}' ['$v_{\mathrm L}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
            end
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == UpperFluid.Velocity/1e3
            output_txt = {['\textbf{',replace(UpperFluid.Name,'_','\_'),'}'] ['$v$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == LowerFluid.Velocity/1e3
            output_txt = {['\textbf{',replace(LowerFluid.Name,'_','\_'),'}'] ['$v$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
            output_txt = 'Marker';
        end
    end
end
function ShowModes_Callback(~,~)
    if  MaterialType == 1
        if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Circumferential')
            ShowModes_Isotropic(Geometry,Symmetric,f.Children(end).Children,SColor,AColor,BColor)
        else
            ShowModes_Isotropic_Rod_Pipe(f.Children(end).Children,SColor,AColor)
        end
    else
        ShowModes_Anisotropic(SuperLayerSize,SymmetricSystem,Symmetric,Decoupled,f.Children(end).Children,SColor,AColor)
    end
end
function Analyze_Callback(~,~)
    f_Analyze = figure('Icon',which('DC_Logo16.png'),'WindowStyle','alwaysontop','NumberTitle','off','Name','Analyze','Visible','off','MenuBar','none','Position',[0 0 210 140],'CloseRequestFcn',@CloseRequest);
    f_Analyze.Units = 'normalized';
    movegui(f_Analyze,'center')
    f_Analyze.Visible = 'on';
    drawnow
    l = line(ax,0,0,'Color','k');
    dt = line(ax,0,0,'Color','k');
    Mode = 1;
    if  XAxisMode == 1 % (kHz)
        Value = .1*XAxis(2);
        String = {'Frequency (kHz)','Phase velocity (m/ms)','Wavelength (mm)','Wavenumber (rad/mm)'};
    elseif XAxisMode == 2 % (MHz)
        Value = .1*XAxis(2)/1e3;
        String = {'Frequency (MHz)','Phase velocity (m/ms)','Wavelength (mm)','Wavenumber (rad/mm)'};
    elseif XAxisMode == 3 % (MHz*mm)
        if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod')
            Value = .1*XAxis(2)*Thickness;
        elseif strcmp(Geometry,'Pipe') || strcmp(Geometry,'Circumferential')
            Value = .1*XAxis(2)/2*(Thickness-ThicknessInner);
        end
        String = {['f',char(8901),'d (MHz',char(8901),'mm)'],'Phase velocity (m/ms)','Wavelength/d (mm/mm)',['Wavenumber',char(8901),'d (rad/mm',char(8901),'mm)']};
    end
    uicontrol('Parent',f_Analyze,'Style','popupmenu','Value',Mode,'Tooltip','Select constant quantity.','String',String,'Position',[10 100 130 23],'Callback',@Mode_Callback);
    ValueUI = uicontrol('Parent',f_Analyze,'Style','edit','String',Value,'Tooltip','Enter a value for the above selected quantity.','Position',[150 100 50 23],'Callback',@Value_Callback);
    function Mode_Callback(source,~)
        Mode = source.Value;
        switch source.Value
        case 1 % frequency
            if  XAxisMode == 1 % (kHz)
                Value = .1*XAxis(2);
            elseif XAxisMode == 2 % (MHz)
                Value = .1*XAxis(2)/1e3;
            elseif XAxisMode == 3 % (MHz*mm)
                if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod')
                    Value = .1*XAxis(2)*Thickness;
                elseif strcmp(Geometry,'Pipe') || strcmp(Geometry,'Circumferential')
                    Value = .1*XAxis(2)/2*(Thickness-ThicknessInner);
                end
            end
            ValueUI.String = Value;
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        case 2 % phase velocity
            Value = .1*YAxis(2);
            ValueUI.String = Value;
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        case 3 % wavelength
            Value = 1;
            ValueUI.String = Value;
            if  XAxisMode == 1 % (kHz)
                XMax = YAxis(2)/Value*1e3;
            else
                XMax = YAxis(2)/Value;
            end
            XData = [0 XMax];
            YData = [0 YAxis(2)];
        case 4 % wavenumber
            Value = 1;
            ValueUI.String = Value;
            if  XAxisMode == 1 % (kHz)
                XMax = YAxis(2)*Value/(2*pi)*1e3;
            else
                XMax = YAxis(2)*Value/(2*pi);
            end
            XData = [0 XMax];
            YData = [0 YAxis(2)];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function Value_Callback(source,~)
        Value = str2double(source.String);
        if  Mode == 1 % frequency
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        elseif Mode == 2 % phase velocity
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        elseif Mode == 3 % wavelength
            if  XAxisMode == 1 % (kHz)
                XMax = YAxis(2)/Value*1e3;
            else
                XMax = YAxis(2)/Value;
            end
            XData = [0 XMax];
            YData = [0 YAxis(2)];
         elseif Mode == 4 % wavenumber
            if  XAxisMode == 1 % (kHz)
                XMax = YAxis(2)*Value/(2*pi)*1e3;
            else
                XMax = YAxis(2)*Value/(2*pi);
            end
            XData = [0 XMax];
            YData = [0 YAxis(2)];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function ListModes
        delete(dt)
        if  Mode == 1 % frequency
            for i = 1:length(Line)
                if  Line(i).XData(1) < Value && Line(i).XData(end) >= Value
                    [~,q] = min(abs(Line(i).XData-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
                end
            end
        elseif Mode == 2 % phase velocity
            for i = 1:length(Line)
                if  min(Line(i).YData) < Value && max(Line(i).YData) > Value
                    [~,q] = min(abs(Line(i).YData-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
                end
            end
        elseif Mode == 3 % wavelength
            for i = 1:length(Line)
                Wavelength = Line(i).YData./Line(i).XData;
                if  XAxisMode == 1 % (kHz)
                    Wavelength = Wavelength*1e3;
                end
                if  min(Wavelength) < Value && max(Wavelength) > Value
                    [~,q] = min(abs(Wavelength-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
                end
            end
        elseif Mode == 4 % wavenumber
            for i = 1:length(Line)
                Wavenumber = 2*pi*Line(i).XData./Line(i).YData;
                if  XAxisMode == 1 % (kHz)
                    Wavenumber = Wavenumber/1e3;
                end
                if  min(Wavenumber) < Value && max(Wavenumber) > Value
                    [~,q] = min(abs(Wavenumber-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
                end
            end
        end
    end
    function CloseRequest(~,~)
        delete(l)
        delete(dt)
        delete(f_Analyze)
    end
end
end