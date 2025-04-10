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
function PhaseVelocity_Anisotropic(Hybrid,Phi,LayupString,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,BulkVelocities,PNGresolution,SColor,AColor,BColor,ALamb,AntisymmetricModes,AShear,AScholte,BLamb,BoxLineWidth,BShear,BScholte,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,HeadLine,HigherOrderModes,LambModes,LineWidth,Material,PDF,FileName,PlateThickness,PNG,PropagationAngle,ShearHorizontalModes,ScholteModes,SLamb,SShear,SScholte,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,XAxis,XAxisMode,YAxis,Decoupled)
%#ok<*FXUP>
%#ok<*AGROW>

% sLamb = SLamb BLamb
% aLamb = ALamb CLamb

% sShear = SShear BShear CShear
% aShear = AShear

% sScholte = SScholte BScholte
% aScholte = AScholte

CLamb{1}=[];CShear{1}=[];
f = figure('Name','Phase velocity','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
hold on
if  Symmetric
    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
        sLamb = plot(SLamb{1}(:,XAxisMode),SLamb{1}(:,4),'LineWidth',LineWidth,'Color',SColor);
        if  ~Decoupled
            sLamb(2) = plot(SLamb{2}(:,XAxisMode),SLamb{2}(:,4),'LineWidth',LineWidth,'Color',SColor);
            N = 3;
        else
            N = 2;
        end
        if  HigherOrderModes
            for i = N:size(SLamb,2)
                sLamb(i) = plot(SLamb{i}(:,XAxisMode),SLamb{i}(:,4),'LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
        aLamb = plot(ALamb{1}(:,XAxisMode),ALamb{1}(:,4),'LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(ALamb,2)
                aLamb(i) = plot(ALamb{i}(:,XAxisMode),ALamb{i}(:,4),'LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
        sShear = plot(SShear{1}(:,XAxisMode),SShear{1}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SShear,2)
                sShear(i) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
        for i = 1:size(AShear,2)
            aShear(i) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
        end
    end
    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
        sScholte = plot(SScholte{1}(:,XAxisMode),SScholte{1}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SScholte,2)
                sScholte(i) = plot(SScholte{i}(:,XAxisMode),SScholte{i}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
        aScholte = plot(AScholte{1}(:,XAxisMode),AScholte{1}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(AScholte,2)
                aScholte(i) = plot(AScholte{i}(:,XAxisMode),AScholte{i}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
else
    if  LambModes && ~isempty(BLamb{1})
        sLamb = plot(BLamb{1}(:,XAxisMode),BLamb{1}(:,4),'LineWidth',LineWidth,'Color',BColor);
        sLamb(2) = plot(BLamb{2}(:,XAxisMode),BLamb{2}(:,4),'LineWidth',LineWidth,'Color',BColor);
        if  ~Decoupled
            sLamb(3) = plot(BLamb{3}(:,XAxisMode),BLamb{3}(:,4),'LineWidth',LineWidth,'Color',BColor);
            N = 4;
        else
            N = 3;
        end
        if  HigherOrderModes
            for i = N:size(BLamb,2)
                sLamb(i) = plot(BLamb{i}(:,XAxisMode),BLamb{i}(:,4),'LineWidth',LineWidth,'Color',BColor);
            end
        end
    end
    if  ShearHorizontalModes
        if  SuperLayerSize > 1 && ~SymmetricSystem && ~isempty(BShear{1})
            sShear = plot(BShear{1}(:,XAxisMode),BShear{1}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
            if  HigherOrderModes
                for i = 2:size(BShear,2)
                    sShear(i) = plot(BShear{i}(:,XAxisMode),BShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                end
            end
        elseif SuperLayerSize == 1 || SymmetricSystem
            if  SymmetricModes && ~isempty(SShear{1})
                sShear = plot(SShear{1}(:,XAxisMode),SShear{1}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                if  HigherOrderModes
                    for i = 2:size(SShear,2)
                        sShear(i) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                    end
                end
            end
            if  AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                for i = 1:size(AShear,2)
                    aShear(i) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
                end
            end
        end
    end
    if  ScholteModes && ~isempty(BScholte{1})
        for i = 1:length(BScholte)
            sScholte(i) = plot(BScholte{i}(:,XAxisMode),BScholte{i}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
        end
    end
end
if  SuperLayerSize == 1 && BulkVelocities
    n = [cosd(Phi) sind(Phi)];
    C = Material{1}.C;
    A1 = (-C(1,1)*n(1)^2-C(2,2)*n(2)^2-C(4,4)*n(2)^2-C(5,5)*n(1)^2-C(6,6)*n(1)^2-C(6,6)*n(2)^2)/Material{1}.Density;
    A2 = (C(1,1)*C(5,5)*n(1)^4+C(2,2)*C(4,4)*n(2)^4+C(1,1)*C(6,6)*n(1)^4+C(2,2)*C(6,6)*n(2)^4+C(4,4)*C(6,6)*n(2)^4+C(5,5)*C(6,6)*n(1)^4-C(1,2)^2*n(1)^2*n(2)^2+C(1,1)*C(2,2)*n(1)^2*n(2)^2+C(1,1)*C(4,4)*n(1)^2*n(2)^2+C(2,2)*C(5,5)*n(1)^2*n(2)^2-2*C(1,2)*C(6,6)*n(1)^2*n(2)^2+C(4,4)*C(6,6)*n(1)^2*n(2)^2+C(5,5)*C(6,6)*n(1)^2*n(2)^2)/Material{1}.Density^2;
    A3 = (C(1,2)^2*C(4,4)*n(1)^2*n(2)^4+C(1,2)^2*C(5,5)*n(1)^4*n(2)^2-C(1,1)*C(5,5)*C(6,6)*n(1)^6-C(2,2)*C(4,4)*C(6,6)*n(2)^6-C(1,1)*C(2,2)*C(4,4)*n(1)^2*n(2)^4-C(1,1)*C(2,2)*C(5,5)*n(1)^4*n(2)^2-C(1,1)*C(4,4)*C(6,6)*n(1)^4*n(2)^2+2*C(1,2)*C(4,4)*C(6,6)*n(1)^2*n(2)^4+2*C(1,2)*C(5,5)*C(6,6)*n(1)^4*n(2)^2-C(2,2)*C(5,5)*C(6,6)*n(1)^2*n(2)^4)/Material{1}.Density^3;
    S_fast = abs(real((-A1^2/9+A2/3)/(2*(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))-(((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)/2-(3^(1/2)*((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)+(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))*1i)/2-A1/3)^(1/2)); % fast quasi shear wave, S_fast
    S_slow = abs(real((-A1^2/9+A2/3)/(2*(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))-(((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)/2+(3^(1/2)*((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)+(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))*1i)/2-A1/3)^(1/2)); % slow shear wave, S_slow
    L = abs(real(-((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)-A1/3-(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))^(1/2))); % quasi longitudinal wave, L
    plot(XAxis,[S_fast S_fast]/1e3,'LineWidth',LineWidth,'Color','c');
    plot(XAxis,[S_slow S_slow]/1e3,'LineWidth',LineWidth,'Color','c');
    plot(XAxis,[L L]/1e3,'LineWidth',LineWidth,'Color','c');
    if  FluidLoading
        if  Symmetric
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
if  Hybrid
    Material{1}.Name = 'hybrid';
end
if  HeadLine == 1
    if  XAxisMode ~= 3
        String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
    else
        String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_')];
    end
elseif ~SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  XAxisMode ~= 3
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
        else
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
        end
    elseif Repetitions > 1
        if  XAxisMode ~= 3
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        else
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        end
    end
elseif SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  XAxisMode ~= 3
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
        else
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
        end
    elseif Repetitions > 1
        if  XAxisMode ~= 3
            String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
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
if  XAxisMode == 1
    ax.XLabel.String = 'Frequency (kHz)';
    ax.XLim = XAxis;
elseif XAxisMode == 2
    ax.XLabel.String = 'Frequency (MHz)';
    ax.XLim = XAxis/1e3;
elseif XAxisMode == 3
    ax.XLabel.String = 'Frequency$\cdot$thickness (MHz$\cdot$mm)';
    ax.XLim = XAxis*PlateThickness;
end
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
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
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj)
    if  Symmetric
        if  ~Decoupled
            if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor)
                for i = 1:length(SLamb)
                    if  i == 1 && event_obj.Target.YData(1) == SLamb{1}(1,4)
                        ModeName = 'S$_0$';
                        SecondaryData = SecondaryDataExtract(SLamb{1});
                        break
                    elseif i == 2 && event_obj.Target.YData(1) == SLamb{2}(1,4)
                        ModeName = 'S$_1$';
                        SecondaryData = SecondaryDataExtract(SLamb{2});
                        break
                    elseif i > 2 && (event_obj.Target.XData(1) == SLamb{i}(1,1) || event_obj.Target.XData(1) == SLamb{i}(1,2) || event_obj.Target.XData(1) == SLamb{i}(1,3))
                        ModeName = ['S$_{',num2str(i-1),'}$'];
                        SecondaryData = SecondaryDataExtract(SLamb{i});
                        break
                    end
                end
            elseif ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor)
                for i = 1:length(ALamb)
                    if  event_obj.Target.XData(1) == ALamb{i}(1,1) || event_obj.Target.XData(1) == ALamb{i}(1,2) || event_obj.Target.XData(1) == ALamb{i}(1,3)
                        ModeName = ['A$_{',num2str(i-1),'}$'];
                        SecondaryData = SecondaryDataExtract(ALamb{i});
                        break
                    end
                end
            end
        else
            if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(SLamb)
                    if  event_obj.Target.XData(1) == SLamb{i}(1,1) || event_obj.Target.XData(1) == SLamb{i}(1,2) || event_obj.Target.XData(1) == SLamb{i}(1,3)
                        ModeName = ['S$_{',num2str(i-1),'}$'];
                        SecondaryData = SecondaryDataExtract(SLamb{i});
                        break
                    end
                end
            elseif ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
                for i = 1:length(SShear)
                    if  event_obj.Target.XData(1) == SShear{i}(1,1) || event_obj.Target.XData(1) == SShear{i}(1,2) || event_obj.Target.XData(1) == SShear{i}(1,3)
                        ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                        SecondaryData = SecondaryDataExtract(SShear{i});
                        break
                    end
                end
            elseif ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(ALamb)
                    if  event_obj.Target.XData(1) == ALamb{i}(1,1) || event_obj.Target.XData(1) == ALamb{i}(1,2) || event_obj.Target.XData(1) == ALamb{i}(1,3)
                        ModeName = ['A$_{',num2str(i-1),'}$'];
                        SecondaryData = SecondaryDataExtract(ALamb{i});
                        break
                    end
                end
            elseif ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
                for i = 1:length(AShear)
                    if  event_obj.Target.XData(1) == AShear{i}(1,1) || event_obj.Target.XData(1) == AShear{i}(1,2) || event_obj.Target.XData(1) == AShear{i}(1,3)
                        ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                        SecondaryData = SecondaryDataExtract(AShear{i});
                        break
                    end
                end
            end
        end
        if  ~isempty(SScholte{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(SScholte)
                if  event_obj.Target.XData(1) == SScholte{i}(1,1) || event_obj.Target.XData(1) == SScholte{i}(1,2) || event_obj.Target.XData(1) == SScholte{i}(1,3)
                    ModeName = ['S$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    SecondaryData = SecondaryDataExtract(SScholte{i});
                    break
                end
            end
        elseif ~isempty(AScholte{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(AScholte)
                if  event_obj.Target.XData(1) == AScholte{i}(1,1) || event_obj.Target.XData(1) == AScholte{i}(1,2) || event_obj.Target.XData(1) == AScholte{i}(1,3)
                    ModeName = ['A$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    SecondaryData = SecondaryDataExtract(AScholte{i});
                    break
                end
            end
        end
    else
        if  ~Decoupled
            if  ~isempty(BLamb{1})
                for i = 1:length(BLamb)
                    if  i == 1 && event_obj.Target.YData(1) == BLamb{1}(1,4)
                        ModeName = 'B$_0$';
                        SecondaryData = SecondaryDataExtract(BLamb{1});
                        break 
                    elseif i == 2 && event_obj.Target.YData(1) == BLamb{2}(1,4)
                        ModeName = 'B$_1$';
                        SecondaryData = SecondaryDataExtract(BLamb{2});
                        break
                    elseif i == 3 && event_obj.Target.YData(1) == BLamb{3}(1,4)
                        ModeName = 'B$_2$';
                        SecondaryData = SecondaryDataExtract(BLamb{3});
                        break
                    elseif i > 3 && (event_obj.Target.XData(1) == BLamb{i}(1,1) || event_obj.Target.XData(1) == BLamb{i}(1,2) || event_obj.Target.XData(1) == BLamb{i}(1,3))
                        ModeName = ['B$_{',num2str(i-1),'}$'];
                        SecondaryData = SecondaryDataExtract(BLamb{i});
                        break
                    end
                end
            end
        else
            if  ~isempty(BLamb{1}) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(BLamb)
                    if  i == 1 && event_obj.Target.YData(1) == BLamb{1}(1,4)
                        ModeName = 'B$_0$';
                        SecondaryData = SecondaryDataExtract(BLamb{1});
                        break
                    elseif i == 2 && event_obj.Target.YData(1) == BLamb{2}(1,4)
                        ModeName = 'B$_1$';
                        SecondaryData = SecondaryDataExtract(BLamb{2});
                        break
                    elseif i > 2 && (event_obj.Target.XData(1) == BLamb{i}(1,1) || event_obj.Target.XData(1) == BLamb{i}(1,2) || event_obj.Target.XData(1) == BLamb{i}(1,3))
                        ModeName = ['B$_{',num2str(i-1),'}$'];
                        SecondaryData = SecondaryDataExtract(BLamb{i});
                        break
                    end
                end
            end
            if  SuperLayerSize > 1 && ~SymmetricSystem
                if  ~isempty(BShear{1}) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(BShear)
                        if  event_obj.Target.XData(1) == BShear{i}(1,1) || event_obj.Target.XData(1) == BShear{i}(1,2) || event_obj.Target.XData(1) == BShear{i}(1,3)
                            ModeName = ['B$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                            SecondaryData = SecondaryDataExtract(BShear{i});
                            break
                        end
                    end
                end
            elseif SuperLayerSize == 1 || SymmetricSystem
                if  ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(SShear)
                        if  event_obj.Target.XData(1) == SShear{i}(1,1) || event_obj.Target.XData(1) == SShear{i}(1,2) || event_obj.Target.XData(1) == SShear{i}(1,3)
                            ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                            SecondaryData = SecondaryDataExtract(SShear{i});
                            break
                        end
                    end
                elseif ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(AShear)
                        if  event_obj.Target.XData(1) == AShear{i}(1,1) || event_obj.Target.XData(1) == AShear{i}(1,2) || event_obj.Target.XData(1) == AShear{i}(1,3)
                            ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                            SecondaryData = SecondaryDataExtract(AShear{i});
                            break
                        end
                    end
                end
            end
        end
        if  ~isempty(BScholte{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(BScholte)
                if  i == 1 && event_obj.Target.YData(1) == BScholte{1}(1,4)
                    ModeName = 'B$^{\mathrm{Scholte}}_0$';
                    SecondaryData = SecondaryDataExtract(BScholte{1});
                    break
                elseif i == 2 && event_obj.Target.YData(1) == BScholte{2}(1,4)
                    ModeName = 'B$^{\mathrm{Scholte}}_1$';
                    SecondaryData = SecondaryDataExtract(BScholte{2});
                    break
                elseif i > 2 && (event_obj.Target.XData(1) == BScholte{i}(1,1) || event_obj.Target.XData(1) == BScholte{i}(1,2) || event_obj.Target.XData(1) == BScholte{i}(1,3))
                    ModeName = ['B$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    SecondaryData = SecondaryDataExtract(BScholte{i});
                    break
                end
            end
        end
    end
    if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        if  ~Decoupled
            output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'] ['$c_{\mathrm{e}1}$: \textbf{',num2str(SecondaryData(1),6),'} m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'} Np/m']};
        else
            output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'] ['$c_{\mathrm{e}}$: \textbf{',num2str(SecondaryData(1),6),'} m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'} Np/m']};
        end
    elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        if  ~Decoupled
            output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} MHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'] ['$c_{\mathrm{e}1}$: \textbf{',num2str(SecondaryData(1),6),'} m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'} Np/m']};
        else
            output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} MHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'] ['$c_{\mathrm{e}}$: \textbf{',num2str(SecondaryData(1),6),'} m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'} Np/m']};
        end
    elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        if  ~Decoupled
            output_txt = {['\textbf{',ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'} MHz$\cdot$mm'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'] ['$c_{\mathrm{e}1}$: \textbf{',num2str(SecondaryData(1),6),'} m/ms'] ['$\alpha d$: \textbf{',num2str(SecondaryData(2),6),'} Np$\cdot$mm/m']};
        else
            output_txt = {['\textbf{',ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'} MHz$\cdot$mm'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'] ['$c_{\mathrm{e}}$: \textbf{',num2str(SecondaryData(1),6),'} m/ms'] ['$\alpha d$: \textbf{',num2str(SecondaryData(2),6),'} Np$\cdot$mm/m']};
        end
    elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == S_fast/1e3
        if  ~Decoupled
            output_txt = {'\textbf{Fast quasi shear velocity}' ['$v_{\mathrm S_{\mathrm{fast}}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
        else
            output_txt = {'\textbf{Fast shear velocity}' ['$v_{\mathrm S_{\mathrm{fast}}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
        end
    elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == S_slow/1e3
        if  Decoupled && PropagationAngle == 0 && strcmp(Material{1}.Class,'Transversely isotropic')
            output_txt = {'\textbf{Fast \& slow shear velocity}' ['$v_{\mathrm SLamb}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
        elseif Decoupled && ((PropagationAngle == 90 && strcmp(Material{1}.Class,'Transversely isotropic')) || strcmp(Material{1}.Class,'Orthotropic'))
            output_txt = {'\textbf{Slow shear velocity}' ['$v_{\mathrm S_{\mathrm{slow}}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
        else
            output_txt = {'\textbf{Slow quasi shear velocity}' ['$v_{\mathrm S_{\mathrm{slow}}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
        end
    elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == L/1e3
        if  ~Decoupled
            output_txt = {'\textbf{Quasi longitudinal velocity}' ['$v_{\mathrm L}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
        else
            output_txt = {'\textbf{Longitudinal velocity}' ['$v_{\mathrm L}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
        end
    elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == UpperFluid.Velocity/1e3
        output_txt = {['\textbf{',replace(UpperFluid.Name,'_','\_'),'}'] ['$v$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == LowerFluid.Velocity/1e3
        output_txt = {['\textbf{',replace(LowerFluid.Name,'_','\_'),'}'] ['$v$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end
    function SecondaryData = SecondaryDataExtract(Data)
        z = find(event_obj.Position(2) == Data(:,4));
        if  ~isscalar(z)
            if  XAxisMode == 1
                z = find(event_obj.Position(1) == Data(:,1));
            elseif XAxisMode == 2
                z = find(event_obj.Position(1) == Data(:,2));
            elseif XAxisMode == 3
                z = find(event_obj.Position(1) == Data(:,3));
            end
        end
        SecondaryData(1) = Data(z,5);
        if  XAxisMode == 3
            SecondaryData(2) = Data(z,7)*PlateThickness*1e3;
        else
            SecondaryData(2) = Data(z,7);
        end
    end
end
function ShowModes_Callback(~,~)
    ShowModes_Anisotropic(SuperLayerSize,SymmetricSystem,Symmetric,Decoupled,f.Children(end).Children,SColor,AColor)
end
function Analyze_Callback(~,~)
    f_Analyze = figure('NumberTitle','off','Name','Analyze','Visible','off','MenuBar','none','Position',[0 0 210 60],'CloseRequestFcn',@CloseRequest);
    f_Analyze.Units = 'normalized';
    movegui(f_Analyze,'center')
    f_Analyze.Visible = 'on';
    drawnow
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.fHG2Client.getWindow.setAlwaysOnTop(true)
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
        Value = .1*XAxis(2)*PlateThickness;
        String = {['f',char(8901),'d (MHz',char(8901),'mm)'],'Phase velocity (m/ms)','Wavelength/d (mm/mm)',['Wavenumber',char(8901),'d (rad/mm',char(8901),'mm)']};
    end
    uicontrol('Parent',f_Analyze,'Style','popupmenu','Value',Mode,'TooltipString','Select constant quantity.','String',String,'Position',[10 20 130 23],'Callback',@Mode_Callback);
    ValueUI = uicontrol('Parent',f_Analyze,'Style','edit','String',Value,'TooltipString','Enter a value for the above selected quantity.','Position',[150 20 50 23],'Callback',@Value_Callback);
    function Mode_Callback(source,~)
        Mode = source.Value;
        switch source.Value
        case 1 % frequency
            if  XAxisMode == 1 % (kHz)
                Value = .1*XAxis(2);
            elseif XAxisMode == 2 % (MHz)
                Value = .1*XAxis(2)/1e3;
            elseif XAxisMode == 3 % (MHz*mm)
                Value = .1*XAxis(2)*PlateThickness; 
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
        if  Mode == 1
            if  LambModes && ((AntisymmetricModes && ~isempty(ALamb{1})) || ~isempty(CLamb{1}))
                DatatipPlacer_Frequency(aLamb)
            end
            if  LambModes && ((SymmetricModes && ~isempty(SLamb{1})) || ~isempty(BLamb{1}))
                DatatipPlacer_Frequency(sLamb)
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                DatatipPlacer_Frequency(aShear)
            end
            if  ShearHorizontalModes && ((SymmetricModes && ~isempty(SShear{1})) || ~isempty(BShear{1}) || ~isempty(CShear{1}))
                DatatipPlacer_Frequency(sShear)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_Frequency(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_Frequency(sScholte)
            end
        elseif Mode == 2
            if  LambModes && ((AntisymmetricModes && ~isempty(ALamb{1})) || ~isempty(CLamb{1}))
                DatatipPlacer_PhaseVelocity(aLamb)
            end
            if  LambModes && ((SymmetricModes && ~isempty(SLamb{1})) || ~isempty(BLamb{1}))
                DatatipPlacer_PhaseVelocity(sLamb)
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                DatatipPlacer_PhaseVelocity(aShear)
            end
            if  ShearHorizontalModes && ((SymmetricModes && ~isempty(SShear{1})) || ~isempty(BShear{1}) || ~isempty(CShear{1}))
                DatatipPlacer_PhaseVelocity(sShear)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_PhaseVelocity(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_PhaseVelocity(sScholte)
            end
        elseif Mode == 3
            if  LambModes && ((AntisymmetricModes && ~isempty(ALamb{1})) || ~isempty(CLamb{1}))
                DatatipPlacer_Wavelength(aLamb)
            end
            if  LambModes && ((SymmetricModes && ~isempty(SLamb{1})) || ~isempty(BLamb{1}))
                DatatipPlacer_Wavelength(sLamb)
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                DatatipPlacer_Wavelength(aShear)
            end
            if  ShearHorizontalModes && ((SymmetricModes && ~isempty(SShear{1})) || ~isempty(BShear{1}) || ~isempty(CShear{1}))
                DatatipPlacer_Wavelength(sShear)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_Wavelength(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_Wavelength(sScholte)
            end
        elseif Mode == 4
            if  LambModes && ((AntisymmetricModes && ~isempty(ALamb{1})) || ~isempty(CLamb{1}))
                DatatipPlacer_Wavenumber(aLamb)
            end
            if  LambModes && ((SymmetricModes && ~isempty(SLamb{1})) || ~isempty(BLamb{1}))
                DatatipPlacer_Wavenumber(sLamb)
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                DatatipPlacer_Wavenumber(aShear)
            end
            if  ShearHorizontalModes && ((SymmetricModes && ~isempty(SShear{1})) || ~isempty(BShear{1}) || ~isempty(CShear{1}))
                DatatipPlacer_Wavenumber(sShear)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_Wavenumber(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_Wavenumber(sScholte)
            end
        end
        function DatatipPlacer_Frequency(Data)
            for i = 1:length(Data)
                if  Data(i).XData(1) < Value
                    if  Data(i).XData(end) >= Value
                        [~,q] = min(abs(Data(i).XData-Value));
                        dt(end+1) = datatip(Data(i),Data(i).XData(q),Data(i).YData(q));
                    end
                else
                    break
                end
            end
        end
        function DatatipPlacer_PhaseVelocity(Data)
            for i = 1:length(Data)
                if  min(Data(i).YData) < Value && max(Data(i).YData) > Value
                    [~,q] = min(abs(Data(i).YData-Value));
                    dt(end+1) = datatip(Data(i),Data(i).XData(q),Data(i).YData(q));
                end
            end
        end
        function DatatipPlacer_Wavelength(Data)
            for i = 1:length(Data)
                Wavelength = Data(i).YData./Data(i).XData;
                if  XAxisMode == 1 % (kHz)
                    Wavelength = Wavelength*1e3;
                end
                if  min(Wavelength) < Value && max(Wavelength) > Value
                    [~,q] = min(abs(Wavelength-Value));
                    dt(end+1) = datatip(Data(i),Data(i).XData(q),Data(i).YData(q));
                end
            end
        end
        function DatatipPlacer_Wavenumber(Data)
            for i = 1:length(Data)
                Wavenumber = 2*pi*Data(i).XData./Data(i).YData;
                if  XAxisMode == 1 % (kHz)
                    Wavenumber = Wavenumber/1e3;
                end
                if  min(Wavenumber) < Value && max(Wavenumber) > Value
                    [~,q] = min(abs(Wavenumber-Value));
                    dt(end+1) = datatip(Data(i),Data(i).XData(q),Data(i).YData(q));
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