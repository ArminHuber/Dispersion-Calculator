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
function ModeShapeLive(Hybrid,LayupString,PropagationAngle,FrequencyLimit,PhaseVelocityLimit,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,ALamb,AntisymmetricModes,AShear,AScholte,BLamb,BShear,BScholte,HigherOrderModes,LambModes,Material,PlateThickness,ShearHorizontalModes,ScholteModes,SLamb,SShear,SScholte,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,Decoupled,c,Delta,Layers,I,I1,Pattern,LayerThicknesses)
%#ok<*FXUP>
%#ok<*AGROW>

Mode = '';
Frequency = 0;
SamplesX3 = round(500/Layers);
ShowHalfSpaces = 0;
HalfSpaces = 1;
Phase = 0;

Plot = true(1,21);
Colors = [1 0 0; % u1 (1)
    .13 .55 .13; % u2 (2)
    0 0 1; % u3 (3)
    1 0 0; % sigma11 (4)
    .13 .55 .13; % sigma22 (5)
    0 0 1; % sigma33 (6)
    1 0 1; % sigma23 (7)
    0 0 0; % sigma13 (8)
    0 1 1; % sigma12 (9)
    1 0 0; % epsilon11 (10)
    .13 .55 .13; % epsilon22 (11)
    0 0 1; % epsilon33 (12)
    1 0 1; % epsilon23 (13)
    0 0 0; % epsilon13 (14)
    0 1 1; % epsilon12 (15)
    0 0 1; % Esrn (16)
    1 0 0; % Ekin (17)
    0 0 0; % Etot (18)
    1 0 0; % p1 (19)
    .13 .55 .13; % p2 (20)
    0 0 1]; % p3 (21)

SColor = [1 0 0];
AColor = [0 0 1];
BColor = [.5 0 1];

HeadLine = 2;

LineWidth = 1;
BoxLineWidth = .5;
FontSizeHeadLine = 24;
FontSizeAxesLabels = 14;
FontSizeAxes = 12;
FontSizeLegend = 12;

f = figure('Name','Live-through-thickness profiles','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Options','MenuSelectedFcn',@Options_Callback)
datacursormode on
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;

t = tiledlayout(3,3);
if  Hybrid
    Material{1}.Name = 'hybrid';
end
if  HeadLine > 0
    if  HeadLine == 1
        String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
    elseif ~SymmetricSystem && HeadLine == 2
        if  Repetitions == 1
            String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
        else
            String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        end
    elseif SymmetricSystem && HeadLine == 2
        if  Repetitions == 1
            String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
        else
            String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
        end
    end
    if  FluidLoading
        if  ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
        end
    end
    t.Title.String = String;
end
t.Title.FontSize = FontSizeHeadLine;
t.Title.Interpreter = 'latex';

ax1 = nexttile;
x = xline(0,'Color',[.6 .6 .6]);
hasbehavior(x,'legend',false);
hold on
u1 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(1,:));
u2 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(2,:));
u3 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(3,:));
ax1.Box = 'on';
ax1.LineWidth = BoxLineWidth;
ax1.FontSize = FontSizeAxes;
ax1.XLabel.FontSize = FontSizeAxesLabels;
ax1.YLabel.FontSize = FontSizeAxesLabels;
ax1.XLabel.String = 'Displacement (nm)';
ax1.YLabel.String = '$x_3$ (mm)';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';
ax1.TickLabelInterpreter = 'latex';
legend({'$u_1$','$u_2$','$u_3$'},'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax2 = nexttile;
x = xline(0,'Color',[.6 .6 .6]);
hasbehavior(x,'legend',false);
hold on
Esrn = plot(0,0,'LineWidth',LineWidth,'Color',Colors(16,:));
Ekin = plot(0,0,'LineWidth',LineWidth,'Color',Colors(17,:));
Etot = plot(0,0,'LineWidth',LineWidth,'Color',Colors(18,:));
ax2.Box = 'on';
ax2.LineWidth = BoxLineWidth;
ax2.FontSize = FontSizeAxes;
ax2.XLabel.FontSize = FontSizeAxesLabels;
ax2.YLabel.FontSize = FontSizeAxesLabels;
ax2.XLabel.String = 'Energy density (J/m$^2$)';
ax2.YLabel.String = '$x_3$ (mm)';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.TickLabelInterpreter = 'latex';
legend({'$E_\mathrm{srn}$','$E_\mathrm{kin}$','$E_\mathrm{tot}$'},'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax3 = nexttile;
x = xline(0,'Color',[.6 .6 .6]);
hasbehavior(x,'legend',false);
hold on
p1 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(19,:));
p2 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(20,:));
p3 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(21,:));
ax3.Box = 'on';
ax3.LineWidth = BoxLineWidth;
ax3.FontSize = FontSizeAxes;
ax3.XLabel.FontSize = FontSizeAxesLabels;
ax3.YLabel.FontSize = FontSizeAxesLabels;
ax3.XLabel.String = 'Power flow density (W/m)';
ax3.YLabel.String = '$x_3$ (mm)';
ax3.XLabel.Interpreter = 'latex';
ax3.YLabel.Interpreter = 'latex';
ax3.TickLabelInterpreter = 'latex';
legend({'$p_1$','$p_2$','$p_3$'},'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax4 = nexttile;
x = xline(0,'Color',[.6 .6 .6]);
hasbehavior(x,'legend',false);
hold on
sigma11 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(4,:));
sigma22 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(5,:));
sigma33 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(6,:));
sigma23 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(7,:));
sigma13 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(8,:));
sigma12 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(9,:));
ax4.Box = 'on';
ax4.LineWidth = BoxLineWidth;
ax4.FontSize = FontSizeAxes;
ax4.XLabel.FontSize = FontSizeAxesLabels;
ax4.YLabel.FontSize = FontSizeAxesLabels;
ax4.XLabel.String = 'Stress (kPa)';
ax4.YLabel.String = '$x_3$ (mm)';
ax4.XLabel.Interpreter = 'latex';
ax4.YLabel.Interpreter = 'latex';
ax4.TickLabelInterpreter = 'latex';
legend({'$\sigma_{11}$','$\sigma_{22}$','$\sigma_{33}$','$\sigma_{23}$','$\sigma_{13}$','$\sigma_{12}$'},'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax7 = nexttile(7);
x = xline(0,'Color',[.6 .6 .6]);
hasbehavior(x,'legend',false);
hold on
epsilon11 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(10,:));
epsilon22 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(11,:));
epsilon33 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(12,:));
epsilon23 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(13,:));
epsilon13 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(14,:));
epsilon12 = plot(0,0,'LineWidth',LineWidth,'Color',Colors(15,:));
ax7.Box = 'on';
ax7.LineWidth = BoxLineWidth;
ax7.FontSize = FontSizeAxes;
ax7.XLabel.FontSize = FontSizeAxesLabels;
ax7.YLabel.FontSize = FontSizeAxesLabels;
ax7.XLabel.String = 'Strain';
ax7.YLabel.String = '$x_3$ (mm)';
ax7.XLabel.Interpreter = 'latex';
ax7.YLabel.Interpreter = 'latex';
ax7.TickLabelInterpreter = 'latex';
legend({'$\varepsilon_{11}$','$\varepsilon_{22}$','$\varepsilon_{33}$','$\varepsilon_{23}$','$\varepsilon_{13}$','$\varepsilon_{12}$'},'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax5 = nexttile(5,[2 2]);
hold on
if  Symmetric
    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
        sLamb = plot(SLamb{1}(:,1),SLamb{1}(:,4),'LineWidth',LineWidth,'Color',SColor);
        if  ~Decoupled
            sLamb(2) = plot(SLamb{2}(:,1),SLamb{2}(:,4),'LineWidth',LineWidth,'Color',SColor);
            N = 3;
        else
            N = 2;
        end
        if  HigherOrderModes
            for i = N:size(SLamb,2)
                sLamb(i) = plot(SLamb{i}(:,1),SLamb{i}(:,4),'LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
        aLamb = plot(ALamb{1}(:,1),ALamb{1}(:,4),'LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(ALamb,2)
                aLamb(i) = plot(ALamb{i}(:,1),ALamb{i}(:,4),'LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
        sShear = plot(SShear{1}(:,1),SShear{1}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SShear,2)
                sShear(i) = plot(SShear{i}(:,1),SShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
        for i = 1:size(AShear,2)
            aShear(i) = plot(AShear{i}(:,1),AShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
        end
    end
    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
        sScholte = plot(SScholte{1}(:,1),SScholte{1}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SScholte,2)
                sScholte(i) = plot(SScholte{i}(:,1),SScholte{i}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
        aScholte = plot(AScholte{1}(:,1),AScholte{1}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(AScholte,2)
                aScholte(i) = plot(AScholte{i}(:,1),AScholte{i}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
else
    if  LambModes && ~isempty(BLamb{1})
        sLamb = plot(BLamb{1}(:,1),BLamb{1}(:,4),'LineWidth',LineWidth,'Color',BColor);
        sLamb(2) = plot(BLamb{2}(:,1),BLamb{2}(:,4),'LineWidth',LineWidth,'Color',BColor);
        if  ~Decoupled
            sLamb(3) = plot(BLamb{3}(:,1),BLamb{3}(:,4),'LineWidth',LineWidth,'Color',BColor);
            N = 4;
        else
            N = 3;
        end
        if  HigherOrderModes
            for i = N:size(BLamb,2)
                sLamb(i) = plot(BLamb{i}(:,1),BLamb{i}(:,4),'LineWidth',LineWidth,'Color',BColor);
            end
        end
    end
    if  ShearHorizontalModes
        if  SuperLayerSize > 1 && ~SymmetricSystem && ~isempty(BShear{1})
            sShear = plot(BShear{1}(:,1),BShear{1}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
            if  HigherOrderModes
                for i = 2:size(BShear,2)
                    sShear(i) = plot(BShear{i}(:,1),BShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                end
            end
        elseif SuperLayerSize == 1 || SymmetricSystem
            if  SymmetricModes && ~isempty(SShear{1})
                sShear = plot(SShear{1}(:,1),SShear{1}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                if  HigherOrderModes
                    for i = 2:size(SShear,2)
                        sShear(i) = plot(SShear{i}(:,1),SShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                    end
                end
            end
            if  AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                for i = 1:size(AShear,2)
                    aShear(i) = plot(AShear{i}(:,1),AShear{i}(:,4),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
                end
            end
        end
    end
    if  ScholteModes && ~isempty(BScholte{1})
        for i = 1:length(BScholte)
            sScholte(i) = plot(BScholte{i}(:,1),BScholte{i}(:,4),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
        end
    end
end
ax5.Box = 'on';
ax5.LineWidth = BoxLineWidth;
ax5.FontSize = FontSizeAxes;
ax5.XLabel.FontSize = FontSizeAxesLabels;
ax5.YLabel.FontSize = FontSizeAxesLabels;
ax5.XLabel.String = 'Frequency (kHz)';
ax5.YLabel.String = 'Phase velocity (m/ms)';
ax5.XLabel.Interpreter = 'latex';
ax5.YLabel.Interpreter = 'latex';
ax5.TickLabelInterpreter = 'latex';
ax5.XLim = [0 FrequencyLimit];
ax5.YLim = [0 PhaseVelocityLimit/1e3];
tb = axtoolbar('default');
tb.Visible = 'on';

Options_Callback

function output_txt = Cursor(~,event_obj)
    if  event_obj.Target.Parent.Layout.Tile == 1
        if  Phase
            if  event_obj.Target.Color == Colors(1,:)
                output_txt = {['$\varphi(u_1)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(2,:)
                output_txt = {['$\varphi(u_2)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(3,:)
                output_txt = {['$\varphi(u_3)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        else
            if  event_obj.Target.Color == Colors(1,:)
                output_txt = {['$u_1$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(2,:)
                output_txt = {['$u_2$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(3,:)
                output_txt = {['$u_3$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        end
    elseif event_obj.Target.Parent.Layout.Tile == 4
        if  Phase
            if  event_obj.Target.Color == Colors(6,:)
                output_txt = {['$\varphi(\sigma_{33})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(8,:)
                output_txt = {['$\varphi(\sigma_{13})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(7,:)
                output_txt = {['$\varphi(\sigma_{23})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(4,:)
                output_txt = {['$\varphi(\sigma_{11})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(5,:)
                output_txt = {['$\varphi(\sigma_{22})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(9,:)
                output_txt = {['$\varphi(\sigma_{12})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        else
            if  event_obj.Target.Color == Colors(6,:)
                output_txt = {['$\sigma_{33}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(8,:)
                output_txt = {['$\sigma_{13}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(7,:)
                output_txt = {['$\sigma_{23}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(4,:)
                output_txt = {['$\sigma_{11}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(5,:)
                output_txt = {['$\sigma_{22}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(9,:)
                output_txt = {['$\sigma_{12}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        end
    elseif event_obj.Target.Parent.Layout.Tile == 7
        if  Phase
            if  event_obj.Target.Color == Colors(12,:)
                output_txt = {['$\varphi(\varepsilon_{33})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(14,:)
                output_txt = {['$\varphi(\varepsilon_{13})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(13,:)
                output_txt = {['$\varphi(\varepsilon_{23})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(10,:)
                output_txt = {['$\varphi(\varepsilon_{11})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(11,:)
                output_txt = {['$\varphi(\varepsilon_{22})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(15,:)
                output_txt = {['$\varphi(\varepsilon_{12})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        else
            if  event_obj.Target.Color == Colors(12,:)
                output_txt = {['$\varepsilon_{33}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(14,:)
                output_txt = {['$\varepsilon_{13}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(13,:)
                output_txt = {['$\varepsilon_{23}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(10,:)
                output_txt = {['$\varepsilon_{11}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(11,:)
                output_txt = {['$\varepsilon_{22}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Colors(15,:)
                output_txt = {['$\varepsilon_{12}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        end
    elseif event_obj.Target.Parent.Layout.Tile == 2
        if  event_obj.Target.Color == Colors(16,:)
            output_txt = {['$E_\mathrm{strain}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
        elseif event_obj.Target.Color == Colors(17,:)
            output_txt = {['$E_\mathrm{kin}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
        elseif event_obj.Target.Color == Colors(18,:)
            output_txt = {['$E_\mathrm{total}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
        end
    elseif event_obj.Target.Parent.Layout.Tile == 3
        if  event_obj.Target.Color == Colors(19,:)
            output_txt = {['$p_1$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
        elseif event_obj.Target.Color == Colors(20,:)
            output_txt = {['$p_2$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
        elseif event_obj.Target.Color == Colors(21,:)
            output_txt = {['$p_3$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
        end
    elseif event_obj.Target.Parent.Layout.Tile == 5
        if  Symmetric
            if  ~Decoupled
                if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor)
                    for i = 1:length(SLamb)
                        if  i == 1 && event_obj.Target.YData(1) == SLamb{1}(1,4)
                            ModeName = 'S$_0$';
                            Mode = 'S0';
                            SecondaryData = SecondaryDataExtract(SLamb{1});
                            break
                        elseif i == 2 && event_obj.Target.YData(1) == SLamb{2}(1,4)
                            ModeName = 'S$_1$';
                            Mode = 'S1';
                            SecondaryData = SecondaryDataExtract(SLamb{2});
                            break
                        elseif i > 2 && (event_obj.Target.XData(1) == SLamb{i}(1,1) || event_obj.Target.XData(1) == SLamb{i}(1,2) || event_obj.Target.XData(1) == SLamb{i}(1,3))
                            ModeName = ['S$_{',num2str(i-1),'}$'];
                            Mode = ['S',num2str(i-1)];
                            SecondaryData = SecondaryDataExtract(SLamb{i});
                            break
                        end
                    end
                elseif ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor)
                    for i = 1:length(ALamb)
                        if  event_obj.Target.XData(1) == ALamb{i}(1,1) || event_obj.Target.XData(1) == ALamb{i}(1,2) || event_obj.Target.XData(1) == ALamb{i}(1,3)
                            ModeName = ['A$_{',num2str(i-1),'}$'];
                            Mode = ['A',num2str(i-1)];
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
                            Mode = ['S',num2str(i-1)];
                            SecondaryData = SecondaryDataExtract(SLamb{i});
                            break
                        end
                    end
                elseif ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(SShear)
                        if  event_obj.Target.XData(1) == SShear{i}(1,1) || event_obj.Target.XData(1) == SShear{i}(1,2) || event_obj.Target.XData(1) == SShear{i}(1,3)
                            ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                            Mode = ['SSH',num2str(i-1)];
                            SecondaryData = SecondaryDataExtract(SShear{i});
                            break
                        end
                    end
                elseif ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-')
                    for i = 1:length(ALamb)
                        if  event_obj.Target.XData(1) == ALamb{i}(1,1) || event_obj.Target.XData(1) == ALamb{i}(1,2) || event_obj.Target.XData(1) == ALamb{i}(1,3)
                            ModeName = ['A$_{',num2str(i-1),'}$'];
                            Mode = ['A',num2str(i-1)];
                            SecondaryData = SecondaryDataExtract(ALamb{i});
                            break
                        end
                    end
                elseif ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(AShear)
                        if  event_obj.Target.XData(1) == AShear{i}(1,1) || event_obj.Target.XData(1) == AShear{i}(1,2) || event_obj.Target.XData(1) == AShear{i}(1,3)
                            ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                            Mode = ['ASH',num2str(i)];
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
                        Mode = ['SScholte',num2str(i-1)];
                        SecondaryData = SecondaryDataExtract(SScholte{i});
                        break
                    end
                end
            elseif ~isempty(AScholte{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-.')
                for i = 1:length(AScholte)
                    if  event_obj.Target.XData(1) == AScholte{i}(1,1) || event_obj.Target.XData(1) == AScholte{i}(1,2) || event_obj.Target.XData(1) == AScholte{i}(1,3)
                        ModeName = ['A$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                        Mode = ['AScholte',num2str(i-1)];
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
                            Mode = 'B0';
                            SecondaryData = SecondaryDataExtract(BLamb{1});
                            break 
                        elseif i == 2 && event_obj.Target.YData(1) == BLamb{2}(1,4)
                            ModeName = 'B$_1$';
                            Mode = 'B1';
                            SecondaryData = SecondaryDataExtract(BLamb{2});
                            break
                        elseif i == 3 && event_obj.Target.YData(1) == BLamb{3}(1,4)
                            ModeName = 'B$_2$';
                            Mode = 'B2';
                            SecondaryData = SecondaryDataExtract(BLamb{3});
                            break
                        elseif i > 3 && (event_obj.Target.XData(1) == BLamb{i}(1,1) || event_obj.Target.XData(1) == BLamb{i}(1,2) || event_obj.Target.XData(1) == BLamb{i}(1,3))
                            ModeName = ['B$_{',num2str(i-1),'}$'];
                            Mode = ['B',num2str(i-1)];
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
                            Mode = 'B0';
                            SecondaryData = SecondaryDataExtract(BLamb{1});
                            break
                        elseif i == 2 && event_obj.Target.YData(1) == BLamb{2}(1,4)
                            ModeName = 'B$_1$';
                            Mode = 'B1';
                            SecondaryData = SecondaryDataExtract(BLamb{2});
                            break
                        elseif i > 2 && (event_obj.Target.XData(1) == BLamb{i}(1,1) || event_obj.Target.XData(1) == BLamb{i}(1,2) || event_obj.Target.XData(1) == BLamb{i}(1,3))
                            ModeName = ['B$_{',num2str(i-1),'}$'];
                            Mode = ['B',num2str(i-1)];
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
                                Mode = ['BSH',num2str(i-1)];
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
                                Mode = ['SSH',num2str(i-1)];
                                SecondaryData = SecondaryDataExtract(SShear{i});
                                break
                            end
                        end
                    elseif ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
                        for i = 1:length(AShear)
                            if  event_obj.Target.XData(1) == AShear{i}(1,1) || event_obj.Target.XData(1) == AShear{i}(1,2) || event_obj.Target.XData(1) == AShear{i}(1,3)
                                ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                                Mode = ['ASH',num2str(i)];
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
                        Mode = 'BScholte0';
                        SecondaryData = SecondaryDataExtract(BScholte{1});
                        break
                    elseif i == 2 && event_obj.Target.YData(1) == BScholte{2}(1,4)
                        ModeName = 'B$^{\mathrm{Scholte}}_1$';
                        Mode = 'BScholte1';
                        SecondaryData = SecondaryDataExtract(BScholte{2});
                        break
                    elseif i > 2 && (event_obj.Target.XData(1) == BScholte{i}(1,1) || event_obj.Target.XData(1) == BScholte{i}(1,2) || event_obj.Target.XData(1) == BScholte{i}(1,3))
                        ModeName = ['B$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                        Mode = ['BScholte',num2str(i-1)];
                        SecondaryData = SecondaryDataExtract(BScholte{i});
                        break
                    end
                end
            end
        end
        if  ~Decoupled
            output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'] ['$c_{\mathrm{e}1}$: \textbf{',num2str(SecondaryData(1),6),'} m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'} Np/m']};
        else
            output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'] ['$c_{\mathrm{e}}$: \textbf{',num2str(SecondaryData(1),6),'} m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'} Np/m']};
        end
        Frequency = event_obj.Position(1);
        UpdateModeShapes
    end
    function SecondaryData = SecondaryDataExtract(Data)
        z = find(event_obj.Position(2) == Data(:,4));
        if  ~isscalar(z)
            z = find(event_obj.Position(1) == Data(:,1));
        end
        SecondaryData = Data(z,[5 7]);
    end
end
function Options_Callback(~,~)
    f_Options = figure('NumberTitle','off','Name','Options','Visible','off','MenuBar','none','Position',[0 0 210 235],'CloseRequestFcn',@CloseRequest);
    f_Options.Units = 'normalized';
    movegui(f_Options,'east')
    f_Options.Visible = 'on';
    drawnow
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.fHG2Client.getWindow.setAlwaysOnTop(true)
    uicontrol('Parent',f_Options,'Style','text','String','Samples x3','Position',[10 205 58 13]);
    uicontrol('Parent',f_Options,'Style','text','String','Half-spaces','Position',[10 175 61 13]);
    uicontrol('Parent',f_Options,'Style','text','String','Phase','Position',[10 145 32 13]);
    uicontrol('Parent',f_Options,'Style','edit','String',SamplesX3,'TooltipString',['Enter the number of sample points over the plate''s thickness (x3)',newline,'at which the selected quantities are calculated.'],'Position',[110 200 50 23],'Callback',@Callback,'Tag','1');
    ShowHalfSpacesUI = uicontrol('Parent',f_Options,'Style','checkbox','Value',ShowHalfSpaces,'TooltipString','Check this to show the quantities in the upper and lower fluid.','Position',[80 170 20 23],'Callback',@Callback,'Tag','2');        
    HalfSpacesUI = uicontrol('Parent',f_Options,'Style','edit','String',HalfSpaces,'TooltipString','Set the height of the half-spaces in plate thicknesses.','Position',[110 170 50 23],'Callback',@Callback,'Tag','3');
    uicontrol('Parent',f_Options,'Style','checkbox','Value',Phase,'TooltipString','Check this to plot the phase of the field components.','Position',[80 140 20 23],'Callback',@Callback,'Tag','4');
    uicontrol('Parent',f_Options,'Style','checkbox','String','u1','Value',Plot(1),'Position',[10 110 35 23],'Callback',@Callback,'Tag','5');
    uicontrol('Parent',f_Options,'Style','checkbox','String','u2','Value',Plot(2),'Position',[10 90 35 23],'Callback',@Callback,'Tag','6');
    uicontrol('Parent',f_Options,'Style','checkbox','String','u3','Value',Plot(3),'Position',[10 70 35 23],'Callback',@Callback,'Tag','7');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(963),'11'],'Value',Plot(4),'Position',[55 110 50 23],'Callback',@Callback,'Tag','8');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(963),'22'],'Value',Plot(5),'Position',[55 90 50 23],'Callback',@Callback,'Tag','9');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(963),'33'],'Value',Plot(6),'Position',[55 70 50 23],'Callback',@Callback,'Tag','10');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(963),'23'],'Value',Plot(7),'Position',[55 50 50 23],'Callback',@Callback,'Tag','11');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(963),'13'],'Value',Plot(8),'Position',[55 30 50 23],'Callback',@Callback,'Tag','12');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(963),'12'],'Value',Plot(9),'Position',[55 10 50 23],'Callback',@Callback,'Tag','13');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(949),'11'],'Value',Plot(10),'Position',[105 110 50 23],'Callback',@Callback,'Tag','14');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(949),'22'],'Value',Plot(11),'Position',[105 90 50 23],'Callback',@Callback,'Tag','15');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(949),'33'],'Value',Plot(12),'Position',[105 70 50 23],'Callback',@Callback,'Tag','16');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(949),'23'],'Value',Plot(13),'Position',[105 50 50 23],'Callback',@Callback,'Tag','17');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(949),'13'],'Value',Plot(14),'Position',[105 30 50 23],'Callback',@Callback,'Tag','18');
    uicontrol('Parent',f_Options,'Style','checkbox','String',[char(949),'12'],'Value',Plot(15),'Position',[105 10 50 23],'Callback',@Callback,'Tag','19');
    uicontrol('Parent',f_Options,'Style','checkbox','String','Esrn','Value',Plot(16),'Position',[155 110 60 23],'Callback',@Callback,'Tag','20');
    uicontrol('Parent',f_Options,'Style','checkbox','String','Ekin','Value',Plot(17),'Position',[155 90 60 23],'Callback',@Callback,'Tag','21');
    uicontrol('Parent',f_Options,'Style','checkbox','String','Etot','Value',Plot(18),'Position',[155 70 60 23],'Callback',@Callback,'Tag','22');
    uicontrol('Parent',f_Options,'Style','checkbox','String','p1','Value',Plot(19),'Position',[155 50 35 23],'Callback',@Callback,'Tag','23');
    uicontrol('Parent',f_Options,'Style','checkbox','String','p2','Value',Plot(20),'Position',[155 30 35 23],'Callback',@Callback,'Tag','24');
    uicontrol('Parent',f_Options,'Style','checkbox','String','p3','Value',Plot(21),'Position',[155 10 35 23],'Callback',@Callback,'Tag','25');
    if  FluidLoading
        if  ShowHalfSpaces
            HalfSpacesUI.Enable = 'on';
        else
            HalfSpacesUI.Enable = 'off';
        end
        ShowHalfSpacesUI.Enable = 'on';
    else
        HalfSpacesUI.Enable = 'off';
        ShowHalfSpacesUI.Enable = 'off';
    end
    function Callback(source,~)
        if  strcmp(source.Tag,'1') % Samples x3
            SamplesX3 = str2double(source.String);
        elseif strcmp(source.Tag,'2') % Half-spaces (on/off)
            ShowHalfSpaces = source.Value;
            if  ShowHalfSpaces
                HalfSpacesUI.Enable = 'on';
            else
                HalfSpacesUI.Enable = 'off';
            end
        elseif strcmp(source.Tag,'3') % Half-spaces (number of)
            HalfSpaces = str2double(source.String);
        elseif strcmp(source.Tag,'4') % Phase
            Phase = source.Value;
            if  Phase
                ax1.XLabel.String = 'Displacement phase ($^\circ$)';
                ax2.XLabel.String = 'Energy density phase ($^\circ$)';
                ax3.XLabel.String = 'Power flow density phase ($^\circ$)';
                ax4.XLabel.String = 'Stress phase ($^\circ$)';
                ax7.XLabel.String = 'Strain phase ($^\circ$)';
            else
                ax1.XLabel.String = 'Displacement (nm)';
                ax2.XLabel.String = 'Energy density (J/m$^2$)';
                ax3.XLabel.String = 'Power flow density (W/m)';
                ax4.XLabel.String = 'Stress (kPa)';
                ax7.XLabel.String = 'Strain';
            end
        elseif strcmp(source.Tag,'5')
            Plot(1) = source.Value;
        elseif strcmp(source.Tag,'6')
            Plot(2) = source.Value;
        elseif strcmp(source.Tag,'7')
            Plot(3) = source.Value;
        elseif strcmp(source.Tag,'8')
            Plot(4) = source.Value;
        elseif strcmp(source.Tag,'9')
            Plot(5) = source.Value;
        elseif strcmp(source.Tag,'10')
            Plot(6) = source.Value;
        elseif strcmp(source.Tag,'11')
            Plot(7) = source.Value;
        elseif strcmp(source.Tag,'12')
            Plot(8) = source.Value;
        elseif strcmp(source.Tag,'13')
            Plot(9) = source.Value;
        elseif strcmp(source.Tag,'14')
            Plot(10) = source.Value;
        elseif strcmp(source.Tag,'15')
            Plot(11) = source.Value;
        elseif strcmp(source.Tag,'16')
            Plot(12) = source.Value;
        elseif strcmp(source.Tag,'17')
            Plot(13) = source.Value;
        elseif strcmp(source.Tag,'18')
            Plot(14) = source.Value;
        elseif strcmp(source.Tag,'19')
            Plot(15) = source.Value;
        elseif strcmp(source.Tag,'20')
            Plot(16) = source.Value;
        elseif strcmp(source.Tag,'21')
            Plot(17) = source.Value;
        elseif strcmp(source.Tag,'22')
            Plot(18) = source.Value;
        elseif strcmp(source.Tag,'23')
            Plot(19) = source.Value;
        elseif strcmp(source.Tag,'24')
            Plot(20) = source.Value;
        elseif strcmp(source.Tag,'25')
            Plot(21) = source.Value;
        end
        UpdateModeShapes
    end
    function CloseRequest(~,~)
        delete(f_Options)
    end
end
function UpdateModeShapes
    [u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,x3Total,~] = ModeShapeLinesComputer_Anisotropic(FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,ALamb,AShear,AScholte,BLamb,BShear,BScholte,SLamb,SShear,SScholte,c,Delta,Material,Layers,Frequency,Mode,PlateThickness,SamplesX3,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,ShowHalfSpaces,HalfSpaces,Phase);
    if  Plot(1)
        u1.YData = -x3Total*1e3;
        if  Phase
            u1.XData = uPhase(:,1);
        else
            u1.XData = real(u(:,1))*1e9;
        end
    else
        u1.XData = 0;
        u1.YData = 0;
    end
    if  Plot(2)
        u2.YData = -x3Total*1e3;
        if  Phase
            u2.XData = uPhase(:,2);
        else
            u2.XData = real(u(:,2))*1e9;
        end
    else
        u2.XData = 0;
        u2.YData = 0;
    end
    if  Plot(3)
        u3.YData = -x3Total*1e3;
        if  Phase
            u3.XData = uPhase(:,3);
        else
            u3.XData = real(u(:,3))*1e9;
        end
    else
        u3.XData = 0;
        u3.YData = 0;
    end
    ax1.XLimMode = 'auto';
    ax1.YLimMode = 'auto';
    ax1.XLim = max(abs(ax1.XLim))*[-1 1];
    ax1.YLim = -1e3*[x3Total(end) x3Total(1)];
    if  Plot(16)
        Esrn.YData = -x3Total*1e3;
        if  Phase
            Esrn.XData = zeros(length(x3Total),1);
        else
            Esrn.XData = StrainEnergyDensity;
        end
    else
        Esrn.XData = 0;
        Esrn.YData = 0;
    end
    if  Plot(17)
        Ekin.YData = -x3Total*1e3;
        if  Phase
            Ekin.XData = zeros(length(x3Total),1);
        else
            Ekin.XData = KineticEnergyDensity;
        end
    else
        Ekin.XData = 0;
        Ekin.YData = 0;
    end
    if  Plot(18)
        Etot.YData = -x3Total*1e3;
        if  Phase
            Etot.XData = zeros(length(x3Total),1);
        else
            Etot.XData = TotalEnergyDensity;
        end
    else
        Etot.XData = 0;
        Etot.YData = 0;
    end
    ax2.XLimMode = 'auto';
    ax2.YLimMode = 'auto';
    ax2.XLim = max(abs(ax2.XLim))*[-1 1];
    ax2.YLim = -1e3*[x3Total(end) x3Total(1)];
    if  Plot(19)
        p1.YData = -x3Total*1e3;
        if  Phase
            p1.XData = zeros(length(x3Total),1);
        else
            p1.XData = PowerFlowDensity(:,1);
        end
    else
        p1.XData = 0;
        p1.YData = 0;
    end
    if  Plot(20)
        p2.YData = -x3Total*1e3;
        if  Phase
            p2.XData = zeros(length(x3Total),1);
        else
            p2.XData = PowerFlowDensity(:,2);
        end
    else
        p2.XData = 0;
        p2.YData = 0;
    end    
    if  Plot(21)
        p3.YData = -x3Total*1e3;
        if  Phase
            p3.XData = zeros(length(x3Total),1);
        else
            p3.XData = PowerFlowDensity(:,3);
        end
    else
        p3.XData = 0;
        p3.YData = 0;
    end
    ax3.XLimMode = 'auto';
    ax3.YLimMode = 'auto';
    ax3.XLim = max(abs(ax3.XLim))*[-1 1];
    ax3.YLim = -1e3*[x3Total(end) x3Total(1)];
    if  Plot(4)
        sigma11.YData = -x3Total*1e3;
        if  Phase
            sigma11.XData = sigmaPhase(:,1);
        else
            sigma11.XData = real(sigma(:,1))/1e3;
        end
    else
        sigma11.XData = 0;
        sigma11.YData = 0;
    end
    if  Plot(5)
        sigma22.YData = -x3Total*1e3;
        if  Phase
            sigma22.XData = sigmaPhase(:,2);
        else
            sigma22.XData = real(sigma(:,2))/1e3;
        end
    else
        sigma22.XData = 0;
        sigma22.YData = 0;
    end
    if  Plot(6)
        sigma33.YData = -x3Total*1e3;
        if  Phase
            sigma33.XData = sigmaPhase(:,3);
        else
            sigma33.XData = real(sigma(:,3))/1e3;
        end
    else
        sigma33.XData = 0;
        sigma33.YData = 0;
    end
    if  Plot(7)
        sigma23.YData = -x3Total*1e3;
        if  Phase
            sigma23.XData = sigmaPhase(:,4);
        else
            sigma23.XData = real(sigma(:,4))/1e3;
        end
    else
        sigma23.XData = 0;
        sigma23.YData = 0;
    end
    if  Plot(8)
        sigma13.YData = -x3Total*1e3;
        if  Phase
            sigma13.XData = sigmaPhase(:,5);
        else
            sigma13.XData = real(sigma(:,5))/1e3;
        end
    else
        sigma13.XData = 0;
        sigma13.YData = 0;
    end
    if  Plot(9)
        sigma12.YData = -x3Total*1e3;
        if  Phase
            sigma12.XData = sigmaPhase(:,6);
        else
            sigma12.XData = real(sigma(:,6))/1e3;
        end
    else
        sigma12.XData = 0;
        sigma12.YData = 0;
    end
    ax4.XLimMode = 'auto';
    ax4.YLimMode = 'auto';
    ax4.XLim = max(abs(ax4.XLim))*[-1 1];
    ax4.YLim = -1e3*[x3Total(end) x3Total(1)];
    if  Plot(10)
        epsilon11.YData = -x3Total*1e3;
        if  Phase
            epsilon11.XData = epsilonPhase(:,1);
        else
            epsilon11.XData = real(epsilon(:,1));
        end
    else
        epsilon11.XData = 0;
        epsilon11.YData = 0;
    end
    if  Plot(11)
        epsilon22.YData = -x3Total*1e3;
        if  Phase
            epsilon22.XData = epsilonPhase(:,2);
        else
            epsilon22.XData = real(epsilon(:,2));
        end
    else
        epsilon22.XData = 0;
        epsilon22.YData = 0;
    end
    if  Plot(12)
        epsilon33.YData = -x3Total*1e3;
        if  Phase
            epsilon33.XData = epsilonPhase(:,3);
        else
            epsilon33.XData = real(epsilon(:,3));
        end
    else
        epsilon33.XData = 0;
        epsilon33.YData = 0;
    end
    if  Plot(13)
        epsilon23.YData = -x3Total*1e3;
        if  Phase
            epsilon23.XData = epsilonPhase(:,4);
        else
            epsilon23.XData = real(epsilon(:,4));
        end
    else
        epsilon23.XData = 0;
        epsilon23.YData = 0;
    end
    if  Plot(14)
        epsilon13.YData = -x3Total*1e3;
        if  Phase
            epsilon13.XData = epsilonPhase(:,5);
        else
            epsilon13.XData = real(epsilon(:,5));
        end
    else
        epsilon13.XData = 0;
        epsilon13.YData = 0;
    end
    if  Plot(15)
        epsilon12.YData = -x3Total*1e3;
        if  Phase
            epsilon12.XData = epsilonPhase(:,6);
        else
            epsilon12.XData = real(epsilon(:,6));
        end
    else
        epsilon12.XData = 0;
        epsilon12.YData = 0;
    end
    ax7.XLimMode = 'auto';
    ax7.YLimMode = 'auto';
    ax7.XLim = max(abs(ax7.XLim))*[-1 1];
    ax7.YLim = -1e3*[x3Total(end) x3Total(1)];
end
end