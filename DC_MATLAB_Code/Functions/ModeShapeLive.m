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
function ModeShapeLive(Geometry,Hybrid,LayupString,PropagationAngle,FrequencyLimit,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Sink,ALamb,AntisymmetricModes,AShear,BLamb,BShear,CLamb,CShear,FlexuralModes,FlexuralModeOrders,LongitudinalModes,TorsionalModes,F,L,T,LineColors,LambModes,Material,Thickness,ThicknessInner,ShearHorizontalModes,SLamb,SShear,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Decoupled,c,Delta,Layers,I,I1,Pattern,LayerThicknesses)
%#ok<*AGROW>
Ro = Thickness/2;
Ri = ThicknessInner/2;
Mode = '';
Frequency = 0;
SamplesX3 = round(500/Layers);
ShowHalfSpaces = false;
HalfSpaces = 1;
Phase = false;
DispersionDiagrams = 1; % 1:cp-ce 2:cp-alpha 3:ce-alpha
YAxesEnable = false;
YAxisText1UIString = 'cp-axis (m/ms)';
YAxisText2UIString = 'ce-axis (m/ms)';
YAxis_cp = 0;
YAxis_ce = 0;
YAxis_alpha = 0;
Quantity = [4 5]; % 4:cp 5:ce1 7:alpha

Plot = true(1,21);
Colors = [1 0 0 % u1 (1)
    .13 .55 .13 % u2 (2)
    0 0 1 % u3 (3)
    1 0 0 % sigma11 (4)
    .13 .55 .13 % sigma22 (5)
    0 0 1 % sigma33 (6)
    1 0 1 % sigma23 (7)
    0 0 0 % sigma13 (8)
    0 1 1 % sigma12 (9)
    1 0 0 % epsilon11 (10)
    .13 .55 .13 % epsilon22 (11)
    0 0 1 % epsilon33 (12)
    1 0 1 % epsilon23 (13)
    0 0 0 % epsilon13 (14)
    0 1 1 % epsilon12 (15)
    0 0 1 % Esrn (16)
    1 0 0 % Ekin (17)
    0 0 0 % Etot (18)
    1 0 0 % p1 (19)
    .13 .55 .13 % p2 (20)
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

if  strcmp(Geometry,'Plate')
    YLabelString = '$x_3$';
elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
    YLabelString = '$r$';
elseif strcmp(Geometry,'Circumferential')
    YLabelString = '$d$';
end
if  strcmp(Geometry,'Plate')
    uString = {'$u_1$' '$u_2$' '$u_3$'};
    sString = {'$\sigma_{11}$' '$\sigma_{22}$' '$\sigma_{33}$' '$\sigma_{23}$' '$\sigma_{13}$' '$\sigma_{12}$'};
    eString = {'$\varepsilon_{11}$' '$\varepsilon_{22}$' '$\varepsilon_{33}$' '$\varepsilon_{23}$' '$\varepsilon_{13}$' '$\varepsilon_{12}$'};
    pString = {'$p_1$' '$p_2$' '$p_3$'};
    uString2 = {'$\varphi(u_1)$' '$\varphi(u_2)$' '$\varphi(u_3)$'};
    sString2 = {'$\varphi(\sigma_{11})$' '$\varphi(\sigma_{22})$' '$\varphi(\sigma_{33})$' '$\varphi(\sigma_{23})$' '$\varphi(\sigma_{13})$' '$\varphi(\sigma_{12})$'};
    eString2 = {'$\varphi(\varepsilon_{11})$' '$\varphi(\varepsilon_{22})$' '$\varphi(\varepsilon_{33}$)' '$\varphi(\varepsilon_{23})$' '$\varphi(\varepsilon_{13})$' '$\varphi(\varepsilon_{12})$'};
    uCheckBoxString = {'u1' 'u2' 'u3'};
    sCheckBoxString = {[char(963) '11'] [char(963) '22'] [char(963) '33'] [char(963) '23'] [char(963) '13'] [char(963) '12']};
    eCheckBoxString = {[char(949) '11'] [char(949) '22'] [char(949) '33'] [char(949) '23'] [char(949) '13'] [char(949) '12']};
    pCheckBoxString = {'p1' 'p2' 'p3'};
else
    uString = {'$u_z$' '$u_\theta$' '$u_r$'};
    sString = {'$\sigma_{zz}$' '$\sigma_{\theta\theta}$' '$\sigma_{rr}$' '$\sigma_{\theta r}$' '$\sigma_{zr}$' '$\sigma_{z\theta}$'};
    eString = {'$\varepsilon_{zz}$' '$\varepsilon_{\theta\theta}$' '$\varepsilon_{rr}$' '$\varepsilon_{\theta r}$' '$\varepsilon_{zr}$' '$\varepsilon_{z\theta}$'};
    pString = {'$p_z$' '$p_\theta$' '$p_r$'};
    uString2 = {'$\varphi(u_z)$' '$\varphi(u_\theta)$' '$\varphi(u_r)$'};
    sString2 = {'$\varphi(\sigma_{zz})$' '$\varphi(\sigma_{\theta\theta})$' '$\varphi(\sigma_{rr})$' '$\varphi(\sigma_{\theta r})$' '$\varphi(\sigma_{zr})$' '$\varphi(\sigma_{z\theta})$'};
    eString2 = {'$\varphi(\varepsilon_{zz})$' '$\varphi(\varepsilon_{\theta\theta})$' '$\varphi(\varepsilon_{rr})$' '$\varphi(\varepsilon_{\theta r})$' '$\varphi(\varepsilon_{zr})$' '$\varphi(\varepsilon_{z\theta})$'};
    uCheckBoxString = {'uz' ['u',char(952)] 'ur'};
    sCheckBoxString = {[char(963),'zz'] [char(963) char(952) char(952)] [char(963) 'rr'] [char(963) char(952) 'r'] [char(963) 'zr'] [char(963) 'z' char(952)]};
    eCheckBoxString = {[char(949),'zz'] [char(949) char(952) char(952)] [char(949) 'rr'] [char(949) char(952) 'r'] [char(949) 'zr'] [char(949) 'z' char(952)]};
    pCheckBoxString = {'pz' ['p',char(952)] 'pr'};
end
EString = {'$E_\mathrm{srn}$','$E_\mathrm{kin}$','$E_\mathrm{tot}$'};

f = figure('Icon',which('DC_Logo16.png'),'Name','Live-through-thickness profiles','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Options','MenuSelectedFcn',@Options_Callback)

t = tiledlayout(3,3);
t.TileSpacing = 'compact'; % loose compact tight none
% t.Padding = 'compact'; % loose compact tight
t.Title.FontSize = FontSizeHeadLine;
t.Title.Interpreter = 'latex';
t1 = tiledlayout(t,2,1);
t1.Layout.Tile = 5;
t1.Layout.TileSpan = [2 2];
t1.TileSpacing = 'none';
if  isempty(LayupString)
    if  HeadLine
        if  strcmp(Geometry,'Plate')
            String = ['Through-thickness profiles for ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' plate'];
        elseif strcmp(Geometry,'Rod')
            String = ['Through-thickness profiles for ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' rod'];
        elseif strcmp(Geometry,'Pipe')
            String = ['Through-thickness profiles for ',num2str(Thickness*1e3),'\,$\times$\,',num2str((Thickness-ThicknessInner)*1e3/2),'\,mm ',replace(Material{1}.Name,'_','\_'),' pipe'];
        elseif strcmp(Geometry,'Circumferential')
            String = ['Through-thickness profiles for ',num2str(Thickness*1e3),'\,$\times$\,',num2str((Thickness-ThicknessInner)*1e3/2),'\,mm ',replace(Material{1}.Name,'_','\_'),' circumference'];
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
        t.Title.String = String;
    end
else
    if  Hybrid
        Material{1}.Name = 'hybrid';
    end
    if  HeadLine > 0
        if  HeadLine == 1
            String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
        elseif ~SymmetricSystem && HeadLine == 2
            if  Repetitions == 1
                String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
            else
                String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
            end
        elseif SymmetricSystem && HeadLine == 2
            if  Repetitions == 1
                String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
            else
                String = ['Through-thickness profiles for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
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
end
if  ~Decoupled
    ce = '$c_{\mathrm e1}$';
else
    ce = '$c_{\mathrm e}$';
end

ax1 = nexttile;
x = xline(0,'Color',[.6 .6 .6],'PickableParts','none');
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
ax1.YLabel.String = [YLabelString,' (mm)'];
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';
ax1.TickLabelInterpreter = 'latex';
legend(uString,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax2 = nexttile;
x = xline(0,'Color',[.6 .6 .6],'PickableParts','none');
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
ax2.YLabel.String = [YLabelString,' (mm)'];
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.TickLabelInterpreter = 'latex';
legend(EString,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax3 = nexttile;
x = xline(0,'Color',[.6 .6 .6],'PickableParts','none');
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
ax3.YLabel.String = [YLabelString,' (mm)'];
ax3.XLabel.Interpreter = 'latex';
ax3.YLabel.Interpreter = 'latex';
ax3.TickLabelInterpreter = 'latex';
legend(pString,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax4 = nexttile;
x = xline(0,'Color',[.6 .6 .6],'PickableParts','none');
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
ax4.YLabel.String = [YLabelString,' (mm)'];
ax4.XLabel.Interpreter = 'latex';
ax4.YLabel.Interpreter = 'latex';
ax4.TickLabelInterpreter = 'latex';
legend(sString,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax7 = nexttile(7);
x = xline(0,'Color',[.6 .6 .6],'PickableParts','none');
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
ax7.YLabel.String = [YLabelString,' (mm)'];
ax7.XLabel.Interpreter = 'latex';
ax7.YLabel.Interpreter = 'latex';
ax7.TickLabelInterpreter = 'latex';
legend(eString,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex');
tb = axtoolbar('default');
tb.Visible = 'on';

ax51 = nexttile(t1);
Plotter(ax51,1,Quantity(1))
ax52 = nexttile(t1);
Plotter(ax52,2,7) % make test plot to obtain YAxis_alpha 
YAxis_alpha = ['[',num2str(ax52.YLim),']'];
delete(ax52.Children)
Plotter(ax52,2,Quantity(2))
Marker1 = scatter(ax51,-1e3,0,'k','PickableParts','none');
Marker2 = scatter(ax52,-1e3,0,'k','PickableParts','none');
YAxis1 = ['[',num2str(ax51.YLim),']'];
YAxis2 = ['[',num2str(ax52.YLim),']'];
YAxis_cp = YAxis1;
YAxis_ce = YAxis2;
Options_Callback
f.WindowState = 'maximized';

function Plotter(ax,m,Q)
    hold on
    Line = line().empty;
    if  strcmp(Geometry,'Plate')
        if  LambModes && SymmetricModes && ~isempty(SLamb{1})
            for i = 1:size(SLamb,2)
                Line(end+1) = plot(ax,SLamb{i}(:,1),SLamb{i}(:,Q),'LineWidth',LineWidth,'Color',SColor);
                Line(end).UserData.Mode = ['S',num2str(i-1)];
                Line(end).UserData.ModeName = ['S$_{',num2str(i-1),'}$'];
                Line(end).UserData.FullData = SLamb{i}(:,[4 5 7]);
            end
        end
        if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
            for i = 1:size(ALamb,2)
                Line(end+1) = plot(ax,ALamb{i}(:,1),ALamb{i}(:,Q),'LineWidth',LineWidth,'Color',AColor);
                Line(end).UserData.Mode = ['A',num2str(i-1)];
                Line(end).UserData.ModeName = ['A$_{',num2str(i-1),'}$'];
                Line(end).UserData.FullData = ALamb{i}(:,[4 5 7]);
            end
        end
        if  LambModes && ~isempty(BLamb{1})
            for i = 1:size(BLamb,2)
                Line(end+1) = plot(ax,BLamb{i}(:,1),BLamb{i}(:,Q),'LineWidth',LineWidth,'Color',BColor);
                Line(end).UserData.Mode = ['B',num2str(i-1)];
                Line(end).UserData.ModeName = ['B$_{',num2str(i-1),'}$'];
                Line(end).UserData.FullData = BLamb{i}(:,[4 5 7]);
            end
        end    
        if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
            for i = 1:size(SShear,2)
                Line(end+1) = plot(ax,SShear{i}(:,1),SShear{i}(:,Q),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                Line(end).UserData.Mode = ['SSH',num2str(i-1)];
                Line(end).UserData.ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                Line(end).UserData.FullData = SShear{i}(:,[4 5 7]);
            end
        end
        if  ShearHorizontalModes && AntisymmetricModes && ~isempty(AShear{1})
            for i = 1:size(AShear,2)
                Line(end+1) = plot(ax,AShear{i}(:,1),AShear{i}(:,Q),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
                Line(end).UserData.Mode = ['ASH',num2str(i)];
                Line(end).UserData.ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                Line(end).UserData.FullData = AShear{i}(:,[4 5 7]);
            end
        end
        if  ShearHorizontalModes && ~isempty(BShear{1})
            for i = 1:size(BShear,2)
                Line(end+1) = plot(ax,BShear{i}(:,1),BShear{i}(:,Q),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                Line(end).UserData.Mode = ['BSH',num2str(i-1)];
                Line(end).UserData.ModeName = ['B$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                Line(end).UserData.FullData = BShear{i}(:,[4 5 7]);
            end
        end
    elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
        if  LongitudinalModes && ~isempty(L{1})
            for i = 1:length(L)
                Line(end+1) = plot(ax,L{i}(:,1),L{i}(:,Q),'LineWidth',LineWidth,'Color',SColor);
                Line(end).UserData.Mode = ['L(0,',num2str(i),')'];
                Line(end).UserData.ModeName = ['L(0,',num2str(i),')'];
                Line(end).UserData.FullData = L{i}(:,[4 5 7]);
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
                    Line(end+1) = plot(ax,F{n}{i}(:,1),F{n}{i}(:,Q),'LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));
                    Line(end).UserData.Mode = ['F(',num2str(n),',',num2str(i),')'];
                    Line(end).UserData.ModeName = ['F(',num2str(n),',',num2str(i),')'];
                    Line(end).UserData.FullData = F{n}{i}(:,[4 5 7]);
                end
                if  LineColorsIndex == height(LineColors)
                    LineColorsIndex = 0;
                end
            end
        end
        if  TorsionalModes && ~isempty(T{1})
            for i = 1:length(T)
                Line(end+1) = plot(ax,T{i}(:,1),T{i}(:,Q),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                Line(end).UserData.Mode = ['T(0,',num2str(i),')'];
                Line(end).UserData.ModeName = ['T(0,',num2str(i),')'];
                Line(end).UserData.FullData = T{i}(:,[4 5 7]);
            end
        end
    elseif strcmp(Geometry,'Circumferential')
        if  LambModes && ~isempty(CLamb{1})
            for i = 1:size(CLamb,2)
                Line(end+1) = plot(ax,CLamb{i}(:,1),CLamb{i}(:,Q),'LineWidth',LineWidth,'Color',BColor);
                Line(end).UserData.Mode = ['C',num2str(i-1)];
                Line(end).UserData.ModeName = ['C$_{',num2str(i-1),'}$'];
                Line(end).UserData.FullData = CLamb{i}(:,[4 5 7]);
            end
        end
        if  ShearHorizontalModes && ~isempty(CShear{1})
            for i = 1:size(CShear,2)
                Line(end+1) = plot(ax,CShear{i}(:,1),CShear{i}(:,Q),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                Line(end).UserData.Mode = ['CSH',num2str(i-1)];
                Line(end).UserData.ModeName = ['C$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                Line(end).UserData.FullData = CShear{i}(:,[4 5 7]);
            end
        end
    end
    delete(ax.Children(end))
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    ax.XLim = [0 FrequencyLimit];
    tb = axtoolbar('default');
    tb.Visible = 'on';
    datacursormode on
    d = datacursormode(f);
    d.Interpreter = 'latex';
    d.UpdateFcn = @Cursor;
    if  m == 1
        if  Q == 4
            ax.YLabel.String = 'Phase velocity (m/ms)';
        elseif Q == 5
            ax.YLabel.String = 'Energy velocity (m/ms)';
        end
        ax.XTickLabel = [];
    elseif m == 2
        if  Q == 5
            ax.YLabel.String = 'Energy velocity (m/ms)';
        elseif Q == 7
            ax.YLabel.String = 'Attenuation (Np/m)';
        end
        ax.XLabel.String = 'Frequency (kHz)';
        ax.XLabel.FontSize = FontSizeAxesLabels;
        ax.XLabel.Interpreter = 'latex';
        ax.YTickMode = 'auto';
        ax.YTick(end) = [];
    end
end
function output_txt = Cursor(~,event_obj)
    if  contains(event_obj.Target.Parent.XLabel.String,'Displacement')
        if  Phase
            if  event_obj.Target.Color == Colors(1,:)
                output_txt = {[uString2{1},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(2,:)
                output_txt = {[uString2{2},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(3,:)
                output_txt = {[uString2{3},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            end
        else
            if  event_obj.Target.Color == Colors(1,:)
                output_txt = {[uString{1},': \textbf{',num2str(event_obj.Position(1),6),'}\,nm'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(2,:)
                output_txt = {[uString{2},': \textbf{',num2str(event_obj.Position(1),6),'}\,nm'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(3,:)
                output_txt = {[uString{3},': \textbf{',num2str(event_obj.Position(1),6),'}\,nm'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            end
        end
    elseif contains(event_obj.Target.Parent.XLabel.String,'Stress')
        if  Phase
            if  event_obj.Target.Color == Colors(4,:)
                output_txt = {[sString2{1},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(5,:)
                output_txt = {[sString2{2},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(6,:)
                output_txt = {[sString2{3},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(7,:)
                output_txt = {[sString2{4},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(8,:)
                output_txt = {[sString2{5},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(9,:)
                output_txt = {[sString2{6},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};            
            end
        else
            if  event_obj.Target.Color == Colors(4,:)
                output_txt = {[sString{1},': \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(5,:)
                output_txt = {[sString{2},': \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(6,:)
                output_txt = {[sString{3},': \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(7,:)
                output_txt = {[sString{4},': \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(8,:)
                output_txt = {[sString{5},': \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(9,:)
                output_txt = {[sString{6},': \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};            
            end
        end
    elseif contains(event_obj.Target.Parent.XLabel.String,'Strain')
        if  Phase
            if  event_obj.Target.Color == Colors(10,:)
                output_txt = {[eString2{1},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(11,:)
                output_txt = {[eString2{2},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(12,:)
                output_txt = {[eString2{3},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(13,:)
                output_txt = {[eString2{4},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(14,:)
                output_txt = {[eString2{5},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(15,:)
                output_txt = {[eString2{6},': \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};            
            end
        else
            if  event_obj.Target.Color == Colors(10,:)
                output_txt = {[eString{1},': \textbf{',num2str(event_obj.Position(1),6),'}'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(11,:)
                output_txt = {[eString{2},': \textbf{',num2str(event_obj.Position(1),6),'}'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(12,:)
                output_txt = {[eString{3},': \textbf{',num2str(event_obj.Position(1),6),'}'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(13,:)
                output_txt = {[eString{4},': \textbf{',num2str(event_obj.Position(1),6),'}'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(14,:)
                output_txt = {[eString{5},': \textbf{',num2str(event_obj.Position(1),6),'}'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Colors(15,:)
                output_txt = {[eString{6},': \textbf{',num2str(event_obj.Position(1),6),'}'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};            
            end
        end
    elseif contains(event_obj.Target.Parent.XLabel.String,'Energy density')
        if  event_obj.Target.Color == Colors(16,:)
            output_txt = {[EString{1},': \textbf{',num2str(event_obj.Position(1),6),'}\,J/m$^2$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
        elseif event_obj.Target.Color == Colors(17,:)
            output_txt = {[EString{2},': \textbf{',num2str(event_obj.Position(1),6),'}\,J/m$^2$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
        elseif event_obj.Target.Color == Colors(18,:)
            output_txt = {[EString{3},': \textbf{',num2str(event_obj.Position(1),6),'}\,J/m$^2$'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
        end
    elseif contains(event_obj.Target.Parent.XLabel.String,'Power flow density')
        if  event_obj.Target.Color == Colors(19,:)
            output_txt = {[pString{1},': \textbf{',num2str(event_obj.Position(1),6),'}\,W/m'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
        elseif event_obj.Target.Color == Colors(20,:)
            output_txt = {[pString{2},': \textbf{',num2str(event_obj.Position(1),6),'}\,W/m'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
        elseif event_obj.Target.Color == Colors(21,:)
            output_txt = {[pString{3},': \textbf{',num2str(event_obj.Position(1),6),'}\,W/m'],[YLabelString,': \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
        end
    elseif strcmp(event_obj.Target.Parent.YLabel.String,'Phase velocity (m/ms)')
        Mode = event_obj.Target.UserData.Mode;
        Frequency = event_obj.Position(1);
        SecondaryData = event_obj.Target.UserData.FullData(event_obj.Target.Children(1).DataIndex,[2 3]);
        Marker1.XData = -1e3;
        Marker1.YData = 0;
        Marker2.XData = Frequency;
        if  DispersionDiagrams == 1
            Marker2.YData = SecondaryData(1);
        elseif DispersionDiagrams == 2
            Marker2.YData = SecondaryData(2);
        end
        output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(Frequency,6),'}\,kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'] [ce,': \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'}\,Np/m']};
        UpdateModeShapes
    elseif strcmp(event_obj.Target.Parent.YLabel.String,'Energy velocity (m/ms)')
        Mode = event_obj.Target.UserData.Mode;
        Frequency = event_obj.Position(1);
        SecondaryData = event_obj.Target.UserData.FullData(event_obj.Target.Children(1).DataIndex,[1 3]);
        if  DispersionDiagrams == 1
            Marker2.XData = -1e3;
            Marker2.YData = 0;
            Marker1.XData = Frequency;
            Marker1.YData = SecondaryData(1);
        elseif DispersionDiagrams == 3
            Marker1.XData = -1e3;
            Marker1.YData = 0;
            Marker2.XData = Frequency;
            Marker2.YData = SecondaryData(2);
        end
        output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(Frequency,6),'}\,kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] [ce,': \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'] ['$\alpha$: \textbf{',num2str(SecondaryData(2),6),'}\,Np/m']};
        UpdateModeShapes
    elseif strcmp(event_obj.Target.Parent.YLabel.String,'Attenuation (Np/m)')
        Mode = event_obj.Target.UserData.Mode;
        Frequency = event_obj.Position(1);
        SecondaryData = event_obj.Target.UserData.FullData(event_obj.Target.Children(1).DataIndex,[1 2]);
        Marker2.XData = -1e3;
        Marker2.YData = 0;
        Marker1.XData = Frequency;
        if  DispersionDiagrams == 2
            Marker1.YData = SecondaryData(1);
        elseif DispersionDiagrams == 3
            Marker1.YData = SecondaryData(2);
        end
        output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(Frequency,6),'}\,kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(SecondaryData(1),6),'}\,m/ms'] [ce,': \textbf{',num2str(SecondaryData(2),6),'}\,m/ms'] ['$\alpha$: \textbf{',num2str(event_obj.Position(2),6),'}\,Np/m']};
        UpdateModeShapes
    end
end
function Options_Callback(~,~)
    f_Options = figure('Icon',which('DC_Logo16.png'),'WindowStyle','alwaysontop','NumberTitle','off','Name','Options','Visible','off','MenuBar','none','Position',[0 0 210 325],'CloseRequestFcn',@CloseRequest);
    f_Options.Units = 'normalized';
    movegui(f_Options,'east')
    f_Options.Visible = 'on';
    drawnow
    uicontrol('Parent',f_Options,'Style','text','HorizontalAlignment','left','String','Dispersion diagrams','Position',[10 295 103 13]);
    YAxisText1UI = uicontrol('Parent',f_Options,'Style','text','HorizontalAlignment','left','String',YAxisText1UIString,'Position',[10 265 70 13]);
    YAxisText2UI = uicontrol('Parent',f_Options,'Style','text','HorizontalAlignment','left','String',YAxisText2UIString,'Position',[10 235 70 13]);
    if  strcmp(Geometry,'Plate')
        uicontrol('Parent',f_Options,'Style','text','HorizontalAlignment','left','String','Samples x3','Position',[10 205 61 13]);
    else
        uicontrol('Parent',f_Options,'Style','text','HorizontalAlignment','left','String','Samples r','Position',[10 205 61 13]);
    end
    uicontrol('Parent',f_Options,'Style','text','HorizontalAlignment','left','String','Half-spaces','Position',[10 175 64 13]);
    uicontrol('Parent',f_Options,'Style','text','HorizontalAlignment','left','String','Phase','Position',[10 145 35 13]);
    uicontrol('Parent',f_Options,'Style','popupmenu','String',{'cp-ce',['cp-',char(945)],['ce-',char(945)]},'Value',DispersionDiagrams,'Tooltip','Select dispersion diagrams.','Position',[115 290 65 23],'Callback',@Callback,'Tag','26');
    uicontrol('Parent',f_Options,'Style','checkbox','Value',YAxesEnable,'Tooltip','Activate to manually set the Y-axes limits.','Position',[85 245 20 23],'Callback',@Callback,'Tag','27');
    YAxis1UI = uicontrol('Parent',f_Options,'Style','edit','String',YAxis1,'Tooltip','Set the Y-axis limits of the upper dispersion diagram.','Position',[115 260 75 23],'Callback',@Callback,'Tag','28');
    YAxis2UI = uicontrol('Parent',f_Options,'Style','edit','String',YAxis2,'Tooltip','Set the Y-axis limits of the lower dispersion diagram.','Position',[115 230 75 23],'Callback',@Callback,'Tag','29');
    uicontrol('Parent',f_Options,'Style','edit','String',SamplesX3,'Tooltip','Enter the number of sample points over the plate''s thickness (x3) at which the selected quantities are calculated.','Position',[115 200 50 23],'Callback',@Callback,'Tag','1');
    ShowHalfSpacesUI = uicontrol('Parent',f_Options,'Style','checkbox','Value',ShowHalfSpaces,'Tooltip','Check this to show the quantities in the upper and lower fluid.','Position',[85 170 20 23],'Callback',@Callback,'Tag','2');        
    HalfSpacesUI = uicontrol('Parent',f_Options,'Style','edit','String',HalfSpaces,'Tooltip','Set the height of the half-spaces in plate thicknesses.','Position',[115 170 50 23],'Callback',@Callback,'Tag','3');
    uicontrol('Parent',f_Options,'Style','checkbox','Value',Phase,'Tooltip','Check this to plot the phase of the field components.','Position',[85 140 20 23],'Callback',@Callback,'Tag','4');
    uicontrol('Parent',f_Options,'Style','checkbox','String',uCheckBoxString{1},'Value',Plot(1),'Position',[10 110 35 23],'Callback',@Callback,'Tag','5');
    uicontrol('Parent',f_Options,'Style','checkbox','String',uCheckBoxString{2},'Value',Plot(2),'Position',[10 90 35 23],'Callback',@Callback,'Tag','6');
    uicontrol('Parent',f_Options,'Style','checkbox','String',uCheckBoxString{3},'Value',Plot(3),'Position',[10 70 35 23],'Callback',@Callback,'Tag','7');
    uicontrol('Parent',f_Options,'Style','checkbox','String',sCheckBoxString{1},'Value',Plot(4),'Position',[55 110 50 23],'Callback',@Callback,'Tag','8');
    uicontrol('Parent',f_Options,'Style','checkbox','String',sCheckBoxString{2},'Value',Plot(5),'Position',[55 90 50 23],'Callback',@Callback,'Tag','9');
    uicontrol('Parent',f_Options,'Style','checkbox','String',sCheckBoxString{3},'Value',Plot(6),'Position',[55 70 50 23],'Callback',@Callback,'Tag','10');
    uicontrol('Parent',f_Options,'Style','checkbox','String',sCheckBoxString{4},'Value',Plot(7),'Position',[55 50 50 23],'Callback',@Callback,'Tag','11');
    uicontrol('Parent',f_Options,'Style','checkbox','String',sCheckBoxString{5},'Value',Plot(8),'Position',[55 30 50 23],'Callback',@Callback,'Tag','12');
    uicontrol('Parent',f_Options,'Style','checkbox','String',sCheckBoxString{6},'Value',Plot(9),'Position',[55 10 50 23],'Callback',@Callback,'Tag','13');
    uicontrol('Parent',f_Options,'Style','checkbox','String',eCheckBoxString{1},'Value',Plot(10),'Position',[105 110 50 23],'Callback',@Callback,'Tag','14');
    uicontrol('Parent',f_Options,'Style','checkbox','String',eCheckBoxString{2},'Value',Plot(11),'Position',[105 90 50 23],'Callback',@Callback,'Tag','15');
    uicontrol('Parent',f_Options,'Style','checkbox','String',eCheckBoxString{3},'Value',Plot(12),'Position',[105 70 50 23],'Callback',@Callback,'Tag','16');
    uicontrol('Parent',f_Options,'Style','checkbox','String',eCheckBoxString{4},'Value',Plot(13),'Position',[105 50 50 23],'Callback',@Callback,'Tag','17');
    uicontrol('Parent',f_Options,'Style','checkbox','String',eCheckBoxString{5},'Value',Plot(14),'Position',[105 30 50 23],'Callback',@Callback,'Tag','18');
    uicontrol('Parent',f_Options,'Style','checkbox','String',eCheckBoxString{6},'Value',Plot(15),'Position',[105 10 50 23],'Callback',@Callback,'Tag','19');
    uicontrol('Parent',f_Options,'Style','checkbox','String','Esrn','Value',Plot(16),'Position',[155 110 60 23],'Callback',@Callback,'Tag','20');
    uicontrol('Parent',f_Options,'Style','checkbox','String','Ekin','Value',Plot(17),'Position',[155 90 60 23],'Callback',@Callback,'Tag','21');
    uicontrol('Parent',f_Options,'Style','checkbox','String','Etot','Value',Plot(18),'Position',[155 70 60 23],'Callback',@Callback,'Tag','22');
    uicontrol('Parent',f_Options,'Style','checkbox','String',pCheckBoxString{1},'Value',Plot(19),'Position',[155 50 35 23],'Callback',@Callback,'Tag','23');
    uicontrol('Parent',f_Options,'Style','checkbox','String',pCheckBoxString{2},'Value',Plot(20),'Position',[155 30 35 23],'Callback',@Callback,'Tag','24');
    uicontrol('Parent',f_Options,'Style','checkbox','String',pCheckBoxString{3},'Value',Plot(21),'Position',[155 10 35 23],'Callback',@Callback,'Tag','25');
    if  FluidLoading
        if  strcmp(Geometry,'Pipe') && ~ToggleUpperFluid
            HalfSpacesUI.Enable = 'off';
            ShowHalfSpacesUI.Enable = 'off';
        else
            if  ShowHalfSpaces
                HalfSpacesUI.Enable = 'on';
            else
                HalfSpacesUI.Enable = 'off';
            end
            ShowHalfSpacesUI.Enable = 'on';
        end
    else
        HalfSpacesUI.Enable = 'off';
        ShowHalfSpacesUI.Enable = 'off';
    end
    if  YAxesEnable
        YAxis1UI.Enable = 'on';
        YAxis2UI.Enable = 'on';
    else
        YAxis1UI.Enable = 'off';
        YAxis2UI.Enable = 'off';
    end
    function Callback(source,~)
        if  strcmp(source.Tag,'26') % Dispersion diagrams
            DispersionDiagrams = source.Value;
            if  source.Value == 1
                Quantity = [4 5]; % cp-ce
            elseif source.Value == 2
                Quantity = [4 7]; % cp-alpha
            elseif source.Value == 3
                Quantity = [5 7]; % ce-alpha
            end
            ax51 = nexttile(t1,1);
            delete(ax51.Children)
            Plotter(ax51,1,Quantity(1))
            ax52 = nexttile(t1,2);
            delete(ax52.Children)
            Plotter(ax52,2,Quantity(2))
            Marker1 = scatter(ax51,-1e3,0,'k','PickableParts','none');
            Marker2 = scatter(ax52,-1e3,0,'k','PickableParts','none');
            if  source.Value == 1
                YAxisText1UIString = 'cp-axis (m/ms)';
                YAxisText2UIString = 'ce-axis (m/ms)';
            elseif source.Value == 2
                YAxisText1UIString = 'cp-axis (m/ms)';
                YAxisText2UIString = [char(945),'-axis (Np/m)'];
            elseif source.Value == 3
                YAxisText1UIString = 'ce-axis (m/ms)';
                YAxisText2UIString = [char(945),'-axis (Np/m)'];
            end
            YAxisText1UI.String = YAxisText1UIString;
            YAxisText2UI.String = YAxisText2UIString;
            if  YAxesEnable
                if  source.Value == 1
                    YAxis1 = YAxis_cp;
                    YAxis2 = YAxis_ce;
                elseif source.Value == 2
                    YAxis1 = YAxis_cp;
                    YAxis2 = YAxis_alpha;
                elseif source.Value == 3
                    YAxis1 = YAxis_ce;
                    YAxis2 = YAxis_alpha;
                end
                ax51.YLim = eval(YAxis1);
                ax52.YLim = eval(YAxis2);
                ax51.YTickMode = 'auto';
                ax52.YTickMode = 'auto';
                ax52.YTick(end) = [];
            else
                YAxis1 = ['[',num2str(ax51.YLim),']'];
                YAxis2 = ['[',num2str(ax52.YLim),']'];
            end
            YAxis1UI.String = YAxis1;
            YAxis2UI.String = YAxis2;
        elseif strcmp(source.Tag,'27') % Y-axes (on/off) 
            YAxesEnable = source.Value;
            if  YAxesEnable
                YAxis1UI.Enable = 'on';
                YAxis2UI.Enable = 'on';
                if  DispersionDiagrams == 1
                    YAxis1 = YAxis_cp;
                    YAxis2 = YAxis_ce;
                elseif DispersionDiagrams == 2
                    YAxis1 = YAxis_cp;
                    YAxis2 = YAxis_alpha;
                elseif DispersionDiagrams == 3
                    YAxis1 = YAxis_ce;
                    YAxis2 = YAxis_alpha;
                end
                ax51.YLim = eval(YAxis1);
                ax52.YLim = eval(YAxis2);
                ax51.YTickMode = 'auto';
                ax52.YTickMode = 'auto';
                ax52.YTick(end) = [];
            else
                YAxis1UI.Enable = 'off';
                YAxis2UI.Enable = 'off';
                ax51.YLimMode = 'auto';
                ax52.YLimMode = 'auto';
                ax52.YTick(end) = [];
                YAxis1 = ['[',num2str(ax51.YLim),']'];
                YAxis2 = ['[',num2str(ax52.YLim),']'];
            end
            YAxis1UI.String = YAxis1;
            YAxis2UI.String = YAxis2;
        elseif strcmp(source.Tag,'28') % Y-axis1
            if  DispersionDiagrams < 3
                YAxis_cp = source.String;
            else
                YAxis_ce = source.String;
            end
            YAxis1 = source.String;
            ax51.YLim = eval(YAxis1);
            ax51.YTickMode = 'auto';
        elseif strcmp(source.Tag,'29') % Y-axis2
            if  DispersionDiagrams == 1
                YAxis_ce = source.String;
            else
                YAxis_alpha = source.String;
            end
            YAxis2 = source.String;
            ax52.YLim = eval(YAxis2);
            ax52.YTickMode = 'auto';
            ax52.YTick(end) = [];
        elseif strcmp(source.Tag,'1') % Samples x3
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
        if  str2double(source.Tag) < 26
            UpdateModeShapes
        end
    end
    function CloseRequest(~,~)
        delete(f_Options)
    end
end
function UpdateModeShapes
    if  strcmp(Geometry,'Plate')
        [u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,x3Total,~,~] = ModeShapeLinesComputer(FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,ALamb,AShear,BLamb,BShear,SLamb,SShear,c,Delta,Material,Layers,Frequency,Mode,Thickness,SamplesX3,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,ShowHalfSpaces,HalfSpaces,Phase);
    elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe') % Material{1}!
        [u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,r,~,~,~] = ModeShapeLinesComputer_Isotropic_Rod_Pipe(Geometry,Material{1},FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Sink,F,L,T,Frequency,Mode,Ro,Ri,SamplesX3,ShowHalfSpaces,HalfSpaces,Phase);
    elseif strcmp(Geometry,'Circumferential') % Material{1}!
        [u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,r,~,~] = ModeShapeLinesComputer_Isotropic_Circumferential(Material{1},CLamb,CShear,Frequency,Mode,Ro,Ri,SamplesX3,Phase);
    end
    if  Plot(1)
        if  strcmp(Geometry,'Plate')
            u1.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            u1.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            u1.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            u2.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            u2.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            u2.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            u3.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            u3.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            u3.YData = (r-Ri)*1e3;
        end
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
    if  strcmp(Geometry,'Plate')
        ax1.YLim = -1e3*[x3Total(end) x3Total(1)];
    elseif strcmp(Geometry,'Rod') || (strcmp(Geometry,'Pipe') && ToggleLowerFluid)
        ax1.YLim = [0 r(end)*1e3];
    elseif strcmp(Geometry,'Pipe') && ~ToggleLowerFluid
        ax1.YLim = [r(1) r(end)]*1e3;
    elseif strcmp(Geometry,'Circumferential')
        ax1.YLim = [0 (r(end)-Ri)*1e3];
    end
    if  Plot(16)
        if  strcmp(Geometry,'Plate')
            Esrn.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            Esrn.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            Esrn.YData = (r-Ri)*1e3;
        end
        if  Phase
            Esrn.XData = zeros(size(u,1),1);
        else
            Esrn.XData = StrainEnergyDensity;
        end
    else
        Esrn.XData = 0;
        Esrn.YData = 0;
    end
    if  Plot(17)
        if  strcmp(Geometry,'Plate')
            Ekin.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            Ekin.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            Ekin.YData = (r-Ri)*1e3;
        end
        if  Phase
            Ekin.XData = zeros(size(u,1),1);
        else
            Ekin.XData = KineticEnergyDensity;
        end
    else
        Ekin.XData = 0;
        Ekin.YData = 0;
    end
    if  Plot(18)
        if  strcmp(Geometry,'Plate')
            Etot.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            Etot.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            Etot.YData = (r-Ri)*1e3;
        end
        if  Phase
            Etot.XData = zeros(size(u,1),1);
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
    if  strcmp(Geometry,'Plate')
        ax2.YLim = -1e3*[x3Total(end) x3Total(1)];
    elseif strcmp(Geometry,'Rod') || (strcmp(Geometry,'Pipe') && ToggleLowerFluid)
        ax2.YLim = [0 r(end)*1e3];
    elseif strcmp(Geometry,'Pipe') && ~ToggleLowerFluid
        ax2.YLim = [r(1) r(end)]*1e3;
    elseif strcmp(Geometry,'Circumferential')
        ax2.YLim = [0 (r(end)-Ri)*1e3];
    end
    if  Plot(19)
        if  strcmp(Geometry,'Plate')
            p1.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            p1.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            p1.YData = (r-Ri)*1e3;
        end
        if  Phase
            p1.XData = zeros(size(u,1),1);
        else
            p1.XData = PowerFlowDensity(:,1);
        end
    else
        p1.XData = 0;
        p1.YData = 0;
    end
    if  Plot(20)
        if  strcmp(Geometry,'Plate')
            p2.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            p2.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            p2.YData = (r-Ri)*1e3;
        end
        if  Phase
            p2.XData = zeros(size(u,1),1);
        else
            p2.XData = PowerFlowDensity(:,2);
        end
    else
        p2.XData = 0;
        p2.YData = 0;
    end    
    if  Plot(21)
        if  strcmp(Geometry,'Plate')
            p3.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            p3.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            p3.YData = (r-Ri)*1e3;
        end
        if  Phase
            p3.XData = zeros(size(u,1),1);
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
    if  strcmp(Geometry,'Plate')
        ax3.YLim = -1e3*[x3Total(end) x3Total(1)];
    elseif strcmp(Geometry,'Rod') || (strcmp(Geometry,'Pipe') && ToggleLowerFluid)
        ax3.YLim = [0 r(end)*1e3];
    elseif strcmp(Geometry,'Pipe') && ~ToggleLowerFluid
        ax3.YLim = [r(1) r(end)]*1e3;
    elseif strcmp(Geometry,'Circumferential')
        ax3.YLim = [0 (r(end)-Ri)*1e3];
    end
    if  Plot(4)
        if  strcmp(Geometry,'Plate')
            sigma11.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            sigma11.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            sigma11.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            sigma22.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            sigma22.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            sigma22.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            sigma33.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            sigma33.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            sigma33.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            sigma23.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            sigma23.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            sigma23.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            sigma13.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            sigma13.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            sigma13.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            sigma12.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            sigma12.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            sigma12.YData = (r-Ri)*1e3;
        end
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
    if  strcmp(Geometry,'Plate')
        ax4.YLim = -1e3*[x3Total(end) x3Total(1)];
    elseif strcmp(Geometry,'Rod') || (strcmp(Geometry,'Pipe') && ToggleLowerFluid)
        ax4.YLim = [0 r(end)*1e3];
    elseif strcmp(Geometry,'Pipe') && ~ToggleLowerFluid
        ax4.YLim = [r(1) r(end)]*1e3;
    elseif strcmp(Geometry,'Circumferential')
        ax4.YLim = [0 (r(end)-Ri)*1e3];
    end
    if  Plot(10)
        if  strcmp(Geometry,'Plate')
            epsilon11.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            epsilon11.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            epsilon11.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            epsilon22.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            epsilon22.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            epsilon22.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            epsilon33.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            epsilon33.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            epsilon33.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            epsilon23.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            epsilon23.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            epsilon23.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            epsilon13.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            epsilon13.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            epsilon13.YData = (r-Ri)*1e3;
        end
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
        if  strcmp(Geometry,'Plate')
            epsilon12.YData = -x3Total*1e3;
        elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            epsilon12.YData = r*1e3;
        elseif strcmp(Geometry,'Circumferential')
            epsilon12.YData = (r-Ri)*1e3;
        end
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
    if  strcmp(Geometry,'Plate')
        ax7.YLim = -1e3*[x3Total(end) x3Total(1)];
    elseif strcmp(Geometry,'Rod') || (strcmp(Geometry,'Pipe') && ToggleLowerFluid)
        ax7.YLim = [0 r(end)*1e3];
    elseif strcmp(Geometry,'Pipe') && ~ToggleLowerFluid
        ax7.YLim = [r(1) r(end)]*1e3;
    elseif strcmp(Geometry,'Circumferential')
        ax7.YLim = [0 (r(end)-Ri)*1e3];
    end
end
end