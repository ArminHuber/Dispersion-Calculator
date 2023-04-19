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
function ModeShapeAnimation_Anisotropic(Hybrid,Layers,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,LayupString,Directory,FontSizeHeadLine,FontSizeAxesLabels,FrameRate,Frequency,GridLine,HeadLine,LineWidth,Material,Mode,Movie,FileName,MovieQuality,PlateThickness,PropagationAngle,SamplesX1,SamplesX3,Scale,Repetitions,SuperLayerSize,SymmetricSystem,Time,u,Undistorted,x1,x3,x3Total,p,ShowHalfSpaces,HalfSpaces)
%#ok<*AGROW>
[X1,X3] = meshgrid(x1,x3Total);
Ratio = (abs(x3Total(1))+abs(x3Total(end)))/x1(end);
for i = 1:size(u,1)
    u1Max(i) = abs(real(u{i,1}(1,1)));
    u3Max(i) = abs(real(u{i,1}(1,2))); % using real instead of imag is correct here; it is a shift in phase
end
u1Max = max(u1Max);
u3Max = max(u3Max);
if  u1Max > u3Max
    for l = 1:length(Time)
        for i = 1:size(X1,2)
            for j = 1:size(X1,1)
                if  FluidLoading && ShowHalfSpaces
                    X1Distorted{l}(j,i) = X1(j,i)+Scale*x1(2)*SamplesX1/(80*(1+2*HalfSpaces))/u1Max*real(u{j,i}(l,1));
                    X3Distorted{l}(j,i) = X3(j,i)+Scale*x1(2)*SamplesX1/(80*(1+2*HalfSpaces))/u1Max*real(u{j,i}(l,2))*Ratio;
                else
                    X1Distorted{l}(j,i) = X1(j,i)+Scale*x1(2)*SamplesX1/80/u1Max*real(u{j,i}(l,1));
                    X3Distorted{l}(j,i) = X3(j,i)+Scale*x1(2)*SamplesX1/80/u1Max*real(u{j,i}(l,2))*Ratio;
                end
            end
        end
    end
else
    for l = 1:length(Time)
        for i = 1:size(X1,2)
            for j = 1:size(X1,1)
                if  FluidLoading && ShowHalfSpaces
                    X1Distorted{l}(j,i) = X1(j,i)+Scale*abs(x3(2)-x3(1))*SamplesX3/(40*(1+2*HalfSpaces))/u3Max*real(u{j,i}(l,1));
                    X3Distorted{l}(j,i) = X3(j,i)+Scale*abs(x3(2)-x3(1))*SamplesX3/(40*(1+2*HalfSpaces))/u3Max*real(u{j,i}(l,2))*Ratio;
                else
                    X1Distorted{l}(j,i) = X1(j,i)+Scale*abs(x3(2)-x3(1))*SamplesX3/40/u3Max*real(u{j,i}(l,1));
                    X3Distorted{l}(j,i) = X3(j,i)+Scale*abs(x3(2)-x3(1))*SamplesX3/40/u3Max*real(u{j,i}(l,2))*Ratio;
                end
            end
        end
    end
end
f = figure('Name','Mode shape','MenuBar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));    
ax = gca;
axis off
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
if  Hybrid
    Material{1}.Name = 'hybrid';
end
if  HeadLine > 0
    if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
        ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
    elseif contains(Mode,'SH')
        ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
    elseif contains(Mode,'Scholte')
        ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
    end
    if  HeadLine == 1
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_'))];
    elseif ~SymmetricSystem && HeadLine == 2
        if  Repetitions == 1
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']'];
        elseif Repetitions > 1
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        end
    elseif SymmetricSystem && HeadLine == 2
        if  Repetitions == 1
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{\mathrm s}$'];
        elseif Repetitions > 1
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
        end
    end
    if  FluidLoading
        if  ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/',char(join(split(LowerFluid.Name,'_'),'\_')));
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/vacuum');
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in vacuum/',char(join(split(LowerFluid.Name,'_'),'\_')));
        end
    end
    ax.Title.String = String;
end
if  ~contains(Mode,'SH')
    text(.5-.155*FontSizeAxesLabels/30,.05,'Propagation direction ($x_1)$','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels)
    text(-.02,.5-.14*FontSizeAxesLabels/30,'Thickness ($x_3$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels,'rotation',90);
else
    ax.Position = [0.13,0.21,0.775,0.615]; % default [.13 .11 .775 .815]
    text(.5-.125*FontSizeAxesLabels/30,-.073,'Shear horizontal ($x_2$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels)
    text(-.02,.5-.18*FontSizeAxesLabels/30,'Thickness ($x_3$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels,'rotation',90);
end
ax.XLim = [-.075*x1(end) x1(end)+.075*x1(end)];
ax.YDir = 'reverse';
hold on
if  FluidLoading && ShowHalfSpaces
	n = HalfSpaces*SamplesX3*Layers;
    if  strcmp(UpperFluid.Name,LowerFluid.Name)
        UpperFluidColor = [0 .5 1];
        LowerFluidColor = [0 .5 1];
    else
        if  UpperFluid.Density*UpperFluid.Velocity > LowerFluid.Density*LowerFluid.Velocity
            UpperFluidColor = [0 .5 1];
            LowerFluidColor = [0 .7 1];
        else
            UpperFluidColor = [0 .7 1];
            LowerFluidColor = [0 .5 1];
        end
    end
end
if  Undistorted
    if  ~FluidLoading || ~ShowHalfSpaces
        line(ax,X1(:,1:GridLine:end),X3(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0]);
        line(ax,X1(1:GridLine:end,:)',X3(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
    elseif FluidLoading && ShowHalfSpaces
        line(ax,X1(1:n+1,1:GridLine:end),X3(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(1:GridLine:n+1,:)',X3(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(end-n:end,1:GridLine:end),X3(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(end-n:GridLine:end,:)',X3(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(n+2:end-n-1,1:GridLine:end),X3(n+2:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(n+2:GridLine:end-n-1,:)',X3(n+2:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
    end
end
i = 0;
if  ~Movie
    while i <= length(Time) && ishghandle(f) == 1
        i = i+1;
        if  ~FluidLoading || ~ShowHalfSpaces
            h1 = line(ax,X1Distorted{i}(:,1:GridLine:end),X3Distorted{i}(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(1:GridLine:end,:)',X3Distorted{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
        elseif FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid && ToggleLowerFluid
                h3 = line(ax,X1Distorted{i}(1:n+1,1:GridLine:end),X3Distorted{i}(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',UpperFluidColor);
                h4 = line(ax,X1Distorted{i}(1:GridLine:n+1,:)',X3Distorted{i}(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',UpperFluidColor);
                h5 = line(ax,X1Distorted{i}(end-n:end,1:GridLine:end),X3Distorted{i}(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',LowerFluidColor);
                h6 = line(ax,X1Distorted{i}(end-n:GridLine:end,:)',X3Distorted{i}(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',LowerFluidColor);
                h1 = line(ax,X1Distorted{i}(n+2:end-n-1,1:GridLine:end),X3Distorted{i}(n+2:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
                h2 = line(ax,X1Distorted{i}(n+2:GridLine:end-n-1,:)',X3Distorted{i}(n+2:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k');
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                h3 = line(ax,X1Distorted{i}(1:n+1,1:GridLine:end),X3Distorted{i}(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
                h4 = line(ax,X1Distorted{i}(1:GridLine:n+1,:)',X3Distorted{i}(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
                h1 = line(ax,X1Distorted{i}(n+2:end,1:GridLine:end),X3Distorted{i}(n+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
                h2 = line(ax,X1Distorted{i}(n+2:GridLine:end,:)',X3Distorted{i}(n+2:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                h5 = line(ax,X1Distorted{i}(end-n:end,1:GridLine:end),X3Distorted{i}(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
                h6 = line(ax,X1Distorted{i}(end-n:GridLine:end,:)',X3Distorted{i}(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
                h1 = line(ax,X1Distorted{i}(1:end-n-1,1:GridLine:end),X3Distorted{i}(1:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
                h2 = line(ax,X1Distorted{i}(1:GridLine:end-n-1,:)',X3Distorted{i}(1:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k');
            end
        end
        drawnow
        delete(h1)
        delete(h2)
        if  FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid
                delete(h3)
                delete(h4)
            end
            if  ToggleLowerFluid
                delete(h5)
                delete(h6)
            end
        end
        if  i == length(Time)
            i = 1;
        end
    end
else
    for i = 1:length(Time)
        if  ~FluidLoading || ~ShowHalfSpaces
            h1 = line(ax,X1Distorted{i}(:,1:GridLine:end),X3Distorted{i}(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(1:GridLine:end,:)',X3Distorted{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
        elseif FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid && ToggleLowerFluid
                h3 = line(ax,X1Distorted{i}(1:n+1,1:GridLine:end),X3Distorted{i}(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',UpperFluidColor);
                h4 = line(ax,X1Distorted{i}(1:GridLine:n+1,:)',X3Distorted{i}(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',UpperFluidColor);
                h5 = line(ax,X1Distorted{i}(end-n:end,1:GridLine:end),X3Distorted{i}(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',LowerFluidColor);
                h6 = line(ax,X1Distorted{i}(end-n:GridLine:end,:)',X3Distorted{i}(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',LowerFluidColor);
                h1 = line(ax,X1Distorted{i}(n+2:end-n-1,1:GridLine:end),X3Distorted{i}(n+2:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
                h2 = line(ax,X1Distorted{i}(n+2:GridLine:end-n-1,:)',X3Distorted{i}(n+2:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k');
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                h3 = line(ax,X1Distorted{i}(1:n+1,1:GridLine:end),X3Distorted{i}(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
                h4 = line(ax,X1Distorted{i}(1:GridLine:n+1,:)',X3Distorted{i}(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
                h1 = line(ax,X1Distorted{i}(n+2:end,1:GridLine:end),X3Distorted{i}(n+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
                h2 = line(ax,X1Distorted{i}(n+2:GridLine:end,:)',X3Distorted{i}(n+2:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                h5 = line(ax,X1Distorted{i}(end-n:end,1:GridLine:end),X3Distorted{i}(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
                h6 = line(ax,X1Distorted{i}(end-n:GridLine:end,:)',X3Distorted{i}(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
                h1 = line(ax,X1Distorted{i}(1:end-n-1,1:GridLine:end),X3Distorted{i}(1:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
                h2 = line(ax,X1Distorted{i}(1:GridLine:end-n-1,:)',X3Distorted{i}(1:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k');
            end
        end
        if  SuperLayerSize == 1 
            Frames(i) = getframe(f);
        else
            Frames{i} = getframe(f);
        end
        delete(h1)
        delete(h2)
        if  FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid
                delete(h3)
                delete(h4)
            end
            if  ToggleLowerFluid
                delete(h5)
                delete(h6)
            end
        end
    end
    if  SuperLayerSize > 1
        Frames = cell2mat(Frames);
    end
    try
        Movie = VideoWriter(fullfile(Directory,FileName));
        Movie.FrameRate = FrameRate;
        Movie.Quality = MovieQuality;
        open(Movie);
        writeVideo(Movie,Frames);
        close(Movie);
        delete(f)
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export movie')
        return
    end
end