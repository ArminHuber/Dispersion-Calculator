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
function ModeShapeGridAnimation(Hybrid,Layers,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,LayupString,Directory,FontSizeHeadLine,FontSizeAxesLabels,FrameRate,Frequency,GridLine,HeadLine,LineWidth,Material,Mode,Movie,FileName,MovieQuality,PlateThickness,PropagationAngle,SamplesX3,Gain,Repetitions,SymmetricSystem,Time,u,Undistorted,x1,x3,p,ShowHalfSpaces,HalfSpaces)
Cycles = 5; % for movie export

%#ok<*AGROW>
[X1,X3] = meshgrid(x1,x3);
Ratio = (abs(x3(1))+abs(x3(end)))/x1(end);
for i = 1:size(u,2)
    u1Max(i) = max(abs(real(u{1,i}(:,1))));
    u3Max(i) = max(abs(real(u{1,i}(:,2))));
end
u1Max = max(u1Max);
u3Max = max(u3Max);
if  FluidLoading && ShowHalfSpaces
    k = 1+2*HalfSpaces;
else
    k = 1;
end
if  u1Max > u3Max
    Compensation = Gain*x1(end)/80/k/u1Max;
else
    Compensation = Gain*(abs(x3(1))+abs(x3(end)))/40/k/u3Max;
end
for l = 1:length(Time)
    for i = 1:size(X1,2)
        for j = 1:size(X1,1)
            X1Distorted{l}(j,i) = X1(j,i)+real(u{j,i}(l,1))*Compensation;
            X3Distorted{l}(j,i) = X3(j,i)+real(u{j,i}(l,2))*Compensation*Ratio;
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
if  isempty(LayupString)
    if  HeadLine
        if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
            ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
        elseif contains(Mode,'SH')
            ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
        elseif contains(Mode,'Scholte')
            ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
        end
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(PlateThickness),'\,mm ',replace(Material.Name,'_','\_'),' plate'];
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
else
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
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
        elseif ~SymmetricSystem && HeadLine == 2
            if  Repetitions == 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
            elseif Repetitions > 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
            end
        elseif SymmetricSystem && HeadLine == 2
            if  Repetitions == 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
            elseif Repetitions > 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
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
        ax.Title.String = String;
    end
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
if  FluidLoading && ShowHalfSpaces
	k = HalfSpaces*SamplesX3*Layers;
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
        line(ax,X1(1:k+1,1:GridLine:end),X3(1:k+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(1:GridLine:k+1,:)',X3(1:GridLine:k+1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(end-k:end,1:GridLine:end),X3(end-k:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(end-k:GridLine:end,:)',X3(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(k+2:end-k-1,1:GridLine:end),X3(k+2:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(ax,X1(k+2:GridLine:end-k-1,:)',X3(k+2:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
    end
end
if  Movie
    Frames(length(Time)*Cycles) = struct('cdata',[],'colormap',[]);
    h = msgbox('Recording cycle. Please wait...');
end
i = 0;
while i <= length(Time) && ishghandle(f) == 1
    i = i+1;
    if  ~FluidLoading || ~ShowHalfSpaces
        h1 = line(ax,X1Distorted{i}(:,1:GridLine:end),X3Distorted{i}(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X1Distorted{i}(1:GridLine:end,:)',X3Distorted{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
    elseif FluidLoading && ShowHalfSpaces
        if  ToggleUpperFluid && ToggleLowerFluid
            h3 = line(ax,X1Distorted{i}(1:k+1,1:GridLine:end),X3Distorted{i}(1:k+1,1:GridLine:end),'LineWidth',LineWidth,'Color',UpperFluidColor);
            h4 = line(ax,X1Distorted{i}(1:GridLine:k+1,:)',X3Distorted{i}(1:GridLine:k+1,:)','LineWidth',LineWidth,'Color',UpperFluidColor);
            h5 = line(ax,X1Distorted{i}(end-k:end,1:GridLine:end),X3Distorted{i}(end-k:end,1:GridLine:end),'LineWidth',LineWidth,'Color',LowerFluidColor);
            h6 = line(ax,X1Distorted{i}(end-k:GridLine:end,:)',X3Distorted{i}(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',LowerFluidColor);
            h1 = line(ax,X1Distorted{i}(k+2:end-k-1,1:GridLine:end),X3Distorted{i}(k+2:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(k+2:GridLine:end-k-1,:)',X3Distorted{i}(k+2:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color','k');
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            h3 = line(ax,X1Distorted{i}(1:k+1,1:GridLine:end),X3Distorted{i}(1:k+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
            h4 = line(ax,X1Distorted{i}(1:GridLine:k+1,:)',X3Distorted{i}(1:GridLine:k+1,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
            h1 = line(ax,X1Distorted{i}(k+2:end,1:GridLine:end),X3Distorted{i}(k+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(k+2:GridLine:end,:)',X3Distorted{i}(k+2:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            h5 = line(ax,X1Distorted{i}(end-k:end,1:GridLine:end),X3Distorted{i}(end-k:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
            h6 = line(ax,X1Distorted{i}(end-k:GridLine:end,:)',X3Distorted{i}(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
            h1 = line(ax,X1Distorted{i}(1:end-k-1,1:GridLine:end),X3Distorted{i}(1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(1:GridLine:end-k-1,:)',X3Distorted{i}(1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color','k');
        end
    end
    drawnow
    if  Movie
        Frames(i) = getframe(f);
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
    if  i == length(Time)
        i = 0;
        if  Movie
            for n = 1:Cycles-1
                for i = 1:length(Time)
                    Frames(i+n*length(Time)).cdata = Frames(i).cdata;
                end
            end
            try
                Movie = VideoWriter(fullfile(Directory,FileName));
                Movie.FrameRate = FrameRate;
                Movie.Quality = MovieQuality;
                close(h)
                h = msgbox('Writing movie. Please wait...');
                open(Movie);
                writeVideo(Movie,Frames);
                close(Movie);
                close([f h])
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export movie')
                return
            end
            break
        end
    end
end