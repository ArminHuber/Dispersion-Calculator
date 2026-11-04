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
function ModeShapeGrid_Plate(Hybrid,Layers,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,LayupString,Directory,FontSizeHeadLine,FontSizeAxesLabels,Frequency,GridLine,HeadLine,LineWidth,Material,Mode,FileName,PlateThickness,PropagationAngle,SamplesX3,Gain,u,Undistorted,x1,x3,p,ShowHalfSpaces,HalfSpaces,PNGresolution,Export,PDF,PNG,FrameRate,Movie,MovieQuality)
Cycles = 5; % for movie export

%#ok<*AGROW>
Length = size(u,3); % number of frames
[X1,X3] = meshgrid(x1,x3);
Ratio = (abs(x3(1))+abs(x3(end)))/x1(end);
u1Max = max(abs(u(:,:,:,1)),[],'all');
u3Max = max(abs(u(:,:,:,2)),[],'all');
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
X1Distorted = X1+u(:,:,:,1)*Compensation;
X3Distorted = X3+u(:,:,:,2)*Compensation*Ratio;
if  FrameRate
    f = figure('Icon',which('DC_Logo16.png'),'Name','Mode shape','MenuBar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
else
    f = figure('Icon',which('DC_Logo16.png'),'Name','Mode shape','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
end
ax = gca;
axis off
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
if  HeadLine
    if  contains(Mode,'SH')
        ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
    else
        ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
    end
    if  Hybrid
        MaterialName = 'hybrid';
    else
        MaterialName = replace(Material{1}.Name,'_','\_');
    end
    if  isempty(LayupString)
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(PlateThickness),'\,mm ',MaterialName,' plate'];
    else
        if  HeadLine == 1
            LayupString = '';
        end
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness),'\,mm ',MaterialName,LayupString,' plate'];
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
if  ~FluidLoading || ~ShowHalfSpaces
    h1 = line(ax,X1Distorted(:,1:GridLine:end,1),X3Distorted(:,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
    h2 = line(ax,X1Distorted(1:GridLine:end,:,1)',X3Distorted(1:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
elseif FluidLoading && ShowHalfSpaces
    if  ToggleUpperFluid && ToggleLowerFluid
        h3 = line(ax,X1Distorted(1:k+1,1:GridLine:end,1),X3Distorted(1:k+1,1:GridLine:end,1),'LineWidth',LineWidth,'Color',UpperFluidColor);
        h4 = line(ax,X1Distorted(1:GridLine:k+1,:,1)',X3Distorted(1:GridLine:k+1,:,1)','LineWidth',LineWidth,'Color',UpperFluidColor);
        h5 = line(ax,X1Distorted(end-k:end,1:GridLine:end,1),X3Distorted(end-k:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color',LowerFluidColor);
        h6 = line(ax,X1Distorted(end-k:GridLine:end,:,1)',X3Distorted(end-k:GridLine:end,:,1)','LineWidth',LineWidth,'Color',LowerFluidColor);
        h1 = line(ax,X1Distorted(k+2:end-k-1,1:GridLine:end,1),X3Distorted(k+2:end-k-1,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X1Distorted(k+2:GridLine:end-k-1,:,1)',X3Distorted(k+2:GridLine:end-k-1,:,1)','LineWidth',LineWidth,'Color','k');
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        h3 = line(ax,X1Distorted(1:k+1,1:GridLine:end,1),X3Distorted(1:k+1,1:GridLine:end,1),'LineWidth',LineWidth,'Color',[0 .5 1]);
        h4 = line(ax,X1Distorted(1:GridLine:k+1,:,1)',X3Distorted(1:GridLine:k+1,:,1)','LineWidth',LineWidth,'Color',[0 .5 1]);
        h1 = line(ax,X1Distorted(k+2:end,1:GridLine:end,1),X3Distorted(k+2:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X1Distorted(k+2:GridLine:end,:,1)',X3Distorted(k+2:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        h5 = line(ax,X1Distorted(end-k:end,1:GridLine:end,1),X3Distorted(end-k:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color',[0 .5 1]);
        h6 = line(ax,X1Distorted(end-k:GridLine:end,:,1)',X3Distorted(end-k:GridLine:end,:,1)','LineWidth',LineWidth,'Color',[0 .5 1]);
        h1 = line(ax,X1Distorted(1:end-k-1,1:GridLine:end,1),X3Distorted(1:end-k-1,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X1Distorted(1:GridLine:end-k-1,:,1)',X3Distorted(1:GridLine:end-k-1,:,1)','LineWidth',LineWidth,'Color','k');
    end
end
if  FrameRate
    if  Movie
        Frames(Length*Cycles) = struct('cdata',[],'colormap',[]);
        h = msgbox('Recording cycle. Please wait...');
    end
    i = 0;
    while i <= Length && ishghandle(f)
        i = i+1;
        if  ~FluidLoading || ~ShowHalfSpaces
            for j = 1:length(h1)
                n = GridLine*(j-1)+1;
                h1(j).XData = X1Distorted(:,n,i);
                h1(j).YData = X3Distorted(:,n,i);
            end
            for j = 1:length(h2)
                n = GridLine*(j-1)+1;
                h2(j).XData = X1Distorted(n,:,i);
                h2(j).YData = X3Distorted(n,:,i);
            end
        elseif FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid && ToggleLowerFluid
                for j = 1:length(h3)
                    n = GridLine*(j-1)+1;
                    h3(j).XData = X1Distorted(1:k+1,n,i);
                    h3(j).YData = X3Distorted(1:k+1,n,i);
                end
                for j = 1:length(h4)
                    n = GridLine*(j-1)+1;
                    h4(j).XData = X1Distorted(n,:,i);
                    h4(j).YData = X3Distorted(n,:,i);
                end
                for j = 1:length(h5)
                    n = GridLine*(j-1)+1;
                    h5(j).XData = X1Distorted(end-k:end,n,i);
                    h5(j).YData = X3Distorted(end-k:end,n,i);
                end
                for j = 1:length(h6)
                    n = GridLine*(j-1)+1;
                    h6(j).XData = X1Distorted(end-k+n-1,:,i);
                    h6(j).YData = X3Distorted(end-k+n-1,:,i);
                end
                for j = 1:length(h1)
                    n = GridLine*(j-1)+1;
                    h1(j).XData = X1Distorted(k+2:end-k-1,n,i);
                    h1(j).YData = X3Distorted(k+2:end-k-1,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X1Distorted(k+1+n,:,i);
                    h2(j).YData = X3Distorted(k+1+n,:,i);
                end
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                for j = 1:length(h3)
                    n = GridLine*(j-1)+1;
                    h3(j).XData = X1Distorted(1:k+1,n,i);
                    h3(j).YData = X3Distorted(1:k+1,n,i);
                end
                for j = 1:length(h4)
                    n = GridLine*(j-1)+1;
                    h4(j).XData = X1Distorted(n,:,i);
                    h4(j).YData = X3Distorted(n,:,i);
                end
                for j = 1:length(h1)
                    n = GridLine*(j-1)+1;
                    h1(j).XData = X1Distorted(k+2:end,n,i);
                    h1(j).YData = X3Distorted(k+2:end,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X1Distorted(k+1+n,:,i);
                    h2(j).YData = X3Distorted(k+1+n,:,i);
                end
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                for j = 1:length(h5)
                    n = GridLine*(j-1)+1;
                    h5(j).XData = X1Distorted(end-k:end,n,i);
                    h5(j).YData = X3Distorted(end-k:end,n,i);
                end
                for j = 1:length(h6)
                    n = GridLine*(j-1)+1;
                    h6(j).XData = X1Distorted(end-k+n-1,:,i);
                    h6(j).YData = X3Distorted(end-k+n-1,:,i);
                end
                for j = 1:length(h1)
                    n = GridLine*(j-1)+1;
                    h1(j).XData = X1Distorted(1:end-k-1,n,i);
                    h1(j).YData = X3Distorted(1:end-k-1,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X1Distorted(n,:,i);
                    h2(j).YData = X3Distorted(n,:,i);
                end
            end
        end
        drawnow
        if  Movie
            Frames(i) = getframe(f);
        end
        if  i == Length
            i = 0;
            if  Movie
                for n = 1:Cycles-1
                    for i = 1:Length
                        Frames(i+n*Length).cdata = Frames(i).cdata;
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
else
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
end
function output_txt = Cursor(~,~)
    output_txt = {[]};
end
end