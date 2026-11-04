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
function ModeShapeGrid_Circumferential(Hybrid,LayupString,Directory,FontSizeHeadLine,FontSizeAxesLabels,Frequency,GridLine,HeadLine,LineWidth,Material,Mode,FileName,Ro,Ri,Gain,u,Undistorted,x1,r,p,PNGresolution,Export,PDF,PNG,FrameRate,Movie,MovieQuality)
Cycles = 5; % for movie export

%#ok<*AGROW>
Length = size(u,3); % number of frames
[X1,X3] = meshgrid(x1,r);
Ratio = (r(end)-r(1))/x1(end);
u1Max = max(abs(u(:,:,:,1)),[],'all');
u3Max = max(abs(u(:,:,:,2)),[],'all');
if  u1Max > u3Max
    Compensation = Gain*x1(end)/40/u1Max;
else
    Compensation = Gain*(r(end)-r(1))/20/u3Max;
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
    if  ~contains(Mode,'SH')
        ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
    else
        ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
    end
    if  Hybrid
        MaterialName = 'hybrid';
    else
        MaterialName = replace(Material{1}.Name,'_','\_');
    end
    if  HeadLine == 1
        LayupString = '';
    end
    String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2),'\,$\times$\,',num2str(Ro-Ri),'\,mm ',MaterialName,LayupString,' circumference'];
    ax.Title.String = String;
end
if  ~contains(Mode,'SH')
    text(.5-.155*FontSizeAxesLabels/30,.05,'Propagation direction ($\theta$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels)
    text(-.02,.5-.14*FontSizeAxesLabels/30,'Thickness ($r$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels,'rotation',90);
else
    text(.5-.125*FontSizeAxesLabels/30,-.073,'Shear horizontal ($z$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels)
    text(-.02,.5-.18*FontSizeAxesLabels/30,'Thickness ($r$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels,'rotation',90);
end
ax.XLim = [-.075*x1(end) x1(end)+.075*x1(end)];
if  Undistorted
    line(ax,X1(:,1:GridLine:end),X3(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0]);
    line(ax,X1(1:GridLine:end,:)',X3(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
end
h1 = line(ax,X1Distorted(:,1:GridLine:end,1),X3Distorted(:,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
h2 = line(ax,X1Distorted(1:GridLine:end,:,1)',X3Distorted(1:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
if  FrameRate
    if  Movie
        Frames(Length*Cycles) = struct('cdata',[],'colormap',[]);
        h = msgbox('Recording cycle. Please wait...');
    end
    i = 0;
    while i <= Length && ishghandle(f)
        i = i+1;
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