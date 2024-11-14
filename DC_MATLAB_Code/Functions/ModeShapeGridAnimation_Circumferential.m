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
function ModeShapeGridAnimation_Circumferential(Hybrid,Layers,LayupString,Directory,FontSizeHeadLine,FontSizeAxesLabels,FrameRate,Frequency,GridLine,HeadLine,LineWidth,Material,Mode,Movie,FileName,MovieQuality,Ro,Ri,PropagationAngle,Gain,Repetitions,SymmetricSystem,Time,u,Undistorted,x1,r,p)
Cycles = 5; % for movie export

%#ok<*AGROW>
[X1,X3] = meshgrid(x1,r);
Ratio = (r(end)-r(1))/x1(end);
for i = 1:size(u,1)
    u1Max(i) = max(abs(real(u{i,1}(:,1))));
    u3Max(i) = max(abs(real(u{i,1}(:,2))));
end
u1Max = max(u1Max);
u3Max = max(u3Max);
if  u1Max > u3Max
    Compensation = Gain*x1(end)/40/u1Max;
else
    Compensation = Gain*(r(end)-r(1))/20/u3Max;
end
for j = 1:length(Time)
    for i = 1:size(X1,2)
        X1Distorted{j}(:,i) = X1(:,i)+real(u{j,i}(:,1))*Compensation;
        X3Distorted{j}(:,i) = X3(:,i)+real(u{j,i}(:,2))*Compensation*Ratio;
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
        if  ~contains(Mode,'SH')
            ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
        else
            ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
        end
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2),'\,$\times$\,',num2str(Ro-Ri),'\,mm ',replace(Material.Name,'_','\_'),' circumference'];
        ax.Title.String = String;
    end
else

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
if  Movie
    Frames(length(Time)*Cycles) = struct('cdata',[],'colormap',[]);
    h = msgbox('Recording cycle. Please wait...');
end
i = 0;
while i <= length(Time) && ishghandle(f) == 1
    i = i+1;
    h1 = line(ax,X1Distorted{i}(:,1:GridLine:end),X3Distorted{i}(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
    h2 = line(ax,X1Distorted{i}(1:GridLine:end,:)',X3Distorted{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
    drawnow
    if  Movie
        Frames(i) = getframe(f);
    end
    delete(h1)
    delete(h2)
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