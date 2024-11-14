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
function ModeShapeGridAnimation_Rod_Pipe(Geometry,Plane,Hybrid,Layers,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,LayupString,Directory,FontSizeHeadLine,FontSizeAxesLabels,FrameRate,Frequency,GridLine,HeadLine,LineWidth,Material,Mode,Movie,FileName,MovieQuality,Ro,Ri,PropagationAngle,SamplesR,Gain,Repetitions,SymmetricSystem,Time,u,ua,ub,Undistorted,x1,r,Thetaa,Thetab,p,ShowHalfSpace,HalfSpaces)
Cycles = 5; % for movie export
Azimuth = -82;
Elevation = 8;

%#ok<*AGROW>
if  strcmp(Plane,'r-z')
    [X1,X3] = meshgrid(x1,r);
    if  ~contains(Mode,'T')
        Ratio = (r(end)-r(1))/x1(end);
        for i = 1:size(u,1)
            u1Max(i) = max(abs(real(u{i,1}(:,1))));
            u3Max(i) = max(abs(real(u{i,1}(:,3))));
        end
        u1Max = max(u1Max);
        u3Max = max(u3Max);
        if  ToggleOuterFluid && ShowHalfSpace
            k = 1+HalfSpaces;
        else
            k = 1;
        end
        if  u1Max > u3Max
            Compensation = Gain*x1(end)/40/k/u1Max;
        else
            Compensation = Gain*(r(end)-r(1))/20/k/u3Max;
        end
        for j = 1:length(Time)
            for i = 1:size(X1,2)
                X1Distorted{j}(:,i) = X1(:,i)+real(u{j,i}(:,1))*Compensation;
                X3Distorted{j}(:,i) = X3(:,i)+real(u{j,i}(:,3))*Compensation*Ratio;
            end
        end
    else
        X1Distorted = repmat({X1},1,length(Time));
        X3Distorted = repmat({X3},1,length(Time));
    end
else
    X1a = repmat(r,[1,length(Thetaa)]);
    X2a = repmat(Thetaa',[length(r),1]);
    X3a = zeros(size(X1a));
    X1b = repmat(r(end),[length(x1),length(Thetab)]);
    X2b = repmat(Thetab',[length(x1),1]);
    X3b = repmat(x1',[1,length(Thetab)]);
    for i = 1:size(ua,1)
        uMax(i) = max(abs(real(ua{i,1})),[],'all');
    end
    uMax = max(uMax);
    Compensation = Gain*r(end)/20/uMax;
    if  contains(Mode,'F')
        for j = 1:length(Time)
            for i = 1:size(X1a,2)
                X1Distorteda{j}(:,i) = X1a(:,i)+real(ua{j,i}(:,3))*Compensation; % r
                X2Distorteda{j}(:,i) = X2a(:,i)+real(ua{j,i}(:,2))./r*Compensation; % theta[rad] = b[m]/r[m]
                X3Distorteda{j}(:,i) = real(ua{j,i}(:,1))*Compensation; % z
            end
            for i = 1:size(X1b,2)
                X1Distortedb{j}(:,i) = X1b(:,i)+real(ub{j,i}(:,3))*Compensation;
                X2Distortedb{j}(:,i) = X2b(:,i)+real(ub{j,i}(:,2))/r(end)*Compensation;
                X3Distortedb{j}(:,i) = X3b(:,i)+real(ub{j,i}(:,1))*Compensation;
            end
        end
    else
        for j = 1:length(Time)
            X1Distorteda{j} = repmat(X1a(:,1)+real(ua{j}(:,3))*Compensation,[1,size(X1a,2)]);
            X2Distorteda{j} = X2a(1,:)+real(ua{j}(:,2))./r*Compensation;
            X3Distorteda{j} = repmat(real(ua{j}(:,1))*Compensation,[1,size(X1a,2)]);
            X1Distortedb{j} = repmat(X1b(:,1)+real(ub{j}(:,3))*Compensation,[1,size(X1b,2)]);
            X2Distortedb{j} = X2b(1,:)+real(ub{j}(:,2))/r(end)*Compensation;
            X3Distortedb{j} = repmat(X3b(:,1)+real(ub{j}(:,1))*Compensation,[1,size(X1b,2)]);
        end
    end
    [X1a,X2a] = pol2cart(X2a,X1a);
    [X1b,X2b] = pol2cart(X2b,X1b);
    for j = 1:length(Time)
        [X1Distorteda{j},X2Distorteda{j}] = pol2cart(X2Distorteda{j},X1Distorteda{j});
        [X1Distortedb{j},X2Distortedb{j}] = pol2cart(X2Distortedb{j},X1Distortedb{j});
    end
end
f = figure('Name','Mode shape','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
ax = gca;
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
if  isempty(LayupString)
    if  HeadLine
        if  ~contains(Mode,'Scholte')
            ModeName = [Mode(1),'(',num2str(p(1)),',',num2str(p(2)),')'];
        else
            ModeName = [Mode(1),'(',num2str(p(1)),',',num2str(p(2)),')$^\mathrm{Scholte}$'];
        end
        if  strcmp(Geometry,'Rod')
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
        elseif strcmp(Geometry,'Pipe')
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2),'\,$\times$\,',num2str(Ro-Ri),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
        end
        if  FluidLoading
            if  strcmp(Geometry,'Rod')
                String = append(String,' in ',replace(OuterFluid.Name,'_','\_'));
            elseif strcmp(Geometry,'Pipe')
                if  ToggleOuterFluid && ToggleInnerFluid
                    String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/',replace(InnerFluid.Name,'_','\_'));
                elseif ToggleOuterFluid && ~ToggleInnerFluid
                    String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/vacuum');
                elseif ~ToggleOuterFluid && ToggleInnerFluid
                    String = append(String,' in vacuum/',replace(InnerFluid.Name,'_','\_'));
                end
            end
        end
        ax.Title.String = String;
    end
else

end
if  strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
    k = HalfSpaces*SamplesR;
    if  strcmp(OuterFluid.Name,InnerFluid.Name)
        OuterFluidColor = [0 .5 1];
        InnerFluidColor = [0 .5 1];
    else
        if  OuterFluid.Density*OuterFluid.Velocity > InnerFluid.Density*InnerFluid.Velocity
            OuterFluidColor = [0 .5 1];
            InnerFluidColor = [0 .7 1];
        else
            OuterFluidColor = [0 .7 1];
            InnerFluidColor = [0 .5 1];
        end
    end
end
if  Undistorted
    if  strcmp(Plane,'r-z')
        if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
            line(ax,X1(:,1:GridLine:end),X3(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(1:GridLine:end,:)',X3(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
            line(ax,X1(1:SamplesR,1:GridLine:end),X3(1:SamplesR,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(1:GridLine:SamplesR,:)',X3(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(end-k:end,1:GridLine:end),X3(end-k:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(end-k:GridLine:end,:)',X3(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(SamplesR+1:end-k-1,1:GridLine:end),X3(SamplesR+1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(SamplesR+1:GridLine:end-k-1,:)',X3(SamplesR+1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
            line(ax,X1(SamplesR+2:end,1:GridLine:end),X3(SamplesR+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(SamplesR+2:GridLine:end,:)',X3(SamplesR+2:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(1:SamplesR+1,1:GridLine:end),X3(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(1:GridLine:SamplesR+1,:)',X3(1:GridLine:SamplesR+1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
            line(ax,X1(1:SamplesR,1:GridLine:end),X3(1:SamplesR,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(1:GridLine:SamplesR,:)',X3(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(SamplesR+1:end,1:GridLine:end),X3(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X1(SamplesR+1:GridLine:end,:)',X3(SamplesR+1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
    else
        if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
            line(ax,X3a(:,1:GridLine:end),X2a(:,1:GridLine:end),X1a(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(1:GridLine:end,:)',X2a(1:GridLine:end,:)',X1a(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3b(:,GridLine:GridLine:end),X2b(:,GridLine:GridLine:end),X1b(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3b(1:GridLine:end,:)',X2b(1:GridLine:end,:)',X1b(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
            line(ax,X3a(1:SamplesR+1,1:GridLine:end),X2a(1:SamplesR+1,1:GridLine:end),X1a(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(1:GridLine:SamplesR,:)',X2a(1:GridLine:SamplesR,:)',X1a(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(end-k-1:end,1:GridLine:end),X2a(end-k-1:end,1:GridLine:end),X1a(end-k-1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(end-k:GridLine:end,:)',X2a(end-k:GridLine:end,:)',X1a(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(SamplesR+1:end-k-1,1:GridLine:end),X2a(SamplesR+1:end-k-1,1:GridLine:end),X1a(SamplesR+1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(SamplesR+1:GridLine:end-k-1,:)',X2a(SamplesR+1:GridLine:end-k-1,:)',X1a(SamplesR+1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3b(:,GridLine:GridLine:end),X2b(:,GridLine:GridLine:end),X1b(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3b(1:GridLine:end,:)',X2b(1:GridLine:end,:)',X1b(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
            line(ax,X3a(SamplesR+1:end,1:GridLine:end),X2a(SamplesR+1:end,1:GridLine:end),X1a(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(SamplesR+2:GridLine:end,:)',X2a(SamplesR+2:GridLine:end,:)',X1a(SamplesR+2:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(1:SamplesR+1,1:GridLine:end),X2a(1:SamplesR+1,1:GridLine:end),X1a(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(1:GridLine:SamplesR+1,:)',X2a(1:GridLine:SamplesR+1,:)',X1a(1:GridLine:SamplesR+1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3b(:,GridLine:GridLine:end),X2b(:,GridLine:GridLine:end),X1b(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3b(1:GridLine:end,:)',X2b(1:GridLine:end,:)',X1b(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
            line(ax,X3a(1:SamplesR+1,1:GridLine:end),X2a(1:SamplesR+1,1:GridLine:end),X1a(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(1:GridLine:SamplesR,:)',X2a(1:GridLine:SamplesR,:)',X1a(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(SamplesR+1:end,1:GridLine:end),X2a(SamplesR+1:end,1:GridLine:end),X1a(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3a(SamplesR+1:GridLine:end,:)',X2a(SamplesR+1:GridLine:end,:)',X1a(SamplesR+1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3b(:,GridLine:GridLine:end),X2b(:,GridLine:GridLine:end),X1b(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
            line(ax,X3b(1:GridLine:end,:)',X2b(1:GridLine:end,:)',X1b(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        end
    end
end
if  strcmp(Plane,'r-z')
    axis off
    text(-.02,.5-.05*FontSizeAxesLabels/30,'$r$','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels,'rotation',90)
    text(.5-.155*FontSizeAxesLabels/30,-.01,'Propagation direction ($z$)','Units','Normalized','Interpreter','latex','FontSize',FontSizeAxesLabels)
    ax.XLim = [-.075*x1(end) x1(end)+.075*x1(end)];
else
    rotate3d on
    axis equal tight off
    view(Azimuth,Elevation)
end
if  Movie
    Frames(length(Time)*Cycles) = struct('cdata',[],'colormap',[]);
    h = msgbox('Recording cycle. Please wait...');
end
i = 0;
while i <= length(Time) && ishghandle(f) == 1
    i = i+1;
    if  strcmp(Plane,'r-z')
        if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
            h1 = line(ax,X1Distorted{i}(:,1:GridLine:end),X3Distorted{i}(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(1:GridLine:end,:)',X3Distorted{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
        elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
            h3 = line(ax,X1Distorted{i}(1:SamplesR,1:GridLine:end),X3Distorted{i}(1:SamplesR,1:GridLine:end),'LineWidth',LineWidth,'Color',InnerFluidColor);
            h4 = line(ax,X1Distorted{i}(1:GridLine:SamplesR,:)',X3Distorted{i}(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',InnerFluidColor);
            h5 = line(ax,X1Distorted{i}(end-k:end,1:GridLine:end),X3Distorted{i}(end-k:end,1:GridLine:end),'LineWidth',LineWidth,'Color',OuterFluidColor);
            h6 = line(ax,X1Distorted{i}(end-k:GridLine:end,:)',X3Distorted{i}(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',OuterFluidColor);
            h1 = line(ax,X1Distorted{i}(SamplesR+1:end-k-1,1:GridLine:end),X3Distorted{i}(SamplesR+1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(SamplesR+1:GridLine:end-k-1,:)',X3Distorted{i}(SamplesR+1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color','k');
        elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
            h5 = line(ax,X1Distorted{i}(SamplesR+2:end,1:GridLine:end),X3Distorted{i}(SamplesR+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
            h6 = line(ax,X1Distorted{i}(SamplesR+2:GridLine:end,:)',X3Distorted{i}(SamplesR+2:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
            h1 = line(ax,X1Distorted{i}(1:SamplesR+1,1:GridLine:end),X3Distorted{i}(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(1:GridLine:SamplesR+1,:)',X3Distorted{i}(1:GridLine:SamplesR+1,:)','LineWidth',LineWidth,'Color','k');
        elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
            h3 = line(ax,X1Distorted{i}(1:SamplesR,1:GridLine:end),X3Distorted{i}(1:SamplesR,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
            h4 = line(ax,X1Distorted{i}(1:GridLine:SamplesR,:)',X3Distorted{i}(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
            h1 = line(ax,X1Distorted{i}(SamplesR+1:end,1:GridLine:end),X3Distorted{i}(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X1Distorted{i}(SamplesR+1:GridLine:end,:)',X3Distorted{i}(SamplesR+1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
        end
    else
        if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
            h1 = line(ax,X3Distorteda{i}(:,1:GridLine:end),X2Distorteda{i}(:,1:GridLine:end),X1Distorteda{i}(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X3Distorteda{i}(1:GridLine:end,:)',X2Distorteda{i}(1:GridLine:end,:)',X1Distorteda{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
            h7 = line(ax,X3Distortedb{i}(:,GridLine:GridLine:end),X2Distortedb{i}(:,GridLine:GridLine:end),X1Distortedb{i}(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h8 = line(ax,X3Distortedb{i}(1:GridLine:end,:)',X2Distortedb{i}(1:GridLine:end,:)',X1Distortedb{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
        elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
            h3 = line(ax,X3Distorteda{i}(1:SamplesR+1,1:GridLine:end),X2Distorteda{i}(1:SamplesR+1,1:GridLine:end),X1Distorteda{i}(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',InnerFluidColor);
            h4 = line(ax,X3Distorteda{i}(1:GridLine:SamplesR,:)',X2Distorteda{i}(1:GridLine:SamplesR,:)',X1Distorteda{i}(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',InnerFluidColor);
            h5 = line(ax,X3Distorteda{i}(end-k-1:end,1:GridLine:end),X2Distorteda{i}(end-k-1:end,1:GridLine:end),X1Distorteda{i}(end-k-1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',OuterFluidColor);
            h6 = line(ax,X3Distorteda{i}(end-k:GridLine:end,:)',X2Distorteda{i}(end-k:GridLine:end,:)',X1Distorteda{i}(end-k:GridLine:end,:)','LineWidth',LineWidth,'Color',OuterFluidColor);
            h1 = line(ax,X3Distorteda{i}(SamplesR+1:end-k-1,1:GridLine:end),X2Distorteda{i}(SamplesR+1:end-k-1,1:GridLine:end),X1Distorteda{i}(SamplesR+1:end-k-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X3Distorteda{i}(SamplesR+1:GridLine:end-k-1,:)',X2Distorteda{i}(SamplesR+1:GridLine:end-k-1,:)',X1Distorteda{i}(SamplesR+1:GridLine:end-k-1,:)','LineWidth',LineWidth,'Color','k');
            h7 = line(ax,X3Distortedb{i}(:,GridLine:GridLine:end),X2Distortedb{i}(:,GridLine:GridLine:end),X1Distortedb{i}(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',OuterFluidColor);
            h8 = line(ax,X3Distortedb{i}(1:GridLine:end,:)',X2Distortedb{i}(1:GridLine:end,:)',X1Distortedb{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color',OuterFluidColor);
        elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
            h5 = line(ax,X3Distorteda{i}(SamplesR+1:end,1:GridLine:end),X2Distorteda{i}(SamplesR+1:end,1:GridLine:end),X1Distorteda{i}(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
            h6 = line(ax,X3Distorteda{i}(SamplesR+2:GridLine:end,:)',X2Distorteda{i}(SamplesR+2:GridLine:end,:)',X1Distorteda{i}(SamplesR+2:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
            h1 = line(ax,X3Distorteda{i}(1:SamplesR+1,1:GridLine:end),X2Distorteda{i}(1:SamplesR+1,1:GridLine:end),X1Distorteda{i}(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X3Distorteda{i}(1:GridLine:SamplesR+1,:)',X2Distorteda{i}(1:GridLine:SamplesR+1,:)',X1Distorteda{i}(1:GridLine:SamplesR+1,:)','LineWidth',LineWidth,'Color','k');
            h7 = line(ax,X3Distortedb{i}(:,GridLine:GridLine:end),X2Distortedb{i}(:,GridLine:GridLine:end),X1Distortedb{i}(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
            h8 = line(ax,X3Distortedb{i}(1:GridLine:end,:)',X2Distortedb{i}(1:GridLine:end,:)',X1Distortedb{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
        elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
            h3 = line(ax,X3Distorteda{i}(1:SamplesR+1,1:GridLine:end),X2Distorteda{i}(1:SamplesR+1,1:GridLine:end),X1Distorteda{i}(1:SamplesR+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1]);
            h4 = line(ax,X3Distorteda{i}(1:GridLine:SamplesR,:)',X2Distorteda{i}(1:GridLine:SamplesR,:)',X1Distorteda{i}(1:GridLine:SamplesR,:)','LineWidth',LineWidth,'Color',[0 .5 1]);
            h1 = line(ax,X3Distorteda{i}(SamplesR+1:end,1:GridLine:end),X2Distorteda{i}(SamplesR+1:end,1:GridLine:end),X1Distorteda{i}(SamplesR+1:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h2 = line(ax,X3Distorteda{i}(SamplesR+1:GridLine:end,:)',X2Distorteda{i}(SamplesR+1:GridLine:end,:)',X1Distorteda{i}(SamplesR+1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
            h7 = line(ax,X3Distortedb{i}(:,GridLine:GridLine:end),X2Distortedb{i}(:,GridLine:GridLine:end),X1Distortedb{i}(:,GridLine:GridLine:end),'LineWidth',LineWidth,'Color','k');
            h8 = line(ax,X3Distortedb{i}(1:GridLine:end,:)',X2Distortedb{i}(1:GridLine:end,:)',X1Distortedb{i}(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k');
        end
    end
    drawnow
    if  Movie
        Frames(i) = getframe(f);
    end
    delete(h1)
    delete(h2)
    if  ToggleOuterFluid && ShowHalfSpace
        delete(h5)
        delete(h6)
    end
    if  strcmp(Geometry,'Pipe') && ToggleInnerFluid
        delete(h3)
        delete(h4)
    end
    if  strcmp(Plane,'r-theta')
        delete(h7)
        delete(h8)
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