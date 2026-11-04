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
function ModeShapeGrid_Rod_Pipe(Hybrid,LayupString,Geometry,Plane,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Directory,FontSizeHeadLine,FontSizeAxesLabels,Frequency,GridLine,HeadLine,LineWidth,Material,Mode,FileName,Ro,Ri,SamplesR,Gain,u,ua,ub,Undistorted,x1,r,Thetaa,Thetab,p,ShowHalfSpace,HalfSpaces,PNGresolution,Export,PDF,PNG,FrameRate,Movie,MovieQuality)
Cycles = 5; % for movie export
Azimuth = -82;
Elevation = 8;

%#ok<*AGROW>
if  strcmp(Plane,'r-z')
    Length = size(u,3); % number of frames
    [X1,X3] = meshgrid(x1,r);
    if  ~contains(Mode,'T')
        Ratio = (r(end)-r(1))/x1(end);
        u1Max = max(abs(u(:,:,:,1)),[],'all');
        u3Max = max(abs(u(:,:,:,3)),[],'all');
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
        X1Distorted = X1+u(:,:,:,1)*Compensation;
        X3Distorted = X3+u(:,:,:,3)*Compensation*Ratio;
    else
        X1Distorted = X1;
        X3Distorted = X3;
    end
else
    Length = size(ua,3);
    X1a = repmat(r,1,length(Thetaa));
    X2a = repmat(Thetaa,length(r),1);
    X3a = zeros(size(X1a));
    X1b = repmat(r(end),length(x1),length(Thetab));
    X2b = repmat(Thetab,length(x1),1);
    X3b = repmat(x1',1,length(Thetab));
    uaMax = max(abs(ua),[],'all');
    ubMax = max(abs(ub),[],'all');
    if  uaMax > ubMax
        Compensation = Gain*r(end)/20/uaMax;
    else
        Compensation = Gain*r(end)/20/ubMax;
    end
    X1Distorteda = X1a+ua(:,:,:,3)*Compensation; % r
    X2Distorteda = X2a+ua(:,:,:,2)*Compensation./r; % theta[rad] = b[m]/r[m]
    X3Distorteda = ua(:,:,:,1)*Compensation; % z
    X1Distortedb = X1b+ub(:,:,:,3)*Compensation;
    X2Distortedb = X2b+ub(:,:,:,2)*Compensation/r(end);
    X3Distortedb = X3b+ub(:,:,:,1)*Compensation;
    [X1a,X2a] = pol2cart(X2a,X1a);
    [X1b,X2b] = pol2cart(X2b,X1b);
    [X1Distorteda,X2Distorteda] = pol2cart(X2Distorteda,X1Distorteda);
    [X1Distortedb,X2Distortedb] = pol2cart(X2Distortedb,X1Distortedb);
end
if  FrameRate
    f = figure('Icon',which('DC_Logo16.png'),'Name','Mode shape','MenuBar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
else
    f = figure('Icon',which('DC_Logo16.png'),'Name','Mode shape','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
end
ax = gca;
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
if  HeadLine
    ModeName = [Mode(1),'(',num2str(p(1)),',',num2str(p(2)),')'];
    if  Hybrid
        MaterialName = 'hybrid';
    else
        MaterialName = replace(Material{1}.Name,'_','\_');
    end
    if  HeadLine == 1
        LayupString = '';
    end
    if  strcmp(Geometry,'Rod')
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2),'\,mm ',MaterialName,LayupString,' rod'];
    elseif strcmp(Geometry,'Pipe')
        String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2),'\,$\times$\,',num2str(Ro-Ri),'\,mm ',MaterialName,LayupString,' pipe'];
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
if  strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
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
k = HalfSpaces*SamplesR;
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
        elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
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
        elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
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
    if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
        h1 = line(ax,X1Distorted(:,1:GridLine:end,1),X3Distorted(:,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X1Distorted(1:GridLine:end,:,1)',X3Distorted(1:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
    elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
        h3 = line(ax,X1Distorted(1:SamplesR,1:GridLine:end,1),X3Distorted(1:SamplesR,1:GridLine:end,1),'LineWidth',LineWidth,'Color',InnerFluidColor);
        h4 = line(ax,X1Distorted(1:GridLine:SamplesR,:,1)',X3Distorted(1:GridLine:SamplesR,:,1)','LineWidth',LineWidth,'Color',InnerFluidColor);
        h5 = line(ax,X1Distorted(end-k:end,1:GridLine:end,1),X3Distorted(end-k:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color',OuterFluidColor);
        h6 = line(ax,X1Distorted(end-k:GridLine:end,:,1)',X3Distorted(end-k:GridLine:end,:,1)','LineWidth',LineWidth,'Color',OuterFluidColor);
        h1 = line(ax,X1Distorted(SamplesR+1:end-k-1,1:GridLine:end,1),X3Distorted(SamplesR+1:end-k-1,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X1Distorted(SamplesR+1:GridLine:end-k-1,:,1)',X3Distorted(SamplesR+1:GridLine:end-k-1,:,1)','LineWidth',LineWidth,'Color','k');
    elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
        h5 = line(ax,X1Distorted(end-k:end,1:GridLine:end,1),X3Distorted(end-k:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color',[0 .5 1]);
        h6 = line(ax,X1Distorted(end-k:GridLine:end,:,1)',X3Distorted(end-k:GridLine:end,:,1)','LineWidth',LineWidth,'Color',[0 .5 1]);
        h1 = line(ax,X1Distorted(1:end-k-1,1:GridLine:end,1),X3Distorted(1:end-k-1,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X1Distorted(1:GridLine:end-k-1,:,1)',X3Distorted(1:GridLine:end-k-1,:,1)','LineWidth',LineWidth,'Color','k');
    elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
        h3 = line(ax,X1Distorted(1:SamplesR,1:GridLine:end,1),X3Distorted(1:SamplesR,1:GridLine:end,1),'LineWidth',LineWidth,'Color',[0 .5 1]);
        h4 = line(ax,X1Distorted(1:GridLine:SamplesR,:,1)',X3Distorted(1:GridLine:SamplesR,:,1)','LineWidth',LineWidth,'Color',[0 .5 1]);
        h1 = line(ax,X1Distorted(SamplesR+1:end,1:GridLine:end,1),X3Distorted(SamplesR+1:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X1Distorted(SamplesR+1:GridLine:end,:,1)',X3Distorted(SamplesR+1:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
    end
else
    if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
        h1 = line(ax,X3Distorteda(:,1:GridLine:end,1),X2Distorteda(:,1:GridLine:end,1),X1Distorteda(:,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X3Distorteda(1:GridLine:end,:,1)',X2Distorteda(1:GridLine:end,:,1)',X1Distorteda(1:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
        h7 = line(ax,X3Distortedb(:,GridLine:GridLine:end,1),X2Distortedb(:,GridLine:GridLine:end,1),X1Distortedb(:,GridLine:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h8 = line(ax,X3Distortedb(1:GridLine:end,:,1)',X2Distortedb(1:GridLine:end,:,1)',X1Distortedb(1:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
    elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
        h3 = line(ax,X3Distorteda(1:SamplesR,1:GridLine:end,1),X2Distorteda(1:SamplesR,1:GridLine:end,1),X1Distorteda(1:SamplesR,1:GridLine:end,1),'LineWidth',LineWidth,'Color',InnerFluidColor);
        h4 = line(ax,X3Distorteda(1:GridLine:SamplesR,:,1)',X2Distorteda(1:GridLine:SamplesR,:,1)',X1Distorteda(1:GridLine:SamplesR,:,1)','LineWidth',LineWidth,'Color',InnerFluidColor);
        h5 = line(ax,X3Distorteda(end-k:end,1:GridLine:end,1),X2Distorteda(end-k:end,1:GridLine:end,1),X1Distorteda(end-k:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color',OuterFluidColor);
        h6 = line(ax,X3Distorteda(end-k:GridLine:end,:,1)',X2Distorteda(end-k:GridLine:end,:,1)',X1Distorteda(end-k:GridLine:end,:,1)','LineWidth',LineWidth,'Color',OuterFluidColor);
        h1 = line(ax,X3Distorteda(SamplesR+1:end-k-1,1:GridLine:end,1),X2Distorteda(SamplesR+1:end-k-1,1:GridLine:end,1),X1Distorteda(SamplesR+1:end-k-1,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X3Distorteda(SamplesR+1:GridLine:end-k-1,:,1)',X2Distorteda(SamplesR+1:GridLine:end-k-1,:,1)',X1Distorteda(SamplesR+1:GridLine:end-k-1,:,1)','LineWidth',LineWidth,'Color','k');
        h7 = line(ax,X3Distortedb(:,GridLine:GridLine:end,1),X2Distortedb(:,GridLine:GridLine:end,1),X1Distortedb(:,GridLine:GridLine:end,1),'LineWidth',LineWidth,'Color',OuterFluidColor);
        h8 = line(ax,X3Distortedb(1:GridLine:end,:,1)',X2Distortedb(1:GridLine:end,:,1)',X1Distortedb(1:GridLine:end,:,1)','LineWidth',LineWidth,'Color',OuterFluidColor);
    elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
        h5 = line(ax,X3Distorteda(end-k:end,1:GridLine:end,1),X2Distorteda(end-k:end,1:GridLine:end,1),X1Distorteda(end-k:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color',[0 .5 1]);
        h6 = line(ax,X3Distorteda(end-k:GridLine:end,:,1)',X2Distorteda(end-k:GridLine:end,:,1)',X1Distorteda(end-k:GridLine:end,:,1)','LineWidth',LineWidth,'Color',[0 .5 1]);
        h1 = line(ax,X3Distorteda(1:end-k-1,1:GridLine:end,1),X2Distorteda(1:end-k-1,1:GridLine:end,1),X1Distorteda(1:end-k-1,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X3Distorteda(1:GridLine:end-k-1,:,1)',X2Distorteda(1:GridLine:end-k-1,:,1)',X1Distorteda(1:GridLine:end-k-1,:,1)','LineWidth',LineWidth,'Color','k');
        h7 = line(ax,X3Distortedb(:,GridLine:GridLine:end,1),X2Distortedb(:,GridLine:GridLine:end,1),X1Distortedb(:,GridLine:GridLine:end,1),'LineWidth',LineWidth,'Color',[0 .5 1]);
        h8 = line(ax,X3Distortedb(1:GridLine:end,:,1)',X2Distortedb(1:GridLine:end,:,1)',X1Distortedb(1:GridLine:end,:,1)','LineWidth',LineWidth,'Color',[0 .5 1]);
    elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
        h3 = line(ax,X3Distorteda(1:SamplesR,1:GridLine:end,1,1),X2Distorteda(1:SamplesR,1:GridLine:end,1),X1Distorteda(1:SamplesR,1:GridLine:end,1),'LineWidth',LineWidth,'Color',[0 .5 1]);
        h4 = line(ax,X3Distorteda(1:GridLine:SamplesR,:,1)',X2Distorteda(1:GridLine:SamplesR,:,1)',X1Distorteda(1:GridLine:SamplesR,:,1)','LineWidth',LineWidth,'Color',[0 .5 1]);
        h1 = line(ax,X3Distorteda(SamplesR+1:end,1:GridLine:end,1),X2Distorteda(SamplesR+1:end,1:GridLine:end,1),X1Distorteda(SamplesR+1:end,1:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h2 = line(ax,X3Distorteda(SamplesR+1:GridLine:end,:,1)',X2Distorteda(SamplesR+1:GridLine:end,:,1)',X1Distorteda(SamplesR+1:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
        h7 = line(ax,X3Distortedb(:,GridLine:GridLine:end,1),X2Distortedb(:,GridLine:GridLine:end,1),X1Distortedb(:,GridLine:GridLine:end,1),'LineWidth',LineWidth,'Color','k');
        h8 = line(ax,X3Distortedb(1:GridLine:end,:,1)',X2Distortedb(1:GridLine:end,:,1)',X1Distortedb(1:GridLine:end,:,1)','LineWidth',LineWidth,'Color','k');
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
if  FrameRate
    if  Movie
        Frames(Length*Cycles) = struct('cdata',[],'colormap',[]);
        h = msgbox('Recording cycle. Please wait...');
    end
    i = 0;
    while i <= Length && ishghandle(f)
        i = i+1;
        if  strcmp(Plane,'r-z')
            if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
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
            elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
                for j = 1:length(h3)
                    n = GridLine*(j-1)+1;
                    h3(j).XData = X1Distorted(1:SamplesR,n,i);
                    h3(j).YData = X3Distorted(1:SamplesR,n,i);
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
                    h1(j).XData = X1Distorted(SamplesR+1:end-k-1,n,i);
                    h1(j).YData = X3Distorted(SamplesR+1:end-k-1,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X1Distorted(SamplesR+n,:,i);
                    h2(j).YData = X3Distorted(SamplesR+n,:,i);
                end
            elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
                for j = 1:length(h5)
                    n = GridLine*(j-1)+1;
                    h5(j).XData = X1Distorted(end-k:end,n,i);
                    h5(j).YData = X3Distorted(end-k:end,n,i);
                end
                for j = 1:length(h6)
                    n = GridLine*(j-1)+1;
                    h6(j).XData = X1Distorted(end-k-1+n,:,i);
                    h6(j).YData = X3Distorted(end-k-1+n,:,i);
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
            elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
                for j = 1:length(h3)
                    n = GridLine*(j-1)+1;
                    h3(j).XData = X1Distorted(1:SamplesR,n,i);
                    h3(j).YData = X3Distorted(1:SamplesR,n,i);
                end
                for j = 1:length(h4)
                    n = GridLine*(j-1)+1;
                    h4(j).XData = X1Distorted(n,:,i);
                    h4(j).YData = X3Distorted(n,:,i);
                end
                for j = 1:length(h1)
                    n = GridLine*(j-1)+1;
                    h1(j).XData = X1Distorted(SamplesR+1:end,n,i);
                    h1(j).YData = X3Distorted(SamplesR+1:end,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X1Distorted(SamplesR+n,:,i);
                    h2(j).YData = X3Distorted(SamplesR+n,:,i);
                end
            end
        else
            if  (strcmp(Geometry,'Pipe') && (~FluidLoading || (ToggleOuterFluid && ~ToggleInnerFluid && ~ShowHalfSpace))) || (strcmp(Geometry,'Rod') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)))
                for j = 1:length(h1)
                    n = GridLine*(j-1)+1;
                    h1(j).XData = X3Distorteda(:,n,i);
                    h1(j).YData = X2Distorteda(:,n,i);
                    h1(j).ZData = X1Distorteda(:,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X3Distorteda(n,:,i);
                    h2(j).YData = X2Distorteda(n,:,i);
                    h2(j).ZData = X1Distorteda(n,:,i);
                end
                for j = 1:length(h7)
                    n = GridLine*(j-1)+1;
                    h7(j).XData = X3Distortedb(:,GridLine-1+n,i);
                    h7(j).YData = X2Distortedb(:,GridLine-1+n,i);
                    h7(j).ZData = X1Distortedb(:,GridLine-1+n,i);
                end
                for j = 1:length(h8)
                    n = GridLine*(j-1)+1;
                    h8(j).XData = X3Distortedb(n,:,i);
                    h8(j).YData = X2Distortedb(n,:,i);
                    h8(j).ZData = X1Distortedb(n,:,i);
                end
            elseif strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
                for j = 1:length(h3)
                    n = GridLine*(j-1)+1;
                    h3(j).XData = X3Distorteda(1:SamplesR,n,i);
                    h3(j).YData = X2Distorteda(1:SamplesR,n,i);
                    h3(j).ZData = X1Distorteda(1:SamplesR,n,i);
                end
                for j = 1:length(h4)
                    n = GridLine*(j-1)+1;
                    h4(j).XData = X3Distorteda(n,:,i);
                    h4(j).YData = X2Distorteda(n,:,i);
                    h4(j).ZData = X1Distorteda(n,:,i);
                end
                for j = 1:length(h5)
                    n = GridLine*(j-1)+1;
                    h5(j).XData = X3Distorteda(end-k:end,n,i);
                    h5(j).YData = X2Distorteda(end-k:end,n,i);
                    h5(j).ZData = X1Distorteda(end-k:end,n,i);
                end
                for j = 1:length(h6)
                    n = GridLine*(j-1)+1;
                    h6(j).XData = X3Distorteda(end-k-1+n,:,i);
                    h6(j).YData = X2Distorteda(end-k-1+n,:,i);
                    h6(j).ZData = X1Distorteda(end-k-1+n,:,i);
                end
                for j = 1:length(h1)
                    n = GridLine*(j-1)+1;
                    h1(j).XData = X3Distorteda(SamplesR+1:end-k-1,n,i);
                    h1(j).YData = X2Distorteda(SamplesR+1:end-k-1,n,i);
                    h1(j).ZData = X1Distorteda(SamplesR+1:end-k-1,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X3Distorteda(SamplesR+n,:,i);
                    h2(j).YData = X2Distorteda(SamplesR+n,:,i);
                    h2(j).ZData = X1Distorteda(SamplesR+n,:,i);
                end
                for j = 1:length(h7)
                    n = GridLine*(j-1)+1;
                    h7(j).XData = X3Distortedb(:,GridLine-1+n,i);
                    h7(j).YData = X2Distortedb(:,GridLine-1+n,i);
                    h7(j).ZData = X1Distortedb(:,GridLine-1+n,i);
                end
                for j = 1:length(h8)
                    n = GridLine*(j-1)+1;
                    h8(j).XData = X3Distortedb(n,:,i);
                    h8(j).YData = X2Distortedb(n,:,i);
                    h8(j).ZData = X1Distortedb(n,:,i);
                end
            elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
                for j = 1:length(h5)
                    n = GridLine*(j-1)+1;
                    h5(j).XData = X3Distorteda(end-k:end,n,i);
                    h5(j).YData = X2Distorteda(end-k:end,n,i);
                    h5(j).ZData = X1Distorteda(end-k:end,n,i);
                end
                for j = 1:length(h6)
                    n = GridLine*(j-1)+1;
                    h6(j).XData = X3Distorteda(end-k-1+n,:,i);
                    h6(j).YData = X2Distorteda(end-k-1+n,:,i);
                    h6(j).ZData = X1Distorteda(end-k-1+n,:,i);
                end
                for j = 1:length(h1)
                    n = GridLine*(j-1)+1;
                    h1(j).XData = X3Distorteda(1:end-k-1,n,i);
                    h1(j).YData = X2Distorteda(1:end-k-1,n,i);
                    h1(j).ZData = X1Distorteda(1:end-k-1,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X3Distorteda(n,:,i);
                    h2(j).YData = X2Distorteda(n,:,i);
                    h2(j).ZData = X1Distorteda(n,:,i);
                end
                for j = 1:length(h7)
                    n = GridLine*(j-1)+1;
                    h7(j).XData = X3Distortedb(:,GridLine-1+n,i);
                    h7(j).YData = X2Distortedb(:,GridLine-1+n,i);
                    h7(j).ZData = X1Distortedb(:,GridLine-1+n,i);
                end
                for j = 1:length(h8)
                    n = GridLine*(j-1)+1;
                    h8(j).XData = X3Distortedb(n,:,i);
                    h8(j).YData = X2Distortedb(n,:,i);
                    h8(j).ZData = X1Distortedb(n,:,i);
                end
            elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid || (ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
                for j = 1:length(h3)
                    n = GridLine*(j-1)+1;
                    h3(j).XData = X3Distorteda(1:SamplesR,n,i);
                    h3(j).YData = X2Distorteda(1:SamplesR,n,i);
                    h3(j).ZData = X1Distorteda(1:SamplesR,n,i);
                end
                for j = 1:length(h4)
                    n = GridLine*(j-1)+1;
                    h4(j).XData = X3Distorteda(n,:,i);
                    h4(j).YData = X2Distorteda(n,:,i);
                    h4(j).ZData = X1Distorteda(n,:,i);
                end
                for j = 1:length(h1)
                    n = GridLine*(j-1)+1;
                    h1(j).XData = X3Distorteda(SamplesR+1:end,n,i);
                    h1(j).YData = X2Distorteda(SamplesR+1:end,n,i);
                    h1(j).ZData = X1Distorteda(SamplesR+1:end,n,i);
                end
                for j = 1:length(h2)
                    n = GridLine*(j-1)+1;
                    h2(j).XData = X3Distorteda(SamplesR+n,:,i);
                    h2(j).YData = X2Distorteda(SamplesR+n,:,i);
                    h2(j).ZData = X1Distorteda(SamplesR+n,:,i);
                end
                for j = 1:length(h7)
                    n = GridLine*(j-1)+1;
                    h7(j).XData = X3Distortedb(:,GridLine-1+n,i);
                    h7(j).YData = X2Distortedb(:,GridLine-1+n,i);
                    h7(j).ZData = X1Distortedb(:,GridLine-1+n,i);
                end
                for j = 1:length(h8)
                    n = GridLine*(j-1)+1;
                    h8(j).XData = X3Distortedb(n,:,i);
                    h8(j).YData = X2Distortedb(n,:,i);
                    h8(j).ZData = X1Distortedb(n,:,i);
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