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
function ModeShapeGrid_Isotropic(PNGresolution,Material,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,ALamb,AShear,AScholte,BLamb,BScholte,Directory,Export,FontSizeHeadLine,FontSizeAxesLabels,Frequency,GridLine,HeadLine,Length,LineWidth,Mode,FileName,PDF,PNG,Thickness,SamplesX1,SamplesX3,Gain,SLamb,SShear,SScholte,Undistorted,ShowHalfSpaces,HalfSpaces)
%#ok<*AGROW>
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
x3 = (0:Thickness/SamplesX3:Thickness)'; % (m)
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'SH')
    if  ~contains(Mode,'Scholte') && Mode(1) == 'S'
        q = find(SLamb{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SLamb{p}(:,1))) || Frequency < min(SLamb{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SLamb{p}(:,1))),' and ',num2str(ceil(max(SLamb{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(SLamb{p}(:,1)-Frequency));
                Frequency = SLamb{p}(q,1);
            end
        end
        PhaseVelocity = SLamb{p}(q,4)*1e3;
        Attenuation = SLamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif ~contains(Mode,'Scholte') && Mode(1) == 'A'
        q = find(ALamb{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(ALamb{p}(:,1))) || Frequency < min(ALamb{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(ALamb{p}(:,1))),' and ',num2str(ceil(max(ALamb{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(ALamb{p}(:,1)-Frequency));
                Frequency = ALamb{p}(q,1);
            end
        end
        PhaseVelocity = ALamb{p}(q,4)*1e3;
        Attenuation = ALamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif ~contains(Mode,'Scholte') && Mode(1) == 'B'
        q = find(BLamb{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(BLamb{p}(:,1))) || Frequency < min(BLamb{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BLamb{p}(:,1))),' and ',num2str(ceil(max(BLamb{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(BLamb{p}(:,1)-Frequency));
                Frequency = BLamb{p}(q,1);
            end
        end
        PhaseVelocity = BLamb{p}(q,4)*1e3;
        Attenuation = BLamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'S' 
        q = find(SScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SScholte{p}(:,1))) || Frequency < min(SScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SScholte{p}(:,1))),' and ',num2str(ceil(max(SScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(SScholte{p}(:,1)-Frequency));
                Frequency = SScholte{p}(q,1);
            end
        end
        PhaseVelocity = SScholte{p}(q,4)*1e3;
        Attenuation = SScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'A'
        q = find(AScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(AScholte{p}(:,1))) || Frequency < min(AScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AScholte{p}(:,1))),' and ',num2str(ceil(max(AScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(AScholte{p}(:,1)-Frequency));
                Frequency = AScholte{p}(q,1);
            end
        end
        PhaseVelocity = AScholte{p}(q,4)*1e3;
        Attenuation = AScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'B'
        q = find(BScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(BScholte{p}(:,1))) || Frequency < min(BScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(BScholte{p}(:,1))),' and ',num2str(ceil(max(BScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(BScholte{p}(:,1)-Frequency));
                Frequency = BScholte{p}(q,1);
            end
        end
        PhaseVelocity = BScholte{p}(q,4)*1e3;
        Attenuation = BScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    end
    if  ~ToggleUpperFluid
        UpperFluid.Velocity = 1e-10;
        UpperFluid.Density = 1e-10;
    end
    if  ~ToggleLowerFluid
        LowerFluid.Velocity = 1e-10;
        LowerFluid.Density = 1e-10;
    end
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
    AngularFrequency2 = AngularFrequency^2;
    k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
    k2 = k^2;
    k3(1) = sqrt(AngularFrequency2/Material.LongitudinalVelocity_complex^2-k2);
    k3(2) = sqrt(AngularFrequency2/Material.TransverseVelocity_complex^2-k2);
    k3UpperFluid = sqrt(AngularFrequency2/UpperFluid.Velocity^2-k2);
    k3LowerFluid = sqrt(AngularFrequency2/LowerFluid.Velocity^2-k2);
    if  contains(Mode,'Scholte') && Attenuation ~= 0
        if  PhaseVelocity < UpperFluid.Velocity
            k3UpperFluid = -k3UpperFluid;
        end
        if  PhaseVelocity < LowerFluid.Velocity
            k3LowerFluid = -k3LowerFluid;
        end
    end
    W = (Material.Density*AngularFrequency2-(Lambda+2*Mu)*k2-Mu*k3.^2)./((Lambda+Mu)*k*k3); 
    WUpperFluid = k3UpperFluid/k;
    WLowerFluid = k3LowerFluid/k;
    D3 = 1i*(k*Lambda+(Lambda+2*Mu)*k3.*W); % sigma33
    D5 = 1i*Mu*(k3+k*W); % sigma13
    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency2/k; % sigma11, sigma22, sigma33 in the upper fluid
    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency2/k; % in the lower fluid
    E = exp(1i*k3*Thickness);
    Z1 = [-W(2) W.*E -WUpperFluid 0;-D3(2) -D3.*E DUpperFluid 0;-D5(2) D5.*E 0 0;-W(2)*E(2) W 0 WLowerFluid;-D3(2)*E(2) -D3 0 DLowerFluid];
    Z2 = [W(1);D3(1);D5(1);W(1)*E(1);D3(1)*E(1)];
    U = Z1\Z2;
    for i = 1:length(x1)
        E = [exp(1i*(k*x1(i)+k3.*x3)) exp(1i*(k*x1(i)+k3.*(Thickness-x3)))];
        u{i,1}(:,1) = E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4);
        u{i,1}(:,2) = W(1)*E(:,1)+W(2)*U(1)*E(:,2)-W(1)*U(2)*E(:,3)-W(2)*U(3)*E(:,4);
    end
elseif contains(Mode,'SH')
    if  Mode(1) == 'S'
        q = find(SShear{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SShear{p}(:,1))) || Frequency < min(SShear{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SShear{p}(:,1))),' and ',num2str(ceil(max(SShear{p}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(SShear{p}(:,1)-Frequency));
                Frequency = SShear{p}(q,1);
            end
        end
        PhaseVelocity = SShear{p}(q,4)*1e3;
        Attenuation = SShear{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
        for i = 1:length(x1)
            E = [exp(1i*(k*x1(i)+k3*x3)) exp(1i*(k*x1(i)+k3*(Thickness-x3)))];
            u{i,1}(:,1) = E(:,1)+E(:,2);
            u{i,1}(:,2) = 0;
        end
    elseif Mode(1) == 'A'
        q = find(AShear{p-1}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(AShear{p-1}(:,1))) || Frequency < min(AShear{p-1}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(AShear{p-1}(:,1))),' and ',num2str(ceil(max(AShear{p-1}(:,1)))),' kHz.'],'Error');
                return
            else
                [~,q] = min(abs(AShear{p-1}(:,1)-Frequency));
                Frequency = AShear{p-1}(q,1);
            end
        end
        PhaseVelocity = AShear{p-1}(q,4)*1e3;
        Attenuation = AShear{p-1}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
        AngularFrequency = 2*pi*Frequency*1e3;
        k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
        k3 = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
        x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
        for i = 1:length(x1)
            E = [exp(1i*(k*x1(i)+k3*x3)) exp(1i*(k*x1(i)+k3*(Thickness-x3)))];
            u{i,1}(:,1) = E(:,1)-E(:,2);
            u{i,1}(:,2) = 0;
        end
    end
end
if  FluidLoading && ShowHalfSpaces
    x3Fluid = (0:Thickness/SamplesX3:HalfSpaces*Thickness)';
    if  ToggleUpperFluid
        if  ~contains(Mode,'SH')
            for i = 1:length(x1)
                EUpperFluid = exp(1i*(k*x1(i)+k3UpperFluid*x3Fluid));
                uUpperFluid{i,1}(:,1) = U(4)*EUpperFluid;
                uUpperFluid{i,1}(:,2) = -WUpperFluid*U(4)*EUpperFluid;
            end
        elseif contains(Mode,'SH')
            for i = 1:length(x1)
                uUpperFluid{i,1}(length(x3Fluid),2) = 0;
            end
        end
        x3 = vertcat(-flipud(x3Fluid),x3);
        for i = 1:length(x1)
            u{i,:} = vertcat(flipud(uUpperFluid{i,:}),u{i,:});
        end
    end
    if  ToggleLowerFluid
        if  ~contains(Mode,'SH')
            for i = 1:length(x1)
                ELowerFluid = exp(1i*(k*x1(i)+k3LowerFluid*x3Fluid));
                uLowerFluid{i,1}(:,1) = U(5)*ELowerFluid;
                uLowerFluid{i,1}(:,2) = WLowerFluid*U(5)*ELowerFluid;
            end
        elseif contains(Mode,'SH')
            for i = 1:length(x1)
                uLowerFluid{i,1}(length(x3Fluid),2) = 0;
            end
        end
        x3 = vertcat(x3,x3(end)+x3Fluid);
        for i = 1:length(x1)
            u{i,:} = vertcat(u{i,:},uLowerFluid{i,:});
        end
    end
end
[X1,X3] = meshgrid(x1,x3);
Ratio = (abs(x3(1))+abs(x3(end)))/x1(end);
u1Max = max(abs(real(u{1}(:,1))));
u3Max = max(abs(real(u{1}(:,2))));
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
for i = 1:size(X1,2)
    X1Distorted(:,i) = X1(:,i)+real(u{i}(:,1))*Compensation;
    X3Distorted(:,i) = X3(:,i)+real(u{i}(:,2))*Compensation*Ratio;
end
f = figure('Name','Mode shape','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
if  ~FluidLoading || ~ShowHalfSpaces
    if  Undistorted
        line(X1(:,1:GridLine:end),X3(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(1:GridLine:end,:)',X3(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
    end
    line(X1Distorted(:,1:GridLine:end),X3Distorted(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
    line(X1Distorted(1:GridLine:end,:)',X3Distorted(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
elseif FluidLoading && ShowHalfSpaces
    n = HalfSpaces*SamplesX3;
    if  Undistorted
        line(X1(1:n+1,1:GridLine:end),X3(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(1:GridLine:n+1,:)',X3(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(end-n:end,1:GridLine:end),X3(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(end-n:GridLine:end,:)',X3(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(n+2:end-n-1,1:GridLine:end),X3(n+2:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
        line(X1(n+2:GridLine:end-n-1,:)',X3(n+2:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color',[1 .7 0])
    end
    if  ToggleUpperFluid && ToggleLowerFluid
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
        line(X1Distorted(1:n+1,1:GridLine:end),X3Distorted(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',UpperFluidColor)
        line(X1Distorted(1:GridLine:n+1,:)',X3Distorted(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',UpperFluidColor)
        line(X1Distorted(end-n:end,1:GridLine:end),X3Distorted(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',LowerFluidColor)
        line(X1Distorted(end-n:GridLine:end,:)',X3Distorted(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',LowerFluidColor)
        line(X1Distorted(n+2:end-n-1,1:GridLine:end),X3Distorted(n+2:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(n+2:GridLine:end-n-1,:)',X3Distorted(n+2:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k')
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        line(X1Distorted(1:n+1,1:GridLine:end),X3Distorted(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(1:GridLine:n+1,:)',X3Distorted(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(n+2:end,1:GridLine:end),X3Distorted(n+2:end,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(n+2:GridLine:end,:)',X3Distorted(n+2:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        line(X1Distorted(end-n:end,1:GridLine:end),X3Distorted(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(end-n:GridLine:end,:)',X3Distorted(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1])
        line(X1Distorted(1:end-n-1,1:GridLine:end),X3Distorted(1:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
        line(X1Distorted(1:GridLine:end-n-1,:)',X3Distorted(1:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k')
    end
end
ax = gca;
axis off
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
if  HeadLine
    if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
        ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
    elseif contains(Mode,'SH')
        ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
    elseif contains(Mode,'Scholte')
        ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
    end
    String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Thickness*1e3),'\,mm ',replace(Material.Name,'_','\_'),' plate'];
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
function output_txt = Cursor(~,~)
    output_txt = {[]};
end
end