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
function ModeShapeImage_Isotropic(Crop,PNGresolution,Material,Viscoelastic,FluidLoading,Fluid,ALamb,AShear,AScholte,Directory,Export,FontSizeHeadLine,FontSizeAxesLabels,Frequency,GridLine,HeadLine,Length,LineWidth,Mode,FileName,PDF,PNG,Thickness,SamplesX1,SamplesX3,Scale,SLamb,SShear,SScholte,Undistorted,ShowHalfSpaces,HalfSpaces)
%#ok<*AGROW>
Lambda = conj(Material.Lambda_complex);
Mu = conj(Material.Mu_complex);
x3 = 0:Thickness/SamplesX3:Thickness;
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'SH')
    if  ~contains(Mode,'Scholte') && Mode(1) == 'S'
        q = find(SLamb{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SLamb{p}(:,1))) || Frequency < min(SLamb{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SLamb{p}(:,1))),' and ',num2str(ceil(max(SLamb{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(SLamb{p}(:,1)-Frequency) == min(abs(SLamb{p}(:,1)-Frequency)));
                q = q(1);
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
                q = find(abs(ALamb{p}(:,1)-Frequency) == min(abs(ALamb{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = ALamb{p}(q,1);
            end
        end
        PhaseVelocity = ALamb{p}(q,4)*1e3;
        Attenuation = ALamb{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    elseif contains(Mode,'Scholte') && Mode(1) == 'S' 
        q = find(SScholte{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SScholte{p}(:,1))) || Frequency < min(SScholte{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SScholte{p}(:,1))),' and ',num2str(ceil(max(SScholte{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(SScholte{p}(:,1)-Frequency) == min(abs(SScholte{p}(:,1)-Frequency)));
                q = q(1);
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
                q = find(abs(AScholte{p}(:,1)-Frequency) == min(abs(AScholte{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = AScholte{p}(q,1);
            end
        end
        PhaseVelocity = AScholte{p}(q,4)*1e3;
        Attenuation = AScholte{p}(q,6)*PhaseVelocity/(Frequency*1e3); % Np/wavelength
    end
    if  ~FluidLoading
        Fluid.Density = 1e-10;
        Fluid.Velocity = 1e-10;
    end
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
    k = AngularFrequency/PhaseVelocity*(1+1i*Attenuation/(2*pi));
    k3(1) = sqrt(AngularFrequency^2/Material.LongitudinalVelocity_complex^2-k^2);
    k3(2) = sqrt(AngularFrequency^2/Material.TransverseVelocity_complex^2-k^2);
    k3Fluid = sqrt(AngularFrequency^2/Fluid.Velocity^2-k^2);
    if  Viscoelastic && contains(Mode,'Scholte')
        k3Fluid = -k3Fluid;
    end
    W = (Material.Density*AngularFrequency^2-(Lambda+2*Mu)*k^2-Mu*k3.^2)./((Lambda+Mu)*k*k3); 
    WFluid = k3Fluid/k;
    D3 = 1i*(k*Lambda+(Lambda+2*Mu)*k3.*W); % sigma33
    D5 = 1i*Mu*(k3+k*W); % sigma13
    DFluid = 1i*Fluid.Density*AngularFrequency^2/k; % sigma11, sigma22, sigma33 in the fluid
    E = exp(1i*k3*Thickness);
    Z1 = [-W(2) W.*E -WFluid 0;-D3(2) -D3.*E DFluid 0;-D5(2) D5.*E 0 0;-W(2)*E(2) W 0 WFluid;-D3(2)*E(2) -D3 0 DFluid];
    Z2 = [W(1);D3(1);D5(1);W(1)*E(1);D3(1)*E(1)];
    U = Z1\Z2;
    for i = 1:length(x1)
        E = [exp(1i*(k*x1(i)+k3.*x3')) exp(1i*(k*x1(i)+k3.*(Thickness-x3)'))];
        u{i,1}(:,1) = E(:,1)+U(1)*E(:,2)+U(2)*E(:,3)+U(3)*E(:,4);
        u{i,1}(:,2) = W(1)*E(:,1)+W(2)*U(1)*E(:,2)-W(1)*U(2)*E(:,3)-W(2)*U(3)*E(:,4);
    end
    if  FluidLoading && ShowHalfSpaces
        x3Fluid = 0:Thickness/SamplesX3:HalfSpaces*Thickness;
        for i = 1:length(x1)
            EFluid = exp(1i*(k*x1(i)+k3Fluid*x3Fluid));
            uFluid0{i,1}(:,1) = U(4)*EFluid;
            uFluid0{i,1}(:,2) = -WFluid*U(4)*EFluid;
            uFluid1{i,1}(:,1) = U(5)*EFluid;
            uFluid1{i,1}(:,2) = WFluid*U(5)*EFluid;
        end
    end
elseif contains(Mode,'SH')
    if  Mode(1) == 'S'
        q = find(SShear{p}(:,1) == Frequency);
        if  isempty(q)
            if  Frequency > ceil(max(SShear{p}(:,1))) || Frequency < min(SShear{p}(:,1))
                errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(SShear{p}(:,1))),' and ',num2str(ceil(max(SShear{p}(:,1)))),' kHz.'],'Error');
                return
            else
                q = find(abs(SShear{p}(:,1)-Frequency) == min(abs(SShear{p}(:,1)-Frequency)));
                q = q(1);
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
            E = [exp(1i*(k*x1(i)+k3*x3')) exp(1i*(k*x1(i)+k3*(Thickness-x3)'))];
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
                q = find(abs(AShear{p-1}(:,1)-Frequency) == min(abs(AShear{p-1}(:,1)-Frequency)));
                q = q(1);
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
            E = [exp(1i*(k*x1(i)+k3*x3')) exp(1i*(k*x1(i)+k3*(Thickness-x3)'))];
            u{i,1}(:,1) = E(:,1)-E(:,2);
            u{i,1}(:,2) = 0;
        end
    end
    Shift = exp(-1i*angle(u{1,1}(2,1)));
    for i = 1:length(x1)
        u{i,1}(:,1) = u{i,1}(:,1)*Shift;
    end
    if  FluidLoading && ShowHalfSpaces
        x3Fluid = 0:Thickness/SamplesX3:HalfSpaces*Thickness;
        for i = 1:length(x1)
            uFluid0{i,1}(length(x3Fluid),2) = 0;
            uFluid1{i,1}(length(x3Fluid),2) = 0;
        end
    end
end
if  FluidLoading && ShowHalfSpaces
    x3 = horzcat(-fliplr(x3Fluid),x3,x3(end)+x3Fluid);
    for i = 1:length(x1)
        u{i,:} = vertcat(flipud(uFluid0{i,:}),u{i,:},uFluid1{i,:});
    end
end
if  ~contains(Mode,'SH')
    Shift = exp(-1i*angle(u{1,1}(2,1)));
    for i = 1:length(x1)
        u{i,1}(:,1) = u{i,1}(:,1)*Shift;
        u{i,1}(:,2) = u{i,1}(:,2)*Shift;
    end
end
[X1,X3] = meshgrid(x1,x3);
Ratio = (abs(x3(1))+abs(x3(end)))/x1(end);
u1Max = max(abs(real(u{1,1}(:,1))));
u3Max = max(abs(real(u{1,1}(:,2)))); % using real instead of imag is correct here; it is only a shift in phase between both
if  u1Max > u3Max
    for i = 1:size(X1,2)
        if  FluidLoading && ShowHalfSpaces
            X1Distorted(:,i) = X1(:,i)+Scale*x1(2)*SamplesX1/(80*(1+2*HalfSpaces))/u1Max*real(u{i}(:,1));
            X3Distorted(:,i) = X3(:,i)+Scale*x1(2)*SamplesX1/(80*(1+2*HalfSpaces))/u1Max*real(u{i}(:,2))*Ratio;
        else
            X1Distorted(:,i) = X1(:,i)+Scale*x1(2)*SamplesX1/80/u1Max*real(u{i}(:,1));
            X3Distorted(:,i) = X3(:,i)+Scale*x1(2)*SamplesX1/80/u1Max*real(u{i}(:,2))*Ratio;
        end
    end
else
    for i = 1:size(X1,2)
        if  FluidLoading && ShowHalfSpaces
            X1Distorted(:,i) = X1(:,i)+Scale*abs(x3(2)-x3(1))*SamplesX3/(40*(1+2*HalfSpaces))/u3Max*real(u{i}(:,1));
            X3Distorted(:,i) = X3(:,i)+Scale*abs(x3(2)-x3(1))*SamplesX3/(40*(1+2*HalfSpaces))/u3Max*real(u{i}(:,2))*Ratio;
        else
            X1Distorted(:,i) = X1(:,i)+Scale*abs(x3(2)-x3(1))*SamplesX3/40/u3Max*real(u{i}(:,1));
            X3Distorted(:,i) = X3(:,i)+Scale*abs(x3(2)-x3(1))*SamplesX3/40/u3Max*real(u{i}(:,2))*Ratio;
        end
    end
end
f = figure('Name','Mode shape','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
hold on
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
    line(X1Distorted(1:n+1,1:GridLine:end),X3Distorted(1:n+1,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
    line(X1Distorted(1:GridLine:n+1,:)',X3Distorted(1:GridLine:n+1,:)','LineWidth',LineWidth,'Color',[0 .5 1])
    line(X1Distorted(end-n:end,1:GridLine:end),X3Distorted(end-n:end,1:GridLine:end),'LineWidth',LineWidth,'Color',[0 .5 1])
    line(X1Distorted(end-n:GridLine:end,:)',X3Distorted(end-n:GridLine:end,:)','LineWidth',LineWidth,'Color',[0 .5 1])
    line(X1Distorted(n+2:end-n-1,1:GridLine:end),X3Distorted(n+2:end-n-1,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
    line(X1Distorted(n+2:GridLine:end-n-1,:)',X3Distorted(n+2:GridLine:end-n-1,:)','LineWidth',LineWidth,'Color','k')
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
    String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Thickness*1e3),'\,mm ',char(join(split(Material.Name,'_'),'\_'))];
    if  FluidLoading
        String = append(String,' in ',char(join(split(Fluid.Name,'_'),'\_')));
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
            if  ~Crop
                if  ~HeadLine
                    set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[43 24])
                else
                    set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[43 30])
                end
                print(f,fullfile(Directory,FileName),'-dpdf');
            else
                exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
            end
        end
        if  PNG
            if  ~Crop
                print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
            else
                exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
            end
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