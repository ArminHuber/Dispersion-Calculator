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
function ModeShapeGrid_Isotropic_Circumferential(PNGresolution,Material,CLamb,CShear,Directory,Export,FontSizeHeadLine,FontSizeAxesLabels,Frequency,GridLine,HeadLine,Length,LineWidth,Mode,FileName,PDF,PNG,Ro,Ri,SamplesX1,SamplesR,Gain,Undistorted)
%#ok<*AGROW>
r = (Ri:(Ro-Ri)/SamplesR:Ro)';
Ro2 = Ro^2;
Ri2 = Ri^2;
p = str2double(regexp(Mode,'\d*','match'))+1;
if  ~contains(Mode,'SH')
    q = find(CLamb{p}(:,1) == Frequency);
    if  isempty(q)
        if  Frequency > ceil(max(CLamb{p}(:,1))) || Frequency < min(CLamb{p}(:,1))
            errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(CLamb{p}(:,1))),' and ',num2str(ceil(max(CLamb{p}(:,1)))),' kHz.'],'Error');
            return
        else
            [~,q] = min(abs(CLamb{p}(:,1)-Frequency));
            Frequency = CLamb{p}(q,1);
        end
    end
    PhaseVelocity = CLamb{p}(q,4)*1e3;
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
    k = AngularFrequency/PhaseVelocity;
    kL = AngularFrequency/Material.LongitudinalVelocity;
    kT = AngularFrequency/Material.TransverseVelocity;
    kL2 = kL^2;
    kT2 = kT^2;
    kLr = kL*r;
    kTr = kT*r;
    kRo = k*Ro;
    kRo2 = kRo^2;
    kRo1 = 1/(kRo+1);
    kRo_1 = 1/(kRo-1);
    JkLr = besselj(kRo,kLr);
    JkTr = besselj(kRo,kTr);
    YkLr = bessely(kRo,kLr);
    YkTr = bessely(kRo,kTr);
    J_2kLr = besselj(kRo-2,kLr);
    J_2kTr = besselj(kRo-2,kTr);
    Y_2kLr = bessely(kRo-2,kLr);
    Y_2kTr = bessely(kRo-2,kTr);
    J2kLr = besselj(kRo+2,kLr);
    J2kTr = besselj(kRo+2,kTr);
    Y2kLr = bessely(kRo+2,kLr);
    Y2kTr = bessely(kRo+2,kTr);
    dJkLr = kL2*r/4.*(kRo_1*(J_2kLr+JkLr)-kRo1*(JkLr+J2kLr));
    dJkTr = kT2*r/4.*(kRo_1*(J_2kTr+JkTr)-kRo1*(JkTr+J2kTr));
    dYkLr = kL2*r/4.*(kRo_1*(Y_2kLr+YkLr)-kRo1*(YkLr+Y2kLr));
    dYkTr = kT2*r/4.*(kRo_1*(Y_2kTr+YkTr)-kRo1*(YkTr+Y2kTr));
    d2JkLr = kL2/4*(J_2kLr-2*JkLr+J2kLr);
    d2JkTr = kT2/4*(J_2kTr-2*JkTr+J2kTr);
    d2YkLr = kL2/4*(Y_2kLr-2*YkLr+Y2kLr);
    d2YkTr = kT2/4*(Y_2kTr-2*YkTr+Y2kTr);
    Z1(1,1) = 2i*Material.Mu*kRo*(dYkLr(1)/Ri-YkLr(1)/Ri2);
    Z1(1,2) = -Material.Mu*(d2JkTr(1)-dJkTr(1)/Ri+kRo2*JkTr(1)/Ri2);
    Z1(1,3) = -Material.Mu*(d2YkTr(1)-dYkTr(1)/Ri+kRo2*YkTr(1)/Ri2);
    Z1(2,1) = (Material.Lambda+2*Material.Mu)*d2YkLr(end)+Material.Lambda*(dYkLr(end)/Ro-kRo2*YkLr(end)/Ro2);
    Z1(2,2) = 2i*Material.Mu*kRo*(dJkTr(end)/Ro-JkTr(end)/Ro2);
    Z1(2,3) = 2i*Material.Mu*kRo*(dYkTr(end)/Ro-YkTr(end)/Ro2);
    Z1(3,1) = 2i*Material.Mu*kRo*(dYkLr(end)/Ro-YkLr(end)/Ro2);
    Z1(3,2) = -Material.Mu*(d2JkTr(end)-dJkTr(end)/Ro+kRo2*JkTr(end)/Ro2);
    Z1(3,3) = -Material.Mu*(d2YkTr(end)-dYkTr(end)/Ro+kRo2*YkTr(end)/Ro2);
    Z2(1,1) = 2i*Material.Mu*kRo*(dJkLr(1)/Ri-JkLr(1)/Ri2);
    Z2(2,1) = (Material.Lambda+2*Material.Mu)*d2JkLr(end)+Material.Lambda*(dJkLr(end)/Ro-kRo2*JkLr(end)/Ro2);
    Z2(3,1) = 2i*Material.Mu*kRo*(dJkLr(end)/Ro-JkLr(end)/Ro2);
    U = Z1\-Z2;
    f = JkLr+U(1)*YkLr; % L_in + L_out
    df = dJkLr+U(1)*dYkLr;
    g = U(2)*JkTr+U(3)*YkTr; % SV_in + SV_out
    dg = U(2)*dJkTr+U(3)*dYkTr;
    u2 = 1i*kRo*f./r-dg; % uq
    u3 = df+1i*kRo*g./r; % ur
    E = exp(1i*k*x1);
    for i = 1:length(x1)
        u{i,1}(:,1) = u2*E(i);
        u{i,1}(:,2) = u3*E(i);
    end
elseif contains(Mode,'SH')
    q = find(CShear{p}(:,1) == Frequency);
    if  isempty(q)
        if  Frequency > ceil(max(CShear{p}(:,1))) || Frequency < min(CShear{p}(:,1))
            errordlg(['Selected frequency outside frequency range! Select between ',num2str(min(CShear{p}(:,1))),' and ',num2str(ceil(max(CShear{p}(:,1)))),' kHz.'],'Error');
            return
        else
            [~,q] = min(abs(CShear{p}(:,1)-Frequency));
            Frequency = CShear{p}(q,1);
        end
    end
    PhaseVelocity = CShear{p}(q,4)*1e3;
    x1 = 0:Length*PhaseVelocity/(1e3*SamplesX1*Frequency):Length*PhaseVelocity/(1e3*Frequency);
    AngularFrequency = 2*pi*Frequency*1e3;
    k = AngularFrequency/PhaseVelocity;
    kT = AngularFrequency/Material.TransverseVelocity;
    kTr = kT*r;
    kRo = k*Ro;
    J = besselj(kRo,kTr);
    Y = bessely(kRo,kTr);
    J_1 = besselj(kRo-1,kTr);
    Y_1 = bessely(kRo-1,kTr);
    J1 = besselj(kRo+1,kTr);
    Y1 = bessely(kRo+1,kTr);
    U = (Y_1(end)-Y1(end))\-(J_1(end)-J1(end));
    u1 = J+U*Y; % uz
    E = exp(1i*k*x1);
    for i = 1:length(x1)
        u{i,1}(:,1) = u1*E(i);
        u{i,1}(:,2) = 0;
    end
end
[X1,X3] = meshgrid(x1,r);
Ratio = (r(end)-r(1))/x1(end);
u1Max = max(abs(real(u{1}(:,1))));
u3Max = max(abs(real(u{1}(:,2))));
if  u1Max > u3Max
    Compensation = Gain*x1(end)/40/u1Max;
else
    Compensation = Gain*(r(end)-r(1))/20/u3Max;
end
for i = 1:size(X1,2)
    X1Distorted(:,i) = X1(:,i)+real(u{i}(:,1))*Compensation;
    X3Distorted(:,i) = X3(:,i)+real(u{i}(:,2))*Compensation*Ratio;
end
f = figure('Name','Mode shape','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
if  Undistorted
    line(X1(:,1:GridLine:end),X3(:,1:GridLine:end),'LineWidth',LineWidth,'Color',[1 .7 0])
    line(X1(1:GridLine:end,:)',X3(1:GridLine:end,:)','LineWidth',LineWidth,'Color',[1 .7 0])
end
line(X1Distorted(:,1:GridLine:end),X3Distorted(:,1:GridLine:end),'LineWidth',LineWidth,'Color','k')
line(X1Distorted(1:GridLine:end,:)',X3Distorted(1:GridLine:end,:)','LineWidth',LineWidth,'Color','k')
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
    String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,$\times$\,',num2str((Ro-Ri)*1e3),'\,mm ',replace(Material.Name,'_','\_'),' circumference'];
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