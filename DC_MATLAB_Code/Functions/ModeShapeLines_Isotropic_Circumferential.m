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
function ModeShapeLines_Isotropic_Circumferential(FunctionMode,MakeFigure,ExportData,XSLX,TXT,MAT,Plot,PNGresolution,Material,Color1,Color2,Color3,Color4,Color5,Color6,CLamb,CShear,BoxLineWidth,Directory,FileName,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeLegend,Frequency,HeadLine,LegendLocation,LineWidth,Mode,PDF,PNG,Ro,Ri,SamplesR,Phase)
%#ok<*FNDSB>
r = (Ri:(Ro-Ri)/SamplesR:Ro)';
r2 = r.^2;
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
    d2f = d2JkLr+U(1)*d2YkLr;
    g = U(2)*JkTr+U(3)*YkTr; % SV_in + SV_out
    dg = U(2)*dJkTr+U(3)*dYkTr;
    d2g = U(2)*d2JkTr+U(3)*d2YkTr;
    u(:,2) = 1i*kRo*f./r-dg;
    u(:,3) = df+1i*kRo*g./r;
    v = -1i*AngularFrequency*u;
    epsilon(:,2) = (1i*kRo*u(:,2)+u(:,3))./r;
    epsilon(:,3) = d2f+1i*kRo*(dg./r-g./r2);
    epsilon(:,4) = 1i*kRo*(df./r-f./r2+u(:,3)./r)-d2g-u(:,2)./r;
    epsilon(:,6) = 0;
    sigma(:,1:3) = Material.Lambda*(d2f+df./r-kRo2*f./r2)+2*Material.Mu*epsilon(:,1:3);
    sigma(:,4) = Material.Mu*epsilon(:,4);
    sigma(:,6) = 0;
    StrainEnergyDensity = .5*real(epsilon(:,2).*conj(sigma(:,2))+epsilon(:,3).*conj(sigma(:,3))+epsilon(:,4).*conj(sigma(:,4)));
    KineticEnergyDensity = .5*Material.Density*(abs(v(:,2)).^2+abs(v(:,3)).^2);
    PowerFlowDensity(:,2) = -.5*real(sigma(:,2).*conj(v(:,2))+sigma(:,4).*conj(v(:,3)));
    if  Phase
        uPhase = rad2deg(angle(u*exp(-1i*angle(u(1,1)))));
        sigmaPhase = rad2deg(angle(sigma*exp(-1i*angle(u(1,1)))));
        epsilonPhase = rad2deg(angle(epsilon*exp(-1i*angle(u(1,1)))));
    end
    u(:,2) = u(:,2)*exp(-1i*angle(u(2,2)));
    u(:,3) = u(:,3)*exp(-1i*angle(u(2,3)));
    epsilon(:,2) = epsilon(:,2)*exp(-1i*angle(epsilon(2,2)));
    epsilon(:,3) = epsilon(:,3)*exp(-1i*angle(epsilon(2,3)));
    epsilon(:,4) = epsilon(:,4)*exp(-1i*angle(epsilon(2,4)));
    sigma(:,1) = sigma(:,1)*exp(-1i*angle(sigma(2,1)));
    sigma(:,2) = sigma(:,2)*exp(-1i*angle(sigma(2,2)));
    sigma(:,3) = sigma(:,3)*exp(-1i*angle(sigma(2,3)));
    sigma(:,4) = sigma(:,4)*exp(-1i*angle(sigma(2,4)));
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
    u(:,1) = J+U*Y;
    u(:,3) = 0;
    v = -1i*AngularFrequency*u;
    epsilon(:,5) = kT/2*(J_1-J1+U*(Y_1-Y1));
    epsilon(:,6) = 1i*kRo./r.*u(:,1);
    sigma(:,5:6) = Material.Mu*epsilon(:,5:6);
    StrainEnergyDensity = .5*real(epsilon(:,5).*conj(sigma(:,5))+epsilon(:,6).*conj(sigma(:,6)));
    KineticEnergyDensity = .5*Material.Density*abs(v(:,1)).^2;
    PowerFlowDensity(:,2) = -.5*real(sigma(:,6).*conj(v(:,1)));
    if  Phase
        uPhase = rad2deg(angle(u));
        sigmaPhase = rad2deg(angle(sigma));
        epsilonPhase = rad2deg(angle(epsilon));
    end
    u(:,1) = u(:,1)*exp(-1i*angle(u(2,1)));
    epsilon(:,5) = epsilon(:,5)*exp(-1i*angle(epsilon(2,5)));
    epsilon(:,6) = epsilon(:,6)*exp(-1i*angle(epsilon(2,6)));
    sigma(:,5) = sigma(:,5)*exp(-1i*angle(sigma(2,5)));
    sigma(:,6) = sigma(:,6)*exp(-1i*angle(sigma(2,6)));
end
PowerFlowDensity(:,3) = 0;
PowerFlow = trapz(r,PowerFlowDensity(:,2));
% disp(['ce: ',num2str(PowerFlow/(.5*trapz(r,StrainEnergyDensity+KineticEnergyDensity))),' m/s'])
u = u/sqrt(PowerFlow);
epsilon = epsilon/sqrt(PowerFlow);
sigma = sigma/sqrt(PowerFlow);
if  real(u(2,1)) < 0
    u = -u;
end
if  real(epsilon(2,1)) < 0
    epsilon = -epsilon;
end
if  real(sigma(2,1)) < 0
    sigma = -sigma;
end
if  Phase
    uPhase(find(round(uPhase) == -180)) = 180;
    sigmaPhase(find(round(sigmaPhase) == -180)) = 180;
    epsilonPhase(find(round(epsilonPhase) == -180)) = 180;
end
StrainEnergyDensity = .5*StrainEnergyDensity/PowerFlow;
KineticEnergyDensity = .5*KineticEnergyDensity/PowerFlow;
TotalEnergyDensity = StrainEnergyDensity+KineticEnergyDensity;
PowerFlowDensity = PowerFlowDensity/PowerFlow;
if  ExportData
    if  FunctionMode == 1
        if  ~Phase
            Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'d (mm)','uz (nm)','uq (nm)','ur (nm)'});
            Table(:,1) = num2cell((r-Ri)*1e3);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
        else
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'d (mm)','uz (nm)','uq (nm)','ur (nm)','Phase uz (deg)','Phase uq (deg)','Phase ur (deg)'});
            Table(:,1) = num2cell((r-Ri)*1e3);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
            Table(:,5) = num2cell(uPhase(:,1));
            Table(:,6) = num2cell(uPhase(:,2));
            Table(:,7) = num2cell(uPhase(:,3));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.txt']));
            end
            if  MAT
                Displacement = Table;
                save(fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.mat']),'Displacement')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end  
    elseif FunctionMode == 2
        if  ~Phase
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'d (mm)','sigmazz (kPa)','sigmaqq (kPa)','sigmarr (kPa)','sigmaqr (kPa)','sigmazr (kPa)','sigmazq (kPa)'});
            Table(:,1) = num2cell((r-Ri)*1e3);
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
        else
            Table = table('Size',[length(r) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'d (mm)','sigmazz (kPa)','sigmaqq (kPa)','sigmarr (kPa)','sigmaqr (kPa)','sigmazr (kPa)','sigmazq (kPa)','Phase sigmazz (deg)','Phase sigmaqq (deg)','Phase sigmarr (deg)','Phase sigmaqr (deg)','Phase sigmazr (deg)','Phase sigmazq (deg)'});
            Table(:,1) = num2cell((r-Ri)*1e3);
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
            Table(:,8) = num2cell(sigmaPhase(:,1));
            Table(:,9) = num2cell(sigmaPhase(:,2));
            Table(:,10) = num2cell(sigmaPhase(:,3));
            Table(:,11) = num2cell(sigmaPhase(:,4));
            Table(:,12) = num2cell(sigmaPhase(:,5));
            Table(:,13) = num2cell(sigmaPhase(:,6));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.txt']));
            end
            if  MAT
                Stress = Table;
                save(fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.mat']),'Stress')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    elseif FunctionMode == 3
        if  ~Phase
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'d (mm)','epsilonzz','epsilonqq','epsilonrr','epsilonqr','epsilonzr','epsilonzq'});
            Table(:,1) = num2cell((r-Ri)*1e3);
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
        else
            Table = table('Size',[length(r) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'d (mm)','epsilonzz','epsilonqq','epsilonrr','epsilonqr','epsilonzr','epsilonzq','Phase epsilonzz (deg)','Phase epsilonqq (deg)','Phase epsilonrr (deg)','Phase epsilonqr (deg)','Phase epsilonzr (deg)','Phase epsilonzq (deg)'});
            Table(:,1) = num2cell((r-Ri)*1e3);
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
            Table(:,8) = num2cell(epsilonPhase(:,1));
            Table(:,9) = num2cell(epsilonPhase(:,2));
            Table(:,10) = num2cell(epsilonPhase(:,3));
            Table(:,11) = num2cell(epsilonPhase(:,4));
            Table(:,12) = num2cell(epsilonPhase(:,5));
            Table(:,13) = num2cell(epsilonPhase(:,6));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.txt']));
            end
            if  MAT
                Strain = Table;
                save(fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.mat']),'Strain')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end   
    elseif FunctionMode == 4
        Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'d (mm)','Estrain (J/m2)','Ekin (J/m2)','Etotal (J/m2)'});
        Table(:,1) = num2cell((r-Ri)*1e3);
        Table(:,2) = num2cell(StrainEnergyDensity);
        Table(:,3) = num2cell(KineticEnergyDensity);
        Table(:,4) = num2cell(TotalEnergyDensity);
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.txt']));
            end
            if  MAT
                EnergyDensity = Table;
                save(fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.mat']),'EnergyDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end        
    elseif FunctionMode == 5
        Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'d (mm)','Pz (W/m)','Pq (W/m)','Pr (W/m)'});
        Table(:,1) = num2cell((r-Ri)*1e3);
        Table(:,2) = num2cell(PowerFlowDensity(:,1));
        Table(:,3) = num2cell(PowerFlowDensity(:,2));
        Table(:,4) = num2cell(PowerFlowDensity(:,3));
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.txt']));
            end
            if  MAT
                PowerFlowDensity = Table;
                save(fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*(Ro-Ri)),'MHzmm.mat']),'PowerFlowDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    end
end
if  MakeFigure
    if  Phase
        if  FunctionMode == 1 
            f = figure('Name','Displacement phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        elseif FunctionMode == 2
            f = figure('Name','Stress phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        elseif FunctionMode == 3
            f = figure('Name','Strain phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        end   
        datacursormode on
        x = xline(0,'Color',[.6 .6 .6]); 
        hasbehavior(x,'legend',false);
        hold on
        if  FunctionMode == 1
            if  Plot(3)
                plot(uPhase(:,3),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end
        end    
        ax = gca;
        ax.Box = 'on';
        ax.LineWidth = BoxLineWidth;
        ax.FontSize = FontSizeAxes;
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
        ax.XLabel.Interpreter = 'latex';
        ax.XLabel.FontSize = FontSizeAxesLabels;
        ax.XLim = max(abs(ax.XLim))*[-1 1];
        ax.YLabel.Interpreter = 'latex';
        ax.YLabel.FontSize = FontSizeAxesLabels;
        ax.YLabel.String = 'Thickness (mm)';
        ax.TickLabelInterpreter = 'latex';
        ax.YLim = [0 (r(end)-Ri)*1e3];
        if  FunctionMode == 1
            ax.XLabel.String = 'Displacement phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$u_r$','$u_z$','$u_\theta$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
            else
                LegendNames = {'Radial ($u_r$)','Axial ($u_z$)','Circumferential ($u_\theta$)'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')
            end
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        elseif FunctionMode == 2
            ax.XLabel.String = 'Stress phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$\sigma_{rr}$','$\sigma_{zz}$','$\sigma_{\theta\theta}$','$\sigma_{zr}$','$\sigma_{\theta r}$','$\sigma_{z\theta}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Radial ($\sigma_{rr}$)','Axial ($\sigma_{zz}$)','Circumferential ($\sigma_{\theta\theta}$)','Shear ($\sigma_{zr}$)','Shear ($\sigma_{\theta r}$)','Shear ($\sigma_{z\theta}$)'}; 
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
            end
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        elseif FunctionMode == 3
            ax.XLabel.String = 'Strain phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$\varepsilon_{rr}$','$\varepsilon_{zz}$','$\varepsilon_{\theta\theta}$','$\varepsilon_{zr}$','$\varepsilon_{\theta r}$','$\varepsilon_{z\theta}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Radial ($\varepsilon_{rr}$)','Axial ($\varepsilon_{zz}$)','Circumferential ($\varepsilon_{\theta\theta}$)','Shear ($\varepsilon_{zr}$)','Shear ($\varepsilon_{\theta r}$)','Shear ($\varepsilon_{z\theta}$)'}; 
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
            end 
            if  Export
                try
                    if  PDF
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.pdf']),'ContentType','vector')
                    end
                    if  PNG
                        exportgraphics(f,fullfile(Directory,[FileName,'_Phase.png']),'Resolution',PNGresolution)
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
                    return
                end
            end
        end
        tb = axtoolbar('default');
        tb.Visible = 'on';
        d = datacursormode(f);
        d.Interpreter = 'latex';
        d.UpdateFcn = @CursorPhase;
    end
    if  FunctionMode == 1
        f = figure('Name','Displacement','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 2
        f = figure('Name','Stress','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 3
        f = figure('Name','Strain','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 4
        f = figure('Name','Energy density','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif FunctionMode == 5 
        f = figure('Name','Power flow density','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    end
    datacursormode on
    x = xline(0,'Color',[.6 .6 .6]); 
    hasbehavior(x,'legend',false);
    hold on
    if  FunctionMode == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),(r-Ri)*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
    end    
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
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
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.XLim = max(abs(ax.XLim))*[-1 1];
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.String = 'Thickness (mm)';    
    ax.TickLabelInterpreter = 'latex';    
    ax.YLim = [0 (r(end)-Ri)*1e3];
    if  FunctionMode == 1
        ax.XLabel.String = 'Displacement (nm)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$u_r$','$u_z$','$u_\theta$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
        else
            LegendNames = {'Radial ($u_r$)','Axial ($u_z$)','Circumferential ($u_\theta$)'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')
        end
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
    elseif FunctionMode == 2
        ax.XLabel.String = 'Stress (kPa)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$\sigma_{rr}$','$\sigma_{zz}$','$\sigma_{\theta\theta}$','$\sigma_{zr}$','$\sigma_{\theta r}$','$\sigma_{z\theta}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Radial ($\sigma_{rr}$)','Axial ($\sigma_{zz}$)','Circumferential ($\sigma_{\theta\theta}$)','Shear ($\sigma_{zr}$)','Shear ($\sigma_{\theta r}$)','Shear ($\sigma_{z\theta}$)'}; 
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
        end
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
    elseif FunctionMode == 3
        ax.XLabel.String = 'Strain';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$\varepsilon_{rr}$','$\varepsilon_{zz}$','$\varepsilon_{\theta\theta}$','$\varepsilon_{zr}$','$\varepsilon_{\theta r}$','$\varepsilon_{z\theta}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Radial ($\varepsilon_{rr}$)','Axial ($\varepsilon_{zz}$)','Circumferential ($\varepsilon_{\theta\theta}$)','Shear ($\varepsilon_{zr}$)','Shear ($\varepsilon_{\theta r}$)','Shear ($\varepsilon_{z\theta}$)'}; 
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')                
        end 
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
    elseif FunctionMode == 4
        ax.XLabel.String = 'Energy density (J/m$^2$)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$E_\mathrm{strain}$','$E_\mathrm{kin}$','$E_\mathrm{total}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Strain energy density','Kinetic energy density','Total energy density'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')        
        end       
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
    elseif FunctionMode == 5
        ax.XLabel.String = 'Power flow density (W/m)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$p_r$','$p_z$','$p_\theta$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Radial ($p_r$)','Axial ($p_z$)','Circumferential ($p_\theta$)'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeastoutside','FontSize',FontSizeLegend,'Interpreter','latex')        
        end
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
    end
    tb = axtoolbar('default');
    tb.Visible = 'on';
    d = datacursormode(f);
    d.Interpreter = 'latex';
    d.UpdateFcn = @Cursor;
end
function output_txt = Cursor(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = replace(InnerFluid.Name,'_','\_');
        else
            output_txt = replace(OuterFluid.Name,'_','\_');
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$u_z$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$u_\theta$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$u_r$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\sigma_{rr}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\sigma_{zr}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\sigma_{\theta r}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\sigma_{zz}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\sigma_{\theta\theta}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\sigma_{z\theta}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varepsilon_{rr}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varepsilon_{zr}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varepsilon_{\theta r}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varepsilon_{zz}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varepsilon_{\theta\theta}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varepsilon_{z\theta}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 4
            if  event_obj.Target.Color == Color1
                output_txt = {['$E_\mathrm{strain}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$E_\mathrm{kin}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$E_\mathrm{total}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 5 
            if  event_obj.Target.Color == Color1
                output_txt = {['$p_z$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$p_\theta$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$p_r$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        end
    end
end
function output_txt = CursorPhase(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = replace(InnerFluid.Name,'_','\_');
        else
            output_txt = replace(OuterFluid.Name,'_','\_');
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$\varphi(u_z)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(u_\theta)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$\varphi(u_r)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\sigma_{rr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\sigma_{zr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\sigma_{\theta r})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\sigma_{zz})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\sigma_{\theta\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\sigma_{z\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\varepsilon_{rr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\varepsilon_{zr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\varepsilon_{\theta r})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\varepsilon_{zz})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\varepsilon_{\theta\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\varepsilon_{z\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$d$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        end
    end
end
end