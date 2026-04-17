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
function ModeShapeLines_Isotropic_Rod_Pipe(Quantity,FunctionMode,DataFormat,Geometry,Plot,PNGresolution,Material,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Color1,Color2,Color3,Color4,Color5,Color6,F,L,T,BoxLineWidth,Directory,FileName,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeLegend,Frequency,HeadLine,LegendLocation,LineWidth,Mode,PDF,PNG,Ro,Ri,SamplesR,ShowHalfSpace,HalfSpaces,Phase)
[u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,r,n,p,Frequency] = ModeShapeLinesComputer_Isotropic_Rod_Pipe(Geometry,Material,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,F,L,T,Frequency,Mode,Ro,Ri,SamplesR,ShowHalfSpace,HalfSpaces,Phase);
if  FunctionMode == 2
    if  strcmp(Geometry,'Rod')
        Product = num2str(Frequency*2*Ro);
    else
        Product = num2str(Frequency*(Ro-Ri));
    end
    if  Quantity == 1
        if  ~Phase
            Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'r (mm)','uz (nm)','uq (nm)','ur (nm)'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
        else
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','uz (nm)','uq (nm)','ur (nm)','Phase uz (deg)','Phase uq (deg)','Phase ur (deg)'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
            Table(:,5) = num2cell(uPhase(:,1));
            Table(:,6) = num2cell(uPhase(:,2));
            Table(:,7) = num2cell(uPhase(:,3));
        end
        try
            if  DataFormat == 1 % mat
                save(fullfile(Directory,[FileName,'_Displacement_',Mode,'@',Product,'MHzmm.mat']),'Table')
            elseif DataFormat == 2 % xlsx
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',Product,'MHzmm.xlsx']));
            elseif DataFormat == 3 % txt
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',Product,'MHzmm.txt']));
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end  
    elseif Quantity == 2
        if  ~Phase
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','sigmazz (kPa)','sigmaqq (kPa)','sigmarr (kPa)','sigmaqr (kPa)','sigmazr (kPa)','sigmazq (kPa)'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
        else
            Table = table('Size',[length(r) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','sigmazz (kPa)','sigmaqq (kPa)','sigmarr (kPa)','sigmaqr (kPa)','sigmazr (kPa)','sigmazq (kPa)','Phase sigmazz (deg)','Phase sigmaqq (deg)','Phase sigmarr (deg)','Phase sigmaqr (deg)','Phase sigmazr (deg)','Phase sigmazq (deg)'});
            Table(:,1) = num2cell(1e3*r);
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
            if  DataFormat == 1 % mat
                save(fullfile(Directory,[FileName,'_Stress_',Mode,'@',Product,'MHzmm.mat']),'Table')
            elseif DataFormat == 2 % xlsx
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',Product,'MHzmm.xlsx']));
            elseif DataFormat == 3 % txt
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',Product,'MHzmm.txt']));
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end   
    elseif Quantity == 3
        if  ~Phase
            Table = table('Size',[length(r) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','epsilonzz','epsilonqq','epsilonrr','epsilonqr','epsilonzr','epsilonzq'});
            Table(:,1) = num2cell(1e3*r);
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
        else
            Table = table('Size',[length(r) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'r (mm)','epsilonzz','epsilonqq','epsilonrr','epsilonqr','epsilonzr','epsilonzq','Phase epsilonzz (deg)','Phase epsilonqq (deg)','Phase epsilonrr (deg)','Phase epsilonqr (deg)','Phase epsilonzr (deg)','Phase epsilonzq (deg)'});
            Table(:,1) = num2cell(1e3*r);
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
            if  DataFormat == 1 % mat
                save(fullfile(Directory,[FileName,'_Strain_',Mode,'@',Product,'MHzmm.mat']),'Table')
            elseif DataFormat == 2 % xlsx
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',Product,'MHzmm.xlsx']));
            elseif DataFormat == 3 % txt
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',Product,'MHzmm.txt']));
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end   
    elseif Quantity == 4
        Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'r (mm)','Estrain (J/m2)','Ekin (J/m2)','Etotal (J/m2)'});
        Table(:,1) = num2cell(1e3*r);
        Table(:,2) = num2cell(StrainEnergyDensity);
        Table(:,3) = num2cell(KineticEnergyDensity);
        Table(:,4) = num2cell(TotalEnergyDensity);
        try
            if  DataFormat == 1 % mat
                save(fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',Product,'MHzmm.mat']),'Table')
            elseif DataFormat == 2 % xlsx
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',Product,'MHzmm.xlsx']));
            elseif DataFormat == 3 % txt
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',Product,'MHzmm.txt']));
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    elseif Quantity == 5
        Table = table('Size',[length(r) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'r (mm)','Pz (W/m)','Pq (W/m)','Pr (W/m)'});
        Table(:,1) = num2cell(1e3*r);
        Table(:,2) = num2cell(PowerFlowDensity(:,1));
        Table(:,3) = num2cell(PowerFlowDensity(:,2));
        Table(:,4) = num2cell(PowerFlowDensity(:,3));
        try
            if  DataFormat == 1 % mat
                save(fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',Product,'MHzmm.mat']),'Table')
            elseif DataFormat == 2 % xlsx
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',Product,'MHzmm.xlsx']));
            elseif DataFormat == 3 % txt
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',Product,'MHzmm.txt']));
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    end
elseif FunctionMode == 1
    if  Phase
        if  Quantity == 1 
            f = figure('Icon',which('DC_Logo16.png'),'Name','Displacement phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        elseif Quantity == 2
            f = figure('Icon',which('DC_Logo16.png'),'Name','Stress phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        elseif Quantity == 3
            f = figure('Icon',which('DC_Logo16.png'),'Name','Strain phase','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
        end   
        x = xline(0,'Color',[.6 .6 .6],'PickableParts','none'); 
        hasbehavior(x,'legend',false);
        hold on
        if  Quantity == 1
            if  Plot(3)
                plot(uPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif Quantity == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),r*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),r*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),r*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif Quantity == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),r*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),r*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),r*1e3,'LineWidth',LineWidth,'Color',Color6);
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
            if  ~contains(Mode,'Scholte')
                ModeName = [Mode(1),'(',num2str(n),',',num2str(p(2)),')'];
            else
                ModeName = [Mode(1),'(',num2str(n),',',num2str(p(2)),')$^\mathrm{Scholte}$'];
            end
            if  strcmp(Geometry,'Rod')
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
            elseif strcmp(Geometry,'Pipe')
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,$\times$\,',num2str((Ro-Ri)*1e3),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
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
        ax.XLabel.Interpreter = 'latex';
        ax.XLabel.FontSize = FontSizeAxesLabels;
        ax.XLim = max(abs(ax.XLim))*[-1 1];
        ax.YLabel.Interpreter = 'latex';
        ax.YLabel.FontSize = FontSizeAxesLabels;
        ax.YLabel.String = '$r$ (mm)';
        ax.TickLabelInterpreter = 'latex';
        if  strcmp(Geometry,'Rod') || (strcmp(Geometry,'pipe') && ToggleInnerFluid)
            ax.YLim = [0 r(end)*1e3];
        else
            ax.YLim = [r(1) r(end)]*1e3;
        end
        if  strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
            if  strcmp(OuterFluid.Name,InnerFluid.Name)
                OuterFaceAlpha = .2;
                InnerFaceAlpha = .2;
            else
                if  OuterFluid.Density*OuterFluid.Velocity > InnerFluid.Density*InnerFluid.Velocity
                    OuterFaceAlpha = .2;
                    InnerFaceAlpha = .1;
                else
                    OuterFaceAlpha = .1;
                    InnerFaceAlpha = .2;
                end
            end
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[Ro*1e3 Ro*1e3 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',OuterFaceAlpha,'EdgeColor','none','PickableParts','none')
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 Ri*1e3 Ri*1e3],'b','FaceAlpha',InnerFaceAlpha,'EdgeColor','none','PickableParts','none')
        elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[Ro*1e3 Ro*1e3 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none','PickableParts','none')
        elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 Ri*1e3 Ri*1e3],'b','FaceAlpha',.2,'EdgeColor','none','PickableParts','none')
        end
        if  Quantity == 1
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
        elseif Quantity == 2
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
        elseif Quantity == 3
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
        datacursormode on
        d = datacursormode(f);
        d.Interpreter = 'latex';
        d.UpdateFcn = @CursorPhase;
    end
    if  Quantity == 1
        f = figure('Icon',which('DC_Logo16.png'),'Name','Displacement','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif Quantity == 2
        f = figure('Icon',which('DC_Logo16.png'),'Name','Stress','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif Quantity == 3
        f = figure('Icon',which('DC_Logo16.png'),'Name','Strain','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif Quantity == 4
        f = figure('Icon',which('DC_Logo16.png'),'Name','Energy density','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    elseif Quantity == 5 
        f = figure('Icon',which('DC_Logo16.png'),'Name','Power flow density','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    end
    x = xline(0,'Color',[.6 .6 .6],'PickableParts','none'); 
    hasbehavior(x,'legend',false);
    hold on
    if  Quantity == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end
    elseif Quantity == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,r*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif Quantity == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),r*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),r*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),r*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif Quantity == 4
        if  Plot(5)
            plot(StrainEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif Quantity == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),r*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),r*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),r*1e3,'LineWidth',LineWidth,'Color',Color2);
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
        if  ~contains(Mode,'Scholte')
            ModeName = [Mode(1),'(',num2str(n),',',num2str(p(2)),')'];
        else
            ModeName = [Mode(1),'(',num2str(n),',',num2str(p(2)),')$^\mathrm{Scholte}$'];
        end
        if  strcmp(Geometry,'Rod')
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
        elseif strcmp(Geometry,'Pipe')
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Ro*2e3),'\,$\times$\,',num2str((Ro-Ri)*1e3),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
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
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.XLim = max(abs(ax.XLim))*[-1 1];
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.String = '$r$ (mm)';    
    ax.TickLabelInterpreter = 'latex';    
    if  strcmp(Geometry,'Rod') || (strcmp(Geometry,'Pipe') && ToggleInnerFluid)
        ax.YLim = [0 r(end)*1e3];
    else
        ax.YLim = [r(1) r(end)]*1e3;
    end
    if  strcmp(Geometry,'Pipe') && ToggleOuterFluid && ToggleInnerFluid && ShowHalfSpace
        if  strcmp(OuterFluid.Name,InnerFluid.Name)
            OuterFaceAlpha = .2;
            InnerFaceAlpha = .2;
        else
            if  OuterFluid.Density*OuterFluid.Velocity > InnerFluid.Density*InnerFluid.Velocity
                OuterFaceAlpha = .2;
                InnerFaceAlpha = .1;
            else
                OuterFaceAlpha = .1;
                InnerFaceAlpha = .2;
            end
        end
        patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[Ro*1e3 Ro*1e3 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',OuterFaceAlpha,'EdgeColor','none','PickableParts','none')
        patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 Ri*1e3 Ri*1e3],'b','FaceAlpha',InnerFaceAlpha,'EdgeColor','none','PickableParts','none')
    elseif (strcmp(Geometry,'Pipe') && ToggleOuterFluid && ~ToggleInnerFluid && ShowHalfSpace) || (strcmp(Geometry,'Rod') && ToggleOuterFluid && ShowHalfSpace)
        patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[Ro*1e3 Ro*1e3 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none','PickableParts','none')
    elseif strcmp(Geometry,'Pipe') && (~ToggleOuterFluid ||(ToggleOuterFluid && ~ShowHalfSpace)) && ToggleInnerFluid
        patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 Ri*1e3 Ri*1e3],'b','FaceAlpha',.2,'EdgeColor','none','PickableParts','none')
    end
    if  Quantity == 1
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
    elseif Quantity == 2
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
    elseif Quantity == 3
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
    elseif Quantity == 4
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
    elseif Quantity == 5
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
    datacursormode on
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
        if  Quantity == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$u_z$: \textbf{',num2str(event_obj.Position(1),6),'}\,nm'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$u_\theta$: \textbf{',num2str(event_obj.Position(1),6),'}\,nm'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$u_r$: \textbf{',num2str(event_obj.Position(1),6),'}\,nm'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            end
        elseif Quantity == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\sigma_{rr}$: \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\sigma_{zr}$: \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\sigma_{\theta r}$: \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\sigma_{zz}$: \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\sigma_{\theta\theta}$: \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\sigma_{z\theta}$: \textbf{',num2str(event_obj.Position(1),6),'}\,kPa'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};            
            end
        elseif Quantity == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varepsilon_{rr}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varepsilon_{zr}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varepsilon_{\theta r}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varepsilon_{zz}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varepsilon_{\theta\theta}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varepsilon_{z\theta}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};            
            end
        elseif Quantity == 4
            if  event_obj.Target.Color == Color1
                output_txt = {['$E_\mathrm{strain}$: \textbf{',num2str(event_obj.Position(1),6),'}\,J/m$^2$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$E_\mathrm{kin}$: \textbf{',num2str(event_obj.Position(1),6),'}\,J/m$^2$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$E_\mathrm{total}$: \textbf{',num2str(event_obj.Position(1),6),'}\,J/m$^2$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            end
        elseif Quantity == 5 
            if  event_obj.Target.Color == Color1
                output_txt = {['$p_z$: \textbf{',num2str(event_obj.Position(1),6),'}\,W/m'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$p_\theta$: \textbf{',num2str(event_obj.Position(1),6),'}\,W/m'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$p_r$: \textbf{',num2str(event_obj.Position(1),6),'}\,W/m'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
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
        if  Quantity == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$\varphi(u_z)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(u_\theta)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$\varphi(u_r)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            end
        elseif Quantity == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\sigma_{rr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\sigma_{zr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\sigma_{\theta r})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\sigma_{zz})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\sigma_{\theta\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\sigma_{z\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};            
            end
        elseif Quantity == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\varepsilon_{rr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\varepsilon_{zr})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\varepsilon_{\theta r})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\varepsilon_{zz})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\varepsilon_{\theta\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\varepsilon_{z\theta})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$r$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};            
            end
        end
    end
end
end