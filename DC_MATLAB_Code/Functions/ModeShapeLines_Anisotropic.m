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
function ModeShapeLines_Anisotropic(FunctionMode,MakeFigure,ExportData,XSLX,TXT,MAT,Hybrid,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,LayupString,Plot,PNGresolution,Color1,Color2,Color3,Color4,Color5,Color6,ALamb,AShear,AScholte,BLamb,BShear,BScholte,SLamb,SShear,SScholte,BoxLineWidth,c,Delta,Material,Layers,Directory,FileName,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeLegend,Frequency,HeadLine,LegendLocation,LineWidth,Mode,PDF,PNG,PlateThickness,PropagationAngle,SamplesX3,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,ShowHalfSpaces,HalfSpaces,Phase)
%#ok<*AGROW>
[u,epsilon,sigma,StrainEnergyDensity,KineticEnergyDensity,TotalEnergyDensity,PowerFlowDensity,uPhase,epsilonPhase,sigmaPhase,x3Total,p] = ModeShapeLinesComputer_Anisotropic(FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,ALamb,AShear,AScholte,BLamb,BShear,BScholte,SLamb,SShear,SScholte,c,Delta,Material,Layers,Frequency,Mode,PlateThickness,SamplesX3,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,ShowHalfSpaces,HalfSpaces,Phase);
if  ExportData
    if  FunctionMode == 1
        if  ~Phase
            Table = table('Size',[length(x3Total) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','u1 (nm)','u2 (nm)','u3 (nm)'});
            Table(:,1) = num2cell(-1e3*x3Total);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
        else
            Table = table('Size',[length(x3Total) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','u1 (nm)','u2 (nm)','u3 (nm)','Phase u1 (deg)','Phase u2 (deg)','Phase u3 (deg)'});
            Table(:,1) = num2cell(-1e3*x3Total);
            Table(:,2) = num2cell(real(u(:,1))*1e9);
            Table(:,3) = num2cell(real(u(:,2))*1e9);
            Table(:,4) = num2cell(real(u(:,3))*1e9);
            Table(:,5) = num2cell(uPhase(:,1));
            Table(:,6) = num2cell(uPhase(:,2));
            Table(:,7) = num2cell(uPhase(:,3));
        end
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                Displacement = Table;
                save(fullfile(Directory,[FileName,'_Displacement_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'Displacement')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end    
    elseif FunctionMode == 2
        if  ~Phase
            Table = table('Size',[length(x3Total) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','sigma11 (kPa)','sigma22 (kPa)','sigma33 (kPa)','sigma23 (kPa)','sigma13 (kPa)','sigma12 (kPa)'});
            Table(:,1) = num2cell(-1e3*x3Total);
            Table(:,2) = num2cell(real(sigma(:,1))/1e3);
            Table(:,3) = num2cell(real(sigma(:,2))/1e3);
            Table(:,4) = num2cell(real(sigma(:,3))/1e3);
            Table(:,5) = num2cell(real(sigma(:,4))/1e3);
            Table(:,6) = num2cell(real(sigma(:,5))/1e3);
            Table(:,7) = num2cell(real(sigma(:,6))/1e3);
        else
            Table = table('Size',[length(x3Total) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','sigma11 (kPa)','sigma22 (kPa)','sigma33 (kPa)','sigma23 (kPa)','sigma13 (kPa)','sigma12 (kPa)','Phase sigma11 (deg)','Phase sigma22 (deg)','Phase sigma33 (deg)','Phase sigma23 (deg)','Phase sigma13 (deg)','Phase sigma12 (deg)'});
            Table(:,1) = num2cell(-1e3*x3Total);
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
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                Stress = Table;
                save(fullfile(Directory,[FileName,'_Stress_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'Stress')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end        
    elseif FunctionMode == 3 
        if  ~Phase
            Table = table('Size',[length(x3Total) 7],'VariableTypes',{'double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','epsilon11','epsilon22','epsilon33','epsilon23','epsilon13','epsilon12'});
            Table(:,1) = num2cell(-1e3*x3Total);
            Table(:,2) = num2cell(real(epsilon(:,1)));
            Table(:,3) = num2cell(real(epsilon(:,2)));
            Table(:,4) = num2cell(real(epsilon(:,3)));
            Table(:,5) = num2cell(real(epsilon(:,4)));
            Table(:,6) = num2cell(real(epsilon(:,5)));
            Table(:,7) = num2cell(real(epsilon(:,6)));
        else
            Table = table('Size',[length(x3Total) 13],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'x3 (mm)','epsilon11','epsilon22','epsilon33','epsilon23','epsilon13','epsilon12','Phase epsilon11 (deg)','Phase epsilon22 (deg)','Phase epsilon33 (deg)','Phase epsilon23 (deg)','Phase epsilon13 (deg)','Phase epsilon12 (deg)'});
            Table(:,1) = num2cell(-1e3*x3Total);
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
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                Strain = Table;
                save(fullfile(Directory,[FileName,'_Strain_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'Strain')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end 
    elseif FunctionMode == 4   
        Table = table('Size',[length(x3Total) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','Estrain (J/m2)','Ekin (J/m2)','Etotal (J/m2)'});
        Table(:,1) = num2cell(-1e3*x3Total);
        Table(:,2) = num2cell(StrainEnergyDensity);
        Table(:,3) = num2cell(KineticEnergyDensity);
        Table(:,4) = num2cell(TotalEnergyDensity);
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                EnergyDensity = Table;
                save(fullfile(Directory,[FileName,'_EnergyDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'EnergyDensity')
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export data')
            return
        end
    elseif FunctionMode == 5
        Table = table('Size',[length(x3Total) 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x3 (mm)','P1 (W/m)','P2 (W/m)','P3 (W/m)'});
        Table(:,1) = num2cell(-1e3*x3Total);
        Table(:,2) = num2cell(PowerFlowDensity(:,1));
        Table(:,3) = num2cell(PowerFlowDensity(:,2));
        Table(:,4) = num2cell(PowerFlowDensity(:,3));
        try
            if  XSLX
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.xlsx']));
            end
            if  TXT
                writetable(Table,fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.txt']));
            end
            if  MAT
                PowerFlowDensity = Table;
                save(fullfile(Directory,[FileName,'_PowerFlowDensity_',Mode,'@',num2str(Frequency*PlateThickness),'MHzmm.mat']),'PowerFlowDensity')
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
                plot(uPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
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
        if  Hybrid
            Material{1}.Name = 'hybrid';
        end
        if  HeadLine > 0
            if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
                ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
            elseif contains(Mode,'SH')
                ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
            elseif contains(Mode,'Scholte')
                ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
            end
            if  HeadLine == 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
            elseif ~SymmetricSystem && HeadLine == 2
                if  Repetitions == 1
                    String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
                else
                    String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
                end
            elseif SymmetricSystem && HeadLine == 2
                if  Repetitions == 1
                    String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
                else
                    String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
                end
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
        ax.XLabel.Interpreter = 'latex';
        ax.XLabel.FontSize = FontSizeAxesLabels;
        ax.XLim = max(abs(ax.XLim))*[-1 1];
        ax.YLabel.Interpreter = 'latex';
        ax.YLabel.FontSize = FontSizeAxesLabels;
        ax.YLabel.String = '$x_3$ (mm)';
        ax.TickLabelInterpreter = 'latex';
        ax.YLim = -1e3*[x3Total(end) x3Total(1)];
        if  FluidLoading && ShowHalfSpaces
            if  ToggleUpperFluid && ToggleLowerFluid
                if  strcmp(UpperFluid.Name,LowerFluid.Name)
                    UpperFaceAlpha = .2;
                    LowerFaceAlpha = .2;
                else
                    if  UpperFluid.Density*UpperFluid.Velocity > LowerFluid.Density*LowerFluid.Velocity
                        UpperFaceAlpha = .2;
                        LowerFaceAlpha = .1;
                    else
                        UpperFaceAlpha = .1;
                        LowerFaceAlpha = .2;
                    end
                end
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',UpperFaceAlpha,'EdgeColor','none')
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -PlateThickness*1e3 -PlateThickness*1e3],'b','FaceAlpha',LowerFaceAlpha,'EdgeColor','none')
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none')
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -PlateThickness*1e3 -PlateThickness*1e3],'b','FaceAlpha',.2,'EdgeColor','none')
            end
        end
        if  FunctionMode == 1
            if  Plot(3)
                plot(uPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;
            end
            if  Plot(5)
                plot(uPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;
            end
            if  Plot(4)
                plot(uPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;
            end
        elseif FunctionMode == 2
            if  Plot(3)
                plot(sigmaPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(sigmaPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(sigmaPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(sigmaPhase(:,5),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(sigmaPhase(:,4),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(sigmaPhase(:,6),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end 
        elseif FunctionMode == 3
            if  Plot(3)
                plot(epsilonPhase(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
                z = 1;
            else
                z = 0;        
            end     
            if  Plot(1)
                plot(epsilonPhase(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
                z(2) = 1;
            else
                z(2) = 0;        
            end 
            if  Plot(2)
                plot(epsilonPhase(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
                z(3) = 1;
            else
                z(3) = 0;        
            end
            if  Plot(5)
                plot(epsilonPhase(:,5),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
                z(4) = 1;
            else
                z(4) = 0;        
            end    
            if  Plot(4)
                plot(epsilonPhase(:,4),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
                z(5) = 1;
            else
                z(5) = 0;
            end
            if  Plot(6)
                plot(epsilonPhase(:,6),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
                z(6) = 1;
            else
                z(6) = 0;        
            end
        end
        if  FunctionMode == 1
            ax.XLabel.String = 'Displacement phase ($^\circ$)';
            if  strcmp(LegendLocation,'in')
                LegendNames = {'$u_3$','$u_1$','$u_2$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
            else
                LegendNames = {'Out-of-plane ($u_3$)','In-plane ($u_1$)','Shear horizontal ($u_2$)'};
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
                LegendNames = {'$\sigma_{33}$','$\sigma_{11}$','$\sigma_{22}$','$\sigma_{13}$','$\sigma_{23}$','$\sigma_{12}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Out-of-plane ($\sigma_{33}$)','In-plane ($\sigma_{11}$)','In-plane ($\sigma_{22}$)','Shear ($\sigma_{13}$)','Shear ($\sigma_{23}$)','Shear ($\sigma_{12}$)'}; 
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
                LegendNames = {'$\varepsilon_{33}$','$\varepsilon_{11}$','$\varepsilon_{22}$','$\varepsilon_{13}$','$\varepsilon_{23}$','$\varepsilon_{12}$'};
                LegendNames = LegendNames(z == 1);
                legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
            else
                LegendNames = {'Out-of-plane ($\varepsilon_{33}$)','In-plane ($\varepsilon_{11}$)','In-plane ($\varepsilon_{22}$)','Shear ($\varepsilon_{13}$)','Shear ($\varepsilon_{23}$)','Shear ($\varepsilon_{12}$)'}; 
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
            plot(real(u(:,3))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
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
    if  Hybrid
        Material{1}.Name = 'hybrid';
    end
    if  HeadLine > 0
        if  ~contains(Mode,'Scholte') && ~contains(Mode,'SH')
            ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
        elseif contains(Mode,'SH')
            ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
        elseif contains(Mode,'Scholte')
            ModeName = [Mode(1),'$^{\mathrm{Scholte}}_{',num2str(p-1),'}$'];
        end
        if  HeadLine == 1
            String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
        elseif ~SymmetricSystem && HeadLine == 2
            if  Repetitions == 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
            else
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
            end
        elseif SymmetricSystem && HeadLine == 2
            if  Repetitions == 1
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
            else
                String = [ModeName,' @ ',num2str(Frequency),'\,kHz for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
            end
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
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.XLim = max(abs(ax.XLim))*[-1 1];
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.String = '$x_3$ (mm)';    
    ax.TickLabelInterpreter = 'latex';    
    ax.YLim = -1e3*[x3Total(end) x3Total(1)];
    if  FluidLoading && ShowHalfSpaces
        if  ToggleUpperFluid && ToggleLowerFluid
            if  strcmp(UpperFluid.Name,LowerFluid.Name)
                UpperFaceAlpha = .2;
                LowerFaceAlpha = .2;
            else
                if  UpperFluid.Density*UpperFluid.Velocity > LowerFluid.Density*LowerFluid.Velocity
                    UpperFaceAlpha = .2;
                    LowerFaceAlpha = .1;
                else
                    UpperFaceAlpha = .1;
                    LowerFaceAlpha = .2;
                end
            end
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',UpperFaceAlpha,'EdgeColor','none')
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -PlateThickness*1e3 -PlateThickness*1e3],'b','FaceAlpha',LowerFaceAlpha,'EdgeColor','none')
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[0 0 ax.YLim(2) ax.YLim(2)],'b','FaceAlpha',.2,'EdgeColor','none')
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            patch([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[ax.YLim(1) ax.YLim(1) -PlateThickness*1e3 -PlateThickness*1e3],'b','FaceAlpha',.2,'EdgeColor','none')
        end
    end
    if  FunctionMode == 1
        if  Plot(3)
            plot(real(u(:,3))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(real(u(:,1))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;
        end
        if  Plot(4)
            plot(real(u(:,2))*1e9,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;
        end 
    elseif FunctionMode == 2
        if  Plot(3)
            plot(real(sigma(:,3))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(sigma(:,1))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(sigma(:,2))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(sigma(:,5))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(sigma(:,4))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(sigma(:,6))/1e3,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end 
    elseif FunctionMode == 3
        if  Plot(3)
            plot(real(epsilon(:,3)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;        
        end     
        if  Plot(1)
            plot(real(epsilon(:,1)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end 
        if  Plot(2)
            plot(real(epsilon(:,2)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
        if  Plot(5)
            plot(real(epsilon(:,5)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color5);
            z(4) = 1;
        else
            z(4) = 0;        
        end    
        if  Plot(4)
            plot(real(epsilon(:,4)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color4);
            z(5) = 1;
        else
            z(5) = 0;
        end
        if  Plot(6)
            plot(real(epsilon(:,6)),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color6);
            z(6) = 1;
        else
            z(6) = 0;        
        end
    elseif FunctionMode == 4
        if  Plot(5)
            plot(StrainEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z = 1;
        else
            z = 0;
        end
        if  Plot(4)
            plot(KineticEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(3)
            plot(TotalEnergyDensity,-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z(3) = 1;
        else
            z(3) = 0;        
        end        
    elseif FunctionMode == 5
        if  Plot(3)
            plot(PowerFlowDensity(:,3),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color3);
            z = 1;
        else
            z = 0;
        end
        if  Plot(5)
            plot(PowerFlowDensity(:,1),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color1);
            z(2) = 1;
        else
            z(2) = 0;        
        end
        if  Plot(4)
            plot(PowerFlowDensity(:,2),-x3Total*1e3,'LineWidth',LineWidth,'Color',Color2);
            z(3) = 1;
        else
            z(3) = 0;        
        end
    end
    if  FunctionMode == 1
        ax.XLabel.String = 'Displacement (nm)';
        if  strcmp(LegendLocation,'in')
            LegendNames = {'$u_3$','$u_1$','$u_2$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')
        else
            LegendNames = {'Out-of-plane ($u_3$)','In-plane ($u_1$)','Shear horizontal ($u_2$)'};
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
            LegendNames = {'$\sigma_{33}$','$\sigma_{11}$','$\sigma_{22}$','$\sigma_{13}$','$\sigma_{23}$','$\sigma_{12}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($\sigma_{33}$)','In-plane ($\sigma_{11}$)','In-plane ($\sigma_{22}$)','Shear ($\sigma_{13}$)','Shear ($\sigma_{23}$)','Shear ($\sigma_{12}$)'}; 
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
            LegendNames = {'$\varepsilon_{33}$','$\varepsilon_{11}$','$\varepsilon_{22}$','$\varepsilon_{13}$','$\varepsilon_{23}$','$\varepsilon_{12}$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($\varepsilon_{33}$)','In-plane ($\varepsilon_{11}$)','In-plane ($\varepsilon_{22}$)','Shear ($\varepsilon_{13}$)','Shear ($\varepsilon_{23}$)','Shear ($\varepsilon_{12}$)'}; 
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
            LegendNames = {'$p_3$','$p_1$','$p_2$'};
            LegendNames = LegendNames(z == 1);
            legend(LegendNames,'Location','northeast','FontSize',FontSizeLegend,'Interpreter','latex')        
        else
            LegendNames = {'Out-of-plane ($p_3$)','In-plane ($p_1$)','Shear horizontal ($p_2$)'};
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
            output_txt = replace(UpperFluid.Name,'_','\_');
        else
            output_txt = replace(LowerFluid.Name,'_','\_');
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$u_1$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$u_2$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$u_3$: \textbf{',num2str(event_obj.Position(1),6),'} nm'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\sigma_{33}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\sigma_{13}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\sigma_{23}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\sigma_{11}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\sigma_{22}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\sigma_{12}$: \textbf{',num2str(event_obj.Position(1),6),'} kPa'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varepsilon_{33}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varepsilon_{13}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varepsilon_{23}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varepsilon_{11}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varepsilon_{22}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varepsilon_{12}$: \textbf{',num2str(event_obj.Position(1),6),'}'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 4
            if  event_obj.Target.Color == Color1
                output_txt = {['$E_\mathrm{strain}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$E_\mathrm{kin}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$E_\mathrm{total}$: \textbf{',num2str(event_obj.Position(1),6),'} J/m$^2$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 5 
            if  event_obj.Target.Color == Color1
                output_txt = {['$p_1$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$p_2$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$p_3$: \textbf{',num2str(event_obj.Position(1),6),'} W/m'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        end
    end
end    
function output_txt = CursorPhase(~,event_obj)
    if  length(event_obj.Target.XData) == 4
        if  event_obj.Target.YData(1) == 0
            output_txt = replace(UpperFluid.Name,'_','\_');
        else
            output_txt = replace(LowerFluid.Name,'_','\_');
        end
    else
        if  FunctionMode == 1 
            if  event_obj.Target.Color == Color1
                output_txt = {['$\varphi(u_1)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(u_2)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color3
                output_txt = {['$\varphi(u_3)$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            end
        elseif FunctionMode == 2
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\sigma_{33})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\sigma_{13})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\sigma_{23})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\sigma_{11})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\sigma_{22})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\sigma_{12})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        elseif FunctionMode == 3
            if  event_obj.Target.Color == Color3
                output_txt = {['$\varphi(\varepsilon_{33})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color5
                output_txt = {['$\varphi(\varepsilon_{13})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color4
                output_txt = {['$\varphi(\varepsilon_{23})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color1
                output_txt = {['$\varphi(\varepsilon_{11})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color2
                output_txt = {['$\varphi(\varepsilon_{22})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
            elseif event_obj.Target.Color == Color6
                output_txt = {['$\varphi(\varepsilon_{12})$: \textbf{',num2str(event_obj.Position(1),6),'}\,$^\circ$'],['$x_3$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};            
            end
        end
    end
end
end