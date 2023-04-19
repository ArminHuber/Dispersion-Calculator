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
function PropagationTime_Anisotropic(Hybrid,Crop,LayupString,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Distance,PNGresolution,SColor,AColor,BColor,A,ALamb,AntisymmetricModes,AShear,AScholte,B,BLamb,BoxLineWidth,BShear,BScholte,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeModeLabels,FileName,HeadLine,HigherOrderModes,LambModes,LineWidth,Material,ModeLabels,ModeLabel1Y,ModeLabel2Y,ModeLabel3Y,PDF,PlateThickness,PNG,PropagationAngle,S,ShearHorizontalModes,ScholteModes,SLamb,SShear,SScholte,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,XAxis,YAxisMode,YAxis,Decoupled)
%#ok<*FXUP>
%#ok<*AGROW>
%#ok<*CHAIN>
f = figure('Name','Propagation time','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
hold on
if  Symmetric
    if  ~Decoupled
        if  SymmetricModes && ~isempty(S{1})
            S{1}(:,5) = Distance./S{1}(:,5); % calculate propagation time
            s = plot(S{1}(:,5),S{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',SColor);
            if  ModeLabels
                z = find(abs(S{1}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(S{1}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(S{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1),'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 2
                    text(S{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 3
                    text(S{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            S{2}(:,5) = Distance./S{2}(:,5); % calculate propagation time
            s(2) = plot(S{2}(:,5),S{2}(:,YAxisMode),'LineWidth',LineWidth,'Color',SColor);
            if  ModeLabels
                z = find(abs(S{2}(:,1)-(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(S{2}(:,1)-(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(S{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1),'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 2
                    text(S{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 3
                    text(S{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 3:size(S,2)
                    S{i}(:,5) = Distance./S{i}(:,5); % calculate propagation time
                    s(i) = plot(S{i}(S{i}(:,5) > 0,5),S{i}(S{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',SColor);
                end
            end
        end
        if  AntisymmetricModes && ~isempty(A{1})
            A{1}(:,5) = Distance./A{1}(:,5);
            a = plot(A{1}(:,5),A{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',AColor);
            if  ModeLabels
                z = find(abs(A{1}(:,1)-(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(A{1}(:,1)-(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(A{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1),'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 2
                    text(A{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 3
                    text(A{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 2:size(A,2)
                    A{i}(:,5) = Distance./A{i}(:,5);
                    a(i) = plot(A{i}(A{i}(:,5) > 0,5),A{i}(A{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',AColor);
                end
            end
        end
    else
        if  LambModes && SymmetricModes && ~isempty(SLamb{1})
            SLamb{1}(:,5) = Distance./SLamb{1}(:,5); % calculate propagation time
            sLamb = plot(SLamb{1}(:,5),SLamb{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',SColor);
            if  ModeLabels
                z = find(abs(SLamb{1}(:,1)-(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(SLamb{1}(:,1)-(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(SLamb{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1),'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 2
                    text(SLamb{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 3
                    text(SLamb{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 2:size(SLamb,2)
                    SLamb{i}(:,5) = Distance./SLamb{i}(:,5); % calculate propagation time
                    sLamb(i) = plot(SLamb{i}(SLamb{i}(:,5) > 0,5),SLamb{i}(SLamb{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',SColor);
                end
            end
        end
        if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
            ALamb{1}(:,5) = Distance./ALamb{1}(:,5);
            aLamb = plot(ALamb{1}(:,5),ALamb{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',AColor);
            if  ModeLabels
                z = find(abs(ALamb{1}(:,1)-(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(ALamb{1}(:,1)-(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(ALamb{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1),'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 2
                    text(ALamb{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 3
                    text(ALamb{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 2:size(ALamb,2)
                    ALamb{i}(:,5) = Distance./ALamb{i}(:,5);
                    aLamb(i) = plot(ALamb{i}(ALamb{i}(:,5) > 0,5),ALamb{i}(ALamb{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',AColor);
                end
            end
        end
        if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
            SShear{1}(:,5) = Distance./SShear{1}(:,5);
            sShear = plot(SShear{1}(:,5),SShear{1}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            if  ModeLabels
                z = find(abs(SShear{1}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(SShear{1}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(SShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1),'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 2
                    text(SShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 3
                    text(SShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 2:size(SShear,2)
                    SShear{i}(:,5) = Distance./SShear{i}(:,5);
                    sShear(i) = plot(SShear{i}(:,5),SShear{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                end
            end
        end
        if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
            for i = 1:size(AShear,2)
                AShear{i}(:,5) = Distance./AShear{i}(:,5);
                aShear(i) = plot(AShear{i}(:,5),AShear{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
        SScholte{1}(:,5) = Distance./SScholte{1}(:,5); % calculate propagation time
        sScholte = plot(SScholte{1}(:,5),SScholte{1}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SScholte,2)
                SScholte{i}(:,5) = Distance./SScholte{i}(:,5); % calculate propagation time
                sScholte(i) = plot(SScholte{i}(:,5),SScholte{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
        AScholte{1}(:,5) = Distance./AScholte{1}(:,5);
        aScholte = plot(AScholte{1}(:,5),AScholte{1}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(AScholte,2)
                AScholte{i}(:,5) = Distance./AScholte{i}(:,5);
                aScholte(i) = plot(AScholte{i}(:,5),AScholte{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end    
else
    if  ~Decoupled
        if  ~isempty(B{1})
            B{1}(:,5) = Distance./B{1}(:,5);
            b = plot(B{1}(:,5),B{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(B{1}(:,1)-(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(B{1}(:,1)-(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(B{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1),'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 2
                    text(B{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 3
                    text(B{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            B{2}(:,5) = Distance./B{2}(:,5); % calculate propagation time
            b(2) = plot(B{2}(:,5),B{2}(:,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(B{2}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(B{2}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(B{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1),'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 2
                    text(B{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 3
                    text(B{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            B{3}(:,5) = Distance./B{3}(:,5);
            b(3) = plot(B{3}(:,5),B{3}(:,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(B{3}(:,1)-(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(B{3}(:,1)-(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(B{3}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1),'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 2
                    text(B{3}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 3
                    text(B{3}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 4:size(B,2)
                    B{i}(:,5) = Distance./B{i}(:,5);
                    b(i) = plot(B{i}(B{i}(:,5) > 0,5),B{i}(B{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
                end
            end
        end
    else
        if  LambModes
            BLamb{1}(:,5) = Distance./BLamb{1}(:,5);
            bLamb = plot(BLamb{1}(:,5),BLamb{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(BLamb{1}(:,1)-(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(BLamb{1}(:,1)-(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(BLamb{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1),'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 2
                    text(BLamb{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 3
                    text(BLamb{1}(z(1),5)+(XAxis(2)-XAxis(1))/80,(ModeLabel3Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            BLamb{2}(:,5) = Distance./BLamb{2}(:,5); % calculate propagation time
            bLamb(2) = plot(BLamb{2}(:,5),BLamb{2}(:,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(BLamb{2}(:,1)-(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(BLamb{2}(:,1)-(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                if  YAxisMode == 1
                    text(BLamb{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1),'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif YAxisMode == 2
                    text(BLamb{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                elseif YAxisMode == 3
                    text(BLamb{2}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel1Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 3:size(BLamb,2)
                    BLamb{i}(:,5) = Distance./BLamb{i}(:,5);
                    bLamb(i) = plot(BLamb{i}(BLamb{i}(:,5) > 0,5),BLamb{i}(BLamb{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
                end
            end
        end
        if  ShearHorizontalModes
            if  SuperLayerSize > 1 && ~SymmetricSystem
                BShear{1}(:,5) = Distance./BShear{1}(:,5);
                bShear = plot(BShear{1}(:,5),BShear{1}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                if  ModeLabels
                    z = find(abs(BShear{1}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(BShear{1}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                    if  YAxisMode == 1
                        text(BShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1),'B$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                    elseif YAxisMode == 2
                        text(BShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'B$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                    elseif YAxisMode == 3
                        text(BShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'B$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                    end
                end
                if  HigherOrderModes
                    for i = 2:size(BShear,2)
                        BShear{i}(:,5) = Distance./BShear{i}(:,5);
                        bShear(i) = plot(BShear{i}(:,5),BShear{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                    end
                end
            elseif SuperLayerSize == 1 || SymmetricSystem
                if  ShearHorizontalModes && ~isempty(SShear{1})
                    SShear{1}(:,5) = Distance./SShear{1}(:,5);
                    sShear = plot(SShear{1}(:,5),SShear{1}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                    if  ModeLabels
                        z = find(abs(SShear{1}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))) == min(abs(SShear{1}(:,1)-(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1)))));
                        if  YAxisMode == 1
                            text(SShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1),'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');                    
                        elseif YAxisMode == 2
                            text(SShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))/1e3,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                        elseif YAxisMode == 3
                            text(SShear{1}(z(1),5)-(XAxis(2)-XAxis(1))/30,(ModeLabel2Y*(YAxis(2)-YAxis(1))+YAxis(1))*PlateThickness,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                        end
                    end
                    if  HigherOrderModes
                        for i = 2:size(SShear,2)
                            SShear{i}(:,5) = Distance./SShear{i}(:,5);
                            sShear(i) = plot(SShear{i}(:,5),SShear{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                        end
                    end
                end
                if  AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                    for i = 1:size(AShear,2)
                        AShear{i}(:,5) = Distance./AShear{i}(:,5);
                        aShear(i) = plot(AShear{i}(:,5),AShear{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
                    end
                end
            end
        end 
    end
    if  ScholteModes && ~isempty(BScholte{1})
        for i = 1:length(BScholte)
            BScholte{i}(:,5) = Distance./BScholte{i}(:,5);
            bScholte(i) = plot(BScholte{i}(:,5),BScholte{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
        end
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
if  HeadLine == 1
    if  YAxisMode ~= 3
       String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' @ ',num2str(Distance),'\,mm'];
    else
       String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' @ ',num2str(Distance),'\,mm'];
    end
elseif ~SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  YAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,'] @ ',num2str(Distance),'\,mm'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,'] @ ',num2str(Distance),'\,mm'];
        end
    elseif Repetitions > 1
        if  YAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'}$ @ ',num2str(Distance),'\,mm'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'}$ @ ',num2str(Distance),'\,mm'];
        end
    end
elseif SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  YAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{\mathrm s}$ @ ',num2str(Distance),'\,mm'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{\mathrm s}$ @ ',num2str(Distance),'\,mm'];
        end
    elseif Repetitions > 1
        if  YAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$ @ ',num2str(Distance),'\,mm'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$ @ ',num2str(Distance),'\,mm'];
        end
    end
end
if  HeadLine > 0
    if  FluidLoading
        if  ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/',char(join(split(LowerFluid.Name,'_'),'\_')));
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/vacuum');
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in vacuum/',char(join(split(LowerFluid.Name,'_'),'\_')));
        end
    end
    ax.Title.String = String;
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
ax.XLabel.String = 'Propagation time ($\mu$s)';
ax.XLim = XAxis;
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSizeAxesLabels;
if  YAxisMode == 1
    ax.YLabel.String = 'Frequency (kHz)';
    ax.YLim = YAxis;
elseif YAxisMode == 2
    ax.YLabel.String = 'Frequency (MHz)';
    ax.YLim = YAxis/1e3;
elseif YAxisMode == 3
    ax.YLabel.String = 'Frequency$\cdot$thickness (MHz$\cdot$mm)';
    ax.YLim = YAxis/1e3*PlateThickness*1e3;
end
ax.TickLabelInterpreter = 'latex';
if  Export
    try
        if  PDF
            if  ~Crop
                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[50 30])
                print(f,fullfile(Directory,FileName),'-dpdf')
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
% c1 = d.UIContextMenu;
function output_txt = Cursor(~,event_obj)
    if  Symmetric
        if  ~Decoupled
            if  ~isempty(S{1}) && all(event_obj.Target.Color == SColor)
                for i = 1:length(S)
                    z = find(S{i}(:,5) > 0);
                    if  i == 1 && event_obj.Target.XData(1) == S{1}(1,5)
                        ModeName = 'S$_0$';
                        break
                    elseif i == 2 && event_obj.Target.XData(1) == S{2}(1,5)
                        ModeName = 'S$_1$';
                        break
                    elseif i > 2 && (event_obj.Target.YData(1) == S{i}(z(1),1) || event_obj.Target.YData(1) == S{i}(z(1),2) || event_obj.Target.YData(1) == S{i}(z(1),3))
                        ModeName = ['S$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  ~isempty(A{1}) && all(event_obj.Target.Color == AColor)
                for i = 1:length(A)
                    z = find(A{i}(:,5) > 0);
                    if  i == 1 && event_obj.Target.XData(1) == A{1}(1,5)
                        ModeName = 'A$_0$';
                        break
                    elseif i > 1 &&  event_obj.Target.YData(1) == A{i}(z(1),1) || event_obj.Target.YData(1) == A{i}(z(1),2) || event_obj.Target.YData(1) == A{i}(z(1),3)
                        ModeName = ['A$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
        else
            if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(SLamb)
                    z = find(SLamb{i}(:,5) > 0);
                    if  event_obj.Target.YData(1) == SLamb{i}(z(1),1) || event_obj.Target.YData(1) == SLamb{i}(z(1),2) || event_obj.Target.YData(1) == SLamb{i}(z(1),3)
                        ModeName = ['S$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
                for i = 1:length(SShear)
                    if  event_obj.Target.YData(1) == SShear{i}(1,1) || event_obj.Target.YData(1) == SShear{i}(1,2) || event_obj.Target.YData(1) == SShear{i}(1,3)
                        ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(ALamb)
                    z = find(ALamb{i}(:,5) > 0);
                    if  event_obj.Target.YData(1) == ALamb{i}(z(1),1) || event_obj.Target.YData(1) == ALamb{i}(z(1),2) || event_obj.Target.YData(1) == ALamb{i}(z(1),3)
                        ModeName = ['A$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
                for i = 1:length(AShear)
                    if  event_obj.Target.YData(1) == AShear{i}(1,1) || event_obj.Target.YData(1) == AShear{i}(1,2) || event_obj.Target.YData(1) == AShear{i}(1,3) 
                        ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                        break
                    end
                end
            end
        end
        if  ~isempty(SScholte{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(SScholte)
                if  event_obj.Target.YData(1) == SScholte{i}(1,1) || event_obj.Target.YData(1) == SScholte{i}(1,2) || event_obj.Target.YData(1) == SScholte{i}(1,3)
                    ModeName = ['S$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    break
                end
            end
        end        
        if  ~isempty(AScholte{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(AScholte)
                if  event_obj.Target.YData(1) == AScholte{i}(1,1) || event_obj.Target.YData(1) == AScholte{i}(1,2) || event_obj.Target.YData(1) == AScholte{i}(1,3)
                    ModeName = ['A$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    break
                end
            end
        end        
    else
        if  ~Decoupled
            if  ~isempty(B{1})
                for i = 1:length(B)
                    z = find(B{i}(:,5) > 0);
                    if  i == 1 && event_obj.Target.XData(1) == B{1}(1,5)
                        ModeName = 'B$_0$';
                        break 
                    elseif i == 2 && event_obj.Target.XData(1) == B{2}(1,5)
                        ModeName = 'B$_1$';
                        break
                    elseif i == 3 && event_obj.Target.XData(1) == B{3}(1,5)
                        ModeName = 'B$_2$';
                        break
                    elseif i > 3 && (event_obj.Target.YData(1) == B{i}(z(1),1) || event_obj.Target.YData(1) == B{i}(z(1),2) || event_obj.Target.YData(1) == B{i}(z(1),3))
                        ModeName = ['B$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end                
        else
            if  ~isempty(BLamb{1}) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(BLamb)
                    z = find(BLamb{i}(:,5) > 0);
                    if  i == 1 && event_obj.Target.XData(1) == BLamb{1}(1,5)
                        ModeName = 'B$_0$';
                        break
                    elseif i == 2 && event_obj.Target.XData(1) == BLamb{2}(1,5)
                        ModeName = 'B$_1$';
                        break
                    elseif i > 2 && (event_obj.Target.YData(1) == BLamb{i}(z(1),1) || event_obj.Target.YData(1) == BLamb{i}(z(1),2) || event_obj.Target.YData(1) == BLamb{i}(z(1),3))
                        ModeName = ['B$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  SuperLayerSize > 1 && ~SymmetricSystem
                if  ~isempty(BShear{1}) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(BShear)
                        if  event_obj.Target.YData(1) == BShear{i}(1,1) || event_obj.Target.YData(1) == BShear{i}(1,2) || event_obj.Target.YData(1) == BShear{i}(1,3)
                            ModeName = ['B$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                            break
                        end
                    end
                end
            elseif SuperLayerSize == 1 || SymmetricSystem
                if  ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(SShear)
                        if  event_obj.Target.YData(1) == SShear{i}(1,1) || event_obj.Target.YData(1) == SShear{i}(1,2) || event_obj.Target.YData(1) == SShear{i}(1,3)
                            ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                            break
                        end
                    end
                end
                if  ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(AShear)
                        if  event_obj.Target.YData(1) == AShear{i}(1,1) || event_obj.Target.YData(1) == AShear{i}(1,2) || event_obj.Target.YData(1) == AShear{i}(1,3)
                            ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                            break
                        end
                    end
                end
            end
        end
        if  ~isempty(BScholte{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(BScholte)
                if  i == 1 && event_obj.Target.XData(1) == BScholte{1}(1,5)
                    ModeName = 'B$^{\mathrm{Scholte}}_0$';
                    break
                elseif i == 2 && event_obj.Target.XData(1) == BScholte{2}(1,5)
                    ModeName = 'B$^{\mathrm{Scholte}}_1$';
                    break
                elseif i > 2 && (event_obj.Target.YData(1) == BScholte{i}(1,1) || event_obj.Target.YData(1) == BScholte{i}(1,2) || event_obj.Target.YData(1) == BScholte{i}(1,3))
                    ModeName = ['B$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    break
                end
            end
        end        
    end
    if  YAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1),6),'} $\mu$s'] ['$f$: \textbf{',num2str(event_obj.Position(2),6),'} kHz']};
    elseif YAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1),6),'} $\mu$s'] ['$f$: \textbf{',num2str(event_obj.Position(2),6),'} MHz']};
    elseif YAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1),6),'} $\mu$s'] ['$f\cdot d$: \textbf{',num2str(event_obj.Position(2),6),'} MHz$\cdot$mm']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end
%     c2 = uicontextmenu;
%     d.UIContextMenu = c2;
%     uimenu(c2,'Label','Hide','Callback',@State_Callback);
%     function State_Callback(~,~)
%         event_obj.Target.LineStyle  = 'none';
%         d.UIContextMenu = c1;
%     end        
end
function ShowModes_Callback(~,~)
    ShowModes_Anisotropic(SuperLayerSize,SymmetricSystem,Symmetric,Decoupled,f.Children(end).Children,SColor,AColor)
end
function Analyze_Callback(~,~)
    f_Analyze = figure('NumberTitle','off','Name','Analyze','Visible','off','MenuBar','none','Position',[0 0 210 60],'CloseRequestFcn',@CloseRequest);
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
    f_Analyze.Units = 'normalized';
    movegui(f_Analyze,'center')
    f_Analyze.Visible = 'on';
    drawnow
    jframe.fHG2Client.getWindow.setAlwaysOnTop(true)
    l = line(ax,0,0,'Color','k');
    dt = line(ax,0,0,'Color','k');
    Mode = 1;
    if  YAxisMode == 1 % (kHz)
        String = {['Propagation time (',char(181),'s)'],'Frequency (kHz)'};
    elseif YAxisMode == 2 % (MHz)
        String = {['Propagation time (',char(181),'s)'],'Frequency (MHz)'};
    elseif YAxisMode == 3 % (MHz*mm)
        String = {['Propagation time (',char(181),'s)'],['f',char(8901),'d (MHz',char(8901),'mm)']};
    end
    Value = .9*XAxis(2);
    uicontrol('Parent',f_Analyze,'Style','popupmenu','Value',Mode,'TooltipString','Select constant quantity.','String',String,'Position',[10 20 130 23],'Callback',@Mode_Callback);
    ValueUI = uicontrol('Parent',f_Analyze,'Style','edit','String',Value,'TooltipString','Enter a value for the above selected quantity.','Position',[150 20 50 23],'Callback',@Value_Callback);
    function Mode_Callback(source,~)
        Mode = source.Value;
        switch source.Value
        case 1 % propagation time
            Value = .9*XAxis(2);
            ValueUI.String = Value;
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        case 2 % frequency
            if  YAxisMode == 1 % (kHz)
                Value = .1*YAxis(2);
            elseif YAxisMode == 2 % (MHz)
                Value = .1*YAxis(2)/1e3;
            elseif YAxisMode == 3 % (MHz*mm)
                Value = .1*YAxis(2)*PlateThickness;
            end
            ValueUI.String = Value;
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function Value_Callback(source,~)
        Value = str2double(source.String);
        if  Mode == 1 % propagation time
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        elseif Mode == 2 % frequency
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function ListModes
        delete(dt)
        if  Symmetric
            if  ~Decoupled
                if  Mode == 1 % propagation time
                    if  AntisymmetricModes && ~isempty(A{1})
                        for i = 1:length(a)
                            if  min(a(i).XData) < Value && max(a(i).XData) > Value
                                z = abs(a(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(a(i),a(i).XData(q),a(i).YData(q));
                            end
                        end
                    end
                    if  SymmetricModes && ~isempty(S{1})
                        for i = 1:length(s)
                            if  min(s(i).XData) < Value && max(s(i).XData) > Value
                                z = abs(s(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(s(i),s(i).XData(q),s(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                        for i = 1:length(aScholte)
                            if  min(aScholte(i).XData) < Value && max(aScholte(i).XData) > Value
                                z = abs(aScholte(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                        for i = 1:length(sScholte)
                            if  min(sScholte(i).XData) < Value && max(sScholte(i).XData) > Value
                                z = abs(sScholte(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                            end
                        end
                    end
                elseif Mode == 2 % frequency
                    if  AntisymmetricModes && ~isempty(A{1})
                        for i = 1:length(a)
                            if  a(i).YData(1) < Value
                                if  a(i).YData(end) >= Value
                                    z = abs(a(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(a(i),a(i).XData(q),a(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  SymmetricModes && ~isempty(S{1})
                        for i = 1:length(s)
                            if  s(i).YData(1) < Value
                                if  s(i).YData(end) >= Value
                                    z = abs(s(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(s(i),s(i).XData(q),s(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                        for i = 1:length(aScholte)
                            if  aScholte(i).YData(1) < Value
                                if  aScholte(i).YData(end) >= Value
                                    z = abs(aScholte(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                        for i = 1:length(sScholte)
                            if  sScholte(i).YData(1) < Value
                                if  sScholte(i).YData(end) >= Value
                                    z = abs(sScholte(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                end
            else
                if  Mode == 1 % propagation time
                    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                        for i = 1:length(aLamb)
                            if  min(aLamb(i).XData) < Value && max(aLamb(i).XData) > Value
                                z = abs(aLamb(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(aLamb(i),aLamb(i).XData(q),aLamb(i).YData(q));
                            end
                        end
                    end
                    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
                        for i = 1:length(sLamb)
                            if  min(sLamb(i).XData) < Value && max(sLamb(i).XData) > Value
                                z = abs(sLamb(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(sLamb(i),sLamb(i).XData(q),sLamb(i).YData(q));
                            end
                        end
                    end
                    if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                        for i = 1:length(aShear)
                            if  min(aShear(i).XData) < Value && max(aShear(i).XData) > Value
                                z = abs(aShear(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                            end
                        end
                    end
                    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                        for i = 1:length(sShear)
                            if  min(sShear(i).XData) < Value && max(sShear(i).XData) > Value
                                z = abs(sShear(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                        for i = 1:length(aScholte)
                            if  min(aScholte(i).XData) < Value && max(aScholte(i).XData) > Value
                                z = abs(aScholte(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                        for i = 1:length(sScholte)
                            if  min(sScholte(i).XData) < Value && max(sScholte(i).XData) > Value
                                z = abs(sScholte(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                            end
                        end
                    end
                elseif Mode == 2 % frequency
                    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                        for i = 1:length(aLamb)
                            if  aLamb(i).YData(1) < Value
                                if  aLamb(i).YData(end) >= Value
                                    z = abs(aLamb(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aLamb(i),aLamb(i).XData(q),aLamb(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
                        for i = 1:length(sLamb)
                            if  sLamb(i).YData(1) < Value
                                if  sLamb(i).YData(end) >= Value
                                    z = abs(sLamb(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sLamb(i),sLamb(i).XData(q),sLamb(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                        for i = 1:length(aShear)
                            if  aShear(i).YData(1) < Value
                                if  aShear(i).YData(end) >= Value
                                    z = abs(aShear(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                        for i = 1:length(sShear)
                            if  sShear(i).YData(1) < Value
                                if  sShear(i).YData(end) >= Value
                                    z = abs(sShear(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                        for i = 1:length(aScholte)
                            if  aScholte(i).YData(1) < Value
                                if  aScholte(i).YData(end) >= Value
                                    z = abs(aScholte(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                        for i = 1:length(sScholte)
                            if  sScholte(i).YData(1) < Value
                                if  sScholte(i).YData(end) >= Value
                                    z = abs(sScholte(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                end
            end
        else
            if  ~Decoupled
                if  Mode == 1 % propagation time
                    if  ~isempty(B{1})
                        for i = 1:length(b)
                            if  min(b(i).XData) < Value && max(b(i).XData) > Value
                                z = abs(b(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(b(i),b(i).XData(q),b(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && ~isempty(BScholte{1})
                        for i = 1:length(bScholte)
                            if  min(bScholte(i).XData) < Value && max(bScholte(i).XData) > Value
                                z = abs(bScholte(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                            end
                        end
                    end
                elseif Mode == 2 % frequency
                    if  ~isempty(B{1})
                        for i = 1:length(b)
                            if  b(i).YData(1) < Value
                                if  b(i).YData(end) >= Value
                                    z = abs(b(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(b(i),b(i).XData(q),b(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && ~isempty(BScholte{1})
                        for i = 1:length(bScholte)
                            if  bScholte(i).YData(1) < Value
                                if  bScholte(i).YData(end) >= Value
                                    z = abs(bScholte(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                end
            else
                if  Mode == 1 % propagation time
                    if  LambModes && ~isempty(BLamb{1})
                        for i = 1:length(bLamb)
                            if  min(bLamb(i).XData) < Value && max(bLamb(i).XData) > Value
                                z = abs(bLamb(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(bLamb(i),bLamb(i).XData(q),bLamb(i).YData(q));
                            end
                        end
                    end
                    if  SuperLayerSize > 1 && ~SymmetricSystem
                        if  ShearHorizontalModes && ~isempty(BShear{1})
                            for i = 1:length(bShear)
                                if  min(bShear(i).XData) < Value && max(bShear(i).XData) > Value
                                    z = abs(bShear(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(bShear(i),bShear(i).XData(q),bShear(i).YData(q));
                                end
                            end
                        end
                    elseif SuperLayerSize == 1 || SymmetricSystem
                        if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                            for i = 1:length(aShear)
                                if  min(aShear(i).XData) < Value && max(aShear(i).XData) > Value
                                    z = abs(aShear(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                                end
                            end
                        end
                        if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                            for i = 1:length(sShear)
                                if  min(sShear(i).XData) < Value && max(sShear(i).XData) > Value
                                    z = abs(sShear(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                                end
                            end
                        end
                    end
                    if  ScholteModes && ~isempty(BScholte{1})
                        for i = 1:length(bScholte)
                            if  min(bScholte(i).XData) < Value && max(bScholte(i).XData) > Value
                                z = abs(bScholte(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                            end
                        end
                    end
                elseif Mode == 2 % frequency
                    if  LambModes && ~isempty(BLamb{1})
                        for i = 1:length(bLamb)
                            if  bLamb(i).YData(1) < Value
                                if  bLamb(i).YData(end) >= Value
                                    z = abs(bLamb(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(bLamb(i),bLamb(i).XData(q),bLamb(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  SuperLayerSize > 1 && ~SymmetricSystem
                        if  ShearHorizontalModes && ~isempty(BShear{1})
                            for i = 1:length(bShear)
                                if  bShear(i).YData(1) < Value
                                    if  bShear(i).YData(end) >= Value
                                        z = abs(bShear(i).YData-Value);
                                        q = find(z == min(z),1);
                                        dt(end+1) = datatip(bShear(i),bShear(i).XData(q),bShear(i).YData(q));
                                    end
                                else
                                    break
                                end
                            end
                        end
                    elseif SuperLayerSize == 1 || SymmetricSystem
                        if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                            for i = 1:length(aShear)
                                if  aShear(i).YData(1) < Value
                                    if  aShear(i).YData(end) >= Value
                                        z = abs(aShear(i).YData-Value);
                                        q = find(z == min(z),1);
                                        dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                                    end
                                else
                                    break
                                end
                            end
                        end
                        if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                            for i = 1:length(sShear)
                                if  sShear(i).YData(1) < Value
                                    if  sShear(i).YData(end) >= Value
                                        z = abs(sShear(i).YData-Value);
                                        q = find(z == min(z),1);
                                        dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                                    end
                                else
                                    break
                                end
                            end
                        end
                    end
                    if  ScholteModes && ~isempty(BScholte{1})
                        for i = 1:length(bScholte)
                            if  bScholte(i).YData(1) < Value
                                if  bScholte(i).YData(end) >= Value
                                    z = abs(bScholte(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    function CloseRequest(~,~)
        delete(l)
        delete(dt)
        delete(f_Analyze)
    end
end
end