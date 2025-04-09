clear
clc

Directory = 'F:\ARHGAP18 OE DMSO ERKi 1st try\Cell 1\New folder';
               %'D:\Akib\RhoA2G 10th try\RhoA2G 1H2C\output\BiosensorsPackage\ratio_tiffs';

MovieWriteDir = 'D:\Akib\LOK mutant FRET baselevel 1st try\LOK WT 1\Cell 1\output';               
DoMovie = false;

dp = 108; % nanometers
dt = 20; % seconds

SmoothPath = true; % (or false)
nPts = 11; % Number of points to smooth path (ideally odd number)

T_treatment = 91; % Timepoint treatment is added
Tc = 100*dt; % Cuttoff period for low pass filtering of cell centroid XY coordinates
PixelPrctileThreshold = 60;
TPtRange = [1,90];
DistPercentile = 75;
%=======================================================================================

[files,~] = FilesAndFolders(Directory);
nT = size(files,1); % Number of time points

info = imfinfo( fullfile(Directory,files{1}) );
ImageData = zeros(info.Height,info.Width,nT);
ImageMask = false(info.Height,info.Width,nT);
[x,y] = meshgrid( dp*(1:size(ImageData,2)), dp*(1:size(ImageData,1)) );
COM = NaN(nT,2);

if DoMovie        
        prompt = {'What would you like to name the movie?'};
        dlgtitle = 'Movie name?';
        fieldsize = [1 45];
        defaultinput = {'Movie1'};
        answer = inputdlg(prompt,dlgtitle,fieldsize,defaultinput);
        if isempty(answer)
            DoMovie = false;
            MovieName = '';
        else
            MovieName = answer{1};
        end
end

% Read in all images into ------------------------------------------------
for t = 1:nT
    ImageData(:,:,t) = double( imread(fullfile(Directory,files{t})) ) ;
    ImageMask(:,:,t) = logical(ImageData(:,:,t));
    COM(t,1) = mean( x(ImageMask(:,:,t)) );
    COM(t,2) = mean( y(ImageMask(:,:,t)) );
end

ThresholdData = ImageData(:,:,TPtRange(1):TPtRange(2));
Threshold = prctile(ThresholdData(ThresholdData>0),PixelPrctileThreshold);

% Smooth path ------------------------------------------------------------
if SmoothPath
    %Method 1
    % Xs = smooth(COM(:,1),nPts);
    % Ys = smooth(COM(:,2),nPts);

    %Method 2
    % fn = 1/(2*dt);
    % fc = 1/Tc;
    % Wn = fc/fn;
    % [b,a] = butter(2,Wn,'low');
    % Xs = filtfilt(b,a,COM(:,1));
    % Ys = filtfilt(b,a,COM(:,2));

    % % Method 3
    Px = polyfit(1:size(COM,1),COM(:,1),1);
    Py = polyfit(1:size(COM,1),COM(:,2),1);
    Xs = polyval(Px,1:size(COM,1));
    Ys = polyval(Py,1:size(COM,1));

    figure(12); clf
    set(gcf,'Color','w')
    subplot(1,2,1); plot(COM(:,1),'b'); ylabel('X'); hold on; plot(Xs,'-r')
    subplot(1,2,2); plot(COM(:,2),'b'); ylabel('Y'); hold on; plot(Ys,'-r')    
else
    Xs = COM(:,1);
    Ys = COM(:,2);
end
%-------------------------------------------------------------------------

% Calculate direction vector ---------------------------------------------
xi = interp1(1.5:1:nT-0.5,diff(Xs)',1:1:nT,'linear','extrap');
yi = interp1(1.5:1:nT-0.5,diff(Ys)',1:1:nT,'linear','extrap');

FH = figure(1); clf; 
set(FH,'Color','w')
set(FH,'Position',[-1853, 275, 1816, 754])
FH.Name = Directory;

TL = tiledlayout(1,3,"TileSpacing","compact","Padding","compact");
P = prctile(ImageData(ImageData>0),[1,99]);

FrontMaskArea   = NaN(nT,1);
BackMaskArea    = NaN(nT,1);
FrontAreaAboveT = NaN(nT,1);
BackAreaAboveT  = NaN(nT,1);
FrontMaskMean   = NaN(nT,1);
BackMaskMean    = NaN(nT,1);
FrontMaskSum    = NaN(nT,1);
BackMaskSum     = NaN(nT,1);
FrontTipCOM     = NaN(nT,2);
BackTipCOM      = NaN(nT,2);
TotalMean       = NaN(nT,1);
TotalSum        = NaN(nT,1);
COMDistFront    = NaN(nT,1);
COMDistBack     = NaN(nT,1);

if DoMovie
    v = VideoWriter( fullfile(MovieWriteDir,MovieName),'MPEG-4');
    v.FrameRate = 5;
    open(v)
end

for t = 1:nT   
    nexttile(TL,1); cla
        imagesc(x(1,:),y(:,1)',ImageData(:,:,t)); axis equal tight; hold on
        %set(gca,'YDir','normal')
        colormap(gca,[[1,1,1];turbo])
        clim([P(1),P(2)])
        plot(Xs(t),Ys(t),'.r','MarkerSize',30)
        plot(COM(:,1),COM(:,2),'-k','LineWidth',1); 
        hold off       
        title(['Timepoint = ',num2str(t)],'FontSize',14,'Interpreter','none','HorizontalAlignment','left')
        
    nexttile(TL,2); cla
        Vector = [xi(t);yi(t)]/sqrt(xi(t)^2 + yi(t)^2);
        m = -xi(t)/yi(t); % Slope of line perpendicular to direction vector (negative inverse)
        b = -m*Xs(t) + Ys(t);
        Xperp = dp*find( any(ImageMask(:,:,t),1) );
        Yperp = m*Xperp + b;
        % Create mask for backside of perpendicular line
        BW1 = false(size(x));
        if Vector(2) < 0
            BW1( (m*x+b-y) < 0) = true;
        else
            BW1( (m*x+b-y) > 0) = true;
        end
        BW_Perp_Line = false(size(x));
        BW_Perp_Line((m*x+b-y) > -1.5*dp & (m*x+b-y) < 1.5*dp) = true;
        BW_Perp_Line = BW_Perp_Line.*ImageMask(:,:,t);
        %------------------------------------------------

        BackOfCellMask  = logical(ImageMask(:,:,t).* BW1);
        BackOfCellDist  = BackOfCellMask.*bwdist(~BW1);
        BackTipThresh = prctile(BackOfCellDist(BackOfCellDist>0),DistPercentile);      %(DistPercentile/100)*max(BackOfCellDist(:));
        BackOfCellTipMask = false(size(BackOfCellMask));
        BackOfCellTipMask(BackOfCellDist >= BackTipThresh) = true;

        FrontOfCellMask = logical(ImageMask(:,:,t).*~BW1);
        FrontOfCellDist = FrontOfCellMask.*bwdist(BW1);
        FrontTipThresh = prctile(FrontOfCellDist(FrontOfCellDist>0),DistPercentile);      %(DistPercentile/100)*max(FrontOfCellDist(:));
        FrontOfCellTipMask = false(size(FrontOfCellMask));
        FrontOfCellTipMask(FrontOfCellDist >= FrontTipThresh) = true;

        FrontTipCOM(t,:) = [mean(x(FrontOfCellTipMask)), mean(y(FrontOfCellTipMask))];
        BackTipCOM(t,:)  = [mean(x(BackOfCellTipMask)), mean(y(BackOfCellTipMask))];

        MaskPerim = bwperim(ImageMask(:,:,t));
        imagesc(x(1,:),y(:,1)',logical(MaskPerim + BW_Perp_Line + BackOfCellTipMask + 2*FrontOfCellTipMask)); axis equal tight; hold on
        %set(gca,'YDir','normal')
        colormap(gca,gray)
        plot(Xs(t),Ys(t),'.r','MarkerSize',30)
        plot(FrontTipCOM(t,1),FrontTipCOM(t,2),'.b','MarkerSize',20)
        plot(BackTipCOM(t,1),BackTipCOM(t,2),'.b','MarkerSize',20)
        plot(Xs,Ys,'-c'); 
        % plot(Xperp,Yperp,':r')        
        L = 300;
        plot([Xs(t);Xs(t)+Vector(1)*L*dt],[Ys(t);Ys(t)+Vector(2)*L*dt],'-r','LineWidth',2)
        hold off
        %title(Directory,'Interpreter','none','FontSize',10)
        
    nexttile(TL,3); cla        
        %Threshold = 1050;
        IMmaskedFront = ImageData(:,:,t).*FrontOfCellTipMask;
        IMmaskedBack  = ImageData(:,:,t).*BackOfCellTipMask;
        FrontBackCombinedTip = IMmaskedFront + IMmaskedBack;
        FrontBackCombinedTip(FrontBackCombinedTip < Threshold) = NaN;
        % IM = IMmasked;
        % IMmasked(IMmasked<Threshold) = NaN;
        imagesc(x(1,:),y(:,1)',FrontBackCombinedTip); axis equal tight; hold on
        %imagesc(x(1,:),y(:,1)',IMmaskedBack);
        %set(gca,'YDir','normal')
        colormap(gca,turbo)
        clim([P(1),P(2)])
        plot(Xs(t),Ys(t),'.r','MarkerSize',30)
        plot(Xs,Ys,'-w'); 

        ImDat = ImageData(:,:,t);

        FrontMaskArea(t,1)   = sum(FrontOfCellTipMask(:));
        BackMaskArea(t,1)    = sum(BackOfCellTipMask(:));
        
        FrontAreaAboveT(t,1) = numel(find(IMmaskedFront > Threshold));
        BackAreaAboveT(t,1)  = numel(find(IMmaskedBack >  Threshold));
        
        FrontMaskMean(t,1)   = mean( ImDat(FrontOfCellTipMask),'omitnan');
        BackMaskMean(t,1)    = mean( ImDat(BackOfCellTipMask),'omitnan');
        
        FrontMaskSum(t,1)    = sum( ImDat(FrontOfCellTipMask),'omitnan');
        BackMaskSum(t,1)     = sum( ImDat(BackOfCellTipMask),'omitnan');
        
        TotalMean(t,1)       = mean( ImDat(ImageMask(:,:,t)) );
        TotalSum(t,1)        = sum( ImDat(ImageMask(:,:,t)) );
        COMDistFront(t,1)    = sqrt(  ( Xs(t)-FrontTipCOM(t,1) ).^2 + ( Ys(t)-FrontTipCOM(t,2) ).^2 );
        COMDistBack(t,1)     = sqrt(  ( Xs(t)-BackTipCOM(t,1)  ).^2 + ( Ys(t)-BackTipCOM(t,2)  ).^2 );

        if DoMovie
            frame = getframe(FH);
            writeVideo(v,frame)
        end
        pause(0.01)     
end

if DoMovie
    close(v)
end

XTickDelta = 20;
XTickValues = 0:XTickDelta:size(ImageData,3); %floor(nT/XTickDelta)*XTickDelta;

%=========================================================================================================================================
%=========================================================================================================================================

FH = figure(2); clf; set(gcf,'Color','w')
set(2,'Position',[ -1251 , 46, 928 , 1053])
FS = 8; % Figure font size
TL2 = tiledlayout(FH,4,3,"TileSpacing","tight","Padding","loose");
nexttile(TL2,1); cla
    Velocity = sqrt(diff(COM(:,1)).^2 + diff(COM(:,2)).^2)/dt;
    %VelocityInterp = interp1(1.5:1:nT-0.5,Velocity',1:nT,'linear','extrap');
    plot(1:nT-1, Velocity  ); hold on
    plot([1,nT],mean(Velocity).*ones(1,2),'-r'); hold off; grid on
    %xlabel('Timepoint')
    ylabel('Velocity (nm/s)')
    title('Cell center of mass velocity','FontSize',FS)
    set(gca,'FontSize',FS,'YLim',[0,round(1.1*max(Velocity))],'XLim',[0,nT],'XTick',XTickValues)

% nexttile(TL2,4); cla
%     P1 = plot(1:nT, FrontMaskArea );  hold on
%     P2 = plot(1:nT, BackMaskArea ); hold off; grid on
%     xlabel('Timepoint')
%     ylabel('Mask area (pixels)')
%     set(gca,'FontSize',FS,'XLim',[0,nT])
%     title('Mask area','FontSize',FS)
%     legend([P1,P2],[{'Front'},{'Back'}])

nexttile(TL2,2); cla
    P1 = plot(1:nT, FrontAreaAboveT );  hold on
    P2 = plot(1:nT, BackAreaAboveT ); hold off; grid on
    %xlabel('Timepoint')
    ylabel('Mask area (pixels)')
    set(gca,'FontSize',FS,'XLim',[0,nT],'XTick',XTickValues)
    title('Mask area above threshold','FontSize',FS)
    legend([P1,P2],[{'Front'},{'Back'}])

nexttile(TL2,3); cla
    P1 = plot(1:nT, FrontMaskMean );  hold on
    P2 = plot(1:nT, BackMaskMean );   hold off; grid on
    %xlabel('Timepoint')
    ylabel('Mean')
    set(gca,'FontSize',FS,'XLim',[0,nT],'XTick',XTickValues)
    title('Mask pixel value mean','FontSize',FS)
    legend([P1,P2],[{'Front'},{'Back'}])

 nexttile(TL2,4); cla
    VelocityFront = sqrt(diff(FrontTipCOM(:,1)).^2 + diff(FrontTipCOM(:,2)).^2)/dt; 
    P1 = plot(1:nT-1, VelocityFront  ); hold on ; grid on
    %xlabel('Timepoint')
    ylabel('Velocity (nm/s)')
    set(gca,'FontSize',FS,'XLim',[0,nT],'XTick',XTickValues)
    ylim([0,40])
    title('Front mask COM velocity','FontSize',FS)

 nexttile(TL2,5); cla
    VelocityBack = sqrt(diff(BackTipCOM(:,1)).^2 + diff(BackTipCOM(:,2)).^2)/dt; 
    VelocityBackSmooth = smoothdata(VelocityBack,1,'gaussian',7);
    P1 = plot(1:nT-1, VelocityBack,'-r','Color',[0.8500    0.3250    0.0980]); hold on ; grid on
    P2 = plot(1:nT-1, VelocityBackSmooth,'--k'); hold off ;
    %xlabel('Timepoint')
    ylabel('Velocity (nm/s)')
    set(gca,'FontSize',FS,'XLim',[0,nT],'XTick',XTickValues)
    ylim([0,40])
    title('Back mask COM velocity','FontSize',FS)

nexttile(TL2,7); cla
    P1 = plot(1:nT, FrontMaskArea );  hold on
    P2 = plot(1:nT, BackMaskArea ); hold off; grid on
    ylabel('Mask area (pixels)')
    set(gca,'FontSize',FS,'XLim',[0,nT],'XTick',XTickValues)
    %ylim([0,40])
    title('Mask Areas','FontSize',FS)
    legend([P1,P2],[{'Front'},{'Back'}])

nexttile(TL2,6); cla
    P1 = plot(1:nT, FrontMaskSum );  hold on
    P2 = plot(1:nT, BackMaskSum );  hold off; grid on
    %xlabel('Timepoint')
    ylabel('Sum')
    set(gca,'FontSize',FS,'XLim',[0,nT],'XTick',XTickValues)
    title('Mask pixel value sum','FontSize',FS)
    legend([P1,P2],[{'Front'},{'Back'}])

 nexttile(TL2,8); cla
    [c,lags] = xcorr( VelocityBack(1:T_treatment), BackMaskMean(1:T_treatment),'normalized');
    P1 = stem(lags,c,'.','MarkerSize',5,'MarkerEdgeColor','r','Color','b','LineWidth',1);
    grid on
    xlabel('Lags')
    ylabel('Xcorr')
    set(gca,'FontSize',FS)
    title('Xcorr: Before treatment','FontSize',FS)

 nexttile(TL2,9); cla
    [c,lags] = xcorr( VelocityBack(end-90:end), BackMaskMean(end-90:end),'normalized');
    P1 = stem(lags,c,'.','MarkerSize',5,'MarkerEdgeColor','r','Color','b','LineWidth',1);
    grid on
    xlabel('Lags')
    ylabel('Xcorr')
    set(gca,'FontSize',FS)
    title('Xcorr: last 90 time points','FontSize',FS)

nexttile(TL2,10); cla
    P1 = plot(1:nT, TotalMean);  grid on
    xlabel('Timepoint')
    ylabel('Mean')
    set(gca,'FontSize',FS,'XLim',[0,nT])
    title('Total mean','FontSize',FS)

nexttile(TL2,11); cla
    P1 = plot(1:nT, TotalSum );  grid on
    xlabel('Timepoint')
    ylabel('Sum')
    set(gca,'FontSize',FS,'XLim',[0,nT])
    title('Total sum','FontSize',FS)

nexttile(TL2,12); cla
    P1 = plot(1:nT, COMDistFront );  hold on
    P2 = plot(1:nT, COMDistBack );  hold off; grid on
    xlabel('Timepoint')
    ylabel('Distance')
    set(gca,'FontSize',FS,'XLim',[0,nT],'XTick',XTickValues)
    title('Center of mass distances','FontSize',FS)
    legend([P1,P2],[{'Front'},{'Back'}])

AXtext = axes(FH,'Visible','off','Units','normalized','Position',[0.1,0.97,0.1,0.1]);    
text(AXtext,0,0,Directory,'FontSize',9,'FontWeight','bold')


