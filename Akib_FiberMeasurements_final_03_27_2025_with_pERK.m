clear
clc
DataDirectory = 'F:\High Low pERk phalloidin 1st replicate - Copy';
[Files,~] = FilesAndFolders(DataDirectory);
idx488 = contains(Files,'488');
Files488 = Files(idx488);
idx561 = contains(Files,'561');
Files561 = Files(idx561);

nCases = 1; %!!!!
BranchDensity  = cell(1,nCases);
FiberThickness = cell(1,nCases);
FiberLengthSum = cell(1,nCases);
FiberLengths   = cell(1,nCases);
Density561     = cell(1,nCases);
Sum561         = cell(1,nCases);
Table = [];

dp = 75.6; % pixel size (nm)
EdgeWidth = 5000; % (nm) width of mask edge (for removing protrusions from analysis)
AbsNum = 0;
for G = 1:nCases
        % Only grab images from the current Group: A,E,or WT ------------------------------------------------------------
        switch G
            case 1
                idx2 = contains(Files488,'Celll');  
            case 2
                idx2 = contains(Files488,'A DMSO'); 
            case 3
                idx2 = contains(Files488,'EZRi'); 
            case 4
                idx2 = contains(Files488,'combined'); 
            case 5
                idx2 = contains(Files488,'EZRi') & contains(Files488,'ERKi'); 
        end

        ImageFiles488 = Files488(idx2);
        ImageFiles561 = Files561(idx2);

        % For each file in this Ggroup ----------------------------------------------------------------------------------
        for F = 1:length(ImageFiles488)
                % Get image info: nZ, Height,Width etc
                    INFO = imfinfo( fullfile(DataDirectory,ImageFiles488{F}) );
                    nZ = length(INFO);
                    ImHeight = INFO.Height;
                    ImWidth  = INFO.Width;    

                % read in image z-stack
                    %ImStack = zeros(ImHeight,ImWidth,nZ,'uint16');
                    
                    %Entropy = zeros(nZ,1);

                    %for z = 1:nZ
                    %    ImStack(:,:,z) = imread( fullfile(DataDirectory,ImageFiles{F}),"Index",z);
                    %    Earray = entropyfilt( squeeze(ImStack(:,:,z)) );
                    %    Entropy(z,1) = sum( Earray(:) );
                    %end
                         
                    %[~,Fidx] = max(Entropy);
                % Z-compress image
                    %IM1 = double(  max( ImStack,[],3) );
                    %MidZ = floor( size(ImStack,3)/2 );
                    IM488 = double( imread( fullfile(DataDirectory,ImageFiles488{F})) );
                    IM561 = double( imread( fullfile(DataDirectory,ImageFiles561{F})) );
                    %IM1 = double(  ImStack(:,:,Fidx) );

                % Crop Image and create cell mask 
                    % Calculate image threshold ---------------------------
                    edges = min(IM488(:)):10:max(IM488(:));
                    [N,edges] = histcounts(IM488(:),edges);
                    Ecent = edges(1:end-1)+ diff(edges(1:2))/2;
                    N = smoothdata(abs(log(N))',1,'gaussian',5);
                    [~,locs1] = findpeaks(N);
                    [~,locs2] = findpeaks(-N);
                    %------------------------------------------------------
                    threshold = Ecent(locs1(1)) + 100;
                    figure(5); clf
                    plot(Ecent,N)
                    title(num2str(Ecent(locs2(1)) - Ecent(locs1(1))))

                    if Ecent(locs2(1)) <= Ecent(locs1(1))
                        threshold = Ecent(locs1(1)) + 100;
                    end

                    BW0 = imfill( imbinarize(IM488, threshold), 'holes');
                    BW0 = bwareafilt(BW0,[1E4,Inf]);
                    [L1,nObj] = bwlabel(BW0); % How many cells are in the image (usually one, but if there are more process them individually)
                    
                    % Loop through each cell in the image
                    for n = 1:nObj
                            AbsNum = AbsNum + 1;
                            BW1 = false(size(BW0));
                            BW1(L1==n) = true; % Create mask for current cell in image
                            BW1(:,[1,end]) = false; % remove mask contact with edge so fill function will work properly
                            BW1([1,end],:) = false;
                            DT  = dp*bwdist(~BW1);
                            BW2 = BW1;
                            BW2(BW1 & (DT>EdgeWidth)) = false;
                            BW3 = bwskel(BW2);
                            BW3 = bwmorph(BW3,'spur',Inf);
                            BW3 = imdilate( imfill(BW3,'holes'), strel('disk',3) ); % Create mask for current cell without protrusions
                        
                            IM2 = IM488.*BW1; %.*~BW2; % Single cell image (if there is more than one cell in image)                            
                            BWperim =  double( bwmorph( bwperim(BW1), 'dilate' ) );

                            % Tophat filtering with line structural element (Hysteresis thresholding)
                            IM2smooth = imgaussfilt(IM2,1);
                            Theta = 180/10;
                            ThetaVector = 0:Theta:(180-Theta);
                            TopHatImages = NaN([size(IM2),length(ThetaVector)]);
                            for v = 1:length(ThetaVector)
                                SE = strel("line",20,ThetaVector(v));
                                TopHatImages(:,:,v) = imtophat(IM2smooth, SE);
                            end
                            TopHatIm = max(TopHatImages,[],3);
                            TopHatIm = BW3.*TopHatIm;

                            % Get only bright fiber structures for width and length measurements ----------------
                            FM2 = ~BW2.*fibermetric(TopHatIm,1:10,"StructureSensitivity",200);
                            BW_fib2 = imbinarize(FM2,0.01);
                            BW_fib2 = bwareafilt(BW_fib2,[10,Inf]);
                            DT2 = bwdist(~BW_fib2);
                            BW_skel2 = bwskel(BW_fib2);
                            FiberThicknessArray = 2*(dp*DT2.*BW_skel2); % Calculate fiber width
                            FiberThicknessValues = FiberThicknessArray(FiberThicknessArray > 0);
                            FiberLengthSumValues = dp*sum(BW_skel2(:));
                            stats = regionprops(BW_skel2,'Area');
                            FiberLengthValues = dp*cell2mat(struct2cell(stats))';

                            % Get all fiber strucures for connection density measurements -----------------------
                            FM1 = ~BW2.*fibermetric(TopHatIm,1,"StructureSensitivity",1);
                            BW_bin1 = imbinarize(FM1,0.01);                  
                            BW_skel1 = bwskel( BW_bin1 );
                            BW_skel1a = bwareafilt(BW_skel1,[10,Inf]);
                            BW_skel1 = ~BW_fib2.*~BW2.*bwareafilt(BW_skel1,[10,Inf]);
                            BW5 = bwmorph(BW_skel1,"branchpoints");
                            [~,nBranchPts] = bwlabel(BW5,8);
                            BD = 1E6*nBranchPts/(sum(BW3(:))*dp^2); 

                            IMtemp = IM561.*BW3;
                            Image561Density = sum(IMtemp(:))/sum(BW3(:));
                            Image561Sum = sum(IMtemp(:));

                            % Record measurements----------------------------------------------------------------
                            BranchDensity{1,G}  = [BranchDensity{1,G}; BD];
                            FiberThickness{1,G} = [FiberThickness{1,G}; mean(FiberThicknessValues/1000)];
                            FiberLengths{1,G}   = [FiberLengths{1,G};   mean(FiberLengthValues/1000)];
                            FiberLengthSum{1,G} = [FiberLengthSum{1,G}, FiberLengthSumValues/1000];
                            Density561{1,G}     = [Density561{1,G}, Image561Density];
                            Sum561{1,G}         = [Sum561{1,G}, Image561Sum];

                            Table(AbsNum,1).FileName = ImageFiles488{F};
                            Table(AbsNum,1).BranchDensity = BD; 
                            Table(AbsNum,1).FiberThickness = mean(FiberThicknessValues/1000); 
                            Table(AbsNum,1).FiberLengths = mean(FiberLengthValues/1000);
                            Table(AbsNum,1).Sum = FiberLengthSumValues/1000;
                            Table(AbsNum,1).Density561 = Image561Density;
                            Table(AbsNum,1).Sum561 = Image561Sum;


                            % Create plots
                            figure(3); clf
                            set(gcf,'Color','w','Position',[-1851 ,144, 1687, 926])
                            TL3 = tiledlayout(2,3,"TileSpacing",'compact',"Padding","tight");
                                AX1 = nexttile(TL3,1); imagesc(IM2); axis equal tight off; title([ImageFiles488{F},'   cell ',num2str(n)],'Interpreter','none')
                                AX2 = nexttile(TL3,4); imagesc( max(TopHatIm(:))*BWperim + TopHatIm); axis equal tight off; title("Smoothed image + tophat filtering")
                                AX3 = nexttile(TL3,2); imagesc( max(FM1(:))*BWperim + FM1); axis equal tight off; title("Fibermetric filter for all fibers (high sensitivity setting)")
                                AX4 = nexttile(TL3,5); imagesc( BWperim + BW_skel1 + BW5); axis equal tight off; title("Mask for all fibers and intersection points")
                                AX5 = nexttile(TL3,3); imagesc( max(FM2(:))*BWperim + FM2); axis equal tight off; title("Fibermetric filter for all fibers (low sensitivity setting)")
                                AX6 = nexttile(TL3,6); imagesc( max(DT2(:))*BWperim + DT2); axis equal tight off; title("Distance transform of high sensitivity fiber mask")
                                linkaxes([AX1, AX2, AX3, AX4, AX5, AX6],'xy')
                            drawnow
                            colormap( [[0,0,0]; turbo(254); [1,1,1]])
                            colormap(AX4, [[0,0,0];[1,1,1]])
                            

                            F2 = getframe(3);
                            if AbsNum == 1
                                imwrite(F2.cdata, fullfile(DataDirectory,'FigureOutputImages.tif'),'WriteMode','overwrite')
                            else
                                imwrite(F2.cdata, fullfile(DataDirectory,'FigureOutputImages.tif'),'WriteMode','append')
                            end
                    end
        end
end

figure(5); set(gcf,'Color','w');clf
TL = tiledlayout(3,2,"TileSpacing","compact","Padding","compact");
% Names = [{'siCtrl DMSO'},{'siCtrl ERKi'},{'siGAP18 DMSO'},{'siGAP18 ERKi'}];
% Names = [{'WT'},{'E'},{'E+ROCKi'}];
% Names = [{'DMSO'},{'ERKi'},{'RhoAi'}];
% Names = [{'DMSO'},{'ERKi'},{'EZRi'},{'ERKi+EZRi'}];
Names = [{'Celll'}];

AX1 = nexttile(TL,1);
    D = CellArray2PaddedNanArray(BranchDensity);
    X = 1:size(D,2);
    Xvals = repmat(X,size(D,1),1);
    SH1 = swarmchart(Xvals,D,'.k','SizeData',50,'XJitterWidth',0.2); hold on;
    BP1 = boxplot(D,'Position',X,'Notch','on','Symbol','','Whisker',0); hold off;
    set(BP1,'LineWidth',2)
    ylabel('Fiber intersection density (\mum^-2)')
    set(gca,'XTickLabel',Names,'FontSize',18,'LineWidth',2)
    title("For dim fibers",'FontSize',12)
    %ylim([0.6,2])

AX2 = nexttile(TL,2);
    D = CellArray2PaddedNanArray(FiberThickness);
    X = 1:size(D,2);
    Xvals = repmat(X,size(D,1),1);
    SH1 = swarmchart(Xvals,D,'.k','SizeData',50,'XJitterWidth',0.2); hold on;
    BP1 = boxplot(D,'Position',X,'Notch','on','Symbol','','Whisker',0); hold off;
    set(BP1,'LineWidth',2)
    ylabel('Mean fiber thickness (\mum)')
    set(gca,'XTickLabel',Names,'FontSize',18,'LineWidth',2)
    title("Bright fibers only",'FontSize',12)
    %ylim([0.24,0.38])

AX3 = nexttile(TL,3);    
    D = CellArray2PaddedNanArray(FiberLengths);
    X = 1:size(D,2);
    Xvals = repmat(X,size(D,1),1);
    SH1 = swarmchart(Xvals,D,'.k','SizeData',50,'XJitterWidth',0.2); hold on;
    BP1 = boxplot(D,'Position',X,'Notch','on','Symbol','','Whisker',0); hold off;
    set(BP1,'LineWidth',2)
    ylabel('Mean fiber lengths (\mum)')
    set(gca,'XTickLabel',Names,'FontSize',18,'LineWidth',2)
    title("Bright fibers only",'FontSize',12)
    %ylim([0.8,2.8])

AX4 = nexttile(TL,4);    
    D = CellArray2PaddedNanArray(FiberLengthSum);
    X = 1:size(D,2);
    Xvals = repmat(X,size(D,1),1);
    SH1 = swarmchart(Xvals,D,'.k','SizeData',50,'XJitterWidth',0.2); hold on;
    BP1 = boxplot(D,'Position',X,'Notch','on','Symbol','','Whisker',0); hold off;
    set(BP1,'LineWidth',2)
    ylabel('Fiber length sum per cell (\mum)')
    set(gca,'XTickLabel',Names,'FontSize',18,'LineWidth',2)  
    title("Bright fibers only",'FontSize',12)
    %ylim([100,1800])

AX5 = nexttile(TL,5);    
    D = CellArray2PaddedNanArray(Density561);
    X = 1:size(D,2);
    Xvals = repmat(X,size(D,1),1);
    SH1 = swarmchart(Xvals,D,'.k','SizeData',50,'XJitterWidth',0.2); hold on;
    BP1 = boxplot(D,'Position',X,'Notch','on','Symbol','','Whisker',0); hold off;
    set(BP1,'LineWidth',2)
    ylabel('561 density (pixel-sum/mask-area)')
    set(gca,'XTickLabel',Names,'FontSize',18,'LineWidth',2)  
    title("Bright fibers only",'FontSize',12)    

AX6 = nexttile(TL,6);    
    D = CellArray2PaddedNanArray(Sum561);
    X = 1:size(D,2);
    Xvals = repmat(X,size(D,1),1);
    SH1 = swarmchart(Xvals,D,'.k','SizeData',50,'XJitterWidth',0.2); hold on;
    BP1 = boxplot(D,'Position',X,'Notch','on','Symbol','','Whisker',0); hold off;
    set(BP1,'LineWidth',2)
    ylabel('561 sum (pixel-sum)')
    set(gca,'XTickLabel',Names,'FontSize',18,'LineWidth',2)  
    title("Bright fibers only",'FontSize',12)    
