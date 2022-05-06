% ModelSingleRun.m Multicellular model of autonomous HES5 noisy/oscillatory
% dynamic coupled between cells by a lateral inhibition Hill function
% interaction. This code runs the model once and generates outputs
% dependent on options in 'Select simulation outputs' section. Full
% description of the model is given in the paper 'Dynamic Switching of 
% Lateral Inhibition Spatial Patterns' (Joshua Hawley, Paul Glendinning, 
% and Nancy Papalopulu) (2022). Correspondence: joshhawley369@gmail.com

%Required toolboxes: Signal analysis toolbox, statistics and machine
%learning toolbox.

clear; clc;
addpath('Functions')

%% Choose the colour of all graphs
graphAppearance = 0;   % 0 = Normal, 1 = Dark grey, 2 = dark green/blue, 3 = light blue

if graphAppearance==0 %Default MATLAB graph appearance
    get(0,'Factory');                           set(0,'defaultfigurecolor',[0.94 0.94 0.94])
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor','w')
    set(0,'defaultAxesXColor','k');             set(0,'defaultAxesYColor','k')
    set(0,'defaultLegendTextColor','k');        set(0,'defaultTextColor','k')
end

if graphAppearance > 0
    if graphAppearance==1   
        backgroundColour=38/255*[1 1 1]; textColour=[1 1 1];
    elseif graphAppearance==2   
        backgroundColour=1/255*[0 28 22]; textColour=[1 1 1];
    elseif graphAppearance==3   
        backgroundColour=1/255*[208 230 229]; textColour=[0 0 0];
    end
    get(0,'Factory');                           set(0,'defaultfigurecolor',backgroundColour)
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor',backgroundColour)
    set(0,'defaultAxesXColor',textColour);      set(0,'defaultAxesYColor',textColour)
    set(0,'defaultLegendTextColor',textColour); set (0,'defaultTextColor',textColour)  
    set(0,'defaultAxesZColor',textColour);
end

%% Define grid size
rows=26;                 % Number of rows of cells in hexagonal grid
cols=1;                  % Number of columns of cells in hexagonal grid
cells=rows*cols;         % Total number of cells

%% Define simulation time, step size, deterministic/stochastic, coupling, and HES5 auotrepression
t0=0;                    % Start time
tf_hours=200;            % Final time (hours)
tf=tf_hours*60;          % Final time (min)

coupledCells=1;          % Simulate with lateral inhibition coupling = 1
boundary=0;              % 0 = 'hard' boundary, 1 = periodic boundary
turnOffAutorepression=0; % Reduce model to just lateral inhibition without autonomous Hes oscillators in each cell
stochastic=1;            % 1 = stochastic, 0 = deterministic

dt=2;                    % Time step size (min) for deterministic simulation (when stochastic = 0) 
if stochastic==1
    dt=1;                % Stochastic uses Euler-Maruyama method rather than Runge-Kutta, so this requires a smaller step size      
end
Nt=(tf-t0)/dt;           % Number of time elements
T=t0:dt:tf;              % Time vector

%% Set up simulation inputs

% For each of the following assingments, 1=yes, 0=no. (Some have additional
% options defined in the comments)

%Differentiation-based perturbation (DBP) algorithm parameters
CrudeDiff = 1;                     % 1 = include DBP, 0 = no DBP
diffSignallingMethod = 1;          % Determines what differentiation does to the system: 0=no change to system, 1=signalling changes in the cell above or below (randomly selected), 2= same as 1 but the cells selected are always the differentiating cell and the one below it, 3= replace diff cell with population mean, 4=local/neighbour mean replacement, 5=replace with random neighbour level,
diffGamma=3;                       % Strength of changed signalling - diffGamma multiplies the Hes5 translation rate in the cells affected by diffSignallingMethod
diffLengthHours=7;                 % Time spent in altered coupling strength (hours)
diffLength=diffLengthHours*60/dt;  % Time spent in altered coupling strength (simulation time steps)

recievingGamma=diffGamma;          % The gamma that a recieving cell experiences if in contact with a differentiating cell
nominalGamma=1;                    % Value of gamma when cells are not differentiating
lowerCoupling=1;                   % Keep as 1

AvM       = 0;                     % Absolute vs moving mean threshold (0=Absolute thresh, 1=Moving mean thresh)
wl_mins   = 100;                   % Window length in mins to take the moving mean from
wl        = wl_mins/dt;            % Window length converted to number of vector elements
S         = 0.003;                 % Rate of differentiation. Nominal value of 0.003

DiffTime_hours=50;                 % Time at which differention can start to occur to allow the system to approach a dynamic equilibrium (hours)
DiffTime=DiffTime_hours*60/dt;     % Diff time in simulation time steps


%Protrusions - distal cell parameters
vertProtStrength=1.5;              % Distal signalling efficiency in the vertical protrusions
numDistalNeigh=4;                  % 4 distal neighbours uses one cell directly above and one below, one  directly left and and one directly right. 6 = one directly above and below, and then two offset vertically from the horizontal on the left and two offset on the right. 
noProximal = 0;                    % 0 means that proximal cells will be included
vertProtrusions=2;                 % Signalling distance of the vertical protrusions (2 means 2 cell diameters away)
horzProtStrength=vertProtStrength; % Distal signalling efficiency in the horizontal protrusions
HorzProtrusions=vertProtrusions;   % Signalling distance of the horizontal protrusions (2 means 2 cell diameters away)
StrengthGrad=0;                    %Gradient in protrusion interactions with distance (leave as 0)         


%Cell-movement/swapping parameters
ImplementSwapping   = 0;   % Cells will randomly swap in the horizontal direction
Pm=0.005;                  % Total probability of movement either left or right (max Pm value is 0.5) Nominal value 0.005
SwapThresh=1-Pm/2;         % Treshold to be used in probability of swapping function within diffSolver.m
Ts=5;                      % Time (in mins) between each swapping event (nominal Ts=5)       

if cols==1
    kymoCols=1;            % Used later for plotting kymographs        
else
    kymoCols=2;
end

%% Select simulation outputs (0 = don't plot, 1 = show plot)

%Various frequency analysis options
TemporalFourier     = 1;   %Gives average Fourier period of the cell population
TemporalWavelet     = 0;   %Preliminary implementation of wavelet to examine whether Hes5 switches between periodic and noisy
SpatialFourier      = 1;   %Detect significant spatial periodic expression of Hes5
PerformMFV          = 1;   % Calculate mean field value (degree of synchronisation)

%Visualising HES5 protein levels
AnimateGrid         = 1;   % Show hexagonal grid with HES5 expression animated
AnimationSpeed      = 0.9; % Multiplier for speed of animation
AnimationSubplots   = 0;   % 0=Just Hes levels plot, 1=show swapping, 2=show crude differentiation, 3=Wave/single column visualisation
ShowRandomCells     = 1;   % Time traces of individual cells
ShowLastFrame       = 0;   % Shows the last time point in the hexagonal lattice arrangement of cells
AnimatePhaseSpace   = 0;   % Animate Hes protein levels vs mRNA levels over time
plotSwitchTime      = 1;   

%Saving animation options
saveAnimation   = 0;       % 1=yes, 0=no. Any animations that run will be made into MPEG4 files.
filename1 = 'Video1';      % Specify the output file name

%% Produce neighbour matrix

if cells==1
    boundary=0;    % Other boundaries don't make sense for 1 cell!
end

[NM,NumNeigh]=neighbours(rows, cols, boundary, vertProtrusions, HorzProtrusions, vertProtStrength, horzProtStrength, StrengthGrad, numDistalNeigh, noProximal); % Defines a matrix of signalling neighbours
NM=sparse(NM);     % Making the Neighbour Matrix (NM) sparse saves a lot on computational cost for large grid sizes
NM=NM./max(NM(:)); 

eps=(1./NumNeigh); % eps defines the contribution of the neighbouring cells
eps=eps*(1/6)/min(eps);
if cells==1
    eps=0;
end

%% Model parameters
%Select a parameter set from the Bayesian inferred single cell parameters 
%from (Manning et al. 2019
%https://www.nature.com/articles/s41467-019-10734-8)

% Columns in accepted_parameters: 1:a_m, 2:a_p, 3:P_H0, 4:TauH, 5:n_H
load accepted_parameters %Selcted for coherence values between 0.05 and 0.015
load summary_stats

%Selecting parameter sets that give the following single cell summary
%stats:
Select(:,1)=accepted_parameters(:,4)>15;                        %TauH selection
Select(:,2)=accepted_parameters(:,5)<3;                         %nH selection
Select(:,3)=summary_stats(:,4)<0.13 & summary_stats(:,4)>0.12;  %Coherence selection
Select(:,4)=summary_stats(:,1)<5e4 & summary_stats(:,1)>3.5e4;  %Mean selection
Select(:,5)=summary_stats(:,2)>0.05 & summary_stats(:,2)<0.1;   %Coefficient of variation                         %Standard deviation selection
Select(:,6)=summary_stats(:,3)<10*60;                           %Peak period oscillation selection 

accept=sum(Select,2);
accept=find(accept==size(Select,2));
fprintf(sprintf('Number of accepted parameter sets = %.f \n',numel(accept)))

PlotSummaryStats=0;
if PlotSummaryStats==1
    
    figure(21)
    clf
    
    subplot(231); histogram(accepted_parameters(accept,4)); title('TauH')
    subplot(232); histogram(accepted_parameters(accept,5)); title('n_H')
    subplot(233); histogram(summary_stats(accept,4));       title('Coherence')
    subplot(234); histogram(summary_stats(accept,1));       title('Mean')
    subplot(235); histogram(summary_stats(accept,2));       title('Coefficient of variation')
    subplot(236); histogram(summary_stats(accept,3)/60);    title('Period')

    figure(22)
    clf
    plot(summary_stats(:,3)/60,summary_stats(:,4),'.','color','w');hold on
    plot(summary_stats(accept,3)/60,summary_stats(accept,4),'.','color',[0.2 0.9 0.4])
    xlabel('Period (hours)')
    ylabel('Coherence')
end

% Choose single cell parameter set from single cell parameterisation
AP=accepted_parameters;
SS=summary_stats;

% s=197; %0.05 Coherence
% s=663; %0.05 Coherence
% s=987; %0.05 Coherence

% s=300; %0.1 Coherence
% s=350; %0.1 Coherence
% s=371; %0.1 Coherence

s=3700; %0.15 Coherence %%%%Nominal param set%%%%
% s=616; %0.15 Coherence
% s=713; %0.15 Coherence

%__________________________________________________________________________
%Single-cell parameters
a_m   = AP(s,1);   % Transcription rate (min^-1)
a_p   = AP(s,2);   % Translation rate (min^-1)
P_H0  = AP(s,3);   % HES5 autorepression repression threshold (# molecules)
TauH  = AP(s,4);   % HES5 autorepression time delay (min)
n_H   = AP(s,5);   % HES5 autorepression Hill coefficient    
u_m   = log(2)/30; % mRNA half-life of 30min
u_p   = log(2)/90; % Protein half-life of 90min

%Multicellular parameters (specifies the lateral inhibition Hill function parameters)
P_ND0 = 3800;      % Lateral inhibition repression threshold (# molecules)
n_ND  = 3;         % Lateral inhibition Hill coefficient
TauND = 0;         % Lateral inhibition time delay (min)

%Additional specification of the min and max values in the Hill function
%(leave these alone!)
Up=1;
Low=0;
a=Low;
b=Up-a;
%__________________________________________________________________________

gamma=1; % Maximum intercellular Hill function value
% gammaRand=rndrng(cells,1,0,1);

TauH_step=round(TauH/dt);   % Conversion to simulation time steps
TauND_step=round(TauND/dt); % Conversion to simulation time steps

if coupledCells==0
    eps=eps.*0; %Set eps=0 to decouple the cells 
end

if CrudeDiff==0
    DiffTime=Nt*2; % Ensure the differentiation algorithm condition never triggers if DBP is turned off
end


%% Delay differential equations (using anonymous functions (@ symbol) rather than standard functions)
dm=@(m,p,p_delay,Psum_delay,gamma,eps) a_m*1./(1+(p_delay./P_H0).^n_H).*gamma.*(  a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND)  ) - u_m.*m; % Describe mRNA change in time                                                                           

if turnOffAutorepression==1
    dm=@(m,p,p_delay,Psum_delay,gamma,eps) a_m.*gamma.* (a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND)) - u_m.*m; % No self repression version for simple Notch Delta modelling
end

dp=@(m,p,p_delay,Psum_delay,gamma,eps) a_p.*m - u_p.*p; % Describes protein change in time

%% Stochastic delay differential equations
dm1=dm;

dm2=@(m,p,p_delay,Psum_delay,gamma,eps) sqrt(a_m.*gamma./(1+(p_delay./P_H0).^n_H).*( a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND) ) + u_m.*m);

if turnOffAutorepression==1
   dm2=@(m,p,p_delay,Psum_delay,gamma,eps) sqrt(a_m.*gamma.*( a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND) ) + u_m.*m);
end

dp1=dp;
dp2=@(m,p,p_delay,Psum_delay,gamma,eps) sqrt(a_p.*m + u_p.*p);

%Initialise vectors and initial values (random)____________________________

% rng(1); % Same seed for random number generator (prevents initial heterogeneities from changing on every simulation, so comment out this line if you want different initial conditions each time)
P=[rndrng(cells,1,0.2*P_H0,1.8*P_H0),zeros(cells,Nt-1),zeros(cells,1)]; %Vector that will store protein values from the simulation
M=[rndrng(cells,1,0,20),zeros(cells,Nt-1),zeros(cells,1)];              %Vector that will store mRNA values from the simulation


%% Main loop (solves differential equations)
PercentCounter=1; %Print progress of the differential solver to the command line (1=yes, 0=no)
[P, M, t_step, DiffYN, DiffYNflash, CT, DiffThresh, M, AbsThresh, GAMMA,countArray,EPS]=diffSolver(P, M, Nt, TauND_step, TauH_step, NM, gamma, dm, dp, dm1, dm2, dp1, dp2, dt, stochastic, rows, cols, DiffTime, S, ImplementSwapping, SwapThresh, Ts, PercentCounter, AvM, diffSignallingMethod, wl, eps, diffLength, diffGamma, recievingGamma, nominalGamma);

%%
% Following code is for visualising differentiation events in the animated output
longExposure=DiffYNflash;       
longExposureFade=longExposure;
DiffElems=find(DiffYNflash==1);
ExpTimeHr=15;
ExposureTime=4*60/dt;
for E=1:ExposureTime
    Elems=DiffElems+cells*E;
    Elems(Elems>numel(longExposure))=[];
    longExposure(Elems)=1;
    longExposureFade(Elems)=1-E/ExposureTime;
end

if eps==0
    fprintf('\nCells are uncoupled!\n')
end
if turnOffAutorepression==1
    fprintf('\nAutorepression is turned off!\n')
end

%==========================================================================
%%                        Animation of Grid
%==========================================================================

if AnimateGrid==1
    n=1; %Length of hexagon side
    [X_vertices, Y_vertices]=hex(rows,cols,n); % Returns hexagonal grid vertices
    
    colour_index1=reshape(flipud(vecTOmat(P(:,1),cols)),[1,cols*rows]);
    
    figure(1)
    clf;
    fig=gcf;
    fig.InvertHardcopy = 'on';

%     map=redblue(1000);
%     map=pmkmp(100,'Swtth');
    map=viridis(1000);
%     map=gray;
    if AnimationSubplots >0
        
        hex1=subplot(121);
        colormap(hex1,map); 
    else
        colormap(map);
    end
    
    if AnimationSubplots == 3
        hex1=subplot(121);
    end
    
    hexagons1 = patch(X_vertices,Y_vertices,colour_index1,'edgecolor','none');
    set(gca,'xtick',[],'ytick',[]); 
  
    
    if AnimationSubplots~=3
        colorbar; 
    end
    axis equal; title('HES5 expression')
    set(gca,'Visible','off')
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'fontsize',15)
    caxis([min(min(P(:,round(0.5*Nt):end))) 0.7*max(max(P(:,round(0.5*Nt):end)))])
    xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
    ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
    
    
    if AnimationSubplots==1
        colour_index2=reshape(flipud(vecTOmat(CT(:,1),cols)),[1,cols*rows]);
        
        hex2=subplot(122);
        hexagons2 = patch(X_vertices,Y_vertices,colour_index2,'edgecolor','none');
        set(gca,'xtick',[],'ytick',[]); 
        colormap(hex2,randColourMapSparse(cells,0.1));
        colorbar; axis equal; title('Cell tracking')
        set(gca,'Visible','off')
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        set(gca,'fontsize',10)
        xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
        ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
        
    elseif AnimationSubplots==2
        exposureColorMap=gray;
        
        if graphAppearance==2
            exposureColorMap=exposureColorMap(60:end,:);
        end
        
        P_exposure=GAMMA(:,1:size(longExposureFade,2));
        colour_index3=reshape(flipud(vecTOmat(P_exposure(:,1),cols)),[1,cols*rows]);
    
        hex3=subplot(122);
        hexagons3 = patch(X_vertices,Y_vertices,colour_index3,'edgecolor','none');
        set(gca,'xtick',[],'ytick',[]); 
        colorbar
        colormap(hex3,exposureColorMap);
        axis equal; title('Crude diff')
        set(gca,'Visible','off')
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        set(gca,'fontsize',15)
        caxis([1 diffGamma])
        xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
        ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
        
        elseif AnimationSubplots==3
            subplot(122)
            xmin=min(min(P(:,Nt/2:Nt)));
            xmax=max(max(P(:,Nt/2:Nt)));
            p1=plot(flipud(P((cols/2-1)*rows+1:cols/2*rows,1)),1:rows,'k'); hold on
            p2=plot(flipud(P((cols/2-1)*rows+1:cols/2*rows,1)),1:rows,'k.'); hold off
            title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));
            ylabel('Cell number')
            xlabel('Hes5 abundance')
            xlim([xmin xmax])
            ylim([1 rows])

    end
      
    startT = Nt/2;
    T_STEP=startT:round(AnimationSpeed*Nt/400):Nt;
    ti=length(T_STEP);
    im=struct([]);
   
    for idx=1:ti
        
        t_step=T_STEP(idx);
        
        colour_index1=reshape(flipud(vecTOmat(P(:,t_step),cols)),[1,cols*rows])';
        set(hexagons1, 'FaceVertexCData',colour_index1); 
        title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));
        
        if AnimationSubplots==2
            subplot(121)
        end
        
        if AnimationSubplots==1
            colour_index2=reshape(flipud(vecTOmat(CT(:,t_step),cols)),[1,cols*rows])';
            set(hexagons2, 'FaceVertexCData',colour_index2); 
        elseif AnimationSubplots==2
            colour_index3=reshape(flipud(vecTOmat(P_exposure(:,t_step),cols)),[1,cols*rows])';
            set(hexagons3, 'FaceVertexCData',colour_index3); 
            subplot(122)
            title('Differentiation events');
        elseif AnimationSubplots==3
            set(p1, 'xData',flipud(P((cols/2-1)*rows+1:cols/2*rows,t_step)));
            set(p2, 'xData',flipud(P((cols/2-1)*rows+1:cols/2*rows,t_step)));
            subplot(122)
            title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));
        end
        drawnow;  

        if saveAnimation==1
            F1=getframe(gcf); %gca makes frame without labels, gcf makes frame including labels
            im{idx}=frame2im(F1);
            F_all{idx}=F1;
        end

    end
    
    if saveAnimation==1
        IDX=idx;
        
        video = VideoWriter(filename1,'MPEG-4'); %create the video object
        video.FrameRate = 30;
        open(video);                             %open the file for writing
        for idx = 1:IDX
            writeVideo(video,F_all{idx});        %write the image to file
        end 
        close(video);                            %close the file
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Calculate population differentiation rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CrudeDiff==1

    %%
    rand_ind=ceil(cells*rand(1,ceil(cells*0.1)));
    
    if cells==2
        rand_ind=[1 2];
    end
    
    T_plot=T(DiffTime+1:end);
    cell_num=1;
    P_plot=P(cell_num,DiffTime+1:end);
    DT=DiffThresh(cell_num,end);
    C=P_plot-DT;
    C(C>0)=0;
    C=C/min(C);
    
    load 'cmap.mat'
    
    figure(102)
    clf
    set(gcf,'renderer','Painters')
    
    plot(T_plot./60,DiffThresh(cell_num,DiffTime+1:end),'--','color','k','linewidth',1.5);hold on
    scatter(T_plot./60,P_plot,12,C,'filled');hold on
    xlim([DiffTime_hours tf_hours])
    ylabel('HES5 abundance')
    xlabel('Time (hours)')
    colormap(cmap)
    h=colorbar;
    ylim([0.95*min(P_plot) 1.05*max(P_plot)])
    set(gca,'FontSize',12)
    legend('Differentiation threshold')
    title('Probability of differentiation')
    
    h=colorbar;
    set(h,'XTickLabel',{'Zero probability of diff', 'Higher probability'}); 
    set(h,'XTick',[0 1]); 
   
    NumberDiff=sum(DiffYN);

    %Estimate half life
    NumberDiffSub=NumberDiff-cells/2;
    [~,I]=min(abs(NumberDiffSub));
    t_half=T(I);

    t_half=t_half(1)-DiffTime; %Half life in minutes
    lam=log(2)/t_half;

    Td=T(DiffTime/dt+1:end);

    D=cells.*(1-exp(-lam.*(Td-DiffTime)));

    DiffSignal=DiffYNflash(:,DiffTime/dt:end);
    WindowSize_Minutes=150;
    WindowSize=WindowSize_Minutes/dt;
    
    fprintf('\nDifferentiation rate over whole time (%% of cell population) = %.1f %%/hour', 100*sum(DiffSignal(:))/((tf_hours-DiffTime_hours)*cells));
    fprintf('\nDifferentiation rate over whole time (absolute number of cells) = %.1f cells/hour', sum(DiffSignal(:))/((tf_hours-DiffTime_hours)));

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Final time point of hexagonal grid plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ShowLastFrame==1
    figure(103)
    set(gcf,'renderer','painters') %For EPS file export
    clf;
    
    fig=gcf;
    if graphAppearance==1
        fig.InvertHardcopy = 'off';
    else
        fig.InvertHardcopy = 'on';
    end
    
    n=1; %Length of hexagon side
    [X_vertices, Y_vertices]=hex(rows,cols,n); % Returns hexagonal grid vertices
    colour_index=reshape(flipud(vecTOmat(P(:,t_step),cols)),[1,cols*rows]);
    map=viridis(1000);
    
    hexagons = patch(X_vertices,Y_vertices,colour_index,'edgecolor','none');
    colormap(map); 
    colorbar
    axis equal;
    set(gca,'Visible','off')
%     caxis([min(min(P(:,round(0.2*Nt):end))) max(max(P(:,round(0.2*Nt):end)))])
    colour_index=reshape(flipud(vecTOmat(P(:,end),cols)),[1,cols*rows])';
    set(hexagons, 'FaceVertexCData',colour_index);
%     title(sprintf('Time: %.0f hours', tf/60));
title({})
    xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
    ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'FontSize',12)

    drawnow
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Plots of HES5 levels in random cells 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ShowRandomCells==1
%     tStart=(tf_hours-50)*60/dt;
%     tStart=Nt/2;
    tStart=1;
    
    polyOrder=1;
    smoothP2=sgolayfilt(P,polyOrder,165,[],2);
    pop=P(:,Nt/2:end); %Population of cells to analyse
    meanPop=mean(pop,'all'); %Mean population expression value
    
    figure(2) 
    set(gcf,'renderer','Painters')
    fig=gcf;
    if graphAppearance==1
        fig.InvertHardcopy = 'off';
    else
        fig.InvertHardcopy = 'on';
    end
    clf
    
    %How many cells to plot
    desired_cells=7;
    if desired_cells>cells
        desired_cells=cells;
    end
    rand_ind=randperm(cells,desired_cells);
    
    num=length(rand_ind);
    subplot(1,4,[1 2 3])
    h=plot(T(tStart:end)/60,P(rand_ind,tStart:end),'linewidth',1); hold on
    plot(T((DiffTime+1)/dt:end)/60,DiffThresh(1,(DiffTime+1)/dt:end),'--','color',[0.5 0.5 0.5],'linewidth',1.5); hold on
    set(h, {'color'}, num2cell(randColourMap(num),2));
    xlabel('Time (hours)')
    ylabel(sprintf('Hes protein count'))
    xlim([tStart*dt/60 tf_hours])
    set(gca,'FontSize',15)
    title('Subset of cell Hes5 time traces')
    ylim([0 1.01*max(max(P(:,Nt/2:Nt)))])
    
    subplot(1,4,4)
    histogram(P(:,Nt/2:Nt),'normalization','probability','edgecolor','none','facecolor',[0.5 0.5 0.5])
    xlim([0 1.01*max(max(P(:,Nt/2:Nt)))])
    ylabel('Density')
    title(sprintf('Median=%.f',median(P(:,Nt/2:Nt),'all')))
    set(gca,'view',[90 -90])
    
    
    %Seperate trajectories on subplots
    numPlots=24;
    ymin=min(min(P(:,Nt/2:Nt)));
    ymax=max(max(P(:,Nt/2:Nt)));
    if numPlots>cells
        numPlots=cells;
    end
    rand_ind=randperm(cells,numPlots);
%     rand_ind=ceil(numPlots*rand(1,numPlots));
    t_start=1;
    t_start=100*60/dt;
    t_end=120*60/dt;

%     figure(8)
%     clf
%     sgtitle('Individual plots of HES5 abundance')
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
%     
%     
%     for n=1:numPlots
%         
%         subplot(4,6,n) 
% %         plot(T(t_start:end)/60,smoothP(rand_ind(n),t_start:end),'linewidth',0.5,'color',[1 1 1]); hold on
% %         plot(T(t_start:end)/60,smoothP2(rand_ind(n),t_start:end),'linewidth',0.5,'color',[1 1 1]); hold on
%         h11=plot(T(t_start:end)/60,P(rand_ind(n),t_start:end),'linewidth',1); hold on
%         plot([t_start/60,tf_hours],meanPop*ones(1,2),'k--')
%         
% %         plot(T/60,DiffThresh(1,:),'w--','linewidth',1.5)
%         set(h11, {'color'}, num2cell(randColourMapBright(1,[0.2 0.5 0.2]),2));
%         xlabel('Time (hours)')
%         ylabel(sprintf('Hes protein'))
%         set(gca,'FontSize',8)
%         xlim([100,115]);
% %         ylim([ymin ymax]);
%         ylim([0 ymax]);
%     end
%     drawnow
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           Switching Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotSwitchTime==1 && CrudeDiff == 1
polyOrder=1;
smoothP2=sgolayfilt(P,polyOrder,165,[],2); %Set to match params used on data (polyorder=1)
tStart_h=100;
tStart = 100*60/dt;
pop=smoothP2(:,tStart:end); %Population of cells to analyse

% %Convert to grid/array format
% lookUp=reshape(1:cells,rows,cols);
% %Extract every other row/column and average
% count=0;
% for r=1:2:rows
%     for c=1:2:cols
%         count=count+1;
%         idx=lookUp(r:r+1,c:c+1);
%         idx=idx(:);
%         popClust(count,:)=mean(pop(idx,:),1);
%     end
% end
    
% pop=popClust;

tNew=T(tStart:end);
popPlot=pop;

if AvM==1
popDiffThresh=DiffThresh(:,tStart:end);
pop=pop-popDiffThresh;
elseif AvM==0
meanPop=mean(pop,'all'); %Mean population expression value
pop=pop-meanPop;
end
pop(pop<0)=-1;
pop(pop>0)=1;
% imagesc(pop)
hi=[];
lo=[];
switchEvents=pop*0;
diffEvents=DiffYNflash(:,tStart:end);
for cell=1:size(pop,1)
    
    [B, N, IB]=RunLength_M(pop(cell,:));%N is the run length, B is the info on directionality, IB gives the element info on where the switch occurs
    
    B=B(2:end-1); %Removes clipped high and low regions at the start and end of signal
    N=N(2:end-1); %Removes clipped high and low regions at the start and end of signal
    IB=IB(2:end-1);
    hi=[hi; N(B==1)];
    lo=[lo; N(B==-1)];
    
    switchEvents(cell,IB(B==1))=1;
    switchEvents(cell,IB(B==-1))=1;
    
end

longExposure=diffEvents;
DiffElems=find(diffEvents==1);
ExpTimeHr=2;
ExposureTime=ExpTimeHr*60/dt;
for E=1:ExposureTime
    Elems=DiffElems+cells*E;
    ElemsNeg=DiffElems-cells*E;
    Elems(Elems>numel(longExposure))=[];
    ElemsNeg(ElemsNeg<1)=[];
    longExposure(Elems)=1;
    longExposure(ElemsNeg)=1;
    
end

switchRateNoProcessing=sum(abs(switchEvents),'all')/(size(switchEvents,2)*dt/60); %Differentiation events per hour

diffSwitchEvent=find(switchEvents+longExposure==2);
switchEventsPlot=switchEvents;
switchEventsPlot(diffSwitchEvent)=2;
switchEvents(diffSwitchEvent)=0;

switchRate=sum(abs(switchEvents),'all')/(size(switchEvents,2)*dt/60); %Differentiation events per hour
fprintf(sprintf('\n Switch rate with not accounting for differentiation reset in single cells = %.2f high/low switch events per hour \n Switch rate with correction = %.2f high/low switch events per hour \n',switchRateNoProcessing, switchRate))
CoV=std(P(:,tStart:end),[],'all')/mean(P(:,tStart:end),'all');

% cell=cells/2;
cell=7;

if CrudeDiff==1
figure(329)
set(gcf,'renderer','Painters')
clf
adjustSize = length(GAMMA(cell, :)) - length(P(cell,:));
normalSignalling = find(GAMMA(cell,Nt/2:end-adjustSize)<diffGamma);
changedSignalling = P(cell,tStart:end);
changedSignalling(normalSignalling)=meanPop;

aboveMean = P(cell,tStart:end);
belowMean = P(cell,tStart:end);
aboveMean(aboveMean<meanPop) = meanPop;
belowMean(belowMean>meanPop) = meanPop;

area(tNew.*dt./60, aboveMean, meanPop, 'facecolor', 1/255*[255 234 236], 'edgecolor', 'none'); hold on; %Fill area above mean in pink
area(tNew.*dt./60, belowMean, meanPop, 'facecolor', 1/255*[219 239 255], 'edgecolor', 'none'); hold on; %Fill area below mean in blue
area(tNew.*dt./60, changedSignalling, meanPop, 'facecolor', [0 0 0], 'facealpha',0.1, 'edgecolor', 'none'); hold on; %Fill areas where signalling is changed

plot(tNew.*dt./60, P(cell,tStart:end),'k', 'linewidth', 0.6);hold on
plot(tNew.*dt./60, popPlot(cell,:), 'color', [1 0.2 0.2], 'linewidth', 1); hold on

pbaspect([1.13 1 1])

title('Switching events due to DBP')
set(gca, 'fontsize', 14)

if AvM==0
plot([tNew(1)*dt/60 tNew(end)*dt/60], [meanPop meanPop],'k--', 'linewidth', 0.8)
elseif AvM==1
plot(T((DiffTime+1)/dt:end)/60,DiffThresh(1,(DiffTime+1)/dt:end),'--','color',[0.5 0.5 0.5],'linewidth',1.5); hold on
end

tSelect1=tNew(abs(switchEvents(cell,:))==1);
switchSelect=tSelect1*0;
tSelect2=tNew(diffEvents(cell,:)==1);
diffSelect=tSelect2*0;
tSelect3=tNew(switchEventsPlot(cell,:)==2);
switchSelect2=tSelect3*0;

xlim([tNew(1)*dt/60 tNew(end)*dt/60])
xlim([100 tf_hours])
max1=P(cell,Nt/2:end);
max2=meanPop*1.2;
ylim([0 max([max1(:);max2])])
ylim([0 4.5e4])

legend('High state','Low state', 'Higher transcription rate due to DBP', 'HES5','Smoothed HES5', 'Mean population abundance ')
xlabel('Time (h)')
ylabel('Hes5 abundance')
end

numPlots=20;
if numPlots>cells
    numPlots=cells;
end

rand_ind=randperm(cells,numPlots);

% figure(332)
% clf
% for n=1:numPlots
%     cell=rand_ind(n);
%   
%     subplot(4,5,n) 
% 
%     plot(tNew.*dt./60, P(cell,tStart:end),'w');hold on
%     plot(tNew.*dt./60, popPlot(cell,:),'color',[1 0.2 0.2]); hold on
%     if AvM==0
%     plot([tNew(1)*dt/60 tNew(end)*dt/60], [meanPop meanPop],'k--')
%     elseif AvM==1
%     plot(T((DiffTime+1)/dt:end)/60,DiffThresh(1,(DiffTime+1)/dt:end),'--','color',[1 1 1],'linewidth',1.5); hold on
%     end
% 
%     tSelect1=tNew(abs(switchEvents(cell,:))==1);
%     switchSelect=tSelect1*0;
%     tSelect2=tNew(diffEvents(cell,:)==1);
%     diffSelect=tSelect2*0;
%     tSelect3=tNew(switchEventsPlot(cell,:)==2);
%     switchSelect2=tSelect3*0;
% 
%     plot(tSelect2*dt/60,diffSelect+meanPop,'*','markersize',10,'color',[0.2 1 0.6])
%     plot(tSelect1*dt/60,switchSelect+meanPop,'.','markersize',20,'color',[0.2 0.6 1])
%     plot(tSelect3*dt/60,switchSelect2+meanPop,'o','markersize',7,'color',[0.3 0.7 1])
%   
%     
%     xlim([tNew(1)*dt/60 tNew(end)*dt/60])
%     max1=P(cell,tStart:end);
%     max2=meanPop*1.2;
%     ylim([0 max([max1(:);max2])])
%     xlabel('Time (h)')
%     ylabel('Hes5 abundance')
%     xlim([100 115])
%     
% end
% sgtitle('Differentiation and switching events')
% drawnow
end

%% Switching Time

if plotSwitchTime==1
polyOrder=1;
smoothP2=sgolayfilt(P,polyOrder,165,[],2);
startT=Nt/2;
Tm=(Nt-startT)*dt/60; %Tm=measurement time
pop=smoothP2(:,startT:end); %Population of cells to analyse
meanPop=mean(pop,'all'); %Mean population expression value

pop=pop-meanPop;
popPlot=pop;
pop(pop<0)=-1;
pop(pop>0)=1;
hi=[];
lo=[];

for i=1:cells
    [B, N, IB]=RunLength_M(pop(i,:));
    hi=[hi; N(B==1)];
    lo=[lo; N(B==-1)];
    
    alphaHi=sum(N(B==1))/(size(pop,2)); %Proportion of the time spent in high state
    alphaLo=sum(N(B==-1))/(size(pop,2)); %Proportion of the time spent in low state
    DcCells(i)=2*(1-max([alphaHi alphaLo]));
end

figure(444)
clf
histo=histogram(DcCells,'normalization','probability');
set(histo,'edgecolor','w');
set(histo,'facecolor',1/255*[50 50 50]);
title(sprintf('Mean D_C = %.2f', mean(DcCells)))
xlabel('Dynamicity')
ylabel('Count (normalised)')
xlim([0 1])
set(gca, 'fontsize', 15)


hi=hi*dt/60;
lo=lo*dt/60;
hi_all=hi;
lo_all=lo;
hi(hi>16)=[];
lo(lo>16)=[];

hiX=hi*0+1;
loX=lo*0+2;


% figure(333)
% clf
% subplot(211)
% histo=histogram(hi,'normalization','probability');
% set(histo,'edgecolor','w');
% set(histo,'facecolor',[0.9 1 0.7]);
% xlabel('Time (h)')
% ylabel('Occurence')
% title(sprintf('Time spent in high state (mean=%.1fh)',mean(hi)))
% % title('Time spent in high state')
% set(gca,'fontsize',12)
% subplot(212)
% histo=histogram(lo,'normalization','probability');
% set(histo,'edgecolor','w');
% set(histo,'facecolor',[1 0.6 0.7]);
% xlabel('Time (h)')
% ylabel('Occurence')
% title(sprintf('Time spent in low state (mean=%.1fh)',mean(lo)))
% set(gca,'fontsize',12)
%%
% figure(334)
% set(gcf,'renderer','Painters')
% clf
% subplot(211)
% histo=histogram(hi_all,'normalization','probability');
% % set(histo,'edgecolor','w');
% % set(histo,'facecolor',[0.9 1 0.7]);
% set(histo,'facecolor', [0 0 0], 'FaceAlpha', 0.6);
% xlabel('Time (h)')
% ylabel('Frequency')
% % xlim([0 170])
% title(sprintf('Time spent in high state (mean=%.1fh)',mean(hi_all)))
% % title('Time spent in high state')
% set(gca,'fontsize',12)
% 
% subplot(212)
% histo=histogram(lo_all,'normalization','probability');
% % set(histo,'edgecolor','w');
% % set(histo,'facecolor',[1 0.6 0.7]);
% set(histo,'facecolor',[0 0 0], 'FaceAlpha', 0.8);
% xlabel('Time (h)')
% ylabel('Frequency')
% title(sprintf('Time spent in low state (mean=%.1fh)',mean(lo_all)))
% set(gca,'fontsize',12)
% % xlim([0 170])

figure(335)
set(gcf,'renderer','Painters')
clf
% subplot(211)
histo=histogram(lo_all,'normalization','probability'); hold on
% set(histo,'facecolor', 1/255*[219 239 255], 'edgecolor', 'none');
set(histo,'facecolor', 1/255*[160 200 255], 'edgecolor', 'none');
histo=histogram(hi_all,'normalization','probability'); hold on
% set(histo,'facecolor', 1/255*[255 234 236], 'edgecolor', 'none');
set(histo,'facecolor', 1/255*[255 170 180], 'edgecolor', 'none');
plot(mean(lo_all), 0, 'd','color', 'b', 'markerfacecolor', 'b', 'markersize', 8)
plot(mean(hi_all), 0, 'd','color', 'r', 'markerfacecolor', 'r', 'markersize', 8)

% set(histo,'edgecolor','w');
% set(histo,'facecolor',[0.9 1 0.7]);
% set(histo,'facecolor', [0 0 0], 'FaceAlpha', 0.6);
xlabel('Persistence time (h)')
ylabel('Normalised count')
% xlim([0 170])
% title(sprintf('Time spent in high state (mean=%.1fh)',mean(hi_all)))
% title('Time spent in high state')
set(gca,'fontsize',14)
pbaspect([1.13 1 1])
% xlim([0 170])

legend('Persistence time low','Persistence time high', sprintf('Mean low = %.1fh',mean(lo_all)), sprintf('Mean high = %.1fh',mean(hi_all)))

end



%==========================================================================
%%                        Temporal Frequency Analysis
%==========================================================================
    
  
    %Changes fraction of signal to use (deterministic uses earlier part of
    %the signal as it is a damped oscillator)
    if stochastic==1
        frac=0.2;
    else   
        frac=0.01;
    end
%     Osc=find(max(P(:,0.7*Nt:Nt),[],2)>1000);
    Osc=find(max(P(:,0.7*Nt:Nt),[],2)>0);
    Y_raw=P(Osc,:);
     
    %Detrending
    polyOrder=2;                  %Order of detrending polynomial in detrend.m function
    frameTime=40;                 %Frame length in hours for detrending window
    frameLength=frameTime*60./dt; %Conversion to window length in elements  
%     frameLength=75;
    [Ydetrend,t,Ysmooth,f,P1,coherence,f_C1,P_C1,avgFourier,ind_per,I]=tempFourier(T,Y_raw,polyOrder,frac,frameLength,Nt,dt);
    
if TemporalFourier==1   
    
    fprintf('\nExpected coherence from single cell posterior data = %.2f\n', SS(s,4))    
    fprintf('Coherence of oscillators in this model = %.2f\n',coherence)
    fprintf('\nExpected peak period from single cell summary stats = %.2f\n', summary_stats(s,3)/60)
    fprintf('Oscillators with periods below 7 hours = %.2f%%\n', 100*sum(ind_per<7)/numel(ind_per))
%__________________________________________________________________________    
%__________________________________Plots___________________________________

    figure(14)
    clf
    fig=gcf;
    fig.InvertHardcopy = 'off';
    txtSize=10;
    sampleCell=size(Y_raw,1);
    
    subplot(321)
    plot(t./60,Ysmooth(sampleCell,:),'--', 'color','k'); hold on
    plot(T(frac*Nt:end)./60,Y_raw(sampleCell,frac*Nt:end),'linewidth',2,'color', [0 0.6 0.6])
    title('Signal with detrend line')
    xlim([dt*Nt*frac/60 dt*Nt/60])
    xlabel('Time (hours)')
    set(gca,'Fontsize',txtSize)
    
    subplot(322)
    plot([t(1)/60 t(end)/60], [0 0],'k--'); hold on
    plot(t./60,Ydetrend(sampleCell,:),'color',[0 0.6 0.6],'linewidth',2); 
    xlim([t(1)./60 t(end)./60])
    title('Detrended signal')
    xlabel('Time (hours)')
    set(gca,'Fontsize',txtSize)
    
    subplot(3,2,3)
    h0=plot(60*f,P1);
    xlim([0 0.5])
    xlabel('Frequency (1/h)')
    set(h0, {'color'}, num2cell(randColourMapBright(size(P1,1),[1 0 0]),2));
    ylim([0 1.2*max(max(P1))])
    title('Individual cell Fourier transforms')
    set(gca,'Fontsize',txtSize)
    
    subplot(324)
    histo=histogram(ind_per,15,'normalization','probability');
    set(histo,'edgecolor','w');
    set(histo,'facecolor',[0.5 0.5 0.5]);
    xlabel('Dominant period (h)')
    ylabel('Num of cells')
    title('Individual cell dominant periods')
    set(gca,'fontsize',txtSize)
    
    subplot(3,2,[5 6])
    h1=area(60*f,avgFourier,'LineStyle','none'); hold on;
    h1.FaceColor=[0 0.6 0.6];
    h2=area(60*f_C1,P_C1,'LineStyle','none'); hold on;
    h2.FaceColor=[0.9 0.4 0.4];
    plot(60*f,avgFourier,'.','linewidth',1,'color','w'); 
%     plot(60*f(I-ten_percent),avgFourier(I-ten_percent),'ro','linewidth',2)
%     plot(60*f(I+ten_percent),avgFourier(I+ten_percent),'ro','linewidth',2)
    plot(60*f(I),avgFourier(I),'.','markersize',20,'color', [1 1 1])
    title('All cells Fourier transform')
    xlabel('Frequency (1/h)')
    xlim([0 0.5])
    ylim([0 1.2*max(avgFourier)])
    title(sprintf('Sum of Fourier Transforms   |   Dominant period: %.2f hours   |   Coherence: %.2f', 1./(60*f(I)), coherence))
    set(gca,'Fontsize',txtSize)
    drawnow
    
    
%__________________________________________________________________________
%Light mode report version

    numPlots=24;
    if numPlots>cells
        numPlots=cells;
    end
    
    rand_ind=randperm(cells,numPlots);
%     rand_ind=rand_ind(1:numPlots);
    
%     t_start=80*60/dt; %fractional start time in elements
    t_start = Nt - 50*60/dt; 
    t_start_h = t_start*dt/60; %fractional start time in hours

    figure(8)
    clf
    sgtitle('Individual cell HES5 plots with detected peak temporal period')
    fig=gcf;
    fig.InvertHardcopy = 'off';
    
    for n=1:numPlots
      
        subplot(4,6,n) 
        plot(t./60, Ysmooth(rand_ind(n),:),'k--'); hold on
        h12=plot(T(t_start:end)/60,Y_raw(rand_ind(n),t_start:end),'linewidth',1);
        set(h12, {'color'}, num2cell(randColourMapBright(1,[0.2 0.5 0.2]),2));
        xlabel('Time (hours)')
        ylabel(sprintf('Hes protein'))
        set(gca,'FontSize',8)
        xlim([t_start_h,tf_hours])
        title(sprintf('T=%.1fh',ind_per(rand_ind(n))))
        
    end
    drawnow
    
end

%==========================================================================
%%                        Synchronisation tests
%==========================================================================

if PerformMFV==1
    
    [PH,PH_unwrap,COP]=phase(t,Ydetrend); % Uses the Hilbert transform to extract phase from oscillatory signals
    Sync_val=mean(abs(COP(round(0.5*length(COP)):round(0.8*length(COP)))));
    fprintf(sprintf('\nMean field/Kuramoto order value = %.2f \n', Sync_val))
    
    mi=size(PH_unwrap,1);
    x_poly    = linspace(t0,tf,Nt+1);
    PH_smooth=zeros(mi,length(x_poly));
    
    for m=1:mi
        [p,S,mu]  = polyfit(t,PH_unwrap(m,:),8);
        PH_smooth(m,:) = polyval(p,x_poly,S,mu);
    end

    elem=1;
    PH1=PH_smooth(:,1:end-elem);
    PH2=PH_smooth(:,1+elem:end);
    inst_freq=((PH2-PH1)./(2*pi))./(elem*dt); %cycles/min    
end



%==========================================================================
%%                          Wavelet Transform
%==========================================================================

if TemporalWavelet==1
    
% wc value is the central frequency, typical values range from 0.5 to 4: 
%                  wc  =  0.5   1   2   3   4    
% Better time resolution <-------------------> Better frequency resolution

    wc=4; %Central frequency
    
    Ntime=length(T)/100; %Reduce time resolution to speed up computation
    Nw=100;
   
    w0h=0.01;       %Lower frequency bound in 1/hours
    wfh=1;          %Upper frequency bound in 1/hours
    w0h_c=w0h/wc;   %Correcting for wc scaling
    wfh_c=wfh/wc;   %Correcting for wc scaling
    w0=w0h_c/3600;  %lower frequency bound in Hz
    wf=wfh_c/3600;  %final frequency bound in Hz
%     dw=(wf-w0)/Nw;

    ts=T*60; %Time in seconds
    
    fh=0.8;
    fs=fh/3600;
    x=P(2,:);
    [x, ~, ~]=detrendSgolay(P(2,:), 1, 0, T, 165);
    
    [WT,W,t]=wavelet(x,Nw,w0,wf,wc,Ntime,ts);

    Wh=W*3600; %Frequency in 1/hours
    Th=t/3600; %Time in hours
   
    [~,I]=max(abs(WT).^2,[],1);
    domPer=1./Wh(I);
    domPer(I==1)=NaN;
    
    if wfh==1
        Wh(end)=1;
    end
    
    %plot
    figure(1212)
    clf
    subplot(311)
    plot(T./60,x,'k')
    ylabel('Hes protein number')
    xlabel('Time (hours)')

    subplot(312)
    load cmap
    surf(Th,Wh,abs(WT).^(2),'EdgeColor', 'None', 'facecolor', 'interp');
    xlabel('Time (hours)')
    ylabel('Frequency (1/hours)')
    xlim([t0 tf_hours])
    ylim([Wh(1) Wh(end)])
    view(2);
    set(gca,'yscale','log')
    colormap(gca,cmap)
    
    subplot(313)
    plot(Th,domPer,'k', 'linewidth', 2)
    xlabel('Time (hours)')
    ylabel('Dominant period (hours)')
    xlim([t0 tf_hours])
end

%==========================================================================
%%                        Spatial Frequency Analysis
%==========================================================================


if SpatialFourier==1 && cols==1
    startTime=tf_hours-100;
    time_frac=startTime/tf_hours;
    
    if time_frac>=1
        error('Increase length of simulation or decrease startTime for spatial frequency analysis.')
    end
    
    t=T(time_frac*Nt:end)./60;
    ti=time_frac*Nt;

    [Pnew]=improvedKymo(P,rows,cols,kymoCols); %Improved selection of kymograph
    
    kymoWidth=1; %Number of cells to use for kymograph selection width
    K=cols/kymoWidth;
    
    t_concat=t;
    Yraw=Pnew(:,ti:end);

    map=viridis(500);
    X=t(1:10:end);
    Y=1:rows;
    plotData=Yraw(:,1:10:end);
    Y_RAW=Yraw;
    
    if CrudeDiff==1
    
        figure(301)
        clf
        subplot(1,4,[1 2 3 ])
        j=imagesc(X,Y,plotData);
        xlim([X(1) tf_hours])
        ylabel('Row index')
        xlabel('Time (h)')
        set(gca,'YTickLabel',[]);
        set(gca,'FontSize',10)
        colormap(gca,map)
        colorbar('Location','westoutside');
        view(0,90)
        set(gca,'fontsize',10)
        
        
        subplot(1,4,4)
        histData=sum(DiffYNflash,2);
        histData= flipud(histData);
        bar(histData,1,'facecolor',1/256.*[80 80 80],'edgecolor','w')
        set(gca,'view',[90 -90])
        xlim([0.5 rows+0.5])
        set(gca,'XTickLabel',[]);
        set(gca,'fontsize',10)
        ylabel('Diff event count')
        
        diffSignal=DiffYNflash(:,DiffTime/dt:end);
        indivDiffFreq = sum(diffSignal,2)/(tf_hours-DiffTime_hours);
        
%         figure(302)
%         clf
%         histogram(indivDiffFreq,'facecolor',1/256.*[80 80 80], 'facealpha',1, 'edgecolor','w')
%         xlabel('Events/hour')
%         xlim([0 0.18])
%         title(sprintf('Diff = %.1f %%/h', 100*sum(DiffSignal(:))/((tf_hours-DiffTime_hours)*cells)))
    end
    

    %%

    %Split Y_RAW into x hour windows of expression and then average
    t_elem=(1-time_frac)*Nt; %-1 is for the fact the first column is not included in this analysis
    split_time=0.1;          %in hours (Nominal value of 2 hours)
    if split_time==0
        split_elem=1;
    else
        split_elem=round(split_time*60/dt);
    end
    
    
    M=1:split_elem:t_elem;
    M=round(M);
    II=length(M);
    
    split_t=linspace(0,tf_hours*K,II);
    
    Y_split=Y_RAW(:,M);
  
    t_split=t_concat(M);
    
    clipPos=t_split(2:end)-t_split(1:end-1);
    clipPos(clipPos<0)=-1;        %End positions marked by -1
    clipPos(clipPos>0)=0;         %Fill rest with zeros
    startPos=find(clipPos==-1)+1; %Find start positions of clips
    clipPos(startPos)=1;
    clipPos=[clipPos, -1];
    clipPos(1)=1;
    clipStart=find(clipPos==1);
    clipEnd=find(clipPos==-1);
    
    YDETREND2=Y_split-repmat(mean(Y_split,1),[rows 1]); %Detrending in spatial direction by subrtracing poulation mean at each time point
    
    Y=fft(YDETREND2,[],1);
    winL=length(Y(:,1));
    P2 = abs(Y/winL).^2;
    P1 = P2(1:winL/2+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:); 

    Fs=1;
    f = Fs*(0:(winL/2))/winL;
 
    KymLong=YDETREND2;
    
    [FisherPer,Occurrence]=Fisher(P1,f);
    
%     %Individual fourier plots
%     figure(60)
%     clf
%     rPlot=7;
%     cPlot=7;
%     yMax=max(P1(:));
%     numPlots=rPlot*cPlot;
%     rand_ind=ceil(numPlots*rand(1,numPlots));
%     for n=1:numPlots
%         
%         subplot(rPlot,cPlot,n)
%         plot(f,P1(:,rand_ind(n)),'k')
%         ylim([0 yMax])
%     end
%         
%     ymin=min(min(P(:,Nt/2:Nt)));
%     ymax=max(max(P(:,Nt/2:Nt)));
%     if numPlots>cells
%         numPlots=cells;
%     end
    
    %%
    figure(61) 
    clf
    set(gcf,'renderer','Painters')
    fig=gcf;
    fontSize=10;
    
    if graphAppearance==1
        fig.InvertHardcopy = 'off';
        
        subplot(2, 6, [1 2 3 4 5])
        load cmap
        colormap(gca,cmap)
        P1=[P1; P1(1,:)];
        f=[f, (2*f(end)-f(end-1))];

        t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));

        P1_lines=P1;
        P1_lines(:,clipEnd)=NaN;

        h=surf(t_plot,f,P1_lines);
        xlabel('Column over time')
        ylabel('Frequency (1/cell)')
        title('Periodogram over time for individual columns of cells')
        xlim([ t_plot(1) t_plot(end)])
        ylim([f(1) 0.5+f(end)-f(end-1)])
        view(0,90)
        colorbar
        set(h,'LineStyle','none')
        set(gca,'FontSize',fontSize)
        pbaspect([2.3 1 1])
        colormap(gca,cmap)

        centre=(t_plot(end)-t_plot(1))./((cols-1)*2); %cols -1 for improved kymo/removed first column
        xticks(linspace(t_plot(1)+centre,t_plot(end)-centre,(cols-1)))
        xticklabels({num2str([2:(cols)]')})

        subplot(266)
        [~,maxIdx]=max(mean(P1,2));
        plot(mean(P1,2),f,'w')
        title({'Mean FT. ',sprintf('Peak period = %.1f', mean(1./f(maxIdx)))})
        ylabel('Frequency (1/cells)')
        xlabel('Power')
        pbaspect([0.5 1.4 1])
        set(gca,'FontSize',fontSize)
        ylim([0 0.5])

        subplot(2,6,[7 8 9 10 11 12 ])
        plot(split_t,FisherPer,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:));hold on

        title(sprintf('Fisher G significant period detected %.0f%% of the time', 100*Occurrence))
        ylabel('Significant period (cells)')
        xlabel('Column over time')
        set(gca,'FontSize',fontSize)
        xlim([split_t(1) split_t(end)])
        ylim([0 5])

        centre=(split_t(end)-split_t(1))./((cols-1)*2);
        xticks(linspace(split_t(1)+centre,split_t(end)-centre,(cols-1)))
        xticklabels({num2str([2:(cols)]')})
        h=gca; h.XAxis.TickLength = [0 0];

        for n=1:numel(clipEnd)
            plot([split_t(clipEnd(n)), split_t(clipEnd(n))],[0,10],'w')
        end

        drawnow

    else

        fig.InvertHardcopy = 'on';
        subplot(2, 6, [1 2 3 4 5])
        load cmap
        colormap(gca,cmap)
        P1=[P1; P1(1,:)];
        f=[f, (2*f(end)-f(end-1))];

        t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));

        P1_lines=P1;
        P1_lines(:,clipEnd)=NaN;

        h=surf(t_plot,f,P1_lines);
        xlabel('Column')
        ylabel('Frequency (1/cell)')
        title('Periodogram over time for individual columns of cells')
        xlim([ t_plot(1) t_plot(end)])
        ylim([f(1) 0.5+f(end)-f(end-1)])
        view(0,90)
        colorbar
        set(h,'LineStyle','none')
        set(gca,'FontSize',fontSize)
        pbaspect([2.3 1 1])
        colormap(gca,cmap)

        centre=(t_plot(end)-t_plot(1))./((cols-1)*2); %cols -1 for improved kymo/removed first column
        xticks(linspace(t_plot(1)+centre,t_plot(end)-centre,(cols-1)))
        xticklabels({num2str([2:(cols)]')})

        subplot(266)
        [~,maxIdx]=max(mean(P1,2));
        plot(mean(P1,2),f,'k')
        title({'Mean FT. ',sprintf('Peak period = %.1f', mean(1./f(maxIdx)))})
        ylabel('Frequency (1/cells)')
        xlabel('Power')
        pbaspect([0.5 1.4 1])
        set(gca,'FontSize',fontSize)
        ylim([0 0.5])

        subplot(2,6,[7 8 9 10 11 12 ])
        plot(split_t,FisherPer,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:));hold on

        title(sprintf('Fisher G significant period detected %.0f%% of the time', 100*Occurrence))
        ylabel('Significant period (cells)')
        xlabel('Column')
        set(gca,'FontSize',fontSize)
        xlim([split_t(1) split_t(end)])
        ylim([0 5])

        centre=(split_t(end)-split_t(1))./((cols-1)*2);
        xticks(linspace(split_t(1)+centre,split_t(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})
        h=gca; h.XAxis.TickLength = [0 0];

        for n=1:numel(clipEnd)
            plot([split_t(clipEnd(n)), split_t(clipEnd(n))],[0,10],'k')
        end

        drawnow
    end 
    
    %% Subplots of all columns/kymos
    clipPos=t_split(2:end)-t_split(1:end-1);
    clipPos(clipPos<0)=-1; %End positions marked by -1
    clipPos(clipPos>0)=0;  %Fill rest with zeros
    startPos=find(clipPos==-1)+1; %Find start positions of clips
    clipPos(startPos)=1;
    clipPos=[clipPos, -1];
    clipPos(1)=1;
    
    clipStart=find(clipPos==1);
    clipEnd=find(clipPos==-1);
    
    plotCols=3;
    plotRows=cols/plotCols;
    plotCols=plotCols*2;
    %% 
    
    
    thresh=0.6;
    kymo=KymLong;

    for i=1:size(kymo,2)
        vect=kymo(:,i);
        [l(i),p(i)]=microClustDetect(vect,thresh);
    end
    clusterMean=nanmean(l);
    clusterOcc=1-sum(isnan(l))/numel(l);
        
end

vid=viridis(500);
if SpatialFourier==1 && cols>1
    startTime=tf_hours-100;
    time_frac=startTime/tf_hours;
    
    if time_frac>=1
        error('Increase length of simulation or decrease startTime for spatial frequency analysis.')
    end
    
    t=T(time_frac*Nt:end)./60;
    ti=time_frac*Nt;
    [Pnew]=improvedKymo(P,rows,cols,kymoCols); %Improved selection of kymograph
    
    kymoWidth=1; %Number of cells to use for kymograph selection width
    K=cols/kymoWidth;
    
    
    if mod(cols,2)==1
        error('Make cols even')
    end
    
    figure(23)
    clf 
    t_concat=t;
    
    for k=2:K %Starts at k=2 because the first column in not included in the analysis (to do with the improvedKymo function)

        Yraw(:,:,k)=Pnew((k-1)*rows+1:k*rows,ti:end);
        
        if k==2
            Y_RAW=Yraw(:,:,2);
            t_concat=t;
        else
            Y_RAW=[Y_RAW Yraw(:,:,k)];
            t_concat=[t_concat t];
        end
            
            
        subplot(3,4,k-1)
        
        X=t;
        Y=1:rows;
        j=surf(X,Y,squeeze(Yraw(:,:,k)));
        title(sprintf('Kymograph %.f',k))
        ylim([Y(1) Y(end)])
        xlim([t(1) t(end)])
        xlabel('Time (hrs)')
        ylabel('Row index')
        set(j, 'LineStyle','none')
        set(gca,'YTickLabel',[]);
        set(gca,'FontSize',10)
        colormap(gca,vid)
        colorbar
        view(0,90)
    end
    
%     %Single Kymo
%     
%     figure(24)
%     set(gcf,'renderer','Painters')
%     clf
%     X=t;
%     Y=1:rows;
%     j=surf(X,Y,squeeze(Yraw(:,:,1)));
% %     title(sprintf('Kymograph %.f',k))
%     ylim([Y(1) Y(end)])
%     xlim([t(1) t(end)])
%     xlabel('Time (hrs)')
%     ylabel('Cell index (row number)')
%     set(j, 'LineStyle','none')
%     set(gca,'YTickLabel',[]);
%     set(gca,'FontSize',12)
%     colormap(gca,vid)
%     colorbar
%     view(0,90)
   
    
    
    %Split Y_RAW into x hour windows of expression and then average
    t_elem=(1-time_frac)*Nt*(K-1); %-1 is for the fact the first column is not included in this analysis
    split_time=2; %in hours (Nominal value of 2 hours)
    if split_time==0
        split_elem=1;
    else
        split_elem=round(split_time*60/dt);
    end
    
    
    M=1:split_elem:t_elem;
    M=round(M);
    II=length(M);
    
    split_t=linspace(0,tf_hours*K,II);
    
    Y_split=Y_RAW(:,M);
  
    t_split=t_concat(M);
    
    clipPos=t_split(2:end)-t_split(1:end-1);
    clipPos(clipPos<0)=-1; %End positions marked by -1
    clipPos(clipPos>0)=0;  %Fill rest with zeros
    startPos=find(clipPos==-1)+1; %Find start positions of clips
    clipPos(startPos)=1;
    clipPos=[clipPos, -1];
    clipPos(1)=1;
    clipStart=find(clipPos==1);
    clipEnd=find(clipPos==-1);
    
    YDETREND2=Y_split-repmat(mean(Y_split,1),[rows 1]); %Detrending in spatial direction by subrtracing poulation mean at each time point
    
    
    Y=fft(YDETREND2,[],1);
    winL=length(Y(:,1));
    P2 = abs(Y/winL).^2;
    P1 = P2(1:winL/2+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);

    Fs=1;
    f = Fs*(0:(winL/2))/winL;
 
    KymLong=YDETREND2;
    
    [FisherPer,Occurrence]=Fisher(P1,f);
    
    
 
    %Individual fourier plots
    figure(60)
    clf
    rPlot=7;
    cPlot=7;
    yMax=max(P1(:));
    numPlots=rPlot*cPlot;
    rand_ind=ceil(numPlots*rand(1,numPlots));
    for n=1:numPlots
        
        subplot(rPlot,cPlot,n)
        plot(f,P1(:,rand_ind(n)),'w')
        ylim([0 yMax])
    end
%         
%     figure(59)
%     clf
%     rPlot=cols/4;
%     cPlot=4;
%     numPlots=rPlot*cPlot;
%     rand_ind=ceil(numPlots*rand(1,numPlots));
%     for n=1:numPlots
%         
%         subplot(rPlot,cPlot,n)
%         plot(1:rows,P((n-1)*rows+1:n*rows,rand_ind(n)),'w')
% %         ylim([0 yMax])
%     end
        
%     ymin=min(min(P(:,Nt/2:Nt)));
%     ymax=max(max(P(:,Nt/2:Nt)));
%     if numPlots>cells
%         numPlots=cells;
%     end
    
    %%
    figure(61) 
    clf
    set(gcf,'renderer','Painters')
    fig=gcf;
    fontSize=10;
    
    if graphAppearance==1
        fig.InvertHardcopy = 'off';
        subplot(2, 6, [1 2 3 4 5])
        load cmap
        colormap(gca,cmap)
    %     P1=[P1(1,:) ;P1];
    %     f=[f, (2*f(end)-f(end-1))];
        P1=[P1; P1(1,:)];
        f=[f, (2*f(end)-f(end-1))];

        t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));

        P1_lines=P1;
        P1_lines(:,clipEnd)=NaN;

        h=surf(t_plot,f,P1_lines);
    %     h=imagesc(t_plot,f,P1_lines);
        xlabel('Column over time')
        ylabel('Frequency (1/cell)')
        title('Periodogram over time for individual columns of cells')
        xlim([ t_plot(1) t_plot(end)])
        ylim([f(1) 0.5+f(end)-f(end-1)])
        view(0,90)
        colorbar
        set(h,'LineStyle','none')
        set(gca,'FontSize',fontSize)
        pbaspect([2.3 1 1])
        colormap(gca,cmap)

        centre=(t_plot(end)-t_plot(1))./((cols-1)*2); %cols -1 for improved kymo/removed first column
        xticks(linspace(t_plot(1)+centre,t_plot(end)-centre,(cols-1)))
        xticklabels({num2str([2:(cols)]')})

        subplot(266)
        [~,maxIdx]=max(mean(P1,2));
        plot(mean(P1,2),f,'w')
    %     plot(f,P1(:,1))
    %     title('Mean Fourier transform')
        title({'Mean FT. ',sprintf('Peak period = %.1f', mean(1./f(maxIdx)))})
        ylabel('Frequency (1/cells)')
        xlabel('Power')
        pbaspect([0.5 1.4 1])
        set(gca,'FontSize',fontSize)
        ylim([0 0.5])

        subplot(2,6,[7 8 9 10 11 12 ])
    %     [~,II]=max(P1,[],1);
    %     plot(split_t,FisherPer,'.','color',cmap(size(cmap,1)/2,:),'linewidth',2);hold on
        plot(split_t,FisherPer,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:));hold on

        title(sprintf('Fisher G significant period detected %.0f%% of the time', 100*Occurrence))
    %     title('Significant period (Fisher G)')
        ylabel('Significant period (cells)')
        xlabel('Column over time')
        set(gca,'FontSize',fontSize)
        xlim([split_t(1) split_t(end)])
        ylim([0 5])

        centre=(split_t(end)-split_t(1))./((cols-1)*2);
        xticks(linspace(split_t(1)+centre,split_t(end)-centre,(cols-1)))
        xticklabels({num2str([2:(cols)]')})
        h=gca; h.XAxis.TickLength = [0 0];

        for n=1:numel(clipEnd)
            plot([split_t(clipEnd(n)), split_t(clipEnd(n))],[0,10],'w')
        end

        drawnow

    else
        fig.InvertHardcopy = 'on';
        subplot(2, 6, [1 2 3 4 5])
        load cmap
        colormap(gca,cmap)
    %     P1=[P1(1,:) ;P1];
    %     f=[f, (2*f(end)-f(end-1))];
        P1=[P1; P1(1,:)];
        f=[f, (2*f(end)-f(end-1))];

        t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));

        P1_lines=P1;
        P1_lines(:,clipEnd)=NaN;

        h=surf(t_plot,f,P1_lines);
    %     h=imagesc(t_plot,f,P1_lines);
        xlabel('Column')
        ylabel('Frequency (1/cell)')
        title('Periodogram over time for individual columns of cells')
        xlim([ t_plot(1) t_plot(end)])
        ylim([f(1) 0.5+f(end)-f(end-1)])
        view(0,90)
        colorbar
        set(h,'LineStyle','none')
        set(gca,'FontSize',fontSize)
        pbaspect([2.3 1 1])
        colormap(gca,cmap)

        centre=(t_plot(end)-t_plot(1))./((cols-1)*2); %cols -1 for improved kymo/removed first column
        xticks(linspace(t_plot(1)+centre,t_plot(end)-centre,(cols-1)))
        xticklabels({num2str([2:(cols)]')})

        subplot(266)
        [~,maxIdx]=max(mean(P1,2));
        plot(mean(P1,2),f,'k')
    %     plot(f,P1(:,1))
    %     title('Mean Fourier transform')
        title({'Mean FT. ',sprintf('Peak period = %.1f', mean(1./f(maxIdx)))})
        ylabel('Frequency (1/cells)')
        xlabel('Power')
        pbaspect([0.5 1.4 1])
        set(gca,'FontSize',fontSize)
        ylim([0 0.5])

        subplot(2,6,[7 8 9 10 11 12 ])
    %     [~,II]=max(P1,[],1);
    %     plot(split_t,FisherPer,'.','color',cmap(size(cmap,1)/2,:),'linewidth',2);hold on
        plot(split_t,FisherPer,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:));hold on

        title(sprintf('Fisher G significant period detected %.0f%% of the time', 100*Occurrence))
    %     title('Significant period (Fisher G)')
        ylabel('Significant period (cells)')
        xlabel('Column')
        set(gca,'FontSize',fontSize)
        xlim([split_t(1) split_t(end)])
        ylim([0 5])

        centre=(split_t(end)-split_t(1))./((cols-1)*2);
        xticks(linspace(split_t(1)+centre,split_t(end)-centre,(cols-1)))
        xticklabels({num2str([2:(cols)]')})
        h=gca; h.XAxis.TickLength = [0 0];

        for n=1:numel(clipEnd)
            plot([split_t(clipEnd(n)), split_t(clipEnd(n))],[0,10],'k')
        end

        drawnow
    end
    
    
    %% Subplots of all columns/kymos
    clipPos=t_split(2:end)-t_split(1:end-1);
    clipPos(clipPos<0)=-1; %End positions marked by -1
    clipPos(clipPos>0)=0;  %Fill rest with zeros
    startPos=find(clipPos==-1)+1; %Find start positions of clips
    clipPos(startPos)=1;
    clipPos=[clipPos, -1];
    clipPos(1)=1;
    
    clipStart=find(clipPos==1);
    clipEnd=find(clipPos==-1);
    
    plotCols=3;
    plotRows=cols/plotCols;
    plotCols=plotCols*2;
    
%     figure(58) 
%     clf
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
%     
%     for n=1:numel(clipStart)
% 
%         tClip=t_split(clipStart(n):clipEnd(n));
%         P1Clip=P1(:,clipStart(n):clipEnd(n));
%         P1Clip=P1Clip./max(P1Clip,[],'all');
% 
%         [FisherPer2,Occurrence2]=Fisher(P1Clip,f);
%         
%         subplot(plotRows, plotCols, 2*n-1)
%         load cmap
%         colormap(gca,cmap)
% 
%         h=surf(tClip,f,P1Clip);
%     %     xlabel('Time (hours)')
%     %     ylabel('Frequency (1/cell)')
% %         title('Spatial Fourier transform over time')
%         title(sprintf('Occ=%.1f',Occurrence2))
%         xlim([t_split(1) tf_hours])
%         ylim([f(1) f(end)])
%         view(0,90)
%         colorbar
%         set(h,'LineStyle','none')
%         set(gca,'FontSize',8)
%     %     pbaspect([2 1 1])
%         colormap(gca,cmap)
% 
%         subplot(plotRows, plotCols, 2*n)
%         [~,maxIdx]=max(mean(P1Clip,2));
%         plot(mean(P1Clip,2),f,'w')
%         ylim([f(1) f(end)])
%     %     title('Mean Fourier transform')
%     %     title(sprintf('Mean FT. Peak period = %.1f', mean(1./f(maxIdx))))
%     %     ylabel('Frequency (1/cells)')
%     %     xlabel('Contribution/power')
%         pbaspect([0.5 1 1])
%         set(gca,'FontSize',8)
%     
%     end
%     
%     figure(57)
%     clf
%     for n=1:numel(clipStart)
% 
%         tClip=t_split(clipStart(n):clipEnd(n));
%         P1Clip=P1(:,clipStart(n):clipEnd(n));
%         P1Clip=P1Clip./max(P1Clip,[],'all');
% 
%         [FisherPer3,Occurrence3]=Fisher(P1Clip,f);
% 
%         subplot(plotRows, plotCols/2, n)
%     
%         plot(tClip,FisherPer3,'.','color',cmap(size(cmap,1)/2,:),'linewidth',2)
%         title(sprintf('Fisher G significance: %.1f%%', 100*Occurrence3))
%     %     title('Significant period (Fisher G)')
%         ylabel('Significant period (cells)')
%         xlabel('Time (hours)')
%         set(gca,'FontSize',8)
%         ylim([0 5])
%         xlim([tClip(1) tClip(end)])
%     end
    
    thresh=0.6;
    kymo=KymLong;

    for i=1:size(kymo,2)
        vect=kymo(:,i);
        [l(i),p(i)]=microClustDetect(vect,thresh);
    end
    clusterMean=nanmean(l)
    clusterOcc=1-sum(isnan(l))/numel(l)
%     figure(13),histogram(l,'normalization','probability'),ylabel('Frequency/occurence'),xlabel('Cluster size (radius)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(56)
% %     set(gcf,'renderer','Painters')
%     clf
%     N=cols-1;
%    
%     for n=1:N
% %     n=2; %column to use for kymo
%     
%     subplot(3,N,n)
%     X=t;
%     Y=1:rows;
%     j=surf(X,Y,squeeze(Yraw(:,:,n+1)));
%     title(sprintf('Kymograph %.f',n+1))
%     ylim([Y(1) Y(end)])
%     xlim([t(1) t(end)])
%     xlabel('Time (hrs)')
%     ylabel('Row index')
%     set(j, 'LineStyle','none')
%     set(gca,'YTickLabel',[]);
%     set(gca,'FontSize',8)
%     colormap(gca,vid)
% %     colorbar
%     view(0,90)
%     
%     tClip=t_split(clipStart(n):clipEnd(n));
%     P1Clip=P1(:,clipStart(n):clipEnd(n));
%     P1Clip=P1Clip./max(P1Clip,[],'all');
% 
%     [FisherPer4,Occurrence4]=Fisher(P1Clip,f);
% 
%     subplot(3,N,n+N)
%     h=surf(tClip,f,P1Clip);
%     %     xlabel('Time (hours)')
%     %     ylabel('Frequency (1/cell)')
% %         title('Spatial Fourier transform over time')
%     title(sprintf('Occ=%.1f',Occurrence4))
%     xlim([t_split(1) tf_hours])
%     ylim([f(1) f(end)])
%     view(0,90)
% %     colorbar
%     set(h,'LineStyle','none')
%     set(gca,'FontSize',8)
%     %     pbaspect([2 1 1])
%     colormap(gca,cmap)
%     
%     
%     subplot(3,N,n+2*N)
% 
%     plot(tClip,FisherPer4,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:))
%     title(sprintf('Significant: %.1f%%', 100*Occurrence4))
%     ylabel('Significant period (cells)')
%     xlabel('Time (hours)')
%     set(gca,'FontSize',8)
%     ylim([0 5])
%     xlim([tClip(1) tClip(end)])
%     
%     end
        

end


%==========================================================================
%%                  Phase space of mRNA and protein
%==========================================================================
if AnimatePhaseSpace==1
    idx=0;
    start=Nt/2;
    cell=1;
    for t=start:round(AnimationSpeed*Nt/400):Nt
        figure(101)
        plot(P(cell,start:t),MM(cell,start:t),'--','color','w');hold on
        plot(P(cell,t),MM(cell,t),'o','linewidth',6,'color',[0.8 0.2 0.2]); hold off
        xlabel('Protein count','color',[64 171 153]./255)
        ylabel('mRNA count','color',[205 102 120]./255)
        xlim([min(P(1,0.01*Nt:end)) max(P(1,0.01*Nt:end))])
        ylim([min(MM(1,0.01*Nt:end)) max(MM(1,0.01*Nt:end))])
        title(sprintf('Time=%.0f hours',T(t)/60))
        set(gca,'fontsize',20)
        axis square
        drawnow

%         if MakeGIF==1
%             open(vidObj);   
%         end

        if saveAnimation==1
            idx=idx+1;
            F1=getframe(gcf); %Make video with title and colour bar included
            im{idx}=frame2im(F1);
%             writeVideo(vidObj,F);
        end
    end
   
    if saveAnimation==1
        IDX=idx;
        for idx = 1:IDX
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,filename3,'gif','LoopCount',Inf,'DelayTime',1/30);
            else
                imwrite(A,map,filename3,'gif','WriteMode','append','DelayTime',1/30);
            end
        end
    end    
end


%% Swapping of neighbours analysis
if ImplementSwapping==1
    add=meshgrid(0:cells:cells*(Nt+1)-cells, 1:cells);
    [~,sortMatrix]=sort(CT);
    recoverCellPos=sortMatrix+add;

    P_recover=P(recoverCellPos); %Recovering cell positions so that individual cell dynamics can be plotted

    cellPos=cellPosition(CT,rows,cols);
    cellPosZero=cellPos-repmat(cellPos(:,1),1,Nt+1);
    
    %Expected movement over 2.5h
    tWin=2.5;
    tEl=tWin*60/dt;
%     size(cellPos(:,1+tEl:end))
%     size(cellPos(:,1:end-tEl))
    mean(abs(cellPos(:,1+tEl:end)-cellPos(:,1:end-tEl)),'all')

    Dist=cellPos(:,end)-cellPos(:,1);
    expectedRootMeanSquare=sqrt((1-SwapThresh)*tf/Ts);
    rootMeanSquare=sqrt(mean(Dist.^2));

    figure(22123)
    clf
    subplot(221)
    plot(T/60,cellPos(1:20:cells,:))
    ylabel('Distance (cells)')
    xlabel('Time (hours)')
    title('Trajectories of cells')
    subplot(222)
    plot(T/60,cellPosZero,'w','linewidth',0.05)
    ylabel('Distance (cells)')
    xlabel('Time (hours)')
    title({'Trajectories of cells','starting from zero'})


    subplot(223)
    histogram(Dist)
    xlim([min(Dist) max(Dist)])
    xlabel('Distance (cells)')
    ylabel('Frequency')
    title({'Distance moved', 'from inital position'})
    subplot(224)
    histogram(sqrt(Dist.^2));hold on
    scatter(rootMeanSquare,0,'r','linewidth',5)
    xlim([0 max(Dist)])
    xlabel('Distance (cells)')
    ylabel('Frequency')
    title(sprintf('Expected RMS distance= %.1f cells \nSimulation RMS distance = %.1f cells',expectedRootMeanSquare, rootMeanSquare))
    drawnow

    avg_dist=zeros(1,cols);
    n_arr=zeros(1,cols);
    
    for n=1:cols
        avg_dist(n)=mean(abs(cellPosZero(cellPos(:,1)==n,end)));
        n_arr(n)=n;
    end
    figure(9786)
    plot(n_arr,avg_dist)
end

