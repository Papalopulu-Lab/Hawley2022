% ParameterSpace_Fig7and8.m Multicellular model of autonomous HES5 noisy/oscillatory
% dynamic coupled between cells by a lateral inhibition Hill function
% interaction. This code generates the parameter space ouputs in Supplementary Fig2&3
% in the paper 'Dynamic Switching of Lateral Inhibition Spatial Patterns' 
% (Joshua Hawley, Paul Glendinning, and Nancy Papalopulu) (2022). 
% Correspondence: joshhawley369@gmail.com

% Required toolboxes: Signal analysis toolbox, statistics and machine
% learning toolbox, Parallel computing toolbox

% HOW TO USE THIS CODE
% To skip generating data (which takes a long time), you can load the
% workspaces in the location: 
% Workspaces used to generate figures in paper/Param space output for time delay exploration - Supp Figure 2 and 3
% Once you have opened one of the workspaces, run this first section only so
% that the Fucntions folder is added, then scroll down to the Plots section
% (line 185) and run that section to generate the graphs.

clear;clc;
addpath('Functions')

%% Define grid size
rows=26;         %Number of rows of cells in hexagonal grid 30
cols=1 ;         %Number of columns of cells in hexagonal grid 12

%% Define simulation time and step size
t0=0;            % Start time
tf_hours=200;    % Final time (hours)

%% Set up simulation inputs

%Core options
Stochastic=1;              %1 = stochastic, 0 = deterministic  (deterministic currently has a bug) 
CoupledCells=1;            %Simulate with Notch-Delta coupling = 1
Boundary=0;                %0 = 'hard' boundary, 1 = periodic boundary (set to 1)
TurnOffAutorepression=0;   %Reduce model to just lateral inhibition without autonomous Hes oscillators in each cell

%Extended neighbour options 
VertProtrusions=2; HorzProtrusions=VertProtrusions; %Increases number of signalling neighbours in the vertical direction 
StrengthGrad=0;
numDistalNeigh=4;          % 4, 6, or 10

%Differentiation options
CrudeDiff            = 1;          %Cells will be marked...
AvM                  = 0;          %Absolute vs moving mean threshold (0=Absolute thresh, 1=Moving mean thresh)
diffSignallingMethod = 1;         %Determines what differentiation does to the system: 0=no change to system, 1=signalling changes in the cell above or below (randomly selected), 2= same as 1 but the cells selected are always the differentiating cell and the one below it, 3= replace diff cell with population mean, 4=local/neighbour mean replacement, 5=replace with random neighbour level,
diffLengthHours      = 0;          %Time spent in altered coupling strength 
S                    = 0.003;      %Probability of differentiation
%% Set up parameter space size/repeats

HH=4;   %Number of points to plot in parameter space - diffLength (Z-axis)
II=20;  %Number of points to plot in parameter space - Repression threshold (y-axis)
JJ=20;  %Number of points to plot in parameter space - Distal signalling efficiency (x-axis)
KK=20;  %Number of times to repeat a simulation per parameter point with random initial conditions

if HH==1
    hh_arr=0;                   % diffLength (h)
else   
    hh_arr = [0 30 60 100];     % diffLength (h)
end

ii_arr=linspace(1000,10000,II); %Repression threshold for NICD
jj_arr=linspace(0,3,JJ);        %Distal signalling efficiency 

%% Model parameters
% Columns in accepted_parameters: 1:a_m, 2:a_p, 3:P_H0, 4:TauH, 5:n_H
load accepted_parameters %Selcted for coherence values between 0.05 and 0.015
load summary_stats

%Select a parameter set from the Bayesian inferred single cell parameters 
%from (Manning et al. 2019) 
AP=accepted_parameters;
SS=summary_stats;

s=3700; %0.15 Coh %%%%Nominal param set%%%%

%% Other parameter set-up
cells=rows*cols; %Total number of cells
tf=tf_hours*60;  % Final time (min)

dt=2;            % Time step size (min)
if Stochastic==1
    dt=1;        % Stochastic uses Euler method rather than Runge-Kutta, so this requires a smaller step size      
end
Nt=(tf-t0)/dt;   % Number of time elements
T=t0:dt:tf;      % Time vector

if cells==1
    Boundary=0; % Other boundaries don't make sense for 1 cell!
end

%% Cell-movement/swapping parameters
ImplementSwapping   = 0;   %Cells will randomly swap in the horizontal direction
Pm=0.005; % Total probability of movement either left or right (max Pm value is 0.5) Nominal value 0.005
SwapThresh=1-Pm/2; % Treshold to be used in probability of swapping function within diffSolver.m
Ts=5; % Time (in mins) between each swapping event (nominal Ts=5)   

%% Differentiation selection and parameters
wl_mins             = 100;        % Window length in mins to take the moving mean from
wl=wl_mins/dt;                    % Window length converted to number of vector elements
diffLength=diffLengthHours*60/dt;
DiffTime_hours=50;                % Time at which differention can start to occur (hours)
DiffTime=DiffTime_hours*60/dt;
nominalGamma=1;                   % Value of gamma when cells are not differentiating
lowerCoupling=1;

%% Cell-movement/swapping parameters
ImplementSwapping   = 0;  %Cells will randomly swap in the horizontal direction
Pm=0.005;       % Total probability of movement either left or right (max Pm value is 0.5)
SwapThresh=1-Pm/2; % Treshold to be used in probability of swapping function within diffSolver.m
Ts=5;              % Time (in mins) between each swapping event (nominal Ts=5)       

%% Use object oriented variable to pass parameters to the processedOutput function 
p.hh_arr                = hh_arr; %differentiation gamma
p.ii_arr                = ii_arr; %P_ND0
p.jj_arr                = jj_arr; %VertProtrusionStrength

p.a_m                   = AP(s,1);
p.a_p                   = AP(s,2);
p.P_H0                  = AP(s,3);
p.TauH                  = AP(s,4);
p.n_H                   = AP(s,5);
p.u_m                   = log(2)/30;
p.u_p                   = log(2)/90;
p.n_ND                  = 3;
% p.TauND                 = hh_arr;

p.VertProtrusions       = VertProtrusions;
p.HorzProtrusions       = HorzProtrusions;
p.StrengthGrad          = StrengthGrad;
p.numDistalNeigh        = numDistalNeigh;

p.rows                  = rows;
p.cols                  = cols;
p.t0                    = t0;
p.tf_hours              = tf_hours;
p.Stochastic            = Stochastic;
p.CoupledCells          = CoupledCells;
p.Boundary              = Boundary;
p.TurnOffAutorepression = TurnOffAutorepression;
p.cells                 = cells;
p.tf                    = tf;
p.dt                    = dt;
p.Nt                    = Nt;
p.T                     = T;
p.CrudeDiff             = CrudeDiff;
p.AvM                   = AvM;
p.wl_mins               = wl_mins;
p.wl                    = wl;
p.diffSignallingMethod  = diffSignallingMethod;
p.diffLength            = diffLength;
p.DiffTime_hours        = DiffTime_hours;
p.DiffTime              = DiffTime;
p.S                     = S;
p.ImplementSwapping     = ImplementSwapping;
p.Pm                    = Pm;
p.SwapThresh            = SwapThresh;
p.Ts                    = Ts;
p.HH                    = HH;
p.II                    = II;
p.JJ                    = JJ;
p.KK                    = KK;
p.S                     = S;
p.nominalGamma          = nominalGamma;
p.diffGamma             = 3;
p.paramSpaceOutput      = 2;

%% Run Model
[output]=processedOutput(p);

%Save workspace
if cols==1
    str1="Single col, ";
else
    str1="Grid, ";
end

str2 = sprintf("Stochastic = %.0f, diffSignallingMethod = %.f, TauND exploration 0 30 60 100, Repeats = %.0f, S = %.4f, diffLength = %.1fh, HH = %.0f, II = %.0f, JJ = %.0f, numDistalNeigh = %.0f, vertProtrusions = %.1f, turnOffAutorepression = %.0f.mat", Stochastic, diffSignallingMethod, KK, S, diffLengthHours, HH, II, JJ, numDistalNeigh, VertProtrusions, TurnOffAutorepression);
STR=append("ParamSpaceOutput - ",str1,str2);
save(STR)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose the colour of all graphs
GraphAppearance = 0;   % 0 = Normal, 1 = Dark grey, 2=dark green/blue

if GraphAppearance==2   
    BackgroundColour=1/255*[0 28 22]; TextColour=[1 1 1];
    get(0,'Factory');                           set(0,'defaultfigurecolor',BackgroundColour)
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor',BackgroundColour)
    set(0,'defaultAxesXColor',TextColour);      set(0,'defaultAxesYColor',TextColour)
    set(0,'defaultLegendTextColor',TextColour); set(0,'defaultTextColor',TextColour)  
    set(0,'defaultAxesZColor',TextColour);
elseif GraphAppearance==1   
    BackgroundColour=38/255*[1 1 1]; TextColour=[1 1 1];
    get(0,'Factory');                           set(0,'defaultfigurecolor',BackgroundColour)
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor',BackgroundColour)
    set(0,'defaultAxesXColor',TextColour);      set(0,'defaultAxesYColor',TextColour)
    set(0,'defaultLegendTextColor',TextColour); set(0,'defaultTextColor',TextColour)  
    set(0,'defaultAxesZColor',TextColour);
elseif GraphAppearance==0
    get(0,'Factory');                           set(0,'defaultfigurecolor',[0.94 0.94 0.94])
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor','w')
    set(0,'defaultAxesXColor','k');             set(0,'defaultAxesYColor','k')
    set(0,'defaultLegendTextColor','k');        set(0,'defaultTextColor','k')
    BackgroundColour=[1 1 1];
end

load 'cmap.mat'


figure(1)
clf
set(gcf,'renderer','Painters') %For EPS file export
colbar=1;
textSize=10;

colourMap = magma(300);

for hh=1:HH
    
    figure(1)

    data=squeeze(output.FisherPer(hh,:,:));

    subplot(3,HH,hh)
    surf(jj_arr,ii_arr,data);
    title({sprintf('Spatial period | Diff length = %.0fh',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Distal signalling efficiency');   
    caxis([2 4])
    if hh==1
        ylabel({'Repression','threshold'});
    end
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,colourMap);
    if colbar==1
        colorbar
    end
    ax = gca;
    ax.YAxis.Exponent = 3;

%_____________________________________________________________________ 


    data=squeeze(output.Occurrence(hh,:,:));
    
    subplot(3,HH,hh+HH)
    surf(jj_arr,ii_arr,data);
    title({sprintf('Occurrence | Diff length = %.0fh',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Distal signalling efficiency');    
    if hh==1
        ylabel({'Repression','threshold'});
    end
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,colourMap);
    if colbar==1
        colorbar
    end
    ax = gca;
    ax.YAxis.Exponent = 3;

%__________________________________________________________________________

    data=squeeze(output.Dc(hh,:,:));

    subplot(3,HH,hh+2*HH)
    surf(jj_arr,ii_arr,data);
    title({sprintf('Dc | Diff length = %.0fh',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Distal signalling efficiency');  
    caxis([0 1])
    if hh==1
        ylabel({'Repression','threshold'});
    end
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,colourMap);
    if colbar==1
        colorbar
    end
    ax = gca;
    ax.YAxis.Exponent = 3;

%__________________________________________________________________________

    figure(2)
    set(gcf,'renderer','Painters') %For EPS file export
    PerLow=3;
    OccLow=0.4;
    DcLow=0.4;
    
    passPer=output.FisherPer(hh,:,:);
    passOcc=output.Occurrence(hh,:,:);
    passDc=output.Dc(hh,:,:);
    
    passPer(passPer<PerLow)=0;
    passPer(passPer>PerLow)=1;
    passOcc(passOcc<OccLow)=0;
    passOcc(passOcc>OccLow)=1;
    passDc(passDc<DcLow)=0;
    passDc(passDc>DcLow)=1;
    
    pass=passPer+passOcc+passDc;
    pass(pass<3)=0;
    pass(pass==3)=1;
    
    
    subplot(1,HH,hh)
    surf(jj_arr,ii_arr,squeeze(pass));
    sgtitle(sprintf('Spatial Per > %.1f  |  Occurrence > %.1f  |  Dc > %.1f', PerLow, OccLow, DcLow))
    title({sprintf('Diff length = %.0fh', hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Distal signalling efficiency');  
    colormap([BackgroundColour; 1.6*[0/256 95/256 0/256]])
    colormap([BackgroundColour; [220/256 60/256 60/256]])
    if hh==1
        ylabel({'Repression','threshold'});
    end
    axis square; view(0,90); set(gca,'fontsize',textSize);
    ax = gca;
    ax.YAxis.Exponent = 3; 
end

%__________________________________________________________________________

figure(3)
clf
set(gcf,'renderer','Painters') %For EPS file export
colbar=1;
textSize=12;

for hh=1:HH
    
    figure(3)

    data=squeeze(output.meanLo(hh,:,:));

    subplot(2,HH,hh)
    surf(jj_arr,ii_arr,data);
    title({sprintf('Mean low persistence (h) | Diff length= %.0fh',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Distal signalling efficiency');   
    if hh==1
        ylabel({'Repression','threshold'});
    end
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
    if colbar==1
        colorbar
    end
    
    data=squeeze(output.meanHi(hh,:,:));
    subplot(2,HH,hh+HH)
    surf(jj_arr,ii_arr,data);
    title({sprintf('Mean high persistence (h) | Diff length= %.0fh',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Distal signalling efficiency');    
    if hh==1
        ylabel({'Repression','threshold'});
    end
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
    if colbar==1
        colorbar
    end
end