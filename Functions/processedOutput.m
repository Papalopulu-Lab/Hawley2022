function [output]=processedOutput(p)

    rows                  = p.rows;
    cols                  = p.cols;
    t0                    = p.t0;
    tf_hours              = p.tf_hours;
    Stochastic            = p.Stochastic;
    CoupledCells          = p.CoupledCells;
    Boundary              = p.Boundary;             
    TurnOffAutorepression = p.TurnOffAutorepression;
    
    cells                 = p.cells;
    tf                    = p.tf;
    dt                    = p.dt;
    Nt                    = p.Nt;
    T                     = p.T;
    VertProtrusions       = p.VertProtrusions;
    HorzProtrusions       = p.HorzProtrusions;
    StrengthGrad          = p.StrengthGrad;
    numDistalNeigh        = p.numDistalNeigh;
%     HorzProtStrength      = p.HorzProtStrength;
    CrudeDiff             = p.CrudeDiff;
    AvM                   = p.AvM;
    wl_mins               = p.wl_mins;
    wl                    = p.wl;
    diffSignallingMethod  = p.diffSignallingMethod;
    diffLength            = p.diffLength;
    DiffTime_hours        = p.DiffTime_hours;
    DiffTime              = p.DiffTime;
    S                     = p.S;
    nominalGamma          = p.nominalGamma;
    ImplementSwapping     = p.ImplementSwapping;
    Pm                    = p.Pm;
    SwapThresh            = p.SwapThresh;
    Ts                    = p.Ts;
    HH                    = p.HH;
    II                    = p.II;
    JJ                    = p.JJ;
    KK                    = p.KK;
    plotSwitchTime        = p.plotSwitchTime;
    SpatialFourier        = p.SpatialFourier;

%% Preallocation for parfor loop
FisherPer=zeros(HH,II,JJ);
Occurrence=zeros(HH,II,JJ);
Dc=zeros(HH,II,JJ);
meanHi=zeros(HH,II,JJ);
meanLo=zeros(HH,II,JJ);
diffRate=zeros(HH,II,JJ);
CoV=zeros(HH,II,JJ);

parfor_progress(HH*II*JJ*KK); 
%% Loops 
for hh=1:HH
    for ii=1:II    
        parfor jj=1:JJ
            
            VertProtStrength=p.jj_arr(jj);
            HorzProtStrength=VertProtStrength;
            
            if p.paramSpaceOutput == 1
                diffGamma=p.hh_arr(hh); 
                TauND = p.TauND;
            elseif p.paramSpaceOutput == 2
                TauND=p.hh_arr(hh); 
                diffGamma = p.diffGamma;
            end
            
            recievingGamma=diffGamma;
            
            %Other parfor loop variables required to be moved inside loop
            S=p.S;
            DiffTime_hours=50;         %Time at which differention can start to occur (hours)
            DiffTime=DiffTime_hours*60/dt;
            

            [NM,NumNeigh]=neighbours(rows,cols,Boundary,VertProtrusions,HorzProtrusions,VertProtStrength,HorzProtStrength,StrengthGrad,numDistalNeigh, 0); 
            %               neighbours(rows,cols,boundary,VertProtrusions,HorzProtrusions,VertProtStrength,HorzProtStrength)
            NM=sparse(NM); % Saves a lot on computational cost for large grid sizes!
            NM=NM./max(NM(:));

            % NM(NM==0.5)=0.2;

            % figure(21)
            % clf
            % imagesc(NM)
            % % load cmap
            % colormap(randColourMapBright(100,[0 0 0]))
            % colorbar
            % axis square
            % % colormap(hsv)

            eps=(1./NumNeigh); % Make this as an output of neighbours.m!
            eps=eps*(1/6)/min(eps);
            if cells==1
                eps=0;
            end


            %__________________________________________________________________________
            %Single-cell parameters

            P_ND0 = p.ii_arr(ii); 
 
            a_m                   = p.a_m;
            a_p                   = p.a_p;
            P_H0                  = p.P_H0;
            TauH                  = p.TauH;
            n_H                   = p.n_H;
            u_m                   = p.u_m;
            u_p                   = p.u_p;
%             P_ND0                 = p.P_ND0;
            n_ND                  = p.n_ND;
%             TauND                 = p.TauND;
            
    
            Up=1;
            Low=0;
            a=Low;
            b=Up-a;
            %__________________________________________________________________________

            gamma=1; % Maximum intercellular Hill function value
%             gammaRand=rnd(cells,1,0,1);

            TauH_step=round(TauH/dt);   % Conversion to simulation time steps
            TauND_step=round(TauND/dt); % Conversion to simulation time steps

            if CoupledCells==0
                eps=eps.*0; %Set eps=0 to decouple the cells 
            end
            
            if eps==0
                fprintf('\nCells are uncoupled!\n')
            end
            if TurnOffAutorepression==1
                fprintf('\nAutorepression is turned off!\n')
            end

            if CrudeDiff==0
                DiffTime=Nt*2;
            end

            %% DDEs (using anonymous functions (@ symbol) rather than standard functions)
            dm=@(m,p,p_delay,Psum_delay,gamma,eps) a_m*1./(1+(p_delay./P_H0).^n_H).*gamma.*(  a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND)  ) - u_m.*m; % Describe mRNA change in time


            if TurnOffAutorepression==1
                dm=@(m,p,p_delay,Psum_delay,gamma,eps) a_m.*gamma.* (a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND)) - u_m.*m; % No self repression version for simple Notch Delta modelling
            end

            dp=@(m,p,p_delay,Psum_delay,gamma,eps) a_p.*m - u_p.*p; % Describes protein change in time

            %% SDDEs
            dm1=dm;

            dm2=@(m,p,p_delay,Psum_delay,gamma,eps) sqrt(a_m.*gamma./(1+(p_delay./P_H0).^n_H).*( a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND) ) + u_m.*m);

            if TurnOffAutorepression==1
               dm2=@(m,p,p_delay,Psum_delay,gamma,eps) sqrt(a_m.*gamma.*( a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND) ) + u_m.*m);
            end

            dp1=dp;
            dp2=@(m,p,p_delay,Psum_delay,gamma,eps) sqrt(a_p.*m + u_p.*p);

            
%==========================================================================
%%                       Main loop (solves DDEs)
%==========================================================================
            Dc_kk=zeros(1,KK);
            meanHi_kk=zeros(1,KK);
            meanLo_kk=zeros(1,KK);  
            FisherPer_kk=zeros(1,KK);
            Occurrence_kk=zeros(1,KK);
            diffRate_kk=zeros(1,KK);
            CoV_kk=zeros(1,KK);

            for kk=1:KK
                if KK==1
                    rng(3); % Same seed for random number generator (prevents initial heterogeneities from changing on every simulation, so comment out this line if you want different initial conditions each time)
                end
                
                P=[rndrng(cells,1,0.2*P_H0,1.8*P_H0),zeros(cells,Nt-1),zeros(cells,1)]; %Vector that will store protein values from the simulation
                M=[rndrng(cells,1,0,20),zeros(cells,Nt-1),zeros(cells,1)];              %Vector that will store mRNA values from the simulation

                P=[rndrng(cells,1,0.2*P_H0,1.8*P_H0),zeros(cells,Nt-1),zeros(cells,1)]; %Vector that will store protein values from the simulation



                PercentCounter=0; %Print progress of the differential solver to the command line (1=yes, 0=no)
                [P, M, t_step, DiffYN, DiffYNflash, CT, DiffThresh, M, AbsThresh, GAMMA,countArray,EPS]=diffSolver(P, M, Nt, TauND_step, TauH_step, NM, gamma, dm, dp, dm1, dm2, dp1, dp2, dt, Stochastic, rows, cols, DiffTime, S, ImplementSwapping, SwapThresh, Ts, PercentCounter, AvM, diffSignallingMethod, wl, eps,diffLength,diffGamma,recievingGamma,nominalGamma);


                %% Differentiation rate 
                if CrudeDiff==1
                    DiffSignal=DiffYNflash(:,DiffTime/dt:end);
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

                    [hi, lo, DcCells]=HiLoLoop(pop); %This replaces code below

                    hi=hi*dt/60;
                    lo=lo*dt/60;
                    hi_all=hi;
                    lo_all=lo;

                    hiX=hi*0+1;
                    loX=lo*0+2;
                end

                %% Spatial Fourier
                if SpatialFourier==1
                    startTime=tf_hours-100;
                    time_frac=startTime/tf_hours;

                    if time_frac>=1
                        error('Increase length of simulation or decrease startTime for spatial frequency analysis.')
                    end

                    t=T(time_frac*Nt:end)./60;
                    ti=time_frac*Nt;
                    
                    if cols==1
                        kymoCols=1;
                    else
                        kymoCols=2; % This value should be 1 or 2. Value of 1 means spatial signal is taken from each individual column, value of 2 means the hexagonal geometry is taken into account and a kymo slice of width 1 is taken but averages the odd rows where two cells are present
                    end
                    
                    [Pnew]=improvedKymo(P,rows,cols,kymoCols); %Improved selection of kymograph (Set to 1 or 2 columns)

%                     figure(98)
%                     clf
%                     subplot(211)
%                     imagesc(Pnew)

                    kymoWidth=1; %Number of cells to use for kymograph selection width
                    K=cols/kymoWidth;

                    if cols==1
                        Y_RAW=P(:,ti:end);
                        t_concat=t;
                        K=2;
                    else

                        if mod(cols,2)==1
                            error('Make cols even')
                        end

                        [Y_RAW,t_concat]=stackCols(Pnew,K,rows,t,ti);
              
                    end

%                     figure(101)
%                     clf
%                     subplot(211)
%                     imagesc(Y_RAW)

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
    
                    Y=fft(YDETREND2,[],1); % Fast Fourier Transform
                    
                    winL=length(Y(:,1));
                    P2 = abs(Y/winL).^2;
                    P1 = P2(1:winL/2+1,:);
                    P1(2:end-1,:) = 2*P1(2:end-1,:);

                    Fs=1;
                    f = Fs*(0:(winL/2))/winL;

                    KymLong=YDETREND2;


                    [FP,Occ]=Fisher(P1,f);

                end
                
                FisherPer_kk(kk)=nanmean(FP);
                Occurrence_kk(kk)=Occ;
                Dc_kk(kk)=mean(DcCells);
                meanHi_kk(kk)=mean(hi);
                meanLo_kk(kk)=mean(lo);  

                diffRate_kk(kk)=100*sum(DiffSignal(:))/((tf_hours-DiffTime_hours)*cells);
                CoV_kk(kk)=std(P(:,Nt/2:end),[],'all')/mean(P(:,Nt/2:end),'all');

                clc
                parfor_progress;
    
    
            end
            Dc(hh,ii,jj)=nanmean(Dc_kk);
            meanHi(hh,ii,jj)=nanmean(meanHi_kk);
            meanLo(hh,ii,jj)=nanmean(meanLo_kk);
            FisherPer(hh,ii,jj)=nanmean(FisherPer_kk);
            Occurrence(hh,ii,jj)=nanmean(Occurrence_kk);
            diffRate(hh,ii,jj)=nanmean(diffRate_kk);
            CoV(hh,ii,jj)=nanmean(CoV_kk);

        end
    end
end

output.Dc=Dc;
output.meanHi=meanHi;
output.meanLo=meanLo;
output.FisherPer=FisherPer;
output.Occurrence=Occurrence;
output.diffRate=diffRate;
output.CoV=CoV;
