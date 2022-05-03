%diffSolver.m Solves the differential equations specified in Model.m either
%deterministically or stochastically

function [P, M, t_step, DiffYN, DiffYNflash, CT,DiffThresh, MM, AbsThresh, GAMMA, countArray,EPS]=diffSolver(P, M, Nt, TauND_step, TauH_step, NM, gamma, dm, dp, dm1, dm2, dp1, dp2, dt, Stochastic, rows, cols, DiffTime, S, ImplementSwapping, SwapThresh, Ts, PercentCounter, AvM, diffSignallingMethod, wl,eps,diffLength,diffGamma,recievingGamma,nominalGamma)
cells=rows*cols;
CT=[(1:cells)' ,zeros(cells,Nt)];
MM=P.*0;
count=0;
initial=1;
GAMMA=P.*0+nominalGamma;
% epsOriginal=eps;
EPS=repmat(eps,1,size(P,2));

countArray=P.*0;
for t_step=1:Nt
%     eps=EPS(:,t_step); 
    %% Percentage counter for the terminal
    if PercentCounter==1
        if ~mod(t_step,floor(Nt*5/100)) 
            clc
            fprintf('Main loop progress: %.0f%% \n', t_step/Nt*100);
        end
    end
    
    %% Cell swapping
    if ImplementSwapping==1
        SwapTheCells=0;
        if mod(t_step-1,Ts)==0
        SwapTheCells=1;
        [P_swap,P_binary]=probSwap(rows,cols,SwapThresh);

        [I,J]=find(P_binary==1); % I and J contains the ith cell to swap with the jth cell
        [~,posI,posJ]=intersect(I,J);

        if ~isempty(posI) %This is a 'if not empty' array condition to deal with a cell that is swapping twice in one move
            remove=0;
            for nv=1:length(posI)
                m1=I(posI(nv));
                n1=J(posI(nv));
                m2=I(posJ(nv));
                n2=J(posJ(nv));
                if P_swap(m1,n1)<P_swap(m2,n2)
                    remove(nv)=posI(nv);
%                     I(posI(nv))=-1;
%                     J(posI(nv))=-1;
                else
                    remove(nv)=posJ(nv);
%                     I(posJ(nv))=-1;
%                     J(posJ(nv))=-1;
                end
            end
            I(remove)=[];
            J(remove)=[];
        end

        ExMat=sparse(eye(cells)); %Exchange matrix

        ExMat=swapRow(ExMat,I,J);
        
%         figure(772)
%         surf(ExMat)
%         view(0,90)
%         drawnow
%         pause(0.1)
        
        
        P(:,t_step)=ExMat*P(:,t_step);
        M(:,t_step)=ExMat*M(:,t_step);
        end
        
        if SwapTheCells==1
%             ExMat
%             cellTracker(:,t_step)
%             ExMat*cellTracker(:,t_step)

%             cellTracker(:,t_step+1)=ExMat*cellTracker(:,t_step);
            CT(:,t_step)=ExMat*CT(:,t_step);
            CT(:,t_step+1)=CT(:,t_step);
            
        else
            CT(:,t_step+1)=CT(:,t_step);
        end
        
%         cellTracker
        
%         if mod(t_step-1,100)==0
%         figure(32347)
%         h=surf(cellTracker);
%         view(0,90)
%         set(h,'LineStyle','none')
%         colormap(jet)
%         drawnow
%         end
    end
    

    %% Crude differentiation method
    
    if t_step<DiffTime/dt
        DiffThresh=zeros(cells,Nt+1);
    end
    
    if t_step>=DiffTime/dt
        if t_step==DiffTime/dt
            DiffYN=zeros(cells,Nt);
            DiffYNflash=zeros(cells,Nt);
            wl_mins=1000;    %Window length in mins to take the moving mean from 
            wl=wl_mins/dt; %Window length converted to number of vector elements 
            
            AbsThresh=mean(mean(P(:,0.5*DiffTime/dt:DiffTime/dt))); %Absolute threshold
            AbsThreshM=mean(mean(M(:,0.5*DiffTime/dt:DiffTime/dt)));
        end
%         [MMThresh]=MovingMean(P,wl,t_step); %Moving mean threshold
        if AvM==0
%         DiffThresh(:,t_step+1)=(1-AvM).*AbsThresh + AvM.*MMThresh;
        DiffThresh(:,t_step+1)=AbsThresh;
        elseif AvM==1
            DiffThresh(:,t_step+1)=1*mean(P(:,t_step-wl:t_step),'all');
        end
            
        
        
        
%         ProbDiff=(DiffThresh(:,t_step)-P(:,t_step)).*S;
        ProbDiff=((DiffThresh(:,t_step)-P(:,t_step))./DiffThresh(:,t_step)).*S;
        ProbVal=ProbDiff-rand(cells,1);

        DiffYN(:,t_step+1)=ProbVal;
        DiffYN(ProbVal>0,t_step+1)=1;
        DiffYN(ProbVal<0,t_step+1)=0;
        DiffYNflash(:,t_step+1)=DiffYN(:,t_step+1);
        DiffYN(DiffYN(:,t_step)==1,t_step+1)=1; %If previously a cell has initiated a diff event, then keep it with a value of 1 for all subsequent time steps.
    else
        DiffYN=0;
        DiffYNflash=0; 
        AbsThresh=1;
    end
    
    %% Moving mean differentiation threshold
%     wl_mins=100;    %Window length in mins to take the moving mean from
%     wl=wl_mins/dt; %Window length converted to number of vector elements

    [MM(:,t_step)]=MovingMean(P,wl,t_step);
    
    %% Conditions for t<Tau (time delays)
    if max(TauND_step)+1>t_step
        Psum_delay=NM*P(:,1);
    else  
        Psum_delay=NM*P(:,t_step-TauND_step); %Average effect of Hes expressed by neighbouring 6 cells
    end

%     if max(TauH_step)+1>t_step
%         p_delay=P(:,1);
%     else  
%         p_delay=P(:,t_step-TauH_step);
%     end
    
    if max(TauH_step)+1>t_step
        p_delay=P(:,1);
    else
        if numel(TauH_step)==1
            p_delay=P(:,t_step-TauH_step);
        else
            p_delay=P(sub2ind(size(P),[1:size(P,1)]',t_step-TauH_step)); %Written this way to deal with the case when time delay is different across cells
        end
    end

%     if t_step<Nt/2
%         gamma=1;
%     else
%         gamma=1-(t_step-Nt/2)/(Nt/2);
%     end
%     GAMMA(t_step)=gamma;
%     gamma=1-t_step/Nt;

% % % % %     switchT_h=8;
% % % % %     switchT=switchT_h*60/dt;
% % % % %     count=count+1;
% % % % %     spatPer=2;
% % % % %     
% % % % %     
% % % % %     if count==switchT
% % % % %         count=0;
% % % % %         initial=initial+1;
% % % % %         if initial==spatPer+1
% % % % %             initial=1;
% % % % %         end
% % % % %     end
% % % % %     
% % % % %     
% % % % %     select=zeros(rows,cols);
% % % % %     select(initial:spatPer:end,:)=1;
% % % % %     select=select(:);
% % % % %     gamma=zeros(cells,1);
% % % % %     gamma(select==1)=0.3;
% % % % %     gamma(select==0)=3;
% gamma=gammaRand;
%     P(select==1,t_step+1)=1;
%     P(select==0,t_step+1)=0;
    if Stochastic==1
        [new1, new2]=EulerStoch(cells,M(:,t_step),P(:,t_step),p_delay,Psum_delay,gamma,EPS(:,t_step),dm1,dm2,dp1,dp2,dt);
    else
        [new1, new2]=RK(M(:,t_step),P(:,t_step),p_delay,Psum_delay,gamma,EPS(:,t_step),dm,dp,dt);                  % 4th order Runge-Kutta solver
    end
    

    M(:,t_step+1)=new1; % Storing new mRNA value
    P(:,t_step+1)=new2; % Storing new protein value
    
%     initial=1;
%     select=zeros(rows,cols);
%     select(initial:4:end,:)=1;
%     select=select(:);
%     P(select==1,t_step+1)=5e3;

    
%     startT=round(Nt/5);
%     startT=20*60/dt;
%     cellsWide=1;
%     if t_step>startT && t_step<startT+5*60
%         P(rows*cols/2+1:(rows)*(cols/2+cellsWide),t_step+1)=0;
%     end
    
%% Replace differentiating cells
    if diffSignallingMethod==1 && t_step>=DiffTime || diffSignallingMethod==2 && t_step>=DiffTime

        t1=t_step+1; %Time to start changed signalling in the gamma parameter
        t2=t_step+1+diffLength; %Time that the altered signalling should end
        
        %Find Cells that need updating on this time step
        diffCells=find(DiffYNflash(:,t_step+1)==1); %These cells will need to be updated with changed signalling in GAMMA

        %CountArray ensures a differentiation event doesnt happen while a
        %cell is already differentiating.
        stopDiffRepeat=find(countArray(diffCells,t_step)==1); %Check one time step behind to see if cell is already differentiating
        diffCells(stopDiffRepeat)=[];
        countArray(diffCells,t_step+1:t_step+1+diffLength)=1; %Diff cells and how long they diff for entries = 1 means they are differentiating
        
        plusMinusOne=[-1 1];        
        idx=randi(2,1);
        plusMinusOne=plusMinusOne(idx);
        
        if diffSignallingMethod==2
            plusMinusOne=1; %Replace cells are always on the same side of a differentiation
        end

        recievingCells=[diffCells+plusMinusOne, diffCells+rows, diffCells+rows+plusMinusOne];
%         recievingCells=diffCells+plusMinusOne; %Only change signal in 2 cells, not 4
            
            
        recievingCells(recievingCells>cells)=[];
        recievingCells(recievingCells<1)=[];
  
        GAMMA(recievingCells, t_step+1:t_step+1+diffLength)=recievingGamma;

        GAMMA(countArray(:,t_step+1)==1,t_step+1)=diffGamma;
        gamma=GAMMA(:,t_step+1);


    elseif diffSignallingMethod==3 && t_step>=DiffTime %Replace differentiating cells with population mean
        
        if sum(DiffYNflash(:,t_step+1))>0
            P(DiffYNflash(:,t_step+1)==1, t_step+1)=AbsThresh;
            M(DiffYNflash(:,t_step+1)==1, t_step+1)=AbsThreshM;
        end
        
    elseif diffSignallingMethod==4 && t_step>=DiffTime %Replace with neighbour/local average
        
        P_localAvg=eps.*(NM*P(:,t_step));
        M_localAvg=eps.*(NM*M(:,t_step));
%             DiffYNflash(:,t_step+1)==1
%             P_localAvg(DiffYNflash(:,t_step+1)==1)
%             P(DiffYNflash(:,t_step+1)==1, t_step+1)
        P(DiffYNflash(:,t_step+1)==1, t_step+1)=P_localAvg(DiffYNflash(:,t_step+1)==1);
        M(DiffYNflash(:,t_step+1)==1, t_step+1)=M_localAvg(DiffYNflash(:,t_step+1)==1);
        
    elseif diffSignallingMethod==5 && t_step>=DiffTime %Replace with a random cell from neighbourhood

        randNeigh=NM.*rand(size(NM)); %Gives random numbers in place of the ones in NM

        [~,randNeigh]=max(randNeigh,[],2);

        randNeigh(randNeigh==0)=1;
        randNeigh(randNeigh<0)=0;
        
        P_randReplace=P(randNeigh,t_step);
        M_randReplace=M(randNeigh,t_step);
        
        P(DiffYNflash(:,t_step+1)==1, t_step+1)=P_randReplace(DiffYNflash(:,t_step+1)==1);
        M(DiffYNflash(:,t_step+1)==1, t_step+1)=M_randReplace(DiffYNflash(:,t_step+1)==1);       


 
    end


end
end


function [P_swap,P_binary]=probSwap(r,c,thresh)
    
    P_vec=rand(r*(c-1),1);
    P_swap=sparse(diag(P_vec,r));
    
    P_binary=P_swap;
    P_binary(P_swap>thresh)=1;
    P_binary(P_swap<thresh)=0;
    P_binary=sparse(P_binary);
    
%     figure(8667)
%     clf
%     h=surf(P_binary)
%     set(h, 'edgecolor','none')
%     view(0,90)
%     caxis([0 1])
%     colorbar
%     drawnow
end

function [new1, new2]=EulerStoch(cells,x1,x2,x3,x4,x5,x6,f11,f12,f21,f22,dt)

    new1=x1 + f11(x1,x2,x3,x4,x5,x6)*dt + f12(x1,x2,x3,x4,x5,x6)*sqrt(dt).*normrnd(0,1,[cells, 1]);
    new2=x2 + f21(x1,x2,x3,x4,x5,x6)*dt + f22(x1,x2,x3,x4,x5,x6)*sqrt(dt).*normrnd(0,1,[cells, 1]);

    if min(new1)<0
        new1(new1<0)=0;
    end

    if min(new2)<0
        new2(new2<0)=0;
    end
end

%RK.m Implementation of the 4th-order Runge Kutta method. Outputs 
%the next step values for t+dt.

function [new1,new2]=RK(x1,x2,x3,x4,x5,x6,f1,f2,dt)

    % x1 = Protein conc
    k1_1=f1(x1,       x2,x3,x4,x5,x6)*dt;
    k2_1=f1(x1+k1_1/2,x2,x3,x4,x5,x6)*dt;
    k3_1=f1(x1+k2_1/2,x2,x3,x4,x5,x6)*dt;
    k4_1=f1(x1+k3_1,  x2,x3,x4,x5,x6)*dt;

    % x2 = mRNA conc
    k1_2=f2(x1,x2,       x3,x4,x5,x6)*dt;
    k2_2=f2(x1,x2+k1_2/2,x3,x4,x5,x6)*dt;
    k3_2=f2(x1,x2+k2_2/2,x3,x4,x5,x6)*dt;
    k4_2=f2(x1,x2+k3_2,  x3,x4,x5,x6)*dt;

    %% New values
    new1=x1+(k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
    new2=x2+(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;

end

%swapRow.m Swaps rows I (vector of row numbers) with rows J (also vector of
%row numbers to be swapped)
function A=swapRow(A,I,J)
A([J,I],:)=A([I,J],:);
end