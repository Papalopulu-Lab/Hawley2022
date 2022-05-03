function [NM,NumNeighbours]=neighbours(rows,cols,boundary,VertProtrusions,HorzProtrusions,VertProtStrength,HorzProtStrength,StrengthGrad,numDistalNeigh,noProximal)
% neighbours.m is a function that outputs a 'neighbour matrix'(NM). 
% Each row of the NM is used to extract the neighbouring cell values from
% the vector it is being multiplied by (e.g. the protein concentration
% vector). Each row corresponds to it's cell number i.e. row 2 describes
% cell 2's neighbours (using 1's and 0's).

% Rows and cols refer to number of rows and columns in the hexagonal grid
% that this function is assuming. Boundary=0 creates a NM with hard
% boundaries, boundary=1 creates a periodic boundary NM, boundary=2 creates
% a 'soft' boundary where the edge cells with missing neighbours will have
% increased influence/weighting from it's 3 or 4 immediate neighbours
% (rather than periodic neighbours).

N=6; % Number of immediate neighbours
cells=rows*cols;
NM=zeros(cells,cells);
NumNeighbours=zeros(cells,1);

n=0; % Keeps track of row number of current cell
m=0; % Keeps track of column number of current cell

cellVec=1:cells;
cellArr=reshape(cellVec,[rows,cols]);

for cell=1:cells
%% Row and column counter
    if n==rows
        n=0;
    end
    
    if n==0
        m=m+1;
    end
    
    n=n+1;

%% Explanation of how neighbours.m function works
    % Cells in a grid are numbered as follows for a n=3,m=3 matrix:
    % 1  4  7
    % 2  5  8
    % 3  6  9
    
    %Notation used in variable 'compass' below refers to the following:
    
    % (NW)                (N)                (NE)
    %     cell-rows-1   cell-1   cell+rows-1
    %                     |
    %                     |
    % (W) cell-rows ---- cell -- cell+rows   (E)  
    %                     |
    %                     |
    %     cell-rows+1   cell+1   cell+rows+1
    % (SW)               (S)                 (SE)
    
%% Compass generates the indices that define the current cell's neighbours
    %         1=N       2=NE         3=E        4=SE       5=S       6=SW         7=W        8=NW
    compass=[cell-1; cell+rows-1; cell+rows; cell+rows+1; cell+1; cell-rows+1; cell-rows; cell-rows-1];

%% Accounting for inccorect/non-existent compass values at the boundaries
    if n==1 %If on upper boundary, remove following entries:
        compass(1)=0; %N
        compass(2)=0; %NE
        compass(8)=0; %NW
    end
    
    if n==rows % If on lower boundary, remove following entries:
        compass(5)=0; %S
        compass(4)=0; %SE
        compass(6)=0; %SW
    end
    
    if m==1 % If on left boundary, remove following entries:
        compass(7)=0; %W
        compass(8)=0; %NW
        compass(6)=0; %SW
    end
    
    if m==cols % If on right boundary, remove following entries:
        compass(3)=0; %E
        compass(2)=0; %NE
        compass(4)=0; %SE
    end
    
%% Accounting for shifted odd and even rows in hexagonal geometry
    if mod(n,2)==0 % If on even row, select approriate neighbours
        compass=[compass(1); compass(2); compass(3); compass(4); compass(5); compass(7)];

    else % If on odd row
        compass=[compass(1); compass(3); compass(5); compass(6); compass(7); compass(8)];

    end
%% Periodic boundary conditions
    if boundary==1
        
        if mod(rows,2)~=0
            error('When using periodic boundary conditions, number of rows must be even! (Due to hexagonal geometry and such)')
        end
        
        p=zeros(6,1);
%________________________________________________________________________

        if n==1 %If on upper boundary, make periodic with 
            p(1)=cell-1;
            p(4)=cell+rows-1;
        end
%________________________________________________________________________

        if n==rows % If on lower boundary, make periodic with 
            p(3)=cell-rows+1;
            p(6)=cell+1;
        end
%________________________________________________________________________

        if m==1 % If on left boundary
            
            if mod(n,2)==1 % If on odd row
                p(2)=cell+(cols-1)*rows;
                p(1)=p(2)-1;
                p(3)=p(2)+1;
                
                if n==1 && m==1 % If at top left
                    p(1)=cols*rows;
                end
                
            else % If on even row
                p(2)=cell+(cols-1)*rows;
            end   
        end
%________________________________________________________________________
        
        if m==cols % If on right boundary, make periodic with 
            
            if mod(n,2)==0 % If on even row
                p(5)=cell-(cols-1)*rows;
                p(4)=p(5)-1;
                p(6)=p(5)+1;

                if n==rows && m==cols % If at bottom right
                    p(6)=1;
                end
                
            else
                p(5)=cell-(cols-1)*rows;
            end
        end 
%________________________________________________________________________
        compass=[compass; p]; % Add periodic boundary conditions into compass array  
    end
    
%% Vertical Protrusions
%     if VertProtrusions>1
%         % VertProtrusions is the number of cells away that a protrusion can reach
% %         cell
%         for pp=2:VertProtrusions
%             
%             j=pp-1;
%             prot_neigh(2*j-1:2*j,1)=[cell-pp; cell+pp];
% %             prot_neigh
% %             
% %             pp
% %             n
%             if n-pp<1
%                 prot_neigh(2*j-1,1)=0;
%             end
%             
%             if pp+n>rows
%                 prot_neigh(2*j,1)=0;
%             end
%             
%         end
% %         prot_neigh
%         compass=[compass; prot_neigh];
%     end
    
if VertProtrusions>1
    vertH0=cellArr(:,m)'; %Cells on the same row as recieving cell
    rankV=1:rows; rankV=abs(rankV-n);
  
    vertH1=[vertH0-rows;vertH0+rows];
    
    vertH0=vertH0(rankV>1 & rankV<=VertProtrusions);
    vertH1=vertH1(:,rankV>1 & rankV<=VertProtrusions);
    rankVW=rankV(:,rankV>1 & rankV<=VertProtrusions);
    
    vertH1(vertH1<1 | vertH1>cells)=[];
%     cell
%     vertH1
end  
    
%% Horizontal Protrusions
if HorzProtrusions>1
    horzV0=cellArr(n,:); %Cells on the same row as recieving cell
    
    rankH=1:cols; rankH=abs(rankH-m); %rankH indexes cells on the same row with the distance it is from the current cell
    rankH1=rankH-m;
    
    if n==1 
        horzV1=horzV0+1;
    elseif n==rows
        horzV1=horzV0-1;
    else
        horzV1=[horzV0-1; horzV0+1];
    end  
    
%     if noProximal == 1
%         proximalV0=horzV0(rankH==1);
%     end
    horzV0=horzV0(rankH>1 & rankH<=HorzProtrusions);    
    horzV1=horzV1(:,rankH>1 & rankH<=HorzProtrusions);
    rankHW=rankH(:,rankH>1 & rankH<=HorzProtrusions);
    
    
% % % %     %Make neighbours assymetric - move right hand side neighbour move up by 1
% % % %     %1
% % % %     if n==1
% % % %         horzV0(horzV0>cell)=[];
% % % %     else
% % % %         horzV0(horzV0>cell)=horzV0(horzV0>cell)-1;
% % % %     end
    
% % %     %Make neighbours assymetric - move left and right hand side neighbour move up by 1
% % % %     if n==1
% % % %         horzV0=[];
% % % %     else
% % % %         horzV0=horzV0-1;
% % % %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%     if mod(n,2)==1
%         right=horzV1>cell;
%         horzV1(right)=horzV1(right)-rows;
%     elseif mod(n,2)==0
%         left=horzV1<cell;
%         horzV1(left)=horzV1(left)+rows;
%     end
  
    %Change top right cell to right (break symmetry)
%     horzV1
%     if mod(n,2)==1
%         horzV1(horzV1>cell+2 & horzV1<cell+rows)=cell+rows;
%     elseif mod(n,2)==0
%         horzV1(horzV1>cell+2 & horzV1<cell+2*rows)=cell+2*rows;
%     end
%     horzV1
end  
%% Final processing of indices
    keep=compass>0 & compass<=cells; % Ensures that index value less than zero or greater than number of cells are not included
    
    NI=compass(keep); %neighbour index for immediate neighbours

    
    NV=zeros(1,cells);
    if noProximal == 0
        NV(NI)=1; % Neighbour Vector (NV) (for this loop) for immediate neighbours
    end
    if HorzProtrusions>1
        if StrengthGrad==1
%             NV(horzV0)=weight(rankHW,HorzProtrusions,HorzProtStrength);
%             NV(horzV1)=repmat(weight(rankHW,HorzProtrusions,HorzProtStrength),[size(horzV1,1),1]);
        elseif StrengthGrad==0
                
            if numDistalNeigh == 4
                
                NV(horzV0)=HorzProtStrength;
%                 if noProximal==1
%                     NV(proximalV0)=0;
%                 end
                
            elseif numDistalNeigh == 6
                
                NV(horzV1)=HorzProtStrength;
%                 if noProximal==1
%                     NV(proximalV0)=0;
%                 end
                
            elseif numDistalNeigh == 12
                
                NV(horzV0)=HorzProtStrength;
                NV(horzV1)=HorzProtStrength;
%                 if noProximal==1
%                     NV(proximalV0)=0;
%                 end
                
            end
        end
         
%         %P1 adds a row above and below

    end
    
    if VertProtrusions>1
        if StrengthGrad==1
%             NV(vertH0)=weight(rankVW,VertProtrusions,VertProtStrength);
%             NV(vertH1)=repmat(weight(rankVW,VertProtrusions,VertProtStrength),[size(vertH1,1),1]);
        elseif StrengthGrad==0
            if numDistalNeigh==2
                
                NV(vertH0)=VertProtStrength;
                
            elseif numDistalNeigh==4
                
                NV(vertH0)=VertProtStrength;
                
            elseif numDistalNeigh==6
                
                NV(vertH0)=VertProtStrength;
                
            elseif numDistalNeigh==12
                
                NV(vertH0)=VertProtStrength;
                NV(vertH1)=VertProtStrength;
                
            end
        end

    end
    
    
    if boundary==2 && sum(NV)<N
        NV=NV.*N./sum(NV);
    end
%     cell
%     find(NV)
    NM(cell,:)=NV; % Neighbour matrix
    NumNeighbours(cell)=sum(NV>0);
end   
end

function [W]=weight(rank,HP,HPS)
    HP=HP;
    m=-1/HP;
    c=1+1/HP;
    W=HPS.*(m*rank+c);
end
