function [hi, lo, DcCells]=HiLoLoop(pop)
    hi=[];
    lo=[];
    cells=size(pop,1);
    for i=1:cells
        [B, N, IB]=RunLength_M(pop(i,:));
    %     B=B(2:end-1); %Removes clipped high and low regions at the start and end of signal
    %     N=N(2:end-1); %Removes clipped high and low regions at the start and end of signal
        hi=[hi; N(B==1)];
        lo=[lo; N(B==-1)];

        alphaHi=sum(N(B==1))/(size(pop,2)); %Proportion of the time spent in high state
        alphaLo=sum(N(B==-1))/(size(pop,2)); %Proportion of the time spent in low state
        DcCells(i)=2*(1-max([alphaHi alphaLo]));
    end