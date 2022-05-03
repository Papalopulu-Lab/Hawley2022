function [Y_RAW,t_concat]=stackCols(Pnew,K,rows,t,ti)

%This takes the m=cells by n=time matrix that stores the protein expression
%over time and splits every full column and places it to the right
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

end
