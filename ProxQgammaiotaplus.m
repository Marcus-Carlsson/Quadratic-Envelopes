
function x = ProxQgammaiotaplus(y,k,gamma,rho) %k is the cardinality, r=rho and g=gamma, y the point in which we compute prox

%this script computes the proximal operator of Q_gamma(iota_k^+) as
%explained in Sec 4.2 of the paper On phase retrieval via matrix completion and the estimation of low rank
%PSD matrices


[n1,n2]=size(y);if n2~=1, y=y';end;%if we get a row vector, make it column vector

n=length(y);

[ysort,id]=sort(y,'descend');
idinv(id)=[1:n];

%we first compute prox of the S transform

%1st step
if ysort(k)<0
    tildek=max(find(ysort>=0));
    x=[ysort(1:tildek);rho/gamma*ysort(tildek+1:n)];
    x=x(idinv);
elseif (rho/gamma)*ysort(k+1)<ysort(k)
    x=[ysort(1:k);rho/gamma*ysort(k+1:n)];
    x=x(idinv);

else    
    %2nd step: We compute the indices called l* and j* in the paper; we call them just l and j.
    xnew=[ysort(1:k);rho/gamma*ysort(k+1:n)];%this vector will be the new radii, to be updated before the end
    temp=find(xnew<=xnew(k+1));j=min(temp);
    temp=find(xnew>=xnew(k));l=max(temp);
    
    z=sort(xnew(j:l),'descend'); %3rd step. We compute the vector z.
    for m=1:length(z)-1
        s=(z(m)+z(m+1))/2;
        %4th step. We compute the indices called l and j in the paper. We call them l1 and j1.
        temp=find(xnew<=s);j1=min(temp);
        temp=find(xnew>=s);l1=max(temp);
        sI=(rho*sum(ysort(j1:l1)))/((k+1-j1)*rho + (l1-k)*gamma);%5th step
        if (sI>=z(m+1) && z(m)>=sI)%6th step
            x=[max(xnew(1:k),sI);min(xnew(k+1:n),sI)];
            x=x(idinv);
            break
        end
    end    
end
    x=(rho*y-gamma*x)/(rho-gamma);%this step finally takes the Sprox and gives the Qprox, via formula (20) in the paper
    
if n2~=1, x=x';end;%switch back to row vector in case this is what came in


