function [it,label,B] = fast_afsc(label, K,S, P,B, p1)

D=diag(sum(S,2));
L=D-S;
n=size(L,1);
Y = label2binary(label);
ca1=diag(Y'*L*Y);
ca2=diag(Y'*D*Y);
ca3=Y'*L;
ca4=sum(Y,1);%cluster size
Ln=(P-2*B)*P';
ca5=diag(Y'*Ln*Y);
ca6=Y'*Ln;
ca7=Ln*Y;

maxit = 1000;
last = zeros(n,1);
it=0;
while any(last~=label)&&it<=maxit
    last=label;
    Y=label2binary(label);

    for i=1:n
        nca1=zeros(K,1);
        nca2=zeros(K,1);
        l=find(Y(i,:)==1);
        nca5=zeros(K,1);
        if (ca4(l)==1)
            continue;
        end
        Y1=Y;
        Y1(i,l)=0;
        Y2=Y1;
        obj=zeros(K,1);
        for j=1:K
            Y2(i,j)=1;
            if(j==l)
                o1=ca1(j)/ca2(j)-(ca1(j)-2*ca3(j,i)+L(i,i))/(ca2(j)-D(i,i));
                o2=ca5(j)/ca4(j)-(ca5(j)-ca6(j,i)-ca7(i,j)+Ln(i,i))/(ca4(j)-1);
                obj(j)=o1+p1*o2/n;
                nca1(j)=ca1(j)-2*ca3(j,i)+L(i,i);
                nca2(j)=ca2(j)-D(i,i);
                nca5(j)=ca5(j)-ca6(j,i)-ca7(i,j)+Ln(i,i);
            end
            if(j~=l)
                o1=(ca1(j)+2*ca3(j,i)+L(i,i))/(ca2(j)+D(i,i))-ca1(j)/ca2(j);
                o2=(ca5(j)+ca6(j,i)+ca7(i,j)+Ln(i,i))/(ca4(j)+1)-ca5(j)/ca4(j);
                obj(j)=o1+p1*o2/n;
                nca1(j)=ca1(j)+2*ca3(j,i)+L(i,i);
                nca2(j)=ca2(j)+D(i,i);
                nca5(j)=ca5(j)+ca6(j,i)+ca7(i,j)+Ln(i,i);
            end


        end
        r=find(obj==min(obj));
%       update
        if(r~=l)
            Y(i,l)=0;
            Y(i,r)=1;

            ca1(r)=nca1(r);
            ca2(r)=nca2(r);
            ca5(r)=nca5(r);
            ca4(r)=ca4(r)+1;
            for e =1:n
                ca3(r,e)=ca3(r,e)+L(i,e);
            end
            for e =1:n
                ca6(r,e)=ca6(r,e)+Ln(i,e);
            end
            for e =1:n
                ca7(e,r)=ca7(e,r)+Ln(e,i);
            end

            ca1(l)=nca1(l);
            ca2(l)=nca2(l);
            ca5(l)=nca5(l);
            ca4(l)=ca4(l)-1;
            for e =1:n
                ca3(l,e)=ca3(l,e)-L(i,e);
            end
            for e =1:n
                ca6(l,e)=ca6(l,e)-Ln(i,e);
            end           
            for e =1:n
                ca7(e,l)=ca7(e,l)-Ln(e,i);
            end

        end
    [~, label] = max(Y, [], 2);
    end
  
    it=it+1; 
end

end
