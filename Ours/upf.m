function [it2,F] = upf(Y,n,K,ca1,Ln)

    last2=zeros(n,K);
    it2=0;
    while any(last2~=Y)
        last2=Y;
        for i=1:n

            l=find(Y(i,:)==1);
            if (ca1(l)==1)
                continue;
            end
            Y1=Y;
            Y1(i,l)=0;
            Y2=Y1;
            obj=zeros(K,1);
            for j=1:K
                Y2(i,j)=1;
                obj(j)=(Y2(:,j)'*Ln*Y2(:,j))/(Y2(:,j)'*Y2(:,j))-(Y1(:,j)'*Ln*Y1(:,j))/(Y1(:,j)'*Y1(:,j));
            end
            r=find(obj==min(obj));
            if(r~=l)
                ca1(l)=ca1(l)-1;
                ca1(r)=ca1(r)+1;
                Y(i,l)=0;
                Y(i,r)=1;
            end
            
        end
        it2=it2+1;
    end
    F=Y*(Y'*Y)^(-1/2);
end