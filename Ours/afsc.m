function [it,label,B] = afsc(label, K,S, P,B, p1)
m=size(P,2);
D=diag(sum(S,2));
L=D-S;
n=size(L,1);
Y = label2binary(label);
ca1=sum(Y,1);%cluster size
maxit = 1000;
last = zeros(n,1);
it=0;
while any(last~=label)&&it<=maxit
    last=label;
    Y=label2binary(label);
    %update F
    Ln=(P-2*B)*P';
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
            o1=(Y2(:,j)'*L*Y2(:,j))/(Y2(:,j)'*D*Y2(:,j))-(Y1(:,j)'*L*Y1(:,j))/(Y1(:,j)'*D*Y1(:,j));
            o2=(Y2(:,j)'*Ln*Y2(:,j))/(Y2(:,j)'*Y2(:,j))-(Y1(:,j)'*Ln*Y1(:,j))/(Y1(:,j)'*Y1(:,j));
            obj(j)=o1+p1*o2/n;
        end
        r=find(obj==min(obj));
        if(r~=l)
            ca1(l)=ca1(l)-1;
            ca1(r)=ca1(r)+1;
            Y(i,l)=0;
            Y(i,r)=1;
        end
    [~, label] = max(Y, [], 2);
    end
    %update B
%    U=F*F'*P;
%     for i=1:n
%         o=ones(m,1);
%         v=U(i,:)+o'/n-(U(i,:)*o*o')/n;
%         d=newton_method2(m,n,v);
%         for j=1:m
%             B(i,j)=max(v(j)-d,0);
%         end
%     end
%     for z=1:n
%         label(z)=find(F(z,:)~=0);
%     end
    it=it+1; 
end

end
