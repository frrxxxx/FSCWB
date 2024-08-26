function a = calRcut(X,K, W,label)
D=diag(sum(W,2));
L=D-W;

a=0;
Y = label2binary(label);
for i=1:K
    tmp=((Y(:,i))'*L*Y(:,i))/((Y(:,i))'*Y(:,i));
    a=a+tmp;
end


end