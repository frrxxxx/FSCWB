close all;
clear
dataname="2d-3c-no0"
X_raw = load("../Datasets/"+dataname+".txt");
Color = load("../Datasets/"+dataname+"_color.txt");
K =3;
sigma=1.0;
n=size(X_raw,1)
[P, v] = SensCNVT(Color);

% Standardization and Normalizationelliptical
B=repmat(v',n,1);
X = normalize(X_raw, 1);
X = X./repmat(sqrt(sum(X.^2,2)),1, size(X,2));
distM=squareform(pdist(X));
tmp = distM.^2 / (2 * sigma );
W = exp(-tmp);
W(logical(eye(size(W))))=0;
%labels = spectralcluster(W,K,'Distance','precomputed','LaplacianNormalization','symmetric');
%labels=kmeans(X,K);
labels = randsrc(n,1,1:K);
%labels=Color;
%Draw(X_raw,labels);

p1 =8;
p2=8;
tic;
t1=clock;
[it,label,B] = afscb(labels, K,W, P,B, p1,p2);
 %[it,label] = facd( labels, K,W, P,v, p1);
t2=clock;
time=etime(t2,t1)
toc;
Draw(X_raw,label);

% filename='elliptical_AFSC.txt';
% dlmwrite("elliptical_AFSC.txt",label)
Y = label2binary(label);

% YY = Y'*Y;
% PYYY = P'*Y/YY;
% AW = 0;
% for i = 1:size(PYYY, 2)
%     tmp = YY(i, i) * ws_distance(PYYY(:, i), v);
%     AW = AW + tmp;
% end
it;
% filename=dataname+'_AFSC.txt';
% dlmwrite(dataname+"_AFSC.txt",label);


BAL = calcBAL(label, Color, K)';
bal=min(BAL)
wbal=calwBal(label,Color,K)
%AW = AW/size(X,1);
%DI = dunns(K, distM, label);
% SSE = calcSSE(X', label)
ncut=calNcut(K,W,label)