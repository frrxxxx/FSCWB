
close all;
clear
dataname="subsampled_adult"
X_raw = load("../Datasets/"+dataname+".txt");
Color = load("../Datasets/"+dataname+"_color.txt");
k =10;
sigma=1.0;
n=size(X_raw,1)
[P, v] = SensCNVT(Color);
% Standardization and Normalizationelliptical
X = normalize(X_raw, 1);
X = X./repmat(sqrt(sum(X.^2,2)),1, size(X,2));
distM=squareform(pdist(X));
tmp = distM.^2 / (2 * sigma );
W = exp(-tmp);
W(logical(eye(size(W))))=0;

% edges = readmatrix('lastfm_clean_edges.csv');
% m = max(edges,[],'all') + 1;
% A = zeros(m,m);
% for k = 1:size(edges,1)
%     i = edges(k,1) + 1;
%     j = edges(k,2) + 1;
%     A(i,j) = 1;
% end
% A = (A + A')/2;
% 
% G = graph(A);
% 
% [bin,binsize] = conncomp(G);
% idx = binsize(bin) == max(binsize);
% SG = subgraph(G, idx);
% 
% 
% W = adjacency(SG);
% n = size(W, 1);
% 
% 
D = diag(W*ones(n,1));
% 
% 
% countries_raw = readmatrix('lastfm_clean_countries.csv');
% 
% countries = countries_raw(idx,:);
% f = countries(:,2);
% uniquef = unique(f);
% h = length(uniquef);
sensitive=Color;
sens_unique=unique(sensitive);
h = length(sens_unique);
sens_unique=reshape(sens_unique,[1,h]);

sensitiveNEW=sensitive;

temp=1;
for ell=sens_unique
    sensitiveNEW(sensitive==ell)=temp;
    temp=temp+1;
end
    
F=zeros(n,h-2);

for ell=1:(h-2)
    temp=(sensitiveNEW == ell);
    F(temp,ell)=1; 
    groupSize = sum(temp);
    F(:,ell) = F(:,ell)-groupSize/n;
end



% F = zeros(n, h-2);
% G = zeros(n, h-1);
% 
% for i = 1:h-1
%     idx = f == uniquef(i);
%     count = sum(idx);
%     fprintf("size of group # %d\n", count);
%     if i < h-1
%         F(:,i) = idx - count/n;
%     end
%     G(:,i) = idx;
% end




% fprintf('----------SC------------\n');
% labels1 = alg1(W, D, k);
% fractions1 = computeFraction(labels1, G);

% fprintf('----------FairSC------------\n');
% labels2 = alg2(W, D, F, k);
% fractions2 = computeFraction(labels2, G);

fprintf('----------s-FairSC------------\n');
tic;
labels3 = alg3(W, D, F, k);
toc;



% filename=dataname+'_sFSC.txt';
% dlmwrite(dataname+"_sFSC.txt",labels3)
% fractions3 = computeFraction(labels3, G);
% Draw(X_raw,labels3);
% BAL = calcBAL(labels3, Color, k)'
% bal=min(BAL)
% ncut=calNcut(k,W,labels3)
% wbal=calwBal(labels3,Color,k)