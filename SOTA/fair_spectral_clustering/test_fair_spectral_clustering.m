clear
X_raw = load("../Datasets/2d-2c-no0.txt");
Color = load("../Datasets/2d-2c-no0_Color.txt");
K =2;

% Normalization and Standardization
X = normalize(X_raw, 1);
X = X./repmat(sqrt(sum(X.^2,2)),1, size(X,2));

% Contruct Similarity Matrix
distM=squareform(pdist(X));
sigma=1.0;
tmp = distM.^2 / (2 * sigma );
W = exp(-tmp);
W(logical(eye(size(W))))=0;

tic;label = Fair_SC_unnormalized(W, K, Color);toc;
% label = load('bank_FSCN.txt');
bal = calcBAL(label, Color, K);
min(bal)
ncut=calNcut(K,W,label)
wbal=calwBal(label,Color,K)
[unique_labels, n_clusters] = unique_count(label);
Draw(X_raw,label)
%rcut=calRcut(X,K,W,label)

%tic;label2 = Fair_SC_normalized(W, K, Color);toc;
% label = load('bank_FSCN.txt');

%bal2 = calcBAL(label2, Color, K);
%min(bal2)
%sse2 = calcSSE(X', label2)
%ncut2=calNcut(X,K,W,label2)

function [unique_a, cnt_unique] = unique_count(a)
unique_a = unique(a);
cnt_unique = [];
for i = 1:length(unique_a)
    IDX = a == unique_a(i);
    cnt_unique = [cnt_unique, sum(IDX)];
end
end

