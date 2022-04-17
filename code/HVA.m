function [result,w] = HVA(signal,param)

[spec,param] = STFT(signal,param);
param.signalSpec = spec;
spec = whitening(spec);

[w,y,X,Xt,M,mu1,mu2,alpha,param] = initialization(spec,param);
for k = 1:param.iterNum
    wOld = w;
    yOld = y;
    w = proxLogDet(w - mu1*mu2*Xt(y), mu1);
    z = y + X(2*w - wOld);
    y = (1 - M(z)) .* z;
    y = alpha*y + (1-alpha)*yOld;
    w = alpha*w + (1-alpha)*wOld;
end

sep = backProjection(X(w),param);
result = iSTFT(sep,param);
result(size(signal,1)+1:end,:) = [];
end



%%% Local Functions %%%

function w = proxLogDet(w,mu)
for f = 1:size(w,3)
    [U,S,V] = svd(w(:,:,f));
    s = diag(S);
    s = (s + sqrt(s.^2 + 4*mu))/2;
    w(:,:,f) = U*diag(s)*V';
end
end

function M = HVAmask(x,lambda,kappa,fftn)
y = log(abs(x) + 1e-3);
meanY = mean(y,3);
y = y - meanY;
z = fft(y,fftn,3)/size(y,3);
M = min(1,abs(z)/lambda);
for n = 1:kappa
    M = (1-cos(pi*M))/2;
end
z = M .* z;
z = ifft(z,[],3,'symmetric')*size(y,3);
y = z(:,:,1:size(y,3));
y = y + meanY;
y = exp(2*y);
M = (y ./ sum(y)).^(1/size(y,1));
end

function y = filt3D(w,x,y)
for f = 1:size(w,3)
    y(:,:,f) = w(:,:,f) * x(:,:,f);
end
end

function w = filt3Dtranspose(y,x,w)
for m = 1:size(x,1)
    w(:,m,:) = sum(conj(x(m,:,:)).*y,2);
end
end



%%% Initialization %%%

function [w1,y1,X,Xt,mask,mu1,mu2,alpha,param] = initialization(x,param)
w1 = [eye(param.sourceNum) zeros(param.sourceNum,size(x,1)-param.sourceNum)];
w1 = repmat(w1,1,1,size(x,3)); % initial w

y1 = zeros(size(x)); % initial y (=0)

zeroW = zeros(size(w1)); zeroY = zeros(size(w1,1),size(x,2),size(x,3));

x = x/operatorNorm(x,w1,zeroW,zeroY); % normalized X

X  = @(w) filt3D(w,x,zeroY);          % multiplication of X
Xt = @(y) filt3Dtranspose(y,x,zeroW); % multiplication of Hermite transpose of X

mask = @(z) HVAmask(z,2*param.lambda,param.kappa,2^nextpow2(size(x,3)));

mu1 = param.mu1;
mu2 = param.mu2;
alpha = param.alpha;
end

function opNorm = operatorNorm(x,w,zw,zy)
opt.issym = true; opt.isreal = true;
XtX = @(z) reshape(filt3Dtranspose(filt3D(reshape(z,size(w)),x,zy),x,zw),[],1);
opNorm = sqrt(eigs(XtX,numel(w),1,'lm',opt));
end



%%% Audio BSS Tools %%%

function [X,param] = STFT(sig,param)
wLen = param.STFTwindowLength;
skip = param.STFTshiftSize;

x = [zeros(wLen,size(sig,2)); sig; zeros(wLen,size(sig,2))];
win = hann(wLen,'periodic');
param.window = win;

idx = (1:wLen)' + (0:skip:length(x));
x = [x; zeros(max(idx(:))-length(x),size(x,2))];
idx = idx + length(x)*reshape((0:size(x,2)-1),1,1,[]);

X = fft(win.*x(idx));
X = X(1:floor(wLen/2)+1,:,:);
X = permute(X,[3 2 1]);
end

function x = iSTFT(X,param)
win = param.window;
skip = param.STFTshiftSize;

X = ipermute(X,[3 2 1]);
win = calcCanonicalDualWindow(win,skip);
X = win.*ifft([X;zeros(size(X)-[2 0 0])],'symmetric'); [X1,X2,X3] = size(X);

vec = @(x) x(:);
X = sparse( ...
    vec((1:X1)' + skip*(0:X2-1) + (skip*X2+X1-skip)*reshape(0:X3-1,1,1,[])), ...
    vec(repmat(1:X2,X1,1,X3)), X(:) );
x = reshape(full(sum(X,2)),[],X3);
x(1:length(win),:) = [];
end

function dualWin = calcCanonicalDualWindow(win,skip)
dualWin = [win; zeros(lcm(skip,length(win))-length(win),1)];
dualWin = reshape(dualWin,skip,[]);
dualWin = dualWin./sum(abs(dualWin).^2,2);
dualWin = reshape(dualWin(1:length(win)),[],1);
end

function X = whitening(X)
for f = 1:size(X,3)
    [U,S,~] = svd(X(:,:,f),'econ');
    X(:,:,f) = U*((1./diag(S)).*(U'*X(:,:,f)));
end
end

function y = backProjection(y,param)
BP = @(y,s) ((s*y')/(y*y')).'.*y;
for f = 1:size(y,3)
    y(:,:,f) = BP(y(:,:,f),param.signalSpec(param.refMicIndex,:,f));
end
end