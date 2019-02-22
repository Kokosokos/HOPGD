% =======================================
% By I. Kriuchevskyi kruchevskyivan@gmail.com,
% reworked algorithm of Ye Lu, ye.lu@insa-lyon.fr 
%
% X: snapshots matrix (any dimensionality >= 2)
% Model: separated modes
% e: error
% =======================================
function [ Model e] = Copy_2_of_HOPGD_6( X, ec,nbrModes)

nmax=nbrModes;
vmax=2400;
nd=length(size(X)); % dimension of X
dimArr=1:nd;

eBc=1e-6;
K=size(X,nd);
f=zeros(size(X));
fna=0;
B.F=cell(1,nd);
idx=cell(1,nd);
Model=cell(1,nd);
bta1=zeros(1,nd);
br1=ones(1,nd);
bt1=zeros(1,nd);
xSize=size(X);

xSizeTot=1;

for m=1:nd
    B.F{m}=ones(xSize(m),1);
    idx{m}=1:xSize(m);
    xSizeTot = xSizeTot*xSize(m);
end

norm_X=sqrt(sum(sum(X.*X)));

%n-number of modes !!!!
for n=1:nmax
    R=X-fna;
    norm_R=sqrt(sum(sum(R.*R)));
    
    e=norm_R./norm_X;
    %Loop over additional dimensions (exclude first two, space and time)
    for m=3:nd
        e=max(e);
    end
    if e<ec
        break;
    end
    
    for v=1:vmax
        e1=1;
        %Loop over dimensions
        for m=1:nd
            permId=dimArr;
            permId(1)=m;
            permId(m)=1;
            RR=permute(R,permId);
            tempkron=1;
            normBF=1;
            for mm=permId
                if mm ~= m
                    tempkron=kron(B.F{mm},tempkron);
                    normBF = normBF*(B.F{mm}'*B.F{mm});
                end
            end
            M=sum(reshape(RR,xSize(m),xSizeTot/xSize(m))'.*tempkron)';
            B.F{m}=M/normBF ;
            bt1(m)=norm(B.F{m});
            %                     norm(bt1(m)-bta1(m))/br1(m)
            if (v==0)
                br1(m)=bt1(m);
            end
            
            e1=e1&&(norm(bt1(m)-bta1(m))/br1(m)<eBc);
        end
        if e1
            f=B.F{1}*B.F{2}';
            for mm=3:nd
                f = bsxfun(@times,f,reshape(B.F{mm},[ones(1,mm-1) numel(B.F{mm})]));
            end
            for mm=1:nd
                Model{mm}=[Model{mm} B.F{mm}];
            end
            fna=fna+f;
            bta1=zeros(1,nd);
            for m=2:nd
                B.F{m}=ones(xSize(m),1);
            end
            break
        else
            bta1=bt1;
        end
    end
    if v==vmax
        disp('no convergence')
    end
end
end



