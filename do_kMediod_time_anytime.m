function [c,itr]= do_kMediod_time_anytime (nor_traj_raw,k,isRand,ctmp,p,varargin)

options = struct('alphabet_size',0,'compression_ratio',0,'weight',[],'rep','RAW');
optionNames = fieldnames(options);
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive
    if any(strmatch(inpName,optionNames))
        options.(inpName) = pair{2};
        % disp(['  ', num2str(inpName),' : ',num2str(options.(inpName))]);
        %    else
        %       error('%s is not a recognized parameter name',inpName)
    end
end
% representation
nor_traj=represent_TS(nor_traj_raw,options.rep,varargin{:});
itr=1;
Dist=Mtx_Distance(nor_traj,nor_traj,'same','Org','dis_method','DTW','dtw_bound',0.1,'alphabet_size',0,'compression_ratio',0);
start=length(Dist)-1;
for i=start:-1:1
    for j=i+1:length(Dist)
        Dist(i,j)=dis_dtw3(nor_traj{i},nor_traj{j},length(nor_traj{i}));
        Dist(j,i)=Dist(i,j);
    end
    Dis2=squareform(Dist);
    [c, cost, medoids] = do_kMedoids_keogh(k,Dis2);
    %evaluation
     for jj=1:k
        l3_mems=find(c==jj);
        for ii=1:length(l3_mems)
            sub_members=find(ctmp(:,1)==l3_mems(ii));
            ctmp(sub_members,2)=jj;
        end
     end
    ccc=ctmp(:,2);
    [SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality]= do_Evaluate(p,ccc,nor_traj,[],[]);
    disp(['  --> i:', num2str(i), ' | ','incremental quality:',num2str(quality)]);
end

end