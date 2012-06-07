function [c,Zh]=do_Hierarchical_time(nor_traj_raw,k,dis_method,linkage_Method,Cutoff,varargin)

options = struct('alphabet_size',0,'compression_ratio',0,'rep','RAW');
optionNames = fieldnames(options);
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1}); %# make case insensitive
   if any(strmatch(inpName,optionNames))
      options.(inpName) = pair{2};
%    else
%       error('%s is not a recognized parameter name',inpName)
    end
end


% representation
nor_traj=represent_TS(nor_traj_raw,options.rep,varargin{:});

if nargin<1,     end
if nargin<2,    k=1;        end

Rows=length(nor_traj);
if k>Rows
    return;
end

dis=Mtx_Distance(nor_traj,nor_traj,'same',dis_method,varargin{:});
dis= squareform(dis+dis');
Zh = linkage(dis,linkage_Method);

if Cutoff==-1
    c = cluster(Zh,'maxclust',k);
else
    c = cluster(Zh,'cutoff',cutoff);
end
end




