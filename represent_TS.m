function nor_traj=represent_TS(nor_traj_raw,rep,varargin)


options = struct('alphabet_size',0,'compression_ratio',0);
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
%         error('%s is not a recognized parameter name',inpName)
     end
end

switch rep
    %-----------------------------------------------------------------
    case 'SAX'
        alphabet_size=options.alphabet_size;
        compression_ratio=options.compression_ratio;
        data_len=(floor(size(nor_traj_raw{1},2)/compression_ratio))*compression_ratio;
        nseg = data_len/compression_ratio;
        for i=1:length(nor_traj_raw)
            sym=rep_SAX(nor_traj_raw{i}, data_len, nseg, alphabet_size);
            nor_traj{i}=sym(1,:) ;
        end
        %-------------------------------------------------------------
    case 'RAW'
        nor_traj = nor_traj_raw;
        %-------------------------------------------------------------
    otherwise
        error(sprintf('DOCLUSTERING - unsupported algorithm "%s"',rep_method))
end %of switch/case
end
