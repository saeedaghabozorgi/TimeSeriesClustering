function dismatrix=Mtx_Distance(x,y,cond,dis_method,varargin)


options = struct('dtw_bound',0,'alphabet_size',0,'compression_ratio',0);
optionNames = fieldnames(options);
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive
    if any(strmatch(inpName,optionNames))
        options.(inpName) = pair{2};
    else
%         error('%s is not a recognized parameter name',inpName)
%     end
end


%compression_ratio: original_data_len / symbolic_len
if strcmp(cond,'same')
    data_n = length(x);
    dismatrix=zeros(data_n,data_n);
    for i = 1:data_n
        a=x{i};
        for j = i:data_n,
            if i ~= j
                b=x{j};
                switch dis_method
                    case 'Euclid'
                        dis=dis_euclidean(a,b);
                    case 'DTW'
                      %  dis=dis_dtw3(a,b,round(length(a)*options.dtw_bound));
                        [dis,~,~,~]=dis_dtw_complete(a,b);
                        
                    case 'LCSS'
                        dis=dis_lcs(a, b, 3, .3);
                    case 'SAXminDis'
                        dis=dis_min_dist(a, b, options.alphabet_size,options.compression_ratio);
                    case 'SAXmaxDis'
                        dis=dis_max_dist(a, b, options.alphabet_size,options.compression_ratio);
                    case 'SAXAPX'
                        dis=dis_SAX_apx(a, b, options.alphabet_size,options.compression_ratio);
                    otherwise
                        error(sprintf('DOCLUSTERING - unsupported algorithm "%s"',dis_method))
                end %of switch/case
                dismatrix(i, j) =dis;
            end
        end
    end
else
    dismatrix=zeros(length(x),length(y));
    for i = 1:length(x)
        a=x{i};
        for j = 1:length(y)
            b=y{j};
            if (length(a) ~= length(b))
                display('error: the strings must have equal length!');
                return;
            end
            switch dis_method
                case 'Euclid'
                    dis=dis_euclidean(a,b);
                case 'DTW'
                    dis=dis_dtw3(a,b,length(a)*options.dtw_bound);
                  %  [dis,~,~,~]=dis_dtw_complete(a,b);
                case 'LCSS'
                    dis=dis_lcs(a, b, 3, .3);
                case 'SAXminDis'
                    dis=dis_min_dist(a, b, options.alphabet_size,options.compression_ratio);
                case 'SAXmaxDis'
                    dis=dis_max_dist(a, b, options.alphabet_size,options.compression_ratio);
                case 'SAXAPX'
                    dis=dis_SAX_apx(a, b, options.alphabet_size,options.compression_ratio);
                otherwise
                    error(sprintf('DOCLUSTERING - unsupported algorithm "%s"',dis_method))
            end %of switch/case
            dismatrix(i, j) =dis;
            %t(j)=toc;

        end
    end
end
Nor = dismatrix - min( dismatrix(:) );
dismatrix = Nor / max( Nor(:) );

end