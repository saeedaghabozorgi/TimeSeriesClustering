function dist = dis_SAX_apx(str1, str2, alphabet_size, compression_ratio)

%-----demo
% str1= [ 1 2 2];
% str2= [1 2 3];
% alphabet_size=4;
% compression_ratio=3;
%---------

% approximate distance using the mean of each word
if (length(str1) ~= length(str2))
    display('error: the strings must have equal length!');
    return;
end

if (any(str1 > alphabet_size) | any(str2 > alphabet_size))
    display('error: some symbol(s) in the string(s) exceed(s) the alphabet size!');
    return;
end

dist_matrix = build_dist_table(alphabet_size);

dist = 0;
tt=dist_matrix(str1,str2);
pp=diag(tt);
dist = sqrt(compression_ratio * sum(pp));

function dist_table = build_dist_table2(alphabet_size)
steps=8/alphabet_size;
cutlines=[((alphabet_size/2)-1)*(-1)*steps:steps:((alphabet_size/2)-1)*steps];

meanline=[((alphabet_size/2)-1)*(-1)*steps-(steps/2):steps:((alphabet_size/2)-1)*steps+(steps/2)];

dist_table=zeros(alphabet_size,alphabet_size);

for i = 1 : alphabet_size
    
    % the min_dist for adjacent symbols are 0, so we start with i+2
    for j = 1 : alphabet_size
        
        % square the distance now for future use
        dist_matrix(i,j)=(meanline(i)-meanline(j))^2;
        
        % the distance matrix is symmetric
        dist_table(j,i) = dist_table(i,j);
    end;
end;

function dist_table = build_dist_table(alphabet_size)

switch alphabet_size
    case 2, cutlines  = [0];
    case 3, cutlines  = [-0.43 0.43];
    case 4, cutlines  = [-0.67 0 0.67];
    case 5, cutlines  = [-0.84 -0.25 0.25 0.84];
    case 6, cutlines  = [-0.97 -0.43 0 0.43 0.97];
    case 7, cutlines  = [-1.07 -0.57 -0.18 0.18 0.57 1.07];
    case 8, cutlines  = [-1.15 -0.67 -0.32 0 0.32 0.67 1.15];
    case 9, cutlines  = [-1.22 -0.76 -0.43 -0.14 0.14 0.43 0.76 1.22];
    case 10, cutlines = [-1.28 -0.84 -0.52 -0.25 0. 0.25 0.52 0.84 1.28];
    otherwise, disp('WARNING:: Alphabet size too big');
end;
switch alphabet_size
    case 1, meanline  = [0];
    case 2, meanline  = [-0.43 0.43];
    case 3, meanline  = [-0.67 0 0.67];
    case 4, meanline  = [-0.84 -0.25 0.25 0.84];
    case 5, meanline  = [-0.97 -0.43 0 0.43 0.97];
    case 6, meanline  = [-1.07 -0.57 -0.18 0.18 0.57 1.07];
    case 7, meanline  = [-1.15 -0.67 -0.32 0 0.32 0.67 1.15];
    case 8, meanline  = [-1.22 -0.76 -0.43 -0.14 0.14 0.43 0.76 1.22];
    case 9, meanline = [-1.28 -0.84 -0.52 -0.25 0. 0.25 0.52 0.84 1.28];
    case 10, meanline = [ -1.34 -0.91 -0.6 -0.35 -0.11 0.11 0.35 0.6 0.91 1.34];
    otherwise, disp('WARNING:: Alphabet size too big');
end;

dist_table=zeros(alphabet_size,alphabet_size);

for i = 1 : alphabet_size
    
    % the min_dist for adjacent symbols are 0, so we start with i+2
    for j = 1 : alphabet_size
        
        % square the distance now for future use
        dist_table(i,j)=(meanline(i)-meanline(j))^2;
        
        % the distance matrix is symmetric
        dist_table(j,i) = dist_table(i,j);
    end;
end;
