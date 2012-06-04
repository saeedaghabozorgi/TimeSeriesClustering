function dis = dis_lcs(a, b, delta, epsilon, transpose)

% some checks
setTranspose = 1;
if nargin == 0, demo; return; end
if nargin < 4 ,
    disp('similarity = lcsMatching(a, b, delta, epsilon, transpose)');
end
if nargin < 5,
    setTranspose = 0;
end


m = length(a);
n = length(b);

% put the shorter first
if n<m,
    temp = a;
    a = b;
    b = temp;
    m = length(a);
    n = length(b);
end

lcstable = zeros(m+1, n+1);
prevx = zeros(m+1, n+1);
prevy = zeros(m+1, n+1);
% Find LCS using dynamic programming
for i=1:m,
    for j = (i-delta):1:(i+delta),
        if (j<=0 || j>n); %#ok<ALIGN>
        elseif ( (b(j)+epsilon)>=a(i) && (b(j)-epsilon)<=a(i))
			lcstable(i+1,j+1) = lcstable(i,j)+1;
			prevx(i+1,j+1) = i;
			prevy(i+1,j+1) = j;
		elseif lcstable(i,j+1) > lcstable(i+1,j),
			lcstable(i+1,j+1) = lcstable(i,j+1);
			prevx(i+1,j+1) = i;
			prevy(i+1,j+1) = j+1;
		else
			lcstable(i+1,j+1) = lcstable(i+1,j);
			prevx(i+1,j+1) = i+1;
			prevy(i+1,j+1) = j;
		end 
	end
end

% Get rid of initial conditions
lcstable = lcstable(2:end, 2:end);
prevx = prevx(2:end, 2:end)-1;
prevy = prevy(2:end, 2:end)-1;

% ====== LCS similarity 
[lcs, pos]= max(lcstable(m, :));
%similarity = lcs / (max(m,n));
similarity = lcs / ((m+n)/2);
dis=1-similarity;

   
