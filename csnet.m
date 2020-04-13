function csn = csnet(data,c,alpha,boxsize,weighted,textdata)

%Construction of cell-specific networks
%The function performs the transformation from gene expression matrix to
%cell-specific network (csn).
%data: Gene expression matrix, rows = genes, columns = cells
%c: Construct the CSNs for all cells, set c = [] (Default);
%   Construct the CSN for cell k, set  c = k
%alpha: Significant level (eg. 0.001, 0.01, 0.05 ...)
%       larger alpha leads to more edges, Default = 0.01
%boxsize: Size of neighborhood, Default = 0.1
%weighted: 1  edge is weighted
%          0  edge is not weighted (Default)
%csn: Cell-specific network, the kth CSN is in csn{k}
%     rows = genes, columns = genes

% Too many cells or genes may lead to out of memory.

if nargin < 5 || isempty(weighted)
    weighted = 0;
end
if nargin < 4 || isempty(boxsize)
    boxsize = 0.1;
end
if nargin <3 || isempty(alpha)
    alpha = 0.01;
end

% devuelve el nÃºmero de filas y columnas cuando A es una matriz.
[n1,n2] = size(data);

%disp(n1)
%disp(n2)

if nargin <2 || isempty(c)
    c = 1 : n2;
end
 
%Define the neighborhood of each plot
upper = zeros(n1,n2);
lower = zeros(n1,n2);
for i = 1 : n1
    [s1,s2] = sort(data(i,:));
    n3 = n2-sum(sign(s1));
    h = round(boxsize/2*sum(sign(s1)));
    k = 1;
    while k <= n2
        s = 0;
        while k+s+1 <= n2 && s1(k+s+1) == s1(k)
            s = s+1;
        end
        if s >= h
            upper(i,s2(k:k+s)) = data(i,s2(k));
            lower(i,s2(k:k+s)) = data(i,s2(k));
        else
            upper(i,s2(k:k+s)) = data(i,s2(min(n2,k+s+h)));
            lower(i,s2(k:k+s)) = data(i,s2(max(n3*(n3>h)+1,k-h)));
        end
        k = k+s+1;
    end
end
 
%Construction of cell-specific network

csn = cell(1,n2);
%  Por ejemplo, zeros(2,3) devuelve una matriz de 2 por 3.
B = zeros(n1,n2);
p = -icdf('norm',alpha,0,1);
for k = c
    for j = 1 : n2
        B(:,j) = data(:,j) <= upper(:,k) & data(:,j) >= lower(:,k);
    end
    a = sum(B,2);
    d = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
    d(1 : n1+1 : end) = 0;
    if weighted
        csn{k} = d.*(d > 0);
        
    else
        csn{k} = sparse(d > p);
        
    end
    
    disp(['Cell ' num2str(k) ' is completed']);
    S = '----Hello World----';
    %disp(textdata)
    nombre_sample = textdata(1,:)

end

dim = n1-1; 
disp(dim)

for l=1: size(csn,1) %scan rows    
    for m=1:size(csn,2) %scan columns
        hij=csn{l,m}(:,:); %convert matrix form
        disp(l)
        disp(m)
        
        formatSpec = "%s_%d";
        A3 = nombre_sample{m+1}; 
        A1 = m;
        str = sprintf(formatSpec,A3,A1);
        %Xfilename = sprintf('Cell_%d.txt',m); 
        % que quede como nombre de archivo el nombre de la muestra
        % generar archivo con las descripciones de los genes
        disp(str)
        writematrix(hij,str,'Delimiter','tab')
        fileID = fopen('Gene_names.txt','w');
        fprintf(fileID,'%s;\n',textdata{:,1});
        
    end 
end

fclose(fileID);
