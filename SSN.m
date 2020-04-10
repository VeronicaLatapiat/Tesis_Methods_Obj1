%%%%%%%% SSN method %%%%%%%%

function [ index_R,p ] = SSN( sample,ref )
%function:construct the SSN
%   Input:
%         sample:calculated sample
%          ref:the reference samples
%   Output:
%         adjacency_matrix:the network structure
%a example
% sample=new_T(:,1);
% ref=new_N;

% [R,P] = corrcoef(___) devuelve la matriz de coeficientes de correlación 
% y la matriz de valores p para probar la hipótesis de que no existe relación
% entre los fenómenos observados (hipótesis nula). Si un elemento 
% OFF-diagonal de P es menor que el nivel de significancia (el valor 
% predeterminado es 0.05), la correlación correspondiente en R se considera
% significativa. Esta sintaxis no es válida si R contiene elementos complejos.

[R,P]=corrcoef(ref');
final_R0=R;
final_R0(isnan(final_R0))=0;

%disp("final_R0_0")
%disp(final_R0)

%NEW_data=[sample ref];
NEW_data=[ref sample];
[R1,P1]=corrcoef(NEW_data');
final_R1=R1;
final_R1(isnan(final_R1))=0;

index_R=final_R1-final_R0;

%disp("final_R1")
%disp(final_R1)
%disp("final_R0")
%disp(final_R0)

%disp(index_R)

[m,n]=size(ref);
Z=index_R./((1-final_R0.^2)/(n-1));
Z(Z==inf)=max(max(Z));
Z(Z==-inf)=-max(max(Z));
Z(isnan(Z))=0;

clear NEW_data final_R1 final_R0 R0 R1 P P1
p=1-normcdf(abs(Z));

end

