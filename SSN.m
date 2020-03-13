#Functions operate on variables within their own workspace, which is also called the local workspace, separate from the workspace you access at the MATLAB command prompt which is called the base workspace.
#Functions can accept more than one input arguments and may return more than one output arguments.
#Syntax of a function statement is −
#function [out1,out2, ..., outN] = myfun(in1,in2,in3, ..., inN)

function [ index_R,p ] = SSN( sample,ref )
%function:construct the SSN
%   Input:
%         sample:calculated sample
%          ref:the reference samples
%   Output:
%         adjacency_matrix:the network structure
%a example
# first dimension is all "row" and the second is "column"
% sample=new_T(:,1); y a mi me interesan todas las columnas
# entonces deberia ser sample=new_T(:,:);
% ref=new_N;

#[R,P] = corrcoef(___) returns the matrix of correlation coefficients and the matrix of p-values 
#for testing the hypothesis that there is no relationship between the observed phenomena (null hypothesis).
#Use this syntax with any of the arguments from the previous syntaxes. If an off-diagonal element of P is 
#smaller than the significance level (default is 0.05), then the corresponding correlation in R is considered significant.
#This syntax is invalid if R contains complex elements.

[R,P]=corrcoef(ref');
final_R0=R; # matriz de correlacion
final_R0(isnan(final_R0))=0; # returns a logical array containing 1 (true) where the elements of A are NaN, and 0 (false) 
                             # where they are not. If A contains complex numbers, isnan(A) contains 1 for elements with 
                             # either real or imaginary part is NaN, and 0 for elements where both real and imaginary parts are not NaN.

NEW_data=[ref sample];
[R1,P1]=corrcoef(NEW_data');
final_R1=R1;
final_R1(isnan(final_R1))=0;

index_R=final_R1-final_R0;
[m,n]=size(ref);
Z=index_R./((1-final_R0.^2)/(n-1));
Z(Z==inf)=max(max(Z));
Z(Z==-inf)=-max(max(Z));
Z(isnan(Z))=0;

clear NEW_data final_R1 final_R0 R0 R1 P P1
p=1-normcdf(abs(Z)); # returns the cumulative distribution function (cdf) of the standard normal distribution, evaluated at the values in x.

end
