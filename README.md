# Tesis_Methods_Obj1

cd /opt/MATLAB/R2019b/bin/glnxa64

./MATLAB

# csnet method

data= importdata('/Users/Vero_Latapiat/Documents/MATLAB/simulated_expr-data_TGFBeta_sc.txt')

csnet(data.data, [], 0.01, 0.1, 0, data.textdata)

