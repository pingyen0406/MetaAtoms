inf = 'D:/Dropbox/MATLAB/MetaAtoms/test/focus_180.txt';
R_list = readmatrix(inf);
R_list = reshape(R_list,36,108);
surface(R_list,'CDataMapping','scaled');
view(3);
