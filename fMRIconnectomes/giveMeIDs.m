load('goc_rfmri_fc.mat')
name = {subs.name}.';
IDs = zeros(length(name),1); 

for i=1:length(name)
A = name{i};
IDs(i) = str2double(A(5:end));
end
IDs = IDs'; 