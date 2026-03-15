function factor = fun_cal_usfactor(data)
%   Calculate undersampled factor of data
data = abs(data(:))>0;
factor = numel(data)/sum(data>0);
end