function [data_out] = normalizeData(data_in)
data_n = size(data_in,1);
min_values = min(data_in);
max_values = max(data_in);
repmat_min_values = repmat(min_values,data_n,1);
repmat_max_values = repmat(max_values,data_n,1);
data_out = (data_in-repmat_min_values)./(repmat_max_values-repmat_min_values);
end