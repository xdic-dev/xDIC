function data_array = Mycell2array(data_cell)
   N = length(data_cell); 
   siz = size(data_cell{1}); 
   data_array = zeros(N,siz(1)); 
   for ii = 1:N
       data_array(ii,:) = data_cell{ii};
   end
end