function imcell = Myarray2cell(imarray)
   N = size(imarray,3); %number of image
   imcell = cell(1,N); 
   for i = 1:N
       imcell{1,i} = imarray(:,:,i);
   end
end