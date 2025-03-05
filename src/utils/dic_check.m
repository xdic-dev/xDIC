% parallel checking
if parallel_processing
    disp('Parallel processing is enabled');
    delete(gcp('nocreate'));
end