function print_data(fname,header,M)

% Description: prints array M to a file

[~,m] = size(M);

fid = fopen(fname,'wt');
fprintf(fid, '%s ', header);
fprintf(fid, '\n');
fmt = repmat('%15.15E ', 1, m);
fmt(end:end+1) = '\n';
fprintf(fid, fmt, M.'); 
fclose(fid);

end