fid = fopen('result_parallel.dat','r');
datacell = textscan(fid, '%f');
fclose(fid);

result_parallel = datacell{1};


fid = fopen('result_single.dat.','r');
datacell = textscan(fid, '%f');
fclose(fid);

result_single = datacell{1};
diff = result_single-result_parallel;

assert(max(diff) == 0 && min(diff) == 0 )

