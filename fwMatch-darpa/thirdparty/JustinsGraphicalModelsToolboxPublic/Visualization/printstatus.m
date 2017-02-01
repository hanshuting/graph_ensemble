function printstatus(num,of,delete);
% num ranges from 0 to 1
% of is resolution printed

num = round(num*of);
if delete
    fprintf(repmat('\b',1,of+2));
end
fprintf([ '[' repmat('-',1,num) repmat(' ',1,of-num) ']']);
