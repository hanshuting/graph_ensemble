function param = getMethodParam(method, freq)
%INPUT:
%method is char
%freq is char

%OUTPUT:
%parameter for method

%The incorporated available parameter decision rules are as of 9 Mar 2015.

param = [];

switch method
    case 'EFI'
        param = getEFIParam(freq);
    case 'EWI'
        param = getEWIParam(freq);
    otherwise
       display('This method currently not available!');
       return
end

end


function paramEFI = getEFIParam(freq) 

paramEFI = {};
erMess = 'This frequency not assigned parameters yet!';
switch freq
    case 'daily'
        paramEFI = [0.5:0.5:1];
    case 'monthly'
        paramEFI = [0.5:0.5:1];
    case 'quarterly'
        paramEFI = [0.5:0.5:1];
    otherwise 
        display(erMess);
        return
end
end


function paramEWI = getEWIParam(freq)
erMess = 'This frequency not assigned parameters yet!';

switch freq
    case 'daily'
        paramEWI = 10;
    case 'monthly'
        paramEWI = 4;
    case 'quarterly'
        paramEWI = 2;
    otherwise 
        display(erMess);
        return
end

end