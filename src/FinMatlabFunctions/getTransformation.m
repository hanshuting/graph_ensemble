function trans_list = getTransformation(indicator_name)
%original: 1
%difference: 2
%percent change: 3

[transInfo,tickerList] = xlsread('transformationList.xlsx');

loc = find(strcmpi(tickerList,indicator_name));

trans_list_raw = transInfo(loc,:);

trans_list = cell(0);

trans_list_raw = trans_list_raw(~isnan(trans_list_raw));

for k =1:length(trans_list_raw)
    currentCode = trans_list_raw(k);
    
    switch currentCode
        case 1
            trans_list = [trans_list;'original'];
        case 2
            trans_list = [trans_list;'difference'];
        case 3
            trans_list = [trans_list;'percent change'];
        %case NaN
        %    break;
        otherwise
            display('No such transformation');
            display(currentCode);
            return
    end
    
end
