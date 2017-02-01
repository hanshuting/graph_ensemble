function trDatesMacroTS = toTrDateMacroTS(tradeDates,macroDates, macroVec)
%an auxilliary function
%tradeDates and macroDates are cells
%macroVec is a column matrix, each value corresponds to macroDates
%trDateMacroTS is a column matrix, each values corresponds to tradeDates

trDatesMacroTS = NaN(length(tradeDates),1);

for k = 1:length(macroDates)
    ind = find(strcmp(tradeDates,macroDates{k}));
    trDatesMacroTS(ind) = macroVec(k);
end


