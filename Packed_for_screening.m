function Packed_for_screening(filename, initialFract, initialRatio, Parameters, iter)

SetParameters
%     save([filename '_parameters.mat'])

filename = [filename '_' num2str(iter,'%02d')];
BranchingColonyMultispecies_Core
SaveFigure
