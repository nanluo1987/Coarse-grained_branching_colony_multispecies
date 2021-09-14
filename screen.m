for P1 = 1.1 : 0.1 : 1.3
    
    for P2 = 0.3
        
        prefix = ['aC2=' num2str(P1) ',mu=' num2str(P2) '_'];
%         disp(prefix(1:end-1))
        aCs_act = [1, P1, 1] * 0.75;  
        mu = P2;
        BranchingColonyMultispecies_analyze
        fprintf([prefix(1:end-1) '  %f\n'],sum(C{2}(:))/(sum(C{2}(:))+sum(C{1}(:))))
        
    end
    
end


