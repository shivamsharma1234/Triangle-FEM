function [ a ] = Assemble_shivam_70503(k, connectivity)
    elements = size(connectivity, 1);
   
    
    kRows = size(k, 1);
    kCols = size(k, 2);
    
    if (kCols > 1)
        % Assemble K
        a = zeros(max(max(connectivity)));

        for e=1:elements
            for i=1:kRows
                for j=1:kCols
                    rowIndex = connectivity(e, i);
                    colIndex = connectivity(e, j);

                    a(rowIndex, colIndex) = a(rowIndex, colIndex) + k(i, j);
                end
            end
        end
    else
        % Assemble F
        a = zeros(max(max(connectivity)), 1);
        
        for e=1:elements
            for i=1:kRows
                a(connectivity(e, i)) = a(connectivity(e, i)) + k(i);
            end
        end
    end
end

