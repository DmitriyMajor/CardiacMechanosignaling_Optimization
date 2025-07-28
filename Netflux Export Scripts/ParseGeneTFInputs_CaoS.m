function [oneTF, twoTF, multiTF] = ParseGeneTFInputs_CaoS(geneIdx)
    % ParseTFCount: Parse the number of upstream TFs from Excel file
    %   oneTF    - indices (from geneIdx) of genes with 1 upstream TF
    %   twoTF    - indices (from geneIdx) of genes with 2 upstream TFs
    %   multiTF  - indices (from geneIdx) of genes with 3 or more upstream TFs

    filename = 'RxnParsed_CaoSMechanosignaling.xlsx';
    sheetname = 'species';

    % Read column 10 from row 3 to end (assuming max 1000 rows or can use detect)
    dataRange = 'J3:J1000';  % adjust max row as needed
    numData = readmatrix(filename, 'Sheet', sheetname, 'Range', dataRange);


    % numData will be a column vector of numbers or NaNs if empty

    % Filter by geneIdx (geneIdx indexes species starting at row 3, so relative to row 3)
    % geneIdx refers to species indices relative to the full species list, 
    % but data starts at row 3, so geneIdx should be indices starting at 1 for row 3

    % Extract TF counts for geneIdx (check indexing carefully)
    tfCounts = numData(geneIdx);

    % Find indices within geneIdx by TF counts
    oneTF = geneIdx(tfCounts == 1);
    twoTF = geneIdx(tfCounts == 2);
    multiTF = geneIdx(tfCounts >= 3);
end
