function [y0_MA, y30m_MA, y30m_MA_trans, y4h_MA, y4h_MA_trans, matchedGeneNames] = ImportRNA_CaoS(speciesNames, speciesTypes, y0)

    % Load RNA-seq data
    RNA = readtable('RNA_Data.xlsx');

    % Extract gene identifiers from RNA data
    RNA_symbols = RNA.Row_names;

    % Parse species that are genes
    geneIdx = find(contains(speciesTypes, 'Gene'));
    geneNames = speciesNames(geneIdx);

    % Initialize output vectors
    y0_MA = zeros(length(speciesNames), 1);
    y30m_MA = y0_MA;
    y30m_MA_trans = y0_MA;
    y4h_MA = y0_MA;
    y4h_MA_trans = y0_MA;

    matchedGeneNames = {};

    % Loop through all species
    for i = 1:length(speciesNames)
        if ismember(i, geneIdx)
            gene = extractAfter(speciesNames{i}, 'gene_');
            rowIdx = find(strcmp(RNA_symbols, gene));
            if ~isempty(rowIdx)
                cc = RNA.Control(rowIdx);
                y0_MA(i) = cc;
                % y30m_MA(i) = 2^RNA.log2FC_1(rowIdx) * cc;
                % y30m_MA_trans(i) = 2^RNA.log2FC(rowIdx) * cc;
                % y4h_MA(i) = 2^RNA.log2FC_3(rowIdx) * cc;
                % y4h_MA_trans(i) = 2^RNA.log2FC_2(rowIdx) * cc;
                y30m_MA(i) = 2^RNA.log2FC_1(rowIdx);
                y30m_MA_trans(i) = 2^RNA.log2FC(rowIdx);
                y4h_MA(i) = 2^RNA.log2FC_3(rowIdx);
                y4h_MA_trans(i) = 2^RNA.log2FC_2(rowIdx);

                matchedGeneNames{end+1} = speciesNames{i};
            else
                % RNA data missing — assign default value
                y0_MA(i) = 20;
                y30m_MA(i) = 20;
                y30m_MA_trans(i) = 20;
                y4h_MA(i) = 20;
                y4h_MA_trans(i) = 20;
            end
        else
            % Non-gene species — retain original initial value
            y0_MA(i) = y0(i);
            y30m_MA(i) = y0(i);
            y30m_MA_trans(i) = y0(i);
            y4h_MA(i) = y0(i);
            y4h_MA_trans(i) = y0(i);
        end
    end
end