function validatedIdx = getValidatedGeneIdx(speciesNames, geneIdx)
% GETVALIDATEDGENEIDX  Return indices of validated gene species using readmatrix/readcell
%   validatedIdx = getValidatedGeneIdx(speciesNames, geneIdx)
%   Reads PubMed IDs via readmatrix and rules via readcell from hard-coded Excel.

    % Hard-coded file and sheet
    filename = 'RxnParsed_CaoSMechanosignaling.xlsx';
    sheet    = 'reactions';

    % Read numeric PMIDs for rows 126:695, column I
    pmids = readmatrix(filename, 'Sheet', sheet, 'Range', 'I126:I695');

    % Read rule strings for rows 126:695, column C
    rawRules = readcell(filename, 'Sheet', sheet, 'Range', 'C126:C695', 'TextType', 'string');
    % Convert to cell array of char
    if isstring(rawRules)
        rules = cellstr(rawRules);
    else
        rules = rawRules;
    end

    validatedIdx = [];
    for k = 1:numel(rules)
        % Proceed only if PMID is non-NaN and non-zero
        val = pmids(k);
        if ~isnan(val) && val~=0
            ruleStr = rules{k};
            if ischar(ruleStr) && contains(ruleStr, '=>')
                % Parse output species (right of '=>')
                parts   = strsplit(ruleStr, '=>');
                outName = strtrim(parts{2});

                % Match to speciesNames (with fallback)
                idx = find(strcmp(speciesNames, outName), 1);
                if isempty(idx)
                    idx = find(strcmp(speciesNames, ['gene_' outName]), 1);
                end
                % Collect only if it's in geneIdx
                if ~isempty(idx) && any(geneIdx == idx)
                    validatedIdx(end+1) = idx; %#ok<AGROW>
                end
            end
        end
    end

    % Return unique sorted indices
    validatedIdx = unique(validatedIdx);
end
