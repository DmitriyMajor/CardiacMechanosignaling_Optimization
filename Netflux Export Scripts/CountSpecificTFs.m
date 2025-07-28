function TFCounts = CountSpecificTFs(filename, geneIdx, speciesNames)
    % CountSpecificTFs_GenesOnly: Count TF occurrences in reactions targeting genes
    %
    % Inputs:
    %   filename      - Excel file name (e.g., 'RxnParsed_CaoSMechanosignaling.xlsx')
    %   geneIdx       - Row indices (relative to row 3) of species that are genes
    %   speciesNames  - Cell array of species names in Excel order (row 3 onward)
    %
    % Optional (commented): validatedIdx - indices of validated genes (subset of geneIdx)

    % Define upstream TFs of interest
    tfList = {'NFAT', 'CREB', 'MEF2', 'cJun', 'cMyc', 'cFos', ...
              'FOXO', 'SRF', 'NFkB', 'STAT', 'GATA4'};

    % Initialize TF count map
    TFCounts = containers.Map(tfList, zeros(1, numel(tfList)));

    % Read entire reactions sheet
    raw = readcell(filename, 'Sheet', 'reactions');

    % Extract reaction rules (starting from row 3, column 3)
    rules_raw = raw(3:end, 3);  % column C, row 3 onward
    rules_raw = rules_raw(~cellfun(@(x) isempty(x) || (isstring(x) && strlength(x) == 0), rules_raw));

    % Convert geneIdx to row indices in spreadsheet (offset by 2 since row 3 is index 1)
    geneNames = speciesNames(geneIdx);


    % (Optional) To use validatedIdx instead, comment above and uncomment this:
    validatedIdx = 85:654;
    genedx = validatedIdx;
    geneNames = speciesNames(validatedIdx);

    % Parse reactions
    for i = 1:length(rules_raw)
        rule = strtrim(rules_raw{i});
        if contains(rule, '=>')
            parts = strsplit(rule, '=>');
            lhs = strtrim(parts{1});
            rhs = strtrim(parts{2});

            if isempty(lhs) || isempty(rhs)
                continue
            end

            % Check that RHS is a gene
            if ~any(strcmp(rhs, geneNames))
                continue
            end

            % Parse LHS regulators
            lhs_clean = strrep(lhs, '!', '');
            regulators = strsplit(lhs_clean, '&');
            regulators = strtrim(regulators);

            % Count TFs if in list
            for r = 1:length(regulators)
                reg = regulators{r};
                for t = 1:length(tfList)
                    if strcmpi(reg, tfList{t})
                        TFCounts(tfList{t}) = TFCounts(tfList{t}) + 1;
                    end
                end
            end
        end
    end

    %% Plot bar chart
    tfVals = cell2mat(values(TFCounts, tfList));

    figure;
    b = bar(tfVals, 'FaceColor', 'flat');
    blueColor = [0.2 0.6 0.8];
    greenColor   = [155, 243, 193] / 255;
barColors = [
    blueColor;
    blueColor;
    blueColor;
    blueColor;
    blueColor;
    blueColor; blueColor; blueColor; greenColor; greenColor; blueColor
];
b.CData = barColors;  
    set(gca, 'XTickLabel', tfList, 'XTick', 1:numel(tfList), ...
        'XTickLabelRotation', 45, 'FontSize', 12)
    ylabel('Number of Gene-Targeting Reactions')

    if geneIdx(end) == 654
        title('Upstream TF Frequency in Gene-Targeting Reactions (Validated Genes)')
        ylim([0 170])
    else
        title('Upstream TF Frequency in Gene-Targeting Reactions (All Genes)')
        ylim([0 350])
    end

    grid on

    % Show count values on top of bars
    xtips = b.XEndPoints;
    ytips = b.YEndPoints;
    labels = string(tfVals);
    text(xtips, ytips, labels, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 11)
end
