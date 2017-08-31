function cac = ReadDefectFile(file)
    %%
        str = fileread(file);
        str = regexprep( str, '(?<=\r?\n)[ ]*(?=\r?\n)' ...
                  , 'nan nan', 'emptymatch' );
        cac = textscan( str, '%f%f', 'CollectOutput', true );
    end 