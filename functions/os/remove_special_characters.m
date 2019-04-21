function string = remove_special_characters(string)
    string = regexprep(string,'[^0-9a-zA-Z -]','');
    string = regexprep(string,'[ ]','_');
end

