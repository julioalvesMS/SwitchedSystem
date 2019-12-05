function uncomment_blocks(simulink_name)
%uncomment_blocks(simulink_name);
% inputs: simulink_name    -> simulink file name

    simulink_name = split(simulink_name, '.');
    simulink_name = simulink_name{1};
    
    blocks = find_system(simulink_name, 'IncludeCommented', 'on');
    
    for id=2:length(blocks)
        set_param(blocks{id},'Commented','off');
    end
    
end