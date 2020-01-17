function comment_blocks(simulink_name, block_list, active)
%comment_blocks(simulink_name, block_name, active);
% inputs: simulink_name    -> simulink file name
%         block_name       -> name of the simulink block
%         active           -> active block or comment
%
%comment_blocks(simulink_name, block_list);
% inputs: simulink_name    -> simulink file name
%         block_list       -> struct array of blocks
%                               block.name
%                               block.active

    simulink_name = split(simulink_name, '.');
    simulink_name = simulink_name{1};
    
    blocks = find_system(simulink_name, 'IncludeCommented', 'on');

    is_single_block = (exist('active','var') && ~isempty(active));
    
    if is_single_block
        comment_block(simulink_name, blocks, block_list, active);
    else
        for block=block_list
            comment_block(simulink_name, blocks, block.name, block.active);
        end
    end
    
end



function comment_block(simulink_name, blocks, block_name, active)
    regex_pattern = strcat(block_name, '$');
    matches = regexp(blocks, regex_pattern);
    index = find(~cellfun(@isempty,matches));

    if isempty(index)
        message = strcat('Não foi possível encontrar o bloco ''', ...
            block_name, ''' na simulação ''', simulink_name, '''.');
        warning(message)
        return
    end
    if length(index)~=1
        message = strcat('Mais de um bloco com o nome ''', ...
            block_name, ''' na simulação ''', simulink_name, '''.');

        conflict_source = sprintf('\\n\\t%s', blocks{index});
        message = strcat(message, conflict_source);

        disp(message)
        ME = MException('CommentException:multipleMatch', ...
            message);
        throw(ME)
    end

    block_id = blocks{index}; 

    if active == 1
        set_param(block_id,'Commented','off');
    else
        set_param(block_id,'Commented','on');
    end
end