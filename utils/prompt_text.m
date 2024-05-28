function response = prompt_text(question, default_response, confirm_response)
    arguments
        question {mustBeText}
        default_response {mustBeText} = ''
        confirm_response logical = true
    end

    % Set the default response
    response = default_response;

    % Format the question
    if isempty(default_response)
        question_full = question;
    else
        question_full = sprintf('%s (default = ''%s'')', question, default_response);
    end
    question_full = append(question_full, '\n>> ');

    % Prompt the user
    response = input(question_full, 's');

    % Parse the response
    if isempty(response)
        % Return the default response if one was provided and response
        % is empty
        if ~isempty(default_response)
            response = default_response;
            return
        % If response is empty, try again
        else
            fprintf('A response is required\n');
            response = prompt_text(question, default_response, confirm_response);
            return
        end
    end

    % Get confirmation
    if confirm_response
        xx = false;
        require_response = true;
        confirmed = prompt_bool( ...
            sprintf('Is ''%s'' correct?', response), xx, require_response ...
        );
        if ~confirmed
            response = prompt_text(question, default_response, confirm_response);
        end
    end
end
