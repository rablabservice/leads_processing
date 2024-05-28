function response = prompt_bool( ...
    yes_no_question, ...
    default_response, ...
    require_response, ...
    ask_twice ...
)
    % Ask the user a yes/no question and return the response
    %
    % Parameters
    % ----------
    % yes_no_question : char or str
    %     The question to ask the user
    % default_response : logical, optional
    %     The default response if the user does not provide one and
    %     require_response is false
    %     - If default_response is true, returns true unless the user
    %       responds 'n'
    %     - If default_response is false, returns false unless the user
    %       responds 'y'
    % require_response : logical, optional
    %     If require_response is true, the user must provide a response
    %     of 'y' or 'n'
    % ask_twice : logical, optional
    %     If ask_twice is true, the user will be prompted to confirm
    %     their response a second time
    % ------------------------------------------------------------------
    arguments
        yes_no_question {mustBeText}
        default_response logical = false
        require_response logical = false
        ask_twice logical = false
    end

    % Some variables we may use
    xx = false;
    dbl_check = 'Are you sure?';

    % Enter the prompt loop
    if require_response
        % Format the question
        yes_no_question_full = append(yes_no_question, ' (y/n)\n>> ');

        % Prompt the user
        user_response = input(yes_no_question_full, 's');

        % If the user didn't respond, try again
        if isempty(user_response)
            fprintf('A response is required\n');
            response = prompt_bool(yes_no_question, xx, require_response, ask_twice);
            return
        end

        % Take the lowercase first letter of the response
        user_response = lower(user_response(1));

        % Parse the response
        if strcmp(user_response, 'y')
            if ask_twice
                ask_twice = false;
                confirmed = prompt_bool(dbl_check, xx, require_response, ask_twice);
                if confirmed
                    response = true;
                else
                    ask_twice = true;
                    response = prompt_bool( ...
                        yes_no_question, xx, require_response, ask_twice ...
                    );
                end
            else
                response = true;
            end
        elseif strcmp(user_response, 'n')
            if ask_twice
                ask_twice = false;
                confirmed = prompt_bool(dbl_check, xx, require_response, ask_twice);
                if confirmed
                    response = false;
                else
                    ask_twice = true;
                    response = prompt_bool( ...
                        yes_no_question, xx, require_response, ask_twice ...
                    );
                end
            else
                response = false;
            end
        else
            fprintf('Response must be ''y'' or ''n''\n');
            response = prompt_bool(yes_no_question, xx, require_response, ask_twice);
        end
    else
        % Set the default response
        response = default_response;

        % Format the question
        if default_response
            yes_no_question_full = append(yes_no_question, ' (default=y)');
        else
            yes_no_question_full = append(yes_no_question, ' (default=n)');
        end
        yes_no_question_full = append(yes_no_question_full, ' (y/n)\n>> ');

        % Prompt the user
        user_response = input(yes_no_question_full, 's');

        % If the user didn't respond, return the default response
        if isempty(user_response)
            return
        end

        % Take the lowercase first letter of the response
        user_response = lower(user_response(1));

        % Parse the response
        if default_response
            if strcmp(user_response, 'y')
                return
            elseif strcmp(user_response, 'n')
                if ask_twice
                    require_response = true;
                    ask_twice = false;
                    response = ~prompt_bool(dbl_check, xx, require_response, ask_twice);
                else
                    response = false;
                end
            else
                fprintf('Response must be ''y'' or ''n''\n');
                response = prompt_bool( ...
                    yes_no_question, default_response, require_response, ask_twice ...
                );
            end
        else
            if strcmp(user_response, 'n')
                return
            elseif strcmp(user_response, 'y')
                if ask_twice
                    require_response = true;
                    ask_twice = false;
                    response = prompt_bool(dbl_check, xx, require_response, ask_twice);
                else
                    response = true;
                end
            else
                fprintf('Response must be ''y'' or ''n''\n');
                response = prompt_bool( ...
                    yes_no_question, default_response, require_response, ask_twice ...
                );
            end
        end
    end
end
