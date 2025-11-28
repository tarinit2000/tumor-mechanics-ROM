function log_debug(level, msg, varargin)
% Simple leveled logger. Levels: DEBUG, INFO, WARN, ERROR
    persistent LOG_LEVEL
    if isempty(LOG_LEVEL)
        LOG_LEVEL = 'INFO'; % change to 'DEBUG' for verbose output
    end
    levels = {'DEBUG','INFO','WARN','ERROR'};
    lvl_idx = find(strcmp(levels, level));
    cur_idx = find(strcmp(levels, LOG_LEVEL));
    if isempty(lvl_idx) || isempty(cur_idx)
        fprintf('[%s] %s\n', datestr(now,'HH:MM:SS'), sprintf(msg,varargin{:}));
        return;
    end
    if lvl_idx >= cur_idx
        fprintf('[%s] %5s: %s\n', datestr(now,'HH:MM:SS'), level, sprintf(msg,varargin{:}));
    end
end
