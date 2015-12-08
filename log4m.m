classdef log4m < handle 
    % It is modified version of LOG4M (http://www.mathworks.com/matlabcentral/fileexchange/37701-log4m-a-powerful-and-simple-logger-for-matlab)
    % which description is as follows:
    % 
    % LOG4M This is a simple logger based on the idea of the popular log4j.
    %
    % Description: Log4m is designed to be relatively fast and very easy to
    % use. It has been designed to work well in a matlab environment.
    % Please contact me (info below) with any questions or suggestions!
    % 
    %
    % Author: 
    %       Luke Winslow <lawinslow@gmail.com>
    % Heavily modified version of 'log4matlab' which can be found here:
    %       http://www.mathworks.com/matlabcentral/fileexchange/33532-log4matlab
    %
    %
    % Copyright (c) 2012, Luke Winslow
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are
    % met:
    % 
    %     * Redistributions of source code must retain the above copyright
    %       notice, this list of conditions and the following disclaimer.
    %     * Redistributions in binary form must reproduce the above copyright
    %       notice, this list of conditions and the following disclaimer in
    %       the documentation and/or other materials provided with the distribution
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.
    %

    properties (Constant)
        ALL = 0;
        TRACE = 1;
        DEBUG = 2;
        INFO = 3;
        WARN = 4;
        ERROR = 5;
        FATAL = 6;
        OFF = 7;
    end
        
    properties(Access = protected)
        logger;
        lFile;
        fid;
        isActive;
        T
    end
    
    properties(SetAccess = protected)
        fullpath = 'log4m.log';  % Default file
        commandWindowLevel = log4m.ALL;
        logLevel = log4m.INFO;
    end
    
    methods (Static)
        function obj = getLogger( logPath, ~ )        
            persistent localObj;
            if(nargin == 0 || nargin == 2) 
                if isempty(localObj) || ~isvalid(localObj)
                    if (nargin == 0)
                        % warning('Logger is not initialized');
                    end
                    localObj = log4m();
                    localObj.isActive = false;            
                end
            else
                localObj = log4m(logPath);
                localObj.isActive = true;
            end            
            obj = localObj;                
        end
        
        function obj = deactivateLogger()
            obj = log4m.getLogger([], []);
        end
    end
    
    
%% Public Methods Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods                
        
        function setCommandWindowLevel(self,loggerIdentifier)
            self.commandWindowLevel = loggerIdentifier;
        end

        function setLogLevel(self,logLevel)
            self.logLevel = logLevel;
        end
        
        function result = isDebugOn(self)
            result = (self.commandWindowLevel <= log4m.DEBUG) || (self.logLevel <= log4m.DEBUG);
        end
        
        function setActive(self, isActive)
            self.isActive = isActive;
        end

        function setTime(self, T)
            self.T = T;
        end
        
        function close(self)
        	fclose(self.fid);
        end
        

%% The public Logging methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function trace(self, funcName, message)
            %TRACE Log a message with the TRACE level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            % 
            self.writeLog(self.TRACE,funcName,message);
        end
        
        function debug(self, funcName, message)
            %TRACE Log a message with the DEBUG level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            % 
            self.writeLog(self.DEBUG,funcName,message);
        end
        
 
        function info(self, funcName, message)
            %TRACE Log a message with the INFO level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            % 
            self.writeLog(self.INFO,funcName,message);
        end
        

        function warn(self, funcName, message)
            %TRACE Log a message with the WARN level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            % 
            self.writeLog(self.WARN,funcName,message);
        end
        

        function error(self, funcName, message)
            %TRACE Log a message with the ERROR level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            % 
            self.writeLog(self.ERROR,funcName,message);
        end
        

        function fatal(self, funcName, message)
            %TRACE Log a message with the FATAL level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            % 
            self.writeLog(self.FATAL,funcName,message);
        end
        
        
    end

%% Private Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Unless you're modifying this, these should be of little concern to you.
    methods (Access = private)
        
        function self = log4m(path)
            if (nargin == 1)
                self.setFilename(path);
            end
        end

        function setFilename(self,logPath)
            %SETFILENAME Change the location of the text log file.
            %
            %   PARAMETERS:
            %       logPath - Name or full path of desired logfile
            %
            
            [self.fid,message] = fopen(logPath, 'a');
            
            if(self.fid < 0)
                error(['Problem with supplied logfile path: ' message]);
            end
            
            self.fullpath = logPath;
        end

        
%% WriteToFile        
        function writeLog(self,level,scriptName,message)
            
            if(self.isActive == false)
                return;
            end
            
            message = sprintf('T[%d] %s', self.T, message);
            
            % If necessary write to command window
            if( self.commandWindowLevel <= level )
                fprintf('%s: %s\n', scriptName, message);
            end
            
            %If currently set log level is too high, just skip this log
            if(self.logLevel > level)
                return;
            end 
            
            % set up our level string
            switch level
                case{self.TRACE}
                    levelStr = 'TRACE';
                case{self.DEBUG}
                    levelStr = 'DEBUG';
                case{self.INFO}
                    levelStr = 'INFO';
                case{self.WARN}
                    levelStr = 'WARN';
                case{self.ERROR}
                    levelStr = 'ERROR';
                case{self.FATAL}
                    levelStr = 'FATAL';
                otherwise
                    levelStr = 'UNKNOWN';
            end

            % Append new log to log file
            try
                fprintf(self.fid,'%s %s %s - %s\r\n' ...
                    , datestr(now,'yyyy-mm-dd HH:MM:SS') ...
                    , levelStr ...
                    , scriptName ... % Have left this one with the '.' if it is passed
                    , message);
            catch ME_1
                display(ME_1);
            end
        end
    end
    
end

