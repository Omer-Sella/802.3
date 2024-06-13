function [config_file,num_fext,num_next,Remember_keyword,OP,varargin]=COM_CommandLine_Parse(OP,varargin)


keywords={'Legacy' 'TD' 'Config2Mat'};
Remember_keyword='Legacy';
OP.TDMODE=false;
OP.GET_FD=true;
OP.CONFIG2MAT_ONLY=false;
config_file='';
num_fext=[];
num_next=[];
if ~isempty(varargin)
    if ~ischar(varargin{1})
        error('First input must be a string');
    end
    keyword_idx=find(strcmpi(keywords,varargin{1}));
    if isempty(keyword_idx)
        %Legacy Mode
        [config_file,varargin]=varargin_extractor(varargin{:});
        [num_fext,varargin]=varargin_extractor(varargin{:});
        [num_next,varargin]=varargin_extractor(varargin{:});
    else
        %Keyword Mode
        my_keyword=varargin{1};
        Remember_keyword=my_keyword;
        varargin(1)=[];
        switch my_keyword
            
            case 'Legacy'
                [config_file,varargin]=varargin_extractor(varargin{:});
                [num_fext,varargin]=varargin_extractor(varargin{:});
                [num_next,varargin]=varargin_extractor(varargin{:});
            case 'TD'
                OP.TDMODE=true;
                OP.GET_FD=false;
                [config_file,varargin]=varargin_extractor(varargin{:});
                [num_fext,varargin]=varargin_extractor(varargin{:});
                [num_next,varargin]=varargin_extractor(varargin{:});
            case 'Config2Mat'
                OP.CONFIG2MAT_ONLY=true;
                [config_file,varargin]=varargin_extractor(varargin{:});
        end
    end
end