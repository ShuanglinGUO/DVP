function options = setOptions(options, options_default)
%Set the parameters of options using the input options and options_default
% Author: Shuanglin GUO @PIMM of ENSAM
% February 22, 2021
% Email: shuanglin.guo@ensam.eu

if ~isstruct(options)
    error('''options'' should be struct class!');
end
if ~isstruct(options_default)
    error('''options_default'' should be struct class!');
end

fieldlist = fieldnames(options);

for i = 1:length(fieldlist)
    if isfield(options_default, fieldlist{i})
        if isempty( options.(fieldlist{i}) )
            error(['''options.',fieldlist{i},''' should not be empty! Please provide a valid value.']);
        end
        
        if isa( options.(fieldlist{i}), class(options_default.(fieldlist{i})) )
            options_default.(fieldlist{i}) = options.(fieldlist{i});
        else
            error(['The class of ''options.',fieldlist{i},''' should be consistent with the class of ''options_default.',fieldlist{i},'''']);
        end
        
    else
        error(['There is no field ''',fieldlist{i},''' in the options_default.''']);
    end
end

options = options_default;

end