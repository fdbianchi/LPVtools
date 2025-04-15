
clc
clearvars

constraint(1) = synConst.Gain({'w1','w2'},'z1','bound',2);
constraint(2) = synConst.Gain({'w3'},'z2','factor',0.5);
constraint(3) = synConst.Gain({'w4','w5'},{'z1','z2'},'factor',0.5);
constraint(4) = synConst.Gain({'w3','t'},'z2','bound',3);
constraint(6) = synConst.Gain(1,[2 3],'bound',3);
constraint(7) = synConst.Poles('MaxFreq',1000);


n = length(constraint); jj = 1;
for ii = 1:n
    if ~isempty(constraint(ii))
        [s1,s2,s3] = char(constraint(ii));
        str1{jj} = s1; str2{jj} = s2; str3{jj} = s3;
        jj = jj + 1;
    end
end
idx_min = find(ismember(str1,'min'));
idx_st  = find(ismember(str1,'s.t.'));

fprintf('Solved the following synthesis problem:\n\n')
if any(idx_min)
    fprintf('\tminimize:\n')
    for ii = idx_min(1:end-1)
        fprintf('           %s + \t\t\t(%s)\n',str2{ii},str3{ii})
    end
    ii = idx_min(end);
    fprintf('           %s  \t\t\t(%s)\n',str2{ii},str3{ii})
end
fprintf('\n')
if any(idx_st)
    fprintf('\tsubject to:\n')
    for ii = idx_st
        if strcmp(str3{ii},'')
        fprintf('           %s\n',str2{ii})
        else
        fprintf('           %s \t(%s)\n',str2{ii},str3{ii})
        end
    end
end
fprintf('\n')


