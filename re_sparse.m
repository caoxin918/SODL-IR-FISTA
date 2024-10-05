function [vector] = re_sparse(vector, sourceNum, sparsity, threshold1, threshold2)
% Re-sparse stage, removing outliers based on the number of light sources,
% sparsity and threshold.
% input:
%   vector: The vector need to re-sparse.
%   sparsity: Number of non-zero elements in ans vector.
%   sourceNum: Number of light sources.
%   threshold1: For single source.
%   threshold2: For double(or more) source.
% output:
%   vector: Sparse vector.
%   outlier_indices: The coordinates of outliers.
% 
% 2024/04/17

numIdx = find(vector>0);
distances = zeros(size(numIdx));
[len, ~] = size(distances);
outlier_indices = [];
if len>sparsity
    for i=1:len
        for j=1:len
            distances(i) = distances(i) + abs(numIdx(i)-numIdx(j));    
        end
    end
    avg_distance = mean(distances);
    outlier_indices = numIdx(distances > threshold1 * avg_distance);    
      
    % 将离群点的值赋为零
    for i=1:length(outlier_indices)
        vector(outlier_indices(i)) = 0;
    end
    
    numIdx = find(vector>0);
    cnt = length(outlier_indices) + 1;
    if sourceNum>=2
        [groups, sortedIdx] = group_close_elements(numIdx, threshold2);
        count = cnt;
        for i=sourceNum+1:length(sortedIdx)
            temp = groups{sortedIdx(i)};
            for j=1:length(temp)
                outlier_indices(count) = temp(j);
                count = count + 1;
            end
        end
            
    end
    
    % 将非最大的两个组的元素值赋为零
    for i=cnt:length(outlier_indices)
        vector(outlier_indices(i)) = 0;
    end
end

end

function [groups, sortedIndices] = group_close_elements(vector, threshold)  
    % input:
    %   vector: Vector that needs to be grouped.
    %   threshold: Grouping criteria.
    % output:
    %   groups: Group based on threshold.
    %   soutedIndices: The number of elements in each group.(in descending order)
    %   

    % 初始化一个空的cell数组来存储组  
    groups = {};  
    % 当前组的起始索引  
    group_start = 1;  
      
    % 遍历向量  
    for i = 2:length(vector)  
        % 计算当前元素与前一个元素的差值  
        diff = abs(vector(i) - vector(i-1));  
          
        % 如果差值大于阈值，则当前元素不属于当前组  
        if diff > threshold  
            % 将当前组添加到groups中  
            groups{end+1} = vector(group_start:i-1);  
            % 更新当前组的起始索引  
            group_start = i;  
        end  
    end  
      
    % 添加最后一个组（如果存在）  
    if group_start <= length(vector)  
        groups{end+1} = vector(group_start:end);  
    end  

    if numel(groups)>=2

        for i=1:numel(groups)  
            lengths(i) = length(groups{i});  
        end  

        % 对向量长度进行降序排序，并获取排序后的索引
        [sortedLengths, sortedIndices] = sort(lengths, 'descend');
        % 提取前两个最大长度的向量的索引
        % maxIndices = sortedIndices(1:2);
    else
        sortedIndices = 1;
    end
end



