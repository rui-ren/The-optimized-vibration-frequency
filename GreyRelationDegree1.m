
function GreyRelationDegree1(stats)    % stats 是一个m×n 的评价矩阵， m 个评价对象、n 个评价指标
% 设置参考数列，即各指标的理想最优组成的行向量，维数与stats相同。
% 原始评价矩阵及样本序号
[r,c]=size(stats);            % stats的行数和列数，即评价对象的个数及评价指标的个数
samNo=1:r;                    % 样本数量
% 数据规范化处理，将各指标数据与参考数列一起规范化到0-1之间
stdMatrix=zeros(r,c);         % 给标准化矩阵分配空间，第一行为参考数列的标准化值，第二行至最后一行为原始评价矩阵的标准化值
maxOfCols=max(stats);         % 包括参数列在内的各列的最大值
minOfCols=min(stats);         % 包括参数列在内的各列的最小值

for j=1:c
    for i=1:r
        stdMatrix(i,j)=(stats(i,j)-minOfCols(j))./(maxOfCols(j)-minOfCols(j)); % 计算标准化
    end
end

R_0=[1 1 1 0 1 0];

% 计算关联系数
absValue=zeros(r,c);                        % 给绝对差值序列分配空间
for i=1:r
    absValue(i,:)=abs(stdMatrix(i,:)-R_0);  % 绝对差值序列计算
end
minAbsValueOfCols=min(absValue,[],1);    % absValue每一列的最小值
maxAbsValueOfCols=max(absValue,[],1);    % absValue每一列的最大值
minAbsValue=min(minAbsValueOfCols);      % absValue的最小值
maxAbsValue=max(maxAbsValueOfCols);      % absValue的最大值
defCoeff=0.5;                            % 设置分辨系数为0.5
relCoeff=(minAbsValue+defCoeff*maxAbsValue)./(absValue+defCoeff*maxAbsValue);  % 关联系数计算
 
% 计算关联度
% 构造A-m
[V,D]=eigs(P);
V=V(:,1);
w=V(:,1)/sum(V(:,1));
R=zeros(r,1);
for i=1:r
        R(i,1)=relCoeff(i,:)*w;  % 关联度计算
end
 
% 权重可视化
[sortW,IXW]=sort(w,'descend');   % 权重降序排序，IXW确保对应的指标名称一致
indexes={};
for i=1:c                        % c为总共的指标参数
    indexes(i)={strcat('指标',num2str(i))}; % 指标名称为“指标1”、指标“2”……
end
sortIndex=indexes(IXW);                % 排序后与权重对应的指标名称
figure;
subplot(1,2,1);
bar(w);
xlim([0 c+1]);   % 设置x轴范围
xlabel('指标名称','FontSize',12,'FontWeight','bold');
set(gca,'xtick',1:c);
set(gca,'XTickLabel',indexes,'FontWeight','light');
ylabel('权重','FontSize',12,'FontWeight','bold');
set(gca,'YGrid','on');
for i=1:c
    text(i-0.35,w(i)+0.005,sprintf('%.3f',w(i)));
end
title('指标权重可视化');
box off;
 
subplot(1,2,2);
bar(sortW);
xlim([0 c+1]);   % 设置x轴范围
xlabel('指标名称','FontSize',12,'FontWeight','bold');
set(gca,'xtick',1:c);
set(gca,'XTickLabel',sortIndex,'FontWeight','light');
ylabel('权重','FontSize',12,'FontWeight','bold');
set(gca,'YGrid','on');
for i=1:c
    text(i-0.35,sortW(i)+0.005,sprintf('%.3f',sortW(i)));
end
title('指标权重可视化（降序排列）');
box off;
 
 
% 关联度分析结果展示
[sortR,IX]=sort(R,'descend');  % 关联度降序排序，IX确保对应的样本序号一致
sortSamNo=samNo(IX);           % 排序后与关联度对应的样本序号
figure;
subplot(2,1,1);
plot(R,'--ro',...
    'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'MarkerSize',10);
xlim([1 r]);   % 设置x轴范围
xlabel('17 1/2 钻头','FontSize',10,'FontWeight','bold');
set(gca,'xtick',1:r);
set(gca,'XTickLabel',samNo,'FontWeight','light');
ylabel('关联度','FontSize',10,'FontWeight','bold');
title('钻头灰色关联度综合评价结果');
grid on;

 
subplot(2,1,2);
plot(sortR,'--ro',...
    'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'MarkerSize',10);
xlim([1 r]);   % 设置x轴范围
xlabel('17 1/2 钻头','FontSize',10,'FontWeight','bold');
set(gca,'xtick',1:r);
set(gca,'XTickLabel',sortSamNo,'FontWeight','light');
ylabel('关联度','FontSize',10,'FontWeight','bold');
title('钻头灰色关联度综合评价结果');
grid on;
hold off;
 
end


   