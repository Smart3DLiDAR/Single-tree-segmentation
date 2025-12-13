function [finalroots,finaltrunk,finalpoints] = ExtrationTreeTrunk2(gp,cloth,gridsize,H_diff,cellsize2)
%%% 通过多联合约束（线性特征、法向量夹角特征、垂直度）将特征相似的点分割特征相似的点簇
%%% 计算点簇的高宽比特征，将线性特征、高宽比特征明显的树干点云提取出来，并找到根节点的位置
%%% input 
% gp                   原始点云
% gridsize             格网大小
%%% output  
% roots                树根节点（树的最低点）
% trunks               树桩点云
% roots_rowvector      树轴方向向量
% delete_toosmall      单木点云（未分割的点云）

%% 树干探测 提取根节点
% % 分块后再图分割 格网化截取部分地面以上点
xgrid=ceil(gp(:,1)/gridsize);
ygrid=ceil(gp(:,2)/gridsize);
XY_grid=[xgrid ygrid];
[C,ia,ic] = unique(XY_grid,'rows');
root0926=[];
trunk0926 = {};
%roots_gv=[];
%trees=[];
remov = [];

for o=1:size(C,1)
    xind=find(xgrid==C(o,1));
    yind=find(ygrid==C(o,2));
    xy_ind=intersect(xind,yind);
    grid_p=gp(xy_ind,:);

    figure(1)
    pcshow(grid_p(:,1:3));

    % 去除坡度影响  截取DTM2米以上的区域进行树干探测，找到
    z_cloud_interp = griddata(cloth(:, 1), cloth(:, 2), cloth(:, 3), grid_p(:, 1), grid_p(:, 2), 'linear');
    ob_data = grid_p((grid_p(:, 3) - z_cloud_interp) <= H_diff,:);
    if isempty(ob_data)
        continue
    end

    figure(2)
    pcshow(ob_data(:,1:3));
    
    above=ob_data(:,1:3);
    above(:,4)=0;
    tic
    ft_threshold = 0.8;paral = 0;plot = 1;
    [linearity,s1,Seg_final] = RecursiveSegmentation_release1(above(:,1:3),ft_threshold,paral,plot);
    for i=1:size(Seg_final,1)
        Seg_final{i,1}(:,4) = i;
        Seg_final{i,1}(:,5) = linearity(i,1);
        Seg_final{i,1}(:,6) = s1(i,1);
    end
    Seg_final2 = cell2mat(Seg_final);
    % 由于点比较稀疏点密度约束设置20(原30)      %%初始为50和0.9
    remove1 = find( Seg_final2(:,6)>30 & Seg_final2(:,5)>0.65 );
    %remove1 = find( Seg_final2(:,5)>0.65 );

    Seg_final2 = Seg_final2(remove1,:);
    ft_threshold = 0.8;%0.9
    ft_k = 40; j = 0;
    
    result2 = [];
    for i = 1:max(Seg_final2(:,4))
        wood_est = [];
        wood_est = Seg_final2(Seg_final2(:,4)==i,:);
        if ~isempty(wood_est)
            j = j+1;
            z_diff = max(wood_est(:,3))-min(wood_est(:,3));
            x_diff = max(wood_est(:,1))-min(wood_est(:,1));
            y_diff = max(wood_est(:,2))-min(wood_est(:,2));
            ratio(j) = z_diff/(sqrt(x_diff^2+y_diff^2));%高宽比
            wood_est(:,6) = j;
            wood_est(:,7) = ratio(j);
            result2 = [result2;wood_est];
        end
    end

    if isempty(result2)
        continue
    end

    if size(result2,1)<ft_k
        normals = pcnormals(pointCloud(result2(:,1:3)),size(result2,1)-1);
    else
        normals = pcnormals(pointCloud(result2(:,1:3)),ft_k);
    end
    dv = 1-abs(normals(:,3));%法向量夹角约束
    result_final = result2(result2(:,7)> 0.85 & dv>=ft_threshold,:);% 树桩部分
% %     tempremov = result2(result2(:,7)<1&dv<ft_threshold&Seg_final2(:,3)<min_h+1,:); 
% %     remov = [remov;tempremov];% 剔除的非树桩部分

    if isempty(result_final)
        continue
    end
    figure(3)
    pcshow(result_final(:,1:3),result_final(:,4),'MarkerSize',10);
    
    [~,ia,ib]=intersect(above(:,1:3),result_final(:,1:3),'rows');
    result_final(ib,8) = above(ia,4);
    roots = [];
    all_class=unique(result_final(:,6));
% %     temproots_gv = [];
    trunk1 = {};
    for i = 1:size(all_class,1)
        idx = [];
        roots_obj = [];
        roots_obj = result_final(result_final(:,6)==all_class(i),:);

        if ~isempty(roots_obj)
            [idx,~] = find(roots_obj(:,3)==min(roots_obj(:,3)));
            roots = [roots;roots_obj(idx(1),:)];
            temptrunk = roots_obj;
            trunk1 = [trunk1;temptrunk];
% %             if size(roots_obj,1)>=2
% %                 [coeff,~,~] = pca(roots_obj(:,1:3));
% %                 [~,maxid]= max(coeff(3,:));
% %                 tempgv = coeff(:,maxid)';
% %             else
% %                 tempgv = [0 0 1];
% %             end
% %             temproots_gv = [temproots_gv;tempgv];
% %             
% %         else
% %             temptrunk = [];
        end
    end
    
% %     % % 找到树根点之后 通过体素化提取与树根点相连的树干点
% %     cellsize = 0.5;%0.5
% % 
% %     points = grid_p(:,1:3);
% %     
% %     roots(:,9) = 0;
% %     
% %     [num,grid] = Grid_num(points,cellsize);
% %     [~,ia,ib] = intersect(roots(:,1:3),points,'rows');
% %     L = bwlabeln(grid,6);
% %     LL = reshape(L,[],1);
% %     tree_grid2 = cell(j,1);
% %     for j=1:size(ib,1)
% %         L3 = [];
% %         L2(j,1) = L(num(ib(j),1),num(ib(j),2),num(ib(j),3));
% %         L3 = intersect(L2(1:j-1,1),L2(j,1)); % 判断是否有重复探测
% %         if isempty(L3)
% %             trees_idx{j,:} =  find(LL==L(num(ib(j),1),num(ib(j),2),num(ib(j),3)));
% %             sz = [size(grid,1) size(grid,2) size(grid,3)];
% %             [I1,I2,I3] = ind2sub(sz,trees_idx{j,:});
% %             tree_grid2{j,1}(:,1:3) = [I1,I2,I3];
% %             clear I1 I2 I3
% %         end
% %     end
% % 
% %     for j=1:size(tree_grid2,1)
% %         if ~isempty(tree_grid2{j,:})
% %             gridnum = tree_grid2{j,:};
% %             tree_ini1 = [];
% %             for k=1:size(tree_grid2{j,:},1)
% %                 idx3 = [];
% %                 idx3 = find(num(:,1)==gridnum(k,1) & num(:,2)==gridnum(k,2) & num(:,3)==gridnum(k,3));
% %                 tree_ini1 = [tree_ini1;points(idx3,1:3)];
% %             end
% %             trees_ini2{j,1}(:,1:3) = tree_ini1;
% %             trees_ini2{j,1}(:,4) = j;
% %         end
% %     end
% %     trees_final2 = cell2mat(trees_ini2);
% %     
% %     figure(4)
% %     pcshow(trees_final2(:,1:3),trees_final2(:,4),'MarkerSize',2);hold on
% %     pcshow(roots(:,1:3),'red','MarkerSize',20);hold off
% %     
% %     toc
% %     if o~=1
% %         trees_final2(:,4)=trees_final2(:,4)+max(trees_final2(:,4));
% %     end

    root0926=[root0926;roots];
    roots=[];

    trunk0926 = [trunk0926;trunk1];
    trunk1 = [];
    clear tree_grid2 trees_ini2 num grid
% %     trees=[trees;trees_final2];
% %     trees_final2=[];
% 
% %     roots_gv =[roots_gv;temproots_gv];
% %     temproots_gv =[];
end

%%  连通性分析
points = gp;
points(:,4) = 0;

% cellsize2 = 0.5;
[num,grid] = Grid_num(points,cellsize2);
L = bwlabeln(grid,18);% LL = reshape(L,[],1);
[C,~, ic] = unique(num, 'rows');% C表示非空体素（体素内有点云）

grid_idx = accumarray(ic, (1:size(points,1)).', [], @(x){x}); 
for i = 1:size(C,1)
    curvoxel = C(i,:);
    label = L(curvoxel(1),curvoxel(2),curvoxel(3));
    points(grid_idx{i,1},4) = label;
end

label = accumarray(points(:,4), (1:size(points,1)).', [], @(x){x});
[~,id1,~] = intersect(points(:,1:3),root0926(:,1:3),'rows');
% tal = tabulate(LL); %计算每个连通分量内的点数(tal(:,2))
LL = ones(length(label),1);LL(:,2) = 0; %第一列表示是否剔除||第二列代表剔除后的索引
id = {};% 储存保留下来的连通分量内各点的类别
for i = 1:length(label)
    n = length(label{i,1});
    c1 = intersect(label{i,1},id1);
    if  n < 100 && isempty(c1)% 当前连通部分内点过少需要剔除，数量较多的保留
        LL(i,1) = 0; 
    elseif n >= 100 && n < 300 && isempty(c1) 
        % 当前连通部分内点数量适中的则需要进一步判断是否有根节点再进行剔除
        LL(i,1) = 0;
    end
    [rows,~] = find(LL(1:i,1)==0);
    temp = {};
    if LL(i,1)~=0
        LL(i,2) = i-length(rows);
        temp = repmat(LL(i,2),n,1);
    end
    id = [id;temp];
end

rows = cell2mat(label(LL(:,1)==1,1));
finalpoints = [points(rows,1:3) cell2mat(id)];

[c1,ia,ib] = intersect(finalpoints(:,1:3),root0926(:,1:3),'rows');
finalroots = [c1 finalpoints(ia,4)];% 第4列表示根节点所属的连通分量
finaltrunk = trunk0926(ib,:);

figure(5);hold on
pcshow(finalpoints(:,1:3),finalpoints(:,4),'MarkerSize',2);
hold on
pcshow(finalroots(:,1:3),'red','MarkerSize',50)
hold off


%% 删除过小的联通部分
% delete_overh = trees;
% label=unique(delete_overh(:,4));
% j=true(size(delete_overh,1),1);
% for i=1:size(label,1)
%     cur_label=label(i,:);
%     if nnz(delete_overh(:,4)==cur_label)<5000
%         delete_ind=find(delete_overh(:,4)==cur_label);
%         j(delete_ind,:)=false;
%     end
% end
% delete_toosmall=delete_overh(j,:);
% figure
% pcshow(delete_toosmall(:,1:3),delete_toosmall(:,4))
% 
% [~,~,ib] = intersect(delete_toosmall(:,1:3),root0926(:,1:3),'rows');
% roots = root0926(ib,:);
% trunks = trunk0926(ib,:);
% roots_rowvector = roots_gv(ib,:);

end



