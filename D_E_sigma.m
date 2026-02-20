function phase_diagram_3D_Gemstone_Final_4Lines
%% 0. 环境初始化
    tic;
    clear; close all; clc;
    
    if isempty(gcp('nocreate'))
        try
            parpool; 
        catch
            fprintf('无法启动并行池，将使用单核运行。\n');
        end
    end
    
    %% 1. 物理参数
    p.k_off   = 1;    
    p.eta     = 0.05; 
    p.alpha   = 0.01; 
    p.sigma_L = 30;    
    p.K       = 50;   
    
    Cs = 1 / (1 + p.k_off); 
    
    %% 2. 扫描参数设置
    Num_Slices = 100;      % E轴切片数量
    Res_D      = 200;      % D轴精度
    x_max_Sigma = 50;      % Sigma 轴最大值
    
    D_vec = linspace(0.00, 5, Res_D); 
    E_vec = linspace(0.1, 5, Num_Slices); 
    
    [D_grid, E_grid] = meshgrid(D_vec, E_vec);
    
    D_flat = D_grid(:);
    E_flat = E_grid(:);
    total_points = numel(D_flat);
    
    Sig_RH_flat    = zeros(total_points, 1); 
    Sig_Delta_flat = zeros(total_points, 1);
    
    k_val_vec = 0.00:0.05:15; 
    k_vec = k_val_vec * pi;   
    K2_vec = k_vec.^2;  
    
    %% 3. 核心物理计算
    fprintf('正在进行物理场计算 (%d 个网格点)...\n', total_points);
    
    p_k_off = p.k_off; p_eta = p.eta; 
    p_alpha = p.alpha; p_sigma_L = p.sigma_L; p_K = p.K;
    
    parfor i = 1:total_points
        current_D = D_flat(i);
        current_E = E_flat(i);
        
        % === Part A: Hopf ===
        Ak_base_all = 1 + p_k_off + current_D * K2_vec;
        Bk_all      = 1 + p_eta * K2_vec;
        m2_vec = -Cs * K2_vec; 
        c2_vec = Bk_all .* (Ak_base_all + p_K - p_alpha * p_sigma_L) + current_E * K2_vec;
        m1_vec = -Cs * p_K * K2_vec;
        c1_vec = current_E * K2_vec .* (Ak_base_all + p_K - p_alpha * p_sigma_L) + Ak_base_all .* Bk_all * p_K;
        c0_vec = current_E * p_K * K2_vec .* Ak_base_all;
        c3_vec = Bk_all;
        
        quad_a = m1_vec .* m2_vec;
        quad_b = m1_vec .* c2_vec + m2_vec .* c1_vec;
        quad_c = c1_vec .* c2_vec - c0_vec .* c3_vec;
        delta_rh = quad_b.^2 - 4 * quad_a .* quad_c;
        valid_mask = delta_rh >= 0;
        
        if any(valid_mask)
            r1 = (-quad_b(valid_mask) + sqrt(delta_rh(valid_mask))) ./ (2 * quad_a(valid_mask));
            r2 = (-quad_b(valid_mask) - sqrt(delta_rh(valid_mask))) ./ (2 * quad_a(valid_mask));
            pos_r = [r1(r1>0), r2(r2>0)];
            if isempty(pos_r), Sig_RH_flat(i) = 10000; else, Sig_RH_flat(i) = min(pos_r); end
        else
            Sig_RH_flat(i) = 10000;
        end
        
        % === Part B: Delta ===
        min_delta_k = 10000;
        m2_sub = m2_vec; c2_sub = c2_vec; m1_sub = m1_vec; c1_sub = c1_vec; c0_sub = c0_vec; c3_sub = c3_vec;
        for j = 1:length(K2_vec)
            pA2 = [m2_sub(j), c2_sub(j)]; pA1 = [m1_sub(j), c1_sub(j)];
            c0 = c0_sub(j); c3 = c3_sub(j);
            pA2_sq = conv(pA2, pA2); pA2_cb = conv(pA2_sq, pA2); 
            pA1_sq = conv(pA1, pA1); pA1_cb = conv(pA1_sq, pA1); 
            T1 = 18 * c3 * c0 * conv(pA2, pA1); T2 = 4 * c0 * pA2_cb; 
            T3 = conv(pA2_sq, pA1_sq); T4 = 4 * c3 * pA1_cb; T5 = 27 * c3^2 * c0^2; 
            coeffs = zeros(1, 5);
            l1=length(T1); coeffs(6-l1:5) = coeffs(6-l1:5) + T1;
            l2=length(T2); coeffs(6-l2:5) = coeffs(6-l2:5) - T2;
            l3=length(T3); coeffs(6-l3:5) = coeffs(6-l3:5) + T3;
            l4=length(T4); coeffs(6-l4:5) = coeffs(6-l4:5) - T4;
            coeffs(5)     = coeffs(5)     - T5; 
            r_d = roots(coeffs);
            valid_r = r_d(imag(r_d)==0 & real(r_d)>0);
            if ~isempty(valid_r), min_delta_k = min(min_delta_k, max(valid_r)); end
        end
        Sig_Delta_flat(i) = min_delta_k;
    end
    
    Sigma_Hopf_Mat = reshape(Sig_RH_flat, size(D_grid));
    Sigma_Delta_Mat = reshape(Sig_Delta_flat, size(D_grid));
    Sigma_Collapse_Mat = max(Sigma_Hopf_Mat, Sigma_Delta_Mat);
    
    % 数据清洗
    Sigma_Hopf_Mat(Sigma_Hopf_Mat > x_max_Sigma) = NaN;
    Sigma_Collapse_Mat(Sigma_Collapse_Mat > x_max_Sigma) = NaN;

    %% 4. 绘图与渲染
    figure('Color', 'white', 'Position', [100, 50, 1000, 800]);
    hold on; grid on; box on;
    
    Color_Stable   = [0.50, 0.76, 1.0];  
    Color_Osc      = [0.0, 0.6, 0.9];    
    Color_Collapse = [0.1, 0.2, 0.5];    
    
    fprintf('正在绘制层切片...\n');
    
    % --- 4a. 绘制常规切片 ---
    for i = 1:Num_Slices
        current_E = E_vec(i);
        Z_plane = ones(size(D_vec)) * current_E;
        
        Curve_Hopf = Sigma_Hopf_Mat(i, :);
        Curve_Col  = Sigma_Collapse_Mat(i, :);
        Hopf_Clean = Curve_Hopf; Hopf_Clean(isnan(Hopf_Clean)) = x_max_Sigma; 
        Col_Clean = Curve_Col; Col_Clean(isnan(Col_Clean)) = x_max_Sigma;
        
        % 绘制 Body
        X_S = [zeros(size(D_vec)), fliplr(Hopf_Clean)];
        Y_S = [D_vec, fliplr(D_vec)]; Z_S = [Z_plane, fliplr(Z_plane)];
        fill3(X_S, Y_S, Z_S, Color_Stable, 'EdgeColor', 'none', 'FaceAlpha', 0.02, 'Tag', 'Body_Patch'); 
        
        Osc_Mask = Col_Clean > Hopf_Clean;
        if any(Osc_Mask)
            X_OL = Hopf_Clean(Osc_Mask); X_OR = Col_Clean(Osc_Mask);
            Y_O = D_vec(Osc_Mask); Z_O = Z_plane(Osc_Mask);
            X_O_Plot = [X_OL, flipud(X_OR')']; Y_O_Plot = [Y_O, flipud(Y_O')']; Z_O_Plot = [Z_O, flipud(Z_O')'];
            fill3(X_O_Plot, Y_O_Plot, Z_O_Plot, Color_Osc, 'EdgeColor', 'none', 'FaceAlpha', 0.08, 'Tag', 'Body_Patch');
        end
        
        Col_Mask = Col_Clean < x_max_Sigma;
        if any(Col_Mask)
            X_CL = Col_Clean(Col_Mask); Y_C = D_vec(Col_Mask); Z_C = Z_plane(Col_Mask);
            X_CR = ones(size(X_CL)) * x_max_Sigma;
            X_C_Plot = [X_CL, flipud(X_CR')']; Y_C_Plot = [Y_C, flipud(Y_C')']; Z_C_Plot = [Z_C, flipud(Z_C')'];
            fill3(X_C_Plot, Y_C_Plot, Z_C_Plot, Color_Collapse, 'EdgeColor', 'none', 'FaceAlpha', 0.04, 'Tag', 'Body_Patch');
        end
        
        % 绘制切片线条 (首尾 + 中间)
        if i == 1 || i == Num_Slices
            lineTag = 'Start_End_Line'; lw = 2.0; 
        elseif mod(i, 20) == 0
            lineTag = 'Intermediate_Line'; lw = 0.5;
        else
            lineTag = 'None';
        end
        if ~strcmp(lineTag, 'None')
            plot3(Hopf_Clean, D_vec, Z_plane, 'Color', [0.0, 0.4, 0.8], 'LineWidth', lw, 'Tag', lineTag);
            plot3(Col_Clean, D_vec, Z_plane, 'Color', [0.0, 0.1, 0.4], 'LineWidth', lw, 'Tag', lineTag);
        end
    end
    
    % --- 4b. 【新增】绘制 4 条纵向脊线 (Corner Lines) ---
    fprintf('正在绘制 4 条纵向脊线...\n');
    % 提取数据：第1列对应 D=0 (min)，最后1列对应 D=5 (max)
    
    % 线条 1: Hopf @ D=0
    Hopf_D0 = Sigma_Hopf_Mat(:, 1);
    Hopf_D0(isnan(Hopf_D0)) = x_max_Sigma; % 同样做清洗以连接到边界
    
    % 线条 2: Hopf @ D=5
    Hopf_Dmax = Sigma_Hopf_Mat(:, end);
    Hopf_Dmax(isnan(Hopf_Dmax)) = x_max_Sigma;
    
    % 线条 3: Collapse @ D=0
    Col_D0 = Sigma_Collapse_Mat(:, 1);
    Col_D0(isnan(Col_D0)) = x_max_Sigma;
    
    % 线条 4: Collapse @ D=5
    Col_Dmax = Sigma_Collapse_Mat(:, end);
    Col_Dmax(isnan(Col_Dmax)) = x_max_Sigma;
    
    % Z 轴数据
    Z_Line = E_vec(:);
    
    % 统一绘制设置
    CornerColor = [0.8, 0.2, 0.2]; % 使用醒目的红褐色，区别于蓝色的切片线
    CornerLW    = 2.5;             % 加粗
    CornerTag   = 'Corner_Line';   % 专用 Tag
    
    % D=0 的两条线 (Y坐标全为 0)
    plot3(Hopf_D0, zeros(size(Z_Line)), Z_Line, 'Color', CornerColor, 'LineWidth', CornerLW, 'Tag', CornerTag);
    plot3(Col_D0,  zeros(size(Z_Line)), Z_Line, 'Color', CornerColor, 'LineWidth', CornerLW, 'Tag', CornerTag);
    
    % D=5 的两条线 (Y坐标全为 5)
    plot3(Hopf_Dmax, ones(size(Z_Line))*5, Z_Line, 'Color', CornerColor, 'LineWidth', CornerLW, 'Tag', CornerTag);
    plot3(Col_Dmax,  ones(size(Z_Line))*5, Z_Line, 'Color', CornerColor, 'LineWidth', CornerLW, 'Tag', CornerTag);
    
    
    %% 5. 坐标轴与外框设置
    box_x = [0 x_max_Sigma]; box_y = [0 5]; box_z = [min(E_vec) max(E_vec)];
    BoxColor = [0.2, 0.2, 0.2]; BoxWidth = 1.0;
    args = {'Color', BoxColor, 'LineWidth', BoxWidth, 'Tag', 'Frame_Line'};
    
    % 绘制盒子
    line([box_x(1) box_x(1)], [box_y(1) box_y(1)], box_z, args{:});
    line([box_x(2) box_x(2)], [box_y(1) box_y(1)], box_z, args{:});
    line([box_x(2) box_x(2)], [box_y(2) box_y(2)], box_z, args{:});
    line([box_x(1) box_x(1)], [box_y(2) box_y(2)], box_z, args{:});
    line([box_x(1) box_x(2)], [box_y(1) box_y(1)], [box_z(1) box_z(1)], args{:});
    line([box_x(2) box_x(2)], [box_y(1) box_y(2)], [box_z(1) box_z(1)], args{:});
    line([box_x(2) box_x(1)], [box_y(2) box_y(2)], [box_z(1) box_z(1)], args{:});
    line([box_x(1) box_x(1)], [box_y(2) box_y(1)], [box_z(1) box_z(1)], args{:});
    line([box_x(1) box_x(2)], [box_y(1) box_y(1)], [box_z(2) box_z(2)], args{:});
    line([box_x(2) box_x(2)], [box_y(1) box_y(2)], [box_z(2) box_z(2)], args{:});
    line([box_x(2) box_x(1)], [box_y(2) box_y(2)], [box_z(2) box_z(2)], args{:});
    line([box_x(1) box_x(1)], [box_y(2) box_y(1)], [box_z(2) box_z(2)], args{:});
    
    %% 6. 视觉美化
    view([-40, 25]); axis tight;
    xlim([0 x_max_Sigma]); ylim([0 5]); zlim([min(E_vec) max(E_vec)]);
    xlabel('Active Stress \sigma_0', 'FontSize', 14, 'Interpreter', 'tex', 'FontWeight', 'bold');
    ylabel('Diffusion D', 'FontSize', 14, 'Interpreter', 'tex', 'FontWeight', 'bold');
    zlabel('Elastic Modulus E', 'FontSize', 14, 'Interpreter', 'tex', 'FontWeight', 'bold');
    set(gca, 'XDir', 'normal');
    material shiny; camlight('headlight'); camlight(180, 45); lighting phong;         
    
    %% 7. 终极分层导出方案 (Layer 1, 2, 3, 4)
    fprintf('正在执行分层导出 (EPS/PNG)...\n');
    
    % --- 准备对象句柄 ---
    h_BodyPatch     = findobj(gca, 'Tag', 'Body_Patch');
    h_FrameLine     = findobj(gca, 'Tag', 'Frame_Line');
    h_StartEndLines = findobj(gca, 'Tag', 'Start_End_Line');
    h_InterLines    = findobj(gca, 'Tag', 'Intermediate_Line');
    h_CornerLines   = findobj(gca, 'Tag', 'Corner_Line'); % 获取新的脊线句柄
    
    h_AllBodyLines = [h_StartEndLines; h_InterLines; h_CornerLines]; % 包含所有实体相关的线
    
    % 全局设置
    set(gcf, 'Color', 'none'); set(gca, 'Color', 'none'); set(gcf, 'InvertHardcopy', 'off');
    
    % === Layer 1: 3D 实体 (PNG) ===
    % 显示: Body, 所有线条(增加质感)
    fprintf('  - 导出 Layer 1: 3D 实体 (PNG)...\n');
    axis off; grid off;
    set(h_FrameLine, 'Visible', 'off'); 
    set(h_BodyPatch, 'Visible', 'on');
    set(h_AllBodyLines, 'Visible', 'on'); 
    set(gcf, 'Renderer', 'opengl');
    print(gcf, 'Layer1_Crystal_Body.png', '-dpng', '-r600', '-opengl');
    
    % === Layer 2: 坐标轴与边框 (EPS) ===
    % 显示: Frame, 坐标轴
    fprintf('  - 导出 Layer 2: 坐标轴与边框 (EPS)...\n');
    axis on; box on; set(gca, 'GridAlpha', 0.15); 
    set(h_BodyPatch, 'Visible', 'off');
    set(h_AllBodyLines, 'Visible', 'off'); 
    set(h_FrameLine, 'Visible', 'on');
    set(gcf, 'Renderer', 'painters'); 
    print(gcf, 'Layer2_Axis_Frame.eps', '-depsc', '-r600', '-painters');
    
    % === Layer 3: 切片层面线条 (EPS) ===
    % 显示: Start_End Lines
    fprintf('  - 导出 Layer 3: 切片层面线条 (EPS)...\n');
    axis off; grid off; box off;
    set(h_FrameLine, 'Visible', 'off');
    set(h_BodyPatch, 'Visible', 'off');
    set(h_CornerLines, 'Visible', 'off'); % 隐藏纵向脊线
    set(h_InterLines, 'Visible', 'off');
    set(h_StartEndLines, 'Visible', 'on'); % 只显示首尾切片线
    set(gcf, 'Renderer', 'painters');
    print(gcf, 'Layer3_Slice_Lines.eps', '-depsc', '-r600', '-painters');
    
    % === 【新增】Layer 4: 纵向脊线 (EPS) ===
    % 显示: Corner Lines
    fprintf('  - 导出 Layer 4: 纵向脊线 (EPS)...\n');
    axis off; grid off; box off;
    % 隐藏其他所有
    set(h_StartEndLines, 'Visible', 'off');
    
    % 只显示脊线
    set(h_CornerLines, 'Visible', 'on');
    
    set(gcf, 'Renderer', 'painters');
    print(gcf, 'Layer4_Corner_Lines.eps', '-depsc', '-r600', '-painters');
    
    % --- 恢复现场 ---
    set(h_BodyPatch, 'Visible', 'on');
    set(h_AllBodyLines, 'Visible', 'on');
    set(h_FrameLine, 'Visible', 'on');
    axis on; box on;
    set(gcf, 'Color', 'white'); 
    set(gcf, 'Renderer', 'opengl'); 
    
    toc;
    fprintf('导出完成！共 4 层文件。\n');
end